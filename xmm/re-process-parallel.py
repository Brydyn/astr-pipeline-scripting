"""
Welcome to the generalized script for XMM-Newton ODF preprocessing.
Author: Brydyn Mac Intyre
Date: Updated [July 2025]

Purpose:
Automate ODF retrieval and reduction using SAS. Handles multiple exposures per instrument (e.g., scheduled “S” and unscheduled “U”)
and keeps all resulting files for downstream selection. Final filtered and flare-cleaned (espfilt) files are output as:
  {OBS}_{CAMERA}_{SUFFIX}_filt_esp.fits
  
THIS VERSION RUNS IN PARALLEL! Make sure you know how many max_workers you can run before starting!

This version uses concurrent.futures.ProcessPoolExecutor to reduce multiple
observations concurrently (default max_workers = 8).

Each OBSID is handled independently with:
- Isolated working directories
- Separate environment variables
- Per-OBSID logs (e.g., 0123456789.log)

Assumptions:
- Internet access for `curl` to fetch ODF tarballs
- SAS setup in environment
- astropy for FITS ONTIME and header parsing
- Final user selection occurs before running specextract_v27.py

Outputs:
- Full reduction per camera/exposure
- All event files retained and uniquely named by exposure suffix
- Summary written to `xmm_reduction.log`

Usage:
    python3 re-process-parallel.py 0123456789 0123456790 ...
"""

from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import subprocess
import sys
from pathlib import Path
from astropy.io import fits

def run_cmd(cmd, env=None, log_file=None):
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, env=env)
    if log_file:
        with open(log_file, "a") as log:
            log.write(f"\n$ {cmd}\n")
            log.write(result.stdout)
            log.write(result.stderr)
    result.check_returncode()

def process_single_obsid(obsid, base_dir_str):
    base_dir = Path(base_dir_str)
    log_path = base_dir / f"{obsid}.log"
    spectra_dir = base_dir / "Spectra"
    spectra_dir.mkdir(exist_ok=True)

    obs_dir = base_dir / obsid
    odf_dir = obs_dir / "ODF"
    proc_dir = obs_dir / "PROC"
    odf_dir.mkdir(parents=True, exist_ok=True)
    proc_dir.mkdir(parents=True, exist_ok=True)

    os.chdir(odf_dir)
    tar_url = f"https://nxsa.esac.esa.int/nxsa-sl/servlet/data-action?OBSERVATION_ID={obsid}&LEVEL=ODF&PROTOCOL=HTTP"
    run_cmd(f"curl -f -L -o ODF.tar \"{tar_url}\"", log_file=log_path)
    odf_tar = odf_dir / "ODF.tar"
    if not odf_tar.exists() or odf_tar.stat().st_size < 1024:
        raise RuntimeError(f"{obsid}: ODF.tar was invalid (empty or corrupted).")
    run_cmd("tar xvf ODF.tar", log_file=log_path)
    odf_tar.unlink()

    inner_tar = next(odf_dir.glob("*.TAR"), None)
    if not inner_tar:
        raise RuntimeError("No inner .TAR file found after initial extraction.")
    run_cmd(f"tar xvf {inner_tar.name}", log_file=log_path)
    inner_tar.unlink()

    env = os.environ.copy()
    env["SAS_ODF"] = str(odf_dir)
    os.chdir(odf_dir)
    run_cmd("cifbuild", env=env, log_file=log_path)
    env["SAS_CCF"] = str(odf_dir / "ccf.cif")
    run_cmd("odfingest", env=env, log_file=log_path)

    sum_file = next(odf_dir.glob("*SCX00000SUM.SAS"), None)
    if not sum_file:
        raise RuntimeError("SCX00000SUM.SAS file not found.")
    env["SAS_ODF"] = str(sum_file)

    os.chdir(proc_dir)
    run_cmd("emproc", env=env, log_file=log_path)
    run_cmd("epproc", env=env, log_file=log_path)

    cam_types = {
        'mos1': '*EMOS1*_ImagingEvts.ds',
        'mos2': '*EMOS2*_ImagingEvts.ds',
        'pn': '*EPN*_ImagingEvts.ds'
    }

    exposure_files = []

    for cam, pattern in cam_types.items():
        candidates = list(proc_dir.glob(pattern))
        if not candidates:
            continue

        for f in candidates:
            with fits.open(f) as hdul:
                hdr = hdul[0].header
                ontime = hdul[1].header.get('ONTIME', 0.0)
                expidstr = hdr.get('EXPIDSTR', 'XXX')
            suffix = expidstr
            newname = f"{obsid}_{cam}_{suffix}.fits"
            run_cmd(f"cp {f} {proc_dir / newname}", log_file=log_path)
            exposure_files.append((cam, suffix, proc_dir / newname, odf_dir / "ccf.cif"))

    os.chdir(proc_dir)

    for cam, suffix, infile, _ in exposure_files:
        expr = {
            'mos1': f"(PATTERN <= 12)&&(PI in [200:12000])&&#XMMEA_EM",
            'mos2': f"(PATTERN <= 12)&&(PI in [200:12000])&&#XMMEA_EM",
            'pn': f"(PATTERN <= 4)&&(PI in [200:15000])&&(FLAG==0)&&#XMMEA_EP"
        }[cam]
        outfile = infile.stem + "_filt.fits"
        run_cmd(f"evselect table={infile} withfilteredset=yes "
                f"filteredset={outfile} filtertype=expression "
                f"expression='{expr}' keepfilteroutput=yes updateexposure=yes filterexposure=yes", log_file=log_path)

    for cam, suffix, infile, _ in exposure_files:
        filt_file = Path(f"{infile.stem}_filt.fits")
        run_cmd(f"espfilt eventfile={filt_file}", log_file=log_path)
        matched = list(Path('.').glob(f"{cam}{suffix}*-allevc.fits"))
        if matched:
            final = f"{obsid}_{cam}_{suffix}_filt_esp.fits"
            run_cmd(f"mv {matched[0]} {final}", log_file=log_path)
        for f in Path('.').glob(f"{cam}{suffix}*"):
            if "-allevc.fits" not in str(f):
                f.unlink()

    spectra_dir = base_dir / "Spectra"
    written = set()
    with open(spectra_dir / "dir.txt", "a") as f1, open(spectra_dir / "cifdir.txt", "a") as f2:
        for _, suffix, _, cif in exposure_files:
            for f in proc_dir.glob(f"{obsid}_*_{suffix}_filt_esp.fits"):
                resolved = str(f.resolve())
                if resolved not in written:
                    f1.write(resolved + "\n")
                    f2.write(str(cif.resolve()) + "\n")
                    written.add(resolved)

    return f"{obsid} complete"

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 re-process-parallel.py OBSID1 OBSID2 ...")
        sys.exit(1)

    base_dir = Path.cwd()
    spectra_dir = base_dir / "Spectra"
    spectra_dir.mkdir(exist_ok=True)

    obsids = sys.argv[1:]
    max_workers = min(8, len(obsids))  # scale up to 8 jobs max

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(process_single_obsid, obsid, str(base_dir)): obsid for obsid in obsids}
        for future in as_completed(futures):
            try:
                print(future.result())
            except Exception as e:
                print(f"Error processing {futures[future]}: {e}")

    os.chdir(spectra_dir)
    print("\nAll outputs complete. Please:")
    print("- Copy the most recent specextract_v27.py into this 'Spectra/' directory.")
    print("- Complete your *_Regions.txt or .reg files for each OBSID.")
    print("- Then run specextract_v27.py to continue with spectral extraction.")