"""
Welcome to the generalized script for XMM-Newton ODF preprocessing.
Author: Brydyn Mac Intyre
Date: Updated [July 2025]

Purpose:
Automate ODF retrieval and reduction using SAS. Handles multiple exposures per instrument (e.g., scheduled “S” and unscheduled “U”)
and keeps all resulting files for downstream selection. Final filtered and flare-cleaned (espfilt) files are output as:
  {OBS}_{CAMERA}_{SUFFIX}_filt_esp.fits

Assumptions:
- Internet access for `curl` to fetch ODF tarballs
- SAS setup in environment
- astropy for FITS ONTIME and header parsing
- Final user selection occurs before running specextract_v27.py

Outputs:
- Full reduction per camera/exposure
- All event files retained and uniquely named by exposure suffix
- Summary written to `xmm_reduction.log`

Key steps per observation:
- Create ODF/ and PROC/ directories
- Download ODF data tarball from the NXSA
- Extract both outer and inner tar archives
- Set SAS_ODF and run cifbuild + odfingest
- Re-set SAS_ODF using the *SCX00000SUM.SAS file
- Run emproc and epproc to create event files
- Rename event files based on camera and suffix
- Filter with evselect using standard expressions
- Run espfilt and retain only the *_allevc.fits outputs
- Rename final files with _filt_esp.fits suffix
- Record paths in dir.txt and cifdir.txt
- Prompt user to proceed with specextract_v27.py

Run with:
    python3 re-process.py 0123456789 0123456790 ...
"""

import os
import subprocess
import sys
from pathlib import Path
from astropy.io import fits

def run_cmd(cmd, env=None):
    print(f"\n$ {cmd}")
    subprocess.run(cmd, shell=True, check=True, env=env)

def download_and_extract_odf(obsid, base_dir):
    obs_dir = base_dir / obsid
    odf_dir = obs_dir / "ODF"
    proc_dir = obs_dir / "PROC"
    odf_dir.mkdir(parents=True, exist_ok=True)
    proc_dir.mkdir(parents=True, exist_ok=True)

    os.chdir(odf_dir)
    tar_url = f"https://nxsa.esac.esa.int/nxsa-sl/servlet/data-action?OBSERVATION_ID={obsid}&LEVEL=ODF&PROTOCOL=HTTP"
    run_cmd(f"curl -f -L -o ODF.tar \"{tar_url}\"")
    odf_tar = odf_dir / "ODF.tar"
    if not odf_tar.exists() or odf_tar.stat().st_size < 1024:
        raise RuntimeError(f"{obsid}: ODF.tar was invalid (empty or corrupted).")
    run_cmd("tar xvf ODF.tar")
    odf_tar.unlink()

    inner_tar = next(odf_dir.glob("*.TAR"), None)
    if not inner_tar:
        raise RuntimeError("No inner .TAR file found after initial extraction.")
    run_cmd(f"tar xvf {inner_tar.name}")
    inner_tar.unlink()

    return odf_dir, proc_dir

def configure_sas(obsid, odf_dir):
    env = os.environ.copy()
    env["SAS_ODF"] = str(odf_dir)
    os.chdir(odf_dir)

    run_cmd("cifbuild", env=env)
    env["SAS_CCF"] = str(odf_dir / "ccf.cif")

    run_cmd("odfingest", env=env)

    sum_file = next(odf_dir.glob("*SCX00000SUM.SAS"), None)
    if not sum_file:
        raise RuntimeError("SCX00000SUM.SAS file not found.")
    env["SAS_ODF"] = str(sum_file)

    return env

def run_sas_pipelines(env, proc_dir):
    os.chdir(proc_dir)
    run_cmd("emproc", env=env)
    run_cmd("epproc", env=env)

def extract_and_rename(obsid, odf_dir, proc_dir, log_lines):
    os.chdir(proc_dir)
    cam_types = {
        'mos1': '*EMOS1*_ImagingEvts.ds',
        'mos2': '*EMOS2*_ImagingEvts.ds',
        'pn': '*EPN*_ImagingEvts.ds'
    }

    exposure_durations = {}
    exposure_files = []

    for cam, pattern in cam_types.items():
        candidates = list(proc_dir.glob(pattern))
        if not candidates:
            continue

        exposure_durations[cam] = []
        for f in candidates:
            with fits.open(f) as hdul:
                hdr = hdul[0].header
                ontime = hdul[1].header.get('ONTIME', 0.0)
                expidstr = hdr.get('EXPIDSTR', 'XXX')
            suffix = expidstr  # e.g., 'S001' or 'U002'
            newname = f"{obsid}_{cam}_{suffix}.fits"
            run_cmd(f"cp {f} {newname}")
            exposure_durations[cam].append((newname, ontime))
            exposure_files.append((cam, suffix, Path(newname), odf_dir / "ccf.cif"))

    for cam, entries in exposure_durations.items():
        if len(entries) > 1:
            entries.sort(key=lambda x: -x[1])
            log_lines.append(
                f"OBS {obsid} - {cam.upper()} multiple exposures:\n  " +
                "\n  ".join([f"{Path(f).name} = {t:.1f} s" for f, t in entries]) +
                f"\n  => Longest: {Path(entries[0][0]).name}"
            )

    return exposure_files

def filter_events(obsid, exposure_files):
    for cam, suffix, infile, _ in exposure_files:
        expr = {
            'mos1': f"(PATTERN <= 12)&&(PI in [200:12000])&&#XMMEA_EM",
            'mos2': f"(PATTERN <= 12)&&(PI in [200:12000])&&#XMMEA_EM",
            'pn': f"(PATTERN <= 4)&&(PI in [200:15000])&&(FLAG==0)&&#XMMEA_EP"
        }[cam]
        outfile = infile.stem + "_filt.fits"
        run_cmd(f"evselect table={infile} withfilteredset=yes "
                f"filteredset={outfile} filtertype=expression "
                f"expression='{expr}' keepfilteroutput=yes updateexposure=yes filterexposure=yes")

def run_espfilt_and_cleanup(obsid, exposure_files):
    for cam, suffix, infile, _ in exposure_files:
        filt_file = Path(f"{infile.stem}_filt.fits")
        run_cmd(f"espfilt eventfile={filt_file}")
        matched = list(Path('.').glob(f"{cam}{suffix}*-allevc.fits"))
        if matched:
            final = f"{obsid}_{cam}_{suffix}_filt_esp.fits"
            run_cmd(f"mv {matched[0]} {final}")
        for f in Path('.').glob(f"{cam}{suffix}*"):
            if "-allevc.fits" not in str(f):
                f.unlink()

def process_obsid(obsid, base_dir, spectra_dir, dir_entries, log_lines):
    print(f"\n=== Processing OBSID {obsid} ===")

    odf_dir, proc_dir = download_and_extract_odf(obsid, base_dir)
    env = configure_sas(obsid, odf_dir)
    run_sas_pipelines(env, proc_dir)
    exposure_files = extract_and_rename(obsid, odf_dir, proc_dir, log_lines)
    os.chdir(proc_dir)
    filter_events(obsid, exposure_files)
    run_espfilt_and_cleanup(obsid, exposure_files)

    for _, _, _, cif_path in exposure_files:
        espfilt_files = sorted(Path('.').glob(f"{obsid}_*_*_filt_esp.fits"))
        for f in espfilt_files:
            dir_entries.append((str(f.resolve()), str(cif_path.resolve())))

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 re-process.py OBSID1 OBSID2 ...")
        sys.exit(1)

    base_dir = Path.cwd()
    spectra_dir = base_dir / "Spectra"
    spectra_dir.mkdir(exist_ok=True)

    log_lines = []
    dir_entries = []

    for obsid in sys.argv[1:]:
        try:
            process_obsid(obsid, base_dir, spectra_dir, dir_entries, log_lines)
        except Exception as e:
            log_lines.append(f"ERROR processing {obsid}: {e}")

    # Write log and dir.txt/cifdir.txt
    os.chdir(spectra_dir)
    if log_lines:
        with open("xmm_reduction.log", "w") as f:
            for line in log_lines:
                f.write(line + "\n")
        print("Log written to xmm_reduction.log")

    with open("dir.txt", "w") as f1, open("cifdir.txt", "w") as f2:
        for d, c in dir_entries:
            f1.write(f"{d}\n")
            f2.write(f"{c}\n")

    print("\nAll outputs complete. Please:")
    print("- Copy the most recent specextract_v27.py into this 'Spectra/' directory.")
    print("- Complete your *_Regions.txt or .reg files for each OBSID.")
    print("- Then run specextract_v27.py to continue with spectral extraction.")