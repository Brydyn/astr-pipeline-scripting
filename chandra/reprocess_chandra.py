#!/usr/bin/env python3
"""
CIAO Chandra Reprocessing Script

Automates reprocessing of multiple OBSIDs with optional flare filtering.

Steps:
1. Download data with `download_chandra_obsid`.
2. Run `chandra_repro` for each OBSID.
3. Apply energy filtering (0.3–10 keV).
4. Prompt user to generate DS9 region file for flare screening.
5. Extract and plot lightcurve.
6. Prompt user for flare filtering threshold.
7. Optionally apply GTI filtering with `dmgti`.
8. Reproject all observations into `reproj/`.
9. Run `flux_obs` to create exposure-corrected images.

Usage:
    python3 reprocess_chandra.py 1625 5796 6257

Requirements:
- CIAO environment must be active
- `pycrates` and `matplotlib` must be installed in Python environment
"""

import sys
import subprocess
import os
from pathlib import Path
import matplotlib.pyplot as plt
from pycrates import read_file

def run_cmd(cmd, cwd=None):
    print(f"\nRunning: {cmd}")
    subprocess.run(cmd, shell=True, check=True, cwd=cwd)

def process_obsid(obsid):
    print(f"\n=== Processing OBSID {obsid} ===")
    repro_dir = Path(f"{obsid}/repro")
    cwd = Path.cwd()
    os.chdir(repro_dir)

    # Set ARDLIB to current observation
    bpix = list(Path('.').glob("*_repro_bpix1.fits"))[0].name
    run_cmd(f"punlearn ardlib && acis_set_ardlib {bpix}")

    # Energy filtering
    evt2 = list(Path('.').glob("*_evt2.fits"))[0].name
    base = evt2.replace("_evt2.fits", "")
    filtered_evt = f"{base}_0.3-10.fits"
    run_cmd(f"punlearn dmcopy && dmcopy \"{evt2}[energy=300:10000]\" {filtered_evt}")

    # Prompt for user-generated region file
    print(f"\nOpen DS9, load {filtered_evt} (bin 4), draw full-region circle,")
    print(f"then save as {obsid}_ltcrv.reg in CIAO (physical) format.")
    input("Press Enter once the region file is saved...")

    # Extract lightcurve
    run_cmd("punlearn dmextract")
    run_cmd(f"pset dmextract infile=\"{filtered_evt}[exclude sky=region({obsid}_ltcrv.reg)][bin time=::3.24104]\"")
    run_cmd("pset dmextract outfile=background_lc.fits")
    run_cmd("pset dmextract opt=ltc1")
    run_cmd("dmextract")

    # Plot lightcurve
    tab = read_file("background_lc.fits")
    time = tab.get_column("time").values
    rate = tab.get_column("count_rate").values
    plt.plot(time, rate, marker="o", linestyle="")
    plt.title(f"Lightcurve for OBSID {obsid}")
    plt.xlabel("Time")
    plt.ylabel("Count Rate")
    plt.show()

    # Ask user whether to apply filtering
    userlimit = input("Enter count_rate threshold for GTI filtering (or press Enter to skip): ").strip()

    # Save original evt2
    run_cmd(f"mv {evt2} {evt2.replace('.fits', '_unfilt.fits')}")
    
    if userlimit:
        try:
            val = float(userlimit)
            run_cmd("punlearn dmgti")
            run_cmd("pset dmgti infile=background_lc.fits")
            run_cmd("pset dmgti outfile=bkg_gti.fits")
            run_cmd(f"pset dmgti userlimit=\"count_rate<={val}\"")
            run_cmd("dmgti")
            run_cmd(f"dmcopy \"{filtered_evt}[@bkg_gti.fits]\" {evt2}")
        except ValueError:
            print("Invalid input. Skipping GTI filtering.")
            run_cmd(f"dmcopy \"{filtered_evt}\" {evt2}")

    os.chdir(cwd)

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 reprocess_chandra.py OBSID1 OBSID2 ...")
        sys.exit(1)

    obsids = [obsid.strip() for obsid in sys.argv[1:]]

    run_cmd(f"download_chandra_obsid {','.join(obsids)}")
    run_cmd(f"punlearn chandra_repro && chandra_repro indir={','.join(obsids)} mode=h")

    for obsid in obsids:
        process_obsid(obsid)
        
    # Reprojection
    if len(obsids) > 1:
        run_cmd("punlearn reproject_obs")
        run_cmd(f"reproject_obs {','.join(obsids)} reproj/")
        flux_input = "reproj/*_reproj_evt.fits"
    else:
        print("Only one OBSID provided — skipping reproject_obs.")
        flux_input = f"{obsids[0]}/repro/*evt2.fits"  # or evt2_unfilt.fits, depending on use case
    
    # Run flux_obs regardless
    run_cmd("punlearn flux_obs")
    run_cmd(f"flux_obs \"{flux_input}\" merged/ "
            "bands=\"broad,0.5:1.0:0.75,1.0:2.0:1.5,2.0:10.0:6.0\" binsize=1")

    run_cmd("punlearn reproject_obs")
    run_cmd(f"reproject_obs {','.join(obsids)} reproj/")

    # Flux-corrected images
    run_cmd("punlearn flux_obs")
    run_cmd("flux_obs \"reproj/*_reproj_evt.fits\" merged/ "
            "bands=\"broad,0.5:1.0:0.75,1.0:2.0:1.5,2.0:10.0:6.0\" binsize=1")

    print("\nAll processing complete. Outputs are in reproj/ and merged/.")
    print("\nSuggested next steps: Chandra_U2+U30_All_Images.py and Chandra_multi-src_specextract+group.py")

if __name__ == "__main__":
    main()
