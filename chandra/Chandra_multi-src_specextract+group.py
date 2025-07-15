"""
Updated July 2025

@author: Brydyn Mac Intyre

Combined Chandra Spectral Extraction and Grouping Script

Performs:
1. Spectral extraction (with CIAO `specextract`) for multiple OBSIDs and multiple source regions
2. Spectral grouping (with HEASoft `ftgrouppha` + `grppha`) using a common background

USAGE:
    Run this script from any directory:
    $ python3 Chandra_multi-src_specextract+group.py

REQUIREMENTS:
- CIAO must be initialized in the shell running this script
- HEASoft tools (`ftgrouppha`, `grppha`) must be available in PATH
- Regions must be named `reg1.reg`, `reg2.reg`, ..., in CIAO physical format
- Region files must contain one inclusion region (optionally with "-" exclusions)
- Any necessary filtering already completed
- Event files for each observation must be reprojected and match region coordinate system
- Input and output paths must be edited in the config section

EXAMPLE FILES:
- Input event:           reproj/1234_reproj_evt.fits
- Region file:           Regions/ciao_indiv/reg1.reg
- Output spectrum:       Spectra/1234_src1.pi
- Grouped spectrum:      Spectra/Grouped/1234_src1_grp3sn.fits

BEHAVIOUR:
- Automatically detects how many regions exist in the specified region directory
- Uses region number `background_region` (e.g., 9) as the background file for all regions

CONFIGURATION:
- Set `obsids`, `evt_dir`, `region_dir`, `output_dir`, and `background_region`
- Grouping S/N ratio is controlled by `snr`
"""

# =================== USER CONFIGURATION ===================
obsids = ["1625", "5796", "6257"]                # OBSIDs to process
evt_dir = "/home/Name/Object/Chandra/reproj"     # Directory with *_reproj_evt.fits
region_dir = "/home/Name/Object/Regions/ciao_indiv"  # Directory with reg{i}.reg files
output_dir = "/home/Name/Object/Spectra"         # Where output spectra go

background_region = 9                            # Region number to use for background
snr = 3.0                                         # S/N ratio for grouping
# ==========================================================

import os
import subprocess
import re
from glob import glob

# Detect number of regions by listing files
region_files = sorted([
    f for f in os.listdir(region_dir)
    if re.fullmatch(r"reg\d+\.reg", f)
], key=lambda x: int(re.findall(r"\d+", x)[0]))
n_regions = len(region_files)
if n_regions == 0:
    raise RuntimeError(f"No region files found in {region_dir}")

grouped_dir = os.path.join(output_dir, "Grouped")
os.makedirs(grouped_dir, exist_ok=True)

# --------- SPECEEXTRACT: Extract ungrouped spectra ---------
print("\n=== STEP 1: Extracting Spectra ===")

for obsid in obsids:
    evt_file = os.path.join(evt_dir, f"{obsid}_reproj_evt.fits")
    if not os.path.exists(evt_file):
        fallback = glob(os.path.join(evt_dir, f"acisf{obsid}_repro_evt2.fits"))
        if len(fallback) == 1:
            evt_file = fallback[0]
            print(f"[{obsid}] No reproj file — using {evt_file}")
        elif len(fallback) == 0:
            raise FileNotFoundError(f"[{obsid}] No event file found in {evt_dir}")
        else:
            raise RuntimeError(f"[{obsid}] Multiple evt2 files found in {evt_dir}; please disambiguate.")
            
    for i in range(1, n_regions + 1):
        reg_file = os.path.join(region_dir, f"reg{i}.reg")
        infile = f"{evt_file}[sky=region({reg_file})]"
        outroot = os.path.join(output_dir, f"{obsid}_src{i}")

        print(f"[{obsid}] Region {i}/{n_regions} → {outroot}.pi")
        subprocess.run(["punlearn", "specextract"], check=True)
        subprocess.run([
            "specextract",
            f"infile={infile}",
            f"outroot={outroot}",
            "grouptype=NONE",
            "mode=h"
        ], check=True)

# --------- GROUPPHA: Group spectra by S/N ---------
print("\n=== STEP 2: Grouping Spectra ===")

snr_str = str(int(round(snr)))
for obsid in obsids:
    bkg_file = os.path.join(output_dir, f"{obsid}_src{background_region}.pi")
    for i in range(1, n_regions + 1):
        pi_file = os.path.join(output_dir, f"{obsid}_src{i}.pi")
        arf_file = os.path.join(output_dir, f"{obsid}_src{i}.arf")
        rmf_file = os.path.join(output_dir, f"{obsid}_src{i}.rmf")
        temp_file = f"{obsid}_src{i}_grp{snr_str}sn_temp.fits"
        grouped_file = os.path.join(grouped_dir, f"{obsid}_src{i}_grp{snr_str}sn.fits")

        print(f"[{obsid}] Grouping region {i} → {grouped_file}")
        subprocess.run([
            "ftgrouppha",
            f"infile={pi_file}",
            f"backfile={bkg_file}",
            f"respfile={rmf_file}",
            f"outfile={temp_file}",
            "grouptype=snmin",
            f"groupscale={snr}"
        ], check=True)

        comm = f"chkey ANCRFILE {arf_file} & chkey RESPFILE {rmf_file} & exit"
        subprocess.run([
            "grppha",
            f"infile={temp_file}",
            f"outfile={grouped_file}",
            f"comm={comm}"
        ], check=True)

        os.remove(temp_file)

print("\nAll spectra extracted and grouped.")
print(f"Grouped output is in: {grouped_dir}")
