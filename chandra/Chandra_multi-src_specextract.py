# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 2025

@author: Brydyn Mac Intyre
"""

"""
Purpose:
    Extract spectra for a list of ObsIDs and multiple regions using CIAO's `specextract`.
    Each region is processed individually for each ObsID.

Requirements:
    - CIAO installed and configured (e.g. `ciao` initialized in the terminal)
    - Region files in ciao format
    - Region files must contain one inclusion region (optionally with "-" exclusions)
    - Event files for each observation must be reprojected and match region coordinate system
    - Any necessary filtering already completed

Inputs:
    - Event files: {evt_dir}/{obsid}_reproj_evt.fits
    - Region files: {region_dir}/reg{i}.reg (for i in 1 to n_regions)

Outputs:
    - Spectra written to: {output_dir}/{obsid}_src{i}.pi (plus .arf and .rmf files)

"""

# ==== FILL IN THESE VALUES ====
obsids = ["1625", "5796", "6257"]         # List of ObsIDs to process
n_regions = 27                            # Total number of regions (1-based index)
evt_dir = "/home/Brydyn/HerA/Chandra/reproj"
region_dir = "/home/Brydyn/HerA/Regions/ciao_indiv"
output_dir = "/home/Brydyn/HerA/Spectra"
# ==============================

import subprocess
import os

for obsid in obsids:
    evt_file = os.path.join(evt_dir, f"{obsid}_reproj_evt.fits")

    for i in range(1, n_regions + 1):
        reg_file = os.path.join(region_dir, f"reg{i}.reg")
        infile = f"{evt_file}[sky=region({reg_file})]"
        outroot = os.path.join(output_dir, f"{obsid}_src{i}")

        print(f"[{obsid}] Extracting region {i} of {n_regions} → {outroot}.pi")

        subprocess.run(["punlearn", "specextract"], check=True)

        subprocess.run([
            "specextract",
            f"infile={infile}",
            f"outroot={outroot}",
            "grouptype=NONE",
            "mode=h"
        ], check=True)
