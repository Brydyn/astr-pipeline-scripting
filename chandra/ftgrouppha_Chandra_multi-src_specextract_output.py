# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 2025

@author: Brydyn Mac Intyre
"""


"""
Purpose:
    Group spectra from multiple ObsIDs and regions by a specified signal-to-noise ratio (SNR).
    Uses `ftgrouppha` to bin the spectrum and `grppha` to write correct header keywords.

Requirements:
    - HEASoft installed and initialized
    - Output from `specextract` must already exist (i.e. .pi, .arf, .rmf files)

Inputs:
    - PI files:       {spectra_dir}/{obsid}_src{i}.pi
    - ARF/RMF files:  {spectra_dir}/{obsid}_src{i}.{arf|rmf}
    - Background PI:  Uses one specified region number for all other regions

    WARNING:
    All grouped spectra will use the same background region for each ObsID, defined below.

Outputs:
    - Grouped FITS:   {grouped_dir}/{obsid}_src{i}_grp{snr}sn.fits

"""

# ==== FILL IN THESE VALUES ====
obsids = ["1625", "5796", "6257"]          # ObsIDs to process
n_regions = 27                             # Number of extracted regions per ObsID
background_region = 9                      # Region number to use as background (used for ALL)
spectra_dir = "/home/Brydyn/HerA/Spectra"  # Directory with extracted spectra
snr = 3.0                                  # Signal-to-noise ratio for binning
# ==============================

import subprocess
import os

grouped_dir = os.path.join(spectra_dir, "Grouped")
os.makedirs(grouped_dir, exist_ok=True)

# Convert SNR to string for use in filenames (e.g., 3 → "3", 3.0 → "3")
snr_str = str(round(snr))

for obsid in obsids:
    bkg_file = f"{spectra_dir}/{obsid}_src{background_region}.pi"

    for i in range(1, n_regions + 1):
        pi_file = f"{spectra_dir}/{obsid}_src{i}.pi"
        arf_file = f"{spectra_dir}/{obsid}_src{i}.arf"
        rmf_file = f"{spectra_dir}/{obsid}_src{i}.rmf"

        temp_file = f"{obsid}_src{i}_grp{snr_str}sn_temp.fits"
        grouped_file = os.path.join(grouped_dir, f"{obsid}_src{i}_grp{snr_str}sn.fits")

        print(f"[{obsid}] Grouping region {i} → {grouped_file}")

        # Step 1: group by SNR using ftgrouppha
        subprocess.run([
            "ftgrouppha",
            f"infile={pi_file}",
            f"backfile={bkg_file}",
            f"respfile={rmf_file}",
            f"outfile={temp_file}",
            f"grouptype=snmin",
            f"groupscale={snr}"
        ], check=True)

        # Step 2: set correct ARF and RMF headers using grppha
        comm = f"chkey ANCRFILE {arf_file} & chkey RESPFILE {rmf_file} & exit"
        subprocess.run([
            "grppha",
            f"infile={temp_file}",
            f"outfile={grouped_file}",
            f"comm={comm}"
        ], check=True)

        # Cleanup: remove temporary file
        os.remove(temp_file)