"""
Welcome to the generalized script for XMM-Newton spectral extraction.
Author: Brydyn Mac Intyre
Email: macintyb@myumanitoba.ca
Date: Updated [June 2025]

Purpose:
Automate XMM-Newton point source spectral extraction with region exclusion, background handling,
and response generation. Final grouping is performed using ftgrouppha + grppha at S/N ≥ 3.0.

Assumptions:
- Point source extraction
- DS9 PHYSICAL region files
- Each source region grouped to all background regions (1:many)

Inputs:
- dir.txt             : one full path to each filtered FITS file
- cifdir.txt          : matching order CIF paths for each FITS file
- *_Regions.txt/.reg  : DS9-style PHYSICAL region files per observation (marked "background" or "-")
                        The script automatically detects .txt or .reg extensions for region files.

Output Naming Convention:
- Intermediate:    {CAM}_{OBS}_pi_reg{i}.fits, _rmf_reg{i}.fits, _arf_reg{i}.fits
- Background:      ..._pi_backreg{j}.fits
- Grouped Final:   ./Grouped/{OBS}_{CAM}_src{i}_grp3sn.fits

    Example dir.txt:
/home/usrNam/XMM/0123456789/PROC/anyFileName.fits
/home/secObs/oddName/1234567890_mos1_filt.fits
/home/usr2/XMM/0123456789/PROC/pn_0123456789_softfilt.fits

    Example cifdir.txt:
/home/usrNam/XMM/0123456789/ODF/ccf.cif
/home/secObs/oddName/huh.cif
/home/usr2/XMM/1234567890/ODF/ccf.cif

    Example 0000000000_Regions.txt:
# Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
physical
polygon(25730.343,27702.546,26724.297,27702.153,26725.732) # color=cyan
-polygon(25828.608,27207.606,26747.072,27587.613) # color=cyan
polygon(25965.214,27702.546,26824.237,27782.553,26555.622) # color=cyan background
"""

import os
import numpy as np
from astropy.io import fits
import subprocess
import sys

# === PARAMETERS ===
group_snr = 3.0
grouped_dir = "Grouped"
os.makedirs(grouped_dir, exist_ok=True)

# === FILE LOADING ===
try:
    filePath = np.loadtxt("dir.txt", dtype=str)
    cifPath = np.loadtxt("cifdir.txt", dtype=str)
except Exception as e:
    print(f"Error loading input files: {e}")
    sys.exit(1)

camNum = len(filePath)
if len(filePath) != len(cifPath):
    print("Mismatch between number of FITS files and CIF paths.")
    sys.exit(1)

camType, camName, obsID = [], [], []
for path in filePath:
    try:
        with fits.open(path) as hdul:
            header = hdul[0].header
            camName.append(header['INSTRUME'])
            camType.append(f"{header['INSTRUME']}_{header['EXPIDSTR']}")
            obsID.append(header['OBS_ID'])
    except Exception as e:
        print(f"Error reading FITS header from {path}: {e}")
        sys.exit(1)

filePath, cifPath = map(np.array, (filePath, cifPath))
camType, camName, obsID = map(np.array, (camType, camName, obsID))

# Sort by OBSID
sorter = np.argsort(obsID)
filePath, cifPath, camName, camType, obsID = [arr[sorter] for arr in [filePath, cifPath, camName, camType, obsID]]

uniObs = np.unique(obsID)
obsNum = len(uniObs)

# === PARSE REGIONS WITH LOCAL NUMBERING ===
include, exclude, backRegs = [], [], []
include_indices_per_obs = [[] for _ in range(obsNum)]
background_indices_per_obs = [[] for _ in range(obsNum)]
obsRegCounterInc = np.zeros(obsNum, dtype=int)
obsRegCounterExc = np.zeros(obsNum, dtype=int)
obsRegCounterBac = np.zeros(obsNum, dtype=int)
removal = []

for obs_idx, obs in enumerate(uniObs):
    reg_file = f"{obs}_Regions.txt" if os.path.exists(f"{obs}_Regions.txt") else f"{obs}_Regions.reg"
    if not os.path.exists(reg_file):
        print(f"No region file found for observation {obs}")
        sys.exit(1)

    with open(reg_file) as f:
        region_lines = f.readlines()[3:]
        local_index = 1
        exclusion_regions = []
    
        for line in region_lines:
            full_line = line.strip()
            is_background = "background" in full_line.lower()
            is_excluded = full_line.strip().startswith("-")
    
            # Remove leading minus sign and strip trailing comments
            clean_line = full_line.lstrip("-").split("#")[0].strip()
    
            if is_background:
                backRegs.append(clean_line)
                background_indices_per_obs[obs_idx].append(local_index)
                obsRegCounterBac[obs_idx] += 1
            elif is_excluded:
                exclude.append(clean_line)
                exclusion_regions.append(clean_line)
                obsRegCounterExc[obs_idx] += 1
            else:
                include.append(clean_line)
                include_indices_per_obs[obs_idx].append(local_index)
                obsRegCounterInc[obs_idx] += 1
    
            local_index += 1

        if exclusion_regions:
            excl_expr = "&&!".join([f"((X,Y) in {e})" for e in exclusion_regions])
            excl_expr = f"&&!{excl_expr}"
        else:
            excl_expr = ""
        removal.append(excl_expr)

# Flatten per-obs indices
include_indices = [idx for sublist in include_indices_per_obs for idx in sublist]
background_indices = [idx for sublist in background_indices_per_obs for idx in sublist]

obsRegCounterIncSum = np.cumsum(obsRegCounterInc)
obsRegCounterBacSum = np.cumsum(obsRegCounterBac)

print(f"Info: Found {len(include)} total source regions and {len(backRegs)} background regions across {obsNum} unique observations and {camNum} cameras.")

# === EXTRACTION + RESPONSE GENERATION ===
srcNames, bkgNames = [], []

obsTracker2 = 0
os.environ['SAS_CCF'] = cifPath[obsTracker2]

for n in range(camNum):
    if obsID[n] != uniObs[obsTracker2]:
        obsTracker2 += 1
        os.environ['SAS_CCF'] = cifPath[obsTracker2]
    start = obsRegCounterIncSum[obsTracker2 - 1] if obsTracker2 > 0 else 0
    stop = obsRegCounterIncSum[obsTracker2]

    for s in range(start, stop):
        reg_number = include_indices[s]
        filtered = f"{camType[n]}_{obsID[n]}_filtered_reg{reg_number}.fits"
        spectrum = f"{camType[n]}_{obsID[n]}_pi_reg{reg_number}.fits"
        rmf = f"{camType[n]}_{obsID[n]}_rmf_reg{reg_number}.fits"
        arf = f"{camType[n]}_{obsID[n]}_arf_reg{reg_number}.fits"
        expr = f"(FLAG==0)&&(PATTERN <= 4)&&((X,Y) in {include[s]}){removal[obsTracker2]}" if camName[n] == "EPN" else f"((X,Y) in {include[s]}){removal[obsTracker2]}"
        chan_max = "20479" if camName[n] == "EPN" else "11999"

        try:
            os.system(f"evselect table='{filePath[n]}' energycolumn='PI' withfilteredset=yes filteredset={filtered} keepfilteroutput=yes filtertype='expression' expression='{expr}' withspectrumset=yes spectrumset='{spectrum}' spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax={chan_max}")
            os.system(f"backscale spectrumset={spectrum} badpixlocation={filePath[n]}")
            os.system(f"rmfgen rmfset={rmf} spectrumset={spectrum}")
            os.system(f"arfgen arfset={arf} spectrumset={spectrum} withrmfset=yes rmfset={rmf} withbadpixcorr=yes badpixlocation={filePath[n]}")
        except Exception as e:
            print(f"Error during extraction for {spectrum}: {e}")
            continue

        srcNames.append((obsID[n], camType[n], reg_number, spectrum, rmf, arf))

# === BACKGROUND EXTRACTION ===
obsTracker3 = 0
if len(backRegs) == 0:
    print("No background regions found. Cannot proceed to grouping.")
    sys.exit(1)

os.environ['SAS_CCF'] = cifPath[obsTracker3]

for j in range(camNum):
    if obsID[j] != uniObs[obsTracker3]:
        obsTracker3 += 1
        os.environ['SAS_CCF'] = cifPath[obsTracker3]

    start = obsRegCounterBacSum[obsTracker3 - 1] if obsTracker3 > 0 else 0
    stop = obsRegCounterBacSum[obsTracker3]

    for q in range(start, stop):
        reg_number = background_indices[q]
        filtered = f"{camType[j]}_{obsID[j]}_filtered_backreg{reg_number}.fits"
        bkg = f"{camType[j]}_{obsID[j]}_pi_backreg{reg_number}.fits"
        expr = f"(FLAG==0)&&(PATTERN <= 4)&&((X,Y) in {backRegs[q]}){removal[obsTracker3]}" if camName[j] == "EPN" else f"((X,Y) in {backRegs[q]}){removal[obsTracker3]}"
        chan_max = "20479" if camName[j] == "EPN" else "11999"

        try:
            os.system(f"evselect table='{filePath[j]}' energycolumn='PI' withfilteredset=yes filteredset={filtered} keepfilteroutput=yes filtertype='expression' expression='{expr}' withspectrumset=yes spectrumset='{bkg}' spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax={chan_max}")
            os.system(f"backscale spectrumset={bkg} badpixlocation={filePath[j]}")
        except Exception as e:
            print(f"Error extracting background spectrum {bkg}: {e}")
            continue

        bkgNames.append((obsID[j], camType[j], reg_number, bkg))

# === GROUPING ===
grouped_count = 0
for obsid, cam, src_idx, pi, rmf, arf in srcNames:
    for bkg_obs, bkg_cam, bkg_idx, bkg in bkgNames:
        if obsid == bkg_obs and cam == bkg_cam:
            grouped_out = os.path.join(grouped_dir, f"{obsid}_{cam}_src{src_idx}_grp3sn.fits")
            temp_grp = f"{obsid}_{cam}_src{src_idx}_back{bkg_idx}_temp_grp.fits"

            # Convert to absolute paths
            pi_abs = os.path.abspath(pi)
            bkg_abs = os.path.abspath(bkg)
            rmf_abs = os.path.abspath(rmf)
            arf_abs = os.path.abspath(arf)

            try:
                subprocess.run([
                    "ftgrouppha",
                    f"infile={pi_abs}",
                    f"backfile={bkg_abs}",
                    f"respfile={rmf_abs}",
                    f"outfile={temp_grp}",
                    f"grouptype=snmin",
                    f"groupscale={group_snr}"
                ], check=True)

                comm = f"chkey ANCRFILE {arf_abs} & chkey RESPFILE {rmf_abs} & exit"
                subprocess.run([
                    "grppha",
                    f"infile={temp_grp}",
                    f"outfile={grouped_out}",
                    f"comm={comm}"
                ], check=True)

                os.remove(temp_grp)
                grouped_count += 1
            except subprocess.CalledProcessError as e:
                print(f"Grouping failed for {pi} with background {bkg}: {e}")

print(f"\nExtraction complete:")
print(f"  Sources extracted:   {len(srcNames)}")
print(f"  Backgrounds used:    {len(bkgNames)}")
print(f"  Grouped products:    {grouped_count}")
