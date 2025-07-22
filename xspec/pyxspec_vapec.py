"""
XSPEC Batch Spectral Fitter for Absorbed and Unabsorbed vapec Model with Abundance Iteration
By Brydyn Mac Intyre — Updated 16/06/2025

Description:
    Fits tbabs*cflux*vapec with:
      - Initial fit with all abundances frozen at 1
      - Then thaw abundances one at a time, refitting
      - If abundance out of [0.1,10], reset to 1 and freeze again
      - Skip abundance thawing if reduced chi-square ≥ 2
      - Outputs final fit parameters only

Requirements / Restrictions:
    - pyXspec must be installed
    - Grouped spectra must be named as follows:
        - XMM: <obsid>_<inst>_src<region>_grp3sn.fits (e.g., 0401730101_EMOS1_src1_grp3sn.fits)
        - Chandra: <obsid>_src<region>_grp3sn.fits (e.g., 1625_src1_grp3sn.fits)
    - Backgrounds must already be accounted for in spectra grouping (no explicit background files used here).
    - Grouping convention: 
        - MOS1 and MOS2 in XSPEC data groups 1 and 2
        - pn in XSPEC data group 3
        - Chandra in XSPEC data group 4

Output:
    - Model fit outputs saved in ./Pyxspec_Output/{Instru}/vapec/
    - Fit summary saved per region to: xaf_<region>_fit_out.txt
    - Master log of all fits written to: ./Pyxspec_Output/{Instru}/vapec/Results.txt
    - Plot images saved for each region
"""

import os
import re
from xspec import *
import subprocess

# Definitions
lowBound = 0.4
upBound = 7.0

# Set starting values and freeze True/False
nHSet = 0.0529
nHFreeze = True
kTSet = 1.0
kTFreeze = False
AbundSet = 1.0
AbundFreeze = False
zSet = 0.1549
zFreeze = True

# Which telescopes are being fit?
PROCESS = "XMM+Chandra"  # Options: "XMM", "Chandra", "XMM+Chandra"

# Make output directories
if nHFreeze == True:
    output_dir = f"./Pyxspec_Output/{PROCESS}/vapec"
else:
    output_dir = f"./Pyxspec_Output/{PROCESS}/vapec/nHThaw"
os.makedirs(output_dir, exist_ok=True)

# Count number of regions based on filenames
reg_pattern = re.compile(r'_src(\d+)_grp3sn.fits')
max_region = -1
for f in os.listdir("."):
    m = reg_pattern.search(f)
    if m:
        r = int(m.group(1))
        if r > max_region:
            max_region = r
regNum = max_region + 1

# Identify all spectrum files
all_files = [f for f in os.listdir(".") if f.endswith("_grp3sn.fits")]

# Correct abundance parameter indices for this model:
# param 1 = const
# params 2-4 = cflux
# param 5 = tbabs nH
# param 6 = vapec kT
# params 7-19 = vapec abundances (He, C, N, O, Ne, Mg, Al, Si, S, Ar, Ca, Fe, Ni)
# param 20 = vapec redshift
ABUNDANCE_PARAMS = list(range(7, 20))  # 7 to 19 inclusive
REDSHIFT_PARAM = 20

ABUNDANCE_LABELS = {
    7: "He", 8: "C", 9: "N", 10: "O", 11: "Ne", 12: "Mg", 13: "Al",
    14: "Si", 15: "S", 16: "Ar", 17: "Ca", 18: "Fe", 19: "Ni"
}

# Open master log file
master_log = open(os.path.join(output_dir, "Results.txt"), "w")

for reg in range(regNum):
    print(f"Processing region {reg}...")
    master_log.write(f"\n%%%%%%%%%% Region {reg} %%%%%%%%%%\n\n")
    
    # Organize spectra by instrument type
    mos1_group = []
    mos2_group = []
    pn_group = []
    cxo_group = []

    for f in all_files:
        if f"_src{reg}_" not in f:
            continue
        if "EMOS1" in f:
            mos1_group.append(f)
        elif "EMOS2" in f:
            mos2_group.append(f)
        elif "EPN" in f:
            pn_group.append(f)
        else:
            cxo_group.append(f)

    k = 1
    group_num = 1
    data_cmd = ""
    
    if PROCESS in ("XMM", "XMM+Chandra"):
        if mos1_group:
            for file in mos1_group:
                data_cmd += f"{group_num}:{k} {file} "
                k += 1
            group_num += 1
        if mos2_group:
            for file in mos2_group:
                data_cmd += f"{group_num}:{k} {file} "
                k += 1
            group_num += 1
        if pn_group:
            for file in pn_group:
                data_cmd += f"{group_num}:{k} {file} "
                k += 1
            group_num += 1
    
    if PROCESS in ("Chandra", "XMM+Chandra"):
        if cxo_group:
            for file in cxo_group:
                data_cmd += f"{group_num}:{k} {file} "
                k += 1
            group_num += 1

    if not (mos1_group or mos2_group or pn_group or cxo_group):
        print(f"No grouped spectra found for region {reg}. Skipping.")
        master_log.write(f"No grouped spectra found for region {reg}. Skipping.\n")
        continue

    # Load into XSPEC
    AllData.clear()
    AllModels.clear()
    AllData(data_cmd)
    AllData.ignore(f"**-{lowBound},{upBound}-**")
    AllData.ignore("bad")

    # ABSORBED model: const*cflux*tbabs*vapec
    m1 = Model("const*cflux*tbabs*vapec")
    
    # Set constants for all data groups
    m1(1).frozen = True  # first constant frozen
    for i in range(2, AllData.nGroups + 1):
        try:
            g = AllModels(i)
            g(1).values = 1.0
            g(1).frozen = False
        except Exception as e:
            print(f"Skipping group {i}: {e}")
    
    # Set cflux bounds
    m1(2).values = lowBound
    m1(2).frozen = True
    m1(3).values = upBound
    m1(3).frozen = True
    m1(4).frozen = False  # log flux
    
    m1(5).values = nHSet  # nH
    m1(5).frozen = nHFreeze  # constant frozen for nH 
    
    m1(6).values = kTSet  # kT in KeV
    m1(6).frozen = kTFreeze

    # Freeze all abundances to 1 initially
    for p in ABUNDANCE_PARAMS:
        m1(p).values = 1.0
        m1(p).frozen = True
        
    m1(REDSHIFT_PARAM).values = zSet  # Redshift
    m1(REDSHIFT_PARAM).frozen = zFreeze

    Fit.query = "no"
    Fit.renorm()
    Fit.nIterations = 1000   # Prevent infinit fitting
    try:
        Fit.perform()
    except Exception as e:
        print(f"Fit failed for region {reg}: {e}")
        master_log.write(f"Fit failed for region {reg}: {e}\n")
        continue
    
    chi2 = Fit.statistic
    dof = Fit.dof
    rchi = chi2 / dof
    
    # If poor fit, skip abundance thawing and write final results
    if rchi >= 2.0:
        master_log.write(f"Initial fit for region {reg} poor (RedChiSq={rchi:.2f}). Skipping abundance thawing.\n")
    else:
        # Iterate thawing abundances one-by-one with rollback if out of bounds
        p_list = {}
        
        for p in ABUNDANCE_PARAMS:
            m1(p).frozen = False
            Fit.renorm()
            try:
                Fit.perform()
            except Exception as e:
                print(f"Fit failed for region {reg}: {e}")
                master_log.write(f"Fit failed for region {reg}: {e}\n")
        
            out_of_bounds = False
            current_val = m1(p).values[0]
        
            if current_val < 0.1 or current_val > 10.0:
                out_of_bounds = True
            else:
                # Check if thawing this parameter caused earlier ones to go out of bounds
                for q in p_list:
                    q_val = m1(q).values[0]
                    if q_val < 0.1 or q_val > 10.0:
                        out_of_bounds = True
                        break
        
            if out_of_bounds:
                # Revert the current one
                m1(p).values = 1.0
                m1(p).frozen = True
        
                # Revert and thaw previous good ones
                for q, val in p_list.items():
                    m1(q).values = val
                    m1(q).frozen = False
        
                # Refit to stabilize
                Fit.renorm()
                try:
                    Fit.perform()
                except Exception as e:
                    print(f"Refit after rollback failed for region {reg}: {e}")
                    master_log.write(f"Refit after rollback failed for region {reg}: {e}\n")
            else:
                # Store accepted value
                p_list[p] = current_val

    master_log.write(f"RedChiSq: {rchi:.3f} ({dof})\n")

    if rchi < 2:
        Fit.error("2.706 4")
        Fit.error("2.706 5")
        Fit.error("2.706 6")
        Fit.error("2.706 7")
        Fit.error("2.706 8")
    for p in ABUNDANCE_PARAMS:
        if not m1(p).frozen:
            try:
                Fit.error(f"2.706 {p}")
            except Exception:
                pass

    # Save absorbed model
    Xset.save(os.path.join(output_dir, f"xcm_all_vapec+cfluxAbs_r{reg}"), info="a")
    Xset.save(os.path.join(output_dir, f"xcm_model_vapec+cfluxAbs_r{reg}"), info="m")

    # Errors
    cFluxAbs = 10**m1(4).values[0]
    if hasattr(m1(4), "error") and m1(4).error:
        cFluxAbsErrorLow = 10**m1(4).error[0] or 0
        cFluxAbsErrorHigh = 10**m1(4).error[1] or 0
    else:
        cFluxAbsErrorLow = 0
        cFluxAbsErrorHigh = 0
    nH = m1(5).values[0]
    if hasattr(m1(5), "error") and m1(5).error:
        nHErrorLow = m1(5).error[0] or 0
        nHErrorHigh = m1(5).error[1] or 0
    else:
        nHErrorLow = 0
        nHErrorHigh = 0
    kT = m1(6).values[0]
    if hasattr(m1(6), "error") and m1(6).error:
        kTErrorLow = m1(6).error[0] or 0
        kTErrorHigh = m1(6).error[1] or 0
    else:
        kTErrorLow = 0
        kTErrorHigh = 0
    zshift = m1(REDSHIFT_PARAM).values[0]
    if hasattr(m1(REDSHIFT_PARAM), "error") and m1(8).error:
        zshiftErrorLow = m1(REDSHIFT_PARAM).error[0] or 0
        zshiftErrorHigh = m1(REDSHIFT_PARAM).error[1] or 0
    else:
        zshiftErrorLow = 0
        zshiftErrorHigh = 0
        
    # Get +/-
    cFluxAbsErrorLow = cFluxAbs - cFluxAbsErrorLow
    cFluxAbsErrorHigh = cFluxAbsErrorHigh - cFluxAbs
    nHErrorLow = nH - nHErrorLow
    nHErrorHigh = nHErrorHigh - nH
    kTErrorLow = kT - kTErrorLow
    kTErrorHigh = kTErrorHigh - kT
    zshiftErrorLow = zshift - zshiftErrorLow
    zshiftErrorHigh = zshiftErrorHigh - zshift

    abundances = {}
    for p in ABUNDANCE_PARAMS:
        val = m1(p).values[0]
        err_low, err_high = (0,0)
        if hasattr(m1(p), "error") and m1(p).error:
            err_low = m1(p).error[0]
            err_high = m1(p).error[1]
            err_low = val - err_low
            err_high = err_high - val
        frozen = m1(p).frozen
        abundances[p] = (val, err_low, err_high, frozen)

    # Write to master log of all region fits
    master_log.write(f"nH: {nH} -{nHErrorLow} +{nHErrorHigh}\n")
    master_log.write(f"cFluxAbs: {cFluxAbs} -{cFluxAbsErrorLow} +{cFluxAbsErrorHigh}\n")
    master_log.write(f"kT: {kT} -{kTErrorLow} +{kTErrorHigh}\n")
    for p, (val, err_lo, err_hi, frozen) in abundances.items():
        label = ABUNDANCE_LABELS.get(p, f"Param{p}")
        status = "frozen" if frozen else "thawed"
        master_log.write(f"{label}: {val:.4f} -{err_lo:.4f} +{err_hi:.4f} ({status})\n")
    master_log.write(f"Redshift: {zshift} -{zshiftErrorLow} +{zshiftErrorHigh}\n")

    # Plot and save absorbed fit
    Plot.xAxis = "keV"
    Plot.addCommand("time off")
    Plot.addCommand("r y 0.00001 1.0")
    Plot.addCommand(f"r x {lowBound} {upBound}")
    Plot.addCommand(f"label Top Reg{reg} vapec")
    Plot.addCommand("plot")
    Plot.addCommand(f"h {output_dir}/r{reg}_vapec.ps/cps")
    Plot("ldata delchi")
    Plot.commands = ()

    # UNABSORBED model: const*tbabs*cflux*vapec
    m2 = Model("const*tbabs*cflux*vapec")
    
    # Set constants for all data groups
    m2(1).frozen = True  # first constant frozen
    for i in range(2, AllData.nGroups + 1):
        try:
            g = AllModels(i)
            g(1).values = 1.0
            g(1).frozen = False
        except Exception as e:
            print(f"Skipping group {i}: {e}")
            
    m2(2).values = nH  # nH
    m2(2).frozen = True  # constant frozen for nH 
    
    # Set cflux bounds
    m2(3).values = lowBound
    m2(3).frozen = True
    m2(4).values = upBound
    m2(4).frozen = True
    m2(5).frozen = False  # log flux
            
    m2(6).values = kT  # kT in KeV
    m2(6).frozen = True

    for p in ABUNDANCE_PARAMS:  # Solar Abundances
        val, _, _, frozen = abundances[p]
        m2(p).values = val
        m2(p).frozen = True
    
    m2(REDSHIFT_PARAM).values = zshift  # Redshift
    m2(REDSHIFT_PARAM).frozen = True

    Fit.renorm()
    Fit.nIterations = 1000   # Prevent infinit fitting
    try:
        Fit.perform()
    except Exception as e:
        print(f"Fit failed for region {reg}: {e}")
        master_log.write(f"Fit failed for region {reg}: {e}\n")
        continue
    if rchi <= 2:
        Fit.error("2.706 5")

    cFluxUnAbs = 10**m2(5).values[0]
    if hasattr(m2(5), "error") and m2(5).error:
        cFluxUnAbsErrorLow = 10**m2(5).error[0]
        cFluxUnAbsErrorHigh = 10**m2(5).error[1]
    else:
        cFluxUnAbsErrorLow = 0
        cFluxUnAbsErrorHigh = 0
    
    cFluxUnAbsErrorLow = cFluxUnAbs - cFluxUnAbsErrorLow
    cFluxUnAbsErrorHigh = cFluxUnAbsErrorHigh - cFluxUnAbs

    master_log.write(f"cFluxUnAbs: {cFluxUnAbs} -{cFluxUnAbsErrorLow} +{cFluxUnAbsErrorHigh}\n")

    # Save contbin-readable file
    with open(os.path.join(output_dir, f"xaf_{reg}_fit_out.txt"), "w") as cf:
        cf.write(f"RedChiSqu {rchi}\n")
        cf.write(f"nH {nH}\n")
        cf.write(f"kT {kT}\n")
        for p, (val, _, _, frozen) in abundances.items():
            label = ABUNDANCE_LABELS.get(p, f"p{p}")
            cf.write(f"{label} {val:.4f}\n")
        cf.write(f"Redshift {zshift}\n")
        cf.write(f"FluxAbs {cFluxAbs}\n")
        cf.write(f"FluxUnAbs {cFluxUnAbs}\n")

# Convert plots and clean .ps files
# Change to output directory and run mogrify
subprocess.run(f"cd {output_dir} && mogrify -rotate 90 -format png *.ps", shell=True)
# Delete .ps files
subprocess.run(f"rm {output_dir}/*.ps", shell=True)


master_log.close()
