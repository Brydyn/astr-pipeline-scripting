"""
Automated XSPEC Fitting and Summary Pipeline
Author: Brydyn Mac Intyre
Date: 19-06-2025

This script automates the XSPEC spectral fitting and summarizing process for three models:
- apec
- powerlaw
- vapec

For each model, it:
1. Sets PROCESS to "XMM", "Chandra", and "both", and for each:
   a. Runs the model script with nHFreeze = True
   b. Runs the model script with nHFreeze = False
2. After all model scripts are run, it runs the summary scripts:
   - Summarize_wNames.py
   - SummarizeCSV.py
3. Creates contbin_colours directories under:
   ./Pyxspec_Output/{tele}/{model}/
   ./Pyxspec_Output/{tele}/{model}/nHThaw/
   where tele = XMM, Chandra, or XMM+Chandra, and model = apec, powerlaw, vapec
4. Copies region_list.txt into each of those model directories before running paint_output_images
5. Runs paint_output_images in each contbin_colours folder.

Assumptions:
- You are in the directory containing the scripts.
- Each model script contains exactly one line starting with `nHFreeze =` and `PROCESS =`.
- Model scripts can be modified and restored during execution.
- `paint_output_images` is in PATH and executable.
- binmap path is fixed at /home/Brydyn/HerA/Artificial_contbin/contbin_binmap_cropped.fits
- region_list.txt is present in the working directory.
"""

import subprocess
import os
import shutil

# Change these lines if you need to
    # All other things should be following previous script structures
binmap_path = "/home/Brydyn/HerA/Artificial_contbin/contbin_binmap_cropped.fits"
region_list_src = "region_list.txt"

# Model scripts
models = ["apec", "powerlaw", "vapec"]
tele_list = ["XMM", "Chandra", "XMM+Chandra"]
nH_options = [True, False]

def run_script(script_name):
    print(f"Running {script_name}...")
    subprocess.run(["python3", script_name], check=True)

def modify_script(script_path, process_value, nH_value):
    with open(script_path, 'r') as f:
        lines = f.readlines()

    with open(script_path, 'w') as f:
        for line in lines:
            if line.strip().startswith("PROCESS"):
                f.write(f'PROCESS = "{process_value}"  # Options: "XMM", "Chandra", "both"\n')
            elif line.strip().startswith("nHFreeze"):
                f.write(f"nHFreeze = {nH_value}\n")
            else:
                f.write(line)

# Step 1: Run pyxspec scripts for all combinations of PROCESS and nHFreeze
for model in models:
    script = f"pyxspec_{model}.py"
    for process in tele_list:
        for nH in nH_options:
            print(f"\n=== Running {script} with PROCESS={process}, nHFreeze={nH} ===")
            modify_script(script, process, nH)
            run_script(script)

# Step 2: Run summary scripts
print("\n=== Running summary scripts ===")
run_script("Summarize_wNames.py")
run_script("SummarizeCSV.py")

# Step 3: Create contbin_colours folders, copy region_list.txt, and run paint_output_images
for tele in tele_list:
    for model in models:
        for sub in ["", "nHThaw"]:
            base_path = f"./Pyxspec_Output/{tele}/{model}"
            if sub:
                base_path = os.path.join(base_path, sub)
            colour_dir = os.path.join(base_path, "contbin_colours")
            
            os.makedirs(colour_dir, exist_ok=True)
            
            # Copy region_list.txt into the model directory (not contbin_colours)
            try:
                shutil.copy(region_list_src, base_path)
                print(f"Copied {region_list_src} to {base_path}")
            except Exception as e:
                print(f"Failed to copy {region_list_src} to {base_path}: {e}")
            
            print(f"\n=== Running paint_output_images in {colour_dir} ===")
            paint_cmd = [
                "paint_output_images",
                f"--binmap={binmap_path}",
                "--input_dir=../"
            ]
            subprocess.run(paint_cmd, cwd=colour_dir, check=True)
