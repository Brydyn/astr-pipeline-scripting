"""
Date Created: 17/06/25
Author: Brydyn Mac Intyre
Description: 
This script processes spectral fitting results stored in text files, extracting parameter values with errors
for different models, datasets, and nH states across multiple regions. It parses the outputs,
splits abundance and flux values into value and error components,
and generates CSV summary files per dataset. The CSV includes detailed parameter info per model
and a summary section highlighting the best-fitting models based on RedChiSq statistics.
"""

import os
import re
import csv
from collections import defaultdict

base_dir = "./Pyxspec_Output"
models = ["powerlaw", "apec", "vapec"]
datasets = ["XMM", "Chandra", "XMM+Chandra"]
nH_states = ["frozen", "nHThaw"]

region_header_re = re.compile(r"^%{10} Region (\d+) %{10}")
rchi_re = re.compile(r"RedChiSq: ([\d.eE+-]+) \((\d+)\)")

# Updated regex to strictly parse value and error terms
triple_val_re = re.compile(
    r"([\-+]?[0-9]*\.?[0-9]+(?:[eE][+\-]?[0-9]+)?)\s+"    # value (float, optional sign and scientific)
    r"([\-][0-9]*\.?[0-9]+(?:[eE][+\-]?[0-9]+)?)\s+"     # error low (must start with '-')
    r"([+][0-9]*\.?[0-9]+(?:[eE][+\-]?[0-9]+)?)"          # error high (must start with '+')
)

def clean_param(name):
    return name.replace(":", "").strip()

def extract_val_err(s):
    m = triple_val_re.match(s)
    if m:
        val = float(m.group(1))
        err_low = float(m.group(2))
        err_high = float(m.group(3))
        return val, err_low, err_high
    try:
        return float(s), "", ""
    except:
        return s.strip(), "", ""

def parse_region(lines):
    result = {}
    for line in lines:
        if line.startswith("RedChiSq:"):
            m = rchi_re.search(line)
            if m:
                result["RedChiSq"] = (float(m.group(1)), "", int(m.group(2)))  # Store dof in Error High
        elif ":" in line:
            key, val = line.split(":", 1)
            key = clean_param(key)
            result[key] = extract_val_err(val.strip())
    return result

# dataset -> region -> model_key -> param dict
all_results = defaultdict(lambda: defaultdict(dict))

for model in models:
    for dataset in datasets:
        for nH_state in nH_states:
            subdir = model if nH_state == "frozen" else f"{model}/nHThaw"
            results_path = os.path.join(base_dir, dataset, subdir, "Results.txt")
            if not os.path.exists(results_path):
                continue

            with open(results_path, "r") as f:
                lines = f.readlines()

            current_region = None
            region_data = []

            for line in lines:
                header_match = region_header_re.match(line.strip())
                if header_match:
                    if current_region is not None:
                        model_key = f"{model} {nH_state}"
                        all_results[dataset][current_region][model_key] = parse_region(region_data)
                    current_region = int(header_match.group(1))
                    region_data = []
                elif current_region is not None:
                    region_data.append(line.strip())

            if current_region is not None and region_data:
                model_key = f"{model} {nH_state}"
                all_results[dataset][current_region][model_key] = parse_region(region_data)

# Write CSV summaries
for dataset, regions in all_results.items():
    summary_file = f"{base_dir}/Summary_{dataset}.csv"
    with open(summary_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)

        for region in sorted(regions.keys()):
            writer.writerow([f"Region {region}"])

            model_keys = [f"{m} {n}" for m in models for n in nH_states]
            header_row = []
            subheader_row = []

            for mk in model_keys:
                header_row += [f"Model: {mk}", "", "", "", ""]
                subheader_row += ["Parameter", "Value", "Error Low", "Error High", ""]

            writer.writerow(header_row)
            writer.writerow(subheader_row)

            # Gather all parameter names
            all_params = set()
            for mk in model_keys:
                all_params.update(regions[region].get(mk, {}).keys())
            all_params = sorted(all_params)

            # Compute best fits
            best_lowest = None
            best_closest = None
            lowest_val = float("inf")
            closest_val = float("inf")
            best_lowest_val = None
            best_closest_val = None

            for mk in model_keys:
                rcs = regions[region].get(mk, {}).get("RedChiSq")
                if rcs:
                    val = rcs[0]
                    if val < lowest_val:
                        lowest_val = val
                        best_lowest = mk
                        best_lowest_val = val
                    if abs(val - 1) < closest_val:
                        closest_val = abs(val - 1)
                        best_closest = mk
                        best_closest_val = val

            for param in all_params:
                row = []
                for mk in model_keys:
                    pdata = regions[region].get(mk, {}).get(param, ("", "", ""))
                    row += [param, *pdata, ""]
                writer.writerow(row)

            # Summary section
            writer.writerow([])
            writer.writerow(["Model Summary"])
            writer.writerow(["RedChiSq", "Value", "Model"])
            writer.writerow(["Lowest", best_lowest_val, best_lowest])
            writer.writerow(["Closest to 1", best_closest_val, best_closest])
            writer.writerow([])

print("CSV summaries written, including per-region model comparison.")
