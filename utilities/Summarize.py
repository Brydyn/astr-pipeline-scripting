# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 2025

@author: Brydyn Mac Intyre

Purpose: Parse all Results.txt files created in subfolders to find the best fits.
"""

import os
import re

# Paths and options
base_dir = "./Pyxspec_Output"
models = ["powerlaw", "apec", "vapec"]
datasets = ["XMM", "Chandra", "XMM+Chandra"]
nH_states = ["frozen", "nHThaw"]

# Regex patterns
region_header_re = re.compile(r"^%{10} Region (\d+) %{10}")
rchi_re = re.compile(r"RedChiSq: ([\d.]+) \((\d+)\)")

# Track best fits by two criteria
region_best_fit_lowest = {}
region_best_fit_closest = {}

# --- Parse a region block ---
def parse_region(lines):
    result = {}
    for line in lines:
        if line.startswith("RedChiSq:"):
            m = rchi_re.search(line)
            if m:
                result["RedChiSq"] = float(m.group(1))
                result["dof"] = int(m.group(2))
        elif ":" in line:
            key, value = line.split(":", 1)
            result[key.strip()] = value.strip()
    return result

# --- Parse Results.txt file ---
def parse_results_file(path, model, dataset, nH_state):
    if not os.path.exists(path):
        return {}

    with open(path, "r") as f:
        lines = f.readlines()

    results = {}
    current_region = None
    current_data = []

    for line in lines:
        header_match = region_header_re.match(line.strip())
        if header_match:
            if current_region is not None:
                results[current_region] = parse_region(current_data)
                results[current_region]["model"] = model
                results[current_region]["dataset"] = dataset
                results[current_region]["nH_state"] = nH_state
            current_region = int(header_match.group(1))
            current_data = []
        elif current_region is not None:
            current_data.append(line.strip())

    # Final region
    if current_region is not None and current_data:
        results[current_region] = parse_region(current_data)
        results[current_region]["model"] = model
        results[current_region]["dataset"] = dataset
        results[current_region]["nH_state"] = nH_state

    return results

# --- Scan and collect best fits ---
for model in models:
    for dataset in datasets:
        for nH_state in nH_states:
            subdir = model if nH_state == "frozen" else f"{model}/nHThaw"
            results_path = os.path.join(base_dir, dataset, subdir, "Results.txt")
            parsed = parse_results_file(results_path, model, dataset, nH_state)

            for region, data in parsed.items():
                if "RedChiSq" not in data:
                    print(f"Skipping region {region} in {results_path}: Missing RedChiSq")
                    continue

                rchi = data["RedChiSq"]

                # For lowest RedChiSq
                best_lowest = region_best_fit_lowest.get(region)
                if best_lowest is None or rchi < best_lowest.get("RedChiSq", float("inf")):
                    region_best_fit_lowest[region] = data

                # For closest-to-1 RedChiSq
                best_closest = region_best_fit_closest.get(region)
                if best_closest is None or abs(rchi - 1) < abs(best_closest.get("RedChiSq", float("inf")) - 1):
                    region_best_fit_closest[region] = data

# --- Output function ---
def write_summary(summary_path, region_dict):
    with open(summary_path, "w") as summary:
        for region in sorted(region_dict.keys()):
            data = region_dict[region]
            summary.write(f"\n{'='*20} Region {region} {'='*20}\n\n")
            summary.write(f"Best Fit Model : {data['model']}\n")
            summary.write(f"Data Used      : {data['dataset']}\n")
            summary.write(f"nH Setting     : {data['nH_state']}\n")
            summary.write(f"RedChiSq       : {data.get('RedChiSq', 'N/A')} ({data.get('dof', '?')} dof)\n")
            summary.write("\nParameters:\n")
            for key, val in data.items():
                if key not in {"model", "dataset", "nH_state", "RedChiSq", "dof"}:
                    summary.write(f"  {key}: {val}\n")

# --- Write both summaries ---
summary0_path = os.path.join(base_dir, "Summary0.txt")
summary1_path = os.path.join(base_dir, "Summary1.txt")

write_summary(summary0_path, region_best_fit_lowest)
write_summary(summary1_path, region_best_fit_closest)

print(f"Summary0 (lowest RCS) written to {summary0_path}")
print(f"Summary1 (closest to 1 RCS) written to {summary1_path}")