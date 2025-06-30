# Astronomical Pipeline Scripting

This repository contains modular Python scripts used for the reduction, analysis, and spectral modelling of X-ray data from *Chandra*, *XMM-Newton*, and *NuSTAR*. These tools were developed to automate large-scale multi-region analysis tasks, including batch spectral extraction, pipeline management, and statistical post-processing.

While originally developed for astrophysical jet and supernova remnant studies, many of the utilities are suitable for general high-energy astrophysical workflows.

## Repository Structure

### `/chandra/`
Scripts for processing Chandra data using CIAO:
- `Chandra_multi-src_specextract.py`: Multi-region spectral extraction with automated grouping and file handling.
- `Chandra_U2+U30_All_Images.py`: Image production across multiple energy bands.
- `ftgrouppha_Chandra_multi-src_specextract_output.py`: Batch grouping of extracted spectra with `ftgrouppha`.

### `/xmm/`
XMM-specific pipelines:
- `specextract_v27.py`: Batch spectral extraction for extended and diffuse sources.
- `specextract_v27_ptsrc.py`: Adapted extraction for point sources.

### `/xspec/`
Automated fitting using XSPEC via `pyXspec`:
- `pyxspec_powerlaw.py`: Fits an absorbed and unabsorbed powerlaw model to each region, saving summary results, log files, and plots.
- `pyxspec_apec.py`: (If included) Thermal plasma modelling using the APEC model.
- `pyxspec_vapec.py`: Vapec model with iterative thawing of elemental abundances. Automatically handles freeze/revert logic and exports full results.

### `/utilities/`
General-purpose tools:
- `crop_fits_to_region.py`: Crops FITS files to a specified DS9 region.
- `fitting_pipeline.py`: Meta-script for chaining spectral extraction and fitting workflows.
- `Handmade_Binmap.py`: Creates binning masks by-hand for spatial mapping.
- `Summarize.py`, `SummarizeCSV.py`: Collates result outputs into tables or CSV for downstream processing.

## Requirements

- Python 3.8+
- CIAO (for Chandra scripts)
- HEASoft / SAS (for XMM scripts)
- XSPEC and `pyXspec` Python bindings
- Standard scientific Python stack (`astropy`, `numpy`, `matplotlib`, `scipy`)

## Usage Notes

Each folder is self-contained. Scripts expect certain filename conventions (described in header comments) and assume grouped spectra have been properly background-subtracted.

The output directories (`Pyxspec_Output/...`) will be created automatically. Results include logs, plots, and structured summary files compatible with visualisation tools (e.g., Contbin, DS9 overlays).

## Author

Brydyn Mac Intyre  
Ph.D. Candidate, Department of Physics and Astronomy  
University of Manitoba  
[LinkedIn](https://www.linkedin.com/in/brydynmacintyre)

## License

This code is released under a custom restrictive license. You are free to view and examine the scripts, but any use, distribution, or derivative work requires explicit written permission. See `LICENSE` for details.
