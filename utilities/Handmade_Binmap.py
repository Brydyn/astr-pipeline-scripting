# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 2025

@author: Brydyn Mac Intyre

Purpose: Create handmade_binmaps when having hand-selected all your regions!

Example usage: 
    python3 Handmade_Binmap.py xray_image.fits my_regions.reg my_binmap.fits
or
    python3 Handmade_Binmap.py xray_image.fits my_regions.reg
    
Note: For any regions that you do not want drawn, replace with small regions 
    away from the source in question. (I used circle(252,4.75,0.1") for example.)
"""

from astropy.io import fits
from astropy.wcs import WCS
from regions import Regions
import numpy as np
import sys

# Corrected import path for DS9ParserError
try:
    from regions.io.ds9.core import DS9ParserError
except ImportError:
    DS9ParserError = Exception  # fallback to generic exception

def create_binmap_from_regions(fits_image_path, region_file, output_binmap_path="handmade_binmap.fits"):
    # Load the base image
    hdu = fits.open(fits_image_path)
    data = hdu[0].data
    if data.ndim > 2:
        data = data[0]  # Handle data cubes by selecting the first slice
    header = hdu[0].header
    wcs = WCS(header)
    
    binmap = np.zeros_like(data, dtype=np.int16)  # Output image

    # Load region definitions
    try:
        regions = Regions.read(region_file, format='ds9')
    except DS9ParserError as e:
        raise RuntimeError(f"Error reading region file: {e}")

    for idx, region in enumerate(regions, start=1):  # 1-based bin ID
        # Convert region from sky to pixel coordinates
        try:
            pix_region = region.to_pixel(wcs)
        except Exception as e:
            print(f"Could not convert region {idx} to pixel: {e}")
            continue

        # Create a mask for the region
        mask = pix_region.to_mask(mode='center')
        if mask is None:
            print(f"No mask produced for region {idx}")
            continue

        region_mask = mask.to_image(data.shape)
        if region_mask is None:
            print(f"Empty mask for region {idx}")
            continue

        mask_bool = region_mask.astype(bool)
        binmap[mask_bool] = idx  # Overwrites pixels regardless of previous assignment

    # Save the binmap as a FITS file
    hdu_out = fits.PrimaryHDU(binmap, header=header)
    hdu_out.writeto(output_binmap_path, overwrite=True)
    print(f"Saved binmap to {output_binmap_path}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 Handmade_Binmap.py <fits_image.fits> <region_file.reg> [output_binmap.fits]")
        sys.exit(1)

    fits_image = sys.argv[1]
    region_file = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) > 3 else "handmade_binmap.fits"

    create_binmap_from_regions(fits_image, region_file, output_file)
