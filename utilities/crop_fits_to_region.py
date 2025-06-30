# -*- coding: utf-8 -*-
"""
Date: 17/06/25

Author: Brydyn Mac Intyre

Crop a FITS image using a rectangular DS9 region.

Usage:
    python3 crop_fits_to_region.py input.fits crop_region.reg

Output:
    input_cropped.fits
"""

import sys
import os
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from regions import Regions

def crop_fits_to_region(fits_file, region_file):
    # Load FITS image and WCS
    with fits.open(fits_file) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        wcs = WCS(header)

    # Read region (assumes one rectangle)
    regions = Regions.read(region_file, format='ds9')
    if len(regions) != 1:
        raise ValueError("Region file must contain exactly one rectangular region.")

    pix_region = regions[0].to_pixel(wcs)
    bbox = pix_region.bounding_box

    # Convert bounding box to integer pixel slice
    x_min = int(np.floor(bbox.ixmin))
    x_max = int(np.ceil(bbox.ixmax))
    y_min = int(np.floor(bbox.iymin))
    y_max = int(np.ceil(bbox.iymax))

    # Handle potential out-of-bounds
    x_min = max(x_min, 0)
    y_min = max(y_min, 0)
    x_max = min(x_max, data.shape[1])
    y_max = min(y_max, data.shape[0])

    # Crop data
    cropped_data = data[y_min:y_max, x_min:x_max]

    # Update header/WCS
    cropped_wcs = wcs.slice((slice(y_min, y_max), slice(x_min, x_max)))
    new_header = cropped_wcs.to_header()
    for key in header:
        if key not in new_header and key not in ['COMMENT', 'HISTORY']:
            new_header[key] = header[key]

    # Save cropped file
    base, ext = os.path.splitext(fits_file)
    output_file = f"{base}_cropped.fits"
    hdu_out = fits.PrimaryHDU(data=cropped_data, header=new_header)
    hdu_out.writeto(output_file, overwrite=True)
    print(f"Cropped file saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 crop_fits_to_region.py <image.fits> <rectangular_region.reg>")
        sys.exit(1)

    crop_fits_to_region(sys.argv[1], sys.argv[2])
