# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 13:29:45 2025

@author: Brydyn Mac Intyre
"""

from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import gaussian_filter
import numpy as np
import os

# === List of input FITS files ===
input_files = [
    "broad_flux.img",
    "0.5-1.0_flux.img",
    "1.0-2.0_flux.img",
    "2.0-10.0_flux.img"
]

for filename in input_files:
    # === Load FITS file ===
    hdul = fits.open(filename)
    data = hdul[0].data
    header = hdul[0].header
    wcs = WCS(header)

    # === Get pixel scale in arcsec/pixel (absolute value of CDELT1 in deg) ===
    pixel_scale = abs(header.get("CDELT1", 1.0)) * 3600

    # === Define Gaussian kernel sizes in pixels ===
    sigma_2arcsec = 2.0 / pixel_scale
    sigma_30arcsec = 30.0 / pixel_scale

    # === Apply Gaussian smoothing ===
    U2 = gaussian_filter(data, sigma=sigma_2arcsec)
    U30 = gaussian_filter(data, sigma=sigma_30arcsec)

    # === Compute unsharp mask ===
    numerator = U2 - U30
    denominator = U2 + U30
    with np.errstate(divide='ignore', invalid='ignore'):
        unsharp = np.true_divide(numerator, denominator)
        unsharp[~np.isfinite(unsharp)] = 0.0

    # === Define output file names ===
    base, _ = os.path.splitext(filename)
    output_U2 = f"{base}_U2.fits"
    output_U30 = f"{base}_U30.fits"
    output_unsharp = f"{base}_unsharp.fits"

    # === Save each result ===
    fits.PrimaryHDU(data=U2, header=header).writeto(output_U2, overwrite=True)
    fits.PrimaryHDU(data=U30, header=header).writeto(output_U30, overwrite=True)
    fits.PrimaryHDU(data=unsharp, header=header).writeto(output_unsharp, overwrite=True)

    print(f"Processed: {filename}")
    print(f"  -> {output_U2}")
    print(f"  -> {output_U30}")
    print(f"  -> {output_unsharp}")
