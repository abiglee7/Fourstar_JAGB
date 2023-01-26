#!/usr/bin/python 

from matplotlib import gridspec
import pandas as pd
from astropy.io import fits
import astropy.wcs as wcs

def read_file(filename, image, column_names=['ID','X', 'Y','J','J_err', 'H','H_err','K','K_err','chi','sharp']):
    """
    Reads in a photometry .raw file and turns it into a pandas dataframe for easier use
    """
    
    with open(filename) as f:
        galaxy = f.readlines()
    galaxy = pd.DataFrame([x.strip().split() for x in galaxy])
    
    galaxy.columns = column_names
    galaxy = galaxy.iloc[3:]
    
    
    hdulist = fits.open(image)
    w = wcs.WCS(hdulist[0].header)
    ra, dec = w.all_pix2world(galaxy['X'].astype(float), galaxy['Y'].astype(float), 1)
    
    galaxy['ra'] = ra
    galaxy['dec'] = dec
    
    return galaxy
    
    


