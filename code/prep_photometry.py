#!/usr/bin/python 

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits
import astropy.wcs as wcs

def read_file(filename):
    """
    Reads in a photometry .txt file and turns it into a pandas dataframe for easier use
    """
    
    with open(filename) as f:
        galaxy = f.readlines()
    galaxy = pd.DataFrame([x.strip().split() for x in galaxy])
    
    galaxy.columns = ['ID','X', 'Y','Mag','Error', 'Ext_err','Num','chi','sharp','var','blunder']
    galaxy = galaxy.iloc[3:]
    
    return galaxy
    
def assign_radec(filename, image):
    
    hdulist = fits.open(image)
    w = wcs.WCS(hdulist[0].header)
    ra, dec = w.all_pix2world(filename['X'].astype(float), filename['Y'].astype(float), 1)
    
    return ra, dec
