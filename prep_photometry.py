"""
Author: Abigail Lee

This script reads in raw .txt photometry files (read_file), cleans the J band magnitudes based on their chi, sharp, and error values, and outputs a cleaned .csv file for each epoch(clean_photometry).

This script also checks for 2MASS zeropoints.
"""

#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
from matplotlib import gridspec
from astroquery.vizier import Vizier
import astropy.units as u
from astropy.coordinates import SkyCoord


def read_file(filename):

    """
    Reads in a photometry .txt file and turns it into a pandas dataframe for easier use

    filename: .txt file
    """

    with open(filename) as f:
        galaxy = f.readlines()
    galaxy = pd.DataFrame([x.strip().split() for x in galaxy])
    return galaxy
    


    
def check_2mass_zps(directory,ra_list, dec_list,maxsep=2.0,jcol=3,max_mag=None,mask=None):

    """
    Compares brightest stars in catalog to 2mass to check for offsets in the J band.

    directory: directory with all raw photometry files
    maxsep: maximum seperation between 2mass catalog and merged Fourstar catalog in arcseconds
    jcol: column in text file with J magnitude (default: 3)
    ra_list, dec_list = list of column indices for ra/dec
    mask: list of indices to remove from directory of photometry files. Only relevant if photometry file lacks J-band data
    """

    # read in Merged photometry catalog
    dir_ = glob.glob(directory)
    print(dir_) # double check for ra_list/dec_list ordering

    max_sep = maxsep * u.arcsec
    
    # mask out files that do not have J-band photometry. Note that the indices shift after each deletion. 
    if mask!=None:
        for j in mask:
            del dir_[j]

    # loop over each individual epoch
    for i in range(len(dir_)):
        
        cat = read_file(dir_[i])
        
        racol = ra_list[i]
        deccol = dec_list[i]
        
        # read in 2mass catalog from Vizier using central coordinates from merged catalog
        v = Vizier(columns=["*","+Jmag"], catalog="II/246")
        center = SkyCoord(ra=np.mean(cat[racol].astype(float)), dec= np.mean(cat[deccol].astype(float)), unit=(u.deg, u.deg), frame='icrs' )
        result = v.query_region(center, radius="20m")
        tmass = result[0]

        if max_mag==None:
            max_mag = np.max(tmass['Jmag'])+.5

        # cut catalog so that it only contains brightest stars to make matching easier
        cat = cat[cat[jcol].astype(float)<(max_mag)]
        
        if len(cat)==0:
            print('No matches found')
            break

        # match two catalogs together
        c = SkyCoord(ra=tmass['RAJ2000'].data * u.deg, dec=tmass['DEJ2000'].data * u.deg ,frame='icrs') # convert to Skycoord
        inst = SkyCoord(cat[racol].astype(float)*u.deg,cat[deccol].astype(float)*u.deg,frame='icrs',unit='deg') # convert to Skycoord
    
        idx, d2d, d3d = c.match_to_catalog_sky(inst)
        sep_con = (d2d < max_sep) # only take matches less than 2 arcseconds away from each other
        inst_matches = cat.iloc[idx[sep_con]]
        xdata= tmass['Jmag'][sep_con].data

        # calculate median offset and its error on the mean
        mu = np.median(xdata - inst_matches[jcol].astype(float).values)
        scatter = np.sqrt(np.sum(((xdata - inst_matches[jcol].astype(float).values)-mu)**2))/(len(xdata))
               
        # remove matches that are more than 2-sigma away (e.g., Hatt+17), but only do this for stars with more than 1 match
        if len(xdata)>1:
            dummy = xdata

            xdata = xdata[((xdata - inst_matches[jcol].astype(float).values)<(2*scatter+mu))&((xdata - inst_matches[jcol].astype(float).values)>(mu-2*scatter))]
            inst_matches = inst_matches[((dummy - inst_matches[jcol].astype(float).values)<(2*scatter+mu))&((dummy - inst_matches[jcol].astype(float).values)>(mu-2*scatter))]

        
            # recalculate median offset and error on the mean
            mu = np.median(xdata - inst_matches[jcol].astype(float).values)
            scatter = np.sqrt(np.sum(((xdata - inst_matches[jcol].astype(float).values)-mu)**2))/(len(xdata))
        if (abs(mu)>.5):
            print('No good matches found with less than 2-sigma')
        if len(xdata)==0:
            print('No matches found after cleaning')
        
        # double check that the matches are working
        plt.figure(figsize=(5,5))
        plt.scatter(cat[racol].astype(float), cat[deccol].astype(float),alpha=.5,label='Fourstar')
        plt.scatter(tmass['RAJ2000'],tmass['DEJ2000'],label='2MASS',color='red',alpha=.5)
        plt.scatter(inst_matches[racol].astype(float), inst_matches[deccol].astype(float), facecolor='none',s=100,edgecolor='black',label='Match') # calculated matches
        plt.xlabel('R.A.',fontsize=20)
        plt.ylabel('Dec.',fontsize=20)
        plt.legend()

        # plot photometric offsets plot
        plt.figure(figsize=(7,4))
        plt.scatter(xdata, xdata - inst_matches[jcol].astype(float).values,color='black')
        plt.axhline(mu,color='black',ls='--')
        plt.ylim(-.5,.5)
        plt.axhspan(mu-scatter, mu+scatter,color='grey',alpha=0.4)
        plt.xlabel('2MASS J',fontsize=20)
        plt.ylabel('2MASS J - Fourstar J',fontsize=20)
        print('Median Offset: '+str(np.round(mu,3))+' +/- ' +str(np.round(scatter,3))+' mag')
