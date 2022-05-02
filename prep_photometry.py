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
        
def merge_photometry_tmp(catalog, racol, deccol):

    """
    Unfortunately, TOPCAT only allows for <=5 matches.
    So if there are 6+ files, I need to do two seperate matches.
    This function reads in initial 5 matches, averages their ra/dec,
    and then matches this file again to export to TOPCAT.

    This function should only be used if there are 5+ epochs for a given galaxy
    (i.e., IC 1613, NGC 300, NGC 3109, M83, Cen A, NGC 247, NGC 6822)

    catalog location: directory with catalog output from TOPCAT (should be 5 matches)
    racol/deccol: list of indices for ra/dec for each catalog. Note: this will correspond to the order in which you input the files into TOPCAT.

    """

    match = pd.read_csv(str(catalog)+'matchJHK.csv')
    match = match.iloc[1:] # delete first line because TOPCAT puts in an unncessary header

    match['ra'] = np.nanmean([ match[col] for col in racol], axis=0)
    match['dec'] = np.nanmean([ match[col] for col in deccol], axis=0)

    match.to_csv(str(catalog)+'matchJHK_test.csv')

def merge_photometry(catalog, jcol, jcolerr, n_epochs,racol, deccol, galaxyname,hcol=None,
                hcolerr=None, kcol=None, kcolerr=None):
    """
    catalog location: directory with catalog output from TOPCAT (should be 5 matches)
    racol/deccol: list of indices for ra/dec for each catalog. Note: this will correspond to the order in which you input the files into TOPCAT.
    jcol, hcol, kcol, etc.: lists of indices for magnitudes
    n_epochs: total number of epochs
    galaxyname: str for naming the final merged catalog
    """

    match = pd.read_csv(str(catalog)+'finalmatch.csv')
    match = match.iloc[1:] # delete first line because TOPCAT puts in an unncessary header
    print('Total Length of Original Catalog: '+str(len(match)))

    # Check that matches were okay. None of these should exceed ~2 mag.
    print('Max Photometeric Offsets:')
    for i in range(1, n_epochs):
        print(np.max(match[jcol[0]]-match[jcol[i]]))

    # remove stars that don't have K (or H) band data to speed up calculations below (there is no point in having single-band photometry of a given star)
    match['dummy'] = np.nanmean([ match[col] for col in kcol], axis=0) # creates a dummy variable that shows where nans are
    match = match[match['dummy'].notna()]
    print('New Length of Original Catalog: '+str(len(match)))


    # create a new dataframe for FINAL clean, merged catalog
    cat = pd.DataFrame({})
    cat['ra'] = np.nanmean([ match[col] for col in racol], axis=0)
    cat['dec'] = np.nanmean([ match[col] for col in deccol], axis=0)

    # calculate weighted average photometry
    J = []; J_err=[]; H = []; H_err=[]; K = []; K_err = []

    for i in range(len(match)):
        
        # J  band
        a_J = [ match[col].iloc[i] for col in jcol]
        weights_J = [ match[col].iloc[i] for col in jcolerr]
        # mask out the dates that don't have values
        ma_J = np.ma.MaskedArray(a_J, mask=np.isnan(a_J)).compressed()
        maw_J = np.ma.MaskedArray(weights_J, mask=np.isnan(weights_J)).compressed()
        # then average
        J_i = np.ma.average(ma_J, weights=1/np.array(maw_J**2))
        uncert_J = np.sqrt(np.sum(np.array(maw_J)**2))/len(maw_J)
        J.append(J_i)
        J_err.append(uncert_J)

        # H band
        a_H = [ match[col].iloc[i] for col in hcol]
        weights_H = [ match[col].iloc[i] for col in hcolerr]
        ma_H = np.ma.MaskedArray(a_H, mask=np.isnan(a_H)).compressed()
        maw_H = np.ma.MaskedArray(weights_H, mask=np.isnan(weights_H)).compressed()
        H_i = np.ma.average(ma_H, weights=1/np.array(maw_H**2))
        uncert_H = np.sqrt(np.sum(np.array(maw_H)**2))/len(maw_H)
        H.append(H_i)
        H_err.append(uncert_H)

        # K band
        a_K = [ match[col].iloc[i] for col in kcol]
        weights_K = [ match[col].iloc[i] for col in kcolerr]
        ma_K = np.ma.MaskedArray(a_K, mask=np.isnan(a_K)).compressed()
        maw_K = np.ma.MaskedArray(weights_K, mask=np.isnan(weights_K)).compressed()
        K_i = np.ma.average(ma_K, weights=1/np.array(maw_K**2))
        uncert_K = np.sqrt(np.sum(np.array(maw_K)**2))/len(maw_K)
        K.append(K_i)
        K_err.append(uncert_K)


    J = np.array(J)
    H = np.array(H)
    K = np.array(K)

    cat['J']=J
    cat['H']=H
    cat['K']=K
    cat['J_err'] = J_err
    cat['H_err'] = H_err
    cat['K_err'] = K_err

    cat.to_csv('/Users/abigaillee/Photometry/JAGB stars MISC/Fourstar/Merged/'+str(galaxyname)+'.csv')

