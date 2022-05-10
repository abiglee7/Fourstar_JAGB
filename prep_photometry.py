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
import math

def read_file(filename):

    """
    Reads in a photometry .txt file and turns it into a pandas dataframe for easier use

    filename: .txt file
    """

    with open(filename) as f:
        galaxy = f.readlines()
    galaxy = pd.DataFrame([x.strip().split() for x in galaxy])
    return galaxy
    


    
        
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


    
def check_2mass_zps(directory,ra_list, dec_list,jcol,hcol, kcol, maxsep=2.0,max_mag=None):

    """
    Compares brightest stars in catalog to 2mass to check for offsets in the J band.

    directory: directory with all raw photometry files
    maxsep: maximum seperation between 2mass catalog and merged Fourstar catalog in arcseconds
    jcol: column in text file with J magnitude (default: 3)
    ra_list, dec_list = list of column indices for ra/dec

    """

    # read in photometry catalog
    dir_ = glob.glob(directory)
    
    max_sep = maxsep * u.arcsec
    

    # loop over each individual epoch
    tmass_mags = ['Jmag','Hmag','Kmag']
    J_offsets = []; H_offsets=[]; K_offsets=[]
    for i in range(len(dir_)):
        
        print(dir_[i]) # double check for ra_list/dec_list ordering
        
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
        cat = cat[cat[jcol[i]].astype(float)<(max_mag)]
        
        if len(cat)==0:
            print('No matches found')
            break

        # match two catalogs together
        c = SkyCoord(ra=tmass['RAJ2000'].data * u.deg, dec=tmass['DEJ2000'].data * u.deg ,frame='icrs') # convert to Skycoord
        inst = SkyCoord(cat[racol].astype(float)*u.deg,cat[deccol].astype(float)*u.deg,frame='icrs',unit='deg') # convert to Skycoord
    
        idx, d2d, d3d = c.match_to_catalog_sky(inst)
        sep_con = (d2d < max_sep) # only take matches less than 2 arcseconds away from each other
        inst_matches = cat.iloc[idx[sep_con]]
        
        
        
        # double check that the matches are working
        plt.figure(figsize=(5,5))
        plt.scatter(cat[racol].astype(float), cat[deccol].astype(float),alpha=.5,label='Fourstar')
        plt.scatter(tmass['RAJ2000'],tmass['DEJ2000'],label='2MASS',color='red',alpha=.5)
        plt.scatter(inst_matches[racol].astype(float), inst_matches[deccol].astype(float), facecolor='none',s=100,edgecolor='black',label='Match') # calculated matches
        plt.xlabel('R.A.',fontsize=20)
        plt.ylabel('Dec.',fontsize=20)
        plt.legend()
        plt.show()

        
        for count, band in enumerate([jcol[i], hcol[i],kcol[i]]): # loop over JHK
            
            if count==0:
                print('J')
                tmass_mag = 'Jmag'
            elif count==1:
                print('H')
                tmass_mag = 'Hmag'
            elif count==2:
                print('K')
                tmass_mag = 'Kmag'
            
            if math.isnan(band)==True: # if no X (jhk) band data for this epoch
                print('No data in this band')
                if count==0:
                    J_offsets.append(0)
                elif count==1:
                    H_offsets.append(0)
                elif count==2:
                    K_offsets.append(0)
            elif math.isnan(band)==False: # if X (jhk) band data for this epoch

                xdata= tmass[tmass_mag][sep_con].data
                inst_matches_tmp = inst_matches # start afresh
                
                # calculate median offset and its error on the mean
                mu = np.median(xdata - inst_matches_tmp[band].astype(float).values)
                scatter = np.sqrt(np.sum(((xdata - inst_matches_tmp[band].astype(float).values)-mu)**2))/np.sqrt(len(xdata))

                # remove matches that are more than 2-sigma away (e.g., Hatt+17), but only do this for stars with more than 1 match
                if len(xdata)>1:
                    dummy = xdata

                    xdata = xdata[((xdata - inst_matches_tmp[band].astype(float).values)<(2*scatter+mu))&((xdata - inst_matches_tmp[band].astype(float).values)>(mu-2*scatter))&(xdata-inst_matches_tmp[band].astype(float).values>-.4)&(xdata-inst_matches_tmp[band].astype(float).values<.4)]
                    inst_matches_tmp = inst_matches_tmp[((dummy - inst_matches_tmp[band].astype(float).values)<(2*scatter+mu))&((dummy - inst_matches_tmp[band].astype(float).values)>(mu-2*scatter))&(dummy-inst_matches_tmp[band].astype(float).values>-.4)&(dummy-inst_matches_tmp[band].astype(float).values<.4)]


                    # recalculate median offset and error on the mean
                    mu = np.median(xdata - inst_matches_tmp[band].astype(float).values)
                    scatter = np.sqrt(np.sum(((xdata - inst_matches_tmp[band].astype(float).values)-mu)**2))/np.sqrt(len(xdata))
                if (abs(mu)>.5):
                    print('No good matches found with less than 2-sigma')
                if len(xdata)==0:
                    print('No matches found after cleaning')
                    
                # we will not use the photometric correction if there are fewer than 10 stars
                if count==0:
                    if len(inst_matches_tmp)<10:
                        J_offsets.append(np.nan)
                    elif len(inst_matches_tmp)>=10:
                        J_offsets.append(mu)
                elif count==1:
                    if len(inst_matches_tmp)<10:
                        H_offsets.append(np.nan)
                    elif len(inst_matches_tmp)>=10:
                        H_offsets.append(mu)
                elif count==2:
                    if len(inst_matches_tmp)<10:
                        K_offsets.append(np.nan)
                    elif len(inst_matches_tmp)>=10:
                        K_offsets.append(mu)
                
                print(len(inst_matches_tmp))
                
                # plot photometric offsets plot
                plt.figure(figsize=(7,4))
                plt.scatter(xdata, xdata - inst_matches_tmp[band].astype(float).values,color='black')
                plt.axhline(mu,color='black',ls='--')
                plt.ylim(-.5,.5)
                plt.axhspan(mu-scatter, mu+scatter,color='grey',alpha=0.4)
                plt.xlabel('2MASS',fontsize=20)
                plt.ylabel('2MASS - Fourstar',fontsize=20)
                plt.show()
                print('Median Offset: '+str(np.round(mu,3))+' +/- ' +str(np.round(scatter/np.sqrt(len(xdata)),3))+' mag')

    # write all offsets to list. Note that 0 means there is no data in that band.
    offset_df = pd.DataFrame({'Date':[date[-10:] for date in dir_],
                              'J':J_offsets,'H':H_offsets,'K':K_offsets})
    
    
    return offset_df

def check_2mass_zps_sparse(fiducial,directory,ra_list, dec_list,jcol,hcol, kcol, non_epochs,df,
                       maxsep=2.0,max_mag=18, fiducial_cols=[], min_mag=13):

    """
    Check the 2MASS ZPs for epochs that fewer than 10 matches with 2MASS. Here, we will match a given epoch with the epoch with the most 2MASS stars.

    fiducial: dataframe of "best epoch" that will be matched against the other epochs
    directory: directory with all raw photometry files
    maxsep: maximum seperation between 2mass catalog and merged Fourstar catalog in arcseconds
    jcol: column in text file with J magnitude (default: 3)
    ra_list, dec_list: list of column indices for ra/dec
    non_epochs: list of epochs that DON'T need correcting
    df: dataframe with offsets from 2MASS
    fiducial_cols: indexes with columns for fiducial dataframe (ra, dec, j, h, k)
    min_mag prevents using saturated stars!
    """

    fiducial = read_file(fiducial)
    fiducial = fiducial[(fiducial[fiducial_cols[2]].astype(float)<max_mag)&(fiducial[fiducial_cols[2]].astype(float)>min_mag)] # make it easier to only match bright unsaturated stars
    print('No. of bright stars = '+str(len(fiducial)))



    catalog = SkyCoord(ra=fiducial[fiducial_cols[0]].astype(float)*u.degree, dec=fiducial[fiducial_cols[1]].astype(float)*u.degree)


    max_sep = maxsep * u.arcsec

    # remove epochs that already have good matches, they are different lengths sometimes so had to do two conditions
    dir_ = glob.glob(directory)
    dir_ = [epoch for epoch in dir_ if str(epoch[-10:]) not in non_epochs]
    dir_ = [epoch for epoch in dir_ if str(epoch[-9:]) not in non_epochs]


    # loop over each individual epoch
    tmass_mags = ['Jmag','Hmag','Kmag']
    J_offsets = []; H_offsets=[]; K_offsets=[]


    for i in range(len(dir_)): # loop over all epoochs
        
        print(dir_[i])
        
        cat = read_file(dir_[i])
        cat = cat[cat[jcol[i]].astype(float)<max_mag]
        
        racol = ra_list[i]
        deccol = dec_list[i]
        
        c = SkyCoord(ra=np.array(cat[racol]).astype(float)*u.degree, dec=np.array(cat[deccol].astype(float))*u.degree) # given epoch
        idx, d2d, d3d = c.match_to_catalog_sky(catalog) # match two catalogs
        sep = (d2d < max_sep) # only use matches within this matching radius
        
        c_matches = c[sep]
        catalog_matches = catalog[idx[sep]]
        

        
        # double check matches are working
        plt.figure(figsize=(5,5))
        plt.scatter(cat[racol].astype(float), cat[deccol].astype(float),alpha=.5,label='Given epoch i')
        plt.scatter(np.array(fiducial[fiducial_cols[0]].astype(float)), np.array(fiducial[fiducial_cols[1]].astype(float)),alpha=.5,label='Fiducial')
        plt.scatter(c_matches.ra.degree, c_matches.dec.degree, facecolor='none',s=100,edgecolor='black',label='Match') # calculated matches
        
        plt.xlabel('R.A.',fontsize=20)
        plt.ylabel('Dec.',fontsize=20)
        plt.legend()
        plt.show()
        
        # loop over each band:
        for count, band in enumerate([jcol[i], hcol[i],kcol[i]]):
            
            if math.isnan(band)==True:
                if count==0:
                    J_offsets.append(np.nan)
                elif count==1:
                    H_offsets.append(np.nan)
                elif count==2:
                    K_offsets.append(np.nan)
                
            elif math.isnan(band)==False:

                # calculate median offset and its error on the mean
                mu = np.median(fiducial.iloc[idx[sep]][band].astype(float) - np.array(cat.iloc[sep][band].astype(float)))
                scatter = np.sqrt(np.sum(((np.array(fiducial.iloc[idx[sep]][band].astype(float))- np.array(cat.iloc[sep][band].astype(float)))**2)))/np.sqrt(len(fiducial.iloc[idx[sep]][band].astype(float)))

                fiducial_data = fiducial.iloc[idx[sep]][band].astype(float)
                epoch_data = np.array(cat.iloc[sep][band].astype(float))
                
                # remove matches greater than 2-sigma from mean. Also there are sometimes randomly saturated stars
                dummy = fiducial_data
                fiducial_data = fiducial_data[((fiducial_data-epoch_data)<(2*scatter+mu))&((fiducial_data-epoch_data)>(-2*scatter+mu))&(fiducial_data-epoch_data>-.4)]
                epoch_data = epoch_data[((dummy-epoch_data)<(2*scatter+mu))&((dummy-epoch_data)>(-2*scatter+mu))&(dummy-epoch_data>-.4)]

                # recalculate median offset and error on the mean
                mu = np.median(fiducial_data-epoch_data)
                scatter = np.sqrt(np.sum(((fiducial_data-epoch_data)-mu)**2))/np.sqrt(len(epoch_data))

                plt.figure(figsize=(6,4))
                plt.scatter( fiducial_data,  fiducial_data- epoch_data,color='black')
                plt.ylim(-.5,.5)
                plt.xlabel('Fiducial',fontsize=12)
                plt.ylabel('Fiducial - Epoch',fontsize=12)
                plt.axhline(mu,color='black',ls='--')
                plt.axhspan(mu-scatter, mu+scatter,color='grey',alpha=0.4)
                plt.show()

                print('Median Offset: '+str(np.round(mu,3))+' +/- ' +str(np.round(scatter/np.sqrt(len(fiducial_data)),3))+' mag')
                if count==0:
                    J_offsets.append(mu)
                elif count==1:
                    H_offsets.append(mu)
                elif count==2:
                    K_offsets.append(mu)
    offset_df = pd.DataFrame({'Date':[date[-10:] for date in dir_], 'J':J_offsets, 'H':H_offsets,'K':K_offsets})
    return offset_df
            
