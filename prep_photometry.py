"""
Author: Abigail Lee

This script reads in raw .txt photometry files, cleans the J band magnitudes based on their chi, sharp, and error values, and outputs a cleaned .csv file for each epoch.
"""

#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
from matplotlib import gridspec

def read_file(filename):

    """
    Reads in a photometry .txt file and turns it into a pandas dataframe for easier use

    filename: .txt file
    """

    with open(filename) as f:
        galaxy = f.readlines()
    galaxy = pd.DataFrame([x.strip().split() for x in galaxy])
    return galaxy
    


def clean_photometry(galaxy, jcol, secondcol=None, secondcolname=None,
                    cleanerr=False, cleanchi=False, cleansharp=False,
                    jcolerr=None, jcolsharp=None, jcolchi=None,
                    const_err=None, A_err=None, p_err=None,
                    const_sharp=None, A_sharp=None, p_sharp=None,
                    const_chi=None, p_chi=None,
                    plot=True,
                    jlower=None, jupper=None,save_fig=False,name=None, date=None):
    """
    Given a Pandas dataframe, clean J-band magnitudes based on error, chi, sharp
    
    galaxy: pandas dataframe
    jcol, jcolerr, jcolsharp, jcolchi: column index (e.g., 1)
    secondcol, secolname: index and name of secondcol (H or K) to make a test CMD
    const_err, A_err, p_err: parameters for cleaning functions (see paper for exact functional form)
    cleanerr, cleanchi, cleansharp: boolean, whether or not you want to clean on this parameter
    plot: do you want to plot it or not?
    jlower, jupper = limits for magnitude plots
    """
        
    def error_cut(mag):
        """
        Cleans photometry based on error value.
        """
        return const_err+ A_err*np.exp(mag-p_err)

    def sharp_cut(mag):
        """
        Cleans photometry based on sharp value.
        """
        return const_sharp + A_sharp*np.exp(mag-p_sharp)

    def chi_cut(mag):
        """
        Cleans photometry based on chi value.
        """
        return const_chi+np.exp(-mag+p_chi)
    
    
    # plot photometric quality cuts
    if plot==True:
        plt.figure(figsize=(7,7))
        
        gs = gridspec.GridSpec(3, 1,height_ratios=[1,1, 1])
        ax0 = plt.subplot(gs[0])
        ax0.set_ylabel('$\sigma$',fontsize=20)
        ax0.set_ylim(0,.2)
        ax0.plot(np.arange(15,23,.1),error_cut(np.arange(15,23,.1)), color='red')
        
        
        
        # error
        bad_error = galaxy.loc[galaxy[jcolerr].astype(float)>error_cut(galaxy[jcol].astype(float))]
        good_error= galaxy.loc[galaxy[jcolerr].astype(float)<error_cut(galaxy[jcol].astype(float))]
        ax0.scatter(bad_error[jcol].astype(float),bad_error[jcolerr].astype(float),alpha=.3,color='grey',s=.5)
        ax0.scatter(good_error[jcol].astype(float),good_error[jcolerr].astype(float),alpha=1,color='black',s=.5)
        
        
        ax1 = plt.subplot(gs[1])
        
        ax1.plot(np.arange(15,23,.1), sharp_cut(np.arange(15,23,.1)),color='red')
        ax1.plot(np.arange(15,23,.1), -sharp_cut(np.arange(15,23,.1)),color='red')
        
        bad_sharp_1 = galaxy.loc[(galaxy[jcolsharp].astype(float)<-sharp_cut(galaxy[jcol].astype(float)))]
        bad_sharp_2 = galaxy.loc[(galaxy[jcolsharp].astype(float)>sharp_cut(galaxy[jcol].astype(float)))]
        good_sharp = galaxy.loc[(galaxy[jcolsharp].astype(float)<sharp_cut(galaxy[jcol].astype(float)))&(galaxy[jcolsharp].astype(float)>-sharp_cut(galaxy[jcol].astype(float)))]



        ax1.scatter(bad_sharp_1[jcol].astype(float),bad_sharp_1[jcolsharp].astype(float),alpha=.3,color='grey',s=.5)
        ax1.scatter(bad_sharp_2[jcol].astype(float),bad_sharp_2[jcolsharp].astype(float),alpha=.3,color='grey',s=.5)
        ax1.scatter(good_sharp[jcol].astype(float),good_sharp[jcolsharp].astype(float),alpha=1,color='black',s=.5)


        ax1.set_ylim(-.9,.9)
        ax1.set_ylabel('Sharp',fontsize=20)

        
        ax1.set_ylim(-2,2)
        
        ax2 = plt.subplot(gs[2])
        ax2.set_ylim(-.04,.3)
        ax2.set_ylabel('$\chi$',fontsize=20)
        ax2.plot(np.arange(15,23,.1), chi_cut(np.arange(15,23,.1)),color='red')

        bad_chi = galaxy.loc[galaxy[jcolchi].astype(float)>chi_cut(galaxy[jcol].astype(float))]
        good_chi= galaxy.loc[galaxy[jcolchi].astype(float)<chi_cut(galaxy[jcol].astype(float))]
        ax2.scatter(bad_chi[jcol].astype(float),bad_chi[jcolchi].astype(float),alpha=.3,color='grey',s=.5)
        ax2.scatter(good_chi[jcol].astype(float),good_chi[jcolchi].astype(float),alpha=1,color='black',s=.5)

        

        ax0.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both maKor and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False)

        ax1.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both maKor and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False)
        
        ax0.set_xlim(jlower,jupper)
        ax1.set_xlim(jlower,jupper)
        ax2.set_xlim(jlower,jupper)
        
        plt.xlabel('J',fontsize=20)
        plt.subplots_adjust(hspace=0)
        ax0.set_title(str(name)+': '+str(date),fontsize=20)
        
        if save_fig==True:
            plt.savefig('/Users/abigaillee/Documents/Research/Fourstar_images/'+str(name)+str(date)+'.png')
    
    if cleanerr==True:
        galaxy = galaxy.loc[galaxy[jcolerr].astype(float)<error_cut(galaxy[jcol].astype(float))]
        
    
    if cleansharp==True:
        galaxy = galaxy.loc[galaxy[jcolsharp].astype(float)<sharp_cut(galaxy[jcol].astype(float))]
        galaxy = galaxy.loc[galaxy[jcolsharp].astype(float)>-sharp_cut(galaxy[jcol].astype(float))]
     
    if cleanchi==True:
        galaxy = galaxy.loc[galaxy[jcolchi].astype(float)<chi_cut(galaxy[jcol].astype(float))]
    
    # plot a CMD to check how well it was cleaned:
    if secondcolname=='H':
        plt.figure(figsize=(4,5))
        plt.scatter(galaxy[jcol].astype(float)-galaxy[secondcol].astype(float), galaxy[jcol].astype(float),s=1,color='black')
        plt.xlim(0.5,1.5)
        plt.ylim(21,16)
        plt.axvline(1.1)
        plt.axvline(1.35)
    if secondcolname=='K':
        plt.figure(figsize=(4,5))
        plt.scatter(galaxy[jcol].astype(float)-galaxy[secondcol].astype(float), galaxy[jcol].astype(float),s=1,color='black')
        plt.xlim(0.5,2.5)
        plt.ylim(21,16)
        plt.axvline(1.4)
        plt.axvline(2)
    

    return galaxy

def combine(jcol, cleanerr=False, cleanchi=False, cleansharp=False,
                    jcolerr=None, jcolsharp=None, jcolchi=None,
                    const_err=None, A_err=None, p_err=None,
                    const_sharp=None, A_sharp=None, p_sharp=None,
                    const_chi=None, p_chi=None,
                    plot=True,
                    jlower=None, jupper=None,save_fig=False,name=None, folder_name=None, date=None, secondcol=None, secondcolname=None):
    """
    Reads in .txt file and then outputs cleaned photometry pandas dataframe
    """

    folder = glob.glob('/Users/abigaillee/Photometry/JAGB stars MISC/Fourstar/Raw/'+str(folder_name)+'/*')
    for i in range(len(folder)):
        file = read_file(folder[i])
        
        file = clean_photometry(file, jcol=jcol,
            jcolsharp=jcolsharp,jcolerr=jcolerr,jcolchi=jcolchi,
            cleansharp=cleansharp, cleanerr=cleanerr,cleanchi=cleanchi,
            const_err=const_err, A_err=A_err, p_err=p_err,
            const_sharp=const_sharp, A_sharp = A_sharp, p_sharp = p_sharp,
            const_chi=const_chi, p_chi=p_chi,
            plot=plot,
            jlower=jlower, jupper=jupper,save_fig=save_fig,
            name=name, date=date, secondcol=secondcol, secondcolname = secondcolname)
            
        file.to_csv('/Users/abigaillee/Photometry/JAGB stars MISC/Fourstar/Cleaned/'+str(name)+'/'+str(name)+'_'+str(date)+'.csv')
