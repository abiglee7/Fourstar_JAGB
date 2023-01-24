#!/usr/bin/python 

from matplotlib import gridspec
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


def clean_photometry(galaxy,
            
                     
                    const_err=None, A_err=None, p_err=None,
                    const_sharp=None, A_sharp=None, p_sharp=None,
                    const_chi=None, p_chi=None,
                     
                    jlower=None, jupper=None):
    """
    Given a Pandas dataframe, clean magnitudes based on error, chi, sharp
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
    plt.figure(figsize=(7,7))

    gs = gridspec.GridSpec(3, 1,height_ratios=[1,1, 1])
    ax0 = plt.subplot(gs[0])
    ax0.set_ylabel('$\sigma$',fontsize=20)
    ax0.set_ylim(0,.2)
    ax0.plot(np.arange(15,23,.1),error_cut(np.arange(15,23,.1)), color='red')



    # error
    bad_error = galaxy.loc[galaxy['J_err'].astype(float)>error_cut(galaxy['J'].astype(float))]
    good_error= galaxy.loc[galaxy['J_err'].astype(float)<error_cut(galaxy['J'].astype(float))]
    ax0.scatter(bad_error['J'].astype(float),bad_error['J_err'].astype(float),alpha=.3,color='grey',s=.5)
    ax0.scatter(good_error['J'].astype(float),good_error['J_err'].astype(float),alpha=1,color='black',s=.5)


    ax1 = plt.subplot(gs[1])

    ax1.plot(np.arange(15,23,.1), sharp_cut(np.arange(15,23,.1)),color='red')
    ax1.plot(np.arange(15,23,.1), -sharp_cut(np.arange(15,23,.1)),color='red')

    bad_sharp_1 = galaxy.loc[(galaxy['sharp'].astype(float)<-sharp_cut(galaxy['J'].astype(float)))]
    bad_sharp_2 = galaxy.loc[(galaxy['sharp'].astype(float)>sharp_cut(galaxy['J'].astype(float)))]
    good_sharp = galaxy.loc[(galaxy['sharp'].astype(float)<sharp_cut(galaxy['J'].astype(float)))&(galaxy['sharp'].astype(float)>-sharp_cut(galaxy['J'].astype(float)))]
    ax1.scatter(bad_sharp_1['J'].astype(float),bad_sharp_1['sharp'].astype(float),alpha=.3,color='grey',s=.5)
    ax1.scatter(bad_sharp_2['J'].astype(float),bad_sharp_2['sharp'].astype(float),alpha=.3,color='grey',s=.5)
    ax1.scatter(good_sharp['J'].astype(float),good_sharp['sharp'].astype(float),alpha=1,color='black',s=.5)


    ax1.set_ylim(-.9,.9)
    ax1.set_ylabel('Sharp',fontsize=20)


    ax1.set_ylim(-2,2)

    ax2 = plt.subplot(gs[2])
    ax2.set_ylim(-.04,2)
    ax2.set_ylabel('$\chi$',fontsize=20)
    ax2.plot(np.arange(15,23,.1), chi_cut(np.arange(15,23,.1)),color='red')

    bad_chi = galaxy.loc[galaxy['chi'].astype(float)>chi_cut(galaxy['J'].astype(float))]
    good_chi= galaxy.loc[galaxy['chi'].astype(float)<chi_cut(galaxy['J'].astype(float))]
    ax2.scatter(bad_chi['J'].astype(float),bad_chi['chi'].astype(float),alpha=.3,color='grey',s=.5)
    ax2.scatter(good_chi['J'].astype(float),good_chi['chi'].astype(float),alpha=1,color='black',s=.5)



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
    
        
    galaxy = galaxy.loc[(galaxy['J_err'].astype(float)<error_cut(galaxy['J'].astype(float)))&
                       (galaxy['sharp'].astype(float)<sharp_cut(galaxy['J'].astype(float)))&
                       (galaxy['sharp'].astype(float)>-sharp_cut(galaxy['J'].astype(float)))&
                       (galaxy['chi'].astype(float)<chi_cut(galaxy['J'].astype(float)))]
        
    
    # plot a CMD to check how well it was cleaned:
    plt.figure(figsize=(4,5))
    plt.scatter(galaxy['J'].astype(float)-galaxy['K'].astype(float), galaxy['J'].astype(float),s=1,color='black')
    plt.xlim(0,1.3)
    plt.ylim(20.5,14)
    plt.axvline(1.4)
    plt.axvline(2)


    return galaxy
