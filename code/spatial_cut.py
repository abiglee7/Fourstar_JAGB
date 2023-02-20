#!/usr/bin/env python

from astropy.coordinates import ICRS, Distance, Angle, SkyCoord
from astropy import units as u
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS

def compute_rgc(ra,
                dec,
                glx_ctr,
                glx_PA,
                glx_incl,
                glx_dist):
    """
    Computes deprojected radial distance "RGC" of a single object in a galaxy.
    Inspired by P.Barmby (https://gist.github.com/PBarmby/9355846)

    coord: coordinates of given point
    ra/dec: array of all ra/dec of objects
    glx_ctr: galaxy center, list, e.g. [0,0]
    glx_PA: galaxy position angle in degrees,e.g. 38
    glx_incl: galaxy inclination angle. Should be <90
    glx_dist: distance to galaxy
    """


    # Put everything into SkyCoord objects
    coord = SkyCoord(ra*u.degree, dec*u.degree, frame='icrs')
    glx_PA = Angle(str(glx_PA)+'d')
    glx_incl =Angle(str(glx_incl)+'d')
    glx_dist  =Distance(glx_dist, unit=u.Mpc)
    

    sky_radius = glx_ctr.separation(coord)
    avg_dec = 0.5 * (glx_ctr.dec + coord.dec).radian

    x = (glx_ctr.ra - coord.ra) * np.cos(avg_dec)
    y = glx_ctr.dec - coord.dec

    phi = glx_PA - Angle('90d') + Angle(np.arctan(y.arcsec / x.arcsec), unit=u.rad)
    xp = (sky_radius * np.cos(phi.radian)).arcmin
    yp = (sky_radius * np.sin(phi.radian)).arcmin

    ypp = yp / np.cos(glx_incl.radian)
    obj_radius = np.sqrt(xp ** 2 + ypp ** 2)  # in arcmin
    obj_dist = Distance(Angle(obj_radius, unit=u.arcmin).radian * glx_dist,
            unit=glx_dist.unit)

    return obj_dist

def compute_rgc_wholegalaxy(catalog,glx_ctr, glx_PA, glx_incl, glx_dist,testing_phase=None):
    
    """
    Returns RGC (in kpc) for every object in a galaxy catalog.
    
    testing_phase: integer of numbers to test the catalog
    
    """
    
    if testing_phase==None:
        rgc=[]
        for i in range(len(catalog)):

            rgc.append(np.array(compute_rgc(ra=np.array(catalog['ra'])[i], dec=np.array(catalog['dec'])[i],
                               glx_ctr=glx_ctr,
                               glx_PA = glx_PA,
                               glx_incl=glx_incl,
                               glx_dist=glx_dist )
                               ))
        return np.array(rgc)*1000

    else: # test to see if PA/incl is right
        test = np.random.randint(len(catalog), size=testing_phase)
        rgc=[]
        for i in test:

            rgc.append(np.array(compute_rgc(ra=np.array(catalog['ra'])[i], dec=np.array(catalog['dec'])[i],
                               glx_ctr=glx_ctr,
                               glx_PA = glx_PA,
                               glx_incl=glx_incl,
                               glx_dist=glx_dist )
                               ))
    
        return np.array(rgc)*1000, test

