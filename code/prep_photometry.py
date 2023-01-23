#!/usr/bin/python 

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def read_file(filename):

    """
    Reads in a photometry .txt file and turns it into a pandas dataframe for easier use
    
    jcol: column index that contains J-band photometry (e.g. 3)
    second_col: column index that contains either J-band or H-band photometry
    second_col_name: Indicate whether second magnitude is H or K, str
    
    """
    
    with open(filename) as f:
        galaxy = f.readlines()
    galaxy = pd.DataFrame([x.strip().split() for x in galaxy])
    
    galaxy.columns = ['ID','X', 'Y','Mag','Error', 'Ext_err','Num','chi','sharp','var','blunder']
    galaxy = galaxy.iloc[3:]


    
    return galaxy
