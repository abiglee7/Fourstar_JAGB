#!/usr/bin/env python

import numpy as np
import pandas as pd


def custom_hist(input_data, bin_size): # bins data
    start = np.min(input_data) -2
    stop = np.max(input_data) + 1
    n_bins = (stop - start) / bin_size
    bins = np.linspace(start, stop, int(n_bins) + 1)
    bin_centers = (bins + bin_size / 2.)[:-1]
    values = np.histogram(input_data, bins=bins)[0]
    return bin_centers, values
    

def simple_gaussian(x_i, center,sigma):

    return np.exp(-.5*((x_i-center )/sigma)**2)


def GLOESS(x,y, s, order=2):

    output = []  # smoothed output, same length as data
    data = pd.DataFrame({'x':x,'y':y})
    #  window size from smoothing parameter
    data_length = len(x)
    window_size = int(s * data_length)

    #  iterate over rows of data and calculate a smoothed y value at each x point
    for idx, row in data.iterrows():
        
        #  for each window, compute weighted regression coefficients
        points_ahead = min(data_length - idx, int(window_size / 2.))
        points_behind = min(idx, int(window_size / 2.))
        window = data[idx - points_behind: idx + points_ahead]

        #  gaussian weights
        weights = simple_gaussian(window['x'], row['x'], s).values
        
        #  local regression and output
        coeffs = np.polyfit(window['x'], window['y'], deg=order, w=weights)
        output.append(np.polyval(coeffs, row['x']))
        
    return output
