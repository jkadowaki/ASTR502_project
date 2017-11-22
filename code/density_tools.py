# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import


import numpy as np
import matplotlib.pyplot as plt


"""
TODO:
----

1. Function to identify the location of the over-dense regions and
   return it's coordinates.

"""

def density_peaks(dens, sigma_t, xmin=0, xmax=1, ymin=0, ymax=1):
    """
    Function to plot contours on densities:

    Parameters:
    -----------
    dens : numpy array
          A 2D array with the densities.
    sigma_t : int
          Upper value in terms of \sigma to plot contours.
    xmin : int
          Minimum value in the x-direction (default = 0). 
    xmax : int
          Maximum value in the x-direction (default = 1). 
    ymin : int
          Minimum value in the y-direction (default = 0). 
    ymax : int
          Maximum value in the y-direction (default = 1). 

    Returns:
    --------
    Figure with contours.

    TODO:
    -----
    1. Fix axes ranges.
    2. Make label in a for loop.
    3. Check how to proplerly return a figure.

    """
    assert xmin < xmax, "xmax should be greater than xmin"
    assert ymin < ymax, "ymax should be greater than ymin"
    assert sigma_t >0, "sigma_t should be larger than 0"
    assert type(sigma_t) == int, "sigma_t should be a integer"

    # Defining grid
    x = np.linspace(xmin, xmax, shape(dens)[0])
    y = np.linspace(xmin, xmax, shape(dens)[1])
    X, Y = np.meshgrid(x, y)

    # Defining sigma as the standar deviation of the data
    sigma = np.std(dens.flatten())
    # Finding the median of the all the data in the field
    dens_median = np.median(dens.flatten())
    # Defining the contours range. 
    overdensities = []
    for i in range(1, sigma_t+1):
        overdensities.append(dens_median + i*sigma)

    # Creating contour plot:
    fig = plt.figure(figsize=(6, 6))
    plt.contourf(X, Y, dens, overdensities)
    cbar = plt.colorbar()
    cbar.ax.set_yticklabels([r'$1\sigma$',r'$2 \sigma$',r'$3 \sigma$',r'$4 \sigma$', r'$5 \sigma$'])

    return fig

