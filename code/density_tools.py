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

def od_coordinates(dens, sigma_t, xmin=0, xmax=1, ymin=0, ymax=1):
    assert xmin < xmax, "xmax should be greater than xmin"
    assert ymin < ymax, "ymax should be greater than ymin"
    assert sigma_max >0, "sigma_t should be larger than 0"
    assert type(sigma_max) == int, "sigma_t should be a integer"

    # Defining grid
    x = np.linspace(xmin, xmax, np.shape(dens)[0])
    y = np.linspace(ymin, ymax, np.shape(dens)[1])
    dx = x[1]-x[0]
    dy = y[1]-y[0]
    # Defining sigma as the standard deviation of the data
    sigma = np.std(dens.flatten())
    # Finding the median of the all the data in the field
    dens_median = np.median(dens.flatten())
    # Defining the contours range.·
    overdensities = []
    color_bar_labels = []

    index = np.where(dens>(dens_median+sigma*sigma_t))

    x_coord = index[0]*dx
    y_coord = index[1]*dy

    x_mean = np.sum(x_coord)/len(x_coord)
    y_mean = np.sum(y_coord)/len(y_coord)

    return x_mean, y_mean

def density_peaks(dens, sigma_min, sigma_max, ax=None,  xmin=0, xmax=1, ymin=0, ymax=1):
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
    1. Allow to use more than 10 ticklabels in the colorbar

    """
    assert xmin < xmax, "xmax should be greater than xmin"
    assert ymin < ymax, "ymax should be greater than ymin"
    assert sigma_max >0, "sigma_t should be larger than 0"
    assert type(sigma_max) == int, "sigma_t should be a integer"

    # Defining grid
    x = np.linspace(xmin, xmax, np.shape(dens)[0])
    y = np.linspace(ymin, ymax, np.shape(dens)[1])
    X, Y = np.meshgrid(x, y)

    # Defining sigma as the standard deviation of the data
    sigma = np.std(dens.flatten())
    # Finding the median of the all the data in the field
    dens_median = np.median(dens.flatten())
    # Defining the contours range. 
    overdensities = []
    color_bar_labels = []

    for i in range(sigma_min, sigma_max+1):
        overdensities.append(dens_median + i*sigma)
        color_bar_labels.append(str(i) + '$\sigma$')
    # Creating contour plot:
    fig = plt.figure(figsize=(6, 6))
    plt.contourf(X, Y, dens, overdensities)
    cbar = plt.colorbar()
    cbar.ax.set_yticks(np.arange(0, len(color_bar_labels)+1, 1))
    cbar.ax.set_yticklabels(color_bar_labels)
    #fig.clf()
    if ax is None:
        ax = plt.gca()
    bp = ax.boxplot(more_data)
    return bp
    #return fig

