# -*- coding: utf-8 -*-
from __future__ import print_function, absolute_import

import numpy as np
from astroML.density_estimation import KDE, KNeighborsDensity
from sklearn.neighbors import KernelDensity



def NN_bayesian_density(x, y, NN, grid_size):
    """
    Function to compute the density of a distribution of particles
    using the K-Nearest Neighboors
    method from:
    http://www.astroml.org/modules/generated/astroML.density_estimation.KNeighborsDensity.html#astroML.density_estimation.KNeighborsDensity
    See Ivezic 10? for the details on how the algorithm works.

    Input:
    ------
    x : 1D numpy.array
        Array with the x-coordinates of the data.
    y : 1D numpy.array
        Array with the y-coordinates of the data.
    NN : int
        Number of neighboors to compute the desnity.
    grid_size : int
        Grid size in which the density is going to be evaluated.

    """
    assert len(x)==len(y), "Input data have different size"
    assert type(NN) == int, "NN should be of type int"
    assert type(grid_size) == int, "grid_zise should be of type int"

    # Grid parameters
    Nx = grid_size
    Ny = grid_size
    xmin, xmax = (min(x), max(x))
    ymin, ymax = (min(y), max(y))

    # Making a grid
    Xgrid = np.vstack(map(np.ravel, np.meshgrid(np.linspace(xmin, xmax, Nx),\
                                                np.linspace(ymin, ymax, Ny)))).T
    # Putting data in 2d-array
    X = np.array([x, y]).T

    # Computing the density
    knn = KNeighborsDensity('bayesian', NN)
    dens_KNN = knn.fit(X).eval(Xgrid).reshape((Ny, Nx))

    return dens_KNN




def KDE(ra, dec, kernel = 'gaussian'):
    """
    Compute the Kernel Density Estimators (KDE) from scikit learn

    Input:
    -----

    x : 1D numpy.array
        Array with the x-coordinates of the data.
    y : 1D numpy.array
        Array with the y-coordinates of the data.
    kernel : str
        Optional kernel: gaussian (default), tophat, exponential

    Output:
    -------

    dens : numpy.array 
        2d array with the density estimation

    """


    ramin = min(ra)
    ramax = max(ra)
    decmin = min(dec)
    decmax = max(dec)

    x, y = np.mgrid[ramin:ramax:300j, decmin:decmax:300j]
    Xgrid = np.vstack([x.ravel(), y.ravel()]).T

    values = np.vstack([ra, dec])
    X = values.T

    kde1 = KernelDensity(0.01, kernel=kernel)
    log_dens1 = kde1.fit(X).score_samples(Xgrid)
    dens1 = X.shape[0] * np.exp(log_dens1).reshape((300,300))

    return (dens1.T)


