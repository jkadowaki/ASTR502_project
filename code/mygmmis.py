
import pygmmis
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as lines
import matplotlib.cm
import datetime
from functools import partial
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astroML.decorators import pickle_results
from astroML.plotting.tools import draw_ellipse
import scipy as scipy

from astroML.density_estimation import XDGMM
hdulist = fits.open('../data/sdss_fields/draco.fits')
hdulist = fits.open('../data/mock_fields/test_field1.fits')


tbdata = hdulist[1].data
cols = hdulist[1].columns
ra=tbdata['ra']
dec=tbdata['dec']
print "column available", cols
# plt.scatter(tbdata['ra'],tbdata['dec'],s=0.01)
meanra=np.mean(ra)
meandec=np.mean(dec)
rar=max(ra)-meanra
decr=max(dec-meandec)
def plotResults(orig, data, gmm, patch=None, description=None, disp=None):
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111, aspect='equal')

    # plot inner and outer points
    ax.plot(orig[:,0], orig[:,1], 'o', mfc='None', mec='r', mew=1)
    missing = np.isnan(data)
    if missing.any():
        data_ = data.copy()
        data_[missing] = -5 # put at limits of plotting range
    else:
        data_ = data
    ax.plot(data_[:,0], data_[:,1], 's', mfc='b', mec='None')#, mew=1)

    # prediction
    B = 100
    x,y = np.meshgrid(np.linspace(-5,15,B), np.linspace(-5,15,B))
    coords = np.dstack((x.flatten(), y.flatten()))[0]

    # compute sum_k(p_k(x)) for all x
    p = gmm(coords).reshape((B,B))
    # for better visibility use arcshinh stretch
    p = np.arcsinh(p/1e-4)
    cs = ax.contourf(p, 10, extent=(-5,15,-5,15), cmap=plt.cm.Greys)
    for c in cs.collections:
        c.set_edgecolor(c.get_facecolor())

    # plot boundary
    if patch is not None:
        import copy
        if hasattr(patch, '__iter__'):
            for p in patch:
                ax.add_artist(copy.copy(p))
        else:
            ax.add_artist(copy.copy(patch))

    # add description and complete data logL to plot
    logL = gmm(orig, as_log=True).mean()
    if description is not None:
        ax.text(0.05, 0.95, r'%s' % description, ha='left', va='top', transform=ax.transAxes, fontsize=20)
        ax.text(0.05, 0.89, '$\log{\mathcal{L}} = %.3f$' % logL, ha='left', va='top', transform=ax.transAxes, fontsize=20)
    else:
        ax.text(0.05, 0.95, '$\log{\mathcal{L}} = %.3f$' % logL, ha='left', va='top', transform=ax.transAxes, fontsize=20)

    # show size of error dispersion as Circle
    if disp is not None:

        circ1 = patches.Circle((12.5, -2.5), radius=disp, fc='b', ec='None', alpha=0.5)
        circ2 = patches.Circle((12.5, -2.5), radius=2*disp, fc='b', ec='None', alpha=0.3)
        circ3 = patches.Circle((12.5, -2.5), radius=3*disp, fc='b', ec='None', alpha=0.1)
        ax.add_artist(circ1)
        ax.add_artist(circ2)
        ax.add_artist(circ3)
        ax.text(12.5, -2.5, r'$\sigma$', color='w', fontsize=20, ha='center', va='center')

    ax.set_xlim(-5, 15)
    ax.set_ylim(-5, 15)
    ax.set_xticks([])
    ax.set_yticks([])
    fig.subplots_adjust(bottom=0.01, top=0.99, left=0.01, right=0.99)
    plt.show()

def getSelection(type="ellipse", rng=np.random):
    if type == "ellipse":
        cb = getellipse
        
    
        
if __name__ == '__main__':




    # set up test

    # sel_type = "ellipse"    # type of selection
    # disp = 0.01          # additive noise dispersion
    # default_covar = np.eye(D) * dispersion**2
    # covar_cb = partial(pygmmis.covar_callback_default, default=default)
    # bg_amp = 0.9        # fraction of background samples
    # w = 0.1             # minimum covariance regularization [data units]
    # cutoff = 5          # cutoff distance between components [sigma]
    # seed = 8366         # seed value
    # pygmmis.VERBOSITY = 1
    # pygmmis.OVERSAMPLING = 10
    # data=np.asarray([ra,dec])
    #
    # # define RNG for run
    # from numpy.random import RandomState
    # rng = RandomState(seed)
    def cb(coords):
        a=rar
        b=decr
        normra=coords[:,0]-meanra
        normdec=coords[:,1]-meandec
        return ((normra/a)**2+(normdec/b)**2 <1 )
        
    K = 2               # number of components
    D=2                 # Data dimensionality
    T=10
    data=np.asarray((ra,dec)).T
    print "data",data.shape, cb(data)
    gmm = pygmmis.GMM(K=K, D=D)
    data = pygmmis.createShared(data)
    
        
    # positional uncertainties    
    dispersion = .001
    default_covar = np.eye(D) * dispersion**2
    covar_cb = partial(pygmmis.covar_callback_default, default=default_covar)
 
        
    #background    
    footprint = data.min(axis=0), data.max(axis=0)
    bg = pygmmis.Background(footprint)
    bg.amp = 0.9
    bg.amp_max = 0.999
    bg.adjust_amp = True
    
    #run the fitter
    w = 0.01    # minimum covariance regularization, same units as data
    cutoff = 5 # segment the data set into neighborhood within 5 sigma around components
    tol = 1e-5 # tolerance on logL to terminate EM
    pygmmis.VERBOSITY = 1      # 0,1,2
    pygmmis.OVERSAMPLING = 100  # number of imputation samples per data sample

    # define RNG for deterministic behavior
    from numpy.random import RandomState
    # seed = 42
    rng = RandomState()

    # run EM
    print "before fit",gmm.mean 
    
    logL, U = pygmmis.fit(gmm, data, init_callback=pygmmis.initFromSimpleGMM,\
                          sel_callback=cb, covar_callback=covar_cb, w=w, cutoff=cutoff,\
                          background=bg, tol=tol, rng=rng)
    

    
    
    
    
    # log of p(x)
    p = gmm(data, as_log=False)
    N_s = 1000
    # draw samples from GMM
    samples = gmm.draw(N_s)

    # draw sample from the model with noise, background, and selection:
    # if you want to get the missing sample, set invert_sel=True.
    # N_orig is the estimated number of samples prior to selection
    obs_size = len(data)
    samples, covar_samples, N_orig = pygmmis.draw(gmm, obs_size, sel_callback=cb,\
                                                  invert_sel=False, orig_size=None,\
                                                  covar_callback=covar_cb,background=bg)
    # print "samples",len(np.where(U[0]==True)[0]),len(np.where(U[1]==True)[0]),len(np.where(U[2]==True)[0])
    # plt.scatter(data[U[0],0],data[U[0],1],s=0.5) 
    print "after fit", gmm.mean 
    size = 100
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)
    from matplotlib.patches import Ellipse
    
    for i in range(K):
        sigma_x = gmm.mean[0,0]
        sigma_y = gmm.mean[0,1]
        
        cov=gmm.covar
        eigval,eigvec=np.linalg.eig(cov[i])
        argsort=np.argsort(eigval)
        eigvec=eigvec[argsort]
        eigval=eigval[argsort]
        a,b=eigval[0],eigval[1]
        center=gmm.mean[i]
        amp=gmm.amp[i]
        
        alpha=np.arctan(eigvec[1,0]/eigvec[0,0])

        print "a,b,center",eigvec.shape
        x = np.linspace(center[0]-2*a, center[0]+2*a, num=size)
        y = np.linspace(center[1]-2*a, center[1]+2*a, num=size)

        x, y = np.meshgrid(x, y)
        z = (1/(2*np.pi*a*b) * np.exp(-(x**2/(2*a**2)
             + y**2/(2*b**2))))
        plt.contourf(x, y, z,alpha=0.0)

        for N in (1,2,3):
            ax.add_patch(Ellipse(center, N * a, N * b,
                                     angle=alpha * 180. / np.pi, lw=1,
                                     ec='k', fc='none'))
        plt.colorbar()
    
    plt.scatter(data[U[1],0],data[U[1],1],s=0.5)
    plt.show()

    










   