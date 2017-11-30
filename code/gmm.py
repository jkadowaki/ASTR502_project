from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astroML.decorators import pickle_results
from astroML.plotting.tools import draw_ellipse

from astroML.density_estimation import XDGMM
hdulist = fits.open('../data/sdss_fields/leo2.fits')
tbdata = hdulist[1].data
cols = hdulist[1].columns
ra=tbdata['ra']
dec=tbdata['dec']
print "column available", cols
plt.scatter(tbdata['ra'],tbdata['dec'],s=0.01)

#plate scale =0.4 arcsec /pixel
#decerr=5e-5
#raerr=
#stack RA,DEC to an array like mag
# X = np.vstack([data_noisy[f + 'RawPSF'] for f in 'ugriz']).T
X=np.vstack((ra,dec)).T


# compute error covariance from mixing matrix
# Xcov[:, range(Xerr.shape[1]), range(Xerr.shape[1])] = Xerr ** 2
npt=len(ra)
Xcov = np.zeros((npt,2,2))
Xcov[:,range(2),range(2)]=5e-5

# @pickle_results("GMM.pkl")
def compute_XD(n_clusters=2, rseed=0, n_iter=10, verbose=True):
    np.random.seed(rseed)
    clf = XDGMM(n_clusters, n_iter=n_iter, tol=1E-5, verbose=verbose)
    clf.fit(X, Xcov)
    return clf

clf = compute_XD(2)
print "Xshape",X.shape
X_sample = clf.sample(X.shape[0])
rasamp=X_sample[:,0]
decsamp=X_sample[:,1]

fig = plt.figure()
xmin=min(ra)
xmax=max(ra)
ymin=min(dec)
ymax=max(dec)
ax1 = fig.add_subplot(131)
ax1.scatter(ra, dec , s=0.5, lw=0, c='k')
ax2 = fig.add_subplot(132)
ax2.scatter(rasamp, decsamp , s=0.5, lw=0, c='k')

ax3 = fig.add_subplot(133)
print "clf.n_components",clf.n_components

for i in range(clf.n_components):
    print clf.mu[i, 0:2], clf.V[i, 0:2, 0:2]
    draw_ellipse(clf.mu[i, 0:2], clf.V[i, 0:2, 0:2], scales=[2],
                 ec='k', fc='gray', alpha=0.2, ax=ax3)
ax1.set_xlim([xmin,xmax])
ax2.set_xlim([xmin,xmax])
ax3.set_xlim([xmin,xmax])
ax1.set_ylim([ymin,ymax])
ax2.set_ylim([ymin,ymax])
ax3.set_ylim([ymin,ymax])
ax1.set_title("Original Data")
ax2.set_title("Extreme Deconvolution Resampling ")
ax3.set_title("Extreme Deconvolution Clusters location")
plt.show()



