from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astroML.decorators import pickle_results
from astroML.plotting.tools import draw_ellipse

from astroML.density_estimation import XDGMM
hdulist = fits.open('../data/sdss_fields/willman1.fits')
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
Z=np.vstack((ra,dec)).T
binnum=100
binz,xedge,yedge=np.histogram2d(ra,dec,bins=binnum)
realfilter=np.where(binz>0)
binz=binz-np.mean(binz[realfilter])
fig=plt.figure()
ax = fig.add_subplot(111)
X, Y = np.meshgrid(xedge, yedge)
cax=ax.pcolormesh(X, Y, binz)
plt.colorbar(cax)
plt.show()

# ragrid=np.linspace(min(ra),max(ra),num=binnum)
# decgrid=np.linspace(min(ra),max(dec),num=binnum)
newZ=[]
print "lenbiz",len(binz[:,0]),len(binz[0,:])
for i in range(len(binz[:,0])):
    for j in range(len(binz[0,:])):
        if int(binz[i,j])>0:
            
            temp=[[xedge[i],yedge[j]] for k in range(int(binz[i,j]))]
            newZ+=temp
# # calculate area of the points

newZ=np.asarray(newZ)

 
# from astropy import units as u
# from astropy.coordinates import SkyCoord
# minpt=SkyCoord(ra=min(ra)*u.degree,dec=min(dec)*u.degree)
# maxpt=SkyCoord(ra=max(ra)*u.degree,dec=max(dec)*u.degree)
# ang=minpt.separation(maxpt)/2. # angular radius

# for i in range(binnum):
    
# compute error covariance from mixing matrix
# Xcov[:, range(Xerr.shape[1]), range(Xerr.shape[1])] = Xerr ** 2
npt=len(ra)
Zcov = np.zeros((len(newZ[:,0]),2,2))
Zcov[:,range(2),range(2)]=5e-5
print "z shape",newZ.shape, binz.shape, Zcov.shape

# @pickle_results("willman1.pkl")
def compute_XD(n_clusters=2, rseed=0, n_iter=30, verbose=True):
    np.random.seed(rseed)
    clf = XDGMM(n_clusters, n_iter=n_iter, tol=1E-5, verbose=verbose)
    clf.fit(newZ, Zcov)
    return clf

clf = compute_XD(7)
print "Zshape",Z.shape
Z_sample = clf.sample(Z.shape[0])
rasamp=Z_sample[:,0]
decsamp=Z_sample[:,1]

fig = plt.figure()
xmin=min(ra)
xmax=max(ra)
ymin=min(dec)
ymax=max(dec)
ax1 = fig.add_subplot(131)
# for i in range(newZ.shape[0]):
ax1.scatter(newZ[:,0], newZ[:,1] , s=0.5, lw=0, c='k')
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
ax1.set_title("Binnws Data willman1")
ax2.set_title("Extreme Deconvolution Resampling ")
ax3.set_title("Extreme Deconvolution Clusters location")
plt.show()



