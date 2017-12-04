# Written by Ben Lew on Dec 2nd, 2017
# wplew@lpl.arizona.edu
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
from matplotlib.patches import Ellipse
from numpy.random import RandomState
from astroML.density_estimation import XDGMM

def in_ellipse(x,y,a,b,alpha,datax,datay):
#https://stackoverflow.com/questions/7946187/point-and-ellipse-rotated-position-test-algorithm
    a=a*3
    b=b*3
    first=  (np.cos(alpha)*(datax-x)+np.sin(alpha)*(datay-y))**2/a**2 
    second= (np.sin(alpha)*(datax-x)+np.cos(alpha)*(datay-y))**2/b**2 
    print  alpha,first +second
    print first+second < 1
    return (first+second) < 1  
hdulist = fits.open('../data/sdss_fields/draco.fits')
# hdulist = fits.open('../data/mock_fields/background/background_0.5_1000_0.15_10000.fits')

tbdata = hdulist[1].data
cols = hdulist[1].columns
ra=tbdata['ra']
dec=tbdata['dec']

meanra=np.mean(ra)
meandec=np.mean(dec)
rar=max(ra)-meanra
decr=max(dec-meandec)
def cb(coords):
    a=rar
    b=decr
    normra=coords[:,0]-meanra
    normdec=coords[:,1]-meandec
    return ((normra/a)**2+(normdec/b)**2 <1 )
    
K = 2               # number of components
D=2                 # Data dimensionality
data=np.asarray((ra,dec)).T
print "data",data.shape, cb(data)
# seed = 42
niter=50
centerlist=[]
alphalist=[]
alist=[]
blist=[]
randcenterlist=[]
realcenter=[meanra,meandec]
maxL=-100

origra=ra
origdec=dec
simamp=1
gaussnum=1000

for j in range(niter):
    # simulate fake gaussian
    rand=int(np.random.rand(1)*len(ra))
    randcenter=[ra[rand],dec[rand]]
    simcov=np.diagflat([0.1,0.1])/simamp
    simgauss= np.random.multivariate_normal(randcenter,simcov,gaussnum)
    randcenterlist.append(randcenter)
    # print "simgauss",simgauss.shape
    # plt.scatter(simgauss[:,0],simgauss[:,1])
    # plt.show()

    ra= np.concatenate((origra,simgauss[:,0]))
    dec=np.concatenate((origdec,simgauss[:,1]))

    # ----------------gmm initialization----------------
    gmm = pygmmis.GMM(K=K, D=D)
    data = pygmmis.createShared(data)

    # positional uncertainties    
    dispersion = .01
    default_covar = np.eye(D) * dispersion**2
    covar_cb = partial(pygmmis.covar_callback_default, default=default_covar)

    
    #background    
    footprint = data.min(axis=0), data.max(axis=0)
    bg = pygmmis.Background(footprint)
    bg.amp = 0.95
    bg.amp_max = 0.999
    bg.adjust_amp = True

    #run the fitter
    w = 0.1    # minimum covariance regularization, same units as data
    cutoff = 1.5 # segment the data set into neighborhood within cutoff* sigma around components
    tol = 1e-9 # tolerance on logL to terminate EM
    pygmmis.VERBOSITY = 2      # 0,1,2
    pygmmis.OVERSAMPLING = 100  # number of imputation samples per data sample

    # define RNG for deterministic behavior
    rng = RandomState()

    # run EM
    # print "before fit",gmm.mean
    logL, U = pygmmis.fit(gmm, data, init_callback=pygmmis.initFromSimpleGMM,\
                          sel_callback=cb, covar_callback=covar_cb, w=w, cutoff=cutoff,\
                          background=bg, tol=tol, rng=rng)
    falselist=np.where(U[0]==False)[0]
    truelist=np.where(U[0]==True)[0]

    # print "num of U==true", len(truelist)*1.0/len(U[0])
    print logL
    if logL > maxL:
        bestgmm=gmm
        maxL=logL
    size = 100

    color=["red","green","magenta","black","blue","orange"]

    for i in range(K):
        centerlist.append(gmm.mean[i])
        cov=gmm.covar
        eigval,eigvec=np.linalg.eig(cov[i])
        argsort=np.argsort(eigval)
        eigvec=eigvec[argsort]
        eigval=eigval[argsort]
        a,b=eigval[0],eigval[1]
        alpha=np.arctan(eigvec[1,0]/eigvec[0,0])
        
        alphalist.append(alpha)
        alist.append(a)
        blist.append(b)
        

# plotgraph=111
fig = plt.figure(figsize=(5, 5))

ax = fig.add_subplot(111)
# ax.patch.set_facecolor('purple')

for i in range(K):
    sigma_x = bestgmm.mean[0,0]
    sigma_y = bestgmm.mean[0,1]

    cov=bestgmm.covar
    eigval,eigvec=np.linalg.eig(cov[i])
    argsort=np.argsort(eigval)
    eigvec=eigvec[argsort]
    eigval=eigval[argsort]
    a,b=eigval[0],eigval[1]
    center=bestgmm.mean[i]
    amp=bestgmm.amp[i]
    # centerlist.append(center)


    alpha=np.arctan(eigvec[1,0]/eigvec[0,0])
    print "amp 1 and 2",color[i],amp
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
                                 ec='red', fc='none'))
# plt.colorbar()
# count how many detect injected gaussian within 1-sigmna#

count=0
centerlist=np.asarray(centerlist)
print "real center",realcenter
for p in range(niter):
    print p, "iter:"
    randcenter=randcenterlist[p]
    if (in_ellipse(randcenter[0],randcenter[1],alist[p*K],blist[p*K],alphalist[p*K],centerlist[p*K,0],centerlist[p*K,1])) or \
        (in_ellipse(randcenter[0],randcenter[1],alist[p*K+1],blist[p*K+1],alphalist[p*K+1],centerlist[p*K+1,0],centerlist[p*K+1,1])):
        count+=1
        
   
plt.title("GMM for Draco +injected gauss pts="+str(gaussnum)+", amp="+str(simamp)+"\n True Positive:"+str(1.0*count/niter))
plt.scatter(data[:,0],data[:,1],s=0.4,label="data")
# plt.scatter(randcenter[0],randcenter[1],marker="x",s=100,color='k',label="injected Gaussian")
# plt.scatter(data[U[1],0],data[U[1],1],s=0.5)
# plt.scatter(data[falselist,0],data[falselist,1],s=0.5,alpha=0.5)
# plt.scatter(data[U[1],0],data[U[1],1],s=0.5,alpha=0.5,color="red")
# fig2=plt.figure()
# ax2 = fig2.add_subplot(111)
#
plt.scatter(centerlist[:,0],centerlist[:,1],marker="o",color="green",alpha=0.4,label="best-fit center")
print "center shape", centerlist.shape
print "centerlist ra",centerlist[:,0]
# H, xedges, yedges= np.histogram2d(centerlist[:,0], centerlist[:,1], bins=6, range=None, normed=False, weights=None)
# X, Y = np.meshgrid(xedges, yedges)
# cap=ax2.pcolormesh(X, Y, H)
# plt.colorbar(cap)
plt.legend()

plt.show()

    
            









   