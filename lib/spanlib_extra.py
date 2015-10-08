"""Extra stuff for instance for running test scripts"""

import sys,  os
from ConfigParser import SafeConfigParser
import numpy as npy

def setup_data1(nx=100, nt=120, xfreq=2, tfreq=3, xyfact=3, masked=True):
    """Generate space-time data sample with a single spatial dimension

    Returned array has masked values if ``masked`` is set to ``True``
    """
    var = npy.zeros((nt,nx))

    xx = npy.arange(nx,dtype='d')
    tt = npy.arange(nt,dtype='d')

    for ifreq in [.5, 3]:
        tvar = (npy.cos(2*npy.pi * tfreq * tt/(nt-1))*tt/(nt-1)).reshape((nt,1))
        xvar = npy.cos(2*npy.pi * xfreq*ifreq * xx/(nx-1)).reshape((1,nx))
        var += npy.multiply(tvar,xvar)
    for i in xrange(nt):
        for j in xrange(nx):
            var[i,j] += xyfact*(i+j)/(nx+nt)
    jj, ii = npy.indices((nt, nx)).astype('d')
    var += npy.exp(-((ii-nx/2.)/(nx/2.))**2 - ((jj-nt/2.)/(nt/2.))**2)
    var = npy.ma.asarray(var)
    if masked and var.shape[1]>1: var[:, 1] = npy.ma.masked

    return var

def setup_data2(nx=30, ny=20, **kwargs):
    """Same as setup_data1 but with two spatial dimensions

    Extra keywords are passed to setup_data1
    """

    data1 = setup_data1(nx=nx*ny, **kwargs)
    return data1.reshape((-1, ny, nx))

def setup_data0(nt = 300):
    """Get a 1D signal suitable for SSA tests"""
    t = npy.arange(nt*1.)
    p0 = 23.
    p1 = p0*3.6
    p2 = p0/10.
    return npy.cos(t/p0*2*npy.pi)+npy.sin(t/p1*2*npy.pi+1.3)+npy.cos(t/p2*2*npy.pi+1.4)


def pca_numpy(var, nmode, cov=None, evsum=False):
    """PCA using numpy library"""
    if cov is None:
        if not npy.ma.isMA(var):
            cov = npy.cov(var, bias=1, rowvar=0)
        else:
            nn = 1-~var.mask.astype('i')
            cov = npy.dot(var.T, var)
            nn2 = npy.dot(nn.T, nn)
            cov /= nn2
            del nn, nn2

    _ev, eof = npy.linalg.eigh(cov)
    isort = npy.argsort(_ev)[::-1]
    ev = _ev[isort][:nmode]
    eof = eof[:, isort][:, :nmode]
    pc = npy.ma.dot(var, eof)
    if evsum: return eof, pc, ev, _ev.sum()
    return eof, pc, ev

def svd_numpy(varl, varr, nmode):
    """SVD of the covariance matrix of two variables"""
    cov = npy.dot(varl, varr.T)/varl.shape[1]
    eofl, ev, eofrt = npy.linalg.svd(cov, full_matrices=False)
    eofr = eofrt.T
    eofl = eofl[::nmode]
    eofr = eofr[::nmode]
    pcl = npy.dot(varl.T, eofl)
    pcr = npy.dot(varr.T, eofr)
    return eofl, eofr, pcl, pcr, ev

def gensin1d(per=50, nper=10, p=0.2):
    """Generate an 1D sinusoid as an numpy array of length ``nt``, period ``T`` and phase ``p``

    :Params:

        - **per**, optional: Period in time steps.
        - **nper**, optional: Sample size in number of periods.
        - **p**, optional: Phase relative to the period (within [-1,1]).
    """
    p *= 2*npy.pi
    return npy.ma.sin(npy.arange(int(round(nper*per)))*2*npy.pi/per+p)

def gensin2d(xper=50, xnper=10, xp=0.2, yper=15, ynper=10, yp=0.6):
    """Generate a 2D sinusoid like :func:`gensin`"""
    x = gensin1d(per=xper, nper=xnper, p=xp)
    y = gensin1d(per=yper, nper=ynper, p=yp)
    xx, yy = npy.meshgrid(x, y)
    return npy.ma.asarray(xx*yy)
