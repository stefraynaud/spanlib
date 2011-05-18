from ConfigParser import SafeConfigParser
import os, sys
import numpy as npy

def insert_local_path():
    if os.path.exists('../config.cfg'):
        cfg = SafeConfigParser()
        cfg.read('../config.cfg')
        if cfg.has_option('paths', 'build_lib'):
            sys.path.insert(0, cfg.get('paths', 'build_lib'))

def setup_data1(nx=100, nt=120, xfreq=2, tfreq=3, xyfact=3):
    """Generate space-time data sample with a single spatial dimension
    
    Returned array has masked values
    """
    var = npy.zeros((nt,nx))
    
    xx = npy.arange(nx,dtype='d')
    tt = npy.arange(nt,dtype='d')
    
    for ifreq in .5,3:
        tvar = (npy.cos(2*npy.pi * tfreq * tt/(nt-1))*tt/(nt-1)).reshape((nt,1))
        xvar = npy.cos(2*npy.pi * xfreq*ifreq * xx/(nx-1)).reshape((1,nx))
        var += npy.multiply(tvar,xvar)
    for i in xrange(nt):
        for j in xrange(nx):
            var[i,j] += xyfact*(i+j)/(nx+nt)
    jj, ii = npy.indices((nt, nx)).astype('d')
    var += npy.exp(-((ii-nx/2.)/(nx/2.))**2 - ((jj-nt/2.)/(nt/2.))**2)
    var = npy.ma.asarray(var)
    var[:, 1] = npy.ma.masked

    return var

def setup_data2(nx=30, ny=20, **kwargs):
    """Same as setup_data1 but with two spatial dimensions
    
    Extra keywords are passed to setup_data1
    """
    
    data1 = setup_data1(nx=nx*ny, **kwargs)
    return data1.reshape((-1, ny, nx))

def pca_numpy(var, nmode):
    """PCA using numpy library"""
    cov = npy.cov(var.T, bias=0)
    ev, eof = npy.linalg.eigh(cov)
    isort = npy.argsort(ev)[::-1]
    ev = ev[isort][:nmode]
    eof = eof[:, isort][:, :nmode]
    pc = npy.dot(var, eof)
    return eof, pc, ev


