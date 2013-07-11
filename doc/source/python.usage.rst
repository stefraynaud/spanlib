Schematic usages
================

Notations:
    
    - ``nt``: number of observations.
    - ``ns``: number of spatial channels.
    - ``nx``, ``ny``: number of spatial channels along X and Y dimensions.
    - ``npca``: number PCA modes.
    - ``nmssa``: number MSSA modes.
    - ``nsvd``: number SVD modes.
    - ``nw``: size of the window in MSSA analysis.


PCA
---

Basic analysis
~~~~~~~~~~~~~~

>>> span = Analyzer(data, npca=npca) # data of shape (nt,ns) or (nt,ny,nx) for example
>>> span.pca() # optional call
>>> pc = span.pca_pc() # (npca,nt)
>>> eof = span.pca_eof() # (npca,ns) or (npca,ny,nx)
>>> ev = span.pca_ev() # (npca)

Expansion coefficients
----------------------

>>> ec = span.pca_ec(data2) # projection of data2 onto EOFs of data
>>> ec = span.pca_ec(xdata=data2) # equivalent
>>> ec = span.pca_ec(xdata=data, xeof=span.pca_eof()) # = pc!

Reconstructions
---------------

>>> rec = span.pca_rec() #  (nt,ns) or (nt,ny,nx)
>>> rec = span.pca_rec(xeof=eof2, xpc=pc2) # with alternate EOFs and PCs
>>> rec = span.pca_rec(xpc= span.pca_ec(xdata=data2)) # filter data2 with data1
>>> rec = span.pca_rec(modes=[-4,7]) # modes 0,1,2,3,4,7 # not 5, 6 and other


SSA or MSSA
-----------

Basic usages
~~~~~~~~~~~~

Usage are the same as for PCA, but shapes are differents.

>>> span = Analyzer(data, nmssa=nmssa)
>>> pc = span.mssa_pc() # (nmssa,nt-nw+1)
>>> eof = span.mssa_eof() # (nmssa,nw,ns) or (nmssa,nw,ny,nx)
>>> ev = span.mssa_ev() # (nmssa)
>>> ec = span.mssa_ec(xdata=xdata, xeof=xeof) # (nmssa,nt-nw+1)
>>> rec = span.mssa_rec(xdata=xdata, xpc=xpc) #  (nt,ns) or (nt,ny,nx)

Monte-Carlo significance test
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

>>> ev,evmin,evmax = span.mssa_ev(mctest=True)
>>> print 'Bad: ', (ev>evmin) & (ev<evmax)

SVD
---

Basic usages
~~~~~~~~~~~~

Usage are the same as for PCA, except that you work two groups of variables (see :ref:`user.schem.py.sev`).

>>> span = DualAnalyzer(left_data, right_data, sequential=True, nsvd=nsvd)
>>> left_pc,right_pc = span.svd_pc() # (nsvd,nt)
>>> left_eof, right_eof = span.svd_eof() # (nsvd,ns) or (nmssa,nw,ny,nx)
>>> ev = span.svd_ev() # (nsvd)


Other cases
-----------

Analyzing several variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~

>>> span = Analyzer((data1,data2), norms=[1,1]) # joint analysis
>>> eof1, eof2 = span.pca_eof()
>>> pc = span.pca_pc()
>>> ev = span.pca_ev()


