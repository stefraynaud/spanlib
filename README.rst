SpanLib - Spectral Analysis Library
===================================

Purpose
-------
SpanLib is a python library that wraps fortran routines
to perform analyses such as Principal Component
Analysis (PCA) and Multichannel Singular
Spectrum Analysis (MSSA).
This packages is suitable for analysis studies of
climate or financial variability.

Requirements
------------
- Python
- Numpy
- Blas/Lapack (at compilation time)

Web site
--------
http://spanlib.sourceforge.net

Authors
-------
Stephane Raynaud and Charles Doutriaux.

Licence
-------
Lesser GNU Public Licence.
See LICENSE file.

Documentation
-------------
- In the package tree: doc/index.html or doc/spanlib.html
- Once the packaged is installed, typically in:
  file://<prefix_dir>/share/spanlib/html
- On the web: http://sourceforge.net/docman/index.php?group_id=168272
  (no html formatting)

Quick installation
------------------

    $ python setup.py install # basic install
    $ python setup.py install --user # install in user home
    $ python setup.py install build_ext --inplace # install inplace for dev


See INSTALL file.


