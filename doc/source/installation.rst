.. _install:

Installation
############

Requirements
============

You need the `BLAS/LAPACK <http://www.netlib.org/lapack>`_ library
and :mod:`numpy` module.

More functionalities are available if you have the `UV-CDAT <http://uv-cdat.llnl.gov>`_ package.


Download
********

You can download the sources (tarball) of the package from the
`Sourceforge repositories <http://sourceforge.net/project/showfiles.php?group_id=168272>`_ of the project.
Then, just type for example::

    $ tar xzf spanlib-X.X.X.tgz
    $ cd spanlib-X.X.X

You can also get the latest version with::

    $ git clone https://github.com/stefraynaud/spanlib.git


Compilation and installation
============================

The :envvar:`BLAS`, :envvar:`BLAS_LIB`, :envvar:`LAPACK` and :envvar:`LAPACK_LIB`
environment variables can be used to setup the installation.

You can also use the template file :file:`setup.cfg.example` to create your own
:file:`setup.cfg` file.

Then, the compilation and installation are performed with::

    $ python setup.py install


