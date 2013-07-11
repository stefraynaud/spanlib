.. _install:

Installation
############

Download
********

You can download the sources (tarball) of the package from the 
`Sourceforge repositories <http://sourceforge.net/project/showfiles.php?group_id=168272>`_ of the project.
Then, just type for example::

    $ tar xzf spanlib-2.0.tgz
    $ cd spanlib-2.0

You can also get the latest version with::
    
    $ svn export https://svn.code.sf.net/p/spanlib/svn/ spanlib-svn


Fortran library
***************

Requirements
============

You need a F90 compiler to compile it, and `BLAS/LAPACK <http://www.netlib.org/lapack>`_
libraries compilated using this F90 compiler to be able to link with library.
If you are using the Intel fortran compiler, you may be able to link the library against
the `Math Kernel Library <http://www.intel.com/cd/software/products/asmo-na/eng/perflib/307757.htm>`_.

To run the F90 examples, you also need the F90 `netcdf library <http://www.unidata.ucar.edu/software/netcdf>`_.


Compilation and installation
============================

First, compile it with::

    $ ./configure
    $ make

Second, install it as root::

    $ sudo make install

Here is an example of :command:`./configure`::

    $ ./configure --with-blas-lib=/usr/local/install/lapack-3.0/lib \
    --with-netcdf-lib=/usr/local/install/netcdf-3.6.1/lib \
    --with-netcdf-inc=/usr/local/install/netcdf-3.6.1/include \
    --prefix=$HOME --with-pythondir=$HOME/python FC=ifort

You can list all :command:`configure` options using::

    $ ./configure --help

You can also list usefull :command:`make` targets with::

    $ make  help

The following variables can be used for the dependencies:
    
.. envvar:: BLAS

    The name of the BLAS library provided to the compiler
    using the ``-l${BLAS}``.
    
.. envvar:: BLAS_LIB

    The name of the BLAS library directory provided to the compiler
    using the ``-L${BLAS_LIB}``.


.. envvar:: LAPACK

    The name of the LAPACK library provided to the compiler
    using the ``-l${LAPACK}``.
    
.. envvar:: LAPACK_LIB

    The name of the LAPACK library directory provided to the compiler
    using the ``-L${LAPACK_LIB}``.

.. envvar:: NETCDF_LIB

    The name of the netcdf library directory provided to the compiler
    using the ``-L${NETCDF_LIB}``.
    
.. envvar:: NETCDF_INC

    The name of the netcdf include directory provided to the compiler
    using the ``-L${NETCDF_INC}``.
    
Python package
**************


Requirements
============

You need the `BLAS/LAPACK <http://www.netlib.org/lapack>`_ library
and :mod:`numpy` module.

More functionalities are available if you have the `UV-CDAT <http://uv-cdat.llnl.gov>`_ package.


Compilation and installation
============================

The :envvar:`BLAS`, :envvar:`BLAS_LIB`, :envvar:`LAPACK` and :envvar:`LAPACK_LIB`
environment variables can be used to setup the installation.

You can also use the template file :file:`setup.cfg.example to create your own
:file:`setup.cfg` file.

Then, the compilation and installation are performed with:
    
    >>> python setup.py install


