######################################################################
## SpanLib, Raynaud 2006-2013
######################################################################

import os
import sys
import re
from ConfigParser import SafeConfigParser
from numpy.distutils.core import setup, Extension
from numpy.f2py import crackfortran

# Revision number
rootdir = os.path.dirname(__file__)
# - from __init__.py
f = open(os.path.join(rootdir, 'lib/spanlib/__init__.py'))
for line in f:
    line = line[:-1].strip()
    if line.startswith('__version__'):
        exec(line)
        release = __version__
        break
f.close()
version_sphinx = release
# - from svn
if os.path.exists(os.path.join(rootdir,'.svn/entries')):
    f = open('.svn/entries')
    line = f.readlines()[3][:-1]
    f.close()
    release += '-svn%i'%eval(line)
release_sphinx = release

# Some info
version = release
description='Python extension to the spanlib fortran library'
author = 'Stephane Raynaud and Charles Doutriaux'
author_email = 'stephane.raynaud@gmail.com'
url="http://spanlib.sf.net"

# Gather up all the files we need
spanlib_files = ['src/spanlib.pyf', 'src/spanlib_pywrap.f90', 'src/spanlib.f90',  'src/anaxv.f90']
anaxv_files = ['src/anaxv.pyf', 'src/anaxv.f90']


# Paths for libs
libs = ['lapack', 'blas']
liddirs = ['/usr/local/lib']
incdirs =  ['/usr/local/include']
cfg = SafeConfigParser()
cfg.read('setup.cfg')
site_libs = []
site_libdirs = []
site_incdirs = []
for libname in libs:
    lib = libdir =  ''
    # - Libraries
    if os.getenv(libname.upper()) is not None:
        lib = os.getenv(libname.upper())
    elif cfg.has_option('blaslapack', libname):
        lib = cfg.get('blaslapack', libname)
    if lib.startswith('-l'): lib = lib[2:]
    if lib: site_libs.append(lib) 
    # - Library dirs
    if os.getenv(libname.upper()+'_LIB') is not None:
        libdir = os.getenv(libname.upper()+'_LIB')
    elif cfg.has_option('blaslapack', libname+'_lib'):
        libdir = cfg.get('blaslapack', libname+'_lib')
    if libdir.startswith('-L'): libdir = libdir[2:]
    if libdir: site_libdirs.append(libdir)      
# - final setup
if len(site_libs): libs = site_libs
if len(site_libdirs): libdirs = site_libdirs
if len(site_incdirs): incdirs = site_incdirs
print 'incdirs',incdirs
print 'libdirs',libdirs    



# Some useful directories.  
## from distutils.sysconfig import get_python_inc, get_python_lib
## python_incdir = os.path.join( get_python_inc(plat_specific=1) )
## python_libdir = os.path.join( get_python_lib(plat_specific=1) )

# Special setups
extra_link_args=[]
if sys.platform=='darwin':
    extra_link_args += ['-bundle','-bundle_loader '+sys.prefix+'/bin/python']
kwext = dict(libraries=libs, library_dirs=libdirs, include_dirs=incdirs, extra_link_args=extra_link_args)

if __name__=='__main__':
    
    # Generate pyf files
    crackfortran.f77modulename = '_fortran'
    pyfcode = crackfortran.crack2fortran(crackfortran.crackfortran(['src/spanlib_pywrap.f90']))
    f = open('src/spanlib.pyf', 'w')
    f.write(pyfcode)
    f.close()
    crackfortran.f77modulename = 'anaxv'
    pyfcode = crackfortran.crack2fortran(crackfortran.crackfortran(['src/anaxv.f90']))
    f = open('src/anaxv.pyf', 'w')
    f.write(pyfcode)
    f.close()

    # Setup the python module
    s = setup(name="spanlib",
        version=version,
        description=description,
        author=author,
        author_email=author_email,
        maintainer=author,
        maintainer_email=author_email,
        license="GNU LGPL",
    
        # Fortran wrapper
        ext_modules = [
            Extension('spanlib._fortran', spanlib_files, **kwext), 
            Extension('spanlib.anaxv', anaxv_files, **kwext)
        ],
          
        # Install these to their own directory
        package_dir={'spanlib':'lib/spanlib'},
        packages = ["spanlib"],
          
    )
    
    # Save info
    if 'build' in s.command_obj:
        cfgfile = 'config.cfg'
        cfg = SafeConfigParser()
        cfg.read(cfgfile)
        if not cfg.has_section('paths'):
            cfg.add_section('paths')
        cfg.set('paths', 'build_lib', os.path.abspath(s.command_obj['build'].build_lib)) 
        if os.path.exists(cfgfile):
            os.remove(cfgfile)
        f = open(cfgfile, 'w')
        cfg.write(f)
        f.close()




