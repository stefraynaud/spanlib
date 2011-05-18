######################################################################
## SpanLib, Raynaud 2006-2011
######################################################################

import os
import sys
import re
from ConfigParser import SafeConfigParser

# Some info
version='0.9' # Like 0.1
description='Python extension to spanlib fortran library'
author = 'Stephane Raynaud and Charles Doutriaux'
author_email = 'stephane.raynaud@gmail.com'
url="http://spanlib.sf.net"
from numpy.distutils.core import setup, Extension


# Gather up all the files we need.
files = ['src/spanlib.pyf', 'src/spanlib_pywrap.f90', 'src/spanlib.f90']

# Generate pyf file
fw = open('src/spanlib_pywrap.f90')
fp = open('src/spanlib.pyf','w')
def add_nextline(ntab=0,nnl=1):
    global i
    fline = ntab*'\t'+lines[i]
    while lines[i].endswith('&'):
        i+=1
        fline += lines[i]
    fline = fline.replace('&','')+'\n'*nnl
    fp.write(fline)
    i+=1
redec = re.compile('(real|logical|integer)')
fp.write("! -*- Mode: f90 -*-\n\npython module spanlib_fort\ninterface\n")
lines = [l[:-1].strip() for l in fw.xreadlines()]
i = 0
while i < len(lines):
    # Find subroutine blocks
    if "subroutine" in lines[i].lower():
        
        # Subroutine declaration
        add_nextline()
        
        # Variales
        while redec.search(lines[i]) is None: i+=1
        while len(lines[i]): add_nextline(ntab=1)
        
        # End subroutine
        while 'end subroutine' not in lines[i].lower(): i+=1
        add_nextline(nnl=2)
        
    i+=1
fw.close()
fp.write("end interface\nend python module\n\n")
fp.close()


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

## Lapack95 modules and functions
#l95mod = {'mod':'f95_lapack', 'pre':'la_'}
#for tt in 'mod', 'pre':
    #if os.getenv('LAPACK95_'+tt.upper()) is not None:
        #l95mod[tt] = os.getenv('LAPACK95_'+tt.upper())
    #elif cfg.has_option('blaslapack', 'lapack95_'+tt):
        #l95mod[tt] = cfg.get('blaslapack', 'lapack95_'+tt)
#f = open('src/template.spanlib_lapack95.f90')
#txt = f.read()
#f.close()
#txt = re.sub(r'\b(?i)f95_lapack\b', l95mod['mod'], txt)
#txt = re.sub(r'\b(?i)la_(\S+)', r'la_\1 => %s\1'%l95mod['pre'], txt)
#f = open('src/spanlib_lapack95.f90', 'w')
#f.write(txt)
#f.close()


# Some useful directories.  
## from distutils.sysconfig import get_python_inc, get_python_lib
## python_incdir = os.path.join( get_python_inc(plat_specific=1) )
## python_libdir = os.path.join( get_python_lib(plat_specific=1) )

# Special setups
extra_link_args=[]
if sys.platform=='darwin':
    extra_link_args += ['-bundle','-bundle_loader '+sys.prefix+'/bin/python']
    
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
        Extension('spanlib.spanlib_fort',
            files,
            libraries=libs,
            library_dirs=libdirs,
            include_dirs=incdirs,
            extra_link_args=extra_link_args,
        ),
    ],
      
    # Install these to their own directory
    package_dir={'spanlib':'lib'},
    packages = ["spanlib"],
      
)

# Save info
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




