######################################################################
## SpanLib, Raynaud 2006
######################################################################

lib_LIBRARIES = libspanlib.a
include_HEADERS=spanlib.mod
libspanlib_a_SOURCES=spanlib.f90
SUFFIXES=.f90 .a .o .mod
AM_FCFLAGS=$(LAPACK_LIBDIR) $(LAPACK) $(BLAS_LIBDIR) $(BLAS)
CLEANFILES=*.o *.a *.mod *.pyf *.pfp

GENPYF=$(top_srcdir)/scripts/genpyf.pl
LAPACKDIR=/usr/local/install/lapack-3.0/lib/
MOD=$(PACKAGE).mod
WRAPPER=$(PACKAGE)_pywrap.f90
PYF=$(PACKAGE).pyf
PFP=$(PACKAGE).pfp
if WITH_PYTHON
PYTHON_MODULE=python_module
endif

######################################################################
######################################################################

all-local: $(LIBRARIES) $(HEADERS) $(PYTHON_MODULE)

if WITH_PYTHON
install-exec-hook:
	@$(NORMAL_INSTALL)
	test -z "$(PYPACKDIR)" || $(mkdir_p) $(PYPACKDIR)/spanlib_dir
	@echo "Installing $(PACKAGE) python module in $(PYPACKDIR)"
	cp `ls -td build/lib.*/spanlib|head -1`/* $(PYPACKDIR)/spanlib_dir
	echo spanlib_dir > $(PYPACKDIR)/spanlib_dir.pth
uninstall-hook:
	@echo "Uninstalling $(PACKAGE) python module"
	rm -rf $(PYPACKDIR)/spanlib_dir $(PYPACKDIR)/spanlib_dir.pth
clean-local:
	-rm -rf build
endif

fortran_library: $(LIBRARIES) $(HEADERS)

python_module: pfp $(WRAPPER) $(PYF) $(LIBRARIES) $(HEADERS)
	$(PYFORT) $(PYFORT_BUILD) $(PFP)

$(PYF): $(WRAPPER)
	@echo Creating $(PYF)
	$(top_srcdir)/scripts/genpyf.pl $(WRAPPER) $(PYF)

pfp:
	@echo Creating $(PFP)
	@echo "pyf('$(PYF)',\
	sources=['$(WRAPPER)'],\
	libraries='$(PACKAGE) $(PYFORT_LIBS)',\
	library_directories='. $(PYFORT_DIRS)',\
	compiler_options='-I.',\
	python_directory='$(top_srcdir)/lib',\
	package_name='$(PACKAGE)',\
	generate_as='$(PACKAGE)_fort',\
	freeform=1)" > $(PFP)

