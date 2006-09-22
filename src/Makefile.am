######################################################################
## SpanLib, Raynaud 2006
######################################################################

lib_LIBRARIES = libspanlib.a
nodist_include_HEADERS=spanlib_precision.mod spanlib.mod
nodist_libspanlib_a_SOURCES=spanlib_precision.f90
dist_libspanlib_a_SOURCES=spanlib.f90
noinst_PROGRAMS = getf90precision
getf90precision_SOURCES = getf90precision.f90
CLEANFILES=*.pyf *.pfp *.mod
EXTRA_DIST=spanlib_pywrap.f90
AM_FCFLAGS = $(LAPACK_INC)

WRAPPER = $(PACKAGE)_pywrap.f90
PYTHON_PACKAGE = $(PACKAGE)_fort
if WITH_PYTHON
PYTHON_LIBRARY_DIR = $(srcdir)/build/$(PYFORT_BUILD_DIR)/$(PACKAGE)
PYTHON_LIBRARY = $(PYTHON_LIBRARY_DIR)/$(PACKAGE)_fort.so
endif

######################################################################
######################################################################

all-local: $(LIBRARIES) $(HEADERS) $(PYTHON_LIBRARY)

spanlib.o: spanlib_precision.mod

spanlib_precision.mod: spanlib_precision.o

spanlib_precision.f90: spanlib_precision.F90
	test -f $@ || cp $^ $@
	@precision=`grep "wp *= *dp" $@ | cut -d\  -f6` ; \
	test "$$precision" = "$(WP)" || \
	( rm -f $@ && sed "s/rp *= *rdouble/wp = $(WP)/" $^ > $@ && \
	echo "Switched to $(PRECISION) precision" )

if WITH_PYTHON
install-exec-hook:
	test -z "$(PYTHONDIR)" || $(mkdir_p) $(PYTHONDIR)
	@echo "Installing $(PACKAGE) python module in $(PYTHONDIR)"
	$(INSTALL) -d -m 755 $(PYTHONDIR)/spanlib
	$(INSTALL_DATA) $(PYTHON_LIBRARY_DIR)/* $(PYTHONDIR)/spanlib

uninstall-hook:
	@echo "Uninstalling $(PACKAGE) python module"
	-rm -rf $(PYTHONDIR)/spanlib

clean-local:
	-rm -rf build $(PACKAGE).pyf $(PACKAGE).pfp

$(PYTHON_LIBRARY): $(HEADERS) $(PACKAGE).pfp $(WRAPPER) $(PACKAGE).pyf $(LIBRARIES) $(top_srcdir)/lib/__init__.py $(top_srcdir)/lib/spanlib_python.py
	$(PYFORT) -c $(FC) $(PYFORT_BUILD_OPT) $(PACKAGE).pfp

$(PACKAGE).pyf: $(top_srcdir)/scripts/genpyf.pl $(WRAPPER) getf90precision
	@echo Creating $@
	$(top_srcdir)/scripts/genpyf.pl $(WRAPPER) $@ $(PRECISION) `getf90precision`

$(PACKAGE).pfp:
	@echo Creating $@
	@echo "pyf('$(PACKAGE).pyf',\
	sources=['$(WRAPPER)'],\
	libraries='$(PACKAGE) $(PYFORT_LIBS)',\
	library_directories='. $(PYFORT_LIBDIRS)',\
	compiler_options='-I../.. $(PYFORT_INCDIRS)',\
	python_directory='$(top_srcdir)/lib',\
	package_name='$(PACKAGE)',\
	generate_as='$(PYTHON_PACKAGE)',\
	freeform=1)" > $(PACKAGE).pfp

endif

