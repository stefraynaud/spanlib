######################################################################
## SpanLib, Raynaud 2006
######################################################################

lib_LIBRARIES = libspanlib.a
nodist_include_HEADERS=spanlib_precision.mod spanlib.mod
nodist_libspanlib_a_SOURCES=spanlib_precision.f90
dist_libspanlib_a_SOURCES=spanlib.f90
noinst_PROGRAMS = getf90precision
getf90precision_SOURCES = getf90precision.f90
CLEANFILES=*.pyf setup.py setup.cfg *.mod spanlib_precision.f90
EXTRA_DIST=spanlib_pywrap.f90
AM_FCFLAGS = $(LAPACK95_INC)

WRAPPER = $(PACKAGE)_pywrap.f90
if WITH_PYTHON
PYTHON_LIBRARY_DIR = $(srcdir)/build/$(F2PY_BUILD_DIR)/$(PACKAGE)
PYTHON_LIBRARY = $(PYTHON_LIBRARY_DIR)/$(PACKAGE)_fort.so
endif

######################################################################
######################################################################

all-local: check-precision $(LIBRARIES) $(HEADERS) $(PYTHON_LIBRARY)

spanlib.o: spanlib_precision.mod

spanlib_precision.mod: spanlib_precision.o

spanlib_precision.f90: check-precision

check-precision: spanlib_precision.F90
	@echo "Checking precision..."
	@test -f spanlib_precision.f90 || cp $^ spanlib_precision.f90
	@precision=`grep "wp *= *sp" spanlib_precision.f90 | cut -d\  -f6` ; \
	test "$$precision" = "$(WP)" || \
	( rm -f spanlib_precision.f90 && sed "s/wp *= *sp/wp = $(WP)/" $^ > spanlib_precision.f90 && \
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
	-rm -rf build $(PACKAGE).pyf setup.py setup.cfg

$(PYTHON_LIBRARY): $(HEADERS) setup.py setup.cfg $(PACKAGE).pyf $(WRAPPER) $(LIBRARIES) $(top_srcdir)/lib/__init__.py $(top_srcdir)/lib/spanlib_python.py
	$(PYTHON) setup.py install

$(PACKAGE).pyf: $(top_srcdir)/scripts/genpyf.pl $(WRAPPER) getf90precision
	@echo Creating $@
	$(top_srcdir)/scripts/genpyf.pl $(WRAPPER) $@ $(PRECISION) `getf90precision`

setup.py: FORCE
	@test -f $@ || ( cp template.$@ $@ && echo "$@ created" )
	@version=`grep "AM_INIT_AUTOMAKE" $(top_srcdir)/configure.in | cut -d"," -f2 | cut -d")" -f1` ; \
	test "x$$version" = "x`grep version=\'.*\' $@ | cut -d\' -f2`" -a \
	"x$(F2PY_LIBS)" = "x`grep libs= $@ | cut -f1 | cut -d= -f2`" -a \
	"x$(F2PY_LIBDIRS)" = "x`grep libDirs= $@ | cut -f1 | cut -d= -f2`" -a \
	"x$(F2PY_INCDIRS)" = "x`grep incDirs= $@ | cut -f1 | cut -d= -f2`"|| \
	( rm -f $@ && sed \
	-e "s/_VERSION_/$$version/" \
	-e "s/_LIBS_/$(F2PY_LIBS)/" \
	-e "s=_LIBDIRS_=$(F2PY_LIBDIRS)=" \
	-e "s=_INCDIRS_=$(F2PY_INCDIRS)=" \
	template.$@ > $@  && \
	echo "$@ updated" )

setup.cfg: FORCE
	@test -f $@ || ( cp template.$@ $@ && echo "$@ created" )
	@test x$(F2PY_FC_VENDOR) = x`grep "^fcompiler=" $@ | cut -f1 | cut -d= -f2` || \
	( rm -f $@ && \
	sed "s/_FCVENDOR_/$(F2PY_FC_VENDOR)/" template.$@ > $@  && \
	echo "$@ updated" )

FORCE:

endif

