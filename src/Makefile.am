######################################################################
## SpanLib, Raynaud 2006-2007
######################################################################

lib_LIBRARIES = libspanlib.a
if UPPER_CASE_MOD
nodist_include_HEADERS=SPANLIB_PRECISION.mod SPANLIB_LAPACK95.mod SPANLIB.mod
else
nodist_include_HEADERS=spanlib_precision.mod spanlib_lapack95.mod spanlib.mod
endif
nodist_libspanlib_a_SOURCES=spanlib_precision.f90 spanlib_lapack95.f90
dist_libspanlib_a_SOURCES=spanlib.f90
noinst_PROGRAMS = getf90precision
getf90precision_SOURCES = getf90precision.f90
CLEANFILES=*.pyf setup.cfg *.mod $nodist_libspanlib_a_SOURCES
AM_FCFLAGS=$(LAPACK95_INC)
AM_F90FLAGS=$(AM_FCFLAGS)

WRAPPER=spanlib_pywrap.f90
if WITH_PYTHON
PYTHON_LIBRARY_DIR = $(srcdir)/build/$(F2PY_BUILD_DIR)/spanlib
PYTHON_LIBRARY = $(PYTHON_LIBRARY_DIR)/spanlib_fort.so
endif

EXTRA_DIST=$(WRAPPER)

######################################################################
######################################################################

all-local: check-precision $(LIBRARIES) $(HEADERS) $(PYTHON_LIBRARY)

if UPPER_CASE_MOD
spanlib.o: SPANLIB_LAPACK95.mod spanlib_precision.mod
SPANLIB_PRECISION.mod: spanlib_precision.o
SPANLIB_LAPACK95.mod: spanlib_lapack95.o
else
spanlib.o: spanlib_lapack95.mod spanlib_precision.mod
spanlib_precision.mod: spanlib_precision.o
spanlib_lapack95.mod: spanlib_lapack95.o
endif
spanlib_precision.f90: check-precision

# spanlib_lapack95.f90: check-lapack95

# check-precision: check-lapack95 spanlib_precision.F90
check-precision: spanlib_precision.F90 spanlib_lapack95.f90
	@echo "Checking precision..."
	test -f spanlib_precision.f90 || cp spanlib_precision.F90 spanlib_precision.f90
	@precision=`grep "wp *= *sp" spanlib_precision.f90 | cut -d\  -f6` ; \
	test "$$precision" = "$(WP)" || \
	( rm -f spanlib_precision.f90 && sed "s/wp = sp/wp = $(WP)/" $^ > spanlib_precision.f90 && \
	echo "Switched to $(PRECISION) precision" )

# check-lapack95: spanlib_lapack95.F90
spanlib_lapack95.f90: spanlib_lapack95.F90
	@echo "Updating lapack95 config..."
	@test -f spanlib_lapack95.f90  || cp $^ spanlib_lapack95.f90
	@if test "$(USE_MKL)" = "yes" ; then \
      rm -f spanlib_lapack95.f90 ; \
      sed -e "s/use la_precision/use mkl95_precision/" \
          -e "s/use f95_lapack/use mkl95_lapack/" \
          $^ > spanlib_lapack95.f90 ; \
      echo "Switched to MKL (intel) lapack95 config" ; \
   else  \
      rm -f spanlib_lapack95.f90 ; \
      sed "s/ => [_a-z]*//g" $^ > spanlib_lapack95.f90 ; \
      echo "Switched back to default lapack95 config" ; \
   fi
 	
# check-lapack95: spanlib_lapack95.F90
# 	@echo "Checking lapack95 config..."
# 	@test -f spanlib_lapack95.f90  || cp $^ spanlib_lapack95.f90
# 	@use_mkl_lapack95=`grep "use .*mkl95_lapack" spanlib_lapack95.f90`  ; \
# 	if test "x$$use_mkl_lapack95" = "x" -a "$(USE_MKL_LAPACK95)" = "yes" ; then \
#       rm -f spanlib_lapack95.f90 ; \
#       sed -e "s/use la_precision/use mkl95_precision/" \
#           -e "s/use f95_lapack/use mkl95_lapack/" \
#           $^ > spanlib_lapack95.f90 ; \
#       echo "Switched to MKL (intel) lapack95 config" ; \
#    elif test "x$$use_mkl_lapack95" != "x" -a "$(USE_MKL_LAPACK95)" = "no" ; then \
#       rm -f spanlib_lapack95.f90 ; \
#       sed "s/ => [_a-z]*//g" $^ > spanlib_lapack95.f90 ; \
#       echo "Switched back to default lapack95 config" ; \
#    fi

if WITH_PYTHON
install-exec-hook:
if LOCAL_PYTHON_INSTALL
	test -z "$(PYTHONDIR)" || $(mkdir_p) $(PYTHONDIR)
	@echo "Installing $(PACKAGE) python module in $(PYTHONDIR)"
	$(INSTALL) -d -m 755 $(PYTHONDIR)/spanlib
	$(INSTALL_DATA) $(PYTHON_LIBRARY_DIR)/* $(PYTHONDIR)/spanlib
else
	$(PYTHON) setup.py install
endif

uninstall-hook:
	@echo "Uninstalling $(PACKAGE) python module"
	-rm -rf $(PYTHONDIR)/spanlib*

clean-local:
	-test ! -f setup.py || $(PYTHON) setup.py clean
	-rm -rf build $(PACKAGE).pyf setup.py 

$(PYTHON_LIBRARY): $(HEADERS) setup.py setup.cfg $(PACKAGE).pyf $(WRAPPER) \
   $(LIBRARIES) $(top_srcdir)/lib/__init__.py $(top_srcdir)/lib/spanlib_python.py
	$(PYTHON) setup.py build
	@echo "Adding tmp lib dir"
	-rm -f build/tmp_lib ; ln -s $(F2PY_BUILD_DIR) build/tmp_lib

$(PACKAGE).pyf: $(top_srcdir)/scripts/genpyf.pl $(WRAPPER) getf90precision
	@echo Creating $@
	$(top_srcdir)/scripts/genpyf.pl $(WRAPPER) $@ $(PRECISION) `./getf90precision`

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

