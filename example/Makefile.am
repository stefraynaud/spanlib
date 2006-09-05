######################################################################
## SpanLib, Raynaud 2006
######################################################################

if WITH_EXAMPLE
noinst_PROGRAMS = example
example_SOURCES = example.f90
endif

LOCAL_LIB=$(top_srcdir)/src/libspanlib.a
AM_FCFLAGS = -I$(top_srcdir)/src $(NETCDF_INC)
example_LDFLAGS = -L$(top_srcdir)/src $(NETCDF_LIB) $(LAPACK_LIB) $(BLAS_LIB)
example_LDADD = -lnetcdf $(LOCAL_LIB) $(LAPACK) $(BLAS)

URL = http://stefdeperou.free.fr/pub
NCINPUT = data2.cdf
F90NCOUTPUT = output.nc
if WITH_PYTHON_EXAMPLE
VARIABLE = pair1
CLEANFILES = $(NCOUTPUT)
endif

EXTRA_DIST = example1.py example2.py

######################################################################
######################################################################

## F90 example

if WITH_EXAMPLE
all-local: $(NCDATA) view
else
all-local:
	@echo "###########################################################"
	@echo "# Your configuration has no f90 netcdf support so you can't
	@echo "# run the example. If you wish to do it, use ./configure
	@echo "# options --with-netcdf-inc and --with-netcdf-lib or
	@echo "# or associated variables."
	@echo "###########################################################"
endif

$(NCINPUT):
if HAS_DOWNLOADER
	@echo "###########################################################"
	@echo "# Download of input netcdf file:"
	@echo "# $(URL)/$@"
	@echo "###########################################################"
if USE_LINKS
	$(DOWNLOADER) -source 1 $(URL)/$@ > $<
else
	$(DOWNLOADER) $(URL)/$@
endif
else
	@echo "###########################################################"
	@echo "# You must first download the netcdf data file at url:"
	@echo "# $(URL)/$(NCINPUT)"
	@echo "###########################################################"
endif

view: $(F90NCOUTPUT)
if HAS_NCVIEWER
	@echo "###########################################################"
	@echo "# Visualisation of output netcdf file $(F90NCOUTPUT)"
	@echo "# This contains the following variables of interest:"
	@echo "# - orig:  original field (before its analysis)"
	@echo "# - reco1: reconstruction of the first pair of mssa modes"
	@echo "# - pair1: phases composites of the first pair"
	@echo "###########################################################"
	$(NCVIEWER) $(F90NCOUTPUT) $(VARIABLE)
else
	@echo "###########################################################"
	@echo "# Results are stored in the following netcdf file:"
	@echo "# $(F90NCOUTPUT)"
	@echo "###########################################################"
endif

$(F90NCOUTPUT): $(NCINPUT) $(LOCAL_LIB) $(PROGRAMS)
	@echo "###########################################################"
	@echo "# Running the F90 example program..."
	@echo "###########################################################"
	./example


$(LOCAL_LIB):
	cd ../src && make


## Python example

if WITH_PYTHON_EXAMPLE
python1 python2: $(NCDATA)
	@echo "###########################################################"
	@echo "# Running the python example `echo $@ |cut -c7` program..."
	@echo "###########################################################"
	$(PYTHON) example`echo $@ |cut -c7`.py
else
python1 python2:
	@echo "###########################################################"
	@echo "# You can't run python examples with your configuration"
	@echo "###########################################################"
endif

