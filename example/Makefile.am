######################################################################
## SpanLib, Raynaud 2006
######################################################################

if WITH_EXAMPLE
noinst_PROGRAMS = example1 example2
example1_SOURCES = example1.f90
example2_SOURCES = example2.f90
endif

LOCAL_LIB=$(top_srcdir)/src/libspanlib.a
AM_FCFLAGS = -I$(top_srcdir)/src $(NETCDF_INC) $(LAPACK_INC)
##example1_LDFLAGS = -L$(top_srcdir)/src $(NETCDF_LIB) $(LAPACK_LIB) $(BLAS_LIB)
##example2_LDFLAGS = $(example1_LDFLAGS)
AM_LDFLAGS = -L$(top_srcdir)/src $(NETCDF_LIB) $(LAPACK_LIB) $(BLAS_LIB)
LDADD = -lnetcdf $(LOCAL_LIB) $(LAPACK95) $(LAPACK) $(BLAS)

URL = http://stefdeperou.free.fr/pub
NCINPUT = data2.cdf
F90NCOUTPUT1 = output1.nc
F90NCOUTPUT2 = output2.nc
if WITH_PYTHON_EXAMPLE
VARIABLE1 = pair1
VARIABLE2 = leof1
CLEANFILES = $(NCOUTPUT)
endif

EXTRA_DIST = example1.py example2.py

######################################################################
######################################################################

###################
## common config ##
###################

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

#################
## F90 example ##
#################
if WITH_EXAMPLE
all-local: $(NCINPUT) fortran1
else
all-local:
	@echo "###########################################################"
	@echo "# Your configuration has no f90 netcdf support so you can't"
	@echo "# run the example. If you wish to do it, use ./configure"
	@echo "# options --with-netcdf-inc and --with-netcdf-lib or"
	@echo "# or associated variables."
	@echo "###########################################################"
endif

fortran1: $(F90NCOUTPUT1)
if HAS_NCVIEWER
	@test -f $^ &&  ( \
	echo "###########################################################" ; \
	echo "# Visualisation of output netcdf file $^"; \
	echo "# This file contains the following variables of interest:"; \
	echo "# - orig:  original field (before its analysis)"; \
	echo "# - reco1: reconstruction of the first pair of mssa modes"; \
	echo "# - pair1: phases composites of the first pair"; \
	echo "###########################################################" )
	test -f $^ && $(NCVIEWER) $(NCVIEWER_ARGS1)
else
	@echo "###########################################################"
	@echo "# Results are stored in the following netcdf file:"
	@echo "# $^"
	@echo "###########################################################"
endif

fortran2: $(F90NCOUTPUT2)
if HAS_NCVIEWER
	@test -f $^ &&  ( \
	echo "###########################################################" ; \
	echo "# Visualisation of output netcdf file $^"; \
	echo "# This contains the following variables of interest:"; \
	echo "# - lmodes1/rmode1: first left/right reconstructed mode"; \
	echo "# - leof1/reof1: first left/right EOF"; \
	echo "# - lpc1/rpc1: first left/right PC"; \
	echo "# - ev: eigen values"; \
	echo "###########################################################" )
	test -f $^ && $(NCVIEWER) $(NCVIEWER_ARGS2)
else
	@echo "###########################################################"
	@echo "# Results are stored in the following netcdf file:"
	@echo "# $^"
	@echo "###########################################################"
endif

$(F90NCOUTPUT1): $(NCINPUT) $(LOCAL_LIB) example1
	@rm -f $(F90NCOUTPUT1)
	@echo "###########################################################"
	@echo "# Running the F90 example1 program..."
	@echo "###########################################################"
	./example1
	@test -f $(F90NCOUTPUT1) || echo "No netcdf file created"

$(F90NCOUTPUT2): $(NCINPUT) $(LOCAL_LIB) example2
	@rm -f $(F90NCOUTPUT2)
	@echo "###########################################################"
	@echo "# Running the F90 example2 program..."
	@echo "###########################################################"
	./example2
	@test -f $(F90NCOUTPUT2) || echo "No netcdf file created"


$(LOCAL_LIB):
	cd ../src && make


####################
## Python example ##
####################

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

