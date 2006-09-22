######################################################################
## SpanLib, Raynaud 2006
######################################################################

if WITH_EXAMPLE
noinst_PROGRAMS = example1 example2 example3
example1_SOURCES = example1.f90
example2_SOURCES = example2.f90
example3_SOURCES = example3.f90
endif

LOCAL_LIB=$(top_srcdir)/src/libspanlib.a
AM_FCFLAGS = -I$(top_srcdir)/src $(NETCDF_INC) $(LAPACK_INC)
##example1_LDFLAGS = -L$(top_srcdir)/src $(NETCDF_LIB) $(LAPACK_LIB) $(BLAS_LIB)
##example2_LDFLAGS = $(example1_LDFLAGS)
AM_LDFLAGS = -L$(top_srcdir)/src $(NETCDF_LIB) $(LAPACK_LIB) $(BLAS_LIB)
LDADD = -lnetcdf $(LOCAL_LIB) $(LAPACK95) $(LAPACK) $(BLAS)

URL = http://stefdeperou.free.fr/pub
NCINPUT = data2.cdf
F90NCOUTPUT1 = output_fortran1.nc
F90NCOUTPUT2 = output_fortran2.nc
F90NCOUTPUT3 = output_fortran3.nc
if WITH_PYTHON_EXAMPLE
VARIABLE1 = pair1
VARIABLE2 = eof_box1
VARIABLE3 = sst_model_box2
CLEANFILES = $(NCOUTPUT)
endif

EXTRA_DIST = example1.py example2.py

VIEW_HEADER=echo "###########################################################" ; \
	echo "# Visualisation of output netcdf file $^"; \
	echo "# This file contains the following variables of interest:"
VIEW_FOOTER=echo "###########################################################"
VIEW_MISSING=echo "###########################################################" ; \
echo "# Results are stored in the following netcdf file:" ;\
echo "# $^" ;\
echo "###########################################################"
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
all-local: fortran1
else
all-local:
	@echo "###########################################################"
	@echo "# Your configuration has no f90 netcdf support so you can't"
	@echo "# run the example. If you wish to do it, use ./configure"
	@echo "# options --with-netcdf-inc and --with-netcdf-lib or"
	@echo "# or associated variables."
	@echo "###########################################################"
endif

fortran1: $(NCINPUT) $(LOCAL_LIB) $(F90NCOUTPUT1)
if HAS_NCVIEWER
	@test -f $^ &&  ( \
	$(VIEW_HEADER) ; \
	echo "# - orig:  original field (before its analysis)"; \
	echo "# - reco1: reconstruction of the first pair of mssa modes"; \
	echo "# - pair1: phases composites of the first pair"; \
	$(VIEW_FOOTER) )
	test -f $^ && $(NCVIEWER) $(NCVIEWER_ARGS1)
else
	@$(VIEW_MISSING)
endif

fortran2: $(NCINPUT) $(LOCAL_LIB) $(F90NCOUTPUT2)
if HAS_NCVIEWER
	@test -f $^ &&  ( \
	$(VIEW_HEADER) ; \
	echo "# - sst_box1/2: left/right original fields"; \
	echo "# - eof_box1/2: left/right SVD EOFs"; \
	echo "# - pc_box1/2:  left/right SVD PCs"; \
	echo "# - ev: eigen values"; \
	$(VIEW_FOOTER) )
	test -f $^ && $(NCVIEWER) $(NCVIEWER_ARGS2)
else
	@$(VIEW_MISSING)
endif

fortran3: $(NCINPUT) $(LOCAL_LIB) $(F90NCOUTPUT3)
if HAS_NCVIEWER
	@test -f $^ &&  ( \
	$(VIEW_HEADER) ; \
	echo "# - sst_box1/2:        left/right original fields"; \
	echo "# - sst_pca_box1/2:    left/right SST from pre-PCA"; \
	echo "# - sst_model_box1/2:  right SST from SVD model"; \
	$(VIEW_FOOTER) )
	test -f $^ && $(NCVIEWER) $(NCVIEWER_ARGS3)
else
	@$(VIEW_MISSING)
endif

output_fortran%.nc: example%
	@rm -f $@
	@echo "###########################################################"
	@echo "# Running the F90 $^ program..."
	@echo "###########################################################"
	./$^
	@test -f $@ || echo "No netcdf file created"

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

