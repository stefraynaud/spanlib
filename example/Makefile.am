######################################################################
## SpanLib, Raynaud 2006
######################################################################

# Fortran examples
FORTRAN_EXAMPLES=fortran_pcamssa fortran_pcasvd fortran_svdmodel
if ENABLE_FORTRAN_EXAMPLES
noinst_PROGRAMS = $(FORTRAN_EXAMPLES)
endif
SRC_F90=$(top_srcdir)/src/spanlib.f90
LOCAL_LIB_F90=$(top_srcdir)/src/libspanlib.a
fortran_pcamssa_SOURCES = fortran_pcamssa.f90
fortran_pcasvd_SOURCES = fortran_pcasvd.f90
fortran_svdmodel_SOURCES = fortran_svdmodel.f90
AM_FCFLAGS = -I$(top_srcdir)/src $(NETCDF_INC) $(LAPACK95_INC)
AM_LDFLAGS = -L$(top_srcdir)/src $(NETCDF_LIB) $(LAPACK95_LIB) $(LAPACK_LIB) $(BLAS_LIB)
LDADD = -lnetcdf $(LOCAL_LIB_F90) $(LAPACK95) $(LAPACK) $(BLAS)

# Python examples
PYTHON_EXAMPLES=python_pcamssa_comp.py python_pcamssa_2regions.py python_mssa_mctest.py
LOCAL_LIB_PY=$(top_srcdir)/src/build/tmp_lib/spanlib/spanlib_fort.so
SRC_PY=$(SRC_F90) $(top_srcdir)/lib/spanlib_python.py $(top_srcdir)/src/spanlib_pywrap.f90

URL = http://stefdeperou.free.fr/pub
NCINPUT = data2.cdf

######################################################################
if ENABLE_FORTRAN_EXAMPLES
define variables_fortran_pcamssa
echo "# - orig:  original field (before its analysis)"; \
echo "# - reco1: reconstruction of the first pair of mssa modes"; \
echo "# - pair1: phases composites of the first pair"
endef
define variables_fortran_pcasvd
echo "# - sst_box1/2:    left/right original fields"; \
echo "# - eof_box1/2:    left/right SVD EOFs"; \
echo "# - pc_box1/2:     left/right SVD PCs"; \
echo "# - sstrec_box1/2: left/right reconstructed fields"; \
echo "# - ev:            eigen values"
endef
define variables_fortran_svdmodel
echo "# - sst_box1/2:      left/right original fields"; \
echo "# - sst_pca_box1/2:  left/right SST from pre-PCA"; \
echo "# - sst_model_box2:  right SST from SVD model"
endef
endif

######################################################################
EXTRA_DIST = $(PYTHON_EXAMPLES)
CLEANFILES = output_*.nc


##################
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
if ENABLE_FORTRAN_EXAMPLES
fortran: $(FORTRAN_EXAMPLES)
else
fortran:
	@echo "###########################################################"
	@echo "# Your configuration has no f90 netcdf support so you can't"
	@echo "# run the example. If you wish to do it, use ./configure"
	@echo "# options --with-netcdf-inc and --with-netcdf-lib or"
	@echo "# or associated variables."
	@echo "###########################################################"
endif

fortran_%: $(NCINPUT) $(LOCAL_LIB_F90) output_fortran_%.nc
	@echo "###########################################################"
	@echo "# Results are stored in the following netcdf file:"
	@echo "# $(word 3,$^)"
	@echo "###########################################################"


output_fortran_%.nc: fortran_%
	@rm -f $@
	@echo "###########################################################"
	@echo "# Running the F90 $^ program..."
	@echo "###########################################################"
	-./$^
	-@test -f $@ || echo "No netcdf file created"

$(LOCAL_LIB_F90): $(SRC_F90)
	cd $(top_srcdir)/src && make


####################
## Python example ##
####################

if ENABLE_PYTHON
python: $(PYTHON_EXAMPLES)
$(PYTHON_EXAMPLES): $(NCINPUT) $(LOCAL_LIB_F90) $(LOCAL_LIB_PY)
	@echo "###########################################################"
	@echo "# Running the python example `echo $@ |cut -c7` program..."
	@echo "###########################################################"
	$(PYTHON) example`echo $@ |cut -c7`.py
$(LOCAL_LIB_PY): $(SRC_PY)
	cd $(top_srcdir)/src && make
else
$(PYTHON_EXAMPLES):
	@echo "###########################################################"
	@echo "# You can't run python examples with your configuration"
	@echo "###########################################################"
endif



