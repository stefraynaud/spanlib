######################################################################
## SpanLib, Raynaud 2006
######################################################################

if WITH_EXAMPLE
noinst_PROGRAMS = example
example_SOURCES = example.f90
endif

AM_FCFLAGS = $(NETCDF_INC) -I$(top_srcdir)/src
AM_LDFLAGS = $(NETCDF_LIB) -L$(top_srcdir)/src $(LAPACK_LIB) $(BLAS_LIB)
LDADD = -lspanlib -lnetcdf -llapack -lblas

URL = http://stefdeperou.free.fr/pub
NCINPUT = data.cdf
NCOUTPUT = output.nc

CLEANFILES = $(NCOUTPUT)

######################################################################
######################################################################

if WITH_EXAMPLE
all-local: $(NCDATA) view
else
all-local:
	@echo "###########################################################"
	@echo "# Your configuration has to f90 netcdf support so you can't
	@echo "# run the example. If you wish to do it, use ./configure
	@echo "# options --with-netcdf-inc and --with-netcdf-lib or
	@echo "# or associated variables."
	@echo "###########################################################"
endif

$(NCINPUT):
if HAS_DOWNLOADER
	@echo "###########################################################"
	@echo "# Download of input netcdf file:"
	@echo "# $(URL)/$(NCINPUT)"
	@echo "###########################################################"
if USE_LINKS
	$(DOWNLOADER) -source 1 $(URL)/$(NCINPUT) > $(NCINPUT)
else
	$(DOWNLOADER) $(URL)/$(NCINPUT)
endif
else
	@echo "###########################################################"
	@echo "# You must first download the netcdf data file at url:"
	@echo "# $(URL)/$(NCINPUT)"
	@echo "###########################################################"
endif

view: $(NCOUTPUT)
if HAS_NCVIEWER
	@echo "###########################################################"
	@echo "# Visualisation of output netcdf file $(NCOUTPUT)"
	@echo "###########################################################"
	$(NCVIEW) $(NCOUTPUT)
else
	@echo "###########################################################"
	@echo "# Results are stored in the following netcdf file:"
	@echo "# $(NCOUTPUT)"
	@echo "###########################################################"
endif

$(NCOUTPUT): $(NCINPUT) $(PROGRAMS)
	@echo "###########################################################"
	@echo "# Running the example program..."
	@echo "###########################################################"
	./example
