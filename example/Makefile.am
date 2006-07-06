######################################################################
## SpanLib, Raynaud 2006
######################################################################

noinst_PROGRAMS = example
example_SOURCES = example.f90

AM_FCFLAGS=-L($top_srcdir)/src -I($top_srcdir)/src -lspanlib $(LAPACK_LIBDIR) $(LAPACK) $(BLAS_LIBDIR) $(BLAS)

URL = http://stefdeperou.free.fr/pub
NCINPUT = data2.cdf
NCOUTPUT = output.nc

######################################################################
######################################################################

if WITH_EXAMPLE
all-local: $(NCDATA) view
else
all-local:
endif

$(NCDATA):
if HAS_DOWNLOADER
	$(WGET) $(URL)/$(NCINPUT)
else
	@echo "###########################################################"
	@echo "# You must first download the netcdf data file at url:"
	@echo "# $(URL)"
	@echo "###########################################################"
endif

view: $(NCOUTPUT)
if HAS_NCVIEWVER
	$(NCVIEWER) $(NCOUTPUT)
else
	@echo "###########################################################"
	@echo "# Results are stored in the following netcdf file:"
	@echo "# $(NCOUTPUT)"
	@echo "###########################################################"
endif

$(NCOUTPUT): $(PROGRAMS)
	./$(PROGRAMS)
