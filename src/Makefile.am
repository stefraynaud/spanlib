######################################################################
## SpanLib, Raynaud 2006
######################################################################

lib_LIBRARIES = libspanlib.a
nodist_include_HEADERS=spanlib.mod
libspanlib_a_SOURCES=spanlib.f90
CLEANFILES=*.pyf *.pfp

WRAPPER = $(PACKAGE)_pywrap.f90
PYTHON_PACKAGE = $(PACKAGE)_fort
if WITH_PYTHON
PYTHON_LIBRARY_DIR = $(srcdir)/build/$(PYFORT_BUILD_DIR)/$(PACKAGE)
PYTHON_LIBRARY = $(PYTHON_LIBRARY_DIR)/$(PACKAGE)_fort.so
endif

######################################################################
######################################################################

all-local: $(LIBRARIES) $(HEADERS) $(PYTHON_LIBRARY)

##fortran_library: $(LIBRARIES) $(HEADERS)

if WITH_PYTHON
install-exec-hook:
	@$(NORMAL_INSTALL)
	test -z "$(PYTHONDIR)" || $(mkdir_p) $(PYTHONDIR)/spanlib_dir
	@echo "Installing $(PACKAGE) python module in $(PYTHONDIR)"
	cp $(PYTHON_LIBRARY_DIR)/* $(PYTHONDIR)/spanlib_dir
	echo spanlib_dir > $(PYTHONDIR)/spanlib_dir.pth

uninstall-hook:
	@echo "Uninstalling $(PACKAGE) python module"
	-rm -rf $(PYTHONDIR)/spanlib_dir $(PYTHONDIR)/spanlib_dir.pth
clean-local:
	-rm -rf build

$(PYTHON_LIBRARY): $(PACKAGE).pfp $(WRAPPER) $(PACKAGE).pyf $(LIBRARIES) $(HEADERS)
	$(PYFORT) $(PYFORT_BUILD_OPT) $(PACKAGE).pfp

$(PACKAGE).pyf: $(WRAPPER)
	@echo Creating $(PACKAGE).pyf
	$(top_srcdir)/scripts/genpyf.pl $(WRAPPER) $(PACKAGE).pyf

$(PACKAGE).pfp:
	@echo Creating $(PACKAGE).pfp
	@echo "pyf('$(PACKAGE).pyf',\
	sources=['$(WRAPPER)'],\
	libraries='$(PACKAGE) $(PYFORT_LIBS)',\
	library_directories='. $(PYFORT_DIRS)',\
	compiler_options='-I.',\
	python_directory='$(top_srcdir)/lib',\
	package_name='$(PACKAGE)',\
	generate_as='$(PYTHON_PACKAGE)',\
	freeform=1)" > $(PACKAGE).pfp

endif

