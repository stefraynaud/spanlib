################################
## MAKEFILE FOR DOCUMENTATION ##
################################
TARGET=pcamssa
TOPDIR=..
SRCDIR=$(TOPDIR)/src
DOCDIR=$(TOPDIR)/doc
EXADIR=$(TOPDIR)/example
SCRIPTDIR=$(TOPDIR)/scripts
INSTALLDIR=$(HOME)
INSTALLDOCDIR=$(INSTALLDIR)/doc/$(TARGET)
GENDOC=$(SCRIPTDIR)/gendoc.pl
DOC_MAIN=$(DOCDIR)/fortran_main.html
DOC_CODE=$(DOCDIR)/fortran_code.html
DOC_CSS=$(DOCDIR)/fortran.css
DOC_TEMPLATES=$(DOCDIR)/fortran_top.html $(DOCDIR)/fortran_mid.html $(DOCDIR)/fortran_bottom.html
SRC_MAIN=$(SRCDIR)/$(TARGET).f90
SRC_EXAMPLE=$(EXADIR)/example.f90


######################################################################

.SUFFIXES: .html

all: install

install: install-doc

install-doc: $(DOC_MAIN) $(DOC_CODE) $(DOC_CSS)
	install -d --mode=755 $(INSTALLDOCDIR)
	install --mode=644 $(DOC_MAIN) $(DOC_CODE) $(DOC_CSS) $(INSTALLDOCDIR)

doc: $(DOC_MAIN) $(DOC_CODE)

$(DOC_MAIN): $(GENDOC) $(DOC_TEMPLATES) $(SRC_MAIN) $(SRC_EXAMPLE)
	$^ $(DOC_MAIN) $(DOC_CODE)

$(DOC_CODE): $(DOC_MAIN)

clear:
	rm -f $(DOC_MAIN) $(DOC_CODE)









