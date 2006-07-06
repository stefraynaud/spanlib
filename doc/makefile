################################
## MAKEFILE FOR DOCUMENTATION ##
################################
TOPDIR=..
SRCDIR=$(TOPDIR)/src
DOCDIR=$(TOPDIR)/doc
EXADIR=$(TOPDIR)/example
SCRIPTDIR=$(TOPDIR)/scripts
HTMLDIR=$(DOCDIR)/html
XMLDIR=$(DOCDIR)/xml
##INSTALLDIR=$(HOME)
##INSTALLDOCDIR=$(INSTALLDIR)/doc/$(TARGET)
GENDOC=$(SCRIPTDIR)/gendoc.pl
XML_LIBRARY=$(XMLDIR)/f90_library.xml
HTML_LIBRARY=$(HTMLDIR)/f90_library.html
F90_LIBRARY=$(SRCDIR)/spanlib.f90
XML_EXAMPLE=$(XMLDIR)/f90_example.xml
HTML_EXAMPLE=$(HTMLDIR)/f90_example.html
F90_EXAMPLE=$(EXADIR)/example.f90
HTML=$(HTMLDIR)/index.html
CSS=$(HTMLDIR)/spanlib.css
XSL=$(XMLDIR)/spanlib.xsl
XML=$(XMLDIR)/spanlib.xml
##XSLTARGS=--xinclude --stringparam  html.stylesheet "$(CSS)"
XSLTARGS=--xinclude
DIRS=$(XMLDIR) $(HTMLDIR)
XMLS=$(XML) $(XML_LIBRARY) $(XML_EXAMPLE)

######################################################################

.SUFFIXES: .html .xml .xsl .css

.PHONY: html

$(HTML): $(XSL) $(CSS) $(XMLS) $(HTMLDIR)
	xsltproc $(XSLTARGS) --stringparam  html.stylesheet "'$(CSS)'" -o $(HTML_LIBRARY) $(XSL) $(XML_LIBRARY)
	xsltproc $(XSLTARGS) --stringparam  html.stylesheet "'$(CSS)'" -o $(HTML_EXAMPLE) $(XSL) $(XML_EXAMPLE)
	xsltproc $(XSLTARGS) -o $(HTML) $(XSL) $(XML)

$(XMLS): $(GENDOC) $(XMLDIR) $(F90_LIBRARY) $(F90_EXAMPLE)
	$^

$(DIRS):
	mkdir -p $@

doc: $(HTML)

all: doc

clear:
	rm -f $(XML) $(DOC_CODE)



######################################################################
##	xsltproc --docbook -o index.html --stringparam  html.stylesheet "$(DOC_CSS)" spanlib.xsl template.xml

##install: install-doc

##install-doc: $(DOC_MAIN) $(DOC_CODE) $(DOC_CSS)
##	install -d --mode=755 $(INSTALLDOCDIR)
##	install --mode=644 $(DOC_MAIN) $(DOC_CODE) $(DOC_CSS) $(INSTALLDOCDIR)







