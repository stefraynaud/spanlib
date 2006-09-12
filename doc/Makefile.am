######################################################################
## SpanLib, Raynaud 2006
######################################################################

## Source files
SOURCES = $(top_srcdir)/src/spanlib.f90 $(top_srcdir)/example/example.f90 $(top_srcdir)/lib/spanlib_python.py $(top_srcdir)/example/example1.py $(top_srcdir)/example/example2.py

## XML translation of source files
XML_SRC = src_f90_lib.xml src_f90_exa.xml src_pyt_mod.xml src_pyt_ex1.xml src_pyt_ex2.xml

## XML translation of source files to be included
XML_SRC_INC = $(XML_SRC:.xml=_inc.xml)

## HTML translation of source files
HTML_SRC = $(XML_SRC:.xml=.html)

## XML list of subroutines
XML_DOC_SUB = doc_f90_sub_inc.xml doc_pyt_sub_inc.xml

## Vurrent version
XML_DOC_VER_INC = doc_ver_inc.xml

## Source parsing command
PARSE = $(top_srcdir)/scripts/gendoc.pl $(SOURCES)

## For the distro
EXTRA_DIST = spanlib.xml spanlib.xsl spanlib_single.xsl src.xsl common.xsl $(XML_SRC)

## For the installation
docdir = $(datadir)/@PACKAGE@/html
dist_doc_DATA = \
	$(HTML_SRC) index.html spanlib.css spanlib.html \
	doc_ins.html doc_f90.html doc_pyt.html doc_lin.html doc_bib.html \
	doc_ins_dow.html doc_ins_ins.html \
	doc_f90_sub.html doc_f90_exa.html \
	doc_pyt_sub.html doc_pyt_exa.html

if HAS_DOC_SUPPORT
CLEANFILES = *_inc.xml
endif

######################################################################
######################################################################

if HAS_DOC_SUPPORT
spanlib.html: spanlib.xml spanlib_single.xsl common.xsl $(XML_SRC_INC) $(XML_DOC_SUB) spanlib.css $(XML_DOC_VER_INC)
	@echo "Generating single page html documentation..."
	@$(XSLTPROC) $(XSLTPROC_FLAGS) -o $@ spanlib_single.xsl $<

index.html: spanlib.xml spanlib.xsl common.xsl $(XML_SRC_INC) $(XML_DOC_SUB) $(XML_DOC_VER_INC)
	@echo "Generating mutiple page html documentation..."
	@$(XSLTPROC) $(XSLTPROC_FLAGS) -o $@ spanlib.xsl $<

src_%.html: src_%.xml src_%_inc.xml src.xsl common.xsl
	@echo "Generating $@..."
	@$(XSLTPROC) $(XSLTPROC_FLAGS) -o $@ src.xsl $<

$(XML_SRC_INC) $(XML_DOC_SUB): $(PARSE)
	@echo "Generating xml documentations from code sources..."
	@$^

$(XML_DOC_VER_INC): $(top_srcdir)/configure.in
	@echo "Checking package version..."
	@test -f $@ || echo "" > $@
	@version=`grep "AM_INIT_AUTOMAKE" $^ | cut -d"," -f2 | cut -d")" -f1` ; \
	test "x$$version" = "x`grep Version $@ | cut -d\  -f2 | cut -d'<' -f1`" ;\
	rm -f $@ && echo '<para>Version '$$version'</para>' > $@ ; \
	echo "Version updated"

else
endif
