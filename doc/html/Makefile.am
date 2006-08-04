######################################################################
## SpanLib, Raynaud 2006
######################################################################

docdir=$(datadir)/@PACKAGE@

HTML_F90_LIB = src_f90_lib.html
HTML_PYT_MOD = src_pyt_mod.html
HTML_F90_EXA = src_f90_exa.html
HTML_PYT_EXA = src_pyt_exa.html

XSL = ../xml/spanlib.xsl
XSLSRC = ../sml/src.xsl

dist_doc_DATA = index.html spanlib.css \
	doc_ins.html doc_f90.html doc_pyt.html doc_bib.html \
	$(HTML_F90_LIB) $(HTML_PYT_MOD) $(HTML_F90_EXA) $(HTML_PYT_EXA)


######################################################################
######################################################################

if HAS_DOC_SUPPORT
else
src_%.html: ../html/src_%.xml $(XSLSRC)
	xsltproc --xinclude -o $@ $(XSLSRC) $<
endif

../html/src_%.xml:
	cd ../html && make
