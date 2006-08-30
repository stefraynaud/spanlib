######################################################################
## SpanLib, Raynaud 2006
######################################################################

EXTRA_DIST = COPYING INSTALL AUTHORS README ChangeLog
SUBDIRS = src lib scripts doc
DIST_SUBDIRS = $(SUBDIRS) example

doc: documentation

documentation:
	cd $(top_srcdir)/doc && make

example:
	cd $(top_srcdir)/example && make


