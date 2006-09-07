dnl AC_SR_DOCBOOK([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl checks if xsltproc can build docbook documentation
dnl (which is possible if the catalog is set up properly
dnl I also tried checking for a specific version and type of docbook
dnl but xsltproc seemed to happily run anyway, so we can t check for that
dnl and version
dnl this macro takes inspiration from
dnl http://www.movement.uklinux.net/docs/docbook-autotools/configure.html
dnl *** Inspired from ax_check_docbook.m4 *** Raynaud, 2006
AC_DEFUN([AC_SR_DOCBOOK],
[
	XSLTPROC_FLAGS="--nonet --xinclude"
	DOCBOOK_ROOT=
	TYPE_LC=xml
	TYPE_UC=XML
	DOCBOOK_VERSION="4.2-Based Variant"

#	for i in /usr/share/sgml/docbook/stylesheet/xsl/nwalsh /usr/share/sgml/docbook/xsl-stylesheets/;
	for i in /usr/share/sgml/docbook/xsl-stylesheets/;
	do
		if test -d "$i"; then
			DOCBOOK_ROOT=$i
		fi
	done

	dnl We need xsltproc to process the test
	AC_CHECK_PROG(XSLTPROC,xsltproc,xsltproc,)

	XSLTPROC_WORKS=no
	if test -n "$XSLTPROC" -a -n "$DOCBOOK_ROOT" ; then
		AC_MSG_CHECKING([whether xsltproc docbook processing works])
		DB_FILE="$DOCBOOK_ROOT/xhtml/docbook.xsl"
		$XSLTPROC $XSLTPROC_FLAGS $DB_FILE >/dev/null 2>&1 << END
<?xml version="1.0" encoding='ISO-8859-1'?>
<!DOCTYPE book PUBLIC "-//KDE//DTD DocBook $TYPE_UC V$DOCBOOK_VERSION//EN" "/usr/share/sgml/docbook/xml-dtd-4.2-1.0-17.2/docbookx.dtd">
<article id="test">
</article>
END
		if test "$?" = 0; then
			XSLTPROC_WORKS=yes
		fi
		AC_MSG_RESULT($XSLTPROC_WORKS)
	fi

	if test "x$XSLTPROC_WORKS" = "xyes"; then
		dnl execute ACTION-IF-FOUND
		ifelse([$1], , :, [$1])
	else
		dnl execute ACTION-IF-NOT-FOUND
		ifelse([$2], , :, [$2])
	fi

	AC_SUBST(XSLTPROC_FLAGS)
	AC_SUBST(DOCBOOK_ROOT)
])