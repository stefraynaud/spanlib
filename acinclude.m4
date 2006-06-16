#####################################################################
# AC macros for spanlib
# Stephane Raynaud 2006
#####################################################################

################################################################################
# Colorisation for fun!
################################################################################
AC_DEFUN([AC_SR_SET_YELLOWINK],
[
	echo -en "\\033\\1331;31m"
])
AC_DEFUN([AC_SR_SET_REDINK],
[
	echo -en "\\033\\1331;31m"
])
AC_DEFUN([AC_SR_SET_GREENINK],
[
	echo -en "\\033\\1331;32m"
])
AC_DEFUN([AC_SR_SET_NORMALINK],
[
	echo -en "\\033\\1330;39m"
])

################################################################################
# Special messages
################################################################################
AC_DEFUN([AC_SR_HEADER],
[
	AC_SR_SET_GREENINK
	echo "################################################################################"
	echo -e "### $1"
	echo "################################################################################"
	AC_SR_SET_NORMALINK
])
AC_DEFUN([AC_SR_WARNING],
[
	AC_SR_SET_YELLOWINK
	echo "################################################################################"
	echo -e "$1"
	echo "################################################################################"
	AC_SR_SET_NORMALINK
])
AC_DEFUN([AC_SR_ERROR],
[
	AC_SR_SET_REDINK
	echo "################################################################################"
	echo -e "$1"
	echo "################################################################################"
	AC_SR_SET_NORMALINK
	exit 1
])
	



################################################################################
# Check for python version (taken from AC_PYTHON_DEVEL)
################################################################################
AC_DEFUN([AC_SR_PYTHON],
[
	AC_ARG_VAR([PYTHON_VERSION],[The installed Python
		version to use, for example '2.3'. This string 
		will be appended to the Python interpreter
		canonical name.])

	AC_PATH_PROG([PYTHON],[python[$PYTHON_VERSION]])
	if test -z "$PYTHON"; then
	   AC_MSG_ERROR([Cannot find python$PYTHON_VERSION in your system path])
	   PYTHON_VERSION=""
	fi

	#
	# Check for a version of Python >= 2.1.0
	#
	AC_MSG_CHECKING([for a version of Python >= '2.1.0'])
	ac_supports_python_ver=`$PYTHON -c "import sys, string; \
		ver = string.split(sys.version)[[0]]; \
		print ver >= '2.1.0'"`
	if test "$ac_supports_python_ver" != "True"; then
		if test -z "$PYTHON_NOVERSIONCHECK"; then
			AC_MSG_RESULT([no])
			AC_MSG_FAILURE([
This version of the AC@&t@_PYTHON_DEVEL macro
doesn't work properly with versions of Python before
2.1.0. You may need to re-run configure, setting the
variables PYTHON_CPPFLAGS, PYTHON_LDFLAGS, PYTHON_SITE_PKG,
PYTHON_EXTRA_LIBS and PYTHON_EXTRA_LDFLAGS by hand.
Moreover, to disable this check, set PYTHON_NOVERSIONCHECK
to something else than an empty string.
])
		else
			AC_MSG_RESULT([skip at user request])
		fi
	else 
		AC_MSG_RESULT([yes])
	fi
	
	#
	# if the macro parameter ``version'' is set, honour it
	#
	if test -n "$1"; then
		AC_MSG_CHECKING([for a version of Python $1])
		ac_supports_python_ver=`$PYTHON -c "import sys, string; \
			ver = string.split(sys.version)[[0]]; \
			print ver $1"`
		if test "$ac_supports_python_ver" = "True"; then
	   	   AC_MSG_RESULT([yes])
		else
			AC_MSG_RESULT([no])
			AC_MSG_ERROR([this package requires Python $1. 
If you have it installed, but it isn't the default Python
interpreter in your system path, please pass the PYTHON_VERSION 
variable to configure. See ``configure --help'' for reference.
])
			PYTHON_VERSION=""
		fi
	fi
	
])

################################################################################
# Check for perl version (http://autoconf-archive.cryp.to/ac_prog_perl_version.m4)
################################################################################
dnl @synopsis AC_PROG_PERL_VERSION(VERSION, [ACTION-IF-TRUE], [ACTION-IF-FALSE])
dnl
dnl Makes sure that perl supports the version indicated. If true the
dnl shell commands in ACTION-IF-TRUE are executed. If not the shell
dnl commands in ACTION-IF-FALSE are run. Note if $PERL is not set (for
dnl example by running AC_CHECK_PROG or AC_PATH_PROG),
dnl AC_CHECK_PROG(PERL, perl, perl) will be run.
dnl
dnl Example:
dnl
dnl   AC_PROG_PERL_VERSION(5.6.0)
dnl
dnl This will check to make sure that the perl you have supports at
dnl least version 5.6.0.
dnl
dnl @category InstalledPackages
dnl @author Dean Povey <povey@wedgetail.com>
dnl @version 2002-09-25
dnl @license AllPermissive

AC_DEFUN([AC_PROG_PERL_VERSION],[dnl
# Make sure we have perl
if test -z "$PERL"; then
AC_CHECK_PROG(PERL,perl,perl)
fi

# Check if version of Perl is sufficient
ac_perl_version="$1"

if test "x$PERL" != "x"; then
  AC_MSG_CHECKING(for perl version greater than or equal to $ac_perl_version)
  # NB: It would be nice to log the error if there is one, but we cannot rely
  # on autoconf internals
  $PERL -e "use $ac_perl_version;" > /dev/null 2>&1
  if test $? -ne 0; then
    AC_MSG_RESULT(no);
    $3
    PERL=""
  else
    AC_MSG_RESULT(ok);
    $2
  fi
else
  AC_MSG_WARN(could not find perl)
fi
])dnl


################################################################################
# Check for docbook
################################################################################
AC_DEFUN([AX_CHECK_DOCBOOK], [
# It's just rude to go over the net to build
	XSLTPROC_FLAGS=--nonet
	DOCBOOK_ROOT=""
	if test ! -f /etc/xml/catalog; then
        for i in /usr/share/sgml/docbook/stylesheet/xsl/nwalsh /usr/share/sgml/docbook/xsl-stylesheets/;
        do
                if test -d "$i"; then
                        DOCBOOK_ROOT=$i
                fi
        done

        # Last resort - try net
        if test -z "$DOCBOOK_ROOT"; then
                XSLTPROC_FLAGS=
        fi
	else
        XML_CATALOG=/etc/xml/catalog
        CAT_ENTRY_START='<!--'
        CAT_ENTRY_END='-->'
	fi

	AC_CHECK_PROG(XSLTPROC,xsltproc,xsltproc,)
	XSLTPROC_WORKS=no
	if test -n "$XSLTPROC"; then
        AC_MSG_CHECKING([whether xsltproc works])

        if test -n "$XML_CATALOG"; then
                DB_FILE="http://docbook.sourceforge.net/release/xsl/current/xhtml/docbook.xsl"
        else
                DB_FILE="$DOCBOOK_ROOT/docbook.xsl"
        fi

#        $XSLTPROC $XSLTPROC_FLAGS $DB_FILE >/dev/null 2&>1 << END
        $XSLTPROC $DB_FILE >/dev/null 2&>1 << END
<?xml version="1.0" encoding='ISO-8859-1'?>
<!DOCTYPE book PUBLIC "-//OASIS//DTD DocBook XML V4.1.2//EN" "http://www.oasis-open.org/docbook/xml/4.1.2/docbookx.dtd">
<book id="test">
</book>
END
echo  yo $?
        if test "$?" = 0; then
                XSLTPROC_WORKS=yes
        fi
        AC_MSG_RESULT($XSLTPROC_WORKS)
	fi
	AM_CONDITIONAL(have_xsltproc, test "$XSLTPROC_WORKS" = "yes")

	AC_SUBST(XML_CATALOG)
	AC_SUBST(XSLTPROC_FLAGS)
	AC_SUBST(DOCBOOK_ROOT)
	AC_SUBST(CAT_ENTRY_START)
	AC_SUBST(CAT_ENTRY_END)
])
################################################################################
# F90 compiler
################################################################################
AC_DEFUN([AC_SR_FORTRAN],
[
AC_LANG(Fortran)
AC_PROG_FC(ifort fort xlf90 pgf90 epcf90 pathf90 ifc efc f90 xlf95 lf95 g95 f95)

])
dnl @synopsis AC_SR_BLAS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the BLAS
dnl linear-algebra interface (see http://www.netlib.org/blas/). On
dnl success, it sets the BLAS_LIBS output variable to hold the
dnl requisite library linkages.
dnl
dnl To link with BLAS, you should link with:
dnl
dnl 	$BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order. FLIBS is the output variable of the
dnl AC_FC_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS), and
dnl is sometimes necessary in order to link with FC libraries. Users
dnl will also need to use AC_FC_DUMMY_MAIN (see the autoconf manual),
dnl for the same reason.
dnl
dnl Many libraries are searched for, from ATLAS to CXML to ESSL. The
dnl user may also use --with-blas=<lib> in order to use some specific
dnl BLAS library <lib>. In order to link successfully, however, be
dnl aware that you will probably need to use the same Fortran compiler
dnl (which can be set via the FC env. var.) as was used to compile the
dnl BLAS library.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a BLAS
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands to
dnl run it if it is not found. If ACTION-IF-FOUND is not specified, the
dnl default action will define HAVE_BLAS.
dnl
dnl This macro requires autoconf 2.50 or later.
dnl
dnl @category InstalledPackages
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl @version 2001-12-13
dnl @license GPLWithACException

AC_DEFUN([AC_SR_BLAS], [
AC_PREREQ(2.50)
AC_REQUIRE([AC_FC_LIBRARY_LDFLAGS])
acx_blas_ok=no

AC_ARG_WITH(blas,
	[AC_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>])])
case $with_blas in
	yes | "") ;;
	no) acx_blas_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) BLAS_LIBS="$with_blas" ;;
	*) BLAS_LIBS="-l$with_blas" ;;
esac

# Get fortran linker names of BLAS functions to check for.
AC_FC_FUNC(sgemm)
AC_FC_FUNC(dgemm)
echo zzzzzzzz $sgemm

acx_blas_save_LIBS="$LIBS"
LIBS="$LIBS $FLIBS"

# First, check BLAS_LIBS environment variable
if test $acx_blas_ok = no; then
if test "x$BLAS_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
	AC_MSG_CHECKING([for $sgemm in $BLAS_LIBS])
	AC_COMPILE_IFELSE(program conftest
	$sgemm
	end, [acx_blas_ok=yes], [BLAS_LIBS=""])
dnl	AC_TRY_LINK_FUNC($sgemm, [acx_blas_ok=yes], [BLAS_LIBS=""])
	AC_MSG_RESULT($acx_blas_ok)
	LIBS="$save_LIBS"
fi
fi

# BLAS linked to by default?  (happens on some supercomputers)
if test $acx_blas_ok = no; then
	save_LIBS="$LIBS"; LIBS="$LIBS"
	AC_CHECK_FUNC($sgemm, [acx_blas_ok=yes])
	LIBS="$save_LIBS"
fi

# BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(atlas, ATL_xerbla,
		[AC_CHECK_LIB(f77blas, $sgemm,
		[AC_CHECK_LIB(cblas, cblas_dgemm,
			[acx_blas_ok=yes
			 BLAS_LIBS="-lcblas -lf77blas -latlas"],
			[], [-lf77blas -latlas])],
			[], [-latlas])])
fi

# BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(dgemm, $dgemm,
		[AC_CHECK_LIB(sgemm, $sgemm,
			[acx_blas_ok=yes; BLAS_LIBS="-lsgemm -ldgemm -lblas"],
			[], [-lblas])],
			[], [-lblas])])
fi

# BLAS in Alpha CXML library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(cxml, $sgemm, [acx_blas_ok=yes;BLAS_LIBS="-lcxml"])
fi

# BLAS in Alpha DXML library? (now called CXML, see above)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(dxml, $sgemm, [acx_blas_ok=yes;BLAS_LIBS="-ldxml"])
fi

# BLAS in Sun Performance library?
if test $acx_blas_ok = no; then
	if test "x$GCC" != xyes; then # only works with Sun CC
		AC_CHECK_LIB(sunmath, acosp,
			[AC_CHECK_LIB(sunperf, $sgemm,
        			[BLAS_LIBS="-xlic_lib=sunperf -lsunmath"
                                 acx_blas_ok=yes],[],[-lsunmath])])
	fi
fi

# BLAS in SCSL library?  (SGI/Cray Scientific Library)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(scs, $sgemm, [acx_blas_ok=yes; BLAS_LIBS="-lscs"])
fi

# BLAS in SGIMATH library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(complib.sgimath, $sgemm,
		     [acx_blas_ok=yes; BLAS_LIBS="-lcomplib.sgimath"])
fi

# BLAS in IBM ESSL library? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(essl, $sgemm,
			[acx_blas_ok=yes; BLAS_LIBS="-lessl -lblas"],
			[], [-lblas $FLIBS])])
fi

# Generic BLAS library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm, [acx_blas_ok=yes; BLAS_LIBS="-lblas"])
fi

AC_SUBST(BLAS_LIBS)

LIBS="$acx_blas_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_blas_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_BLAS,1,[Define if you have a BLAS library.]),[$1])
        :
else
        acx_blas_ok=no
        $2
fi
])dnl ACX_BLAS
