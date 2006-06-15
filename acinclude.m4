#####################################################################
# AC macros for spanlib
# Stephane Raynaud 2006
#####################################################################

################################################################################
# Colorisation for fun!
################################################################################
AC_DEFUN([AC_SR_HEADER],
[
        echo -en "\\033\\1331;32m"
        echo "################################################################################"
        echo "### $1"
        echo "################################################################################"
        echo -en "\\033\\1330;39m"
])
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
