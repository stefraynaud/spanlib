################################################################################
# BLAS and LAPACK
# Raynaud 2006
################################################################################

AC_DEFUN([AC_SR_BLASLAPACK],
[dnl

#################################
# Initialisations
#################################
AS_VAR_SET(save_libs,$LIBS)
AS_VAR_SET(save_fcflags,$FCFLAGS)

LIBS=$FLIBS


#################################
# BLAS
#################################

# Library name or LIBS
AC_ARG_VAR(BLAS,[Library name or LIBS flag(s)])
AC_ARG_WITH(blas, dnl
 	AS_HELP_STRING(--with-blas=LIB, dnl
 		[Library name or LIBS flag(s)]),
 		[case AS_VAR_GET(with_blas) in
			no)AC_SR_ERROR([You cant disable blas]);;
			yes)AS_VAR_SET(BLAS,-lblas);;
			*)AS_VAR_SET(BLAS,-l$with_blas);;
		esac]
)
AS_VAR_SET_IF([BLAS],,[AS_VAR_SET(BLAS,-lblas)])
case AS_VAR_GET(BLAS) in
	-l* | */* | *.a | *.so | *.so.* | *.o):;;
	*)AS_VAR_SET(BLAS,"-l$BLAS");;
esac
AS_VAR_SET(LIBS,"$BLAS $LIBS")

# Directory where to find the library
AC_ARG_VAR(BLAS_LIBDIR,[Location of the BLAS library (compile-time)])
AC_ARG_WITH(blas-libdir, dnl
	AS_HELP_STRING(--with-blas-libdir=DIR, dnl
		[Location of the BLAS library (compile-time)]), dnl
		AS_IF([test $with_blas_libdir != yes -a $with_blas_libdir != no],
				AS_VAR_SET(BLAS_LIBDIR,$with_blas_libdir))
)
AS_VAR_SET_IF([BLAS_LIBDIR],[
	case AS_VAR_GET(BLAS_LIBDIR) in
		-L*):;;
		*)AS_VAR_SET(BLAS_LIBDIR,"-L$BLAS_LIBDIR");;
	esac
	AS_VAR_SET(FCFLAGS,"$BLAS_LIBDIR $FCFLAGS")
])

# Try sgemm with blas
AC_CACHE_CHECK([for sgemm of the blas library],ac_cv_blasok,
[AC_TRY_LINK_FUNC([sgemm],
                 [AS_VAR_SET(ac_cv_blasok,yes)],
                 [AS_VAR_SET(ac_cv_blasok,no)])
])

# Error
AS_IF([test AS_VAR_GET(ac_cv_blasok) = yes],,
		[AC_SR_ERROR([Impossible to find blas library. Try with --with-blas and --with-blas-libdir])])


#################################
# LAPACK
#################################

# Library name or LIBS
AC_ARG_VAR(LAPACK,Library name or LIBS flag(s))
AC_ARG_WITH(lapack, dnl
 	AS_HELP_STRING(--with-lapack=LIB, dnl
 		[Library name or LIBS flag(s)]),
		[case AS_VAR_GET(with_lapack) in
			no)AC_SR_ERROR([You cant disable lapack]);;
			yes)AS_VAR_SET(LAPACK,-llapack);;
			*)AS_VAR_SET(LAPACK,-l$with_lapack);;
		esac]
)
AS_VAR_SET_IF([LAPACK],,[AS_VAR_SET(LAPACK,-llapack)])
case AS_VAR_GET(LAPACK) in
	-l* | */* | *.a | *.so | *.so.* | *.o):;;
	*)AS_VAR_SET(LAPACK,"-l$LAPACK");;
esac
AS_VAR_SET(LIBS,"$LAPACK $LIBS")

# Library dir name or FCFLAGS
AC_ARG_VAR(LAPACK_LIBDIR,Location of the LAPACK library (compile-time))
AC_ARG_WITH(lapack-libdir, dnl
	AS_HELP_STRING(--with-lapack-libdir=DIR, dnl
		[Location of the LAPACK library (compile-time)]), dnl
		AS_IF([test $with_lapack_libdir != yes -a $with_lapack_libdir != no],
				AS_VAR_SET(LAPACK_LIBDIR,$with_lapack_libdir))
)
AS_VAR_SET_IF([LAPACK_LIBDIR],[
	case AS_VAR_GET(LAPACK_LIBDIR) in
		-L*):;;
		*)AS_VAR_SET(LAPACK_LIBDIR,"-L$LAPACK_LIBDIR");;
	esac
	AS_VAR_SET(FCFLAGS,"$LAPACK_LIBDIR $FCFLAGS")
])

# Try ssyev with lapack
AC_CACHE_CHECK([for ssyev of the lapack library],ac_cv_lapackok,
[AC_TRY_LINK_FUNC([ssyev],
                 [AS_VAR_SET(ac_cv_lapackok,yes)],
                 [AS_VAR_SET(ac_cv_lapackok,no)])
])

# Error
AS_IF([test AS_VAR_GET(ac_cv_lapackok) = yes],,
		[AC_SR_ERROR([Impossible to find lapack library. Try with --with-lapack and --witch-lapack-libdir])])



#################################
# Ending
#################################

AS_VAR_SET(LIBS,$save_libs)
AS_VAR_SET(FCFLAGS,$save_fcflags)
AC_SUBST(BLAS)
AC_SUBST(BLAS_LIBDIR)
AC_SUBST(LAPACK)
AC_SUBST(LAPACK_LIBDIR)

])dnl AC_SR_LAPACK
################################################################################
################################################################################
# F90 compiler
################################################################################
AC_DEFUN([AC_SR_FORTRAN],
[
## OPTIMIZATION LEVELS
##
AC_ARG_ENABLE(
    optimization,
    [  --enable-optimization=level - Control optimization level.
                             The following levels are supported.
       debug     - debugging compiler options will be selected.
       normal    - soft optimization (default).
       aggressive - aggressive optimization (YOU HAVE TO VERIFY YOUR RESULTS!).],
    ,
    enable_optimization=normal
             )
#
## PROFILING (if it is available)
##
AC_ARG_ENABLE(
    profiling,
    [  --enable-profiling - Turn on profiling compiler options.],
    enable_prof=yes,
    enable_prof=no
             )
# Main program
AC_LANG(Fortran)
AC_PROG_FC(ifort fort xlf90 pgf90 epcf90 pathf90 ifc efc f90 xlf95 lf95 g95 f95 sxf90)
if test ! -n "$FC" ; then
	AC_SR_ERROR([No Fortran 90 compiler available on this machine.
               Please use FC to specify it or
               update your environnement variable PATH or
               install a Fortran 90 compiler.])
fi

# LD FLAGS
  case "$FC" in
##    GENERIC FORTRAN COMPILER (SGI-IRIX, HP-TRUE64, NEC-SX )
      f90)
        case "$host" in
        *-sgi-irix*)
          case "$enable_optimization" in
	    debug)
              AC_MSG_NOTICE([**** DEBUGGING OPTIONS are SELECTED *****])
              FCFLAGS="-g -O0 -C -fullwarn -DEBUG:trap_uninitialized=ON:subscript_check=ON"
              LDFLAGS="-g"
            ;;
            aggressive)
              AC_MSG_NOTICE([**** AGGRESSIVE COMPILER OPTIONS are SELECTED *****])
    	      FCFLAGS="-g3 -O3 -ipa -listing"
              LDFLAGS="-g3"
            ;;
            normal|*)
              AC_MSG_NOTICE([**** NORMAL MODE *****])
	      FCFLAGS="-g3 -O2 -listing"
              LDFLAGS="-g3"
	    ;;
          esac
	  if test "$enable_prof" = "yes" ; then
            AC_SR_WARNING([!!! NO PROFILING COMPILER OPTIONS ON IRIX SYSTEM !!!])
            AC_SR_WARNING([!!!        PLEASE READ SPEEDSHOP MANUAL          !!!])
          fi
        ;;
        alpha*-dec-osf*)
          case "$enable_optimization" in
	    debug)
              AC_MSG_NOTICE([**** DEBUGGING OPTIONS are SELECTED *****])
              FCFLAGS="-V -ladebug -g -O0 -C -check overflow -check underflow -warn nouninitialized -warn argument_checking"
              LDFLAGS="-ladebug -g"
            ;;
            aggressive)
              AC_MSG_NOTICE([**** AGGRESSIVE COMPILER OPTIONS are SELECTED *****])
    	      FCFLAGS="-V -g3 -fast -math_library fast"
              LDFLAGS="-g3 -fast -math_library fast"
            ;;
            normal|*)
              AC_MSG_NOTICE([**** NORMAL MODE *****])
	      FCFLAGS="-V -g3 -O"
              LDFLAGS=""
	    ;;
          esac
	  if test "$enable_prof" = "yes" ; then
            AC_MSG_NOTICE([**** PROFILING is SELECTED (gprof) *****])
            FCFLAGS="-pg $FCFLAGS"
            LDFLAGS="-pg $LDFLAGS"
          fi
        ;;
        *nec*superux*)
          case "$enable_optimization" in
	    debug)
              AC_MSG_NOTICE([**** DEBUGGING OPTIONS are SELECTED *****])
              FCFLAGS='-C debug -eR -eP -R1 -R5 -Wf"-L nostdout" -Wf"-L source mrgmsg" -Wf"-L summary" -Wf"-init stack=nan" -Wf"-init heap=nan" -Wl"-f nan" Wf"-msg d" -Wf"-msg o"'
              LDFLAGS="-C debug"
            ;;
            aggressive)
              AC_MSG_NOTICE([**** AGGRESSIVE COMPILER OPTIONS are SELECTED *****])
    	      FCFLAGS='-C hopt -R1 -R5 -Wf"-L nostdout"  -Wf"-L summary" -Wf"-pvctl fullmsg" -Wf"-O infomsg"'
              LDFLAGS="-C hopt"
            ;;
            normal|*)
              AC_MSG_NOTICE([**** NORMAL MODE *****])
	      FCFLAGS='-R1 -R5 -Wf"-L nostdout"  -Wf"-L summary" -Wf"-pvctl fullmsg" -Wf"-O infomsg"'
              LDFLAGS=""
	    ;;
          esac
	  if test "$enable_prof" = "yes" ; then
            AC_MSG_NOTICE([**** PROFILING is SELECTED (gprof) *****])
            FCFLAGS="-ftrace $FCFLAGS"
            LDFLAGS="-ftrace $LDFLAGS"
          fi
        ;;
        *)
          AC_MSG_WARN([!!! HOST and/or SYSTEM is UNKNOWN : $host !!!])
          exit
        ;;
        esac
        ;;
##    INTEL FORTRAN COMPILER on LINUX OPERATING SYSTEM
      ifort|efc|ifc)
	case "$enable_optimization" in
	  debug)
            AC_MSG_NOTICE([**** DEBUGGING OPTIONS are SELECTED *****])
            FCFLAGS="-g -O0 -no_cpprt -check all -traceback -auto -warn all -warn unused -debug variable_locations -inline_debug_info"
	    LDFLAGS="-g -O0 -no_cpprt -check all -traceback -auto -inline_debug_info"
## if idb bugs use          FCFLAGS="-g -O0 "
## if idb bugs use 	    LDFLAGS="-g -O0 "
          ;;
          aggressive)
            AC_MSG_NOTICE([**** AGGRESSIVE COMPILER OPTIONS are SELECTED *****])
    	    FCFLAGS="-fast"
            LDFLAGS="-fast"
          ;;
          normal|*)
            AC_MSG_NOTICE([**** NORMAL MODE *****])
	    FCFLAGS="-g -O3 -132 -check bounds"
	  ;;
	esac
	if test "$enable_prof" = "yes" ; then
          AC_MSG_NOTICE([**** PROFILING is SELECTED (gprof) *****])
          FCFLAGS="-pg $FCFLAGS"
          LDFLAGS="-pg $LDFLAGS"
        fi
        ;;
##    IBM FORTRAN COMPILER on AIX OPERATING SYSTEM
      xlf90|xlf95)
	case "$enable_optimization" in
	  debug)
            FCFLAGS="-qsuffix=f=f90 -qfree=f90 -g -qnooptimize -C -qinitauto=7FBFFFFF -qflttrap=overflow:underflow:zerodivide:invalid:enable -qfloat=nans -qsigtrap -qextchk"
          ;;
          aggressive)
            FCFLAGS="-qsuffix=f=f90 -qfree=f90 -O3 -qstrict"
          ;;
          normal|*)
            FCFLAGS="-qsuffix=f=f90 -qfree=f90 -O5 -qipa=level=2 -qessl -qhot=vector -qunroll"
            LDFLAGS="-qessl"
	  ;;
	esac
	if test "$enable_prof" = "yes" ; then
          AC_MSG_NOTICE([**** PROFILING is SELECTED (gprof) *****])
          FCFLAGS="-pg $FCFLAGS"
          LDFLAGS="-pg $LDFLAGS"
        fi
        ;;
##    PORTLAND GROUP FORTRAN COMPILER
      pgf90)
	FCFLAGS="-g"
        ;;
##    GENERIC Fortran 95 compiler (not tested)
      f95)
	FCFLAGS="-g"
        ;;
##    HP_COMPAQ ALPHASERVER FORTRAN COMPILER (LINUX OPERATING SYSTEM)
      fort)
	FCFLAGS="-g"
        ;;
##    Lahey-Fujitsu compiler
      lf95)
	FCFLAGS="-g"
        ;;
##    PATHSCALE FORTRAN COMPILER (AMD-OPTERON) (Not Tested)
      pathf90)
	FCFLAGS="-g"
        ;;
##    GNU FORTRAN 90/95 COMPILER (Tested on Intel-PC and Mac OS X)
      g95)
	case "$enable_optimization" in
	  debug)
            AC_MSG_NOTICE([**** DEBUGGING OPTIONS are SELECTED *****])
            FCFLAGS="-g -O0 -fno-second-underscore -Wall -Wunset-vars -Wunused-vars -fbounds-check "
	    LDFLAGS="-g -O0-fno-second-underscore"
          ;;
          aggressive)
            AC_MSG_NOTICE([**** AGGRESSIVE COMPILER OPTIONS are SELECTED *****])
    	    FCFLAGS="-g -O3 -fno-second-underscore"
            LDFLAGS="-g -O3 -fno-second-underscore"
          ;;
          normal|*)
            AC_MSG_NOTICE([**** NORMAL MODE *****])
	    FCFLAGS="-g -O -fno-second-underscore"
            LDFLAGS="-g -O -fno-second-underscore"
	  ;;
	esac
	if test "$enable_prof" = "yes" ; then
          AC_MSG_NOTICE([**** PROFILING is SELECTED (gprof) *****])
          FCFLAGS="-pg $FCFLAGS"
          LDFLAGS="-pg $LDFLAGS"
        fi
        ;;
  esac
	AC_FC_SRCEXT(f90)
	AC_FC_FREEFORM()
])
################################################################################
# F90 netcdf
################################################################################
AC_DEFUN([AC_SR_NETCDF],
[

# Save variables
AS_VAR_SET(save_libs,$LIBS)
AS_VAR_SET(save_fcflags,$FCFLAGS)


# Module location
AC_ARG_VAR(NETCDF_INC,Location of netCDF module (compile-time))



AC_ARG_WITH(netcdf-inc, dnl
	AS_HELP_STRING(--with-netcdf-inc=DIR, dnl
		[Location of netCDF module (compile-time)]), dnl
	[
		AS_IF([test $with_netcdf_inc != yes -a $with_netcdf_inc != no],
			AS_VAR_SET(NETCDF_INC,$with_netcdf_inc))
	]
)
AS_VAR_SET_IF(NETCDF_INC,[
	case AS_VAR_GET(NETCDF_INC) in
		-I*);;
		*)AS_VAR_SET(NETCDF_INC,"-I$NETCDF_INC");;
	esac
	AS_VAR_SET(FCFLAGS,"$NETCDF_INC $FCFLAGS")
])

# Library
AC_ARG_VAR(NETCDF_LIB,Location of netCDF library (compile-time))
AC_ARG_WITH(netcdf-lib, dnl
	AS_HELP_STRING(--with-netcdf-lib=DIR, dnl
		[Location of netCDF library (compile-time)]), dnl
	[
		AS_IF([test $with_netcdf_lib != yes -a $with_netcdf_lib != no],
			AS_VAR_SET(NETCDF_LIB,$with_netcdf_lib))
	]
)
AS_VAR_SET_IF(NETCDF_LIB,[
	case AS_VAR_GET(NETCDF_LIB) in
		-L*);;
		*)AS_VAR_SET(NETCDF_LIB,"-L$NETCDF_LIB");;
	esac
	AS_VAR_SET(FCFLAGS,"$NETCDF_LIB $FCFLAGS")
])
AS_VAR_SET(FCFLAGS,"$FCFLAGS -lnetcdf")

# Check
AC_MSG_CHECKING([for f90 netcdf support])
AC_COMPILE_IFELSE([subroutine conftest_routine
	use netcdf
	integer :: n
	n = nf90_close(1)
end subroutine conftest_routine
],AS_VAR_SET(HAS_F90NETCDF,yes),AS_VAR_SET(HAS_F90NETCDF,no))
AC_MSG_RESULT([AS_VAR_GET(HAS_F90NETCDF)])

# # End
AS_VAR_SET(LIBS,$save_libs)
AS_VAR_SET(FCFLAGS,$save_fcflags)
AC_SUBST(NETCDF_LIB)
AC_SUBST(NETCDF_INC)
AC_SUBST(HAS_F90NETCDF)

])
################################################################################
################################################################################
# Spanlib, Pyfort with BLAS/LAPACK
# Raynaud 2006
################################################################################
# AC_SR_STRIPFLAGS([FLAGS],[OUTPUT_VARIABLE])
# Remove -L and -l
AC_DEFUN([AC_SR_PYFORT_STRIPFLAGS],
[
	AS_VAR_SET($2,`echo $1 | sed -r 's/(^|\s)-(L|l)/\1/g'`)
])
# Setup library dir and name variables PYFORT_DIRS and PYFORT_LIBS
AC_DEFUN([AC_SR_PYFORT],
[
	# Full path to pyfort
	AC_PATH_PROG(PYFORT,pyfort,pyfort)
	AM_CONDITIONAL([HAS_PYFORT],[test AS_VAR_GET(PYFORT) != "no"])

	# Blas and Lapack libraries and directories
	AC_SR_PYFORT_STRIPFLAGS("AS_VAR_GET(BLAS) AS_VAR_GET(LAPACK)",PYFORT_LIBS)
	AC_SUBST(PYFORT_LIBS)
	AC_SR_PYFORT_STRIPFLAGS("AS_VAR_GET(BLAS_LIBDIR) AS_VAR_GET(LAPACK_LIBDIR)",PYFORT_DIRS)
	AC_SUBST(PYFORT_DIRS)

	# Set where to install the python module
	AC_ARG_VAR(PYPACKDIR,[Directory where to install the spanlib python module])
	AC_ARG_WITH(pypackdir, dnl
		AS_HELP_STRING(--with-pypackdir=DIR, dnl
				[Directory where to install the spanlib python module]),
				[case AS_VAR_GET(with_pypackdir) in
					no|yes);;
					*)AS_VAR_SET(PYPACKDIR,AS_VAR_GET(with_pypackdir));;
				esac]
	)
	AS_VAR_SET_IF(PYPACKDIR,
		[AS_VAR_SET(PYFORT_BUILD,[-b])],
		[
#			AS_VAR_SET(PYFORT_BUILD,[-i])
			AS_IF([test AS_VAR_GET(PYFORT) != "no"],
				[AS_VAR_SET(PYPACKDIR,
					[`AS_VAR_GET(PYTHON) -c ["import sys; print sys.prefix+\"/lib/python\"+sys.version[:3]+\"/site-packages\""]`])],
				[AS_VAR_SET(PYPACKDIR,"")]
			)
		]
	)
	AC_SUBST(PYFORT_BUILD)

])
################################################################################
################################################################################
### BLAS and LAPACK
### Raynaud 2006
################################################################################
AC_DEFUN([AC_SR_SPANLIB_DOC],[

	# docbook (xslt) processor support
	AS_DOCBOOK()

	# We need perl and xslt processor
	AM_CONDITIONAL([HAS_DOC_SUPPORT],
		[test AS_VAR_GET(PERL) != "no" -a AS_VAR_GET(XSLTPROC) != "no"])

])
################################################################################
################################################################################
################################################################################

################################################################################
### BLAS and LAPACK
### Raynaud 2006
################################################################################
AC_DEFUN([AC_SR_SPANLIB_EXAMPLE],[

	# We need netcdf
	AC_SR_NETCDF()
	AM_CONDITIONAL([WITH_EXAMPLE],[test AS_VAR_GET(HAS_F90NETCDF) = "yes"])
	AS_IF([test AS_VAR_GET(HAS_F90NETCDF) = "yes"],[
		AS_VAR_SET(ac_cv_text_example,["Finally, to run the f90 example, go to 'example', and type 'make'"])
	],[
		AC_SR_WARNING([Without f90 netcdf support, you wont be able to run the example])
	])

	# A commandline downloader may be useful to get the data
	AC_CHECK_PROG(WGET,wget,wget)
	AC_CHECK_PROG(LINKS,links,links)
	AS_VAR_SET(ac_cv_dl,
		[`test AS_VAR_GET(WGET) != "no" -o AS_VAR_GET(LINKS) != "no"`])
	AM_CONDITIONAL([HAS_DOWNLOADER],[AS_VAR_GET(ac_cv_dl)])
	AS_IF([AS_VAR_GET(ac_cv_dl)],,
		[AC_SR_WARNING([No commandline downloader found:
you will have to download yourself the input data file to run the example])])

	# A viewver may be useful to visualise the output netcdf file
AC_SR_WARNING([YOOOOO1])
	AC_CHECK_PROG(NCVIEW,ncview,ncview)
AC_SR_WARNING([YOOOOO2])
	AS_VAR_SET(ac_cv_vw,
		[`test AS_VAR_GET(CDAT) != "no" -o AS_VAR_GET(NCVIEW) != "no"`])
AC_SR_WARNING([YOOOOO3])
	AM_CONDITIONAL([HAS_NCVIEWVER],[AS_VAR_GET(ac_cv_vw)])
AC_SR_WARNING([YOOOOO4])
	AS_IF([AS_VAR_GET(ac_cv_vw)],,
		[AC_SR_WARNING([No netcdf viewer available:
you will have to visualise yhe output netcdf file by own])])
AC_SR_WARNING([YOOOOO5])

])
################################################################################
################################################################################
################################################################################

################################################################################
### BLAS and LAPACK
### Raynaud 2006
################################################################################
AC_DEFUN([AC_SR_SPANLIB_FORTRAN],[

	# Fortran lib
	AC_CONFIG_SRCDIR([src/spanlib.f90])
	AC_SR_FORTRAN()

	# Tools
	AC_CHECK_TOOL(SED, sed, sed)
	AS_IF([test AS_VAR_GET(SED) = "no"],
		[AC_SR_ERROR([You need sed to build the library])])
	AC_SR_BLASLAPACK()
	AC_PROG_RANLIB()
	AC_CHECK_TOOL(AR, ar, ar)
	AS_IF([test AS_VAR_GET(AR) = "no"],
		[AC_SR_ERROR([You need ar to build the library])])

])
################################################################################
################################################################################
################################################################################
################################################################################
### BLAS and LAPACK
### Raynaud 2006
################################################################################
AC_DEFUN([AC_SR_SPANLIB_PYTHON],[

	AC_CHECK_PROG(PERL,perl,perl)
	AC_PROG_CC()
	AX_WITH_PYTHON(2.4)
	AC_SR_PYFORT()

	AS_VAR_SET(ac_cv_py,
		[`test AS_VAR_GET(PYFORT) != "no" -a AS_VAR_GET(PERL) != "no"`])
	AM_CONDITIONAL([WITH_PYTHON],[AS_VAR_GET(ac_cv_py)])
	AS_IF([AS_VAR_GET(ac_cv_py)],,
		[AC_SR_WARNING([You wont be able to build the python module.
You need perl and pyfort (cdat from python).])])


])
################################################################################
################################################################################
################################################################################
################################################################################
### SpanLib
### Raynaud 2006
################################################################################
# Colorisation for fun!
################################################################################
AC_DEFUN([AC_SR_SET_YELLOWINK],
[
	echo -en "\\033\\1331;33m"
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
################################################################################
################################################################################

dnl AS_DOCBOOK([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]])
dnl checks if xsltproc can build docbook documentation
dnl (which is possible if the catalog is set up properly
dnl I also tried checking for a specific version and type of docbook
dnl but xsltproc seemed to happily run anyway, so we can t check for that
dnl and version
dnl this macro takes inspiration from
dnl http://www.movement.uklinux.net/docs/docbook-autotools/configure.html
AC_DEFUN([AS_DOCBOOK],
[
  XSLTPROC_FLAGS=--nonet
  DOCBOOK_ROOT=
  TYPE_LC=xml
  TYPE_UC=XML
  DOCBOOK_VERSION=4.1.2

  if test ! -f /etc/xml/catalog; then
    for i in /usr/share/sgml/docbook/stylesheet/xsl/nwalsh /usr/share/sgml/docbook/xsl-stylesheets/;
    do
      if test -d "$i"; then
        DOCBOOK_ROOT=$i
      fi
    done
  else
    XML_CATALOG=/etc/xml/catalog
    CAT_ENTRY_START='<!--'
    CAT_ENTRY_END='-->'
  fi

  dnl We need xsltproc to process the test
  AC_CHECK_PROG(XSLTPROC,xsltproc,xsltproc,)
  XSLTPROC_WORKS=no
  if test -n "$XSLTPROC"; then
    AC_MSG_CHECKING([whether xsltproc docbook processing works])

    if test -n "$XML_CATALOG"; then
      DB_FILE="http://docbook.sourceforge.net/release/xsl/current/xhtml/docbook.xsl"
    else
      DB_FILE="$DOCBOOK_ROOT/docbook.xsl"
    fi
    $XSLTPROC $XSLTPROC_FLAGS $DB_FILE >/dev/null 2>&1 << END
<?xml version="1.0" encoding='ISO-8859-1'?>
<!DOCTYPE book PUBLIC "-//OASIS//DTD DocBook $TYPE_UC V$DOCBOOK_VERSION//EN" "http://www.oasis-open.org/docbook/$TYPE_LC/$DOCBOOK_VERSION/docbookx.dtd">
<book id="test">
</book>
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

  AC_SUBST(XML_CATALOG)
  AC_SUBST(XSLTPROC_FLAGS)
  AC_SUBST(DOCBOOK_ROOT)
  AC_SUBST(CAT_ENTRY_START)
  AC_SUBST(CAT_ENTRY_END)
])dnl @synopsis AX_WITH_PYTHON([minimum-version], [value-if-not-found], [path])
dnl
dnl Locates an installed Python binary, placing the result in the
dnl precious variable $PYTHON. Accepts a present $PYTHON, then
dnl --with-python, and failing that searches for python in the given
dnl path (which defaults to the system path). If python is found,
dnl $PYTHON is set to the full path of the binary; if it is not found,
dnl $PYTHON is set to VALUE-IF-NOT-FOUND, which defaults to 'python'.
dnl
dnl Example:
dnl
dnl   AX_WITH_PYTHON(2.2, missing)
dnl
dnl @category InstalledPackages
dnl @author Dustin Mitchell <dustin@cs.uchicago.edu>
dnl @version 2005-01-22
dnl @license GPLWithACException

dnl Raynaud 2006

AC_DEFUN([AX_WITH_PYTHON],
[
  AC_ARG_VAR([PYTHON],[absolute path name of Python executable])

  dnl unless PYTHON was supplied to us (as a precious variable)
dnl  if test -z "$PYTHON"
  AS_VAR_SET_IF(PYTHON,,[
    AC_MSG_CHECKING(for --with-python)
    AC_ARG_WITH(python,
                AC_HELP_STRING([--with-python=PYTHON],
                               [absolute path name of Python executable]),
                [ if test "$withval" != "yes"
                  then
                    PYTHON="$withval"
                    AC_MSG_RESULT($withval)
                  else
                    AC_MSG_RESULT(no)
                  fi
                ],
                [ AC_MSG_RESULT(no)
                ])
  ])

  dnl if it's still not found, check the paths, or use the fallback
dnl  if test -z "$PYTHON"
  AS_VAR_SET_IF(PYTHON,,[
    AC_PATH_PROG([PYTHON], python, m4_ifval([$2],[$2],[python]), $3)
  ])

  dnl check version if required
  m4_ifvaln([$1], [
    dnl do this only if we didn't fall back
    if test "$PYTHON" != "m4_ifval([$2],[$2],[python])"
    then
      AC_MSG_CHECKING($PYTHON version >= $1)
      if test `$PYTHON -c ["import sys; print sys.version[:3] >= \"$1\" and \"OK\" or \"OLD\""]` = "OK"
      then
        AC_MSG_RESULT(ok)
      else
        AC_MSG_RESULT(no)
        PYTHON="$2"
      fi
    fi])
])
