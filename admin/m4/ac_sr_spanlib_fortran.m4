################################################################################
### FORTRAN/BLAS/LAPACK
### Raynaud 2006-2010
################################################################################
AC_DEFUN([AC_SR_SPANLIB_FORTRAN],[

    # Fortran lib
    AC_CONFIG_SRCDIR([src/spanlib.f90])
    AC_SR_FORTRAN()

    # Tools
    AC_CHECK_TOOL(SED, sed, sed)
    AS_IF([test "AS_VAR_GET(SED)" = "no"],
        [AC_SR_ERROR([You need sed to build the library])])
    AC_SR_BLASLAPACK()
    AC_PROG_RANLIB()
    AC_CHECK_TOOL(AR, ar, ar)
    AS_IF([test "AS_VAR_GET(AR)" = "no"],
        [AC_SR_ERROR([You need ar to build the library])])

#     # Setup the sources for lapack95 usage
#     AC_MSG_NOTICE(Generating sources for lapack95 usage...)
#     sed -r -e 's/f95_lapack/AS_VAR_GET(LAPACK95_MOD)/' -e 's/\bla_(\S+)/la_\1 => (AS_VAR_GET(LAPACK95_PRE)\1/g' src/template.spanlib_lapack95.f90 > src/spanlib_lapack95.f90
 
    # Install fortran library?
    AC_MSG_CHECKING([whether the pure fortran library should be compilated and installed])
    AC_ARG_ENABLE(
        fortran,
        AS_HELP_STRING([--disable-fortran],[Turn off compilation and installation of the pure fortran library.]),
        [ENABLE_FORTRAN=no],
        [ENABLE_FORTRAN=yes]
    )
    #dnl AC_SUBST([ENABLE_FORTRAN],AS_VAR_GET(enable_fortran))
    AM_CONDITIONAL([ENABLE_FORTRAN],[test "AS_VAR_GET(ENABLE_FORTRAN)" = "yes"])
    AC_MSG_RESULT(AS_VAR_GET(ENABLE_FORTRAN))
        
 
   

])
################################################################################
################################################################################
################################################################################
