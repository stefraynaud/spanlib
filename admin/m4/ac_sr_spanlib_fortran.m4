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
