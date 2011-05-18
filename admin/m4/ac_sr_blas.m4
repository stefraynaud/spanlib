################################################################################
# BLAS and LAPACK
# SpanLib, Raynaud 2006-2011
################################################################################

AC_DEFUN([AC_SR_BLASLAPACK],
[

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
AC_MSG_CHECKING([for blas library flag])
AC_ARG_VAR(BLAS,[Library name or LIBS flag(s)])
AC_ARG_WITH(blas,
    AC_HELP_STRING([--with-blas=LIBNAME],
        [Library name or LIBS flag(s)]),
    [case AS_VAR_GET(with_blas) in
        no)AC_SR_ERROR([You cant disable blas]);;
        yes)AS_VAR_SET(BLAS,"-lblas");;
        *)AS_VAR_SET(BLAS,$with_blas);;
    esac]
)
AS_VAR_SET_IF([BLAS],,[AS_VAR_SET(BLAS,"-lblas")])
case AS_VAR_GET(BLAS) in
    -l* | */* | *.a | *.so | *.so.* | *.o):;;
    *)AS_VAR_SET(BLAS,"-l$BLAS");;
esac
AS_VAR_SET(LIBS,"$BLAS $LIBS")
AC_MSG_RESULT(AS_VAR_GET(BLAS))    

# Directory where to find the library
AC_MSG_CHECKING([for blas library location flag])
AC_ARG_VAR(BLAS_LIB,[Location of the BLAS library (compile-time)])
AC_ARG_WITH(blas-lib,
    AC_HELP_STRING([--with-blas-lib=DIR],
        [Location of the BLAS library (compile-time)]), 
    AS_IF([test "$with_blas_lib" != "yes" -a "$with_blas_lib" != "no"],
            AS_VAR_SET(BLAS_LIB,$with_blas_lib))
)
AS_VAR_SET_IF(BLAS_LIB,,[
    AS_IF(test AS_VAR_GET(prefix) != "None/include",
        AS_VAR_SET(BLAS_LIB,AS_VAR_GET(prefix)/lib))
])
AS_VAR_SET_IF([BLAS_LIB],[
    case AS_VAR_GET(BLAS_LIB) in
        -L*):;;
        *)AS_VAR_SET(BLAS_LIB,"-L$BLAS_LIB");;
    esac
    AS_VAR_SET(FCFLAGS,"$BLAS_LIB $FCFLAGS")
])
AC_MSG_RESULT(AS_VAR_GET(BLAS_LIB))

# Try sgemm with blas
AC_MSG_CHECKING([for dgemm of the blas library])
AC_TRY_LINK_FUNC(dgemm,
                 [AS_VAR_SET(ac_cv_blasok,yes)],
                 [AS_VAR_SET(ac_cv_blasok,no)])
AC_MSG_RESULT(AS_VAR_GET(ac_cv_blasok))

#################################
# LAPACK
#################################


# F77 library name or LIBS
AC_MSG_CHECKING([for lapack f77 library flag])
AC_ARG_VAR(LAPACK,F77 library name or LIBS flag(s))
AC_ARG_WITH(lapack,
    AC_HELP_STRING([--with-lapack=LIBNAME],
        [F77 library name or LIBS flag(s)]),
    [case AS_VAR_GET(with_lapack) in
        no)AC_SR_ERROR([You cant disable lapack]);;
        yes)AS_VAR_SET(LAPACK,"-llapack");;
        *)AS_VAR_SET(LAPACK,$with_lapack);;
    esac] 
)
AS_VAR_SET_IF([LAPACK],,[AS_VAR_SET(LAPACK,"-llapack")])
case AS_VAR_GET(LAPACK) in
    -l* | */* | *.a | *.so | *.so.* | *.o):;;
    *)AS_VAR_SET(LAPACK,"-l$LAPACK");;
esac
AS_VAR_SET(LIBS,"$LAPACK $LIBS")
AC_MSG_RESULT(AS_VAR_GET(LAPACK))

# Library dir name or FCFLAGS for Lapack
AC_MSG_CHECKING([for lapack library location flag])
AC_ARG_VAR(LAPACK_LIB,Location of the LAPACK library (compile-time))
AC_ARG_WITH(lapack-lib,
    AC_HELP_STRING([--with-lapack-lib=DIR],
        [Location of the LAPACK library (compile-time)]),
    AS_IF([test "$with_lapack_lib" != "yes" -a "$with_lapack_lib" != "no"],
        AS_VAR_SET(LAPACK_LIB,$with_lapack_lib))
)
AS_VAR_SET_IF([LAPACK_LIB],[
    case AS_VAR_GET(LAPACK_LIB) in
        -L*):;;
        *)AS_VAR_SET(LAPACK_LIB,"-L$LAPACK_LIB");;
    esac
    AS_VAR_SET(FCFLAGS,"$LAPACK_LIB $FCFLAGS")
])
AC_MSG_RESULT(AS_VAR_GET(LAPACK_LIB))


# Try sdyev with blas
AC_MSG_CHECKING([for dsyev of the lapack77 library])
AC_TRY_LINK_FUNC(dsyev,
                 [AS_VAR_SET(ac_cv_lapackok,yes)],
                 [AS_VAR_SET(ac_cv_lapackok,no)])
AC_MSG_RESULT(AS_VAR_GET(ac_cv_lapackok))


#################################
# Ending
#################################

# Warning
AS_IF([test "AS_VAR_GET(ac_cv_blasok)" != "yes" -o "AS_VAR_GET(ac_cv_lapackok)" != "yes"],
    [
        AS_VAR_SET(HAS_BLASLAPACK,no)
        AC_SR_WARNING([It seems you have no Blas/Lapack developpment support.
Try with switches --with-blas/lapack and/or --with-blas/lapack-lib.
Without it, you wont be able to run the example and use the python interface.])
    ],
    AS_VAR_SET(HAS_BLASLAPACK,yes)
)

# Variables
AS_VAR_SET(LIBS,$save_libs)
AS_VAR_SET(FCFLAGS,$save_fcflags)


])
################################################################################
