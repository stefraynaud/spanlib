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
