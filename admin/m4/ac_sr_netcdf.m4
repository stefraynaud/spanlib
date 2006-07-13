################################################################################
# F90 netcdf
################################################################################
AC_DEFUN([AC_SR_NETCDF],
[

# Save variables
AS_VAR_SET(save_libs,AS_VAR_GET(LIBS))
AS_VAR_SET(save_fcflags,AS_VAR_GET(FCFLAGS))


# Module location
AC_ARG_VAR([NETCDF_INC],[Locationzz of netCDF module (compile-time)])
AC_ARG_WITH(netcdf-inc, dnl
	AS_HELP_STRING(--with-netcdf-inc=DIR, dnl
		[Locationxx of netCDF module (compile-time)]), dnl
	[
		AS_IF([test "$with_netcdf_inc" != "yes" -a "$with_netcdf_inc" != "no"],
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
AC_ARG_VAR([NETCDF_LIB],[Location of netCDF library (compile-time)])
AC_ARG_WITH(netcdf-lib, dnl
	AS_HELP_STRING(--with-netcdf-lib=DIR, dnl
		[Location of netCDF library (compile-time)]), dnl
	[
		AS_IF([test "$with_netcdf_lib" != "yes" -a "$with_netcdf_lib" != "no"],
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
AS_VAR_SET(LIBS,"$LIBS -lnetcdf")

# Check
AC_MSG_CHECKING([for f90 netcdf support])
###dnl AC_COMPILE_IFELSE([subroutine conftest_routine])
AC_LINK_IFELSE([program conftest_routine
	use netcdf
	integer :: n
	n = nf90_close(1)
end program conftest_routine],
	AS_VAR_SET(HAS_F90NETCDF,yes),AS_VAR_SET(HAS_F90NETCDF,no))
AC_MSG_RESULT([AS_VAR_GET(HAS_F90NETCDF)])

# # End
AS_VAR_SET(LIBS,$save_libs)
AS_VAR_SET(FCFLAGS,$save_fcflags)
AC_SUBST(NETCDF_LIB)
AC_SUBST(NETCDF_INC)
AC_SUBST(HAS_F90NETCDF)

])
################################################################################
