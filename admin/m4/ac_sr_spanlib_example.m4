################################################################################
### BLAS and LAPACK
### Raynaud 2006
################################################################################
AC_DEFUN([AC_SR_SPANLIB_EXAMPLE],[

	# We need netcdf
	AC_SR_NETCDF()
	AM_CONDITIONAL([WITH_EXAMPLE],[test "AS_VAR_GET(HAS_F90NETCDF)" = "yes"])
	AS_IF([test "AS_VAR_GET(HAS_F90NETCDF)" = "yes"],[
		AS_VAR_SET(ac_cv_text_example,["Finally, to run the f90 example, go to 'example', and type 'make'"])
	],[
		AC_SR_WARNING([Without f90 netcdf support, you wont be able to run the example])
	])

	# A commandline downloader may be useful to get the data
	AC_CHECK_PROG(WGET,wget,wget)
	AC_CHECK_PROG(LINKS,links,links)
	AS_VAR_SET(ac_cv_dl,
		[`test "AS_VAR_GET(WGET)" != "no" -o "AS_VAR_GET(LINKS)" != "no"`])
	AM_CONDITIONAL([HAS_DOWNLOADER],[AS_VAR_GET(ac_cv_dl)])
	AS_IF([AS_VAR_GET(ac_cv_dl)],,
		[AC_SR_WARNING([No commandline downloader found:
you will have to download yourself the input data file to run the example])])

	# A viewver may be useful to visualise the output netcdf file
	AC_CHECK_PROG(NCVIEW,ncview,ncview)
	AS_VAR_SET(ac_cv_vw,
		[`test "AS_VAR_GET(PYFORT)" != "no" -o "AS_VAR_GET(NCVIEW)" != "no"`])
	AM_CONDITIONAL([HAS_NCVIEWVER],[AS_VAR_GET(ac_cv_vw)])
	AS_IF([AS_VAR_GET(ac_cv_vw)],,
		[AC_SR_WARNING([No netcdf viewer available:
you will have to visualise yhe output netcdf file by own])])

])
################################################################################
################################################################################
################################################################################

