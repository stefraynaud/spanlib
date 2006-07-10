################################################################################
### BLAS and LAPACK
### Raynaud 2006
################################################################################
AC_DEFUN([AC_SR_SPANLIB_EXAMPLE],[

	#################################################################
	# We need netcdf
	#################################################################

	AC_SR_NETCDF()
	AS_VAR_SET(ac_cv_example,[`test "AS_VAR_GET(HAS_F90NETCDF)" = "yes"`])
	AM_CONDITIONAL(WITH_EXAMPLE,AS_VAR_GET(ac_cv_example))
	AS_IF(AS_VAR_GET(ac_cv_example),[
		AS_VAR_SET(ac_cv_text_example,["Finally, to run the f90 example, go to 'example', and type 'make'"])
	],[
		AC_SR_WARNING([Without f90 netcdf support, you wont be able to run the example])
	])

	#################################################################
	# Blas/Lapack (checked before)
	#################################################################
	AS_IF(AS_VAR_GET(ac_cv_example),
		AM_CONDITIONAL([WITH_EXAMPLE],
			[test "AS_VAR_GET(HAS_BLASLAPACK)" != "no"]),,)



	#################################################################
	# A commandline downloader may be useful to get the data
	#################################################################

	AC_CHECK_PROG(WGET,wget,wget)
	AS_IF([test "AS_VAR_GET(WGET)" != "no"],
		AS_VAR_SET(DOWNLOADER,AS_VAR_GET(WGET)))
	AS_VAR_SET_IF(DOWNLOADER,,[
			AC_CHECK_PROG(LINKS,links,links)
			AS_IF([test "AS_VAR_GET(LINKS)" != "no"],
				AS_VAR_SET(DOWNLOADER,AS_VAR_GET(LINKS))) ])
	AS_VAR_SET_IF(DOWNLOADER,,
		AC_SR_WARNING([No commandline downloader found:
you will have to download yourself the input data file to run the example]))
	AC_SUBST(DOWNLOADER)
	AM_CONDITIONAL(HAS_DOWNLOADER,AS_VAR_TEST_SET(DOWNLOADER))
	AM_CONDITIONAL(USE_LINKS,
		[test "AS_VAR_GET(DOWNLOADER)" = "AS_VAR_GET(LINKS)"])



	#################################################################
	# A viewver may be useful to visualise the output netcdf file
	#################################################################

	AS_IF([test "AS_VAR_GET(PYFORT)" != "no"],
		AS_VAR_SET(NCVIEWER,AS_VAR_GET(PYFORT)))
	AS_VAR_SET_IF(NCVIEWER,,[
			AC_CHECK_PROG(NCVIEW,ncview,ncview)
			AS_IF([test "AS_VAR_GET(NCVIEW)" != "no"],
				AS_VAR_SET(NCVIEWER,AS_VAR_GET(NCVIEW))) ])
	AS_VAR_SET_IF(NCVIEWER,,
		AC_SR_WARNING([No netcdf viewer available:
you will have to visualise yhe output netcdf file by own]))
	AC_SUBST(NCVIEWER)
	AM_CONDITIONAL(HAS_NCVIEWER,AS_VAR_TEST_SET(NCVIEWER))

])
################################################################################
################################################################################
################################################################################

