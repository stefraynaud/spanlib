################################################################################
### BLAS and LAPACK
### Raynaud 2006
################################################################################
AC_DEFUN([AC_SR_SPANLIB_EXAMPLE],[

	#################################################################
	# We need netcdf
	#################################################################

	AC_SR_NETCDF()
	AS_IF([test "AS_VAR_GET(HAS_F90NETCDF)" = "yes"],[
		AS_VAR_SET(F90_EXAMPLE_TEXT,["To run the f90 example, type from here '"AS_VAR_GET([NORMAL])"make example"AS_VAR_GET(GREEN)"'."])
	],[
		AC_SR_WARNING([Without f90 netcdf support, you wont be able to run the f90 example])
	])

	#################################################################
	# Blas/Lapack (checked before)
	#################################################################
	AM_CONDITIONAL([WITH_EXAMPLE],
		[test "AS_VAR_GET(HAS_BLASLAPACK)" != "no" -a "AS_VAR_GET(HAS_F90NETCDF)" = "yes"])

	#################################################################
	# Python example
	#################################################################
	AS_VAR_SET_IF(WITH_PYTHON_EXAMPLE,[
		AS_VAR_SET(PYTHON_EXAMPLE_TEXT,["To run the python examples, type from here '"AS_VAR_GET([NORMAL])"cd example && make python1 && make python2"AS_VAR_GET(GREEN)"'."])
	])


	#################################################################
	# A commandline downloader may be useful to get the data
	#################################################################

	AC_CHECK_PROG(WGET,wget,wget,no)
	AS_IF([test "AS_VAR_GET(WGET)" != "no"],
		AS_VAR_SET(DOWNLOADER,AS_VAR_GET(WGET)))
	AS_VAR_SET_IF(DOWNLOADER,,[
			AC_CHECK_PROG(LINKS,links,links,no)
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

	# Ncview
	AC_CHECK_PROG(NCVIEW,ncview,ncview,no)
	AS_IF([test "AS_VAR_GET(NCVIEW)" != "no"],[
		AS_VAR_SET(NCVIEWER,AS_VAR_GET(NCVIEW))
		AS_VAR_SET(NCVIEWER_ARGS,['$(F90NCOUTPUT)'])
	])

	# VCDAT
	AS_VAR_SET_IF(NCVIEWER,,[
		AS_IF([test "AS_VAR_GET(VCDAT)" != "no"],[
			AS_VAR_SET(NCVIEWER,AS_VAR_GET(VCDAT))
			AS_VAR_SET(NCVIEWER_ARGS,[''])
		])])

	# Quick python viewer
	AS_VAR_SET_IF(NCVIEWER,,[
		AS_IF(AS_VAR_GET(WITH_PYTHON_EXAMPLE),[
			AS_VAR_SET(NCVIEWER,[AS_VAR_GET(PYTHON) ../scripts/quickplot.py])
			AS_VAR_SET(NCVIEWER_ARGS,['$(F90NCOUTPUT) $(VARIABLE)'])
		])])

	# So...
	AS_VAR_SET_IF(NCVIEWER,,
		AC_SR_WARNING([No netcdf viewer available:
you will have to visualise the output netcdf file by your own]))
	AC_SUBST(NCVIEWER)
	AC_SUBST(NCVIEWER_ARGS)
	AM_CONDITIONAL(HAS_NCVIEWER,AS_VAR_TEST_SET(NCVIEWER))

])
################################################################################
################################################################################
################################################################################

