################################################################################
### BLAS and LAPACK
### Raynaud 2006-2010
################################################################################
AC_DEFUN([AC_SR_SPANLIB_EXAMPLE],[

    #################################################################
    # We need netcdf
    #################################################################
    AS_IF([test "AS_VAR_GET(ENABLE_FORTRAN)" = "yes"],[
        AC_SR_NETCDF()
        AS_IF([test "AS_VAR_GET(HAS_F90NETCDF)" = "yes"],[
            AS_VAR_SET(F90_EXAMPLE_TEXT,['To run fortran examples, go to the "example" directory and type "make fortran"'])
        ],[
            AC_SR_WARNING([Without f90 netcdf support, you wont be able to run the f90 example])
        ])
        AS_VAR_SET(ENABLE_FORTRAN_EXAMPLES,"yes")
    ],[
        AS_VAR_SET(ENABLE_FORTRAN_EXAMPLES,"no")
    ])
    AM_CONDITIONAL([ENABLE_FORTRAN_EXAMPLES],[test AS_VAR_GET(ENABLE_FORTRAN_EXAMPLES) = "yes"])

    #################################################################
    # Blas/Lapack (checked before)
    #################################################################
    AM_CONDITIONAL([ENABLE_EXAMPLE],
        [test "AS_VAR_GET(HAS_BLASLAPACK)" != "no" -a "AS_VAR_GET(HAS_F90NETCDF)" = "yes"])

    #################################################################
    # Python example help
    #################################################################
    AS_IF([test "AS_VAR_GET(ENABLE_PYTHON)" = "yes" -a "AS_VAR_GET(HAS_BLASLAPACK)" != "no"],
        [AS_VAR_SET(ENABLE_PYTHON_EXAMPLES,"yes")],
        [AS_VAR_SET(ENABLE_PYTHON_EXAMPLES,"no")])
    AM_CONDITIONAL([ENABLE_PYTHON_EXAMPLES],[test AS_VAR_GET(ENABLE_PYTHON_EXAMPLES) = "yes"])


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
you will have to download yourself the input data file to run some examples]))
    AC_SUBST(DOWNLOADER)
    AM_CONDITIONAL(HAS_DOWNLOADER,AS_VAR_TEST_SET(DOWNLOADER))
    AM_CONDITIONAL(USE_LINKS,
        [test "AS_VAR_GET(DOWNLOADER)" = "AS_VAR_GET(LINKS)"])


])
################################################################################
################################################################################
################################################################################

