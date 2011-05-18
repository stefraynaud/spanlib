################################################################################
### SpanLib, Raynaud 2006
################################################################################
AC_DEFUN([AC_SR_SPANLIB_PYTHON],[

    # Tools
#     AC_CHECK_PROG(PERL,perl,perl,no)
    AC_PROG_CC()

    # Default conditionals
    AM_CONDITIONAL(HAS_NUMPY,false)
    AM_CONDITIONAL(LOCAL_PYTHON_INSTALL,false)

    # Install python package?
    AC_MSG_CHECKING([whether the python package should be compilated and installed])
    AC_ARG_ENABLE(
        python,
        AS_HELP_STRING([--disable-python],[Turn off compilation and installation of the python package.]),
        [],
        enable_python=yes
    )
    #dnl AS_VAR_COPY([ENABLE_PYTHON],enable_python)
    AC_SUBST([ENABLE_PYTHON],AS_VAR_GET(enable_python))
    AC_MSG_RESULT(AS_VAR_GET(enable_python))

    # Basic python
    AS_IF([test "$enable_python" = "no" ], [AS_VAR_SET(PYTHON,no)])
    AX_WITH_PYTHON(2.5,no)
    
    # We have python
    AS_IF([test "AS_VAR_GET(PYTHON)" != "no" -a "AS_VAR_GET(HAS_BLASLAPACK)" != "no"],[

        # Current python
        AS_VAR_SET(MYPYTHONPATH,[`AS_DIRNAME(AS_VAR_GET(PYTHON))`])
        
        # Check numpy
        AC_MSG_CHECKING([for numpy support])
        AS_VAR_GET(PYTHON) -c ["import numpy"] 2> /dev/null
        AS_IF([test "$?" = "0"],[
        
            AC_MSG_RESULT("yes")
            
            # Ok, now check install dir
            AC_ARG_VAR(PYTHONDIR,[Directory where to install the spanlib python pure module])
            AC_ARG_WITH(pythondir,
                AC_HELP_STRING(--with-pythondir=DIR,
                        [Directory where to install the spanlib python module]),
                        [case AS_VAR_GET(with_pythondir) in
                            no|yes);;
                            *)AS_VAR_SET(PYTHONDIR,AS_VAR_GET(with_pythondir))
                                AM_CONDITIONAL(LOCAL_PYTHON_INSTALL,true);;
                        esac]
            )
        ],[
            # numpy not available => no python module
            AC_MSG_RESULT("no")
            AC_SR_WARNING([You wont be able to build the python module. You need BLAS/LAPACK and numpy.])
            AS_VAR_SET(ENABLE_PYTHON,no)
        ])
    ])


    # So, for python...
    AM_CONDITIONAL([ENABLE_PYTHON],[test "AS_VAR_GET(ENABLE_PYTHON)" = "yes"])

AS_IF([test "AS_VAR_GET(ENABLE_PYTHON)" == "yes"],[

# Python site.cfg
cat > setup.cfg << EOF
######################################################################
# Python setup config file generated with ./configure 
######################################################################

[[build_ext]]
library_dirs=/usr/lib:/usr/local/lib

[[config_fc]]
f90exec=AS_VAR_GET(FC)
f77exec=AS_VAR_GET(FC)

[[blaslapack]]
blas=AS_VAR_GET(BLAS)
blas_lib=AS_VAR_GET(BLAS_LIB)
lapack=AS_VAR_GET(LAPACK)
lapack_lib=AS_VAR_GET(LAPACK_LIB)
lapack_inc=AS_VAR_GET(LAPACK_INC)

EOF

dnl build_lib=build/lib
dnl lapack95=AS_VAR_GET(LAPACK95)
dnl lapack95_lib=AS_VAR_GET(LAPACK95_LIB)
dnl lapack95_inc=AS_VAR_GET(LAPACK95_INC)
dnl lapack95_mod=AS_VAR_GET(LAPACK95_MOD)
dnl lapack95_pre=AS_VAR_GET(LAPACK95_PRE)

AC_MSG_NOTICE(Created setup.cfg)
])


])

################################################################################
################################################################################
################################################################################
