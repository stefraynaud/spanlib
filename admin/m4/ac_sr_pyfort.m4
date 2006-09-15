################################################################################
# Pyfort with BLAS/LAPACK
# SpanLib, Raynaud 2006
################################################################################
# AC_SR_STRIPFLAGS([FLAGS],[OUTPUT_VARIABLE])
# Remove -L and -l
AC_DEFUN([AC_SR_PYFORT_STRIPFLAGS],
[
	AS_VAR_SET($2,`echo $1 | sed -r 's/(^| )-(L|l)/\1/g'`)
])
# Setup library dir and name variables PYFORT_DIRS and PYFORT_LIBS
AC_DEFUN([AC_SR_PYFORT],
[
	# Full path to pyfort
	AC_ARG_VAR(PYFORT,[Absolute path to pyfort executable])
	AS_VAR_SET_IF([PYFORT],[:],
		[AC_PATH_PROG(PYFORT,pyfort,no,AS_VAR_GET(MYPYTHONPATH))])
	AM_CONDITIONAL([HAS_PYFORT],[`test "AS_VAR_GET(PYFORT)" != "no"`])

	# F90 compiler id
	AS_IF([`test "AS_VAR_GET(PYFORT)" != "no"`],
		AS_IF([scripts/check_fortran.py AS_VAR_GET(FC)],
			AS_VAR_SET(ac_cv_goodfcid,yes),[
			AC_SR_WARNING([You f90 compiler (AS_VAR_GET(FC)) is not in the available list with pyfort.
You wont be able to build the python package.])
		AM_CONDITIONAL([HAS_PYFORT],[`false`])
	]))

	# Blas and Lapack libraries and directories
	AC_SR_PYFORT_STRIPFLAGS("AS_VAR_GET(LAPACK) AS_VAR_GET(BLAS)",PYFORT_LIBS)
	AC_SUBST(PYFORT_LIBS)
	AC_SR_PYFORT_STRIPFLAGS("AS_VAR_GET(LAPACK_LIB) AS_VAR_GET(BLAS_LIB)",PYFORT_LIBDIRS)
	AC_SUBST(PYFORT_LIBDIRS)
	AC_SR_PYFORT_STRIPFLAGS("AS_VAR_GET(LAPACK_INC)",PYFORT_INCDIRS)
	AC_SUBST(PYFORT_INCDIRS)

	# Set manual path to install the python module
	AC_ARG_VAR(PYTHONDIR,[Directory where to install the spanlib python module])
	AC_ARG_WITH(pythondir,
		AC_HELP_STRING(--with-pythondir=DIR,
				[Directory where to install the spanlib python module]),
				[case AS_VAR_GET(with_pythondir) in
					no|yes);;
					*)AS_VAR_SET(PYTHONDIR,AS_VAR_GET(with_pythondir));;
				esac]
	)

	# Build directory
	AS_IF([`test "AS_VAR_GET(PYFORT)" != "no"`],
	[
		AC_MSG_CHECKING([the generic directory name for building python libraries])
		AS_VAR_SET(PYFORT_BUILD_DIR,
			[`AS_VAR_GET(PYTHON) -c ["import sys;from distutils.util import get_platform ; print \"lib.\"+get_platform()+\"-\"+sys.version[0:3]"]`])
		AC_MSG_RESULT(AS_VAR_GET(PYFORT_BUILD_DIR))
	])

	# Build option and default install path
	AS_VAR_SET_IF(PYTHONDIR,
		AS_VAR_SET(PYFORT_BUILD_OPT,[-b]),
		[
			AS_VAR_SET(PYFORT_BUILD_OPT,[-b])
			AS_IF([`test "AS_VAR_GET(PYFORT)" != "no"`],[
					AC_MSG_CHECKING([where is the default place for python packages])
					AS_VAR_SET(PYTHONDIR,
						[`AS_VAR_GET(PYTHON) -c ["from distutils import sysconfig; print sysconfig.get_python_lib(1,0)"]`])
					AC_MSG_RESULT(AS_VAR_GET(PYTHONDIR))
				],
				AS_VAR_SET(PYTHONDIR,"")
			)
		]
	)

	AC_SUBST(PYFORT_BUILD_DIR)
	AC_SUBST(PYFORT_BUILD_OPT)

])
################################################################################
