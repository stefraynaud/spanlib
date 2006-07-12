################################################################################
# Spanlib, Pyfort with BLAS/LAPACK
# Raynaud 2006
################################################################################
# AC_SR_STRIPFLAGS([FLAGS],[OUTPUT_VARIABLE])
# Remove -L and -l
AC_DEFUN([AC_SR_PYFORT_STRIPFLAGS],
[
	AS_VAR_SET($2,`echo $1 | sed -r 's/(^|\s)-(L|l)/\1/g'`)
])
# Setup library dir and name variables PYFORT_DIRS and PYFORT_LIBS
AC_DEFUN([AC_SR_PYFORT],
[
	# Current python
	AS_VAR_SET_IF(MYPYTHON,
		AS_VAR_SET(MYPYTHONPATH,[AS_DIRNAME(AS_VAR_GET(MYPYTHON))]))

	# Full path to pyfort
	AS_VAR_SET_IF(MYPYTHONPATH,
		[AC_PATH_PROG(PYFORT,pyfort,no,AS_VAR_GET(MYPYTHONPATH))],
		[AC_PATH_PROG(PYFORT,pyfort,no)])
	AM_CONDITIONAL([HAS_PYFORT],[test "AS_VAR_GET(PYFORT)" != "no"])

	# F90 compiler id


	# Blas and Lapack libraries and directories
	AC_SR_PYFORT_STRIPFLAGS("AS_VAR_GET(BLAS) AS_VAR_GET(LAPACK)",PYFORT_LIBS)
	AC_SUBST(PYFORT_LIBS)
	AC_SR_PYFORT_STRIPFLAGS("AS_VAR_GET(BLAS_LIBDIR) AS_VAR_GET(LAPACK_LIBDIR)",PYFORT_DIRS)
	AC_SUBST(PYFORT_DIRS)

	# Set manual path to install the python module
	AC_ARG_VAR(PYTHONDIR,[Directory where to install the spanlib python module])
	AC_ARG_WITH(pythondir, dnl
		AS_HELP_STRING(--with-pythondir=DIR, dnl
				[Directory where to install the spanlib python module]),
				[case AS_VAR_GET(with_pythondir) in
					no|yes);;
					*)AS_VAR_SET(PYTHONDIR,AS_VAR_GET(with_pythondir));;
				esac]
	)

	# Build directory
	AS_IF([test "AS_VAR_GET(PYFORT)" != "no"],
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
			AS_IF([test AS_VAR_GET(PYFORT) != "no"],[
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
