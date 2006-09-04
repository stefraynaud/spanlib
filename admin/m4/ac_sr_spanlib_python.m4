################################################################################
### SpanLib, Raynaud 2006
################################################################################
AC_DEFUN([AC_SR_SPANLIB_PYTHON],[

	# Tools
	AC_CHECK_PROG(PERL,perl,perl,no)
	AC_PROG_CC()

	# Basic python
	AX_WITH_PYTHON(2.4,no)

	# Current python
	AS_IF([test "AS_VAR_GET(PYTHON)" != "no"],
		AS_VAR_SET(MYPYTHONPATH,[`AS_DIRNAME(AS_VAR_GET(PYTHON))`]))

	# Pyfort support
	AC_SR_PYFORT()

	# Python modules for spanlib
	AC_MSG_CHECKING([for Numeric and MV module support])
	AS_VAR_GET(PYTHON) -c ["import Numeric;import MV"] 2> /dev/null
	AS_IF([test "$?" = "0"],
		[AC_MSG_RESULT([yes])
			AS_VAR_SET(HAS_PYT_MOD,[true])],
		[AC_MSG_RESULT([no])
			AS_VAR_SET(HAS_PYT_MOD,[false])]
	)

	# So..
	AS_IF([`AS_VAR_SET(HAS_PYT_MOD) && test "AS_VAR_GET(PYFORT)" != "no" -a \
			"AS_VAR_GET(PERL)" != "no" -a "AS_VAR_GET(HAS_BLASLAPACK)" != "no"`],
			[AS_VAR_SET(WITH_PYTHON,"yes")])
	AM_CONDITIONAL([WITH_PYTHON],AS_VAR_TEST_SET(WITH_PYTHON))
	AS_IF(AS_VAR_TEST_SET(WITH_PYTHON),,
		[AC_SR_WARNING([You wont be able to build the python module.
You need perl, BLAS/LAPACK, pyfort and Numeric, MV python modules from CDAT.])]
	)

	# Python modules for python examples
	AC_MSG_CHECKING([for cdms and cvs module support])
	AS_IF(AS_VAR_TEST_SET(WITH_PYTHON),[
		AS_VAR_GET(PYTHON) -c ["import cdms;import vcs"] 2> /dev/null
		AS_IF([test "$?" = "0"],
			[AS_VAR_SET(WITH_PYTHON_EXAMPLE,"yes")])
	])
	AS_IF(AS_VAR_TEST_SET(WITH_PYTHON_EXAMPLE),
		[AC_MSG_RESULT([yes])],
		[AC_MSG_RESULT([no])
			AC_SR_WARNING([You wont be able to run the python example.
You need cdms and vcs python modules from CDAT.])])
	AM_CONDITIONAL([WITH_PYTHON_EXAMPLE],AS_VAR_TEST_SET(WITH_PYTHON_EXAMPLE))

	# Vcdat
	AS_VAR_SET_IF(MYPYTHONPATH,
		AC_CHECK_PROG(VCDAT,vcdat,vcdat,no,AS_VAR_GET(MYPYTHONPATH)),
		[AC_CHECK_PROG(VCDAT,vcdat,vcdat,no)])
	AS_VAR_SET(HAS_VCDAT,[`test "AS_VAR_GET(VCDAT)" != "no"`])


])
################################################################################
################################################################################
################################################################################
