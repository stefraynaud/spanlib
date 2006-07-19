################################################################################
### SpanLib, Raynaud 2006
################################################################################
AC_DEFUN([AC_SR_SPANLIB_PYTHON],[

	# Tools
	AC_CHECK_PROG(PERL,perl,no)
	AC_PROG_CC()

	# Basic python
	AX_WITH_PYTHON(2.4,no)

	# Current python
	AS_IF([test "AS_VAR_GET(PYTHON)" != "no"],
		AS_VAR_SET(MYPYTHONPATH,[`AS_DIRNAME(AS_VAR_GET(PYTHON))`]))

	# Pyfort
	AC_SR_PYFORT()
	AS_VAR_SET(ac_cv_py,
		[`test "AS_VAR_GET(PYFORT)" != "no" -a "AS_VAR_GET(PERL)" != "no"`])
	AM_CONDITIONAL([WITH_PYTHON],AS_VAR_GET(ac_cv_py))
	AS_IF(AS_VAR_GET(ac_cv_py),,
		[AC_SR_WARNING([You wont be able to build the python module.
You need perl and pyfort (cdat from python).])])

	# Blas/Lapack dependency
	AS_IF(AS_VAR_GET(ac_cv_py),
		AM_CONDITIONAL([WITH_PYTHON],
			[test "AS_VAR_GET(HAS_BLASLAPACK)" != "no"]),,)

	# Cdms, cvs
	AC_MSG_CHECKING([for cdms and cvs support])
	AS_VAR_GET(PYTHON) -c ["import cdms;import vcs"] 2> /dev/null
	AS_IF([test "$?" = "0"],
		[AC_MSG_RESULT([yes])
		AS_VAR_SET(HAS_CDMSCVS,`true`)],
		[AC_MSG_RESULT([no])
		AS_VAR_SET(HAS_CDMSCVS,`false`)]
	)

	# Vcdat
	AS_VAR_SET_IF(MYPYTHONPATH,
		[AC_CHECK_PROG(VCDAT,vcdat,vcdat,no,AS_VAR_GET(MYPYTHONPATH))
		AC_SR_WARNING([with path AS_VAR_GET(MYPYTHONPATH)])],
		[AC_CHECK_PROG(VCDAT,vcdat,vcdat,no)])
	AS_VAR_SET(HAS_VCDAT,[`test "AS_VAR_GET(VCDAT)" != "no"`])


])
################################################################################
################################################################################
################################################################################
