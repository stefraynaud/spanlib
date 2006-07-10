################################################################################
### BLAS and LAPACK
### Raynaud 2006
################################################################################
AC_DEFUN([AC_SR_SPANLIB_PYTHON],[

	# Perl & Pyfort
	AC_CHECK_PROG(PERL,perl,perl)
	AC_PROG_CC()
	AX_WITH_PYTHON(2.4)
##	AM_PATH_PYTHON(2.4,,:)
	AC_SR_PYFORT()
	AS_VAR_SET(ac_cv_py,
		[`test "AS_VAR_GET(PYFORT)" != "no" -a "AS_VAR_GET(PERL)" != "no"`])
	AM_CONDITIONAL([WITH_PYTHON],AS_VAR_GET(ac_cv_py))
	AS_IF(AS_VAR_GET(ac_cv_py),,
		[AC_SR_WARNING([You wont be able to build the python module.
You need perl and pyfort (cdat from python).])])

	## Blas/Lapack
	AS_IF(AS_VAR_GET(ac_cv_py),
		AM_CONDITIONAL([WITH_PYTHON],
			[test "AS_VAR_GET(HAS_BLASLAPACK)" != "no"]),,)

])
################################################################################
################################################################################
################################################################################
