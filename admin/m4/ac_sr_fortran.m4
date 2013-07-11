################################################################################
# F90 compiler
################################################################################
AC_DEFUN([AC_SR_FORTRAN],
[


	# Main program
AC_FC_WRAPPERS
AC_LANG(Fortran)
	AC_PROG_FC
	AS_VAR_SET_IF(FC,,
		AC_SR_ERROR([No Fortran 90 compiler available on this machine.
Please use FC to specify it or update your environnement variable PATH or
install a Fortran 90 compiler.])
	)
AC_FC_SRCEXT(f90)
AC_FC_FREEFORM()
	
	

])
