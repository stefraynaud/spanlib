################################################################################
### BLAS and LAPACK
### Raynaud 2006
################################################################################
AC_DEFUN([AC_SR_SPANLIB_DOC],[

	AC_ARG_ENABLE(
		doc-update,
		[  --enable-doc-update - Enable documentation update],
		,[AS_VAR_SET(enableval,"no")]
	)

	AS_IF([test "AS_VAR_GET(enableval)" == "yes"],
		[
			#docbook (xslt) processor support
			AC_SR_DOCBOOK(:,:)
			# We need perl and xslt processor
			AM_CONDITIONAL([HAS_DOC_SUPPORT],
				[test "AS_VAR_GET(PERL)" != "no" -a "AS_VAR_GET(XSLTPROC_WORKS)" == "yes"])
		],[AM_CONDITIONAL([HAS_DOC_SUPPORT],[false])]
	)

])
################################################################################
################################################################################
################################################################################

