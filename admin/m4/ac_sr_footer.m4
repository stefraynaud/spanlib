################################################################################
### SpanLib
### Raynaud 2006
################################################################################
# Configure footer
################################################################################
AC_DEFUN([AC_SR_FOOTER],
[
AC_SR_SET_GREENINK
echo "################################################################################"
echo -en "Now type 'make' to build the library.\n\
Then, to install it, type 'make install'\n\
(you may need to do it as root if you didnt use --prefix)."
AS_VAR_SET_IF(ac_cv_text_example,[echo -e [\nAS_VAR_GET(ac_cv_text_example)]],[echo -e \n])
echo "################################################################################"
AC_SR_SET_NORMALINK
])


