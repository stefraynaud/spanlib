################################################################################
### SpanLib
### Raynaud 2006
################################################################################
# Configure footer
################################################################################
AC_DEFUN([AC_SR_SPANLIB_FOOTER],
[
AC_SR_SET_GREENINK
echo "################################################################################"
echo -e "Now type '"AS_VAR_GET([NORMAL])"make"AS_VAR_GET([GREEN])"' to build the library.\nThen, to install it, type '"AS_VAR_GET([NORMAL])"make install"AS_VAR_GET([GREEN])"'.\nyou may need to do it as "AS_VAR_GET([NORMAL])"root"AS_VAR_GET([GREEN])" if you didnt use "AS_VAR_GET([NORMAL])"--prefix"AS_VAR_GET([GREEN])"."
AS_VAR_SET_IF(F90_EXAMPLE_TEXT,echo -e "AS_VAR_GET(F90_EXAMPLE_TEXT)")
AS_VAR_SET_IF(PYTHON_EXAMPLE_TEXT,echo -e "AS_VAR_GET(PYTHON_EXAMPLE_TEXT)")
echo "################################################################################"
AC_SR_SET_NORMALINK
])


