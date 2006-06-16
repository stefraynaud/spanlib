#! /bin/bash

echo "###########################################"

echo "Generating m4 include file..."
MY_M4="my ac_sr_fortran ac_sr_blas"
# acx_lapack"
rm -f acinclude.m4
for m4 in $MY_M4 ; do
	cat admin/m4/$m4.m4 >> acinclude.m4
done

echo "Aclocal..."
aclocal

echo "Automake..."
automake --add-missing --force-missing

echo "Autoconf..."
autoconf


echo "###########################################"
echo "Now you can use ./configure"
echo "###########################################"
