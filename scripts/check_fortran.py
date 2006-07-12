#!/usr/bin/env python

from Pyfort import fortran_compiler
import sys

if sys.argv[1] in fortran_compiler.compiler_ids:
	sys.exit(0)
else:
	sys.exit(1)
