#!/usr/bin/env python

id2fc =  {	'gfortran':	('gfortran'),
		'intel':	('ifort','ifc'),
		'intelv':	('ifl'),
		'intele':	('efort','efc','ifort'),
		'intelev':	('efl'),
		'pg':		('pgf90'),
		'lahey':	('lf95'),
		'ibm':	('xlf90'),
		'compaqv':	('DF')
		'nag':	('f95'),
		'compaq':	('fort','f90'),
		'absoft':	('f90'),
		'mips':	('f90'),
		'sun':	('f90'),
		'vast':	('f90'),
		'hpux':	('f90')
 }

myId = ""
for id in id2fc.keys():
	if fc in id2fc[id]:
		myId = id
		break

print myId	

#import scipy_distutils.fcompiler
#import sys
#if sys.argv[1].lower in scipy_distutils.fcompiler.fcompiler_class.keys():
#	sys.exit(0)
#else:
#	sys.exit(1)

