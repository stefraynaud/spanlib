import sys
sys.path.append('./build/lib.linux-i686-2.4')
import cdms
import spanlib
import vcs
import MV
import Numeric
import cdutil
import genutil

cdms.axis.latitude_aliases.append('Y')
cdms.axis.longitude_aliases.append('X')
cdms.axis.time_aliases.append('T')

f=cdms.open('../example/data2.cdf')

s=f('ssta')

SP=spanlib.SpAn(s)

eof,pc,ev = SP.pca()

print 'Done PCA, doing mssa'

steof,stpc,stev = SP.mssa(pca=True)

print 'Reconstructing'

ffrec = SP.reconstruct()

phases = spanlib.phases(ffrec)

print phases
#aaa
x=vcs.init()
x.plot(phases)
raw_input('ok')
x.clear()
