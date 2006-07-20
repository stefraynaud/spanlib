import cdms
import spanlib
import MV


cdms.axis.latitude_aliases.append('Y')
cdms.axis.longitude_aliases.append('X')
cdms.axis.time_aliases.append('T')

f=cdms.open('example/data2.cdf')

s2=f('ssta',latitude=(-10,10),longitude=(110,180))
s1=f('ssta',latitude=(-15,15),longitude=(210,250))


print 'Data  in:',s1.shape,s2.shape

res = spanlib.stackData(s1,s2)

print res[0].shape

SP=spanlib.SpAn(MV.array(res[0]),weights=MV.array(res[1]))
eof,pc,ev = SP.pca()

res2 = spanlib.unStackData(eof,res[1],res[2],res[3])



print res2

import vcs

x=vcs.init()
x.plot(res2[1])
raw_input('ok?')
