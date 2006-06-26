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
axes=s.getAxisList()
## s=f('ssta',longitude=(210,215),latitude=(-2,2))
## axes=s.getAxisList()
## cdutil.setTimeBoundsMonthly(s)

## mean=cdutil.averager(s,axis='t')
## s,mean=genutil.grower(s,mean)
## s=s-mean
## s.setAxisList(axes)
print s.shape
sh=list(s.shape)

ns1=sh[-1]
ns2=sh[-2]
nt=sh[0]

## pca param
nkeep=10

## mssa param
nwindow=84
nkeep2=2

## Phases param
nphases = 8
offset = .5
firstphase = 0

mask=s.mask()

if mask is None: # If no mask creates a "dummy" mask
    mask=MV.zeros(s.shape)

M=1.-MV.greater_equal(MV.sum(mask,axis=0),1)

ns=MV.sum(MV.ravel(M))
print ns,ns1*ns2
print spanlib.pack3d.__doc__

s=MV.transpose(s)
M=MV.transpose(M)

s2=spanlib.pack3d(s.filled(1.e20),M.filled(),ns1,ns2,nt,ns)
s=MV.transpose(s)

print s2.shape

w=MV.ones((int(ns)),'f')
print w.shape


eof,pc,ev = spanlib.pca(s2,ns,nt,nkeep,w,1)
f=cdms.open('tmp_python.nc','w')
f.write(eof,id='eof')
f.write(pc,id='pc')
f.write(ev,id='ev')
f.close()

print eof.shape
print pc.shape
print ev.shape

print spanlib.unpack3d.__doc__

## print 'EOF:',eof.shape
## eof=spanlib.unpack3d(M.filled(),ns1,ns2,nkeep,eof,ns,1.e20)

#ffrec = spanlib.pcarec(eof, pc, ns, nt, nkeep, 1,nkeep)

steof, stpc, ev = spanlib.mssa(Numeric.transpose(pc), nkeep, nt, nwindow, nkeep2)


ffrec = spanlib.mssarec(steof, stpc, nkeep, nt, nkeep2, nwindow, 1,nkeep2)
print 'Mssa rec:',ffrec.shape

ffrec = spanlib.pcarec(eof, Numeric.transpose(ffrec), ns, nt, nkeep, 1,nkeep)
print 'pca rec:',ffrec.shape

#ffrec = spanlib.unpack3d(M.filled(),ns1,ns2,nt,ffrec,ns,1.e20)

phases = spanlib.phasecomp(ffrec, ns, nt, nphases, w, offset, firstphase)
phases = spanlib.unpack3d(M.filled(),ns1,ns2,nphases,phases,ns,1.e20)

import EzTemplate

x=vcs.init()
T=EzTemplate.Multi(rows=4,columns=2)

for i in range(8):
    t=T.get()
    x.plot(MV.transpose(phases)[i],t)

raw_input()

## f=cdms.open('out.nc','w')

## f.write(steof,id='steof')
## f.write(stpc,id='stpc')
## f.write(ev,id='ev')

## f.close()

