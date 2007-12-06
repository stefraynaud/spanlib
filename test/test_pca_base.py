import os,sys
if os.path.exists('../src/build/tmp_lib'):
	sys.path.insert(0,'../src/build/tmp_lib')
import spanlib,cdms2 as cdms,MV2 as MV,numpy as npy,pylab as P



var = npy.zeros((50,100),'f')
nt,nx = var.shape

xx = npy.arange(nx,dtype='f')
tt = npy.arange(nt,dtype='f')

xfreq = 3
tfreq = 2
for ifreq in 1,2:
	tvar = npy.cos(npy.pi * xfreq*ifreq * tt/(nt-1)).reshape((nt,1))
	xvar = npy.cos(npy.pi * xfreq*ifreq * xx/(nx-1)).reshape((1,nx))
	var += npy.multiply(tvar,xvar)


lon = cdms.createAxis(xx,id='lon')
lon.units = 'degrees_east'
lon.designateLongitude()

time = cdms.createAxis(tt,id='time')
time.units = 'months since 2000'
time.designateTime()

var = MV.array(var,copy=0,id='var')
var.units = 'm'
var.long_name = 'My variable'
var.setAxisList([time,lon])


span = spanlib.SpAn(var,npca=2)
eof = span.pca_eof()
print 'SHAPE EOF',eof.shape
eof.info()

#pc = span.pca_pc()

f=cdms.open('out.nc','w')
f.write(eof)
f.close()
P.pcolor(eof)
P.savefig('out.png')


