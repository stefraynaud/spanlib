import os,sys
if os.path.exists('../src/build/tmp_lib'):
	sys.path.insert(0,'../src/build/tmp_lib')
import spanlib,cdms2 as cdms,MV2 as MV,numpy as npy,pylab as P


var = npy.zeros((50,100),'f')
nt,nx = var.shape

xx = npy.arange(nx,dtype='f')
tt = npy.arange(nt,dtype='f')

xfreq = 2
tfreq = 3
for ifreq in .5,3:
	tvar = (npy.cos(2*npy.pi * tfreq * tt/(nt-1))*tt/(nt-1)).reshape((nt,1))
	xvar = npy.cos(2*npy.pi * xfreq*ifreq * xx/(nx-1)).reshape((1,nx))
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

f=cdms.open('var.nc','w')
f.write(var)
f.close()

nmode = 6
span = spanlib.SpAn(var,npca=nmode)
eof = span.pca_eof()
print 'SHAPE EOF',eof.shape
pc = span.pca_pc()


pc.setAxis(1,var.getAxis(0))
print 'VAR',var.getTime()
print 'PC',pc.getTime()

f=cdms.open('out.nc','w')
f.write(var,extend=False)
f.write(eof)
f.write(pc)
f.close()

P.figure(figsize=(8,10))
for im in xrange(nmode):
	P.subplot(nmode*2,2,im*2+1)
	P.plot(eof[im].filled())
	P.title('EOF %i'%(im+1))
	P.subplot(nmode*2,2,im*2+2)
	P.plot(pc[im].filled())
	P.title('PC %i'%(im+1))
P.savefig('out.png')



