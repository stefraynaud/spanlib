import os,sys
if os.path.exists('../src/build/tmp_lib'):
	sys.path.insert(0,'../src/build/tmp_lib')
import spanlib,cdms2 as cdms,MV2 as MV,numpy as npy,pylab as P

##var = npy.random.rand(50,100).astype('f')
var = npy.zeros((150,100),'f')
nt,nx = var.shape

xx = npy.arange(nx,dtype='f')
tt = npy.arange(nt,dtype='f')

xfreq = 2
tfreq = 3
for ifreq in .5,3:
	tvar = (npy.cos(2*npy.pi * tfreq * tt/(nt-1))*tt/(nt-1)).reshape((nt,1))
	xvar = npy.cos(2*npy.pi * xfreq*ifreq * xx/(nx-1)).reshape((1,nx))
	var += npy.multiply(tvar,xvar)
for i in xrange(nt):
	for j in xrange(nx):
		var[i,j] += 3*(i+j)/(nx+nt)


vv = var.copy()

cov = npy.cov(var.transpose())
print '- xcov',cov[0,:10]
xev,xeof = npy.linalg.eig(cov)
print '- xev',xev[:6].real,max(xev).real
print '- xeof**2',[(xeof[:,i]**2).sum().real for i in xrange(6)]
print '- ',npy.max(xeof[:,0].real),npy.min(xeof[:,0].real)
print '- ev sums',xev.real.sum(),cov.trace()

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

nmode = 10
span = spanlib.SpAn(var,npca=nmode)
eof = span.pca_eof()
pc = span.pca_pc()
ev = span.pca_ev()
rec= span.pca_rec(imode=-2)

print '+ ev:',ev
print '+ eof sums:',[(thiseof.filled()**2).sum() for thiseof in eof]
print '+ ',npy.max(eof[0]),npy.min(eof[0])
print '+ ev sum:',span.pca_ev(sum=True)


f=cdms.open('out.nc','w')
f.write(var,extend=False)
f.write(eof)
f.write(pc)
f.write(ev)
f.write(rec)
f.close()

P.figure(figsize=(6,8))
P.subplots_adjust(hspace=.4)
for im in xrange(nmode):
	P.subplot(nmode,2,im*2+1)
	P.plot(eof[im].filled())
	P.title('EOF %i'%(im+1))
	P.subplot(nmode,2,im*2+2)
	P.plot(pc[im].filled())
	P.title('PC %i'%(im+1))
P.savefig('out.png')

P.figure(figsize=(6,8))
P.subplots_adjust()
P.subplot(211)
P.pcolor(var.filled()-var.filled().mean(0))
P.colorbar()
P.title('Original')
P.subplot(212)
P.pcolor(rec.filled())
P.colorbar()
P.title('Rec%i'%nmode)
P.savefig('outrec.png')






