# File: example1.py
#
# This file is part of the SpanLib library.
# Copyright (C) 2006  Charles Doutiraux, Stephane Raynaud
# Contact: stephane dot raynaud at gmail dot com
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

###################################################################
# In this example, we perform a pre-PCA to reduce the number of
# d-o-f, then we perform an MSSA to extract the first oscillation
# (first par of modes). Finally, we compute phase composites
# to represent the oscillation over its cycle.
###################################################################
print "#################################################"
print "# PCA+MSSA+Phase_composites of 1st oscillation. #"
print "# Then reconstructions, and plots:              #"
print "# - 1 phase over 2                              #"
print "# - 1 time series                               #"
print "#################################################"


# Needed modules
print 'Importing needed modules...'
import sys
import cdms
import vcs
import MV
import Numeric
import cdutil
import genutil

# Current version of spanlib is prioritary
sys.path.insert(0,'../src/build/tmp_lib')
import spanlib

# We tell cdms that we have longitude, latitude and time
cdms.axis.latitude_aliases.append('Y')
cdms.axis.longitude_aliases.append('X')
cdms.axis.time_aliases.append('T')

# Simply open the netcdf file
print "Open file"
f=cdms.open('data2.cdf')

# Retrieve data
print "Read the whole dataset"
s=f('ssta')

# Create the analysis object
print "Creating SpAn object"
SP=spanlib.SpAn(s)

# Perform a preliminary PCA+MSSA
# (equivalent to simple use SP.mssa(pca=True) later)
print "PCA..."
eof,pc,ev = SP.pca()

# MSSA on PCA results
print 'MSSA...'
steof,stpc,stev = SP.mssa()

# Phase composites of first two MSSA modes
print 'Phase composites...'
out = SP.reconstruct(phases=True,nphases=6,imode=-2)

# Plot 1 phase over two, then a time series
print "Now, plot!"
x=vcs.init()
#x.open()
slices = range(0,out.shape[0])
nslices = len(slices)
import EzTemplate
columns = 3
rows = (nslices-1)/columns+1
print 'R,C',rows,columns
T=EzTemplate.Multi(rows=rows,columns=columns) # Nrow added 1 for original data row
#templ = T.get()
mn,mx=-1,1
print 'Min, max:',mn,mx,vcs.minmax(out)
levels = vcs.mkscale(mn,mx)
levels.insert(0,-1.e20) # extension left side
levels.append(1.e20) # extension right side
colors = vcs.getcolors(levels)
iso = x.createisofill('spanlib')
iso.levels = levels
iso.fillareacolors = colors
#iso.list()
#x.plot(s,templ,iso,ratio='autot')
#templ=T.get() # dummy space
#templ=T.get() # dummy space
f=cdms.open('tmp.nc','w')
f.write(out,id='out',typecode='f')
f.close()
for i in slices:
    print i
    templ = T.get(font=0)
    #templ.data.list()
    x.plot(out[i],templ,iso,ratio='autot',bg=1,title="Phase composites of the first MSSA oscillation")
##     raw_input('map out %i/%i ok?' % ( i+1 , out.shape[0]))
##     x.clear()
x.postscript('crap')
x.showbg()
raw_input('map out ok?')
x.clear()
x.plot(out[:,30,80],title="Cycle of the ocillation")
raw_input('Time series at center of bassin ok?')
x.clear()
