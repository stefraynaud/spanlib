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
import sys
import cdms
import spanlib
import vcs
import MV
import Numeric
import cdutil
import genutil

# We tell cdms that we have longitude, latitude and time
cdms.axis.latitude_aliases.append('Y')
cdms.axis.longitude_aliases.append('X')
cdms.axis.time_aliases.append('T')

# Simply open the netcdf file
print "Open file"
f=cdms.open('../example/data2.cdf')

# Retrieve data
print "Read two different regions"
s=f('ssta',time=slice(0,120))

S# Create the analysis object
print "Creating SpAn object"
P=spanlib.SpAn(s)

# Perform a preliminary PCA+MSSA
# (equivalent to simple use SP.mssa(pca=True) later)
print "PCA..."
eof,pc,ev = SP.pca()

# MSSA on PCA results
print 'MSSA...'
steof,stpc,stev = SP.mssa()

# Phase composites of first two MSSA modes
print 'Phase composites...'
out = SP.reconstruct(phases=True,nphases=16,end=2)

# Plot 1 phase over two, then a time series
print "Now, plot!"
x=vcs.init()
for i in range(0,out.shape[0],out.shape[0]/10):
    x.plot(out[i])
    raw_input('map out %i/%i ok?' % ( i+1 , out.shape[0]))
    x.clear()
for i in range(0,out.shape[0],out.shape[0]/10):
    x.plot(s[i]-out[i])
    raw_input('residual noise map out %i/%i ok?' % ( i+1 , out.shape[0]))
    x.clear()
x.plot(out[:,30,80])
raw_input('Time series at center of bassin ok?')
x.clear()
