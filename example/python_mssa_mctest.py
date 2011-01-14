## Read SST
import cdms2
f = cdms2.open('../data/pacific_sst.nc')
sst = f('ssta')
f.close()


# Analysis with Spanlib
# - init
print "Creating SpAn object"
import sys, os
if os.path.exists('../src/build/lib'):sys.path.insert(0,'../src/build/lib')
import spanlib
span = spanlib.SpAn(sst)
# - MSSA
print "MSSA..."
ev = span.mssa_ev()
print ev
# - Monte-Carlo test
print 'MC test'
mcev = span.mssa_mctest(nens=20)
