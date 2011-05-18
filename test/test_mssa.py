import unittest
import numpy as npy
import os, sys
from util import insert_local_path, setup_data2
insert_local_path()
sys.path.insert(0, '../lib')

from spanlib_python import SpAn

import pylab as P

class TestSequenceFunctions(unittest.TestCase):

#    def test_ssa(self):
#        nt = 300
#        t = npy.arange(nt*1.)
#        p0 = 23.
#        p1 = p0*3.6
#        p2 = p0/10.
#        data = npy.cos(t/p0*2*npy.pi)+npy.sin(t/p1*2*npy.pi+1.3)+npy.cos(t/p2*2*npy.pi+1.4)
#        span = SpAn(data)
#        span.mssa(nmssa=4)
#        self.assertTrue(npy.allclose(
#            span._mssa_raw_ev[0], 
#            npy.array([ 20.09992241,  16.45012377,  16.39094369,  16.18919368])))
#        self.assertTrue(npy.allclose(
#            span._mssa_raw_pc[0][:2,:2], 
#            npy.array([[ 6.91052733,  0.32026904],
#                [ 5.77278847,  2.92765239]])))
#        self.assertTrue(npy.allclose(
#            span._mssa_raw_eof[0][:2,:2], 
#            npy.array([[ 0.13631628, -0.12773358],
#                [ 0.13412104,  0.13178561]])))
        
#    def test_mssa(self):
#        data = setup_data2(nx=3, ny=2)
#        span = SpAn(data)
#        span.mssa(nmssa=4)
#        self.assertTrue(npy.allclose(
#            span._mssa_raw_ev[0], 
#            npy.array([ 83.10160414,  39.19375277,  36.03284783,  21.17232915])))
#        self.assertTrue(npy.allclose(
#            span._mssa_raw_pc[0][:2,:2], 
#            npy.array([[-13.66449141,   2.6225741 ],
#                [-13.28713796,   2.7381359 ]])))
#        self.assertTrue(npy.allclose(
#            span._mssa_raw_eof[0][:2,:2], 
#            npy.array([[ 0.05062128,  0.10591617],
#                [ 0.0529649 ,  0.1076283 ]])))

    def test_mssa_eof(self):
        data = setup_data2(nx=3, ny=2)
        span = SpAn(data)
        eof = span.mssa_eof(nmssa=4)
        self.assertTrue(npy.allclose(
            eof[1,1].compressed(), 
            npy.array([ 0.00751432,  0.03178418,  0.03079323,  0.01421789,  0.00330587])))

    def test_mssa_rec(self):
        data = setup_data2(nx=3, ny=2)
        span = SpAn(data)
        rec = span.mssa_rec(nmssa=4)
        self.assertTrue(npy.allclose(
            rec[1].compressed(), 
            npy.array([ 0.56983066,  0.17809844,  0.25733394,  0.66705892,  0.81728021])))
        

#    def test_mssa_parallel(self):
#        data1 = setup_data2(nx=3, ny=2)
#        data2 = data1**2
#        span = SpAn([data1, data2])
#        span.mssa(nmssa=4)
#        self.assertTrue(npy.allclose(
#            span._mssa_raw_eof[0][:2,:2], 
#            npy.array([[ 0.04292406,  0.06860169],
#                [ 0.04443118,  0.06903215]])))
  
#    def test_mssa_parallel_eof(self):
#        data1 = setup_data2(nx=3, ny=2)
#        data2 = data1**2
#        span = SpAn([data1, data2])
#        span.mssa(nmssa=4)
#        eof1, eof2 = span.mssa_eof(nmssa=4)
#        self.assertTrue(npy.allclose(
#            eof1[1,1].compressed(), 
#            npy.array([-0.00484883,  0.0204779 ,  0.01975063,  0.00289329, -0.00804804])))
#        self.assertTrue(npy.allclose(
#            eof2[1,1].compressed(), 
#            npy.array([  1.01368878e-05,   2.09268367e-02,   2.17408030e-02, 6.91380266e-03,  -3.95801005e-03])))
#            
#    def test_mssa_parallel_pc(self):
#        data1 = setup_data2(nx=3, ny=2)
#        data2 = data1**2
#        span = SpAn([data1, data2])
#        span.mssa(nmssa=4)
#        pc = span.mssa_pc(nmssa=4)
#        self.assertTrue(npy.allclose(pc[1,:2], npy.array([3.44024313,  3.60591829])))
           
        
    
if __name__ == '__main__':
    unittest.main()
