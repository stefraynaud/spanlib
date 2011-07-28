import unittest
import numpy as npy
import os, sys
from util import insert_local_path, setup_data2, setup_data1, setup_data0
insert_local_path()
sys.path.insert(0, '../lib')

from spanlib_python import SpAn

import pylab as P

class TestSequenceFunctions(unittest.TestCase):

    def test_ssa(self):
        data = setup_data0()
        span = SpAn(data)
        span.mssa(nmssa=4)
        self.assertTrue(npy.allclose(
            span._mssa_raw_ev[0], 
            npy.array([ 20.25648796,  16.57547325,  16.51572589,  16.31088181])))
        self.assertTrue(npy.allclose(
            span._mssa_raw_pc[0][:2,:2], 
            npy.array([[ 6.9105607 ,  0.32054539],
       [ 5.77266597,  2.92740569]])))
        self.assertTrue(npy.allclose(
            span._mssa_raw_eof[0][:2,:2], 
            npy.array([[ 0.13637961, -0.12778613],
       [ 0.13417614,  0.13187513]])))
        
    def test_mssa(self):
        data = setup_data2(nx=3, ny=2)
        span = SpAn(data)
        span.mssa(nmssa=4)
        self.assertTrue(npy.allclose(
            span._mssa_raw_ev[0], 
            npy.array([ 84.65284564,  39.94531834,  36.70936929,  21.49023578])))
        self.assertTrue(npy.allclose(
            span._mssa_raw_pc[0][:2,:2], 
            npy.array([[-13.66377591,   2.62887011],
       [-13.28622293,   2.74443558]])))
        self.assertTrue(npy.allclose(
            span._mssa_raw_eof[0][:2,:2], 
            npy.array([[ 0.05067246,  0.10596479],
       [ 0.05300994,  0.10764546]])))

    def test_mssa_eof(self):
        data = setup_data2(nx=3, ny=2)
        span = SpAn(data)
        eof = span.mssa_eof(nmssa=4)
        self.assertTrue(npy.allclose(
            eof[1,1].compressed(), 
            npy.array([ 0.10764546, -0.0763087 , -0.07211439,  0.05060925,  0.128224  ])))

    def test_mssa_eof_restack(self):
        data = setup_data2(nx=3, ny=2)
        span = SpAn(data)
        steof = span.mssa_eof(nmssa=4)
        stacked_eof = npy.ascontiguousarray(span[0].restack(steof, scale=False))
        self.assertTrue(npy.allclose(steof[1,1].compressed(), stacked_eof[:,1,1]))
        raw_eof = span._mssa_raw_eof[0].reshape(stacked_eof.shape)
        self.assertTrue(npy.allclose(stacked_eof, raw_eof))
        
    def test_mssa_ec(self):
        data = setup_data2(nx=3, ny=2)
        span = SpAn(data)
        span.mssa(nmssa=4)
        ec = span.mssa_ec(xeof=span.mssa_eof(), xdata=data)
        pc = span.mssa_pc()
        self.assertTrue(npy.allclose(ec, pc))
#
    def test_pca_mssa_ec(self):
        data = setup_data2(nx=30, ny=20)
        span = SpAn(data)
        span.mssa(nmssa=4)
        ec = span.mssa_ec(xeof=span.mssa_eof(), xdata=data)
        pc = span.mssa_pc()
        self.assertTrue(npy.allclose(ec, pc))

    def test_ssa_rec(self):
        data = setup_data0()
        span = SpAn(data)
        rec = span.mssa_rec(nmssa=6)
        self.assertTrue(npy.allclose(rec[:4], 
            npy.array([  2.14699609,  1.3861899 ,  2.70336958,  0.70994501])))
            
    def test_ssa_xrec(self):
        data = setup_data0()
        span = SpAn(data, nmssa=6)
        rec = span.mssa_rec()
        eof = span.mssa_eof()
        pc = span.mssa_pc()
        xrec = span.mssa_rec(xeof=eof, xpc=pc)
        self.assertTrue(npy.allclose(rec,xrec))

    def test_mssa_rec(self):
        data = setup_data2(nx=3, ny=2)
        span = SpAn(data)
        rec = span.mssa_rec(nmssa=4)
        self.assertTrue(npy.allclose(
            rec[1].compressed(), 
            npy.array([ 0.56972385,  0.17638098,  0.2556532 ,  0.66645171,  0.81734079])))
        
    def test_mssa_xrec(self):
        data = setup_data2(nx=3, ny=2)
        span = SpAn(data)
        rec = span.mssa_rec(nmssa=4)
        xeof = span.mssa_eof()
        xpc = span.mssa_pc()
        xrec = span.mssa_rec(xeof=xeof, xpc=xpc)
        self.assertTrue(npy.allclose(rec, xrec))
        
    def test_mssa_phases(self):
        nt = 100.
        nx = 5
        data = npy.cos(npy.linspace(0, 5*6.2, nt*nx)*npy.pi*2).reshape((nx, nt)).T
        data += npy.cos(npy.linspace(0, 5*3.4, nt*nx)*npy.pi*2).reshape((nx, nt)).T
        span = SpAn(data)
        rec = span.mssa_rec(imodes=[-1], phases=8)
        self.assertTrue(npy.allclose(
            rec[:2,:2].ravel(), 
            npy.array([ 0.5525298 ,  0.32848071,  0.15829024,  1.19338772])))        

    def test_mssa_parallel(self):
        data1 = setup_data2(nx=3, ny=2)
        data2 = data1**2
        span = SpAn([data1, data2])
        span.mssa(nmssa=4)
        self.assertTrue(npy.allclose(
            span._mssa_raw_eof[0][:2,:2], 
            npy.array([[ 0.04298026,  0.06864925],
       [ 0.04448111,  0.06906854]])))
  
    def test_mssa_parallel_eof(self):
        data1 = setup_data2(nx=3, ny=2)
        data2 = data1**2
        span = SpAn([data1, data2])
        span.mssa(nmssa=4)
        eof1, eof2 = span.mssa_eof(nmssa=4)
        self.assertTrue(npy.allclose(
            eof1[1,1].compressed(), 
            npy.array([ 0.06906854, -0.05710377, -0.05434486,  0.02958842,  0.08265851])))
        self.assertTrue(npy.allclose(
            eof2[1,1].compressed(), 
            npy.array([ 0.05833788, -0.06212552, -0.06297591,  0.03415673,  0.08703455])))
            
    def test_mssa_parallel_pc(self):
        data1 = setup_data2(nx=3, ny=2)
        data2 = data1**2
        span = SpAn([data1, data2])
        span.mssa(nmssa=4)
        pc = span.mssa_pc(nmssa=4)
        self.assertTrue(npy.allclose(pc[1,:2], npy.array([ 3.44881773,  3.61397627])))
 

    def test_pca_mssa(self):
        data = setup_data2(nx=30, ny=20)
        span = SpAn(data, nmssa=2)
        self.assertTrue(npy.allclose(
            span.mssa_ev(), 
            npy.array([ 2482.65882191,  2355.53413819])))
        pc = span.mssa_pc()
        self.assertTrue(npy.allclose(
            span.mssa_pc()[:,:2].ravel(), 
            npy.array([-19.47701848, -20.97653195,   2.99480916,  -0.20444392])))
        self.assertAlmostEqual(span.mssa_eof()[1,2,3,4], 0.0041007559488027926)
        
    def test_pca_mssa_rec(self):
        data = setup_data2(nx=30, ny=20)
        span = SpAn(data)
        rec = span.mssa_rec(nmssa=10)
        self.assertTrue(npy.allclose(
            rec[2,:2,:2].compressed(), 
            npy.array([ 0.36712853,  0.28310489,  0.27881545])))

    def test_pca_mssa_xrec(self):
        data = setup_data2(nx=30, ny=20)
        span = SpAn(data)
        rec = span.mssa_rec(nmssa=4)
        xeof = span.mssa_eof()
        xpc = span.mssa_pc()
        xrec = span.mssa_rec(xeof=xeof, xpc=xpc)
        self.assertTrue(npy.allclose(rec, xrec))
        
        
#    def test_mssa_mctest(self):
#        data = setup_data1(nx=5, masked=False)
#        span = SpAn(data)
#        ev, evmin, evmax = span.mssa_ev(nmssa=4, mctest=True, tfreq=4)
   
if __name__ == '__main__':
    unittest.main()
