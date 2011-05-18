import unittest
import numpy as npy
import os, sys
from util import insert_local_path, pca_numpy, setup_data1, setup_data2
insert_local_path()
sys.path.insert(0, '../lib')

from spanlib_python import SpAn



class TestSequenceFunctions(unittest.TestCase):
    
        
    def setup_span3(self):
        var1 = setup_data1()
        var2 = setup_data1(nx=45, tfreq=4, xfreq=5)
        var3 = setup_data1(nt=120, nx=20)
        return SpAn([[var1, var2], var3])
        
        
#    def test_init(self):
#        span = self.setup_span3()
#        input = span[0][0].rescale(span[0].stacked_data[:99], copy=True).T
#        var1 = setup_data1()
#        output = var1.compress(~var1[0].mask, axis=1)
#        self.assertTrue(npy.allclose(input, output))
#
#    def test_pca(self):
#        # spanlib
#        span = SpAn(setup_data1(nt=70, nx=50))
#        span.pca()
#        sp_eof = span._pca_raw_eof[0]
#        sp_pc = span._pca_raw_pc[0]
#        sp_ev = span._pca_raw_ev[0]
#        # numpy
#        np_eof, np_pc, np_ev = pca_numpy(span[0].stacked_data.T, span.npca())
#        # checks
#        self.assertTrue(npy.allclose(sp_ev, np_ev))
#        signs = npy.sign(sp_eof[0])*npy.sign(np_eof[0])
#        np_eof *= npy.resize(signs, np_eof.shape)  
#        np_pc *= npy.resize(signs, np_pc.shape)
#        self.assertTrue(npy.allclose(sp_eof, np_eof))
#        self.assertTrue(npy.allclose(sp_pc, np_pc))
#    
#    def test_pca_ev(self):
#        span = SpAn(setup_data1(nt=70, nx=150))
#        span.pca()
#        ev = span.pca_ev()
#        self.assertTrue(npy.allclose(ev, npy.array(
#        [ 16.41185619,  11.849788  ,   2.36332927,   0.94659187,
#         0.58360009,   0.29865376,   0.16189405,   0.1470843 ,
#         0.11100314,   0.0948228 ])))
#         
#    def test_pca_eof(self):
#        span = SpAn(setup_data1(nt=70, nx=150))
#        span.pca_eof()
#        raw_eof = span.pca_eof(raw=True)
#        self.assertTrue(npy.allclose(raw_eof[:2,:2], npy.array(
#            [[ 0.15702435,  0.14647224],
#            [ 0.05325682,  0.05082613]])))
#        eof = span.pca_eof()
#        self.assertTrue(npy.allclose(eof[:2,:3].filled(0.),  npy.array(
#            [[ 0.15702435,  0.        ,  0.14647224],
#            [ 0.05325682,  0.        ,  0.05082613]])))
#            
#    def test_pca_pc(self):
#        span = SpAn(setup_data1(nt=70, nx=150))
#        span.pca_eof()
#        pc = span.pca_pc()
#        self.assertTrue(npy.allclose(pc[:2,:2],  npy.array(
#            [[  3.24705271,   3.49697418],
#            [-10.34830743, -10.01240894]])))
#        
#    def test_pca_ec(self):
#        data = setup_data1(nt=70, nx=150)
#        span = SpAn(data)
#        span.pca_eof()
#        ec = span.pca_ec()
#        pc = span.pca_pc()
#        self.assertTrue(npy.allclose(ec, pc))
#        ece = span.pca_ec(xeof=span.pca_eof())
#        self.assertTrue(npy.allclose(ece, pc))
#        ecd = span.pca_ec(xdata=data)
#        self.assertTrue(npy.allclose(ecd, pc))
#        ecrd = span.pca_ec(xdata=span[0].stacked_data,raw=True)
#        self.assertTrue(npy.allclose(ecrd, pc))
#        ecre = span.pca_ec(xeof=span._pca_raw_eof[0],raw=True)
#        self.assertTrue(npy.allclose(ecrd, pc))
#     
#    def test_pca_parallel(self):
#        span = self.setup_span3()
#        span.pca()
#        # eof
#        eof = npy.array([[ 0.05771058,  0.20549028],
#       [ 0.09507791, -0.02198945]])
#        self.assertTrue(npy.allclose(span._pca_raw_eof[0][:110:109,:2], eof))
#        # pc
#        pc = npy.array([[-15.94523611,   1.55413017],
#       [  8.18719366,  -1.64042587]])
#        self.assertTrue(npy.allclose(span._pca_raw_pc[0][:110:109,:2], pc))
#        # ev
#        ev = npy.array([ 51.9014321 ,  23.96282796])
#        self.assertTrue(npy.allclose(span._pca_raw_ev[0][:2], ev))
#        # ev sum
#        self.assertAlmostEqual(span._pca_ev_sum[0], 95.226299860954697)
#  
#
#    def test_pca_rec(self):
#        span = self.setup_span3()
#        (rec1, rec2), rec3 = span.pca_rec()
#        self.assertAlmostEqual(rec1.sum(), 18464.31988812385)
#        self.assertAlmostEqual(rec2.sum(), 8253.6128613672063)
#        self.assertAlmostEqual(rec3.sum(), 3557.322411248088)
      
#    def test_pca_ndim3(self):
#        span = SpAn(setup_data2())
#        span.pca()
#        self.assertTrue(npy.allclose(
#            span._pca_raw_eof[0][190:193,5], 
#            npy.array([-0.03623882, -0.02312414, -0.00942417])))
    
    
if __name__ == '__main__':
    unittest.main()
