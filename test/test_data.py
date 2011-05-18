import unittest
import numpy as N
import os, sys
from util import insert_local_path
insert_local_path()
sys.path.insert(0, '../lib')

from spanlib_python import Data


class TestSequenceFunctions(unittest.TestCase):
    
    def setup_data(self):
        self.d = Data(self.data, nvalid=20, weights=self.weights, norm=self.norm)

    def setUp(self):
        nt, ny, nx = 100, 5, 8
        self.shape = nt, ny, nx
        self.data = N.ma.arange(nt*ny*nx, dtype='d').reshape(self.shape)
        self.data[:, 4, 6:] = N.ma.masked # test fixed mask
        self.data[:95, 1, 0] = N.ma.masked # test nvalid
        self.weights = N.arange(ny*nx).reshape((ny, nx))/10.
        self.norm = 100.

    def test_init(self):
        self.setup_data()

    def test_nstot(self):
        self.setup_data()
        self.assertEqual(self.d.ns, 36)
        
    def test_pack(self):
        self.setup_data()
        self.assertAlmostEqual(self.d.packed_data[:, 1].sum(), -698.4)
        
    def test_repack(self):
        self.setup_data()
        self.assertTrue(N.ma.allclose(self.d.repack(self.data), self.d.packed_data))
        
    def test_repack_notime(self):
        self.setup_data()
        self.assertTrue(N.ma.allclose(self.d.repack(self.data[1]), self.d.packed_data[:, 1]))
        
    def test_weights(self):
        self.setup_data()
        self.assertAlmostEqual(self.d.packed_weights.sum(), 69.5)
    
    def test_unpack(self):
        self.setup_data()
        self.assertTrue(N.ma.allclose(self.d.unpack(self.d.packed_data), self.data))
        
    def test_unpack_notime(self):
        self.setup_data()
        self.assertTrue(N.ma.allclose(self.d.unpack(self.d.packed_data[:, 1]), self.data[1]))
        
         
#    def test_invalid(self):
        

if __name__ == '__main__':
    unittest.main()
