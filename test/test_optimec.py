
import unittest
import numpy as npy
import os, sys
sys.path.insert(0, '../lib')
from spanlib_extra import pca_numpy, setup_data1, setup_data2
from spanlib.analyzer import Analyzer, default_missing_value
from spanlib._optimec import optimec_pca



class TSF(unittest.TestCase):


    def test_fortran_pca_filled(self):

        # Direct PC
        A = Analyzer(setup_data1(nt=30, nx=100))
        direct_pc = A.pca_pc(raw=True)
        eof = A.pca_eof(raw=True)

        # PC from optimec
        optim_pc = npy.asfortranarray(direct_pc.copy()*2+npy.random.normal(0,1.,direct_pc.shape))
        optimec_pca(A.stacked_data, eof, optim_pc, default_missing_value)

        # Check
        self.assertTrue(npy.allclose(optim_pc, direct_pc))

if __name__ == '__main__':
    unittest.main()

