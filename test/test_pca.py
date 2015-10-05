import unittest
import numpy as npy
import os, sys
sys.path.insert(0, '../lib')
from spanlib.analyzer import Analyzer
from spanlib_extra import pca_numpy, setup_data1, setup_data2


class TSF(unittest.TestCase):


    def setup_Analyzer3(self):
        var1 = setup_data1()
        var2 = setup_data1(nx=45, tfreq=4, xfreq=5)
        var3 = setup_data1(nt=120, nx=20)
        return Analyzer([[var1, var2], var3])


    def test_pca(self):

        # spanlib
        A = Analyzer(setup_data1(nt=70, nx=50))
        sp_eof = A.pca_eof(raw=True)
        sp_pc = A.pca_pc(raw=True)
        sp_ev = A.pca_ev(raw=True)

        # numpy
        np_eof, np_pc, np_ev = pca_numpy(A.stacked_data.T, A.npca)

        # checks
        self.assertTrue(npy.allclose(sp_ev, np_ev))
        signs = npy.sign(sp_eof[0])*npy.sign(np_eof[0])
        np_eof *= npy.resize(signs, np_eof.shape)
        np_pc *= npy.resize(signs, np_pc.shape)
        self.assertTrue(npy.allclose(sp_eof, np_eof))
        self.assertTrue(npy.allclose(sp_pc, np_pc))
        self.assertTrue(npy.allclose(sp_ev, np_ev))

    def test_pca_ev(self):
        A = Analyzer(setup_data1(nt=70, nx=150))
        ev = A.pca_ev()
        self.assertTrue(npy.allclose(ev, npy.array(
            [ 24.98560392,  18.04025736,   3.59796042,   1.44110266,
            0.88847967,   0.45467402,   0.24646941,   0.22392288,
            0.16899249,   0.14435935]
        )))

    def test_pca_eof(self):
        A = Analyzer(setup_data1(nt=70, nx=150))
        raw_eof = A.pca_eof(raw=True)
        np_eof, np_pc, np_ev = pca_numpy(A.stacked_data.T, A.npca)
        self.assertTrue(npy.allclose(npy.abs(raw_eof[:2,:2]), npy.abs(npy.array(
            [[ 0.15702435,  0.05325682],
            [ 0.14647224,  0.05082613]]))))
        eof = A.pca_eof()
        self.assertTrue(npy.allclose(npy.abs(eof[:2,:3].filled(0.)),  npy.abs(npy.array(
            [[ 0.15702435,  0., 0.14647224],
            [ 0.05325682,  0., 0.05082613]]))))

    def test_pca_pc(self):
        A = Analyzer(setup_data1(nt=70, nx=150))
        pc = A.pca_pc()
        self.assertTrue(npy.allclose(pc[:2,:2],  npy.array(
            [[ 2.75532971,   2.96740389],
            [-8.78119373,  -8.49616261]])))
#            [[  3.24705271,   3.49697418],
#            [-10.34830743, -10.01240894]])))

    def test_pca_ec(self):
        data = setup_data1(nt=70, nx=150)
        A = Analyzer(data)
        A.pca_eof()
        ec = A.pca_ec()
        pc = A.pca_pc()
        self.assertTrue(npy.allclose(ec, pc))
        ece = A.pca_ec(xeof=A.pca_eof(raw=True), xraw=True)
        self.assertTrue(npy.allclose(ece, pc))
        ecd = A.pca_ec(xdata=data)
        self.assertTrue(npy.allclose(ecd, pc))
        ecrd = A.pca_ec(xdata=A.stacked_data, xraw=True)
        self.assertTrue(npy.allclose(ecrd, pc))
        ecre = A.pca_ec(xeof=A.pca_eof(raw=True), xraw=True)
        self.assertTrue(npy.allclose(ecrd, pc))


    def test_pca_rec(self):
        A = Analyzer(setup_data1())
        rec1 = A.pca_rec()
        self.assertAlmostEqual(rec1.sum(), 18464.31988812385)

    def test_pca_xrec(self):
        var = setup_data1()
        A = Analyzer(var)
        rec = A.pca_rec()
        xeof = A.pca_eof()
        xpc = A.pca_pc()
        xrec = A.pca_rec(xeof=xeof, xpc=xpc)
        self.assertTrue(npy.allclose(rec, xrec))

    def test_pca_ndim3(self):
        A = Analyzer(setup_data2())
        A.pca()
        self.assertTrue(npy.allclose(
            npy.abs(A.pca_eof(raw=True)[190:193,5]),
            npy.abs(npy.array([-0.03623882, -0.02312414, -0.00942417]))))


if __name__ == '__main__':
    unittest.main()
