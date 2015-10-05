import unittest
import numpy as N
import os, sys
sys.path.insert(0, '../lib')

from spanlib.data import Dataset


class TestSequenceFunctions(unittest.TestCase):

    def setup_data(self):
        self.d = Dataset((self.data1, self.data2, self.data3),
            weights=(self.weights1, ))

    def setUp(self):
        # data1
        nt, ny, nx = 100, 5, 8
        self.shape = nt, ny, nx
        self.data1 = N.ma.arange(nt*ny*nx, dtype='d').reshape(self.shape)
        self.data1[:, 4, 6:] = N.ma.masked # test fixed mask
        self.weights1 = N.arange(ny*nx).reshape((ny, nx))/10.
        self.data1 -= self.data1.mean(axis=0)
        self.data1 /= self.data1.std()
        # data2
        nx = 70
        self.shape = nt, nx
        self.data2 = N.ma.arange(nt*nx, dtype='d').reshape(self.shape)
        self.data2[:, 0] = N.ma.masked
        self.data2 -= self.data2.mean(axis=0)
        self.data2 /= self.data2.std()
        # data3
        self.shape = nt,
        self.data3 = N.arange(nt, dtype='d')
        self.data3 -= self.data3.mean(axis=0)
        self.data3 /= self.data3.std()

    def test_init(self):
        self.setup_data()

    def test_stack(self):
        self.setup_data()
        self.data1[:, 0, 0] = N.ma.masked
        sum1 = self.data1[1].sum()
        sum2 = self.data2[1].sum()
        sum3 = self.data3[1]
        self.assertAlmostEqual(self.d.stacked_data[:, 1].sum(), sum1+sum2+sum3)

    def test_restack(self):
        self.setup_data()
        self.assertTrue(N.allclose(
            self.d.restack((self.data1, self.data2, self.data3)),
            self.d.stacked_data))

    def test_restack_notime(self):
        self.setup_data()
        self.assertTrue(N.allclose(
            self.d.restack((self.data1[0], self.data2[0], self.data3[0:1])),
            self.d.stacked_data[:, 0]))

    def test_unstack(self):
        self.setup_data()
        data = self.d.unstack(self.d.stacked_data)
        self.assertTrue(N.ma.allclose(self.data1, data[0]))
        self.assertTrue(N.ma.allclose(self.data2, data[1]))
        self.assertTrue(N.ma.allclose(self.data3, data[2]))

    def test_unstack_notime(self):
        self.setup_data()
        data = self.d.unstack(self.d.stacked_data[:, 1])
        self.assertTrue(N.ma.allclose(self.data1[1], data[0]))
        self.assertTrue(N.ma.allclose(self.data2[1], data[1]))
        self.assertTrue(N.ma.allclose(self.data3[1], data[2]))

    def test_toto(self):
        dataset = Dataset((self.data1, self.data1))
        print 1

if __name__ == '__main__':
    unittest.main()
