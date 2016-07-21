import os
import sys
import numpy as npy

TEST_DIR = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(TEST_DIR, '..', 'lib'))


npy.random.seed(0)
