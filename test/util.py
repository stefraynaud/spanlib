from ConfigParser import SafeConfigParser
import os, sys
import numpy as npy
config_file = os.path.join(os.path.dirname(__file__), '..', 'config.cfg')
def insert_local_path():
    if os.path.exists(config_file):
        cfg = SafeConfigParser()
        cfg.read(config_file)
        if cfg.has_option('paths', 'build_lib'):
            sys.path.insert(0, cfg.get('paths', 'build_lib'))
            print 'insert', cfg.get('paths', 'build_lib')
            print sys.path[0]

