import os
import configparser
import getpass

thispath=os.path.dirname(os.path.abspath(__file__))+'/'

class configuration(object):
    def __init__(self):
        self.config_ = configparser.ConfigParser()
        self.config_.read_file(open(thispath+'config.cfg'))
    def __getitem__(self, l):
        return self.config_['general'][l]
    def __setitem__(self, l, v):
        self.config_['general'][l]=str(v)
