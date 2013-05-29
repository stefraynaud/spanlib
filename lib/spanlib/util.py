#################################################################################
# File: util.py
#
# This file is part of the SpanLib library.
# Copyright (C) 2006-2011  Stephane Raynaud
# Contact: stephane dot raynaud at gmail dot com
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#################################################################################

import logging
import logging.handlers
import os
import numpy as N
import re

re_errmsg_split = re.compile(r'[\[\]]').split

def broadcast(input, mylen, fillvalue=None):
    """Make sure that input has the right length"""
    # A single value
    if mylen == 0:
        if isinstance(input,(list,tuple)):
            if not input: return None
            return input[0]
        return input
        
    # Multiple values as a list (or tuple)
    if not isinstance(input,(list,tuple)):
        fillvalue = input
        input = [input]
    if isinstance(input,tuple):
        input = list(input)
    dlen = mylen-len(input)
    if dlen < 0:
        input = input[:mylen]
    elif dlen > 0:
        input.extend([fillvalue]*dlen)
    return input

class SpanlibError(Exception):
    """Reporter of exceptions
    
    :Params:
    
        - **what**: What is the error.
        - **where**, optional: Where the error occured.
        
    :Example:
    
        >>> SpanlibError('Bad number of channels', 'pca')
    """
    def __init__(self, what, where=None):
        Exception.__init__(self)
        self._where = where
        self._what = what
    def __str__(self):
        if self._where is None:
            return 'SpanlibError: %s' %  self._what
        return 'SpanlibError: [%s] %s' % (self._where, self._what)

class SpanlibIter:
    """Iterator over datasets"""
    def __init__(self, span):
        self.span = span
        self.iset = 0
        span._iter = self
    def __iter__(self):
        return self
    def next(self):
        if self.iset<self.span.nd:
            self.iset += 1
            return self.span[self.iset-1]
        del self.span._iter
        raise StopIteration
 

class Logger(object):
    """Class for logging facilities when subclassing.
    Logging may be sent to the console and/or a log file
    
    :Params:
    
        - **name**: Name of the logger.
        - **logfile**, optional: Log file.
        - **console**, optional: Log to the console.
        - **maxlogsize**, optional: Maximal size of log file before rotating it.
        - **maxbackup**, optional: Maximal number of rotated files.
        - **sfmt**, optional: Format of log messages in log file.
        - **cfmt**, optional: Format of log message in console.
        - **asctime**, optional: Time format.
        - **level**, optional: Initialize logging level (see :meth:`set_loglevel`).
    
    :See also: :mod:`logging` module
    
    .. note:: Inspired from :class:`vacumm.misc.io.Logger`
    """
    def __init__(self, name=None, logfile=None, console=True, maxlogsize=0, maxbackup=0, 
        cfmt='%(name)s [%(levelname)-8s] %(message)s', 
        ffmt='%(asctime)s: %(name)s [%(levelname)-8s] %(message)s', 
        asctime='%Y-%m-%d %H:%M', 
        loglevel='warning', logger=None, 
        parser=None, announce=False):
        
        # Setup the logger
        if logger is None:
            # Handlers
            # - file
            if logfile is not None and logfile != '':
                logdir = os.path.dirname(logfile)
                if logdir != '' and not os.path.exists(logdir): os.makedirs(logdir)
                file =  logging.handlers.RotatingFileHandler(logfile, 
                    maxBytes=maxlogsize*1000, backupCount=maxbackup)
                file.setFormatter(logging.Formatter(ffmt, asctime))
            # - console
            if console:
                console = logging.StreamHandler()
                console.setFormatter(logging.Formatter(cfmt))
        
            # Logger
            if name is None: name = self.__class__.__name__.upper()
            logger = logging.getLogger(name)
            if logfile is not None: logger.addHandler(file)
            if console is not None: logger.addHandler(console)
        self.logger = logger
        self.parser = parser
            
        # Set level
        self.set_loglevel(loglevel)
        
        # Announce
        if announce: logger.debug('*** Start log session ***')
        
    def debug(self, text, *args, **kwargs):
        """Send a debug message"""
        self.logger.debug(text, *args, **kwargs)
        
    def info(self, text, *args, **kwargs):
        """Send a info message"""
        self.logger.info(text, *args, **kwargs)
        
    def warning(self, text, *args, **kwargs):
        """Send a warning message"""
        self.logger.warning(text, *args, **kwargs)
        
    def error(self, text, *args, **kwargs):
        """Send an error message and raise :class:`SpanLibError`
        
        :Params:
        
            - **text**: Text to display.
            - **errfunc**, optional: Callable function to display error once logged.
              If ``errfunc`` is not passed and :attr:`parser` has been passed 
              during initialization, :meth:`argparse.ArgumentParser.error` method is used.
            - **errxargs**: List of extra arguments for ``errfunc``.
            
        :Raise:
        
            If ``errfunc`` is not provided or does not exit, :class:`SpanlibError` 
            is raised.
            
        :Examples:
        
            >>> logger.Logger('mylogger')
            >>> logger.error('my error')
            
            >>> myparser = ArgumentParser()
            >>> logger.Logger('mylogger', parser=myparser)
            >>> logger.error('my error') # implicit call to myparser.error
            
            >>> myparser = ArgumentParser()
            >>> logger.Logger('mylogger')
            >>> logger.error('--myoption', errfunc=parser.ArgumentError,
            ... errxargs='error calling this option')

        """
        self.logger.error(text, *args, **kwargs)
        errfunc = kwargs.pop('errfunc', None)
        errargs = kwargs.pop('errxargs', [])
        if not callable(errfunc) and hasattr(self.parser, 'error'):
           errfunc = self.parser.error 
           errxargs = []
        if callable(errfunc):
            if not isinstance(errxargs, (list, tuple)):
                errargs = [errxargs]
            errargs = [text]+list(errxargs)
            errfunc(*errargs)
        raise SpanlibError(text)
        
    def set_loglevel(self, level=None, console=None, file=None):
        """Set the log level (DEBUG, INFO, WARNING, ERROR)
        
        :Example:
        
            >>> logger.set_loglevel('DEBUG', console='INFO')
        """
        if level is not None:
            self.logger.setLevel(self._get_loglevel_(level))
        for handler in self.logger.handlers:
            if isinstance(handler, logging.handlers.RotatingFileHandler):
                if file is not None: handler.setLevel(self._get_loglevel_(file))
            elif console is not None and \
                isinstance(handler, logging.StreamHandler):
                handler.setLevel(self._get_loglevel_(console))
                
    def get_loglevel(self, asstring=False):
        """Get the log level as an integer or a string"""
        if asstring:
            for label in 'NOTSET', 'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL', 'FATAL':
                if self.logger.level == getattr(logging, label):
                    return label
            return 'NOTSET'
        return self.logger.level
        
    def _get_loglevel_(self, level):
        if level is None: level = 'debug'
        if isinstance(level, str): 
            level = getattr(logging, level.upper(), 'DEBUG')
        return level

    def check_fortran_errmsg(self, errmsg):
        errmsg = errmsg.strip()
        if not errmsg: return
        _, loglevel, errmsg = re_errmsg_split(errmsg)
        errmsg = 'FORTRAN - '+errmsg.strip()
        getattr(self, loglevel)(errmsg)


def dict_filter(kwargs, prefix, pop=True):
    """Filter a dictionary by keeping only items starting with a given prefix
    
    :Example:
    
        >>> kwargs = dict(a=3, logfile='mylog.log', loglevel=0) 
        >>> dict_filter(kwargs, 'log')
        {'file':'mylog.log','level':0}
        >>> kwargs
        {'a': 3}
    
    """
    kwfilt = {}
    for key, val in kwargs.items():
        if key.startswith(prefix):
            kwfilt[key[len(prefix):]] = val
            if pop: kwargs.pop(key)
    return kwfilt
    


