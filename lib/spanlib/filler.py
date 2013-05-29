#################################################################################
# File: spanlib_python.py
#
# This file is part of the SpanLib library.
# Copyright (C) 2006-2013  Stephane Raynaud
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

import numpy as npy
from spanlib.util import Logger, broadcast, SpanlibIter, dict_filter
from spanlib_python import Analyzer, default_missing_value
import pylab as P

class Filler(Logger):
    """Class to fill missing value with a MSSA filtered version of the data
    
    The initialization automatically call the :meth:`fill` method.
    """
    
    def __init__(self, data, run=True, logger=None, loglevel=None, **kwargs):
                
        # Logger
        Logger.__init__(self, logger=logger, loglevel=loglevel, **dict_filter(kwargs, 'log_'))
        self._kwargs = kwargs
        self._kwargs['keep_invalids'] = True
        self._kwargs.setdefault('nvalid', 1)
        self._kwargs.setdefault('quiet', True)
        self._data = data
        self.filled = None
        self.filtered = None
        self.nstep = 0
        self.cv = 100
        
        
        # Setup analyzer
        span = self.analyze(**kwargs)
##            kwargs['keep_invalids'] = False
##            kwargs['nvalid'] = None
#        kwargs['keep_invalids'] = True
#        kwargs.setdefault('nvalid', 1)
#        kwargs.setdefault('quiet', True)
#        print '+span'
#        print 'len(data)', len(data)
#        print 'args', npcamax, nmssamax, kwargs
#        kwargs['npca'] = npcamax
#        kwargs['nmssa'] = nmssamax
#        span = self.span = SpAn(data, sequential=False, logger=self.logger, **kwargs)
##        span.update_params()
##        if mssa: span.update_params(anatype='mssa')
#        print '-span'
        
        # Checks
#        if not span.tree_type<2:
#            self.error('Input data must be a single variable or a list of variables')
        for invalid in span.invalids:
            if invalid is not None and invalid.any():
                break
        else:
            self.warning("No gap to fill")
            
        # Keep original data safe
        self.set_field(span.stacked_data, 'orig')
        
        # Start filling?
        if run: self.fill(**kwargs)
        
    def analyze(self, data=None, npcamax=20, nmssamax=15, **kwargs):
        # Setup analyzer
#            kwargs['keep_invalids'] = False
#            kwargs['nvalid'] = None
        if data is None: data = self._data
        kwargs['keep_invalids'] = True
        kwargs.setdefault('nvalid', 1)
        kwargs.setdefault('quiet', True)
        print '+span'
        print 'len(data)', len(data)
        print 'args', npcamax, nmssamax, kwargs
        kwargs['npca'] = npcamax
        kwargs['nmssa'] = nmssamax
        span = self.span = Analyzer(data, sequential=False, logger=self.logger, **kwargs)
        span.update_params()
#        if mssa: span.update_params(anatype='mssa')
        print '-span'
        return span
        

    def fill(self, nitermax=20, errchmax=-0.01, fillmode='normal', testmode='crossvalid', 
        mssa=True, full=True, cvregen=False, skipat=2, nreana=3, **kwargs):
        """Run the filler with a convergence loop
        
        Results are accessible in the following attributes:
        
        .. attribute:: filtered
        
            Filtered data, result from the convergence loop.
            
        .. attribute:: filled
        
            Data filled with :attr:`filtered`
            
        :Parameters:
        
            - **fillmode**: "zeros" or "normal"
            - **nitermax**: Maximal number of iterations
            - **cvmax**: Convergence criterion (%)
            - **npca**: Number of PCA modes (see :class:`Analyzer`)
            - **nmssa**: Number of MSSA modes (see :class:`Analyzer`)
            - Other parameters are passed to :class:`Analyzer`
            
        :Returns:
        
            - :attr:`filled`
        
        """
    
        
        # Parameters
        kwgencv = dict_filter(kwargs, 'cvfield_')
        span = self.span
        if fillmode==0:
            fillmode = "none"
        fillmode = str(fillmode).lower()
        kwargs['zerofill'] = 0 if fillmode.startswith('n') else 2
        self.debug('Filling mode: %s (%i)'%(fillmode, kwargs['zerofill']))
        span.update_params(**kwargs)
        
        # which analyzes types?
        analyzes = []
        if not mssa or span.prepca: analyzes.append('pca')
        if mssa: analyzes.append('mssa')
        
        import pylab as P
        # Which modes?
        testmode = str(testmode).lower()
        if testmode.startswith('c'):
            kwgencv['regen'] = cvregen
            self.gen_cvfield(**kwgencv)
            anamodes = [2, 0]
        else:
            anamodes = [1]
            
        # Loop on analysis modes
        print 'anamodes', anamodes
        modes = {}
        for anamode in anamodes:
            self.debug('Analysis mode: '+['NORMAL', 'SELF-VALIDATION', 'CROSS-VALIDATION'][anamode])
            
            # Initialization according to analysis mode
            if anamode!=2:
                self.load_field('orig')
            if anamode!=0:
                if anamode==2:
                    self.load_field('cvfield')
#                modes = {}
#                if 'pca' in analyzes: modes['pca'] = range(span.npca)
#                if 'mssa' in analyzes: modes['mssa'] = range(span.nmssa)
        
            # Loop on analysis types (PCA, MSSA)
            print '---------------------\n modes', modes
            for self._ana in analyzes:
                self.debug('Running '+self._ana.upper())
                
                # Update the number of pre-PCA modes for MSSA
                if self._ana=='mssa' and span.prepca:
                    span.update_params(prepca = imode+1)
                    
                # Reanalyses loop
                for irf in xrange(nreana):
                    
                    self.debug('Analysis (%i/%i)'%(i+1, nrena))                        
                
                    # Run analysis to get EOFs
                    self.get_func()(force=True)
                    
                    # CV loop on EC estimation (not for PCA with T-EOF)?
                    ecloop = self._ana!= 'pca' or not self.span.useteof
                    niterec = nitermax if ecloop else 1
                    
                    # Number of modes to retain
                    modes.setdefault(self._ana, [])
                    amodes = modes[self._ana]
                    if len(amodes)<irf+1:
                        amodes.append(range(getattr(self.span, 'n'+self._ana))) # test all
                    
                    # Loop on the number of modes 
                    last_mode_err = None
                    self._last_pcs_mode = {}
                    for imode in amodes[irf]:
                        verb = 'Reconstructing' if anamode==0 else 'Trying'
                        self.debug(verb+' with %i mode%s'%(imode+1, 's'*(imode>0)))
                        
                        # Inits
                        self.iter_init()
                        last_iter_err = None
                        skiplast = False
                        if anamode==1 and imode>0:
                            self._last_pcs_mode[self._ana] = \
                                getattr(self.span, '_%s_raw_pc'%self._ana)
                        self._last_pcs_iter = {}
                        
                        # Convergence loop for expansion coefficients
                        for istep in xrange(niterec):
                            if ecloop: self.debug(' EC convergence step: %i'%istep)
                            
                            # Reconstruct
                            if anamode==1 and istep>0:
                                self._last_pcs_iter[self._ana] = \
                                    getattr(self.span, '_%s_raw_pc'%self._ana)
                            self.rec(imode)
        
                            # Current error
                            err = self.get_error(anamode)
                            
                            # Check MSSA full filling
                            if self._ana=='mssa' and full:
                                nem, nemch = self.get_necmiss()
                                if nem and (nemch or nemch is None):
                                    self.debug('  Still %i missing values in %s PCs'%
                                        (nem, self._ana.upper()))
                                    last_iter_err = err
                                    continue
        
                            # Check convergence error for EC
                            if ecloop:
                                self.debug('  Current error: %.1f%%'%err)
                                if istep>0 and last_iter_err is not None:
                                    errch = err - last_iter_err
                                    self.debug('  Error change: %g%%'%errch)
                                    if errch>=errchmax:
                                        if errch>0:
                                            self.debug('  Error change > 0: unstable mode -> step skipped')
                                            err = last_iter_err
                                            if anamode==1:
                                                setattr(self.span, '_%s_raw_pc'%self._ana, 
                                                    self._last_pcs_iter[self._ana])
                                            self.debug('Recovered PCs from last step')
        #                                    skiplast = True
                                        else:
                                            self.debug('  Error change > threshold -> stopping convergence'%errch)
                                        break
                            else:
                                self.debug('  Error: %.1f%%'%err)
                                break
                            last_iter_err = err
                        
                        else:
                            self.debug('  Reached max number of iterations for EC convergence loop')
                                        
                        # Check mode truncature error
                        if anamode!=0 and imode>0:
                            if skiplast: 
                                errch = 1
                            else:
                                errch = err-last_mode_err
                                self.debug(' Error change between %i and %i modes: %g'%
                                    (imode, imode+1, errch))
                            if errch>errchmax:
                                imode -= 1
                                self.debug(' Best number of %s modes: %i'%(self._ana.upper(), imode+1))
    #                            getattr(self.span, '_%s_raw_pc'%self._ana)[0] = self._last_pcs
                                if anamode==1:
                                    setattr(self.span, '_%s_raw_pc'%self._ana, 
                                        self._last_pcs_mode[self._ana])
                                    self.debug('Recovered PCs from last mode')
                                break    
                        last_mode_err = err
                    else:
                        self.debug(' Reached max number of %s modes (%i)'%(self._ana.upper(), imode+1))
                 
                    # Refill
                    if nreana>0:
                        self.span.stacked_data = npy.where(self._orig_mask, 
                            self._rec.filled(default_missing_value), self._orig)
       
                    # Store optimal number of modes for normal analysis after cross-validation
                    if anamode==2:
                        modes[self._ana][irf] = [imode]
#                        print 'store modes[%s][%i] = %s'%(self._ana, irf, modes[self._ana])
                        
        # Final reconstruction
        rec_meth = getattr(self.span, self._ana+'_rec')
        self.filtered = rec_meth(modes=-imode)

    def get_error(self, anamode):
        diff = self._origm-self._rec
        field = self._origm
        if anamode==2:
            diff = npy.ma.masked_where(self._cvfield_kept, diff)
            field = npy.ma.masked_where(self._cvfield_kept, field)
            goodcv = ~diff.mask
        print 'get_error: anamode', anamode, npy.ma.count(diff)*1./diff.size
        return 100*diff.std()/field.std()
        
     
    def iter_init(self):
        self._rec = default_missing_value
        if hasattr(self, '_necmiss'): del self._necmiss
        if hasattr(self, '_necmiss_old'): del self._necmiss_old
        
    def rec(self, imode):
        # Merge
        data = npy.where(self._current_mask, self._rec, self._current)
        
        # Expansion coefficients
        self._last_pcs = {} #getattr(self.span, '_%s_raw_pc'%self._ana)[0]
        if (self._ana=='pca' or self.span.prepca) and not self.span.useteof:
            xxx
            self.span.pca_ec(raw=True, xraw=True, xdata=data, replace=True, demean=False)
        if self._ana=='mssa':
            kw = dict(xdata=data, xraw=True) if not self.span.prepca else {}
            self.span.mssa_ec(raw=True, replace=True, demean=False, **kw)
        del data
    
        # Reconstruction
        recfunc = self.get_func('rec')
        self._rec = recfunc(modes=-imode, raw=2, rescale=False, unmap=False)
        if hasattr(self, '_current_mean'):
            mean = npy.repeat(self._current_mean, self._rec.shape[1], axis=1)
#            mean = mean.filled(default_missing_value)
            self._rec += mean
            del mean
        
#        if self._ana=="mssa":
#            P.plot(self._currentm[0], 'g')
#            P.plot(self._rec[0], 'r')
#            P.show()
    
    def get_func(self, suf=None):
        ana = self._ana
        if suf is not None:
            ana += '_'+suf
        return getattr(self.span, ana)
        
    def get_necmiss(self):
        """Get the number of missing values of current expansion coefficents"""
        if not hasattr(self,  '_ana'): return None, None
        if hasattr(self, '_necmiss'): self._necmiss_old = self._necmiss
        pc0 = npy.ma.masked_values(getattr(self.span, '_%s_raw_pc'%self._ana)[:, 0], 
            default_missing_value)
        self._necmiss =  pc0.size-npy.ma.count(pc0)
        self._necmiss =  pc0.size-npy.ma.count(pc0)
        del pc0
        if not hasattr(self, '_necmiss_old'):  
            necmissch = None
        else:
            necmissch = self._necmiss-self._necmiss_old
        return self._necmiss, necmissch

    def gen_cvfield(self, mask=None, level=15, sscale=1, tscale=1, minrelscale=0.2, regen=False):
        """Generate false missing values
        
        Field used for analysis is stored using `set_field('cvfield')`.
        Retained points are saved as logical array in :attr:`_cvfield_kept`.
        
        :Params:
        
            - **mask**, optional: Mask set to True where values are not retained for
              analysis. These rejected values are used to cross-validation.
              If not set, it is randomly generated.
        """
        if not regen and hasattr(self, '_cvfield'): return

        # Generate the mask
        if mask is None:
            
            # Smooth scale
            ns, nt = self._orig.shape
            if sscale<0:
                sscale = -int(ns*sscale/100.)
            if sscale > minrelscale*ns:
                sscale = 1
            else:
                sscale = sscale/2*2+1
            if tscale<0:
                tscale = -int(nt*tscale/100.)
            if tscale > minrelscale*nt:
                tscale = 1
            else:
                tscale = tscale/2*2+1
            
            # Random values
            tmp = npy.random.random((ns, nt)).astype('f')
            
            # Apply smoothing
            if tscale>1 or sscale>1:
                ww = npy.ones((ns, nt), 'f')
                from numpy import convolve
                kernel = npy.ones(sscale, 'f')
                if sscale>1:
                    for i in xrange(nt):
                        tmp[:, i] = convolve(tmp[:, i], kernel, 'same') / convolve(ww[:, i], kernel, 'same')
                kernel = npy.ones(tscale, 'f')
                if tscale>1:
                    for i in xrange(ns):
                        tmp[i] = convolve(tmp[i], kernel, 'same') / convolve(ww[i], kernel, 'same')
                del ww
            
            # Get max value for noisy field
            sorted = npy.sort(tmp.ravel())
            nlevel = int(npy.round(tmp.size*level/100.))
            maxval = sorted[max(0, nlevel-1)]
            del sorted
            
            # Mask
            mask = tmp<maxval
            del tmp
        
        # Apply mask
        cvfield = npy.where(mask, default_missing_value, self._orig)
        self.set_field(cvfield, 'cvfield', mean=True)
        self._cvfield_kept = ~mask
        
    def set_field(self, field, name, mean=False):
        setattr(self, '_%s'%name, field)
        fieldm = npy.ma.masked_values(field, default_missing_value)
        setattr(self, '_%sm'%name, fieldm)
        setattr(self, '_%s_mask'%name, npy.ma.getmaskarray(fieldm))
        setattr(self, '_%s_std'%name, fieldm.std())
        if mean:
            setattr(self, '_%s_mean'%name, fieldm.mean(axis=1).reshape(field.shape[0], -1))

    def load_field(self, name):
        """Load a stored field to the current field into 
        self.span.stacked_data and self._current
        """
        for att in '', 'm', '_mask', '_std', '_mean':
            if hasattr(self, '_%s'%name+att):
                setattr(self, '_current'+att, getattr(self, '_%s'%name+att))
        self.span.stacked_data = getattr(self, '_current')
    
        
        
#                    # check convergence
#                    self.cv = (energy-last_energy)/energy
#                    if npy.isnan(self.cv): xxxxx
#    #                self.debug('Change rate: %s%%'%(self.cv*100))
#                    print '>>>>>>>>>>> Change rate: %s%%'%(self.cv*100)
#                    print '1'
#                    if self.cv < 0 or npy.abs(self.cv) < cvmax/100.:
#                        self.info('Convergence reached: %.2f%% (%i iterations)'%(100*self.cv, self.nstep))
#                        break
#                    last_energy = energy
#                
#                # Reconstruction
#                print '+fillerrec'
#    #            drec = span.mssa_rec()
#    #            print '>'*10, 'pca_ec'
#    #            ec = span.pca_ec(xdata=data, replace=True)
#                print '>'*10, 'mssa_ec'
#                ecm = span.mssa_ec(xdata=data, replace=True)
#                print '>'*10, 'mssa_rec'
#                datarec = span.mssa_rec()
#    #            datarec = span.pca_rec()
#    #            print '-fillerrec', datarec.sum()
#    #            import pylab as P
#    ##            P.plot(span.mssa_pc())
#    #            P.plot(datarec)
#    #            P.show()
#    #            xxxxx
#                mdatarec = span.remap(datarec)
#    #            print '-rec'
#    
#                # Replace invalid data with reconstructed field
#    #            print '+replace'
#                for ivar in xrange(len(mdatarec)):
#                    if invalids[ivar] is not None:
#                        if callable(filter):
#                            mdatarec[ivar][:] = filter(mdatarec[ivar])
#                        mdata[ivar][:] = eval(span[ivar].array_type).where(invalids[ivar], 
#                            mdatarec[ivar], mdata[ivar])
#    #            print '-replace'
#    
#                # Save current filtered state
#                if getfiltsteps: self.filtsteps.append(datarec)
#                            
#    #            print '+unmap'
#                data = span.unmap(mdata)
#    #            data = span._unmap_([mdata])
#    #            old_span = span                    
#    #            print '-unmap'
#        
#                # Check iteration limit
#                if self.nstep >= nitermax:
#    #                self.warning('Convergence not reached %i%% (%i iterations)'%(100*self.cv, self.nstep))
#                    break
#                
#        # Output
#        span = old_span if old_span is not None else span
#        self.span = span
#        self.nmssa = span.nmssa()
#        self.npca = span.npca()
#        self.filled = data
#        self.filtered = datarec
##        print '+ids'
#        if span[0].has_cdat(): # Ids for CDAT variables
#            if isinstance(self.filled, list):
#                for i, data in enumerate(span[0]):
#                    if not data.has_cdat(): continue
#                    self.filled[i].id +='_filled'
#                    self.filtered[i].id += '_filtered'
#            else:
#                self.filled.id += '_filled'
#                self.filled.id += '_filtered'
##        print '-ids'
#        return data

      
