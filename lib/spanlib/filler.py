#################################################################################
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
from util import Logger, broadcast, SpanlibIter, dict_filter, SpanlibError
from analyzer import Analyzer, default_missing_value
#import pylab as P

class FillerError(SpanlibError):
    pass

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
#        self.filled = None
#        self.filtered = None
        self.nstep = 0
        self.cv = 100
        self._analyzes = []   
        self._kwfill = {}
        self._ana = 'pca'
        
        # Setup analyzer
        span = self.analyze(**kwargs)
        
        # Checks
        if not span.invalids.any():
            self.warning("No gap to fill")
            
        # Keep original data safe
        self._set_field_(span.stacked_data, 'orig')
        
        # Start filling?
        if run: self.fill(**kwargs)
        
    def analyze(self, data=None, npcamax=15, nmssamax=15, **kwargs):
        # Setup analyzer
        if data is None: data = self._data
        kwargs['keep_invalids'] = True
        kwargs.setdefault('nvalid', 1)
        kwargs.setdefault('quiet', True)
        kwargs.setdefault('npca', kwargs.pop('npca', npcamax))
        kwargs.setdefault('nmssa', kwargs.pop('nmssa', nmssamax))
        kwargs.setdefault('window', kwargs.pop('win', None))
        span = self.span = Analyzer(data, logger=self.logger, **kwargs)
        span.update_params()
        self.mask = span.invalids
        return span
        

    def fill(self, nitermax=20, errchmax=-0.01, fillmode='masked', testmode='crossvalid', 
        mssa=True, full=True, cvregen=False, nreanapca=3, nreanamssa=2, errchmaxreana=-1, 
        remode=False, **kwargs):
        """Run the filler with a convergence loop
        
        Results are accessible in the following attributes:
        
        .. attribute:: filtered
        
            Filtered data, result from the convergence loop.
            
        .. attribute:: filled
        
            Data filled with :attr:`filtered`
            
        :Parameters:
        
            - **fillmode**: "zeros" or "masked"
            - **nitermax**: Maximal number of iterations
            - **cvmax**: Convergence criterion (%)
            - **npca**: Number of PCA modes (see :class:`Analyzer`)
            - **nmssa**: Number of MSSA modes (see :class:`Analyzer`)
            - **cvfield_level**: Percent of data used for cross-validation.
            - Other parameters are passed to :class:`Analyzer`
            
        :Returns:
        
            - :attr:`filled`
        
        """
    
        
        # Parameters
        self._kwfill.update(nitermax=nitermax, errchmax=errchmax, fillmode=fillmode, 
            testmode=testmode, mssa=mssa, full=full, cvregen=cvregen, 
            nreanapca=nreanapca, nreanamssa=nreanamssa, **kwargs)
        kwgencv = dict_filter(kwargs, 'cvfield_')
        span = self.span
        if fillmode==0:
            fillmode = "none"
        fillmode = str(fillmode).lower()
        if fillmode.startswith('n') or fillmode.startswith('m'): 
            fillmode = "masked"
        elif fillmode.startswith('z'):
            fillmode = "zeros"
        kwargs['zerofill'] = 0 if fillmode=="masked" else 2
        kwargs['prepca'] = True # always PCA
        self.debug('Filling mode: %s (%i)'%(fillmode, kwargs['zerofill']))
        span.update_params(**kwargs)
        nreana = dict(pca=nreanapca, mssa=nreanamssa)
        
        # which analyzes types?
        analyzes = []
#        if not mssa or span.prepca: 
        analyzes.append('pca')
        if mssa: analyzes.append('mssa')
        self._analyzes = analyzes
        self._nomssaneeded = False
        
        # Which modes?
        testmode = str(testmode).lower()
        if testmode.startswith('c'):
            kwgencv['regen'] = cvregen
            anamodes = [2, 0]
        else:
            anamodes = [1]
            
        # Loop on analysis modes
        self._nmodes = {}
        self._errors = {}
        for anamode in anamodes:
            self.debug('Analysis mode: '+['NORMAL', 'SELF-VALIDATION', 'CROSS-VALIDATION'][anamode])
            
            # Loop on analysis types (PCA, MSSA)
            self._errors[anamode] = {}
            for self._ana in analyzes:
                self.debug(' Running '+self._ana.upper())
                
                # Update the number of pre-PCA modes for MSSA
                if self._ana=='mssa':
                    span.update_params(prepca = imode+1)
                    
                # Link to appropriate data
                
                # - reference data
                self._link_field_('orig' if self._ana=='pca' else 'pcs', 'ref') 
                rmask = self._refm.mask
                saxis = int(self._ana=='mssa')
                if rmask is npy.ma.nomask or not (rmask.all(axis=saxis)|rmask.any(axis=saxis)).any():
                    self.warning('%s: No gap to fill -> skipping'%self._ana.upper())
                    analyzes.remove('mssa')
                    self._nomssaneeded = True
                    self._nmodes.setdefault(self._ana, [[self.span.nmssa]])
                    break
#                    self.warning('%s: No gap to fill -> just analyzing with all modes'%self._ana.upper())
#                    self._get_func_()() # exec
#                    self._nmodes.setdefault(self._ana, [[self.span.nmssa]])
#                    break
                
                # - data to fill
                if anamode==2: # cross-validation
                    self._gen_cvfield_(**kwgencv)
                    self._link_field_('cvfield', 'current')
                else: # normal, self-validation
                    self._link_field_('ref', 'current')
                    
                # Initialize raw data
                self._set_raw_(self._currentm.data)
                
                # Reanalyses loop
                self._errors[anamode][self._ana] = []
                last_reana_err = None
                for ira in xrange(nreana[self._ana]):
                    
                    self.debug('  Analysis (%i/%i)'%(ira+1, nreana[self._ana]))                        
                
                    # Run analysis to get EOFs
                    self._get_func_()(force=True)
                    
                    # CV loop on EC estimation (not for PCA with T-EOF)?
                    ecloop = self._ana!= 'pca' or not self.span.useteof
                    niterec = nitermax if ecloop else 1
                    
                    # Number of modes to retain
                    self._nmodes.setdefault(self._ana, [])
                    amodes = self._nmodes[self._ana]
                    trymodes = anamode!=0
                    if len(amodes)<ira+1:
                        if remode or ira==0:
                            amodes.append(range(getattr(self.span, 'n'+self._ana))) # test all
                        else:
                            trymodes = False
                            amodes.append(amodes[0]) # same as for first analysis
                    
                    # Loop on the number of modes 
                    last_mode_err = None
                    self._last_pcs_mode = {}
                    self._errors[anamode][self._ana].append([])
                    for im, imode in enumerate(amodes[ira]):
                        verb = '   Reconstructing' if not trymodes else '   Trying'
                        self.debug(verb+' with %i mode%s'%(imode+1, 's'*(imode>0)))
                        
                        # Inits
                        self._recm = default_missing_value
                        if hasattr(self, '_necmiss'): del self._necmiss
                        if hasattr(self, '_necmiss_old'): del self._necmiss_old
                        last_iter_err = None
                        skiplast = False
                        if anamode==1 and im>0: # save previous pcs
                            self._last_pcs_mode[self._ana] = \
                                getattr(self.span, '_%s_raw_pc'%self._ana)
                        self._last_pcs_iter = {}
                        
                        # Convergence loop for expansion coefficients
                        self._errors[anamode][self._ana][-1].append([])
                        for istep in xrange(niterec):
                            if ecloop: self.debug('    EC convergence step: %i'%istep)
                            
                            # Reconstruct
                            if anamode==1 and istep>0:
                                self._last_pcs_iter[self._ana] = \
                                    getattr(self.span, '_%s_raw_pc'%self._ana)
                            self._rec_(imode)
        
                            # Current error
                            err = self._get_error_(anamode)
                            self._errors[anamode][self._ana][-1][-1].append(err)
                            
                            # Check MSSA full filling
                            if self._ana=='mssa' and full:
                                nem, nemch = self._get_necmiss_()
                                if nem and (nemch or nemch is None):
                                    self.debug('    Still %i missing values in %s MSSA PCs'%
                                        (nem, self._ana.upper()))
                                    last_iter_err = err
                                    continue
        
                            # Check convergence error for EC
                            if ecloop:
                                self.debug('    Current error: %.1f%%'%err)
                                if istep>0 and last_iter_err is not None:
                                    errch = err - last_iter_err
                                    self.debug('    Error change: %g%%'%errch)
                                    if errch>=errchmax:
                                        if errch>0:
                                            self.debug('  Error change > 0: unstable mode -> step skipped')
                                            err = last_iter_err
                                            if anamode==1:
                                                setattr(self.span, '_%s_raw_pc'%self._ana, 
                                                    self._last_pcs_iter[self._ana])
                                            self.debug('    Recovered PCs from last step')
                                            self._errors[anamode][self._ana][-1][-1] *= -1
        #                                    skiplast = True
                                        else:
                                            self.debug('    Error change > threshold -> stopping convergence'%errch)
                                        break
                            else:
                                self.debug('    Error: %.1f%%'%err)
                                break
                            last_iter_err = err
                        
                        else:
                            self.debug('  Reached max number of iterations for EC convergence loop')
                                        
                        # Check mode truncature error
                        if anamode!=0 and im:
                            if skiplast: 
                                errch = 1
                            else:
                                errch = err-last_mode_err
                                self.debug('   Error change between %i and %i modes: %g'%
                                    (imode, imode+1, errch))
                            if errch>errchmax:
                                imode -= 1
                                self.debug('   Best number of %s modes: %i'%(self._ana.upper(), imode+1))
                                if anamode==1:
                                    setattr(self.span, '_%s_raw_pc'%self._ana, 
                                        self._last_pcs_mode[self._ana])
                                    self.debug('   Recovered PCs from last mode')
                                    if self._errors[anamode][self._ana][-1][-1] >0:
                                        self._errors[anamode][self._ana][-1][-1] *= -1
                                break    
                        last_mode_err = err
                    else:
                        if trymodes:
                            self.debug('   Reached max number of %s modes (%i)'%(self._ana.upper(), imode+1))
                            
                    # -> NMODES
                 
                    # Refill for next reana
                    if ira<nreana[self._ana]-1:
                        self.debug('  Filling for next reanalysis')
                        self._set_raw_(self._currentm.filled(self._recm.data))
                        
                    # Store optimal number of modes info for normal analysis after cross-validation
                    if anamode==2:
                        self._nmodes[self._ana][ira] = [imode]
                        
                    # Error check
                    if ira>0:
                        
                        errch = err-last_reana_err
                        self.debug('  Error change since last analysis: %.2g'%errch)
                        if errch>errchmaxreana:
                            self.debug('  Stopping reanalyzes')
                            break
                            
                    last_reana_err = err
                
                else:
                    self.debug('  Reached max number of reanalyzes for %s (%i)'%(self._ana.upper(), ira+1))
       
                # -> REANA
                
                # Store number of reanalyzes for normal analysis after cross-validation
                if anamode==2:
                    nreana[self._ana] = ira+1

                
                # Store PCA pcs for MSSA
                if self._ana=='pca' and 'mssa' in analyzes:
                    self.span.prepca = imode+1
                    self._set_field_(self.span._pca_raw_pc[:, :imode+1].T, 'pcs', mean=True, std=False)
                
            # -> PCA/MSSA
            
        # -> NORM/SELF/CV
               
        
    def get_filtered(self, mssa=None, **kwargs):
        """Get a filtered version of the original dataset"""
        if not self._analyzes:
            kw = self._kwfill.copy()
            kw.update(kwargs, mssa=mssa is not False)
            self.fill(**kw)
        if mssa is False:
            ana = 'pca'
        elif mssa is None:
            ana = self._ana
        else:
            ana = 'mssa'
            if 'mssa' not in self._analyzes:
                if self._nomssaneeded: # MSSA not needed to fill but we run it
                    self.span.mssa()
                else: # Refill with MSSA
                    kw = self._kwfill.copy()
                    kw.update(kwargs, mssa=True)
                    self.fill(mssa=True, **kwargs)            
        rec_meth = getattr(self.span, ana+'_rec')
        return rec_meth(modes=-self._nmodes[ana][-1][0])
        
    filtered = property(fget=get_filtered, doc='Filtered data')
        
    def get_filled(self, mode="mssa", **kwargs):
        """Get a version of the original dataset where missing data are
        filled with the filtered version
        
        :Params:
        
            - **mode**: Filling mode.
            
                - ``"best"``: Keep original values where they are present,
                  then fill with PCA reconstruction where possible, 
                  finally fill with MSSA reconstruction.
                - ``"pca[+]"``: Filled only with PCA reconstruction. 
                  If ``mode`` has a ``+``, the filtered version  is returned
                  instead of the filled version.
                - ``"mssa[+]"``: Filled with MSSA reconstruction.
                  If ``mode`` has a ``+``, the filtered version  is returned
                  instead of the filled version.
                - ``"both"`` or ``"best+"``: Use the PCA reconstruction as a basis,
                  and fill it with MSSA reconstruction.
                - ``None`` or ``"auto[+]"``: It depends on current analyzes, and try to use the 
                  ``"mssa"`` mode.
                  If ``mode`` has a ``+``, the filtered version  is returned
                  instead of the filled version.

        """ 
        mode = str(mode).lower()
        
        if mode=='best':
            pcafiltered = self.get_filtered(mssa=False)
            out = self.span.fill_invalids(self._data, pcafiltered, copy=True)
            del pcafiltered
            mssafiltered = self.get_filtered(mssa=True)
            out =  self.span.fill_invalids(out, mssafiltered, copy=False, missing=True)   
            return out
            
        if mode=='both' or ('best' in mode and '+' in mode):
            pcafiltered = self.get_filtered(mssa=False, unmap=False)
            mssafiltered = self.get_filtered(mssa=True, unmap=False)
            out = []
            for vp, vm in zip(pcafiltered, mssafiltered):
                out.append(vp)
                vp[:] = npy.ma.where(vp.mask, vm, vp)
            return self.span.unmap(out)
            
        if 'none' in mode or 'auto' in mode or mode=='+':
            mssa = None
        elif 'pca' in mode:
            mssa = False
        elif 'mssa' not in mode:
            raise FillerError('Unknown filling mode: '+mode)
        filtered = self.get_filtered(mssa=mssa, **kwargs)
        if '+' in mode: return filtered
        return self.span.fill_invalids(self._data, filtered)
        
    filled = property(fget=get_filled, doc='Data filled with filtered data where missing')


    def _set_field_(self, field, name, mean=True, std=False):
        """Put a field in self._<name>*"""
        fieldm = npy.ma.masked_values(field, default_missing_value, copy=False)
        setattr(self, '_%sm'%name, fieldm)
        if std:
            setattr(self, '_%s_std'%name, fieldm.std())
        if mean:
            setattr(self, '_%s_mean'%name, fieldm.mean(axis=1).reshape(-1, 1))

    def _link_field_(self, name, to='current'):
        """
        1.) Put a stored field and its attributes in self._current
        2.) Put data it in the current field to be analyzed (stacked data or PCs)
        """
        # Put it in self._current*
        if self._ana=='mssa':
            pass
        for att in 'm', '_mask', '_std', '_mean':
            if hasattr(self, '_%s'%name+att):
                setattr(self, '_'+to+att, getattr(self, '_'+name+att))

    def _set_raw_(self, data):
        if self._ana=='mssa':
            data = data.T
            raw_name = '_pca_raw_pc'
        else:
            raw_name = 'stacked_data'
        if npy.ma.isMA(data): data = data.data
        setattr(self.span, raw_name, data)
        

#        # Put it in analyzed field
#        if self._ana=="mssa":
#            self.span._pca_raw_pc = self._current
#        else:
#            self.span.stacked_data = self._current
    
        
    def _get_error_(self, anamode):
        """Get current reconstruction error"""
        diffm = self._refm-self._recm
        fieldm = self._refm
        if anamode==2:
            diffm = npy.ma.masked_where(self._cvfield_kept, diffm, copy=False)
            fieldm = npy.ma.masked_where(self._cvfield_kept, fieldm, copy=False)
        err = 100*diffm.compressed().std()/fieldm.compressed().std()
        if npy.isnan(err):
            raise FillerError('Error is NaN. Stop.')

        return err
        
     
    def _rec_(self, imode):
        """Get PCA or MSSA recontruction
        
        For PCA, output is same as raw input.
        For MSSA, output is same as PCA PCs input
        """
        
        # Expansion coefficients
        self._last_pcs = {} #getattr(self.span, '_%s_raw_pc'%self._ana)[0]
        if self._ana=='mssa' or not self.span.useteof: 
            
            # Fill masked point with reconstructed data
            rec = self._recm.data if hasattr(self._recm, 'data') else self._recm
            data =  self._currentm.filled(rec)
            
            if self._ana=='mssa': # MSSA rec
            
                kw = dict(xdata=data, xraw=True) if not self.span.prepca else {}
                self.span.mssa_ec(raw=True, replace=True, demean=False, **kw)
                
            else: # PCA rec
            
                self.span.pca_ec(raw=True, xraw=True, xdata=data, replace=True, demean=False)
                
            del data
    
        # Reconstruction (masked)
        recfunc = self._get_func_('rec')
        self._recm = recfunc(modes=-imode, raw=1, rescale=False, unmap=False)
        if hasattr(self, '_current_mean'):
            self._recm += self._current_mean
    
    def _get_func_(self, suf=None):
        """Get a PCA or MSSA related generic function (method)"""
        ana = self._ana
        if suf is not None:
            ana += '_'+suf
        return getattr(self.span, ana)
        
    def _get_necmiss_(self):
        """Get the number of missing values of current expansion coefficents"""
        if not hasattr(self,  '_ana'): return None, None
        if hasattr(self, '_necmiss'): self._necmiss_old = self._necmiss
        pc0 = npy.ma.masked_values(getattr(self.span, '_%s_raw_pc'%self._ana)[:, 0], 
            default_missing_value)
        self._necmiss =  pc0.size-npy.ma.count(pc0)
        del pc0
        if not hasattr(self, '_necmiss_old'):  
            necmissch = None
        else:
            necmissch = self._necmiss-self._necmiss_old
        return self._necmiss, necmissch
        
    def _gen_cvfield_(self, mask=None, level=1.0, regen=False, ntries=50, ncvmin=5):
        """Generate false missing values
        
        Field used for analysis is stored using `set_field('cvfield')`.
        Retained points are saved as logical array in :attr:`_cvfield_kept`.
        
        :Params:
        
            - **mask**, optional: Mask set to True where values are not retained for
              analysis. These rejected values are used to cross-validation.
              If not set, it is randomly generated.
            - **level**, optional: Approximate percentile of cross validation points. 
            - **mincv**, optional: Minimal number of cross-validation points.
        """
        if not regen and hasattr(self, '_cvfield') and self._cvfield_ana==self._ana: return

        # Generate the mask
        if mask is None:
            
            # Data mask
            data = self._refm
            ns, nt = self._refm.shape
            transp = self._ana=='pca'
            
            # Generate new mask
            if ntries>=0:
                dmask = npy.ma.getmaskarray(data)
                da0 = dmask.all(axis=0)
                da1 = dmask.all(axis=1)
                del dmask
            for reduc in npy.arange(1, 0, -0.1):
                for itry in xrange(min(1, ntries)):
                    
                    # Get new mask
                    mask = gen_cv_mask(data, level*reduc, merged=True, nmin=ncvmin)
                    
                    # Check consistency
                    if ntries<0 or (
                        npy.ma.allclose(mask.all(axis=0), da0) and 
                        npy.ma.allclose(mask.all(axis=1), da1)):
                        break
                else:
                    msg = 'Mask does not preserve minimal availability along axes'
                    if reduc!=0.1:
                        self.warning(msg+
                            'trying with  %i%% of level of cross validation points'%(reduc*100))
                    else:
                        self.warning(msg+'skipping.')
                    continue
                break
            
        # Apply mask
        cvfield = npy.where(mask, default_missing_value, self._refm)
        self._set_field_(cvfield, 'cvfield', mean=True, std=False)
        self._cvfield_kept = ~mask
        self._cvfield_ana = self._ana
        
        
def gen_cv_mask(data, level, merged=True, nmin=10):
    """Generate a cross validation mask with density depending on time availability
    
    :Params:
    
        - **dmask**: Channel/record (2D) data mask: True if masked.
        - **level**: Percent of cross-validation points in valid data.
        - **merged**, optional: Merge data mask and cv mask? Some of the data points
          chosen of cross-validation becomes masked.
        
    :Return: Depends on ``merge``:
        
        - True: ``mask``: False at cross-validation points only.
        - False: ``dmask``: Data mask with cross-validation points marked as False.
    """
    # Data mask
    dmask = npy.ma.getmaskarray(data)
    good = ~dmask 
    
    # Weights
    weights = npy.resize((data**2).sum(axis=1), data.shape[::-1]).T
    weights /= weights.max()

    # Max probability
    problim = weights[good] 
    ngood = problim.shape[0]
    probbad = problim==0.
    problim[probbad] = 1.
    prand = npy.random.uniform(0, 1, size=ngood)/problim
    prand[probbad] = 1
    psort = npy.argsort(prand)
    ncv = max(nmin, int(ngood*0.01*level))
    ccv = npy.zeros(problim.shape, '?')
    ccv[psort[:ncv]] = True
    cv = npy.zeros(dmask.shape, '?')
    cv[good] = ccv #; del good

    # Merge
    if merged:
        mask = dmask.copy()
        mask[cv] = True
    else:
        mask = ~cv
    
    return mask



if __name__=='__main__':
    N=npy
    s=N.sin(N.arange(500))
    c = N.linspace(5, 0.01, 6)
    data = N.dot(c.reshape(-1, 1), s.reshape(1, -1))
    mask =  gen_cv_mask(data, 5.)
    P.subplot(211)
    P.plot((mask).sum(axis=1));P.title('mask');
    P.subplot(212)
    P.plot((data**2).sum(axis=1));P.show()
#    

    
