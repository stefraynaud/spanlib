#################################################################################
# File: dual.py
#
# This file is part of the SpanLib library.
# Copyright (C) 2006-2012  Stephane Raynaud, Charles Doutriaux
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


import copy
import gc
from warnings import warn
import _fortran
import numpy as N
npy = N
from data import has_cdat_support, cdms2_isVariable, Data, Dataset, default_missing_value
if has_cdat_support: import MV2, cdms2
import pylab as P
from spanlib.util import Logger, broadcast, SpanlibIter, dict_filter
from analyzer import _BasicAnalyzer_, Analyzer, docs, _filldocs_

class DualAnalyzer(_BasicAnalyzer_, Logger):
    
    _nsvd_max = 100 
    _svd_params = ['nsvd']
    _params = dict(svd=_svd_params)
    _all_params = _svd_params
    _int_params = ['nsvd']
    
    def __init__(self, ldataset, rdataset, lweights=None, rweights=None, lnorms=None, rnorms=None,  
        lminvalid=None, rminvalid=None, logger=None, loglevel=None, zerofill=0, **kwargs):
            
        # Loggers        
        Logger.__init__(self, logger=logger, loglevel=loglevel, **dict_filter(kwargs, 'log_'))
        self._quiet=False
        
        # Left and right Analyzer instances
        if zerofill==2:
            kwargs['zerofill'] = 2
        kwleft, kwright = self._dict_filter_lr_(kwargs)
        kwargs['logger'] = self.logger
        self.lspan = Analyzer(ldataset, weights=lweights, norms=lnorms, 
            minvalid=lminvalid, **kwargs)
        self.rspan = Analyzer(rdataset, weights=rweights, norms=rnorms, 
            minvalid=rminvalid, **kwargs)
            
        # Check and save parameters
        self.update_params(None, **kwargs)
        
    
    @staticmethod
    def _dict_filter_lr_(kwargs):
        return dict_filter(kwargs, 'l', short=True), dict_filter(kwargs, 'r', short=True)  
    
    def update_params(self, anatype=None, **kwargs):
        """Initialize, update and check statistical paremeters.
        A value of  ``None`` is converted to an optimal value.
        Analyses are re-ran if needed by checking dependencies.
        
        :Params:
        
            - **anatype**, optional: current analysis type.
            
                - ``None``: Simple initialization with ``None``
                - ``"svd"``: Check SVD parameters.
               
              If different from ``None``, analysis may be ran again
              if parameters are changed.
              
            - Other keywords are interpreted as analysis parameters.

        :Output: A dictionary of (param name, change status) items.
        """
        
        # Update left and right params first
        kwleft, kwright = self._dict_filter_lr_(kwargs)
        reran = (self.lspan.update_params(**kwleft), 
            self.rspan.update_params(**kwright))
            
        # Initialize old values and defaults changed to False
        old = {}
        init_all = [None]*self.nd
        for param in self._all_params:
            old[param] = getattr(self,'_'+param,None)
            setattr(self, '_'+param, kwargs.pop(param, old[param], None))
            
        # Number of SVD modes
        if self._nsvd is None: # Initialization
            self._nsvd = self._nsvd_default # Default value
        for span in self.lspan, self.rspan: # Check values
            if span._prepca:
                nchanmax = span._prepca # Input channels are from pre-PCA
            else:
                nchanmax = self.ns # Input channels are from real space
            self._nsvd = npy.clip(self._nsvd, 1, min(SpAn._nsvd_max, nchanmax)) # Max
          
        # Re-run analyses when needed
        rerun = dict(
            svd = self._has_run_('svd') and (self._has_changed_(old, 'nsvd') or \
            (self.lspan._prepca and reran[0]['pca']) or \
            (self.rspan._prepca and reran[1]['pca'])), 
        )
        # - SVD
        if rerun['svd']:
            self.debug('Re-running SVD because some parameters changed')
            self.svd()

        # Inform what has reran
        return rerun
            
    def _has_run_(self, anatype=None, iset=0):
        """Check if an analysis has already run"""
        if anatype is None: return False
        return iset in getattr(self, '_%s_raw_eof'%anatype)
    
    def _svd_channel_axis_(self,iset):
        """Get the SVD channel axis for one dataset"""
        return _channel_axis_(iset, 'svd', svd=True)

    def get_nsvd(self):
        """Get :attr:`nsvd`        
        
        :Returns:
            integer
        """
        return self._nsvd
    def set_nsvd(self, nsvd):
        """Set :attr:`nsvd` and update analyzes"""
        self.update(nsvd=nsvd)
    nsvd = property(get_nsvd, set_nsvd, doc="Number of SVD modes")
    
    #################################################################
    ## SVD
    #################################################################
    @_filldocs_
    def svd(self, usecorr=False, largematrix=False, **kwargs):
        """ Singular Value Decomposition (SVD)

        It is called everytime needed by :meth:`svd_eof`, :meth:`svd_pc`, :meth:`svd_ev` and :meth:`svd_rec`.
        Thus, since results are stored in cache, it not necessary call it explicitly.

        :Parameters:
            %(nsvd)s
            %(prepca)s
            - *usecorr*: bool
                Use correlations instead of covariances.
        """

        if self.nd<2:
            self.error('Error you need at least (most) 2 datasets to run svd, otherwise use pca and mssa',  'svd')
        
        # Parameters
        self.update_params('svd', **kwargs)

        # Loop on first two datasets (left and right)
        for iset in xrange(2):

            # Check if old results can be used when nsvd is lower
            if iset in self._svd_raw_pc and \
                self._svd_raw_pc[iset].shape[-1] > self._nsvd:
                continue
            
            # Remove old results
            for att in 'raw_eof','raw_pc','raw_ev','ev_sum':
                dic = getattr(self,'_svd_'+att)
                if isinstance(dic, dict) and iset in dic: 
                    del dic[iset]

            # Prepare input to SVD
            if self._prepca[iset]: # Pre-PCA case
            
                # Pre-PCA
                if not self._pca_raw_pc.has_key(iset): self.pca(iset=iset)
                
                # svd
                data = self._pca_raw_pc[iset][:, :self._prepca[iset]].transpose()
                if iset == 0:
                    left = data
                    lweights = data[:, 0]*0.+1
                else:
                    right = data
                    rweights = data[:, 0]*0.+1

            else: # Direct svd case
                weights = self[iset].stacked_weights
                pdata = self[iset].stacked_data
                if iset == 0:
                    left = pdata
                    lweights = weights

                else:
                    right = pdata
                    rweights = weights

        # Compute SVD
        left = npy.asfortranarray(left, 'd')
        right = npy.asfortranarray(right, 'd')
        raw_eof_left, raw_eof_right, raw_pc_left, raw_pc_right, raw_ev, ev_sum, errmsg = \
            _fortran.svd(left, right, self._nsvd, lweights, rweights, int(usecorr))
        self.check_fortran_errmsg(errmsg)
        if info != 0:
            self.error('Error when running fortran SVD',  'svd')
            
        # Save results
        self._svd_raw_pc[0] = npy.ascontiguousarray(raw_pc_left)
        self._svd_raw_eof[0] = npy.ascontiguousarray(raw_eof_left)
        self._svd_raw_pc[1] = npy.ascontiguousarray(raw_pc_right)
        self._svd_raw_eof[1] = npy.ascontiguousarray(raw_eof_right)
        self._svd_raw_ev = raw_ev
        self._svd_ev_sum = ev_sum


        self._last_anatype = 'svd'
        gc.collect()

    @_filldocs_
    def svd_eof(self,iset=None,scale=False, raw=False, unmap=True, format=True, **kwargs):
        """Get EOFs from SVD analysis

        If SVD was not performed, it is done with all parameters sent to :meth:`svd`

        :Parameters:
            %(scale)s
            %(raw)s
            
        :SVD parameters:
            %(nsvd)s
            %(prepca)s
            
        :Returns:
            Arrays with shape ``(nsvd,...)``
        """

        # Update params
        self.update_params('svd', **kwargs)

        # Of, let's format the variable
        fmt_eof = {}
        for iset in xrange(self.nd): # (window*nchan,nsvd)
        
            # Operate only on selected datasets
            if isets is not None and iset not in isets: continue
        
            # EOF already available 
            if not raw and iset in self._svd_fmt_eof:
                fmt_eof[iset] = self._svd_fmt_eof[iset]
                continue
                
            # No analyses performed?
            if not iset in self._svd_raw_eof: self.svd(iset=iset)
            raw_eof = self._svd_raw_eof[iset]
            nm = self._nsvd
            nl = self._svd_raw_eof[iset].shape[0]
            raw_eof = self._svd_raw_eof[iset].reshape((nl, nm))
            
            if raw: # Do not go back to physical space
                self._svd_fmt_eof[iset] = [npy.ascontiguousarray(raw_eof.T)]
                if format and self[iset].has_cdat():
                    self._mssa_fmt_eof[iset][0] = cdms2.createVariable(self._mssa_fmt_eof[iset])
                    self._mssa_fmt_eof[iset][0].setAxisList(
                        [self._mode_axis_('svd',iset), self._svd_channel_axis_(iset)])
                    
            else: # Get raw data back to physical space
                raw_eof = npy.ascontiguousarray(raw_eof)
                raw_eof.shape = (nl, nm)
                
                firstaxes = self._mode_axis_('svd', iset)
                if not self._prepca[iset]: # No pre-PCA performed
                    self._svd_fmt_eof[iset] = self[iset].unstack(raw_eof,
                        rescale=False, format=1, firstaxes=firstaxes)
                        
                else:
                    proj_eof, smodes = self._raw_rec_(self._pca_raw_eof[iset], raw_eof.T)
                    self._svd_fmt_eof[iset] = self[iset].unstack(proj_eof, 
                        rescale=False, format=1, firstaxes = firstaxes)
                    
            # Set attributes
            for idata,eof in enumerate(self._svd_fmt_eof[iset]):
                
                # Attributes
                if format and cdms2_isVariable(eof):
                    if not raw and not self[iset][idata].id.find('variable_'):
                        eof.id = self[iset][idata].id+'_svd_eof'
                    else:
                        eof.id = 'svd_eof'
                    eof.name = eof.id
                    eof.standard_name = 'empirical_orthogonal_functions_of_svd'
                    eof.long_name = 'SVD empirical orthogonal functions'
                    if not raw:
                        atts = self[iset][idata].atts
                        if atts.has_key('long_name'):
                            eof.long_name += ' of '+atts['long_name']
                        if scale is True and atts.has_key('units'):
                            eof.units = atts['units']
                    
                # Scaling
                if scale:
                    if scale is True:
                        scale = npy.sqrt(self._mssa_raw_ev[iset])*nl
                    eof[:] *= scale
                    
            fmt_eof[iset] = self._svd_fmt_eof[iset]
            
        if not unmap: return fmt_eof
        return self._unmap_(fmt_eof, grouped=raw)

    @_filldocs_
    def svd_pc(self, iset=None, raw=False, unmap=True, format=True, **kwargs):
        """Get PCs from SVD analysis

        If SVD was not performed, it is done with all parameters sent to :meth:`svd`
        
        :Parameters:
            %(raw)s
            
        :SVD parameters:
            %(nsvd)s
            %(prepca)s
            
        :Returns:
            Arrays with shape ``(nsvd,nt)``
        """

        # Dataset selection
        isets = self._check_isets_(iset)

        # Update params
        self.update_params('svd', **kwargs)

        # Of, let's format the variable
        fmt_pc = {}
        for iset in xrange(self.nd):
            
            # PC already available 
            if not raw and iset in self._svd_fmt_pc:
                fmt_pc[iset] = self._svd_fmt_pc[iset]
                continue
                
            # No analyses performed?
            if not iset in self._svd_raw_pc: self.svd(iset=iset)
            
            # Raw?
            raw_pc = self._svd_raw_pc[iset][:,:self._nsvd]
            if raw:
                fmt_pc[iset] = raw_pc
                continue
                        
            # Format the variable
            idata = 0 # Reference is first data
            pc = npy.ascontiguousarray(raw_pc).T
            if self._cdat_inside_(iset):
                pc = cdms2.createVariable(pc)
                pc.setAxis(0,self._mode_axis_('svd',iset))
                pc.setAxis(1,self._time_axis_(iset, idata))
                pc.id = pc.name = 'svd_pc'
#                pc.standard_name = 'principal_components_of_svd'
                pc.long_name = 'SVD principal components'
                atts = self[iset][idata].atts
                if atts.has_key('long_name'): pc.long_name += ' of '+atts['long_name']
    #           if atts.has_key('units'):     pc.units = atts['units']

            fmt_pc[iset] = self._svd_fmt_pc[iset] = pc
            self._check_dataset_tag_('_svd_fmt_pc', iset, svd=True)

        if not unmap: return fmt_pc
        return self._unmap_(fmt_pc, grouped=True)       

    def svd_ec(self, iset=None, xdata=None, xeof=None, xraw=False, 
        raw=False, unmap=True, format=True, **kwargs):
        """Get expansion coefficients from SVD analysis
        
        :Parameters:
            %(raw)s

        :SVD parameters:
            %(nsvd)s
            %(window)s
            %(prepca)s
            
        :Returns:
            Arrays with the shape ``(nsvd,nt)``
        """

        # Dataset selection
        isets = self._check_isets_(iset)

        # Update params
        self.update_params('svd', **kwargs)

        # Remap
        xeof = self._remap_(xeof, grouped=xraw)
        xdata = self._remap_(xdata, grouped=xraw)
            
        fmt_ec = {}
        for iset in isets:
            
            # SVD EOFs used for projection
            if xeof[iset] is None:
                if not iset in self._svd_raw_eof: self.svd(iset=iset)
                raw_eof = self._svd_raw_eof[iset]
            elif xraw:
                raw_eof = xeof[iset]
            elif self._prepca[iset]: # After PCA
                raw_eof = self.pca_ec(iset=iset, xdata=xeof, raw=True, xscale=False, unmap=False)[iset].T
                raw_eof.shape = -1, self._nsvd
            else:
                eofs = self[iset].remap(xeof[iset])
                raw_eof = npy.ascontiguousarray(self[iset].restack(eofs, scale=False))
                del eofs
                
            # Data to project on SVD EOFs
            if xdata[iset] is None: # From input
                if not self._prepca[iset]: # Input data
                    raw_data = self[iset].stacked_data
                else: # After PCA
                    raw_data = self._pca_raw_pc[iset].T
            elif xraw: # Direct use
                raw_data = xdata[iset]
            elif self._prepca[iset]: # After PCA
                raw_data = self.pca_ec(iset=iset, xdata=xdata, raw=True, unmap=False)[iset].T
            else:
                data = self[iset].remap(xdata[iset])
                raw_data = self[iset].restack(data, scale=True)
                    
            
            # Projection
            raw_data = npy.asfortranarray(raw_data)
            raw_eof = npy.asfortranarray(raw_eof[:, :self._nsvd])
            raw_weights = npy.ones(raw_data.shape[0])
            raw_ec = _fortran.pca_getec(raw_data, raw_eof, raw_weights, default_missing_value)
            if raw:
                fmt_ec[iset] = raw_ec
                continue
            ec = npy.ascontiguousarray(raw_ec.T)
            ec = npy.ma.masked_values(ec, default_missing_value, copy=False)
            
            # Format the variable
            if format and self[iset].has_cdat():
                
                ec = cdms2.createVariable(ec)
                ec.setAxis(0,self._mode_axis_('svd',iset))
#                ec.setAxis(1,self._svd_pctime_axis_(iset))
                ec.id = ec.name = 'svd_ec'
#                ec.standard_name = 'expansion coefficient_of_svd of '
                ec.long_name = 'SVD principal components'
                atts = self[iset][0].atts
                if len(self[iset])==1 and atts.has_key('long_name'): 
                    ec.long_name += atts['long_name']
    #           else:
    #               ec.long_name += 'of dataset %i'%iset
                #if (len(self[iset])==1 or npy.allclose(self.norms(iset), 1.)) and atts.has_key('units'):
                #   ec.units = atts['units']
            fmt_ec[iset] = ec
            # FIXME: _check_dataset_tag_
#            self._check_dataset_tag_('_svd_fmt_ec', iset)

        if not unmap: return fmt_ec
        return self._unmap_(fmt_ec, grouped=True)


    @_filldocs_
    def svd_ev(self,relative=False,sum=False,cumsum=False,**kwargs):
        """Get eigen values from SVD analysis

        :Parameters:
          %(relative)s
          %(sum)s
          %(cumsum)s
        
            
        :SVD parameters:
            %(nsvd)s
            %(prepca)s
            
        :Returns:
            Array with shape ``(nsvd,)`` or a float
        """

        # Update params
        self.update_params('svd', **kwargs)

        # No analyses performed?
        if not self._svd_raw_eof.has_key(0): self.svd()

        # We only want the sum
        if sum: return self._svd_ev_sum

        # Format the variable
        id = 'svd_ev'
        long_name = []
        ev = self._svd_raw_ev[:self._nsvd]
        if cumsum:
            ev = raw_ev.cumsum()
            id += '_cumsum'
            long_name.append('cumulative')
        if relative: 
            ev = 100.*raw_ev/self._svd_ev_sum
            id += '_rel'
            long_name.append('relative')
            
        # Format the variables   
        if self._cdat_inside_(0) or self._cdat_inside_(1):
            ev = cdms2.createVariable(ev)
            ev.id = ev.name = id
            long_name = ['SVD eigen values']
            ev.long_name = ' '.join(long_name).title()
            ev.setAxis(0, self._mode_axis_('svd'))
#            ev.standard_name = 'eigen_values_of_svd'
            if cumsum:
                self.id += '_cumsum'
                long_name.append('cumulative')
            if relative:
                ev.id += '_rel'
                ev.units = '% of total variance'
                long_name.append('cumulative')
            ev.long_name = ' '.join(long_name).title()
        return ev

    @_filldocs_
    def svd_rec(self, iset=None, modes=None, raw=False, xpc=None, xeof=None, xraw=False, 
        rescale=True, unmap=True, format=2, **kwargs):
        """Reconstruction of SVD modes
        
        :Parameters:
            %(modes)s
            %(raw)s
            
        :SVD parameters:
            %(nsvd)s
            %(prepca)s
            
        :Returns:
            Arrays with the same shape as input arrays.
        """
        
        # Dataset selection
        isets = self._check_isets_(iset)

        # Update params
        self.update_params('svd', **kwargs)
        
        # Remap alternate arrays
        xeof = self._remap_(xeof, grouped=xraw)
        xpc = self._remap_(xpc)
            
        # Loop on datasets
        fmt_rec = {}
        for iset in isets:

            # No analyses performed?
            if not iset in self._svd_raw_eof: self.svd(iset=iset)
            
            # ST-EOFs used for projection
            if xeof[iset] is None:
                if not iset in self._svd_raw_eof: self.svd(iset=iset)
                raw_eof = self._svd_raw_eof[iset]
            elif xraw:
                raw_eof = xeof[iset]
            elif self._prepca[iset]: # After PCA
                raw_eof = self.pca_ec(iset=iset, xdata=xeof, raw=True, xscale=False, unmap=False)[iset].T
            else:
                eofs = self[iset].remap(xeof[iset])
                raw_eof = npy.ascontiguousarray(self[iset].restack(eofs, scale=False))
                del eofs
            raw_eof = raw_eof[:, :self._nsvd]
                
            # PCs
            if xpc[iset] is None: # From PCA
                if not iset in self._svd_raw_pc: self.svd(iset=iset)
                raw_pc = self._svd_raw_pc[iset]
            elif xraw: # Direct use
                raw_pc = xpc[iset]
            elif npy.ma.isMA(xpc[iset]): # formatted -> pure numpy
                raw_pc = xpc[iset].filled(default_missing_value).T
            else:
                raw_pc = xpc[iset].T
            raw_pc = raw_pc[:, :self._nsvd]
            
            # Projection
            raw_rec,smodes = self._raw_rec_(raw_eof, raw_pc, modes)
                
            # Get raw data back to physical space (nchan,nt)
            taxis = self._time_axis_(iset)
            if not self._prepca[iset]: # No pre-PCA performed
                fmt_rec[iset] = self[iset].unstack(raw_rec, rescale=rescale, format=format)
                
            elif raw: # Force direct result from svd
                if format and self[iset].has_cdat():
                    fmt_rec[iset] = [cdms2.createVariable(raw_rec.T)]
                    fmt_rec[iset][0].setAxisList(0, 
                        [taxis,self._mode_axis_('svd',iset)])
                else:
                    fmt_rec[iset] = raw_rec.T
                    
            else: # With pre-pca
                proj_rec, spcamodes = self._raw_rec_(self._pca_raw_eof[iset], raw_rec.T)
                fmt_rec[iset] = self[iset].unstack(proj_rec, rescale=rescale, format=format)
            del  raw_rec
            
            # Set attributes
            if format:
                for idata,rec in enumerate(fmt_rec[iset]):
                    if not cdms2_isVariable(rec): continue
                    if not self[iset][idata].id.startswith('variable_'):
                        rec.id = self[iset][idata].id+'_svd_rec'
                    else:
                        rec.id = 'svd_rec'
    #               if modes is not None:
                    rec.id += smodes #FIXME: do we keep it?
                    rec.modes = smodes
    #               else:
    #                   rec.modes = '1-%i'%self._nsvd[iset]
                    rec.standard_name = 'recontruction_of_svd_modes'
                    rec.long_name = 'Reconstruction of SVD modes'
                    atts = self[iset][idata].atts
                    if atts.has_key('long_name'):
                        rec.long_name += ' of '+atts['long_name']
                   
        if not unmap: return fmt_rec
        return self._unmap_(fmt_rec, grouped=raw)   
    

class SVDModel(DualAnalyzer):
    """
    :Params:
        - **method**, optional: Method of reconstruction [default: 'std']. 
        
            - 'std' assumes that left and normalized expansion coefficients are equal [Syu_and_Neelin_1995]_.
            - 'regre' does not use right EOFs but regression coefficients [Harrisson_et al_2002]_.
    """

    def __init__(self,predictor, predictand, method='std', **kwargs):

        DualAnalyzer.__init__(self, predictor, predictand ,**kwargs)

        # Perform an SVD between the first two datasets
        self.learn(nsvd=None, npca=None)

        # Setup conversion between the two datasets
        if method.startswith('s'): # standard deviation ratio
            self._l2r = self._svd_raw_pc[1].std(axis=0)/self._svd_raw_pc[0].std(axis=0)
        elif self.method.startswith('r'): # regression
            from npy.linalg import lstsq
            a = self._svd_raw_pc[0].T
            b = self._svd_raw_pc[1].T
            self._l2r, res, rank, s = lstsq(a, b)
        else:
            self.error('Wrong method for right coefficients from left coefficients. '
                'It should one of : "regre",  "std"')
        self.method = method
        

    def _left2right_(self, leftcoefs):
        if self.method.startswith('s'):
            return self._l2r*leftcoefs
        if self.method.startswith('r'):
            if leftcoefs.ndim==1:
                leftcoefs = leftcoefs.reshape(-1, leftcoefs.shape[-1])
            else:
                leftcoefs = leftcoefs.T
            leftcoefs = npy.atleast_2d(leftcoefs)
            return npy.dot(leftcoefs.T, self._l2r)
            
    def learn(self, **kwargs):
        """Learning phase"""
#       self.clean()
        self.svd(**kwargs)
        
    def run(self, predictor, method='pcs', modes=None):
        """Run the SVD model 
        
        
        :Example:
            
            # Init
            >>> model = SVDModel(predictor,predictand)
            # Run method
            >>> predicted1 = model.run(new_predictor1)
            # Call method
            >>> predicted2 = model(new_predictor2)
        """
        
        
        # Predictor as SpAn datasets
        predictor = self._remap_([predictor])[0]
        
        # Get expansion coefficients for predictor only
        lec = self.svd_ec(iset=0, xdata=predictor, raw=True, format=False, unmap=False)[0]
        
        # Convert to predictand coeffisients
        rec = self._left2right_(lec)
            
        # Reconstruct
        rrec = self.svd_rec(iset=1, xpc=rec, xraw=True)[0]
        return rrec
        

    __call__ = run

    
    def clean(self):
        SpAn.clean(self)
        self._regre_ = None
        gc.collect()



   
