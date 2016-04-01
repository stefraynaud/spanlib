# File: spanlib_python.py
#
# This file is part of the SpanLib library.
# Copyright (C) 2006-2014  Stephane Raynaud, Charles Doutriaux
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

import gc
import numpy as N
npy = N
from data import (has_cdat_support, cdms2_isVariable, Data, Dataset,
    default_missing_value)
if has_cdat_support:
    import MV2, cdms2
#from .util import Logger, broadcast, SpanlibIter, dict_filter
#from spanlib.util import Logger, broadcast, SpanlibIter, dict_filter
# import _core
from .util import Logger, broadcast, SpanlibIter, dict_filter

docs = dict(
    npca="""- *npca*: int | ``None``
                Number of PCA modes to keep in analysis (defaults to 10).""",
    prepca="""- *prepca*: int | bool | ``None``
                Number of pre-PCA modes to keep before MSSA and SVD analyses (defaults to ``npca`` if ``True``).
                If the number of input channels is greater than 30, it is automatically switch to ``True``.""",
    nmssa="""- *nmssa*: int | ``None``
                Number of MSSA modes to keep in analysis (defaults to 10).""",
    window="""- *window*: int | ``None``
                Size of the MSSA window parameter (defaults to 1/3 the time length).""",
    nsvd="""- *nsvd*: int | ``None``
                Number of SVD modes to keep in analysis (defaults to 10).""",
    modes="""- *modes*: int | list | tuple
                If ``None``, all modes are summed.
                Example of other usages:

                    - ``4`` or ``[4]`` or ``(4,)``: only mode 4
                    - ``-4``: modes 1 to 4
                    - ``(1,3,-5)``: modes 1, 3, 4 and 5 (``-`` means "until")""",
    raw="""- *raw*: bool
                When pre-PCA is used and ``raw`` is ``True``, it prevents from
                going back to physical space (expansion to PCA EOFs space).""",
    scale="""- *scale*: bool | float
                Apply a factor to EOFs or PCs. This is essentially useful to add
                a quantitative meaning. If ``True``, ``scale`` is chosen so that
                the standard deviation of the mode (EOF or PC) is the square of
                the associated eigen value.""",
    relative="""- *relative*: bool
                Return percentage of variance instead of its absolute value.""",
    sum="""- *sum*: bool
                Return the sum of ALL (not only the selected modes) eigen values (total variance).""",
    cumsum="""- *cumsum*: bool
                Return the cumulated sum of eigen values.""",
)

def _filldocs_(func):
    func.__doc__ = func.__doc__ % docs
    return func


class _BasicAnalyzer_(object):
    """Define fundamental methods for all analysers"""

    @staticmethod
    def _get_imodes_(imode, nmode):
        """First mode at 0"""

        # Which modes
        if imode is None:
            imode = range(0, nmode)#+1)
        elif isinstance(imode, slice):
            imode = range(*imode.indices(nmode))
        else:
            if isinstance(imode,int):
                imode = [imode,]

        # Rearrange modes (imode=[0,3,4,5,9] -> [0,0],[3,5],[9,9])
        imode = [im for im in imode if im < nmode]
        imode.sort(cmp=lambda x,y: cmp(abs(x),abs(y)))
        if imode[0]<0: imode.insert(0,0)
        imodes = []
        im = 0
        while im < len(imode):
            imode1 = imode2 = imode[im]
            for imt in xrange(im+1,len(imode)):
                if imode[imt] > 0  and (imode[imt]-imode2) > 1: # end of group
                    im = imt-1
                    break
                imode2 = abs(imode[imt]) # increase group
                continue
            else:
                if im < len(imode)-1: im = imt
            im += 1
            imodes.append([imode1,imode2])
        return imodes

#    @staticmethod
#    def _check_length_(input, mylen, fillvalue):
#        return broadcast(input, mylen, fillvalue)

    def _has_changed_(self, old, param):
        """Check if a parameter has changed

        :Returns:

            - ``None`` if the parameter is new
            - ``False`` if not integer and not changed
            - ``True`` if not intager and changed
            - ``-1`` if integer, changed and lower
            - ``1`` if integer, changed and upper
            - ``0`` if integer and not changed
        """
        oldv = old[param]
        newv = getattr(self, param)
        if oldv is None: return
        if param in self._int_params: return newv > oldv
        return oldv != newv


    def _mode_axis_(self, analysis_type):
        """Get a mode axis according to the type of modes (pca, mssa, svd)

        If CDAT is not used, length of axis is returned.
        """
        if analysis_type not in self._mode_axes:
            self._mode_axes[analysis_type] = {}
        nn = getattr(self,'_n'+analysis_type)
        if not self.has_cdat():
                self._mode_axes[analysis_type] = nn
        elif self._mode_axes[analysis_type] is None or \
            len(self._mode_axes[analysis_type]) != nn:
            self._mode_axes[analysis_type] = cdms2.createAxis(npy.arange(1,nn+1))
            self._mode_axes[analysis_type].id = analysis_type+'_mode'
            self._mode_axes[analysis_type].long_name = analysis_type.upper()+' modes in decreasing order'
        if analysis_type == 'svd': return self._mode_axes[analysis_type]
        return self._mode_axes[analysis_type]

    def _channel_axis_(self, name, **kwargs):
        """Get the channel axis for one dataset (MSSA or SVD)

        If CDAT is not used, length of axis is returned.
        """
        if not self._prepca:
            nchan = self.ns
        else:
            nchan = self._prepca
        channel_axes = getattr(self,'_%s_channel_axes'%name)
        if not self.has_cdat():
            channel_axes = nchan
        if channel_axes is None or len(channel_axes) != nchan:
            channel_axes = cdms2.createAxis(npy.arange(nchan))
            channel_axes.id = '%s_channel' % name
            channel_axes.long_name = '%s channels'% name.upper()
        return channel_axes

    def replace_invalids(self, data, raw=False, copy=False):
        self.clean()
        return Dataset.replace_invalids(self, data, raw=raw, copy=copy)


    def _cleanattr_(self, att, init=None):
        """Delete an attribute and set it to None"""
        if not hasattr(self, att) or \
            ((callable(init) and not isinstance(getattr(self, att), init))) or \
            getattr(self, att) is not init:
            if callable(init): init = init()
            if hasattr(self, att): delattr(self, att)
            setattr(self, att, init)


class Analyzer(_BasicAnalyzer_, Dataset):
    """ Prepare the Spectral Analysis Object

    Description:
      This function creates an object for future analyses.
      It optionally initializes some parameters.


      data  :: List of data on which to run the PC Analysis
        Last dimensions must represent the spatial dimensions.
        Analysis will be run on the first dimension.
      weights :: If you which to apply weights on some points,
        set weights to "0" where you wish to mask.
        The input data mask will be applied,
        using the union of all none spacial dimension mask.
        If the data are on a regular grid, area weights
        will be generated, if the cdutil (CDAT) module is available.
        default: 1. everywhere]
      npca  :: Number of principal components to return [default: 10]
      nmssa   :: Number of MSSA modes retained [default: 4]
      nsvd  :: Number of SVD modes retained [default: 10]
      window  :: MSSA window parameter [default: time_length/3.]


    Output:
      analysis_object :: SpAn object created for further analysis

    """

    # pylint: disable=too-many-instance-attributes
    _npca_default = 10
    _npca_max = 200
    _nprepca_max = 30
    _useteof_default = -1
    _notpc_default = 0
    _minecvalid_default = 0
    _zerofill_default = 0
    _nmssa_default = _nsvd_default = 8
    _nmssa_max = 20 #_nsvd_max = 20
    _window_default = 1/3. # Relative to time length
    _pca_params = ['npca', 'prepca', 'minecvalid', 'zerofill', 'useteof',
        'notpc', 'pcapf']
    _mssa_params = _pca_params+['nmssa', 'prepca', 'window']
#    _svd_params = _pca_params+['nsvd']
    _params = dict(pca=_pca_params, mssa=_pca_params+_mssa_params)#, svd=_pca_params+_svd_params)
    _all_params = list(set(_pca_params+_mssa_params))#+_svd_params
#    _common_params = ['nsvd'] # common to all datasets
    _int_params = ['npca', 'nmssa', 'prepca', 'window', 'nsvd']

    def __init__(self, dataset, weights=None, norms=None,
            minvalid=None, clean_weights=True, keep_invalids=False, zerofill=0,
            logger=None, loglevel=None, **kwargs):

        # Create Dataset instance
        Dataset.__init__(self, dataset, weights=weights, norms=norms, zerofill=zerofill,
            minvalid=minvalid, clean_weights=clean_weights, keep_invalids=keep_invalids)
        self._quiet=False

        # Init results
        self.clean()

        # Check and save parameters
        if zerofill==2:
            kwargs['zerofill'] = 2
        self.update_params(None, **kwargs)


    #################################################################
    ## Get datasets info
    #################################################################

    def _mssa_channel_axis_(self):
        """Get the MSSA channel axis

        If CDAT is not used, length of axis is returned.
        """
        return self._channel_axis_('mssa')

    def _mssa_window_axis_(self, update=False):
        """Get the MSSA window axis for one dataset

        If CDAT is not used, length of axis is returned.
        """
        if not self.has_cdat():
            self._mssa_window_axes = self._window
        elif self._mssa_window_axes is None or \
            len(self._mssa_window_axes) != self._window:
            self._mssa_window_axes = cdms2.createAxis(npy.arange(self._window))
            self._mssa_window_axes.id = 'mssa_window'
            self._mssa_window_axes.long_name = 'MSSA window time'
        return self._mssa_window_axes

    def _mssa_pctime_axis_(self, idata=0):
        """Get the MSSA PCs time axis for one dataset

        If CDAT is not used, length of axis is returned.
        """
        nt = self.nt - self._window + 1
        if not self.has_cdat():
            self._mssa_pctime_axes = nt
        elif self._mssa_pctime_axes is None or \
            len(self._mssa_pctime_axes) != nt:
            self._mssa_pctime_axes = cdms2.createAxis(npy.arange(nt).astype('d'))
            self._mssa_pctime_axes.id = 'mssa_pctime'
            self._mssa_pctime_axes.long_name = 'MSSA PC time'
            taxis = self[idata].get_time()
            if hasattr(taxis,'units') and taxis.units.split()[0].lower() in \
                ['seconds','minutes','hours','days','months','years']:
                self._mssa_pctime_axes.units = taxis.units.split()[0].lower() + ' since 0001-01-01'
                self._mssa_pctime_axes.designateTime()
        return self._mssa_pctime_axes

#    def _check_dataset_tag_(self, name, key=None, long_name=True, id=True, svd=False):
#        """Mark some attributes as specific to a dataset (only if there are more then one dataset)
#            iset:: ID of the dataset
#            key:: A dictionary key to select the dataset [default: None]
#            long_name:: Mark the long name [defualt: True]
#            id:: Mark the id [default: True]
#            svd: Mark using 'left' or 'right'
#        """
#        if self.nd > 1:
#            targetset = getattr(self, name)
#            if key is not None:
#                targetset = targetset[key]
#            target = targetset[iset]
#            if not cdms2_isVariable(target): return
#            if svd:
#                svdtag = ['left', 'right'][iset]
#            if id:
#                if svd:
#                    target.id  = '%s_%s'%(svdtag, target.id)
#                else:
#                    target.id += '_set%i'%iset
#            if long_name:
#                if svd:
#                    target.long_name += ' for %s dataset'%svdtag
#                else:
#                    target.long_name += ' for dataset #%i'%iset

    def _cdat_ev_(self, ev, analysis_type, relative, cumsum,
        id=None, long_name=None, standard_name=None,  atts=None):
        """Format eigenvalues as CDAT variable"""
        if not self.has_cdat(): return ev
        ev = cdms2.createVariable(ev)
        if id is None:
            id = analysis_type.lower()+'_ev'
        if long_name is None:
            long_name = []
            if cumsum:
                id += '_cumsum'
                long_name.append('cumulative')
            if relative:
                id += '_rel'
                long_name.append('relative')
        ev.id = ev.name = id
        if isinstance(long_name, list):
            long_name.append(analysis_type.upper()+' eigen values')
            ev.long_name = ' '.join(long_name).title()
            atts = self[0].atts
            if atts.has_key('long_name'):
                ev.long_name += ' of '+atts['long_name']
        else:
            ev.long_name = long_name
        ev.setAxisList([self._mode_axis_(analysis_type.lower())])
        if standard_name is not False:
            if standard_name is None:
                ev.standard_name = 'eigen_values_of_'+analysis_type.lower()
            else:
                ev.standard_name = standard_name
        if relative:
            ev.units = '% of total variance'
        return ev



    def _cdat_inside_(self, idata=None):
        """Check if a data var has CDAT support"""
        # Not at all
        if not has_cdat_support: return False

        # Single var
        if idata is not None:
            return self[idata].has_cdat()

        # Full dataset == at least one var
        return self.has_cdat()

    def _has_run_(self, anatype=None):
        """Check if an analysis has already run"""
        if anatype is None: return False
        return getattr(self, '_%s_raw_eof'%anatype) is not None

    def _has_changed_(self, old, param):
        """Check if a parameter has changed

        :Returns:

            - ``None`` if the parameter is new
            - ``False`` if not integer and not changed
            - ``True`` if not intager and changed
            - ``-1`` if integer, changed and lower
            - ``1`` if integer, changed and upper
            - ``0`` if integer and not changed
        """
        oldv = old[param]
        newv = getattr(self, '_'+param)
        if oldv is None: return
        if param in self._int_params: return newv > oldv
        return oldv != newv

    def update_params(self, anatype=None, checkprepca=False, **kwargs):
        """Initialize, update and check statistical paremeters.
        A value of  ``None`` is converted to an optimal value.
        Analyses are re-ran if needed by checking dependencies.

        :Params:

            - **anatype**, optional: current analysis type.

                - ``None``: Simple initialization with ``None``
                - ``"pca"``: Check PCA parameters.
                - ``"mssa"``: Check MSSA parameters.
                - ``"svd"``: Check SVD parameters.

              If different from ``None``, analysis may be ran again
              if parameters are changed.

            - Other keywords are interpreted as analysis parameters.

        :Output: A dictionary of (param name, change status) items.
        """

        # Initialize old values and defaults changed to False
        if 'win' in kwargs:
            kwargs['window'] = kwargs.pop('win')
        old = {}
        for param in self._all_params:
            old[param] = getattr(self, '_'+param, None)
            setattr(self, '_'+param, kwargs.pop(param, old[param]))

        # Number of PCA modes
        # - Guess a value
        if self._npca is None:
            if self._prepca is not None:
                self._npca = self._prepca
            else:
                self._npca = SpAn._npca_default # Default value
        # - Min
        if self._prepca is not None:
            self._npca = max(self._npca, self._prepca)
        # - max
        self._npca = npy.clip(self._npca, 1, min(SpAn._npca_max, self.ns, self.nt))

        # Number of pre-PCA modes before MSSA and SVD
        if self._prepca is None: # Default: pre-PCA needed over max (for MSSA and SVD)
            self._prepca = self.ns > SpAn._nprepca_max
            if not hasattr(self, '_warned_prepca'): self._warned_prepca = {}
            if not self._quiet and self._prepca and \
                self._warned_prepca is False:
                self._warned_prepca = True
                self.warning('[mssa] The number of valid points (%i) '
                    ' is greater than %i, so we perform a pre-PCA'%
                    (self.ns, SpAn._nprepca_max))
        if self._prepca is True: # Defaults to the number of PCA modes
            self._prepca = self._npca
        elif self._prepca: # Max number of prepca modes is number of points
            self._prepca = min(self._prepca, self.ns, self.nt)
        if self._prepca == 0:
            self._prepca = False

        # Dependency rules between prepca and npca
        if self._prepca and self._npca < self._prepca:
            if not self._quiet  and self._prepca:
                    self.warning('The number of pre-PCA modes (%i) is lower' \
                        ' than the number of PCA modes (%i), so we adjust the latter.' \
                        % (self._prepca, self._npca))
            self._npca = self._prepca

        # Use T-EOF?
        if self._useteof is None:
#            self._useteof = SpAn._useteof_default
            self._useteof = self.ns > self.nt
        if self._notpc is None:
            self._notpc = SpAn._notpc_default

        # Fill with zeros?
        if self._zerofill is None:
            self._zerofill = SpAn._zerofill_default

        # Min number of values for expansion coefficients?
        if self._minecvalid is None:
            self._minecvalid = SpAn._minecvalid_default

        # Window extension of MSSA
        if self._window is None: # Initialization
            self._window = SpAn._window_default
        if self.window <0:
            self._window = int(npy.round(npy.clip(-self._window, 0., 100)*self.nt/100))
        elif self._window<1:
            self._window = int(npy.round(self._window*self.nt))
        self._window = int(npy.clip(self._window, 1, max(1, self.nt)))

        # Number of MSSA modes
        if self._nmssa is None: # Initialization
            # Guess a value
            self._nmssa = SpAn._nmssa_default # Default value
        if self._prepca:
            nchanmax = self._prepca # Input channels are from pre-PCA
        else:
            nchanmax = self.ns # Input channels are from real space
        self._nmssa = int(npy.clip(self._nmssa,1,
            min(SpAn._nmssa_max,nchanmax*self._window))) # Max

        # Re-run analyses when needed
        rerun = {}
        rerun['pca'] = self._has_run_('pca') and (
            self._has_changed_(old, 'npca')>0 or
            self._has_changed_(old, 'useteof') or
            self._has_changed_(old, 'notpc')
        )
        rerun['mssa'] = self._has_run_('mssa') and (
            self._has_changed_(old, 'nmssa')<0 or
            self._has_changed_(old, 'window') or
            (self._prepca and rerun['pca'])
        )
        # - PCA
        if rerun['pca']:
            self.debug('Re-running PCA because some important parameters changed')
            self.pca()
        # - MSSA
        if rerun['mssa']:
            self.debug('Re-running MSSA because some important parameters changed')
            self.mssa()

        # Inform what has reran
        return rerun

#    def _check_isets_(self,iset):
#        """Check if an iset is associated to a a valid dataset.
#
#        It can be a list, and it is returned as a list.
#        If an iset is invalid, it is removed from the output list.
#        """
#        if iset is None: return range(self.nd)
#        if iset == 'left':
#            iset = 0
#        elif iset == 'right':
#            iset = 1
#        if iset < 0 or iset >= self.nd:
#            warn('Invalid dataset id: %i. Valid id are < %i'%(iset,self.nd))
#        else:
#            return [iset]

#    def _check_shape_(self, inputs, fillvalue):
#        """Return input as datasets (tree) *shape*"""
#        imap = self._input_map
#        if isinstance(imap, int):
#            return [SpAn._check_length_(inputs, max(1, imap), fillvalue), ]
#        inputs = SpAn._check_length_(inputs,len(imap),fillvalue)
#        for iset, im in enumerate(imap):
#            inputs[iset] = SpAn._check_length_(inputs[iset], max(1, im), fillvalue)
#        return inputs


    #################################################################
    ## PCA
    #################################################################
    @_filldocs_
    def pca(self, force=False, **kwargs):
        """
        Principal Components Analysis (PCA)

        It is called everytime needed by :meth:`pca_eof`, :meth:`pca_pc`,
        :meth:`pca_ev` and :meth:`pca_rec`.
        Thus, since results are stored in cache, it not necessary call it explicitly.

        :Parameters:
            %(npca)s
        """

        # Update params
        self.update_params('pca', **kwargs)

        # Check if old results can be used when npca is lower
        if not force and self._pca_raw_pc is not None and \
            self._pca_raw_pc.shape[-1] >= self._npca:
            return

        # Remove old results
        for att in 'raw_eof','raw_pc','raw_ev','ev_sum':
            self._cleanattr_('_pca_'+att)

        # Compute PCA
        pdata = self.stacked_data
        if pdata.shape[1] == 1: # One single channel, so result is itself
            raw_eof = npy.ones(1, dtype=pdata.dtype)
            raw_pc = pdata
            ev_sum = raw_pc.var()
            raw_ev = npy.atleast_1d(ev_sum)

        else: # Several channels
#            weights = npy.asfortranarray(self.stacked_weights)
#            pdata = npy.asfortranarray(pdata)
            raw_eof, raw_pc, raw_ev, ev_sum, errmsg = \
                _core.pca(pdata, self._npca, default_missing_value,
                useteof=self._useteof, notpc=self._notpc, minecvalid=self._minecvalid,
                zerofill=self._zerofill)
            self.check_fortran_errmsg(errmsg)

        # Post filtering
        if callable(self._pcapf):
            self._pcapf(pdata, raw_eof, raw_pc, raw_ev, default_missing_value)

        # Save results
        self._pca_raw_pc = raw_pc
        self._pca_raw_eof = raw_eof
        self._pca_raw_ev = raw_ev
        self._pca_ev_sum = ev_sum

        # Delete formatted variables
        for vtype in 'pc', 'eof':
            self._cleanattr_('_pca_fmt_'+vtype)
        gc.collect()

        self._last_anatype = 'pca'

    def pca_has_run(self):
        """Check if PCA has already run"""
        return self._has_run_('pca')


    @_filldocs_
    def pca_eof(self, scale=False, raw=False, unmap=True, format=True, **kwargs):
        """Get EOFs from PCA analysis

        :Parameters:
            %(scale)s
            %(raw)s

        :PCA parameters:
            %(npca)s

        :Returns:
            Arrays with shape ``(npca,...)``
        """

        # Update params
        self.update_params('pca', **kwargs)


        # EOF already available
        if not raw and self._pca_fmt_eof is not None:
            return self._pca_fmt_eof

        # First PCA analysis?
        if self._pca_raw_eof is None: self.pca()

        # Raw
        eof = self._pca_raw_eof
        if raw:
            return eof

        # Mask
        eof = npy.ascontiguousarray(eof)
        eof = npy.ma.masked_values(eof, default_missing_value, copy=False)

        # Back to physical space
        self._pca_fmt_eof = self.unstack(eof, rescale=False,
            firstaxes=self._mode_axis_('pca'))

        # Set attributes and scale
        for idata,eof in enumerate(self._pca_fmt_eof):

            # Scaling
            if scale:
                if scale is True: # Std dev of EOF is sqrt(ev)
                    scale = npy.sqrt(self._pca_raw_ev*(self.ns-1))
                    for imode in xrange(eof.shape[0]):
                        eof[imode] *= scale[imode]
                else:
                    eof *= scale

            # Attributes (CDAT)
            if format and cdms2_isVariable(eof):
                if not self[idata].id.startswith('variable_'):
                    eof.id = self[idata].id+'_pca_eof'
                else:
                    eof.id = 'pca_eof'
                eof.name = eof.id
                eof.standard_name = 'empirical_orthogonal_functions_of_pca'
                eof.long_name = 'PCA empirical orthogonal functions'
                atts = self[idata].atts
                if atts.has_key('long_name'):
                    eof.long_name += ' of '+atts['long_name']
                if scale and atts.has_key('units'):
                    eof.units = atts['units']


        if unmap: return self.unmap(self._pca_fmt_eof)
        return self._pca_fmt_eof


    @_filldocs_
    def pca_pc(self, scale=False, raw=False, format=True, **kwargs):
        """Get principal components (PCs) from current PCA decomposition

        :Parameters:

        :PCA parameters:
            %(npca)s

        :Returns:
            Arrays with the shape ``(npca,nt)``
        """
        # Update params
        self.update_params('pca', **kwargs)

        # PC already available
        if self._pca_fmt_pc is not None:
            return self._pca_fmt_pc

        # First PCA analysis?
        if self._pca_raw_pc is None: self.pca()

        # Raw ?
        if raw:
            return self._pca_raw_pc[:,:self._npca]

        # Mask
        pc = npy.ascontiguousarray(self._pca_raw_pc[:,:self._npca].T)
        pc = npy.ma.masked_values(pc, default_missing_value, copy=False)

        # Format (CDAT)
        #TODO: scale pca_pc
        if format and self.has_cdat():
            pc = cdms2.createVariable(pc)
            pc.setAxis(0, self._mode_axis_('pca'))
            pc.setAxis(1, self[0].get_time())
            pc.id = pc.name = 'pca_pc'
            pc.standard_name = 'principal_components_of_pca'
            pc.long_name = 'PCA principal components of '
            atts = self[0].atts
            if len(self)==1 and  atts.has_key('long_name'):
                pc.long_name += atts['long_name']
            if scale and (len(self) == 1 or npy.allclose(self.norms, 1.)) and atts.has_key('units'):
                pc.units = atts['units']


        self._pca_fmt_pc = pc

        return pc

    def pca_ec(self, xdata=None, xeof=None, scale=False, ev=None,
        xraw=False, xscale=True, raw=False, unmap=True, format=True,
        replace=False, **kwargs):
        """Get expansion coefficient using current PCA decomposition

        Expansion coefficients are the projection of data onto EOFs.
        By default, it uses input data and computed EOFs, and thus
        returns principal components.
        You can bypass this default behaviour using ``xeof`` and ``xdata``
        keyword parameters.

        :Parameters:

        :PCA parameters:
            %(npca)s

        :Returns:
            Arrays with the shape ``(npca,nt)``
        """
        # Remap
        if not xraw:
            if xeof is not None: xeof = self.remap(xeof)
            if xdata is not None: xdata = self.remap(xdata)

        # EOFs used for projection
        if xeof is None: # From PCA
            if self._pca_raw_eof is None: self.pca()
            raw_eof = self._pca_raw_eof
        elif xraw: # Direct use
            raw_eof = xeof
        else: # We format then use
            eofs = self.remap(xeof)
            raw_eof = self.restack(eofs, scale=False)

        # Data to project on EOFs
        if xdata is None: # From input
            raw_data = self.stacked_data
            ndim = raw_data.ndim
        elif xraw: # Direct use
            raw_data = xdata
            ndim = raw_data.ndim
        else: # We format then use
            data = self.remap(xdata)
            raw_data = self.restack(data, scale=xscale)
            ndim = raw_data.ndim
            if raw_data.ndim>2:
                raw_data = npy.reshape(raw_data, (raw_data.shape[0], -1))

        # Eigenvalues for normalisation
        if ev is None: ev = self._pca_raw_ev
        if ev is False: ev = ev*0-1

        # Projection
        raw_data = npy.asfortranarray(raw_data)
        raw_eof = npy.asfortranarray(raw_eof)
        raw_ec = _core.pca_getec(raw_data, raw_eof, mv=default_missing_value,
            minvalid=self._minecvalid, zerofill=self._zerofill)

        # Replace current pc with computed ec
        if replace:
            self._pca_raw_pc = raw_ec
            if self._pca_fmt_pc is None:
                self._cleanattr_('_pca_fmt_pc')

        # Raw?
        if raw:
            return raw_ec

        # Mask
        ec = npy.ascontiguousarray(raw_ec.T)
        ec = npy.ma.masked_values(ec, default_missing_value, copy=False)

        # Format (CDAT)
        if format and self.has_cdat():
            ec = cdms2.createVariable(ec)
            ec.setAxis(0, self._mode_axis_('pca'))
            if ndim==2 and ec.shape[1]==self.nt:
                ec.setAxis(1, self[0].get_time())
            ec.id = ec.name = 'pca_ec'
            ec.long_name = 'PCA expansion coefficients'
            atts = self[0].atts
            if len(self)==1 and atts.has_key('long_name'):
                ec.long_name += ' of '+atts['long_name']
            if scale and (len(self)==1 or npy.allclose(self.norms, 1.)) and atts.has_key('units'):
                ec.units = atts['units']

        return ec


    @_filldocs_
    def pca_ev(self, relative=False, sum=False, cumsum=False, format=True, **kwargs):
        """Get eigen values from current PCA decomposition

        :Parameters:
          %(relative)s
          %(sum)s
          %(cumsum)s


        :PCA parameters:
            %(npca)s

        :Returns:
            Arrays with shape ``(npca,)`` or a float
        """

        # Update params
        self.update_params('pca', **kwargs)

        # First PCA analysis?
        if self._pca_raw_eof is None: self.pca()

        # We only want the sum
        if sum:
            return self._pca_ev_sum

        # Data
        ev = self._pca_raw_ev[:self._npca]
        if cumsum:
            ev = ev.cumsum()
        if relative:
            ev = 100.*ev/self._pca_ev_sum

        # Format (CDAT)
        if format and self.has_cdat():

            id = 'pca_ev'
            long_name = []
            if cumsum:
                id += '_cumsum'
                long_name.append('cumulative')
            if relative:
                id += '_rel'
                long_name.append('relative')
            ev = cdms2.createVariable(ev)
            ev.id = ev.name = id
            long_name.append('PCA eigen values')
            ev.long_name = ' '.join(long_name).title()+' of '
            ev.setAxisList([self._mode_axis_('pca')])
            ev.standard_name = 'eigen_values_of_pca'
            atts = self[0].atts
            if len(self)==1 and atts.has_key('long_name'):
                ev.long_name += atts['long_name']
            if relative:
                ev.units = '% of total variance'
            elif (len(self)==1 or npy.allclose(self.norms, 1.)) and atts.has_key('units'):
                ev.units = atts['units']
                for ss in ['^','**',' ']:
                    if ev.units.find(ss) != -1:
                        ev.units = '(%s)^2' % ev.units
                        break

        return ev


    @_filldocs_
    def pca_rec(self, modes=None, raw=False, xpc=None, xeof=None, xraw=False,
        rescale=True, format=2, unmap=True, **kwargs):
        """Reconstruct a set of modes from PCA decomposition

        :Parameters:
            %(modes)s
            %(raw)s

        :PCA parameters:
            %(npca)s

        :Returns:
            Arrays with the same shape as input arrays.
        """
        # Update params
        self.update_params('pca', **kwargs)

        # Remap alternate arrays
        if not xraw:
            if xeof is not None: xeof = self.remap(xeof)
#        xpc = self.remap(xpc)

        # First PCA analysis?
        if self._pca_raw_pc is None: self.pca()

        # EOFs
        if xeof is None: # From PCA
            raw_eof = self._pca_raw_eof
        elif xraw: # Direct use
            raw_eof = xeof
        else: # We format then use
            eofs = self.remap(xeof)
            raw_eof = self.restack(eofs, scale=False)

        # PCs
        if xpc is None: # From PCA
            raw_pc = self._pca_raw_pc
        elif xraw: # Direct use
            raw_pc = xpc
        else:
            raw_pc = xpc.T

        # Reconstruction
        reof = raw_eof[:,:self._npca]
        rpc = raw_pc[:,:self._npca]
        raw_rec, smodes = self._raw_rec_(reof, rpc, modes)

        # Raw?
        if raw:
            return raw_rec


        # Back to physical space
        if rescale: rescale = 2
        pca_fmt_rec = self.unstack(raw_rec, rescale=rescale, format=format)
        del  raw_rec

        # Format (CDAT)
        if format:
            for idata,rec in enumerate(pca_fmt_rec):
                if not cdms2_isVariable(rec): continue
    #               rec[:] *= self._norm_(idata) # Scale
                if not self[idata].id.startswith('variable_'):
                    rec.id = self[idata].id+'_pca_rec'
                else:
                    rec.id = 'pca_rec'
                rec.name = rec.id
    #               if modes is not None:
                rec.id += smodes
                rec.modes = smodes
                rec.standard_name = 'recontruction_of_pca_modes'
                rec.long_name = 'Reconstruction of PCA modes: '+smodes
                atts = self[idata].atts
                if atts.has_key('long_name'):
                    rec.long_name += ' of '+atts['long_name']
                if atts.has_key('units'):
                    rec.units = atts['units']

        if not unmap: return pca_fmt_rec
        return self.unmap(pca_fmt_rec)


    #################################################################
    # MSSA
    #################################################################



    def preproc_raw_output(self, force=False):
        """Get preprocessing raw output

        It is either the original data (:attr:`stacked_data`)
        if not pre-PCA must be performed or the first attr:`prepca`
        PCA PCs :attr:`_pca_raw_pcs` with mean removed.
        """
        # Pre-PCA case
        if self._prepca:

            # PCA
            self.pca(force=int(force)==2)

            # Compute the pre-PCs mean (not always zero!) for future reconstructions
            pca_raw_pc = self._pca_raw_pc[:, :self._prepca]
            pca_raw_pc_masked = npy.ma.masked_values(pca_raw_pc,
                default_missing_value, copy=False)
            self._pca_raw_pc_mean = pca_raw_pc_masked.mean(axis=0).filled(0.)
            del pca_raw_pc_masked
            pca_raw_pc -= self._pca_raw_pc_mean
            self._pca_raw_pc_mean.shape = -1, 1
            return npy.asfortranarray(pca_raw_pc.T)


        # Direct MSSA case
        return npy.asfortranarray(self.stacked_data)



    @_filldocs_
    def mssa(self, force=False, **kwargs):
        """ MultiChannel Singular Spectrum Analysis (MSSA)

        It is called everytime needed by :meth:`mssa_eof`, :meth:`mssa_pc`, :meth:`mssa_ev` and :meth:`mssa_rec`.
        Thus, since results are stored in cache, it not necessary call it explicitly.

        :Parameters:
            %(nmssa)s
            %(window)s
            %(prepca)s
        """

        # Parameters
        self.update_params('mssa', **kwargs)

        # Check if old results can be used when nmssa is lower
        if not force and self._mssa_raw_pc is not None and self._mssa_raw_pc.shape[-1] >= self._nmssa:
            return

        # Remove old results
        for att in 'raw_eof','raw_pc','raw_ev','ev_sum':
            self._cleanattr_('_mssa_'+att)

        # Get input to MSSA
        raw_input = self.preproc_raw_output(force=force)

        # Run MSSA
        raw_eof, raw_pc, raw_ev, ev_sum, errmsg = \
            _core.mssa(raw_input, self._window, self._nmssa,
                default_missing_value, minecvalid=self._minecvalid,
                zerofill=self._zerofill)
        self.check_fortran_errmsg(errmsg)

        # Save results
        self._mssa_raw_pc = raw_pc
        self._mssa_raw_eof = raw_eof
        self._mssa_raw_ev = raw_ev
        self._mssa_ev_sum = ev_sum

        # Delete formmated variables
        for vtype in 'pc', 'eof':
            self._cleanattr_('_mssa_fmt_'+vtype)

        self._last_anatype = 'mssa'
        gc.collect()

    def mssa_has_run(self):
        """Check if MSSA has already run"""
        return self._has_run_('mssa')

    @_filldocs_
    def mssa_eof(self, scale=False, raw=False, format=True, unmap=True, **kwargs):
        """Get EOFs from MSSA analysis

        Shape: (window*nchan,nmssa)

        :Parameters:
            %(scale)s
            %(raw)s

        :MSSA parameters:
            %(nmssa)s
            %(window)s
            %(prepca)s

        :Returns:
            Arrays with shape ``(nmssa,nt-window+1,...)``
        """

        # Update params
        self.update_params('mssa', **kwargs)

#        # EOF already available
#        if not raw and self._mssa_fmt_eof is not None:
#            return self._mssa_fmt_eof

        # No analyses performed?
        if self._mssa_raw_eof is None: self.mssa()

        # Raw eof
        raw_eof = self._mssa_raw_eof[:, :self._nmssa]
        nlxnw, nm = raw_eof.shape
        nw = self._window
        nl = nlxnw/nw
        raw_eof = raw_eof.reshape((nl, nw, nm))

        if raw: # Do not go back to physical space

            return raw_eof
#            self._mssa_fmt_eof = [npy.ascontiguousarray(raw_eof.T)]

#            if format and self.has_cdat(): # Fromat (CDAT)
#                self._mssa_fmt_eof[0] = cdms2.createVariable(self._mssa_fmt_eof)
#                self._mssa_fmt_eof[0].setAxisList(
#                    [self._mode_axis_('mssa'),self._mssa_channel_axis_()])

        else: # Get raw data back to physical space

            raw_eof = npy.ascontiguousarray(raw_eof) # (nl, nw, nm)
            firstaxes = (self._mode_axis_('mssa'), self._mssa_window_axis_())
            if not self._prepca: # No pre-PCA performed
                self._mssa_fmt_eof = self.unstack(raw_eof, rescale=False,
                    format=format, firstaxes=firstaxes)

            else: # With pre-PCA

                proj_eof, smodes = self._raw_rec_(self._pca_raw_eof,
                    raw_eof.T.reshape((nm*nw, nl)))
                self._mssa_fmt_eof = self.unstack(proj_eof,
                    rescale=False, format=format, firstaxes = firstaxes)

        # Scaling
        if scale:
            if scale is True:
                scale = npy.ma.sqrt(self._mssa_raw_ev)*nl*nw
            for idata, eof in enumerate(self._mssa_fmt_eof):
                eof[:] *= scale

        # Format (CDAT)
        if format :
            for idata, eof in enumerate(self._mssa_fmt_eof):
                if cdms2_isVariable(eof):
                    if not raw and not self[idata].id.find('variable_'):
                        eof.id = self[idata].id+'_mssa_eof'
                    else:
                        eof.id = 'mssa_eof'
                    eof.name = eof.id
                    eof.standard_name = 'empirical_orthogonal_functions_of_mssa'
                    eof.long_name = 'MSSA empirical orthogonal functions'
                    if not raw:
                        atts = self[idata].atts
                        if atts.has_key('long_name'):
                            eof.long_name += ' of '+atts['long_name']
                        if scale and atts.has_key('units'):
                            eof.units = atts['units']


        gc.collect()
        if not unmap or raw: return self._mssa_fmt_eof
        return self.unmap(self._mssa_fmt_eof)

    @_filldocs_
    def mssa_pc(self, raw=False, unmap=True, format=True, **kwargs):
        """Get PCs from MSSA analysis

        :Parameters:

        :MSSA parameters:
            %(nmssa)s
            %(window)s
            %(prepca)s

        :Returns:
            Arrays with the shape ``(nmssa,nt)``
        """

        # Update params
        self.update_params('mssa', **kwargs)

#        # PC already available
#        if not raw and self._mssa_fmt_pc is not None:
#            return self._mssa_fmt_pc

        # No analyses performed?
        if self._mssa_raw_pc is None: self.mssa()

        # Raw?
        raw_pc = self._mssa_raw_pc[:,:self._nmssa]
        if raw:
            return raw_pc

        # Mask
        pc = npy.ascontiguousarray(raw_pc.T)
        pc = npy.ma.masked_values(pc, default_missing_value, copy=False)

        # Format (CDAT)
        if format and self.has_cdat():

            pc = cdms2.createVariable(pc)
            pc.setAxis(0,self._mode_axis_('mssa'))
            pc.setAxis(1,self._mssa_pctime_axis_())
            pc.id = pc.name = 'mssa_pc'
            pc.standard_name = 'principal_components_of_mssa of '
            pc.long_name = 'MSSA principal components'
            atts = self[0].atts
            if len(self)==1 and atts.has_key('long_name'):
                pc.long_name += atts['long_name']
            #if (len(self)==1 or npy.allclose(self.norms, 1.)) and atts.has_key('units'):
            #   pc.units = atts['units']

        self._mssa_fmt_pc = pc
        return pc

    def mssa_ec(self, xdata=None, xeof=None, xraw=False,
        raw=False, unmap=True, format=True, replace=False, demean=True, **kwargs):
        """Get expansion coefficients from MSSA analysis

        :Parameters:

        :MSSA parameters:
            %(nmssa)s
            %(window)s
            %(prepca)s

        :Returns:
            Arrays with the shape ``(nmssa,nt)``
        """

        # Update params
        self.update_params('mssa', **kwargs)

        # Remap
        if not xraw:
            if xeof is not None: xeof = self.remap(xeof)
            if xdata is not None: xdata = self.remap(xdata)
        xraw = int(xraw)

        # ST-EOFs used for projection
        if xeof is None:
            if self._mssa_raw_eof is None: self.mssa()
            raw_eof = self._mssa_raw_eof
        elif xraw:
            raw_eof = xeof
        elif self._prepca: # After PCA
            raw_eof = self.pca_ec(xdata=xeof, xscale=False, raw=True,
                unmap=False, demean=0).T
            raw_eof.shape = -1, self._nmssa
        else:
            eofs = self.remap(xeof)
            raw_eof = npy.ascontiguousarray(self.restack(eofs, scale=False))
            raw_eof.shape = -1, raw_eof.shape[-1] # (nc*nw, nm)
            del eofs
        if self._prepca:
            nw = self._window
            nc = raw_eof.shape[0]/nw
            if nc > self._prepca:
                raw_eof = raw_eof.reshape(nc, nw, -1)
                raw_eof = raw_eof[:self._prepca]
                raw_eof = raw_eof.reshape(self._prepca*nw, -1)

        # Data to project on ST-EOFs
        if xdata is None: # From input
            if not self._prepca: # Input data
                raw_data = self.stacked_data
            else: # After PCA
                raw_data = self._pca_raw_pc.T[:self._prepca]
        elif int(xraw)==1: # Direct use
            raw_data = xdata
        elif self._prepca: # After PCA
            raw_data = self.pca_ec(xdata=xdata, raw=True, unmap=False,
                xraw=xraw==2, demean=int(demean)).T[:self._prepca]
        else:
            data = self.remap(xdata)
            raw_data = self.restack(data, scale=True)


        # Projection
        raw_data = npy.asfortranarray(raw_data)
        raw_eof = npy.asfortranarray(raw_eof[:, :self._nmssa])
        raw_ec = _core.mssa_getec(raw_data, raw_eof, self.window,
            default_missing_value,
            minvalid=self._minecvalid, zerofill=self._zerofill)
        if replace:
            self._cleanattr_('_mssa_raw_pc', raw_ec)
            self._cleanattr_('_mssa_fmt_pc')
        if raw:
            return raw_ec

        # Mask
        ec = npy.ascontiguousarray(raw_ec.T) ; del raw_ec
        ec = npy.ma.masked_values(ec, default_missing_value, copy=False)

        # Format
        if format and self.has_cdat():

            ec = cdms2.createVariable(ec)
            ec.setAxis(0,self._mode_axis_('mssa'))
#                ec.setAxis(1,self._mssa_pctime_axis_())
            ec.id = ec.name = 'mssa_ec'
#                ec.standard_name = 'expansion coefficient_of_mssa of '
            ec.long_name = 'MSSA principal components'
            atts = self[0].atts
            if len(self)==1 and atts.has_key('long_name'):
                ec.long_name += atts['long_name']
            #if (len(self)==1 or npy.allclose(self.norms, 1.)) and atts.has_key('units'):
            #   ec.units = atts['units']


        return ec


    @_filldocs_
    def mssa_ev(self, relative=False, sum=False, cumsum=False,
        mctest=False, mcnens=100, mcqt=90, format=True, unmap=True, **kwargs):
        """Get eigen values from current MSSA decomposition

        :Options:

          %(relative)s
          %(sum)s
          %(cumsum)s


        :MSSA options:

            %(nmssa)s
            %(window)s
            %(prepca)s
            - *mctest*: bool
                If ``True``, perfom a Monte-Carlo test if significance
                eigenvalues against red noise, following [Allen_and_Smith_96]_.
            - *mcnens*: int
                Size of the ensemble fo compute quantiles.
            - *mcqt*: float|int
                Value of the higher quantile in %% (the lower is ``100-mcqt``).

        :Returns:

            - Arrays with shape ``(nmssa,)`` or a float
            - Optionally, ``(ev, evmin, evmax)`` if ``mctest is True``

        :Basic example:

        >>> span = SpAn(data)
        >>> ev = span.mssa_ev()
        >>> ev_sum = span.mssa_ev(sum=True)
        >>> ev_cumsum = span.mssa_ev(cumsum=True)
        >>> ev_rel = span.mssa_ev(relative=True)

        :Monte-Carlo test example:

        >>> span = SpAn(data)
        >>> ev, mc_evmin, mc_evmax = span.mssa_ev(mctest=True, mcqt=95, mcnens=100)
        >>> print (ev<evmin)|(ev>evmax)
        """

        # Update params
        self.update_params('mssa', **kwargs)

        # No analyses performed?
        if self._mssa_raw_eof is None: self.mssa()

        # We only want the sum
        if sum:
            return self._mssa_ev_sum

        # Data
        ev = [self._mssa_raw_ev[:self._nmssa], ]

        # Mont-Carlo test
        if mctest:

            # Get reference data
            if self._prepca:
                data = self._pca_raw_pc
            else:
                data = self.stacked_data

            # Inits
            rn = RedNoise(data.T) # red noise generator
            mcev = npy.zeros((mcnens, self.nmssa))

            # Generate and ensemble of surrogate data
            for iens in xrange(mcnens):

                # Create a sample red noise (nt,nchan)
                red_noise = rn.sample().T

                # Block-covariance matrix (nchan*nwindow,nchan*nwindow,)
                cov = _core.stcov(npy.asfortranarray(red_noise),
                    self.window, default_missing_value)
                cov = npy.ascontiguousarray(cov)
                del red_noise

                # Fake eigen values (EOFt.COV.EOF)
                ce = npy.dot(cov, self._mssa_raw_eof) #; del cov
                evmat = npy.dot(self._mssa_raw_eof.T, ce) #; del ce
                mcev[iens] = npy.diag(evmat)[:nmssa] ; del evmat

            mcev.sort(axis=0) # Sort by value inside ensemble

            # Confidence interval
            # - min
            imini, iminf = divmod((1-mcqt/100.)*(mcnens-1), 1)
            evmin = mcev[int(imini)]
            if int(imini) != mcnens-1:
                evmin += (1-iminf)*mcev[int(imini)] + iminf*mcev[int(imini)+1]
            # - max
            imaxi, imaxf = divmod(mcqt/100.*(mcnens-1), 1)
            evmax = mcev[int(imaxi)]
            if int(imaxi) != mcnens-1:
                evmax += (1-imaxf)*mcev[int(imaxi)] + imaxf*mcev[int(imaxi)+1]
            ev.extend([evmin, evmax])

        # Special outputs
        if cumsum:
            ev = [e.cumsum() for e in ev]
        if relative:
            ev = [(100.*e/self._mssa_ev_sum) for e in ev]

        # Format the variables
        if format and self.has_cdat():

            for i, e in enumerate(ev):
                ev[i] = self._cdat_ev_(e,  'mssa', relative, cumsum)
            if mctest:
                ev[1].id += '_mclow'
                ev[1].long_name += ': lower bound of MC confidence interval (%g%%)'%mcqt
                ev[2].id += '_mchigh'
                ev[2].long_name += ': upper bound of MC confidence interval (%g%%)'%(100-mcqt)

        return ev[0] if not mctest else ev




    @_filldocs_
    def mssa_rec(self, modes=None, raw=False,
        xpc=None, xeof=None, xev=None, xraw=False,
        phases=False, rescale=True, format=2, unmap=True, evrenorm=False, **kwargs):
        """Reconstruction of MSSA modes

        :Parameters:
            %(modes)s
            %(raw)s
            *phases* : ``False`` | ``True`` | int
                Return phases composites of the reconstructed field.
                By default, 8 phases are computed. You can psecify
                the number of phases by passing an integer.

        :MSSA parameters:
            %(nmssa)s
            %(window)s
            %(prepca)s

        :Returns:
            Arrays with the same shape as input arrays.
        """
        # Update params
        self.update_params('mssa', **kwargs)
        if isinstance(phases, int):
            self._nphase = phases
        elif phases is not False:
            phases = self._nphase

        # Remap alternate arrays
        if not xraw and xeof is not None:
            xeof = self.remap(xeof)

        # Loop on datasets
        mssa_fmt_rec = {}
        raw = int(raw)

        # ST-EOFs used for projection
        if xeof is None:
            if self._mssa_raw_eof is None: self.mssa()
            raw_eof = self._mssa_raw_eof
        elif xraw:
            raw_eof = xeof
        elif self._prepca: # After PCA
            raw_eof = self.pca_ec(xdata=xeof, raw=True, xscale=False, unmap=False).T
            raw_eof.shape = -1, self.nmssa
        else:
            eofs = self.remap(xeof)
            raw_eof = npy.ascontiguousarray(self.restack(eofs, scale=False))
            raw_eof.shape = -1, raw_eof.shape[-1] # (nl*nw, nm)
            del eofs
        if npy.ma.isMA(raw_eof):
            raw_eof = raw_eof.filled(default_missing_value)
        raw_eof = raw_eof[:, :self.nmssa]
        ncxnw, nm = raw_eof.shape
        nw = self.window
        nc = ncxnw/nw
        raw_eof = raw_eof.reshape((nc, nw, nm))
        if self._prepca and nc > self.prepca:
            raw_eof = raw_eof[:self.prepca]

        # PCs
        if xpc is None: # From PCA
            if self._mssa_raw_pc is None: self.mssa()
            raw_pc = self._mssa_raw_pc
        elif xraw: # Direct use
            raw_pc = xpc
        else:
            raw_pc = xpc.T
        if npy.ma.isMA(raw_pc):
            raw_pc = raw_pc.filled(default_missing_value)
        raw_pc = raw_pc[:, :self.nmssa]

        # Eigenvalues
        if xev is None:
            raw_ev = self._mssa_raw_ev
        else:
            raw_ev = npy.asarray(xev)

        # Projection
        kw = {} if not evrenorm else {'ev':raw_ev}
        raw_rec, smodes = self._raw_rec_(raw_eof, raw_pc, modes, **kw)

        # Back to physical space and format
        return self._mssa_rec_format_(raw_rec, raw, phases, rescale, format,
            unmap, xpc, smodes)


    def _mssa_rec_format_(self, raw_rec, raw=False, phases=False, rescale=True,
        format=2, unmap=True, xpc=None, smodes=''):
        """Format data from an MSSA reconstruction"""


        # Mask
        if not npy.ma.isMA(raw_rec):
            raw_rec = npy.ma.masked_values(raw_rec,  default_missing_value, copy=False)

        # Phases composites
        if phases:
            raw_rec_t = phase_composites(raw_rec.T, phases, format=format and self.has_cdat)
            nt = raw_rec.shape[0]
            if cdms2_isVariable(raw_rec_t):
                raw_rec = MV2.transpose(raw_rec_t)
                taxis = raw_rec.getAxis(1)
            else:
                raw_rec = raw_rec_t.asma().T
                taxis = nt
            del raw_rec_t
        else:
            if raw_rec.shape[1]==self.nt and (xpc is None or xpc is self._mssa_raw_pc):
                taxis = self.get_time()
            else:
                taxis = self.get_time(nt=raw_rec.shape[1])

        # Get raw data back to physical space (nchan,nt)
        if not self.prepca: # No pre-PCA performed
            if raw:

                mssa_fmt_rec = raw_rec

                if format and self.has_cdat() and int(raw)==12: # Formatted

                    mssa_fmt_rec = [cdms2.createVariable(raw_rec.T)]
                    mssa_fmt_rec[0].setAxisList([taxis, self._mode_axis_('mssa')])

                return mssa_fmt_rec

            else:
                mssa_fmt_rec = self.unstack(raw_rec, rescale=rescale, format=format,
                    firstaxes=[taxis])

        else: # With pre-pca

            # Add back the pre-PCs mean
            if rescale:
                pca_raw_pc_mean = npy.repeat(self._pca_raw_pc_mean, raw_rec.shape[1], 1)
                raw_rec += pca_raw_pc_mean
                del pca_raw_pc_mean

            # No PCA reprojection
            if int(raw)==1 or int(raw)==12: # Force direct result from MSSA

                if format and self.has_cdat() and int(raw)==12: # Formatted

                    mssa_fmt_rec = [cdms2.createVariable(raw_rec.T)]
                    mssa_fmt_rec[0].setAxisList([taxis,self._mode_axis_('mssa')])

                else: # Pure raw

                    mssa_fmt_rec = raw_rec

                return raw_rec

            else: # Reprojection

                # Project
                proj_rec, spcamodes = \
                    self._raw_rec_(self._pca_raw_eof[:, :self.prepca], raw_rec.T)

                if int(raw)==2: # Raw results from PCA rec
                    mssa_fmt_rec = proj_rec

                else: # Back to original format
                    mssa_fmt_rec = self.unstack(proj_rec, rescale=rescale, format=format,
                        firstaxes=[taxis])

        del  raw_rec

        # Remove the mean for phases
        if phases:
            if raw:
                mssa_fmt_rec -= mssa_fmt_rec.mean(axis=0)
            else:
                for rec in mssa_fmt_rec:
                    rec[:] -= rec.mean(axis=0)

        # Set attributes
        if format:
            for idata,rec in enumerate(mssa_fmt_rec):
                if not cdms2_isVariable(rec): continue
                if not self[idata].id.startswith('variable_'):
                    rec.id = self[idata].id+'_mssa_rec'
                else:
                    rec.id = 'mssa_rec'
#               if modes is not None:
                rec.id += smodes #FIXME: do we keep it?
                rec.modes = smodes
#               else:
#                   rec.modes = '1-%i'%self.nmssa
                rec.standard_name = 'recontruction_of_mssa_modes'
                rec.long_name = 'Reconstruction of MSSA modes'
                atts = self[idata].atts
                if atts.has_key('long_name'):
                    rec.long_name += ' of '+atts['long_name']

        if unmap: return self.unmap(mssa_fmt_rec)
        return mssa_fmt_rec

    def mssa_phases(self, pair, nphase=8, format=True, unmap=True, **kwargs):
        """Build phase composites using an MSSA pair

        .. note::

            This method is special call to :meth:`mssa_rec`
            where ``nphase`` is not zero and ``modes`` are set
            to ``pair``.

        :Parameters:

            - *pair*: A tuple designating the pair of mode to
              reconstruct. If a integer is given, this mode and the
              the following are used.

        :Options:

            - *iset*: the dataset to work on.
            - *nphase*: The number of phase composites.
            - Other options are passed to :meth:`mssa_rec`

        """
        if isinstance(pair, int):
            pair = (pair, pair+1)

        return self.mssa_rec(modes=pair, phases=nphase, format=format, remap=remap)

#    def rec(self, anatype=None, *args, **kwargs):
#        """Generic method for reconstruction"""
#        if anatype is None:
#            anatype = self._last_anatype
#        else:
#            valid = ['pca','mssa']
#            if anatype not in valid:
#                raise SpanlibException('rec','anatype must be one of '+str(valid))
#        if anatype is None:
#            warnings.warn('Yet no statistics performed, so nothing to reconstruct!')
#        else:
#            return getattr(self,self._last_anatype+'_rec')(*args,**kwargs)

    def clean(self, pca=True, mssa=True):
        """(Re-)Initialization"""
        anatypes = []
        if pca:
            anatypes.append('pca')
            mssa = True
        if mssa:
            anatypes.append('mssa')

        # Register what to clean
        nones = []
        for aa in anatypes:
            nones.append('_%s_ev_sum'%aa)
            for bb in 'raw','fmt':
                for cc in 'eof','pc','ev':
                    nones.append('_%s_%s_%s'%(aa,bb,cc))
        if pca:
            nones.extend(['_pca_raw_pc_mean'])
        if mssa:
            nones.extend(['_mssa_window_axes','_mssa_pctime_axes', '_mssa_channel_axes'])
#        lists = ['_mssa_pairs','_nt','_ns','_ndata','_pdata']
        dicts = ['_mode_axes']

        # Clean now using values or functions
        for ll,init in [(nones, None), (dicts, dict)]:#,(lists,list):
            for att in ll:
                self._cleanattr_(att, init)

        # Integers
#        self.nd = 0
        self._nphase = 8
        self._last_anatype = None
        gc.collect()


    def get_npca(self):
        """Get :attr:`npca`

        :Returns:
            integer
        """
        return self._npca
    def set_npca(self, npca):
        """Set :attr:`npca`"""
        self.update_params(npca=npca)
    npca = property(get_npca, set_npca, doc="Number of PCA modes")

    def get_prepca(self):
        """Get :attr:`prepca`

        :Returns:
            integer
        """
        return self._prepca
    def set_prepca(self, prepca):
        """Set :attr:`prepca`"""
        self.update_params(prepca=prepca)
    prepca = property(get_prepca, set_prepca,
        doc="Number of PCA modes used before MSSA or SVD for d.o.f reduction")

    def get_useteof(self):
        return self._useteof
    def set_useteof(self, value):
        self.update_params(useteof=value)
        return self._useteof
    useteof = property(fget=get_useteof, fset=set_useteof, doc="Use T-EOFs for PCA?")


    def get_nmssa(self):
        """Get :attr:`nmssa`

        :Returns:
            integer or tuple
        """
        return self._nmssa
    def set_nmssa(self, nmssa):
        """Set :attr:`nmssa`"""
        self.update_params(nmssa=nmssa)
    nmssa = property(get_nmssa, set_nmssa, doc="Number of MSSA modes")

    def get_window(self, absolute=False):
        """Get :attr:`window`

        :Options:

            - *absolute*: if False, return the window relative to time length,
              else return the effective window (multiplied by nt)

        :Returns:
            integer or tuple
        """
        win = self._window
        if absolute:
            return int(win*self.nt)
        return win
    def set_window(self, win):
        """Set :attr:`window`"""
        self.update_params(window=window)
    win = window = property(get_window, set_window, doc="Size of MSSA window")

    def get_nc(self):
        """Get :attr:`nc`"""
        if self.prepca:
            return self.prepca
        return self.ns
    nc = property(get_nc,
        doc="Number of channels used for MSSA: either :attr:`prepca` or :attr:`ns`")

#    #################################################################
#    ## Organize datasets
#    #################################################################
#
#    def _map_(self, datasets, sequential):
#        """Make datasets in the right form (a list of variables or lists of variables)"""
#
#         # Input data tree type
#        if not isinstance(datasets, (tuple, list)):
#            self.tree_type = 0
#        else:
#            self.tree_type = 1
#            for dataset in datasets:
#                if isinstance(dataset, (tuple, list)):
#                    self.tree_type = 2
#                    break
#            else:
#                if sequential: self.tree_type = 2
#
#        # Group variables if needed
#        if self.tree_type < 2:
#            datasets = [datasets]
#
#        self._datasets = datasets
#        self.nd = len(datasets)

#    def _remap_(self, values, reshape=True, fill_value=None, grouped=False):
#        """Makes sure that values is a list (or tuple) of length :attr:`ndatasets`
#
#        :func:`broadcast` is called when reshaping.
#        """
#        # We always need a sequence
#        if not isinstance(values, (tuple, list)):
#            values = [values]
#
#        # Check length
#        if not grouped and values[0] is not None and not isinstance(values[0], (tuple, list)):
#            values = [values]
#        if len(values)!=len(self):
#            if not reshape:
#                self.error('Wrong number of input items (%i instead of %i)'
#                    %(len(values), len(self)))
#            values = broadcast(values, len(self), fill_value)
#
#        return values

#    def _unmap_(self, values, grouped=None):
#        """Return values as input dataset (depth and shapes)
#
#        :Params:
#
#            - **values**: List, list of lists or dictionary with integers as keys.
#            - **grouped**, optional: Do not unmap as the :class:`Dataset` level.
#        """
#
#        # Convert list to dictionary
#        if isinstance(values, list):
#            if grouped is None: grouped = not isinstance(values[0], list)
#            ff = {}
#            for i, f in enumerate(values):
#                ff[i] = f
#            values = ff
#        else:
#            if grouped is None: grouped = False
#            values = copy.copy(values)
#
#        # Loop on datasets
#        ret = ()
#        for iset in sorted(values.keys()):
#            val = values[iset]
#            if not grouped:
#                val = self[iset].unmap(val)
#            ret += val,
#
#        if self.tree_type<2:
#            return ret[0]
#        return ret




#    def _input_raw_(self, input, alt):
#        if input is None:
#            return alt
#        if isinstance(input, dict): # Explicit dictionary
#            for i in xrange(self.nd):
#                if not input.has_key(i) and alt.has_key(i):
#                    input[i] = alt[i]
#            return input
#        out = {}
#        if isinstance(input, (list, tuple)): # Sequence
#            for i, data in enumerate(input):
#                out[i] = data
#        else: # Array (=> same for all dataset)
#            for i in xrange(self.nd):
#                out[i] = input
#        return out


    def _raw_rec_(self, raw_eof, raw_pc, imodes=None, ev=None):
        """Generic raw reconstruction of modes for pure PCA, MSSA or SVD, according to EOFs and PCs, for ONE DATASET

        raw_eof: (nspace,nmode)
        raw_pc: (nt,nmode)
        """

#        # Get EOFs and PCs for one dataset
#        if isinstance(raw_eof, (list, tuple, dict)): raw_eof = raw_eof[iset]
#        if isinstance(raw_pc, (list, tuple, dict)): raw_pc = raw_pc[iset]
#
        # Sizes
        ns = raw_eof.shape[0]
        nt = raw_pc.shape[0]
        if raw_eof.ndim==3:
            nw = raw_eof.shape[1]
            raw_eof = raw_eof.reshape((ns*nw, -1))
            nt += nw-1
        else:
            nw = 0

        # Which modes
        nmode = raw_eof.shape[-1]
        imodes = self._get_imodes_(imodes, nmode)

        # Function of reconstruction
        if nw:
            function = _core.mssa_rec # MSSA
        else:
            function = _core.pca_rec  # PCA/SVD

        # Arguments
        if npy.ma.isMA(raw_pc):
            raw_pc = raw_pc.filled(default_missing_value)
        if npy.ma.isMA(raw_eof):
            raw_eof = raw_eof.filled(default_missing_value)
        args = [npy.asfortranarray(var) for var in [raw_eof, raw_pc]]
        if nw:
            args.extend([ns, nt, nw])

        # Loop no modes
        smodes = []
        ffrec = 0.
        for j, ims in enumerate(imodes):

            # python -> fortran index
            ims[0] += 1
            ims[1] += 1

            # check modes
            if ims[0] > nmode: break
            if ims[1] > nmode: ims = (ims[0], nmode)

            # add modes to args
            targs = args+list(ims)
            targs += [default_missing_value]
            if ev is not None:
                targs += [ev]

            # fortran call
            try:
                raw_rec, errmsg = function(*targs)
            except:
                pass
            self.check_fortran_errmsg(errmsg)
            raw_rec = npy.ascontiguousarray(raw_rec) # (nc,nt)
            ffrec += npy.ma.masked_values(raw_rec, default_missing_value, copy=False)

            # mode specs as a string
            if ims[0] == ims[1]:
                smode = str(ims[0])
            else:
                smode = '%i-%i'%tuple(ims)
            smodes.append(smode)

        return ffrec,'+'.join(smodes)


    def __iter__(self):
        return SpanlibIter(self)
    def next(self):
        if not hasattr(self, '_iter'):
            self.error('No iterator have been initialized')
        return span._iter.next()


SpAn = Analyzer


def phase_composites(data,  nphase=8, minamp=.5, firstphase=0, index=None, format=True):
    """ Phase composites for oscillatory fields

      This computes temporal phase composites of a spatio-temporal
      dataset. The dataset is expected to be oscillatory in time.
      It corresponds to a reorganisation of the time axis
      to represents the dataset over its cycle in a arbitrary
      number of phases. It is useful, for example, to have a
      synthetic view of a reconstructed MSSA oscillation.

    :Parameters:

        - *data*: A variable with a time as the first dimension.

    :Options:

        - *nphase*: Number of phases (divisions of the cycle)
        - *minamp*: Minimal value of retained data, relative to standard deviation.
        - *firstphase*: Position of the first phase in the 360 degree cycle.
        - *index*: Index to identify phases. If ``None``, the first PC is used
           (see :meth:`SpAn.pca_pc`).
        - *format: If ``data`` is not a :class:`MV2.array` (CDAT) array,
          and CDAT is supported, output is a class:`MV2.array` array.

    :Returns:
        A ``numpy.ma.array`` (masked array).
    """

    # Get the first PC and its smoothed derivative
    if index is None:
        pc = SpAn(data).pca_pc(npca=1, quiet=True, format=False)[0]
    else:
        pc = npy.array(index)
    pc /= pc.std() # normalization
    if len(pc) > 2: # 1,2,1 smooth
        pc[1:-1] = npy.convolve(pc, [.25, .5, .25], 'valid')
    dpc = npy.gradient(pc)
    if len(pc) > 2: # 1,2,1 smooth
        dpc[1:-1] = npy.convolve(dpc, [.25, .5, .25], 'valid')
    dpc[:] = dpc/dpc.std() # normalization

    # Get amplitude and phase indexes
    amplitudes = npy.hypot(pc, dpc)
    angles = npy.arctan2(dpc, pc)
    dphase = 2*npy.pi/nphase
    angles[:] = npy.where(angles >= 2*npy.pi-dphase*.5, angles-dphase*.5, angles)
    marks = dphase * (npy.arange(nphase+1) - .5) + firstphase*npy.pi/360.
    angles = (angles-marks[0])%(npy.pi*2)+marks[0]

    # Count
    indexes = npy.digitize(angles,marks)-1

    # Initialize output variable
    if format or cdms2_isVariable(data):
        if not cdms2_isVariable(data): data = MV2.asarray(data)
        order = data.getOrder()
        itaxis = npy.clip(order.find('t'), 0, data.ndim-1)
        phases = MV2.repeat(MV2.take(data, (0, ), itaxis), nphase, itaxis)
        paxis = cdms2.createAxis(npy.arange(nphase)*dphase, id='phases')
        paxis.long_name = 'Circular phases'
        paxis.units = 'degrees'
        axes = data.getAxisList()
        axes[itaxis] = paxis
        phases.setAxisList(axes)
        phases = phases
    else:
        data = npy.ma.asarray(data)
        phases = npy.resize(data[0], (nphase, )+data.shape[1:])
    phases[:] = MV2.masked

    # Loop on circular bins to make composites
    slices = [slice(None), ]*phases.ndim
    idx = npy.arange(data.shape[itaxis])
    for iphase in xrange(len(marks)-1):
        slices[itaxis] = iphase
        inbin = amplitudes > minamp
        inbin &= angles >= marks[iphase]
        inbin &= angles < marks[iphase+1]
        phases[tuple(slices)] = data.compress(inbin, axis=itaxis).mean(axis=itaxis)

    return phases



class RedNoise(object):
    """Create a red noise generated based on lag-0 and lag-1 autocovariances of a variable

    :Algorithmm: The algorithm method follows [Allen_and_Smith_1996] use the following formula.

        - The red noise is an autoregressive process of order 1 (AR(1)):

          .. math:: u_t = u_0 + \gamma(u_{t-1}-u_0) + \alpha z_t

          where :math:`z` is a white noise of unit variance.

        - Its Lag-l autocovariance is:

          .. math:: c_l = \frac{\alpha^2\gamma^l}{1-\gamma^2}

        - Estimator of the lag-l autocovariance of a sample:

          .. math:: \hat{c}_l = \frac{1}{N-l} \sum_{i=1}^{N-l} (d_i-\bar{d})(d_{i+l}-\bar{d})

        - Bias of the variance of surrogates generated using :math:`\gamma` and :math:`\alpha`:

          .. math::

            \mu^2(\gamma) = \frac{1}{N} + \frac{1}{N^2} \left[ \frac{N-\gamma^N}{1-\gamma}-
            \frac{\gamma(1-\gamma^{N-1})}{(1-\gamma)^2} \right]

        - Corrected estimator of the lag-l covariance:

          .. math:: \tilde{c}_l \equiv \hat{c}_l + \hat{c}_0 \gamma^2

        - Corrected :math:`\tilde{\gamma}` is the solution of :

          .. math::

            \frac{\hat{c}_1}{\hat{c}_0}  =
            \frac{\tilde{\gamma}-\mu^2(\tilde{\gamma})}{1-\mu^2(\tilde{\gamma})}

          using a Newton-Raphson algorithm.

        - We now have : :math:`\tilde{c}_0 = \hat{c}_0/(1-\mu^2(\tilde{\gamma}))` .
        - :math:`\tilde{\alpha}` is estimated using the second equation.
        - The red noise is then generated using the first formula above.

    :Parameters: *data*: array with time as first dimension

    :Example:

    >>> noiser = RedNoise(data)
    >>> noise1 = noiser.sample()
    >>> noise2 = noiser.sample()
    """

    # Bias
    mu2 = classmethod(lambda cls,  g, N: -1./N + 2./N**2 * ( (N-g**N)/(1-g) + (g**N-g)/(1-g)**2 ))

    # Derivative of the bias with respect to gamma
    dmu2dg = classmethod(lambda cls, g, N: 2.*( -(N+1)*g**N + (N-1)*g**(N+1) + (N+1)*g -N+1 ) / ( N**2*(g-1)**3 ))

    # Newtown iterations
    nitermax = 6

    # Unbiased gamma
    ubg = classmethod(lambda cls, g, N: (g-RedNoise.mu2(g, N))/(1-RedNoise.mu2(g, N)))

    # Derivative of the unbiased gamma with respect to gamma
    dubgdg = classmethod(lambda cls, g, N: (1-RedNoise.mu2(g, N)-RedNoise.dmu2dg(g, N)*(1-g)) / (1-RedNoise.mu2(g, N))**2)

    def __init__(self, data):

        # Get biased statistics (nchan)
        data = data-data.mean(axis=0)
        nt = data.shape[0]
        c0 = data.var(axis=0, ddof=0)
        c1 = (data[:-1]*data[1:]).sum(axis=0)
        c1 /= (nt-1)
        self.shape = data.shape
        del data

        # Find the bias and gamma (Newton-Raphson algo)
        gamma0 = c1/c0
        self.gamma = c1/c0
        for i in xrange(0, self.nitermax):
            self.gamma -= (self.ubg(self.gamma, nt)-gamma0)/self.dubgdg(self.gamma, nt)
        self.bias = self.mu2(self.gamma, nt)

        # Corrections
        c0 /= (1-self.bias)

        # Setup other red noise coefficients
        alpha2 = c0
        alpha2 *= 1-self.gamma**2
        self.alpha = npy.sqrt(alpha2) ; del alpha2
        self.gg = npy.repeat(self.gamma[npy.newaxis, :], nt, axis=0) # (nt,nchan)
        self.gg = self.gg.cumprod(axis=0, out=self.gg)


    def sample(self):
        """Get a red noise sample fitted to input data"""
        white_noise = npy.random.randn(*self.shape) # (nt,nchan)
        white_noise /= white_noise.std(axis=0, ddof=1)
        white_noise *= self.alpha
        white_noise *= self.gg[::-1]
        red_noise = white_noise.cumsum(axis=0, out=white_noise)
        red_noise /= self.gg[::-1]
        red_noise -= red_noise.mean(axis=0)
        return red_noise

def freqfilter(data, low_freq, high_freq, **kwargs):
    """Filter out frequencies using FFT applied to MSSA PCs"""
    span = SpAn(data, **kwargs)
    stpcs = span.mssa_pc()
    if npy.ma.isMA(stpcs): stpcs = stpcs.filled(0.)
    nt = stpcs.shape[1]
    modes = []
    for mode, pc in enumerate(stpcs):
        result = npy.fft.rfft(pc)
        imax = npy.argmax(result)
        freq = npy.fft.fftfreq(nt)[:nt/2]
        freq_max = freq[imax]
        if freq_max < low_freq or freq_max > high_freq:
            modes.append(mode)
    return span.mssa_rec(modes=modes)

