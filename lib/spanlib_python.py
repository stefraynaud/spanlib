#################################################################################
# File: spanlib_python.py
#
# This file is part of the SpanLib library.
# Copyright (C) 2006-2011  Stephane Raynaud, Charles Doutriaux
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

import spanlib_fort,numpy
npy = numpy
try:
    import cdms2 ,MV2
    MV = MV2
    cdms = cdms2
    cdms2_isVariable = cdms2.isVariable
    has_cdat_support = True
except:
    cdms2_isVariable = lambda var: False
    has_cdat_support = False
import copy,gc
from warnings import warn
#print spanlib_fort.__file__

docs = dict(
    npca = """- *npca*: int | ``None``
                Number of PCA modes to keep in analysis (defaults to 10).""", 
    prepca = """- *prepca*: int | bool | ``None``
                Number of pre-PCA modes to keep before MSSA and SVD analyses (defaults to ``npca`` if ``True``).
                If the number of input channels is greater than 30, it is automatically switch to ``True``.""", 
    nmssa = """- *nmssa*: int | ``None``
                Number of MSSA modes to keep in analysis (defaults to 10).""", 
    window = """- *window*: int | ``None``
                Size of the MSSA window parameter (defaults to 1/3 the time length).""", 
    nsvd = """- *nsvd*: int | ``None``
                Number of SVD modes to keep in analysis (defaults to 10).""", 
    iset = """- *iset*: int | ``None``
                Set it to an integer to perform reconstruction on one dataset.
                If ``None``, all dataset are reconstructed.""", 
    modes = """- *modes*: int | list | tuple
                If ``None``, all modes are summed. 
                Example of other usages:
                
                    - ``4`` or ``[4]`` or ``(4,)``: only mode 4
                    - ``-4``: modes 1 to 4
                    - ``(1,3,-5)``: modes 1, 3, 4 and 5 (``-`` means "until")""", 
    raw = """- *raw*: bool
                When pre-PCA is used and ``raw`` is ``True``, it prevents from
                going back to physical space (expansion to PCA EOFs space).""", 
    scale = """- *scale*: bool | float
                Apply a factor to EOFs or PCs. This is essentially useful to add
                a quantitative meaning. If ``True``, ``scale`` is chosen so that
                the standard deviation of the mode (EOF or PC) is the square of
                the associated eigen value.""", 
    relative = """- *relative*: bool
                Return percentage of variance instead of its absolute value.""", 
    sum = """- *sum*: bool
                Return the sum of ALL (not only the selected modes) eigen values (total variance).""", 
    cumsum = """- *cumsum*: bool
                Return the cumulated sum of eigen values.""", 
)

def _filldocs_(func):
    func.__doc__ = func.__doc__ % docs
    return func

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

class Data(object):
    """Class to handle a single variable
    
        This class packs a variable in 2D space by removing
        all masked points and storing a space-time array.
        It performs this operation also on the weights.
        It is used for removing unnecessary points and
        simplifying the input format for analysis functions.
            
        :Parameters:
        
            - **data**: masked array
                Flatten in space an [x,y,t] array by 
                removing its masked point
        
        :Options:
        
            - **weights**: Weights to be flatten also
            - **norm**: Normalisation coefficients
            - **mask**: Integer mask where valid data = 1
    """
    def __init__(self, data, weights=None, norm=None, keep_invalids=False, nvalid=None, clean_weights=True):
        
        # Guess data type and copy
        if cdms2_isVariable(data):
            self.array_type = 'MV2'
            data = data.clone()
        elif npy.ma.isMA(data):
            self.array_type = 'numpy.ma'
            data = data.copy()
        else:
            self.array_type = 'numpy'
            data = data.copy()
        self.data = data
        self.dtype = data.dtype
        data = data.astype('d')
            
         # Shape
        self.shape = data.shape
        self.ndim = data.ndim
        self.nt = self.shape[0]
        self.nstot = data.size/self.nt
        self.nsdim = data.ndim-1
        
        # Check time axis
        if cdms2_isVariable(data) and data.getTime() is not None: 
            order = data.getOrder()
            if not order.startswith('t'):
                warn('Time axis is not the first axis of input variable (order="%s")'%order)

        # Weights ?
        if weights is None or weights is False:
            if weights is not False and data.ndim == 3 and \
                cdms2_isVariable(data) and \
                'x' in data.getOrder() and 'y' in data.getOrder():
                import cdutil
                weights = cdutil.area_weights(data[0]).data.astype('d') # Geographic weights
            elif self.nstot==1:
                weights = npy.ones(1)
            else:
                weights = npy.ones(self.shape[1:])
        elif npy.ma.isMA(weights):
            weights = weight.filled(0.).astype('d')
        else:
            weights = npy.asarray(weights, dtype='d')
        if data.ndim>1 and self.shape[1:] != weights.shape:
            raise SpanlibError('Weights must be of shape %s (instead of %s)'
                %(self.shape[1:],  weights.shape))
            
        # Store some info
        # - time
        if not cdms2_isVariable(data):
            self.taxis = data.shape[0]
        else:
            self.taxis = data.getAxis(0)
        # - others axes and attributes
        if cdms2_isVariable(data): # cdms -> ids
            self.saxes = data.getAxisList()[1:]
            self.id = data.id
            self.atts =  {}
            for att in data.listattributes():
                self.atts[att] = data.attributes[att]
            self.grid = data.getGrid()
            data = data.asma()
        else: # numpy -> length
            self.saxes = data.shape[1:]
            self.id = None
            self.atts = None
            self.grid = None
        # - missing value
        if npy.ma.isMA(data):
            self.missing_value = data.get_fill_value()
        else:
            self.missing_value = 1.e20
        # - special cases
        for att in 'long_name', 'units':
            if hasattr(data, att):
                setattr(self, att, data.attributes[att])
        

        # Mask (1 means good)
        # - real good values
        bmask = npy.ma.getmaskarray(data)
        good = 1-bmask.astype('l')
        # - first from data (integrate) => 1D
        count = npy.atleast_1d(good.sum(axis=0))
        del good
        # - now remove channels where weight is zero
        if clean_weights:
            count[npy.atleast_1d(weights==0.)] = 0
        # - check number of valid data along time
        nvalid = npy.clip(int(nvalid), 1, self.nt) if nvalid is not None else 2
        count[count<nvalid] = 0
        count = npy.clip(count, 0, 1)
        # - save as 0/1
        self.ns = long(count.sum())
        self.compress = count.size != self.ns
        self.good = count>0
        self.nvalid = nvalid
        
        # Scale unpacked data
        # - mean
        self.mean = data.mean(axis=0)
        # - normalisation factor
        if norm is True or norm is None:
            norm = self.data.std() # Standard norm
        elif norm is not False:
            if norm <0: # Relative norm, else strict norm
                norm = abs(norm)*self.data.std()
        else:
            norm = 1.
        self.norm = norm
        # - apply
        self.scale(data)
            
        # Fill data
        # - fill with mean (0.) where possible
        if nvalid != self.nt: 
            invalids = bmask & self.good
            data[invalids] = 0.
            if keep_invalids: 
                self.invalids = invalids
            else: 
                self.invalids = None
                del invalids
        else:
            self.invalids = None
        # - finally fill with missing values at zero
        if npy.ma.isMA(data):
            data_num = data.filled(0.)
        else:
            data_num = data
    
        # Pack 
        # - data
        self.packed_data = self.core_pack(data_num, force2d=True)
        # - weights
        self.packed_weights = self.core_pack(weights)
        

        
    def core_pack(self, data_num, force2d=False):
        """Compress data along space if needed and put time dimension first
        
        :Parameters:
        
            - **data_num**: Pure numpy array.
        """    
        # Remove bad channels ?
        nxdim = data_num.ndim-self.nsdim # dims other than space (2 for steofs)
        if self.compress: # Pack  
            sl = [slice(None)]*nxdim+[self.good]
            pdata = data_num[sl].T
        else: # Don't pack
            if data_num.ndim>2:
                data_num.shape = data_num.shape[:nxdim]+(-1, )
            pdata = data_num.T
            
        # At least 1D
        pdata = npy.atleast_1d(pdata)
        
        # 2D?
        if force2d and pdata.ndim==1:
            if self.nsdim==0:
                pdata = npy.atleast_2d(pdata)
            else:
                pdata.shape = pdata.size, 1
        return pdata
        
    def scale(self, data, copy=False, mean=None, norm=None):
        """Remove mean and normalize unpacked data"""
        if copy: 
            data = cdms2_isVariable(data) and data.clone() or data.copy()
        if mean is not False:
            data[:] -= self.mean if mean is True or mean is None else mean
        if norm is not False:
            data[:] /= self.norm if norm is True or norm is None else norm
        return data
        
    def rescale(self, data, copy=False, mean=None, norm=None):
        """Re-add mean and unnormalize unpacked data"""
        if copy: 
            data = cdms2_isVariable(data) and data.clone() or data.copy()
        if norm is not False:
            data[:] *= self.norm if norm is True or norm is None else norm
        if mean is not False:
            data[:] += self.mean if mean is True or mean is None else mean
        return data
   
#    def scale(self, pdata, copy=False):
#        """Demean and normalize packed data"""
#        if copy: pdata = pdata.copy()
#        if pdata.ndim==1:
#            mean = self.packed_mean[:,0]
#        else:
#            mean = npy.repeat(self.packed_mean, pdata.shape[1], axis=1)
#        pdata[:] -= mean ; del mean
#        pdata[:] /= self.norm
#        return 
        
#    def rescale(self, pdata, copy=False):
#        """Remean and unnormalize packed data"""
#        if copy: pdata = pdata.copy()
#        if pdata.ndim==1:
#            mean = self.packed_mean[:,0]
#        else:
#            mean = npy.repeat(self.packed_mean, pdata.shape[1], axis=1)
#        pdata[:] *= self.norm
#        pdata[:] += mean ; del mean
#        return pdata
   
    def repack(self, data, scale=True, force2d=False):
        """Pack a variable using previously computed mask array"""
        
        # Scale
        if scale:
            data = self.scale(data, copy=True)
        
        # Numpy
        if npy.ma.isMA(data): data = data.filled()
        
        # Check shape
        nsdim = self.ndim-1
        if nsdim!=0: 
            if data.shape[-nsdim:] != self.shape[-nsdim:]:
                raise SpanlibError('Incompatible shape of channels (%s instead of %s)'
                    %(data.shape[-nsdim:], self.shape[-nsdim:]))                
        # Pack it
        return self.core_pack(data, force2d=force2d)
       
    def _get_firstdims_(self, firstdims, firstaxes):
        shape = self.shape[1:]        
        if firstdims is not None:
            if not firstdims:
                firstdims = False
            elif not isinstance(firstdims, tuple):
                firstdims = (firstdims, )
        if firstdims is not False:
            if firstaxes is None:
                firstaxes = [self.taxis]
            else:
                if isinstance(firstaxes, tuple):
                    firstaxes = list(firstaxes)
                elif not isinstance(firstaxes, list):
                    firstaxes = [firstaxes]
            if firstdims is None and firstaxes is not None: 
                firstdims = tuple([(isinstance(a, (int, long))  and a or len(a)) for a in firstaxes])
            shape = firstdims + shape
            if firstaxes and isinstance(firstaxes[0], int): # real axes, not ints
                firstaxes = None
        elif len(shape)==0:
            shape = (1, )
        return shape, firstdims, firstaxes

    def create_array(self, firstdims=None, format=1, firstaxes=None):
        """Initialize an array similar to input array
        
        Type of array is guessed from attribute :attr:`array_type`.
        
        The array is filled with attribute :attr:`missing_value` if
        pure numpy, else with masked values.
        
        :Params:
        
            - **firstdims**, optional: Size of the first dimension(s).
              Defaults to attribute :attr:`nt`.
            - **format**, optional: To format output array as a CDAT variable
              using information from analyzed array.
              
                - ``1``: Add all initial axes but the first one.
                - ``2`` or ``"full"``: Add attributes. 
                - else, does not format.
                
            - **firstaxes**, optional: Set thes axes as the first ones.
              If ``firstaxes is None`` and ``firstdims is None``, it defaults
              to the first axis of analyzed array.
              You may also provide integers instead of axes, which are then
              used to set length of first axes if ``firstdims`` is not
              set.
              

        """
        # Get base array
        # - shape
        shape, firstdims, firstaxes = self._get_firstdims_(firstdims, firstaxes)
        # - create
        MM = eval(self.array_type)
        data = MM.zeros(shape, self.dtype)
        # - missing values
        if self.array_type != 'numpy':
            data[:] = npy.ma.masked
            data.set_fill_value(self.missing_value)
        else:
            data[:] = self.missing_value
            
        # Format CDAT variables
        if self.array_type=='MV2' and format:
            
            # Axes
            for i, axis in enumerate(self.saxes):
                data.setAxis(i+len(firstdims), axis)
            if firstdims is not False and firstaxes is not None:
                for i, a in enumerate(firstaxes):
                    if not isinstance(a, (int, long)):
                        data.setAxis(i, a)
             
            # Attributes
            if format=='full' or format>1:
                data.id = self.id
                for att, val in self.atts.items():
                    setattr(data, att, val)
        return data
        
        
    def unpack(self, pdata, rescale=True, format=1, firstdims=None, firstaxes=None):
        """Unpack data along space, reshape, and optionally unnormalize and remean.
        
        Input is sub_space:other, output is other:split_space.
        For MSSA: pdata(nchan,window,nmssa) (other = window:nmssa)
        
        :Params:
            
            - **pdata**: Packed data.
            - **rescale**, optional: Rescale the variable (mean and norm).
            - **format**, optional: Format the variable (see :meth:`create_array`).
            - **firstaxes**, optional: First axis (see :meth:`create_array`).
        """
        
            
        # Unpack
        # - space is last
        pdata = npy.ascontiguousarray(pdata.T).copy()
        # - first dimensions
        if firstdims is None and firstaxes is None:
            firstdims = 0 if pdata.ndim==1 else pdata.shape[0]
        shape, firstdims, firstaxes = self._get_firstdims_(firstdims, firstaxes) #FIXME:remove _get_firstdims_?
        # - create variable
        data = self.create_array(firstdims=firstdims,  format=format, firstaxes=firstaxes)
        # - check packed data shape
        firstdims = data.shape[:len(data.shape)-self.nsdim]
        if len(firstdims) > 1:
            pdata.shape = firstdims + (-1, )
        # - uncompress
        first_slices = (slice(None), )*max(1, len(firstdims))
        if self.compress:
            mdata = data.asma() if self.array_type=='MV2' else data
            slices = first_slices+(self.good, )
            mdata[slices] = pdata
            data[:] = mdata # just to be sure
        else:
            if self.nsdim==0 and pdata.ndim>len(firstdims): pdata = pdata[first_slices+(0, )]
            data[:] = pdata
        del pdata
        
        # Rescale
        if rescale:
            self.rescale(data)

        return data

    def has_cdat(self):
        """Was input array of CDAT type (:mod:`MV2`)?"""
        return self.array_type == "MV2"

          
class Dataset(object):
    """Class to handle one or a list of variables
    
    This fonction concatenates several dataset that have the
    same time axis. It is useful for analysing for example
    several variables at the same time.
    It takes into account weights, masks and axes.

    :Params:
    
        - *dataset*: A list of data objects.
            They must all have the same time length.
            
    :Options:
    
        - *weights*: Associated weights.    
    """
    
    def __init__(self, dataset, weights=None, norms=None, 
        keep_invalids=False, nvalid=None, clean_weights=True):
        
        if isinstance(dataset, (list, tuple)):
            dataset = list(dataset)
            self.map = len(dataset)
        else:
            dataset = [dataset]
            self.map = 0
        self.ndataset = len(dataset)

        # Inits
        self.data = []
        self.nt = None
        weights = self.remap(weights, reshape=True)
        norms = self.remap(norms, reshape=True)
        self.invalids = []
    
        # Loop on datasets
        for idata,data in enumerate(dataset):
    
            # Create the Data instance and pack array
            dd = Data(data, norm=norms[idata], weights=weights[idata], 
                keep_invalids=keep_invalids, nvalid=nvalid, clean_weights=clean_weights)
            self.data.append(dd)
            self.invalids.append(dd.invalids)
            
            # Check nt
            if self.nt is None:
                self.nt = dd.nt
            elif self.nt != dd.nt:
                raise SpanlibError('Time dimension of variable %i must have length %i (not %i)'%(idata, self.nt, dd.nt))
    
        # Merge
        self.stacked_data = npy.asfortranarray(npy.vstack([d.packed_data for d in self.data]))
        self.splits = npy.cumsum([d.packed_data.shape[0] for d in self.data[:-1]])
        self.stacked_weights = npy.hstack([d.packed_weights for d in self.data])
        self.ns = self.stacked_data.shape[0]
        self.norms = [d.norm for d in self.data]
    
    def restack(self, dataset, scale=True):
        """Stack new variables as a fortran array"""
        
        # Check input data
        if isinstance(dataset, (list, tuple)):
            dataset = list(dataset)
            dmap = len(dataset)
        else:
            dataset = [dataset]
            dmap = 0
        if  len(dataset)!=self.ndataset:
            raise SpanlibError('You must provide %i variable(s) to stack'%len(dataset))
            
        # Pack
        packs = [self[idata].repack(data, scale=scale, force2d=True)
            for idata, data in enumerate(dataset)]
            
        # Check time length (first axis)
        nt1 = dataset[0].size/self[0].nstot
        if len(dataset)>1:
            for i, d in enumerate(dataset[1:]):
                i += 1
                nt2 = d.size/self[i].nstot
                if packs[0].ndim!= packs[i].ndim:
                    if packs[0].ndim==1:
                        SpanlibError('Variable %i must not have a time dimension '% (i+1, ))
                    else:
                        SpanlibError('Variable %i must have a time dimension '% (i+1, ))
                elif nt1!=nt2:
                    raise SpanlibError('Time length of variable %i (%i) different from that of first variable (%i)'
                        % (i+1, nt2, nt1))
            
        # Stack
        stacker = npy.vstack if packs[0].ndim==2 else npy.hstack
        return npy.asfortranarray(stacker(packs))
    
 
    def unstack(self, sdata, rescale=True, format=1, firstdims=None, firstaxes=None):
        """Unstack and unpack data
        
        :Params:
        
            - **sdata**: Stacked data (as computed by :meth:`__init__` or :meth:`restack`).
            - **rescale**, optional: Rescale the variable (mean and norm).
            - **format**, optional: Format the variable (see :meth:`Data.create_array`).
            - **firstaxes**, optional: First axis (see :meth:`Data.create_array`).
            
        :Seel also: :meth:`Data.unpack`
        """
        
        # Unstack
        spliter = npy.hsplit if sdata.ndim==1 else npy.vsplit
        packs = spliter(sdata, self.splits)
        
        # Unpack
        return [self[i].unpack(pdata, rescale=rescale, format=format, firstdims=firstdims, firstaxes=firstaxes)
            for i, pdata in enumerate(packs)]
        
    
#    def get_time(self, idata=0):
#        """Get time axis"""
#        return self.data[idata].get_time()

    def __len__(self):
        return self.ndataset
        
    def __getitem__(self, i):
        return self.data[i]

    def unmap(self, out):
        """Remap out so as to be in the same form as input data"""
        if self.map==0:
            return out[0]
        return out

    def remap(self, values, reshape=True, fill_value=None):
        """Makes sure that values is a list (or tuple) of length :attr:`ndataset`
        
        :func:`broadcast` is called when reshaping.
        """ 
        # We always need a list
        if not isinstance(values, (list, tuple)):
            values = [values]
            
        # Check length
        n = len(self)
        if len(values)!=n:
            if not reshape:
                raise SpanlibError('Wrong number of input items (%i instead of %i)'
                    %(len(values), n))
            values = broadcast(values, n, fill_value)
        return values

    def has_cdat(self):
        """Check if there were at least one input array of CDAT type (:mod:`MV2`)?
        
        :Sea also: :meth:`Data.has_cdat`
        """
        return len(filter(lambda d: d.has_cdat(), self.data))

class SpAn(object):
    
    _npca_default = 10
    _npca_max = 30
    _nmssa_default = _nsvd_default = 8
    _nmssa_max = _nsvd_max = 20
    _window_default = 1/3. # Relative to time length
    _pca_params = ['npca', 'prepca']
    _mssa_params = ['nmssa', 'window']
    _svd_params = ['nsvd']
    _params = dict(pca=_pca_params, mssa=_pca_params+_mssa_params, svd=_pca_params+_svd_params)
    _all_params = _pca_params+_mssa_params+_svd_params

    def __init__(self, datasets, sequential=False, weights=None, norms=True, quiet=False, 
        nvalid=None, clean_weights=True, keep_invalids=False, **kwargs):
        """ Prepare the Spectral Analysis Object

        Description:::
          This function creates an object for future analyses.
          It optionally initializes some parameters.
        :::

        Usage:::
        analysis_object = SpAn(datasets,weights=None,npca=None,window=None,nmssa=None)

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
          sequential :: If we have a list (or tuple) of variables as "datasets", they are analysed independantly (sequential analyses in opposition to parallel analyses) if sequential is True [default: False]. If False, they are packed before being analysed.
        :::

        Output:::
          analysis_object :: SpAn object created for further analysis
        :::
        """

        self._quiet = quiet
        self.clean()
        
        # Form data, weights and norms for :class:`Dataset` 
        self._map_(datasets, sequential)
        weights = self._remap_(weights)
        norms = self._remap_(norms)
        
        # Create Dataset instances
        self.datasets = [Dataset(d, weights=w, norms=n, 
            nvalid=nvalid, clean_weights=clean_weights, keep_invalids=keep_invalids)
            for d, w, n in zip(self._datasets, weights, norms)]
                
#        # Weights and norms in the same form as datasets
#        norms = kwargs.pop('norm', norms)
#        weights = kwargs.pop('weight', weights)
#        if norms is None:   norms = 1.
#        norms = self._check_shape_(norms, True)
#        for iset,d in enumerate(self._datasets):
#            if len(d) == 1: norms[iset] = [1.]
#        weights = self._check_shape_(weights, None)
#
#        # We stack and pack data
#        keep_invalids = kwargs.pop('keep_invalids', None)
#        nvalid = kwargs.pop('nvalid', None)
#        for iset,d in enumerate(self._datasets):
#            self._stack_(d, weights[iset], norms[iset], keep_invalids=keep_invalids, nvalid=nvalid)
#        self.ndataset = len(self._datasets)
#        self._ndata= [len(dataset) for dataset in self._datasets]


        # Check and save parameters
        self._update_(None, **kwargs)
#       self._update_(npca=npca,window=window,nmssa=nmssa,nsvd=nsvd,prepca=prepca)


    #################################################################
    ## Get datasets info
    #################################################################
    
    def __getitem__(self, i):
        return self.datasets[i]
    
    def __len__(self):
        return self.nd

    def _time_axis_(self, iset, idata=0):
        """Get the time axis of data variable of a dataset"""
        return self[iset][idata].taxis

    def _mode_axis_(self, analysis_type, isets=None):
        """Get a mode axis according to the type of modes (pca, mssa, svd)  
        
        If CDAT is not used, length of axis is returned.
        """
        if not self._mode_axes.has_key(analysis_type):
            self._mode_axes[analysis_type] = {}
        single = False
        if analysis_type == 'svd':
            isets = [0]
        elif isets is None: 
            isets = xrange(self.nd)
        elif not isinstance(isets,(list,tuple)):
            isets = [isets]
            single = True
        out = []
        for iset in isets:
            nn = getattr(self,'_n'+analysis_type)
            if isinstance(nn, list): nn = nn[iset]
            if not self[iset].has_cdat():
                    self._mode_axes[analysis_type][iset] = nn
            elif not self._mode_axes[analysis_type].has_key(iset) or \
                len(self._mode_axes[analysis_type][iset]) != nn:
                self._mode_axes[analysis_type][iset] = cdms2.createAxis(npy.arange(1,nn+1))
                self._mode_axes[analysis_type][iset].id = analysis_type+'_mode'
                self._mode_axes[analysis_type][iset].long_name = analysis_type.upper()+' modes in decreasing order'
                if analysis_type!='svd':
                    self._check_dataset_tag_('_mode_axes',iset,analysis_type, svd=analysis_type=='svd')
            if analysis_type == 'svd': return self._mode_axes[analysis_type][iset]
            out.append(self._mode_axes[analysis_type][iset])
        if single: return out[0]
        return out

    def _channel_axis_(self,iset, name, **kwargs):
        """Get the channel axis for one dataset (MSSA or SVD)  
        
        If CDAT is not used, length of axis is returned.
        """
        if not self._prepca[iset]:
            nchan = self[iset].ns
        else:
            nchan = self._prepca[iset]
        channel_axes = getattr(self,'_%s_channel_axes'%name)
        if not self[iset].has_cdat():
            channel_axes[iset] = nchan
        if not channel_axes.has_key(iset) or len(channel_axes[iset]) != nchan:
            channel_axes[iset] = cdms2.createAxis(npy.arange(nchan))
            channel_axes[iset].id = '%s_channel' % name
            channel_axes[iset].long_name = '%s channels'% name.upper()
            self._check_dataset_tag_('_%s_channel_axes'% name, iset, **kwargs)
        return channel_axes[iset]

    def _svd_channel_axis_(self,iset):
        """Get the SVD channel axis for one dataset"""
        return _channel_axis_(iset, 'svd', svd=True)

    def _mssa_channel_axis_(self,iset):
        """Get the MSSA channel axis for one dataset"""
        return _channel_axis_(iset, 'mssa')

    def _mssa_window_axis_(self,iset, update=False):
        """Get the MSSA window axis for one dataset  
        
        If CDAT is not used, length of axis is returned.
        """
        if not self[iset].has_cdat():
            self._mssa_window_axes[iset] = self._window[iset]
        elif not self._mssa_window_axes.has_key(iset) or \
            len(self._mssa_window_axes[iset]) != self._window[iset]:
            self._mssa_window_axes[iset] = cdms2.createAxis(npy.arange(self._window[iset]))
            self._mssa_window_axes[iset].id = 'mssa_window'
            self._mssa_window_axes[iset].long_name = 'MSSA window time'
            self._check_dataset_tag_('_mssa_window_axes',iset)
        return self._mssa_window_axes[iset]

    def _mssa_pctime_axis_(self, iset, idata=0):
        """Get the MSSA PCs time axis for one dataset
        
        If CDAT is not used, length of axis is returned.
        """
        nt = self[iset].nt - self._window[iset] + 1
        if not self[iset].has_cdat():
            self._mssa_pctime_axes[iset] = nt
        elif not self._mssa_pctime_axes.has_key(iset) or \
            len(self._mssa_pctime_axes[iset]) != nt:
            self._mssa_pctime_axes[iset] = cdms2.createAxis(npy.arange(nt).astype('d'))
            self._mssa_pctime_axes[iset].id = 'mssa_pctime'
            self._mssa_pctime_axes[iset].long_name = 'MSSA PC time'
            self._check_dataset_tag_('_mssa_pctime_axes',iset)
            taxis = self._time_axis_(iset,idata)
            if hasattr(taxis,'units') and taxis.units.split()[0].lower() in \
                ['seconds','minutes','hours','days','months','years']:
                self._mssa_pctime_axes[iset].units = taxis.units.split()[0].lower() + ' since 0001-01-01'
                self._mssa_pctime_axes[iset].designateTime()
        return self._mssa_pctime_axes[iset]
        
    def _check_dataset_tag_(self, name, iset, key=None, long_name=True, id=True, svd=False):
        """Mark some attributes as specific to a dataset (only if there are more then one dataset)
            iset:: ID of the dataset
            key:: A dictionary key to select the dataset [default: None]
            long_name:: Mark the long name [defualt: True]
            id:: Mark the id [default: True]
            svd: Mark using 'left' or 'right'
        """
        if self.nd > 1:
            targetset = getattr(self, name)
            if key is not None:
                targetset = targetset[key]
            target = targetset[iset]
            if not cdms2_isVariable(target): return
            if svd:
                svdtag = ['left', 'right'][iset]
            if id: 
                if svd:
                    target.id  = '%s_%s'%(svdtag, target.id)
                else:
                    target.id += '_set%i'%iset
            if long_name:
                if svd:
                    target.long_name += ' for %s dataset'%svdtag
                else:
                    target.long_name += ' for dataset #%i'%iset
                    
    def _cdat_ev_(self, iset,  ev, analysis_type, relative, cumsum, 
        id=None, long_name=None, standard_name=None,  atts=None):
        """Format eigenvalues as CDAT variable"""
        if not self[iset].has_cdat(): return ev
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
            atts = self[iset][0].atts
            if atts.has_key('long_name'):
                ev.long_name += ' of '+atts['long_name']
        else:
            ev.long_name = long_name
        ev.setAxisList([self._mode_axis_(analysis_type.lower(), iset)])
        if standard_name is not False:
            if standard_name is None:
                ev.standard_name = 'eigen_values_of_'+analysis_type.lower()
            else:
                ev.standard_name = standard_name
        if relative:
            ev.units = '% of total variance'
        return ev
                    


    def norms(self, iset, idata=None):
        """Get the normalization factor of a dataset"""
        if idata is None:
            return [d.norm for d in self[iset].data]
        return self[iset][idata].norm
        
    def _cdat_inside_(self, iset, idata=None):
        """Check if a data var or a dataset has CDAT support"""
        # Not at all
        if not has_cdat_support: return False
        
        # Single var
        if idata is not None:
            return self[iset][idata].has_cdat()

        # Full dataset == at least one var
        return self[iset].has_cdat()

        
    def _update_(self, analysis_type=None, **kwargs):
        """Initialize, update and check statistical paremeters.
        A value of None is converted to an optimal value.
        Analyses are re-ran if needed by checking dependencies.
        """
        # Filter parameter list according to analysis_type
        running = isinstance(analysis_type, str)
        if running:
            for param in kwargs.keys():
                if param not in self._params[analysis_type]:
                    del kwargs[param]
                
        req_params = kwargs.keys()
            
        # Initialize old values and defaults changed to False
        old = {}
        changed = {}
        init_all = [None]*self.nd
        for param in self._all_params:
            #print 'check', param
#           if not req_params: continue
            # Get old values , defaults to None and set new value
#           if kwargs.has_key(param):
            if param == 'nsvd':  # Single value for all datasets
                changed[param] = False
                old[param] = getattr(self,'_'+param,None)
                setattr(self,'_'+param,SpAn._check_length_(kwargs.pop(param, old[param]),0,None))
            else:
                changed[param] = [False]*self.nd
                old[param] = getattr(self,'_'+param,list(init_all))
                setattr(self,'_'+param,SpAn._check_length_(kwargs.pop(param, old[param]),self.nd,None))
#               print 'cur', param, getattr(self, '_'+param), changed[param]
#       if not req_params: return changed
        if not running: return changed

        # Number of PCA modes       
#       if 'npca' in req_params:
        for iset in xrange(self.nd):
            if self._npca[iset] is None:
                # Guess a value
                if self._prepca[iset] is not None:
                    self._npca[iset] = self._prepca[iset]
                elif iset:
                    self._npca[iset] = self._npca[iset-1] # From last dataset
                else:
                    self._npca[iset] = SpAn._npca_default # Default value
            if self._prepca[iset] is not None:
                self._npca[iset] = max(self._npca[iset], self._prepca[iset]) # Min
            self._npca[iset] = npy.clip(self._npca[iset], 1, 
                min(SpAn._npca_max, self[iset].ns, self[iset].nt)) # Max
            
        # Number of pre-PCA modes before MSSA and SVD
#       if 'prepca' in req_params:
        for iset in xrange(self.nd):
            if self._prepca[iset] is None: # Default: pre-PCA needed over max (for MSSA and SVD)
                self._prepca[iset] = min(self[iset].ns, self[iset].ns) > SpAn._npca_max
                if not self._quiet and self._prepca[iset] and analysis_type in ['mssa', 'svd'] :
                    print '[mssa/svd] The number of valid points of one of the datasets is greater than %i, so we perform a pre-PCA'%SpAn._npca_max
            if self._prepca[iset] is True: # Defaults to the number of PCA modes
                self._prepca[iset] = self._npca[iset]
            elif self._prepca[iset]: # Max number of prepca modes is number of points
                self._prepca[iset] = min(self._prepca[iset], self[iset].ns, self[iset].nt)
            if self._prepca[iset] == 0:
                self._prepca[iset] = False
            
        # Dependency rules between prepca and npca
        for iset in xrange(self.nd):
            if self._prepca[iset] and self._npca[iset] < self._prepca[iset]:
                if not self._quiet  and self._prepca[iset]:
                        print 'The number of pre-PCA modes (%i) for dataset #%i is lower' \
                            ' than the number of PCA modes (%i), so we adjust the latter.' \
                            % (self._prepca[iset], iset, self._npca[iset])
                self._npca[iset] = self._prepca[iset]
            
        # Window extension of MSSA
#       if 'window' in req_params:
        for iset in xrange(self.nd):
            if self._window[iset] is None: # Initialization
                self._window[iset] = int(self[iset].nt*SpAn._window_default)
            self._window[iset] = npy.clip(self._window[iset], 1, max(1,self[iset].nt))
            
        # Number of MSSA modes
        for iset in xrange(self.nd):
#           if 'nmssa' not in req_params and not changed['prepca'][iset]: continue
#           if not changed['prepca'][iset]: continue
            if self._nmssa[iset] is None: # Initialization
                # Guess a value
                if iset:
                    self._nmssa[iset] = self._nmssa[iset-1] # From last dataset
                else:
                    self._nmssa[iset] = SpAn._nmssa_default # Default value
            if self._prepca[iset]:
                nchanmax = self._prepca[iset] # Input channels are from pre-PCA
            else:
                nchanmax = self[iset].ns # Input channels are from real space
            self._nmssa[iset] = npy.clip(self._nmssa[iset],1,
                min(SpAn._nmssa_max,nchanmax*self._window[iset])) # Max
            
        # Number of SVD modes (special case)
        if self._nsvd is None: # Initialization
            self._nsvd = SpAn._nsvd_default # Default value
        for iset in xrange(self.nd): # Check values
#           if 'nsvd' not in req_params and not changed['prepca'][iset]: continue
#            if not changed['prepca'][iset]: continue #?
            if self._prepca[iset]:
                nchanmax = self._prepca[iset] # Input channels are from pre-PCA
            else:
                nchanmax = self[iset].ns # Input channels are from real space
            self._nsvd = npy.clip(self._nsvd, 1, min(SpAn._nsvd_max, nchanmax)) # Max
            
#       # Check what changed
#       for param in self._all_params:
#           changed[param] = self._changed_param_(old,param)
            
        # Re-run analyses when needed
#       if not kwargs: return changed # Just initializations (dry run to prevent not ending loop)
        changed['nsvd'] = old['nsvd'] != self._nsvd
        runsvd = False
        for iset in xrange(self.nd):
            
            # Check what changed
            for param in self._all_params:
                if param !=  'nsvd':
                    changed[param][iset] = old[param][iset] != getattr(self,'_'+param)[iset]
            
            # Analyses
            # - PCA
            if (analysis_type == 'pca' or self._prepca[iset]) and \
                (self._pca_raw_eof.has_key(iset) and changed['npca'][iset]):
                print 'Rerunning PCA'
                self.pca(iset=iset)
                    
            # - MSSA
            if analysis_type == 'mssa' and \
                (self._mssa_raw_eof.has_key(iset) and
                    (changed['nmssa'][iset] or changed['window'][iset] or 
                    (self._prepca[iset] and changed['prepca'][iset]))):
                print 'Rerunning MSSA'
                self.mssa(iset=iset)
            
            # - SVD
            if not runsvd and analysis_type == 'svd' and (changed['nsvd'] or \
                (self._svd_raw_eof.has_key(iset) and
                    (self._prepca[iset] and changed['prepca'][iset]))):
                runsvd = True
        if runsvd:
            #FIXME: MUST NOT RERUN SVD
#           print 'Rerunning SVD'
            self.svd()
                
        # Inform about which params have been modified for each dataset
        return changed

    def _check_isets_(self,iset):
        """Check if an iset is associated to a a valid dataset.
        
        It can be a list, and it is returned as a list.
        If an iset is invalid, it is removed from the output list.
        """
        if iset is None: return range(self.nd)
        if iset == 'left':
            iset = 0
        elif iset == 'right':
            iset = 1
        if iset < 0 or iset >= self.nd:
            warn('Invalid dataset id: %i. Valid id are < %i'%(iset,self.nd))
        else:
            return [iset]

    def _check_shape_(self, inputs, fillvalue):
        """Return input as datasets (tree) *shape*"""
        imap = self._input_map
        if isinstance(imap, int):
            return [SpAn._check_length_(inputs, max(1, imap), fillvalue), ]
        inputs = SpAn._check_length_(inputs,len(imap),fillvalue)
        for iset, im in enumerate(imap):
            inputs[iset] = SpAn._check_length_(inputs[iset], max(1, im), fillvalue)
        return inputs


    #################################################################
    ## PCA
    #################################################################
    @_filldocs_
    def pca(self, iset=None, **kwargs):
        """ 
        Principal Components Analysis (PCA)
        
        It is called everytime needed by :meth:`pca_eof`, :meth:`pca_pc`, :meth:`pca_ev` and :meth:`pca_rec`.
        Thus, since results are stored in cache, it not necessary call it explicitly.

        :Parameters:
            %(npca)s
            %(iset)s
        """

        # Check on which dataset to operate
        isets = self._check_isets_(iset)

        # Update params
        self._update_('pca', **kwargs)

        # Loop on datasets
        for iset in isets:
            
            # Check if old results can be used when npca is lower
            if iset in self._pca_raw_pc and \
                self._pca_raw_pc[iset].shape[-1] > self._npca[iset]:
                continue

            # Remove old results
            for att in 'raw_eof','raw_pc','raw_ev','ev_sum':
                dic = getattr(self,'_pca_'+att)
                if dic.has_key(iset): del dic[iset]

            # Compute PCA
            pdata = self[iset].stacked_data
            if pdata.shape[1] == 1: # One single channel, so result is itself
                raw_eof = npy.ones(1, dtype=pdata.dtype)
                raw_pc = pdata
                ev_sum = raw_pc.var()
                raw_ev = npy.atleast_1d(ev_sum)

            else: # Several channels
                weights = npy.asfortranarray(self[iset].stacked_weights)
                pdata = npy.asfortranarray(pdata)
                raw_eof, raw_pc, raw_ev, ev_sum = \
                    spanlib_fort.pca(pdata, self._npca[iset], weights, -1)  
                raw_pc = npy.ascontiguousarray(raw_pc)
                raw_eof = npy.ascontiguousarray(raw_eof)

            # Save results
            self._pca_raw_pc[iset] = raw_pc
            self._pca_raw_eof[iset] = raw_eof
            self._pca_raw_ev[iset] = raw_ev
            self._pca_ev_sum[iset] = ev_sum
            
            # Delete formatted variables
            for vtype in 'pc', 'eof':
                vfmt = getattr(self, '_pca_fmt_'+vtype)
                if vfmt.has_key(iset): del vfmt[iset]
            gc.collect()

        self._last_analysis_type = 'pca'


    @_filldocs_
    def pca_eof(self, iset=None, scale=False, raw=False, unmap=True, format=True, **kwargs):
        """Get EOFs from PCA analysis

        :Parameters:
            %(scale)s
            %(raw)s
            %(iset)s
            
        :PCA parameters:
            %(npca)s
            
        :Returns:
            Arrays with shape ``(npca,...)``
        """
    
        # Dataset selection
        isets = self._check_isets_(iset)

        # Update params
        changed =  self._update_('pca', **kwargs)

        # Of, let's format the variables
        fmt_eof = {}
        for iset in isets:
            
            # Operate only on selected datasets
            if isets is not None and iset not in isets: continue
            
            # EOF already available 
            if not raw and iset in self._pca_fmt_eof:
                fmt_eof[iset] = self._pca_fmt_eof[iset]
                continue
                
            # First PCA analysis?
            if not iset in self._pca_raw_eof: self.pca(iset=iset)
                
            # Get raw data back to physical space
            if raw:
                fmt_eof[iset] = self._pca_raw_eof[iset].T
                continue
            self._pca_fmt_eof[iset] = self[iset].unstack(self._pca_raw_eof[iset], rescale=False, 
                firstaxes=self._mode_axis_('pca', iset))
            
            # Set attributes and scale
            for idata,eof in enumerate(self._pca_fmt_eof[iset]):
                
                # Attributes
                if format and cdms2_isVariable(eof):
                    if not self[iset][idata].id.startswith('variable_'):
                        eof.id = self[iset][idata].id+'_pca_eof'
                    else:
                        eof.id = 'pca_eof'
                    eof.name = eof.id
                    eof.standard_name = 'empirical_orthogonal_functions_of_pca'
                    eof.long_name = 'PCA empirical orthogonal functions'
                    atts = self[iset][idata].atts
                    if atts.has_key('long_name'):
                        eof.long_name += ' of '+atts['long_name']
                    print 'scale', scale
                    if scale and atts.has_key('units'):
                        eof.units = atts['units']
                    
                # Scaling
                if scale:
                    if scale is True: # Std dev of EOF is sqrt(ev)
                        scale = npy.sqrt(self._pca_raw_ev[iset]*(self.ns(iset)-1))
                        for imode in xrange(eof.shape[0]):
                            eof[imode] *= scale[imode]
                    else:
                        eof *= scale
                    
            fmt_eof[iset] = self._pca_fmt_eof[iset]

        if not unmap: return fmt_eof
        return self._unmap_(fmt_eof, grouped=raw)       


    @_filldocs_
    def pca_pc(self, iset=None, scale=False, raw=False, unmap=True, format=True, **kwargs):
        """Get principal components (PCs) from current PCA decomposition
        
        :Parameters:
            %(iset)s

        :PCA parameters:
            %(npca)s
            
        :Returns:
            Arrays with the shape ``(npca,nt)``
        """
        # Check on which dataset to operate
        isets = self._check_isets_(iset)

        # Update params
        self._update_('pca', **kwargs)
        
        # Of, let's format the variable
        fmt_pc = {}
        for iset in isets:
            
            # PC already available 
            if iset in self._pca_fmt_pc:
                fmt_pc[iset] = self._pca_fmt_pc[iset]
                continue
            
            # First PCA analysis?
            if not iset in self._pca_raw_pc: self.pca(iset=iset)

            # Raw ?
            if raw:
                fmt_pc[iset] = self._pca_raw_pc[iset][:,:self._npca[iset]]
                continue

            # Format the variable
            pc = npy.ascontiguousarray(self._pca_raw_pc[iset][:,:self._npca[iset]].T)
            #TODO: scale pca_pc
            if format and self[iset].has_cdat():
                pc = cdms2.createVariable(pc)
                pc.setAxis(0, self._mode_axis_('pca',iset))
                pc.setAxis(1, self._time_axis_(iset,0))
                pc.id = pc.name = 'pca_pc'
                pc.standard_name = 'principal_components_of_pca'
                pc.long_name = 'PCA principal components of '
                atts = self[iset][0].atts
                if len(self[iset])==1 and  atts.has_key('long_name'): 
                    pc.long_name += atts['long_name']
    #           else:
    #               pc.long_name += 'dataset %i'%iset
                if scale and (len(self[iset]) == 1 or npy.allclose(self[iset].norms, 1.)) and atts.has_key('units'):
                    pc.units = atts['units']
                
            
            fmt_pc[iset] = self._pca_fmt_pc[iset] = pc
            self._check_dataset_tag_('_pca_fmt_pc',iset)

        if not unmap: return fmt_pc
        return self._unmap_(fmt_pc, grouped=True)       

    def pca_ec(self, xeof=None, xdata=None, iset=None, scale=False, 
        xraw=False, xscale=True, raw=False, unmap=True, format=True, **kwargs):
        """Get expansion coefficient using current PCA decomposition
        
        Expansion coefficients are the projection of data onto EOFs.
        By default, it uses input data and computed EOFs, and thus
        returns principal components. 
        You can bypass this default behaviour using ``xeof`` and ``xdata``
        keyword parameters.
        
        :Parameters:
            %(iset)s

        :PCA parameters:
            %(npca)s
            
        :Returns:
            Arrays with the shape ``(npca,nt)``
        """
        # Check on which dataset to operate
        isets = self._check_isets_(iset)
                
        # Remap
        xeof = self._remap_(xeof)
        xdata = self._remap_(xdata)
            
        # Loop on datasets
        fmt_ec = {}
        for iset in isets:
            
            # EOFs used for projection
            if xeof[iset] is None: # From PCA
                if not iset in self._pca_raw_eof: self.pca(iset=iset)
                raw_eof = self._pca_raw_eof[iset]
            elif xraw: # Direct use
                raw_eof = xeof[iset]
            else: # We format then use
                eofs = self[iset].remap(xeof[iset])
                raw_eof = self[iset].restack(eofs, scale=False)
            
            # Data to project on EOFs
            if xdata[iset] is None: # From input
                raw_data = self[iset].stacked_data
                ndim = raw_data.ndim
            elif xraw: # Direct use
                raw_data = xdata[iset]
                ndim = raw_data.ndim
            else: # We format then use
                data = self[iset].remap(xdata[iset])
                raw_data = self[iset].restack(data, scale=xscale)
                ndim = raw_data.ndim
                if raw_data.ndim>2:
                    raw_data = npy.reshape(raw_data, (raw_data.shape[0], -1))

            # Weights
            weights = self[iset].stacked_weights
            
            # Projection
            raw_data = npy.asfortranarray(raw_data)
            raw_eof = npy.asfortranarray(raw_eof)
            weights = npy.asfortranarray(weights)
            raw_ec = spanlib_fort.pca_getec(raw_data, raw_eof, weights)
            if raw:
                fmt_ec[iset] = raw_ec
                continue
            ec = npy.ascontiguousarray(raw_ec.T)
            
            # Format
            if not format and self._cdat_inside_(iset):
                ec = cdms2.createVariable(ec)
                ec.setAxis(0, self._mode_axis_('pca',iset))
                if ndim==2 and ec.shape[1]==self[iset].nt:
                    ec.setAxis(1, self._time_axis_(iset,0))
                ec.id = ec.name = 'pca_ec'
                ec.long_name = 'PCA expansion coefficients'
                atts = self[iset][0].atts
                if len(self[iset])==1 and atts.has_key('long_name'): 
                    ec.long_name += ' of '+atts['long_name']
    #           else:
    #               ec.long_name += 'dataset %i'%iset
                if scale and (len(self[iset])==1 or npy.allclose(self.norms(iset), 1.)) and atts.has_key('units'):
                    ec.units = atts['units']
            fmt_ec[iset] = ec
            
        if not unmap: return fmt_ec
        return self._unmap_(fmt_ec, grouped=True)       
            

    @_filldocs_
    def pca_ev(self,iset=None,relative=False,sum=False,cumsum=False,**kwargs):
        """Get eigen values from current PCA decomposition

        :Parameters:
          %(relative)s
          %(sum)s
          %(cumsum)s
          %(iset)s
        
            
        :PCA parameters:
            %(npca)s
            
        :Returns:
            Arrays with shape ``(npca,)`` or a float
        """

        # Check on which dataset to operate
        isets = self._check_isets_(iset)

        # Update params
        self._update_('pca', **kwargs)

        # Loop on dataset
        res = {}
        for iset in isets:
            
            # First PCA analysis?
            if not self._pca_raw_eof.has_key(iset): self.pca(iset=iset)
                
            # We only want the sum
            if sum:
                res[iset] = self._pca_ev_sum[iset]
                continue
                
            # Data
            ev = self._pca_raw_ev[iset][:self._npca[iset]]
            if cumsum:
                ev = raw_ev.cumsum()
            if relative: 
                ev = 100.*ev/self._pca_ev_sum[iset]

            # Format the variable
            if self._cdat_inside_(iset):
                
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
                ev.setAxisList([self._mode_axis_('pca',iset)])
                ev.standard_name = 'eigen_values_of_pca'
                atts = self[iset][0].atts
                if len(self[iset])==1 and atts.has_key('long_name'):
                    ev.long_name += atts['long_name']
                else:
                    ev.long_name += 'dataset %i'%iset
                if relative:
                    ev.units = '% of total variance'
                elif (len(self[iset])==1 or npy.allclose(self.norms(iset), 1.)) and atts.has_key('units'):
                    ev.units = atts['units']
                    for ss in ['^','**',' ']:
                        if ev.units.find(ss) != -1:
                            ev.units = '(%s)^2' % ev.units
                            break
                            
            res[iset] = ev

        return self._unmap_(res, grouped=True)      

    @_filldocs_
    def pca_rec(self, iset=None, modes=None, raw=False, xpc=None, xeof=None, xraw=False, rescale=True, format=2, **kwargs):
        """Reconstruct a set of modes from PCA decomposition

        :Parameters:
            %(modes)s
            %(raw)s
            %(iset)s
            
        :PCA parameters:
            %(npca)s
            
        :Returns:
            Arrays with the same shape as input arrays.
        """
        
        # Check on which dataset to operate
        isets = self._check_isets_(iset)
            
        # Update params
        self._update_('pca', **kwargs)
        
        # Remap alternate arrays
        xeof = self._remap_(xeof)
        xpc = self._remap_(xpc)
            
        # Loop on datasets
        pca_fmt_rec = {}
        for iset in isets:
            
            # First PCA analysis?
            if not iset in self._pca_raw_pc: self.pca(iset=iset)

            # EOFs
            if xeof[iset] is None: # From PCA
                raw_eof = self._pca_raw_eof[iset]
            elif xraw: # Direct use
                raw_eof = xeof[iset]
            else: # We format then use
                eofs = self[iset].remap(xeof[iset])
                raw_eof = self[iset].restack(eofs, scale=False)
                
            # PCs
            if xpc[iset] is None: # From PCA
                raw_pc = self._pca_raw_pc[iset]
            elif xraw: # Direct use
                raw_pc = xpc[iset]
            elif cdms2_isVariable(xpc[iset]): # formatted -> pure numpy
                raw_pc = xpc[iset].filled().T
            else:
                raw_pc = xpc[iset].T
                
            # Get raw data back to physical space
            reof = raw_eof[:,:self._npca[iset]]
            rpc = raw_pc[:,:self._npca[iset]]
            raw_rec, smodes = self._raw_rec_(reof, rpc, modes)
            if raw:
                pca_fmt_rec[iset] = raw_rec.T
            else:
                pca_fmt_rec[iset] = self[iset].unstack(raw_rec, rescale=rescale, format=format)
                pass
            del  raw_rec

            # Set attributes and scale
            for idata,rec in enumerate(pca_fmt_rec[iset]):
                if not cdms2_isVariable(rec): continue
#               rec[:] *= self._norm_(iset,idata) # Scale
                if not self[iset][idata].id.startswith('variable_'):
                    rec.id = self[iset][idata].id+'_pca_rec'
                else:
                    rec.id = 'pca_rec'
                rec.name = rec.id
#               if modes is not None:
                rec.id += smodes
#               else:
#                   rec.modes = '1-%i'%self._npca[iset]
                rec.modes = smodes
                rec.standard_name = 'recontruction_of_pca_modes'
                rec.long_name = 'Reconstruction of PCA modes: '+smodes
                atts = self[iset][idata].atts
                if atts.has_key('long_name'):
                    rec.long_name += ' of '+atts['long_name']
                if atts.has_key('units'):
                    rec.units = atts['units']
                    
        return self._unmap_(pca_fmt_rec, grouped=raw)   
    

    #################################################################
    # MSSA
    #################################################################

    @_filldocs_
    def mssa(self,iset=None, **kwargs):
        """ MultiChannel Singular Spectrum Analysis (MSSA)

        It is called everytime needed by :meth:`mssa_eof`, :meth:`mssa_pc`, :meth:`mssa_ev` and :meth:`mssa_rec`.
        Thus, since results are stored in cache, it not necessary call it explicitly.

        :Parameters:
            %(nmssa)s
            %(window)s
            %(prepca)s
            %(iset)s
        """
        # Check on which dataset to operate
        isets = self._check_isets_(iset)

        # Parameters
        self._update_('mssa', **kwargs)

        # Loop on datasets
        for iset in isets:

            # Operate only on selected datasets
            if isets is not None and iset not in isets: continue
            
            # Check if old results can be used when nmssa is lower
            if getattr(self,'_mssa_raw_pc').has_key(iset) and \
                getattr(self,'_mssa_raw_pc')[iset].shape[-1] > self._nmssa[iset]:
                continue
            
            # Remove old results
            for att in 'raw_eof','raw_pc','raw_ev','ev_sum':
                dic = getattr(self,'_mssa_'+att)
                if dic.has_key(iset): del dic[iset]

            # Compute MSSA
            if self._prepca[iset]: # Pre-PCA case
            
                # PCA
                if not self._pca_raw_pc.has_key(iset): self.pca(iset=iset)
                
                # MSSA
                pca_raw_pc = npy.asfortranarray(self._pca_raw_pc[iset][:, :self._prepca[iset]].T)
                raw_eof, raw_pc, raw_ev, ev_sum = \
                  spanlib_fort.mssa(pca_raw_pc, self._window[iset], self._nmssa[iset])
                  
            else: # Direct MSSA case
                pdata = npy.asfortranarray(self[iset].stacked_data)
                raw_eof, raw_pc, raw_ev, ev_sum = \
                  spanlib_fort.mssa(pdata, self._window[iset], self._nmssa[iset])

            # Save results
            self._mssa_raw_pc[iset] = npy.ascontiguousarray(raw_pc)
            self._mssa_raw_eof[iset] = npy.ascontiguousarray(raw_eof)
            self._mssa_raw_ev[iset] = raw_ev
            self._mssa_ev_sum[iset] = ev_sum

            # Delete formmated variables
            for vtype in 'pc', 'eof':
                vfmt = getattr(self, '_mssa_fmt_'+vtype)
                if vfmt.has_key(iset): del vfmt[iset]
                
        self._last_analysis_type = 'mssa'
        gc.collect()

    @_filldocs_
    def mssa_eof(self, iset=None, scale=False, raw=False, format=True, unmap=True, **kwargs):
        """Get EOFs from MSSA analysis

        :Parameters:
            %(scale)s
            %(raw)s
            %(iset)s
            
        :MSSA parameters:
            %(nmssa)s
            %(window)s
            %(prepca)s
            
        :Returns:
            Arrays with shape ``(nmssa,nt-window+1,...)``
        """

        # Dataset selection
        isets = self._check_isets_(iset)

        # Update params
        self._update_('mssa', **kwargs)

        # Of, let's format the variable
        fmt_eof = {}
        for iset in isets: # (window*nchan,nmssa)
        
            # EOF already available 
            if not raw and iset in self._mssa_fmt_eof:
                fmt_eof[iset] = self._mssa_fmt_eof[iset]
                continue
                
            # No analyses performed?
            if not iset in self._mssa_raw_eof: self.mssa(iset=iset)
            
            # Raw eof
            raw_eof = self._mssa_raw_eof[iset][:, :self._nmssa[iset]]
            nlxnw, nm = raw_eof.shape
            nw = self._window[iset]
            nl = nlxnw/nw
            raw_eof = raw_eof.reshape((nl, nw, nm))
            
            if raw: # Do not go back to physical space
                self._mssa_fmt_eof[iset] = [npy.ascontiguousarray(raw_eof.T)]
                if self[iset].has_cdat():
                    self._mssa_fmt_eof[iset][0] = cdms2.createVariable(self._mssa_fmt_eof[iset])
                    self._mssa_fmt_eof[iset][0].setAxisList(
                        [self._mode_axis_('mssa',iset),self._mssa_channel_axis_(iset)])
                    
            else: # Get raw data back to physical space
                raw_eof = npy.ascontiguousarray(raw_eof) # (nl, nw, nm)
                firstaxes = (self._mode_axis_('mssa', iset), self._mssa_window_axis_(iset))
                if not self._prepca[iset]: # No pre-PCA performed
                    self._mssa_fmt_eof[iset] = self[iset].unstack(raw_eof, rescale=False, 
                        format=format, firstaxes=firstaxes)
                        
                else: # With pre-PCA
                    proj_eof, smodes = self._raw_rec_(self._pca_raw_eof[iset], 
                        raw_eof.T.reshape((nm*nw, nl)))
                    self._mssa_fmt_eof[iset] = self[iset].unstack(proj_eof, 
                        rescale=False, format=format, firstaxes = firstaxes)
                    
            # Set attributes
            for idata, eof in enumerate(self._mssa_fmt_eof[iset]):
                
                # Attributes
                if format and cdms2_isVariable(eof):
                    if not raw and not self[iset][idata].id.find('variable_'):
                        eof.id = self[iset][idata].id+'_mssa_eof'
                    else:
                        eof.id = 'mssa_eof'
                    eof.name = eof.id
                    eof.standard_name = 'empirical_orthogonal_functions_of_mssa'
                    eof.long_name = 'MSSA empirical orthogonal functions'
                    if not raw:
                        atts = self[iset][idata].atts
                        if atts.has_key('long_name'):
                            eof.long_name += ' of '+atts['long_name']
                        if scale and atts.has_key('units'):
                            eof.units = atts['units']
                    
                # Scaling
                if scale:
                    if scale is True:
                        scale = npy.sqrt(self._mssa_raw_ev[iset])*nl*nw
                    eof[:] *= scale
                    
            fmt_eof[iset] = self._mssa_fmt_eof[iset]
            
        gc.collect()
        if not unmap: return fmt_eof
        return self._unmap_(fmt_eof, grouped=raw)

    @_filldocs_
    def mssa_pc(self, iset=None, raw=False, unmap=True, format=True, **kwargs):
        """Get PCs from MSSA analysis
        
        :Parameters:
            %(iset)s

        :MSSA parameters:
            %(nmssa)s
            %(window)s
            %(prepca)s
            
        :Returns:
            Arrays with the shape ``(nmssa,nt)``
        """

        # Dataset selection
        isets = self._check_isets_(iset)

        # Update params
        self._update_('mssa', **kwargs)

        # Of, let's format the variable
        fmt_pc = {}
        for iset in isets:
            
            # PC already available 
            if not raw and self._mssa_fmt_pc.has_key(iset):
                fmt_pc[iset] = self._mssa_fmt_pc[iset]
                continue
                
            # No analyses performed?
            if not iset in self._mssa_raw_pc: self.mssa(iset=iset)
            
            # Raw?
            raw_pc = self._mssa_raw_pc[iset][:,:self._nmssa[iset]]
            if raw:
                fmt_pc[iset] = raw_pc
                        
            # Format the variable
            pc = npy.ascontiguousarray(raw_pc.T)
            if format and self._cdat_inside_(iset):
                
                pc = cdms2.createVariable(pc)
                pc.setAxis(0,self._mode_axis_('mssa',iset))
                pc.setAxis(1,self._mssa_pctime_axis_(iset))
                pc.id = pc.name = 'mssa_pc'
                pc.standard_name = 'principal_components_of_mssa of '
                pc.long_name = 'MSSA principal components'
                atts = self[iset][0].atts
                if len(self[iset])==1 and atts.has_key('long_name'): 
                    pc.long_name += atts['long_name']
    #           else:
    #               pc.long_name += 'of dataset %i'%iset
                #if (len(self[iset])==1 or npy.allclose(self.norms(iset), 1.)) and atts.has_key('units'):
                #   pc.units = atts['units']
            fmt_pc[iset] = self._mssa_fmt_pc[iset] = pc
            self._check_dataset_tag_('_mssa_fmt_pc',iset)

        if not unmap: return fmt_pc
        return self._unmap_(fmt_pc, grouped=True)
            
    def mssa_ec(self, iset=None, xdata=None, xeof=None, xraw=False, 
        raw=False, unmap=True, format=True, **kwargs):
        """Get expansion coefficients from MSSA analysis
        
        :Parameters:
            %(iset)s

        :MSSA parameters:
            %(nmssa)s
            %(window)s
            %(prepca)s
            
        :Returns:
            Arrays with the shape ``(nmssa,nt)``
        """

        # Dataset selection
        isets = self._check_isets_(iset)

        # Update params
        self._update_('mssa', **kwargs)

        # Remap
        xeof = self._remap_(xeof)
        xdata = self._remap_(xdata)
            
        fmt_ec = {}
        for iset in isets:
            
            # ST-EOFs used for projection
            if xeof[iset] is None:
                if not iset in self._mssa_raw_eof: self.mssa(iset=iset)
                raw_eof = self._mssa_raw_eof[iset]
            elif xraw:
                raw_eof = xeof[iset]
            elif self._prepca[iset]: # After PCA
                raw_eof = self.pca_ec(iset=iset, xdata=xeof, xscale=False, raw=True, unmap=False)[iset].T
                raw_eof.shape = -1, self._nmssa[iset]
            else:
                eofs = self[iset].remap(xeof[iset])
                raw_eof = npy.ascontiguousarray(self[iset].restack(eofs, scale=False))
                raw_eof.shape = -1, raw_eof.shape[-1] # (nl*nw, nm)
                del eofs
                
            # Data to project on ST-EOFs
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
            raw_eof = npy.asfortranarray(raw_eof[:, :self._nmssa[iset]])
            raw_ec = spanlib_fort.mssa_getec(raw_data, raw_eof, self.window(iset))
            if raw:
                fmt_ec[iset] = raw_ec
                continue
            ec = npy.ascontiguousarray(raw_ec.T)
            
            # Format the variable
            if format and self[iset].has_cdat():
                
                ec = cdms2.createVariable(ec)
                ec.setAxis(0,self._mode_axis_('mssa',iset))
#                ec.setAxis(1,self._mssa_pctime_axis_(iset))
                ec.id = ec.name = 'mssa_ec'
#                ec.standard_name = 'expansion coefficient_of_mssa of '
                ec.long_name = 'MSSA principal components'
                atts = self[iset][0].atts
                if len(self[iset])==1 and atts.has_key('long_name'): 
                    ec.long_name += atts['long_name']
    #           else:
    #               ec.long_name += 'of dataset %i'%iset
                #if (len(self[iset])==1 or npy.allclose(self.norms(iset), 1.)) and atts.has_key('units'):
                #   ec.units = atts['units']
            fmt_ec[iset] = ec
            # FIXME: _check_dataset_tag_
#            self._check_dataset_tag_('_mssa_fmt_ec', iset)

        if not unmap: return fmt_ec
        return self._unmap_(fmt_ec, grouped=True)
#            


    @_filldocs_
    def mssa_ev(self, iset=None, relative=False, sum=False, cumsum=False, 
        mctest=False, mcnens=100, mcqt=90, format=True, unmap=True, **kwargs):
        """Get eigen values from current MSSA decomposition

        :Options:
        
          %(relative)s
          %(sum)s
          %(cumsum)s
          %(iset)s
        
            
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

        # Check on which dataset to operate
        isets = self._check_isets_(iset)

        # Update params
        self._update_('mssa', **kwargs)

        # Loop on dataset
        res = {}
        for iset in isets:
            
            # No analyses performed?
            if not iset in self._mssa_raw_eof: self.mssa(iset=iset)

            # We only want the sum
            if sum:
                res[iset] = self._mssa_ev_sum[iset]
                continue
                
            # Data
            ev = [self._mssa_raw_ev[iset][:self._nmssa[iset]], ]
            
            # Mont-Carlo test
            if mctest:
                
                # Get reference data
                if self._prepca[iset]:
                    data = self._pca_raw_pc[iset]
                else:
                    data = self[iset].stacked_data
                
                # Inits
                rn = RedNoise(data.T) # red noise generator
                nmssa = self._nmssa[iset]
                mcev = npy.zeros((mcnens, nmssa))
                
                # Generate and ensemble of surrogate data
                for iens in xrange(mcnens):
                    
                    # Create a sample red noise (nt,nchan)
                    red_noise = rn.sample().T
                    
                    # Block-covariance matrix (nchan*nwindow,nchan*nwindow,)
                    cov = spanlib_fort.stcov(npy.asfortranarray(red_noise), self.window(iset))
                    cov = npy.ascontiguousarray(cov)
                    del red_noise
                    
                    # Fake eigen values (EOFt.COV.EOF)
                    ce = npy.dot(cov, self._mssa_raw_eof[iset]) #; del cov
                    evmat = npy.dot(self._mssa_raw_eof[iset].T, ce) #; del ce
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
                ev = [(100.*e/self._mssa_ev_sum[iset]) for e in ev]

            # Format the variables
            if format and self[iset].has_cdat():
                
                for i, e in enumerate(ev):
                    ev[i] = self._cdat_ev_(iset, e,  'mssa', relative, cumsum)
                if mctest:
                    ev[1].id += '_mclow'
                    ev[1].long_name += ': lower bound of MC confidence interval (%g%%)'%mcqt 
                    ev[2].id += '_mchigh'
                    ev[2].long_name += ': upper bound of MC confidence interval (%g%%)'%(100-mcqt)
                    
            res[iset] = ev[0] if not mctest else ev

        if not unmap: return res
        return self._unmap_(res, grouped=True)      
            


    @_filldocs_
    def mssa_rec(self, iset=None, modes=None, raw=False, xpc=None, xeof=None, xraw=False, 
        phases=False, rescale=True, format=2, unmap=True, **kwargs):
        """Reconstruction of MSSA modes
        
        :Parameters:
            %(iset)s
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
        
        # Dataset selection
        isets = self._check_isets_(iset)

        # Update params
        self._update_('mssa', **kwargs)
        if isinstance(phases, int):
            self._nphase = phases
        elif phases is not False:
            phases = self._nphase
        
        # Remap alternate arrays
        xeof = self._remap_(xeof)
        xpc = self._remap_(xpc)
            
        # Loop on datasets
        mssa_fmt_rec = {}
        for iset in isets:

            # No analyses performed?
            if not iset in self._mssa_raw_eof: self.mssa(iset=iset)
            
            # ST-EOFs used for projection
            if xeof[iset] is None:
                if not iset in self._mssa_raw_eof: self.mssa(iset=iset)
                raw_eof = self._mssa_raw_eof[iset]
            elif xraw:
                raw_eof = xeof[iset]
            elif self._prepca[iset]: # After PCA
                raw_eof = self.pca_ec(iset=iset, xdata=xeof, raw=True, xscale=False, unmap=False)[iset].T
                raw_eof.shape = -1, self._nmssa[iset]
            else:
                eofs = self[iset].remap(xeof[iset])
                raw_eof = npy.ascontiguousarray(self[iset].restack(eofs, scale=False))
                raw_eof.shape = -1, raw_eof.shape[-1] # (nl*nw, nm)
                del eofs
            raw_eof = raw_eof[:, :self._nmssa[iset]]
            nlxnw, nm = raw_eof.shape
            nw = self._window[iset]
            nl = nlxnw/nw
            raw_eof = raw_eof.reshape((nl, nw, nm))
                
            # PCs
            if xpc[iset] is None: # From PCA
                if not iset in self._mssa_raw_pc: self.mssa(iset=iset)
                raw_pc = self._mssa_raw_pc[iset]
            elif xraw: # Direct use
                raw_pc = xpc[iset]
            elif cdms2_isVariable(xpc[iset]): # formatted -> pure numpy
                raw_pc = xpc[iset].filled().T
            else:
                raw_pc = xpc[iset].T
            raw_pc = raw_pc[:, :self._nmssa[iset]]
            
            # Projection
            raw_rec, smodes = self._raw_rec_(raw_eof, raw_pc, modes)
                
            # Phases composites
            if phases:
                raw_rec_t = phase_composites(raw_rec.T, phases, format=format and self[iset].has_cdat)
                nt = raw_rec.shape[0]
                if cdms2_isVariable(raw_rec_t):
                    raw_rec = MV2.transpose(raw_rec_t)
                    taxis = raw_rec.getAxis(1)
                else:
                    raw_rec = raw_rec_t.asma().T
                    taxis = nt
                del raw_rec_t
            else:
#                taxis = self._time_axis_(iset)
#                nt = len(taxis)
                nt = self[iset].nt
                
            # Get raw data back to physical space (nchan,nt)
            if not self._prepca[iset]: # No pre-PCA performed
                mssa_fmt_rec[iset] = self[iset].unstack(raw_rec, rescale=rescale, format=format)
                
            elif raw: # Force direct result from MSSA
                if format and self[iset].has_cdat():
                    mssa_fmt_rec[iset] = [cdms2.createVariable(raw_rec.T)]
                    mssa_fmt_rec[iset][0].setAxisList(0,[taxis,self._mode_axis_('mssa',iset)])
                else:
                    mssa_fmt_rec[iset] = raw_rec.T
                    
            else: # With pre-pca
                proj_rec, spcamodes = \
                    self._raw_rec_(self._pca_raw_eof[iset], raw_rec.T)
                mssa_fmt_rec[iset] = self[iset].unstack(proj_rec, rescale=rescale, format=format)
            del  raw_rec
            
            # Remove the mean for phases
            if phases:
                for rec in mssa_fmt_rec[iset]:
                    rec[:] -= rec.mean(axis=0)
            
            # Set attributes
            if format:
                for idata,rec in enumerate(mssa_fmt_rec[iset]):
                    if not cdms2_isVariable(rec): continue
                    if not self[iset][idata].id.startswith('variable_'):
                        rec.id = self[iset][idata].id+'_mssa_rec'
                    else:
                        rec.id = 'mssa_rec'
    #               if modes is not None:
                    rec.id += smodes #FIXME: do we keep it?
                    rec.modes = smodes
    #               else:
    #                   rec.modes = '1-%i'%self._nmssa[iset]
                    rec.standard_name = 'recontruction_of_mssa_modes'
                    rec.long_name = 'Reconstruction of MSSA modes'
                    atts = self[iset][idata].atts
                    if atts.has_key('long_name'):
                        rec.long_name += ' of '+atts['long_name']
                  
        if not unmap: return mssa_fmt_rec
        return self._unmap_(mssa_fmt_rec, grouped=raw)
    
    def mssa_phases(self, pair, iset=0, nphase=8, format=True, unmap=True, **kwargs):
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
        
        # Dataset selection
        isets = self._check_isets_(iset)
        return self.mssa_rec(iset=iset, modes=pair, phases=nphase, format=format, remap=remap)

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
            %(iset)s
        """

        if self.nd<2:
            raise SpanlibError('Error you need at least (most) 2 datasets to run svd, otherwise use pca and mssa',  'svd')
        
        # Parameters
        self._update_('svd', **kwargs)

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
        raw_eof_left, raw_eof_right, raw_pc_left, raw_pc_right, raw_ev, ev_sum, info = \
            spanlib_fort.svd(left, right, self._nsvd, lweights, rweights, int(usecorr))
        if info != 0:
            raise SpanlibError('Error when running fortran SVD',  'svd')
            
        # Save results
        self._svd_raw_pc[0] = npy.ascontiguousarray(raw_pc_left)
        self._svd_raw_eof[0] = npy.ascontiguousarray(raw_eof_left)
        self._svd_raw_pc[1] = npy.ascontiguousarray(raw_pc_right)
        self._svd_raw_eof[1] = npy.ascontiguousarray(raw_eof_right)
        self._svd_raw_ev = raw_ev
        self._svd_ev_sum = ev_sum


        self._last_analysis_type = 'svd'
        gc.collect()

    @_filldocs_
    def svd_eof(self,iset=None,scale=False, raw=False, unmap=True, format=True, **kwargs):
        """Get EOFs from SVD analysis

        If SVD was not performed, it is done with all parameters sent to :meth:`svd`

        :Parameters:
            %(scale)s
            %(raw)s
            %(iset)s
            
        :SVD parameters:
            %(nsvd)s
            %(prepca)s
            
        :Returns:
            Arrays with shape ``(nsvd,...)``
        """

        # Dataset selection
        isets = self._check_isets_(iset)

        # Update params
        self._update_('svd', **kwargs)

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
            %(iset)s
            
        :SVD parameters:
            %(nsvd)s
            %(prepca)s
            
        :Returns:
            Arrays with shape ``(nsvd,nt)``
        """

        # Dataset selection
        isets = self._check_isets_(iset)

        # Update params
        self._update_('svd', **kwargs)

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
            %(iset)s

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
        self._update_('svd', **kwargs)

        # Remap
        xeof = self._remap_(xeof)
        xdata = self._remap_(xdata)
            
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
            raw_ec = spanlib_fort.pca_getec(raw_data, raw_eof, raw_weights)
            if raw:
                fmt_ec[iset] = raw_ec
                continue
            ec = npy.ascontiguousarray(raw_ec.T)
            
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
          %(iset)s
        
            
        :SVD parameters:
            %(nsvd)s
            %(prepca)s
            
        :Returns:
            Array with shape ``(nsvd,)`` or a float
        """

        # Update params
        self._update_('svd', **kwargs)

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
            %(iset)s
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
        self._update_('svd', **kwargs)
        
        # Remap alternate arrays
        xeof = self._remap_(xeof)
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
            elif cdms2_isVariable(xpc[iset]): # formatted -> pure numpy
                raw_pc = xpc[iset].filled().T
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
    
    def rec(self, analysis_type=None, *args, **kwargs):
        """Generic method for reconstruction"""
        if analysis_type is None:
            analysis_type = self._last_analysis_type
        else:
            valid = ['pca','mssa','svd']
            if analysis_type not in valid:
                raise SpanlibException('rec','analysis_type must be one of '+str(valid))
        if analysis_type is None:
            warnings.warn('Yet no statistics performed, so nothing to reconstruct!')
        else:
            return getattr(self,self._last_analysis_type+'_rec')(*args,**kwargs)
            
    def clean(self):
        """(Re-)Initialization"""
        dicts = []
        for aa in 'pca','mssa','svd':
            dicts.append('_%s_ev_sum'%aa)
            for bb in 'raw','fmt':
                for cc in 'eof','pc','ev':
                    dicts.append('_%s_%s_%s'%(aa,bb,cc))
        dicts.extend(['_mode_axes','_mssa_window_axes','_mssa_pctime_axes','_mssa_channel_axes','_svd_channel_axes'])
#        lists = ['_mssa_pairs','_nt','_ns','_ndata','_pdata']
        self._nphase = 8
        for ll,func in [(dicts,dict)]:#,(lists,list):
            for att in ll:
                if hasattr(self,att):
                    obj = getattr(self,att)
                    del obj
                setattr(self,att,func())
        self.nd = 0
        self._last_analysis_type = None
        gc.collect()


    def nmssa(self):
        """
        Number of MSSA modes
        
        :Returns:
            integer or tuple
        """
        return self._unmap_(self._nmssa, grouped=True)
        
    def npca(self):
        """
        Number of PCA modes
        
        :Returns:
            integer or tuple
        """ 
        return self._unmap_(self._npca, grouped=True)

    def ns(self):
        """
        Length of channel axis (unmasked input points)
        
        :Returns:
            integer or tuple
        """
        return self._unmap_([d.ns for d in self], grouped=True)
        
    def nt(self):
        """
        Length of time axis
        
        :Returns:
            integer or tuple
        """
        return self._unmap_([d.nt for d in self], grouped=True)
        
    def window(self, absolute=False):
        """
        MSSA window parameter
        
        :Options: 
        
            - *absolute*: if False, return the window relative to time length, 
              else return the effective window (multiplied by nt)
        
        :Returns:
            integer or tuple
        """
        win = self._window
        if absolute:
            win = [int(w*nt) for nt in zip(win, self._nt)]
        return self._unmap_(win, grouped=True)

    def nsvd(self):
        """
        Number of SVD modes
        
        :Returns:
            integer or tuple
        """
        return self._nsvd



    #################################################################
    ## Organize datasets
    #################################################################

    def _map_(self, datasets, sequential):
        """Make datasets in the right form (a list of variables or lists of variables)"""
        
         # Input data tree type
        if not isinstance(datasets, (tuple, list)):
            self.tree_type = 0
        else:
            self.tree_type = 1
            for dataset in datasets:
                if isinstance(dataset, (tuple, list)):
                    self.tree_type = 2
                    break
            else:
                if sequential: self.tree_type = 2
                
        # Group variables if needed
        if self.tree_type < 2:
            datasets = [datasets]
            
        self._datasets = datasets
        self.nd = len(datasets)
            
    def _remap_(self, values, reshape=True, fill_value=None):
        """Makes sure that values is a list (or tuple) of length :attr:`ndatasets`
        
        :func:`broadcast` is called when reshaping.
        """ 
        # We always need a sequence
        if not isinstance(values, (tuple, list)):
            values = [values]
           
        # Check length
        if self.tree_type==2 and len(values)!=len(self): 
            if not reshape:
                raise SpanlibError('Wrong number of input items (%i instead of %i)'
                    %(len(values), len(self)))
            values = broadcast(values, len(self), fill_value)
            
        return values
        
    def _unmap_(self, values, grouped=None):
        """Return values as input dataset (depth and shapes)
        
        :Params:
        
            - **values**: List, list of lists or dictionary with integers as keys.
            - **grouped**, optional: Do not unmap as the :class:`Dataset` level.
        """
        
        # Convert list to dictionary
        if isinstance(values, list):
            if grouped is None: grouped = not isinstance(values[0], list)
            ff = {}
            for i, f in enumerate(values):
                ff[i] = f
            values = ff
        else:
            if grouped is None: grouped = False
            values = copy.copy(values)
        
        # Loop on datasets  
        ret = ()
        for iset in sorted(values.keys()):
            val = values[iset]
            if not grouped:
                val = self[iset].unmap(val)
            ret += val, 
            
        if self.tree_type<2:
            return ret[0]
        return ret

        
        
    
    def _input_raw_(self, input, alt):
        if input is None:
            return alt
        if isinstance(input, dict): # Explicit dictionary
            for i in xrange(self.nd):
                if not input.has_key(i) and alt.has_key(i):
                    input[i] = alt[i]
            return input
        out = {}
        if isinstance(input, (list, tuple)): # Sequence
            for i, data in enumerate(input):
                out[i] = data
        else: # Array (=> same for all dataset)
            for i in xrange(self.nd):
                out[i] = input
        return out

        
    @classmethod
    def _get_imodes_(cls, imode, nmode):
        
        # Which modes
        if imode is None:
            imode = range(0, nmode+1)
        elif isinstance(imode,slice):
            imode = range(imode.start, imode.stop, imode.step)
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
    
    @classmethod
    def _check_length_(cls, input, mylen, fillvalue):
        return broadcast(input, mylen, fillvalue)
    
    def _raw_rec_(self, raw_eof, raw_pc, imodes=None):
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
        imodes = SpAn._get_imodes_(imodes, nmode)

        # Function of reconstruction
        if nw:
            function = spanlib_fort.mssa_rec # MSSA
        else:
            function = spanlib_fort.pca_rec  # PCA/SVD

        # Arguments
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
            
            # fortran call
            ffrec += npy.ascontiguousarray(function(*targs)) # (nc,nt)
            
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
            raise SpanlibError('No iterator have been initialized')
        return span._iter.next()

class SVDModel(SpAn):
    """
    :Params:
        - **method**, optional: Method of reconstruction [default: 'std']. 
        
            - 'std' assumes that left and normalized expansion coefficients are equal [Syu_and_Neelin_1995]_.
            - 'regre' does not use right EOFs but regression coefficients [Harrisson_et al_2002]_.
    """

    def __init__(self,predictor, predictand, method='std', **kwargs):

        kwargs['sequential'] = True
        SpAn.__init__(self,(predictor,predictand) ,**kwargs)

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
            raise SpanlibError('Wrong method for right coefficients from left coefficients. '
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
            leftcoefs = N.atleast_2d(leftcoefs)
            return N.dot(leftcoefs.T, self._l2r)
            
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


class Filler(object):
    """Class to fill missing value with a MSSA filtered version of the data
    
    The initialization automatically call the :meth:`fill` method.
    """
    
    def __init__(self, data, run=True, **kwargs):
        
        self._kwargs = kwargs
        self._kwargs['keep_invalids'] = True
        self._kwargs.setdefault('nvalid', 1)
        self._kwargs.setdefault('quiet', True)
        self._data = data
        self.filled = None
        self.filtered = None
        self.nstep = 0
        self.cv = 100
        
        if run: self.fill(**kwargs)

    def fill(self, nitermax=50, cvmax=2, npca=20, nmssa=15, filter=None, getfiltsteps=False, **kwargs):
        """Run the filler with a convergence loop
        
        Results are accessible in the following attributes:
        
        .. attribute:: filtered
        
            Filtered data, result from the convergence loop.
            
        .. attribute:: filled
        
            Data filled with :attr:`filtered`
            
        :Parameters:
        
            - **nitermax**: Maximal number of iterations
            - **cvmax**: Convergence criterion (%)
            - **npca**: Number of PCA modes (see :class:`SpAn`)
            - **nmssa**: Number of MSSA modes (see :class:`SpAn`)
            - Other parameters are passed to :class:`SpAn`
            
        :Returns:
        
            - :attr:`filled`
        
        """
    
        kwargs['keep_invalids'] = True
        kwargs.setdefault('nvalid', 1)
        kwargs.setdefault('quiet', True)
        data = self._data
        if getfiltsteps: 
            self.filtsteps = []
        else:
            self.filtsteps = None
        old_span = None
        for self.nstep in xrange(nitermax):
            
            # Setup data
            # - run analysis
            if self.nstep:
                kwargs['keep_invalids'] = False
                kwargs['nvalid'] = None
            span = SpAn(data, sequential=False, npca=npca, nmssa=nmssa, **kwargs)
            
            # - first run
            if self.nstep == 0:
                
                # Check input
                if not span.tree_type<2:
                    raise TypeError('Input data must be a single variable or a list of variables')
                
                # Keep original data safe
                mdata = []
                norms = []
                for ivar in xrange(len(span[0])):
                    var = span[0][ivar].data
                    if cdms2_isVariable(var): var = var.clone()
                    else: var = var.copy()  
                    mdata.append(var)
                    norms.append(span[0][ivar].norm)
                    kwargs['norms'] = norms
            
            # Check convergence
            if self.nstep == 0: # reference
            
                self.invalids = invalids = span[0].invalids
                for invalid in invalids:
                    if invalid is not None and invalid.any():
                        break
                else:
                    print "Nothing to fill in"
                    datarec = data
                    if getfiltsteps: self.filtsteps.append(datarec)
                    break
                self.cv = 1.
                last_energy = 0.
                    
            else: # next steps
            
                # new energy
                energy = 0.
                for i, invalid in enumerate(invalids):
                    if invalid is None: 
                            print 'invalid is none'
                            continue
                    basevar = span[0][i].data
                    if span[0][i].has_cdat():
                        basevar = basevar.asma()
                    basevar = span[0][i].scale(basevar, copy=True)#, norm=True)
                    basevar[:] *= basevar
                    energy += basevar[invalid].sum()
                    del basevar
                    
                # check convergence
                self.cv = (energy-last_energy)/energy
                print ' cv', self.cv*100
                if self.cv < 0 or npy.abs(self.cv) < cvmax/100.:
                    print 'Convergence reached: %.2f%% (%i iterations)'%(100*self.cv, self.nstep)
                    break
                last_energy = energy
            
            # Reconstruction
            datarec = span.mssa_rec()
            mdatarec = span[0].remap(datarec)

            # Replace invalid data with reconstructed field
            import pylab as P
            for ivar in xrange(len(mdatarec)):
                if invalids[ivar] is not None:
                    if callable(filter):
                        mdatarec[ivar][:] = filter(mdatarec[ivar])
                    mdata[ivar][:] = eval(span[0][ivar].array_type).where(invalids[ivar], 
                        mdatarec[ivar], mdata[ivar])

            # Save current filtered state
            if getfiltsteps: self.filtsteps.append(datarec)
                        
            data = span._unmap_([mdata])
            old_span = span                    
    
            # Check iteration limit
            if self.nstep >= nitermax:
                print 'Convergence not reached %i%% (%i iterations)'%(100*self.cv, self.nstep)
                break
                
        # Output
        span = old_span if old_span is not None else span
        self.span = span
        self.nmssa = span.nmssa()
        self.npca = span.npca()
        self.filled = data
        self.filtered = datarec
        #FIXME: USE REMAP
        if span[0].has_cdat(): # Ids for CDAT variables
            if isinstance(self.filled, list):
                for i, data in enumerate(span[0]):
                    if not data.has_cdat(): continue
                    self.filled[i].id +='_filled'
                    self.filtered[i].id += '_filtered'
            else:
                self.filled.id += '_filled'
                self.filled.id += '_filtered'
        return data

       

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
    if npy.ma.isMA(stpcs): stpcs = stpcs.filled()
    nt = stpcs.shape[1]
    modes = []
    for mode, pc in enumerate(stpcs):
        result = npy.fft.rfft(pc)
        imax = npy.argmax(result)
        freq = npy.fft.fftfreq(nt)[:nt/2]
        freq_max = freq[imax]
        print mode, 1/freq_max
        if freq_max < low_freq or freq_max > high_freq:
            modes.append(mode)
    return span.mssa_rec(modes=modes)
        
