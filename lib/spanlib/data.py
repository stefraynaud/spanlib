#################################################################################
# File: data.py
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
import numpy
npy = numpy
import _fortran
try:
    import cdms2 ,MV2
    MV = MV2
    cdms = cdms2
    cdms2_isVariable = cdms2.isVariable
    has_cdat_support = True
except:
    cdms2_isVariable = lambda var: False
    has_cdat_support = False
    
default_missing_value = npy.ma.default_fill_value(0.)
from util import Logger, dict_filter, broadcast


class Data(Logger):
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
            - **weights**, optional: masked array
                Flatten in space an [x,y,t] array by 
                removing its masked point
        
        :Options:
        
            - **weights**: Weights to be flatten also
            - **norm**: Normalisation coefficients
            - **mask**: Integer mask where valid data = 1
    """
    def __init__(self, data, weights=None, norm=None, keep_invalids=False, 
        minvalid=None, clean_weights=True, 
        logger=None, loglevel=None, zerofill=False, **kwargs):
        
        # Logger
        Logger.__init__(self, logger=logger, loglevel=loglevel, **dict_filter(kwargs, 'log_'))
        
        # Guess data type and copy
        if cdms2_isVariable(data):
            self.array_type = 'MV2'
            self.array_mod = MV2
            data = data.clone()
        elif npy.ma.isMA(data):
            self.array_type = 'numpy.ma'
            self.array_mod = numpy.ma
            data = data.copy()
        else:
            self.array_type = 'numpy'
            data = data.copy()
            self.array_mod = numpy
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
            if False and weights is not False and data.ndim == 3 and \
                cdms2_isVariable(data) and \
                'x' in data.getOrder() and 'y' in data.getOrder():
                import cdutil# FIXME: WARNING FALSE
                weights = cdutil.area_weights(data[0]).data.astype('d') # Geographic weights
            elif self.nstot==1:
                weights = npy.ones(1)
            else:
                weights = npy.ones(self.shape[1:])
        elif npy.ma.isMA(weights):
            weights = weight.astype('d').filled(0.)
        else:
            weights = npy.asarray(weights, dtype='d')
        if data.ndim>1 and self.shape[1:] != weights.shape:
            self.error('Weights must be of shape %s (instead of %s)'
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
        

        # Masking nans
        nans = npy.isnan(data)
        if nans.any():
            self.warning("Masking %i NaNs"%nans.sum())
            if self.array_type == 'numpy':
                self.array_type = 'numpy.ma'
                data = npy.ma.array(data, mask=nans, copy=False)
            else:
                data[nans] = npy.ma.masked
            self.data = data
            
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
        minvalid = kwargs.pop('nvalid', minvalid)
        if minvalid is not None and minvalid < 0:
            minvalid = -int(round(npy.clip(minvalid, -100, 0)*self.nt/100))
        minvalid = npy.clip(int(minvalid), 1, self.nt) if minvalid is not None else 1
        count[count<minvalid] = 0 # <minvalid -> 0
        count = npy.clip(count, 0, 1)
        # - save as 0/1
        self.ns = long(count.sum())
        self.compress = count.size != self.ns
        self.good = count>0 # points in space where there are enough data in time
        self.minvalid = self.nvalid = minvalid
        
        # Scale unpacked data
        if not self.good.any():
            self.warning('No valid data')
            self.norm = 1.
            self.mean = 0
        else:
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
        # - fill with missing value or mean (0.) where possible
        if minvalid != self.nt: 
#            invalids = bmask & self.good # invalids = masked data that will be analyzed
#            data[invalids] = 0. if zerofill else default_missing_value
#            data[invalids] = default_missing_value
            data[:, ~self.good] = default_missing_value
            if keep_invalids: 
                self.invalids = bmask & self.good # invalids = masked data that will be analyzed
            else: 
                self.invalids = None
                #del invalids
        else:
            self.invalids = None
        # - finally fill with missing values at zero
        if npy.ma.isMA(data):
            data_num = data.filled(default_missing_value)
        else:
            data_num = data

        # Pack 
        # - data
        self.packed_data = self.core_pack(data_num, force2d=True)
        # - weights
        self.packed_weights = self.core_pack(weights)
        

        
    def core_pack(self, data_num, force2d=False):
        """Compress data along space if needed
        
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
        return npy.asfortranarray(pdata)
        
    def scale(self, data, copy=False, mean=None, norm=None, mode=None):
        """Remove mean and normalize unpacked data"""
        if mode=='mean':
            if norm is None:
                norm = False
            if mean is None:
                mean = True
        elif mode=='norm':
            if norm is None:
                norm = True
            if mean is None:
                mean = False            
        if copy: 
            data = data.clone() if cdms2_isVariable(data) else data.copy()
        if mean is not False:
            data[:] -= self.mean if mean is True or mean is None else mean
        if norm is not False:
            data[:] /= self.norm if norm is True or norm is None else norm
        return data
        
    def rescale(self, data, copy=False, mean=None, norm=None, mode=None):
        """Re-add mean and unnormalize unpacked data"""
        if mode=='mean':
            if norm is None:
                norm = False
            if mean is None:
                mean = True
        elif mode=='norm':
            if norm is None:
                norm = True
            if mean is None:
                mean = False            
        if copy: 
            data = data.clone() if cdms2_isVariable(data) else data.copy()
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
        if npy.ma.isMA(data): data = data.filled(default_missing_value)
        
        # Check shape
        nsdim = self.ndim-1
        if nsdim!=0: 
            if data.shape[-nsdim:] != self.shape[-nsdim:]:
                self.error('Incompatible shape of channels (%s instead of %s)'
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
                
            - **firstaxes**, optional: Set the axes as the first ones.
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
                        try:   
                            data.setAxis(i, a)
                        except:
                            pass
             
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
            data[:] = pdata.reshape(data.shape)
        del pdata
        # - mask
        data[:] = npy.ma.masked_values(data, default_missing_value,  copy=0) 
        
        # Rescale
        if rescale:
            self.rescale(data, mode=rescale)

        return data

    def has_cdat(self):
        """Was input array of CDAT type (:mod:`MV2`)?"""
        return self.array_type == "MV2"

    def get_time(self, nt=None, offset=0):
        """Get the time axis or dimension of an input variable
        
        If CDAT is not used, length of axis is returned or nt if different.
        """
        axis = self.taxis
        if not isinstance(axis, int) and nt is not None and nt!=self.nt:
            dt = npy.median(npy.diff(axis[:]))
            axiso = cdms2.createAxis(npy.arange(nt)*dt+axis[0])
            for att, val in axis.attributes.items():
                setattr(axiso, att, val)
            axiso.id = axis.id
            axiso[:] += offset
            return axiso
        if not isinstance(axis, int) and offset:
            axis = axis.clone()
            axis[:] += offset
        elif isinstance(axis, int) and nt is not None and axis!=nt: 
            return nt
        return axis
          
class Dataset(Logger):
    """Class to handle one variable or a list of variables
    
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
        keep_invalids=False, minvalid=None, clean_weights=True, 
        logger=None, loglevel=None, zerofill=False, **kwargs):
        
        # Logger
        Logger.__init__(self, logger=logger, loglevel=loglevel, **dict_filter(kwargs, 'log_'))
        if isinstance(dataset, (list, tuple)):
            dataset = list(dataset)
            self.map = len(dataset)
        else:
            dataset = [dataset]
            self.map = 0
        self.ndataset = self.nd = len(dataset)
        self.dataset = dataset

        # Inits
        self.data = []
        self.nt = None
        weights = self.remap(weights, reshape=True)
        norms = self.remap(norms, reshape=True)
        if self.ndataset==1 and norms[0] is None: norms = [False]
        self._invalids = []
    
        # Loop on datasets
        for idata,data in enumerate(dataset):
    
            # Create the Data instance and pack array
            dd = Data(data, norm=norms[idata], weights=weights[idata], 
                keep_invalids=keep_invalids, minvalid=minvalid, clean_weights=clean_weights, zerofill=zerofill)
            self.data.append(dd)
            self._invalids.append(dd.invalids)
            
            # Check nt
            if self.nt is None:
                self.nt = dd.nt
            elif self.nt != dd.nt:
                self.error('Time dimension of variable %i must have length %i (not %i)'%(idata, self.nt, dd.nt))
    
        # Merge
        self.stacked_data = npy.asfortranarray(npy.vstack([d.packed_data for d in self.data]))
        self.splits = npy.cumsum([d.packed_data.shape[0] for d in self.data[:-1]])
        self.stacked_weights = npy.hstack([d.packed_weights for d in self.data])
        self.ns = self.stacked_data.shape[0]
        self.ntv = (self.stacked_data!=default_missing_value).any(axis=0).sum()
        
    def get_norms(self, idata=None):
        """Get :attr:`norms` for one or all input variables"""
        if idata is None:
            return [d.norm for d in self.dataset.data]
        return self[idata].norm
    norms = property(get_norms, doc="Normalization factors")
    
    def get_invalids(self, stacked=True):
        if self._invalids[0] is None: return
        if not stacked: return self._invalids
        return self.restack(self._invalids, scale=False)
    invalids = property(fget=get_invalids, doc='Final mask of stacked data')
        
    def restack(self, dataset, scale=True):
        """Stack new variables as a fortran array
        
        It has the opposite effect of :meth:`unstack`.  
        
        :Params:
        
            - **dataset**: Argument in the same form as initialization data.
            - **scale**, optional: Scale the variable (mean and norm), and optionally 
              remove spatial mean if == 2.
        
        :Seel also: :meth:`Data.repack`
        """
        
        # Check input data
        if isinstance(dataset, (list, tuple)):
            dataset = list(dataset)
            dmap = len(dataset)
        else:
            dataset = [dataset]
            dmap = 0
        if  len(dataset)!=self.ndataset:
            self.error('You must provide %i variable(s) to stack'%len(dataset))
            
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
                    self.error('Time length of variable %i (%i) different from that of first variable (%i)'
                        % (i+1, nt2, nt1))
            
        # Stack
        stacker = npy.vstack if packs[0].ndim==2 else npy.hstack
        sdata = npy.asfortranarray(stacker(packs))
        
        return sdata
    
 
    def unstack(self, sdata, rescale=True, format=1, firstdims=None, firstaxes=None):
        """Unstack and unpack data
        
        It has the opposite effect of :meth:`restack`.

        
        :Params:
        
            - **sdata**: Stacked data (as computed by :meth:`__init__` or :meth:`restack`).
            - **rescale**, optional: Rescale the variable (mean and norm) and 
              add spatial mean if == 2.
            - **format**, optional: Format the variable (see :meth:`Data.create_array`).
            - **firstaxes**, optional: First axis (see :meth:`Data.create_array`).
            
            
        :Seel also: :meth:`Data.unpack`
        """
        
        # Unstack
        spliter = npy.hsplit if sdata.ndim==1 else npy.vsplit
        packs = spliter(sdata, self.splits)
        
        # Unpack
        return [self[i].unpack(pdata, rescale=rescale, format=format, 
            firstdims=firstdims, firstaxes=firstaxes)
            for i, pdata in enumerate(packs)]
            
    def fill_invalids(self, dataref, datafill, raw=False,  copy=False, unmap=True, 
        missing=False):
        """Fill ``dataref`` with ``datafill`` at registered invalid points 
        or current missing points 
        
        
        :Params:
        
            - **dataref**: Reference dataset than must be filled.
            - **datafill**: Dataset used to fill reference.
            - **raw**, optional: Input data (and mask through missing) are aleardy packed.
            - **missing**, optional: If True, gaps are defined by ``dataref``
              missing values. If False, gaps are those of original data. If an array,
              it is used as the mask to define gaps.
        """
#        if self.invalids is None: 
#            return dataref.clone() if copy else dataref # FIXME: 
        
        # Re stack to have raw data
        if not raw:
            if dataref!='stacked_data': 
                dataref = self.restack(dataref)
            datafill = self.restack(datafill)
            copy = False
        if dataref=='stacked_data': 
             dataref = self.stacked_data
             
        # Check raw shape
        if dataref.shape!=datafill.shape:
            self.error("Can't replace with raw data because of wrong shape (%s!=%s)"%
                (data.shape, datafill.shape))
        
        # Copy for safety
        if copy: dataref = dataref.copy()
        
        # Put it at invalid points
        if missing is not False and missing is not True:
            if raw:
                mask = missing
            else:
                mask = self.restack(missing)
        elif missing:
            mask = npy.ma.masked_values(dataref, default_missing_value, shrink=False).mask
        else:
            mask = self.invalids
        dataref[:] = npy.where(mask, datafill, dataref)            
        del mask
        
        # Unstack ?
        if int(raw)>0:
            data = npy.asfortranarray(dataref)
        else:
            data = self.unstack(dataref)
            if unmap: data = self.unmap(data)
            
        return data
        
                    
    def replace_invalids(self, datafill, raw=False, copy=False):
        if self.invalids is None: return 
        
        # Re stack to have raw data if needed
        if not raw:
            datafill = self.restack(datafill)
        
        # Fill
        return self.fill_invalids(self.stacked_data, datafill, raw=True, copy=False)
      
       
        
    
#    def get_time(self, idata=0):
#        """Get time axis"""
#        return self.data[idata].get_time()



    def __len__(self):
        return self.ndataset
        
    def __getitem__(self, i):
        return self.data[i]

    def unmap(self, out):
        """Remap out so as to be in the same form as input data
        
        It has the opposite effect of :meth:`remap`.
        """
        if self.map==0:
            return out[0]
        return out

    def remap(self, values, reshape=True, fill_value=None):
        """Makes sure that values is a list (or tuple) of length :attr:`ndataset`
        
        :func:`broadcast` is called when reshaping.
        
        It has the opposite effect of :meth:`unmap`.
        """         
        # We always need a list
        if not isinstance(values, (list, tuple)):
            values = [values]
            
        # Check length
        n = len(self)
        if len(values)!=n:
            if not reshape:
                self.error('Wrong number of input items (%i instead of %i)'
                    %(len(values), n))
            values = broadcast(values, n, fill_value)
        return values

    def has_cdat(self, idata=None):
        """Check if there were at least one input array of CDAT type (:mod:`MV2`)
        
        :Sea also: :meth:`Data.has_cdat`
        """
        # Not at all
        if not has_cdat_support: return False
        
        # Single var
        if idata is not None:
            return self[idata].has_cdat()
        
        # At least one
        for d in self.data:
            if d.has_cdat(): return True
        return False

    def get_time(self, nt=None, offset=0, idata=0):
        """Get the time axis of an input variable  
        
        If CDAT is not used, length of axis is returned.
        """
        return self[idata].get_time(nt, offset)
            


