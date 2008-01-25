#################################################################################
# File: spanlib_python.py
#
# This file is part of the SpanLib library.
# Copyright (C) 2006-2008  Charles Doutriaux, Stephane Raynaud
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

import spanlib_fort,numpy as npy,genutil.statistics,genutil.filters,cdms2 as cdms,MV2 as MV
import copy,gc
from warnings import warn

def _pack_(data,weights=None,norm=None):
	""" Pack a dataset and its weights according to its mask

	Description:::
	  This function packs a dataset in 2D space by removing
	  all masked points and returning a space-time array.
	  It performs this operation also on the weights.
	  It is used for removing unnecessary points and
	  simplifying the input format for analysis functions.
	:::

	Usage:::
	packed_data, packed_weights, mask = _pack_(data,weights)

	  data	  :: Flatten in space an [x,y,t] array by removing
	    its masked point
	  weights :: Weights to be flatten also
	:::

	Output:::
	  packed_data	 :: Space-time packed array
	  packed_weights :: Packed weights that were guessed or used
	  mask		     :: Mask that were guessed or used
	:::
	"""

	# Total number of channels
	packed = {}
	if not cdms.isVariable(data):
		data = cdms.createVariable(data,copy=0)
		if npy.core.ma.isMa(data):
			packed['type'] = 'ma'
		else:
			packed['type'] = 'numpy'
	else:
		packed['type'] = 'cdms'
	nt = data.shape[0]
	nstot = npy.multiply.reduce(data.shape[1:])
	print 'pack data shape',data.shape
	sh=list(data.shape)
	packed['norm'] = norm

	# Is it already packed? Check the mask...
	if data.mask is False:
		if weights is None:
			weights = npy.ones(nstot,dtype='f')
		else:
			weights = npy.asarray(weights)
		packed['data'] = data.filled().reshape((nt,nstot)).transpose()
		packed['weights'] = weights.ravel()
		packed['mask'] = npy.ones(nstot)
		check_norm(packed)
		return packed

	# Weights ?
	if weights is None:
		if data.ndim == 3 and \
			data.getAxis(-1).isLongitude() and data.getAxis(-2).isLatitude():
			import cdutil
			weights = cdutil.area_weights(data[0]).raw_data() # Geographic weights
		else:
			weights = npy.ones(nstot,dtype='f').reshape(sh[1:])
	else:
		if cdms.isVariable(weights):
			weights = weights.filled(0.)

	# Mask
	# - First from data
	mask = npy.asarray(npy.sum(MV.getmaskarray(data),axis=0),dtype='i') # no time
	# - Now add the ones from the weights
	mask[:] = mask+npy.equal(weights,0.)
	# - >=1 means masked, Fortran "mask": 1 means data ==> 1-mask
	mask[:] = 1-npy.clip(mask,0,1)
	packed['mask'] = mask

	# Number of valid points in spatial dimension
	ns = long(mask.sum())

	# Pack space
	# - Pack numeric data array
	data_to_pack = data.filled(1.e20).reshape((nt,nstot)) 
	packed['data'] = npy.asarray(spanlib_fort.chan_pack(data_to_pack,mask.flat,ns),
		dtype='f',order='F')
	del data_to_pack
	# - Pack weights
	weights = weights.reshape((1,nstot))
	packed['weights'] = npy.asarray(spanlib_fort.chan_pack(weights,mask.flat,ns)[:,0],
		dtype='f')

	_check_norm_(packed)
	return packed

def _sort_modes_(mode1,mode2):
	return npy.sign(abs(mode2)-abs(mode1))
	
def _check_shape_(inputs,datasets,fillvalue):
	"""Return input as datasets *shape*"""
	inputs = _check_length_(inputs,len(datasets),fillvalue)
	for iset,dataset in enumerate(datasets):
		inputs[iset] = _check_length_(inputs[iset],len(dataset),fillvalue)
	return inputs

def _check_length_(input,mylen,fillvalue):
	if not isinstance(input,(list,tuple)):
		input = [input]
	dlen = mylen-len(input)
	if dlen < 0:
		input = input[:mylen]
	elif dlen > 0:
			input.extend([fillvalue]*dlen)
	return input

def _check_norm_(packed):
	if packed['norm'] in [True,None]:
		packed['norm'] = packed['data'].std() # Standard norm
		packed['data'] /= packed['norm']
	elif packed['norm'] is not False:
		if packed['norm'] <0: # Relative norm, else strict norm
			packed['norm'] = packed['norm']*data.std()
		packed['data'] /= packed['norm']
	else:
		packed['norm'] = 1.
			
	

def get_pairs(pcs,mincorr=0.9,maxdistance=2,smooth=5):
	""" Get pairs in MSSA results

	Description:::
	  This function detects pairs of mode in principal
	  components (typical from MSSA). Modes of a pair
	  have their PCs (and EOFs) in phase quadrature in time.
	  The algorithm uses correlations between a pc and the 
	  derivative of the other pc.
	:::

	Usage:::
	pairs = get_pairs(pcs)

	  pcs         :: Principal components [modes,time]
	  mincorr     :: Minimal correlation [default:0.95]
	  maxdistance :: Maximal difference of indices between a pair of modes [default: 3]
	  smooth      :: Apply 'smooth'-running averages during scan [default: 5].
	                 Not applied if smooth < 1.
	:::

	Output:::
	  pairs :: List of two-element tuples of indices
	:::
	"""

	pcs = MV.filled(pcs)
	nst = pcs.shape[1]
	if nst < 4:
		return []
	smooth = npy.clip(smooth,1,max([1,nst-3]))
	pairs = []
	found = []
	minCorr = npy.clip(minCorr,0.,1.)
	maxDistance = npy.clip(maxDistance,1,len(pcs)-1)
	for ip1 in xrange(len(pcs)-1):
##		print 'ipc',i
##		print 'found',found
		# Mode already found
		if ip1 in found:
			continue
		# Scan next maxDistance modes to find a possible pair
		for ip2 in xrange(ip1+1,min([ip1+maxDistance+1,len(pcs)])):
##			print 'TRY',ip1,ip2
			ip = [ip1,ip2]
			corr = 0.
			for k in 0,1:
				# PC without extreme bounds
				pc = pcs[ip[k],1:-1]
				if smooth > 1:
					pc = genutil.filters.runningaverage(pc,smooth)
				# Derivative
				dpc = pcs[ip[1-k],2:] - pcs[ip[1-k],:-2]
				if smooth > 1:
					dpc = genutil.filters.runningaverage(dpc,smooth)				
				# Increment correlation
				corr += float(npy.absolute(genutil.statistics.correlation(pc,dpc,axis=0)))
##			print 'ssssCORRELATION',float(genutil.statistics.correlation(pc,dpc,axis=0))
##			import vcs
##			v1,v2=vcs.init(),vcs.init()
##			v1.plot(pc)
##			v2.plot(dpc)
##			raw_input('map out ok?')
			if corr >= 2*minCorr:
				# We got a pair!
				print 'We got a pair!'
				break
		else:
			# No pair found
			continue
		found.extend(ip)
		pairs.append(tuple(ip))

	return pairs


def compute_phases(data,nphases=8,offset=.5,firstphase=0):
	""" Phase composites for oscillatory fields

	Description:::
	  This computes temporal phase composites of a spatio-temporal
	  dataset. The dataset is expected to be oscillatory in time.
	  It corresponds to a reoganisation of the time axis to
	  to represents the dataset over its cycle in a arbitrary
	  number of phases. It is useful, for example, to have a
	  synthetic view of an reconstructed MSSA oscillation.
	:::

	Usage:::
	phases = computePhases(data,nphases,offset,firstphase)

	  data	   :: Space-time data oscillatory in time data.shape is rank 2 and dim 0 is space
	  nphases	:: Number of phases (divisions of the cycle)
	  offset	 :: Normalised offset to keep higher values only [default:
	  firstphase :: Position of the first phase in the 360 degree cycle
	:::

	Output:::
	  phases :: Space-phase array
	:::
	"""

	ns=data.shape[0]
	nt=data.shape[1]
	w = MV.ones((ns),dtype='f')
	phases = MV.array(spanlib_fort.phasecomp(data, nphases, w, offset, firstphase))
	axes = MV.array(data).getAxisList()
	phases.id = 'Phase composites'
	ax = phases.getAxis(1)
	ax[:]=ax[:]*360./nphases+firstphase
	ax.id = 'phases'
	ax.long_name = 'Phases'
	ax.units = 'degrees'
	axes[1] = ax
	phases.setAxisList(axes)
	return phases


class SpAn(object):
	
	_npca_default = 10
	_npca_max = 30
	_nmssa_default = _nsvd_default = 8
	_nmssa_max = _nsvd_max = 20
	_window_default = 1/3. # Relative to time length

	def __init__(self, datasets, serial=False, weights=None, norms=None, npca=None, prepca=None, window=None, nmssa=None, nsvd=None):
		""" Prepare the Spectral Analysis Object

		Description:::
		  This function creates an object for future analyses.
		  It optionally initializes some parameters.
		:::

		Usage:::
		analysis_object = SpAn(datasets,weights=None,npca=None,window=None,nmssa=None)

		  data	:: List of data on which to run the PC Analysis
		    Last dimensions must represent the spatial dimensions.
		    Analysis will be run on the first dimension.
		  weights :: If you which to apply weights on some points,
		    set weights to "0" where you wish to mask.
		    The input data mask will be applied,
		    using the union of all none spacial dimension mask.
		    If the data are on a regular grid, area weights
		    will be generated, if the cdutil (CDAT) module is available.
		    default: 1. everywhere]
		  npca	:: Number of principal components to return [default: 10]
		  nmssa   :: Number of MSSA modes retained [default: 4]
		  nsvd	:: Number of SVD modes retained [default: 10]
		  window  :: MSSA window parameter [default: time_length/3.]
		:::

		Output:::
		  analysis_object :: SpAn object created for further analysis
		:::
		"""

		self.clean()

		# We start from a list
		if not isinstance(datasets,(list,tuple)):
			datasets = [datasets]
			self._input_map = 0
		else:
			self._input_map = len(datasets)

		# Check if we are in a serial analysis mode
		for d in datasets:
			if isinstance(d,(list,tuple)):
				serial = True
				break
		if not serial: # Convert to serial mode
			datasets = [datasets]
		else:
			input_map = []
			for iset in xrange(len(datasets)):
				if isinstance(datasets[iset],(list,tuple)):
					input_map[iset] = len(datasets[iset])
				else:
					input_map[iset] = 0
			self._input_map = input_map

		# Weights and norms in the same form as datasets
		if norms is None:	norms = 1.
		norms = _check_shape_(norms,datasets,1.) 
		weights = _check_shape_(weights,datasets,1.)

		# We stack and pack data
		for iset,d in enumerate(datasets):
			self._stack_(d,weights[iset],norms[iset])
		self._ndataset = len(datasets)

		# Check and save parameters
		self._update_params_(npca=npca,window=window,nmssa=nmssa,nsvd=nsvd,prepca=prepca)


	#################################################################
	## Get datasets info
	#################################################################

	def _time_axis_(self,iset,idata=None):
		"""Get the time axis of data variable of a dataset"""
		if idata is None:
			return tuple(self._stack_info[iset]['taxes'])
		return self._stack_info[iset]['taxes'][idata]

	def _space_axis_(self,iset,idata,iaxis=None):
		"""Get a sp_ace axis of data variable of a dataset"""
		if iaxis is None:
			return self._stack_info[iset]['saxes'][idata]
		else:
			return self._stack_info[iset]['saxes'][idata][iaxis]

	def _mode_axis_(self,analysis_type,isets=None):
		"""Get a mode axis according to the type of modes (pca, mssa, svd)"""
		if not self._mode_axes.has_key(analysis_type):
			self._mode_axes[analysis_type] = {}
		if isets is None: 
			isets = xrange(self._ndataset)
		else:
			isets = [isets]
		for iset in isets:
			if not self._mode_axes[analysis_type].has_key(iset):
				self._mode_axes[analysis_type][iset] = cdms.createAxis(npy.arange(1,getattr(self,'_n'+analysis_type)[iset]+1))
				self._mode_axes[analysis_type][iset].id = analysis_type+'_mode'
				self._mode_axes[analysis_type][iset].long_name = analysis_type.upper()+' modes in decreasing order'
				self._check_dataset_tag_('_mssa_mode_axes',iset,analysis_type)
		if isinstance(isets,list):
			return self._mode_axes[analysis_type][iset]
		return self._mode_axes[analysis_type]

	def _mssa_window_axis_(self,iset):
		"""Get the window axis for mssa for one dataset"""
		if not self._mssa_window_axes.has_key(iset) or len(self._mssa_window_axes[iset]) != self._window[iset]:
			self._mssa_window_axes[iset] = cdms.createAxis(npy.arange(self._window[iset]))
			self._mssa_window_axes[iset].id = 'mssa_window'
			self._mssa_window_axes[iset].long_name = 'MSSA window time'
			self._check_dataset_tag_('_mssa_window_axes',iset)
		return self._mssa_window_axes[iset]

	def _mssa_channel_axis_(self,iset):
		"""Get the channel axis for mssa for one dataset"""
		if not self._channel_axes.has_key(iset) or len(self._mssa_channel_axes[iset]) != self._nmssa[iset]:
			self._mssa_channel_axes[iset] = cdms.createAxis(npy.arange(self._nmssa[iset]))
			self._mssa_channel_axes[iset].id = 'mssa_channel'
			self._mssa_channel_axes[iset].long_name = 'MSSA channels'
			self._check_dataset_tag_('_mssa_channel_axes',iset)
		return self._mssa_channel_axes[iset]

	def _mssa_pctime_axis_(self,iset,idata=0):
		"""Get the channel axis for mssa for one dataset"""
		nt = self._nt[iset] - self._window[iset] + 1
		if not self._mssa_pctime_axes.has_key(iset) or len(self._mssa_pctime_axes[iset]) != nt:
			self._mssa_pctime_axes[iset] = cdms.createAxis(npy.arange(nt))
			self._mssa_pctime_axes[iset].id = 'mssa_pctime'
			self._mssa_pctime_axes[iset].long_name = 'MSSA PC time'
			self._check_dataset_tag_('_mssa_pctime_axes',iset)
			taxis = self._time_axis_(iset,idata)
			if hasattr(taxis,'units') and taxis.units.split()[0].lower() in \
				['seconds','minutes','hours','days','months','years']:
				self._mssa_pctime_axes[iset].units = taxis.units.split()[0].lower() + ' since 0001-01-01'
				self._mssa_pctime_axes[iset].designateTime()
		return self._mssa_pctime_axes[iset]
		
	def _check_dataset_tag_(self,name,iset,key=None,long_name=True,id=True):
		"""Mark some attributes as specific to a dataset (only if there are more then one dataset)
			iset:: ID of the dataset
			key:: A dictionary key to select the dataset [default: None]
			long_name:: Mark the long name [defualt: True]
			id:: Mark the id [default: True]
		"""
		if self._ndataset > 1:
			targetset = getattr(self,name)
			if key is not None:
				targetset = targetset[key]
			target = targetset[iset]
			if id: 
				target.id += str(iset)
			if long_name:
				target.long_name += ' for dataset #%i'%iset

	def _norm_(self,iset,idata):
		"""Get the normalization factor of one subdataset of one dataset"""
		return self._stack_info[iset]['norms'][idata]
		
	def _type_(self,iset,idata):
		"""Get 'numpy', 'MA' or 'cdms' for one subdataset of one dataset"""
		return self._stack_info[iset]['types'][idata]

	def _changed_param_(self,old,param):
		"""Check if a parameter is changed for all datasets.
		Return a list of booleans.
		"""
		return [old[param][iset]==getattr(self,'_'+name)[iset] for iset in xrange(self._ndataset)]
		
	def _update_params_(self,**kwargs):
		"""Update and check requested statistical paremeters.
		A value of None is converted to an optimal value.
		"""
		old = {}
		changed = {}
		for param in 'npca','nmssa','window','prepca','nsvd':
			changed[param] = [False]*self._ndataset # Init
			if kwargs.has_key(param):
				old[param] = getattr(self,'_'+param,None)
				if param != 'nsvd':
					setattr(self,'_'+param,_check_length_(kwargs.pop(param),self._ndataset,None))
				else:
					setattr(self,'_'+param,kwargs.pop(param))
					
		# Number of PCA modes		
		if kwargs.has_key('npca'):
			for iset in xrange(self._ndataset):
				if self._npca[iset] is None: self._npca[iset] = SpAn._npca_default # Default value
				self._npca[iset] = npy.clip(self._npca[iset],1,
					min(SpAn._npca_max,self._ns[iset])) # Max
			changed['npca'] = self._changed_param_(old,'npca')
			
		# Should we perform pre-PCA before MSSA or SVD?
		if kwargs.has_key('prepca'):
			for iset in xrange(self._ndataset):
				if self._prepca[iset] is None: # Default: pre-pca needed over max (for MSSA and SVD)
					self._prepca[iset] = self._ns[iset] > SpAn._npca_max
					if self._prepca[iset]:
						print '[mssa] The number of valid points of one of the datasets is greater than %i, so we perform a pre-PCA'%SpAn._npca_max
			changed['prepca'] = self._changed_param_(old,'npca')
			
		# Number of MSSA modes
		if kwargs.has_key('nmssa'):
			for iset in xrange(self._ndataset):
				if self._nmssa[iset] is None: self._nmssa[iset] = SpAn._nmssa_default # Default value
				if self._prepca[iset]:
					nchanmax = self._npca[iset] # Input channels are from PCA
				else:
					nchanmax = self._ns[iset] # Input channels are from real space
				self._nmssa[iset] = npy.clip(self._nmssa[iset],1,
					min(SpAn._nmssa_max,nchanmax)) # Max
			changed['nmssa'] = self._changed_param_(old,'nmssa')
			
		# Window extension of MSSA
		if kwargs.has_key('window'):
			for iset in xrange(self._ndataset):
				if self._window[iset] is None:self._window[iset] = int(self._nt[iset]*SpAn._window_default)
				self._window[iset] = npy.clip(self._window,1,max(1,self._nt[iset]))
			changed['window'] = self._changed_param_(old,'window')
			
		# Number of SVD modes
		if kwargs.has_key('nsvd'):
			if self._nsvd is None: self._nsvd = SpAn._nsvd_default # Default value
			for iset in xrange(self._ndataset):
				if self._prepca[iset]:
					nchanmax = self._npca[iset] # Input channels are from PCA
				else:
					nchanmax = self._ns[iset] # Input channels are from real space
				self._nsvd = npy.clip(self._nsvd,1,
					min(SpAn._nsvd_max,nchanmax)) # Max
			changed['nsvd'] = self._changed_param_(old,'nsvd')
			
		# Inform which params have been modified for each dataset
		return changed

	def _check_isets_(self,iset):
		"""Check if an iset is a valid dataset.
		It can be a list, and it is returned as a list.
		if an iset is invalid, it is removed from the output list.
		"""
		if iset is not None and not isinstance(iset,(list,tuple)):
			isets = [iset]
		else:
			isets = iset
		if isets is not None:
			for i,iset in enumerate(isets):
				if iset < 0 or iset > self._ndataset:
					warn('Invalid dataset id: %i. Valid id are < %i'%(iset,self._ndataset))
					del isets[i]
		return isets


	#################################################################
	## PCA
	#################################################################

	def pca(self,get_ev_sum=False,relative=False,iset=None,**kwargs):
		""" Principal Components Analysis (PCA)

		Descriptions:::
		  This function performs a PCA on the analysis objects
		  and returns EOF, PC and eigen values.
		  EOF are automatically unpacked.
		:::

		Usage:::
		  pca(npca=None,weights=None,relative=False)

		  iset	:: Dataset selection. If None, all are analysed [default: None]
		  npca	:: Number of principal components to return [default: 10]
		  get_ev_sum  :: Also return sum of all eigen values [default: False]
		  relative :: Egein values are normalized to their sum (%) [default: False]

		:::
		"""

		# Check on which dataset to operate
		isets = self._check_isets_(iset)

		# Update params
		self._update_params_(**kwargs)
		print 'doing pca'

		# Loop on datasets
		for iset,pdata in enumerate(self._pdata):
			print '   ',iset
			
			# Operate only on selected datasets
			if isets is not None and iset not in isets: continue
			
			# Remove old results
			for att in 'raw_eof','raw_pc','raw_ev','ev_sum':
				dic = getattr(self,'_pca_'+att)
				if dic.has_key(iset): del dic[iset]

			# Compute PCA
			if pdata.ndim == 1: # One single channel, so result is itself
				raw_eof = npy.ones(1,dtype=pdata.dtype)
				raw_pc = pdata
				raw_ev = raw_pc.var()
				ev_sum = ev

			else: # Several channels
				weights = self._stack_info[iset]['weights']
				raw_eof,raw_pc,raw_ev,ev_sum = spanlib_fort.pca(pdata,self._npca[iset],weights,-1)	

			# Append results
			self._pca_raw_pc[iset] = raw_pc
			self._pca_raw_eof[iset] = raw_eof
			self._pca_raw_ev[iset] = raw_ev
			self._pca_ev_sum[iset] = ev_sum
			gc.collect()

		self._last_analysis_type = 'pca'


	def pca_eof(self,update=False,iset=None,*args,**kwargs):
		"""Get EOFs from current PCA decomposition

		If PCA was not performed ro if update is True, it is done with all parameters sent to pca()
		"""
	
		# Dataset selection
		isets = self._check_isets_(iset)

		# Update params
		changed =  self._update_params_(**kwargs)['npca']
		print 'eof changed',changed

		# Of, let's format the variable
		fmt_eof = {}
		for iset in xrange(self._ndataset):
			
			# Check dataset
			if isets is not None and iset not in isets:
				continue
				
			# EOF already available 
			if self._pca_fmt_eof.has_key(iset) and not update:
				fmt_eof[iset] = self._pca_fmt_eof[iset]
				continue
				
			# Should we rerun PCA for this dataset?
			if not self._pca_raw_eof.has_key(iset) or update or changed[iset]:
				self.pca(iset=iset,*args,**kwargs)
			if self._pca_fmt_eof.has_key(iset): del self._pca_fmt_eof[iset]
			
			# Get raw data back to physical space
			self._pca_fmt_eof[iset] = \
				self._unstack_(iset,self._pca_raw_eof[iset],self._mode_axis_('pca'))
			
			# Set attributes
			for idata,eof in enumerate(self._pca_fmt_eof[iset]):
				if not self._stack_info[iset]['ids'][idata].startswith('variable_'):
					eof.id = self._stack_info[iset]['ids'][idata]+'_pca_eof'
				else:
					eof.id = 'pca_eof'
				eof.name = eof.id
				eof.standard_name = 'empirical_orthogonal_functions_of_pca'
				eof.long_name = 'PCA empirical orthogonal functions'
				atts = self._stack_info[iset]['atts'][idata]
				if atts.has_key('long_name'):
					eof.long_name += ' of '+atts['long_name']
				if atts.has_key('units'):
					del eof.units
			fmt_eof[iset] = self._pca_fmt_eof[iset]

		return self._return_(fmt_eof)		


	def pca_pc(self,update=False,*args,**kwargs):
		"""Get PCs from current PCA decomposition

		If PCA was not performed, it is done with all parameters sent to pca()
		"""
		# Dataset selection
		if not isinstance(iset,(list,tuple)):
			isets = [iset]
		else:
			isets = iset

		# Of, let's format the variable
		fmt_pc = {}
		for iset,raw_pc in self._pca_raw_pc.items():
			
			# Check dataset
			if isets is not None and iset not in isets:
				continue
				
			# PC already available 
			if self._pca_fmt_pc.has_key(iset) and not update:
				fmt_pc[iset] = self._pca_fmt_pc[iset]
				continue
			
			# Should we rerun PCA for this dataset?
			if not self._pca_raw_pc.has_key(iset) or update:
				self.pca(iset=iset,*args,**kwargs)
			if self._pca_fmt_eof.has_key(iset): del self._pca_fmt_eof[iset]

			# Format the variable
			idata = 0 # Reference is first data
			pc = cdms.createVariable(npy.asarray(raw_pc.transpose(),order='C'))
			pc.setAxis(0,self._mode_axis_('pca',iset))
			pc.setAxis(1,self._time_axis_(iset,idata))
			pc.id = 'pca_pc'
			pc.standard_name = 'principal_components_of_pca'
			
			pc.name = pc.id
			pc.long_name = 'PCA principal components'
			atts = self._stack_info[iset]['atts'][idata]
			if atts.has_key('long_name'): pc.long_name += ' of '+atts['long_name']
			if atts.has_key('units'):     pc.units = atts['units']
			
			fmt_pc[iset] = self._pca_fmt_pc[iset] = pc
			self._check_dataset_tag_('_pca_fmt_pc',iset)

		return self._return_(fmt_pc)		


	def pca_ev(self,relative=False,sum=False,cumsum=False,update=False,*args,**kwargs):
		"""Get eigen values from current PCA decomposition

		Inputs:
		  relative :: Return percentage of variance
		  sum :: Return the sum of eigen values (total variance)
		"""

		# Dataset selection
		if not isinstance(iset,(list,tuple)):
			isets = [iset]
		else:
			isets = iset
			
		# Loop on dataset
		res = {}
		for iset in xrange(self._ndataset):
			
			# Check dataset
			if isets is not None and iset not in isets:
				continue
				
			# Should we rerun PCA for this dataset?
			if not self._pca_ev.has_key(iset) or update:
				self.pca(iset=iset,*args,**kwargs)

			# We only want the sum
			if sum:
				res[iset] = self._pca_ev_sum[iset]
				continue

			# Format the variable
			id = 'pca_ev'
			long_name = []
			raw_ev = self._pca_raw_ev[iset]
			if cumsum:
				raw_ev = raw_ev.cumsum()
				id += '_cumsum'
				long_name.append('cumulative')
			if relative: 
				raw_ev = 100.*raw_ev/self._pca_ev_sum[iset]
				id += '_rel'
				long_name.append('relative')
			ev = cdms.createVariable(raw_ev)
			ev.id = ev.name = id
			long_name.append('PCA eigen values')
			ev.long_name = ' '.join(long_name).title()
			ev.setAxisList([self._mode_axis_('pca',iset)])
			ev.standard_name = 'eigen_values_of_pca'
			atts = self._stack_info[iset]['atts'][0]
			if atts.has_key('long_name'):
				ev.long_name += ' of '+atts['long_name']
			if relative:
				ev.units = '%'
			elif atts.has_key('units'):
				ev.units = atts['units']
				for ss in ['^','**',' ']:
					if ev.units.find(ss) != -1:
						ev.units = '(%s)^2' % ev.units
						break
			res[iset] = ev

		return self._return_(res)		

	def pca_rec(self,imode=None,update=False,iset=None):
		"""Reconstruct a set of modes from PCA decomposition

		Inputs:
		  imode :: Selection of Modes. Can be like [1,3,-5] -> modes [1,3,4,5], or None -> all [default: None].
		  update :: Rerun the PCA.
		"""
		
		# Dataset selection
		if not isinstance(iset,(list,tuple)):
			isets = [iset]
		else:
			isets = iset
			
		# Loop on datasets
		pca_fmt_rec = {}
		for iset in xrange(self._ndataset):
			
			# Check dataset
			if isets is not None and iset not in isets:
				continue
				
			# Should we rerun PCA for this dataset?
			if not self._pca_raw_eof.has_key(iset) or update:
				self.pca(iset=iset,*args,**kwargs)
			
			# Get raw data back to physical space
			raw_rec,smodes = self._project_(iset,self._pca_raw_eof,self._pca_raw_pc,imode)
			pca_fmt_rec[iset] = self._unstack_(iset,raw_rec,self._time_axis_(iset))
			del  raw_rec
			
			# Set attributes and scale
			for idata,rec in enumerate(pca_fmt_rec[iset]):
				rec[:] *= self._norm_(iset,idata) # Scale
				if not self._stack_info[iset]['ids'][idata].startswith('variable_'):
					rec.id = self._stack_info[iset]['ids'][idata]+'_pca_rec'
				else:
					rec.id = 'pca_rec'
				rec.name = rec.id
				if imode is not None:
					rec.id += smodes
				rec.standard_name = 'recontruction_of_pca_modes'
				rec.long_name = 'Reconstruction of PCA modes: '+smodes
				atts = self._stack_info[iset]['atts'][idata]
				if atts.has_key('long_name'):
					rec.long_name += ' of '+atts['long_name']
					
		return self._return_(pca_fmt_rec)	
	
	

	#################################################################
	# MSSA
	#################################################################

	def mssa(self,**kwargs):
		""" MultiChannel Singular Spectrum Analysis (MSSA)

		Description:::
		  This function performs a MSSA on the analysis objects
		  and returns EOF, PC and eigen values.
		  Unless pca parameter is set to false, a pre
		  PCA is performed to reduced the number of d-o-f
		  if already done and if the number of channels is
		  greater than 30.
		:::

		Usage:::
		eof, pc, ev = mssa(nmssa,pca,relative=False)

		OR

		eof, pc, ev, ev_sum = mssa(nmssa,pca,get_ev_sum=True,relative=False)

		  nmssa  :: Number of MSSA modes retained
		  window :: MSSA window parameter
		  prepca	:: If True, performs a preliminary PCA
		  get_ev_sum  :: Also return sum of all eigen values (default: False)
		  relative :: Egein values are normalized to their sum (%) [default: False]

		Output:::
		  eof :: EOF array
		  pc  :: Principal Components array
		  ev  :: Eigen Values  array
		  ev_sum :: Sum of all eigen values (even thoses not returned).
		    Returned ONLY if get_ev_sum is True.
		    It can also be retreived with <SpAn_object>.stev_sum.
		:::
		"""

		# Parameters
		self._update_params_(**kwargs)
		print 'mssa params',self._nmssa,self._window,self._prepca

		# Check if we need PCA. If needed for one dataset, it is performed for ALL!
		if self._pre_pca_() and self._pca_raw_pc ==[]: # Still no PCA done
				self.pca(**kwargs)

		# Initializations of data
		for att in 'raw_eof','raw_pc','raw_ev','pairs','ev_sum':
			setattr(self,'_mssa_'+att,[])

		# Loop on datasets
		for iset,pdata in enumerate(self._pdata):

			if self._prepca[iset]: # Pre-PCA case
				if not self._pca_raw_pc.has_key(iset):
					self.pca(iset=iset,**kwargs)
				raw_eof, raw_pc, raw_ev, ev_sum = \
				  spanlib_fort.mssa(self._pca_raw_pc[iset].transpose(), self._window[iset], self._nmssa[iset])
			else: # Direct MSSA case
				raw_eof, raw_pc, raw_ev, ev_sum = \
				  spanlib_fort.mssa(pdata, self._window[iset], self._nmssa[iset])

			# Append results
			self._mssa_raw_pc.append(raw_pc)
			self._mssa_raw_eof.append(raw_eof)
			self._mssa_raw_ev.append(raw_ev)
			self._mssa_ev_sum.append(ev_sum)

			#FIXME:
			"""
			teof = MV.transpose(MV.reshape(ntsteof,(self.window,nspace[i],self._nmssa)))
			teof.id=self.varname[i]+'_mssa_steof'
			teof.standard_name='Empirical Orthogonal Functions'

			ax0=teof.getAxis(0)
			ax0.id='mssa_chan'
			ax0.standard_name='MSSA Channels'

			ax1=teof.getAxis(1)
			ax1.id='mssa_mode'
			ax1.standard_name='MSSA Modes in decreasing order'

			ax2=teof.getAxis(2)
			ax2.id='mssa_window'
			ax2.standard_name='MSSA Window Axis'

			tpc=MV.transpose(MV.array(ntstpc))
			tpc.id=self.varname[i]+'_stpc'
			tpc.standard_name='MSSA Principal Components'
			tpc.setAxis(0,ax0)

			ax3 = tpc.getAxis(1)
			ax3.id='time'
			ax3.standard_name = 'Time'

			tev = MV.array(ntstev,id=self.varname[i]+'_mssa_ev',axes=[ax0])
			tev.standard_name = 'MSSA Eigen Values'
			if relative:
				tev[:] = tev[:] * 100. / ntev_sum
				tev.units = '%'

			for var in teof,tpc,tev,ax0,ax1,ax2,ax3:
				var.name = var.id
				var.long_name = var.standard_name

			self._mssa_pc.append(ntstpc)
			self._mssa_eof.append(ntsteof)

			#self.pairs.append(getPairs(tpc))
			eof.append(teof)
			pc.append(tpc)
			ev.append(tev)
			self._mssa_ev_sum.append(ntev_sum)

		if len(eof)==1:
			#self.pairs = self.pairs[0]
			ret = [eof[0],pc[0],ev[0]]
			self._mssa_ev_sum = self._mssa_ev_sum[0]
		else:
			ret = [eof,pc,ev]

		if get_ev_sum:
			ret.append(self._pca_ev_sum)

		return ret
		"""
		self._last_analysis_type = 'mssa'
		gc.collect()



	def mssa_eof(self,update=False,pure=False,*args,**kwargs):
		"""Get EOFs from current MSSA decomposition

		If MSSA was not performed, it is done with all parameters sent to mssa()
		"""

		# Already available
		if self._mssa_fmt_pc != [] and not update:
			return self._mssa_fmt_eof

		# MSSA still not performed
		if self._mssa_raw_eof == [] or update:
			self.mssa(*args,**kwargs)
		del self._mssa_fmt_eof
		self._mssa_fmt_eof = []

		# Of, let's format the variable
		for iset,raw_eof in enumerate(self._mssa_raw_eof): # (window*nchan,nmssa)
			# Get raw data back to physical space
			if not self._prepca[iset]:
				self._mssa_fmt_eof.append(self._unstack_(iset,raw_eof.transpose(),
					(self._mode_axis_('mssa',iset),self._mssa_window_axis_(iset))))
			elif pure: # No pre-PCA performed
				self._mssa_fmt_eof.append([cdms.createVariable(raw_eof.transpose())])
				self._mssa_fmt_eof[-1][0].setAxisList(
					[self._mode_axis_('mssa',iset),self._mssa_channel_axis_(iset)])
			else:
				nm = self._nmssa[iset] ; nw = self._window[iset] ; nc = self._npca[iset]
				proj_eof,smodes = self._project_(iset,self._pca_raw_eof[iset],
					npy.swapaxes(raw_eof,0,1).reshape((nw*nc,nm),order='F'),nt=nw*nm)
				self._mssa_fmt_eof.append(self._unstack_(iset,proj_eof,
					(self._mode_axis_('mssa',iset),self._mssa_window_axis_(iset))))
			# Set attributes
			for idata,eof in enumerate(self._mssa_fmt_eof[-1]):
				if not self._stack_info[iset]['ids'][idata].find('variable_'):
					eof.id = self._stack_info[iset]['ids'][idata]+'_mssa_eof'
				else:
					eof.id = 'mssa_eof'
				eof.name = eof.id
				eof.standard_name = 'empirical_orthogonal_functions_of_mssa'
				eof.long_name = 'MSSA empirical orthogonal functions'
				atts = self._stack_info[iset]['atts'][idata]
				if atts.has_key('long_name'):
					eof.long_name += ' of '+atts['long_name']
				if atts.has_key('units'):
					del eof.units
		gc.collect()
		return self._return_(self._mssa_fmt_eof)

	def mssa_pc(self,update=False,*args,**kwargs):
		"""Get PCs from current MSSA decomposition

		If MSSA was not performed, it is done with all parameters sent to mssa()
		"""

		# Already available
		if self._mssa_fmt_pc != [] and not update:
			return self._mssa_fmt_pc

		# MSSA still not performed
		if self._mssa_raw_pc == [] or update:
			self.mssa(*args,**kwargs)
		del self._mssa_fmt_eof
		self._mssa_fmt_eof = []

		# Of, let's format the variable
		for iset,raw_pc in enumerate(self._mssa_raw_pc):
			idata = 0 # Reference is first data
			pc = cdms.createVariable(npy.asarray(raw_pc.transpose(),order='C'))
			pc.setAxis(0,self._mode_axis_('mssa',iset))
			print 'mssa raw pc',raw_pc.shape
			pc.setAxis(1,self._mssa_pctime_axis_(iset))
			if not self._stack_info[iset]['ids'][idata].startswith('variable_'):
				pc.id = self._stack_info[iset]['ids'][idata]+'_mssa_pc'
			else:
				pc.id = 'mssa_pc'
			pc.name = pc.id
			pc.standard_name = 'principal_components_of_mssa'
			pc.long_name = 'MSSA principal components'
			atts = self._stack_info[iset]['atts'][idata]
			if atts.has_key('long_name'): pc.long_name += ' of '+atts['long_name']
			if atts.has_key('units'):     pc.units = atts['units']
			self._mssa_fmt_pc.append(pc)

		return self._return_(self._mssa_fmt_pc)		
			

	def mssa_ev(self,relative=False,sum=False,update=False):
		"""Get eigen values from current MSSA decomposition

		Inputs:
		  relative :: Return percentage of variance
		  sum :: Return the sum of eigen values (total variance)
		"""

		# PCA still not performed
		if self._mssa_raw_ev == [] or update: self.mssa()

		if sum:
			res = self._ev_sum
		else:
			res = []
			for iset,raw_ev in enumerate(self._ev):
				if relative: 
					ev = cdms.createVariable(100.*raw_ev/self._ev_sum[iset])
				else:
					ev = cdms.createVariable(raw_ev)
				ev.setAxisList([self.mode_axis('mssa')])
				ev.name = ev.id
				ev.standard_name = 'eigen_values_of_mssa'
				ev.long_name = 'MSSA eigen values'
				atts = self._stack_info[iset]['atts'][idata]
				if atts.has_key('long_name'):
					ev.long_name += ' of '+atts['long_name']
				if relative:
					ev.units = '%'
				elif atts.has_key('units'):
					ev.units = atts['units']
					for ss in ['^','**',' ']:
						if ev.units.find(ss) != -1:
							ev.units = '(%s)2' % ev.units
							break
				res.append(ev)

		return self._return_(res)		


	def mssa_rec(self,imode=None,pure=False):
		
		# MSSA still not performed
		if self._mssa_raw_pc == [] or update: self.mssa(*args,**kwargs)
		del self._mssa_fmt_rec ; gc.collect()
		self._mssa_fmt_rec = []
		
		# Loop on datasets
		for iset,(raw_eof,raw_pc) in enumerate(zip(self._mssa_raw_eof,self._mssa_raw_pc)):
			raw_rec,smodes = self._project_(iset,raw_eof,raw_pc,imode).transpose() # (nt,nchan)
			if not self._prepca[iset]: # No pre-PCA performed
				self._mssa_fmt_rec.append(self._unstack_(iset,raw_rec,self._time_axis_(iset)))
			elif pure: # Force direct result from MSSA
				self._mssa_fmt_rec.append([cdms.createVariable(raw_rec)])
				self._mssa_fmt_rec[-1][0].setAxisList(0,
					[self._time_axis_(iset),self._mode_axis_('mssa',iset)])
			else: # With pre-pca
				proj_rec = self._project_(iset,self._pca_raw_eof[iset], raw_rec, 
					nt=self._window[iset]*self._nmssa[iset])
				self._mssa_fmt_eof.append(self._unstack_(iset,proj_rec,
					(self._mode_axis_('pca',iset),self._mssa_window_axis_(iset))))
			del  raw_rec
			# Set attributes
			for idata,rec in enumerate(self._mssa_fmt_rec[-1]):
				if not self._stack_info[iset]['ids'][idata].startswith('variable_'):
					rec.id = self._stack_info[iset]['ids'][idata]+'_mssa_rec'
				else:
					rec.id = 'mssa_rec'
				if imode is not None:
					rec.id += smodes #FIXME: do we keep it?
				rec.standard_name = 'recontruction_of_mssa_modes'
				rec.long_name = 'Reconstruction of MSSA modes'
				atts = self._stack_info[iset]['atts'][idata]
				if atts.has_key('long_name'):
					rec.long_name += ' of '+atts['long_name']
					
		return self._return_(self._mssa_fmt_rec,grouped=pure)	
	

	#################################################################
	## SVD
	#################################################################

	def svd(self,nsvd=None,pca=None):
		""" Singular Value Decomposition (SVD)

		Descriptions:::
		  This function performs a SVD
		  and returns EOF, PC and eigen values.
		  Unless pca parameter is set to false, a pre
		  PCA is performed to reduced the number of d-o-f
		  if already done and if the number of channels is
		  greater than 30.
		:::

		Usage:::
		eof, pc, ev = svd(nsvd,pca)

		  nsvd  :: Number of SVD modes retained
		  window :: MSSA window parameter
		  pca	:: If True, performs a preliminary PCA

		Output:::
		  eof :: EOF array
		  pc  :: Principal Components array
		  ev  :: Eigen Values  array

		:::
		"""

		# Check we have at least 2 variables!!
		# At the moment we will not use any more variable
		if len(self._pdata)<2:
			raise SpanlibError('svd','Error you need at least (most) 2 datasets to run svd, otherwise use pca and mssa')

		# Check for default values for mssa and pca if not passed by user
		if pca is None:
			if self._pca_raw_pc ==[] and max(self.ns) > 30: # Pre-PCA needed
				print '[svd] The number of valid points is greater than',30,' so we perform a pre-PCA'
				pca = True
			elif self._pca_raw_pc is not None:
				pca = True
			else:
				pca = False

		if pca is True: # From PCA to MSSA
			nspace = [self._npca,]*len(self._pdata)
			if self._pca_raw_pc ==[]: # Still no PCA done
				self.pca()
		else:
			nspace = self.ns

		if nsvd is not None:
			self._nsvd = nsvd


		if pca is True: # Pre-PCA case
			lneof, rneof, lnpc, rnpc, nev = \
			  spanlib_fort.svd(npy.transpose(self._pca_raw_pc[0]), 
			    npy.transpose(self._pca_raw_pc[1]), self._nsvd)
		else: # Direct SVD case
			lneof, rneof, lnpc, rnpc, nev = \
			  spanlib_fort.svd(self._pdata[0], self._pdata[1], self._nsvd)

		self._svd_eof = [lneof,rneof]
		self._svd_pc = [lnpc,rnpc]

		eof=[]
		pc=[]

		for i in range(2):
			teof = MV.transpose(self._svd_eof[i])
			teof.id=self.varname[i]+'_svd_eof'
			teof.standard_name='SVD Empirical Orthogonal Functions'

			ax0=teof.getAxis(0)
			ax0.id='svd_chan'
			ax0.standard_name='Channels'

			ax1=teof.getAxis(1)
			ax1.id=self.varname[i]+'_svd_mode'
			ax1.standard_name='SVD Modes in decreasing order'

			tpc=MV.transpose(MV.array(self._svd_pc[i]))
			tpc.id=self.varname[i]+'_svd_pc'
			tpc.standard_name='SVD Principal Components'
			tpc.setAxis(0,ax0)

			ax3 = tpc.getAxis(1)
			ax3.id='time'

			tev=MV.array(ntstev,id=elf.varname[i]+'_svd_ev',axes=[ax0])
			tev.standard_name='SVD Eigen Values'

			eof.append(teof)
			pc.append(tpc)

		return eof[0],pc[0],eof[1],pc[1],ev


	def _project_(self,iset,reof,rpc,imode=None,ns=None,nt=None,nw=None):
		"""Generic raw construction of modes for pure PCA, MSSA or SVD, according to EOFs and PCs, for ONE DATASET"""

		# Which modes
		imodemax = reof[iset].shape[-1]
		if imode is None:
			imode = range(1,imodemax+1)
		elif isinstance(imode,slice):
			imode = range(imode.start,imode.stop,imode.step)
		else:
			if isinstance(imode,int):
				if imode < 0:
					imode = range(-imode)
				else:
					imode = [imode-1,]
			imode = [im+1 for im in imode]

		# Rearrange modes (imode=[1,3,4,5,9] -> [1,1],[3,5],[9,9])
		imode = [im for im in imode if im > 0 and im < imodemax]
		imode.sort(_sort_modes_)
		imodes = []
		im = 0
		while im < len(imode):
			imode1 = imode2 = last_imode = imode[im]
			for imt in xrange(im+1,len(imode)):
				if imode[imt] > 0  and (abs(imode[imt])-imode2) > 1: # Group
					im = imt-1
					break
				imode2 = abs(imode[imt])
				continue
			else:
				im = imt
			im += 1
			imodes.append((imode1,imode2))

		# Function of reconstruction
		if nw is not None:
			function = spanlib_fort.mssa_rec # MSSA
		else:
			function = spanlib_fort.pca_rec  # PCA/SVD

		print 'imodes',imodes

		# Loop on datasets
##		ffrec=[]
##		for iset,pdata in enumerate(self._pdata):
		# Arguments
		args = [reof[iset],rpc[iset]]
		nmode = len(reof[iset])
		print '_proj',nmode
		kwargs = {}
		if nw is not None:
			# MSSA
			args.extend([ns[iset],nt,nw])
##		else:
			pass
##			# PCA
##			if ns is not None:
##				kwargs['ns'] = ns[iset]
##			else:
##				kwargs['ns'] = 1
			#FIXME: check that ok for ns=1
##			if kwargs['ns'] == 1: # Single space point
##				print '_reco',single
##				ffrec = self._pdata[iset]
##				ffrec.smodes = ''
##				return ffrec
			if nt is not None:
				kwargs['nt'] = nt
		# Loop modes
		smodes = []
		ffrec = 0.
		for j,ims in enumerate(imodes):
			if ims[0] > nmode: break
			if ims[1] > nmode: ims[1] = nmode
			print 'rec',ims
			args.extend(ims)
			print kwargs,'ims',ims,'reof',reof[iset].shape,'rpc',rpc[iset].shape
			ffrec += function(*args,**kwargs) # (nc,nt)
			if ims[0] == ims[1]:
				smode = str(ims[0])
			else:
				smode = '%i-%i'%tuple(ims)
			smodes.append(smode)
		print 'rec shape',ffrec[-1].shape
		return ffrec,'+'.join(smodes)



	def reconstruct(self,imode=None,mssa=None,pca=None,phases=False,nphases=8,offset=.5,firstphase=0,svd=None,ipair=None):
		""" Reconstruct results from mssa or pca

		Description:::
		  This function performs recontructions to retreive the
		  the contribution of a selection of modes to the original field.
		  By default, it recontructs from available PCA and MSSA
		  results. Recontruction of MSSA modes also calls recontruction
		  from of pre-PCA to get back to the original space.
		  This function can optionally performs phase composites
		  (useful for pairs of MSSA modes = oscillations) on MSSA
		  recontructions.
		:::

		Usage:::
		ffrec = reconstruct(imode,mssa,pca)

		  imode  :: Selection of modes [default: None]. If:
		    - None: all modes
		    - > 0: only this mode
		    - < 0: all modes until -imode
		    - list of modes: use it directly
		  ipair  :: Reconstruct pairs from MSSA if available. If:
		    - > 0: only this pair
		    - < 0: all pairs until ipair
		    - list of pairs: use it directly
		    It takes precedence over imode.
		  mssa   :: Reconstruct MSSA if True
		  pca    :: Reconstruct PCA if True
		  phases :: Operate phase reconstruction True/False (default is False)
		:::

		Output:::
		  ffrec :: Reconstructed field
		:::
		"""
		#print 'self._mssa_eof[0].shape',self._mssa_eof[0].shape

		# Which modes
		if imode is not None:
			if type(imode) is type(1):
				if imode < 0:
					imode = range(-imode)
				else:
					imode = [imode-1,]
			imode = (npy.array(imode)+1).tolist()

		# Which pairs (for MSSA)
		if mssa is True and ipair is not None:
			if ipair is type(1):
				if ipair < 0:
					ipair = range(-ipair)
				else:
					ipair = [ipair-1,]
			ipair = (npy.array(ipair)+1).tolist()

		ntimes=self._nt
		comments = 'Reconstructed from'
		axes=list(self.axes)

		# What we need explicitly
		#  - MSSA
		if mssa is True and self._mssa_eof == []:
			self.mssa()
		# - SVD
		if svd is True and self._svd_eof == []:
			self.svd()

		# What we have
		#  - SVD
		if svd is None:
			if self._svd_eof == []:
				svd = False
			elif self._mssa_eof==[]:
				svd = True
		#  - MSSA
		if mssa is None:
			if self._mssa_eof ==[]:
				mssa = False
			elif svd is False:
				mssa = True
			else:
				mssa = False
		# - PCA	
		if pca is None:
			pca = self._pca_raw_pc != []


		# Phase reconstruction for pca or mssa
		if phases and not pca and not mssa:
			raise 'Error you did not do any PCA or MSSA!\n To do a phases analysis only use the function %s in this module.\n%s' % ('computePhases',computePhases.__doc__)

		# MSSA reconstruction
		if mssa:
			comments+=' MSSA '
			# Space dimension
			if pca: # PCA (indirect MSSA)
				nspace = [self._npca,]*len(self._pdata[0])
			else:   # Physical (direct MSSA)
				nspace = [pdata.shape[0] for pdata in self._pdata]

			# Reconstructions
			if ipair is not None:
				# Reconstruct pairs
				if self.pairs == []:
					for i in xrange(len(self._pdata)):
						self.pairs.append(getPairs(MV.transpose(self._mssa_pc[i])))	
			else:
				# Reconstruct other modes
				if imode is None:
					imode=range(1,self._nmssa+1)
				print 'self._mssa_eof[0].shape',self._mssa_eof[0].shape
				ffrec = self._reconstruct(imode,self._mssa_eof,self._mssa_pc,
					ns=nspace,nt=self.nt,nw=self.window)

		# SVD Reconstruction
		if svd:
			comments+=' SVD '
			# Space dimension
			if pca: # PCA
				nspace = [self._npca,self._npca]
			else:   # Physical
				nspace = [pdata.shape[0] for pdata in self._pdata]

			# Reconstructions
			if imode is None:
				imode = range(1,self._nsvd+1)
			ffrec = self._reconstruct(imode,self._svd_eof,self._svd_pc,
				ns=nspace,nt=self.nt)

		# Phase composites reconstuction
		if phases:
			comments+=' Phases'
			if mssa:
				for i in xrange(len(self._pdata)):
					print 'yo',ffrec,computePhases(ffrec[i],nphases,offset,firstphase)
					ffrec[i] = computePhases(ffrec[i],nphases,offset,firstphase)
			else:
				ffrec=[]
				for i in xrange(len(self._pdata)):
					ffrec.append(computePhases(npy.transpose(self._pca_raw_pc[i]),
						nphases,offset,firstphase))

			# Replace time axis with phases axis
			print 'aaaaa',ffrec
			ntimes = nphases
			for iset in xrange(self._ndataset):
				for idata in xrange(self._ndata[iset]):
					if axes[j][i].isTime():
						axes[j][i]=ffrec[j].getAxis(1)
						break


		# PCA reconstruction (alone or for SVD or MSSA)
		if svd:
			ndataset = 2
		else:
			ndataset = self._ndataset
		if pca:
			comments+=' PCA'
			if True in [mssa, phases, svd]: # PCs are reconstructions themselves
				this_pc = [npy.transpose(ffrec[i]) for i in xrange(ndataset)]
				del(ffrec)
			else: # Classic pure PCA case
				this_pc = self._pca_raw_pc
			if mssa or imode is None: imode = range(1,self._npca+1)
			ffrec_raw = [self._reconstruct_modes(imode,self._pca_raw_eof[i],this_pc[i])
				for iset in xrange(ndataset)]


		# Format ouput data
		ffrec = []
		for iset in xrange(ndataset):

			# Unstack data for this dataset
			ffrec[iset] = self._unstack_(iset,ffrec_raw[iset],'pca')

			# Check some properties
			for this_rec in ffrec[iset]:
				if not ffrec[i].id.find('variable_'):
					this_rec.id = 'rec'
				else:
					this_rec.id += '_rec'
				this_rec.name = this_rec.id
				if self._stack_info[iset]['atts'].has_key('long_name'):
					this_rec.long_name = 'Reconstruction'
				else:
					this_rec.long_name = 'Reconstruction of '+this_rec.long_name
				this_rec.comments = comments

			# Back to physical space
			sh = [ffrec[i].shape[1],]
			sh.extend(self.shapes[i][1:])
			if self.mask[i] is not False:
				print ffrec[i].shape
				ffrec[i] = MV.transpose(spanlib_fort.chan_unpack(self.mask[i],ffrec[i],1.e20))
				ffrec.setMissing(1.e20)
				ffrec[i] = MV.reshape(ffrec[i],sh)
				ffrec[i] = MV.masked_object(ffrec[i],1.e20)
			else:
				ffrec[i] = MV.transpose(ffrec[i])
				ffrec[i] = MV.reshape(ffrec[i],sh)
			ffrec[i].setAxisList(axes[i])
			ffrec[i].id=self.varname[i]+'_rec'
			ffrec[i].name = ffrec[i].id
			for att in 'units','long_name':
				if att in self.attributes[i].keys():
					if att is 'long_name':
						ffrec[i].attributes[att] = \
						  'Reconstruction of '+self.attributes[i][att]
					else:
						ffrec[i].attributes[att] = self.attributes[i][att]
			ffrec[i].comment=comments
			ffrec[i].modes = imode
			if not svd:
				ffrec[i].setGrid(self.grids[i])

		return self._return_(ffrec)

	def _return_(self,dataset,grouped=False):
		"""Return dataset as input dataset"""
		# Single variable
		if self._input_map == 0:
#			if cdms.isVariable(dataset[0]):
			print dataset
			if not isinstance(dataset[0],list):
				return dataset[0]
			return dataset[0][0]
		# A single list of stacked variable (not serial mode)
		if isinstance(self._input_map,int):
			return dataset[0]
		# Full case
		for map,iset in enumerate(self._input_map):
#			if (map== 0 and not cdms.isVariable(dataset[iset])) or grouped:
			print isinstance(dataset[iset],list)
			if (map== 0 and isinstance(dataset[iset],list)) and not grouped:
				dataset[iset] = dataset[iset][0]
		gc.collect()
		return dataset

	def clean(self):
		"""(Re-)Initialization"""
		dicts = []
		for aa in 'pca','mssa','svd':
			dicts.append('_%s_ev_sum'%aa)
			for bb in 'raw','fmt':
				for cc in 'eof','pc','ev':
					dicts.append('_%s_%s_%s'%(aa,bb,cc))
		dicts.extend(['_mode_axes','_mssa_window_axes','_mssa_pctime_axes','_mssa_channel_axes','_svd_channel_axes'])
		lists = ['_mssa_pairs','_stack_info','_svd_l2r','_nt','_ns','_ndata','_pdata']
		for ll,func in (dicts,dict),(lists,list):
			for att in ll:
				if hasattr(self,att):
					obj = getattr(self,att)
					del obj
				setattr(self,att,func())
		self._ndataset = 0
		self._last_analysis_type = None
		gc.collect()


	def _stack_(self,dataset,dweights,dnorms):
		""" Takes several data files, of same time and stacks them up together

		Description:::
		This fonction concatenates several dataset that have the
		same time axis. It is useful for analysing for example
		several variables at the same time.
		It takes into account weights, masks and axes.
		:::

		Inputs:::
		dataset   :: A list of data objects.
			They must all have the same time length.
		dweights :: Associated weights.
		:::
		"""

		# Inits
		len_time=None
		taxes = []
		saxes = []
		atts = []
		shapes = []
		ids = []
		grids = []
		mvs = []
		masks = []
		norms = []
		types = []

		# Loop on datasets
		for idata,data in enumerate(dataset):

			# We must work on an cdms variable
			if not cdms.isVariable(data):
				data = cdms.createVariable(data,copy=0)
				if npy.core.ma.isMa(data):
					types.append('ma')
				else:
					types.append('numpy')
			else:
				types.append('cdms')

			# Check time
			if data.getTime() is not None:
				# If a proper time axis is found, bring it to front
				data = data(order='t...')
			if len_time is None:
				len_time = data.shape[0]
			elif len_time != len(data.getAxis(0)):
				raise 'Error all datasets must have the same time length!!!!'

			# Append info
			taxes.append(data.getAxis(0))
			saxes.append(data.getAxisList()[1:])
			shapes.append(data.shape)
			ids.append(data.id)
			atts.append({})
			for att in data.listattributes():
				atts[-1][att] = data.attributes[att]
			grids.append(data.getGrid())
			mvs.append(data.missing_value)

			# Pack 
			packed = _pack_(data,dweights[idata],dnorms[idata])

			# Create or concatenate
			if not idata:
				stacked = packed['data']
				weights = packed['weights']
			else:
				stacked = npy.concatenate((stacked,packed['data']))
				weights = npy.concatenate((weights,packed['weights']))
			norms.append(packed)
			masks.append(packed['mask'])
			del packed
			gc.collect()

		# Store data and information
		self._stack_info.append(dict(ids=ids,taxes=taxes,saxes=saxes,masks=masks,
			weights=weights,atts=atts,shapes=shapes,grids=grids,
			mvs=mvs,types=types,norms=norms))
		self._ndata.append(len(ids))
		self._pdata.append(stacked)
		self._nt.append(stacked.shape[0])
		if len(stacked.shape) == 2:
			self._ns.append(stacked.shape[1])
		else:
			self._ns.append(1)


	def _unstack_(self,iset,pdata,firstaxes):
		"""Return a list of variables in the physical space, for ONE dataset.

		firstaxes: MUST be a tuple
		"""

		## Expansion axis (first axis)
		#if pack_type == 'rec': # Normal time reconstruction
			#first_axis = self._stack_info[iset]['taxes']
		#elif pack_type == 'phase': # Time is phases
			#first_axis = pdata[0].getAxis(0)
		#else: # Reconstruction along pca, mssa or svd modes
			#first_axis = self._mode_axis(pack_type)

		# Loop on stacked data
		istart = 0
		unstacked = []
		if not isinstance(firstaxes,list):
			firstaxes = [firstaxes]
		for idata in xrange(len(self._stack_info[iset]['ids'])):

			# Get needed stuff
			for vname in self._stack_info[iset].keys():
				if vname[-1] == 's':
					exec "%s = self._stack_info[iset]['%s'][idata]" % (vname[:-1],vname)

			# Unpack data
			mlen = int(mask.sum())
			iend = istart + mlen
			unpacked = spanlib_fort.chan_unpack(mask.flat,pdata[istart:iend,:],mv)
			unpacked = npy.asarray(unpacked,order='C')

			# Check axes and shape
			axes = []
			for fa in firstaxes:
				if isinstance(fa,tuple):
					axes.append(fa[idata])
				elif isinstance(fa,dict):
					axes.append(fa[iset])
				else:
					axes.append(fa)
			if len(shape) > 1: axes.extend(saxe)
			shape = tuple([len(axis) for axis in axes])
			if unpacked.shape != shape:
				unpacked = unpacked.reshape(shape,order='C')
				
			# Mask and set attributes
			masked = MV.masked_object(unpacked,mv) ; del unpacked
			masked.setMissing(mv)
			masked.setAxisList(axes)
			masked.setGrid(grid)
			for attn,attv in att.items():
				setattr(masked,attn,attv)

			# Append to output
			unstacked.append(masked)
			istart += mlen

		gc.collect()
		return unstacked

	def rec(self,analysis_type=None,*args,**kwargs):
		if analysis_type is None:
			analysis_type = self._last_analysis_type
		else:
			valid = ['pca','mssa','svd']
			if analysis_type not in valid:
				raise SpanlibException('rec','analysis_type must be one of '+valid)
		if analysis_type is None:
			warnings.warn('Yet no statistics performed, so nothing to reconstruct!')
		else:
			return getattr(self,self._last_analysis_type+'_rec')(*args,**kwargs)

class SVDModel(SpAn):

	def __init__(self,data,**kwargs):

		SpAn.__init__(self,data,**kwargs)

		# Perform an SVD between the first two datasets
		self.svd(nsvd=None,pca=None)

		# Compute the scale factors between the two datasets
		self.scale_factors = npy.sqrt((npy.average(self._svd_pc[0]**2)/nsr - \
									 (npy.average(self._svd_pc[0])/nsr)**2) / \
									(npy.average(self._svd_pc[1]**2)/nsr - \
									 (npy.average(self._svd_pc[1])/nsl)**2))

	def __call__(self,data,nsvdrun=None,method='regre'):
		"""Run the SVD model 
		
		@keyparam method: Method of reconstruction [default: 'regre']. 'direct' assumes that left and normalized expansion coefficients are equal (Syu and Neelin 1995). 'regre' does not use right EOFs but regression coefficients (Harrisson et al 2002)
		"""

		if nsvdrun is not None:
			self._nsvdrun = nsvdrun

		#TODO: finish the svd model man !
		print 'missing code'


	def clean(self):
		SpAn.clean(self)
		self._regre_ = None
		gc.collect()



class SpanlibError(Exception):
	def __init__(self,where,what):
		Exception.__init__(self)
		self._where = where
		self._what = what
	def __str__(self):
		return 'SpanlibError: [%s] %s' % (self._where,self._what)



##	def __old_init__(self,data,weights=None,npca=10,window=None, nmssa=4,nsvd=10,relative=False):
##
##		## Sets all values to None
##		self.clean()
##
##		## Before all makes sure data is list of data
##		if not isinstance(data,(list,tuple)):
##			data=[data,]
##		if weights is None:
##			weights=[None,] * len(data)
##		elif not isinstance(weights,(list,tuple)):
##			weights = [weights,]
##
##		## First pack our data, prepare the weights and mask for PCA
##		self._pdata=[]
##		self.weights=[]
##		self.mask=[]
##		self.attributes=[]
##		self.shapes=[]
##		self.axes=[]
##		self.varname=[]
##		self.grids=[]
##		for i,d in enumerate(data):
##			d=MV.array(d,copy=0)
##			w=weights[i]
##			tmp = pack(d,w)
##			tmpdata,tmpweights,tmpmask = tmp
##			self._pdata.append(tmpdata)
##			self.weights.append(tmpweights)
##			self.mask.append(tmpmask)
##			self.attributes.append(d.attributes)
##			self.shapes.append(d.shape)
##			self.axes.append(d.getAxisList())
##			self.varname.append(d.id)
##			self.grids.append(d.getGrid())
##
##		# Space and time dimensions
##		self.nt = data[0].shape[0]
##		for d in data:
##			if d.shape[0] != self.nt:
##				raise Exception, 'Error your dataset are not all consistent in time length'
##		self.ns = [pdata.shape[0] for pdata in self._pdata]
##
##		# Number of modes
##		self._npca = npca
##		self._nmssa = nmssa
##		self._nsvd = nsvd
##
##		# MSSA window
##		if window is None:
##			self.window = int(self.nt/3.)
##		else:
##			self.window = window



#def stackData(*data):
	#""" Takes several data files, of same time and stacks them up together

	#Description:::
	  #This fonction concatenates several dataset that have the
	  #same time axis. It is useful for analysing for example
	  #several variables at the same time.
	  #It takes into account weights, masks and axes.
	#:::

	#Usage:::
	#dout, weights, mask, axes = stackData(data1[, data2...])

	  #*data   :: One or more data objects to stack.
	    #They must all have the same time length.
	#:::

	#Output:::
	  #dout	:: Stacked data
	  #weights :: Associated stacked weights
	  #masks   :: Associated stacked masks
	  #axes	:: Associated stacked axes
	#:::
	#"""
	#len_time=None
	#axes=[]
	#dout=None # data output
	#for d in data:
		#d = MV.array(d)
		#t=d.getTime()
		#if t is None:
			#t = d.getAxis(0)
			#t.designateTime()
			#d.setAxis(0,t)
		#if len_time is None:
			#len_time=len(t)
		#elif len_time!=len(t):
			#raise 'Error all datasets must have the same time length!!!!'

		#if d.getAxis(0)!=t:
			#d=d(order='t...')

		#axes.append(d.getAxisList())
		#tdata,w,m=pack(d)
		#if dout is None: # Create
			#dout=tdata
			#weights=w
			#masks=[m]
		#else: # Append
			#dout=npy.concatenate((dout,tdata))
			#weights=npy.concatenate((weights,w))
			#masks.append(m)
	#return npy.transpose(dout),weights,masks,axes

#def unStackData(din,weights,masks,axes):
	#""" Unstack data in the form returned from stackData

	#Description:::
	  #This function is the reverse operation of stakData.
	  #It splits stacked datasets into a list.
	#:::

	#Usage:::
	#dout = unStackData(din,weights,mask,axes)

	  #din	 :: Stacked data (see stackData function)
	  #weights :: Associated stacked weights
	  #masks   :: Associated stacked masks
	  #axes	:: Associated stacked axes
	#:::

	#Output:::
	  #dout	:: List of unstacked data
	#:::
	#"""
	#nvar=len(axes)

	#if nvar!=len(masks):
		#raise 'Error masks and input data length not compatible'

	#totsize=0
	#for m in masks:
		#totsize+=int(npy.sum(npy.ravel(m)))
	#if totsize!=din.shape[1]:
		#raise 'Error data and masks are not compatible in length!!!! (%s) and (%s)' \
		      #% (totsize,din.shape[1])

	#istart=0
	#dout=[]
	#missing_value = 1.e20
	#for i in xrange(nvar):
		#m=masks[i]
		#mlen=int(npy.sum(m.flat))
		#iend=istart+mlen
		#data=npy.transpose(din[:,istart:iend])
		#w=weights[istart:iend]
		##FIXME: missing value
		#up=spanlib_fort.chan_unpack(m,data,missing_value)
		#unpacked = MV.transpose(MV.masked_array(up))
		#sh = []
		#for ax in axes[i]:
			#sh.append(len(ax))
		#unpacked = MV.reshape(unpacked,sh)
		#unpacked.setMissing(missing_value)

		##unpacked = MV.masked_where(npy.equal(npy.resize(npy.transpose(m),unpacked.shape),0),unpacked,copy=0)
		#unpacked.setAxisList(axes[i])
		#istart+=mlen
		#dout.append(unpacked)
	#return dout

