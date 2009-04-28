#################################################################################
# File: spanlib_python.py
#
# This file is part of the SpanLib library.
# Copyright (C) 2006-2008  Stephane Raynaud, Charles Doutriaux
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

docs = dict(
	npca = """*npca* : int | ``None``
				Number of PCA modes to keep in analysis (defaults to 10).""", 
	prepca = """*prepca* : int | bool | ``None``
				Number of pre-PCA modes to keep before MSSA and SVD analyses (defaults to ``npca`` if ``True``).
				If the number of input channels is greater than 30, it is automatically switch to ``True``.""", 
	nmssa = """*nmssa* : int | ``None``
				Number of MSSA modes to keep in analysis (defaults to 10).""", 
	window = """*window* : int | ``None``
				Size of the MSSA window parameter (defaults to 1/3 the time length).""", 
	nsvd = """*nsvd* : int | ``None``
				Number of SVD modes to keep in analysis (defaults to 10).""", 
	iset = """*iset* : int | ``None``
				Set it to an integer to perform reconstruction on one dataset.
				If ``None``, all dataset are reconstructed.""", 
	modes = """*modes* : int | list | tuple
				If ``None``, all modes are summed. 
				Example of other usages:
				
					- ``4`` or ``[4]`` or ``(4,)``: only mode 4
					- ``-4``: modes 1 to 4
					- ``(1,3,-5)``: modes 1, 3, 4 and 5 (``-`` means "until")""", 
	raw = """*raw* : bool
				When pre-PCA is used and ``raw`` is ``True``, it prevents from
				going back to physical space (expansion to PCA EOFs space).""", 
	scale = """*scale* : bool | float
				Apply a factor to EOFs or PCs. This is essentially useful to add
				a quantitative meaning. If ``True``, ``scale`` is chosen so that
				the standard deviation of the mode (EOF or PC) is the square of
				the associated eigen value.""", 
	relative = """*relative* : bool
				Return percentage of variance instead of its absolute value.""", 
	sum = """*sum* : bool
				Return the sum of ALL (not only the selected modes) eigen values (total variance).""", 
	cumsum = """*cumsum* : bool
				Return the cumulated sum of eigen values.""", 
)

def _filldocs_(func):
	func.__doc__ = func.__doc__ % docs
	return func

def _pack_(data,weights=None,norm=None):
	""" Pack a dataset and its weights according to its mask

	:Description:
		This function packs a dataset in 2D space by removing
		all masked points and returning a space-time array.
		It performs this operation also on the weights.
		It is used for removing unnecessary points and
		simplifying the input format for analysis functions.


	:Parameters:
		data	: masked array
			Flatten in space an [x,y,t] array by 
	              removingits masked point
	========================================================
	
	========================================================
	Keyrwords     Description
	========================================================
	weights       Weights to be flatten also
	norm          Normalisation coefficients
	========================================================

	:Returns:
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
	if data.ndim == 1:
		nstot = 1
	else:
		nstot = npy.multiply.reduce(data.shape[1:])
	sh=list(data.shape)
	packed['norm'] = norm

	# Is it already packed? Check the mask...
	if data.mask is MV.nomask:
		if weights is None:
			weights = npy.ones(nstot,dtype='f')
		else:
			weights = npy.asarray(weights,dtype='f')
		packed['weights'] = weights.ravel()
		if len(packed['weights']) != nstot:
			packed['weights'] = npy.resize(packed['weights'], (nstot,))
		packed['data'] = data.filled().reshape((nt,nstot)).transpose()
		packed['mask'] = npy.ones(nstot)
		_check_norm_(packed)
		return packed

	# Weights ?
	if weights is None:
		if data.ndim == 3 and \
			data.getAxis(-1).isLongitude() and data.getAxis(-2).isLatitude():
			import cdutil
			weights = cdutil.area_weights(data[0]).raw_data() # Geographic weights
		else:
			weights = npy.ones(nstot,dtype='f').reshape(sh[1:])
	elif cdms.isVariable(weights):
		weights = weights.filled(0.)
	else:
		weights = npy.asarray(weights,dtype='f')

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
	packed['data'] = npy.asarray(spanlib_fort.chan_pack(data_to_pack,mask.flat,ns),dtype='f',order='F')
	del data_to_pack
	# - Pack weights
	weights = weights.reshape((1,nstot))
	packed['weights'] = npy.asarray(spanlib_fort.chan_pack(weights,mask.flat,ns)[:,0],dtype='f')

	_check_norm_(packed)
	return packed

def _get_imodes_(imode, nmode):
	
	# Which modes
	if imode is None:
		imode = range(1,nmode+1)
	elif isinstance(imode,slice):
		imode = range(imode.start,imode.stop,imode.step)
	else:
		if isinstance(imode,int):
			#if imode < 0:
				#imode = range(-imode)
			#else:
			imode = [imode,]
				#imode = [imode-1,]
		#imode = [im+1 for im in imode]

	# Rearrange modes (imode=[1,3,4,5,9] -> [1,1],[3,5],[9,9])
	imode = [im for im in imode if im <= nmode]
	imode.sort(cmp=lambda x,y: cmp(abs(x),abs(y)))
	if imode[0]<0: imode.insert(0,1)
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
		imodes.append((imode1,imode2))
	return imodes

def _check_length_(input,mylen,fillvalue):
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
#		for iadd in xrange(dlen):
			input.extend([fillvalue]*dlen)
	return input

def _check_norm_(packed):
	"""Setup the normalisation factor of a packed dataset"""
	if packed['norm'] is True or packed['norm'] is None:
		packed['norm'] = packed['data'].std() # Standard norm
		packed['data'] = packed['data']/packed['norm'] # Norm copy
	elif packed['norm'] is not False:
		if packed['norm'] <0: # Relative norm, else strict norm
			packed['norm'] = abs(packed['norm'])*packed['data'].std()
		packed['data'] = packed['data']/packed['norm'] # Norm copy
	else:
		packed['norm'] = 1.

def _get_pairs_(pcs, mincorr=0.85, idxtol=2, pertol=.1, smooth=5, evs=None, evtol=.1):
	""" Get pairs in MSSA results

	  This function detects pairs of mode in principal
	  components (typical from MSSA). Modes of a pair
	  have their PCs (and EOFs) in phase quadrature in time.
	  The algorithm detects quadratures using lag correlations.

	:Parameters:
		pcs : 
			Principal components [modes,time]
	  mincorr     :: Minimal correlation [default:0.95]
	  idxtol :: Maximal difference of indices between a pair of modes [default: 3]
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
	if evs is not None and len(evs) != pcs.shape[0]: evs = None
		
#	smooth = npy.clip(smooth,1,max([1,nst-3]))
	pairs = []
	found = []
	mincorr = npy.clip(mincorr,0.,1.)
	nt = len(pcs[0])
	idxtol = npy.clip(idxtol,1,len(pcs)-1)

	# Compute lagged autocorrelations to find dominant periods
	autocorrs = []
	idx = npy.arange(nt-2)
	for ip in xrange(len(pcs)):
		lcs = zip(*_next_max_corrs_(pcs[ip]))
		print lcs
		if not len(lcs):
			autocorrs.append(None)
		else:
			imax = npy.argmax(lcs[1])
		autocorrs.append((float(abs(lcs[0][imax])),lcs[1][imax]))
		print '%i : p=%g c=%g'%(ip, autocorrs[-1][0],autocorrs[-1][1])
			
	# Now the lagged cross correlations to check phase quadrature	
	for ip1 in xrange(len(pcs)-1):
		
		# Trend
		if autocorrs[ip1] is None: continue
		period = autocorrs[ip1][0]
		autocorr = autocorrs[ip1][1]
		
		# Mode already found
		if ip1 in found: continue
		
		# Scan next idxtol modes to find a possible pair
		for ip2 in xrange(ip1+1,min([ip1+idxtol+1,len(pcs)])):
			print ' %i-%i'%( ip1, ip2),

			# Lag correlations
			lcs = _next_max_corrs_(pcs[ip1], pcs[ip2])
			if len(lcs) < 2: continue
			(l1, c1), (l2, c2) = lcs
			
			# Check correlation strengths
			if max([c1, c2])/autocorr < mincorr:
				print ', bad max corr', max([c1, c2])/autocorr
				continue
			
			# Check closeness of periods
			if abs((l2-l1)-period)/period > pertol:
				print ', bad per', abs((l2-l1)-period)/period
				continue
			
			# Check closeness of eigen values 
			if evs is not None and 2*(evs[ip1]-evs[ip2])/(evs[ip1]+evs[ip2]) > evtol: 
				print ', bad ev', 2*(evs[ip1]-evs[ip2])/(evs[ip1]+evs[ip2])
				continue
			
			print 'We got a pair:', ip1, ip2
			break

		else:
			# No pair found
			continue
		found.extend([ip1, ip2])
		pairs.append((ip1, ip2))

	return pairs

def _next_max_corrs_(pc1, pc2=None, maxlag=None, getfirsts=True, decorr=False):
	"""left and right (-lag) max correlations"""
	if maxlag is None:
		maxlag = len(pc1)-2
	right_lags = range(0, maxlag)
	left_lags = range(0, -maxlag-1, -1)
	if pc2 is None:
		pc2 = pc1
	right_corrs = genutil.statistics.laggedcorrelation(pc1[:], pc2[:], right_lags)
	left_corrs = genutil.statistics.laggedcorrelation(pc1[:], pc2[:], left_lags)
	maxima = ()
	for lags, corrs in (left_lags, left_corrs), (right_lags, right_corrs):
		dcorr = npy.diff(corrs)
		mask = (dcorr[:-1]>0) & (dcorr[1:]<0)
		if npy.any(mask):
			mylags = npy.compress(mask, lags[1:])
			mycorrs = corrs[1:].compress(mask)
			if getfirsts:
				maxima += ((mylags[0], mycorrs[0]), )
			else:
				maxima += ((mylags, mycorrs), )
	return maxima



def get_phases(data,nphase=8,minamp=.5,firstphase=0):
	""" Phase composites for oscillatory fields

	  This computes temporal phase composites of a spatio-temporal
	  dataset. The dataset is expected to be oscillatory in time.
	  It corresponds to a reoganisation of the time axis to
	  to represents the dataset over its cycle in a arbitrary
	  number of phases. It is useful, for example, to have a
	  synthetic view of an reconstructed MSSA oscillation.

	:Parameters:
		data	 : array
			A ``cdms2`` variable with a time axis over which
			composites are computed.
		nphase : int
			Number of phases (divisions of the cycle)
		minamp : float 
			Minimal value of retained data, relative to standard deviation.
		firstphase : float
			Position of the first phase in the 360 degree cycle.

	:Returns:
		A ``cdms2`` variable with phase axis (of length ``nphase``)
		instead of a time axis.
	"""
	
	# Get the first PC and its smoothed derivative
	pc = SpAn(data).pca_pc(npca=1, quiet=True)[0]
	pc[:] /= pc.std()
	dpc = npy.gradient(pc)
	if len(dpc) > 2: # 1,2,1 smooth
		dpc[1:-1] = npy.convolve(dpc, [.25, .5, .25])
	dpc[:] = dpc/dpc.std()
	
	# Get amplitude and phase indexes
	amplitudes = npy.hypot(pc, dpc)
	angles = npy.arctan2(dpc, pc)
	dphase = 2*npy.pi/nphase
	angles[:] = npy.where(npy.greater_equal(angles, 2*npy.pi-dphase*.5), angles-dphase*.5, angles)
	
	# Selection according to amplitudes 
	good = npy.greater_equal(amplitudes, minamp)
	pc = npy.compress(good, pc, axis=0)
	dpc = npy.compress(good, dpc, axis=0)
	itaxis = npy.clip(data.getOrder().find('t'), 0, data.ndim-1)
	cdata = data.compress(good, axis=itaxis)
	
	# Initialize output variable
	sl = [slice(None)]*data.ndim
	sl[itaxis] = slice(0, 1)
	phases = MV.repeat(MV.take(data, (0, ), itaxis), nphase, itaxis)
	phases[:] = MV.masked
	paxis = cdms.createAxis(npy.arange(nphase)*dphase, id='phases')
	paxis.long_name = 'Circular phases'
	paxis.units = 'degrees'
	axes = data.getAxisList()
	axes[itaxis] = paxis
	phases.setAxisList(axes)
	
	# Loop on circular bins to make composites
	marks = dphase * (npy.arange(nphase) - .5) + firstphase*npy.pi/360.
	slices = [slice(None), ]*phases.ndim
	idx = npy.arange(data.shape[itaxis])
	for iphase in xrange(len(marks-1)):
		slices[itaxis] = iphase
		inbin = npy.logical_and(npy.geater_equal(angles, marks[iphase]), npy.less(angles, marks[iphase+1]))
		phases[tuple(slices)] = MV.average(data.compress(inbin, axis=itaxis))

	return phases


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

	def __init__(self, datasets, serial=False, weights=None, norms=True, quiet=False, **kwargs):
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
		  serial :: If we have a list (or tuple) of variables as "datasets", they are analysed independantly (serial analyses in opposition to parallel analyses) if serial is True [default: False]. If False, they are packed before being analysed.
		:::

		Output:::
		  analysis_object :: SpAn object created for further analysis
		:::
		"""

		self._quiet = quiet
		self.clean()

		# We start from a list
		if not isinstance(datasets,(list,tuple)):
			datasets = [datasets]
			self._input_map = 0
		else:
			if isinstance(datasets,tuple):
				datasets = list(datasets)
			self._input_map = len(datasets)

		# Check if we are forced to be in serial analysis mode
		for d in datasets:
			if isinstance(d,(list,tuple)):
				serial = True
				break
		if not serial: # Convert to serial like mode
			datasets = [datasets]
		else:
			input_map = []
			for iset in xrange(len(datasets)):
				if isinstance(datasets[iset],(list,tuple)):
					if isinstance(datasets[iset],tuple):
						datasets[iset] = list(datasets[iset])
					input_map.append(len(datasets[iset]))
				else:
					datasets[iset] = [datasets[iset]]
					input_map.append(0)
			self._input_map = input_map

		# Weights and norms in the same form as datasets
		norms = kwargs.pop('norm', norms)
		weights = kwargs.pop('weight', weights)
		if norms is None:	norms = 1.
		norms = self._check_shape_(norms, True)
		for iset,d in enumerate(datasets):
			if len(d) == 1: norms[iset] = [1.]
		weights = self._check_shape_(weights, None)

		# We stack and pack data
		for iset,d in enumerate(datasets):
			self._stack_(d,weights[iset],norms[iset])
		self._ndataset = len(datasets)
		self._ndata= [len(dataset) for dataset in datasets]


		# Check and save parameters
		self._update_(None, **kwargs)
#		self._update_(npca=npca,window=window,nmssa=nmssa,nsvd=nsvd,prepca=prepca)


	#################################################################
	## Get datasets info
	#################################################################

	def _time_axis_(self,iset,idata=None):
		"""Get the time axis of data variable of a dataset"""
		if idata is None:
			return self._stack_info[iset]['taxes']
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
		single = False
		if analysis_type == 'svd':
			isets = [0]
		elif isets is None: 
			isets = xrange(self._ndataset)
		elif not isinstance(isets,(list,tuple)):
			isets = [isets]
			single = True
		out = []
		for iset in isets:
			nn = getattr(self,'_n'+analysis_type)
			if isinstance(nn, list): nn = nn[iset]
			if not self._mode_axes[analysis_type].has_key(iset) or \
				len(self._mode_axes[analysis_type][iset]) != nn:
				self._mode_axes[analysis_type][iset] = cdms.createAxis(npy.arange(1,nn+1))
				self._mode_axes[analysis_type][iset].id = analysis_type+'_mode'
				self._mode_axes[analysis_type][iset].long_name = analysis_type.upper()+' modes in decreasing order'
				self._check_dataset_tag_('_mode_axes',iset,analysis_type, svd=analysis_type=='svd')
			if analysis_type == 'svd': return self._mode_axes[analysis_type][iset]
			out.append(self._mode_axes[analysis_type][iset])
		if single: return out[0]
		return out

	def _channel_axis_(self,iset, name, **kwargs):
		"""Get the channel axis for one dataset (MSSA or SVD)"""
		if not self._prepca[iset]:
			nchan = self._ns[iset]
		else:
			nchan = self._prepca[iset]
		channel_axes = getattr(self,'_%s_channel_axes'%name)
		if not channel_axes.has_key(iset) or len(channel_axes[iset]) != nchan:
			channel_axes[iset] = cdms.createAxis(npy.arange(nchan))
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
		"""Get the MSSA window axis for one dataset"""
		if not self._mssa_window_axes.has_key(iset) or len(self._mssa_window_axes[iset]) != self._window[iset]:
			self._mssa_window_axes[iset] = cdms.createAxis(npy.arange(self._window[iset]))
			self._mssa_window_axes[iset].id = 'mssa_window'
			self._mssa_window_axes[iset].long_name = 'MSSA window time'
			self._check_dataset_tag_('_mssa_window_axes',iset)
		return self._mssa_window_axes[iset]

	def _mssa_pctime_axis_(self,iset,idata=0):
		"""Get the MSSA PCs time axis for one dataset"""
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
		
	def _check_dataset_tag_(self,name,iset,key=None,long_name=True,id=True, svd=False):
		"""Mark some attributes as specific to a dataset (only if there are more then one dataset)
			iset:: ID of the dataset
			key:: A dictionary key to select the dataset [default: None]
			long_name:: Mark the long name [defualt: True]
			id:: Mark the id [default: True]
			svd: Mark using 'left' or 'right'
		"""
		if self._ndataset > 1:
			targetset = getattr(self,name)
			if key is not None:
				targetset = targetset[key]
			target = targetset[iset]
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

	def _norm_(self,iset,idata):
		"""Get the normalization factor of one subdataset of one dataset"""
		return self._stack_info[iset]['norms'][idata]
		
	def norms(self, iset, idata=None):
		"""Get the normalization factor of a dataset"""
		if idata is None:
			return self._stack_info[iset]['norms']
		return self._stack_info[iset]['norms'][idata]
		
	def _type_(self,iset,idata):
		"""Get 'numpy', 'MA' or 'cdms' for one subdataset of one dataset"""
		return self._stack_info[iset]['types'][idata]

#	def _changed_param_(self,old,param):
#		"""Check if a parameter is changed for all datasets.
#		Return a list of booleans.
#		@param old: Dictionary of name and value of parameters
#		@param param: Parameter name.
#		"""
#		if isinstance(old[param],(list,tuple)):
#			return [old[param][iset] != getattr(self,'_'+param)[iset] for iset in xrange(self._ndataset)]
#		else:
#			return old[param] != getattr(self,'_'+param)
		
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
		init_all = [None]*self._ndataset
		for param in self._all_params:
			#print 'check', param
#			if not req_params: continue
			# Get old values , defaults to None and set new value
#			if kwargs.has_key(param):
			if param == 'nsvd':  # Single value for all datasets
				changed[param] = False
				old[param] = getattr(self,'_'+param,None)
				setattr(self,'_'+param,_check_length_(kwargs.pop(param, old[param]),0,None))
			else:
				changed[param] = [False]*self._ndataset
				old[param] = getattr(self,'_'+param,list(init_all))
				setattr(self,'_'+param,_check_length_(kwargs.pop(param, old[param]),self._ndataset,None))
#				print 'cur', param, getattr(self, '_'+param), changed[param]
#		if not req_params: return changed
		if not running: return changed

		# Number of PCA modes		
#		if 'npca' in req_params:
		for iset in xrange(self._ndataset):
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
			self._npca[iset] = npy.clip(self._npca[iset],1,min(SpAn._npca_max,self._ns[iset],self._nt[iset])) # Max
			
		# Number of pre-PCA modes before MSSA and SVD
#		if 'prepca' in req_params:
		for iset in xrange(self._ndataset):
			if self._prepca[iset] is None: # Default: pre-PCA needed over max (for MSSA and SVD)
				self._prepca[iset] = min(self._ns[iset], self._ns[iset]) > SpAn._npca_max
				if not self._quiet and self._prepca[iset] and analysis_type in ['mssa', 'pca'] :
					print '[mssa/svd] The number of valid points of one of the datasets is greater than %i, so we perform a pre-PCA'%SpAn._npca_max
			if self._prepca[iset] is True: # Defaults to the number of PCA modes
				self._prepca[iset] = self._npca[iset]
			elif self._prepca[iset]: # Max number of prepca modes is number of points
				self._prepca[iset] = min(self._prepca[iset], self._ns[iset], self._nt[iset])
			if self._prepca[iset] == 0:
				self._prepca[iset] = False
			
		# Dependency rules between prepca and npca
		for iset in xrange(self._ndataset):
			if self._prepca[iset] and self._npca[iset] < self._prepca[iset]:
				if not self._quiet  and self._prepca[iset]:
						print 'The number of pre-PCA modes (%i) for dataset #%iis lower than the number of PCA modes (%i), so we adjust the latter.' % (self._prepca[iset],iset,self._npca[iset])
				self._npca[iset] = self._prepca[iset]
			
		# Window extension of MSSA
#		if 'window' in req_params:
		for iset in xrange(self._ndataset):
			if self._window[iset] is None: # Initialization
				self._window[iset] = int(self._nt[iset]*SpAn._window_default)
			self._window[iset] = npy.clip(self._window[iset],1,max(1,self._nt[iset]))
			
		# Number of MSSA modes
		for iset in xrange(self._ndataset):
#			if 'nmssa' not in req_params and not changed['prepca'][iset]: continue
#			if not changed['prepca'][iset]: continue
			if self._nmssa[iset] is None: # Initialization
				# Guess a value
				if iset:
					self._nmssa[iset] = self._nmssa[iset-1] # From last dataset
				else:
					self._nmssa[iset] = SpAn._nmssa_default # Default value
			if self._prepca[iset]:
				nchanmax = self._prepca[iset] # Input channels are from pre-PCA
			else:
				nchanmax = self._ns[iset] # Input channels are from real space
			self._nmssa[iset] = npy.clip(self._nmssa[iset],1,
				min(SpAn._nmssa_max,nchanmax*self._window[iset])) # Max
			
		# Number of SVD modes (special case)
		if self._nsvd is None: # Initialization
			self._nsvd = SpAn._nsvd_default # Default value
		for iset in xrange(self._ndataset): # Check values
#			if 'nsvd' not in req_params and not changed['prepca'][iset]: continue
			if not changed['prepca'][iset]: continue
			if self._prepca[iset]:
				nchanmax = self._prepca[iset] # Input channels are from pre-PCA
			else:
				nchanmax = self._ns[iset] # Input channels are from real space
			self._nsvd = npy.clip(self._nsvd,1, min(SpAn._nsvd_max,nchanmax)) # Max
			
#		# Check what changed
#		for param in self._all_params:
#			changed[param] = self._changed_param_(old,param)
			
		# Re-run analyses when needed
#		if not kwargs: return changed # Just initializations (dry run to prevent not ending loop)
		changed['nsvd'] = old['nsvd'] != self._nsvd
		runsvd = False
		for iset in xrange(self._ndataset):
			
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
#			print 'Rerunning SVD'
			self.svd()
				
		# Inform about which params have been modified for each dataset
		return changed

	def _check_isets_(self,iset):
		"""Check if an iset is a valid dataset.
		It can be a list, and it is returned as a list.
		if an iset is invalid, it is removed from the output list.
		"""
		if iset is None: return range(self._ndataset)
		if iset == 'left':
			iset = 0
		elif iset == 'right':
			iset = 1
		if iset < 0 or iset >= self._ndataset:
			warn('Invalid dataset id: %i. Valid id are < %i'%(iset,self._ndataset))
		else:
			return [iset]

	def _check_shape_(self, inputs, fillvalue):
		"""Return input as datasets (tree) *shape*"""
		imap = self._input_map
		if isinstance(imap, int):
			return [_check_length_(inputs, max(1, imap), fillvalue), ]
		inputs = _check_length_(inputs,len(imap),fillvalue)
		for iset, im in enumerate(imap):
			inputs[iset] = _check_length_(inputs[iset], max(1, im), fillvalue)
		return inputs


	#################################################################
	## PCA
	#################################################################
	@_filldocs_
	def pca(self,iset=None,**kwargs):
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
			if getattr(self,'_pca_raw_pc').has_key(iset) and \
				getattr(self,'_pca_raw_pc')[iset].shape[-1] > self._npca[iset]:
				continue
			
			# Remove old results
			for att in 'raw_eof','raw_pc','raw_ev','ev_sum':
				dic = getattr(self,'_pca_'+att)
				if dic.has_key(iset): del dic[iset]

			# Compute PCA
			pdata = self._pdata[iset]
			if pdata.ndim == 1: # One single channel, so result is itself
				raw_eof = npy.ones(1,dtype=pdata.dtype)
				raw_pc = pdata
				raw_ev = raw_pc.var()
				ev_sum = ev

			else: # Several channels
				weights = self._stack_info[iset]['weights']
				raw_eof,raw_pc,raw_ev,ev_sum = spanlib_fort.pca(pdata,self._npca[iset],weights,-1)	

			# Save results
			self._pca_raw_pc[iset] = raw_pc
			self._pca_raw_eof[iset] = raw_eof
			self._pca_raw_ev[iset] = raw_ev
			self._pca_ev_sum[iset] = ev_sum
			
			# Delete formmated variables
			for vtype in 'pc', 'eof':
				vfmt = getattr(self, '_pca_fmt_'+vtype)
				if vfmt.has_key(iset): del vfmt[iset]
			gc.collect()

		self._last_analysis_type = 'pca'


	@_filldocs_
	def pca_eof(self, iset=None, scale=False, **kwargs):
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
			if self._pca_fmt_eof.has_key(iset):
				fmt_eof[iset] = self._pca_fmt_eof[iset]
				continue
				
			# First PCA analysis?
			if not self._pca_raw_eof.has_key(iset): self.pca(iset=iset)
				
			# Get raw data back to physical space
			self._pca_fmt_eof[iset] = \
				self._unstack_(iset,self._pca_raw_eof[iset][:, :self._npca[iset]],self._mode_axis_('pca',iset))
			
			# Set attributes and scale
			for idata,eof in enumerate(self._pca_fmt_eof[iset]):
				
				# Attributes
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

		return self._return_(fmt_eof)		


	@_filldocs_
	def pca_pc(self,iset=None, **kwargs):
		"""Get PCs from current PCA decomposition
		
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
			if self._pca_fmt_pc.has_key(iset):
				fmt_pc[iset] = self._pca_fmt_pc[iset]
				continue
			
			# First PCA analysis?
			if not self._pca_raw_pc.has_key(iset): self.pca(iset=iset)

			# Format the variable
			pc = cdms.createVariable(npy.asarray(self._pca_raw_pc[iset][:,:self._npca[iset]].transpose(),order='C'))
			pc.setAxis(0,self._mode_axis_('pca',iset))
			pc.setAxis(1,self._time_axis_(iset,0))
			pc.id = pc.name = 'pca_pc'
			pc.standard_name = 'principal_components_of_pca'
			pc.long_name = 'PCA principal components of '
			atts = self._stack_info[iset]['atts'][0]
			if self._ndata[iset] == 1 and  atts.has_key('long_name'): 
				pc.long_name += atts['long_name']
#			else:
#				pc.long_name += 'dataset %i'%iset
			if (self._ndata[iset] == 1 or npy.allclose(self.norms(iset), 1.)) and atts.has_key('units'):
				pc.units = atts['units']
				
			
			fmt_pc[iset] = self._pca_fmt_pc[iset] = pc
			self._check_dataset_tag_('_pca_fmt_pc',iset)

		return self._return_(fmt_pc)		


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

			# Format the variable
			id = 'pca_ev'
			long_name = []
			raw_ev = self._pca_raw_ev[iset][:self._npca[iset]]
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
			ev.long_name = ' '.join(long_name).title()+' of '
			ev.setAxisList([self._mode_axis_('pca',iset)])
			ev.standard_name = 'eigen_values_of_pca'
			atts = self._stack_info[iset]['atts'][0]
			if self._ndata[iset] == 1 and atts.has_key('long_name'):
				ev.long_name += atts['long_name']
			else:
				ev.long_name += 'dataset %i'%iset
			if relative:
				ev.units = '% of total variance'
			elif (self._ndata[iset] == 1 or npy.allclose(self.norms(iset), 1.)) and atts.has_key('units'):
				ev.units = atts['units']
				for ss in ['^','**',' ']:
					if ev.units.find(ss) != -1:
						ev.units = '(%s)^2' % ev.units
						break
			res[iset] = ev

		return self._return_(res)		

	@_filldocs_
	def pca_rec(self,iset=None,modes=None,**kwargs):
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

		# Loop on datasets
		pca_fmt_rec = {}
		for iset in isets:
			
			# First PCA analysis?
			if not self._pca_raw_pc.has_key(iset): self.pca(iset=iset)

			# Get raw data back to physical space
			reof = self._pca_raw_eof[iset][:,:self._npca[iset]]
			rpc = self._pca_raw_pc[iset][:,:self._npca[iset]]
			raw_rec,smodes = self._project_(reof,rpc,iset,modes)
			pca_fmt_rec[iset] = self._unstack_(iset,raw_rec,self._time_axis_(iset))
			del  raw_rec
			
			# Set attributes and scale
			for idata,rec in enumerate(pca_fmt_rec[iset]):
#				rec[:] *= self._norm_(iset,idata) # Scale
				if not self._stack_info[iset]['ids'][idata].startswith('variable_'):
					rec.id = self._stack_info[iset]['ids'][idata]+'_pca_rec'
				else:
					rec.id = 'pca_rec'
				rec.name = rec.id
#				if modes is not None:
				rec.id += smodes
#				else:
#					rec.modes = '1-%i'%self._npca[iset]
				rec.modes = smodes
				rec.standard_name = 'recontruction_of_pca_modes'
				rec.long_name = 'Reconstruction of PCA modes: '+smodes
				atts = self._stack_info[iset]['atts'][idata]
				if atts.has_key('long_name'):
					rec.long_name += ' of '+atts['long_name']
				if atts.has_key('units'):
					rec.units = atts['units']
					
		return self._return_(pca_fmt_rec)	
	
	

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
		for iset,pdata in enumerate(self._pdata):

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
				raw_eof, raw_pc, raw_ev, ev_sum = \
				  spanlib_fort.mssa(self._pca_raw_pc[iset][:, :self._prepca[iset]].transpose(), 
				  self._window[iset], self._nmssa[iset])
				  
			else: # Direct MSSA case
				raw_eof, raw_pc, raw_ev, ev_sum = \
				  spanlib_fort.mssa(pdata, self._window[iset], self._nmssa[iset])

			# Save results
			self._mssa_raw_pc[iset] = raw_pc
			self._mssa_raw_eof[iset] = raw_eof
			self._mssa_raw_ev[iset] = raw_ev
			self._mssa_ev_sum[iset] = ev_sum

			# Delete formmated variables
			for vtype in 'pc', 'eof':
				vfmt = getattr(self, '_mssa_fmt_'+vtype)
				if vfmt.has_key(iset): del vfmt[iset]
				
		self._last_analysis_type = 'mssa'
		gc.collect()

	@_filldocs_
	def mssa_eof(self, iset=None, scale=False, raw=False,**kwargs):
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
			if self._mssa_fmt_eof.has_key(iset):
				fmt_eof[iset] = self._mssa_fmt_eof[iset]
				continue
				
			# No analyses performed?
			if not self._mssa_raw_eof.has_key(iset): self.mssa(iset=iset)
			raw_eof = self._mssa_raw_eof[iset]
			nm = self._nmssa[iset] ; nw = self._window[iset]
			nl = self._mssa_raw_eof[iset].shape[0]/nw
			raw_eof = self._mssa_raw_eof[iset].reshape((nl, nw, nm))
			
			if raw: # Do not go back to physical space
				self._mssa_fmt_eof[iset] = [cdms.createVariable(raw_eof.transpose())]
				self._mssa_fmt_eof[iset][0].setAxisList(
					[self._mode_axis_('mssa',iset),self._mssa_channel_axis_(iset)])
					
			else: # Get raw data back to physical space
#				raw_eof = raw_eof.swapaxes(1,2).reshape((nl, nm*nw),order='F')
				raw_eof = raw_eof.reshape((nl, nm*nw),order='F')
				
				if not self._prepca[iset]: # No pre-PCA performed
					self._mssa_fmt_eof[iset] = self._unstack_(iset,raw_eof, 
						(self._mode_axis_('mssa',iset),self._mssa_window_axis_(iset)))
						
				else:
					proj_eof,smodes = self._project_(self._pca_raw_eof[iset],raw_eof.transpose(),iset, nt=nm*nw)
	#					npy.swapaxes(raw_eof,0,1).reshape((nw*nc,nm),order='F'),iset, nt=nw*nm)
					self._mssa_fmt_eof[iset] = self._unstack_(iset,proj_eof,
						(self._mode_axis_('mssa',iset),self._mssa_window_axis_(iset)))
					
			# Set attributes
			for idata,eof in enumerate(self._mssa_fmt_eof[iset]):
				if not raw and not self._stack_info[iset]['ids'][idata].find('variable_'):
					eof.id = self._stack_info[iset]['ids'][idata]+'_mssa_eof'
				else:
					eof.id = 'mssa_eof'
				eof.name = eof.id
				eof.standard_name = 'empirical_orthogonal_functions_of_mssa'
				eof.long_name = 'MSSA empirical orthogonal functions'
				if not raw:
					atts = self._stack_info[iset]['atts'][idata]
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
		return self._return_(fmt_eof)

	@_filldocs_
	def mssa_pc(self, iset=None, **kwargs):
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
			if self._mssa_fmt_pc.has_key(iset):
				fmt_pc[iset] = self._mssa_fmt_pc[iset]
				continue
				
			# No analyses performed?
			if not self._mssa_raw_pc.has_key(iset): self.mssa(iset=iset)
						
			# Format the variable
			pc = cdms.createVariable(npy.asarray(self._mssa_raw_pc[iset][:,:self._nmssa[iset]].transpose(),order='C'))
			pc.setAxis(0,self._mode_axis_('mssa',iset))
			pc.setAxis(1,self._mssa_pctime_axis_(iset))
			pc.id = pc.name = 'mssa_pc'
			pc.standard_name = 'principal_components_of_mssa of '
			pc.long_name = 'MSSA principal components'
			atts = self._stack_info[iset]['atts'][0]
			if self._ndata[iset] == 1 and atts.has_key('long_name'): 
				pc.long_name += atts['long_name']
#			else:
#				pc.long_name += 'of dataset %i'%iset
			if (self._ndata[iset] == 1 or npy.allclose(self.norms(iset), 1.)) and atts.has_key('units'):
				pc.units = atts['units']
			fmt_pc[iset] = self._mssa_fmt_pc[iset] = pc
			self._check_dataset_tag_('_mssa_fmt_pc',iset)

		return self._return_(fmt_pc,grouped=True)
			

	@_filldocs_
	def mssa_ev(self,iset=None,relative=False,sum=False,cumsum=False,**kwargs):
		"""Get eigen values from current MSSA decomposition

		:Parameters:
		  %(relative)s
		  %(sum)s
		  %(cumsum)s
		  %(iset)s
		
			
		:MSSA parameters:
			%(nmssa)s
			%(window)s
			%(prepca)s
			
		:Returns:
			Arrays with shape ``(nmssa,)`` or a float
		"""

		# Check on which dataset to operate
		isets = self._check_isets_(iset)

		# Update params
		self._update_('mssa', **kwargs)

		# Loop on dataset
		res = {}
		for iset in isets:
			
			# No analyses performed?
#			if not self._pca_raw_eof.has_key(iset): self.pca(iset=iset)
			if not self._mssa_raw_eof.has_key(iset): self.mssa(iset=iset)

			# We only want the sum
			if sum:
				res[iset] = self._mssa_ev_sum[iset]
				continue

			# Format the variable
			id = 'mssa_ev'
			long_name = []
			raw_ev = self._mssa_raw_ev[iset][:self._nmssa[iset]]
			if cumsum:
				raw_ev = raw_ev.cumsum()
				id += '_cumsum'
				long_name.append('cumulative')
			if relative: 
				raw_ev = 100.*raw_ev/self._mssa_ev_sum[iset]
				id += '_rel'
				long_name.append('relative')
			ev = cdms.createVariable(raw_ev)
			ev.id = ev.name = id
			long_name.append('MSSA eigen values')
			ev.long_name = ' '.join(long_name).title()
			ev.setAxisList([self._mode_axis_('mssa',iset)])
			ev.standard_name = 'eigen_values_of_mssa'
			atts = self._stack_info[iset]['atts'][0]
			if atts.has_key('long_name'):
				ev.long_name += ' of '+atts['long_name']
			if relative:
				ev.units = '% of total variance'
			res[iset] = ev

		return self._return_(res)		


	@_filldocs_
	def mssa_rec(self, iset=None, modes=None, raw=False, phases=False, **kwargs):
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
		
		# Loop on datasets
		mssa_fmt_rec = {}
		for iset in isets:

			# No analyses performed?
			if not self._mssa_raw_eof.has_key(iset): self.mssa(iset=iset)
			
			# Projection
			raw_rec,smodes = self._project_(self._mssa_raw_eof[iset],
				self._mssa_raw_pc[iset],iset,modes,nw=self._window[iset])
				
			# Phases composites
			if phases:
				raw_rec = get_phases(raw_rec,phases)
				taxis = raw_rec.getAxis(1)
			else:
				taxis = self._time_axis_(iset)
				
			# Get raw data back to physical space (nchan,nt)
			if not self._prepca[iset]: # No pre-PCA performed
				mssa_fmt_rec[iset] = self._unstack_(iset,raw_rec,taxis)
				
			elif raw: # Force direct result from MSSA
				mssa_fmt_rec[iset] = [cdms.createVariable(raw_rec.transpose())]
				mssa_fmt_rec[iset][0].setAxisList(0,[taxis,self._mode_axis_('mssa',iset)])
					
			else: # With pre-pca
				proj_rec, spcamodes = self._project_(self._pca_raw_eof[iset],raw_rec.transpose(),iset,nt=len(taxis))
				mssa_fmt_rec[iset] = self._unstack_(iset,proj_rec,taxis)
			del  raw_rec
			
			# Set attributes
			for idata,rec in enumerate(mssa_fmt_rec[iset]):
				if not self._stack_info[iset]['ids'][idata].startswith('variable_'):
					rec.id = self._stack_info[iset]['ids'][idata]+'_mssa_rec'
				else:
					rec.id = 'mssa_rec'
#				if modes is not None:
				rec.id += smodes #FIXME: do we keep it?
				rec.modes = smodes
#				else:
#					rec.modes = '1-%i'%self._nmssa[iset]
				rec.standard_name = 'recontruction_of_mssa_modes'
				rec.long_name = 'Reconstruction of MSSA modes'
				atts = self._stack_info[iset]['atts'][idata]
				if atts.has_key('long_name'):
					rec.long_name += ' of '+atts['long_name']
					
		return self._return_(mssa_fmt_rec,grouped=raw)	

	#################################################################
	## SVD
	#################################################################
	@_filldocs_
	def svd(self,usecorr=True, **kwargs):
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

		if self._ndataset<2:
			raise SpanlibError('svd','Error you need at least (most) 2 datasets to run svd, otherwise use pca and mssa')
		
		# Parameters
		self._update_('svd', **kwargs)
		print 'xxx running svd'
		# Loop on datasets
		for iset,pdata in enumerate(self._pdata[:2]):

			# Check if old results can be used when nsvd is lower
			if getattr(self,'_svd_raw_pc').has_key(iset) and \
				getattr(self,'_svd_raw_pc')[iset].shape[-1] > self._nsvd[iset]:
				continue
			
			# Remove old results
			for att in 'raw_eof','raw_pc','raw_ev','ev_sum':
				dic = getattr(self,'_svd_'+att)
				if dic.has_key(iset): del dic[iset]

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
				weights = self._stack_info[iset]['weights']
				if iset == 0:
					left = pdata
					lweights = weights

				else:
					right = pdata
					rweights = weights

		# Compute SVD
		raw_eof_left, raw_eof_right, raw_pc_left, raw_pc_right, raw_ev, ev_sum = \
			spanlib_fort.svd(left, right, self._nsvd, lweights, rweights, usecorr)
			
		# Save results
		self._svd_raw_pc[0] = raw_pc_left
		self._svd_raw_eof[0] = raw_eof_left
		self._svd_raw_pc[1] = raw_pc_right
		self._svd_raw_eof[1] = raw_eof_right
		self._svd_raw_ev = raw_ev
		self._svd_ev_sum = ev_sum


		self._last_analysis_type = 'svd'
		gc.collect()

	@_filldocs_
	def svd_eof(self,iset=None,raw=False,**kwargs):
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
		for iset in xrange(self._ndataset): # (window*nchan,nsvd)
		
			# Operate only on selected datasets
			if isets is not None and iset not in isets: continue
		
			# EOF already available 
			if self._svd_fmt_eof.has_key(iset):
				fmt_eof[iset] = self._svd_fmt_eof[iset]
				continue
				
			# No analyses performed?
			if not self._svd_raw_eof.has_key(iset): self.svd(iset=iset)
			raw_eof = self._svd_raw_eof[iset]
			nm = self._nsvd
			nl = self._svd_raw_eof[iset].shape[0]
			raw_eof = self._svd_raw_eof[iset].reshape((nl, nm))
			
			if raw: # Do not go back to physical space
				self._svd_fmt_eof[iset] = [cdms.createVariable(raw_eof.transpose())]
				self._svd_fmt_eof[iset][0].setAxisList(
					[self._mode_axis_('svd',iset),self._svd_channel_axis_(iset)])
					
			else: # Get raw data back to physical space
				raw_eof = raw_eof.reshape((nl, nm),order='F')
				
				if not self._prepca[iset]: # No pre-PCA performed
					self._svd_fmt_eof[iset] = self._unstack_(iset,raw_eof,self._mode_axis_('svd',iset))
						
				else:
					proj_eof,smodes = self._project_(self._pca_raw_eof[iset],raw_eof.transpose(),iset, nt=nm)
					self._svd_fmt_eof[iset] = self._unstack_(iset,proj_eof,self._mode_axis_('svd',iset))
					
			# Set attributes
			for idata,eof in enumerate(self._svd_fmt_eof[iset]):
				if not raw and not self._stack_info[iset]['ids'][idata].find('variable_'):
					eof.id = self._stack_info[iset]['ids'][idata]+'_svd_eof'
				else:
					eof.id = 'svd_eof'
				eof.name = eof.id
				eof.standard_name = 'empirical_orthogonal_functions_of_svd'
				eof.long_name = 'SVD empirical orthogonal functions'
				if not raw:
					atts = self._stack_info[iset]['atts'][idata]
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
			
		gc.collect()
		return self._return_(fmt_eof)

	@_filldocs_
	def svd_pc(self,iset=None,**kwargs):
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
		for iset in xrange(self._ndataset):
			
			# PC already available 
			if self._svd_fmt_pc.has_key(iset):
				fmt_pc[iset] = self._svd_fmt_pc[iset]
				continue
				
			# No analyses performed?
#			if not self._pca_raw_eof.has_key(iset): self.pca(iset=iset)
			if not self._svd_raw_pc.has_key(iset): self.svd(iset=iset)
						
			# Format the variable
			idata = 0 # Reference is first data
			pc = cdms.createVariable(npy.asarray(self._svd_raw_pc[iset][:,:self._nsvd].transpose(),order='C'))
			pc.setAxis(0,self._mode_axis_('svd',iset))
			pc.setAxis(1,self._time_axis_(iset, idata))
			pc.id = pc.name = 'svd_pc'
			pc.standard_name = 'principal_components_of_svd'
			pc.long_name = 'SVD principal components'
			atts = self._stack_info[iset]['atts'][idata]
			if atts.has_key('long_name'): pc.long_name += ' of '+atts['long_name']
#			if atts.has_key('units'):     pc.units = atts['units']

			fmt_pc[iset] = self._svd_fmt_pc[iset] = pc
			self._check_dataset_tag_('_svd_fmt_pc', iset, svd=True)

		return self._return_(fmt_pc,grouped=True)		
			
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
		raw_ev = self._svd_raw_ev[:self._nsvd]
		if cumsum:
			raw_ev = raw_ev.cumsum()
			id += '_cumsum'
			long_name.append('cumulative')
		if relative: 
			raw_ev = 100.*raw_ev/self._svd_ev_sum
			id += '_rel'
			long_name.append('relative')
		ev = cdms.createVariable(raw_ev)
		ev.id = ev.name = id
		long_name.append('SVD eigen values')
		ev.long_name = ' '.join(long_name).title()
		ev.setAxis(0, self._mode_axis_('svd'))
		ev.standard_name = 'eigen_values_of_svd'
#		atts = self._stack_info['atts'][0]
#		if atts.has_key('long_name'):
#			ev.long_name += ' of '+atts['long_name']
		if relative:
			ev.units = '% of total variance'
		return ev

	@_filldocs_
	def svd_rec(self,iset=None,modes=None,raw=False, **kwargs):
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
		if phases in [True,None]:
			phases = self._nphase
		elif phases is not False:
			self._nphase = phases
		
		# Loop on datasets
		svd_fmt_rec = {}
		for iset in isets:

			# No analyses performed?
			if not self._svd_raw_eof.has_key(iset): self.svd(iset=iset)
			
			# Projection
			raw_rec,smodes = self._project_(self._svd_raw_eof[iset],
				self._svd_raw_pc[iset],iset,modes)
				
			# Phases composites
			if phases:
				raw_rec = get_phases(raw_rec,phases)
				taxis = raw_rec.getAxis(1)
			else:
				taxis = self._time_axis_(iset)
				
			# Get raw data back to physical space (nchan,nt)
			if not self._prepca[iset]: # No pre-PCA performed
				svd_fmt_rec[iset] = self._unstack_(iset,raw_rec,taxis)
				
			elif raw: # Force direct result from svd
				svd_fmt_rec[iset] = [cdms.createVariable(raw_rec.transpose())]
				svd_fmt_rec[iset][0].setAxisList(0,[taxis,self._mode_axis_('svd',iset)])
					
			else: # With pre-pca
				proj_rec, spcamodes = self._project_(self._pca_raw_eof[iset],raw_rec.transpose(),iset,nt=len(taxis))
				svd_fmt_rec[iset] = self._unstack_(iset,proj_rec,taxis)
			del  raw_rec
			
			# Set attributes
			for idata,rec in enumerate(svd_fmt_rec[iset]):
				if not self._stack_info[iset]['ids'][idata].startswith('variable_'):
					rec.id = self._stack_info[iset]['ids'][idata]+'_svd_rec'
				else:
					rec.id = 'svd_rec'
#				if modes is not None:
				rec.id += smodes #FIXME: do we keep it?
				rec.modes = smodes
#				else:
#					rec.modes = '1-%i'%self._nsvd[iset]
				rec.standard_name = 'recontruction_of_svd_modes'
				rec.long_name = 'Reconstruction of SVD modes'
				atts = self._stack_info[iset]['atts'][idata]
				if atts.has_key('long_name'):
					rec.long_name += ' of '+atts['long_name']
					
		return self._return_(svd_fmt_rec,grouped=raw)	
	


	def _old_svd_(self,nsvd=None,pca=None):

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



	def _project_(self,reof,rpc,iset=0,imodes=None,ns=None,nt=None,nw=None):
		"""Generic raw construction of modes for pure PCA, MSSA or SVD, according to EOFs and PCs, for ONE DATASET
		
		reof: (nspace,nmode)
		rpc: (nt,nmode)
		"""

		# Get EOFs and PCs for one dataset
		if isinstance(reof, (list, tuple, dict)): reof = reof[iset]
		if isinstance(rpc, (list, tuple, dict)): rpc = rpc[iset]
		if ns is None: 
			ns = reof.shape[0]
			if nw: ns /= nw
		if nt is None: nt = self._nt[iset]

		# Which modes
		nmode = reof.shape[-1]
		imodes = _get_imodes_(imodes, nmode)

		# Function of reconstruction
		if nw is not None:
			function = spanlib_fort.mssa_rec # MSSA
		else:
			function = spanlib_fort.pca_rec  # PCA/SVD


		# Arguments
		args = [reof,rpc]
		kwargs = {}
		if nw is not None:
			# MSSA
			args.extend([ns,nt,nw])
		# Loop modes
		smodes = []
		ffrec = 0.
		for j,ims in enumerate(imodes):
			if ims[0] > nmode: break
			if ims[1] > nmode: ims[1] = nmode
			targs = args+list(ims)
			ffrec += function(*targs,**kwargs) # (nc,nt)
			if ims[0] == ims[1]:
				smode = str(ims[0])
			else:
				smode = '%i-%i'%tuple(ims)
			smodes.append(smode)
		return ffrec,'+'.join(smodes)



	def _old_reconstruct_(self,imode=None,mssa=None,pca=None,phases=False,nphases=8,offset=.5,firstphase=0,svd=None,ipair=None):
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
			raise 'Error you did not do any PCA or MSSA!\n To do a phases analysis only use the function %s in this module.\n%s' % ('get_phases',get_phases.__doc__)

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
					ffrec[i] = get_phases(ffrec[i],nphases,offset,firstphase)
			else:
				ffrec=[]
				for i in xrange(len(self._pdata)):
					ffrec.append(get_phases(npy.transpose(self._pca_raw_pc[i]),
						nphases,offset,firstphase))

			# Replace time axis with phases axis
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
		"""Return dataset as input dataset (depth and shapes)"""
		
		# Grouped case
		imap = self._input_map
		if grouped:
			if isinstance(imap, list):
				imap = [0]*len(imap)
#				imap = [m[0] for m in imap]
			else:
				imap = 0
		
		# Single variable or a single grouped list
		dataset = copy.copy(dataset)
		if imap == 0 :
			while isinstance(dataset, (list, tuple, dict)):
				dataset = dataset[0]
			return dataset
			
		# A single list of stacked variables (not serial mode)
		if isinstance(imap, int):
			return dataset[0]
		
		# Full case
		# - check groups
		for iset in dataset.keys():
			map = imap[iset]
			if imap[iset] == 0 and isinstance(dataset[iset], (list, tuple)):
				dataset[iset] = dataset[iset][0]
		# - single element required
		if len(dataset) != self._ndataset:
			return dataset.values()[0]
		# - convert to tuple
		ret = ()
		for iset in dataset.keys():
			ret += (dataset[iset], )
		gc.collect()
		return ret

	def nmssa(self, iset):
		"""
		Number of MSSA modes
		
		:Returns:
			integer or tuple
		"""
		return self._return_(self._nmssa[iset],grouped=True)
		
	def npca(self, iset):
		"""
		Number of PCA modes
		
		:Returns:
			integer or tuple
		""" 
		return self._return_(self._npca[iset],grouped=True)

	def ns(self, iset):
		"""
		Length of channel axis (unmasked input points)
		
		:Returns:
			integer or tuple
		"""
		return self._return_(self._ns[iset],grouped=True)
		
	def nt(self):
		"""
		Length of time axis
		
		:Returns:
			integer or tuple
		"""
		return self._return_(self._nmssa,grouped=True)
		
	def window(self):
		"""
		MSSA window parameter
		
		:Returns:
			integer or tuple
		"""
		return self._return_(self._nmssa,grouped=True)

	def nsvd(self):
		"""
		Number of SVD modes
		
		:Returns:
			integer or tuple
		"""
		return self._nsvd

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
		self._nphase = 8
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
				if npy.ma.isMA(data):
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
			norms.append(packed['norm'])
			masks.append(packed['mask'])
			del packed
			gc.collect()

		# Store data and information
		self._stack_info.append(dict(ids=ids,taxes=taxes,saxes=saxes,masks=masks,
			weights=weights,atts=atts,shapes=shapes,grids=grids,
			mvs=mvs,types=types,norms=norms))
		self._ndata.append(len(ids))
		self._pdata.append(stacked)
		self._nt.append(stacked.shape[1])
		if len(stacked.shape) == 2:
			self._ns.append(stacked.shape[0])
		else:
			self._ns.append(1)


	def _unstack_(self,iset,pdata,firstaxes):
		"""Return a list of variables in the physical space, for ONE dataset.
		
		pdata: (~ns,nt)
		firstaxes: MUST be a tuple
		"""

		# Loop on stacked data
		istart = 0
		unstacked = []
		#if not isinstance(firstaxes,list):
			#firstaxes = [firstaxes]
		if not isinstance(firstaxes,tuple):
			firstaxes = firstaxes,
		for idata in xrange(len(self._stack_info[iset]['ids'])):

			# Get needed stuff
			for vname in self._stack_info[iset].keys():
				if vname[-1] == 's':
					exec "%s = self._stack_info[iset]['%s'][idata]" % (vname[:-1],vname)

			# Unpack data: input is sub_space:other, output is other:split_space
			# For MSSA: pdata(nchan,window,nmssa) (other = window:nmssa)
			mlen = int(mask.sum())
			iend = istart + mlen
			indata = pdata[istart:iend,:]
			if hasattr(pdata,'filled'):
				indata = idata.filled(mv)
			unpacked = spanlib_fort.chan_unpack(mask.flat,indata,mv)
			unpacked = npy.asarray(unpacked,order='C')

			# Check axes and shape
			axes = []
			for fa in firstaxes:
				if isinstance(fa,(list,tuple)): # Time
					axes.append(fa[idata])
				elif isinstance(fa,dict): # MODE
					axes.append(fa[iset])
				else:
					axes.append(fa) # WINDOW OR OTHERS
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

			# Unnormalize
			if norm != 1.:
				masked[:] *= norm

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
				raise SpanlibException('rec','analysis_type must be one of '+str(valid))
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

