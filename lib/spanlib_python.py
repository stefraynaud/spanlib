#################################################################################
# File: spanlib_python.py
#
# This file is part of the SpanLib library.
# Copyright (C) 2006  Charles Doutiraux, Stephane Raynaud
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

import spanlib_fort
import MV
import Numeric


def stackData(*data):
    """ Takes several data files, of same time and stacks them up together

    Description:::
      This fonction concatenates several dataset that have the
      same time axis. It is useful for analysing for example
      several variables at the same time.
      It takes into account weights, masks and axes.
    :::

    Usage:::
    dout, weights, mask, axes = stackData(data1[, data2...])

      *data   :: One or more data objects to stack.
                 They must all have the same time length.
    :::

    Output:::
      dout    :: Stacked data
      weights :: Associated stacked weights
      masks   :: Associated stacked masks
      axes    :: Associated stacked axes
    :::
    """
    len_time=None
    axes=[]
    dout=None # data output
    for d in data:
        t=d.getTime()
        if t is None:
            raise 'Error, all data muist have a time dimension'
        if len_time is None:
            len_time=len(t)
        elif len_time!=len(t):
            raise 'Error all datasets must have the same time length!!!!'

        if d.getAxis(0)!=t:
            d=d(order='t...')

        axes.append(d.getAxisList())
        tdata,w,m=pack(d)
        if dout is None:
            dout=tdata
            weights=w
            masks=[m]
        else:
##             print dout.shape,tdata.shape
            dout=Numeric.concatenate((dout,tdata))
            weights=Numeric.concatenate((weights,w))
            masks.append(m)

    return Numeric.transpose(dout),weights,masks,axes

def unStackData(din,weights,masks,axes):
    """ Unstack data in the form returned from stackData

    Description:::
      This function is the reverse operation of stakData.
      It splits stacked datasets into a list.
    :::

    Usage:::
    dout = unStackData(din,weights,mask,axes)

      din     :: Stacked data (see stackData function)
      weights :: Associated stacked weights
      masks   :: Associated stacked masks
      axes    :: Associated stacked axes
    :::

    Output:::
      dout    :: List of unstacked data
    :::
    """
    nvar=len(axes)

    if nvar!=len(masks):
        raise 'Error masks and input data length not compatible'

    totsize=0
    for m in masks:
        totsize+=int(MV.sum(MV.ravel(m)))
    if totsize!=din.shape[1]:
        raise 'Error data and masks are not compatible in length!!!! (%s) and (%s)' % (totsize,din.shape[1])

    istart=0
    dout=[]
    for i in range(nvar):
        m=masks[i]
        mlen=int(MV.sum(MV.ravel(m)))
        iend=istart+mlen
        data=Numeric.transpose(din[:,istart:iend])
        w=weights[istart:iend]
        ns1=len(axes[i][-1])
        ns2=len(axes[i][-2])
        up=spanlib_fort.unpack3d(m,ns1,ns2,data.shape[1],data,mlen,1.e20)
        unpacked = MV.transpose(MV.array(up))
        unpacked.setAxisList(axes[i])
        istart+=mlen
        dout.append(unpacked)
    return dout


def pack(data,weights=None):
    """ Pack a dataset and its weights according to its mask

    Description:::
      This function packs a dataset in 2D space by removing
      all masked points and returning a space-time array.
      It performs this operation also on the weights.
      It is used for removing unnecessary points and
      simplifying the input format for analysis functions.
    :::

    Usage:::
    packed_data, packed_weights, mask = pack(data,weights)

      data    :: Flatten in space an [x,y,t] array by removing
                 its masked point
      weights :: Weights to be flatten also
    :::

    Output:::
      packed_data    :: Space-time packed array
      packed_weights :: Packed weights that were guessed or used
      mask           :: Mask that were guessed or used
    :::
    """
    if Numeric.rank(data)==2: # Already packed but then needs weights!
        if weights is None:
            raise 'Error packed data must be sent with weights!'
        else:
            packed_data = Numeric.transpose(data)
            packed_weights = weights
            mask = True
            return packed_data,packed_weights,mask
    sh=list(data.shape)
    ns1=sh[-1]
    ns2=sh[-2]
    nt=sh[0]
    ## Weights ?
    if weights is None:
        try:
            import cdutil
            tmp=data
            while tmp.rank()>2:
                tmp=tmp[0]
            weights=cdutil.area_weights(tmp).raw_data()
            del(tmp)
        except:
            weights=MV.ones((sh[-2],sh[-1]),typecode='f')

    ## Now masking part
    ## First from data
    mask = data.mask()
    if mask is None:
        mask=MV.zeros(sh,typecode='f')

    if data.rank()>2:
        mask=MV.sum(mask,axis=0)

    ## Now add the ones from the weights
    mask=mask.filled()+MV.equal(weights,0.).filled(1)

    ## >=1 means masked, Fortran "mask": 1 means data ==> 1-mask
    mask=1.-MV.greater_equal(mask,1).filled()
    mask=Numeric.transpose(mask)

    ## Number of points in spatial dimension
    ns=int(MV.sum(Numeric.ravel(mask)))

    ## Ok now calls fortran, but need to transpose first
    packed_data=Numeric.transpose(data.filled(1.e20))

    ## Dummy 1D for time for tmask
    ## Pack data
    packed_data = spanlib_fort.pack3d(packed_data,mask,ns1,ns2,nt,ns)

    weights=MV.reshape(weights,(1,sh[-2],sh[-1]))
    ## Pack weights
    tweights=Numeric.transpose(weights.filled(0))
    tweights=Numeric.ones(tweights.shape,'f')
    packed_weights = spanlib_fort.pack3d(tweights,mask,ns1,ns2,1,ns)[:,0].astype('f')

    return packed_data,packed_weights,mask


def computePhases(data,nphases=8,offset=.5,firstphase=0):
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

      data       :: Space-time data oscillatory in time data.shape is rank 2 and dim 0 is space
      nphases    :: Number of phases (divisions of the cycle)
      offset     :: Normalised offset to keep higher values only [default:
      firstphase :: Position of the first phase in the 360 degree cycle
    :::

    Output:::
      phases :: Space-phase array
    :::
    """

    ns=data.shape[0]
    nt=data.shape[1]
    w = MV.ones((ns),typecode='f')
    phases = MV.array(spanlib_fort.phasecomp(data, ns, nt, nphases, w, offset, firstphase))
    axes = MV.array(data).getAxisList()
    phases.id = 'phases'
    ax = phases.getAxis(1)
    ax[:]=ax[:]*360./nphases+firstphase
    ax.id = 'phases'
##     print axes,ax,phases.shape
    axes[1] = ax
    phases.setAxisList(axes)
    return phases

class SpAn(object):
    def __init__(self,data,weights=None,npca=None,window=None,nmssa=None):
        """ Prepare the Spectral Analysis Object

        Description:::
          This function creates an object for future analyses.
          It optionally initializes some parameters.
        :::

        Usage:::
        analysis_object = SpAn(data,weights=None,npca=None,window=None,nmssa=None)

          data    :: Data on which to run the PC Analysis
                     Last dimensions must represent the spatial dimensions.
                     Analysis will be run on the first dimension.
          weights :: If you which to apply weights on some points.
                     Set weights to "0" where you wish to mask.
                     The input data mask will be applied,
                     using the union of all none spacial dimension mask.
                     If the data are on a regular grid, area weights
                     will be generated, if the cdutil (CDAT) module is available.
                     [default: 1. everywhere]
          npca    :: Number of principal components to return [default: 10]
          nmssa   :: Number of MSSA modes retained [default: 4]
          window  :: MSSA window parameter [default: time_length/3.]
        :::

        Output:::
          analysis_object :: Object created for further analysis
        :::
        """
        ## Sets all values to None
        self.clean()
        ## First pack our data, prepare the weights and mask for PCA
        self.pdata,self.weights,self.mask = pack(data,weights)


        ## Store axes for later
        self.axes=data.getAxisList()
        self.varname=data.id

        ## Figures out dimenssions
        sh=list(data.shape)

        self.ns1=sh[-1]
        self.ns2=sh[-2]
        self.nt=sh[0]

        if npca is None:
            self.npca=10

        self.ns = self.pdata.shape[0]

        if nmssa is None:
            self.nmssa = 4

        if window is None:
            self.window = int(self.nt/3.)

##         print 'At the end:',self.pdata.shape,self.ns,self.nt

    def pca(self,npca=None):
        """ Principal Components Analysis (PCA)

        Descriptions:::
          This function performs a PCA on the analysis objects
          and returns EOF, PC and eigen values.
          EOF are automatically unpacked.
        :::

        Usage:::
        eof, pc, ev = pca(data,npca=None,weights=None)

          npca    :: Number of principal components to return, default will be 10
        :::

        Output:::
          eof :: EOF array
          pc  :: Principal Components array
          ev  :: Eigein Values array
        :::
        """

        if npca is None:
            npca=self.npca

        ## Calls Fortran pca
        self.eof,self.pc,self.ev = spanlib_fort.pca(self.pdata,self.ns,self.nt,self.npca,self.weights,1)

        Axes=self.axes[1:]

##         print 'SEF.mask is:',self.mask
        if self.mask is not None:
            eof = MV.transpose(MV.array(spanlib_fort.unpack3d(self.mask,self.ns1,self.ns2,npca,self.eof,self.ns,1.e20)))
            eof.id='EOF'
            eof.standard_name='Empirical Orthogonal Function'
        else:
            eof = MV.transpose(self.eof)
            pc  = MV.array(self.pc)

        ax=eof.getAxis(0)
        ax.id='pc'
        ax.standard_name='Principal Components Axis'
        Axes.insert(0,ax)
        eof.setAxisList(Axes)

        pc=MV.transpose(MV.array(self.pc,axes=[self.axes[0],ax]))
        pc.id='PC'
        pc.standard_name='Principal Components'

        ev=MV.array(self.ev,id='EV',axes=[ax])
        ev.standard_name='Eigen Values'


        return eof,pc,ev


    def mssa(self,nmssa=None,pca=None,window=None):
        """ MultiChannel Singular Spectrum Analysis (MSSA)

        Descriptions:::
          This function performs a MSSA on the analysis objects
          and returns EOF, PC and eigen values.
          Unless pca parameter is set to false, a pre
          PCA is performed to reduced the number of d-o-f
          if already done and if the number of channels is
          greater than 30.
        :::

        Usage:::
        eof, pc, ev = mssa(nmssa,pca)

          nmssa  :: Number of MSSA modes retained
          window :: MSSA window parameter
          pca    :: If True, performs a preliminary PCA

        Output:::
          eof :: EOF array
          pc  :: Principal Components array
          ev  :: Eigen Values  array
        :::
        """

        ## Check for default values for mssa and pca if not passed by user
        if pca is None:
            if self.pc is None and self.ns > 30: # Pre-PCA needed
                print '[mssa] The number of valid points is greater than',30,' so we perform a pre-PCA'
                pca = True
            elif self.pc is not None:
                pca = True
            else:
                pca = False

        if pca is True: # From PCA to MSSA
            nspace = self.npca
            if self.pc is None: # Still no PCA done
                self.pca()
        else:
            nspace = self.ns

        if nmssa is not None:
            self.nmssa = nmssa

        if window is not None:
            self.window = window

        if self.stpc is None: # Still no MSSA
            if pca is True: # Pre-PCA case
                self.steof, self.stpc, self.stev = spanlib_fort.mssa(Numeric.transpose(self.pc), nspace, self.nt, self.window, self.nmssa)
            else: # Direct MSSA case
                self.steof, self.stpc, self.stev = spanlib_fort.mssa(self.pdata, nspace, self.nt, self.window, self.nmssa)

        eof = MV.transpose(MV.reshape(self.steof,(self.window,nspace,self.nmssa)))
        eof.id='EOF'
        eof.standard_name='Empirical Orthogonal Function'

        ax0=eof.getAxis(0)
        ax0.id='mssa'
        ax0.standard_name='MSSA Axis'

        ax1=eof.getAxis(1)
        ax1.id='pc'
        ax1.standard_name='Principal Components Axis'

        ax2=eof.getAxis(2)
        ax2.id='window'
        ax2.standard_name='Window Axis'

        pc=MV.transpose(MV.array(self.stpc))
        pc.id='PC'
        pc.standard_name='Principal Components'
        pc.setAxis(0,ax0)
        ax3 = pc.getAxis(1)
        ax3.id='time'

        ev=MV.array(self.stev,id='EV',axes=[ax0])
        ev.standard_name='Eigen Values'
        return eof,pc,ev


    def reconstruct(self,start=1,end=None,mssa=None,pca=None,phases=False,nphases=8,offset=.5,firstphase=0):
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
        ffrec = reconstruct(start,end,mssa,pca)

          start :: First mode
          end   :: Last mode
          mssa  :: Reconstruct MSSA if True
          pca   :: Reconstruct PCA if True
          phases :: Operate phase reconstruction True/False (default is False)
        :::

        Output:::
          ffrec :: Reconstructed field
        :::
        """
        n1=start
        n2=end
        # FIXME: needs MSSA alone
        axes=self.axes
        ntimes=self.nt
        comments = 'Reconstructed from'


        ## Check for default values for mssa and pca if not passed by user
        if mssa is None:
            if self.steof is None:
                mssa = False
            else:
                mssa = True

        if pca is None:
            if self.pc is None:
                pca = False
            else:
                pca = True


##         print 'mssa,pca,phases:',mssa,pca,phases
        if phases and not pca and not mssa:
            raise 'Error you did not do any PCA or MSSA!\n To do a phases analysis only use the function %s in this module.\n%s' % ('computePhases',computePhases.__doc__)
        ## Now do the actual reconstruct job

        if mssa:
            if pca:
                nspace=self.npca
            else:
                nspace=self.pdata.shape[0]

            if n2 is None:
                n2=self.nmssa

            ffrec = spanlib_fort.mssarec(self.steof, self.stpc, nspace, self.nt, self.nmssa, self.window, n1, n2)
##             print 'Ok did mssa',ffrec.shape
            comments+=' MSSA '

        if phases:
            if mssa:
##                 print 'phase+mssa reconst'
                ffrec = computePhases(ffrec,nphases,offset,firstphase)
            else:
                ffrec = computePhases(Numeric.transpose(self.pc),nphases,offset,firstphase)

            ## Replace time axis with phases axis
            ntimes=nphases
            for i in range(len(axes)):
                if axes[i].isTime():
                    axes[i]=ffrec.getAxis(1)
            ## Attributes
                    comments+=' Phases'


        if pca:
            comments+=' PCA'
            axes=self.axes
            if mssa is True or phases is True:
                pcreconstruct = Numeric.transpose(ffrec) ; del(ffrec)
            else:
                pcreconstruct = self.pc

            if mssa:
               n1 = 1
               n2 = self.npca
            elif n2 is None:
               n2 = self.npca


##             print pcreconstruct.shape,self.ns,ntimes
            ffrec = spanlib_fort.pcarec(self.eof, pcreconstruct, self.ns, ntimes, self.npca, n1, n2)

##         print 'SEF.mask is:',self.mask
        if self.mask is not None:
            ffrec = MV.transpose(spanlib_fort.unpack3d(self.mask,self.ns1,self.ns2,ntimes,ffrec,self.ns,1.e20))
            ffrec.setAxisList(axes)
        else:
            ffrec = MV.transpose(ffrec)
            ffrec.id=self.varname
            ffrec.comment=comments

        return ffrec

    def clean(self):
        self.pc=None
        self.stpc=None
        self.steof=None
        self.eof=None
        self.stev=None
        self.ev=None







