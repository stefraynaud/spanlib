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
        if dout is None: # Create
            dout=tdata
            weights=w
            masks=[m]
        else: # Append
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
        up=spanlib_fort.unpack3d(m,data,1.e20)
        unpacked = MV.transpose(MV.array(up))
        unpacked = MV.masked_where(Numeric.equal(Numeric.resize(Numeric.transpose(m),unpacked.shape),0),unpacked)
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
            mask = None
            return packed_data,packed_weights,mask
    
    sh=list(data.shape)
##     ns1=sh[-1]
##     ns2=sh[-2]
##     nt=sh[0]
    ## Weights ?
    if weights is None:
        try:
            import cdutil
            tmp=data
            while tmp.rank()>2:
                tmp=tmp[0]
            weights=cdutil.area_weights(tmp).raw_data()
            del(tmp)
        except Exception,err:
            weights=MV.ones((sh[-2],sh[-1]),typecode='f')

    ## Now masking part
    ## First from data
    mask = data.mask()
    if mask is None:
        mask=MV.zeros(sh,typecode='f')
    else:
        mask=mask.astype('i')
        
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
    packed_data = spanlib_fort.pack3d(packed_data,mask,ns)
    
    weights=MV.reshape(weights,(1,sh[-2],sh[-1]))
    ## Pack weights
    tweights=Numeric.transpose(weights.filled(0))
    tweights=Numeric.ones(tweights.shape,'f')
    packed_weights = spanlib_fort.pack3d(tweights,mask,ns)[:,0].astype('f')
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
    phases = MV.array(spanlib_fort.phasecomp(data, nphases, w, offset, firstphase))
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
    def __init__(self,data,weights=None,npca=None,window=None,nmssa=None,nsvd=None,relative=False):
        """ Prepare the Spectral Analysis Object

        Description:::
          This function creates an object for future analyses.
          It optionally initializes some parameters.
        :::

        Usage:::
        analysis_object = SpAn(data,weights=None,npca=None,window=None,nmssa=None)

          data    :: List of data on which to run the PC Analysis
                     Last dimensions must represent the spatial dimensions.
                     Analysis will be run on the first dimension.
          weights :: If you which to apply weights on some points,
                     set weights to "0" where you wish to mask.
                     The input data mask will be applied,
                     using the union of all none spacial dimension mask.
                     If the data are on a regular grid, area weights
                     will be generated, if the cdutil (CDAT) module is available.
                     [default: 1. everywhere]
          npca    :: Number of principal components to return [default: 10]
          nmssa   :: Number of MSSA modes retained [default: 4]
          nsvd    :: Number of SVD modes retained [default: 10]
          window  :: MSSA window parameter [default: time_length/3.]
        :::

        Output:::
          analysis_object :: SpAn object created for further analysis
        :::
        """
        ## Sets all values to None
        self.clean()

        ## Before all makes sure data is list of data
        if not isinstance(data,(list,tuple)):
            data=[data,]
        if weights is None:
            weights=[None,] * len(data)
        elif not isinstance(weights,(list,tuple)):
            weights = [weights,]
            
        ## First pack our data, prepare the weights and mask for PCA
        self.pdata=[]
        self.weights=[]
        self.mask=[]
        for i in range(len(data)):
            d=data[i]
            w=weights[i]
            tmp = pack(d,w)
            tmpdata,tmpweights,tmpmask = tmp
            self.pdata.append(tmpdata)
            self.weights.append(tmpweights)
            self.mask.append(tmpmask)

        ## Store axes for later
        self.axes=[]
        self.varname=[]
        self.grids=[]
        for d in data:
            self.axes.append(d.getAxisList())
            self.varname.append(d.id)
            self.grids.append(d.getGrid())

        ## Figures out length of time dimension
        self.nt = data[0].shape[0]
        
        for d in data:
            if d.shape[0] != self.nt:
                raise Exception, 'Error your dataset are not all consistent in time length'

        if npca is None:
            self.npca=10
        else:
            self.npca=npca

        self.ns=[]
        for p in self.pdata:
            self.ns.append(p.shape[0])
            
        if nmssa is None:
            self.nmssa = 4
        else:
            self.nmssa = nmssa
            
        if nsvd is None:
            self.nsvd = 10
        else:
            self.nsvd = nsvd

        if window is None:
            self.window = int(self.nt/3.)
        else:
            self.window = window


    def pca(self,npca=None,get_ev_sum=False,relative=False):
        """ Principal Components Analysis (PCA)

        Descriptions:::
          This function performs a PCA on the analysis objects
          and returns EOF, PC and eigen values.
          EOF are automatically unpacked.
        :::

        Usage:::
        eof, pc, ev = pca(npca=None,weights=None,relative=False)

        OR

        eof, pc, ev, ev_sum = pca(npca=None,weights=None,get_ev_sum=True,relative=False)

          npca    :: Number of principal components to return [default: 10]
          get_ev_sum  :: Also return sum of all eigen values [default: False]
          relative :: Egein values are normalized to their sum (%) [default: False]
        :::

        Output:::
          eof    :: List EOF array (one per data input when created SpAn object)
          pc     :: List of Principal Components array
          ev     :: List of Eigein Values array
          ev_sum :: Sum of all eigen values (even thoses not returned).
                    Returned ONLY if get_ev_sum is True.
                    It can also be retreived with <SpAn_object>.ev_sum.
        :::
        """

        if npca is not None:
            self.npca = npca
            

        ## Calls Fortran pca
        self.eof=[]
        self.pc=[]
        self.ev=[]
        eof=[]
        pc=[]
        ev=[]
        self.ev_sum=[]

        for i in range(len(self.pdata)):
            pdat=self.pdata[i]
            w=self.weights[i]
            nteof,ntpc,ntev,ntev_sum = spanlib_fort.pca(pdat,self.npca,w,1)

            Axes=list(self.axes[i][1:])

            if self.mask[i] is not None:
                teof = MV.transpose(MV.array(spanlib_fort.unpack3d(self.mask[i],nteof,1.e20)))
                teof = MV.masked_where(Numeric.equal(Numeric.resize(Numeric.transpose(self.mask[i]),teof.shape),0),teof)
                teof.id='eof'
                teof.name = teof.id
                teof.standard_name='Empirical Orthogonal Functions'
                teof.long_name = teof.standard_name
            else:
                teof = MV.transpose(nteof)
                tpc  = MV.array(ntpc)

            ax=teof.getAxis(0)
            ax.id='mode'
            ax.standard_name='Modes in decreasing order'
            
            Axes.insert(0,ax)
            teof.setAxisList(Axes)
            teof.setGrid(self.grids[i])

            tpc=MV.transpose(MV.array(ntpc,axes=[self.axes[i][0],ax]))
            tpc.id='pc'
            tpc.standard_name='Principal Components'

            tev=MV.array(ntev,id='ev',axes=[ax])
            tev.standard_name='Eigen Values'
            if relative:
                tev[:] = tev[:] * 100. / ntev_sum
                tev.units = '%'

            for var in tpc,tev,ax:
            	var.name = var.id
            	var.long_name = var.standard_name


            self.pc.append(ntpc)
            self.eof.append(nteof)
            eof.append(teof)
            pc.append(tpc)
            ev.append(tev)
            self.ev_sum.append(ntev_sum)

        if len(eof)==1:
            ret =  [eof[0],pc[0],ev[0]]
            self.ev_sum = self.ev_sum[0]
        else:
            ret =  [eof,pc,ev]

        if get_ev_sum:
            ret.append(self.ev_sum)

        return ret
    

    def mssa(self,nmssa=None,pca=None,window=None,get_ev_sum=False,relative=False):
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
        eof, pc, ev = mssa(nmssa,pca,relative=False)

        OR

        eof, pc, ev, ev_sum = mssa(nmssa,pca,get_ev_sum=True,relative=False)

          nmssa  :: Number of MSSA modes retained
          window :: MSSA window parameter
          pca    :: If True, performs a preliminary PCA
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

        ## Check for default values for mssa and pca if not passed by user
        if pca is None:
            if self.pc ==[] and max(0,self.ns) > 30: # Pre-PCA needed
                print '[mssa] The number of valid points is greater than',30,' so we perform a pre-PCA'
                pca = True
            elif self.pc is not None:
                pca = True
            else:
                pca = False

        if pca is True: # From PCA to MSSA
            nspace = [self.npca,]*len(self.pdata)
            if self.pc ==[]: # Still no PCA done
                self.pca()
        else:
            nspace = self.ns

        if nmssa is not None:
            self.nmssa = nmssa

        if window is not None:
            self.window = window

        self.steof = []
        self.stpc  = []

        eof=[]
        ev=[]
        pc=[]
        self.stev_sum=[]
        
        for i in range(len(self.pdata)):
            if pca is True: # Pre-PCA case
                ntsteof, ntstpc, ntstev, ntev_sum = spanlib_fort.mssa(Numeric.transpose(self.pc[i]), self.window, self.nmssa)
            else: # Direct MSSA case
                ntsteof, ntstpc, ntstev, ntev_sum = spanlib_fort.mssa(self.pdata[i], self.window, self.nmssa)


            teof = MV.transpose(MV.reshape(ntsteof,(self.window,nspace[i],self.nmssa)))
            teof.id='eof'
            teof.standard_name='Empirical Orthogonal Functions'

            ax0=teof.getAxis(0)
            ax0.id='space'
            ax0.standard_name='Space'

            ax1=teof.getAxis(1)
            ax1.id='mode'
            ax1.standard_name='Modes in decreasing order'

            ax2=teof.getAxis(2)
            ax2.id='window'
            ax2.standard_name='Window Axis'

            tpc=MV.transpose(MV.array(ntstpc))
            tpc.id='pc'
            tpc.standard_name='Principal Components'
            tpc.setAxis(0,ax0)

            ax3 = tpc.getAxis(1)
            ax3.id='time'
            ax3.standard_name = 'Time'
            
            tev = MV.array(ntstev,id='ev',axes=[ax0])
            tev.standard_name = 'Eigen Values'
            tev.id = 'ev'
            if relative:
                tev[:] = tev[:] * 100. / ntev_sum
                tev.units = '%'
            
            for var in teof,tpc,tev,ax0,ax1,ax2,ax3:
            	var.name = var.id
            	var.long_name = var.standard_name

            self.stpc.append(ntstpc)
            self.steof.append(ntsteof)
            eof.append(teof)
            pc.append(tpc)
            ev.append(tev)
            self.stev_sum.append(ntev_sum)

        if len(eof)==1:
            ret = [eof[0],pc[0],ev[0]]
            self.stev_sum = self.stev_sum[0]
        else:
            ret = [eof,pc,ev]

        if get_ev_sum:
            ret.append(self.ev_sum)

        return ret
    

    def svd(self,nsvd=None,pca=None):
        """ Singular Value Decomposition (SVD)

        Descriptions:::
          This function performs a SVD
          ---blabla---
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

        ---blabla---
        :::
        """

        ## Check we have at least 2 variables!!
        ## At the moment we will not use any more variable
        if len(self.pdata)<2:
            raise Exception,'Error you need at least (most) 2 datasets to run svd, otherwise use pca and mssa'
        
        ## Check for default values for mssa and pca if not passed by user
        if pca is None:
            if self.pc ==[] and max(0,self.ns) > 30: # Pre-PCA needed
                print '[svd] The number of valid points is greater than',30,' so we perform a pre-PCA'
                pca = True
            elif self.pc is not None:
                pca = True
            else:
                pca = False

        if pca is True: # From PCA to MSSA
            nspace = [self.npca,]*len(self.pdata)
            if self.pc ==[]: # Still no PCA done
                self.pca()
        else:
            nspace = self.ns

        if nsvd is not None:
            self.nsvd = nsvd


        if pca is True: # Pre-PCA case
            lneof, rneof, lnpc, rnpc, nev = spanlib_fort.svd(Numeric.transpose(self.pc[0]), Numeric.transpose(self.pc[1]), self.nsvd)
        else: # Direct SVD case
            lneof, rneof, lnpc, rnpc, nev = spanlib_fort.svd(self.pdata[0], self.pdata[1], self.nsvd)

        self.svdeof = [lneof,rneof]
        self.svdpc = [lnpc,rnpc]

        eof=[]
        pc=[]

        for i in range(2):
            teof = MV.transpose(self.svdeof[i])
            teof.id='EOF'
            teof.standard_name='Empirical Orthogonal Function'

            ax0=teof.getAxis(0)
            ax0.id='svd'
            ax0.standard_name='SVD Axis'

            ax1=teof.getAxis(1)
            ax1.id='pc'
            ax1.standard_name='Principal Components Axis'


            tpc=MV.transpose(MV.array(self.svdpc[i]))
            tpc.id='PC'
            tpc.standard_name='Principal Components'
            tpc.setAxis(0,ax0)
            
            ax3 = tpc.getAxis(1)
            ax3.id='time'

            tev=MV.array(ntstev,id='EV',axes=[ax0])
            tev.standard_name='Eigen Values'

            eof.append(teof)
            pc.append(tpc)

        return eof[0],pc[0],eof[1],pc[1],ev


    def reconstruct(self,start=1,end=None,mssa=None,pca=None,phases=False,nphases=8,offset=.5,firstphase=0,svd=None):
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
        ntimes=self.nt
        comments = 'Reconstructed from'
        axes=list(self.axes)

        if mssa is True and self.steof == []: # Want MSSA and didn't run it!
            raise Exception, 'Error you did not run MSSA yet'

        if svd is True and self.svdeof == []:
            raise Exception, 'Error you did not run SVD yet'
            
        ## Check fr svd
        if svd is None:
            if self.svdeof == []:
                svd = False
            elif self.steof==[]:
                svd = True
                
        ## Check for default values for mssa and pca if not passed by user
        if mssa is None:
            if self.steof ==[]:
                mssa = False
            elif svd is False:
                mssa = True
            else:
                mssa = False

        if pca is None:
            if self.pc == []:
                pca = False
            else:
                pca = True


##         print 'mssa,pca,phases:',mssa,pca,phases
        if phases and not pca and not mssa:
            raise 'Error you did not do any PCA or MSSA!\n To do a phases analysis only use the function %s in this module.\n%s' % ('computePhases',computePhases.__doc__)
        ## Now do the actual reconstruct job

        if mssa:
            comments+=' MSSA '
            if pca:
                nspace=[self.npca,]*len(self.pdata[0])
            else:
                nspace=[]
                for i in range(len(self.pdata)):
                    nspace.append(self.pdata[i].shape[0])

            if n2 is None:
                n2=self.nmssa

            ffrec=[]
            for i in range(len(self.pdata)):
                ffrec.append(spanlib_fort.mssa_rec(self.steof[i], self.stpc[i], nspace[i], self.nt, self.window, n1, n2))
##             print 'Ok did mssa',ffrec.shape

        if svd:
            comments+=' SVD '
            if pca:
                nspace=[self.npca,self.npca]
            else:
                nspace=[]
                for i in range(2):
                    nspace.append(self.pdata[i].shape[0])
                
            if n2 is None:
                n2=self.nsvd

            ffrec=[]
            for i in range(2):
                ffrec.append(spanlib_fort.pca_rec(self.svdeof[i], self.svdpc[i], nspace[i], self.nt, n1, n2))

        if phases:
            comments+=' Phases'
            if mssa:
##                 print 'phase+mssa reconst'
                
                for i in range(len(self.pdata)):
                    ffrec[i] = computePhases(ffrec[i],nphases,offset,firstphase)
            else:
                ffrec=[]
                for i in range(len(self.pdata)):
                    ffrec.append(computePhases(Numeric.transpose(self.pc[i]),nphases,offset,firstphase))

            ## Replace time axis with phases axis
            ntimes=nphases
            for j in range(len(self.pdata)):
                for i in range(len(self.axes[0])):
                    if axes[j][i].isTime():
                        axes[j][i]=ffrec[j].getAxis(1)


        if svd:
            nloop=2
        else:
            nloop = len(self.pdata)
            
        if pca:
            comments+=' PCA'
            if mssa is True or phases is True or svd is True:
                pcreconstruct=[]
                for i in range(nloop):
                    pcreconstruct.append(Numeric.transpose(ffrec[i]))
                del(ffrec)
            else:
                pcreconstruct = self.pc

            if mssa:
               n1 = 1
               n2 = self.npca
            elif n2 is None:
               n2 = self.npca


##             print pcreconstruct.shape,self.ns,ntimes
            ffrec=[]
            for i in range(nloop):
                ffrec.append(spanlib_fort.pca_rec(self.eof[i], pcreconstruct[i], n1, n2))

##         print 'SEF.mask is:',self.mask
                
        for i in range(nloop):
            if self.mask[i] is not None:
                ffrec[i] = MV.transpose(spanlib_fort.unpack3d(self.mask[i],ffrec[i],1.e20))
                ffrec[i] = MV.masked_where(Numeric.equal(Numeric.resize(Numeric.transpose(self.mask[i]),ffrec[i].shape),0),ffrec[i])
            else:
                ffrec[i] = MV.transpose(ffrec[i])
            ffrec[i].setAxisList(axes[i])
            ffrec[i].id=self.varname[i]
            ffrec[i].comment=comments
            if not svd:
                ffrec[i].setGrid(self.grids[i])

        if len(ffrec)==1:
            return ffrec[0]
        else:
            return ffrec

    def clean(self):
        self.pc=[]
        self.eof=[]
        self.stpc=[]
        self.steof=[]
        self.svdpc=[]
        self.svdeof=[]







