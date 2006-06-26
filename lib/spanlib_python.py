import spanlib
import MV
import Numeric


def phases(data,nphases=8,offset=.5,firstphase=0):
    """ phases analysis """
    tdata,w,mask = pack(data)
    ns = tdata.shape[0]
    sh=data.shape
    ns1=sh[-1]
    ns2=sh[-2]
    nt=data.shape[0]
    phases = spanlib.phasecomp(tdata, ns, nt, nphases, w, offset, firstphase)
    phases = MV.transpose(spanlib.unpack3d(mask,ns1,ns2,nphases,phases,ns,1.e20))
    axes = data.getAxisList()
    phases.id='phases'
    ax=phases.getAxis(0)
    ax.id='phases'
    axes[0]=ax
    phases.setAxisList(axes)
    return phases
    
def pack(data,weights=None):
    """ Computes weights and mask"""
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
    mask=mask+MV.equal(weights,0.)

    ## >=1 means masked, Fortran "mask": 1 means data ==> 1-mask
    mask=1.-MV.greater_equal(mask,1).filled()
    mask=Numeric.transpose(mask)

    ## Number of points in spatial dimension
    ns=int(MV.sum(Numeric.ravel(mask)))

    ## Ok now calls fortran, but need to transpose first
    tdata=Numeric.transpose(data.filled(1.e20))

    ## Dummy 1D for time for tmask
    ## Pack data
    tdata = spanlib.pack3d(tdata,mask,ns1,ns2,nt,ns)

    weights=MV.reshape(weights,(1,sh[-2],sh[-1]))
    ## Pack weights
    tweights=Numeric.transpose(weights.filled(0))
    tweights=Numeric.ones(tweights.shape,'f')
    weights = spanlib.pack3d(tweights,mask,ns1,ns2,1,ns)[:,0]

    return tdata,weights.astype('f'),mask


class SpAn(object):
    def __init__(self,data,weights=None,npca=None,window=None,nmssa=None):
        """ Prepare the Spectral Analysis Object"""
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

    def pca(self,npca=None):
        """ Principal Components Analysis tool

        Usage:::
        eof, pc, ev = pca(data,npca=None,weights=None)

          data    :: Data on which to run the PC Analysis
                     Last dimensions must represent the spatial dimensions. Analysis will be run on the first dimension.
          weights :: If you which to apply weights on some points.
                    Set weights to "0" where you wish to mask.
                    The input data mask will be applied, using the union of all none spacial dimension mask.
                    If the data are on a regular grid, area weights will be generated, if the cdutil (CDAT) module is available
          npca    :: Number of principal components to return, default will be 10
        :::
        Output:::
          eof:: EOF array
          pc:: Principal Components Array
        :::
        """

        if npca is None:
            npca=self.npca

        ## Calls Fortran pca
        self.eof,self.pc,self.ev = spanlib.pca(self.pdata,self.ns,self.nt,self.npca,self.weights,1)

        eof = MV.transpose(MV.array(spanlib.unpack3d(self.mask,self.ns1,self.ns2,npca,self.eof,self.ns,1.e20)))
        eof.id='EOF'
        eof.standard_name='Empirical Orthogonal Function'

        ax=eof.getAxis(0)
        ax.id='pc'
        ax.standard_name='Principal Components Axis'
        Axes=self.axes[1:]
        Axes.insert(0,ax)
        eof.setAxisList(Axes)

        pc=MV.transpose(MV.array(self.pc,axes=[self.axes[0],ax]))
        pc.id='PC'
        pc.standard_name='Principal Components'

        ev=MV.array(self.ev,id='EV',axes=[ax])
        ev.standard_name='Eigen Values'


        return eof,pc,ev
    

    def mssa(self,nmssa=None,pca=False):
        """ Principal Components Analysis tool

        Output:::
          eof:: EOF array
          pc:: Principal Components Array
          ev:: Eigen Values
        :::
        """

        if pca is True: # runs the pre PCA
            if self.pc is None:
                self.pca()

        if self.steof is None:
            if pca is True:
                self.steof, self.stpc, self.stev = spanlib.mssa(Numeric.transpose(self.pc), self.npca, self.nt, self.window, self.nmssa)
            else:
                self.steof, self.stpc, self.stev = spanlib.mssa(self.pdata, self.ns, self.nt, self.window, self.nmssa)

            eof = MV.transpose(MV.reshape(self.steof,(self.window,self.npca,self.nmssa)))
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
    

    def reconstruct(self,start=1,end=None,mssa=True,pca=True):
        """ Reconstruct results from mssa or pca"""
        n1=start
        n2=end
        if mssa:
            if n2 is None:
                n2=self.nmssa
            if pca:
                ffrec = spanlib.mssarec(self.steof, self.stpc, self.npca, self.nt, self.nmssa, self.window, n1, n2)
                ffrec = spanlib.pcarec(self.eof, Numeric.transpose(ffrec), self.ns, self.nt, self.npca, 1,self.npca)
                ffrec = MV.transpose(spanlib.unpack3d(self.mask,self.ns1,self.ns2,self.nt,ffrec,self.ns,1.e20))
                ffrec.setAxisList(self.axes)
                ffrec.id=self.varname
                ffrec.comment='Reconstructed from MSSA and PCA'

        return ffrec

    def clean(self):
        self.pc=None
        self.stpc=None
        self.steof=None
        self.eof=None
        self.stev=None
        self.ev=None







