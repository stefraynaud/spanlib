import spanlib_fort
import MV
import Numeric


def stackData(*data):
    """ Takes several data files, of same time and stacks them up together """
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
            mask=[m]
        else:
            print dout.shape,tdata.shape
            dout=Numeric.concatenate((dout,tdata))
            weights=Numeric.concatenate((weights,w))
            mask.append(m)

    return Numeric.transpose(dout),weights,mask,axes

def unStackData(din,weights,mask,axes):
    """Unstack data in the form returned from stackData"""
    nvar=len(axes)

    if nvar!=len(mask):
        raise 'Error mask and var length not compatible'

    totsize=0
    for m in mask:
        totsize+=int(MV.sum(MV.ravel(m)))
    if totsize!=din.shape[1]:
        raise 'Error data and masks are not compatible in length!!!! (%s) and (%s)' % (totsize,din.shape[0])

    istart=0
    out=[]
    for i in range(nvar):
        m=mask[i]
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
        out.append(unpacked)
    return out


def pack(data,weights=None):
    """ Computes weights and mask"""
    if Numeric.rank(data)==2: # Already packed but then vneeds weights!
        if weights is None:
            raise 'Error packed data must be sent with weights!'
        else:
            return Numeric.transpose(data),weights,True
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
    print mask.shape,weights.shape
    mask=mask.filled()+MV.equal(weights,0.).filled(1)
    print 'After:',mask.shape

    ## >=1 means masked, Fortran "mask": 1 means data ==> 1-mask
    mask=1.-MV.greater_equal(mask,1).filled()
    mask=Numeric.transpose(mask)

    ## Number of points in spatial dimension
    ns=int(MV.sum(Numeric.ravel(mask)))

    ## Ok now calls fortran, but need to transpose first
    tdata=Numeric.transpose(data.filled(1.e20))

    ## Dummy 1D for time for tmask
    ## Pack data
    tdata = spanlib_fort.pack3d(tdata,mask,ns1,ns2,nt,ns)

    weights=MV.reshape(weights,(1,sh[-2],sh[-1]))
    ## Pack weights
    tweights=Numeric.transpose(weights.filled(0))
    tweights=Numeric.ones(tweights.shape,'f')
    weights = spanlib_fort.pack3d(tweights,mask,ns1,ns2,1,ns)[:,0]

    return tdata,weights.astype('f'),mask


def phases(data,nphases=8,offset=.5,firstphase=0):
    """ phases analysis """
    tdata,w,mask = pack(data)
    ns = tdata.shape[0]
    sh=data.shape
    ns1=sh[-1]
    ns2=sh[-2]
    nt=data.shape[0]
    phases = spanlib_fort.phasecomp(tdata, ns, nt, nphases, w, offset, firstphase)
    phases = MV.transpose(spanlib_fort.unpack3d(mask,ns1,ns2,nphases,phases,ns,1.e20))
    axes = data.getAxisList()
    phases.id='phases'
    ax=phases.getAxis(0)
    ax.id='phases'
    axes[0]=ax
    phases.setAxisList(axes)
    return phases

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

        print 'At the end:',self.pdata.shape,self.ns,self.nt

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
        self.eof,self.pc,self.ev = spanlib_fort.pca(self.pdata,self.ns,self.nt,self.npca,self.weights,1)

        if self.mask is not True:
            eof = MV.transpose(MV.array(spanlib_fort.unpack3d(self.mask,self.ns1,self.ns2,npca,self.eof,self.ns,1.e20)))
            eof.id='EOF'
            eof.standard_name='Empirical Orthogonal Function'
        else:
            eof=self.eof

        ax=eof.getAxis(0)
        ax.id='pc'
        ax.standard_name='Principal Components Axis'
        Axes=self.axes[1:]
        Axes.insert(0,ax)
        if self.mask is not True:
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
                self.steof, self.stpc, self.stev = spanlib_fort.mssa(Numeric.transpose(self.pc), self.npca, self.nt, self.window, self.nmssa)
            else:
                self.steof, self.stpc, self.stev = spanlib_fort.mssa(self.pdata, self.ns, self.nt, self.window, self.nmssa)

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
                ffrec = spanlib_fort.mssarec(self.steof, self.stpc, self.npca, self.nt, self.nmssa, self.window, n1, n2)
                ffrec = spanlib_fort.pcarec(self.eof, Numeric.transpose(ffrec), self.ns, self.nt, self.npca, 1,self.npca)
                ffrec = MV.transpose(spanlib_fort.unpack3d(self.mask,self.ns1,self.ns2,self.nt,ffrec,self.ns,1.e20))
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







