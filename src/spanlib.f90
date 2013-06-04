! File: spanlib.f90
!
! This file is part of the SpanLib library.
! Copyright (C) 2006-2009  Stephane Raynaud
! Contact: stephane dot raynaud at gmail dot com
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.

! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

module spanlib

    implicit none
    
    integer, parameter :: ierr_warning=0, ierr_error=1
    real(8), parameter :: default_missing_value=1d20, &
        & mvtol=epsilon(1.)
    
    
contains

    ! ############################################################
    ! ############################################################
    ! ## PCA PART ################################################
    ! ############################################################
    ! ############################################################
!    call sl_pca(var, nkeep, xeof=xeof, pc=pc, ev=ev, ev_sum=ev_sum,&
!     & useteof=useteof, mv=mv, &
!     & minecvalid=minecvalid, zerofill=zerofill, errmsg=errmsg)

subroutine sl_pca(var, nkeep, xeof, pc, ev, ev_sum, mv, useteof, &
    notpc, minecvalid, zerofill, errmsg)
    ! **Principal Component Analysis**
    !
    ! :Description:
    !
    !    Perform a decomposition of space-time field in a set of
    !    Empirical Orthogonal Functions (EOFs) and Principal components (PCs).
    !    By default, the analysis computes  "temporal" (T) or classical
    !    spatial (S) EOFs depending on if the space dimension is greater
    !    than the time dimension. This default behavior can be overridden.
    
    !
    ! :Necessary arguments:
    !
    !    - *var (ns,nt)*: Data
    !    - *nkeep*: Maximum number of modes to keep in outputs
    !
    ! :Optional arguments:
    !
    !    - *xeof (ns,nkeep)*: Space-mode array of EOFs
    !    - *pc (nt,nkeep)*: Time-mode array of PCs
    !    - *ev (nkeep)*: Mode array of eigen values (variances)
    !    - *ev_sum*: Sum of all egein values (even thoses not returned)
    !    - *useteof*: To force the use of T or S EOFs [0 = T, 1 = S, -1 = default]
    !    - *mv**: Missing value
    !
    ! :Dependencies:
    !    :func:`dgemm` (BLAS) :func:`dsyrk` (BLAS) :func:`dsyev` (LAPACK) 


    ! Declarations
    ! ============

    implicit none

    ! External
    ! --------
    real(8), intent(in)            :: var(:,:) ! Input data
    integer,  intent(in)            :: nkeep ! Number of modes retained
    real(8), intent(out), optional :: pc(size(var,2),nkeep), & ! Principal components
    &                                xeof(size(var,1),nkeep), & ! Empirical orthogonal functions
    &   ev(nkeep) ! Eigen values (variances)
    real(8), intent(in),  optional :: mv ! Missing value
    integer,  intent(in),  optional :: useteof, & ! Use Spatial or Temporal EOFs  [0 = T, 1 = S, -1 = default]?
                                    notpc
    real(8), intent(out), optional :: ev_sum ! Sum of eigen values (total variance)
    integer, intent(in), optional :: zerofill, minecvalid
    character(len=120), intent(out), optional :: errmsg ! Logging message (len=120)

    ! Internal
    ! --------
    integer               :: ns, nt, odim, cdim
    real(8), allocatable :: cov(:,:), subcov(:,:), nn(:,:)
    real(8), allocatable :: zeof(:,:), zvar(:,:), zmean(:)
    real(8), allocatable :: zev(:), work(:)
    integer, allocatable :: valid(:,:)!, cvalid(:)
    integer :: zuseteof, znkeepmax, i, la_info, lwork, im, nc, no, &
        & nsv, ntv, it, ncf, ic, io
    integer, allocatable :: iselect(:), iselect2(:)
    character(len=120) :: msg
    character(len=1) :: trflag
    real(8) :: zmv, zdmv, zevsumt, zevsums, w0
    logical :: zusetpc

    ! Setups
    ! ======

    ! Sizes
    ! -----
    ns = size(var,1)
    nt = size(var,2)
    znkeepmax = 100
    if(nkeep>znkeepmax)then
        if(present(errmsg))then
            write(msg,*)'You want to keep a number of PCs '//&
                & 'greater than ',znkeepmax
            errmsg = sl_errmsg(ierr_error, 'pca', msg)
        endif
        return
    end if

    ! What does the user want?
    ! ------------------------
    if(.not.present(xeof).and..not.present(pc)&
      &.and..not.present(ev))then
        if(present(errmsg))&
            & errmsg = sl_errmsg(ierr_warning, 'pca', 'Nothing to do. Quit.')
        return
    end if


    ! Valid data in 2D space
    ! ----------------------
    
    ! Missing value
    if(present(mv))then
        zmv = mv
    else
        zmv = default_missing_value
    endif
    if(zmv==0d0)then
        zdmv = 1d0
    else
        zdmv = 1d0/zmv
    endif
    
    ! Valid in the 2D space
    allocate(valid(ns, nt))
    if(present(zerofill).and.zerofill==1)then
        valid = 1
    else
        valid = merge(0, 1, abs((var-zmv)*zdmv)<=mvtol)
    endif

    ! By default, T-EOF decompostion if ns > nt
    ! -----------------------------------------
    zuseteof = -1
    if(present(useteof))zuseteof = useteof
    ntv = count(any(valid==1, dim=1))
    nsv = count(any(valid==1, dim=2))
    if(zuseteof<0)then     
        if(nsv>ntv)then
            zuseteof=1
        else
            zuseteof=0
        endif
    endif
    zusetpc = present(pc) .and. zuseteof==1 .and. (.not. present(notpc) .or. notpc==0)
    znkeepmax=100
    if(zuseteof==1)then
        nc = ntv
        ncf = nt
        no = ns
        nsv = ns
        trflag = 'T'
        odim = 1
        cdim = 2
    else ! classic
        nc = nsv
        ncf = ns
        no = nt
        ntv = nt
        trflag = 'N'
        odim = 2
        cdim = 1
    endif
  
    ! Case where data are always missing
    ! ----------------------------------
    
    ! Valid channels
    allocate(iselect(nc))
    iselect = pack((/(ic,ic=1,ncf)/), any(valid==1, dim=odim))
    if(zuseteof==1)then
        allocate(iselect2(ns))
        iselect2 = pack((/(ic,ic=1,ns)/), any(valid==1, dim=2))
    endif
    iselect = pack((/(ic,ic=1,ncf)/), any(valid==1, dim=odim))
    if(nkeep>nc)then
        if(present(errmsg))then
            if(zuseteof==1)then
                write(msg, 100)'time steps',nc
            else
                write(msg, 100)'channels',nc
            endif
            errmsg = sl_errmsg(ierr_error, 'pca', msg)
        endif
        return
    end if
100     format('You want to keep a number of PCs greater than the number of valid', &
        & A, X, I2)
        
    ! Working var array
    allocate(zvar(nsv, ntv))
    zvar = 0d0
    if(zuseteof==1)then
        zvar = var(:, iselect)
    else
        zvar = var(iselect, :) ! classic
    endif
    zvar = merge(0d0, zvar, abs((zvar-zmv)*zdmv)<=mvtol)
  
    if(zuseteof==0)deallocate(zvar)

    ! Remove the mean along T
    ! -----------------------
    allocate(zmean(nsv))
    zmean = sum(zvar, dim=2)
    if(zuseteof==1)then
        zmean = zmean / dble(sum(valid(:, iselect), dim=2))
    else ! classic
        zmean = zmean / dble(sum(valid(iselect, :), dim=2))
    endif
    do i = 1, ntv
        zvar(:, i) = zvar(:, i) - zmean
    enddo
    deallocate(zmean)

    ! EOF decomposition
    ! =================

    ! Covariances and variances
    ! -------------------------
    ! Covariances
    allocate(nn(nc,nc))
    allocate(cov(nc,nc))
    cov = 0d0
    nn = 0d0
    call dsyrk('U', trflag, nc, no, 1d0, zvar, nsv, 0d0, cov, nc)
    call dsyrk('U', trflag, nc, no, 1d0, dble(valid), nsv, 0d0, nn, nc)
    where(nn>0d0) cov = cov / nn
    deallocate(nn)
    ! Variances
    zevsums = 0d0
    do i = 1, ns
        if(any(valid(i, :)==1))then
            zevsums = zevsums + sum(zvar(i, :)**2)/ dble(sum(valid(i, :)))
        endif
    enddo
    if(zuseteof==1)then
        zevsumt = 0d0
        do i = 1, nt
            if(any(valid(:, i)==1))then
                zevsumt = zevsumt + sum(zvar(:, i)**2)/ dble(sum(valid(:, i)))
            endif
        enddo
    endif
    if(.not.present(pc))deallocate(valid)
   
    ! Diagonalization (cov: input=cov, output=eof)
    ! --------------------------------------------
    allocate(zev(nc))
    allocate(work(1))
    call dsyev('V', 'U', nc, cov, nc, zev, work, -1, la_info)
    if(la_info/=0)then
        if(present(errmsg))&
            & errmsg = sl_errmsg(ierr_error, 'pca', &
            &   la_info=la_info, la_fname='DSYEV (init phase)')
        return
    endif
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))
    call dsyev('V', 'U', nc, cov, nc, zev, work, lwork, la_info)
    if(la_info/=0)then
        if(present(errmsg))&
            & errmsg = sl_errmsg(ierr_error, 'pca', &
            &   la_info=la_info, la_fname='DSYEV')
        return
    endif
    deallocate(work) 

    ! Eigenvalues
    ! -----------
!     if(zuseteof==1) zev = zev * dble(ntv) / dble(nsv)
    zev = merge(0d0, zev, zev<0)
    if(zuseteof==1)zev = zev * zevsums / zevsumt
    if(present(ev))ev = zev(nc:nc-nkeep+1:-1)
    
    ! EOFs
    ! ----
    if(present(pc).or.present(xeof))then
        
        allocate(zeof(ns,nkeep))
        zeof = zmv
        
        if(zuseteof==1)then ! T-EOF
        
            ! PC->EOF: ZVAR*PC=EOF (USE PCA_GETEC ?)
            allocate(subcov(nc,nkeep))
            subcov = cov(:,nc:nc-nkeep+1:-1)
            deallocate(cov)
            call dgemm('N', 'N', ns, nkeep, nc, 1d0, &
                & zvar, ns, subcov, nc, 0d0, zeof, ns)
            deallocate(zvar)
            
            ! Direct PC
            if(zusetpc)then
            
                if(present(zerofill).and.zerofill>=1)then
                    pc = 0d0
                else
                    pc = zmv
                endif
                pc(iselect,:) = subcov
                do im = 1, nkeep
                    pc(iselect, im) = pc(iselect, im)*sqrt(dble(nc)*zev(nc-im+1))
                end do
                
            end if
            
            ! Effect of mask
            do i = 1, ns
                do im = 1, nkeep
                    w0 = sum(subcov(:, im)**2*dble(valid(i, iselect)))
                    if(w0>tiny(0d0))zeof(i, im) = zeof(i, im) / w0
                enddo
            enddo
            deallocate(subcov)
            
            ! Norm
            do im = 1, nkeep
                zeof(:,im) = zeof(:,im) / sqrt(sum(zeof(:,im)**2 ))
            end do
            
        else ! S-EOF
        
            zeof(iselect,:) = cov(:,nc:nc-nkeep+1:-1)
            
        end if
    
    else
        deallocate(cov)
    endif
    
    ! First channel of an EOF is >= 0
    do im = 1, nkeep 
        if(zuseteof==1 .and. zeof(iselect2(1), im)<0)then
            zeof(iselect2, im) = -zeof(iselect2, im)
            if(zusetpc) pc(iselect, im) = -pc(iselect, im)
        else if(zeof(iselect(1), im)<0)then
            zeof(iselect, im) = -zeof(iselect, im)
        endif        
    enddo
    if(zuseteof==1)deallocate(iselect2)

       
    ! Sum of all eigenvalues (useful for percentils)
    ! ----------------------------------------------
    if(present(ev_sum)) ev_sum = zevsums

    ! Free eof array
    ! --------------
    if(present(xeof))then
        xeof = zeof
        if(.not.present(pc) .or. zusetpc) deallocate(zeof)
    end if

    ! Finally get PCs
    ! ===============
    if(present(pc) .and. .not.zusetpc)then
        call sl_pca_getec(var, zeof, pc, mv=zmv, &
            & ev=zev(nc:nc-nkeep+1:-1), minvalid=minecvalid, &
            & zerofill=merge(1,0,present(zerofill).and.zerofill==2), &
            & demean=1)
    end if

end subroutine sl_pca


!############################################################
!############################################################
!############################################################


subroutine sl_pca_getec(var, xeof, ec, mv, ev, minvalid, zerofill, demean)
    ! **Compute PCA expansion coefficients**
    !
    ! :Description:
    !
    !    Get an expansion coefficients from a space-time field
    !    and a set of EOFs computed by PCA. If the input
    !    space-time field is the same as the one used to computes
    !    input ST-EOFs, these expansion coincide with the
    !    associated principal components.
    !
    ! :Necessary arguments:
    !
    !    - *var (ns, nt)*: Data, supposed to be centered by default)
    !    - *xeof (ns, nkeep)*: EOFs
    !    - *ec (nt, nmode)*: Expansion coefficients
    !

    implicit none

    ! Declarations
    ! ============

    ! External
    ! --------
    real(8), intent(in)           :: var(:,:), xeof(:,:)
    real(8), intent(out)          :: ec(size(var,2),size(xeof,2))
    real(8), intent(in), optional :: mv, ev(size(xeof,2))
    integer, intent(in), optional :: minvalid, zerofill
    integer, intent(in), optional :: demean

    ! Internal
    ! --------
    real(8), allocatable :: zvar(:,:), norm(:), zeof(:,:), zmean(:)
    real(8) :: zmv, zdmv
    integer :: ns, nt, nkeep, im, it, nc, ic
    integer,allocatable :: valid(:,:), zcount(:)
    integer :: zminvalid
    logical :: zf
    logical, allocatable :: cvalid(:)
    integer, allocatable :: iselect(:)

    ! Computations
    ! ============

    ! Sizes
    ! -----
    ns = size(var,1)
    nt = size(var,2)
    nkeep = size(xeof,2)

    ! Missing values
    ! --------------
    
    ! Missing value
    if(present(mv))then
        zmv = mv
    else
        zmv = default_missing_value
    endif
    if(zmv==0d0)then
        zdmv = 1d0
    else
        zdmv = 1d0/zmv
    endif
    
    ! Valid data and compression
    allocate(valid(ns,nt), cvalid(ns))
    cvalid = abs((xeof(ns,1)-zmv)*zdmv)>mvtol
    nc = count(cvalid)
    allocate(iselect(nc))
    iselect = pack((/(im,im=1,nc)/), cvalid) ! spatial selection
    deallocate(cvalid)
    valid = merge(1, 0, abs((var-zmv)*zdmv)>mvtol)
    allocate(zvar(nc,nt), zeof(nc,nkeep))
    zvar = merge(var(iselect,:), 0d0, valid(iselect,:)==1)
    zeof = xeof(iselect, :)
    where(abs((zeof-zmv)*zdmv)<=mvtol)zeof = 0d0
    
    ! Min number of valid values for projections
    if(present(minvalid).and.minvalid/=0)then
        zminvalid = minvalid
        
    else
        zminvalid = -50
    endif
    zminvalid = max(zminvalid, -100)
    if(zminvalid<0)zminvalid = -nc*zminvalid/100
    zminvalid = max(1, zminvalid)
    zf = present(zerofill).and.zerofill/=0


    ! Remove the mean along T
    ! -----------------------
    if(present(demean).and.demean==1)then
        allocate(zmean(nc),zcount(nc))
        zcount = sum(valid(iselect,:), dim=2)
        zmean = sum(zvar, dim=2) / dble(merge(1, zcount, zcount==0))
        do it = 1, nt
            zvar(:,it) = zvar(:,it) - zmean
        enddo
        deallocate(zmean,zcount)
    endif

    ! Compute EC
    ! ----------
    call dgemm('T', 'N', nt, nkeep, nc, 1d0, zvar, nc, zeof, nc, 0d0, ec, nt)
    deallocate(zvar)
    
    ! Normalisations
    ! --------------
    allocate(norm(nt))
    do im = 1, nkeep
    
        ! Scale each time step depending on mask
        do it=1, nt
            norm(it) =  sum(zeof(:,im)**2 * dble(merge(1,valid(iselect,it),zf)))
        enddo       
        ec(:,im) = ec(:,im) / merge(norm, 1d0, norm/=0d0)
        
        ! Norm using eigenvalues
        if(present(ev) .and. any(ec>tiny(1d0)))then
            if(ev(im)>tiny(1d0))then
                ec(:,im) = ec(:,im) * sqrt(dble(nt)*ev(im)/sum(ec(:,im)**2))
            else
                ec(:,im) = 0d0
            endif
        endif
        
    end do
    deallocate(norm,zeof)

    ! Mask insignificant value
    ! ------------------------
    do it=1,nt
        if(sum(valid(iselect,it))<zminvalid) ec(it,:)=zmv
    enddo
    deallocate(valid)

end subroutine sl_pca_getec


!############################################################
!############################################################
!############################################################


subroutine sl_pca_rec(xeof, pc, varrec, istart, iend, mv, errmsg)

    ! **Reconstruction of a set of PCA components**
    !
    ! :Description:
    !
    !    Perform a reconstruction using a set of components previously
    !    computed with a PCA. All the reconstructed components are summed.
    !    A reconstructed component is simply the "product" of an EOF
    !    by its PC. The sum of all reconstructed component is the original field.
    !
    ! :Necessary arguments:
    !
    !    - *xeof (ns,nkeep)*: Space-mode array of EOFs
    !    - *pc (nt,nkeep)*: Time-mode array of PCs
    !    - *varrec (ns,nt)*: Space-time array of the reconstructed field
    !
    ! :Optional arguments:
    !
    !    - *istart*: Index of the first component to use
    !    - *iend*: Index of the last component to use


    ! Declarations
    ! ============

    implicit none

    ! External
    ! --------
    real(8),    intent(in)     :: xeof(:,:), pc(:,:)
    real(8),    intent(out)    :: varrec(size(xeof,1),size(pc,1))
    integer,intent(in),    optional    :: istart, iend
    real(8), intent(in), optional :: mv
    character(len=120), optional :: errmsg

    ! Internal
    ! --------
    integer           :: nkept, itmp, zistart, ziend, nt, ns, i, istatus
    real(8), allocatable    :: zpc(:,:), zeof(:,:)
    character(len=120) :: msg
    real(8) :: zmv, zdmv

    ! Setup
    ! =====
    nkept = size(xeof,2)
    zistart=1
    if(present(istart))zistart=istart
    if(present(iend))then
        ziend=iend
    else
        ziend=nkept
    end if
    if(present(errmsg))errmsg = ''
    if(zistart.lt.1 .or. zistart.gt.nkept)then
        zistart=1
        if(present(errmsg)) &
            & errmsg = 'istart lower than 1 => set to 1 | '
    end if
    if(ziend.lt.1 .or. ziend.gt.nkept)then
        ziend=nkept
        if(present(errmsg))then
            write(msg,*)'iend greater than the number '//&
                &'of avalaible modes => reduced to ',ziend
            errmsg = errmsg // msg // ' | '
        endif
    end if
    if(zistart>ziend)then
        itmp=ziend
        ziend=zistart
        zistart=itmp
        if(present(errmsg)) &
            & errmsg = errmsg // 'istart > iend => inversion | '
    end if
    if(present(errmsg).and.len_trim(errmsg)>0)then
        if(len_trim(errmsg)>3) errmsg = errmsg(:len_trim(errmsg)-3)
        errmsg = sl_errmsg(ierr_warning, 'pcarec', errmsg)
    endif
    ns = size(xeof,1)
    nt = size(pc,1)
    if(present(mv))then
        zmv = mv
    else
        zmv = default_missing_value
    endif
    if(zmv==0d0)then
        zdmv = 1d0
    else
        zdmv = 1d0/zmv
    endif


    ! Computation
    ! ===========
    
    ! Product
    ! -------
    allocate(zpc(nt, ziend-zistart+1), zeof(ns, ziend-zistart+1))  
    zpc = pc(:, zistart:ziend)
    zpc = merge(0d0,zpc, abs((zpc-zmv)*zdmv)<=mvtol)
    zeof = xeof(:, zistart:ziend)
    zeof = merge(0d0, zeof, abs((zeof-zmv)*zdmv)<=mvtol)
    call dgemm('N', 'T', ns, nt, ziend-zistart+1, 1d0, &
        & xeof(:, zistart:ziend), ns, zpc, nt, 0d0, varrec, ns)
    deallocate(zpc,zeof)
    
    ! Masking
    ! -------
    do i = 1, nt
        if(abs((pc(i,1)-zmv)*zdmv)<=mvtol) varrec(:,i) = zmv
    enddo
    do i = 1, ns
        if(abs((xeof(i,1)-zmv)*zdmv)<=mvtol) varrec(i,:) = zmv
    enddo
    
end subroutine sl_pca_rec




!############################################################
!############################################################
!## MSSA PART ###############################################
!############################################################
!############################################################

subroutine sl_mssa(var, nwindow, nkeep, steof, stpc, ev, ev_sum, mv, &
    & minecvalid, zerofill, errmsg)

    ! **Multi-channel Singular Spectrum Analysis**
    !
    ! :Description:
    !
    !    Perform a decomposition of space-time field in a set of
    !    space-time Empirical Orthogonal Functions (EOFs) and
    !    time Principal components (PCs), according to a window
    !    parameter.
    !
    ! :Necessary arguments:
    !
    !    - *var*:      Space-time array
    !    - *nwindow*: Window size
    !    - *nkeep*:   Maximum number of modes to keep in outputs
    !
    ! :Optional arguments:
    !
    !    - *steof*: Space-window-mode array of EOFs
    !    - *stpc*: Time-mode array of PCs
    !    - *ev*: Mode array of eigen values (variances)
    !    - *ev_sum*: Sum of all eigen values (even thoses not returned)
    !
    ! :Dependencies:
    !    :f:func:`sl_stcov` :f:func:`dsyev` (LAPACK)


    ! Declarations
    ! ============

!         use spanlib_lapack95, only: la_syevd, la_syev

    implicit none

    ! External
    ! --------
    real(8), intent(in)            :: var(:,:)
    integer,  intent(in)            :: nwindow, nkeep
    real(8), intent(out), optional :: &
        & steof(size(var,1)*nwindow, nkeep), &
        & stpc(size(var,2)-nwindow+1, nkeep), ev(nkeep)
    real(8), intent(out), optional :: ev_sum
    real(8), intent(in), optional :: mv
    integer, intent(in), optional :: zerofill, minecvalid
    character(len=120), optional :: errmsg

    ! Internal
    ! --------
    real(8), allocatable :: cov(:,:), zev(:), &
        & zvar(:,:), zsteof(:,:), work(:), zmean(:)
    integer, allocatable :: valid(: ,:)
    integer :: nchan, nsteof, nt, znkeepmax, la_info, lwork, istatus, im, it
    real(8) :: zmv
    character(len=120) :: msg

    
    ! Setup
    ! =====

    ! Sizes
    ! -----
    nchan = size(var,1)
    nsteof = nchan * nwindow
    nt = size(var,2)
    znkeepmax = 100
    if(nkeep>znkeepmax)then
        if(present(errmsg))then
            write(msg, 100) ':', znkeepmax
            errmsg = sl_errmsg(ierr_error, 'mssa', msg)
        endif
        return
    else if(nkeep>nsteof) then
        if(present(errmsg))then
            write(msg, 100) ' the number of ST-EOFs:', nsteof
            errmsg = sl_errmsg(ierr_error, 'mssa', msg)
        endif
        return
    end if
100     format('You want to keep a number of PCs greater than', &
        & A, X, I2)

    ! Missing values
    ! --------------
    if(present(mv))then
        zmv = mv
    else
        zmv = default_missing_value
    endif
    allocate(valid(nchan, nt), zvar(nchan, nt))
    if(present(zerofill).and.zerofill==1)then
        valid = 1
    else
        valid = merge(0, 1, abs((var-zmv)/zmv)<=mvtol)
    endif
    zvar = merge(0d0, var, abs((var-zmv)/zmv)<=mvtol)

    ! Remove the mean
    ! ---------------
    allocate(zmean(nchan))
    zmean = sum(zvar, dim=2)/dble(sum(valid, dim=2))
    do it = 1, nt
        zvar(:, it) = zvar(:, it) - zmean
    enddo
    deallocate(zmean)
    if(.not.present(stpc)) deallocate(valid)

    ! Set the block-Toeplitz covariance matrix
    ! ========================================
    allocate(cov(nsteof, nsteof))
    call sl_stcov(merge(zmv, zvar, var==zmv), cov, zmv)

    ! Diagonalisation
    ! ===============
    allocate(zev(nsteof))
    allocate(work(1))
    call dsyev('V', 'U', nsteof, cov, nsteof, zev, work, -1, la_info)
    if(la_info/=0)then
        if(present(errmsg))&
            & errmsg = sl_errmsg(ierr_error, 'mssa', &
            &   la_info=la_info, la_fname='DSYEV (init phase)')
        return
    endif
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))
    call dsyev('V', 'U', nsteof, cov, nsteof, zev, work, lwork, la_info)
    if(la_info/=0)then
        if(present(errmsg))&
            & errmsg = sl_errmsg(ierr_error, 'mssa', &
            &   la_info=la_info, la_fname='DSYEV')
        return
    endif
    deallocate(work)            


    ! Get ST-EOFs and eigenvalues
    ! ===========================
    if(present(steof).or.present(stpc))then
        allocate(zsteof(nsteof, nkeep))
        zsteof = cov(:, nsteof : nsteof-nkeep+1 : -1)
        deallocate(cov)
        do im = 1, nkeep ! First point of an EOF is >= 0
            if(zsteof(1, im)<0.)zsteof(:, im) = -zsteof(:, im)
        enddo
        if(present(steof))then
            steof = zsteof
            deallocate(zsteof)
        end if
    end if

    ! Eigen values
    ! ------------
    if(present(ev))then
        ev = zev(nsteof : nsteof-nkeep+1 : -1)
    end if
    if(present(ev_sum)) ev_sum = sum(zev)
    deallocate(zev)


    ! Get ST-PCs
    ! ==========
    if(present(stpc)) call sl_mssa_getec(merge(zvar, zmv, valid==1), &
        & steof, nwindow, stpc, zmv, &
        & minvalid=minecvalid, &
        & zerofill=merge(1,0,present(zerofill).and.zerofill==2))
    deallocate(zvar, valid)

end subroutine sl_mssa




!############################################################
!############################################################
!############################################################

subroutine sl_stcov(var, cov, mv)
! Compute the Block-Toeplitz covariance matrix for MSSA analysis
!
! .. note:: ``var`` does not need to be centered


    implicit none
    
    ! Declarations
    ! ------------
    
    real(8), intent(in)  :: var(:, :)
    real(8), intent(out) :: cov(:, :)
    real(8), intent(in), optional :: mv
    
    real(8), allocatable :: zvar(:,:)
    integer, allocatable :: valid(:,:), nn(:,:)
    integer ::  nchan, nt, nsteof, nwindow
    integer :: iw, iw1, iw2, i1, i2, ic1, ic2
    real(8) :: zmv

    ! Sizes
    ! -----
    nchan = size(var, 1)
    nt = size(var, 2)
    nsteof = size(cov, 1)
    nwindow = nsteof/nchan
    
    ! Missing values
    ! --------------
    if(present(mv))then
        zmv = mv
    else
        zmv = default_missing_value
    endif
    allocate(valid(nchan,nt), zvar(nchan,nt))
    valid = merge(0, 1, abs((var-zmv)/zmv)<=mvtol)
    zvar = merge(0d0, var, valid==0)
   
    ! Anomaly
    ! -------
    zvar = zvar-spread(sum(zvar, dim=2)/sum(dble(valid), dim=2), ncopies=nt, dim=2)

    ! Covariances
    ! -----------
    allocate(nn(size(cov,1),size(cov,2)))
    
    !$OMP PARALLEL &
    !$OMP SHARED(zvar,cov,nchan,nwindow,valid,nn,nt) &
    !$OMP PRIVATE(iw, iw1, iw2, i1, i2, ic1, ic2)
    !$OMP DO 
    do ic1 = 1, nchan
        do ic2 = 1, nchan
            do iw2 = 1, nwindow
                do iw1 = 1, iw2
                    i1 = (ic1-1) * nwindow + iw1
                    i2 = (ic2-1) * nwindow + iw2
                    iw = iw2 - iw1 + 1
                    cov(i1,i2) = &
                        & dot_product(zvar(ic1, 1  : nt-iw+1),  &
                        &             zvar(ic2, iw : nt     ))
                    cov(i2,i1) = cov(i1,i2)
                    nn(i1,i2) = dot_product(valid(ic1, 1  : nt-iw+1), &
                                            valid(ic2, iw : nt     ))
                    nn(i2,i1) = nn(i1,i2)
                end do
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    deallocate(zvar,valid)
    cov = merge(cov/dble(nn), 0d0, nn>0)
    deallocate(nn)

    
!            do i2 = 1, nwindow
!                do i1 = 1, iw2
!                    iw = i2 - i1 + 1
!                    cov(i1,i2) = &
!                        & dot_product(zvar(1  : nt-iw+1),  &
!                        &             zvar(iw : nt     ))! 
!                    cov(i2,i1) = cov(i1,i2)
!                    nn(i1,i2) = dot_product(valid(1  : nt-iw+1), &
!                                            valid(iw : nt     ))
!                    nn(i2,i1) = nn(i1,i2)
!                end do
!            end do
   
end subroutine sl_stcov


!############################################################
!############################################################
!############################################################


subroutine sl_mssa_getec(var, steof, nwindow, stec, mv, minvalid, zerofill)

    ! Computes MSSA expansion coefficients
    !
    ! :Description:
    !    Get an expansion coefficients from a space-time field
    !    and a set of ST-EOFs computed by MSSA. If the input
    !    space-time field is the same as the one used to computes
    !    input ST-EOFs, these expansion coincide with the
    !    associated principal components.
    !
    ! :Necessary arguments:
    !    - var (nc, nt):   Space-time field
    !    - steof (nw*nc, nm): Space-window-mode EOFs
    !    - stec:  Time-mode array of expansion coefficients
    !
    ! :Dependencies:
    !    [sd]gemm(BLAS)
    
    implicit none
    
    ! Declarations
    ! ============
    
    ! External
    ! --------
    real(8), intent(in)  :: var(:,:), steof(:,:)
    real(8), intent(out) :: stec(size(var,2)-nwindow+1,&
        & size(steof,2))
    integer,       intent(in)  :: nwindow
    real(8), intent(in), optional :: mv
    integer, intent(in), optional :: zerofill, & ! Fill var missing values with zeros?
        & minvalid ! Minimal number of data available at one time step during projection
    
    ! Internal
    ! --------
    integer :: nt, nkeep, im, iw, nchan, ntpc, it, zminvalid
    real(8), allocatable :: wpc(:), substeof(:), subvar(:,:), norm(:,:), &
        & zvalid(:,:)
    real(8) :: zmv

    ! Computations
    ! ============
  
    ! Initialisations
    ! ---------------
    stec = 0d0
    nchan = size(var,1)
    nt = size(var,2)
    nkeep = size(steof,2)
    ntpc = size(stec, 1)
    if(present(mv))then
        zmv = mv
    else
        zmv = default_missing_value
    endif
    allocate(wpc(ntpc), substeof(nchan), subvar(nchan, ntpc), norm(ntpc, nkeep))
    allocate(zvalid(nchan, nwindow))
    norm = 0d0
    if(present(minvalid).and.minvalid/=0)then
        zminvalid = minvalid
    else
        zminvalid = -50
    endif
    zminvalid = max(zminvalid, -100)
    if(zminvalid<0)zminvalid = -nchan*nwindow*zminvalid/100
    zminvalid = max(1, zminvalid)
    
    ! Main stuff
    ! ----------
    do im = 1, nkeep
        do iw = 1, nwindow
            subvar = merge(0d0, var(:,iw:iw+ntpc-1), &
                & abs((var(:,iw:iw+ntpc-1)-zmv)/zmv)<=mvtol)
            substeof = steof(iw:iw+(nchan-1)*nwindow:nwindow, im)
            call dgemm('T', 'N', nt-nwindow+1, 1, nchan, 1d0,&
                & subvar, nchan, substeof, nchan, 0d0, wpc, ntpc)
            stec(:, im)  =  stec(:, im) + wpc                    
        end do
        do it=1, ntpc
            zvalid = merge(0d0, 1d0, abs((var(:, it:it+nwindow-1)-zmv)/zmv)<=mvtol)
            if(sum(zvalid)>=zminvalid)then
                if(present(zerofill).and.zerofill/=0)then
                    norm(it, im) = 1d0
                else
                    norm(it, im) = sum(transpose(reshape(steof(:,im)**2, &
                        &(/nwindow, nchan/))) * zvalid) 
                endif 
            endif
        end do
        !stec(:, im) = stec(:, im) / sum(steof(:,im)**2)
    end do
    deallocate(subvar, substeof, zvalid, wpc)
    stec = merge(stec/norm, zmv, norm/=0d0)
    
end subroutine sl_mssa_getec


!############################################################
!############################################################
!############################################################


subroutine sl_mssa_rec(steof, stpc, nwindow, varrec, istart, &
    & iend, mv, ev, errmsg)
    ! Reconstruction of a set of MSSA components
    !
    ! Description:
    !    Same as for the reconstruction of PCA components, but for MSSA.
    !
    ! :Necessary arguments:
    !    - steof:   SpaceXwindow-mode array of EOFs
    !    - stpc:    Time-mode array of PCs
    !    - nwindow: Window size
    !    - varrec:   Space-time array of the reconstructed field
    !
    ! :Optional arguments:
    !    - istart: Index of the first component to use
    !    - iend:   Index of the last component to use
    
    implicit none
    
    
    ! Declarations
    ! ============
    
    ! External
    ! --------
    real(8),   intent(in)  :: steof(:,:), stpc(:,:)
    real(8),   intent(out) :: varrec(size(steof, 1)/nwindow,&
     &                                    size(stpc, 1)+nwindow-1)
    integer,intent(in)           :: nwindow
    integer,intent(in), optional :: istart, iend
    real(8), intent(in), optional :: mv, ev(size(steof, 2))
    character(len=120), intent(out), optional :: errmsg ! Logging message
        
    ! Internal
    ! --------
    integer :: ntpc, nchan, nt, ic, im, iw, nkept, &
     &         itmp, zistart, ziend, istatus, nv, ddof
    character(len=200) :: msg
    real(8), allocatable :: zeof(:), epc(:,:), zpc(:,:), zev(:)
    real(8) :: zmv
    logical, allocatable :: valid(:,:), pvalid(:,:), pcvalid(:,:)
    
    
    ! Setup
    ! =====
    ddof = 0
  
    ! Sizes
    ! -----
    ntpc  = size(stpc, 1)
    nt    = ntpc+nwindow-1
    nchan = size(steof, 1)/nwindow
    nkept = size(steof, 2)
    allocate(zeof(nwindow), zpc(ntpc,nkept))
    allocate(epc(nwindow, ntpc-nwindow+1))
    varrec = 0d0

    ! Missing values
    ! --------------
    if(present(mv))then
        zmv = mv
    else
        zmv = default_missing_value
    endif
    allocate(pvalid(nwindow, ntpc-nwindow+1))
    allocate(valid(size(varrec, 1), size(varrec, 2)))
    allocate(pcvalid(ntpc,nkept))
    pcvalid = abs((stpc-zmv)/zmv)>mvtol
    zpc = merge(stpc, 0d0, pcvalid)
    
    ! Renormalization of ST-PCs using eigenvalues
    ! -------------------------------------------
    if(present(ev).and..not.all(ev==0d0))then
    stop
        allocate(zev(nkept))
!        zev = merge(0d0, ev, ev<0.)
!        do ic = 1, nkept
!            zpc(:, ic) = zpc(:, ic) * sqrt(zev(ic) * &
!                & dble(count(pcvalid(:,1))) / sum(zpc(:, ic)**2))
!            
!        end do
        deallocate(zev)
    end if
    
    ! Range
    ! -----
    if(present(iend))then
        ziend=iend
    else
        ziend=nkept
    end if
    if(present(istart))then
        zistart=istart
    else
        zistart = 1
    endif
    if(present(errmsg))errmsg = ''
    if(zistart.lt.1.or.zistart.gt.nkept)then
        zistart = 1
        if(present(errmsg))&
            & errmsg = 'istart lower than 1 => set to 1 | '
    end if
    if(ziend.lt.1.or.ziend.gt.nkept)then
        ziend = nkept
        if(present(errmsg))then
            write(msg,*)'iend greater than the number of '// &
                &     'avalaible modes => reduced to',iend
            errmsg = errmsg // msg // ' | '
        endif
    end if
    if(zistart>ziend)then
        itmp    = ziend
        ziend   = zistart
        zistart = itmp
        if(present(errmsg))&
            & errmsg = errmsg // 'istart > iend => inversion | '
    end if
    if(present(errmsg).and.len_trim(errmsg)>0)then
        if(len_trim(errmsg)>3) errmsg = errmsg(:len_trim(errmsg)-3)        
        errmsg = sl_errmsg(ierr_warning, errmsg)
    endif
    
    ! Computation
    ! ===========
    varrec = 0d0
    valid = .false.
    do im = zistart, ziend ! sum over the selection of modes
    
        ! (ntpc-nwindow+1) length slices
        do iw = 1, nwindow
            epc(iw,:) = zpc(iw : iw+ntpc-nwindow, im)
            pvalid(iw,:) = pcvalid(iw : iw+ntpc-nwindow, im)
        end do
    
        do ic = 1, nchan ! sum over the channels (= space, or PCs from PCA)
    
            ! reversed eof
            zeof = steof(nwindow+(ic-1)*nwindow : 1+(ic-1)*nwindow : -1, im)
    
            ! * middle * [nwindow length projections]
            varrec(ic, nwindow : ntpc) =  varrec(ic, nwindow : ntpc) + &
                & matmul(zeof, epc) / dble(max(count(pvalid, dim=1),1))! dble(nwindow)
            valid(ic, nwindow : ntpc) = valid(ic, nwindow : ntpc) .or. &
                & any(pvalid, dim=1)
                    
            do iw = 1, nwindow-1
!    
                ! * beginning * [iw length projections]
                varrec(ic, iw) = varrec(ic, iw) + &
                    & dot_product(zeof(nwindow-iw+1:nwindow), &
                    & zpc               (1:iw, im)) / &
                    & dble(max(count(pcvalid(1:iw, im)), 1)) !dble(iw)
                valid(ic, iw) = valid(ic, iw) .or. any(pcvalid(1:iw, im))
                
                ! * end * [iw length projections]                      
                varrec(ic, nt-iw+1) = varrec(ic, nt-iw+1) + &
                    & dot_product(zeof(1:iw), &
                    & zpc               (ntpc-iw+1:ntpc, im)) / &
                    & dble(max(count(pcvalid(ntpc-iw+1:ntpc, im)),1)) !dble(iw)
                valid(ic, nt-iw+1) = valid(ic, nt-iw+1) .or. any(pcvalid(ntpc-iw+1:ntpc, im))
                
            end do

    
        end do
    
    end do
    
    deallocate(zeof, zpc, epc, pvalid, pcvalid)
    where(.not.valid) varrec = zmv
    deallocate(valid)

end subroutine sl_mssa_rec



!############################################################
!############################################################
!## SVD PART ################################################
!############################################################
!############################################################

subroutine sl_svd(ll, rr, nkeep, leof, reof, lpc, rpc, &
    & ev, ev_sum, usecorr, mv, minecvalid, errmsg)
    ! Title:
    !    Singular Value Decomposition
    !
    ! Description:
    !    Singular value decomposition between two datasets having
    !    the same length in time.
    !
    ! :Necessary arguments:
    !    - ll:    Left space-time array
    !    - rr:    Right space-time array
    !    - nkeep: Maximum number of modes to keep in outputs
    !
    ! :Optional arguments:
    !    - leof:  Left EOFs
    !    - reof:  Right EOFs
    !    - lpc:   Left PCs
    !    - rpc:   Right PCs
    !    - ev:    Eigen values
    !    - usecorr:  Use correlations instead of covariances
    !
    ! :Dependencies:
    !    :func:`sdgemm` (BLAS) :func:`dgesvd` (LAPACK) :func:`dgesdd` (LAPACK)
    
    
    ! Declarations
    ! ============
    
    !     use spanlib_lapack95, only: la_gesdd, la_gesvd
    
    implicit none
    
    ! External
    ! --------
    real(8), intent(in)           :: ll(:,:),rr(:,:)
    integer,  intent(in)           :: nkeep
    real(8), intent(out),optional :: lpc(size(ll,2),nkeep), &
        & leof(size(ll,1),nkeep), rpc(size(rr,2),nkeep), &
        & reof(size(rr,1),nkeep),ev(nkeep)
    integer, intent(in),  optional :: usecorr
    real(8), intent(in),  optional :: mv ! Missing value
    real(8), intent(out), optional :: ev_sum
    integer, intent(in), optional :: minecvalid
    character(len=120), intent(out), optional :: errmsg ! Logging message
    
    ! Internal
    ! --------
    integer               :: ns,nsl,nsr,nt,nslv,nsrv
    real(8), allocatable :: zll(:,:), zrr(:,:), cov(:,:), &
        &                    zls(:), zrs(:)
    real(8), allocatable :: zev(:), zleof(:,:), work(:), nn(:,:), zlmean(:), zrmean(:)
    real(8)              :: zvt(1, 1)
    integer               :: znkeepmax, i, it, la_info, lwork, istatus
    integer, allocatable :: lvalid(:,:),slvalid(:),rvalid(:,:),srvalid(:)
    integer               :: zbcorr
    character(len=120) :: msg
    real(8) :: zmv, zdmv
    integer, allocatable :: ilselect(:),irselect(:)
    
    
    ! Sizes
    ! -----
    la_info = 0
    nsl = size(ll,1)
    nsr = size(rr,1)
    nt = size(ll,2)
    !     if(nsl/=nsr.or.nt/=size(rr,2))then
    if(nt/=size(rr,2))then
        if(present(errmsg))&
            & errmsg = sl_errmsg(ierr_error, 'svd', &
                & 'Left and right arrays have incompatible sizes')
        return
    end if
    
    ! What does the user want?
    ! ------------------------
    if(.not.present(leof).and..not.present(lpc).and.&
      &.not.present(reof).and..not.present(rpc).and.&
      &.not.present(ev))then
         if(present(errmsg))&
            & errmsg = sl_errmsg(ierr_warning, 'svd', &
                & 'Nothing to do => quit.')
        return
    end if
    if(present(usecorr))then
        zbcorr = usecorr
    else
        zbcorr = 0
    endif
    
    ! Valid data in 2D space
    ! ----------------------
    
    ! Missing value
    if(present(mv))then
        zmv = mv
    else
        zmv = default_missing_value
    endif
    if(zmv==0d0)then
        zdmv = 1d0
    else
        zdmv = 1d0/zmv
    endif
    
    ! Selection
    allocate(lvalid(nsl,nt),rvalid(nsr,nt))
    lvalid = merge(0, 1, abs((ll-zmv)*zdmv)<=mvtol)
    rvalid = merge(0, 1, abs((rr-zmv)*zdmv)<=mvtol)
    allocate(slvalid(nsl),srvalid(nsr))
    slvalid = sum(lvalid, dim=2)
    srvalid = sum(rvalid, dim=2)
    nslv = count(slvalid/=0)
    nsrv = count(srvalid/=0)
    allocate(zll(nslv,nt),ilselect(nslv))
    allocate(zrr(nsrv,nt),irselect(nsrv))
    ilselect = pack((/(i,i=1,nsl)/), slvalid/=0)
    irselect = pack((/(i,i=1,nsr)/), srvalid/=0)
    zll = ll(ilselect, :)
    zrr = rr(irselect, :)
    zll = merge(0d0, zll, abs((zll-zmv)*zdmv)<=mvtol)
    zrr = merge(0d0, zrr, abs((zrr-zmv)*zdmv)<=mvtol)
        
    ! Channels versus number of modes
    ns = min(nslv,nsrv)
    znkeepmax = 100
    if(nkeep>znkeepmax)then
        if(present(errmsg))then
            write(msg,*)'You want to keep a number of PCs '//&
             & 'greater than ',znkeepmax
            errmsg = sl_errmsg(ierr_error, 'svd', msg)
        endif
        return
    end if
    if(nkeep>ns)then
        if(present(errmsg))then
            write(msg,*)'You want to keep a number of PCs '//&
                &'greater than the number of EOF:',ns
            errmsg = sl_errmsg(ierr_error, 'svd', msg)
        endif
        return
    end if
            
    ! Remove the mean
    ! ---------------
    allocate(zlmean(nslv))
    allocate(zrmean(nsrv))
    zll = ll
    zrr = rr
    zlmean = sum(ll,dim=2)/dble(slvalid(ilselect)) 
    zrmean = sum(rr,dim=2)/dble(srvalid(irselect))
    do it = 1, nt
        zll(:,it) = zll(:,it) - zlmean
        zrr(:,it) = zrr(:,it) - zrmean
    end do
    deallocate(zlmean,zrmean)
    
    ! Standard deviation for correlations
    ! -----------------------------------
    allocate(zls(nslv),zrs(nsrv))
    if(zbcorr==1)then
        zls = sqrt(sum(zll**2,dim=2)/dble(slvalid(ilselect)))
        zrs = sqrt(sum(zrr**2,dim=2)/dble(srvalid(irselect)))
        where(zls==0d0) zls = 1d0
        where(zrs==0d0) zrs = 1d0
    else
        zls = 1d0
        zrs = 1d0
    end if
    
    ! Computations
    ! ============
    
    ! Correlation
    ! -----------
    do i = 1, nt
        zll(:,i) = zll(:,i) / zls
        zrr(:,i) = zrr(:,i) / zrs
    end do
    
    ! Cross-covariances
    ! -----------------
    allocate(cov(nslv,nsrv),nn(nslv,nsrv))
    cov = 0d0
    nn = 0d0
    call dgemm('N', 'T', nslv, nsrv, nt, 1d0, &
        zll, nslv, zrr, nsrv, 0d0, cov, nslv)
    call dgemm('N', 'T', nslv, nsrv, nt, 1d0, &
        dble(lvalid), nslv, dble(rvalid), nsrv, 0d0, cov, nslv)
    where(nn>0d0) cov = cov / nn
    deallocate(nn)
    if(.not.present(lpc)) deallocate(zll, zls)
    if(.not.present(rpc)) deallocate(zrr, zrs)
    
    ! SVD
    ! ---
    allocate(zleof(nslv, ns), zev(ns))
    allocate(work(1))
    call dgesvd('S', 'O', nslv, nsrv, cov, nslv, zev, &
        zleof, nslv, zvt, 1, work, -1, la_info)
    if(la_info/=0)then
        if(present(errmsg))&
            & errmsg = sl_errmsg(ierr_error, 'svd', la_info=la_info, &
                & la_fname='DSYEVD (init phase)')
        return
    endif
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))
    call dgesvd('S', 'O', nslv, nsrv, cov, nslv, zev, zleof, &
        nslv, zvt, 1, work, lwork, la_info)
    if(la_info/=0)then
        if(present(errmsg))&
            & errmsg = sl_errmsg(ierr_error, 'svd', la_info=la_info, &
                & la_fname='DSYEVD')
        return
    endif
    deallocate(work)
    
    ! Get output arrays
    ! =================
    
    ! Sum of all eigenvalues (useful for percentils)
    ! ----------------------------------------------
    if(present(ev_sum)) ev_sum = sum(zev)
    
    ! Eigen values
    ! ------------
    where(zev<0d0) zev=0d0
    if(present(ev)) ev = zev(1:nkeep)
    deallocate(zev)
    
    ! EOFs
    ! ----
    if(present(leof).or.present(lpc))then
        leof(ilselect,:) = zleof(:,1:nkeep)
        do i = 1, nkeep ! First channel of an EOF is >= 0
            if(leof(ilselect(1), i)<0)leof(:, i) = -leof(:, i)
        enddo
        deallocate(zleof)
    end if
    if(present(reof).or.present(rpc))then
        do i = 1, nkeep
            reof(irselect, i) = cov(i,:)
            if(reof(irselect(1), i)<0)reof(:, i) = -reof(:, i)
        end do
        deallocate(cov)
    end if
        
    ! PCs
    ! ---
    if(present(lpc))then
        do i = 1, nt
            zll(:,i) = zll(:,i) * zls ! Correlation case
        end do
        deallocate(zls)
        call sl_pca_getec(ll, leof, lpc, mv=zmv, &
            & minvalid=minecvalid, demean=1)
        deallocate(zll)
    end if
    if(present(rpc))then
        do i = 1, nt
            zrr(:,i) = zrr(:,i) * zrs ! Correlation case
        end do
        deallocate(zrs)
        call sl_pca_getec(rr, reof, rpc, mv=zmv, &
            & minvalid=minecvalid, demean=1)
        deallocate(zrr)
    end if

end subroutine sl_svd


!############################################################
!############################################################
!############################################################


subroutine sl_svd_model_setup(ll, rr, lPcaEof, rPcaEof,&
    & lSvdEof,rSvdEof,l2r,lPcaPc,rPcaPc,lSvdPc,rSvdPc, mv, &
    & minecvalid, errmsg)
    ! Title:
    !    SVD statistical model - Setup part
    !
    ! Description:
    !    Build a SVD-based statistical model to deduce right field
    !    from left field. First, it performs pre-PCA on both dataset,
    !    then it decomposes resulting PCs using a SVD.
    !    Outputs EOFs and PCs from PCA and SVD can further be used
    !    by the model part.
    !
    ! :Necessary arguments:
    !    - ll:    Left space-time array
    !    - rr:    Right space-time array
    !    - mv:    Missing value
    !    - lPcaEof:  Left pre-PCA EOFs
    !    - rPcaEof:  Right pre-PCA EOFs
    !    - lSvdEof:  Left SVD EOFs
    !    - rSvdEof:  Right SVD EOFs
    !    - l2r:      Scale factors to convert from left to right
    !    - lPcsPc:   Left pre-PCA PCs
    !    - rPcsPc:   Right pre-PCA PCs
    !    - lPcsPc:   Left SVD PCs
    !    - rPcsPc:   Right SVD PCs
    
    
    ! Declarations
    ! ============
    
    implicit none
    
    ! External
    ! --------
    real(8), intent(in) :: ll(:,:), rr(:,:)
    real(8), intent(out) ::lPcaEof(:,:), rPcaEof(:,:),&
     & lsvdEof(:,:), rSvdEof(:,:), l2r(:)
    real(8), intent(out), optional :: &
     & lPcaPc(size(ll,2),size(lPcaEof,2)), &
     & rPcaPc(size(ll,2),size(rPcaEof,2)), &
     & lSvdPc(size(ll,2),size(lSvdEof,2)), &
     & rSvdPc(size(ll,2),size(rSvdEof,2)), mv
    integer, intent(in), optional :: minecvalid
    character(len=120), intent(out), optional :: errmsg ! Logging message (len=120)
    
    ! Internal
    ! --------
    integer :: i,nt,nkeepPca, nkeepSvd, nsl, nsr
    
    ! Sizes
    ! -----
    nt = size(ll,2)
    nkeepPca = size(lPcaEof,2)
    nkeepSvd = size(lSvdEof,2)
    nsl = size(ll,1)
    nsr = size(rr,1)
    
    !    ! Missing value
    !    ! -------------
    !    if(present(mv))then
    !        zmv = mv
    !    else
    !        zmv = default_missing_value
    !    endif
    
    ! Computations
    ! ============
    
    ! Pre-PCA
    ! -------
    call sl_pca(ll, nkeepPca, xeof=lPcaEof, pc=lPcaPc, mv=mv, &
        & minecvalid=minecvalid, errmsg=errmsg)
    call sl_pca(rr, nkeepPca, xeof=rPcaEof, pc=rPcaPc, mv=mv, &
        & minecvalid=minecvalid, errmsg=errmsg)
    
    ! SVD
    ! ---
    call sl_svd(transpose(lPcaPc),transpose(rPcaPc),nkeepSvd, &
        & leof=lSvdEof, reof=rSvdEof, lpc=lSvdPc, rpc=rSvdPc, &
        & minecvalid=minecvalid, errmsg=errmsg)
    
    ! Scale factors based on standard deviations
    ! -----------------------------------------
    do i = 1, nkeepSVD
        l2r(i) = sqrt( (sum(rSvdPc(:,i)**2)/dble(nsr) - &
         &             (sum(rSvdPc(:,i))/dble(nsr))**2) / &
         &            (sum(lSvdPc(:,i)**2)/dble(nsr) - &
         &             (sum(lSvdPc(:,i))/dble(nsl))**2))
    end do

end subroutine sl_svd_model_setup


!############################################################
!############################################################
!############################################################


subroutine sl_svd_model_run(ll,rr,&
    & lPcaEof,rPcaEof,lSvdEof,rSvdEof,l2r, mv, minecvalid, errmsg)
    ! Title:
    !    SVD statistical model - Run part
    !
    ! Description:
    !    SVD-based statistical model to deduce right field
    !    from left field. It uses results from pre-PCA
    !    and SVD decompositions performed by sl_svdmodel_build.
    !
    ! :Necessary arguments:
    !    - ll:    Left space array
    !    - rr:    Right space array
    !    - lPcaEof:  Left pre-PCA EOFs
    !    - rPcaEof:  Right pre-PCA EOFs
    !    - lSvdEof:  Left pre-SVD EOFs
    !    - rSvdEof:  Right pre-SVD EOFs
    !    - l2r:      Scale factors to convert from left to right
    !
    ! :Dependencies:
    !    sl_pca_getec sl_pca_rec
    
    
    ! Declarations
    ! ============
    
    implicit none
    
    ! External
    ! --------
    real(8), intent(in) :: ll(:), lPcaEof(:,:), rPcaEof(:,:),&
     & lSvdEof(:,:), rSvdEof(:,:), l2r(:), mv
    real(8), intent(out) :: rr(:)
    integer, intent(in), optional :: minecvalid
    character(len=120), intent(out), optional :: errmsg ! Logging message (len=120)
    
    ! Internal
    ! --------
    integer :: i,nt,nkeepPca, nkeepSvd
    real(8), allocatable :: zlPcaEc(:,:),zlSvdEc(:,:),&
        & zrSvdEc(:,:),zrPcaPc(:,:),zll(:,:),zrr(:,:)
    
    ! Computations
    ! ============
    
    ! Size
    ! ----
    nt = 1
    nkeepPca = size(lPcaEof,2)
    nkeepSvd = size(rSvdEof,2)
    
    ! Get expansion coefficients from re-PCA
    ! --------------------------------------
    allocate(zll(size(ll,1),1),zlPcaEc(nt,nkeepPca))
    zll(:,1) = ll
    call sl_pca_getec(zll,lPcaEof,zlPcaEc, &
        & mv=mv, minvalid=minecvalid, demean=0)
    allocate(zlSvdEc(nt,nkeepSvd))
    call sl_pca_getec(transpose(zlPcaEc),lSvdEof,zlSvdEc, &
        & mv=mv, minvalid=minecvalid, demean=0)
    deallocate(zll,zlPcaEc)
    
    ! Scale factorisation from left to right
    ! --------------------------------------
    allocate(zrSvdEc(nt,nkeepSvd))
    do i = 1, nkeepSvd
        zrSvdEc(:,i) = l2r(i) * zlSvdEc(:,i)
    end do
    deallocate(zlSvdEc)
    
    ! Reconstructions
    ! ---------------
    allocate(zrPcaPc(nkeepPca,nt))
    call sl_pca_rec(rSvdEof,zrSvdEc,zrPcaPc, mv=mv, errmsg=errmsg)
    deallocate(zrSvdEc)
    allocate(zrr(size(rr,1),1))
    call sl_pca_rec(rPcaEof,transpose(zrPcaPc),zrr, mv=mv, errmsg=errmsg)
    rr = zrr(:,1)
    deallocate(zrr,zrPcaPc)

end subroutine sl_svd_model_run


! ############################################################
! ############################################################
! ## TOOLS PART ##############################################
! ############################################################
! ############################################################

function sl_errmsg(ierr, fname, msg, la_info, la_fname)
    ! Generate an error message
    
    ! External
    integer, intent(in), optional :: ierr, & ! Error id
        la_info ! Lapack error id
    character(len=*), intent(in), optional :: fname , &! Routine name
        msg, & ! Error message
        la_fname ! Lapack routine name
    character(len=120) :: sl_errmsg ! Output error message
    
    ! Internal
    character(len=120) :: zla_fname, zmsg, prefix, zla_msg
    
    ! Optional arguments
    if(.not.present(msg))then
        zmsg = ''
    else
        zmsg = ': '//trim(msg)
    endif
    if(.not.present(la_fname))then
        zla_fname = ''
    else
        zla_fname = la_fname
    endif
    sl_errmsg = ''
    
    ! Prefix
    select case(ierr)
        case(ierr_warning)
            sl_errmsg = '[warning] '
        case(ierr_error)
            sl_errmsg = '[error] '
    end select
    
    ! Routine name
    if(present(fname))sl_errmsg = sl_errmsg // fname // ': '
    
    ! Main message
    if(present(msg))sl_errmsg = sl_errmsg // msg
    
    ! Lapack
    zla_msg = ''
    if(present(la_info).and.la_info/=0 .and. &
        & present(la_fname).and.len_trim(la_fname)>0)then
        if(len_trim(sl_errmsg)>0)sl_errmsg = sl_errmsg // ': '
        write(zla_msg,*)'Error running LAPACK routine '//&
            & trim(la_fname)//' (info=', la_info, ')'
        sl_errmsg = sl_errmsg // zla_msg
    endif
                
end function sl_errmsg
        

subroutine sl_phasecomp(varrec, np, phases, offset, firstphase, mv)
    ! Title:
    !    Phase composites
    !
    ! Description:
    !    Performs phase composites of S-T oscillatory field.
    !    This field is typically a reconstructed pair of MSSA modes.
    !    Composites are evaluated according to an index defined by the
    !    first PC of the input field and its derivative.
    !    A minimal normalized amplitude can be also used: when the
    !    index is under value, data are not used to compute phases.
    !    It is also possible so specify the angle of the first phase
    !    in the 360 degrees phase diagram circle: zero means the
    !    the first phase conincides with the maximmum.
    !
    !
    ! :Necessary arguments:
    !    - varrec: Space-time array
    !    - mv: Missing value
    !    - np:    Number of requested phases over the 360 degrees cycle [default:8]
    !
    ! :Optional arguments:
    !    - offset:     Minimal normalized amplitude of the index [default:0.]
    !    - firstphase: Value in degrees of the first phase [default:0.]
    !
    ! :Dependencies:
    !    sl_pca
    
    
    implicit none
    
    ! Declarations
    ! ============
    
    ! External
    ! --------
    integer,       intent(in)           :: np
    real(8), intent(in)           :: varrec(:,:)
    real(8), intent(in), optional :: offset, firstphase, mv
    real(8), intent(out)          :: phases(size(varrec, 1),np)
    
    ! Internal
    ! --------
    real(8), allocatable :: pc(:,:)
    real(8) :: dpc(size(varrec,2)), amp(size(varrec,2))
    integer :: nt, iphase
    real(8) :: angles(np), projection(size(varrec,2))
    real(8) :: pi, deltarad, pcos, psin, zoffset, zfirstphase!, zmv
    logical :: select_amplitude(size(varrec,2)), &
     &         select_phase(size(varrec,2))
    integer :: itime(size(varrec,2)), nsel, i, ns
    integer, allocatable :: isel(:)
    
    
    ! Setup
    ! =====
    nt = size(varrec,2)
    pi = acos(-1d0)
    itime = (/ (i, i=1, nt) /)
    ns = size(varrec, 1)
    if(present(offset))then
        zoffset=offset
    else
        zoffset=0d0
    end if
    !    if(present(mv))then
    !        zmv = mv
    !    else
    !        zmv = default_missing_value
    !    endif
    
    ! Find the first PC and its derivative
    ! ====================================
    allocate(pc(nt,1))
    call sl_pca(varrec, 1, pc=pc, mv=mv)
    pc = pc * sqrt(dble(nt)/sum(pc**2))
    dpc = 0.5d0 * (eoshift(pc(:,1),  1, pc(nt,1)) - &
     &           eoshift(pc(:,1), -1, pc(1,1)))
    dpc((/1,nt/)) = dpc((/1,nt/)) * 2d0
    dpc = dpc * sqrt(dble(nt)/sum(dpc**2))
    amp = sqrt(pc(:,1)**2 + dpc**2)
    
    
    ! Compute the maps
    ! ================
    
    ! Define the marks
    ! ----------------
    deltarad = 2d0 * pi / dble(np)
    if(present(firstphase))then
        zfirstphase = modulo(firstphase * 2d0 * pi / 360d0, 2d0 * pi)
    else
       zfirstphase = 0d0
    end if
    angles = (/ (dble(iphase), iphase=0,np-1) /) * deltarad + &
     &       zfirstphase
    
    ! Compute the phase maps
    ! ----------------------
    phases = 0d0
    select_amplitude = amp >= zoffset
    do iphase = 1, np
        pcos = cos(angles(iphase))
        psin = sin(angles(iphase))
        projection =  (pc(:,1)*pcos+dpc*psin) / amp
        select_phase = ( projection >= cos(0.5d0*deltarad) ) &
         &             .and. select_amplitude
        if(any(select_phase))then
            nsel = count(select_phase)
            allocate(isel(nsel))
            isel = pack(itime, select_phase)
            phases(:,iphase) = sum(varrec(:,isel), dim=2) / &
                & dble(nsel)
            deallocate(isel)
        end if
    end do

end subroutine sl_phasecomp

!subroutine sl_report(text)
!    ! Simply report text to the standard output
!    character(len=*) :: text
!    write(*,*) text
!end subroutine sl_report

end module spanlib
