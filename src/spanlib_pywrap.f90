! File: spanlib_pywrap.f90
!
! This file is part of the SpanLib library.
! Copyright (C) 2006-2013  Stephane Raynaud
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

! Interface to f90
! ================

subroutine pca(var, ns, nt, nkeep, xeof, pc, ev, ev_sum, &
    & mv, useteof, notpc, minecvalid, zerofill, errmsg)

    use spanlib, only: sl_pca

    implicit none

    ! External
    ! --------
    integer, intent(in)  :: ns,nt
    real(8),    intent(in)  :: var(ns,nt)
    integer, intent(in)  :: nkeep
    real(8),    intent(out) :: pc(nt,nkeep), xeof(ns,nkeep), &
        & ev(nkeep)
    real(8),    intent(in)  :: mv
    real(8),    intent(out) :: ev_sum
    integer, intent(in), optional  :: useteof, notpc
    character(len=120), intent(out), optional :: errmsg
    integer, intent(in), optional :: zerofill, minecvalid
        
    ! Call to original subroutine
    ! ---------------------------
    call sl_pca(var, nkeep, xeof=xeof, pc=pc, ev=ev, ev_sum=ev_sum,&
     & mv=mv, useteof=useteof, notpc=notpc, &
     & minecvalid=minecvalid, zerofill=zerofill, errmsg=errmsg)

end subroutine pca

subroutine pca_getec(var, xeof, ns, nt, nkept, ec, mv, &
    & ev, minvalid, zerofill, demean)

    use spanlib, only: sl_pca_getec

    implicit none

    ! External
    ! --------
    integer, intent(in)  :: ns,nt,nkept
    real(8),    intent(in)  :: var(ns,nt), xeof(ns,nkept), mv
    real(8),    intent(out) :: ec(nt,nkept)
    real(8), intent(in), optional :: ev(nkept)
    integer, intent(in), optional :: zerofill, minvalid, demean


    ! Call to original subroutine
    ! ---------------------------
    call sl_pca_getec(var, xeof, ec, mv=mv, ev=ev, &
        & minvalid=minvalid, zerofill=zerofill, demean=demean)

end subroutine pca_getec


subroutine pca_rec(xeof, pc, ns, nt, nkept, varrec, istart, iend, &
    & mv, errmsg)

    use spanlib, only: sl_pca_rec

    implicit none

    ! External
    ! --------
    integer, intent(in)  :: ns, nt, nkept, istart, iend
    real(8),    intent(in)  :: xeof(ns,nkept), pc(nt,nkept), mv
    real(8),   intent(out) :: varrec(ns,nt)
    character(len=120), intent(out), optional :: errmsg

    ! Call to original subroutine
    ! ---------------------------
    call sl_pca_rec(xeof, pc, varrec, istart=istart, iend=iend, mv=mv, errmsg=errmsg)

end subroutine pca_rec



subroutine mssa(var, nchan, nt, nwindow, nkeep, steof, &
    & stpc, ev, ev_sum, mv, minecvalid, zerofill, errmsg)

    use spanlib, only: sl_mssa

    implicit none

    ! External
    ! --------
    integer, intent(in)  :: nchan, nt, nwindow, nkeep
    real(8),    intent(in)  :: var(nchan,nt), mv
    real(8),    intent(out) :: steof(nchan*nwindow,nkeep), &
     & stpc(nt-nwindow+1,nkeep), ev(nkeep)
    real(8),    intent(out) :: ev_sum
    character(len=120), intent(out), optional :: errmsg
    integer, intent(in), optional :: zerofill, minecvalid

    ! Call to original subroutine
    ! ---------------------------
    call sl_mssa(var, nwindow, nkeep, steof=steof, stpc=stpc, &
        & ev=ev, ev_sum=ev_sum, mv=mv, minecvalid=minecvalid, zerofill=zerofill, &
        & errmsg=errmsg)

end subroutine mssa

subroutine stcov(var, cov, nchan, nt, nwindow, mv)

    use spanlib, only: sl_stcov
    
    implicit none
    
    ! External
    ! --------
    integer, intent(in)  :: nchan, nt, nwindow
    real(8),    intent(in)  :: var(nchan,nt), mv
    real(8),    intent(out) :: cov(nchan*nwindow,nchan*nwindow)

    ! Call to original subroutine
    ! ---------------------------
    call sl_stcov(var, cov, mv=mv)

end subroutine stcov


subroutine mssa_getec(var, steof, nchan, nt, nkept, nwindow, stec, &
    & mv, minvalid, zerofill)

    use spanlib, only: sl_mssa_getec

    implicit none

    ! External
    ! --------
    integer, intent(in) :: nchan, nt, nwindow, nkept
    real(8),    intent(in) :: var(nchan,nt), &
     & steof(nchan*nwindow,nkept), mv
    real(8),    intent(out)  :: stec(nt-nwindow+1, nkept)
    integer, intent(in), optional :: zerofill, minvalid

    ! Call to original subroutine
    ! ---------------------------
    call sl_mssa_getec(var, steof, nwindow, stec, mv=mv, &
        & minvalid=minvalid, zerofill=zerofill)

end subroutine mssa_getec

subroutine mssa_rec(steof, stpc, nchan, nt, nkeep, nwindow, &
  & varrec, istart, iend, mv, ev, errmsg)

    use spanlib, only: sl_mssa_rec

    implicit none

    ! External
    ! --------
    integer, intent(in)  :: nchan, nt, nwindow, nkeep
    integer, intent(in)  :: istart, iend
    real(8),    intent(in)  :: steof(nchan*nwindow,nkeep), &
     & stpc(nt-nwindow+1,nkeep), mv
    real(8),    intent(in), optional  :: ev(nkeep)
    real(8),    intent(out) :: varrec(nchan,nt)
    character(len=120), intent(out), optional :: errmsg

    ! Call to original subroutine
    ! ---------------------------
    call sl_mssa_rec(steof, stpc, nwindow, varrec, &
     & istart=istart, iend=iend, mv=mv, ev=ev, errmsg=errmsg)

end subroutine mssa_rec



subroutine phasecomp(varrec, ns, nt, np,  phases, offset, firstphase)

    use spanlib, only: sl_phasecomp

    implicit none

    ! External
    ! --------
    integer, intent(in)  :: ns, nt, np
    real(8),    intent(in)  :: varrec(ns,nt)
    real(8),    intent(in)  :: offset, firstphase
    real(8)    ,intent(out) :: phases(ns, np)

    ! Call to original subroutine
    ! ---------------------------
    call sl_phasecomp(varrec, np, phases, offset=offset, firstphase=firstphase)

end subroutine phasecomp



subroutine svd(ll, nsl, rr, nsr, nt, nkeep, leof, reof, lpc, rpc, ev, &
    & ev_sum, usecorr, mv, minecvalid, errmsg)

    use spanlib, only: sl_svd

    implicit none

    ! External
    ! --------
    integer, intent(in)  :: nsl,nsr,nt
    real(8),    intent(in)  :: ll(nsl,nt),rr(nsr,nt)
    integer, intent(in)  :: nkeep
    real(8),    intent(out) :: lpc(nt,nkeep), leof(nsl,nkeep)
    real(8),    intent(out) :: rpc(nt,nkeep), reof(nsr,nkeep), ev(nkeep)
    real(8),    intent(in)  :: mv
    integer, intent(in)  :: usecorr
    real(8),    intent(out) :: ev_sum
    integer, intent(in), optional :: minecvalid
    character(len=120), intent(out), optional :: errmsg

    ! Internal
    ! --------

    ! Call to original subroutine
    ! ---------------------------
    call sl_svd(ll, rr, nkeep, leof, reof, lpc, rpc, ev, &
        & ev_sum, usecorr, mv, minecvalid, errmsg)

end subroutine svd


! Utilities
! =========

subroutine chan_pack(varNd, mask, nstot, nt, var2d, ns)
    
    implicit none

    ! External
    ! --------
    integer, intent(in)  :: nstot, nt, ns
    real(8),    intent(in)  :: varNd(nt,nstot)
    integer, intent(in)  :: mask(nstot)
    real(8),    intent(out) :: var2d(ns,nt)

    ! Internal
    ! --------
    integer :: it

    ! Call to pack
    ! ------------
    do it = 1, nt
        var2d(:,it) = pack(varNd(it,:), mask==1)
    end do

end subroutine chan_pack


subroutine chan_unpack(varNd, mask, nstot, nt, var2d, ns, &
  & missing_value)

    implicit none
    
    ! External
    ! --------
    integer, intent(in)  :: nstot, nt, ns
    real(8),    intent(out) :: varNd(nt,nstot)
    integer, intent(in)  :: mask(nstot)
    real(8),    intent(in)  :: var2d(ns,nt),missing_value

    ! Internal
    ! --------
    integer :: it

    ! Call to pack
    ! ------------
    do it = 1, nt
        varNd(it,:) = unpack(var2d(:,it), mask==1, missing_value)
    end do

end subroutine chan_unpack

