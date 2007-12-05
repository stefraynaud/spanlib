! File: spanlib_pywrap.f90
!
! This file is part of the SpanLib library.
! Copyright (C) 2006  Stephane Raynaud
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

subroutine pca(ff, ns, nt, nkeep, xeof, pc, ev, ev_sum, weights, useteof)

	use spanlib, only: sl_pca
	use spanlib_precision

	implicit none

	! External
	! --------
	integer, intent(in)  :: ns,nt
	real(wp),intent(in)  :: ff(ns,nt)
	integer, intent(in)  :: nkeep
	real(wp),intent(out) :: pc(nt,nkeep), xeof(ns,nkeep), ev(nkeep)
	real(wp),intent(in)  :: weights(ns)
	integer, intent(in)  :: useteof
	real,    intent(out) :: ev_sum

	! Call to original subroutine
	! ---------------------------
	call sl_pca(ff, nkeep, xeof=xeof, pc=pc, ev=ev, ev_sum=ev_sum,&
	 & weights=weights, useteof=useteof)

end subroutine pca

subroutine pca_getec(ff, xeof, ns, nt, nkept, ec, weights)

	use spanlib, only: sl_pca_getec
	use spanlib_precision

	implicit none

	! External
	! --------
	integer, intent(in)  :: ns,nt,nkept
	real(wp),intent(in)  :: ff(ns,nt), xeof(ns,nkept),weights(ns)
	real(wp),intent(out) :: ec(nt,nkept)


	! Call to original subroutine
	! ---------------------------
	call sl_pca_getec(ff, xeof, ec, weights=weights)

end subroutine pca_getec


subroutine pca_rec(xeof, pc, ns, nt, nkept, ffrec, istart, iend)

	use spanlib, only: sl_pca_rec
	use spanlib_precision

	implicit none

	! External
	! --------
	integer, intent(in)  :: ns, nt, nkept, istart, iend
	real(wp),intent(in)  :: xeof(ns,nkept), pc(nt,nkept)
	real(wp),intent(out) :: ffrec(ns,nt)

	! Call to original subroutine
	! ---------------------------
	call sl_pca_rec(xeof, pc, ffrec=ffrec, istart=istart, iend=iend)

end subroutine pca_rec



subroutine mssa(ff, nchan, nt, nwindow, nkeep, steof, stpc, ev, ev_sum)

	use spanlib, only: sl_mssa
	use spanlib_precision

	implicit none

	! External
	! --------
	integer, intent(in)  :: nchan, nt, nwindow, nkeep
	real(wp),intent(in)  :: ff(nchan,nt)
	real(wp),intent(out) :: steof(nchan*nwindow,nkeep), &
	 & stpc(nt-nwindow+1,nkeep), ev(nkeep)
	real,    intent(out) :: ev_sum

	! Call to original subroutine
	! ---------------------------
	call sl_mssa(ff, nwindow, nkeep, steof=steof, stpc=stpc, ev=ev, &
	 & ev_sum=ev_sum)

end subroutine mssa

subroutine mssa_getec(ff, steof, nchan, nt, nkept, nwindow, stec)

	use spanlib, only: sl_mssa_getec
	use spanlib_precision

	implicit none

	! External
	! --------
	integer, intent(in) :: nchan, nt, nwindow, nkept
	real(wp),intent(in) :: ff(nchan,nt), &
	 & steof(nchan*nwindow,nkept)
	real(wp),intent(out)  :: stec(nt-nwindow+1, nkept)

	! Call to original subroutine
	! ---------------------------
	call sl_mssa_getec(ff, steof, nwindow, stec)

end subroutine mssa_getec

subroutine mssa_rec(steof, stpc, nchan, nt, nkeep, nwindow, &
  & ffrec, istart, iend)

	use spanlib, only: sl_mssa_rec
	use spanlib_precision

	implicit none

	! External
	! --------
	integer, intent(in)  :: nchan, nt, nwindow, nkeep
	integer, intent(in)  :: istart, iend
	real(wp),intent(in)  :: steof(nchan*nwindow,nkeep), &
	 & stpc(nt-nwindow+1,nkeep)
	real(wp),intent(out) :: ffrec(nchan,nt)

	! Call to original subroutine
	! ---------------------------
	call sl_mssa_rec(steof, stpc, nwindow, ffrec=ffrec, &
	 & istart=istart, iend=iend)

end subroutine mssa_rec



subroutine phasecomp(ffrec, ns, nt, np,  phases, weights, &
  & offset, firstphase)

	use spanlib, only: sl_phasecomp
	use spanlib_precision

	implicit none

	! External
	! --------
	integer, intent(in)  :: ns, nt, np
	real(wp),intent(in)  :: ffrec(ns,nt)
	real(wp),intent(in)  :: weights(ns)
	real(wp),intent(in)  :: offset, firstphase
	real(wp),intent(out) :: phases(ns, np)

	! Call to original subroutine
	! ---------------------------
	call sl_phasecomp(ffrec, np, phases, weights=weights, &
	 & offset=offset, firstphase=firstphase)

end subroutine phasecomp



subroutine svd(ll, nsl, rr, nsr, nt, nkeep, leof, reof, lpc, rpc, ev, lweights, rweights)

	use spanlib, only: sl_svd
	use spanlib_precision

	implicit none

	! External
	! --------
	integer, intent(in)  :: nsl,nsr,nt
	real(wp),intent(in)  :: ll(nsl,nt),rr(nsr,nt)
	integer, intent(in)  :: nkeep
	real(wp),intent(out) :: lpc(nt,nkeep), leof(nsl,nkeep)
	real(wp),intent(out) :: rpc(nt,nkeep), reof(nsr,nkeep), ev(nkeep)
	real(wp),intent(in)  :: lweights(nsl),rweights(nsr)

	! Call to original subroutine
	! ---------------------------
	call sl_svd(ll,rr,nkeep,leof,reof,lpc,rpc,ev,lweights,rweights)

end subroutine svd


! Utilities
! =========

subroutine chan_pack(ff3d, mask, nstot, nt, ff2d, ns)

	use spanlib_precision

	implicit none

	! External
	! --------
	integer, intent(in)  :: nstot, nt, ns
	real(wp),intent(in)  :: ff3d(nt,nstot)
	integer, intent(in)  :: mask(nstot)
	real(wp),intent(out) :: ff2d(nt,ns)

	! Internal
	! --------
	integer :: it

	! Call to pack
	! ------------
	do it = 1, nt
		ff2d(it,:) = pack(ff3d(it,:), mask==1)
	end do

end subroutine chan_pack


subroutine chan_unpack(ff3d, mask, nstot, nt, ff2d, ns, &
  & missing_value)

	use spanlib_precision

	implicit none

	! External
	! --------
	integer, intent(in)  :: nstot, nt, ns
	real(wp),intent(out) :: ff3d(nt,nstot)
	integer, intent(in)  :: mask(nstot)
	real(wp),intent(in)  :: ff2d(nt,ns),missing_value

	! Internal
	! --------
	integer :: it

	! Call to pack
	! ------------
	do it = 1, nt
		ff3d(it,:) = unpack(ff2d(it,:), mask==1, missing_value)
	end do

end subroutine chan_unpack

