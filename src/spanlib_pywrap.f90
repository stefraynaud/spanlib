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

subroutine pca(ff, ns, nt, nkeep, xeof, pc, ev, weights, useteof)

	use spanlib, only: sl_pca

	implicit none

	! External
	! --------
	integer, intent(in)  :: ns,nt
	real,    intent(in)  :: ff(ns,nt)
	integer, intent(in)  :: nkeep
	real,    intent(out) :: pc(nt,nkeep), xeof(ns,nkeep), ev(nkeep)
	real,    intent(in)  :: weights(ns)
	integer, intent(in)  :: useteof

	! Call to original subroutine
	! ---------------------------
	call sl_pca(ff, nkeep, xeof=xeof, pc=pc, ev=ev, &
	 & weights=weights, useteof=useteof)

end subroutine pca

subroutine pca_getec(ff, xeof, ns, nt, nkept, ec, weights)

	use spanlib, only: sl_pca_getec

	implicit none

	! External
	! --------
	integer, intent(in)  :: ns,nt,nkept
	real,    intent(in)  :: ff(ns,nt), xeof(ns,nkept),weights(ns)
	real,    intent(out) :: ec(nt,nkept)


	! Call to original subroutine
	! ---------------------------
	call sl_pca_getec(ff, xeof, ec, weights=weights)

end subroutine pca_getec


subroutine pca_rec(xeof, pc, ns, nt, nkept, ffrec, istart, iend)

	use spanlib, only: sl_pca_rec

	implicit none

	! External
	! --------
	integer, intent(in)  :: ns, nt, nkept, istart, iend
	real,    intent(in)  :: xeof(ns,nkept), pc(nt,nkept)
	real,    intent(out) :: ffrec(ns,nt)

	! Call to original subroutine
	! ---------------------------
	call sl_pca_rec(xeof, pc, ffrec=ffrec, istart=istart, iend=iend)

end subroutine pca_rec



subroutine mssa(ff, nchan, nt, nwindow, nkeep, steof, stpc, ev)

	use spanlib, only: sl_mssa

	implicit none

	! External
	! --------
	integer,intent(in)  :: nchan, nt, nwindow, nkeep
	real,   intent(in)  :: ff(nchan,nt)
	real,	  intent(out) :: steof(nchan*nwindow,nkeep), &
	 & stpc(nt-nwindow+1,nkeep), ev(nkeep)

	! Call to original subroutine
	! ---------------------------
	call sl_mssa(ff, nwindow, nkeep, steof=steof, stpc=stpc, ev=ev)

end subroutine mssa

subroutine mssa_getec(ff, steof, nchan, nt, nkept, nwindow, stec)

	use spanlib, only: sl_mssa_getec

	implicit none

	! External
	! --------
	integer,intent(in) :: nchan, nt, nwindow, nkept
	real,	  intent(in) :: ff(nchan,nt), &
	 & steof(nchan*nwindow,nkept)
	real, intent(out)  :: stec(nt-nwindow+1, nkept)

	! Call to original subroutine
	! ---------------------------
	call sl_mssa_getec(ff, steof, nwindow, stec)

end subroutine mssa_getec

subroutine mssa_rec(steof, stpc, nchan, nt, nkeep, nwindow, &
  & ffrec, istart, iend)

	use spanlib, only: sl_mssa_rec

	implicit none

	! External
	! --------
	integer,intent(in)  :: nchan, nt, nwindow, nkeep
	integer,intent(in)  :: istart, iend
	real,	  intent(in)  :: steof(nchan*nwindow,nkeep), &
	 & stpc(nt-nwindow+1,nkeep)
	real,   intent(out) :: ffrec(nchan,nt)

	! Call to original subroutine
	! ---------------------------
	call sl_mssa_rec(steof, stpc, nwindow, ffrec=ffrec, &
	 & istart=istart, iend=iend)

end subroutine mssa_rec



subroutine phasecomp(ffrec, ns, nt, np,  phases, weights, &
  & offset, firstphase)

	use spanlib, only: sl_phasecomp

	implicit none

	! External
	! --------
	integer, intent(in)  :: ns, nt, np
	real,    intent(in)  :: ffrec(ns,nt)
	real,    intent(in)  :: weights(ns)
	real,    intent(in)  :: offset, firstphase
	real,    intent(out) :: phases(ns, np)

	! Call to original subroutine
	! ---------------------------
	call sl_phasecomp(ffrec, np, phases, weights=weights, &
	 & offset=offset, firstphase=firstphase)

end subroutine phasecomp



subroutine pack3d(ff3d, mask, ns1, ns2, nt, ff2d, ns)

	implicit none

	! External
	! --------
	integer, intent(in)  :: ns1, ns2, nt, ns
	real,    intent(in)  :: ff3d(ns1,ns2,nt)
	integer, intent(in)  :: mask(ns1,ns2)
	real,    intent(out) :: ff2d(ns,nt)

	! Internal
	! --------
	integer :: it

	! Call to pack
	! ------------
	do it = 1, nt
		ff2d(:,it) = pack(ff3d(:,:,it), mask==1)
	end do

end subroutine pack3d


subroutine unpack3d(ff3d, mask, ns1, ns2, nt, ff2d, ns, &
  & missing_value)

	implicit none

	! External
	! --------
	integer, intent(in)  :: ns1, ns2, nt, ns
	real,    intent(out) :: ff3d(ns1,ns2,nt)
	integer, intent(in)  :: mask(ns1,ns2)
	real,    intent(in)  :: ff2d(ns,nt)
	real                 :: missing_value

	! Internal
	! --------
	integer :: it

	! Call to pack
	! ------------
	do it = 1, nt
		ff3d(:,:,it) = unpack(ff2d(:,it), mask==1, missing_value)
	end do

end subroutine unpack3d
