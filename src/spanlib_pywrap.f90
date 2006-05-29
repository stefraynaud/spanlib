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

subroutine pca(ff, ns, nt, nkeep, eof, pc, ev, weights, useteof)

	use spanlib, only: sl_pca

	implicit none

	! External
	! --------
	integer, intent(in)           :: ns,nt
	real,    intent(in)           :: ff(ns,nt)
	integer, intent(in)           :: nkeep
	real,    intent(out)          :: pc(nt,nkeep), eof(ns,nkeep), ev(nkeep)
	real,    intent(in)           :: weights(ns)
	integer, intent(in), optional :: useteof

	! Internal
	! --------
	real, allocatable :: mypc(:,:), myeof(:,:), myev(:)

	! Call to original subroutine
	! ---------------------------
	call sl_pca(ff, nkeep=nkeep, eof=myeof, pc=mypc, ev=myev, weights=weights, useteof=useteof)
	eof=myeof
	pc=mypc
	ev=myev

end subroutine pca


subroutine pcarec(eof, pc, ns, nt, nkept, ffrec, istart, iend)

	use spanlib, only: sl_pcarec
	
	implicit none

	! External
	! --------
	integer, intent(in)  :: ns, nt, nkept, istart, iend
	real,    intent(in)  :: eof(ns,nkept), pc(nt,nkept)
	real,    intent(out) :: ffrec(ns,nt)

	! Internal
	! --------
	real, allocatable :: myffrec(:,:)

	! Call to original subroutine
	! ---------------------------
	call sl_pcarec(eof, pc, ffrec=myffrec, istart=istart, iend=iend)

	! Recover results
	! ---------------
	ffrec = myffrec

end subroutine pcarec



subroutine mssa_python(ff, nchan, nt, nwindow, nkeep, steof, stpc, ev)
	
	use spanlib, only: sl_mssa
	
	implicit none

	! External
	! --------
	integer,intent(in)  :: nchan, nt, nwindow, nkeep
	real,   intent(in)  :: ff(nchan,nt)
	real,	  intent(out) :: steof(nchan*nwindow,nkeep), stpc(nt-nwindow+1,nkeep), ev(nkeep)

	! Internal
	! --------
	real, allocatable :: mysteof(:,:), mystpc(:,:), myev(:)

	! Call to original subroutine
	! ---------------------------
	call sl_mssa(ff, nwindow, nkeep, steof=mysteof, stpc=mystpc, ev=myev)
		
	! Recover results
	! ---------------
	steof = mysteof
	stpc  = mystpc
	ev    = myev

end subroutine mssa




subroutine mssarec(steof, stpc, nchan, nt, nkeep, nwindow, ffrec, istart, iend)

	use spanlib, only: sl_mssarec
	
	implicit none

	! External
	! --------
	integer,intent(in)           :: nchan, nt, nwindow, nkeep
	integer,intent(in), optional :: istart, iend
	real,	  intent(in)           :: steof(nchan*nwindow,nkeep), stpc(nt-nwindow+1,nkeep)
	real,   intent(out)          :: ffrec(nchan,nt)

	! Internal
	! --------
	real, allocatable :: myffrec(:,:)

	! Call to original subroutine
	! ---------------------------
	call sl_mssarec(steof, stpc, nwindow, ffrec=myffrec, istart=istart, iend=iend)

	! Recover results
	! ---------------
	ffrec = myffrec

end subroutine mssarec




subroutine phasecomp(ffrec, ns, nt, nphases,  phases, weights, offset, firstphase)

	use spanlib, only: sl_phasecomp

	implicit none

	! External
	! --------
	integer, intent(in)           :: ns, nt, nphases
	real,    intent(in)           :: ffrec(ns,nt)
	real,    intent(in)           :: weights(ns)
	real,    intent(in), optional :: offset, firstphase
	real,    intent(out)          :: phases(ns, nphases)

	! Internal
	! --------
	real, allocatable :: myphases(:,:)

	! Call to original subroutine
	! ---------------------------
	call sl_phasecomp(ffrec, nphases, myphases, weights=weights, offset=offset, firstphase=firstphase)

	! Recover results
	! ---------------
	phases = myphases

end subroutine phasecomp

_phyton

subroutine pack3d(ff3d, mask, ns1, ns2, nt, ff2d, ns)

		implicit none

		! External
		! --------
		integer, intent(in)  :: ns1, ns2, nt, ns
		real,    intent(in)  :: ff3d(ns1,ns2,nt)
		logical, intent(in)  :: mask(ns1,ns2)
		real,    intent(out) :: ff2d(ns,nt)

		! Internal
		! --------
		integer :: it

		! Call to pack
		! ------------
		do it = 1, nt
			ff2d(:,it) = pack(ff3d(:,:,it), mask)
		end do

	end subroutine pack3d


	subroutine unpack3d(ff3d, mask, ns1, ns2, nt, ff2d, ns, missing_value)

		implicit none

		! External
		! --------
		integer, intent(in)  :: ns1, ns2, nt, ns
		real,    intent(out) :: ff3d(ns1,ns2,nt)
		logical, intent(in)  :: mask(ns1,ns2)
		real,    intent(in)  :: ff2d(ns,nt)
		real                 :: missing_value

		! Internal
		! --------
		integer :: it

		! Call to pack
		! ------------
		do it = 1, nt
			ff3d(:,:,it) = unpack(ff2d(:,it), mask, missing_value)
		end do

	end subroutine unpack3d
