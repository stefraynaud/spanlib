! File: anaxv.f90
!
! This file is part of the SpanLib library.
! Copyright (C) 2012-2015  Stephane Raynaud
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

! This module provides subroutines to read and write anaxv files



subroutine write_anaxv_field(fname, t, x, ns, nt, mv)

    implicit none

    ! Declarations
    ! - external
    integer, intent(in) ::  ns, &       ! Spatial (channel) dim
                            nt, &       ! Temporal dim
                            t(nt)       ! Time (like 2012010100)
    character(len=*) :: fname           ! File name
    real, intent(in) :: x(ns, nt), &    ! Input data
                            mv          ! Missing value
    ! - local
    integer :: it, igap
    real :: zx(ns),zn
    logical :: bad(ns)

    ! Write
    zn = -1.
    open(11,file=trim(fname),form='unformatted')
    do it = 1, nt
        zx = x(:, it)
        bad = zx==mv
        where(bad)zx = sqrt(zn)
        if(all(bad))then
            igap = 0
        elseif(any(bad))then
            igap = 10
        else
            igap = 1
        end if
        write(11) t(it), igap, zx
    end do
    close(11)

end subroutine write_anaxv_field


subroutine read_anaxv_field(fname, t, x, ns, nt, mv)
    ! t,x = read_anaxv_field(fname,ns,nt,mv)

    implicit none

    ! Declarations
    ! - external
    integer, intent(in) ::  ns, &       ! Spatial (channel) dim
                            nt          ! Temporal dim
    character(len=*) :: fname           ! File name
    real, intent(out) :: x(ns, nt)      ! Output data
    integer, intent(out) :: t(nt)       ! Spatial (channel) dim
    real, intent(in) :: mv              ! Missing value
    ! - local
    integer :: it, igap


    ! Read
    open(11,file=trim(fname),form='unformatted',status='old')
    do it = 1, nt
        read(11) t(it), igap, x(:, it)
        if(igap==0)then
            x(:, it) = mv
        elseif(igap==10)then
            where(isnan(x(:, it))) x(:, it)=mv
        end if
    end do
    close(11)

end subroutine read_anaxv_field


subroutine read_anaxv_stcov(fname, cov, nc, nw)
    ! Read a square bloc

    implicit none

    ! Declarations
    ! - external
    integer, intent(in) ::  nc, &           ! Number of channels
                            nw              ! Size of MSSA window
    character(len=*) :: fname               ! File name
    real, intent(out) :: cov(nc*nw, nc*nw)  ! Covariance matrix
    ! - local
    integer :: ic, iw
    real :: covar(nc, nc, nw)

    ! Read
    open(11,file=trim(fname),form='unformatted',status='old')
    do iw = 1, nw
        do ic = 1, nc
            read(11) covar(ic, :, iw)
        end do
    end do
    close(11)

    ! Form
    call stcovar2cov(covar, cov, nc, nw)

end subroutine read_anaxv_stcov


subroutine stcovar2cov(covar,cov,nchan,nw)

    real covar(nchan,nchan,nw)
    integer, intent(in) :: nchan,nw
    real, intent(out) :: cov(nchan*nw,nchan*nw)

    integer :: ilign,icolo,nc1,nc2,nw2,nw1,nc2m,nc1m

    nsteof = nchan*nw
    do nc1=1, nchan
     do nc2=1, nchan
        nc1m=(nc1-1)*nw
        nc2m=(nc2-1)*nw
        do nw2=1, nw
           do nw1=1, nw2
              ilign = nc1m + nw1
              icolo = nc2m + nw2
              idiag = nw2 - nw1 + 1
              cov(ilign, icolo) = covar(nc1, nc2, idiag)
              cov(icolo, ilign) = cov(ilign, icolo)
           enddo
        enddo
     enddo
    enddo
end subroutine stcovar2cov

subroutine read_anaxv_eof(fname, xeof, neof)
    ! EOF(ns,neof=ns) or EOF(nchan*nw,neof=nchan*nw)

    implicit none

    ! Declarations
    ! - external
    integer, intent(in) ::  neof          ! Number of channels and EOFs
    character(len=*) :: fname             ! File name
    real, intent(out) :: xeof(neof, neof) ! Output data


    ! Read
    open(11,file=trim(fname),form='unformatted',status='old')
    read(11)xeof
    close(11)

end subroutine read_anaxv_eof


subroutine eof2steof(xeof, steof, nsteof, nw, nkeep)
    ! EOF(nsteof,nsteof) -> STEOF(neof,nw,nchan)
    ! nsteof = nchan*nw

    implicit none

    ! Declarations
    ! - external
    integer, intent(in) ::  nsteof, &          ! Number EOFs
        &   nw, & ! MSSA window size
        &   nkeep  ! Number of MSSA first EOF to keep
!        &   nchan, & ! NUmber of input channels to MSSA
    real, intent(in) :: xeof(nsteof, nsteof) ! Input eof
    real, intent(out) :: steof(nkeep, nw, nsteof/nw) ! Output ST-EOF
    ! - internal
    integer :: n1, n2, n, ilag, ipc, i, nc

    ! Reshape
    nc = nsteof/nw
    do i=1, nkeep
        do ilag=1, nw
            do ipc=1, nc
              n1 = ilag!nw-ilag+1
              n2 = (ipc-1)*nw
              n = n2+n1
              steof(i, ilag, ipc) = xeof(n, i)
            enddo
        enddo
    enddo

end subroutine eof2steof
