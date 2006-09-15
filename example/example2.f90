! File: example2.f90
!
! This file is part of the SpanLib library.
! Copyright (C) 2006  Stephane Raynaud
! Contact: stephane dot raynaud at gmail dot com
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

program example2

	! This simple example shows how to use SVD from this package.
	! Warning: it requires netcdf for in/outputs.
	!
	! We start from longitude/latitude/time value of Pacific Sea Surface Temperature
	! that include the El Nino Southern Oscillation signal.
	! Input is the netcdf file data.cdf.
	! We remove land points from the initial array according
	! to the netcdf missing_value attribute of the analysed variable (data are "packed").
	!
	! Initial data set (data2.cdf):
	!  see example1.f90


	use spanlib
	use netcdf

	implicit none

	! Parameters
	! ----------
	integer,parameter :: pcaNkeep=10, svdNkeep=5,&
		& lons1(2)=(/10,40/),lats1(2)=(/12,49/),&
		& lons2(2)=(/70,150/),lats2(2)=(/16,45/)
	real, parameter ::new_missing_value=-999.
	character(len=20), parameter :: input_nc_file="data2.cdf", &
		& output_nc_file="output2.nc", var_name='ssta'

	! Other declarations
	! ------------------
	real, allocatable :: lon1(:), lat1(:), &
		& lon2(:), lat2(:), time(:), &
		& sst1(:,:,:), sst2(:,:,:), sst(:,:,:)
	logical, allocatable :: mask(:,:),mask1(:,:),mask2(:,:)
	real, allocatable :: packed_sst1(:,:),packed_sst2(:,:)
	real, allocatable :: svdEv(:), &
		& pcaEofs1(:,:),pcaEofs2(:,:),pcaPcs1(:,:),pcaPcs2(:,:), &
		& svdEofs1(:,:),svdEofs2(:,:),svdPcs1(:,:),svdPcs2(:,:), &
		& svdEofsRec1(:,:,:), svdEofsRec2(:,:,:), &
		& packed_svdEofsRec1(:,:), packed_svdEofsRec2(:,:)
	character(len=20) :: dim_names(3), dim_name, &
		& lon_units, lat_units, var_units, &
		&	lon_name, lat_name, time_name, time_units
	integer :: ncid, dimid, dimids(4), varid, dims(3), thisdim, &
		& lonid, latid, phaseid, timeid, phcoid, recoid, origid
	integer(kind=4) :: i,nspace,nlon1,nlat1,nlon2,nlat2,ntime,ns1,ns2
	real :: missing_value

	! Get the initial sst field from the netcdf file
	! ----------------------------------------------
	print*,'Reading inputs...'
	call err(nf90_open(input_nc_file, nf90_nowrite, ncid))
	call err(nf90_inq_varid(ncid, var_name, varid))
	! Dimensions
	nlon1 = lons1(2)-lons1(1)+1 ;	nlat1 = lats1(2)-lats1(1)+1
	nlon2 = lons2(2)-lons2(1)+1 ; nlat2 = lats2(2)-lats2(1)+1
	call err(nf90_inquire_variable(ncid, varid, dimids=dimids(1:3)))
	call err(nf90_inquire_dimension(ncid,dimids(1),name=lon_name))
	call err(nf90_inquire_dimension(ncid,dimids(2),name=lat_name))
	call err(nf90_inquire_dimension(ncid, dimids(3), &
		&	name=time_name, len=ntime))
	! Allocations
	allocate(sst1(nlon1,nlat1,ntime),sst2(nlon2,nlat2,ntime))
	allocate(mask1(nlon1,nlat1),mask2(nlon2,nlat2))
	allocate(lon1(nlon1),lat1(nlat1))
	allocate(lon2(nlon2),lat2(nlat2))
	allocate(time(ntime))
	! SST boxes and attributes
	call err(nf90_get_var(ncid, varid, sst1,&
		& start=(/lons1(1),lats1(1),1/), &
		& count=(/lons1(2)-lons1(1)+1,lats1(2)-lats1(1)+1,ntime/)))
	call err(nf90_get_var(ncid, varid, sst2,&
		& start=(/lons2(1),lats2(1)/), &
		& count=(/lons2(2)-lons2(1)+1,lats2(2)-lats2(1)+1/)))
	call err(nf90_get_att(ncid,varid,'missing_value',missing_value))
	call err(nf90_get_att(ncid,varid,'units',var_units))
	! Longitudes
	call err(nf90_inq_varid(ncid, lon_name, varid))
	call err(nf90_get_var(ncid, varid, lon1, &
		& start=(/lons1(1)/), count=(/lons1(2)-lons1(1)+1/)))
	call err(nf90_get_att(ncid, varid, 'units', lon_units))
	call err(nf90_get_var(ncid, varid, lon2, &
		& start=(/lons2(1)/), count=(/lons2(2)-lons1(1)+1/)))
	! Latitudes
	call err(nf90_inq_varid(ncid, lat_name, varid))
	call err(nf90_get_var(ncid, varid, lat1, &
		& start=(/lats1(1)/), count=(/lats1(2)-lats1(1)+1/)))
	call err(nf90_get_att(ncid, varid, 'units', lat_units))
	call err(nf90_get_var(ncid, varid, lat2, &
		& start=(/lats2(1)/), count=(/lats2(2)-lats1(1)+1/)))
	! Time
	call err(nf90_inq_varid(ncid, time_name, varid))
	call err(nf90_get_var(ncid, varid, time))
	call err(nf90_get_att(ncid, varid, 'units', time_units))
	call err(nf90_close(ncid))


	! Format (pack) data to have only one space dimension
	! ---------------------------------------------------
	print*,'Packing...'

	! Now pack
	mask1 = (sst1(:,:,1) /= missing_value)
	mask2 = (sst1(:,:,1) /= missing_value)
	ns1 = count(mask1) ; ns2 = count(mask2)
	allocate(packed_sst1(ns1, ntime))
	allocate(packed_sst2(ns2, ntime))
	do i=1, ntime
		packed_sst1(:,i) = pack(sst1(:,:,i), mask)
		packed_sst2(:,i) = pack(sst2(:,:,i), mask)
	end do

	! Perform a PCA to reduce the d.o.f
	! ---------------------------------
	print*,'Pre-PCA (sl_pca)...'
	allocate(pcaEofs1(ns1, pcaNkeep))
	allocate(pcaPcs1(ntime,pcaNkeep))
	call sl_pca(packed_sst1, pcaNkeep, xeof=pcaEofs1, pc=pcaPcs1)
	allocate(pcaEofs2(ns1, pcaNkeep))
	allocate(pcaPcs2(ntime,pcaNkeep))
	call sl_pca(packed_sst2, pcaNkeep, xeof=pcaEofs2, pc=pcaPcs2)
	deallocate(packed_sst1,packed_sst2)


	! Perform a SVD on previous PCs
	! -----------------------------
	print*,'SVD (sl_svd)...'
	allocate(svdEofs1(pcaNkeep,svdNkeep),svdPcs1(ntime,pcaNkeep))
	allocate(svdEofs2(pcaNkeep,svdNkeep),svdPcs2(ntime,pcaNkeep))
	call sl_svd(pcaPcs1,pcaPcs2,svdNkeep,svdEofs1,svdEofs2,&
		& svdPcs1,svdPcs2,svdEv)
	deallocate(pcaPcs1,pcaPcs2)

	! Swicth SVD EOFs to the physical space (!)
	! -----------------------------------------
	print*,'Back to the physical space (sl_pcarec)...'
	allocate(packed_svdEofsRec1(ns1,svdNkeep))
	allocate(packed_svdEofsRec2(ns2,svdNkeep))
	call sl_pcarec(pcaEofs1, svdEofs1, packed_svdEofsRec1)
	call sl_pcarec(pcaEofs2, svdEofs2, packed_svdEofsRec2)
	deallocate(pcaEofs1,svdEofs1,pcaEofs2,svdEofs2)

	! Unpacking
	! ---------
	print*,'Unpacking...'
	allocate(svdEofsRec1(nlon1,nlat1,svdNkeep))
	allocate(svdEofsRec2(nlon2,nlat2,svdNkeep))
	do i=1, ntime
		svdEofsRec1(:,:,i) = unpack(packed_svdEofsRec1(:,i), &
			& mask1, new_missing_value)
		svdEofsRec2(:,:,i) = unpack(packed_svdEofsRec2(:,i), &
			& mask2, new_missing_value)
		where(mask1 .eq. .false.)svdEofsRec1(:,:,i) = new_missing_value
		where(mask2 .eq. .false.)svdEofsRec2(:,:,i) = new_missing_value
	end do



end program example2

subroutine err(jstatus)

	use netcdf

	integer :: jstatus

	if (jstatus .ne. nf90_noerr) then
		print *, trim(nf90_strerror(jstatus))
		stop
	end if

end subroutine err

