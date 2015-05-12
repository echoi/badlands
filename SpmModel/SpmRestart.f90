! =====================================================================================
! BADLANDS (BAsin anD LANdscape DynamicS)
!
! Copyright (C) 2015 Tristan Salles 
!
! This program is free software; you can redistribute it and/or modify it under 
! the terms of the GNU General Public License as published by the Free Software 
! Foundation; either version 2 of the License, or (at your option) any later 
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT 
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
! more details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 59 Temple 
! Place, Suite 330, Boston, MA 02111-1307 USA
! ===================================================================================== 

! =====================================================================================
!
!       Filename:  SpmRestart.f90
!
!    Description:  Write SPM mesh HDF5
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles     
!
! =====================================================================================

module spm_restart

  use hdf5
  use parameters
  use topology
  use hydroUtil
  use hydrology

  implicit none

  character(len=128)::frspm

contains

  ! =====================================================================================

  subroutine getSPM_hdf5topography

    logical::compression,simple

    integer::id,i,rank,rpet,localNd,p

    real(kind=8),dimension(:),allocatable::nodes,nID

    character(len=128)::text

    integer(hid_t)::file_id,d_spc
    integer(hid_t)::dset_id,dtype_id
    integer(hsize_t),dimension(2)::dims,maxdims

    do rpet=0,restartPet-1
        frspm=rstfolder
        call noblnk(frspm)
        text='/outputs/SurfaceSPM.'
        call append_str(frspm,text)
        call append_zero(frspm,restartStep)
        text='.p'
        call append_str(frspm,text)
        call append_nb(frspm,rpet)
        text='.h5'
        call append_str(frspm,text)

        ! Initialize predefined datatypes
        call h5open_f(rc)
        call h5zfilter_avail_f(h5z_filter_deflate_f,compression,rc)

        ! Open the file collectively.
        call h5fopen_f(frspm,h5f_acc_rdonly_f,file_id,rc)
        
        ! The node global ID
        text="/nID"
        call h5dopen_f(file_id,trim(text),dset_id,rc)
        call h5dget_type_f(dset_id,dtype_id,rc)
        call h5dget_space_f(dset_id,d_spc,rc)
        call h5sis_simple_f(d_spc,simple,rc)
        call h5sget_simple_extent_ndims_f(d_spc,rank,rc)
        call h5sget_simple_extent_dims_f(d_spc,dims,maxdims,rc)
        allocate(nID(dims(1)*dims(2)))
        localNd=dims(2)
        call h5dread_f(dset_id,h5t_native_double,nID,dims,rc)
        
        ! The Coordinates - vertices
        text="/vertices"
        call h5dopen_f(file_id,trim(text),dset_id,rc)
        call h5dget_type_f(dset_id,dtype_id,rc)
        call h5dget_space_f(dset_id,d_spc,rc)
        call h5sis_simple_f(d_spc,simple,rc)
        call h5sget_simple_extent_ndims_f(d_spc,rank,rc)
        call h5sget_simple_extent_dims_f(d_spc,dims,maxdims,rc)
        allocate(nodes(3*localNd))

        call h5dread_f(dset_id,h5t_native_double,nodes,dims,rc)
        id=1
        do p=1,localNd
            i=int(nID(p))
            if(abs(nodes(id)-tcoordX(i))>0.01.or.abs(nodes(id+1)-tcoordY(i))>0.01)then
                print*,'Problem during restart',i,nodes(id),tcoordX(i),nodes(id+1),tcoordY(i)
                call mpi_finalize(rc)
            endif
            tcoordZ(i)=nodes(id+2)
            id=id+3
        enddo

        ! Close the dataset
        call h5sclose_f(d_spc,rc)
        call h5dclose_f(dset_id,rc)
        
        ! Close the file.
        call h5fclose_f(file_id,rc)

        ! Close interface
        call h5close_f(rc)

        deallocate(nodes,nID)
    enddo

    return

  end subroutine getSPM_hdf5topography
  ! =====================================================================================
  
end module spm_restart
