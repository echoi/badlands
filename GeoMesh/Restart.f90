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
!       Filename:  Restart.f90
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

module restart

  use hdf5
  use parameters
  use topology
  use hydroUtil
  use parallel
  use external_forces
  use kdtree2_module
  use kdtree2_precision_module

  implicit none


  integer::newnds,rstnodes

  character(len=128)::frspm

  real(kind=8),dimension(:,:),allocatable::newcoord,rstXYZ

contains

  ! =====================================================================================

  subroutine getSPM_hdf5topography3D

    logical::compression,simple

    integer::id,i,rank,rpet,localNd,p

    real(kind=8),dimension(:),allocatable::nodes,nID,sedID
    real(kind=8),dimension(:,:),allocatable::prevNd
    real(kind=8)::x,y

    character(len=128)::text

    integer(hid_t)::file_id,d_spc
    integer(hid_t)::dset_id,dtype_id
    integer(hsize_t),dimension(2)::dims,maxdims

    rstnodes=0

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
        if(rpet==0)then
            allocate(prevNd((restartPet+1)*localNd,4))
        endif
        call h5dread_f(dset_id,h5t_native_double,nodes,dims,rc)

        ! The sediment pile thickness
        text="/cumdz"
        call h5dopen_f(file_id,trim(text),dset_id,rc)
        call h5dget_type_f(dset_id,dtype_id,rc)
        call h5dget_space_f(dset_id,d_spc,rc)
        call h5sis_simple_f(d_spc,simple,rc)
        call h5sget_simple_extent_ndims_f(d_spc,rank,rc)
        call h5sget_simple_extent_dims_f(d_spc,dims,maxdims,rc)
        allocate(sedID(dims(1)*dims(2)))
        call h5dread_f(dset_id,h5t_native_double,sedID,dims,rc)

        ! The last sediment pile thickness for flexure
        if(flexure)then
          text="/flex"
          call h5dopen_f(file_id,trim(text),dset_id,rc)
          call h5dget_type_f(dset_id,dtype_id,rc)
          call h5dget_space_f(dset_id,d_spc,rc)
          call h5sis_simple_f(d_spc,simple,rc)
          call h5sget_simple_extent_ndims_f(d_spc,rank,rc)
          call h5sget_simple_extent_dims_f(d_spc,dims,maxdims,rc)
          allocate(sedloader(dims(1)*dims(2)))
          call h5dread_f(dset_id,h5t_native_double,sedloader,dims,rc)
        endif

        id=1
        do p=1,localNd
            i=int(nID(p))
            rstnodes=max(rstnodes,i)
            prevNd(i,1)=nodes(id)
            prevNd(i,2)=nodes(id+1)
            prevNd(i,3)=nodes(id+2)
            prevNd(i,4)=sedID(p)
            id=id+3
        enddo

        ! Close the dataset
        call h5sclose_f(d_spc,rc)
        call h5dclose_f(dset_id,rc)

        ! Close the file.
        call h5fclose_f(file_id,rc)

        ! Close interface
        call h5close_f(rc)

        deallocate(nodes,nID,sedID)
    enddo

    if(allocated(rstXYZ)) deallocate(rstXYZ)
    if(allocated(newcoord)) deallocate(newcoord)
    allocate(rstXYZ(rstnodes,4))
    allocate(newcoord(rstnodes,2))

    newnds=0
    do id=1,rstnodes
        x=prevNd(id,1)
        y=prevNd(id,2)
        rstXYZ(id,1:4)=prevNd(id,1:4)
        if(x>minx+dx.and.x<maxx-dx.and.y>miny+dx.and.y<maxy-dx)then
            newnds=newnds+1
            newcoord(newnds,1)=x
            newcoord(newnds,2)=y
        endif
    enddo

    deallocate(prevNd)

    return

  end subroutine getSPM_hdf5topography3D
  ! =====================================================================================

  subroutine getRestart_topography

    integer::id,k,rc

    real(kind=8),dimension(2)::txy
    real(kind=8),dimension(2,rstnodes)::Fd

    type(kdtree2),pointer::Ftree
    type(kdtree2_result),dimension(1)::FRslt

    if(pet_id==0)then
        do k=1,rstnodes
            Fd(1,k)=rstXYZ(k,1)
            Fd(2,k)=rstXYZ(k,2)
        enddo
        Ftree=>kdtree2_create(Fd,sort=.true.,rearrange=.true.)

        do k=1,dnodes
            txy(1)=tcoordX(k)
            txy(2)=tcoordY(k)
            call kdtree2_n_nearest(Ftree,txy,nn=1,results=FRslt)
            tcoordZ(k)=rstXYZ(FRslt(1)%idx,3)
            sedthick(k)=rstXYZ(FRslt(1)%idx,4)
            if(tcoordX(k)>minx.and.tcoordX(k)<maxx &
              .and.tcoordY(k)>miny.and.tcoordY(k)<maxy)then
              sedthick(k)=sedthick(k)+100000.
            else
              sedthick(k)=100000.
            endif
        enddo
        call kdtree2_destroy(Ftree)
    endif


    if(allocated(bilinearX)) deallocate(bilinearX)
    if(allocated(bilinearY)) deallocate(bilinearY)
    if(allocated(bilinearV)) deallocate(bilinearV)
    allocate(bilinearX(nx+2))
    allocate(bilinearY(ny+2))
    allocate(bilinearV(nx+2,ny+2))
    if(disp3d)then
        if(allocated(bilinearHx)) deallocate(bilinearHx)
        if(allocated(bilinearHy)) deallocate(bilinearHy)
        allocate(bilinearHx(nx+2,ny+2))
        allocate(bilinearHy(nx+2,ny+2))
    endif
    bilinearX(1)=real(minx-dx)
    do id=2,nx+2
        bilinearX(id)=bilinearX(id-1)+real(dx)
    enddo
    bilinearY(1)=real(miny-dx)
    do id=2,ny+2
        bilinearY(id)=bilinearY(id-1)+real(dx)
    enddo

    if(allocated(rstXYZ)) deallocate(rstXYZ)

    call mpi_bcast(tcoordZ,dnodes,mpi_double_precision,0,badlands_world,rc)

    return

  end subroutine getRestart_topography
  ! =====================================================================================

  subroutine getSPM_hdf5topography

    logical::compression,simple

    integer::id,i,rank,rpet,localNd,p

    real(kind=8),dimension(:),allocatable::nodes,nID,sedID,sedlID

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

        ! The sediment pile thickness
        text="/cumdz"
        call h5dopen_f(file_id,trim(text),dset_id,rc)
        call h5dget_type_f(dset_id,dtype_id,rc)
        call h5dget_space_f(dset_id,d_spc,rc)
        call h5sis_simple_f(d_spc,simple,rc)
        call h5sget_simple_extent_ndims_f(d_spc,rank,rc)
        call h5sget_simple_extent_dims_f(d_spc,dims,maxdims,rc)
        allocate(sedID(dims(1)*dims(2)))
        call h5dread_f(dset_id,h5t_native_double,sedID,dims,rc)

        ! The last sediment pile thickness for flexure
        if(flexure)then
          text="/flex"
          call h5dopen_f(file_id,trim(text),dset_id,rc)
          call h5dget_type_f(dset_id,dtype_id,rc)
          call h5dget_space_f(dset_id,d_spc,rc)
          call h5sis_simple_f(d_spc,simple,rc)
          call h5sget_simple_extent_ndims_f(d_spc,rank,rc)
          call h5sget_simple_extent_dims_f(d_spc,dims,maxdims,rc)
          allocate(sedloader(dims(1)*dims(2)))
          call h5dread_f(dset_id,h5t_native_double,sedloader,dims,rc)
        endif

        id=1
        do p=1,localNd
            i=int(nID(p))
            if(abs(nodes(id)-tcoordX(i))>0.01.or.abs(nodes(id+1)-tcoordY(i))>0.01)then
                print*,'Problem during restart',i,nodes(id),tcoordX(i),nodes(id+1),tcoordY(i)
                call mpi_finalize(rc)
            endif
            tcoordZ(i)=nodes(id+2)
            if(tcoordX(i)>minx.and.tcoordX(i)<maxx &
              .and.tcoordY(i)>miny.and.tcoordY(i)<maxy)then
              sedthick(i)=sedID(p)+100000.
            else
              sedthick(i)=100000.
            endif
            id=id+3
        enddo

        ! Close the dataset
        call h5sclose_f(d_spc,rc)
        call h5dclose_f(dset_id,rc)

        ! Close the file.
        call h5fclose_f(file_id,rc)

        ! Close interface
        call h5close_f(rc)

        deallocate(nodes,nID,sedlID,sedID)
    enddo

    if(allocated(bilinearX)) deallocate(bilinearX)
    if(allocated(bilinearY)) deallocate(bilinearY)
    if(allocated(bilinearV)) deallocate(bilinearV)
    allocate(bilinearX(nx+2))
    allocate(bilinearY(ny+2))
    allocate(bilinearV(nx+2,ny+2))
    if(disp3d)then
        if(allocated(bilinearHx)) deallocate(bilinearHx)
        if(allocated(bilinearHy)) deallocate(bilinearHy)
        allocate(bilinearHx(nx+2,ny+2))
        allocate(bilinearHy(nx+2,ny+2))
    endif
    bilinearX(1)=real(minx-dx)
    do id=2,nx+2
        bilinearX(id)=bilinearX(id-1)+real(dx)
    enddo
    bilinearY(1)=real(miny-dx)
    do id=2,ny+2
        bilinearY(id)=bilinearY(id-1)+real(dx)
    enddo

    return

  end subroutine getSPM_hdf5topography
  ! =====================================================================================

end module restart
