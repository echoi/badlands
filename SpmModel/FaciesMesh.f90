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
!       Filename:  FaciesMesh.f90
!
!    Description:  Records facies information through time.
!
!        Version:  1.0
!        Created:  16/02/15 11:54:53
!        Revision:  none
!
!        Author:  Tristan Salles     
!
! =====================================================================================

module facies

  use hdf5
  use ESMF
  use topology
  use parameters
  use hydrology
  use hydroUtil
  use FoX_wxml
  use kdtree2_module
  use kdtree2_precision_module

  implicit none

  character(len=128)::fwedge,ffacies
  integer::wedgeNodeNb,wedgeElmtNb

  integer,dimension(:),allocatable::wedgeEVor,outsnode,outselem
  integer,dimension(:,:),allocatable::wedgeElem

  real(ESMF_KIND_R8),dimension(:,:),allocatable::wedgePt

  real(ESMF_KIND_R8),dimension(:),allocatable::interpX,interpY
  real(ESMF_KIND_R8),dimension(:,:),allocatable::interpZ

contains

  ! =====================================================================================
  subroutine interpolate_voronoi_point

    integer::p,k,id1,id2,id3

    do k=1,dnodes
        if(voronoiCell(k)%border==1)then
            stratalZ(1:layerID,k)=spmZ(voronoiCell(k)%bpoint)
        endif
    enddo

    do p=1,layerID
        do k=1,wedgeNodeNb
            ! In case we are on a voronoi point
            if(vorDel(k,1)>0)then
                id1=vorDel(k,1)
                id2=vorDel(k,2)
                id3=vorDel(k,3)
                call DeriveTrianglePlanes(p,interpX(k),interpY(k),id1,id2,id3,interpZ(p,k))
            ! In case we are on a delaunay point
            elseif(vorDel(k,2)==0)then
                id1=vorDel(k,3)
                interpZ(p,k)=stratalZ(p,id1)
            endif
        enddo
    enddo
  
  end subroutine interpolate_voronoi_point
  ! =====================================================================================
  subroutine voronoi_grid_order

    integer::k,l,i,p,n,id,ide,vID,id1,id2,temp(vnodes),totNds
    integer,dimension(:),allocatable::vtemp,etemp
    
    real(ESMF_KIND_R8)::xy(2),xb(3),yb(3)

    type(kdtree2),pointer::Ftree
    type(kdtree2_result),dimension(12)::FRslt
    
    wedgeElmtNb=0
    wedgeNodeNb=0

    ! Define number of wedge nodes and elements
    temp=-1
    totNds=0
    do k=1,upartN
        i=unodeID(k)
        if(voronoiCell(i)%border==0.and.tcoordX(i)>minx.and.tcoordX(i)<maxx.and. &
            tcoordY(i)>miny.and.tcoordY(i)<maxy)then
            do p=1,voronoiCell(i)%vertexNb
                vID=voronoiCell(i)%vertexID(p)
                if(temp(vID)==-1)then
                    wedgeNodeNb=wedgeNodeNb+1
                    temp(vID)=0
                endif
            enddo
            wedgeNodeNb=wedgeNodeNb+1
            totNds=totNds+voronoiCell(i)%vertexNb+1
            wedgeElmtNb=wedgeElmtNb+voronoiCell(i)%vertexNb
        endif
    enddo

    ! Allocate wedges points and local arrays
    if(allocated(wedgePt)) deallocate(wedgePt)
    allocate(wedgePt(wedgeNodeNb,2))
    if(allocated(interpX)) deallocate(interpX)
    allocate(interpX(wedgeNodeNb))
    if(allocated(interpY)) deallocate(interpY)
    allocate(interpY(wedgeNodeNb))
    if(allocated(interpZ)) deallocate(interpZ)
    allocate(interpZ(layerNb,wedgeNodeNb))
    if(allocated(wedgeElem)) deallocate(wedgeElem)
    allocate(wedgeElem(wedgeElmtNb,3))
    if(allocated(wedgeEVor)) deallocate(wedgeEVor)
    allocate(wedgeEVor(wedgeElmtNb))
    if(allocated(vorDel)) deallocate(vorDel)
    allocate(vorDel(wedgeNodeNb,3))
    allocate(vtemp(totNds))
    allocate(etemp(totNds))
    temp=-1
    vtemp=-1
    etemp=-1
    wedgePt=0
    vorDel=-1
    wedgeEVor=-1
    wedgeElem=-1
    interpZ=-10000.

    ! Define voronoi points ID
    id=0
    id1=0
    id2=0
    do k=1,upartN
        i=unodeID(k)
        if(voronoiCell(i)%border==0.and.tcoordX(i)>minx.and.tcoordX(i)<maxx.and. &
            tcoordY(i)>miny.and.tcoordY(i)<maxy)then
            id=id+1
            id1=id1+1
            vtemp(id)=1
            etemp(id)=id1
            do p=1,voronoiCell(i)%vertexNb
                id=id+1
                vID=voronoiCell(i)%vertexID(p)
                if(temp(vID)==-1)then
                    id1=id1+1
                    vtemp(id)=1
                    etemp(id)=id1
                    temp(vID)=id1
                else
                    etemp(id)=temp(vID)
                endif
            enddo
        endif
    enddo

    ! Find position of nodes
    id=0
    ide=0
    id2=0
    do l=1,upartN
        i=unodeID(l)
        if(voronoiCell(i)%border==0.and.tcoordX(i)>minx.and.tcoordX(i)<maxx.and. &
            tcoordY(i)>miny.and.tcoordY(i)<maxy)then
            ! Get the cell points
            ide=ide+1
            id1=ide
            id2=id2+1
            if(vtemp(id2)>0)then 
                id=id+1
                wedgePt(id,1)=tcoordX(i)
                wedgePt(id,2)=tcoordY(i)
                interpX(id)=tcoordX(i)
                interpY(id)=tcoordY(i)
            endif
            wedgeEVor(ide)=i
            k=1
            wedgeElem(ide,k)=etemp(id2)
            ! Get the voronoi points
            do p=1,voronoiCell(i)%vertexNb
                id2=id2+1
                vID=voronoiCell(i)%vertexID(p)
                if(vtemp(id2)>0)then 
                    id=id+1
                    wedgePt(id,1)=vcoordX(vID)
                    wedgePt(id,2)=vcoordY(vID)
                    interpX(id)=vcoordX(vID)
                    interpY(id)=vcoordY(vID)
                endif
                k=k+1
                wedgeElem(ide,k)=etemp(id2)
                wedgeEVor(ide)=i
                if(k==3.and.p<voronoiCell(i)%vertexNb)then
                    k=2
                    ide=ide+1
                    wedgeElem(ide,1)=wedgeElem(ide-1,1)
                    wedgeElem(ide,2)=wedgeElem(ide-1,3)
                    wedgeEVor(ide)=i
                elseif(k==3.and.p==voronoiCell(i)%vertexNb)then
                    ide=ide+1
                    wedgeElem(ide,1)=wedgeElem(ide-1,1)
                    wedgeElem(ide,2)=wedgeElem(ide-1,3)
                    wedgeElem(ide,3)=wedgeElem(id1,2)
                    wedgeEVor(ide)=i
                endif
            enddo
        endif
    enddo

    Ftree=>kdtree2_create(Fdata,sort=.true.,rearrange=.true.)

    id1=0
    temp=-1
    do k=1,upartN
        i=unodeID(k)
        if(voronoiCell(i)%border==0.and.tcoordX(i)>minx.and.tcoordX(i)<maxx.and. &
            tcoordY(i)>miny.and.tcoordY(i)<maxy)then
            id1=id1+1
            vorDel(id1,1:2)=0
            vorDel(id1,3)=i
            do p=1,voronoiCell(i)%vertexNb
                vID=voronoiCell(i)%vertexID(p)
                if(temp(vID)==-1)then
                    id1=id1+1
                    temp(vID)=id1
                    ! Search the tree
                    xy(1)=vcoordX(vID)
                    xy(2)=vcoordY(vID)
                    call kdtree2_n_nearest(Ftree,xy,nn=12,results=FRslt)
                    lo:do n=1,12
                        vorDel(id1,1:3)=delmt(FRslt(n)%idx,1:3)
                        xb(1)=tcoordX(vorDel(id1,1))
                        xb(2)=tcoordX(vorDel(id1,2))
                        xb(3)=tcoordX(vorDel(id1,3))
                        yb(1)=tcoordY(vorDel(id1,1))
                        yb(2)=tcoordY(vorDel(id1,2))
                        yb(3)=tcoordY(vorDel(id1,3))
                        call inside_triangle(xy(1),xy(2),xb,yb,l)
                        if(l==1)then
                            exit lo
                        elseif(n==12)then
                            print*,'Problem when finding voronoi points within Delaunay cell.'
                        endif
                    enddo lo
                endif
            enddo
        endif
    enddo

    call kdtree2_destroy(Ftree)
    deallocate(vtemp,etemp,Fdata)

  end subroutine voronoi_grid_order
  ! =====================================================================================
  subroutine update_facies

    integer::id,p,maxtype

    real(ESMF_KIND_R8)::tstep,maxfac

    ! Store layers information
    do id=1,dnodes
      
      ! Update facies layers
      if(simulation_time>time_start.and.faciesOn==1)then
        tstep=time_step
        if(simulation_time+tstep>cpl2_time) tstep=cpl2_time-simulation_time
        do p=1,layerID
          stratalZ(p,id)=stratalZ(p,id)+tvertDisp(id)*tstep
        enddo      
        ! In case there is deposition add it to top layer 
        if(spmZ(id)>=stratalZ(layerID,id))then
          
          !----------------------------------
          ! Marine and beach facies
          ! Marine facies
          if(gsea%actual_sea-spmZ(id)>=10.)then
            facType(id,1)=facType(id,1)+spmZ(id)-stratalZ(layerID,id)
          ! Deltaic channel
          elseif(gsea%actual_sea-spmZ(id)>=-0.5.and.strahler(id)>3)then
            facType(id,2)=facType(id,2)+spmZ(id)-stratalZ(layerID,id)
          ! Beach facies
          elseif(gsea%actual_sea-spmZ(id)>=-0.5.and.strahler(id)<=3)then
            facType(id,3)=facType(id,3)+spmZ(id)-stratalZ(layerID,id)
          !----------------------------------
          ! Lacustrine facies
          ! Lacustrine facies
          elseif(watercell(id)-spmZ(id)>=5)then
            facType(id,4)=facType(id,4)+spmZ(id)-stratalZ(layerID,id)
          ! Lacustrine channels
          elseif(watercell(id)-spmZ(id)>=-0.5.and.strahler(id)>3)then
            facType(id,5)=facType(id,5)+spmZ(id)-stratalZ(layerID,id)
          ! Lacustrine shore
          elseif(watercell(id)-spmZ(id)>=-0.5.and.strahler(id)<=3)then
            facType(id,6)=facType(id,6)+spmZ(id)-stratalZ(layerID,id)
          !----------------------------------
          ! Stream facies
          ! River facies
          elseif(strahler(id)>3)then
            facType(id,7)=facType(id,7)+spmZ(id)-stratalZ(layerID,id)
          ! Stream
          elseif(strahler(id)>1)then
            facType(id,8)=facType(id,8)+spmZ(id)-stratalZ(layerID,id)
          ! Hillslope
          else
            facType(id,9)=facType(id,9)+spmZ(id)-stratalZ(layerID,id)
          endif

          maxtype=0
          maxfac=0.
          do p=1,faciestype
            if(facType(id,p)>maxfac)then
              maxfac=facType(id,p)
              maxtype=p
            endif
          enddo
          stratalFacies(layerID,id)=maxtype
          stratalZ(layerID,id)=spmZ(id)

        ! In case there is erosion check how much of the stratigraphy is eroded
        elseif(spmZ(id)<stratalZ(layerID,id))then
          ! If all the layer is eroded
          if(spmZ(id)<stratalZ(layerID-1,id))then
            facType(id,1:faciestype)=0.
            stratalFacies(layerID,id)=0
            !----------------------------------
            ! Erosional facies
            ! Erosional river
            if(strahler(id)>3)then
              stratalFacies(layerID,id)=10
            ! Ephemeral stream
            elseif(strahler(id)>1)then
              stratalFacies(layerID,id)=11
            ! Hillslope 
            else
              stratalFacies(layerID,id)=12
            endif
          endif

          stratalZ(layerID,id)=spmZ(id)
          lpp:do p=layerID-1,1,-1
            if(stratalZ(p,id)>stratalZ(p+1,id))then 
              stratalZ(p,id)=stratalZ(p+1,id)
            else
              exit lpp
            endif
          enddo lpp

        endif
      endif
    enddo

  end subroutine update_facies
  ! =====================================================================================
  subroutine wedge_hdf5(iter)

    logical::compression

    integer(ESMF_KIND_I4)::id,i,k,rank,iter,p,id1
    integer(ESMF_KIND_I4),dimension(:),allocatable::connect

    real(ESMF_KIND_R8),dimension(:),allocatable::nodes,cellID,LayID,facID

    character(len=128)::text,file

    integer(hid_t)::file_id,plist_id
    integer(hid_t)::filespace,dset_id
    integer(hsize_t),dimension(2)::dims

    if(iter==0) call voronoi_grid_order
    if(iter==0.and..not.restartFlag) return
    call interpolate_voronoi_point

    fwedge='Facies'

    file=''
    file=fwedge
    call noblnk(file)
    text='.'
    call append_str(file,text)
    call append_zero(file,iter)
    text='.p'
    call append_str(file,text)
    call append_nb(file,pet_id)
    text='.h5'
    call append_str(file,text)
    call addpath1(file)
    allocate(nodes(3*wedgeNodeNb*layerID))
    allocate(connect(6*wedgeElmtNb*(layerID-1)))
    allocate(cellID(wedgeElmtNb*(layerID-1)))
    allocate(LayID(wedgeElmtNb*(layerID-1)))
    allocate(facID(wedgeElmtNb*(layerID-1)))

    ! Create nodes arrays
    id=1
    do p=1,layerID
        do i=1,wedgeNodeNb
           nodes(id)=wedgePt(i,1)
           nodes(id+1)=wedgePt(i,2)
           nodes(id+2)=interpZ(p,i)
           id=id+3
        enddo
    enddo

    ! Initialize predefined datatypes
    call h5open_f(rc)
    call h5zfilter_avail_f(h5z_filter_deflate_f,compression,rc)

    ! Setup file access property list for MPI-IO access.
    call h5pcreate_f( h5p_file_access_f, plist_id,rc)

    ! Create the file collectively.
    call h5fcreate_f(file,h5f_acc_trunc_f,file_id,rc,access_prp=plist_id)

    ! Surface - connectivity
    id=1
    id1=1
    do p=1,layerID-1
        do i=1,wedgeElmtNb
          connect(id)=wedgeElem(i,1)+wedgeNodeNb*(p-1)
          connect(id+1)=wedgeElem(i,2)+wedgeNodeNb*(p-1)
          connect(id+2)=wedgeElem(i,3)+wedgeNodeNb*(p-1)
          connect(id+3)=wedgeElem(i,1)+wedgeNodeNb*p
          connect(id+4)=wedgeElem(i,2)+wedgeNodeNb*p
          connect(id+5)=wedgeElem(i,3)+wedgeNodeNb*p
          id=id+6
          cellID(id1)=dble(wedgeEVor(i))
          LayID(id1)=(p-1)*layer_interval
          facID(id1)=stratalFacies(p+1,wedgeEVor(i))
          id1=id1+1
        enddo
    enddo
    dims(1)=6
    dims(2)=wedgeElmtNb*(layerID-1)
    rank=2
    call h5screate_simple_f(rank,dims,filespace,rc)
    text=''
    text="/connectivity"

    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_create_f,plist_id,rc)
    call h5pset_chunk_f(plist_id,rank,dims,rc)
    call h5pset_deflate_f(plist_id,9,rc)
    dims(1)=1
    dims(2)=wedgeElmtNb*6*(layerID-1)

    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_integer,filespace,dset_id,rc,plist_id)
    
    ! Write the dataset collectively
    call h5dwrite_f(dset_id,h5t_native_integer,connect,dims,rc)
    call h5pclose_f(plist_id,rc)

    ! Close the dataset
    call h5dclose_f(dset_id,rc)
    call h5sclose_f(filespace,rc)

    ! The Coordinates - vertices
    dims(1)=3
    dims(2)=wedgeNodeNb*layerID
    rank=2
    call h5screate_simple_f(rank,dims,filespace,rc)
    text=''
    text="/vertices"

    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_create_f,plist_id,rc)
    call h5pset_deflate_f(plist_id,9,rc )
    call h5pset_chunk_f(plist_id,rank,dims,rc)
    dims(1)=1
    dims(2)=wedgeElmtNb*(layerID-1)

    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,rc,plist_id)

    ! Write the dataset collectively
    call h5dwrite_f(dset_id,h5t_native_double,nodes,dims,rc)
    call h5pclose_f(plist_id,rc)

    ! Close the dataset
    call h5dclose_f(dset_id,rc)
    call h5sclose_f(filespace,rc)

    dims(1)=1
    dims(2)=wedgeElmtNb*(layerID-1)
    rank=2
    call h5screate_simple_f(rank,dims,filespace,rc)
    text=''
    text="/vorID"
    
    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_create_f,plist_id,rc)
    call h5pset_deflate_f(plist_id,9,rc)
    call h5pset_chunk_f(plist_id,rank,dims,rc)
    
    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,rc,plist_id)
    
    ! Write the datase
    call h5dwrite_f(dset_id,h5t_native_double,cellID,dims,rc)
    call h5pclose_f(plist_id,rc)
    
    ! Close the dataset
    call h5dclose_f(dset_id,rc)
    call h5sclose_f(filespace,rc)

    dims(1)=1
    dims(2)=wedgeElmtNb*(layerID-1)
    rank=2
    call h5screate_simple_f(rank,dims,filespace,rc)
    text=''
    text="/LayID"
    
    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_create_f,plist_id,rc)
    call h5pset_deflate_f(plist_id,9,rc)
    call h5pset_chunk_f(plist_id,rank,dims,rc)
    
    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,rc,plist_id)
    
    ! Write the datase
    call h5dwrite_f(dset_id,h5t_native_double,LayID,dims,rc)
    call h5pclose_f(plist_id,rc)
    
    ! Close the dataset
    call h5dclose_f(dset_id,rc)
    call h5sclose_f(filespace,rc)

    dims(1)=1
    dims(2)=wedgeElmtNb*(layerID-1)
    rank=2
    call h5screate_simple_f(rank,dims,filespace,rc)
    text=''
    text="/facID"
    
    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_create_f,plist_id,rc)
    call h5pset_deflate_f(plist_id,9,rc)
    call h5pset_chunk_f(plist_id,rank,dims,rc)
    
    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,rc,plist_id)
    
    ! Write the datase
    call h5dwrite_f(dset_id,h5t_native_double,facID,dims,rc)
    call h5pclose_f(plist_id,rc)
    
    ! Close the dataset
    call h5dclose_f(dset_id,rc)
    call h5sclose_f(filespace,rc)

    ! Close the file.
    call h5fclose_f(file_id,rc)

    ! Close interface
    call h5close_f(rc)

    deallocate(nodes,connect,cellID,facID)

    return

  end subroutine wedge_hdf5
  ! =====================================================================================
  subroutine wedge_xdmf(iter)

    type(xmlf_t)::xf

    integer::ierr
    integer(ESMF_KIND_I4)::iter,totNdsodes,totelems,k

    character(len=128)::str,stg,filename,filename1,filename2,filename3,filename4,file

    call wedge_hdf5(iter)

    if(.not.allocated(outsnode))then
        allocate(outsnode(npets),outselem(npets))
        call mpi_gather(wedgeNodeNb,1,mpi_integer,outsnode,1,mpi_integer,0,badlands_world,ierr)
        call mpi_gather(wedgeElmtNb,1,mpi_integer,outselem,1,mpi_integer,0,badlands_world,ierr)
    endif
    
    if(pet_id==0)then
        fwedge='Facies'
        file=''
        file=fwedge
        call noblnk(file)
        str='.'
        call append_str(file,str)
        call append_zero(file,iter)
        str='.xmf'
        call append_str(file,str)
        call addpath1(file)

        call xml_OpenFile(file,xf)

        ! Header
        call xml_AddDOCTYPE(xf,"Xdmf","Xdmf.dtd")
        call xml_DeclareNamespace(xf,"http://www.w3.org/2001/XInclude","xi")
        call xml_NewElement(xf,"Xdmf")
        call xml_AddAttribute(xf,"Version","2.0")
        call xml_NewElement(xf,"Domain")
        call xml_NewElement(xf,"Grid")
        call xml_AddAttribute(xf,"GridType","Collection")
        call xml_AddAttribute(xf,"CollectionType","Spatial")
        call xml_NewElement(xf,"Time")
        call xml_AddAttribute(xf,"Type","Single")
        call xml_AddAttribute(xf,"Value",time_display)
        call xml_EndElement(xf,"Time")

        do k=1,npets
            totNdsodes=outsnode(k)*layerID
            totelems=outselem(k)*(layerID-1)

            filename=''
            filename=fwedge
            call noblnk(filename)
            str='.'
            call append_str(filename,str)
            call append_zero(filename,iter)
            str='.p'
            call append_str(filename,str)
            call append_nb(filename,k-1)
            str='.h5'
            call append_str(filename,str)
            filename1=filename
            filename2=filename
            filename3=filename
            filename4=filename
            str=':/connectivity'
            call append_str(filename,str)
            str=':/vertices'
            call append_str(filename1,str)
            str=':/vorID'
            call append_str(filename2,str)
            str=':/LayID'
            call append_str(filename3,str)
            str=':/facID'
            call append_str(filename4,str)

            ! Block begin
            call xml_NewElement(xf,"Grid")
            str='StratBlock.'
            call append_zero(str,iter)
            stg='.p'
            call append_str(str,stg)
            call append_zero(str,k-1)
            call xml_AddAttribute(xf,"Name",trim(str))
            call xml_NewElement(xf,"Topology")
            call xml_AddAttribute(xf,"Type","Wedge")
            call xml_AddAttribute(xf,"NumberOfElements",totelems)
            call xml_AddAttribute(xf,"BaseOffset","1")
            call xml_NewElement(xf,"DataItem")
            call xml_AddAttribute(xf,"Format","HDF")
            call xml_AddAttribute(xf,"DataType","Int")
            str=' '
            call append_nb2(str,totelems)
            call append_nb2(str,6)
            call xml_AddAttribute(xf,"Dimensions",trim(str))
            call xml_AddCharacters(xf,trim(filename))
            call xml_EndElement(xf,"DataItem")
            call xml_EndElement(xf,"Topology")

            ! Geometry
            call xml_NewElement(xf,"Geometry")
            call xml_AddAttribute(xf,"Type","XYZ")
            call xml_NewElement(xf,"DataItem")
            call xml_AddAttribute(xf,"Format","HDF")
            call xml_AddAttribute(xf,"NumberType","Float")
            call xml_AddAttribute(xf,"Precision","8")
            str=' '
            call append_nb2(str,3)
            call append_nb2(str,totNdsodes)
            call xml_AddAttribute(xf,"Dimensions",trim(str))
            call xml_AddCharacters(xf,trim(filename1))
            call xml_EndElement(xf,"DataItem")
            call xml_EndElement(xf,"Geometry")

            ! Cell voronoi ID
            str=' '
            call append_nb2(str,totelems)
            call xml_NewElement(xf,"Attribute")
            call xml_AddAttribute(xf,"Type","Scalar")
            call xml_AddAttribute(xf,"Center","Cell")
            call xml_AddAttribute(xf,"Name","Voronoi ID")
            call xml_NewElement(xf,"DataItem")
            call xml_AddAttribute(xf,"Format","HDF")
            call xml_AddAttribute(xf,"NumberType","Float")
            call xml_AddAttribute(xf,"Precision","8")
            call xml_AddAttribute(xf,"Dimensions",trim(str))
            call xml_AddCharacters(xf,trim(filename2))
            call xml_EndElement(xf,"DataItem")
            call xml_EndElement(xf,"Attribute")

            ! Cell voronoi ID
            str=' '
            call append_nb2(str,totelems)
            call xml_NewElement(xf,"Attribute")
            call xml_AddAttribute(xf,"Type","Scalar")
            call xml_AddAttribute(xf,"Center","Cell")
            call xml_AddAttribute(xf,"Name","Layer Time")
            call xml_NewElement(xf,"DataItem")
            call xml_AddAttribute(xf,"Format","HDF")
            call xml_AddAttribute(xf,"NumberType","Float")
            call xml_AddAttribute(xf,"Precision","8")
            call xml_AddAttribute(xf,"Dimensions",trim(str))
            call xml_AddCharacters(xf,trim(filename3))
            call xml_EndElement(xf,"DataItem")
            call xml_EndElement(xf,"Attribute")

            ! Cell facies ID
            str=' '
            call append_nb2(str,totelems)
            call xml_NewElement(xf,"Attribute")
            call xml_AddAttribute(xf,"Type","Scalar")
            call xml_AddAttribute(xf,"Center","Cell")
            call xml_AddAttribute(xf,"Name","Facies Type")
            call xml_NewElement(xf,"DataItem")
            call xml_AddAttribute(xf,"Format","HDF")
            call xml_AddAttribute(xf,"NumberType","Float")
            call xml_AddAttribute(xf,"Precision","8")
            call xml_AddAttribute(xf,"Dimensions",trim(str))
            call xml_AddCharacters(xf,trim(filename4))
            call xml_EndElement(xf,"DataItem")
            call xml_EndElement(xf,"Attribute")

            call xml_EndElement(xf,"Grid")
        enddo

        ! Footer
        call xml_EndElement(xf,"Grid")
        call xml_EndElement(xf,"Domain")
        call xml_EndElement(xf,"Xdmf")
        call xml_Close(xf)
    endif

    return

  end subroutine wedge_xdmf
  ! =====================================================================================
  subroutine visualise_facies_change(iter)

    ! Parameters Declaration
    type(xmlf_t)::xf

    integer::i,iter,it0
    character(len=128)::filename,str,fname

    call wedge_xdmf(iter)

    if(pet_id==0)then
        filename='Facies_series.xdmf'
        call addpath1(filename)
        
        ! Header
        call xml_OpenFile(filename,xf)
        call xml_AddDOCTYPE(xf,"Xdmf","Xdmf.dtd")
        call xml_DeclareNamespace(xf,"http://www.w3.org/2001/XInclude","xi")
        call xml_NewElement(xf,"Xdmf")
        call xml_AddAttribute(xf,"Version","2.0")
        call xml_NewElement(xf,"Domain")
        call xml_NewElement(xf,"Grid")
        call xml_AddAttribute(xf,"GridType","Collection")
        call xml_AddAttribute(xf,"CollectionType","Temporal")
        
        it0=1
        ! Loop over time step
        do i=it0,iter+1
            fname=''
            ! Grid name
            if(i==1.and..not.restartFlag)then
                fname='SurfaceSPM'
            else
                fname=fwedge
            endif
            call noblnk(fname)
            str='.'
            call append_str(fname,str)
            call append_zero(fname,i-1)
            str='.xmf'
            call append_str(fname,str)
            call xml_NewElement(xf,"xi:include")
            call xml_AddAttribute(xf,"href",trim(fname))
            call xml_AddAttribute(xf,"xpointer","xpointer(//Xdmf/Domain/Grid)")
            call xml_EndElement(xf,"xi:include")
        enddo
        
        ! Footer
        call xml_EndElement(xf,"Grid")
        call xml_EndElement(xf,"Domain")
        call xml_EndElement(xf,"Xdmf")
        call xml_Close(xf)
    endif

    ! Record facies for restart
    call record_facies_h5(iter)

    return

  end subroutine visualise_facies_change
  ! =====================================================================================
  subroutine record_facies_h5(iter)

    logical::compression

    integer(ESMF_KIND_I4)::id,id1,i,k,p,rank,iter,totnodes,ierr
    real(ESMF_KIND_R8),dimension(:),allocatable::zstrat,cID,faciesT

    character(len=128)::text,file

    integer(hid_t)::file_id,plist_id
    integer(hid_t)::filespace,dset_id
    integer(hsize_t),dimension(2)::dims

    ffacies='FaciesData'
    file=''
    file=ffacies
    call noblnk(file)
    text='.'
    call append_str(file,text)
    call append_zero(file,iter)
    text='.p'
    call append_str(file,text)
    call append_nb(file,pet_id)
    text='.h5'
    call append_str(file,text)
    call addpath1(file)
    totnodes=upartN

    allocate(zstrat(totnodes*layerID),faciesT(totnodes*layerID))
    allocate(cID(totnodes))

    id=0
    do k=1,upartN
        i=unodeID(k)
        cID(k)=i
        do p=1,layerID
            id=id+1
            zstrat(id)=stratalZ(p,i)
            faciesT(id)=stratalFacies(p,i)
        enddo
    enddo

    ! Initialize predefined datatypes
    call h5open_f(rc)
    call h5zfilter_avail_f(h5z_filter_deflate_f,compression,rc)

    ! Setup file access property list for MPI-IO access.
    call h5pcreate_f(h5p_file_access_f,plist_id,rc)

    ! Create the file collectively.
    call h5fcreate_f(file,h5f_acc_trunc_f,file_id,rc,access_prp=plist_id)
    
    ! Nodes ID
    dims(1)=1
    dims(2)=totnodes
    rank=2
    call h5screate_simple_f(rank,dims,filespace,ierr)
    text=''
    text="/nID"

    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
    call h5pset_deflate_f(plist_id,9,ierr)
    call h5pset_chunk_f(plist_id,rank,dims,ierr)

    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,ierr,plist_id)

    ! Write the dataset collectively
    call h5dwrite_f(dset_id,h5t_native_double,cID,dims,ierr)
    call h5pclose_f(plist_id,ierr)

    ! Close the dataset
    call h5dclose_f(dset_id,ierr)
    call h5sclose_f(filespace,ierr)

    ! Stratal elevation
    dims(1)=layerID
    dims(2)=upartN
    rank=2
    call h5screate_simple_f(rank,dims,filespace,ierr)
    text=''
    text="/zstrat"
    
    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
    call h5pset_deflate_f(plist_id,9,ierr)
    call h5pset_chunk_f(plist_id,rank,dims,ierr)

    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,ierr,plist_id)

    ! Write the dataset collectively
    call h5dwrite_f(dset_id,h5t_native_double,zstrat,dims,ierr)
    call h5pclose_f(plist_id,ierr)

    ! Close the dataset
    call h5dclose_f(dset_id,ierr)
    call h5sclose_f(filespace,ierr)

    ! Facies type
    dims(1)=layerID
    dims(2)=upartN
    rank=2
    call h5screate_simple_f(rank,dims,filespace,ierr)
    text=''
    text="/facies"
    
    ! Create property list for collective dataset write
    call h5pcreate_f(h5p_dataset_create_f,plist_id,ierr)
    call h5pset_deflate_f(plist_id,9,ierr)
    call h5pset_chunk_f(plist_id,rank,dims,ierr)

    ! Create the dataset with default properties
    call h5dcreate_f(file_id,trim(text),h5t_native_double,filespace,dset_id,ierr,plist_id)

    ! Write the dataset collectively
    call h5dwrite_f(dset_id,h5t_native_double,faciesT,dims,ierr)
    call h5pclose_f(plist_id,ierr)

    ! Close the dataset
    call h5dclose_f(dset_id,ierr)
    call h5sclose_f(filespace,ierr)

    ! Close the file.
    call h5fclose_f(file_id,rc)

    ! Close interface
    call h5close_f(rc)

    deallocate(cID,zstrat,faciesT)

  end subroutine record_facies_h5
  ! =====================================================================================

end module facies