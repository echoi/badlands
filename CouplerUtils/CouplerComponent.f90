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
!       Filename:  CouplerComponent.f90
!
!    Description:  Defines basic functions for communications between different model outputs.
!
!        Version:  1.0
!        Created:  08/05/15 13:05:05
!        Revision:  none
!
!        Author:  Tristan Salles     
!
! =====================================================================================

module coupling

  use geomesh
  use bilinear
  use parallel
  use topology
  use parameters
  use hydroUtil
  use watershed
  use hydrology
  use underworld
  use earthforces
  !use stratal_class
  use external_forces
  use kdtree2_module
  use kdtree2_precision_module

  implicit none

  integer::new_nodes

  real(kind=8),dimension(:,:),allocatable::record
  real(kind=8),dimension(:),allocatable::nhx,nhy,nhz,uzz
  real(kind=8),dimension(:),allocatable::disp_hx,disp_hy 
  
contains

  ! =====================================================================================

  subroutine bilinearTopo

    integer::k,p,id
    real(kind=8),dimension(dnodes)::ubZ
    real,dimension(upartN)::uxpart,uypart,uval

    if(simulation_time==time_start)then
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
      do k=2,nx+2
        bilinearX(k)=bilinearX(k-1)+real(dx)
      enddo
      bilinearY(1)=real(miny-dx)
      do k=2,ny+2
        bilinearY(k)=bilinearY(k-1)+real(dx)
      enddo
    endif

    id=0
    do k=1,ny+2
      do p=1,nx+2
        id=id+1
        bilinearV(p,k)=real(rcoordZ(id))
      enddo
    enddo

    do k=1,upartN
      id=unodeID(k)
      uxpart(k)=real(tcoordX(id))
      uypart(k)=real(tcoordY(id))
    enddo

    ubZ=-1.e6
    call interpolate_grid_bilinear(nx+2,bilinearX,ny+2,bilinearY,bilinearV,upartN,uxpart,uypart,uval)
    
    do k=1,upartN
      id=unodeID(k)
      ubZ(id)=uval(k)
    enddo

    call mpi_allreduce(ubZ,tcoordZ,dnodes,mpi_double_precision,mpi_max,badlands_world,rc)
    
    return

  end subroutine bilinearTopo
  ! =====================================================================================

  subroutine bilinearRain

    integer::k,p,id
    real(kind=8),dimension(dnodes)::ubZ
    real,dimension(upartN)::uxpart,uypart,uval

    id=0
    do k=1,ny+2
      do p=1,nx+2
        id=id+1
        bilinearV(p,k)=real(rainVal(id))
      enddo
    enddo

    do k=1,upartN
      id=unodeID(k)
      uxpart(k)=real(tcoordX(id))
      uypart(k)=real(tcoordY(id))
    enddo

    ubZ=-1.e6

    call interpolate_grid_bilinear(nx+2,bilinearX,ny+2,bilinearY,bilinearV,upartN,uxpart,uypart,uval)
    
    do k=1,upartN
      id=unodeID(k)
      ubZ(id)=uval(k)
    enddo

    call mpi_allreduce(ubZ,precipitation,dnodes,mpi_double_precision,mpi_max,badlands_world,rc)
    
    return

  end subroutine bilinearRain
  ! =====================================================================================

  subroutine bilinearDisp

    integer::k,p,id
    real(kind=8),dimension(dnodes)::ubZ,ub1,ub2
    real,dimension(upartN)::uxpart,uypart,uval,uv1,uv2

    id=0
    do k=1,ny+2
      do p=1,nx+2
        id=id+1
        bilinearV(p,k)=real(rvertDisp(id))
        if(disp3d)then
          bilinearHx(p,k)=real(rhxDisp(id))
          bilinearHy(p,k)=real(rhyDisp(id))
        endif
      enddo
    enddo

    do k=1,upartN
      id=unodeID(k)
      uxpart(k)=real(tcoordX(id))
      uypart(k)=real(tcoordY(id))
    enddo

    ubZ=-1.e6
    if(disp3d)then
      if(allocated(disp_hx)) deallocate(disp_hx)
      if(allocated(disp_hy)) deallocate(disp_hy)
      allocate(disp_hx(dnodes),disp_hy(dnodes))
      ub1=-1.e6
      ub2=-1.e6
    endif

    if(disp3d)then
      call interpolate_grid_bilinear3(nx+2,bilinearX,ny+2,bilinearY,bilinearV,bilinearHx,bilinearHy, &
        upartN,uxpart,uypart,uval,uv1,uv2)    
    else
      call interpolate_grid_bilinear(nx+2,bilinearX,ny+2,bilinearY,bilinearV,upartN,uxpart,uypart,uval)    
    endif

    do k=1,upartN
      id=unodeID(k)
      ubZ(id)=uval(k)
      ub1(id)=uv1(k)
      ub2(id)=uv2(k)
    enddo

    call mpi_allreduce(ubZ,tvertDisp,dnodes,mpi_double_precision,mpi_max,badlands_world,rc)

    if(disp3d)then
      call mpi_allreduce(ub1,disp_hx,dnodes,mpi_double_precision,mpi_max,badlands_world,rc)
      call mpi_allreduce(ub2,disp_hy,dnodes,mpi_double_precision,mpi_max,badlands_world,rc)
    endif
    
    return

  end subroutine bilinearDisp
  ! =====================================================================================

  subroutine bilinearDispS

    integer::k,p,id
    !real,dimension(lsnb)::uxpart,uypart,uval

    id=0
    do k=1,ny+2
      do p=1,nx+2
        id=id+1
        bilinearV(p,k)=real(rvertDisp(id))
      enddo
    enddo

    !do k=1,lsnb
    !  id=lnID(k)
    !  uxpart(k)=real(stcoord(id,1))
    !  uypart(k)=real(stcoord(id,2))
    !enddo

    !call interpolate_grid_bilinear(nx+2,bilinearX,ny+2,bilinearY,bilinearV,lsnb,uxpart,uypart,uval)    
    !svertDisp=uval
    
    return

  end subroutine bilinearDispS
  ! =====================================================================================

  subroutine bilinearGrid

    ! Get rainfall and displacement
    if(rain_event>0) call rainfall
    if(rain_event>0) call bilinearRain
    ! Define coupling time for rain
    if(rain_event==0) cpl1_time=time_end+1000.

    if(disp%event>0.and..not.udwFlag)then 
      call displacement
    elseif(udwFlag)then
      call SurfaceVTK
      call WaitStepCompletion
      call displacement
    endif    

    if(disp%event>0)then 
      call bilinearDisp
      if(totgrn>0) call bilinearDispS
    endif

    ! Define coupling time for geodynamic
    if(disp%event==0) cpl2_time=time_end+1000.

    return

  end subroutine bilinearGrid
  ! =====================================================================================

  subroutine delaunayInterpolant

    integer::gid,k,n,l,pts(3)
    
    real(kind=8)::xy(2),xb(3),yb(3)
    real(kind=8),dimension(bnbnodes)::nzz

    type(kdtree2),pointer::Ftree
    type(kdtree2_result),dimension(12)::FRslt

    nzz=-1.e6

    Ftree=>kdtree2_create(Fdata,sort=.true.,rearrange=.true.)

    do k=1,spartN
    
      gid=snodeID(k)
      xy(1)=rcoordX(gid)
      xy(2)=rcoordY(gid)
      call kdtree2_n_nearest(Ftree,xy,nn=12,results=FRslt)
    
      lo:do n=1,12 
        pts(1:3)=delmt(FRslt(n)%idx,1:3)
        xb(1)=tcoordX(pts(1))
        xb(2)=tcoordX(pts(2))
        xb(3)=tcoordX(pts(3))
        yb(1)=tcoordY(pts(1))
        yb(2)=tcoordY(pts(2))
        yb(3)=tcoordY(pts(3))
        call insideTriangle(xy,xb,yb,l)
        if(l==1)then
          ! derive elevation from triangle plane equation
          call DeriveTrianglePlanes(xy(1),xy(2),pts(1),pts(2),pts(3),nzz(gid))
          exit lo
        elseif(n==12)then
          print*,'Problem when finding points within Delaunay cell.'
        endif
      enddo lo

    enddo
    
    call mpi_allreduce(nzz,rcoordZ,bnbnodes,mpi_double_precision,mpi_max,badlands_world,rc)

    call kdtree2_destroy(Ftree)
 
    return

  end subroutine delaunayInterpolant
  ! =====================================================================================

  subroutine getEarthData

    if(rain_event>0.and.cpl1_time<=simulation_time)then 
      call rainfall
      call bilinearRain
    endif

    if(disp%event>0.and.cpl2_time<=simulation_time)then 
      if(udwFlag)then
        call delaunayInterpolant
        call SurfaceVTK
        call WaitStepCompletion
        call displacement
        call bilinearDisp
      else
        call displacement
        call bilinearDisp
        if(totgrn>0) call bilinearDispS
      endif
    endif

    return

  end subroutine getEarthData
  ! =====================================================================================

  subroutine mvSpmGrid

    integer::k,id,rmvID(dnodes),n,l,h,pts(3),onodes,tt
    real(kind=8)::nx,ny,d1,d2,d3,w1,w2,w3

    real(kind=8),dimension(2)::txy
    real(kind=8),dimension(3)::txb,tyb,tzb
    real(kind=8),dimension(2,dnodes)::Fd1

    type(kdtree2),pointer::Ftree1
    type(kdtree2_result),dimension(2)::FRslt1

    type(kdtree2),pointer::Ftree2
    type(kdtree2_result),dimension(12)::FRslt2

    if(cpl2_time<=simulation_time)then 

      if(allocated(nhx)) deallocate(nhx)
      if(allocated(nhy)) deallocate(nhy)
      if(allocated(nhz)) deallocate(nhz)
      if(allocated(record)) deallocate(record)
      allocate(nhx(dnodes),nhy(dnodes),nhz(dnodes),record(dnodes,2))
 
      ! Apply the displacement to delaunay points
      do k=1,dnodes
        if(tcoordX(k)>minx.and.tcoordX(k)<maxx.and. &
          tcoordY(k)>miny.and.tcoordY(k)<maxy)then
          nhx(k)=tcoordX(k)+disp_hx(k)
          nhy(k)=tcoordY(k)+disp_hy(k)
          nhz(k)=tcoordZ(k)+tvertDisp(k)
        else
          nhx(k)=tcoordX(k)
          nhy(k)=tcoordY(k)
        endif
        Fd1(1,k)=nhx(k)
        Fd1(2,k)=nhy(k)
      enddo

      tt=0
      rmvID=-1
      do k=1,dnodes
        if(tcoordX(k)==minx.or.tcoordX(k)==maxx.or. &
          tcoordY(k)==miny.or.tcoordY(k)==maxy)then
          n=0
          d1=0.
          do id=1,delaunayVertex(k)%ngbNb
            l=delaunayVertex(k)%ngbID(id)
            if(l>0)then
              if(tcoordX(l)>minx.and.tcoordX(l)<maxx.and. &
                tcoordY(l)>miny.and.tcoordY(l)<maxy)then
                n=n+1
                d1=d1+nhz(l)
              endif
            endif
          enddo
          if(n>=1)then
            nhz(k)=d1/n
          else
            tt=tt+1
            rmvID(tt)=k
          endif
        endif
      enddo

      do h=1,tt
        k=rmvID(h)
        n=0
        d1=0.
        do id=1,delaunayVertex(k)%ngbNb
          l=delaunayVertex(k)%ngbID(id)
          if(l>0)then
            if(tcoordX(l)>minx-dx.and.tcoordX(l)<maxx+dx.and. &
              tcoordY(l)>miny-dx.and.tcoordY(l)<maxy+dx)then
              n=n+1
              d1=d1+nhz(l)
            endif
          endif
        enddo
        if(n>=1)then
          nhz(k)=d1/n
        else
          print*,'problem de merde1',k,delaunayVertex(k)%ngbNb
          print*,tcoordX(k),tcoordY(k)
        endif
      enddo
      
      do k=1,dnodes
        if(tcoordX(k)==minx-dx.and.tcoordX(k)==maxx+dx.and. &
          tcoordY(k)==miny-dx.and.tcoordY(k)==maxy+dx)then
          n=0
          d1=0.
          do id=1,delaunayVertex(k)%ngbNb
            l=delaunayVertex(k)%ngbID(id)
            if(l>0)then
              if(tcoordX(l)>minx-dx.and.tcoordX(l)<maxx+dx.and. &
                tcoordY(l)>miny-dx.and.tcoordY(l)<maxy+dx)then
                n=n+1
                d1=d1+nhz(l)
              endif
            endif
          enddo
          if(n>=1)then
            nhz(k)=d1/n
          else
            print*,'problem de merde2',k,delaunayVertex(k)%ngbNb
          endif
        endif
      enddo
      
      Ftree1=>kdtree2_create(Fd1,sort=.true.,rearrange=.true.)

      onodes=dnodes
      tvertDisp=0.0

      if(allocated(delmt2)) deallocate(delmt2)
      if(allocated(Fdata2)) deallocate(Fdata2)
      allocate(delmt2(delem,3))
      allocate(Fdata2(2,delem))
      do k=1,delem
        delmt2(k,1:3)=delmt(k,1:3)
        Fdata2(1,k)=Fdata(1,k)
        Fdata2(2,k)=Fdata(2,k)
      enddo
      Ftree2=>kdtree2_create(Fdata2,sort=.true.,rearrange=.true.)

      if(disp%mindist==0.0)disp%mindist=dx/2.

      ! Check if some points need to be merged
      rmvID=0
      new_nodes=0
      do k=1,dnodes
        if(rmvID(k)==0)then
          if(nhx(k)>minx+dx.and.nhx(k)<maxx-dx.and. &
            nhy(k)>miny+dx.and.nhy(k)<maxy-dx)then
            txy(1)=nhx(k)
            txy(2)=nhy(k)
            call kdtree2_n_nearest(Ftree1,txy,nn=2,results=FRslt1)
            id=FRslt1(2)%idx
            if(sqrt(FRslt1(2)%dis)>disp%mindist)then
              new_nodes=new_nodes+1
              record(new_nodes,1)=nhx(k)
              record(new_nodes,2)=nhy(k)
            elseif(sqrt(FRslt1(2)%dis)<disp%mindist)then              
              nx=nhx(k)+0.5*(nhx(k)-nhx(id))
              ny=nhy(k)+0.5*(nhy(k)-nhy(id))
              if(nx>minx+dx.and.nx<maxx-dx.and.ny>miny+dx &
                .and.ny<maxy-dx)then
                new_nodes=new_nodes+1
                record(new_nodes,1)=nx
                record(new_nodes,2)=ny
                rmvID(k)=1
                rmvID(id)=1
              endif
            endif
          endif
        endif 
      enddo
      call kdtree2_destroy(Ftree1)

      ! Delete the previous delaunay grid
      call UnstructuredMeshDestroy
      ! Clean geomorphic arrays
      call cleanGeomorpho
      ! Create the new grid
      call remesher

      ! Interpolate new elevation
      if(allocated(uzz))deallocate(uzz)
      allocate(uzz(dnodes))
      uzz=-1.e6
      do k=1,upartN
        id=unodeID(k)
        txy(1)=tcoordX(id)
        txy(2)=tcoordY(id)
        call kdtree2_n_nearest(Ftree2,txy,nn=12,results=FRslt2)
        lo1:do n=1,12 
          pts(1:3)=delmt2(FRslt2(n)%idx,1:3)
          txb(1)=nhx(pts(1))
          txb(2)=nhx(pts(2))
          txb(3)=nhx(pts(3))
          tyb(1)=nhy(pts(1))
          tyb(2)=nhy(pts(2))
          tyb(3)=nhy(pts(3))
          tzb(1)=nhz(pts(1))
          tzb(2)=nhz(pts(2))
          tzb(3)=nhz(pts(3))
          call insideTriangle(txy,txb,tyb,l)
          if(l==1)then
            d1=sqrt((txy(1)-txb(1))**2.+(txy(2)-tyb(1))**2.)
            if(d1<1.0e-2)then 
              uzz(id)=tzb(1)
              goto 10
            endif
            d2=sqrt((txy(1)-txb(2))**2.+(txy(2)-tyb(2))**2.)
            if(d2<1.0e-2)then 
              uzz(id)=tzb(2)
              goto 10
            endif
            d3=sqrt((txy(1)-txb(3))**2.+(txy(2)-tyb(3))**2.)
            if(d3<1.0e-2)then 
              uzz(id)=tzb(3)
              goto 10
            endif
            ! derive elevation from triangle plane equation
            call DeriveTrianglePlanes2(txy,txb,tyb,tzb,uzz(id))
            exit lo1
          elseif(n==12)then
            pts(1:3)=delmt2(FRslt2(1)%idx,1:3)
            txb(1)=nhx(pts(1))
            txb(2)=nhx(pts(2))
            txb(3)=nhx(pts(3))
            tyb(1)=nhy(pts(1))
            tyb(2)=nhy(pts(2))
            tyb(3)=nhy(pts(3))
            tzb(1)=nhz(pts(1))
            tzb(2)=nhz(pts(2))
            tzb(3)=nhz(pts(3))
            d1=sqrt((txy(1)-txb(1))**2.+(txy(2)-tyb(1))**2.)
            if(d1<1.0e-2)then 
              uzz(id)=tzb(1)
              goto 10
            endif
            d2=sqrt((txy(1)-txb(2))**2.+(txy(2)-tyb(2))**2.)
            if(d2<1.0e-2)then 
              uzz(id)=tzb(2)
              goto 10
            endif
            d3=sqrt((txy(1)-txb(3))**2.+(txy(2)-tyb(3))**2.)
            if(d3<1.0e-2)then 
              uzz(id)=tzb(3)
              goto 10
            endif
            ! Compute elevation based on inverse weighted averaged distance
            ! Shepard's method adapted by Franke and Nielson, 1980
            w1=(1/d1)**2.
            w2=(1/d2)**2.
            w3=(1/d3)**2.
            uzz(id)=w1*tzb(1)+w2*tzb(2)+w3*tzb(3)
            uzz(id)=(uzz(id))/(w1+w2+w3)
          endif          
        enddo lo1
10 continue        
      enddo
      call mpi_allreduce(uzz,tcoordZ,dnodes,mpi_double_precision,mpi_max,badlands_world,rc)

      ! Update the rain values
      call bilinearRain

      ! Get next displacement timestep
      if(disp%actual<disp%event)then
        cpl2_time=disp_time(disp%actual+1,1)
      else
        cpl2_time=time_end+1000.
      endif

      update3d=.true.

      call kdtree2_destroy(Ftree2)
    endif

    return

  end subroutine mvSpmGrid
  ! =====================================================================================

  subroutine remesher

    ! Get the new TIN
    if(pet_id==0) call delaunayRebuilt
    call mpi_barrier(badlands_world,rc)

    ! GeoMesh main utilities
    call ReadTriangle
    call DelaunayVoronoiDuality
    call DelaunayBorders

    ! Unstructured Grid Partitioning
    call UnstructureGridPart
    ! Structured Grid Partitioning
    call StructureGridPart

    return

  end subroutine remesher
  ! =====================================================================================

  subroutine delaunayRebuilt

    character(len=128)::TINCfile,OPTC,stg

    integer::iu,n,p,tnb,l

    iu=65
    TINfile='TIN.poly'
    call addpath2(TINfile)
    TINCfile='TIN.poly'
    call addpath2(TINCfile)
    TINCfile(len(TINfile):len(TINfile))=CHAR(0)
    open(iu,file=TINfile,status="replace",action="write",iostat=rc)
    if(rc/=0)then
      print*,'Failed to open Triangle Polygon Shape File'
      call mpi_finalize(rc)
    endif    
    rewind(iu)

    ! A box with four points in 2D, no attributes, one boundary marker.
    tnb=4+2*nx+2*(ny-2)+new_nodes
    write(iu,'(I10,1X,3(I2,1X))') tnb,2,0,1
    write(iu,'(I10,1X,2(F16.3,1X),I2)') 1,minx-dx,miny-dx,1
    write(iu,'(I10,1X,2(F16.3,1X),I2)') 2,minx-dx,maxy+dx,1
    write(iu,'(I10,1X,2(F16.3,1X),I2)') 3,maxx+dx,maxy+dx,1
    write(iu,'(I10,1X,2(F16.3,1X),I2)') 4,maxx+dx,miny-dx,1

    ! Write simulation area boundary
    p=4
    do n=1,nx
       p=p+1
       write(iu,'(I10,1X,2(F16.3,1X),I2)') p,minx+(n-1)*dx,miny,0
       p=p+1
       write(iu,'(I10,1X,2(F16.3,1X),I2)') p,minx+(n-1)*dx,maxy,0
    enddo
    do n=2,ny-1
       p=p+1
       write(iu,'(I10,1X,2(F16.3,1X),I2)') p,minx,miny+(n-1)*dx,0
       p=p+1
       write(iu,'(I10,1X,2(F16.3,1X),I2)') p,maxx,miny+(n-1)*dx,0
    enddo

    ! Nodes from previous TIN after displcament 
    do l=1,new_nodes
      p=p+1
      write(iu,'(I10,1X,2(F16.3,1X),I2)') p,record(l,1:2),0
    enddo

    ! Segments of the outer box
    tnb=4
    write(iu,'(I10,1X,I2)') tnb,1
    do n=1,4
       if(n<4)then
          write(iu,'(I10,1X,I10,1X,I10,1X,I2)') n,n,n+1,1
       else
          write(iu,'(I10,1X,I10,1X,I10,1X,I2)') n,n,1,1
       endif
    enddo

    ! Hole not present
    write(iu,'(I10)') 0
    close(iu)

    ! Conforming Delaunay Triangulation options
    OPTC="-Dq"
    call append_nbreal(OPTC,del_angle)
    stg="a"
    call append_str(OPTC,stg)
    call append_nbreal(OPTC,del_area)
    OPTC(len(OPTC):len(OPTC))=CHAR(0)

    ! C-function Triangle call (Shewchuck's algorithm)
    call trianglegen(OPTC,TINCfile)

    return

  end subroutine delaunayRebuilt
  ! =====================================================================================

  subroutine cleanGeomorpho

    if(allocated(nZ)) deallocate(nZ)
    if(allocated(nH)) deallocate(nH)
    if(allocated(newZ)) deallocate(newZ)
    if(allocated(spmZ)) deallocate(spmZ)
    if(allocated(spmH)) deallocate(spmH)
    if(allocated(Qs_in)) deallocate(Qs_in)
    if(allocated(rcvNb)) deallocate(rcvNb)
    if(allocated(dglbID)) deallocate(dglbID)
    if(allocated(rcvIDs)) deallocate(rcvIDs)
    if(allocated(allocs)) deallocate(allocs)
    if(allocated(delemID)) deallocate(delemID)
    if(allocated(cumDisp)) deallocate(cumDisp)
    if(allocated(filldem)) deallocate(filldem)
    if(allocated(strahler)) deallocate(strahler)
    if(allocated(baselist)) deallocate(baselist)
    if(allocated(intArray)) deallocate(intArray)
    if(allocated(rcvprocNb)) deallocate(rcvprocNb)
    if(allocated(rcvprocID)) deallocate(rcvprocID)
    if(allocated(rcvsendID)) deallocate(rcvsendID)
    if(allocated(watercell)) deallocate(watercell)
    if(allocated(receivers)) deallocate(receivers)
    if(allocated(discharge)) deallocate(discharge)
    if(allocated(indexArray)) deallocate(indexArray)
    if(allocated(donorCount)) deallocate(donorCount)
    if(allocated(stackOrder)) deallocate(stackOrder)
    if(allocated(sendprocID)) deallocate(sendprocID)
    if(allocated(donorsList)) deallocate(donorsList)
    if(allocated(lstackOrder)) deallocate(lstackOrder)
    if(allocated(junctionIDs)) deallocate(junctionIDs)
    if(allocated(ljunctionIDs)) deallocate(ljunctionIDs)
    if(allocated(change_local)) deallocate(change_local)
    if(allocated(precipitation)) deallocate(precipitation)
    if(allocated(localNodesGID)) deallocate(localNodesGID)
    if(allocated(subcatchmentID)) deallocate(subcatchmentID)
    if(allocated(subcatchmentProc)) deallocate(subcatchmentProc)

  end subroutine cleanGeomorpho
  ! =====================================================================================

end module coupling
