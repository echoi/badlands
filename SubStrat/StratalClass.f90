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
!       Filename:  StratalClass.f90
!
!    Description:  Implements the geometry and topology functions of the stratigraphic grid
!
!        Version:  1.0 
!        Created:  02/06/15 07:56:05
!        Revision:  none
!
!        Author:  Tristan Salles     
!
! =====================================================================================

module stratal_class

  use bilinear
  use topology
  use parallel
  use hydroUtil
  use parameters
  use kdtree2_module
  use kdtree2_precision_module

  implicit none

  ! Stratal grid resolution
  integer::snx,sny,snodes,slay,ilay,layNb
  real(kind=8)::sdx
  real(kind=8)::time_layer

  real(kind=8),dimension(:,:),allocatable::stratadata
  type(kdtree2),pointer::stratatree
  type(kdtree2_result),dimension(10)::strataRslt

  type sediment_parameters
     ! Diameter in millimetres
     real(kind=8)::dia
     ! Density in kg.m-3
     real(kind=8)::dens
     ! Fall velocity
     real(kind=8)::vfall
     ! Marine diffusivity
     real(kind=8)::diffm
     ! Aerial diffusivity
     real(kind=8)::diffa
  end type sediment_parameters
  type(sediment_parameters),dimension(:),allocatable::sediments

  ! Stratigraphic layer
  character(len=128),dimension(:),allocatable::flayh 

  integer,dimension(:),allocatable::lnID
  integer,dimension(:),allocatable::delPt
  integer,dimension(:),allocatable::regPtNb
  integer,dimension(:,:),allocatable::regPtID

  real(kind=8),dimension(:),allocatable::ilayh
!   real(kind=8),dimension(:,:),allocatable::lay_coord
  real(kind=8),dimension(:),allocatable::lay_base
  real(kind=8),dimension(:,:),allocatable::lay_thick
  real(kind=8),dimension(:,:),allocatable::lay_porosity
  real(kind=8),dimension(:,:,:),allocatable::lay_sed
  real(kind=8),dimension(:,:),allocatable::nlay_sed

  ! Active layer parameters for stratigraphic mesh
  real(kind=8)::active_thick
  real(kind=8),dimension(:),allocatable::alay_thick
  real(kind=8),dimension(:),allocatable::alay_porosity
  real(kind=8),dimension(:,:),allocatable::alay_sed

  ! Active layer parameters for delaunay grid
  real(kind=8),dimension(:),allocatable::alay_dthick
  real(kind=8),dimension(:),allocatable::alay_dporosity
  real(kind=8),dimension(:,:),allocatable::alay_dsed

  real,dimension(:),allocatable::biX,biY,interpVal
  real,dimension(:,:),allocatable::biZ
  real,dimension(:,:,:,:),allocatable::biSed

  ! Sediment influx
  real(kind=8),dimension(:,:),allocatable::Qs_inS
  real(kind=8),dimension(:,:),allocatable::change_localS

contains

  ! =====================================================================================

  subroutine constructStrata

    logical::found

    integer::i,j,p,n,k,iu

    real(kind=8)::s
    real(kind=8),dimension(totgrn)::sth
    real(kind=8),dimension(nbnodes,ilay,totgrn)::lay_sed1

    ! Define the stratal dimension
    !snx=int((maxx-minx)/sdx)+1
    !sny=int((maxy-miny)/sdx)+1
    !snodes=snx*sny
    slay=int((time_end-time_start)/time_layer)+ilay

!     if(allocated(stcoord)) deallocate(stcoord)
!     allocate(stcoord(snodes,2))
!     p=0
!     do j=1,sny
!       do i=1,snx
!         p=p+1
!         stcoord(p,1)=minx+(i-1)*sdx
!         stcoord(p,2)=miny+(j-1)*sdx
!       enddo
!     enddo

    ! Read initial deposit
    do k=1,ilay
      if(ilayh(k)==0)then
        inquire(file=flayh(k),exist=found)
        if(.not.found)then
          if(pet_id==0) print*,'Cannot find XmL stratal file ',k
          call mpi_finalize(rc)
        endif
        iu=79
        open(iu,file=flayh(k),status="old",action="read",iostat=rc)
        rewind(iu)
        sth=0.
        do n=1,nbnodes
          read(iu,*)sth(1:totgrn)
          s=sum(sth)
          do p=1,totgrn
            lay_sed1(n,k,p)=sth(p)/s
          enddo
        enddo
      else
        s=ilayh(k)
        do n=1,nbnodes
          lay_sed1(n,k,1:totgrn)=s/totgrn
        enddo
      endif
    enddo

    ! Define coarse regular mesh for interpolation
    if(allocated(biX)) deallocate(biX)
    if(allocated(biY)) deallocate(biY)
    if(allocated(biZ)) deallocate(biZ)
    if(allocated(biSed)) deallocate(biSed)
    allocate(biX(nx+2))
    allocate(biY(ny+2))
    allocate(biZ(nx+2,ny+2))
    allocate(biSed(ilay,totgrn,nx+2,ny+2))

    biX(1)=real(minx-dx)
    do k=2,nx+2
      biX(k)=biX(k-1)+real(dx)
    enddo
    biY(1)=real(miny-dx)
    do k=2,ny+2
      biY(k)=biY(k-1)+real(dx)
    enddo

    p=0
    do j=2,ny+1
      do i=2,nx+1
        p=p+1
        biZ(i,j)=real(coordZ(p))
        do k=1,ilay
          do n=1,totgrn
            biSed(k,n,i,j)=real(lay_sed1(p,k,n))
          enddo
        enddo
      enddo
    enddo
    biZ(1,1)=biZ(2,2)
    biZ(nx+2,1)=biZ(nx+1,2)
    biZ(1,ny+2)=biZ(2,ny+1)

    biZ(nx+2,ny+2)=biZ(nx+1,ny+1)
    biZ(2:nx+1,1)=biZ(2:nx+1,2)
    biZ(2:nx+1,ny+2)=biZ(2:nx+1,ny+1)
    biZ(1,2:nx+1)=biZ(2,2:nx+1)
    biZ(nx+2,2:nx+1)=biZ(nx+1,2:nx+1)

    do k=1,ilay
      do n=1,totgrn
        biSed(k,n,1,1)=biSed(k,n,2,2)
        biSed(k,n,nx+2,1)=biSed(k,n,nx+1,2)
        biSed(k,n,1,ny+2)=biSed(k,n,2,ny+1)
        biSed(k,n,nx+2,ny+2)=biSed(k,n,nx+1,ny+1)
        biSed(k,n,2:nx+1,1)=biSed(k,n,2:nx+1,2)
        biSed(k,n,2:nx+1,ny+2)=biSed(k,n,2:nx+1,ny+1)
        biSed(k,n,1,2:nx+1)=biSed(k,n,2,2:nx+1)
        biSed(k,n,nx+2,2:nx+1)=biSed(k,n,nx+1,2:nx+1)
      enddo
    enddo

    ! Partition grid by rows
!     call partitionStrata 

    ! Allocate the array dimension
    if(allocated(lay_coord)) deallocate(lay_coord)
    if(allocated(lay_base)) deallocate(lay_base)
    allocate(lay_base(upartN))
    allocate(lay_coord(upartN,2))
    if(allocated(lay_thick)) deallocate(lay_thick)
    if(allocated(lay_sed)) deallocate(lay_sed)
    if(allocated(nlay_sed)) deallocate(nlay_sed)
    allocate(lay_thick(upartN,slay))    
    allocate(lay_sed(lsnb,slay,totgrn))
    allocate(nlay_sed(snodes,totgrn))
    if(allocated(interpVal)) deallocate(interpVal)
    allocate(interpVal(lsnb))
   
    ! Build grid X-Y coordinates for fine stratal mesh
    do k=1,upartN
      n=lnID(k)
      lay_coord(k,1)=stcoord(n,1)
      lay_coord(k,2)=stcoord(n,2)
    enddo

    ! Interpolate values from coarse regular grid to fine stratal grid
    call interpolate_grid_bilinear(nx+2,biX,ny+2,biY,biZ,lsnb,real(lay_coord(1:lsnb,1)),real(lay_coord(1:lsnb,2)),interpVal)
    lay_base=interpVal

    do k=1,ilay
      do n=1,totgrn
        call interpolate_grid_bilinear(nx+2,biX,ny+2,biY,biSed(k,n,1:nx+2,1:ny+2),lsnb,real(lay_coord(1:lsnb,1)),real(lay_coord(1:lsnb,2)),interpVal)
        lay_sed(1:lsnb,k,n)=interpVal(1:lsnb)
      enddo
    enddo

    do n=1,lsnb
      lay_thick(n,1:slay)=0.0
      do k=1,ilay
        s=0.0_8
        do p=1,totgrn
          s=s+lay_sed(n,k,p)
        enddo
        lay_thick(n,k)=s
        lay_base(n)=lay_base(n)-s
      enddo
    enddo

    layNb=ilay

    return

  end subroutine constructStrata
  ! =====================================================================================

  subroutine constructStrataReg

    logical::found

    integer::i,j,p,n,k,iu

    real(kind=8)::s
    real(kind=8),dimension(totgrn)::sth
    real(kind=8),dimension(nbnodes,ilay,totgrn)::lay_sed1

    ! Define the stratal dimension
    snx=int((maxx-minx)/sdx)+1
    sny=int((maxy-miny)/sdx)+1
    snodes=snx*sny
    slay=int((time_end-time_start)/time_layer)+ilay

    if(allocated(stcoord)) deallocate(stcoord)
    allocate(stcoord(snodes,2))
    p=0
    do j=1,sny
      do i=1,snx
        p=p+1
        stcoord(p,1)=minx+(i-1)*sdx
        stcoord(p,2)=miny+(j-1)*sdx
      enddo
    enddo

    ! Read initial deposit
    do k=1,ilay
      if(ilayh(k)==0)then
        inquire(file=flayh(k),exist=found)
        if(.not.found)then
          if(pet_id==0) print*,'Cannot find XmL stratal file ',k
          call mpi_finalize(rc)
        endif
        iu=79
        open(iu,file=flayh(k),status="old",action="read",iostat=rc)
        rewind(iu)
        sth=0.
        do n=1,nbnodes
          read(iu,*)sth(1:totgrn)
          s=sum(sth)
          do p=1,totgrn
            lay_sed1(n,k,p)=sth(p)/s
          enddo
        enddo
      else
        s=ilayh(k)
        do n=1,nbnodes
          lay_sed1(n,k,1:totgrn)=s/totgrn
        enddo
      endif
    enddo

    ! Define coarse regular mesh for interpolation
    if(allocated(biX)) deallocate(biX)
    if(allocated(biY)) deallocate(biY)
    if(allocated(biZ)) deallocate(biZ)
    if(allocated(biSed)) deallocate(biSed)
    allocate(biX(nx+2))
    allocate(biY(ny+2))
    allocate(biZ(nx+2,ny+2))
    allocate(biSed(ilay,totgrn,nx+2,ny+2))

    biX(1)=real(minx-dx)
    do k=2,nx+2
      biX(k)=biX(k-1)+real(dx)
    enddo
    biY(1)=real(miny-dx)
    do k=2,ny+2
      biY(k)=biY(k-1)+real(dx)
    enddo

    p=0
    do j=2,ny+1
      do i=2,nx+1
        p=p+1
        biZ(i,j)=real(coordZ(p))
        do k=1,ilay
          do n=1,totgrn
            biSed(k,n,i,j)=real(lay_sed1(p,k,n))
          enddo
        enddo
      enddo
    enddo
    biZ(1,1)=biZ(2,2)
    biZ(nx+2,1)=biZ(nx+1,2)
    biZ(1,ny+2)=biZ(2,ny+1)

    biZ(nx+2,ny+2)=biZ(nx+1,ny+1)
    biZ(2:nx+1,1)=biZ(2:nx+1,2)
    biZ(2:nx+1,ny+2)=biZ(2:nx+1,ny+1)
    biZ(1,2:nx+1)=biZ(2,2:nx+1)
    biZ(nx+2,2:nx+1)=biZ(nx+1,2:nx+1)

    do k=1,ilay
      do n=1,totgrn
        biSed(k,n,1,1)=biSed(k,n,2,2)
        biSed(k,n,nx+2,1)=biSed(k,n,nx+1,2)
        biSed(k,n,1,ny+2)=biSed(k,n,2,ny+1)
        biSed(k,n,nx+2,ny+2)=biSed(k,n,nx+1,ny+1)
        biSed(k,n,2:nx+1,1)=biSed(k,n,2:nx+1,2)
        biSed(k,n,2:nx+1,ny+2)=biSed(k,n,2:nx+1,ny+1)
        biSed(k,n,1,2:nx+1)=biSed(k,n,2,2:nx+1)
        biSed(k,n,nx+2,2:nx+1)=biSed(k,n,nx+1,2:nx+1)
      enddo
    enddo

    ! Partition grid by rows
    call partitionStrata 

    ! Allocate the array dimension
    if(allocated(lay_coord)) deallocate(lay_coord)
    if(allocated(lay_base)) deallocate(lay_base)
    if(allocated(svertDisp)) deallocate(svertDisp)
    allocate(lay_base(lsnb))
    allocate(svertDisp(snodes))
    allocate(lay_coord(lsnb,2))
    if(allocated(lay_thick)) deallocate(lay_thick)
    if(allocated(lay_sed)) deallocate(lay_sed)
    if(allocated(nlay_sed)) deallocate(nlay_sed)
    allocate(lay_thick(lsnb,slay))    
    allocate(lay_sed(lsnb,slay,totgrn))
    allocate(nlay_sed(snodes,totgrn))
    if(allocated(interpVal)) deallocate(interpVal)
    allocate(interpVal(lsnb))
   
    ! Build grid X-Y coordinates for fine stratal mesh
    do k=1,lsnb
      n=lnID(k)
      lay_coord(k,1)=stcoord(n,1)
      lay_coord(k,2)=stcoord(n,2)
    enddo

    ! Interpolate values from coarse regular grid to fine stratal grid
    call interpolate_grid_bilinear(nx+2,biX,ny+2,biY,biZ,lsnb,real(lay_coord(1:lsnb,1)),real(lay_coord(1:lsnb,2)),interpVal)
    lay_base=interpVal

    do k=1,ilay
      do n=1,totgrn
        call interpolate_grid_bilinear(nx+2,biX,ny+2,biY,biSed(k,n,1:nx+2,1:ny+2),lsnb,real(lay_coord(1:lsnb,1)),real(lay_coord(1:lsnb,2)),interpVal)
        lay_sed(1:lsnb,k,n)=interpVal(1:lsnb)
      enddo
    enddo

    do n=1,lsnb
      lay_thick(n,1:slay)=0.0
      do k=1,ilay
        s=0.0_8
        do p=1,totgrn
          s=s+lay_sed(n,k,p)
        enddo
        lay_thick(n,k)=s
        lay_base(n)=lay_base(n)-s
      enddo
    enddo

    layNb=ilay

    return

  end subroutine constructStrataReg
  ! =====================================================================================

  subroutine partitionStrata

    integer::averow,extras,k,rows,i

    integer,dimension(npets)::strat_row,rs,re

    averow=int((sny-1)/npets)
    extras=mod((sny-1),npets)

    ! Get the number of rows for each processor
    do k=0,npets-1
       rows=averow
       if(k<extras) rows=averow+1
       strat_row(k+1)=rows+1
    enddo

    ! Get number of nodes per processor
    lsnb=(strat_row(pet_id+1))*snx
    rs=0
    re=0
    rs(1)=1
    do k=2,npets
       re(k-1)=strat_row(k-1)*snx+rs(k-1)-1
       rs(k)=re(k-1)-snx+1
    enddo
    re(npets)=snodes
    if(allocated(lnID)) deallocate(lnID)
    allocate(lnID(lsnb))

    i=0
    lnID=-1
    do k=1,snodes
       if(k>=rs(pet_id+1).and.k<=re(pet_id+1))then
          i=i+1
          lnID(i)=k
       endif
    enddo
    
    return

  end subroutine partitionStrata
  ! =====================================================================================

  subroutine connectStrata2TIN

    integer::cell,inpoly,gid,k,p,maxnb,n
    real(kind=8),dimension(2)::xy

    if(allocated(delPt)) deallocate(delPt)
    allocate(delPt(snodes))

    if(allocated(regPtNb)) deallocate(regPtNb)
    allocate(regPtNb(dnodes))

    ! Build the stratal tree search
    if(allocated(stratadata))deallocate(stratadata)
    allocate(stratadata(2,dnodes))
    do n=1,dnodes
      stratadata(1,n)=tcoordX(n)
      stratadata(2,n)=tcoordY(n)
    enddo
    stratatree=>kdtree2_create(stratadata,sort=.true.,rearrange=.true.)

    delPt=-1
    
    do k=1,lsnb
      gid=lnID(k)
      xy(1)=stcoord(gid,1)
      xy(2)=stcoord(gid,2)

      ! kdtree
      call kdtree2_n_nearest(stratatree,xy,nn=10,results=strataRslt)
      lp: do p=1,10
        inpoly=-1
        cell=strataRslt(p)%idx
        call pnpoly(gid,cell,inpoly)
        if(inpoly>=0)then 
          delPt(gid)=cell
          exit lp
        endif
        if(p==10)print*,'Problem finding stratal point in TIN'
      enddo lp
    enddo

    call mpi_allreduce(mpi_in_place,delPt,snodes,mpi_integer,mpi_max,badlands_world,rc)

    call kdtree2_destroy(stratatree)

    maxnb=0
    regPtNb=0
    do k=1,snodes
      gid=delPt(k)
      regPtNb(gid)=regPtNb(gid)+1
      maxnb=max(maxnb,regPtNb(gid))
    enddo

    if(allocated(regPtID)) deallocate(regPtID)
    allocate(regPtID(dnodes,maxnb))
    if(allocated(alay_dthick)) deallocate(alay_dthick)
    if(allocated(alay_dsed)) deallocate(alay_dsed)
    allocate(alay_dthick(dnodes))    
    allocate(alay_dsed(dnodes,totgrn))

    regPtNb=0
    do k=1,snodes
      gid=delPt(k)
      regPtNb(gid)=regPtNb(gid)+1
      regPtID(gid,regPtNb(gid))=k
    enddo

    return

  end subroutine connectStrata2TIN
  ! =====================================================================================

  subroutine pnpoly(pt,cell,flag)     
                                            
    integer::np,flag,i,j,cell,pt,id 
    real(kind=8)::x(20),y(20)                                  
    logical::mxf,myf,nxf,nyf           
 
    np=voronoiCell(cell)%vertexNb         
    
    do 1 i=1,np
      id=voronoiCell(cell)%vertexID(i)                                                   
      x(i)=vcoordX(id)-stcoord(pt,1)                                                   
1     y(i)=vcoordY(id)-stcoord(pt,2) 

    flag=-1                                                          
    do 2 i=1,np                                                        
      j=1+mod(i,np)                                                      
      mxf=x(i)>=0.0                                                    
      nxf=x(j)>=0.0                                                    
      myf=y(i)>=0.0                                                    
      nyf=y(j)>=0.0                                                    
      if(.not.((myf.or.nyf).and.(mxf.or.nxf)).or.(mxf.and.nxf)) goto 2       
      if(.not.(myf.and.nyf.and.(mxf.or.nxf).and..not.(mxf.and.nxf))) goto 3  
      flag=-flag  
    goto 2                                                            
3   if((y(i)*x(j)-x(i)*y(j))/(x(j)-x(i))) 2,4,5                       
4   flag=0                                                           
    return                                                            
5   flag=-flag                                                      
2   continue
  
    return                                                            
      
  end subroutine pnpoly         
  ! =====================================================================================                                                     

end module stratal_class
