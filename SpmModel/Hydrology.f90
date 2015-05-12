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
!       Filename:  Hydrology.f90
!
!    Description:  Compute morphometrics based on hydrological features
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles     
!
! =====================================================================================

module hydrology

  use sorting
  use topology
  use parameters
  use hydroUtil
  use external_forces

  implicit none
   
  ! Depression algorithm
  logical::depressionAlgo

  ! Local minima number
  integer::minimaNb,nonIsoNb

  ! Local arrays
  integer,dimension(:),allocatable::intArray,allocs,donorCount,trueCatch
  integer,dimension(:),allocatable::catchmentID,sillNb
  integer,dimension(:,:),allocatable::sillGID 
  
  real(kind=8),dimension(:),allocatable::cumDisp,watercell
  real(kind=8),dimension(:),allocatable::nZ,nH,change_local,filldem,sillMax
  
contains

  ! =====================================================================================

  subroutine define_landscape_network

    integer::k,p,j,lowestID,maxtree,success

    if(.not.allocated(receivers)) allocate(receivers(dnodes))
    if(.not.allocated(trueCatch)) allocate(trueCatch(dnodes))
    if(.not.allocated(baselist)) allocate(baselist(dnodes))
    if(.not.allocated(donorCount)) allocate(donorCount(dnodes))
    donorCount=0
    receivers=-1
    baselist=-1
    trueCatch=-1
    baseNb=0

    do k=1,dnodes
       receivers(k)=k
       lowestID=k
       do p=1,delaunayVertex(k)%ngbNb
          if(delaunayVertex(k)%ngbID(p)>0)then
            if(filldem(delaunayVertex(k)%ngbID(p))<filldem(lowestID))lowestID=delaunayVertex(k)%ngbID(p)
          endif
       enddo
       receivers(k)=lowestID
       ! Baselevel
       if(receivers(k)==k)then
          baseNb=baseNb+1
          baselist(baseNb)=k
       endif
    enddo

    do k=1,dnodes
      donorCount(receivers(k))=donorCount(receivers(k))+1
    enddo

    ! Index of donors number
    if(.not.allocated(indexArray)) allocate(indexArray(dnodes+1))
    indexArray=0
    maxtree=0
    indexArray(dnodes+1)=dnodes+1
    do k=dnodes,1,-1
       indexArray(k)=indexArray(k+1)-donorCount(k) 
       maxtree=max(maxtree,indexArray(k+1)-indexArray(k))       
    enddo
    maxrcvs=maxtree

    ! List of donors
    if(.not.allocated(intArray)) allocate(intArray(dnodes))
    if(.not.allocated(donorsList)) allocate(donorsList(dnodes))
    intArray=0
    do k=1,dnodes
       donorsList(indexArray(receivers(k))+intArray(receivers(k)))=k
       intArray(receivers(k))=intArray(receivers(k))+1
    enddo

    ! Get the donor information array
    if(allocated(donorInfo))deallocate(donorInfo)
    allocate(donorInfo(dnodes,maxtree))
    donorInfo=0
    do k=1,dnodes
       do p=1,indexArray(k+1)-indexArray(k)
          donorInfo(k,p)=donorsList(indexArray(k)+p-1)
       enddo
    enddo

    ! Build the ordering stack
    if(.not.allocated(catchmentID)) allocate(catchmentID(dnodes))
    if(.not.allocated(stackOrder)) allocate(stackOrder(dnodes))
    if(.not.allocated(allocs)) allocate(allocs(dnodes))

    if(allocated(sillGID)) deallocate(sillGID)
    if(allocated(sillMax)) deallocate(sillMax)
    if(allocated(sillNb)) deallocate(sillNb)
    allocate(sillGID(baseNb,dnodes))
    allocate(sillNb(baseNb))
    allocate(sillMax(baseNb))

    allocs=-1
    j=0
    stackOrder=-1
    catchmentID=0
    minimaNb=0
    nonIsoNb=0
    sillNb=0
!     sillGID=0
    sillMax=-1.e6

    do p=1,baseNb
      allocs=-1
      k=baselist(p)
      j=j+1
      stackOrder(j)=k
      allocs(k)=0
      catchID=k
      sillNb(p)=0
      minimaNb=minimaNb+1 
      catchID=-k
      if(voronoiCell(k)%border<1) depressionAlgo=.true.
      if(catchID>0)then
        nonIsoNb=nonIsoNb+1
        trueCatch(nonIsoNb)=catchID
      endif
      catchmentID(k)=catchID
      if(catchID<0)then
        sillMax(p)=max(sillMax(p),filldem(p))
        sillNb(p)=sillNb(p)+1
        sillGID(p,sillNb(p))=k
      endif
      success=addtostack(p,k,j)
    enddo
    
  end subroutine define_landscape_network
  ! =====================================================================================

  recursive function addtostack(base,donor,stackID) result(success)

    integer::base,donor,stackID,n,success

    success=1
    do n=indexArray(donor),indexArray(donor+1)-1
      if((allocs(donorsList(n))==-1))then
        donor=donorsList(n)
        stackID=stackID+1
        stackOrder(stackID)=donor
        allocs(donor)=0
        catchmentID(donor)=catchID
        if(catchID<0)then
          sillNb(base)=sillNb(base)+1
          sillGID(base,sillNb(base))=donor
          sillMax(base)=max(sillMax(base),filldem(base))
        endif
        success=addtostack(base,donor,stackID)
      endif
    enddo 

    success=0

  end function addtostack
  ! =====================================================================================

  subroutine define_drained_water_thickness

    integer::c,q,cID,n,p,no,k,nid,qo
    integer,dimension(:),allocatable::ptOrder,minimaOrder,maximaOrder

    real(kind=8)::vol,totvol,hstep,height,maxheight,dh
    real(kind=8),dimension(:),allocatable::fz
    ! Order minima from top to bottom
    if(allocated(maximaOrder)) deallocate(maximaOrder)
    if(allocated(minimaOrder)) deallocate(minimaOrder)
    if(allocated(fz)) deallocate(fz)
    allocate(maximaOrder(baseNb))
    allocate(minimaOrder(baseNb))
    allocate(fz(baseNb))

    maximaOrder=-1
    do q=1,baseNb
      c=baselist(q)
      minimaOrder(q)=q
      if(catchmentID(c)<0)then
        fz(q)=spmZ(-catchmentID(c))
      else
        fz(q)=-1.e5-real(q)
      endif
    enddo
    call quick_sort(fz,minimaOrder)
    p=0
    do q=baseNb,1,-1
      p=p+1
      maximaOrder(p)=minimaOrder(q)
    enddo
    do qo=1,baseNb
      ! Get the highest isolated catchment
      q=maximaOrder(qo)
      c=baselist(q)
      if(catchmentID(c)<0)then
        ! Define catchment parameters
        cID=catchmentID(c)
        vol=0.
        hstep=0.1
        height=spmZ(-cID)
        maxheight=sillMax(q)-spmZ(-cID)
        totvol=discharge(-cID)*(1.0-infiltration_evaporation) 

        ! Sort catchment points by increasing elevations
        if(allocated(ptOrder)) deallocate(ptOrder)
        if(allocated(fz)) deallocate(fz)
        allocate(ptOrder(sillNb(q)))
        allocate(fz(sillNb(q)))
        ptOrder=-1
        if(sillNb(q)>1)then
          do n=1,sillNb(q)
            p=sillGID(q,n)
            fz(n)=spmZ(p)
            ptOrder(n)=n
          enddo
          call quick_sort(fz,ptOrder)
        else
          ptOrder(1)=1
        endif
        ! Loop through the nodes to distribute volume of water until:
        ! 1- the volume has been totally distributed
        ! 2- the water flows to another catchment
        do while(vol<totvol.and.height<maxheight) 
          height=height+hstep
          watercell(-cID)=watercell(-cID)+hstep
          vol=vol+hstep*voronoiCell(-cID)%area

          do no=2,sillNb(q)
            n=ptOrder(no)
            p=sillGID(q,n)
            ! Add water to cell
            if(spmZ(p)+watercell(p)<=height)then
              dh=height-spmZ(p)-watercell(p) 
              watercell(p)=watercell(p)+dh
              vol=vol+dh*voronoiCell(p)%area
            endif

            if(watercell(p)>0.)then
              do k=1,delaunayVertex(p)%ngbNb
                nid=delaunayVertex(p)%ngbID(k)
                ! Neighbor nodes in a different catchment
                if(nid>0)then
                  if(spmZ(nid)<spmZ(p)+watercell(p).and. &
                    catchmentID(nid)/=cID) goto 18
                endif
              enddo
            endif
          enddo
        enddo
18 continue      
      endif
    enddo

  end subroutine define_drained_water_thickness
  ! =====================================================================================

  subroutine planchon_dem_fill_algorithm

    logical::flag
    integer::p,k,cID 

    real(kind=8)::step 

    ! Update DEM borders elevation values
    step=0.0001
    filldem=1.e6_8
    do k=1,dnodes
      cID=catchmentID(k)
      if(cID<0)then
        if(watercell(-cID)<fh.and.Cerodibility>0.)then
          watercell(k)=spmZ(-cID)+fh-spmZ(k)
        endif
      endif
      
      if(watercell(k)<0.) watercell(k)=0.
      
      ! On border fix DEM elevation for filling calculation
      if(tcoordX(k)==minx-dx.or.tcoordX(k)==maxx+dx.or.tcoordY(k)==miny-dx &
        .or.tcoordY(k)==maxy+dx)then 
        filldem(k)=spmZ(k)+watercell(k)
      endif
    enddo

    ! Now find the sinks and fill them using Planchon's method
    flag=.true.
    do while(flag)
      flag=.false.
      do k=1,dnodes
        if(tcoordX(k)>minx-dx.and.tcoordX(k)<maxx+dx.and.tcoordY(k)>miny-dx &
          .and.tcoordY(k)<maxy+dx.and.filldem(k)>spmZ(k))then
          ! Loop through the neighbors
          do p=1,delaunayVertex(k)%ngbNb
            if(delaunayVertex(k)%ngbID(p)>0)then
              ! In case delaunay elevation greater than neighbor's ones
              if(spmZ(k)>=filldem(delaunayVertex(k)%ngbID(p))+step)then
                filldem(k)=spmZ(k)
              ! Otherwise this is a sink and we perform sink filling
              else
                if(filldem(k)>filldem(delaunayVertex(k)%ngbID(p))+step)then
                  filldem(k)=filldem(delaunayVertex(k)%ngbID(p))+step
                  if(filldem(k)-spmZ(k)>fh.and.Cerodibility>0.)then 
                    filldem(k)=spmZ(k)+fh
                  else
                    flag=.true.
                  endif
                endif
              endif
            endif
          enddo
        endif
      enddo
    enddo

    if(Cerodibility==0.)then
      watercell=filldem-spmZ
    endif

    return

  end subroutine planchon_dem_fill_algorithm
  ! =====================================================================================

  subroutine compute_vertical_displacement

    integer::k,lid,id

    do lid=1,localNodes
      id=localNodesGID(lid)
      k=stackOrder(id) 
      nZ(k)=nZ(k)+tvertDisp(k)*time_step
    enddo

    do k=1,upartN
      id=unodeID(k)
      cumDisp(k)=cumDisp(k)+tvertDisp(id)*time_step
    enddo

  end subroutine compute_vertical_displacement
  ! =====================================================================================

  subroutine update_grid_borders

    integer::k,p

    if(simulation_time==time_start.or.update3d)then
      if(.not.allocated(delemID)) allocate(delemID(delem))
      delemoo=0
      delemID=-1
      p=0
      do k=1,delem
         if(elemtmask(k)==0.and.uownEID(k)==pet_id)then 
           delemoo=delemoo+1
           p=p+1
           delemID(p)=k
         endif
      enddo
    endif

    ! Define the borders according to boundary definition
    do k=1,dnodes
      ! On border fix DEM elevation
      if(voronoiCell(k)%border==1)then 

        spmH(k)=0.0_8

        ! In case there is an outlet
        if(outlet/=0)then
          if(voronoiCell(k)%btype==outlet)then
            spmZ(k)=spmZ(voronoiCell(k)%bpoint)-1.0_8
          else
             spmZ(k)=1.e6_8 
          endif
        else
          !--------- WEST
          ! South West corner
          if(voronoiCell(k)%btype==1)then
            ! Fixed
            if(bounds(3)==0)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)
            ! Wall
            elseif(bounds(3)==1)then
              spmZ(k)=1.e6_8 
            ! Fall
            elseif(bounds(3)==2)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)-1.0_8
            endif
          endif
          ! North West corner
          if(voronoiCell(k)%btype==3)then
            ! Fixed
            if(bounds(3)==0)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)
            ! Wall
            elseif(bounds(3)==1)then
              spmZ(k)=1.e6_8
            ! Fall
            elseif(bounds(3)==2)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)-1.0_8
            endif
          endif
          ! West border
          if(voronoiCell(k)%btype==5)then
            ! Fixed
            if(bounds(3)==0)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)
            ! Wall
            elseif(bounds(3)==1)then
              spmZ(k)=1.e6_8 
            ! Fall
            elseif(bounds(3)==2)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)-1.0_8
            endif
          endif
          !--------- EAST
          ! South East corner
          if(voronoiCell(k)%btype==2)then
            ! Fixed
            if(bounds(4)==0)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)
            ! Wall
            elseif(bounds(4)==1)then
              spmZ(k)=1.e6_8 
            ! Fall
            elseif(bounds(4)==2)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)-1.0_8
            endif
          endif
          ! North East corner
          if(voronoiCell(k)%btype==4)then
            ! Fixed
            if(bounds(4)==0)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)
            ! Wall
            elseif(bounds(4)==1)then
              spmZ(k)=1.e6_8
            ! Fall
            elseif(bounds(4)==2)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)-1.0_8
            endif
          endif
          ! East border
          if(voronoiCell(k)%btype==6)then
            ! Fixed
            if(bounds(4)==0)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)
            ! Wall
            elseif(bounds(4)==1)then
              spmZ(k)=1.e6_8 
            ! Fall
            elseif(bounds(4)==2)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)-1.0_8
            endif
          endif
          !--------- SOUTH
          ! South East corner
          if(voronoiCell(k)%btype==2)then
            ! Fixed
            if(bounds(2)==0)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)
            ! Wall
            elseif(bounds(2)==1)then
              spmZ(k)=1.e6_8 
            ! Fall
            elseif(bounds(2)==2)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)-1.0_8
            endif
          endif
          ! South West corner
          if(voronoiCell(k)%btype==1)then
            ! Fixed
            if(bounds(2)==0)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)
            ! Wall
            elseif(bounds(2)==1)then
              spmZ(k)=1.e6_8 
            ! Fall
            elseif(bounds(2)==2)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)-1.0_8
            endif
          endif
          ! South border
          if(voronoiCell(k)%btype==7)then
            ! Fixed
            if(bounds(2)==0)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)
            ! Wall
            elseif(bounds(2)==1)then
              spmZ(k)=1.e6_8
            ! Fall
            elseif(bounds(2)==2)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)-1.0_8
            endif
          endif
          !--------- NORTH
          ! North East corner
          if(voronoiCell(k)%btype==4)then
            ! Fixed
            if(bounds(1)==0)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)
            ! Wall
            elseif(bounds(1)==1)then
              spmZ(k)=1.e6_8 
            ! Fall
            elseif(bounds(1)==2)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)-1.0_8
            endif
          endif
          ! North West corner
          if(voronoiCell(k)%btype==3)then
            ! Fixed
            if(bounds(1)==0)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)
            ! Wall
            elseif(bounds(1)==1)then
              spmZ(k)=1.e6_8 
            ! Fall
            elseif(bounds(1)==2)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)-1.0_8
            endif
          endif
          ! North border
          if(voronoiCell(k)%btype==8)then
            ! Fixed
            if(bounds(1)==0)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)
            ! Wall
            elseif(bounds(1)==1)then
              spmZ(k)=1.e6_8 
            ! Fall
            elseif(bounds(1)==2)then
              spmZ(k)=spmZ(voronoiCell(k)%bpoint)-1.0_8
            endif
          endif
        endif

      endif
    enddo

    return

  end subroutine update_grid_borders
  ! =====================================================================================

  subroutine DeriveTrianglePlanes(x,y,id1,id2,id3,z)

    integer::id1,id2,id3
    real(kind=8)::s,x,y,z

    real(kind=8),dimension(3)::d1,d2,n
    real(kind=8),dimension(4)::plane

    d1(1)=tcoordX(id2)-tcoordX(id1)
    d1(2)=tcoordY(id2)-tcoordY(id1)
    d1(3)=tcoordZ(id2)-tcoordZ(id1)

    d2(1)=tcoordX(id3)-tcoordX(id1)
    d2(2)=tcoordY(id3)-tcoordY(id1)
    d2(3)=tcoordZ(id3)-tcoordZ(id1)

    n(1)=d2(2)*d1(3)-d2(3)*d1(2)
    n(2)=d2(3)*d1(1)-d2(1)*d1(3)
    n(3)=d2(1)*d1(2)-d2(2)*d1(1)

    s=1.0/sqrt(n(1)*n(1)+n(2)*n(2)+n(3)*n(3))

    plane(1)=n(1)*s
    plane(2)=n(2)*s
    plane(3)=n(3)*s
    plane(4)= -(plane(1)*tcoordX(id1)+plane(2)*tcoordY(id1)+plane(3)*tcoordZ(id1))

    z=-(plane(1)*x+plane(2)*y+plane(4))/plane(3)

  end subroutine DeriveTrianglePlanes
  ! =====================================================================================

  subroutine DeriveTrianglePlanes2(xy,xa,ya,za,z)

    real(kind=8),dimension(2)::xy
    real(kind=8),dimension(3)::xa,ya,za
    real(kind=8)::s,z

    real(kind=8),dimension(3)::d1,d2,n
    real(kind=8),dimension(4)::plane

    d1(1)=xa(2)-xa(1)
    d1(2)=ya(2)-ya(1)
    d1(3)=za(2)-za(1)

    d2(1)=xa(3)-xa(1)
    d2(2)=ya(3)-ya(1)
    d2(3)=za(3)-za(1)

    n(1)=d2(2)*d1(3)-d2(3)*d1(2)
    n(2)=d2(3)*d1(1)-d2(1)*d1(3)
    n(3)=d2(1)*d1(2)-d2(2)*d1(1)

    s=1.0/sqrt(n(1)*n(1)+n(2)*n(2)+n(3)*n(3))

    plane(1)=n(1)*s
    plane(2)=n(2)*s
    plane(3)=n(3)*s
    plane(4)= -(plane(1)*xa(1)+plane(2)*ya(1)+plane(3)*za(1))

    z=-(plane(1)*xy(1)+plane(2)*xy(2)+plane(4))/plane(3)

  end subroutine DeriveTrianglePlanes2
  ! =====================================================================================
  
  subroutine inside_triangle(xa,ya,xb,yb,l)

    integer,intent(out)::l

    real(kind=8),intent(in)::xa,ya,xb(3),yb(3)
    real(kind=8)::det0,det1,det2

    l=-1

    det0=(xb(2)-xb(1))*(ya-yb(1))-(yb(2)-yb(1))*(xa-xb(1))
    det1=(xb(3)-xb(2))*(ya-yb(2))-(yb(3)-yb(2))*(xa-xb(2))
    det2=(xb(1)-xb(3))*(ya-yb(3))-(yb(1)-yb(3))*(xa-xb(3))

    if(det0>=0.and.det1>=0.and.det2>=0)then
      l=1
    elseif(det0<=0.and.det1<=0.and.det2<=0)then
      l=1
    endif

    return

  end subroutine inside_triangle
  ! =====================================================================================

end module hydrology 