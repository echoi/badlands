! =====================================================================================
! BADLANDS (BAsin anD LANdscape DynamicS)
!
! Copyright (c) Tristan Salles (The University of Sydney)
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by the Free Software
! Foundation; either version 3.0 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
! more details.
!
! You should have received a copy of the GNU Lesser General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 59 Temple
! Place, Suite 330, Boston, MA 02111-1307 USA
! =====================================================================================

! =====================================================================================
!
!       Filename:  IsostaticFlexure.f90
!
!    Description:  Flexural isostasy computation (Li et al. 2004 CompGeo) and compaction
!
!        Version:  1.0
!        Created:  11/07/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles
!
! =====================================================================================


module isoflex

  use parallel
  use topology
  use parameters
  use hydroUtil
  use external_forces

  implicit none

  integer::n2x,n2y,n4x,n4y

  ! Fourth order schema parameters
  real(kind=8),parameter::q1=216.,q2=-18.,q3=-84.
  real(kind=8),parameter::q4=-3.,q5=792.
  real(kind=8),parameter::r1=144.,r2=18.,r3=-48.
  real(kind=8),parameter::r4=-6.,r5=1./240.

  real(kind=8),dimension(:,:),allocatable::w,w1,wx,wy,ld

contains

  ! =====================================================================================
  subroutine porosity_compaction

    integer::k,p

    real(kind=8)::phi,subs,pressure_lithos,mass

    do k=1,dnodes
      pressure_lithos=0.
      if(spmZ(k)<gsea%actual_sea) pressure_lithos=9.81*sea_water_density*(gsea%actual_sea-spmZ(k))
      subs=0.
      do p=flex_lay-1,2,-1
        ! Mass due to water in pore space
        mass=mean_sediment_density*ulay_th(k,p+1)*(1-ulay_phi(k,p+1))
        mass=mass+sea_water_density*ulay_th(k,p+1)*ulay_phi(k,p+1)
        ! There is some sediment in this layer
        if(mass>0.)then
          pressure_lithos=pressure_lithos+mass*9.81
          ! Calculate new porosity with this pressure
          call get_porosity(pressure_lithos,ulay_phi(k,p),phi)
          ! Subsidence due to porosity change
          subs=subs+ulay_th(k,p)*(ulay_phi(k,p)-phi)
          ! Update layer thickness
          ulay_th(k,p)=ulay_th(k,p)-ulay_th(k,p)*(ulay_phi(k,p)-phi)
          ! Update porosity
          ulay_phi(k,p)=phi
        endif
      enddo
      ! Correct the topographic elevation due to compactional subsidence
      spmZ(k)=spmZ(k)-subs
      sedthick(k)=sedthick(k)-subs
    enddo

  end subroutine porosity_compaction
  ! =====================================================================================

  subroutine get_porosity(Plith,in_phi,out_phi)

    integer:: p,fd
    real(kind=8)::in_phi,out_phi,Plith

    ! Go through the pressure field and find the fitting interval
    if(Plith<pressTable(1))then
      out_phi=poroTable(1)
    endif
    fd=0
    loop: do p=2,pressureFields
      if(Plith>=pressTable(p-1).and.Plith<pressTable(p))then
        fd=p
        out_phi=(Plith-pressTable(p-1))/(pressTable(p)-pressTable(p-1))
        out_phi=out_phi*(poroTable(p)-poroTable(p-1))+poroTable(p-1)
        exit loop
      endif
    enddo loop
    if(fd==0) out_phi=poroTable(pressureFields)

    out_phi=min(out_phi,in_phi)

  end subroutine get_porosity
  ! =====================================================================================

  subroutine update_Flex_array

    integer::step,i,j,ic,jc,p,m
    real(kind=8),dimension(nbnodes)::tempZ,tempSh

    step=int(flex_dx/dx)

    ! Read changes in topographic regular grid.
    p=0
    m=0
    do j=1,ny+2
      do i=1,nx+2
        p=p+1
        if(i>1.and.i<nx+2.and.j>1.and.j<ny+2)then
          m=m+1
          tempZ(m)=rcoordZ(p)
          tempSh(m)=rsedload(p)
        endif
      enddo
    enddo

    j=1
    do jc=2,nbfy+1
      i=1
      p=(j-1)*nx+i
      do ic=2,nbfx+1
        flexZ(ic,jc)=tempZ(p)
        flexSed(ic,jc)=tempSh(p)
        i=i+step
        p=(j-1)*nx+i
      enddo
      j=j+step
    enddo

    ! Update border
    flexZ(2:nbfx+1,1)=flexZ(2:nbfx+1,2)
    flexZ(2:nbfx+1,nbfy+2)=flexZ(2:nbfx+1,nbfy+1)
    flexZ(1,2:nbfy+1)=flexZ(2,2:nbfy+1)
    flexZ(nbfx+2,2:nbfy+1)=flexZ(nbfx+1,2:nbfy+1)

    flexSed(2:nbfx+1,1)=flexSed(2:nbfx+1,2)
    flexSed(2:nbfx+1,nbfy+2)=flexSed(2:nbfx+1,nbfy+1)
    flexSed(1,2:nbfy+1)=flexSed(2,2:nbfy+1)
    flexSed(nbfx+2,2:nbfy+1)=flexSed(nbfx+1,2:nbfy+1)

    ! Update corner
    flexZ(1,1)=flexZ(2,2)
    flexZ(1,nbfy+2)=flexZ(2,nbfy+1)
    flexZ(nbfx+2,1)=flexZ(nbfx+1,2)
    flexZ(nbfx+2,nbfy+2)=flexZ(nbfx+1,nbfy+1)

    flexSed(1,1)=flexSed(2,2)
    flexSed(1,nbfy+2)=flexSed(2,nbfy+1)
    flexSed(nbfx+2,1)=flexSed(nbfx+1,2)
    flexSed(nbfx+2,nbfy+2)=flexSed(nbfx+1,nbfy+1)

  end subroutine update_Flex_array
  ! =====================================================================================

  subroutine built_initial_load

    integer::step,i,j,ic,jc,p,m

    real(kind=8)::tmp
    real(kind=8),dimension(nbnodes)::tempZ

    step=int(flex_dx/dx)

    ! Read changes in topographic regular grid.
    p=0
    m=0
    do j=1,ny+2
      do i=1,nx+2
        p=p+1
        if(i>1.and.i<nx+2.and.j>1.and.j<ny+2)then
          m=m+1
          tempZ(m)=rcoordZ(p)
        endif
      enddo
    enddo

    j=1
    do jc=2,nbfy+1
      i=1
      p=(j-1)*nx+i
      do ic=2,nbfx+1
        flexZ(ic,jc)=tempZ(p)
        i=i+step
        p=(j-1)*nx+i
      enddo
      j=j+step
    enddo

    ! Update border
    flexZ(2:nbfx+1,1)=flexZ(2:nbfx+1,2)
    flexZ(2:nbfx+1,nbfy+2)=flexZ(2:nbfx+1,nbfy+1)
    flexZ(1,2:nbfy+1)=flexZ(2,2:nbfy+1)
    flexZ(nbfx+2,2:nbfy+1)=flexZ(nbfx+1,2:nbfy+1)

    ! Update corner
    flexZ(1,1)=flexZ(2,2)
    flexZ(1,nbfy+2)=flexZ(2,nbfy+1)
    flexZ(nbfx+2,1)=flexZ(nbfx+1,2)
    flexZ(nbfx+2,nbfy+2)=flexZ(nbfx+1,nbfy+1)

    if(allocated(sedloader))then
      p=0
      do j=1,nbfy+2
        do i=1,nbfx+2
          p=p+1
          prevload(i,j)=sedloader(p)
        enddo
      enddo
    else
      do j=1,nbfy+2
        do i=1,nbfx+2
          flexSed(i,j)=mean_sediment_density*100000.0*(1-poroTable(pressureFields))+100000.0*poroTable(pressureFields)*sea_water_density
          tmp=flexSed(i,j) !mean_sediment_density*flexSed(i,j) !*(1-flexPor(i,j))+flexPor(i,j)*flexSed(i,j)*sea_water_density
          if(flexZ(i,j)<gsea%actual_sea) tmp=tmp+(gsea%actual_sea-flexZ(i,j))*sea_water_density
          prevload(i,j)=tmp*cst1
        enddo
      enddo
    endif



  end subroutine built_initial_load
  ! =====================================================================================

  subroutine isostatic_flexure

    integer::i,j,p,m
    real(kind=8)::tmp,oldload,dtot,wdiff,wtot

    load=0.

    dtot=0.0
    do j=1,nbfy+2
      do i=1,nbfx+2
        if(j>1.and.j<nbfy+2.and.i>1.and.i<nbfx+2)then
          tmp=flexSed(i,j) !*(1-flexPor(i,j))+flexPor(i,j)*flexSed(i,j)*sea_water_density
          if(flexZ(i,j)<gsea%actual_sea) tmp=tmp+(gsea%actual_sea-flexZ(i,j))*sea_water_density
          oldload=prevload(i,j)
          prevload(i,j)=cst1*tmp
          load(i,j)=prevload(i,j)-oldload
          dtot=dtot+load(i,j)
        endif
      enddo
      if(j>1.and.j<nbfy+2)then
        prevload(1,j)=prevload(2,j)
        prevload(nbfx+2,j)=prevload(nbfx+1,j)
        load(1,j)=load(2,j)
        load(nbfx+2,j)=load(nbfx+1,j)
      endif
    enddo

    ! Tiny added load, ignore isostatic changes
    if(abs(dtot)<cst1*mean_sediment_density)return

    prevload(1:nbfx+2,1)=prevload(1:nbfx+2,2)
    prevload(1:nbfx+2,nbfy+2)=prevload(1:nbfx+2,nbfy+1)
    load(1:nbfx+2,1)=load(1:nbfx+2,2)
    load(1:nbfx+2,nbfy+2)=load(1:nbfx+2,nbfy+1)

    ! Initialise flexure arrays
    if(.not.allocated(w)) allocate(w(nbfx+2,nbfy+2))
    if(.not.allocated(wx)) allocate(wx(nbfx+2,nbfy+2))
    if(.not.allocated(wy)) allocate(wy(nbfx+2,nbfy+2))
    if(.not.allocated(w1)) allocate(w1(nbfx+2,nbfy+2))
    if(.not.allocated(ld)) allocate(ld(nbfx+2,nbfy+2))

    ! Update flexure parameters
    w=0.0
    w1=0.0
    wx=0.0
    wy=0.0
    ld=load

    n2x=int(nbfx/2)
    n2x=n2x*2
    n4x=int(nbfx/4)
    n4x=n4x*4
    n2y=int(nbfy/2)
    n2y=n2y*2
    n4y=int(nbfy/4)
    n4y=n4y*4

    p=0
    ! Iterate until a solution is reached
    do
      p=p+1
      ! Fine grid first
      m=1
      call solve_flexure(m,nbfx+1,nbfy+1,4)
      wdiff=0.
      wtot=0.
      ! Compare with the last finest grid results
      do j=2,nbfy+2-m
        do i=2,nbfx+2-m
          wdiff=wdiff+abs(w1(i,j)-w(i,j))
          wtot=wtot+abs(w(i,j))
        enddo
      enddo

      w1=w
      if(wdiff<torb*wtot.or.p>10000)exit
      if(wtot==0.0)exit

      ! Coarser grid 2dx, full weighting operator
      m=2
      do j=2+m,n2y+2-m,m
        do i=2+m,n2x+2-m,m
          w(i,j)=0.0625_8*(w(i-1,j-1)+w(i+1,j-1)+w(i-1,j+1)+w(i+1,j+1)+ &
            2.0_8*(w(i,j-1)+w(i,j+1)+w(i-1,j)+w(i+1,j))+4.0_8*w(i,j))
          wx(i,j)=0.0625_8*(wx(i-1,j-1)+wx(i+1,j-1)+wx(i-1,j+1)+wx(i+1,j+1)+ &
            2.0_8*(wx(i,j-1)+wx(i,j+1)+wx(i-1,j)+wx(i+1,j))+4.0_8*wx(i,j))
          wy(i,j)=0.0625_8*(wy(i-1,j-1)+wy(i+1,j-1)+wy(i-1,j+1)+wy(i+1,j+1)+ &
            2.0_8*(wy(i,j-1)+wy(i,j+1)+wy(i-1,j)+wy(i+1,j))+4.0_8*wy(i,j))
        enddo
      enddo
      call solve_flexure(m,n2x+2,n2y+2,4)

      ! Coarser grid 4dx, full weighting operator
      m=4
      do j=2+m,n4y+2-m,m
        do i=2+m,n4x+2-m,m
          w(i,j)=0.0625_8*(w(i-2,j-2)+w(i+2,j-2)+w(i-2,j+2)+w(i+2,j+2)+ &
            2.0_8*(w(i,j-2)+w(i,j+2)+w(i-2,j)+w(i+2,j))+4.0_8*w(i,j))
          wx(i,j)=0.0625_8*(wx(i-2,j-2)+wx(i+2,j-2)+wx(i-2,j+2)+wx(i+2,j+2)+ &
            2.0_8*(wx(i,j-2)+wx(i,j+2)+wx(i-2,j)+wx(i+2,j))+4.0_8*wx(i,j))
          wy(i,j)=0.0625_8*(wy(i-2,j-2)+wy(i+2,j-2)+wy(i-2,j+2)+wy(i+2,j+2)+ &
            2.0_8*(wy(i,j-2)+wy(i,j+2)+wy(i-2,j)+w(i+2,j))+4.0_8*wy(i,j))
        enddo
      enddo
      call boundary_flexure(m,n4x+2-m,n4y+2-m,w)
      call boundary_flexure(m,n4x+2-m,n4y+2-m,wx)
      call boundary_flexure(m,n4x+2-m,n4y+2-m,wy)
      call solve_flexure(m,n4x+2,n4y+2,16)

      ! Interpolate to finer grid (2dx)
      do j=2,n4y+2-m,4
        do i=2,n4x+2-m,4
          w(i,j+2)=0.5_8*(w(i,j)+w(i,j+m))
          w(i+2,j)=0.5_8*(w(i,j)+w(i+m,j))
          w(i+2,j+2)=0.25_8*(w(i,j)+w(i,j+m)+w(i+m,j)+w(i+m,j+m))
          wx(i,j+2)=0.5_8*(wx(i,j)+wx(i,j+m))
          wx(i+2,j)=0.5_8*(wx(i,j)+wx(i+m,j))
          wx(i+2,j+2)=0.25_8*(wx(i,j)+wx(i,j+m)+wx(i+m,j)+wx(i+m,j+m))
          wy(i,j+2)=0.5_8*(wy(i,j)+wy(i,j+m))
          wy(i+2,j)=0.5_8*(wy(i,j)+wy(i+m,j))
          wy(i+2,j+2)=0.25_8*(wy(i,j)+wy(i,j+m)+wy(i+m,j)+wy(i+m,j+m))
        enddo
      enddo

      ! At boundaries of finer grid (2dx)
      if(n2y>n4y)then
        do i=2,n4x+2-m,4
          w(i,n2y+2)=w(i,n4y+2)
          w(i+2,n2y+2)=0.5*(w(i,n4y+2)+w(i+m,n4y+2))
          wx(i,n2y+2)=wx(i,n4y+2)
          wx(i+2,n2y+2)=0.5*(wx(i,n4y+2)+wx(i+m,n4y+2))
          wy(i,n2y+2)=wy(i,n4y+2)
          wy(i+2,n2y+2)=0.5*(wy(i,n4y+2)+wy(i+m,n4y+2))
        enddo
      endif

      if(n2x>n4x)then
        do j=2,n4y+2-m,4
          w(n2x+2,j)=w(n4x+2,j)
          w(n2x+2,j+2)=0.5*(w(n4x+2,j)+w(n4x+2,j+m))
          wx(n2x+2,j)=wx(n4x+2,j)
          wx(n2x+2,j+2)=0.5*(wx(n4x+2,j)+wx(n4x+2,j+m))
          wy(n2x+2,j)=wy(n4x+2,j)
          wy(n2x+2,j+2)=0.5*(wy(n4x+2,j)+wy(n4x+2,j+m))
        enddo
      endif
      w(n2x+2,n2y+2)=w(n4x+2,n4y+2)
      wx(n2x+2,n2y+2)=wx(n4x+2,n4y+2)
      wy(n2x+2,n2y+2)=wy(n4x+2,n4y+2)
      m=2
      call solve_flexure(m,n2x+2,n2y+2,4)

      ! Finer grid dx
      do j=2,n2y+2-m,2
        do i=2,n2x+2-m,2
          w(i,j+1)=0.5*(w(i,j)+w(i,j+m))
          w(i+1,j)=0.5*(w(i,j)+w(i+m,j))
          w(i+1,j+1)=0.25*(w(i,j)+w(i,j+m)+w(i+m,j)+w(i+m,j+m))
          wx(i,j+1)=0.5*(wx(i,j)+wx(i,j+m))
          wx(i+1,j)=0.5*(wx(i,j)+wx(i+m,j))
          wx(i+1,j+1)=0.25*(wx(i,j)+wx(i,j+m)+wx(i+m,j)+wx(i+m,j+m))
          wy(i,j+1)=0.5*(wy(i,j)+wy(i,j+m))
          wy(i+1,j)=0.5*(wy(i,j)+wy(i+m,j))
          wy(i+1,j+1)=0.25*(wy(i,j)+wy(i,j+m)+wy(i+m,j)+wy(i+m,j+m))
        enddo
      enddo

      if(nbfy+1>n2y+2)then
        do i=2,n2x+2-m,2
          w(i,nbfy+1)=w(i,n2y+2)
          w(i+1,nbfy+1) =0.5*(w(i,n2y+2)+w(i+m,n2y+2))
          wx(i,nbfy+1) =wx(i,n2y+2)
          wx(i+1,nbfy+1)=0.5*(wx(i,n2y+2)+wx(i+m,n2y+2))
          wy(i,nbfy+1) =wy(i,n2y+2)
          wy(i+1,nbfy+1)=0.5*(wy(i,n2y+2)+wy(i+m,n2y+2))
        enddo
      endif
      if(nbfx+1>n2x+2)then
        do j=2,n2y+2-m,2
          w(nbfx+1,j)=w(n2x+2,j)
          w(nbfx+1,j+1)=0.5*(w(n2x+2,j)+w(n2x+2,j+m))
          wx(nbfx+1,j)=wx(n2x+2,j)
          wx(nbfx+1,j+1)=0.5*(wx(n2x+2,j)+wx(n2x+2,j+m))
          wy(nbfx+1,j)=wy(n2x+2,j)
          wy(nbfx+1,j+1)=0.5*(wy(n2x+2,j)+wy(n2x+2,j+m))
        enddo
      endif
      w(nbfx+1,nbfy+1)=w(n2x+2,n2y+2)
      wx(nbfx+1,nbfy+1)=wx(n2x+2,n2y+2)
      wy(nbfx+1,nbfy+1)=wy(n2x+2,n2y+2)
      do i=1,nbfx+2
        w(i,1)=w(i,2)
        w(i,nbfy+2)=w(i,nbfy+1)
        wx(i,1)=wx(i,2)
        wx(i,nbfy+2)=wx(i,nbfy+1)
        wy(i,1)=wy(i,2)
        wy(i,nbfy+2)=wy(i,nbfy+1)
      enddo
      do j=1,nbfy+2
        w(1,j)=w(2,j)
        w(nbfx+2,j)=w(nbfx+1,j)
        wx(1,j)=wx(2,j)
        wx(nbfx+2,j)=wx(nbfx+1,j)
        wy(1,j)=wy(2,j)
        wy(nbfx+2,j)=wy(nbfx+1,j)
      enddo
    enddo

    if(p>=10000)then
      print*,'Issue: isostacy did not converge to a solution'
      print*,wdiff,1.0E-18_8*wtot,wtot
    endif
    flexDisp=w

  end subroutine isostatic_flexure
  ! =====================================================================================

  subroutine boundary_flexure(ks,ncl,nrw,temp)

    integer::ks,ncl,nrw,i,j,k1
    real(kind=8),dimension(nbfx+2,nbfy+2)::temp

    k1=2+ks
    do j=1,nbfy+2
      do i=1,k1-1
        temp(i,j)=temp(k1,j)
      enddo
      do i=ncl+1,nbfx+2
        temp(i,j)=temp(ncl,j)
      enddo
    enddo

    do i=1,nbfx+2
      do j=1,k1-1
        temp(i,j)=temp(i,k1)
      enddo
      do j=nrw+1,nbfy+2
        temp(i,j)=temp(i,nrw)
      enddo
    enddo

  end subroutine boundary_flexure
  ! =====================================================================================

  subroutine solve_flexure(m,ncl,nrw,nloop)

    integer::m,ncl,nrw,nloop,i,j,isw,jsw,ks,ipass,n

    real(kind=8)::resw,w3,wtot,wdiff,piv,dxm4,el1

    ks=0
    if(m>1)ks=m
    dxm4=(flex_dx**4)*m**4

    do j=2,nrw,m
      do i=2,ncl,m
        el1=flexZ(i,j)
        if(el1<=gsea%actual_sea)then
          ld(i,j)=(load(i,j)-cst2*w(i,j))*dxm4
        else
          ld(i,j)=(load(i,j)-cst3*w(i,j))*dxm4
        endif
      enddo
    enddo
    call boundary_flexure(0,ncl,nrw,ld)

    n=0
    do
      n=n+1
      wtot=0.0
      wdiff=0.0
      jsw=1
      do ipass=1,2
        isw=jsw
        do j=2+ks,nrw-ks,m
          do i=isw+1+ks,ncl-ks,2*m
            el1=flexZ(i,j)
            if(el1<=gsea%actual_sea)then
              piv=q5+11.0*cst2*dxm4
            else
              piv=q5+11.0*cst3*dxm4
            endif

            wx(i,j)=(r1*(w(i+m,j)-w(i-m,j)) &
              +r2*(w(i+m,j+m)-w(i-m,j+m)+w(i+m,j-m)-w(i-m,j-m))  &
              +r3*(wx(i+m,j)+wx(i-m,j)) &
              +r4*(wx(i+m,j+m)+wx(i-m,j+m)+wx(i+m,j-m)+wx(i-m,j-m)) &
              +r4*(wy(i+m,j+m)-wy(i-m,j+m)-wy(i+m,j-m)+wy(i-m,j-m)) &
              +(ld(i+m,j)-ld(i-m,j)))*r5

            wy(i,j)=(r1*(w(i,j+m)-w(i,j-m)) &
              +r2*(w(i+m,j+m)+w(i-m,j+m)-w(i+m,j-m)-w(i-m,j-m))  &
              +r3*(wy(i,j+m)+wy(i,j-m)) &
              +r4*(wy(i+m,j+m)+wy(i-m,j+m)+wy(i+m,j-m)+wy(i-m,j-m)) &
              +r4*(wx(i+m,j+m)-wx(i-m,j+m)-wx(i+m,j-m)+wx(i-m,j-m)) &
              +(ld(i,j+m)-ld(i,j-m)))*r5

            w3=(q1*(w(i+m,j)+w(i-m,j)+w(i,j+m)+w(i,j-m)) &
              +q2*(w(i+m,j+m)+w(i-m,j+m)+w(i+m,j-m)+w(i-m,j-m)) &
              +q3*(wx(i+m,j)-wx(i-m,j)+wy(i,j+m)-wy(i,j-m)) &
              +q4*(wx(i+m,j+m)-wx(i-m,j+m)+wx(i+m,j-m)-wx(i-m,j-m)) &
              +q4*(wy(i+m,j+m)+wy(i-m,j+m)-wy(i+m,j-m)-wy(i-m,j-m)) &
              +(ld(i+m,j)+ld(i-m,j)+ld(i,j+m)+ld(i,j-m)) &
              +11.0*load(i,j)*dxm4)/piv

            if(el1<=gsea%actual_sea)then
              ld(i,j)=(load(i,j)-cst2*w3)*dxm4
            else
              ld(i,j)=(load(i,j)-cst3*w3)*dxm4
            endif

            resw=w3-w(i,j)
            wdiff=wdiff+abs(w3-w(i,j))
            wtot=wtot+abs(w3)
            w(i,j)=w3
          enddo
          isw=m+2-isw
        enddo
        jsw=m+2-jsw
        call boundary_flexure(ks,ncl-ks,nrw-ks,w)
        call boundary_flexure(ks,ncl-ks,nrw-ks,wx)
        call boundary_flexure(ks,ncl-ks,nrw-ks,wy)
        call boundary_flexure(ks,ncl-ks,nrw-ks,ld)
      enddo

      if(wdiff<1.0E-10*wtot.or.n>nloop)exit
      if(wtot==0.0)exit
    enddo

  end subroutine solve_flexure
  ! =====================================================================================

end module isoflex
! =====================================================================================
