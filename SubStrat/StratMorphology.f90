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
!       Filename:  StratMorphology.f90
!
!    Description:  Defines morphological formalism used in BADLANDS for stratigraphic grid.
!
!        Version:  1.0
!        Created:  11/06/15 09:11:25
!        Revision:  none
!
!        Author:  Tristan Salles     
!
! =====================================================================================

module stratmorph

  use bilinear
  use parallel
  use topology
  use parameters
  use hydrology
  use hydroUtil
  use watershed
  use out_stratal
  use strata_evol
  use stratal_class
  use strat_hillslope
  use outspm_surface
  use outspm_drainage
  use external_forces

  implicit none
  
  integer::iter,idiff

  real(kind=8)::cpl_time,max_time
  real(kind=8),dimension(:),allocatable::nth
  real(kind=8),dimension(:,:),allocatable::nsed
!   real(kind=8)::time1,time2,tt1,tt2

contains

  ! =====================================================================================

  subroutine stratgeomorph

    integer::k,ks

    real(kind=8),dimension(dnodes)::tsed,gsed

    if(.not.allocated(nZ)) allocate(nZ(dnodes))
    if(.not.allocated(nth)) allocate(nth(dnodes))
    if(.not.allocated(newZ)) allocate(newZ(dnodes))
    if(.not.allocated(spmZ)) allocate(spmZ(dnodes))
    if(.not.allocated(spmH)) allocate(spmH(dnodes))
    if(.not.allocated(cumDisp)) allocate(cumDisp(upartN))
    if(.not.allocated(nsed)) allocate(nsed(dnodes,totgrn))
    if(.not.allocated(Qs_inS)) allocate(Qs_inS(dnodes,totgrn))
    if(.not.allocated(change_localS)) allocate(change_localS(dnodes,totgrn))
   
    ! Define paramaters
    if(simulation_time==time_start.or.update3d)then
      if(simulation_time==time_start) iter=0
      call connectStrata2TIN
      ! Build active layer
      call buildActiveLayer
      spmZ=tcoordZ
      cumDisp=0.
      CFL_diffusion=display_interval
      idiff=0
      if(simulation_time==time_start) time_step=0. 
      if(simulation_time==time_start) time_display=time_start
      if(Cefficiency>0..or.stream_ero>0.) Tforce=1
    endif

    cpl_time=min(cpl1_time,cpl2_time)
    cpl_time=min(cpl_time,time_end)

    do while(simulation_time<cpl_time+0.001)

      if(.not.allocated(filldem)) allocate(filldem(dnodes))
      if(.not.allocated(watercell)) allocate(watercell(dnodes))

      ! Update borders
      call update_grid_borders

      ! Perform depressionless water filling algo
      call planchon_dem_fill_algorithm

      ! Find network tree based on Braun & Willet 2013
      call define_landscape_network
      
      ! Define subcathcment partitioning
      call compute_subcatchment

      ! Define load balancing
      call bcast_loadbalancing

      if(simulation_time==time_start.or.update3d) newZ=spmZ
      if(pet_id==0)print*,'Current time:',simulation_time
      
      ! Visualisation surface
      if(simulation_time>=time_display)then
        call visualise_surface_changes(iter)
        call visualise_drainage_changes(iter)
        layNb=layNb+1
        call visualise_strata(iter)
        call mpi_barrier(badlands_world,rc)
        if(pet_id==0)print*,'Creating output: ',int(simulation_time)
        time_display=time_display+display_interval
        iter=iter+1
        newZ=spmZ
        cumDisp=0.
      endif
      ! Get time step size for hillslope process and stream power law
      call CFL_conditionS 

      ! Geomorphological evolution
      call geomorphic_evolutionS
      
      ! Advance time
      simulation_time=simulation_time+time_step
      ! Update sea-level
      if(gsea%sealevel) call eustatism
      ! Apply displacement
      if(disp%event>0.and..not.disp3d)then 
        call compute_vertical_displacement
!         call compute_stratal_displacement
      endif

      ! Merge local geomorphic evolution      
      call mpi_allreduce(nZ,spmZ,dnodes,mpi_double_precision,mpi_max,badlands_world,rc)

      ! Merge local stratigraphic evolution     
      call mpi_allreduce(nth,alay_dthick,dnodes,mpi_double_precision,mpi_max,badlands_world,rc) 
      
      do ks=1,totgrn
        tsed(1:dnodes)=nsed(1:dnodes,ks)
        call mpi_allreduce(tsed,gsed,dnodes,mpi_double_precision,mpi_max,badlands_world,rc) 
        alay_dsed(1:dnodes,ks)=gsed(1:dnodes)
      enddo    

      ! Update stratigraphic layer
      call update_stratigraphy_layer
      update3d=.false.
!       stop

    enddo

    if(simulation_time>=time_end)then
      if(time_display<=time_end)then
        call visualise_surface_changes(iter)
        call visualise_drainage_changes(iter)
        layNb=layNb+1
        call visualise_strata(iter)
      endif
      if(pet_id==0)print*,'simulation time: ',int(time_end)
    endif

    ! Update elevation for coupling component assuming infinite boundary
    do k=1,dnodes
      if(voronoiCell(k)%border==0) tcoordZ(k)=spmZ(k)
    enddo
    call update_grid_borders

  end subroutine stratgeomorph
  ! =====================================================================================

  subroutine CFL_conditionS

    integer::k,p,id,lid
    real(kind=8)::denom,distance,dh,ldt,dt

    ! Hillslope CFL conditions
    if(simulation_time==time_start.or.idiff==0.or.sediments(1)%diffa>0.)then 
      call CFLdiffusionS
      CFL_diffusion=display_interval
      idiff=1
      if(CFL_diffusion==0.) time_step=display_interval
    else
      time_step=CFL_diffusion
    endif
    
    ! CFL factor for stream power law (bedrock incision)
    if(Cerodibility>0.)then
      do lid=1,localNodes
        p=localNodesGID(lid)
        k=stackOrder(p) 
        if(voronoiCell(k)%border==0)then
          id=receivers(k)
          if(k/=id)then
            dh=spmZ(k)-spmZ(id)
            distance=sqrt((tcoordX(k)-tcoordX(id))**2.0+(tcoordY(k)-tcoordY(id))**2.0)
            if(spmZ(k)>spmZ(id).and.discharge(k)>0.)then
              denom=Cerodibility*discharge(k)**spl_m*(dh/distance)**(spl_n-1.0)
              time_step=min(time_step,distance/denom)
            endif
          endif
        endif
      enddo
    endif

    if(time_step>display_interval) time_step=display_interval
    if(time_step<force_time) time_step=force_time
    
    ! Get maximum time step for stability
    ldt=time_step    
    call mpi_allreduce(ldt,dt,1,mpi_double_precision,mpi_min,badlands_world,rc)
    time_step=dt

    ! Check time-step in relation to coupling time
    if(simulation_time+time_step>cpl_time) time_step=cpl_time-simulation_time+1.e-4
    if(simulation_time+time_step>time_display) time_step=time_display-simulation_time+1.e-4
    if(simulation_time+time_step>time_end) time_step=time_end-simulation_time+1.e-4
      
  end subroutine CFL_conditionS
  ! =====================================================================================

  subroutine geomorphic_evolutionS

    integer::k,id,rcv,lid,p,q,m,ks
    integer::stat(mpi_status_size),ierr,req(localNodes),r
    
    real(kind=8),dimension(totgrn)::LDL,SPL,Qs1,Qsr

    real(kind=8)::mtime,dt,maxtime,th
    real(kind=8)::distance,diffH,maxh

    max_time=time_step
    Qs_inS=0.0
    nZ=-1.e6
    nth=0.0
    nsed=0.0
    change_localS=-1.e6
    r=0

    do lid=localNodes,1,-1
      k=localNodesGID(lid)
      id=stackOrder(k) 
    
      ! Receive child influx
      if(rcvprocNb(id)>0)then
        do p=1,rcvprocNb(id)
          call mpi_recv(Qsr,totgrn,mpi_double_precision,rcvprocID(id,p),rcvsendID(id,p),badlands_world,stat,ierr)
          do ks=1,totgrn
            Qs_inS(id,ks)=Qs_inS(id,ks)+Qsr(ks) 
          enddo
        enddo
      endif

      if(voronoiCell(id)%btype<0)then
        
        ! Define local parameters
        rcv=receivers(id)
        distance=sqrt((tcoordX(id)-tcoordX(rcv))**2.0+(tcoordY(id)-tcoordY(rcv))**2.0)

        ! Creep processes
        call hillslope_fluxS(id,LDL,diffH)
    
        ! Stream Power Law (detachment-limited) - bedrock incision
        SPL=0.
        Qs1=0.
!         if(Cerodibility>0.) call detachmentlimitedS(id,rcv,distance,diffH,SPL,Qs1)

        ! Detachment-limited condition
        if(Cerodibility>0.)then
          do ks=1,totgrn
            Qs_inS(rcv,ks)=Qs_inS(rcv,ks)+Qs1(ks)
          enddo
        endif
        do ks=1,totgrn
          if(perosive==1.and.SPL(ks)>0.) SPL(ks)=0.
          change_localS(id,ks)=SPL(ks)+LDL(ks)
!           if(abs(change_localS(id,ks))<3)then
!           else
!             print*,ks,id,change_localS(id,ks)
!             print*,sediments(1)%diffm,sediments(1)%diffa
!             stop
!           endif
        enddo

      else
        change_localS(id,1:totgrn)=0.
      endif

      ! Send parent outflux
      if(sendprocID(id)>=0)then
        r=r+1
        call mpi_isend(Qs_inS(rcv,1:totgrn),totgrn,mpi_double_precision,sendprocID(id),id,badlands_world,req(r),ierr)
      endif
    enddo

    ! Wait for each process to finish sending information
    do k=1,r
      call mpi_wait(req(k),stat,ierr)
    enddo

    ! Adapt time step to ensure stability
    if(Tforce==0)then
      maxtime=max_time
    else
      maxtime=time_step
    endif

    if(Tforce==0)then
      do lid=1,localNodes
        k=localNodesGID(lid)
        id=stackOrder(k) 
        rcv=receivers(id)
        th=0.0
        do ks=1,totgrn
          th=th+change_localS(id,ks)
        enddo
        if(tcoordX(id)<minx.or.tcoordX(id)>maxx.or.tcoordY(id)<miny &
            .or.tcoordY(id)>maxy) change_localS(id,1:totgrn)=0.

        if(th<0.)then
          mtime=-active_thick/th
          maxtime=min(maxtime,mtime)
        endif
!         if(lid==44)print*,change_localS(id,1),th,-active_thick/th
        
        if(rcv==id.and.th>0..and.voronoiCell(id)%btype<0)then 
          if(watercell(id)==0.)then
            maxh=0.
            do q=1,delaunayVertex(id)%ngbNb
              m=delaunayVertex(id)%ngbID(q)
              maxh=max(maxh,spmZ(m)-spmZ(id))
            enddo
            if(maxh>0.)then
              mtime=maxh/th 
              maxtime=min(maxtime,mtime)
            endif
          else
            mtime=watercell(id)/th
            maxtime=min(maxtime,mtime)
          endif
        endif

      enddo
    endif
    call mpi_allreduce(maxtime,dt,1,mpi_double_precision,mpi_min,badlands_world,rc)
    time_step=dt   
!     print*,time_step
    if(time_step<force_time) time_step=force_time

    ! Perform elevation and regolith evolution
    do lid=1,localNodes
      id=localNodesGID(lid)
      k=stackOrder(id) 
      rcv=receivers(k)
      if(voronoiCell(k)%border==0)then 
        if(tcoordX(k)==minx.and.bounds(3)==0)change_local(k)=0. 
        if(tcoordX(k)==maxx.and.bounds(4)==0)change_local(k)=0.
        if(tcoordY(k)==miny.and.bounds(2)==0)change_local(k)=0.
        if(tcoordY(k)==maxy.and.bounds(1)==0)change_local(k)=0.
        ! Update surface
        th=0.0
        do ks=1,totgrn
          th=th+time_step*change_localS(k,ks)
          nsed(k,ks)=alay_dsed(k,ks)+time_step*change_localS(k,ks)
          if(nsed(k,ks)<0.0)then
            print*,'fefefe',lid,k,ks,alay_dsed(k,ks),time_step,change_localS(k,ks),nsed(k,ks)
          endif
        enddo
        nZ(k)=spmZ(k)+th
        nth(k)=alay_dthick(k)+th
        if(nth(k)<0.0)then
          print*,'fafafa',lid,k,nth(k),alay_dthick(k),th
          stop
        endif
      endif
    enddo


  end subroutine geomorphic_evolutionS
  ! =====================================================================================

!   subroutine detachmentlimitedS(id,rcv,distance,maxh,ST,Qs)

!     integer::id,rcv
!     real(kind=8)::ST,dh,distance,Qs,maxh

    ! Stream Power Law (detachment-limited) - bedrock incision
!     dh=0.95*(spmZ(id)-spmZ(rcv))
!     if(dh<0.001)dh=0.
!     ST=0.
!     if(rcv/=id.and.dh>0.)then
!       if(watercell(id)<0.001.and.spmZ(id)>=gsea%actual_sea)&    
!         ST=-Cerodibility*(discharge(id))**spl_m*(dh/distance)**spl_n
!     endif

!     ! Stability criteria for deposition solutions
!     if(ST==0..and.perosive==0)then
!       if(maxh>0..and.spmZ(id)<gsea%actual_sea)then
!         if(Tforce==0)then
!           if(Qs_in(id)*max_time/voronoiCell(id)%area<maxh)then
!             ST=Qs_in(id)/voronoiCell(id)%area
!             Qs=0.
!           else
!             ST=maxh/max_time
!             max_time=min(max_time,maxh/ST)
!             Qs=Qs_in(id)-ST*voronoiCell(id)%area
!           endif
!         else
!           if(Qs_in(id)*time_step/voronoiCell(id)%area<maxh)then
!             ST=Qs_in(id)/voronoiCell(id)%area
!             Qs=0.
!           else
!             ST=maxh/time_step
!             Qs=Qs_in(id)-ST*voronoiCell(id)%area
!           endif
!         endif
!       ! Fill depression 
!       elseif(watercell(id)>0.0001.and.rcv/=id)then
!         dh=0.95*watercell(id)
!         if(Tforce==0)then
!           if(Qs_in(id)*max_time/voronoiCell(id)%area<dh)then
!             ST=Qs_in(id)/voronoiCell(id)%area
!             Qs=0.
!           else
!             ST=dh/max_time
!             max_time=min(max_time,dh/ST)
!             Qs=Qs_in(id)-ST*voronoiCell(id)%area
!           endif
!         else
!           if(Qs_in(id)*time_step/voronoiCell(id)%area<dh)then
!             ST=Qs_in(id)/voronoiCell(id)%area
!             Qs=0.
!           else
!             ST=dh/time_step 
!             Qs=Qs_in(id)-ST*voronoiCell(id)%area
!           endif
!         endif
!       elseif(rcv==id)then
!         ST=Qs_in(id)/voronoiCell(id)%area
!         Qs=0.
!       else
!         Qs=Qs_in(id)
!       endif
!     ! For purely erosive case
!     elseif(ST==0..and.perosive==1)then
!       Qs=Qs_in(id)
!     ! Stability criteria for erosion solutions
!     elseif(ST<0.)then
!       if(Tforce==0)then
!         if(-ST*max_time>spmZ(id)-spmZ(rcv)) max_time=min(max_time,-0.99*(spmZ(id)-spmZ(rcv))/ST)
!       else
!         if(-ST*time_step>spmZ(id)-spmZ(rcv)) ST=-0.95*(spmZ(id)-spmZ(rcv))/time_step
!       endif
!       Qs=-ST*voronoiCell(id)%area+Qs_in(id)
!     endif

!   end subroutine detachmentlimitedS
  ! =====================================================================================

end module stratmorph