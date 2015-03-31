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
!       Filename:  Geomorphology.f90
!
!    Description:  Defines morphological formalism used in BADLANDS.
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles     
!
! =====================================================================================

module geomorpho

  use ESMF
  use topology
  use parameters
  use hydrology
  use hydroUtil
  use hillslope
  use watershed
  use facies
  use outspm_surface
  use outspm_drainage
  use external_forces

  implicit none
  
  integer::iter

  real(ESMF_KIND_R8)::cpl_time,max_time
  real(ESMF_KIND_R8)::time1,time2,wt1,wt2

contains

  ! =====================================================================================
  subroutine geomorphology

    integer::k

    if(.not.allocated(nZ)) allocate(nZ(dnodes))
    if(.not.allocated(nH)) allocate(nH(dnodes))
    if(.not.allocated(newZ)) allocate(newZ(dnodes))
    if(.not.allocated(spmZ)) allocate(spmZ(dnodes))
    if(.not.allocated(spmH)) allocate(spmH(dnodes))
    if(.not.allocated(Qs_in)) allocate(Qs_in(dnodes))
    if(.not.allocated(cumDisp)) allocate(cumDisp(upartN))
    if(.not.allocated(change_local)) allocate(change_local(dnodes))
    
    ! Define paramaters
    if(simulation_time==time_start)then
      iter=0
      spmH=0.
      if(stream_ero>0.) spmH=bed_sed_interface
      if(faciesOn==1)then
        if(allocated(facType)) deallocate(facType)
        allocate(facType(dnodes,faciestype))
        facType=0.
      endif
      spmZ=tcoordZ
      if(faciesOn==1.and..not.restartFlag)then
        if(allocated(stratalZ)) deallocate(stratalZ)
        if(allocated(stratalFacies)) deallocate(stratalFacies)
        layerNb=int((time_end-time_start)/layer_interval)+1
        allocate(stratalZ(layerNb,dnodes))
        allocate(stratalFacies(layerNb,dnodes))
        stratalZ=0.
        stratalFacies=0.
        layerID=1
        stratalZ(layerID,:)=spmZ
      endif
      cumDisp=0.
      CFL_diffusion=layer_interval
      time_step=0. 
      time_display=time_start
      layer_time=time_start 
      if(Cefficiency>0..or.stream_ero>0.) Tforce=1
    endif
    cpl_time=min(cpl1_time,cpl2_time)
    cpl_time=min(cpl_time,time_end)

    do while(simulation_time<cpl_time)

!       call ESMF_VMWtime(wt1,rc=rc)
!       call ESMF_VMWtime(time1,rc=rc)
      ! From rough surface find the stream network 
      depressionAlgo=.false.

      ! Update borders
      call update_grid_borders
      if(.not.allocated(filldem)) allocate(filldem(dnodes))
      if(.not.allocated(watercell)) allocate(watercell(dnodes))
      watercell=0.
      filldem=spmZ

      ! Find network tree based on Braun & Willet 2013
      call define_landscape_network
!       call ESMF_VMWtime(time2,rc=rc)
!       if(pet_id==0) print*,'define_landscape_network1',time2-time1
!       call ESMF_VMWtime(time1,rc=rc)

      ! Find drainage area
      call compute_cumulative_discharge
!       call ESMF_VMWtime(time2,rc=rc)
!       if(pet_id==0) print*,'compute_cumulative_discharge1',time2-time1
!       call ESMF_VMWtime(time1,rc=rc)

     ! Depressionless DEM required
      if(depressionAlgo)then 
        ! Define lake height for each local minima
        call define_drained_water_thickness
!       call ESMF_VMWtime(time2,rc=rc)
!       if(pet_id==0) print*,'define_drained_water_thickness',time2-time1
!       call ESMF_VMWtime(time1,rc=rc)
        ! Perform depressionless filling algo
        call planchon_dem_fill_algorithm  
!       call ESMF_VMWtime(time2,rc=rc)
!       if(pet_id==0) print*,'planchon_dem_fill_algorithm',time2-time1
!       call ESMF_VMWtime(time1,rc=rc)
        ! Find network tree based on Braun & Willet 2013
        call define_landscape_network
!       call ESMF_VMWtime(time2,rc=rc)
!       if(pet_id==0) print*,'define_landscape_network2',time2-time1
!       call ESMF_VMWtime(time1,rc=rc)
        ! Find drainage area
        call compute_cumulative_discharge
!       call ESMF_VMWtime(time2,rc=rc)
!       if(pet_id==0) print*,'compute_cumulative_discharge2',time2-time1
!       call ESMF_VMWtime(time1,rc=rc)
      endif

      ! Define subcathcment partitioning
      call compute_subcatchment
!       call ESMF_VMWtime(time2,rc=rc)
!       if(pet_id==0) print*,'compute_subcatchment',time2-time1
!       call ESMF_VMWtime(time1,rc=rc)

      if(simulation_time==time_start) newZ=spmZ
      if(pet_id==0)print*,'Current time:',simulation_time

      ! Visualisation surface
      if(simulation_time>=time_display)then
        call visualise_surface_changes(iter)
        call visualise_drainage_changes(iter)
        if(faciesOn==1)call visualise_facies_change(iter)
        call ESMF_VMBarrier(vm=vm,rc=rc)
        if(pet_id==0)print*,'Creating output: ',int(simulation_time)
        time_display=time_display+display_interval
        iter=iter+1
        newZ=spmZ
        cumDisp=0.
      endif
!       call ESMF_VMWtime(time2,rc=rc)
!       if(pet_id==0) print*,'visualise',time2-time1
!       call ESMF_VMWtime(time1,rc=rc)
      
      ! Add facies layer
      if(simulation_time>=layer_time.and.faciesOn==1)then
        layer_time=layer_time+layer_interval
        layerID=layerID+1
        stratalZ(layerID,:)=stratalZ(layerID-1,:)
        facType=0.0
      endif

      ! Get time step size for hillslope process and stream power law
      call CFL_condition 
!       call ESMF_VMWtime(time2,rc=rc)
!       if(pet_id==0) print*,'CFL_condition',time2-time1
!       call ESMF_VMWtime(time1,rc=rc)
      
      ! Geomorphological evolution
      call geomorphic_evolution
!       call ESMF_VMWtime(time2,rc=rc)
!       if(pet_id==0) print*,'geomorphic_evolution',time2-time1
!       call ESMF_VMWtime(time1,rc=rc)

      ! Advance time
      simulation_time=simulation_time+time_step
      ! Update sea-level
      if(gsea%sealevel) call eustatism
      ! Apply displacement
      if(disp%event>0) call compute_vertical_displacement

!       call ESMF_VMWtime(time2,rc=rc)
!       if(pet_id==0) print*,'compute_vertical_displacement',time2-time1
!       call ESMF_VMWtime(time1,rc=rc)
      ! Merge local geomorphic evolution
      call ESMF_VMAllReduce(vm=vm,sendData=nZ,recvData=spmZ,count=dnodes,&
        reduceflag=ESMF_REDUCE_MAX,rc=rc)
      if(stream_ero>0..or.regoProd>0.) &
        call ESMF_VMAllReduce(vm=vm,sendData=nH,recvData=spmH,count=dnodes,&
        reduceflag=ESMF_REDUCE_MAX,rc=rc)

      if(faciesOn==1) call update_facies
!       call ESMF_VMWtime(time2,rc=rc)
!       if(pet_id==0) print*,'reduce',time2-time1
!       call ESMF_VMWtime(time1,rc=rc)
!       call ESMF_VMWtime(wt2,rc=rc)
!       call ESMF_VMBarrier(vm=vm,rc=rc)
!       if(pet_id==0) print*,'-------- total loop',wt2-wt1
!       if(pet_id==0) stop

    enddo

    if(simulation_time>=time_end)then
      call visualise_surface_changes(iter)
      call visualise_drainage_changes(iter)
      if(faciesOn==1)call visualise_facies_change(iter)
      call ESMF_VMBarrier(vm=vm,rc=rc)
      if(pet_id==0)print*,'simulation time: ',int(time_display)
    endif

    ! Update elevation for coupling component assuming infinite boundary
    do k=1,dnodes
      if(voronoiCell(k)%border==0) tcoordZ(k)=spmZ(k)
    enddo
    call update_grid_borders

  end subroutine geomorphology
  ! =====================================================================================
  subroutine CFL_condition

    integer::k,p,id,lid
    real(ESMF_KIND_R8)::denom,distance,dh,ldt(1),dt(1)

    ! Hillslope CFL conditions
    if(simulation_time==time_start.or.Cdiffusion_d(1)>0..or.Cdiffusion_nl(1)>0..or.CFL_diffusion<1.)then 
      call CFLdiffusion
      CFL_diffusion=display_interval
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
    if(time_step>force_time.and.force_time>1.) time_step=force_time
    if(time_step<force_time) time_step=force_time
    
    ! Get maximum time step for stability
    ldt(1)=time_step
    call ESMF_VMAllReduce(vm=vm,sendData=ldt,recvData=dt,count=1,&
      reduceflag=ESMF_REDUCE_MIN,rc=rc)
    time_step=dt(1)

    ! Check time-step in relation to coupling time
    if(simulation_time+time_step>cpl_time) time_step=cpl_time-simulation_time+1.e-4
    if(simulation_time+time_step>layer_time.and.faciesOn==1) time_step=layer_time-simulation_time+1.e-4
    if(simulation_time+time_step>time_display) time_step=time_display-simulation_time+1.e-4
    if(simulation_time+time_step>time_end) time_step=time_end-simulation_time+1.e-4
      
  end subroutine CFL_condition
  ! =====================================================================================
  subroutine geomorphic_evolution

    integer::k,id,rcv,lid,p,q,m
    integer::stat(mpi_status_size),ierr,req(localNodes),r
    real(ESMF_KIND_R8)::Ps,STL,SPL,STC,LDL,NDL,DDD,ST,Qsr,mtime,dt(1),maxtime(1)
    real(ESMF_KIND_R8)::distance,diffH,Qs1,Qs2,Qs3,S_d,S_t,N_dt,maxh

    max_time=time_step
    Qs_in=0.0
    nZ=-1.e6
    if(stream_ero>0..or.regoProd>0.) nH=-1.e6
    change_local=-1.e6
    r=0
    do lid=localNodes,1,-1
      k=localNodesGID(lid)
      id=stackOrder(k) 

      ! Receive child influx
      if(rcvprocNb(id)>0)then
        do p=1,rcvprocNb(id)
          call mpi_recv(Qsr,1,mpi_double_precision,rcvprocID(id,p),rcvsendID(id,p),badlands_world,stat,ierr)
          Qs_in(id)=Qs_in(id)+Qsr 
        enddo
      endif

      if(voronoiCell(id)%btype<0)then

        ! Define local parameters
        rcv=receivers(id)
        distance=sqrt((tcoordX(id)-tcoordX(rcv))**2.0+(tcoordY(id)-tcoordY(rcv))**2.0)

        ! Creep processes
        call hillslope_flux(id,LDL,DDD,NDL,diffH)
        ! Stream Power Law (detachment-limited) - bedrock incision
        SPL=0.
        Qs1=0.
        if(Cerodibility>0.) call detachmentlimited(id,rcv,distance,diffH,SPL,Qs1)
        
        ! Sediment Transport Law (transport-limited)
        STL=0.
        Qs2=0.
        if(Cefficiency>0.) call transportlimited(id,rcv,distance,diffH,STL,Qs2) 
        
        ! Sediment Transport Capacity model
        STC=0.
        Qs3=0.
        if(stream_ero>0.) call streamcapacity(id,rcv,distance,diffH,STC,Qs3)

        ! Hybrid sediment flux 
        ST=0
        if(Cefficiency>0..and.Cerodibility>0..and.regoProd==0.)then
          N_dt=2. 
          if(discharge(id)>0.)then
            ! Detachment-limited channel slope
            S_d=discharge(id)**(-spl_m/spl_n)*abs(tvertDisp(id)/Cerodibility)**(1/spl_n)
            ! Transport-limited channel slope
            S_t=discharge(id)**(-(stl_m-1)/stl_n)*abs(Fracbed*tvertDisp(id)/Cefficiency)**(1/stl_n)
            ! Transition number
            if(S_t/=0) N_dt=S_d/S_t
          endif
          ! Transport-limited 
          if(N_dt<=1.)then
            ST=STL
            Qs_in(rcv)=Qs_in(rcv)+Qs2
          ! Detachment-limited 
          else
            ST=SPL
            Qs_in(rcv)=Qs_in(rcv)+Qs1
          endif
        ! Hybrid with regolith limitation
        elseif(Cefficiency>0..and.Cerodibility>0..and.regoProd>0.)then
          if(spmH(id)==0.)then
            ST=SPL
            Qs_in(rcv)=Qs_in(rcv)+Qs1
          ! Fluvial processes are limited by the slower rate
          elseif(SPL<0..and.STL<0.)then
            if(SPL>STL)then
              ST=SPL
              Qs_in(rcv)=Qs_in(rcv)+Qs1
            else
              ST=STL
              Qs_in(rcv)=Qs_in(rcv)+Qs2
            endif
          elseif(STL>=0.and.SPL<0.)then
            ST=STL
            Qs_in(rcv)=Qs_in(rcv)+Qs2
          else
            ST=SPL
            Qs_in(rcv)=Qs_in(rcv)+Qs1
          endif
        ! Detachment-limited condition
        elseif(Cerodibility>0.)then
          ST=SPL
          Qs_in(rcv)=Qs_in(rcv)+Qs1
        ! Transport-limited condition
        elseif(Cefficiency>0.)then
          ST=STL
          Qs_in(rcv)=Qs_in(rcv)+Qs2
        ! Capacity model
        elseif(stream_ero>0.)then
          ST=STC
          Qs_in(rcv)=Qs_in(rcv)+Qs3
        endif
        if(perosive==1.and.ST>0.) ST=0.
        change_local(id)=ST+LDL+DDD+NDL 
      else
        change_local(id)=0.
      endif

      ! Send parent outflux
      if(sendprocID(id)>=0)then
        r=r+1
        call mpi_isend(Qs_in(rcv),1,mpi_double_precision,sendprocID(id),id,badlands_world,req(r),ierr)
      endif
    enddo

    ! Wait for each process to finish sending information
    do k=1,r
      call mpi_wait(req(k),stat,ierr)
    enddo

    ! Adapt time step to ensure stability
    if(Tforce==0)then
      maxtime(1)=max_time
    else
      maxtime(1)=time_step
    endif
    if(Tforce==0)then
      do lid=1,localNodes
        k=localNodesGID(lid)
        id=stackOrder(k) 
        rcv=receivers(id)
        if(tcoordX(id)<minx.or.tcoordX(id)>maxx.or.tcoordY(id)<miny &
            .or.tcoordY(id)>maxy) change_local(id)=0.
        if(rcv==id.and.change_local(id)>0..and.voronoiCell(id)%btype<0)then 
          if(watercell(id)==0.)then
            maxh=0.
            do q=1,delaunayVertex(id)%ngbNb
              m=delaunayVertex(id)%ngbID(q)
              maxh=max(maxh,spmZ(m)-spmZ(id))
            enddo
            if(maxh>0.)then
              mtime=maxh/change_local(id)
              maxtime(1)=min(maxtime(1),mtime)
            endif
          else
            mtime=watercell(id)/change_local(id)
            maxtime(1)=min(maxtime(1),mtime)
          endif
        endif
      enddo
    endif
    call ESMF_VMAllReduce(vm=vm,sendData=maxtime,recvData=dt,count=1,&
      reduceflag=ESMF_REDUCE_MIN,rc=rc)
    time_step=dt(1)    

    if(time_step<force_time) time_step=force_time

    ! Perform elevation and regolith evolution
    do lid=1,localNodes
      id=localNodesGID(lid)
      k=stackOrder(id) 
      if(voronoiCell(k)%border==0)then 
        if(tcoordX(k)==minx.and.bounds(3)==0)change_local(k)=0. 
        if(tcoordX(k)==maxx.and.bounds(4)==0)change_local(k)=0.
        if(tcoordY(k)==miny.and.bounds(2)==0)change_local(k)=0.
        if(tcoordY(k)==maxy.and.bounds(1)==0)change_local(k)=0.
        ! Regolith formation
        if(regoProd>0.)then
          ! Soil production
          Ps=regoProd*exp(-spmH(k)/regoDepth) 
          ! Soil production is produced from exposed bedrock 
          nH(k)=spmH(k)+time_step*(change_local(k)+Ps*rock_density/soil_density)
          if(nH(k)<0.)nH(k)=0.
        elseif(stream_ero>0.)then
          nH(k)=spmH(k)+time_step*change_local(k)
          if(nH(k)<0.)nH(k)=0.
        endif
        ! Update surface
        nZ(k)=spmZ(k)+time_step*change_local(k)
      
      endif
    enddo

  end subroutine geomorphic_evolution
  ! =====================================================================================
  subroutine detachmentlimited(id,rcv,distance,maxh,ST,Qs)

    integer::id,rcv
    real(ESMF_KIND_R8)::ST,dh,distance,Qs,maxh

    ! Stream Power Law (detachment-limited) - bedrock incision
    dh=0.95*(spmZ(id)-spmZ(rcv))
    if(dh<0.001)dh=0.
    ST=0.
    if(rcv/=id.and.dh>0.)then
      if(watercell(id)<0.001.and.spmZ(id)>=gsea%actual_sea)&    
        ST=-Cerodibility*(discharge(id))**spl_m*(dh/distance)**spl_n
    endif

    ! Stability criteria for deposition solutions
    if(ST==0..and.perosive==0)then
      if(maxh>0..and.spmZ(id)<gsea%actual_sea)then
        if(Tforce==0)then
          if(Qs_in(id)*max_time/voronoiCell(id)%area<maxh)then
            ST=Qs_in(id)/voronoiCell(id)%area
            Qs=0.
          else
            ST=maxh/max_time
            max_time=min(max_time,maxh/ST)
            Qs=Qs_in(id)-ST*voronoiCell(id)%area
          endif
        else
          if(Qs_in(id)*time_step/voronoiCell(id)%area<maxh)then
            ST=Qs_in(id)/voronoiCell(id)%area
            Qs=0.
          else
            ST=maxh/time_step
            Qs=Qs_in(id)-ST*voronoiCell(id)%area
          endif
        endif
      ! Fill depression 
      elseif(watercell(id)>0.0001.and.rcv/=id)then
        dh=0.95*watercell(id)
        if(Tforce==0)then
          if(Qs_in(id)*max_time/voronoiCell(id)%area<dh)then
            ST=Qs_in(id)/voronoiCell(id)%area
            Qs=0.
          else
            ST=dh/max_time
            max_time=min(max_time,dh/ST)
            Qs=Qs_in(id)-ST*voronoiCell(id)%area
          endif
        else
          if(Qs_in(id)*time_step/voronoiCell(id)%area<dh)then
            ST=Qs_in(id)/voronoiCell(id)%area
            Qs=0.
          else
            ST=dh/time_step 
            Qs=Qs_in(id)-ST*voronoiCell(id)%area
          endif
        endif
      elseif(rcv==id)then
        ST=Qs_in(id)/voronoiCell(id)%area
        Qs=0.
      else
        Qs=Qs_in(id)
      endif
    ! For purely erosive case
    elseif(ST==0..and.perosive==1)then
      Qs=Qs_in(id)
    ! Stability criteria for erosion solutions
    elseif(ST<0.)then
      if(Tforce==0)then
        if(-ST*max_time>spmZ(id)-spmZ(rcv)) max_time=min(max_time,-0.99*(spmZ(id)-spmZ(rcv))/ST)
      else
        if(-ST*time_step>spmZ(id)-spmZ(rcv)) ST=-0.95*(spmZ(id)-spmZ(rcv))/time_step
      endif
      Qs=-ST*voronoiCell(id)%area+Qs_in(id)
    endif

  end subroutine detachmentlimited
  ! =====================================================================================
  subroutine streamcapacity(id,rcv,distance,maxh,ST,Qs)

    integer::id,rcv
    real(ESMF_KIND_R8)::ST,dh,width,distance,Qs,maxh,Qe

    dh=0.95*(spmZ(id)-spmZ(rcv))
    if(dh<0.001)dh=0.
    ST=0.

    if(watercell(id)<0.001.and.spmZ(id)>gsea%actual_sea.and.distance>0.)then
      ! Carrying capacity
      Qe=stream_ero*dh*discharge(id)/distance
      ! Compare sediment flux to carrying capacity
      if(Qs_in(id)>=Qe)then
        ST=(Qs_in(id)-Qe)/voronoiCell(id)%area
        Qs=Qe
        if(maxh>0..and.ST*time_step>maxh)then
          ST=maxh/time_step
          Qs=Qs_in(id)-ST*voronoiCell(id)%area
        elseif(maxh<=0.)then
          ST=0.
          Qs=Qs_in(id)
        endif
      else
        ! Channel width
        width=chan_width*(discharge(id))**chan_exp
        ! For alluvial channel
        if(spmH(id)>0.)then
!           ST=(Qs_in(id)-Qe)*(distance/sed_length)/voronoiCell(id)%area
          ST=(Qs_in(id)-Qe)/(width*sed_length)
          if(-ST*time_step>spmZ(id)-spmZ(rcv)) ST=-0.95*(spmZ(id)-spmZ(rcv))/time_step
          Qs=Qs_in(id)-ST*voronoiCell(id)%area
        ! For bedrock channel
        else
!           ST=(Qs_in(id)-Qe)*(distance/bed_length)/voronoiCell(id)%area
          ST=(Qs_in(id)-Qe)/(width*bed_length)
          if(-ST*time_step>spmZ(id)-spmZ(rcv)) ST=-0.95*(spmZ(id)-spmZ(rcv))/time_step
          Qs=Qs_in(id)-ST*voronoiCell(id)%area
        endif
      endif

    ! Force deposition underwater & in lake
    elseif((watercell(id)>0.001.or.spmZ(id)<gsea%actual_sea).and.Qs_in(id)>0..and.distance>0.)then
      if(maxh>0..and.spmZ(id)<gsea%actual_sea)then
        if(Qs_in(id)*time_step/voronoiCell(id)%area<maxh)then
          ST=Qs_in(id)/voronoiCell(id)%area
          Qs=0.
        else
          ST=maxh/time_step
          Qs=Qs_in(id)-ST*voronoiCell(id)%area
        endif
      ! Fill depression
      elseif(watercell(id)>0.001.and.rcv/=id)then
        dh=0.95*watercell(id)
        if(Qs_in(id)*time_step/voronoiCell(id)%area<dh)then
          ST=Qs_in(id)/voronoiCell(id)%area
          Qs=0.
        else
          ST=dh/time_step 
          Qs=Qs_in(id)-ST*voronoiCell(id)%area
        endif
      elseif(rcv==id)then
        ST=Qs_in(id)/voronoiCell(id)%area
        Qs=0.
      else
        ST=0.
        Qs=Qs_in(id)
      endif
    ! Force deposition at local minima
    elseif(Qs_in(id)>0..and.distance==0.)then
      ST=Qs_in(id)/voronoiCell(id)%area
      Qs=0.
    ! Otherwise propagate the sediment flux downstream
    else
      ST=0.
      Qs=Qs_in(id)
    endif

  end subroutine streamcapacity
  ! =====================================================================================
  subroutine transportlimited(id,rcv,distance,maxh,ST,Qs) 

    integer::id,rcv
    real(ESMF_KIND_R8)::ST,distance,dh,QsOut,Qs,maxh,area

    ST=0.
    QsOut=0.
    Qs=0.
    dh=0.95*(spmZ(id)-spmZ(rcv))
    if(dh<0.001)dh=0.
    area=voronoiCell(id)%area
   
    if(dh>0..and.spmZ(id)>=gsea%actual_sea.and.watercell(id)<0.001.and.distance>0.)then
      QsOut=Cefficiency*(discharge(id))**stl_m*(dh/distance)**stl_n
      ! In regolith production is considered check that there is a regolith thickness
      if(regoProd>0..and.spmH(id)==0.) Qsout=0.
      ST=(Qs_in(id)-QsOut)/area

      ! In case of deposition limit sediment transport by maximum neighbor elevation
      if(ST>0.)then
        if(maxh>0..and.ST*time_step>maxh)then
          ST=maxh/time_step
          Qs=Qs_in(id)-ST*area
        elseif(maxh<=0..and.ST>0.)then
          ST=0.
          Qs=Qs_in(id)
        else
          Qs=QsOut
        endif
      
      ! Limit erosion from relief elevation
      elseif(ST<0..and.regoProd==0..and.dh>=0.)then
        if(-ST*time_step>dh)then 
          ST=-dh/time_step
          Qs=Qs_in(id)-ST*area
        else
          Qs=QsOut
        endif
    
      ! Limit erosion from regolith thickness or relief elevation
      elseif(ST<0..and.regoProd>0..and.dh>=0.)then
        if(spmH(id)==0.)then 
          ST=0.
          Qs=Qs_in(id)
        elseif(maxh>0.)then
          if(-ST*time_step>spmH(id))then 
            ST=-spmH(id)/time_step
            Qs=Qs_in(id)-ST*area
          else
            Qs=QsOut
          endif
          if(-ST*time_step>spmZ(id)-spmZ(rcv))then 
            ST=-(spmZ(id)-spmZ(rcv))/time_step
            Qs=Qs_in(id)-ST*area
          else
            Qs=QsOut
          endif
        else
          Qs=QsOut
        endif

      elseif(ST==0.)then
        Qs=Qs_in(id)
      endif

    ! Force deposition below sea-level & in lake
    elseif((watercell(id)>0.001.or.spmZ(id)<gsea%actual_sea).and.Qs_in(id)>0..and.distance>0.)then
      if(maxh>0..and.spmZ(id)<gsea%actual_sea)then
        if(Qs_in(id)*time_step/area<maxh)then
          ST=Qs_in(id)/area
          Qs=0.
        else
          ST=maxh/time_step
          Qs=Qs_in(id)-ST*area
        endif
      ! Fill depression
      elseif(watercell(id)>0.001.and.rcv/=id)then
        dh=0.95*watercell(id)
        if(Qs_in(id)*time_step/area<dh)then
          ST=Qs_in(id)/area
          Qs=0.
        else
          ST=dh/time_step 
          Qs=Qs_in(id)-ST*area
        endif
      elseif(rcv==id)then
        ST=Qs_in(id)/area
        Qs=0.
      else
        ST=0.
        Qs=Qs_in(id)
      endif

    ! Force deposition on local minima
    elseif(Qs_in(id)>0..and.rcv==id)then
      ST=Qs_in(id)/area
      Qs=0.

    ! Otherwise propagate the sediment flux downstream
    elseif(Qs_in(id)>0.)then
      Qs=Qs_in(id)
      ST=0.
    endif

  end subroutine transportlimited
  ! =====================================================================================

end module geomorpho