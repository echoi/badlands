! =====================================================================================
! =====================================================================================
! BADLANDS (BAsin and LAndscape Dynamics)
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
!       Filename:  badlands_App.f90
!
!    Description:  Top level BADLANDS Application Driver
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles     
!
! =====================================================================================
program BADLANDS_Application

  use geomesh
  use parameters
  use spm_restart

  ! Model Component Registration Routines
  use esmfSPM,only:spm_register
  use esmfOcean,only:ocean_register
  use esmfEarth,only:earth_register
  use esmfCoupler,only:cpl_register

  implicit none
  
  integer::urc,localrc
  integer(ESMF_KIND_I4)::opt

  real(ESMF_KIND_R8)::time1,time2

  ! set rc = ESMF_SUCCES
  rc=ESMF_SUCCESS

  ! Initialize framework and get back default global VM
  call ESMF_Initialize(vm=vm,defaultlogfilename="badlands.Log",rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,&
      file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  call ESMF_LogWrite("BAsin & LAndscape DynamicS ESMF Framework",ESMF_LOGMSG_INFO,rc=localrc)

  ! Get number of PETs we are running with
  call ESMF_VMGet(vm,petCount=npets,localPet=pet_id,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,&
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  if(pet_id==0)print*,'-------------------------'

  ! Get experiment file name
  xmlfile=''
  opt=iargc()
  if(opt<1)then
    call ESMF_LogSetError(rcToCheck=ESMF_RC_VAL_WRONG, &
      msg="Wrong command line, use: mpirun -n X ./badlands <xml-file-name>", &
      line=__LINE__,file=__FILE__) 
    call ESMF_Finalize(endflag=ESMF_END_ABORT)
  elseif(opt==1)then
    call getarg(1,xmlfile)
  endif

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!    Mesh & Grid generation section
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
  ! Define simulation meshes
  call GeoMesher
  if(pet_id==0)print*,'-------------------------'

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!    Grid component creation section
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
  ! Create model components and coupler
  ocean="Ocean Model"
  oceanComp=ESMF_GridCompCreate(name=ocean,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,&
    file=__FILE__,rcToReturn=rc))call ESMF_Finalize(endflag=ESMF_END_ABORT)

  earth="Geodynamic Model"
  earthComp=ESMF_GridCompCreate(name=earth,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,&
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  spm="Surface Process Model"
  spmComp=ESMF_GridCompCreate(name=spm,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,&
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
  
  coupler1="BADLANDS Coupler Ocean 2 SPM"
  cpl1=ESMF_CplCompCreate(name=coupler1,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,&
    file=__FILE__,rcToReturn=rc))call ESMF_Finalize(endflag=ESMF_END_ABORT)

  coupler2="BADLANDS Coupler SPM 2 Ocean"
  cpl2=ESMF_CplCompCreate(name=coupler2,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,&
    file=__FILE__,rcToReturn=rc))call ESMF_Finalize(endflag=ESMF_END_ABORT)

  coupler3="BADLANDS Coupler Earth 2 SPM"
  cpl3=ESMF_CplCompCreate(name=coupler3,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,&
    file=__FILE__,rcToReturn=rc))call ESMF_Finalize(endflag=ESMF_END_ABORT)

  coupler4="BADLANDS Coupler SPM 2 Earth"
  cpl4=ESMF_CplCompCreate(name=coupler4,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__,&
    file=__FILE__,rcToReturn=rc))call ESMF_Finalize(endflag=ESMF_END_ABORT)

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!  Register services section
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
  call ESMF_GridCompSetServices(oceanComp,userRoutine=ocean_register,userRc=urc,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
  if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 

  call ESMF_GridCompSetServices(spmComp,spm_register,userRc=urc,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
  if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 

  call ESMF_GridCompSetServices(earthComp,earth_register,userRc=urc,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
  if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 

  call ESMF_CplCompSetServices(cpl1,cpl_register,userRc=urc,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
  if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 

  call ESMF_CplCompSetServices(cpl2,cpl_register,userRc=urc,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
  if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
  call ESMF_VMWtime(time2,rc=rc)

  call ESMF_CplCompSetServices(cpl3,cpl_register,userRc=urc,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
  if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 

  call ESMF_CplCompSetServices(cpl4,cpl_register,userRc=urc,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
  if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
  call ESMF_VMWtime(time2,rc=rc)

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!  Create and initialize a clock.
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
  ! Initialize geological calendar 
  ! The geological calendar uses: 
  ! - seconds for years
  ! - days for thousands years
  ! - years for millions years
  GeologyCalendar=ESMF_CalendarCreate(secondsPerDay=1000,daysPerYear=1000,name="GeoCal",rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 

  ! Initialize time interval to 100 years
  call ESMF_TimeIntervalSet(timeStep,d=500,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 

  ! Initialize start time to 0
  call ESMF_TimeSet(startTime,yy=0000,s=00,calendar=GeologyCalendar,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 

  ! We run the simulation for 1 millenium
  call ESMF_TimeSet(stopTime,yy=0001,d=00,s=01,calendar=GeologyCalendar,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 

  ! Initialize the clock with the above values
  clock=ESMF_ClockCreate(timeStep,startTime,stopTime=stopTime,name="GeoClock",rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!  Init section
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

  ! Initialise Import/Export states for the Ocean Model 
  call ESMF_VMWtime(time1,rc=rc)
  OstateExp=ESMF_StateCreate(name="ocean export", &
      stateintent=ESMF_STATEINTENT_EXPORT,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 

  OstateImp=ESMF_StateCreate(name="ocean import", &
      stateintent=ESMF_STATEINTENT_EXPORT,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 

  call ESMF_GridCompInitialize(oceanComp,importState=OstateImp,exportState=OstateExp,clock=clock,userRc=urc,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
  if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
 
  ! Initialise Import/Export states for the Earth Model
  EstateImp=ESMF_StateCreate(name="earth import", &
      stateintent=ESMF_STATEINTENT_EXPORT,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 

  EstateExp=ESMF_StateCreate(name="earth export", &
      stateintent=ESMF_STATEINTENT_EXPORT,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 

  call ESMF_GridCompInitialize(earthComp,importState=EstateImp,exportState=EstateExp,clock=clock,userRc=urc,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
  if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 

  ! Initialise Import/Export states for the SPM Model
  SstateImp=ESMF_StateCreate(name="spm import", &
      stateintent=ESMF_STATEINTENT_EXPORT,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 

  SstateExp=ESMF_StateCreate(name="spm export", &
      stateintent=ESMF_STATEINTENT_EXPORT,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 

  call ESMF_GridCompInitialize(spmComp,importState=SstateImp,exportState=SstateExp,clock=clock,userRc=urc,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
  if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

  ! First coupling route import Ocean State (OstateExp) to SPM State (SstateImp)
  ! Note that the coupler's import is oceanComp's export
  call ESMF_CplCompInitialize(cpl1,importState=OstateExp,exportState=SstateImp,phase=1,clock=clock,userRc=urc,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
  if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 

  ! Second coupling route import SPM State (SstateExp) to Ocean State (OstateImp)
  ! Note that the coupler's import is spmComp's export
  call ESMF_CplCompInitialize(cpl2,importState=SstateExp,exportState=OstateImp,phase=2,clock=clock, &
      userRc=urc,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
  if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
  call ESMF_VMWtime(time2,rc=rc)

  ! Third coupling route import Earth State (EstateExp) to SPM State (SstateImp)
  ! Note that the coupler's import is earthComp's export
  call ESMF_CplCompInitialize(cpl3,importState=EstateExp,exportState=SstateImp,phase=3,clock=clock,userRc=urc,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
  if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 

  ! Fourth coupling route import SPM State (SstateExp) to Earth State (EstateImp)
  ! Note that the coupler's import is spmComp's export
  call ESMF_CplCompInitialize(cpl4,importState=SstateExp,exportState=EstateImp,phase=4,clock=clock, &
      userRc=urc,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
  if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
  call ESMF_VMWtime(time2,rc=rc)

  if(pet_id==0) print*,'BADLANDS ESMF Components Initialized (s) ',time2-time1
  if(pet_id==0)print*,'-------------------------'

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!     Run section
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
  oceanFlag=.false.
  geodynamicFlag=.true.
  if(restartFlag)then
    call getSPM_hdf5topography
    if(faciesOn==1) call getSPM_hdf5strata
    updateSPM_elevation=.true.
  else
    updateSPM_elevation=.false.
  endif
  
  do while(simulation_time<time_end)

    ! Get the current topographic state, performs Ocean velocity fields calculation and
    ! prepare export arrays for the SPM model
    ! Note that first pass is not taking into account changes due to SPM as nothing has 
    ! been done with surface processes
    if(oceanFlag)then
      call ESMF_GridCompRun(oceanComp,importState=OstateImp,exportState=OstateExp,clock=clock,userRc=urc,rc=localrc)
      if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
        file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
      if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
        file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT)  
      ! First coupling route is regridding the Ocean velocity field to the SPM model mesh
      call ESMF_CplCompRun(cpl1,importState=OstateExp,exportState=SstateImp,phase=1,clock=clock,userRc=urc,rc=localrc)
      if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
        file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
      if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
        file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
    endif
    
    ! Get the current topographic state, performs Geodynamic Evolution calculation and
    ! prepare export vertical displacement arrays for the SPM model
    ! Note that first pass is not taking into account changes due to SPM as nothing has 
    ! been done with surface processes
    if(geodynamicFlag)then
      call ESMF_GridCompRun(earthComp,importState=EstateImp,exportState=EstateExp,clock=clock,userRc=urc,rc=localrc)
      if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
        file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
      if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
        file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT)

      ! Second coupling route is regridding the vertical displacement field to the SPM model mesh
      call ESMF_CplCompRun(cpl3,importState=EstateExp,exportState=SstateImp,phase=3,clock=clock,userRc=urc,rc=localrc)
      if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
        file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
      if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
        file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
    endif
    
    ! Get the regridded ocean velocity field and perform stratigraphy and geomorphology changes
    ! according to Surface Processes. Then prepare the new topography export arrays for the Ocean model
    call ESMF_GridCompRun(spmComp,importState=SstateImp,exportState=SstateExp,clock=clock,userRc=urc,rc=localrc)
    if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
      file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
    if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
      file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 

    ! Third coupling route is regridding the SPM new topography field to the Ocean model mesh
    if(oceanFlag)then
      call ESMF_CplCompRun(cpl2,importState=SstateExp,exportState=OstateImp,phase=2,clock=clock,userRc=urc,rc=localrc)
      if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
        file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
      if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
        file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
    endif

    ! Fourth coupling route is regridding the SPM new topography field to the Earth model mesh
    if(geodynamicFlag)then
      call ESMF_CplCompRun(cpl4,importState=SstateExp,exportState=EstateImp,phase=4,clock=clock,userRc=urc,rc=localrc)
      if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
        file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
      if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
        file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
    endif

    updateSPM_elevation=.true.

  enddo

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!     Finalize section
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

  ! Finalize Ocean Component
  call ESMF_VMWtime(time1,rc=rc)
  call ESMF_GridCompFinalize(oceanComp,importState=OstateImp,exportState=OstateExp,clock=clock,userRc=urc,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
  if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 

  ! Finalize Earth Component
  call ESMF_GridCompFinalize(earthComp,importState=EstateImp,exportState=EstateExp,clock=clock,userRc=urc,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
  if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 

  ! Finalize SPM Component
  call ESMF_GridCompFinalize(spmComp,importState=SstateImp,exportState=SstateExp,clock=clock,userRc=urc,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
  if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 

  ! Finalize Coupling Component Route 1
  call ESMF_CplCompFinalize(cpl1,importState=OstateExp,exportState=SstateImp,phase=1,clock=clock,userRc=urc,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
  if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 

  ! Finalize Coupling Component Route 2
  call ESMF_CplCompFinalize(cpl2,importState=SstateExp,exportState=OstateImp,phase=2,clock=clock,userRc=urc,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
  if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 

  ! Finalize Coupling Component Route 3
  call ESMF_CplCompFinalize(cpl3,importState=EstateExp,exportState=SstateImp,phase=3,clock=clock,userRc=urc,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
  if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 

  ! Finalize Coupling Component Route 4
  call ESMF_CplCompFinalize(cpl4,importState=SstateExp,exportState=EstateImp,phase=4,clock=clock,userRc=urc,rc=localrc)
  if(ESMF_LogFoundError(rcToCheck=localrc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 
  if(ESMF_LogFoundError(rcToCheck=urc,msg=ESMF_LOGERR_PASSTHRU,line=__LINE__, &
    file=__FILE__,rcToReturn=rc)) call ESMF_Finalize(endflag=ESMF_END_ABORT) 

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!     Destroy section
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

  ! Clean Up
  call ESMF_StateDestroy(OstateImp,rc=rc)
  call ESMF_StateDestroy(OstateExp,rc=rc)
  call ESMF_StateDestroy(SstateExp,rc=rc)
  call ESMF_StateDestroy(SstateImp,rc=rc)
  call ESMF_StateDestroy(EstateExp,rc=rc)
  call ESMF_StateDestroy(EstateImp,rc=rc)
  call ESMF_LogWrite("Topo State Import/Export Destroyed",ESMF_LOGMSG_INFO,rc=localrc)

  ! Destroy Clock
  call ESMF_ClockDestroy(clock,rc=rc)
  call ESMF_CalendarDestroy(GeologyCalendar,rc=rc)
  call ESMF_LogWrite("Geological Calendar Destroyed",ESMF_LOGMSG_INFO,rc=localrc)

  ! Destroy Grids and Coupling Components
  call ESMF_GridCompDestroy(oceanComp,rc=rc)
  call ESMF_GridCompDestroy(spmComp,rc=rc)
  call ESMF_CplCompDestroy(cpl1,rc=rc)
  call ESMF_CplCompDestroy(cpl2,rc=rc)
  call ESMF_CplCompDestroy(cpl3,rc=rc)
  call ESMF_CplCompDestroy(cpl4,rc=rc)
  call ESMF_LogWrite("ESMF Meshes Destroyed",ESMF_LOGMSG_INFO,rc=localrc)

  ! Destroy Geometries
  call UnstructuredMeshDestroy
  call ESMF_LogWrite("Unstructured Mesh Arrays Destroyed",ESMF_LOGMSG_INFO,rc=localrc)
  call RegularGridDestroy
  call ESMF_LogWrite("Regular Mesh Arrays Destroyed",ESMF_LOGMSG_INFO,rc=localrc)
  call ESMF_VMWtime(time2,rc=rc)

  if(pet_id==0) print*,'BADLANDS ESMF Finalized (s) ',time2-time1
  if(pet_id==0)print*,'-------------------------'

  call ESMF_Finalize(rc=rc)
  
end program BADLANDS_Application
