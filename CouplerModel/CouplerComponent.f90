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
!    Description:  Defines the Coupler ESMF Component
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles     
!
! =====================================================================================

module esmfCoupler

  use topology
  use parameters

  implicit none

  public cpl_register

  ! Global data
  type(ESMF_RouteHandle),save::ocean2spm
  type(ESMF_RouteHandle),save::spm2ocean
  type(ESMF_RouteHandle),save::earth2spm
  type(ESMF_RouteHandle),save::spm2earth

contains

  ! =====================================================================================
  subroutine cpl_register(comp,rclocal)

    type(ESMF_CplComp)::comp
    integer,intent(out)::rclocal

    rclocal=ESMF_SUCCESS

    ! Register the callback routines
    call ESMF_CplCompSetEntryPoint(comp,ESMF_METHOD_INITIALIZE,cpl_init1,phase=1,rc=rclocal)
    if(rclocal/=ESMF_SUCCESS) return
	  call ESMF_LogWrite("Coupling Initialize Method 1 Done",ESMF_LOGMSG_INFO,rc=rclocal)
    
    call ESMF_CplCompSetEntryPoint(comp,ESMF_METHOD_INITIALIZE,cpl_init2,phase=2,rc=rclocal)
    if(rclocal/=ESMF_SUCCESS) return
    call ESMF_LogWrite("Coupling Initialize Method 2 Done",ESMF_LOGMSG_INFO,rc=rclocal)

    call ESMF_CplCompSetEntryPoint(comp,ESMF_METHOD_INITIALIZE,cpl_init3,phase=3,rc=rclocal)
    if(rclocal/=ESMF_SUCCESS) return
    call ESMF_LogWrite("Coupling Initialize Method 3 Done",ESMF_LOGMSG_INFO,rc=rclocal)
    
    call ESMF_CplCompSetEntryPoint(comp,ESMF_METHOD_INITIALIZE,cpl_init4,phase=4,rc=rclocal)
    if(rclocal/=ESMF_SUCCESS) return
    call ESMF_LogWrite("Coupling Initialize Method 4 Done",ESMF_LOGMSG_INFO,rc=rclocal)
    
    call ESMF_CplCompSetEntryPoint(comp,ESMF_METHOD_RUN,cpl_run1,phase=1,rc=rclocal)
    if(rclocal/=ESMF_SUCCESS) return
	  call ESMF_LogWrite("Coupling Run Method 1 Done",ESMF_LOGMSG_INFO,rc=rclocal)

    call ESMF_CplCompSetEntryPoint(comp,ESMF_METHOD_RUN,cpl_run2,phase=2,rc=rclocal)
    if(rclocal/=ESMF_SUCCESS) return
    call ESMF_LogWrite("Coupling Run Method 2 Done",ESMF_LOGMSG_INFO,rc=rclocal)
    
    call ESMF_CplCompSetEntryPoint(comp,ESMF_METHOD_RUN,cpl_run3,phase=3,rc=rclocal)
    if(rclocal/=ESMF_SUCCESS) return
    call ESMF_LogWrite("Coupling Run Method 3 Done",ESMF_LOGMSG_INFO,rc=rclocal)

    call ESMF_CplCompSetEntryPoint(comp,ESMF_METHOD_RUN,cpl_run4,phase=4,rc=rclocal)
    if(rclocal/=ESMF_SUCCESS) return
    call ESMF_LogWrite("Coupling Run Method 4 Done",ESMF_LOGMSG_INFO,rc=rclocal)
    
    call ESMF_CplCompSetEntryPoint(comp,ESMF_METHOD_FINALIZE,cpl_final1,phase=1,rc=rclocal)
    if(rclocal/=ESMF_SUCCESS) return
	  call ESMF_LogWrite("Coupling Finalize Method 1 Done",ESMF_LOGMSG_INFO,rc=rclocal)

    call ESMF_CplCompSetEntryPoint(comp,ESMF_METHOD_FINALIZE,cpl_final2,phase=2,rc=rclocal)
    if(rclocal/=ESMF_SUCCESS) return
    call ESMF_LogWrite("Coupling Finalize Method 2 Done",ESMF_LOGMSG_INFO,rc=rclocal)

    call ESMF_CplCompSetEntryPoint(comp,ESMF_METHOD_FINALIZE,cpl_final3,phase=3,rc=rclocal)
    if(rclocal/=ESMF_SUCCESS) return
    call ESMF_LogWrite("Coupling Finalize Method 3 Done",ESMF_LOGMSG_INFO,rc=rclocal)

    call ESMF_CplCompSetEntryPoint(comp,ESMF_METHOD_FINALIZE,cpl_final4,phase=4,rc=rclocal)
    if(rclocal/=ESMF_SUCCESS) return
    call ESMF_LogWrite("Coupling Finalize Method 4 Done",ESMF_LOGMSG_INFO,rc=rclocal)

	  call ESMF_LogWrite("Coupling Registration Completed",ESMF_LOGMSG_INFO,rc=rclocal)

    rclocal=ESMF_SUCCESS

  end subroutine cpl_register
  ! =====================================================================================
  subroutine cpl_init1(comp,importState,exportState,clock,rclocal)

      type(ESMF_CplComp)::comp
      type(ESMF_State)::importState,exportState
      type(ESMF_Clock)::clock
      integer,intent(out)::rclocal

      ! Local variables
      type(ESMF_VM)::vm
      type(ESMF_Field)::oceanH,spmH,oceanU,oceanV,spmU,spmV
      integer::localPET,petCount

      rclocal=ESMF_SUCCESS

      ! Need to reconcile import and export states
      call ESMF_CplCompGet(comp,vm=vm,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return   
	    call ESMF_LogWrite("Coupling Initialize 1: Get Component",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_VMGet(vm,localPet=localPet,petCount=petCount,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
	    call ESMF_LogWrite("Coupling Initialize 1: Get VM",ESMF_LOGMSG_INFO,rc=rclocal)
      
      call ESMF_StateReconcile(importState,vm=vm,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
	    call ESMF_LogWrite("Coupling Initialize 1: Reconcile Import State",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_StateReconcile(exportState,vm=vm,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
	    call ESMF_LogWrite("Coupling Initialize 1: Reconcile Export State",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Get input data
      call ESMF_StateGet(importState,"oceanH",oceanH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
	    call ESMF_LogWrite("Coupling Initialize 1: Get Import State oceanH",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_StateGet(importState,"oceanU",oceanU,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 1: Get Import State oceanU",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_StateGet(importState,"oceanV",oceanV,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 1: Get Import State oceanV",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Get location of output data
      call ESMF_StateGet(exportState,"spmH",spmH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
	    call ESMF_LogWrite("Coupling Initialize 1: Get Export State Interpolate spmH",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_StateGet(exportState,"spmU",spmU,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 1: Get Export State Interpolate spmU",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_StateGet(exportState,"spmV",spmV,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 1: Get Export State Interpolate spmV",ESMF_LOGMSG_INFO,rc=rclocal)

      ! These are fields on different Grids - call RegridStore to set
      ! up the Regrid structure
      ! ESMF_REGRIDMETHOD_PATCH
      call ESMF_FieldRegridStore(srcField=oceanH,dstField=spmH, &
      	routeHandle=ocean2spm,regridmethod=ESMF_REGRIDMETHOD_BILINEAR,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
	    call ESMF_LogWrite("Coupling Initialize 1: Regrid Field oceanH to spmH",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_FieldRegridStore(srcField=oceanU,dstField=spmU, &
        routeHandle=ocean2spm,regridmethod=ESMF_REGRIDMETHOD_BILINEAR,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 1: Regrid Field oceanU to spmU",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_FieldRegridStore(srcField=oceanV,dstField=spmV, &
        routeHandle=ocean2spm,regridmethod=ESMF_REGRIDMETHOD_BILINEAR,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 1: Regrid Field oceanV to spmV",ESMF_LOGMSG_INFO,rc=rclocal)

      rclocal=ESMF_SUCCESS

  end subroutine cpl_init1
  ! =====================================================================================
  subroutine cpl_init2(comp,importState,exportState,clock,rclocal)

      type(ESMF_CplComp)::comp
      type(ESMF_State)::importState,exportState
      type(ESMF_Clock)::clock
      integer,intent(out)::rclocal

      ! Local variables
      type(ESMF_VM)::vm
      type(ESMF_Field)::newoceanH,newspmH
      integer::localPET,petCount

      rclocal=ESMF_SUCCESS

      ! Need to reconcile import and export states
      call ESMF_CplCompGet(comp,vm=vm,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return   
      call ESMF_LogWrite("Coupling Initialize 2: Get Component",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_VMGet(vm,localPet=localPet,petCount=petCount,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 2: Get VM",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_StateReconcile(importState,vm=vm,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 2: Reconcile Import State",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_StateReconcile(exportState,vm=vm,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 2: Reconcile Export State",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Get input data
      call ESMF_StateGet(importState,"newspmH",newspmH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 2: Get Import State spmH",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Get location of output data
      call ESMF_StateGet(exportState,"newoceanH",newoceanH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 2: Get Export State oceanH",ESMF_LOGMSG_INFO,rc=rclocal)

      ! These are fields on different Grids - call RegridStore to set
      ! up the Regrid structure
      ! ESMF_REGRIDMETHOD_PATCH
      call ESMF_FieldRegridStore(srcField=newspmH,dstField=newoceanH, &
        routeHandle=spm2ocean,regridmethod=ESMF_REGRIDMETHOD_BILINEAR,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 2: Regrid Field spmH to oceanH",ESMF_LOGMSG_INFO,rc=rclocal)

  end subroutine cpl_init2
  ! =====================================================================================
  subroutine cpl_init3(comp,importState,exportState,clock,rclocal)

      type(ESMF_CplComp)::comp
      type(ESMF_State)::importState,exportState
      type(ESMF_Clock)::clock
      integer,intent(out)::rclocal

      ! Local variables
      type(ESMF_VM)::vm
      type(ESMF_Field)::earthvD,earthrainH,spmvD,spmrainH
      integer::localPET,petCount

      rclocal=ESMF_SUCCESS

      ! Need to reconcile import and export states
      call ESMF_CplCompGet(comp,vm=vm,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return   
      call ESMF_LogWrite("Coupling Initialize 3: Get Component",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_VMGet(vm,localPet=localPet,petCount=petCount,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 3: Get VM",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_StateReconcile(importState,vm=vm,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 3: Reconcile Import State",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_StateReconcile(exportState,vm=vm,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 3: Reconcile Export State",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Get input data
      call ESMF_StateGet(importState,"earthvD",earthvD,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 3: Get Import State earthvD",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Get input data
      call ESMF_StateGet(importState,"earthrainH",earthrainH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 3: Get Import State earthrainH",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Get location of output data
      call ESMF_StateGet(exportState,"spmvD",spmvD,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 3: Get Export State spmvD",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Get location of output data
      call ESMF_StateGet(exportState,"spmrainH",spmrainH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 3: Get Export State spmrainH",ESMF_LOGMSG_INFO,rc=rclocal)

      ! These are fields on different Grids - call RegridStore to set
      ! up the Regrid structure
      ! ESMF_REGRIDMETHOD_PATCH
      call ESMF_FieldRegridStore(srcField=earthvD,dstField=spmvD, &
        routeHandle=earth2spm,regridmethod=ESMF_REGRIDMETHOD_BILINEAR,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 3: Regrid Field earthvD to spmvD",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_FieldRegridStore(srcField=earthrainH,dstField=spmrainH, &
        routeHandle=earth2spm,regridmethod=ESMF_REGRIDMETHOD_BILINEAR,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 3: Regrid Field earthrainH to spmrainH",ESMF_LOGMSG_INFO,rc=rclocal)

  end subroutine cpl_init3
  ! =====================================================================================
  subroutine cpl_init4(comp,importState,exportState,clock,rclocal)

      type(ESMF_CplComp)::comp
      type(ESMF_State)::importState,exportState
      type(ESMF_Clock)::clock
      integer,intent(out)::rclocal

      ! Local variables
      type(ESMF_VM)::vm
      type(ESMF_Field)::newearthH,newspmH
      integer::localPET,petCount

      rclocal=ESMF_SUCCESS

      ! Need to reconcile import and export states
      call ESMF_CplCompGet(comp,vm=vm,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return   
      call ESMF_LogWrite("Coupling Initialize 4: Get Component",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_VMGet(vm,localPet=localPet,petCount=petCount,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 4: Get VM",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_StateReconcile(importState,vm=vm,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 4: Reconcile Import State",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_StateReconcile(exportState,vm=vm,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 4: Reconcile Export State",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Get input data
      call ESMF_StateGet(importState,"newspmH",newspmH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 4: Get Import State spmH",ESMF_LOGMSG_INFO,rc=rclocal)


      ! Get location of output data
      call ESMF_StateGet(exportState,"newearthH",newearthH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 4: Get Export State newearthH",ESMF_LOGMSG_INFO,rc=rclocal)

      ! These are fields on different Grids - call RegridStore to set
      ! up the Regrid structure
      ! ESMF_REGRIDMETHOD_PATCH
      call ESMF_FieldRegridStore(srcField=newspmH,dstField=newearthH, &
        routeHandle=spm2earth,regridmethod=ESMF_REGRIDMETHOD_BILINEAR,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Initialize 4: Regrid Field spmH to earthH",ESMF_LOGMSG_INFO,rc=rclocal)

  end subroutine cpl_init4
  ! =====================================================================================
  subroutine cpl_run1(comp,importState,exportState,clock,rclocal)

      type(ESMF_CplComp)::comp
      type(ESMF_State)::importState,exportState
      type(ESMF_Clock)::clock
      integer,intent(out)::rclocal

      ! Local variables
      type(ESMF_Field)::oceanH,spmH,oceanU,oceanV,spmU,spmV
      integer::status

      rclocal=ESMF_SUCCESS

      ! Get input data
      if(.not.updateSPM_elevation)then
      
        call ESMF_StateGet(importState,"oceanH",oceanH,rc=rclocal)
        if(rclocal/=ESMF_SUCCESS) return
	      call ESMF_LogWrite("Coupling Run 1: Get Import State oceanH",ESMF_LOGMSG_INFO,rc=rclocal)
        
        ! Get location of output data
        call ESMF_StateGet(exportState,"spmH",spmH,rc=rclocal)
        if(rclocal/=ESMF_SUCCESS) return
        call ESMF_LogWrite("Coupling Run 1: Get Export State Interpolate spmH",ESMF_LOGMSG_INFO,rc=rclocal)

        call ESMF_FieldRegrid(oceanH,spmH,ocean2spm,rc=status)
        if(rclocal/=ESMF_SUCCESS) return
        call ESMF_LogWrite("Coupling Run 1: Regrid Field oceanH to spmH",ESMF_LOGMSG_INFO,rc=rclocal)

      endif
      
      call ESMF_StateGet(importState,"oceanU",oceanU,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Run 1: Get Import State oceanU",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_StateGet(importState,"oceanV",oceanV,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Run 1: Get Import State oceanV",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Get location of output data
      call ESMF_StateGet(exportState,"spmU",spmU,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Run 1: Get Export State Interpolate spmU",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_StateGet(exportState,"spmV",spmV,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Run 1: Get Export State Interpolate spmV",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_FieldRegrid(oceanU,spmU,ocean2spm,rc=status)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Run 1: Regrid Field oceanU to spmU",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_FieldRegrid(oceanV,spmV,ocean2spm,rc=status)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Run 1: Regrid Field oceanV to spmV",ESMF_LOGMSG_INFO,rc=rclocal)

      rclocal=ESMF_SUCCESS

  end subroutine cpl_run1
  ! =====================================================================================
  subroutine cpl_run2(comp,importState,exportState,clock,rclocal)

      type(ESMF_CplComp)::comp
      type(ESMF_State)::importState,exportState
      type(ESMF_Clock)::clock
      integer,intent(out)::rclocal

      ! Local variables
      type(ESMF_Field)::newoceanH,newspmH
      integer::status

      rclocal=ESMF_SUCCESS

      ! Get input data
      call ESMF_StateGet(importState,"newspmH",newspmH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Run 2: Get Import State spmH",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Get location of output data
      call ESMF_StateGet(exportState,"newoceanH",newoceanH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Run 2: Get Export State oceanH",ESMF_LOGMSG_INFO,rc=rclocal)

      ! These are fields on different Grids - call Regrid to rearrange
      ! the data. The communication pattern was computed at init,
      ! this simply has to execute the send and receive equivalents.
      call ESMF_FieldRegrid(newspmH,newoceanH,spm2ocean,rc=status)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Run 2: Regrid Field spmH to oceanH",ESMF_LOGMSG_INFO,rc=rclocal)

      rclocal=ESMF_SUCCESS

  end subroutine cpl_run2
  ! =====================================================================================
  subroutine cpl_run3(comp,importState,exportState,clock,rclocal)

      type(ESMF_CplComp)::comp
      type(ESMF_State)::importState,exportState
      type(ESMF_Clock)::clock
      integer,intent(out)::rclocal

      ! Local variables
      type(ESMF_Field)::earthvD,earthrainH,spmvD,spmrainH
      integer::status

      rclocal=ESMF_SUCCESS
      
      call ESMF_StateGet(importState,"earthvD",earthvD,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Run 3: Get Import State earthvD",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_StateGet(importState,"earthrainH",earthrainH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Run 3: Get Import State earthrainH",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Get location of output data
      call ESMF_StateGet(exportState,"spmvD",spmvD,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Run 3: Get Export State Interpolate spmvD",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_FieldRegrid(earthvD,spmvD,earth2spm,rc=status)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Run 3: Regrid Field earthvD to spmvD",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_StateGet(exportState,"spmrainH",spmrainH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Run 3: Get Export State Interpolate spmrainH",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_FieldRegrid(earthrainH,spmrainH,earth2spm,rc=status)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Run 3: Regrid Field earthrainH to spmrainH",ESMF_LOGMSG_INFO,rc=rclocal)

      rclocal=ESMF_SUCCESS

  end subroutine cpl_run3
  ! =====================================================================================
  subroutine cpl_run4(comp,importState,exportState,clock,rclocal)

      type(ESMF_CplComp)::comp
      type(ESMF_State)::importState,exportState
      type(ESMF_Clock)::clock
      integer,intent(out)::rclocal

      ! Local variables
      type(ESMF_Field)::newearthH,newspmH
      integer::status

      rclocal=ESMF_SUCCESS

      ! Get input data
      call ESMF_StateGet(importState,"newspmH",newspmH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Run 4: Get Import State spmH",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Get location of output data
      call ESMF_StateGet(exportState,"newearthH",newearthH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Run 4: Get Export State newearthH",ESMF_LOGMSG_INFO,rc=rclocal)

      ! These are fields on different Grids - call Regrid to rearrange
      ! the data. The communication pattern was computed at init,
      ! this simply has to execute the send and receive equivalents.
      call ESMF_FieldRegrid(newspmH,newearthH,spm2earth,rc=status)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Run 4: Regrid Field spmH to earthH",ESMF_LOGMSG_INFO,rc=rclocal)

      rclocal=ESMF_SUCCESS

  end subroutine cpl_run4
  ! =====================================================================================
  subroutine cpl_final1(comp,importState,exportState,clock,rclocal)

      type(ESMF_CplComp)::comp
      type(ESMF_State)::importState,exportState
      type(ESMF_Clock)::clock
      integer,intent(out)::rclocal

      rclocal=ESMF_SUCCESS

      ! Release resources stored for the Regridding.
      call ESMF_FieldRegridRelease(ocean2spm,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Finalize: Regrid Field Release Route 1",ESMF_LOGMSG_INFO,rc=rclocal)
      
      ! Return success
      rclocal=ESMF_SUCCESS 

  end subroutine cpl_final1
  ! =====================================================================================
  subroutine cpl_final2(comp,importState,exportState,clock,rclocal)

      type(ESMF_CplComp)::comp
      type(ESMF_State)::importState,exportState
      type(ESMF_Clock)::clock
      integer,intent(out)::rclocal

      rclocal=ESMF_SUCCESS

      ! Release resources stored for the Regridding.
      call ESMF_FieldRegridRelease(spm2ocean,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Finalize: Regrid Field Release Route 2",ESMF_LOGMSG_INFO,rc=rclocal)
      
      ! Return success
      rclocal=ESMF_SUCCESS 

  end subroutine cpl_final2
  ! =====================================================================================
  subroutine cpl_final3(comp,importState,exportState,clock,rclocal)

      type(ESMF_CplComp)::comp
      type(ESMF_State)::importState,exportState
      type(ESMF_Clock)::clock
      integer,intent(out)::rclocal

      rclocal=ESMF_SUCCESS

      ! Release resources stored for the Regridding.
      call ESMF_FieldRegridRelease(earth2spm,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Finalize: Regrid Field Release Route 3",ESMF_LOGMSG_INFO,rc=rclocal)
      
      ! Return success
      rclocal=ESMF_SUCCESS 

  end subroutine cpl_final3
  ! =====================================================================================
  subroutine cpl_final4(comp,importState,exportState,clock,rclocal)

      type(ESMF_CplComp)::comp
      type(ESMF_State)::importState,exportState
      type(ESMF_Clock)::clock
      integer,intent(out)::rclocal

      rclocal=ESMF_SUCCESS

      ! Release resources stored for the Regridding.
      call ESMF_FieldRegridRelease(spm2earth,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Coupling Finalize: Regrid Field Release Route 4",ESMF_LOGMSG_INFO,rc=rclocal)
      
      ! Return success
      rclocal=ESMF_SUCCESS 

  end subroutine cpl_final4
  ! =====================================================================================
 
end module esmfCoupler