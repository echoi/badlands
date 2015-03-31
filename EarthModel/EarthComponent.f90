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
!       Filename:  EarthComponent.f90
!
!    Description:  Defines the Earth Geodyanmic ESMF Component
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles     
!
! =====================================================================================

module esmfEarth

  use topology
  use earthforces
  use underworld
  use parameters

  implicit none

  public earth_register

contains

  ! =====================================================================================
  subroutine earth_register(comp,rclocal)

    type(ESMF_GridComp)::comp
    integer,intent(out)::rclocal

    rclocal=ESMF_SUCCESS

    ! Register the callback routines
    call ESMF_GridCompSetEntryPoint(comp,ESMF_METHOD_INITIALIZE,earth_init,rc=rclocal)
    if(rclocal/=ESMF_SUCCESS) return
    call ESMF_LogWrite("SPM Initialize Method Done",ESMF_LOGMSG_INFO,rc=rclocal)

    call ESMF_GridCompSetEntryPoint(comp,ESMF_METHOD_RUN,earth_run,rc=rclocal)
    if(rclocal/=ESMF_SUCCESS) return
    call ESMF_LogWrite("SPM Run Method Done",ESMF_LOGMSG_INFO,rc=rclocal)

    call ESMF_GridCompSetEntryPoint(comp,ESMF_METHOD_FINALIZE,earth_final,rc=rclocal)
    if(rclocal/=ESMF_SUCCESS) return
    call ESMF_LogWrite("SPM Finalize Method Done",ESMF_LOGMSG_INFO,rc=rclocal)

    call ESMF_LogWrite("SPM Model Component Registration Completed",ESMF_LOGMSG_INFO,rc=rclocal)

  end subroutine earth_register
  ! =====================================================================================
  subroutine earth_init(comp,importState,exportState,clock,rclocal)

      type(ESMF_GridComp)::comp
      type(ESMF_State)::importState,exportState
      type(ESMF_Clock)::clock
      integer,intent(out)::rclocal

      ! Local variables
      type(ESMF_VM)::vm
      type(ESMF_Field)::earthvD,newearthH,earthrainH
      type(ESMF_Mesh)::earthMesh
      type(ESMF_ArraySpec)::arrayspec
      integer::localPET,petCount,localrc,i,k,p
      real(ESMF_KIND_R8),pointer::fptr1D(:),fptr1D2(:),fptr1D3(:)
      integer,pointer::nodeIds(:),nodeOwners(:)
      real(ESMF_KIND_R8),pointer::nodeCoords(:)
      integer,pointer::elemIds(:),elemTypes(:),elemConn(:)
      integer::numNodes,numElems
      integer::numOwnedNodes

      rclocal=ESMF_SUCCESS

      ! Query component for VM and create a layout with the right breakdown
      call ESMF_GridCompGet(comp,vm=vm,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Earth Initialize: Get Component",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_VMGet(vm,localPet=localPet,petCount=petCount,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Earth Initialize: Get VM",ESMF_LOGMSG_INFO,rc=rclocal)

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!     Build Earth mesh
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

      ! Set number of nodes for the considered PET
      numNodes=spartN

      ! Allocate and fill the node id array
      ! Allocate and fill node coordinate array.
      ! Since this is a 2D Mesh the size is 2x the
      ! number of nodes.
      allocate(nodeIds(numNodes),nodeCoords(2*numNodes))
      if(allocated(earthIDs)) deallocate(earthIDs)
      allocate(earthIDs(sOwnedNode))
      allocate(nodeOwners(numNodes))

      k=1
      p=1
      do i=1,numNodes 
        nodeIds(i)=snodeId(i)
        nodeCoords(k)=rcoordX(nodeIds(i))
        k=k+1
        nodeCoords(k)=rcoordY(nodeIds(i))
        k=k+1
        nodeOwners(i)=sownedID(i)
        if(nodeOwners(i)==localPet)then 
          earthIDs(p)=snodeID(i)
          p=p+1
        endif
      enddo

      ! Set the number of each type of element, plus the total number.
      numElems=spartE

      ! Allocate and fill the element id array.
      ! Allocate and fill the element topology type array.
      ! Allocate and fill the element connection type array.
      allocate(elemIds(numElems))
      allocate(elemTypes(numElems))
      allocate(elemConn(4*numElems))
      k=1
      do i=1,numElems 
        elemIds(i)=selemID(i)
        elemTypes(i)=ESMF_MESHELEMTYPE_QUAD
        elemConn(k)=snodeLID(relmt(elemIds(i),1))
        k=k+1
        elemConn(k)=snodeLID(relmt(elemIds(i),2))
        k=k+1
        elemConn(k)=snodeLID(relmt(elemIds(i),3))
        k=k+1
        elemConn(k)=snodeLID(relmt(elemIds(i),4))
        k=k+1
      enddo

      ! Create Mesh structure in 1 step
      earthMesh=ESMF_MeshCreate(parametricDim=2,spatialDim=2,nodeIds=nodeIds, &
        nodeCoords=nodeCoords,nodeOwners=nodeOwners,elementIds=elemIds,&
        elementTypes=elemTypes,elementConn=elemConn,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("Earth Initialize: Create Mesh",ESMF_LOGMSG_INFO,rc=rclocal)

      ! deallocate node data
      deallocate(nodeIds)
      deallocate(nodeCoords)
      deallocate(nodeOwners)

      ! deallocate elem data
      deallocate(elemIds)
      deallocate(elemTypes)
      deallocate(elemConn)

      ! Create dem topography field
      call ESMF_ArraySpecSet(arrayspec,1,ESMF_TYPEKIND_R8,rc=rclocal)

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!     Export states
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

      ! Get rainfall and displacement
      if(rain_event>0) call rainfall
      if(disp%event>0.and..not.udwFlag)then 
        call displacement
      elseif(udwFlag)then
        call SurfaceVTK
        call WaitStepCompletion
        call displacement
      endif

      ! Define coupling time for rain
      if(rain_event==0) cpl1_time=time_end+1000.

      ! Define coupling time for geodynamic
      if(disp%event==0) cpl2_time=time_end+1000.

      ! Create Export Fields
      earthvD=ESMF_FieldCreate(earthMesh,arrayspec,name="earthvD",rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("Earth Initialize: Create Export Field Vertical Displacement",ESMF_LOGMSG_INFO,rc=rclocal)
      
      ! Load data into the source Field
      call ESMF_FieldGet(earthvD,farrayPtr=fptr1D,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
          rclocal=ESMF_FAILURE
          return
      endif      
      call ESMF_LogWrite("Earth Initialize: Get Export Field Vertical Displacement",ESMF_LOGMSG_INFO,rc=rclocal)
      
      ! Create Export Fields
      earthrainH=ESMF_FieldCreate(earthMesh,arrayspec,name="earthrainH",rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("Earth Initialize: Create Export Field Rain",ESMF_LOGMSG_INFO,rc=rclocal)
      
      ! Load data into the source Field
      call ESMF_FieldGet(earthrainH,farrayPtr=fptr1D3,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
          rclocal=ESMF_FAILURE
          return
      endif      
      call ESMF_LogWrite("Earth Initialize: Get Export Field Rain",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Set vertical displacement field
      call ESMF_MeshGet(earthMesh,numOwnedNodes=numOwnedNodes,rc=rclocal)  
      do i=1,numOwnedNodes
        if(disp%event>0)then 
          fptr1D(i)=0.0
        else
          fptr1D(i)=rvertDisp(earthIDs(i))
        endif
        if(rain_event>0)then
          fptr1D3(i)=rainVal(earthIDs(i))
        else
          fptr1D3(i)=1.0
        endif
      enddo
      
      ! Add Fields to State
      call ESMF_StateAdd(exportState,(/earthvD,earthrainH/),rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Earth Initialize: Add Export States",ESMF_LOGMSG_INFO,rc=rclocal)

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!     Import states
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

      newearthH=ESMF_FieldCreate(earthMesh,arrayspec,name="newearthH",rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("Earth Initialize: Create Import Field newearthH",ESMF_LOGMSG_INFO,rc=rclocal)
      
      call ESMF_FieldGet(newearthH,farrayPtr=fptr1D2,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
          rclocal=ESMF_FAILURE
          return
      endif      
      call ESMF_LogWrite("Earth Initialize: Get Import Field newearthH",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Clear interpolated topography field
      fptr1D2=0.0

      ! Add Field to State
      call ESMF_StateAdd(importState,(/newearthH/),rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Earth Initialize: Add Import State",ESMF_LOGMSG_INFO,rc=rclocal)

      rclocal=ESMF_SUCCESS

  end subroutine earth_init
  ! =====================================================================================
  subroutine earth_run(comp,importState,exportState,clock,rclocal)

      type(ESMF_GridComp)::comp
      type(ESMF_State)::importState,exportState
      type(ESMF_Clock)::clock
      integer,intent(out)::rclocal
      real(ESMF_KIND_R8)::time1,time2
      real(ESMF_KIND_R8),dimension(:),allocatable::newelev

      ! Local variables
      type(ESMF_VM)::vm
      type(ESMF_Field)::earthvD,newearthH,earthrainH
      type(ESMF_Mesh)::earthMesh
      integer::localrc,i
      integer::localPet,petCount
      real(ESMF_KIND_R8),pointer::fptr1D(:),fptr1D2(:),fptr1D3(:)
      integer::numOwnedNodes
      real(ESMF_KIND_R8),pointer::ownedNodeCoords(:)

      rclocal=ESMF_SUCCESS

      call ESMF_VMWtime(time1,rc=rc)
      ! Query component for VM and create a layout with the right breakdown
      call ESMF_GridCompGet(comp,vm=vm,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Run: Get Grid Component",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_VMGet(vm,localPet=localPet,petCount=petCount,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Run: Get VM",ESMF_LOGMSG_INFO,rc=rclocal)

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!     Import states
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

      if(updateSPM_elevation)then

        ! Get information from the component.
        call ESMF_StateGet(importState,"newearthH",newearthH,rc=rclocal)
        if(rclocal/=ESMF_SUCCESS) return
        call ESMF_LogWrite("Earth Run: Get Import State",ESMF_LOGMSG_INFO,rc=rclocal)

        ! Get Grid from field
        call ESMF_FieldGet(newearthH,mesh=earthMesh,rc=localrc)
        if(localrc/=ESMF_SUCCESS)then
          rclocal=ESMF_FAILURE
          return
        endif
        call ESMF_LogWrite("Earth Run: Get Grid from Field",ESMF_LOGMSG_INFO,rc=rclocal)

        ! Get number of local nodes to allocate space
        ! to hold local node coords
        call ESMF_MeshGet(earthMesh,numOwnedNodes=numOwnedNodes,rc=rclocal)  
        if(rclocal/=ESMF_SUCCESS) return
        call ESMF_LogWrite("Earth Run: Get Structured Mesh",ESMF_LOGMSG_INFO,rc=rclocal)    

        ! Allocate space to hold local node coordinates
        ! (spatial dimension of Mesh*number of local nodes)
        allocate(ownedNodeCoords(2*numOwnedNodes))

        ! Get local node coordinates
        call ESMF_MeshGet(earthMesh,ownedNodeCoords=ownedNodeCoords,rc=rclocal)
        if(rclocal/=ESMF_SUCCESS) return
        call ESMF_LogWrite("Earth Run: Get Mesh Coordinates",ESMF_LOGMSG_INFO,rc=rclocal)   

        ! Check destination field
        ! Should only be 1 localDE
        call ESMF_FieldGet(newearthH,farrayPtr=fptr1D2,rc=localrc)
        if(localrc/=ESMF_SUCCESS)then
          rclocal=ESMF_FAILURE
          return
        endif
        call ESMF_LogWrite("Earth Run: Get Field",ESMF_LOGMSG_INFO,rc=rclocal)

        ! Loop through nodes and defined the interpolated topographic value
        allocate(newelev(bnbnodes))
        newelev=-1.e6
        do i=1,numOwnedNodes
          newelev(earthIDs(i))=fptr1D2(i)
        enddo
        call ESMF_VMAllReduce(vm=vm,sendData=newelev,recvData=rtectoZ,count=bnbnodes,&
          reduceflag=ESMF_REDUCE_MAX,rc=rc)
        deallocate(newelev)
        
        ! Reset destination to 0
        fptr1D2=0.0

        rclocal=ESMF_SUCCESS

      endif
      
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!     Export states
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

      ! This is where we compute Earth geodynamic field and rain climate
      if(updateSPM_elevation)then
        if(rain_event>0.and.cpl1_time<=simulation_time) call rainfall
        if(disp%event>0.and.cpl2_time<=simulation_time)then 
          if(udwFlag)then
            call SurfaceVTK
            call WaitStepCompletion
            call displacement
          else
            call displacement
          endif
        endif
      endif
      call ESMF_VMWtime(time2,rc=rc)
      if(pet_id==0) print*,'BADLANDS Displacement & Rainfall step Completed (s) ',time2-time1
      if(pet_id==0)print*,'-------------------------'

      ! Prepare the Export State
      call ESMF_StateGet(exportState,"earthvD",earthvD,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Earth Run: Get Export State earthvD",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_FieldGet(earthvD,mesh=earthMesh,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("Earth Run: Get Grid from Export Field earthvD",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_StateGet(exportState,"earthrainH",earthrainH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Earth Run: Get Export State earthrainH",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_FieldGet(earthrainH,mesh=earthMesh,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("Earth Run: Get Grid from Export Field earthrainH",ESMF_LOGMSG_INFO,rc=rclocal)

      if(.not.updateSPM_elevation)then

        ! Get number of local nodes to allocate space
        ! to hold local node coords
        call ESMF_MeshGet(earthMesh,numOwnedNodes=numOwnedNodes,rc=rclocal)  
        if(rclocal/=ESMF_SUCCESS) return
        call ESMF_LogWrite("Earth Run: Get Structured Mesh",ESMF_LOGMSG_INFO,rc=rclocal)    

        ! Allocate space to hold local node coordinates
        ! (spatial dimension of Mesh*number of local nodes)
        allocate(ownedNodeCoords(2*numOwnedNodes))

        ! Get local node coordinates
        call ESMF_MeshGet(earthMesh,ownedNodeCoords=ownedNodeCoords,rc=rclocal)
        if(rclocal/=ESMF_SUCCESS) return
        call ESMF_LogWrite("Earth Run: Get Mesh Coordinates",ESMF_LOGMSG_INFO,rc=rclocal)   

      endif

      call ESMF_FieldGet(earthvD,farrayPtr=fptr1D,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("Earth Run: Get Export Field earthvD",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_FieldGet(earthrainH,farrayPtr=fptr1D3,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("Earth Run: Get Export Field earthrainH",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Loop through nodes and get the displacement and rain values
      do i=1,numOwnedNodes
        fptr1D(i)=rvertDisp(earthIDs(i)) 
        fptr1D3(i)=rainVal(earthIDs(i))
      enddo

      ! Deallocate space to hold local node coordinates
      deallocate(ownedNodeCoords)

      rclocal=ESMF_SUCCESS

  end subroutine earth_run
  ! =====================================================================================
  subroutine earth_final(comp,importState,exportState,clock,rclocal)

      type(ESMF_GridComp)::comp
      type(ESMF_State)::importState,exportState
      type(ESMF_Clock)::clock
      integer,intent(out)::rclocal

      ! Local variables
      type(ESMF_Field)::earthvD,newearthH,earthrainH
      type(ESMF_Mesh)::earthMesh

      rclocal=ESMF_SUCCESS

      ! Check validity of results
      ! Get Fields from import state
      call ESMF_StateGet(importState,"newearthH",newearthH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Earth Finalize: Get Import State newearthH",ESMF_LOGMSG_INFO,rc=rclocal)    

      ! Get Fields from export state
      call ESMF_StateGet(exportState,"earthvD",earthvD,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Earth Finalize: Get Export State earthvD",ESMF_LOGMSG_INFO,rc=rclocal)    

      ! Get Fields from export state
      call ESMF_StateGet(exportState,"earthrainH",earthrainH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Earth Finalize: Get Export State earthrainH",ESMF_LOGMSG_INFO,rc=rclocal)  
      
      ! Garbage collection
      call ESMF_FieldGet(newearthH,mesh=earthMesh,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Earth Finalize: Get Import Fiel newearthH",ESMF_LOGMSG_INFO,rc=rclocal)    

      call ESMF_FieldGet(earthvD,mesh=earthMesh,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Earth Finalize: Get Export Field earthvD",ESMF_LOGMSG_INFO,rc=rclocal)    

      call ESMF_FieldGet(earthrainH,mesh=earthMesh,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Earth Finalize: Get Export Field earthrainH",ESMF_LOGMSG_INFO,rc=rclocal)   

      call ESMF_FieldDestroy(newearthH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Earth Finalize: Import Field newearthH Destroyed",ESMF_LOGMSG_INFO,rc=rclocal)    

      call ESMF_FieldDestroy(earthvD,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Earth Finalize: Export Field earthvD Destroyed",ESMF_LOGMSG_INFO,rc=rclocal)    

      call ESMF_FieldDestroy(earthrainH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Earth Finalize: Export Field earthrainH Destroyed",ESMF_LOGMSG_INFO,rc=rclocal)   

      call ESMF_MeshDestroy(earthMesh,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Earth Finalize: Mesh Destroyed",ESMF_LOGMSG_INFO,rc=rclocal)    
      
      ! Return success
      rclocal=ESMF_SUCCESS 

  end subroutine earth_final
  ! =====================================================================================
 
end module esmfEarth
    
