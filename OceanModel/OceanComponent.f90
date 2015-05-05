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
!       Filename:  OceanComponent.f90
!
!    Description:  Defines the Ocean Model ESMF Component
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles     
!
! =====================================================================================

module esmfOcean

  use topology
  use parameters

  implicit none

  public ocean_register

contains

  ! =====================================================================================
  subroutine ocean_register(comp,rclocal)

    type(ESMF_GridComp)::comp
    integer,intent(out)::rclocal

    rclocal=ESMF_SUCCESS

    ! Register the callback routines
    call ESMF_GridCompSetEntryPoint(comp,ESMF_METHOD_INITIALIZE,ocean_init,rc=rclocal)
    if(rclocal/=ESMF_SUCCESS) return
    call ESMF_LogWrite("Ocean Initialize Method Done",ESMF_LOGMSG_INFO,rc=rclocal)

    call ESMF_GridCompSetEntryPoint(comp,ESMF_METHOD_RUN,ocean_run,rc=rclocal)
    if(rclocal/=ESMF_SUCCESS) return
    call ESMF_LogWrite("Ocean Run Method Done",ESMF_LOGMSG_INFO,rc=rclocal)
    
    call ESMF_GridCompSetEntryPoint(comp,ESMF_METHOD_FINALIZE,ocean_final,rc=rclocal)
    if(rclocal/=ESMF_SUCCESS) return
    call ESMF_LogWrite("Ocean Finalize Method Done",ESMF_LOGMSG_INFO,rc=rclocal)
 
    call ESMF_LogWrite("Ocean Model Component Registration Completed",ESMF_LOGMSG_INFO,rc=rclocal)

  end subroutine ocean_register
  ! =====================================================================================
  subroutine ocean_init(comp,importState,exportState,clock,rclocal)

      type(ESMF_GridComp)::comp
      type(ESMF_State)::importState,exportState
      type(ESMF_Clock)::clock
      integer,intent(out)::rclocal

      ! Local variables
      type(ESMF_VM)::vm
      type(ESMF_Field)::oceanH,newoceanH,oceanU,oceanV
      type(ESMF_Mesh)::oceanMesh
      type(ESMF_ArraySpec)::arrayspec
      integer::localPET,petCount,localrc,i,k,p
      real(ESMF_KIND_R8),pointer::fptr1D(:),fptr1D2(:),fptr1D3(:),fptr1D4(:)
      integer,pointer::nodeIds(:),nodeOwners(:)
      real(ESMF_KIND_R8),pointer::nodeCoords(:)
      integer,pointer::elemIds(:),elemTypes(:),elemConn(:)
      integer::numNodes,numElems
      integer::numOwnedNodes

      rclocal=ESMF_SUCCESS

      ! Query component for VM and create a layout with the right breakdown
      call ESMF_GridCompGet(comp,vm=vm,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Ocean Initialize: Get Component",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_VMGet(vm,localPet=localPet,petCount=petCount,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Ocean Initialize: Get VM",ESMF_LOGMSG_INFO,rc=rclocal)
     
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!     Build Ocean mesh
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

      ! Set number of nodes for the considered PET
      numNodes=spartN

      ! Allocate and fill the node id array
      ! Allocate and fill node coordinate array.
      ! Since this is a 2D Mesh the size is 2x the
      ! number of nodes.
      allocate(nodeIds(numNodes),nodeCoords(2*numNodes))
      allocate(nodeOwners(numNodes))
      if(allocated(oceanIDs)) deallocate(oceanIDs)
      allocate(oceanIDs(sOwnedNode))

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
          oceanIDs(p)=snodeID(i)
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
      oceanMesh=ESMF_MeshCreate(parametricDim=2,spatialDim=2,nodeIds=nodeIds, &
        nodeCoords=nodeCoords,nodeOwners=nodeOwners,elementIds=elemIds,&
        elementTypes=elemTypes,elementConn=elemConn,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("Ocean Initialize: Create Mesh",ESMF_LOGMSG_INFO,rc=rclocal)
      
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

      ! Create Export Fields
      oceanH=ESMF_FieldCreate(oceanMesh,arrayspec,name="oceanH",rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("Ocean Initialize: Create Export Field oceanH",ESMF_LOGMSG_INFO,rc=rclocal)
      
      ! Load data into the source Field
      call ESMF_FieldGet(oceanH,farrayPtr=fptr1D,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
          rclocal=ESMF_FAILURE
          return
      endif      
      call ESMF_LogWrite("Ocean Initialize: Get Export Field oceanH",ESMF_LOGMSG_INFO,rc=rclocal)

      oceanU=ESMF_FieldCreate(oceanMesh,arrayspec,name="oceanU",rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("Ocean Initialize: Create Export Field oceanU",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Load data into the source Field
      call ESMF_FieldGet(oceanU,farrayPtr=fptr1D3,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
          rclocal=ESMF_FAILURE
          return
      endif      
      call ESMF_LogWrite("Ocean Initialize: Get Export Field oceanU",ESMF_LOGMSG_INFO,rc=rclocal)

      oceanV=ESMF_FieldCreate(oceanMesh,arrayspec,name="oceanV",rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("Ocean Initialize: Create Export Field oceanV",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Load data into the source Field
      call ESMF_FieldGet(oceanV,farrayPtr=fptr1D4,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
          rclocal=ESMF_FAILURE
          return
      endif      
      call ESMF_LogWrite("Ocean Initialize: Get Export Field oceanV",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Set ocean velocities field
      call ESMF_MeshGet(oceanMesh,numOwnedNodes=numOwnedNodes,rc=rclocal)
      do i=1,numOwnedNodes
        fptr1D(i)=rcoordZ(oceanIDs(i))
      enddo
      fptr1D3=0.0
      fptr1D4=0.0

      ! Add Export Field to State
      call ESMF_StateAdd(exportState,(/oceanH,oceanU,oceanV/),rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Ocean Initialize: Add Export States",ESMF_LOGMSG_INFO,rc=rclocal)

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!     Import states
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

      newoceanH=ESMF_FieldCreate(oceanMesh,arrayspec,name="newoceanH",rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("Ocean Initialize: Create Import Field newoceanH",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_FieldGet(newoceanH,farrayPtr=fptr1D2,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
          rclocal=ESMF_FAILURE
          return
      endif      
      call ESMF_LogWrite("Ocean Initialize: Get Import Field newoceanH",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Clear interpolated topography field
      fptr1D2=0.0

      ! Add Field to State
      call ESMF_StateAdd(importState,(/newoceanH/),rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Ocean Initialize: Add Import State",ESMF_LOGMSG_INFO,rc=rclocal)

      rclocal=ESMF_SUCCESS

  end subroutine ocean_init
  ! =====================================================================================
  subroutine ocean_run(comp,importState,exportState,clock,rclocal)

      type(ESMF_GridComp)::comp
      type(ESMF_State)::importState,exportState
      type(ESMF_Clock)::clock
      integer,intent(out)::rclocal

      ! Local variables
      type(ESMF_VM)::vm
      type(ESMF_Field)::newoceanH,oceanU,oceanV
      type(ESMF_Mesh)::oceanMesh
      integer::localrc,i
      integer::localPet,petCount
      real(ESMF_KIND_R8),pointer::fptr1D2(:),fptr1D3(:),fptr1D4(:)
      integer::numOwnedNodes
      real(ESMF_KIND_R8),pointer::ownedNodeCoords(:)

      rclocal=ESMF_SUCCESS

      ! Query component for VM and create a layout with the right breakdown
      call ESMF_GridCompGet(comp,vm=vm,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Ocean Run: Get Grid Component",ESMF_LOGMSG_INFO,rc=rclocal)
      
      call ESMF_VMGet(vm,localPet=localPet,petCount=petCount,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Ocean Run: Get VM",ESMF_LOGMSG_INFO,rc=rclocal) 

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!     Import state
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

      if(updateSPM_elevation)then

        ! Get information from the component.
        call ESMF_StateGet(importState,"newoceanH",newoceanH,rc=rclocal)
        if(rclocal/=ESMF_SUCCESS) return
        call ESMF_LogWrite("Ocean Run: Get State",ESMF_LOGMSG_INFO,rc=rclocal)

        ! Get Grid from field
        call ESMF_FieldGet(newoceanH,mesh=oceanMesh,rc=localrc)
        if(localrc/=ESMF_SUCCESS)then
          rclocal=ESMF_FAILURE
          return
        endif
        call ESMF_LogWrite("Ocean Run: Get Grid from Field",ESMF_LOGMSG_INFO,rc=rclocal)

        ! Get number of local nodes to allocate space
        ! to hold local node coords
        call ESMF_MeshGet(oceanMesh,numOwnedNodes=numOwnedNodes,rc=rclocal)  
        if(rclocal/=ESMF_SUCCESS) return
        call ESMF_LogWrite("Ocean Run: Get Structured Mesh",ESMF_LOGMSG_INFO,rc=rclocal)    

        ! Allocate space to hold local node coordinates
        ! (spatial dimension of Mesh*number of local nodes)
        allocate(ownedNodeCoords(2*numOwnedNodes))

        ! Get local node coordinates
        call ESMF_MeshGet(oceanMesh,ownedNodeCoords=ownedNodeCoords,rc=rclocal)
        if(rclocal/=ESMF_SUCCESS) return
        call ESMF_LogWrite("Ocean Run: Get Mesh Coordinates",ESMF_LOGMSG_INFO,rc=rclocal)   

        ! Check destination field
        ! Should only be 1 localDE
        call ESMF_FieldGet(newoceanH,farrayPtr=fptr1D2,rc=localrc)
        if(localrc/=ESMF_SUCCESS)then
          rclocal=ESMF_FAILURE
          return
        endif
        call ESMF_LogWrite("Ocean Run: Get Field",ESMF_LOGMSG_INFO,rc=rclocal)

        ! Loop through nodes and defined the interpolated topographic value
        do i=1,numOwnedNodes
          rcoordZ(oceanIDs(i))=fptr1D2(i)
!           if(i==2)print*,'frfr',rcoordX(i),rcoordY(i),rcoordZ(i)
        enddo

        ! Reset destination to 0
        fptr1D2=0.0

        rclocal=ESMF_SUCCESS

      endif

      ! This is where we compute Ocean velocity field
      

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!     Export states
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

      ! Prepare the Export State
      call ESMF_StateGet(exportState,"oceanU",oceanU,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Ocean Run: Get Export State oceanU",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_StateGet(exportState,"oceanV",oceanV,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Ocean Run: Get Export State oceanV",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_FieldGet(oceanU,mesh=oceanMesh,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("Ocean Run: Get Grid from Export Field oceanU",ESMF_LOGMSG_INFO,rc=rclocal)

      if(.not.updateSPM_elevation)then

        ! Get number of local nodes to allocate space
        ! to hold local node coords
        call ESMF_MeshGet(oceanMesh,numOwnedNodes=numOwnedNodes,rc=rclocal)  
        if(rclocal/=ESMF_SUCCESS) return
        call ESMF_LogWrite("Ocean Run: Get Structured Mesh",ESMF_LOGMSG_INFO,rc=rclocal)    

        ! Allocate space to hold local node coordinates
        ! (spatial dimension of Mesh*number of local nodes)
        allocate(ownedNodeCoords(2*numOwnedNodes))

        ! Get local node coordinates
        call ESMF_MeshGet(oceanMesh,ownedNodeCoords=ownedNodeCoords,rc=rclocal)
        if(rclocal/=ESMF_SUCCESS) return
        call ESMF_LogWrite("Ocean Run: Get Mesh Coordinates",ESMF_LOGMSG_INFO,rc=rclocal)   

      endif

      call ESMF_FieldGet(oceanV,mesh=oceanMesh,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("Ocean Run: Get Grid from Export Field oceanV",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_FieldGet(oceanU,farrayPtr=fptr1D3,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("Ocean Run: Get Export Field oceanU",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_FieldGet(oceanV,farrayPtr=fptr1D4,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("Ocean Run: Get Export Field oceanV",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Loop through nodes and get the ocean Velocity along X/Y directions
      do i=1,numOwnedNodes
        fptr1D3(i)=-1.0 ! velocityX(oceanIDs(i))
        fptr1D4(i)=-0.34 ! velocityY(oceanIDs(i))
      enddo

      ! Deallocate space to hold local node coordinates
      deallocate(ownedNodeCoords)

      rclocal=ESMF_SUCCESS

  end subroutine ocean_run
  ! =====================================================================================
  subroutine ocean_final(comp,importState,exportState,clock,rclocal)

      type(ESMF_GridComp)::comp
      type(ESMF_State)::importState,exportState
      type(ESMF_Clock)::clock
      integer,intent(out)::rclocal

      ! Local variables
      type(ESMF_Field)::oceanH,newoceanH,oceanU,oceanV
      type(ESMF_Mesh)::oceanMesh

      rclocal=ESMF_SUCCESS

      ! Check validity of results
      ! Get Fields from export state
      call ESMF_StateGet(exportState,"oceanH",oceanH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Ocean Finalize: Get Export State oceanH",ESMF_LOGMSG_INFO,rc=rclocal)    

      call ESMF_StateGet(exportState,"oceanU",oceanU,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Ocean Finalize: Get Export State oceanU",ESMF_LOGMSG_INFO,rc=rclocal)    

      call ESMF_StateGet(exportState,"oceanV",oceanV,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Ocean Finalize: Get Export State oceanV",ESMF_LOGMSG_INFO,rc=rclocal)    
      
      ! Get Fields from import state
      call ESMF_StateGet(importState,"newoceanH",newoceanH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Ocean Finalize: Get Import State",ESMF_LOGMSG_INFO,rc=rclocal)    
      
      ! Garbage collection
      call ESMF_FieldGet(oceanH,mesh=oceanMesh,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Ocean Finalize: Get Export Field oceanH",ESMF_LOGMSG_INFO,rc=rclocal)    

      call ESMF_FieldGet(oceanU,mesh=oceanMesh,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Ocean Finalize: Get Export Field oceanU",ESMF_LOGMSG_INFO,rc=rclocal)  

      call ESMF_FieldGet(oceanV,mesh=oceanMesh,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Ocean Finalize: Get Export Field oceanV",ESMF_LOGMSG_INFO,rc=rclocal)  

      call ESMF_FieldGet(newoceanH,mesh=oceanMesh,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Ocean Finalize: Get Import Field",ESMF_LOGMSG_INFO,rc=rclocal)    

      call ESMF_FieldDestroy(oceanH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Ocean Finalize: Export Field Destroyed oceanH",ESMF_LOGMSG_INFO,rc=rclocal)   

      call ESMF_FieldDestroy(oceanU,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Ocean Finalize: Export Field Destroyed oceanU",ESMF_LOGMSG_INFO,rc=rclocal)   

      call ESMF_FieldDestroy(oceanV,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Ocean Finalize: Export Field Destroyed oceanV",ESMF_LOGMSG_INFO,rc=rclocal)   

      call ESMF_FieldDestroy(newoceanH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Ocean Finalize: Import Field Destroyed",ESMF_LOGMSG_INFO,rc=rclocal)   

      call ESMF_MeshDestroy(oceanMesh,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("Ocean Finalize: Mesh Destroyed",ESMF_LOGMSG_INFO,rc=rclocal)    
      
      ! Return success
      rclocal=ESMF_SUCCESS 

  end subroutine ocean_final
  ! =====================================================================================
 
end module esmfOcean
    
