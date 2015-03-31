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
!       Filename:  SpmComponent.f90
!
!    Description:  Defines the Surface Process Model ESMF Component
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles     
!
! =====================================================================================

module esmfSPM

  use topology
  use hydrology
  use geomorpho
  use parameters

  implicit none

  public spm_register

contains

  ! =====================================================================================
  subroutine spm_register(comp,rclocal)

    type(ESMF_GridComp)::comp
    integer,intent(out)::rclocal

    rclocal=ESMF_SUCCESS

    ! Register the callback routines
    call ESMF_GridCompSetEntryPoint(comp,ESMF_METHOD_INITIALIZE,spm_init,rc=rclocal)
    if(rclocal/=ESMF_SUCCESS) return
    call ESMF_LogWrite("SPM Initialize Method Done",ESMF_LOGMSG_INFO,rc=rclocal)

    call ESMF_GridCompSetEntryPoint(comp,ESMF_METHOD_RUN,spm_run,rc=rclocal)
    if(rclocal/=ESMF_SUCCESS) return
    call ESMF_LogWrite("SPM Run Method Done",ESMF_LOGMSG_INFO,rc=rclocal)

    call ESMF_GridCompSetEntryPoint(comp,ESMF_METHOD_FINALIZE,spm_final,rc=rclocal)
    if(rclocal/=ESMF_SUCCESS) return
    call ESMF_LogWrite("SPM Finalize Method Done",ESMF_LOGMSG_INFO,rc=rclocal)

    call ESMF_LogWrite("SPM Model Component Registration Completed",ESMF_LOGMSG_INFO,rc=rclocal)

  end subroutine spm_register
  ! =====================================================================================
  subroutine spm_init(comp,importState,exportState,clock,rclocal)

      type(ESMF_GridComp)::comp
      type(ESMF_State)::importState,exportState
      type(ESMF_Clock)::clock
      integer,intent(out)::rclocal

      ! Local variables
      type(ESMF_VM)::vm
      type(ESMF_Field)::spmH,newspmH,spmU,spmV,spmvD,spmrainH
      type(ESMF_Mesh)::spmMesh
      type(ESMF_ArraySpec)::arrayspec
      integer::localPET,petCount,localrc,i,k,p
      real(ESMF_KIND_R8),pointer::fptr1D(:),fptr1D2(:),fptr1D3(:)
      real(ESMF_KIND_R8),pointer::fptr1D4(:),fptr1D5(:),fptr1D6(:)
      integer,pointer::nodeIds(:),nodeOwners(:)
      real(ESMF_KIND_R8),pointer::nodeCoords(:)
      integer,pointer::elemIds(:),elemTypes(:),elemConn(:)
      integer::numNodes,numElems

      rclocal=ESMF_SUCCESS

      ! Query component for VM and create a layout with the right breakdown
      call ESMF_GridCompGet(comp,vm=vm,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Initialize: Get Component",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_VMGet(vm,localPet=localPet,petCount=petCount,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Initialize: Get VM",ESMF_LOGMSG_INFO,rc=rclocal)

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!     Build SPM mesh
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

      ! Set number of nodes
      numNodes=upartN

      ! Allocate and fill the node id array
      ! Allocate and fill node coordinate array.
      ! Since this is a 2D Mesh the size is 2x the
      ! number of nodes.
      allocate(nodeIds(numNodes),nodeCoords(2*numNodes))
      if(allocated(spmIDs)) deallocate(spmIDs)
      allocate(spmIDs(uOwnedNode))
      allocate(nodeOwners(numNodes))
      k=1
      p=1
      do i=1,numNodes 
        nodeIds(i)=unodeId(i)
        nodeCoords(k)=tcoordX(nodeIds(i))
        k=k+1
        nodeCoords(k)=tcoordY(nodeIds(i))
        k=k+1
        nodeOwners(i)=uownedID(i)
        if(nodeOwners(i)==localPet)then 
          spmIDs(p)=unodeID(i)
          p=p+1
        endif
      enddo

      ! Set the number of each type of element, plus the total number.
      numElems=upartE

      ! Allocate and fill the element id array.
      ! Allocate and fill the element topology type array.
      ! Allocate and fill the element connection type array.
      allocate(elemIds(numElems))
      allocate(elemTypes(numElems))
      allocate(elemConn(3*numElems))
      k=1
      do i=1,numElems 
        elemIds(i)=uelemID(i)
        elemTypes(i)=ESMF_MESHELEMTYPE_TRI
        elemConn(k)=unodeLID(delmt(elemIds(i),1))
        k=k+1
        elemConn(k)=unodeLID(delmt(elemIds(i),2))
        k=k+1
        elemConn(k)=unodeLID(delmt(elemIds(i),3))
        k=k+1
      enddo

      ! Create Mesh structure in 1 step
      spmMesh=ESMF_MeshCreate(parametricDim=2,spatialDim=2,nodeIds=nodeIds, &
        nodeCoords=nodeCoords,nodeOwners=nodeOwners,elementIds=elemIds,&
        elementTypes=elemTypes,elementConn=elemConn,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("SPM Initialize: Create Mesh",ESMF_LOGMSG_INFO,rc=rclocal)

      ! deallocate node data
      deallocate(nodeIds)
      deallocate(nodeCoords)
      deallocate(nodeOwners)

      ! deallocate elem data
      deallocate(elemIds)
      deallocate(elemTypes)
      deallocate(elemConn)
    
      ! Create interpolated topography field
      call ESMF_ArraySpecSet(arrayspec,1,ESMF_TYPEKIND_R8,rc=rclocal)

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!     Import states
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

      ! Create Import Fields
      spmH=ESMF_FieldCreate(spmMesh,arrayspec,name="spmH",rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("SPM Initialize: Create Import Field spmH",ESMF_LOGMSG_INFO,rc=rclocal)
      
      ! Load data into the source Field
      call ESMF_FieldGet(spmH,farrayPtr=fptr1D,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
          rclocal=ESMF_FAILURE
          return
      endif      
      call ESMF_LogWrite("SPM Initialize: Get Import Field spmH",ESMF_LOGMSG_INFO,rc=rclocal)
      
      spmU=ESMF_FieldCreate(spmMesh,arrayspec,name="spmU",rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("SPM Initialize: Create Import Field spmU",ESMF_LOGMSG_INFO,rc=rclocal)
      
      ! Load data into the source Field
      call ESMF_FieldGet(spmU,farrayPtr=fptr1D3,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
          rclocal=ESMF_FAILURE
          return
      endif      
      call ESMF_LogWrite("SPM Initialize: Get Import Field spmU",ESMF_LOGMSG_INFO,rc=rclocal)
   
      spmV=ESMF_FieldCreate(spmMesh,arrayspec,name="spmV",rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("SPM Initialize: Create Import Field spmV",ESMF_LOGMSG_INFO,rc=rclocal)
   
      ! Load data into the source Field
      call ESMF_FieldGet(spmV,farrayPtr=fptr1D4,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
          rclocal=ESMF_FAILURE
          return
      endif      
      call ESMF_LogWrite("SPM Initialize: Get Import Field spmV",ESMF_LOGMSG_INFO,rc=rclocal)
   
      spmvD=ESMF_FieldCreate(spmMesh,arrayspec,name="spmvD",rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("SPM Initialize: Create Import Field spmvD",ESMF_LOGMSG_INFO,rc=rclocal)
   
      call ESMF_FieldGet(spmvD,farrayPtr=fptr1D5,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
          rclocal=ESMF_FAILURE
          return
      endif      
      call ESMF_LogWrite("SPM Initialize: Get Import Field spmvD",ESMF_LOGMSG_INFO,rc=rclocal)
   
      spmrainH=ESMF_FieldCreate(spmMesh,arrayspec,name="spmrainH",rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("SPM Initialize: Create Import Field spmrainH",ESMF_LOGMSG_INFO,rc=rclocal)
   
      ! Load data into the source Field
      call ESMF_FieldGet(spmrainH,farrayPtr=fptr1D6,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
          rclocal=ESMF_FAILURE
          return
      endif      
      call ESMF_LogWrite("SPM Initialize: Get Import Field spmrainH",ESMF_LOGMSG_INFO,rc=rclocal)
   
      ! Clear interpolated topography field
      fptr1D=0.0
      fptr1D3=0.0
      fptr1D4=0.0
      fptr1D5=0.0
      fptr1D6=0.0
     
      ! Add Fields to State
      call ESMF_StateAdd(importState,(/spmH,spmU,spmV,spmvD,spmrainH/),rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Initialize: Add Import States",ESMF_LOGMSG_INFO,rc=rclocal)

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!     Export states
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

      newspmH=ESMF_FieldCreate(spmMesh,arrayspec,name="newspmH",rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("SPM Initialize: Create Export Field newspmH",ESMF_LOGMSG_INFO,rc=rclocal)
      
      call ESMF_FieldGet(newspmH,farrayPtr=fptr1D2,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
          rclocal=ESMF_FAILURE
          return
      endif      
      call ESMF_LogWrite("SPM Initialize: Get Export Field newspmH",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Clear interpolated topography field
      fptr1D2=0.0

      ! Add Field to State
      call ESMF_StateAdd(exportState,(/newspmH/),rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Initialize: Add Export State",ESMF_LOGMSG_INFO,rc=rclocal)
      
      rclocal=ESMF_SUCCESS

  end subroutine spm_init
  ! =====================================================================================
  subroutine spm_run(comp,importState,exportState,clock,rclocal)

      type(ESMF_GridComp)::comp
      type(ESMF_State)::importState,exportState
      type(ESMF_Clock)::clock
      integer,intent(out)::rclocal

      ! Local variables
      type(ESMF_VM)::vm
      type(ESMF_Field)::spmH,newspmH,spmU,spmV,spmvD,spmrainH
      type(ESMF_Mesh)::spmMesh
      real(ESMF_KIND_R8)::time1,time2
      integer::localrc,i,id
      integer::localPet,petCount
      real(ESMF_KIND_R8),pointer::fptr1D(:),fptr1D2(:),fptr1D3(:),fptr1D4(:),fptr1D5(:),fptr1D6(:)
      real(ESMF_KIND_R8),pointer:: uval(:),vval(:),vDval(:),dem(:),rainval(:)
      integer::numOwnedNodes
      real(ESMF_KIND_R8),pointer::ownedNodeCoords(:),newelev(:)

      rclocal=ESMF_SUCCESS

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

      if(.not.updateSPM_elevation)then

        ! Get information from the component.
        call ESMF_StateGet(importState,"spmH",spmH,rc=rclocal)
        if(rclocal/=ESMF_SUCCESS) return
        call ESMF_LogWrite("SPM Run: Get Import State spmH",ESMF_LOGMSG_INFO,rc=rclocal)

        ! Get Grid from field
        call ESMF_FieldGet(spmH,mesh=spmMesh,rc=localrc)
        if(localrc/=ESMF_SUCCESS)then
          rclocal=ESMF_FAILURE
          return
        endif
        call ESMF_LogWrite("SPM Run: Get Grid from Import Field spmH",ESMF_LOGMSG_INFO,rc=rclocal)
       
        ! Get number of local nodes to allocate space
        ! to hold local node coords
        call ESMF_MeshGet(spmMesh,numOwnedNodes=numOwnedNodes,rc=rclocal)  
        call ESMF_LogWrite("SPM Run: Get Unstructured Mesh",ESMF_LOGMSG_INFO,rc=rclocal)    

        ! Allocate space to hold local node coordinates
        ! (spatial dimension of Mesh*number of local nodes)
        allocate(ownedNodeCoords(2*numOwnedNodes))

        ! Get local node coordinates
        call ESMF_MeshGet(spmMesh,ownedNodeCoords=ownedNodeCoords,rc=rclocal)
        call ESMF_LogWrite("SPM Run: Get Mesh Coordinates",ESMF_LOGMSG_INFO,rc=rclocal)    

        allocate(dem(numOwnedNodes))

        ! Check destination field
        call ESMF_FieldGet(spmH,farrayPtr=fptr1D,rc=localrc)
        if(localrc/=ESMF_SUCCESS)then
          rclocal=ESMF_FAILURE
          return
        endif
        call ESMF_LogWrite("SPM Run: Get Import Field spmH",ESMF_LOGMSG_INFO,rc=rclocal)

      endif
      
      ! Get information from the component.
      call ESMF_StateGet(importState,"spmU",spmU,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Run: Get Import State spmU",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Get information from the component.
      call ESMF_StateGet(importState,"spmV",spmV,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Run: Get Import State spmV",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Get information from the component.
      call ESMF_StateGet(importState,"spmvD",spmvD,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Run: Get Import State spmvD",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Get information from the component.
      call ESMF_StateGet(importState,"spmrainH",spmrainH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Run: Get Import State spmrainH",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Get Grid from field
      call ESMF_FieldGet(spmU,mesh=spmMesh,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("SPM Run: Get Grid from Import Field spmU",ESMF_LOGMSG_INFO,rc=rclocal)
      
      if(updateSPM_elevation)then

        ! Get number of local nodes to allocate space
        ! to hold local node coords
        call ESMF_MeshGet(spmMesh,numOwnedNodes=numOwnedNodes,rc=rclocal)  
        call ESMF_LogWrite("SPM Run: Get Unstructured Mesh",ESMF_LOGMSG_INFO,rc=rclocal)    

        ! Allocate space to hold local node coordinates
        ! (spatial dimension of Mesh*number of local nodes)
        allocate(ownedNodeCoords(2*numOwnedNodes))

        ! Get local node coordinates
        call ESMF_MeshGet(spmMesh,ownedNodeCoords=ownedNodeCoords,rc=rclocal)
        call ESMF_LogWrite("SPM Run: Get Mesh Coordinates",ESMF_LOGMSG_INFO,rc=rclocal)    

        allocate(dem(numOwnedNodes))

        ! Get information from the component.
        call ESMF_StateGet(importState,"spmH",spmH,rc=rclocal)
        if(rclocal/=ESMF_SUCCESS) return
        call ESMF_LogWrite("SPM Run: Get Import State spmH",ESMF_LOGMSG_INFO,rc=rclocal)

        ! Check destination field
        call ESMF_FieldGet(spmH,farrayPtr=fptr1D,rc=localrc)
        if(localrc/=ESMF_SUCCESS)then
          rclocal=ESMF_FAILURE
          return
        endif
        call ESMF_LogWrite("SPM Run: Get Import Field spmH",ESMF_LOGMSG_INFO,rc=rclocal)

      endif
     
      ! Get Grid from field
      call ESMF_FieldGet(spmV,mesh=spmMesh,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("SPM Run: Get Grid from Import Field spmV",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Get Grid from field
      call ESMF_FieldGet(spmvD,mesh=spmMesh,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("SPM Run: Get Grid from Import Field spmvD",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Get Grid from field
      call ESMF_FieldGet(spmrainH,mesh=spmMesh,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("SPM Run: Get Grid from Import Field spmrainH",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Check destination field
      call ESMF_FieldGet(spmU,farrayPtr=fptr1D3,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("SPM Run: Get Import Field spmU",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Check destination field
      call ESMF_FieldGet(spmV,farrayPtr=fptr1D4,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("SPM Run: Get Import Field spmV",ESMF_LOGMSG_INFO,rc=rclocal) 

      ! Check destination field
      call ESMF_FieldGet(spmvD,farrayPtr=fptr1D5,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("SPM Run: Get Import Field spmvD",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Check destination field
      call ESMF_FieldGet(spmrainH,farrayPtr=fptr1D6,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("SPM Run: Get Import Field spmrainH",ESMF_LOGMSG_INFO,rc=rclocal)

      ! Loop through nodes and defined the interpolated topographic value
      if(.not.updateSPM_elevation)then
        allocate(newelev(dnodes))
        newelev=-1.e6
      endif
      allocate(uval(dnodes),vval(dnodes))
      allocate(vDval(dnodes),rainval(dnodes))
      uval=-1.e6
      vval=-1.e6
      vDval=-1.e6
      rainval=-1.e6
 
      do i=1,numOwnedNodes
        id=spmIDs(i)
        if(.not.updateSPM_elevation) newelev(id)=fptr1D(i)
        uval(id)=fptr1D3(i)
        vval(id)=fptr1D4(i)
        vDval(id)=fptr1D5(i)
        rainval(id)=fptr1D6(i)
      enddo
      
      if(.not.updateSPM_elevation)then
        call ESMF_VMAllReduce(vm=vm,sendData=newelev,recvData=tcoordZ,count=dnodes,&
          reduceflag=ESMF_REDUCE_MAX,rc=rc)
        deallocate(newelev)
      endif
      
      call ESMF_VMAllReduce(vm=vm,sendData=vDval,recvData=tvertDisp,count=dnodes,&
          reduceflag=ESMF_REDUCE_MAX,rc=rc)

      call ESMF_VMAllReduce(vm=vm,sendData=rainval,recvData=precipitation,count=dnodes,&
          reduceflag=ESMF_REDUCE_MAX,rc=rc)

      ! Reset destination to 0
      if(updateSPM_elevation) fptr1D=0.0
      fptr1D3=0.0
      fptr1D4=0.0
      fptr1D5=0.0
      
      ! This is where we do SPM run
      call ESMF_VMWtime(time1,rc=rc)
      
      ! Gormorphology model
      call geomorphology
      if(pet_id==0)print*,'-------------------------'
      call ESMF_VMWtime(time2,rc=rc)
      if(pet_id==0) print*,'BADLANDS SPM step Completed (s) ',time2-time1
      if(pet_id==0)print*,'-------------------------'

      ! Update elevation due to surface processes
      do i=1,numOwnedNodes
        dem(i)=tcoordZ(spmIDs(i))
      enddo

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!     Export state
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

      call ESMF_StateGet(exportState,"newspmH",newspmH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Run: Get Export State",ESMF_LOGMSG_INFO,rc=rclocal)

      call ESMF_FieldGet(newspmH,mesh=spmMesh,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("SPM Run: Get Grid from Export Field",ESMF_LOGMSG_INFO,rc=rclocal)


      call ESMF_FieldGet(newspmH,farrayPtr=fptr1D2,rc=localrc)
      if(localrc/=ESMF_SUCCESS)then
        rclocal=ESMF_FAILURE
        return
      endif
      call ESMF_LogWrite("SPM Run: Get Export Field",ESMF_LOGMSG_INFO,rc=rclocal)

      fptr1D2=dem

      rclocal=ESMF_SUCCESS

      ! Deallocate space to hold local node coordinates
      deallocate(ownedNodeCoords,uval,vval,dem,vDval,rainval)

  end subroutine spm_run
  ! =====================================================================================
  subroutine spm_final(comp,importState,exportState,clock,rclocal)

      type(ESMF_GridComp)::comp
      type(ESMF_State)::importState,exportState
      type(ESMF_Clock)::clock
      integer,intent(out)::rclocal

      ! Local variables
      type(ESMF_Field)::spmH,newspmH,spmU,spmV,spmvD,spmrainH
      type(ESMF_Mesh)::spmMesh

      rclocal=ESMF_SUCCESS

      ! Clean SPM model 
      call clean_SPM_model

      ! Check validity of results
      ! Get Fields from import state
      call ESMF_StateGet(importState,"spmH",spmH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Finalize: Get Import State spmH",ESMF_LOGMSG_INFO,rc=rclocal)    

      call ESMF_StateGet(importState,"spmU",spmU,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Finalize: Get Import State spmU",ESMF_LOGMSG_INFO,rc=rclocal)    

      call ESMF_StateGet(importState,"spmV",spmV,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Finalize: Get Import State spmV",ESMF_LOGMSG_INFO,rc=rclocal)    
          
      call ESMF_StateGet(importState,"spmvD",spmvD,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Finalize: Get Import State spmvD",ESMF_LOGMSG_INFO,rc=rclocal)    

      call ESMF_StateGet(importState,"spmrainH",spmrainH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Finalize: Get Import State spmrainH",ESMF_LOGMSG_INFO,rc=rclocal)   
  
      ! Get Fields from export state
      call ESMF_StateGet(exportState,"newspmH",newspmH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Finalize: Get Export State",ESMF_LOGMSG_INFO,rc=rclocal)    
      
      ! Garbage collection
      call ESMF_FieldGet(spmH,mesh=spmMesh,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Finalize: Get Import Fiel spmH",ESMF_LOGMSG_INFO,rc=rclocal)    

      call ESMF_FieldGet(spmU,mesh=spmMesh,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Finalize: Get Import Field spmU",ESMF_LOGMSG_INFO,rc=rclocal)    
     
      call ESMF_FieldGet(spmV,mesh=spmMesh,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Finalize: Get Import Field spmV",ESMF_LOGMSG_INFO,rc=rclocal)    

      call ESMF_FieldGet(spmvD,mesh=spmMesh,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Finalize: Get Import Field spmvD",ESMF_LOGMSG_INFO,rc=rclocal)    

      call ESMF_FieldGet(spmrainH,mesh=spmMesh,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Finalize: Get Import Field spmrainH",ESMF_LOGMSG_INFO,rc=rclocal)    

      call ESMF_FieldGet(newspmH,mesh=spmMesh,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Finalize: Get Export Field",ESMF_LOGMSG_INFO,rc=rclocal)    

      call ESMF_FieldDestroy(spmH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Finalize: Import Field spmH Destroyed",ESMF_LOGMSG_INFO,rc=rclocal)    

      call ESMF_FieldDestroy(spmU,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Finalize: Import Field spmU Destroyed",ESMF_LOGMSG_INFO,rc=rclocal)  

      call ESMF_FieldDestroy(spmV,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Finalize: Import Field spmV Destroyed",ESMF_LOGMSG_INFO,rc=rclocal) 

      call ESMF_FieldDestroy(spmvD,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Finalize: Import Field spmvD Destroyed",ESMF_LOGMSG_INFO,rc=rclocal) 

      call ESMF_FieldDestroy(spmrainH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Finalize: Import Field spmrainH Destroyed",ESMF_LOGMSG_INFO,rc=rclocal) 

      call ESMF_FieldDestroy(newspmH,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Finalize: Export Field Destroyed",ESMF_LOGMSG_INFO,rc=rclocal)    

      call ESMF_MeshDestroy(spmMesh,rc=rclocal)
      if(rclocal/=ESMF_SUCCESS) return
      call ESMF_LogWrite("SPM Finalize: Mesh Destroyed",ESMF_LOGMSG_INFO,rc=rclocal)    
      
      ! Return success
      rclocal=ESMF_SUCCESS 

  end subroutine spm_final
  ! =====================================================================================
 
end module esmfSPM
    
