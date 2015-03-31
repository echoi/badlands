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
!       Filename:  WatershedDelineation.f90
!
!    Description:  Sub-basin partitioning based upon channel network junctions
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles     
!
! =====================================================================================

module watershed

  use topology
  use parameters
  use hydroUtil
  use hydrology

  implicit none

  integer::junctionsNb
  integer,dimension(:),allocatable::rcvNb
  integer,dimension(:,:),allocatable::rcvIDs
  integer,dimension(:),allocatable::junctionIDs

  integer,dimension(:),allocatable::options
  integer,pointer::vsize,adjwgt 
  real(ESMF_KIND_R8),pointer::tpwgts,ubvec

contains

  ! =====================================================================================
  subroutine compute_cumulative_discharge

    integer::k,rcv 
    
    if(.not.allocated(discharge)) allocate(discharge(dnodes))

    discharge=0.

    if(allocated(rcvNb)) deallocate(rcvNb)
    if(allocated(rcvIDs)) deallocate(rcvIDs)
    allocate(rcvNb(dnodes),rcvIDs(dnodes,maxrcvs*2))

    rcvNb=0
    rcvIDs=0

    ! Compute local discharge 
    do k=1,dnodes
!       if(watercell(k)==0.)then 
        discharge(k)=precipitation(k)*voronoiCell(k)%area 
!       else
!         discharge(k)=0.
!       endif
    enddo

    ! Compute drainage area and receivers IDs
    do k=dnodes,1,-1
      rcv=receivers(stackOrder(k))
      ! Drainage
      if(rcv/=stackOrder(k))then
        discharge(rcv)=discharge(rcv)+discharge(stackOrder(k))
      endif
      ! Receivers
      rcvNb(rcv)=rcvNb(rcv)+1
      rcvIDs(rcv,rcvNb(rcv))=stackOrder(k)
    enddo

  end subroutine compute_cumulative_discharge
  ! =====================================================================================
  subroutine compute_subcatchment
  
    logical::inTree
    integer::id,k,n,p,s,maxs,subcatchID
    integer::localjctns(npets+1),lsubID(dnodes)
    real(ESMF_KIND_R8)::accu_limit

    accu_limit=exp(accu_thres)
    if(.not.allocated(strahler)) allocate(strahler(dnodes))
    if(.not.allocated(subcatchmentID)) allocate(subcatchmentID(dnodes))
    if(.not.allocated(junctionIDs)) allocate(junctionIDs(dnodes))

    strahler=0
    junctionsNb=0
    
    subcatchID=0
    junctionIDs=0

    ! Compute Strahler stream order
    do id=dnodes,1,-1
      maxs=0
      k=stackOrder(id)
      ! outlet
      if(receivers(k)==k)then 
        junctionsNb=junctionsNb+1 
        junctionIDs(junctionsNb)=id
      else
        if(discharge(k)>accu_limit)then
          n=0
          do p=1,rcvNb(k)
            s=strahler(rcvIDs(k,p))
            if(s>=maxs)then
              if(s==maxs)then
                n=n+1
              else
                n=1
              endif
              maxs=s
            endif
          enddo
          if(n==1)strahler(k)=maxs
          if(n>1)strahler(k)=maxs+1
          if(receivers(k)==k)then
            junctionsNb=junctionsNb+1 
            junctionIDs(junctionsNb)=id
          else
            loop_stream:do p=1,rcvNb(k)
              s=strahler(rcvIDs(k,p))
              if(s<strahler(k).and.s>0)then
                junctionsNb=junctionsNb+1 
                junctionIDs(junctionsNb)=id
                exit loop_stream
              endif
            enddo loop_stream
          endif
        endif
      endif
    enddo

    if(allocated(subcatchNb)) deallocate(subcatchNb)
    allocate(subcatchNb(junctionsNb))

    ! Distribute junctions to pets
    localjctns=0
    localjctns(1)=1
    p=junctionsNb/npets
    n=mod(junctionsNb,npets)
    do s=0,npets-1
      if(s<n)then
        k=p+1
      else 
        k=p
      endif
      localjctns(s+2)=k+localjctns(s+1)
    enddo

    ! Define subcatchments parameters
    lsubID=junctionsNb+1
    do s=localjctns(pet_id+1),localjctns(pet_id+2)-1 
      id=junctionIDs(s)
      subcatchID=s 
      loop_stack:do p=id,dnodes 
        n=stackOrder(p)
        k=receivers(n)
        inTree=.false.
        loop_rcv:do while(filldem(k)>=filldem(stackOrder(id)))
          if(k==stackOrder(id))then
            inTree=.true.
            exit loop_rcv
          endif
          if(k==receivers(k)) exit loop_rcv
          k=receivers(k)
          if(filldem(k)<filldem(stackOrder(id)))exit loop_stack
        enddo loop_rcv
        if(lsubID(n)==junctionsNb+1.and.inTree) & 
          lsubID(n)=subcatchID
        if(s<junctionsNb)then
          if(stackOrder(p)==stackOrder(junctionIDs(s+1)))exit loop_stack
        endif
      enddo loop_stack
      ! Add the subcatchment outlet
      lsubID(stackOrder(id))=subcatchID
    enddo

    ! Merge subcatchment parameters 
    call ESMF_VMAllReduce(vm=vm,sendData=lsubID,recvData=subcatchmentID,count=dnodes,&
      reduceflag=ESMF_REDUCE_MIN,rc=rc)

    ! Count nodes per subcatchment and store layers information
    subcatchNb=0
    do id=1,dnodes
      subcatchNb(subcatchmentID(id))=subcatchNb(subcatchmentID(id))+1
    enddo

    ! Perform subcatchment load-balancing and partitioning
    call metis_loadbalancing

  end subroutine compute_subcatchment
  ! =====================================================================================
  subroutine metis_loadbalancing

    integer::l,s,id,k,p,objval

    integer,dimension(:),allocatable::partition,vwgt,CRSadj,CRSadjncy,connectsNb
    integer,dimension(:,:),allocatable::connects

    if(allocated(partition)) deallocate(partition)
    allocate(partition(junctionsNb+1))

    if(pet_id==0)then
      if(allocated(vwgt)) deallocate(vwgt,CRSadj,CRSadjncy)
      if(allocated(connects)) deallocate(connects,connectsNb)
      allocate(vwgt(junctionsNb+1))
      allocate(CRSadj(junctionsNb+2),CRSadjncy(2*(junctionsNb)))
      allocate(connectsNb(junctionsNb+1),connects(junctionsNb+1,junctionsNb+1))

      connects=0
      connectsNb=0
      do s=1,junctionsNb+1
        if(s<=junctionsNb)then
          id=stackOrder(junctionIDs(s))
          k=receivers(id)
          l=subcatchmentID(k)
          vwgt(s)=subcatchNb(s)
        else
          l=s
          vwgt(s)=1
        endif
        ! Declare tree connection parameters
        if(l/=s)then
          connectsNb(s)=connectsNb(s)+1
          connectsNb(l)=connectsNb(l)+1      
          connects(s,connectsNb(s))=l
          connects(l,connectsNb(l))=s
        ! Link to the top root 
        elseif(s<=junctionsNb)then
          l=junctionsNb+1
          connectsNb(s)=connectsNb(s)+1
          connectsNb(l)=connectsNb(l)+1      
          connects(s,connectsNb(s))=l
          connects(l,connectsNb(l))=s
        endif
      enddo

      ! Define the metis graph data structure CRS
      CRSadj=1
      do s=2,junctionsNb+2
        CRSadj(s)=CRSadj(s-1)+connectsNb(s-1)
      enddo

      CRSadjncy=0
      k=1
      do p=1,junctionsNb+1
        do l=1,connectsNb(p)
          CRSadjncy(k)=connects(p,l)
          k=k+1
        enddo
      enddo
      
      ! METIS_OPTION_NUMBERING
      allocate(options(0:40))
      call METIS_SetDefaultOptions(options)
      options(17)=1

      ! Nullify pointers
      vsize=>null()
      adjwgt=>null() 
      tpwgts=>null()
      ubvec=>null()

      ! K-way metis partitioning
      partition=1
      if(npets>1) call METIS_PartGraphKway(junctionsNb+1,1,CRSadj,CRSadjncy,vwgt,vsize,adjwgt,npets,tpwgts,ubvec,options,objval,partition) 
      deallocate(options)
    endif

    ! Broadcast partitioning
    call ESMF_VMBroadcast(vm=vm,bcstData=partition,count=junctionsNb+1,rootPet=0,rc=rc)
    
    if(.not.allocated(subcatchmentProc)) allocate(subcatchmentProc(dnodes))
    subcatchmentProc=-1

    ! Define for each partition local nodes and their global IDs
    localNodes=0
    do k=1,dnodes
      subcatchmentProc(k)=partition(subcatchmentID(k))-1
      if(pet_id==subcatchmentProc(k))localNodes=localNodes+1
    enddo

    ! Define local and inter-processor communications
    if(.not.allocated(sendprocID)) allocate(sendprocID(dnodes))
    sendprocID=-1
    if(.not.allocated(rcvprocNb)) allocate(rcvprocNb(dnodes))
    rcvprocNb=0
    if(.not.allocated(rcvprocID)) allocate(rcvprocID(dnodes,maxrcvs))
    rcvprocID=-1
    if(.not.allocated(rcvsendID)) allocate(rcvsendID(dnodes,maxrcvs))
    rcvsendID=-1
    if(allocated(localNodesGID)) deallocate(localNodesGID)
    allocate(localNodesGID(localNodes))

    s=1
    drainOde=localNodes
    do k=1,dnodes
      p=stackOrder(k) 
      id=receivers(p)
      ! Inter-partition nodes variables
      if(subcatchmentProc(p)/=subcatchmentProc(id))then
        sendprocID(p)=subcatchmentProc(id)
        rcvprocNb(id)=rcvprocNb(id)+1
        rcvprocID(id,rcvprocNb(id))=subcatchmentProc(p)
        rcvsendID(id,rcvprocNb(id))=p
        drainOde=drainOde+1
      endif
      if(pet_id==subcatchmentProc(p))then
        localNodesGID(s)=k
        s=s+1
      endif
    enddo

    ! Cleaning process
    nullify(vsize,adjwgt,tpwgts,ubvec)

  end subroutine metis_loadbalancing
  ! =====================================================================================

end module watershed 