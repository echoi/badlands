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
!       Filename:  StratalFunctions.f90
!
!    Description:  Stratigraphic grid evolution functions.
!
!        Version:  1.0 
!        Created:  03/06/15 09:05:27
!        Revision:  none
!
!        Author:  Tristan Salles     
!
! =====================================================================================

module strata_evol

  use parallel
  use parameters
  use out_stratal
  use stratal_read
  use stratal_class

  implicit none

contains

  ! =====================================================================================

  subroutine StrataGen
    
    ! Get stratal definition
    call stratal_parser

    ! Built stratigraphy grid
    if(totgrn>0) call constructStrata

    return

  end subroutine StrataGen
  ! =====================================================================================

  subroutine compute_stratal_displacement

    integer::k

    do k=1,lsnb
      lay_base(k)=lay_base(k)+svertDisp(k)*time_step
    enddo

  end subroutine compute_stratal_displacement
  ! =====================================================================================

  subroutine update_stratigraphy_layer
    
    integer::k,l,gid,ks,nb,ll

    real(kind=8)::remain,sed

    ! Loop through the stratigraphic nodes associated to the vertex and 
    ! distribute the new active layer composition to the regular nodes
    do k=1,dnodes
      nb=regPtNb(k)
      do l=1,nb
        gid=regPtID(k,l)
        do ks=1,totgrn
          nlay_sed(gid,ks)=alay_dsed(k,ks)
        enddo
      enddo
    enddo

    ! Record the active layer changes in the stratigraphic mesh
    do k=1,lsnb
      gid=lnID(k)
      do ks=1,totgrn
        sed=nlay_sed(gid,ks)-alay_sed(gid,ks)
        ! Deposition
        if(sed>=0)then
          lay_sed(k,layNb,ks)=lay_sed(k,layNb,ks)+sed
          lay_thick(k,layNb)=lay_thick(k,layNb)+sed
        else
          remain=-sed
          lpp: do ll=layNb,1,-1
            if(lay_sed(k,ll,ks)>remain)then
              lay_sed(k,ll,ks)=lay_sed(k,ll,ks)-remain
              lay_thick(k,ll)=lay_thick(k,ll)-remain
              exit lpp
            else
              remain=remain-lay_sed(k,ll,ks)
              lay_thick(k,ll)=lay_thick(k,ll)-lay_sed(k,ll,ks)
              lay_sed(k,ll,ks)=0.0
            endif
          enddo lpp
        endif
      enddo
    enddo

    ! Build new active layer
    call buildActiveLayer

    return

  end subroutine update_stratigraphy_layer
  ! =====================================================================================

  subroutine buildActiveLayer

    integer::k,l,ks,gid,nb

    real(kind=8)::th,sed(totgrn),prop(totgrn)
    real(kind=8),dimension(snodes)::tsed,gsed    

    if(.not.allocated(alay_thick))then
      allocate(alay_thick(snodes))
      allocate(alay_sed(snodes,totgrn))
    endif

    ! Get the composition and thickness of the active layer
    alay_thick=0.0
    alay_sed=0.0
    do k=1,lsnb
      th=0.0
      sed=0.0
      gid=lnID(k)
      lp:do l=layNb,1,-1
        th=th+lay_thick(k,l)
        if(th>active_thick)then
          th=active_thick
          do ks=1,totgrn
            prop(ks)=lay_sed(k,l,ks)/lay_thick(k,l)
            sed(ks)=sed(ks)+prop(ks)*(th-alay_thick(gid))
          enddo
          alay_thick(gid)=th
        else
          alay_thick(gid)=th
          do ks=1,totgrn
            sed(ks)=sed(ks)+lay_sed(k,l,ks) !*lay_thick(k,l)
          enddo
        endif
        if(th>=active_thick)then 
          alay_sed(gid,1:totgrn)=sed(1:totgrn)
          exit lp
        endif
      enddo lp
    enddo

    ! Broadcast local stratigraphic layer globally
    call mpi_allreduce(mpi_in_place,alay_thick,snodes,mpi_double_precision,mpi_max,badlands_world,rc)
    do ks=1,totgrn
      tsed(1:snodes)=alay_sed(1:snodes,ks)
      call mpi_allreduce(tsed,gsed,snodes,mpi_double_precision,mpi_max,badlands_world,rc)
      alay_sed(1:snodes,ks)=gsed(1:snodes)
    enddo

    ! Define delaunay vertex associated active layer
    do k=1,dnodes
      nb=regPtNb(k)
      th=0.0
      sed=0.0
      ! Loop through the stratigraphic nodes associated to the vertex
      do l=1,nb
        gid=regPtID(k,l)
        th=th+alay_thick(gid)
        do ks=1,totgrn
          sed(ks)=sed(ks)+alay_sed(gid,ks)
        enddo
      enddo
      ! Define combined active layer for the vertex node
      if(nb>0)then
        alay_dthick(k)=th/nb
        do ks=1,totgrn
          alay_dsed(k,ks)=sed(ks)/nb
        enddo
      else
        alay_dthick(k)=0.0
        alay_dsed(k,1:totgrn)=0.0
      endif
    enddo

    return

  end subroutine buildActiveLayer
  ! =====================================================================================

end module strata_evol