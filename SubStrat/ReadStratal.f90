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
!       Filename:  ReadStratal.f90
!
!    Description:  Reads stratigraphic grid and sediment parameters
!
!        Version:  1.0 
!        Created:  03/06/15 07:50:05
!        Revision:  none
!
!        Author:  Tristan Salles     
!
! =====================================================================================

module stratal_read

  use parallel
  use FoX_sax
  use topology
  use parameters
  use hydroUtil
  use FoX_common
  use stratal_class

  implicit none

  integer::sedn,stratn

  logical,save::in_Strata=.false.
  logical,save::in_strat=.false.
  logical,save::in_sdx=.false.
  logical,save::in_Sedi=.false.
  logical,save::in_sed=.false.
  logical,save::in_sedNb=.false.
  logical,save::in_dia=.false.
  logical,save::in_dens=.false.
  logical,save::in_diffa=.false.
  logical,save::in_diffm=.false.
  logical,save::in_strath=.false.
  logical,save::in_stratname=.false.
  logical,save::in_laynb=.false.
  logical,save::in_laythick=.false.
  logical,save::in_laytime=.false.

contains

  ! ============================================================================

  subroutine startDocument_handler

  end subroutine startDocument_handler
  ! =====================================================================================

  subroutine endDocument_handler

  end subroutine endDocument_handler
  ! =====================================================================================

  subroutine startElement_handler(namespaceURI,localname,name,atts)

    character(len=*),intent(in)::namespaceURI
    character(len=*),intent(in)::localname
    character(len=*),intent(in)::name
    type(dictionary_t),intent(in)::atts

    ! Sediment element
    if(name=='struct_sediment')in_Sedi=.true.
    if(in_Sedi) call SsedElement_handler(name)
    if (name=='sed')then
       in_sed=.true.
       sedn=sedn+1
    endif
    if(in_sed) call SsedfieldElement_handler(name)

    ! Rainfall element
    if(name=='struct_stratal')in_Strata=.true.
    if(in_Strata) call SstrataElement_handler(name)
    if (name=='strata')then
       in_strat=.true.
       stratn=stratn+1
    endif
    if(in_strat) call SstratafieldElement_handler(name)

  end subroutine startElement_handler
  ! =====================================================================================

  subroutine endElement_handler(namespaceURI,localname,name)

    character(len=*),intent(in)::namespaceURI
    character(len=*),intent(in)::localname
    character(len=*),intent(in)::name

    call EsedElement_handler(name)
    call EsedfieldElement_handler(name)
    call EstrataElement_handler(name)
    call EstratafieldElement_handler(name)

  end subroutine endElement_handler
  ! =====================================================================================

  subroutine characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_Sedi) call sed_characters_handler(chars)
    if(in_sed) call sedfield_characters_handler(chars)
    if(in_Strata) call strata_characters_handler(chars)
    if(in_strat) call stratafield_characters_handler(chars)

  end subroutine characters_handler
  ! =====================================================================================

  subroutine SsedElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='sed_nb') in_sedNb=.true.

  end subroutine SsedElement_handler
  ! =====================================================================================

  subroutine SsedfieldElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='dia') in_dia=.true.
    if(name=='dens') in_dens=.true.
    if(name=='diffa') in_diffa=.true.
    if(name=='diffm') in_diffm=.true.

  end subroutine SsedfieldElement_handler
  ! =====================================================================================

  subroutine SstrataElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='strat_dx') in_sdx=.true.
    if(name=='layer_time') in_laytime=.true.
    if(name=='active_lay_thick') in_laythick=.true.
    if(name=='layer_nb') in_laynb=.true.

  end subroutine SstrataElement_handler
  ! =====================================================================================

  subroutine SstratafieldElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='strat_name') in_stratname=.true.
    if(name=='strat_thick') in_strath=.true.

  end subroutine SstratafieldElement_handler
  ! =====================================================================================

  subroutine EsedElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='sed_nb') in_sedNb=.false.

  end subroutine EsedElement_handler
  ! =====================================================================================

  subroutine EsedfieldElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='dia') in_dia=.false.
    if(name=='dens') in_dens=.false.
    if(name=='diffa') in_diffa=.false.
    if(name=='diffm') in_diffm=.false.

  end subroutine EsedfieldElement_handler
  ! =====================================================================================

  subroutine EstrataElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='strat_dx') in_sdx=.false.
    if(name=='layer_time') in_laytime=.false.
    if(name=='active_lay_thick') in_laythick=.false.
    if(name=='layer_nb') in_laynb=.false.

  end subroutine EstrataElement_handler
  ! =====================================================================================

  subroutine EstratafieldElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='strat_name') in_stratname=.false.
    if(name=='strat_thick') in_strath=.false.

  end subroutine EstratafieldElement_handler
  ! =====================================================================================

  subroutine sed_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_sedNb)then
       call rts(chars,totgrn)
       allocate(sediments(totgrn))
    endif

  end subroutine sed_characters_handler
  ! =====================================================================================

  subroutine sedfield_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_dia) then
       call rts(chars,sediments(sedn)%dia)
    endif
    if(in_dens)then
       call rts(chars,sediments(sedn)%dens)
    endif
    if(in_diffa) then
       call rts(chars,sediments(sedn)%diffa)
    endif
    if(in_diffm)then
       call rts(chars,sediments(sedn)%diffm) 
    endif

  end subroutine sedfield_characters_handler
  ! =====================================================================================

  subroutine strata_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_sdx)then
       call rts(chars,sdx)
    endif
    if(in_laytime)then
       call rts(chars,time_layer)
    endif
    if(in_laythick)then
       call rts(chars,active_thick)
    endif
    if(in_laynb)then
       call rts(chars,ilay)
       allocate(ilayh(ilay))
       allocate(flayh(ilay))
       ilayh=0.0_8
       flayh=''
    endif

  end subroutine strata_characters_handler
  ! =====================================================================================

  subroutine stratafield_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_stratname)then
       flayh(stratn)=chars
    endif
    if(in_strath)then
       call rts(chars,ilayh(stratn))
    endif

  end subroutine stratafield_characters_handler
  ! =====================================================================================

  subroutine stratal_parser

    type(xml_t)::xf

    sedn=0
    ilay=0
    stratn=0
    totgrn=0
    sdx=dx*0.1_8
    active_thick=1.0_8
    time_layer=display_interval

    ! Open file
    call open_xml_file(xf,xmlfile,rc)
    if(rc/=0)then
       print*,'Failed to open XmL input file'
       call mpi_finalize(rc)
    endif

    ! Parser 
    call parse(xf, &
         startDocument_handler=startDocument_handler, &
         startElement_handler=startElement_handler, &
         endElement_handler=endElement_handler, &
         characters_handler=characters_handler, &
         endDocument_handler=endDocument_handler)

    ! Close file
    call close_xml_t(xf)

    if(time_layer>display_interval) time_layer=display_interval

  end subroutine stratal_parser
  ! =====================================================================================

end module stratal_read
