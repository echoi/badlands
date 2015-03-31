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
!       Filename:  ReadForces.f90
!
!    Description:  Read the XmL input file
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles     
!
! =====================================================================================

module readforces

  use FoX_sax
  use topology
  use parameters
  use hydroUtil
  use FoX_common
  use external_forces

  implicit none
  
  integer::dispn,rainn

  logical,save::in_Sea=.false.
  logical,save::in_SeaFile=.false.
  logical,save::in_DispF=.false.
  logical,save::in_DispNb=.false.
  logical,save::in_DispFile=.false.
  logical,save::in_DispET=.false.
  logical,save::in_DispST=.false.
  logical,save::in_disp=.false.
  logical,save::in_RainGrid=.false.
  logical,save::in_RainNb=.false.
  logical,save::in_RainMFile=.false.
  logical,save::in_RainET=.false.
  logical,save::in_RainST=.false.
  logical,save::in_rain=.false.

contains

  ! =====================================================================================
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

    ! Ocean element
    if(name=='struct_ocean')in_Sea=.true.
    if(in_Sea) call SseaElement_handler(name)

    ! Geodynamic element
    if(name=='struct_geodyn') in_DispF=.true.
    if(in_DispF) call SvdispElement_handler(name)
    if(name=='disp')then
      in_disp=.true.
      dispn=dispn+1
    endif
    if(in_disp) call SdispElement_handler(name)

    ! Rainfall element
    if(name=='struct_rainfall')in_RainGrid=.true.
    if(in_RainGrid) call SrainmapElement_handler(name)
    if (name=='rain')then
      in_rain=.true.
      rainn=rainn+1
    endif
    if(in_rain) call SrainfieldElement_handler(name)

  end subroutine startElement_handler
  ! =====================================================================================
  subroutine endElement_handler(namespaceURI,localname,name)

    character(len=*),intent(in)::namespaceURI
    character(len=*),intent(in)::localname
    character(len=*),intent(in)::name

    call EseaElement_handler(name)
    call EvdispElement_handler(name)
    call EdispElement_handler(name)
    call ErainmapElement_handler(name)
    call ErainfieldElement_handler(name)

  end subroutine endElement_handler
  ! =====================================================================================
  subroutine characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_Sea) call sea_characters_handler(chars)
    if(in_DispF) call vdisp_characters_handler(chars)
    if(in_disp) call disp_characters_handler(chars)
    if(in_RainGrid) call rainmap_characters_handler(chars)
    if(in_rain) call rainfield_characters_handler(chars)
    
  end subroutine characters_handler
  ! =====================================================================================
  subroutine SrainmapElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='rain_nb') in_RainNb=.true.

  end subroutine SrainmapElement_handler
  ! =====================================================================================
  subroutine SrainfieldElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='rain_file') in_RainMFile=.true.
    if(name=='rain_start') in_RainST=.true.
    if(name=='rain_end') in_RainET=.true.

  end subroutine SrainfieldElement_handler
  ! =====================================================================================
  subroutine SseaElement_handler(name)

    character(len=*), intent(in) :: name

    if (name=='ocean_file') in_SeaFile=.true.

  end subroutine SseaElement_handler
  ! =====================================================================================
  subroutine SvdispElement_handler(name)

    character(len=*),intent(in)::name

    if (name=='interval_nb') in_DispNb=.true.

  end subroutine SvdispElement_handler
  ! =====================================================================================
  subroutine SdispElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='disp_file') in_DispFile=.true.
    if(name=='disp_start') in_DispST=.true.
    if(name=='disp_end') in_DispET=.true.

  end subroutine SdispElement_handler
  ! =====================================================================================
  subroutine ErainmapElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='rain_nb') in_RainNb=.false.

  end subroutine ErainmapElement_handler
  ! =====================================================================================
  subroutine ErainfieldElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='rain_file') in_RainMFile=.false.
    if(name=='rain_start') in_RainST=.false.
    if(name=='rain_end') in_RainET=.false.

  end subroutine ErainfieldElement_handler
  ! =====================================================================================
  subroutine EseaElement_handler(name)

    character(len=*), intent(in) :: name

    if (name=='ocean_file') in_SeaFile=.false.

  end subroutine EseaElement_handler
  ! =====================================================================================
  subroutine EvdispElement_handler(name)

    character(len=*),intent(in)::name

    if (name=='interval_nb') in_DispNb=.false.

  end subroutine EvdispElement_handler
  ! =====================================================================================
  subroutine EdispElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='disp_file') in_DispFile=.false.
    if(name=='disp_start') in_DispST=.false.
    if(name=='disp_end') in_DispET=.false.

  end subroutine EdispElement_handler
  ! =====================================================================================
  subroutine rainmap_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_RainNb)then
      call rts(chars,rain_event)
      allocate(frainmap(rain_event))
      allocate(rain_tend(rain_event))
      allocate(rain_tstart(rain_event))
    endif

  end subroutine rainmap_characters_handler
  ! =====================================================================================
  subroutine rainfield_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_RainMFile)then
      frainmap(rainn)=chars
    endif
    if(in_RainST) then
      call rts(chars,rain_tstart(rainn))
    endif
    if(in_RainET)then
      call rts(chars,rain_tend(rainn))
    endif

  end subroutine rainfield_characters_handler
  ! =====================================================================================
  subroutine sea_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_SeaFile)then
      seafile=chars
      gsea%sealevel=.true.
    endif

    end subroutine sea_characters_handler
  ! =====================================================================================
  subroutine vdisp_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_DispNb)then
      call rts(chars,disp%event)
      allocate(fgeodyn(disp%event))
      allocate(disp_time(disp%event,2))
    endif

  end subroutine vdisp_characters_handler
  ! =====================================================================================
  subroutine disp_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_DispFile) fgeodyn(dispn)=chars
    if(in_DispST)then
      call rts(chars,disp_time(dispn,1))
    endif
    if(in_DispET)then
      call rts(chars,disp_time(dispn,2))
    endif

  end subroutine disp_characters_handler
  ! =====================================================================================
  subroutine forces_parser

    type(xml_t)::xf

    rainn=0
    dispn=0
    rain_event=0
    if(.not.udwFlag)disp%event=0
    gsea%sealevel=.false.
    gsea%actual_sea=0.

    ! Open file
    call open_xml_file(xf,xmlfile,rc)
    if(rc/=0)then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_FILE_OPEN, &
        msg="Failed to open namelist file 'XmL_input_file'", &
        line=__LINE__,file=__FILE__)
      call ESMF_Finalize(endflag=ESMF_END_ABORT)
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

    if(gsea%sealevel) call read_sealevel_file

  end subroutine forces_parser
  ! =====================================================================================

end module readforces
! =====================================================================================
