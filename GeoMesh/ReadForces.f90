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

  use parallel
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
  logical,save::in_3Dfield=.false.
  logical,save::in_mrgd=.false.
  logical,save::in_DispNb=.false.
  logical,save::in_DispFile=.false.
  logical,save::in_DispET=.false.
  logical,save::in_DispST=.false.
  logical,save::in_disp=.false.
  logical,save::in_RainGrid=.false.
  logical,save::in_RainNb=.false.
  logical,save::in_RainMFile=.false.
  logical,save::in_RainMVal=.false.
  logical,save::in_RainET=.false.
  logical,save::in_RainST=.false.
  logical,save::in_rain=.false.
  logical,save::in_Ice=.false.
  logical,save::in_icedx=.false.
  logical,save::in_kernel=.false.
  logical,save::in_icedt=.false.
  logical,save::in_icedeform=.false.
  logical,save::in_iceslide=.false.
  logical,save::in_iceZsld=.false.
  logical,save::in_iceTstep=.false.
  logical,save::in_iceero=.false.
  logical,save::in_m1=.false.
  logical,save::in_m2=.false.
  logical,save::in_ELAfile=.false.
  logical,save::in_flex=.false.
  logical,save::in_flexdt=.false.
  logical,save::in_flexdx=.false.
  logical,save::in_flexsedth=.false.
  logical,save::in_flexsedporo=.false.
  logical,save::in_flexseddens=.false.
  logical,save::in_flexdecay=.false.
  logical,save::in_flexthick=.false.
  logical,save::in_mantledens=.false.

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

    ! Ice element
    if(name=='struct_ice')in_Ice=.true.
    if(in_Ice) call SiceElement_handler(name)

    ! Flexure element
    if(name=='struct_flex')in_flex=.true.
    if(in_flex) call SflexElement_handler(name)

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
    call EiceElement_handler(name)
    call EflexElement_handler(name)

  end subroutine endElement_handler
  ! =====================================================================================

  subroutine characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_Sea) call sea_characters_handler(chars)
    if(in_DispF) call vdisp_characters_handler(chars)
    if(in_disp) call disp_characters_handler(chars)
    if(in_RainGrid) call rainmap_characters_handler(chars)
    if(in_rain) call rainfield_characters_handler(chars)
    if(in_Ice) call ice_characters_handler(chars)
    if(in_flex) call flex_characters_handler(chars)

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
    if(name=='rain_value') in_RainMVal=.true.
    if(name=='rain_start') in_RainST=.true.
    if(name=='rain_end') in_RainET=.true.

  end subroutine SrainfieldElement_handler
  ! =====================================================================================

  subroutine SseaElement_handler(name)

    character(len=*), intent(in) :: name

    if (name=='ocean_file') in_SeaFile=.true.

  end subroutine SseaElement_handler
  ! =====================================================================================

  subroutine SiceElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='ice_dx') in_icedx=.true.
    if(name=='ice_dt') in_icedt=.true.
    if(name=='gauss_kernel') in_kernel=.true.
    if(name=='deformCst') in_icedeform=.true.
    if(name=='slideCst') in_iceslide=.true.
    if(name=='slideZ_ELA') in_iceZsld=.true.
    if(name=='ice_Tstep') in_iceTstep=.true.
    if(name=='ELA_file') in_ELAfile=.true.
    if(name=='melt1_ELA') in_m1=.true.
    if(name=='melt2_ELA') in_m2=.true.
    if(name=='ice_ero') in_iceero=.true.

  end subroutine SiceElement_handler
  ! =====================================================================================

  subroutine SflexElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='flex_dx') in_flexdx=.true.
    if(name=='flex_dt') in_flexdt=.true.
    if(name=='sed_dens') in_flexseddens=.true.
    if(name=='sed_poro') in_flexsedporo=.true.
    if(name=='comp_decay') in_flexdecay=.true.
    if(name=='mantle_dens') in_mantledens=.true.
    if(name=='flex_thick') in_flexthick=.true.

  end subroutine SflexElement_handler
  ! =====================================================================================

  subroutine SvdispElement_handler(name)

    character(len=*),intent(in)::name

    if (name=='fields_3D') in_3Dfield=.true.
    if (name=='merge_dist') in_mrgd=.true.
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
    if(name=='rain_value') in_RainMVal=.false.
    if(name=='rain_start') in_RainST=.false.
    if(name=='rain_end') in_RainET=.false.

  end subroutine ErainfieldElement_handler
  ! =====================================================================================

  subroutine EseaElement_handler(name)

    character(len=*), intent(in) :: name

    if (name=='ocean_file') in_SeaFile=.false.

  end subroutine EseaElement_handler
  ! =====================================================================================

  subroutine EiceElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='ice_dx') in_icedx=.false.
    if(name=='ice_dt') in_icedt=.false.
    if(name=='gauss_kernel') in_kernel=.false.
    if(name=='deformCst') in_icedeform=.false.
    if(name=='slideCst') in_iceslide=.false.
    if(name=='slideZ_ELA') in_iceZsld=.false.
    if(name=='ice_Tstep') in_iceTstep=.false.
    if(name=='ELA_file') in_ELAfile=.false.
    if(name=='melt1_ELA') in_m1=.false.
    if(name=='melt2_ELA') in_m2=.false.
    if(name=='ice_ero') in_iceero=.false.

  end subroutine EiceElement_handler
  ! =====================================================================================

  subroutine EflexElement_handler(name)

    character(len=*),intent(in)::name

    if(name=='flex_dx') in_flexdx=.false.
    if(name=='flex_dt') in_flexdt=.false.
    if(name=='sed_dens') in_flexseddens=.false.
    if(name=='comp_decay') in_flexdecay=.false.
    if(name=='sed_poro') in_flexsedporo=.false.
    if(name=='mantle_dens') in_mantledens=.false.
    if(name=='flex_thick') in_flexthick=.false.

  end subroutine EflexElement_handler
  ! =====================================================================================

  subroutine EvdispElement_handler(name)

    character(len=*),intent(in)::name

    if (name=='fields_3D') in_3Dfield=.false.
    if (name=='merge_dist') in_mrgd=.false.
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
      allocate(frainval(rain_event))
      frainval=-1
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
    if(in_RainMVal)then
      call rts(chars,frainval(rainn))
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

  subroutine ice_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_icedx)then
      call rts(chars,ice_dx)
    elseif(in_icedt)then
      call rts(chars,ice_dt)
    elseif(in_kernel)then
      call rts(chars,gausskernelSize)
    elseif(in_icedeform)then
      call rts(chars,ice_deform)
    elseif(in_iceslide)then
      call rts(chars,ice_slide)
    elseif(in_iceZsld)then
      call rts(chars,ice_Zsld)
    elseif(in_iceTstep)then
      call rts(chars,ice_Tstep)
    elseif(in_ELAfile)then
      elafile=chars
      gela%ela=.true.
    elseif(in_m1)then
      call rts(chars,ice_m1)
    elseif(in_m2)then
      call rts(chars,ice_m2)
    elseif(in_iceero)then
      call rts(chars,IceEro)
    endif

  end subroutine ice_characters_handler
  ! =====================================================================================

  subroutine flex_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_flexdx)then
      flexure=.true.
      call rts(chars,flex_dx)
    elseif(in_flexdt)then
      call rts(chars,flex_dt)
    elseif(in_flexseddens)then
      call rts(chars,mean_sediment_density)
    elseif(in_flexdecay)then
      call rts(chars,comp_decay)
    elseif(in_flexsedporo)then
      call rts(chars,comp_poro)
    elseif(in_mantledens)then
      call rts(chars,mean_mantle_density)
    elseif(in_flexthick)then
      call rts(chars,flex_thick)
    endif

  end subroutine flex_characters_handler
  ! =====================================================================================

  subroutine vdisp_characters_handler(chars)

    character(len=*),intent(in)::chars

    if(in_3Dfield)then
      disp3d=.true.
    elseif(in_mrgd)then
      call rts(chars,disp%mindist)
    elseif(in_DispNb)then
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
    real(kind=8)::temp

    disp3d=.false.
    rainn=0
    dispn=0
    rain_event=0
    if(.not.udwFlag)disp%event=0
    gsea%sealevel=.false.
    gsea%actual_sea=0.
    disp%mindist=0.0
    gausskernelSize=2
    ice_dx=0.0
    gela%ela=.false.
    gela%actual_ela=2000.
    ice_m1=0.01
    ice_m2=0.01
    ice_dt=0.1
    ice_Tstep=0.
    ice_Zsld=-200.
    IceEro=0.001
    flex_dx=0.0
    flex_thick=10000.
    comp_decay=4.0
    comp_poro=0.
    flexure=.false.

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
    if(gsea%sealevel) call read_sealevel_file
    if(gela%ela) call read_ELA_file
    if(ice_m1<ice_m2)then
      temp=ice_m1
      ice_m1=ice_m2
      ice_m2=temp
    endif
    if(ice_dx>0.)then
      cpl3_time=time_start
    else
      cpl3_time=time_end+1000.0
    endif

    if(flexure)then
      cpl4_time=time_start
      flex_rigid=70.e9/(12.0*0.75**2.)
      flex_rigid=flex_rigid*flex_thick**3.
    else
      cpl4_time=time_end+1000.0
    endif

  end subroutine forces_parser
  ! =====================================================================================

end module readforces
! =====================================================================================
