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
!       Filename:  ForcesClass.f90
!
!    Description:  Define forces parameters such as rain, sea-level and tectonic evolution
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles     
!
! =====================================================================================

module external_forces

  use parallel  
  use hydroUtil
  use parameters

  implicit none

  logical::disp3d

  ! Number of rain events 
  integer::rain_event

  ! Rain event files
  character(len=128),dimension(:),allocatable::frainmap

  ! Rain event starting and ending time
  real(kind=8),dimension(:),allocatable::rain_tstart,rain_tend

  ! Precipitation
  real(kind=8),dimension(:),allocatable::precipitation

  ! Sea-level fluctuation file
  character(len=128)::seafile

  !> Sea level fluctuations type
  type sea_fluc
    ! Flag for sea level fluctuations
    logical::sealevel
    ! Number of sea fluctuation events
    integer::event
    ! Old sea level value
    real(kind=8)::last_sea
    ! Actual sea level value
    real(kind=8)::actual_sea
  end type sea_fluc
  type(sea_fluc)::gsea

  ! Sea-level elevation data
  real(kind=8),dimension(:,:),allocatable::sealvl

  !> Sea-level parameters
  type sl_par
    ! Lower sea-level elevation
    real(kind=8)::sea1
    ! Upper sea-level elevation
    real(kind=8)::sea2
    ! Lower sea-level time
    real(kind=8)::time1
    ! Upper sea-level time
    real(kind=8)::time2
    ! Simulation current time
    real(kind=8)::tsim
  end type sl_par

  !> Geodynamic type
  type geodyn
    ! Actual displacement event
    integer::actual
    ! Displacements event numbers
    integer::event
    ! Minimum distance between points for merging
    real(kind=8)::mindist
    ! Time from last displacement call
    real(kind=8)::lastdisp
    ! Time elapsed between current and previous displacement time
    real(kind=8)::disptime
  end type geodyn
  type(geodyn)::disp,vdisp

  ! Geodynamic file names
  character(len=128),dimension(:),allocatable::fgeodyn

  ! Geodynamic event time
  real(kind=8),dimension(:,:),allocatable::disp_time,vdisp_time

  ! Geodynamic event time
  integer,dimension(:),allocatable::disp_fill,vdisp_fill

  ! Regular grid field arrays
  real(kind=8),dimension(:),allocatable::rtectoZ,rvertDisp,rainVal,rhxDisp,rhyDisp
  real(kind=8),dimension(:,:),allocatable::rDisp

  ! Vertical displacement rate on unstructured grid
  real(kind=8),dimension(:),allocatable::tvertDisp

contains

  ! =====================================================================================

  subroutine read_sealevel_file

    integer::iu,nbsea,k,i,i2,ios
    
    character(len=128)::line

    iu=42
    open(iu,file=seafile,status="old",action="read",iostat=rc)
    if(rc/=0)then
      if(pet_id==0)print*,'ERROR: failed to open sealevel file'
      call mpi_finalize(rc)
    endif

    ! Determine total number of sea level records
    nbsea=0
    do 
      read(iu,*,iostat=rc,end=30)line
       nbsea=nbsea+1
    enddo 
30 continue    
    rewind(iu)

    ! Define the number of events and allocate array
    gsea%event=nbsea
    if(allocated(sealvl)) deallocate(sealvl)
    allocate(sealvl(2,nbsea))
    ! Read sea-level
    do k=1,gsea%event
      read(iu,'(a128)',iostat=ios,end=60)line
      i2=len_trim(line)
      if(line(i2:i2)==char(13))then
        i2=i2-1
      endif
      read(line(1:i2),*,iostat=ios)(sealvl(i,k),i=1,2)
      if(ios/=0)then
        if(pet_id==0)print*,'ERROR: reading sea level file at line:',k
        call mpi_finalize(rc)
      endif
      if(k>1)then
        if(sealvl(1,k)<=sealvl(1,k-1))then
          if(pet_id==0)print*,'ERROR: in sea level file value at line',k
          if(pet_id==0)print*,'the time is not defined in increasing order'
          call mpi_finalize(rc)
        endif
      endif
    enddo
60 close(iu)

    gsea%last_sea=0.0
    gsea%actual_sea=0.0

    ! Update the sea level value
    call eustatism

  end subroutine read_sealevel_file
  ! =====================================================================================

  subroutine eustatism

    integer::i
    type(sl_par)::sl_param

    gsea%last_sea=gsea%actual_sea
    i=1
    ! Find time for the considered event
    do
      if(i>gsea%event)exit
      if(simulation_time<sealvl(1,i))exit
      sl_param%time1=sealvl(1,i)
      sl_param%sea1=sealvl(2,i)
      i=i+1
    enddo

    ! Linear interpolation of sea level elevation based on current time
    if(i<=gsea%event)then
      sl_param%time2=sealvl(1,i)
      sl_param%sea2=sealvl(2,i)
      sl_param%tsim=simulation_time
      gsea%actual_sea=sealvl_interpolation(sl_param)
    ! If event greater than max defined allocate with last value
    else
      gsea%actual_sea=sl_param%sea1
      sl_param%time2=sealvl(1,i-1)
    endif

  end subroutine eustatism
  ! =====================================================================================
  
  function sealvl_interpolation(sl_param) result(slinterpol)

    type(sl_par)::sl_param
    real(kind=8)::slinterpol

    slinterpol=(sl_param%tsim-sl_param%time1)/(sl_param%time2-sl_param%time1)
    slinterpol=slinterpol*(sl_param%sea2-sl_param%sea1)+sl_param%sea1

  end function sealvl_interpolation
  ! =====================================================================================

end module external_forces