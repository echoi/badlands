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
!       Filename:  HydrologyClass.f90
!
!    Description:  Store hydrological parameters
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles     
!
! =====================================================================================

module hydroUtil

  use ESMF

  implicit none

  ! Logical flag for regolith file
  logical::regofileflag

  ! Receivers nodes
  integer,dimension(:),allocatable::receivers

  ! Number of donors
  integer,dimension(:),allocatable::donorsList

  ! Index array
  integer,dimension(:),allocatable::indexArray

  ! Donors information array
  integer,dimension(:,:),allocatable::donorInfo

  ! Nodes ordering stack
  integer,dimension(:),allocatable::stackOrder

  ! List of baselevel nodes
  integer::baseNb
  integer,dimension(:),allocatable::baselist

  ! Discharge
  real(ESMF_KIND_R8),dimension(:),allocatable::discharge

  ! Maximum number of receiver nodes
  integer::maxrcvs

  ! Stream network accumulation minimal value
  real(ESMF_KIND_R8)::accu_thres

  ! Check if facies is required
  real(ESMF_KIND_R8)::faciesOn
  
  ! Strahler stream order
  integer,dimension(:),allocatable::strahler

  ! Cathcment ID
  integer::catchID
!   integer,dimension(:),allocatable::catchmentID
  integer,dimension(:),allocatable::subcatchmentID
  integer,dimension(:),allocatable::subcatchmentProc

  ! Simulation time parameters
  integer::Tforce
  real(ESMF_KIND_R8)::time_step,simulation_time,layer_time,udw_time,force_time
  real(ESMF_KIND_R8)::time_start,time_end,time_display,cpl1_time,cpl2_time
  real(ESMF_KIND_R8)::display_interval,layer_interval,CFL_diffusion
  
  ! Diffusion coefficient
  real(ESMF_KIND_R8),dimension(2)::Cdiffusion ! linear
  real(ESMF_KIND_R8),dimension(2)::Cdiffusion_d ! depth-dependent
  real(ESMF_KIND_R8)::Cdiff_m ! depth-dependent 
  real(ESMF_KIND_R8)::Cdiff_n ! depth-dependent 
  real(ESMF_KIND_R8),dimension(2)::Cdiffusion_nl ! non-linear
  real(ESMF_KIND_R8)::slope_critical 
  
  ! Surface erodibility coefficient
  integer::perosive
  real(ESMF_KIND_R8)::Cerodibility
  
  ! Stream Power Law coefficient
  real(ESMF_KIND_R8)::spl_m
  real(ESMF_KIND_R8)::spl_n

  ! Sefiment transport efficiency coefficient
  real(ESMF_KIND_R8)::Cefficiency

  ! Fraction of total load delivered to channels as bedload
  real(ESMF_KIND_R8)::Fracbed
  
  ! Stream Power Law coefficient
  real(ESMF_KIND_R8)::stl_m
  real(ESMF_KIND_R8)::stl_n

  ! Regolith formation
  real(ESMF_KIND_R8)::regoProd
  real(ESMF_KIND_R8)::regoDepth
  real(ESMF_KIND_R8)::soil_density
  real(ESMF_KIND_R8)::rock_density

  ! Capacity model coefficient
  real(ESMF_KIND_R8)::chan_exp
  real(ESMF_KIND_R8)::chan_width
  real(ESMF_KIND_R8)::bed_length
  real(ESMF_KIND_R8)::sed_length
  real(ESMF_KIND_R8)::stream_ero
  real(ESMF_KIND_R8)::bed_sed_interface

  ! New elevation change
  real(ESMF_KIND_R8),dimension(:),allocatable::newZ,spmZ,spmH

  ! Stratigraphic record
  integer::layerNb,layerID
  integer,parameter::faciestype=12
  real(ESMF_KIND_R8),dimension(:,:),allocatable::stratalZ
  real(ESMF_KIND_R8),dimension(:,:),allocatable::stratalFacies

  ! Sediment influx
  real(ESMF_KIND_R8),dimension(:),allocatable::Qs_in

  ! Subcatchement declaration type
  integer(ESMF_KIND_I4),dimension(:),allocatable::subcatchNb

  ! Network parallelisation stream declaration
  integer::localNodes
  integer,dimension(:),allocatable::sendprocID,rcvprocNb
  integer,dimension(:,:),allocatable::rcvsendID,rcvprocID
  integer,dimension(:),allocatable::localNodesGID

end module hydroUtil