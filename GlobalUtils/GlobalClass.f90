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
!       Filename:  GlobalClass.f90
!
!    Description:  Store global application parameters
!
!        Version:  1.0
!        Created:  11/02/15 05:05:05
!        Revision:  none
!
!        Author:  Tristan Salles     
!
! =====================================================================================

module parallel

#include "mpif.h"

  integer,parameter::badlands_world=mpi_comm_world

end module parallel

module parameters

  use ESMF

  implicit none

  logical::updateSPM_elevation,oceanFlag,geodynamicFlag,restartFlag,udwFlag

  ! Step for restarting simulation
  integer::restartStep,restartPet

  ! ESMF Virtual Machine Class
  type(ESMF_VM),save::vm

  ! Persistent Execution Threads (ID and total number)
  integer::pet_id,npets,rc 

  ! Grid / Coupler Component Names
  character(len=128)::ocean,spm,earth,coupler1,coupler2,coupler3,coupler4

  ! Grid Component Declaration
  type(ESMF_GridComp),save::oceanComp,spmComp,earthComp

  ! Coupler Component Declaration
  type(ESMF_CplComp),save::cpl1,cpl2,cpl3,cpl4

  ! Import and Export States Declaration
  type(ESMF_State),save::SstateExp,SstateImp
  type(ESMF_State),save::OstateExp,OstateImp
  type(ESMF_State),save::EstateExp,EstateImp

  ! Instantiate a clock, a calendar, and timesteps
  type(ESMF_Clock),save::clock
  type(ESMF_Calendar),save::GeologyCalendar
  type(ESMF_TimeInterval),save::timeStep
  type(ESMF_Time),save::startTime,stopTime,currentTime

  ! Simulation directories and file names
  character(len=128)::regularfile,outdir,xmlfile
  character(len=128)::outdir1,outdir2,runfiles,outputs

  ! Underworld files and folders
  character(len=128)::outdir3,fudw,fudisp,maestro

  ! Regolith file name
  character(len=128)::regofile

  ! Regular structured grid nodes on X/Y directions 
  integer(ESMF_KIND_I4)::nx,ny

  ! Regular structured grid resolution 
  real(ESMF_KIND_R8)::dx

  ! Regular structured grid region extensions 
  real(ESMF_KIND_R8)::minx,miny,maxx,maxy

  ! Regular structured grid X,Y,Z coordinates 
  real(ESMF_KIND_R8),dimension(:),allocatable::coordX,coordY,coordZ

  ! Regular structured grid X,Y,Z coordinates with added cells 
  ! on the edges and vertical displacement
  real(ESMF_KIND_R8),dimension(:),allocatable::rcoordX,rcoordY,rcoordZ

  ! Conforming Delaunay Triangulation (CDT) parameters
  integer(ESMF_KIND_I4)::tnodes,dnodes,delem,dedge

  ! Output parameters for CDT and drainage
  integer(ESMF_KIND_I4)::delemo,delemoo,drainOde
  integer(ESMF_KIND_I4),dimension(:),allocatable::outelem,outnode,dglbID,doutelem,doutnode

  ! Voronoi Diagram of the CDT parameters
  integer(ESMF_KIND_I4)::vnodes,vedge,vcellIN,velemIN

  ! Conforming Delaunay Triangulation X,Y,Z coordinates and vertical displacement
  real(ESMF_KIND_R8),dimension(:),allocatable::tcoordX,tcoordY,tcoordZ

  ! Voronoi point triangle face pt ID
  integer(ESMF_KIND_I4),dimension(:,:),allocatable::vorDel

  ! Voronoi Diagram X,Y,Z coordinates 
  real(ESMF_KIND_R8),dimension(:),allocatable::vcoordX,vcoordY 

  ! Regular structured elements (square cells) with added cells on the edges 
  integer(ESMF_KIND_I4),dimension(:,:),allocatable::relmt

  ! Number of grid nodes owned by the processor on a given pet 
  integer(ESMF_KIND_I4)::sOwnedNode,uOwnedNode

  ! Partition processor ID for unstructured grid
  integer(ESMF_KIND_I4),dimension(:),allocatable::snodeID,selemID,sownedID,snodeLID

  ! Partition processor ID for unstructured grid
  integer(ESMF_KIND_I4),dimension(:),allocatable::unodeID,unodeIDI,uelemID,uownEID,uownedID,unodeLID,delemID

  ! Number of nodes and elements on each partition mesh
  integer(ESMF_KIND_I4)::upartN,upartE,spartN,spartE !,upartI

  ! Redistribruted array of arbitrary distrinuted ESMF mesh indexes
  integer(ESMF_KIND_I4),dimension(:),allocatable::earthIDs,oceanIDs,spmIDs

  ! Voronoi Cell Declaration Type
  type voronoi
     integer(ESMF_KIND_I4)::vertexNb
     integer(ESMF_KIND_I4)::border
     integer(ESMF_KIND_I4)::btype
     integer(ESMF_KIND_I4)::bpoint
     integer(ESMF_KIND_I4),dimension(20)::vertexID
     real(ESMF_KIND_R8)::perimeter
     real(ESMF_KIND_R8)::area
  end type voronoi
  type(voronoi),dimension(:),allocatable::voronoiCell

  ! Conforming Delaunay Triangulation Declaration Type
  type delaunay
     integer(ESMF_KIND_I4)::ngbNb
     integer(ESMF_KIND_I4)::sortedHullNb
     integer(ESMF_KIND_I4),dimension(20)::ngbID
     integer(ESMF_KIND_I4),dimension(20)::sortedHull
     real(ESMF_KIND_R8)::slope
     real(ESMF_KIND_R8),dimension(20)::voronoi_edge
     real(ESMF_KIND_R8),dimension(20)::distance
     real(ESMF_KIND_R8),dimension(20)::weight
  end type delaunay
  type(delaunay),dimension(:),allocatable::delaunayVertex

  ! Number of seconds per year
  real(ESMF_KIND_R8),parameter::secyear=31536000.0

  ! Maximum filling algorithm height
  real(ESMF_KIND_R8)::fh

  ! Infiltration evaporation percentage for water in lakes
  real(ESMF_KIND_R8)::infiltration_evaporation

contains
  
  ! =====================================================================================
  subroutine term_command(cmds)

    logical(4)::result

    character(len=128)::cmds

    result = .false.

    ! INTEL FORTRAN COMPILER
    !        result = systemqq( cmds )
    ! GNU FORTRAN COMPILER
    call system(cmds)

    return

  end subroutine term_command
  ! =====================================================================================

end module parameters