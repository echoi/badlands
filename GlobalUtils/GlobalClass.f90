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

  use parallel

  implicit none

  logical::updateSPM_elevation,oceanFlag,geodynamicFlag,restartFlag,udwFlag,update3d

  ! Step for restarting simulation
  integer::restartStep,restartPet

  ! Persistent Execution Threads (ID and total number)
  integer::pet_id,npets,rc 

  ! Grid / Coupler Component Names
!   character(len=128)::ocean,spm,earth

  ! Number of total grains
  integer::totgrn

  ! Simulation directories and file names
  character(len=128)::regularfile,outdir,xmlfile
  character(len=128)::outdir1,outdir2,runfiles,outputs

  ! Underworld files and folders
  character(len=128)::outdir3,fudw,fudisp,maestro

  ! Regolith file name
  character(len=128)::regofile

  ! Regular structured grid nodes on X/Y directions 
  integer::nx,ny

  ! Regular structured grid resolution 
  real(kind=8)::dx

  ! Regular structured grid region extensions 
  real(kind=8)::minx,miny,maxx,maxy

  ! Regular structured grid X,Y,Z coordinates 
  real(kind=8),dimension(:),allocatable::coordX,coordY,coordZ

  ! Regular structured grid X,Y,Z coordinates with added cells 
  ! on the edges and vertical displacement
  real(kind=8),dimension(:),allocatable::rcoordX,rcoordY,rcoordZ

  ! Regular structured grid bilinear
  real,dimension(:),allocatable::bilinearX,bilinearY
  real,dimension(:,:),allocatable::bilinearV,bilinearHx,bilinearHy

  ! Conforming Delaunay Triangulation (CDT) parameters
  integer::tnodes,dnodes,delem,dedge

  ! Output parameters for CDT and drainage
  integer::delemo,delemoo,drainOde
  integer,dimension(:),allocatable::outelem,outnode,dglbID,doutelem,doutnode

  ! Voronoi Diagram of the CDT parameters
  integer::vnodes,vedge,vcellIN,velemIN

  ! Conforming Delaunay Triangulation X,Y,Z coordinates and vertical displacement
  real(kind=8),dimension(:),allocatable::tcoordX,tcoordY,tcoordZ

  ! Voronoi point triangle face pt ID
  integer,dimension(:,:),allocatable::vorDel

  ! Voronoi Diagram X,Y,Z coordinates 
  real(kind=8),dimension(:),allocatable::vcoordX,vcoordY 

  ! Regular structured elements (square cells) with added cells on the edges 
  integer,dimension(:,:),allocatable::relmt

  ! Number of grid nodes owned by the processor on a given pet 
  integer::sOwnedNode,uOwnedNode

  ! Partition processor ID for unstructured grid
  integer,dimension(:),allocatable::snodeID,selemID,sownedID,snodeLID

  ! Partition processor ID for unstructured grid
  integer,dimension(:),allocatable::unodeID,unodeIDI,uelemID,uownEID,uownedID,unodeLID,delemID

  ! Number of nodes and elements on each partition mesh
  integer::upartN,upartE,spartN,spartE !,upartI

  ! Redistribruted array of arbitrary distrinuted mesh indexes
  integer,dimension(:),allocatable::earthIDs,oceanIDs,spmIDs

  ! Voronoi Cell Declaration Type
  type voronoi
     integer::vertexNb
     integer::border
     integer::btype
     integer::bpoint
     integer,dimension(20)::vertexID
     real(kind=8)::perimeter
     real(kind=8)::area
  end type voronoi
  type(voronoi),dimension(:),allocatable::voronoiCell

  ! Conforming Delaunay Triangulation Declaration Type
  type delaunay
     integer::ngbNb
     integer::sortedHullNb
     integer,dimension(20)::ngbID
     integer,dimension(20)::sortedHull
     real(kind=8)::slope
     real(kind=8),dimension(20)::voronoi_edge
     real(kind=8),dimension(20)::distance
     real(kind=8),dimension(20)::weight
  end type delaunay
  type(delaunay),dimension(:),allocatable::delaunayVertex

  ! Number of seconds per year
  real(kind=8),parameter::secyear=31536000.0

  ! Maximum filling algorithm height
  real(kind=8)::fh

  ! Infiltration evaporation percentage for water in lakes
  real(kind=8)::infiltration_evaporation

contains
  
  ! =====================================================================================
  
  subroutine term_command(cmds)

    logical(4)::result

    character(len=128)::cmds

    result = .false.

    ! INTEL FORTRAN COMPILER
    ! result = systemqq( cmds )
    ! GNU FORTRAN COMPILER
    call system(cmds)

    return

  end subroutine term_command
  ! =====================================================================================

end module parameters