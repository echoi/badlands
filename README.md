Badlands
=====

<div align="center">
    <img width=400 src="https://github.com/badlands-model/Badlands-doc/blob/master/riseofthephoenix.png" alt="Rise of the Phoenix" title="Example of Landscape evolution with Badlands"</img>
</div>

## Overview

**Ba**sin an**d** **Lan**dscape **D**ynamic**s** (**Badlands**) is a parallel TIN-based landscape evolution model, built to simulate topography development at various space and time scales. The model is presently capable of simulating hillslope processes (**linear** diffusion), fluvial incision (*'modified'* **SPL**:  erosion/transport/deposition), spatially and temporally varying geodynamic (horizontal + vertical displacements) and climatic forces which can be used to simulate changes in base level, as well as effects of climate changes or sea-level fluctuations. 

## The specs...

The model is mainly written in fortran and is based on the following characteristics: 
* The finite volume approach from Tucker et al. (2001) based on the dual Delaunay-Voronoi framework is used to solve the continuity equation explicitly, 
* Node ordering is perform efficiently based on the work from Braun & Willett (2013),
* A Hilbert Space-Filling Curve method algorithm (Zoltan) is used to partition the TIN-based surface into subdomains,
* Drainage network partitioning is generated through METIS library. 

## Community driven

This program is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with this program.  If not, see <http://www.gnu.org/licenses/lgpl-3.0.en.html>.

The code is conceived as an open-source project, and is an ideal tool for both **Research** and **Learning** purposes.

## Going further

A compilation of examples and some documentation related to the physics, assumptions and installation of Badlands can be found in the <a href='https://github.com/badlands-model/Badlands-doc'>Badlands-doc repository<a/>

