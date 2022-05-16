!   ////////////////////////////////////////////////////////////////////////////////////////////////
!   //                                                                                            //
!   // Copyright (2022) Patrick A. Bonnaud                                                        //
!   //                                                                                            //
!   // This file is part of SLAB (Build Molecular Models for Molecular Simulations using          //
!   // Models of Molecules taken from a library).                                                 //
!   //                                                                                            //
!   // SLAB is free software; you can redistribute it and/or modify it under the terms            //
!   // of the GNU General Public License as published by the Free Software Foundation; either     //
!   // version 2 of the License, or (at your option) any later version.                           //
!   //                                                                                            //
!   // SLAB is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;          //
!   // without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  //
!   // See the GNU General Public License for more details.                                       //
!   //                                                                                            //
!   // You should have received a copy of the GNU General Public License along with this program. //
!   // If not, see <http://www.gnu.org/licenses/>.                                                //
!   //                                                                                            //
!   ////////////////////////////////////////////////////////////////////////////////////////////////


module module_simulation_box


!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    real (kind=8), dimension(1:3) :: LAMMPS_CELL_AXIS, LAMMPS_CELL_ANGDEG;  

    real (kind=8), dimension(1:6) :: LAMMPS_MATA;

    real (kind=8), dimension(1:3,1:3) :: LAMMPS_PASSA, LAMMPS_PASSB;

!   ************************************************************************************************

    character (len=250) :: LAMMPS_CREATE_BOX_RNAME;

!   character (len=250), allocatable, dimension(:) :: LAMMPS_REGION_NAME,  &
!                                                     LAMMPS_REGION_STYLE;


end module module_simulation_box
