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


module module_regions


!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4) :: NLAMMPS_REGIONS;

!   ************************************************************************************************


!   integer (kind=4), allocatable, dimension(:) :: CONFIG_ATOM_TYPE_DUPLICATE;  


!   ************************************************************************************************

    real (kind=8), allocatable, dimension(:,:) :: LAMMPS_REGION_ARGS, LAMMPS_REGION_DIMENSIONS;

!   real (kind=8), allocatable, dimension(:) :: ATOM_MASSES_DUPLICATE;

!   real (kind=8), allocatable, dimension(:,:) :: CONFIG_RI_DUPLICATE;

!   ************************************************************************************************

!   real (kind=8), dimension(1:3) :: CELL_AXIS_DUPLICATE, CELL_ANGDEG_DUPLICATE;

!   real (kind=8), dimension(1:3,1:3) :: PASSA_DUPLICATE, PASSB_DUPLICATE;

!   ************************************************************************************************

!   character (len=50) :: DUPLICATE_CHFILE_FORMAT;

!   character (len=250) :: DUPLICATE_CHNAME_FILE;

!   ************************************************************************************************

    character (len=250), allocatable, dimension(:) :: LAMMPS_REGION_NAME,  &
                                                      LAMMPS_REGION_STYLE;

!   character (len=20), allocatable, dimension(:) :: ATOM_LABEL_DUPLICATE;

!   character (len=250) :: CH_UNITS_STYLE_DUPLICATE, &
!                          CH_ATOM_STYLE_DUPLICATE;

end module module_regions
