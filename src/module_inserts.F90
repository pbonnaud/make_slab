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

module module_inserts


!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4) :: NLAMMPS_INSERTS, ILAMMPS_INSERTS;

    integer (kind=4) :: NFILE_EXT;

    integer (kind=4) :: iglobal_file, ifile_info;

    integer (kind=4) :: LAMMPS_INSERT_IPOTENTIAL_CLASS2, LAMMPS_INSERT_ITETHER;

!   ************************************************************************************************

    integer (kind=4), allocatable, dimension(:) :: LAMMPS_INSERT_FLAG;

    integer (kind=4), allocatable, dimension(:) :: LAMMPS_INSERT_RANDOM_SEED;

    integer (kind=4), allocatable, dimension(:) :: LAMMPS_INSERT_REGION_IRANK;

    integer (kind=4), allocatable, dimension(:) :: LAMMPS_INSERT_NMOLECULE;

    integer (kind=4), allocatable, dimension(:) :: LAMMPS_INSERT_NMODIFIED;

    integer (kind=4), allocatable, dimension(:) :: LAMMPS_INSERT_NTYPE_ATOM_MAX;

    integer (kind=4), allocatable, dimension(:) :: LAMMPS_INSERT_NTYPE_ATOM_MAX_OLD;

!   ************************************************************************************************

    real (kind=8), allocatable, dimension(:,:) :: LAMMPS_INSERT_COM_RG;

    real (kind=8), allocatable, dimension(:,:) :: LAMMPS_INSERT_TETHER_COEFFS;

    real (kind=8), allocatable, dimension(:,:,:) :: LAMMPS_INSERT_POTENTIAL_CLASS2;

!   ************************************************************************************************

    character (len=250), dimension(1:3) :: LAMMPS_INSERT_CHNAME, LAMMPS_INSERT_CHFILE_EXT;

!   ************************************************************************************************

    character (len=250), allocatable, dimension(:) :: LAMMPS_INSERT_NAME,          &
                                                      LAMMPS_INSERT_FILE_OPTION,   &
                                                      LAMMPS_INSERT_LIB_NAME,      &
                                                      LAMMPS_INSERT_CHFILE_FORMAT, &
                                                      LAMMPS_INSERT_CHNAME_FILE,   &
                                                      LAMMPS_INSERT_METHOD,        &
                                                      LAMMPS_INSERT_REGION_NAME;

    character (len=250), allocatable, dimension(:) :: LAMMPS_INSERT_NAME_OLD;

    character (len=250), allocatable, dimension(:) :: LAMMPS_INSERT_POTENTIAL_CLASS2_CHTYPE, &
                                                      LAMMPS_INSERT_TETHER_CHTYPE;

    character (len=250), allocatable, dimension(:,:) :: LAMMPS_INSERT_MODIFY_STYLE;

    character (len=250), allocatable, dimension(:,:,:) :: LAMMPS_INSERT_POTENTIAL_CLASS2_NAT;

end module module_inserts
