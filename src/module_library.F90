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


module module_library

    use module_size_arrays;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4) :: NPATH_TO_LIBRARY;

    integer (kind=4) :: NFILE_LIBRARY,             &
                        NFILE_LIBRARY_REMOVE,      &
                        NTYPE_DUMMY_PARTICLE,      &
                        NGenerate_polymer_species, &
                        NTotal_monomers;

    integer (kind=4) :: IMAX_NATOM,     &
                        IMAX_NBOND,     &
                        IMAX_NANGLE,    &
                        IMAX_NDIHEDRAL, & 
                        IMAX_NIMPROPER;

    integer (kind=4) :: IMAX_NTYPE_NATOM,      &
                        IMAX_NTYPE_BOND,       &
                        IMAX_NTYPE_ANGLE,      &
                        IMAX_NTYPE_DIHEDRAL,   &
                        IMAX_NTYPE_IMPROPER,   &
                        IMAX_NPAIR_COEFF_CROSS;

    integer (kind=4), allocatable, dimension(:) :: NMOLECULE_INSERTION,         &
                                                   NMOLECULE_DELETION,          &
                                                   NDUMMY_PARTICLE,             &
                                                   NPOLYMER_CHAINS,             &
                                                   NSPECIES_PER_POLYMER_CHAINS, &
                                                   IPOTENTIAL_CLASS2_LIBRARY,   &
                                                   NPAIR_COEFF_CROSS_LIBRARY,   &
                                                   NATOM_LIBRARY,               &
                                                   NBOND_LIBRARY,               &
                                                   NANGLE_LIBRARY,              &
                                                   NDIHEDRAL_LIBRARY,           &
                                                   NIMPROPER_LIBRARY,           &
                                                   NTYPE_ATOM_LIBRARY,          &
                                                   NTYPE_BOND_LIBRARY,          &
                                                   NTYPE_ANGLE_LIBRARY,         &
                                                   NTYPE_DIHEDRAL_LIBRARY,      &
                                                   NTYPE_IMPROPER_LIBRARY,      &
                                                   NPARAM_BONDS_LIBRARY,        &
                                                   NPARAM_ANGLES_LIBRARY,       &
                                                   NPARAM_DIHEDRALS_LIBRARY,    &
                                                   ILIBRARY_FILE_INFO,          &
                                                   INSERTION_METHOD_LIBRARY,    &
                                                   ILOCR1_LIBRARY,              &
                                                   ILOCR2_LIBRARY;

    integer (kind=4), allocatable, dimension(:) :: LIBRARY_RANDOM_SEED;

    integer (kind=4), allocatable, dimension(:) :: LIBRARY_REGION_IRANK;

    integer (kind=4), allocatable, dimension(:,:) :: NMONOMER_UNITS_PER_CHAINS,   &
                                                     CONFIG_ATOM_TYPE_LIBRARY,    &
                                                     BOND_TYPE_LIBRARY,           & 
                                                     ANGLE_TYPE_LIBRARY,          &
                                                     DIHEDRAL_TYPE_LIBRARY,       &
                                                     IMPROPER_TYPE_LIBRARY;

    integer (kind=4), allocatable, dimension(:,:) :: CONFIG_ATOMID_LIBRARY;

    integer (kind=4), allocatable, dimension(:,:,:) :: BOND_ATOMID_LIBRARY,      &
                                                       ANGLE_ATOMID_LIBRARY,     &
                                                       DIHEDRAL_ATOMID_LIBRARY,  &
                                                       IMPROPER_ATOMID_LIBRARY,  &
                                                       PAIR_ATOMID_CROSS_LIBRARY;

!   ************************************************************************************************

    real (kind=8), allocatable, dimension(:,:) :: MONOMER_INSERTION_PARAM;

!   ************************************************************************************************

    real (kind=8), allocatable, dimension(:) :: RADIUS_DUMMY_PARTICLE, &
                                                MASSE_DUMMY_PARTICLE;

    real (kind=8), allocatable, dimension(:,:) :: ATOM_MASSES_LIBRARY, &
                                                  CELL_AXIS_LIBRARY;

    real (kind=8), allocatable, dimension(:,:,:) :: POTENTIAL_CLASS2_LIBRARY, &
                                                    PAIR_COEFF_CROSS_LIBRARY, &
                                                    MOLECULE_COM_LIBRARY;

    real (kind=8), allocatable, dimension(:,:) :: CONFIG_QI_LIBRARY;

    real (kind=8), allocatable, dimension(:,:,:) :: IMPROPER_COEFFS_LIBRARY;

    real (kind=8), allocatable, dimension(:,:,:) :: BONDBOND_COEFFS_LIBRARY,           &
                                                    ANGLEANGLETORSION_COEFFS_LIBRARY,  &
                                                    BONDBOND13_COEFFS_LIBRARY;

    real (kind=8), allocatable, dimension(:,:,:) :: CONFIG_RI_LIBRARY, &
                                                    CONFIG_VI_LIBRARY;

    real (kind=8), allocatable, dimension(:,:,:) :: BOND_COEFFS_LIBRARY,               &
                                                    ANGLE_COEFFS_LIBRARY,              &
                                                    BONDANGLE_COEFFS_LIBRARY,          &
                                                    MIDDLEBONDTORSION_COEFFS_LIBRARY;

    real (kind=8), allocatable, dimension(:,:,:) :: ANGLEANGLE_COEFFS_LIBRARY,         &
                                                    DIHEDRAL_COEFFS_LIBRARY;

    real (kind=8), allocatable, dimension(:,:,:) :: ENDBONDTORSION_COEFFS_LIBRARY,     &
                                                    ANGLETORSION_COEFFS_LIBRARY;

!   ************************************************************************************************

    character (len=50) :: CHFILE_FORMAT,        &
                          CHFILE_FORMAT_REMOVE, &
                          Chfile_format_polymer;

    character (len=150), allocatable, dimension(:) :: CHNAME_FILE_LIBRARY_REMOVE, &
                                                      CHNAME_DUMMY_PARTICLE,      &
                                                      Chname_file_library_monomer;

    character (len=250),allocatable, dimension(:) :: CH_UNITS_STYLE_LIBRARY,          &
                                                     CH_SPECIAL_BONDS_LIBRARY,        &
                                                     CH_PAIR_MODIFY_LIBRARY,          &
                                                     CH_KSPACE_STYLE_LIBRARY,         &
                                                     CH_ATOM_STYLE_LIBRARY,           &
                                                     CH_BOND_STYLE_LIBRARY,           &
                                                     CH_ANGLE_STYLE_LIBRARY,          &
                                                     CH_DIHEDRAL_STYLE_LIBRARY,       &
                                                     POTENTIAL_CLASS2_CHTYPE_LIBRARY;

    character (len=250), allocatable, dimension(:) :: LIBRARY_REGION_NAME_INSERTION;

    character (len=250), allocatable, dimension(:) :: LIB_PATH_NAME, PATH_TO_LIBRARY;

    character (len=150), allocatable, dimension(:,:) :: CHNAME_FILE_LIBRARY;

    character (len=20), allocatable, dimension(:,:) :: ATOM_LABEL_LIBRARY;

    character (len=20), allocatable, dimension(:,:) :: CONFIG_NAT_LIBRARY;

end module module_library
