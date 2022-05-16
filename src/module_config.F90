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



module module_config

    use module_size_arrays;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4) :: IPOTENTIAL_CLASS2, ATOMS_FLAG;

    integer (kind=4) :: NATOM,            &
                        NBOND,            &
                        NANGLE,           &
                        NDIHEDRAL,        &
                        NIMPROPER,        &
                        NTYPE_ATOM,       &
                        NTYPE_BOND,       &
                        NTYPE_ANGLE,      &
                        NTYPE_DIHEDRAL,   &
                        NTYPE_IMPROPER,   &
                        NPARAM_BONDS,     &
                        NPARAM_ANGLES,    &
                        NPARAM_DIHEDRALS, &
                        NPAIR_COEFF_CROSS;

    integer (kind=4), allocatable, dimension(:) :: CONFIG_ATOM_TYPE,      &
                                                   BOND_TYPE,             &
                                                   ANGLE_TYPE,            &
                                                   IMPROPER_TYPE,         &
                                                   DIHEDRAL_TYPE,         &
                                                   CONFIG_ATOMID,         &
                                                   CONFIG_MOLECULEID;

    integer (kind=4), allocatable, dimension(:) :: CONFIG_ATOM_TYPE_OLD,  &
                                                   BOND_TYPE_OLD,         &
                                                   ANGLE_TYPE_OLD,        &
                                                   IMPROPER_TYPE_OLD,     &
                                                   DIHEDRAL_TYPE_OLD,     &
                                                   CONFIG_ATOMID_OLD,     &
                                                   CONFIG_MOLECULEID_OLD;

    integer (kind=4), allocatable, dimension(:,:) :: BOND_ATOMID,     &
                                                     ANGLE_ATOMID,    &
                                                     IMPROPER_ATOMID, &
                                                     DIHEDRAL_ATOMID;

    integer (kind=4), allocatable, dimension(:,:) :: BOND_ATOMID_OLD,       &
                                                     ANGLE_ATOMID_OLD,      &
                                                     IMPROPER_ATOMID_OLD,   &
                                                     DIHEDRAL_ATOMID_OLD,   &
                                                     PAIR_ATOMID_CROSS,     &
                                                     PAIR_ATOMID_CROSS_OLD;

!   ************************************************************************************************

    real (kind=8), allocatable, dimension(:) :: ATOM_MASSE, ATOM_MASSE_OLD;

    real (kind=8), allocatable, dimension(:) :: CONFIG_QI, CONFIG_QI_OLD;

    real (kind=8), allocatable, dimension(:,:) :: POTENTIAL_CLASS2,      &
                                                  POTENTIAL_CLASS2_OLD,  &
                                                  PAIR_COEFF_CROSS,      &
                                                  PAIR_COEFF_CROSS_OLD;

    real (kind=8), allocatable, dimension(:,:) :: CONFIG_RI,     &
                                                  CONFIG_VI,     &
                                                  CONFIG_RI_OLD, &
                                                  CONFIG_VI_OLD;

    real (kind=8), allocatable, dimension(:,:) :: IMPROPER_COEFFS;

    real (kind=8), allocatable, dimension(:,:) :: BONDBOND_COEFFS,           &
                                                  ANGLEANGLETORSION_COEFFS,  &
                                                  BONDBOND13_COEFFS;

    real (kind=8), allocatable, dimension(:,:) :: BOND_COEFFS,               &
                                                  ANGLE_COEFFS,              &
                                                  BONDANGLE_COEFFS,          &
                                                  MIDDLEBONDTORSION_COEFFS;

    real (kind=8), allocatable, dimension(:,:) :: ANGLEANGLE_COEFFS,         &
                                                  DIHEDRAL_COEFFS;

    real (kind=8), allocatable, dimension(:,:) :: ENDBONDTORSION_COEFFS,     &
                                                  ANGLETORSION_COEFFS;


    real (kind=8), allocatable, dimension(:,:) :: BOND_COEFFS_OLD,              &
                                                  ANGLE_COEFFS_OLD,             &
                                                  BONDBOND_COEFFS_OLD,          &
                                                  BONDANGLE_COEFFS_OLD,         &
                                                  DIHEDRAL_COEFFS_OLD,          &
                                                  MIDDLEBONDTORSION_COEFFS_OLD, &
                                                  ENDBONDTORSION_COEFFS_OLD,    &
                                                  ANGLETORSION_COEFFS_OLD,      &
                                                  ANGLEANGLETORSION_COEFFS_OLD, &
                                                  BONDBOND13_COEFFS_OLD,        &
                                                  IMPROPER_COEFFS_OLD,          &
                                                  ANGLEANGLE_COEFFS_OLD;

!   ************************************************************************************************

    character (len=250) :: CHWORKING_FILE_REMOVE, &
                           CHWORKING_FILE_INSERTION;

    character (len=250) :: CH_UNITS_STYLE,          &
                           CH_ATOM_STYLE,           &
                           POTENTIAL_CLASS2_CHTYPE, &
                           CH_SPECIAL_BONDS,        &
                           CH_PAIR_MODIFY,          &
                           CH_KSPACE_STYLE,         &
                           CH_BOND_STYLE,           &
                           CH_ANGLE_STYLE,          &
                           CH_DIHEDRAL_STYLE,       &
                           CH_PAIR_STYLE;

    character (len=20), allocatable, dimension(:) :: ATOM_LABEL, ATOM_LABEL_OLD;

    character (len=20), allocatable, dimension(:) :: CONFIG_NAT, CONFIG_NAT_OLD;

end module module_config
