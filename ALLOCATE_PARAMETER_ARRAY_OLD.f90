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

subroutine ALLOCATE_PARAMETER_ARRAY_OLD(icanal,               &
                                        NTYPE_ATOM_TMP,       &
                                        NTYPE_BOND_TMP,       &
                                        NTYPE_ANGLE_TMP,      &
                                        NTYPE_DIHEDRAL_TMP,   &
                                        NTYPE_IMPROPER_TMP,   &
                                        NPAIR_COEFF_CROSS_TMP)

!   ************************************************************************************************
!   **                           ALLOCATE OLD PARAMETER ARRAYS                                     **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                    : CANAL ON WHICH OUTPUT DATA ARE WRITTEN                         **
!   **                                                                                            **
!   ************************************************************************************************

    use module_size_arrays;

    use module_library;

    use module_config;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    integer (kind=4), intent(in) :: NTYPE_ATOM_TMP,       &
                                    NTYPE_BOND_TMP,       &
                                    NTYPE_ANGLE_TMP,      &
                                    NTYPE_DIHEDRAL_TMP,   &
                                    NTYPE_IMPROPER_TMP,   &
                                    NPAIR_COEFF_CROSS_TMP;

!   ************************************************************************************************
                                                                                         !
    if ( NTYPE_ATOM_TMP > 0 ) then;                                                      !
                                                                                         !
        allocate(ATOM_MASSE_OLD(1:NTYPE_ATOM_TMP));                                      !
                                                                                         !
        allocate(ATOM_LABEL_OLD(1:NTYPE_ATOM_TMP));                                      !
                                                                                         !
        allocate(POTENTIAL_CLASS2_OLD(1:2,1:NTYPE_ATOM_TMP));                            !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_BOND_TMP > 0 )  then;                                                     !
                                                                                         !
        allocate(BOND_COEFFS_OLD(1:4,NTYPE_BOND_TMP));                                   !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_ANGLE_TMP > 0 ) then;                                                     !
                                                                                         !
        allocate(ANGLE_COEFFS_OLD(1:4,1:NTYPE_ANGLE_TMP));                               !
                                                                                         !
        allocate(BONDBOND_COEFFS_OLD(1:3,1:NTYPE_ANGLE_TMP));                            !
                                                                                         !
        allocate(BONDANGLE_COEFFS_OLD(1:4,1:NTYPE_ANGLE_TMP));                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_DIHEDRAL_TMP > 0 ) then;                                                  !
                                                                                         !
        allocate(DIHEDRAL_COEFFS_OLD(1:6,1:NTYPE_DIHEDRAL_TMP));                         !
                                                                                         !
        allocate(MIDDLEBONDTORSION_COEFFS_OLD(1:4,1:NTYPE_DIHEDRAL_TMP));                !
                                                                                         !
        allocate(ENDBONDTORSION_COEFFS_OLD(1:8,1:NTYPE_DIHEDRAL_TMP));                   !
                                                                                         !
        allocate(ANGLETORSION_COEFFS_OLD(1:8,1:NTYPE_DIHEDRAL_TMP));                     !
                                                                                         !
        allocate(ANGLEANGLETORSION_COEFFS_OLD(1:3,1:NTYPE_DIHEDRAL_TMP));                !
                                                                                         !
        allocate(BONDBOND13_COEFFS_OLD(1:3,1:NTYPE_DIHEDRAL_TMP));                       !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_IMPROPER_TMP > 0 ) then;                                                  !
                                                                                         !
        allocate(IMPROPER_COEFFS_OLD(1:2,1:NTYPE_IMPROPER_TMP));                         !
                                                                                         !
        allocate(ANGLEANGLE_COEFFS_OLD(1:6,1:NTYPE_IMPROPER_TMP));                       !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    if ( NPAIR_COEFF_CROSS_TMP > 0 ) then;                                               !
                                                                                         !
        allocate(PAIR_COEFF_CROSS_OLD(1:2,1:NPAIR_COEFF_CROSS_TMP));                     !
                                                                                         !
        allocate(PAIR_ATOMID_CROSS_OLD(1:2,1:NPAIR_COEFF_CROSS_TMP));                    !
                                                                                         !
    end if                                                                               !
                                                                                         !
end subroutine ALLOCATE_PARAMETER_ARRAY_OLD











