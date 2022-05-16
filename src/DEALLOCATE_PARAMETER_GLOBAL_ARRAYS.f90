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

subroutine DEALLOCATE_PARAMETER_GLOBAL_ARRAYS(icanal,                &
                                              NTYPE_ATOM_TMP,        &
                                              NTYPE_BOND_TMP,        &
                                              NTYPE_ANGLE_TMP,       &
                                              NTYPE_DIHEDRAL_TMP,    &
                                              NTYPE_IMPROPER_TMP,    &
                                              NPAIR_COEFF_CROSS_TMP);

!   ************************************************************************************************
!   **                           DEALLOCATE CONFIG GLOBAL ARRAYS                                  **
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
        deallocate(ATOM_MASSE);                                                          !
                                                                                         !
        deallocate(ATOM_LABEL);                                                          ! 
                                                                                         !
        deallocate(POTENTIAL_CLASS2);                                                    !
                                                                                         !
    end if                                                                               !
                                                                                         ! 
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    if ( NTYPE_BOND_TMP > 0 ) then;                                                      !
                                                                                         !
        deallocate(BOND_COEFFS);                                                         !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    if ( NTYPE_ANGLE_TMP > 0 ) then;                                                     !
                                                                                         !
        deallocate(ANGLE_COEFFS);                                                        !
                                                                                         !
        deallocate(BONDBOND_COEFFS);                                                     !
                                                                                         !
        deallocate(BONDANGLE_COEFFS);                                                    !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    if ( NTYPE_DIHEDRAL_TMP > 0 ) then;                                                  !
                                                                                         !
        deallocate(DIHEDRAL_COEFFS);                                                     !
                                                                                         !
        deallocate(MIDDLEBONDTORSION_COEFFS);                                            !
                                                                                         !
        deallocate(ENDBONDTORSION_COEFFS);                                               !
                                                                                         !
        deallocate(ANGLETORSION_COEFFS);                                                 !
                                                                                         !
        deallocate(ANGLEANGLETORSION_COEFFS);                                            !
                                                                                         !  
        deallocate(BONDBOND13_COEFFS);                                                   !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    if ( NTYPE_IMPROPER_TMP > 0 ) then;                                                  !
                                                                                         !
        deallocate(IMPROPER_COEFFS);                                                     !
                                                                                         !
        deallocate(ANGLEANGLE_COEFFS);                                                   !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    if ( NPAIR_COEFF_CROSS_TMP > 0 ) then;                                               !
                                                                                         !
        deallocate(PAIR_COEFF_CROSS);                                                    !
                                                                                         !
        deallocate(PAIR_ATOMID_CROSS);                                                   !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!

end subroutine DEALLOCATE_PARAMETER_GLOBAL_ARRAYS











