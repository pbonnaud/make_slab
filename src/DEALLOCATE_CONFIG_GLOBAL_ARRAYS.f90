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

subroutine DEALLOCATE_CONFIG_GLOBAL_ARRAYS(icanal,         &
                                           NATOM_TMP,      &
                                           NBOND_TMP,      &
                                           NANGLE_TMP,     &
                                           NDIHEDRAL_TMP,  &
                                           NIMPROPER_TMP)

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

    integer (kind=4), intent(in) :: NATOM_TMP,     &
                                    NBOND_TMP,     &
                                    NANGLE_TMP,    &
                                    NDIHEDRAL_TMP, &
                                    NIMPROPER_TMP;

!   ************************************************************************************************
                                                                                         !
    if ( NATOM_TMP > 0 ) then;                                                           !
                                                                                         !
        deallocate(CONFIG_MOLECULEID);                                                   !
                                                                                         !
        deallocate(CONFIG_ATOMID);                                                       !
                                                                                         !
        deallocate(CONFIG_RI);                                                           !
                                                                                         !
        deallocate(CONFIG_VI);                                                           !
                                                                                         !
        deallocate(CONFIG_ATOM_TYPE);                                                    !
                                                                                         !
        deallocate(CONFIG_QI);                                                           !
                                                                                         !
        deallocate(CONFIG_NAT);                                                          !
                                                                                         !
    end if;                                                                              !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    if ( NBOND_TMP > 0 ) then;                                                           !
                                                                                         !
        deallocate(BOND_ATOMID);                                                         !
                                                                                         !
        deallocate(BOND_TYPE);                                                           !
                                                                                         !
    end if;                                                                              !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    if ( NANGLE_TMP > 0 ) then;                                                          !
                                                                                         !
        deallocate(ANGLE_ATOMID);                                                        !
                                                                                         !
        deallocate(ANGLE_TYPE);                                                          !
                                                                                         !
    end if;                                                                              !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    if ( NDIHEDRAL_TMP > 0 ) then;                                                       !
                                                                                         !
        deallocate(DIHEDRAL_ATOMID);                                                     !
                                                                                         !
        deallocate(DIHEDRAL_TYPE);                                                       !
                                                                                         !
    end if;                                                                              !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    if ( NIMPROPER_TMP > 0 ) then;                                                       !
                                                                                         !
        deallocate(IMPROPER_ATOMID);                                                     !
                                                                                         !
        deallocate(IMPROPER_TYPE);                                                       !
                                                                                         !
    end if;                                                                              !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine DEALLOCATE_CONFIG_GLOBAL_ARRAYS











