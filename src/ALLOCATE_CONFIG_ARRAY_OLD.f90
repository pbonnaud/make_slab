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



subroutine ALLOCATE_CONFIG_ARRAY_OLD(icanal,         &
                                     NATOM_TMP,      &
                                     NBOND_TMP,      &
                                     NANGLE_TMP,     &
                                     NDIHEDRAL_TMP,  &
                                     NIMPROPER_TMP)


!   ************************************************************************************************
!   **                             ALLOCATE OLD CONFIG ARRAYS                                     **
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
        allocate(CONFIG_MOLECULEID_OLD(1:NATOM_TMP));                                    !
                                                                                         !
        allocate(CONFIG_ATOMID_OLD(1:NATOM_TMP));                                        !
                                                                                         !
        allocate(CONFIG_RI_OLD(1:3,1:NATOM_TMP));                                        !
                                                                                         !
        allocate(CONFIG_VI_OLD(1:3,1:NATOM_TMP));                                        !
                                                                                         !
        allocate(CONFIG_ATOM_TYPE_OLD(1:NATOM_TMP));                                     !
                                                                                         !
        allocate(CONFIG_QI_OLD(1:NATOM_TMP));                                            !
                                                                                         !
        allocate(CONFIG_NAT_OLD(1:NATOM_TMP));                                           !
                                                                                         !
    end if;                                                                              !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    if ( NBOND_TMP > 0 ) then;                                                           !
                                                                                         ! 
        allocate(BOND_ATOMID_OLD(1:2,1:NBOND_TMP));                                      !
                                                                                         !
        allocate(BOND_TYPE_OLD(1:NBOND_TMP));                                            !
                                                                                         !
    end if;                                                                              !
                                                                                         !
    if ( NANGLE_TMP > 0 ) then;                                                          !
                                                                                         !
        allocate(ANGLE_ATOMID_OLD(1:3,1:NANGLE_TMP));                                    !
                                                                                         !
        allocate(ANGLE_TYPE_OLD(1:NANGLE_TMP));                                          !
                                                                                         !
    end if;                                                                              !
                                                                                         !
    if ( NDIHEDRAL_TMP > 0 ) then;                                                       !
                                                                                         !
        allocate(DIHEDRAL_ATOMID_OLD(1:4,1:NDIHEDRAL_TMP));                              !
                                                                                         !
        allocate(DIHEDRAL_TYPE_OLD(1:NDIHEDRAL_TMP));                                    !
                                                                                         !
    end if;                                                                              !
                                                                                         !
    if ( NIMPROPER_TMP > 0 ) then;                                                       !
                                                                                         !
        allocate(IMPROPER_ATOMID_OLD(1:4,1:NIMPROPER_TMP));                              !
                                                                                         !
        allocate(IMPROPER_TYPE_OLD(1:NIMPROPER_TMP));                                    !
                                                                                         !
    end if;                                                                              !
                                                                                         !
end subroutine ALLOCATE_CONFIG_ARRAY_OLD











