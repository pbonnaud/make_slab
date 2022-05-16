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


subroutine UPDATE_CONFIG_ARRAY(icanal,         &
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
        CONFIG_MOLECULEID(1:NATOM_TMP) = CONFIG_MOLECULEID_OLD(1:NATOM_TMP);             !
                                                                                         !
        CONFIG_ATOMID(1:NATOM_TMP)     = CONFIG_ATOMID_OLD(1:NATOM_TMP);                 !
                                                                                         !
        CONFIG_ATOM_TYPE(1:NATOM_TMP)  = CONFIG_ATOM_TYPE_OLD(1:NATOM_TMP);              !
                                                                                         !
        CONFIG_RI(1:3,1:NATOM_TMP)     = CONFIG_RI_OLD(1:3,1:NATOM_TMP);                 !
                                                                                         !   
        CONFIG_VI(1:3,1:NATOM_TMP)     = CONFIG_VI_OLD(1:3,1:NATOM_TMP);                 !
                                                                                         !
        CONFIG_QI(1:NATOM_TMP)         = CONFIG_QI_OLD(1:NATOM_TMP);                     !
                                                                                         !
        CONFIG_NAT(1:NATOM_TMP)        = CONFIG_NAT_OLD(1:NATOM_TMP);                    !
                                                                                         !
    end if;                                                                              !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    if ( NBOND_TMP > 0 ) then;                                                           !
                                                                                         !
        BOND_ATOMID(1:2,1:NBOND_TMP) = BOND_ATOMID_OLD(1:2,1:NBOND_TMP);                 !
                                                                                         !
        BOND_TYPE(1:NBOND_TMP)       = BOND_TYPE_OLD(1:NBOND_TMP);                       !
                                                                                         !
    end if;                                                                              !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    if ( NANGLE_TMP > 0 ) then;                                                          !
                                                                                         !
        ANGLE_ATOMID(1:3,1:NANGLE_TMP) = ANGLE_ATOMID_OLD(1:3,1:NANGLE_TMP);             !
                                                                                         !  
        ANGLE_TYPE(1:NANGLE_TMP)       = ANGLE_TYPE_OLD(1:NANGLE_TMP);                   !
                                                                                         !
    end if;                                                                              !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    if ( NDIHEDRAL_TMP > 0 ) then;                                                       !
                                                                                         !
        DIHEDRAL_ATOMID(1:4,1:NDIHEDRAL_TMP) = DIHEDRAL_ATOMID_OLD(1:4,1:NDIHEDRAL_TMP); !
                                                                                         !
        DIHEDRAL_TYPE(1:NDIHEDRAL_TMP)       = DIHEDRAL_TYPE_OLD(1:NDIHEDRAL_TMP);       !
                                                                                         !
    end if;                                                                              !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    if ( NIMPROPER_TMP > 0 ) then;                                                       !
                                                                                         !
        IMPROPER_ATOMID(1:4,1:NIMPROPER_TMP) = IMPROPER_ATOMID_OLD(1:4,1:NIMPROPER_TMP); !
                                                                                         !
        IMPROPER_TYPE(1:NIMPROPER_TMP)       = IMPROPER_TYPE_OLD(1:NIMPROPER_TMP);       !
                                                                                         !
    end if;                                                                              !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine UPDATE_CONFIG_ARRAY











