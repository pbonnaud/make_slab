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



subroutine UPDATE_SIZE_OF_ARRAYS(icanal,                   &
                                   NATOM_TMP,              &
                                   NBOND_TMP,              &
                                   NANGLE_TMP,             &
                                   NDIHEDRAL_TMP,          &
                                   NIMPROPER_TMP,          &
                                   NTYPE_ATOM_TMP,         &
                                   NTYPE_BOND_TMP,         &
                                   NTYPE_ANGLE_TMP,        &
                                   NTYPE_DIHEDRAL_TMP,     &
                                   NTYPE_IMPROPER_TMP,     &
                                   NPAIR_COEFF_CROSS_TMP);



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

    integer (kind=4), intent(inout) :: NATOM_TMP,            &
                                       NBOND_TMP,            &
                                       NANGLE_TMP,           &
                                       NDIHEDRAL_TMP,        &
                                       NIMPROPER_TMP,        &
                                       NTYPE_ATOM_TMP,       &
                                       NTYPE_BOND_TMP,       &
                                       NTYPE_ANGLE_TMP,      &
                                       NTYPE_DIHEDRAL_TMP,   &
                                       NTYPE_IMPROPER_TMP,   &
                                       NPAIR_COEFF_CROSS_TMP;

!   ************************************************************************************************

    integer (kind=4) :: i;

!   ************************************************************************************************
                                                                                         !
    do i = 1, NFILE_LIBRARY;                                                             ! Loop over the molecule read in the library
                                                                                         !
        NATOM_TMP = NATOM_TMP + NATOM_LIBRARY(i) * NMOLECULE_INSERTION(i);               !
                                                                                         !
        if ( NBOND_LIBRARY(i) > 0 ) then;                                                !
                                                                                         !
            NBOND_TMP = NBOND_TMP + NBOND_LIBRARY(i) * NMOLECULE_INSERTION(i);           !
                                                                                         !
        end if                                                                           !
                                                                                         !
        if ( NANGLE_LIBRARY(i) > 0 ) then;                                               !
                                                                                         !
            NANGLE_TMP = NANGLE_TMP + NANGLE_LIBRARY(i) * NMOLECULE_INSERTION(i);        !
                                                                                         !
        end if                                                                           !
                                                                                         !
        if ( NDIHEDRAL_LIBRARY(i) > 0 ) then;                                            !
                                                                                         !
            NDIHEDRAL_TMP          = &                                                   !
            NDIHEDRAL_TMP          + &                                                   !
            NDIHEDRAL_LIBRARY(i)   * &                                                   !
            NMOLECULE_INSERTION(i);                                                      !
                                                                                         !
        end if                                                                           !
                                                                                         ! 
        if ( NIMPROPER_LIBRARY(i) > 0 ) then;                                            ! 
                                                                                         !
            NIMPROPER_TMP          = &                                                   !
            NIMPROPER_TMP          + &                                                   !
            NIMPROPER_LIBRARY(i)   * &                                                   ! 
            NMOLECULE_INSERTION(i);                                                      !
                                                                                         !
        end if                                                                           !
                                                                                         !
        NTYPE_ATOM_TMP = NTYPE_ATOM_TMP + NTYPE_ATOM_LIBRARY(i);                         !
                                                                                         !
        if ( NTYPE_BOND_LIBRARY(i) > 0 ) then;                                           !
                                                                                         !
            NTYPE_BOND_TMP = NTYPE_BOND_TMP + NTYPE_BOND_LIBRARY(i);                     !
                                                                                         !
        end if                                                                           !
                                                                                         !
        if ( NTYPE_ANGLE_LIBRARY(i) > 0 ) then;                                          !
                                                                                         !
            NTYPE_ANGLE_TMP       = &                                                    !
            NTYPE_ANGLE_TMP       + &                                                    !
            NTYPE_ANGLE_LIBRARY(i);                                                      !
                                                                                         !
        end if                                                                           !
                                                                                         !
        if ( NTYPE_DIHEDRAL_LIBRARY(i) > 0 ) then;                                       !
                                                                                         !
            NTYPE_DIHEDRAL_TMP       = &                                                 !
            NTYPE_DIHEDRAL_TMP       + &                                                 !
            NTYPE_DIHEDRAL_LIBRARY(i);                                                   !
                                                                                         !
        end if                                                                           !
                                                                                         !
        if ( NTYPE_IMPROPER_LIBRARY(i) > 0 ) then;                                       !
                                                                                         !
            NTYPE_IMPROPER_TMP       = &                                                 !
            NTYPE_IMPROPER_TMP       + &                                                 !
            NTYPE_IMPROPER_LIBRARY(i);                                                   !
                                                                                         !
        end if                                                                           !
                                                                                         !
        if ( NPAIR_COEFF_CROSS_LIBRARY(i) > 0 ) then;                                    !
                                                                                         !
            NPAIR_COEFF_CROSS_TMP       = &                                              !
            NPAIR_COEFF_CROSS_TMP       + &                                              !
            NPAIR_COEFF_CROSS_LIBRARY(i);                                                !
                                                                                         !
        end if                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   ### Write the new size of configuration arrays #################################################
                                                                                         !
    write(icanal,'(a23,i8,a39)') '| NATOM(UPDATE)      : ', &                            !
                                 NATOM_TMP,                 &                            !
                                 REPEAT(' ',38)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    if ( NBOND_TMP > 0 ) then;                                                           !
                                                                                         !
        write(icanal,'(a23,i8,a39)') '| NBOND(UPDATE)      : ', &                        !
                                     NBOND_TMP,                 &                        !
                                     REPEAT(' ',38)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NANGLE_TMP > 0 ) then;                                                          !
                                                                                         !
        write(icanal,'(a23,i8,a39)') '| NANGLE(UPDATE)     : ', &                        !
                                     NANGLE_TMP,                &                        !
                                     REPEAT(' ',38)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NDIHEDRAL_TMP > 0 ) then;                                                       !
                                                                                         !
        write(icanal,'(a23,i8,a39)') '| NDIHEDRAL(UPDATE)  : ', &                        !
                                     NDIHEDRAL_TMP,             &                        !
                                     REPEAT(' ',38)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NIMPROPER_TMP > 0 ) then;                                                       !
                                                                                         !
        write(icanal,'(a23,i8,a39)') '| NIMPROPER(UPDATE)  : ', &                        !
                                     NIMPROPER_TMP,             &                        !
                                     REPEAT(' ',38)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_ATOM_TMP > 0 ) then;                                                      !
                                                                                         !
       write(icanal,'(a23,i8,a39)') '| NTYPE_ATOM(UPDATE) : ', &                         !
                                    NTYPE_ATOM_TMP,            &                         !
                                    REPEAT(' ',38)//'|';                                 !
                                                                                         !
       write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                   !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_BOND_TMP > 0 ) then;                                                      !
                                                                                         !
        write(icanal,'(a23,i8,a39)') '| NTYPE_BOND(UPDATE) : ', &                        !
                                     NTYPE_BOND_TMP,            &                        !
                                     REPEAT(' ',38)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_ANGLE_TMP > 0 ) then;                                                     !
                                                                                         !
       write(icanal,'(a24,i8,a38)') '| NTYPE_ANGLE(UPDATE) : ', &                        !
                                    NTYPE_ANGLE_TMP,            &                        !
                                    REPEAT(' ',37)//'|';                                 !
                                                                                         !
       write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                   !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_DIHEDRAL_TMP > 0 ) then;                                                  !
                                                                                         !
        write(icanal,'(a27,i8,a35)') '| NTYPE_DIHEDRAL(UPDATE) : ', &                    !
                                     NTYPE_DIHEDRAL_TMP,            &                    !
                                     REPEAT(' ',34)//'|';                                !
                                                                                         ! 
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_IMPROPER_TMP > 0 ) then;                                                  !
                                                                                         !
        write(icanal,'(a27,i8,a35)') '| NTYPE_IMPROPER(UPDATE) : ', &                    !
                                     NTYPE_IMPROPER_TMP,            &                    !
                                     REPEAT(' ',34)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NPAIR_COEFF_CROSS_TMP > 0 ) then;                                               !
                                                                                         !
        write(icanal,'(a30,i8,a32)') '| NPAIR_COEFF_CROSS(UPDATE) : ', &                 !
                                     NPAIR_COEFF_CROSS_TMP,            &                 !
                                     REPEAT(' ',31)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine UPDATE_SIZE_OF_ARRAYS











