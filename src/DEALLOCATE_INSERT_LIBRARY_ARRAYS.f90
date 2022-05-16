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

subroutine DEALLOCATE_INSERT_LIBRARY_ARRAYS(icanal,ILOCAL_INSERTS)

!   ************************************************************************************************
!   **                      Allocate library arrays for the insert command                        **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                           **
!   **                                                                                            **
!   ************************************************************************************************

    use module_size_arrays;

    use module_library;

    use module_osef;

    use module_inserts;

!   ************************************************************************************************

    implicit none;

!   ***********************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    integer (kind=4), intent(in) :: ILOCAL_INSERTS;

!   ************************************************************************************************

    integer (kind=4) :: i, j, k;

    integer (kind=4) :: IFOUND;

    integer (kind=4) :: ilocal_file

    character (len=20), dimension(1:100) :: LOCAL_ATOM_LABEL;

!   ************************************************************************************************

    integer (kind=4) :: NFILE_LIBRARY_TMP;    

!   ************************************************************************************************

    integer (kind=4) :: IOSEF_NATOM,     &
                        IOSEF_NBOND,     &
                        IOSEF_NANGLE,    &
                        IOSEF_NDIHEDRAL, &
                        IOSEF_NIMPROPER;

    integer (kind=4) :: IOSEF_NTYPE_NATOM,       &
                        IOSEF_NTYPE_BOND,        &
                        IOSEF_NTYPE_ANGLE,       &
                        IOSEF_NTYPE_DIHEDRAL,    &
                        IOSEF_NTYPE_IMPROPER,    &
                        IOSEF_NPAIR_COEFF_CROSS;

    character (len=250) :: CHEXT;

!   ************************************************************************************************

    character (len=250) :: CHARLINE;

    integer (kind=4) :: EOF;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the routine #############################################################
                                                                                         !
    CHTITLE = 'Dellocate library arrays for the insert command';                         !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate arrays for insertion files ########################################################
                                                                                         !
    if ( NFILE_LIBRARY > 0 ) then;                                                       !
                                                                                         !
        deallocate(NMOLECULE_INSERTION);                                                 !
                                                                                         !
        deallocate(CHNAME_FILE_LIBRARY);                                                 !
                                                                                         !
        deallocate(INSERTION_METHOD_LIBRARY);                                            !
                                                                                         !
        deallocate(ILIBRARY_FILE_INFO);                                                  !
                                                                                         !
        write(icanal,'(a70)') '| Arrays for insertion files were deallocated'// &        !
                              REPEAT(' ',24)//'|';                                       ! 
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set properties depending on the insertion method ###########################################
                                                                                         !
    select case(TRIM(LAMMPS_INSERT_METHOD(ILOCAL_INSERTS)));                             !
                                                                                         !
        case('com');                                                                     !
                                                                                         !
            deallocate(MOLECULE_COM_LIBRARY);                                            !
                                                                                         !
        case('rand/region');                                                             !
                                                                                         !
            deallocate(LIBRARY_REGION_NAME_INSERTION);                                   !
                                                                                         !
            deallocate(LIBRARY_RANDOM_SEED);                                             !
                                                                                         !
            deallocate(LIBRARY_REGION_IRANK);                                            !
                                                                                         !
        case default;                                                                    !
                                                                                         !
            write(icanal,'(a70)') '| Not implemented - stop '// &                        !
                                  REPEAT(' ',44)//'|';                                   !
                                                                                         !
            call CLOSING_PROGRAM(icanal,1);                                              !
                                                                                         !
    end select                                                                           !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate library arrays for cell dimensions ################################################
                                                                                         !
    deallocate(CELL_AXIS_LIBRARY);                                                       !
                                                                                         !
!   ### Allocate library arrays ####################################################################
                                                                                         !
    deallocate(NATOM_LIBRARY);                                                           !
                                                                                         !
    deallocate(NBOND_LIBRARY);                                                           !
                                                                                         !
    deallocate(NANGLE_LIBRARY);                                                          !
                                                                                         !
    deallocate(NDIHEDRAL_LIBRARY);                                                       !
                                                                                         !
    deallocate(NIMPROPER_LIBRARY);                                                       !
                                                                                         !
    deallocate(NTYPE_ATOM_LIBRARY);                                                      !
                                                                                         !
    deallocate(NTYPE_BOND_LIBRARY);                                                      !
                                                                                         !
    deallocate(NTYPE_ANGLE_LIBRARY);                                                     !
                                                                                         !
    deallocate(NTYPE_DIHEDRAL_LIBRARY);                                                  !
                                                                                         ! 
    deallocate(NTYPE_IMPROPER_LIBRARY);                                                  !
                                                                                         !
    deallocate(CH_UNITS_STYLE_LIBRARY);                                                  !
                                                                                         !
    deallocate(CH_SPECIAL_BONDS_LIBRARY);                                                !
                                                                                         !
    deallocate(CH_PAIR_MODIFY_LIBRARY);                                                  !
                                                                                         !
    deallocate(CH_KSPACE_STYLE_LIBRARY);                                                 !
                                                                                         !
    deallocate(IPOTENTIAL_CLASS2_LIBRARY);                                               !
                                                                                         !
    deallocate(NPAIR_COEFF_CROSS_LIBRARY);                                               ! 
                                                                                         !
    deallocate(CH_ATOM_STYLE_LIBRARY);                                                   !
                                                                                         !
    deallocate(CH_BOND_STYLE_LIBRARY);                                                   ! 
                                                                                         !
    deallocate(CH_ANGLE_STYLE_LIBRARY);                                                  !
                                                                                         !
    deallocate(CH_DIHEDRAL_STYLE_LIBRARY);                                               !
                                                                                         !
    deallocate(POTENTIAL_CLASS2_CHTYPE_LIBRARY);                                         !
                                                                                         !
    deallocate(NPARAM_BONDS_LIBRARY);                                                    !
                                                                                         !
    deallocate(NPARAM_ANGLES_LIBRARY);                                                   !
                                                                                         !
    deallocate(NPARAM_DIHEDRALS_LIBRARY);                                                !
                                                                                         !
    if ( NGenerate_polymer_species > 0 ) then;                                           !
                                                                                         !
        deallocate(ILOCR1_LIBRARY);                                                      !
                                                                                         !
        deallocate(ILOCR2_LIBRARY);                                                      !
                                                                                         !
    end if                                                                               !
                                                                                         !
    write(icanal,'(a70)') '| Library arrays were deallocated'//REPEAT(' ',36)//'|';      !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Deallocation of library arrays #############################################################
                                                                                         !
    if ( ( IMAX_NATOM > 0 ) .AND. ( NFILE_LIBRARY > 0 ) ) then;                          !
                                                                                         !
        deallocate(CONFIG_ATOMID_LIBRARY);                                               !
                                                                                         !
        deallocate(CONFIG_ATOM_TYPE_LIBRARY);                                            !
                                                                                         !
        deallocate(CONFIG_QI_LIBRARY);                                                   !
                                                                                         !
        deallocate(CONFIG_VI_LIBRARY);                                                   !
                                                                                         !
        deallocate(CONFIG_RI_LIBRARY);                                                   !
                                                                                         !
        deallocate(CONFIG_NAT_LIBRARY);                                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( ( IMAX_NBOND > 0 ) .AND. ( NFILE_LIBRARY > 0 ) ) then;                          !
                                                                                         !
        deallocate(BOND_TYPE_LIBRARY);                                                   !
                                                                                         !
        deallocate(BOND_ATOMID_LIBRARY);                                                 !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( ( IMAX_NANGLE > 0 ) .AND. ( NFILE_LIBRARY > 0 ) ) then;                         !
                                                                                         !
        deallocate(ANGLE_TYPE_LIBRARY);                                                  !
                                                                                         !
        deallocate(ANGLE_ATOMID_LIBRARY);                                                !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( ( IMAX_NDIHEDRAL > 0 ) .AND. ( NFILE_LIBRARY > 0 ) ) then;                      !
                                                                                         !
        deallocate(DIHEDRAL_TYPE_LIBRARY);                                               !
                                                                                         !
        deallocate(DIHEDRAL_ATOMID_LIBRARY);                                             !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( ( IMAX_NIMPROPER > 0 ) .AND. ( NFILE_LIBRARY > 0 ) ) then;                      ! 
                                                                                         !
        deallocate(IMPROPER_TYPE_LIBRARY);                                               !
                                                                                         !
        deallocate(IMPROPER_ATOMID_LIBRARY);                                             !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( ( IMAX_NTYPE_NATOM > 0 ) .AND. ( NFILE_LIBRARY > 0 ) ) then;                    !
                                                                                         !
        deallocate(ATOM_MASSES_LIBRARY);                                                 !
                                                                                         !
        deallocate(ATOM_LABEL_LIBRARY);                                                  !
                                                                                         !
        deallocate(POTENTIAL_CLASS2_LIBRARY);                                            !
                                                                                         !
        deallocate(PAIR_COEFF_CROSS_LIBRARY);                                            !
                                                                                         !
        deallocate(PAIR_ATOMID_CROSS_LIBRARY);                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( ( IMAX_NTYPE_BOND > 0 ) .AND. ( NFILE_LIBRARY > 0 ) ) then;                     !
                                                                                         !
        deallocate(BOND_COEFFS_LIBRARY);                                                 !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( ( IMAX_NTYPE_ANGLE > 0 ) .AND. ( NFILE_LIBRARY > 0 ) ) then;                    !
                                                                                         !
        deallocate(ANGLE_COEFFS_LIBRARY);                                                !
                                                                                         !
        deallocate(BONDBOND_COEFFS_LIBRARY);                                             !
                                                                                         !
        deallocate(BONDANGLE_COEFFS_LIBRARY);                                            !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( ( IMAX_NTYPE_DIHEDRAL > 0 ) .AND. ( NFILE_LIBRARY > 0 ) ) then;                 !              
                                                                                         !
        deallocate(DIHEDRAL_COEFFS_LIBRARY);                                             !
                                                                                         !
        deallocate(MIDDLEBONDTORSION_COEFFS_LIBRARY);                                    !
                                                                                         !
        deallocate(ENDBONDTORSION_COEFFS_LIBRARY);                                       !
                                                                                         !
        deallocate(ANGLETORSION_COEFFS_LIBRARY);                                         !
                                                                                         !
        deallocate(ANGLEANGLETORSION_COEFFS_LIBRARY);                                    !
                                                                                         !
        deallocate(BONDBOND13_COEFFS_LIBRARY);                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( ( IMAX_NTYPE_IMPROPER > 0 ) .AND. ( NFILE_LIBRARY > 0 ) ) then;                 !
                                                                                         !
        deallocate(IMPROPER_COEFFS_LIBRARY);                                             !
                                                                                         !
        deallocate(ANGLEANGLE_COEFFS_LIBRARY);                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine DEALLOCATE_INSERT_LIBRARY_ARRAYS
