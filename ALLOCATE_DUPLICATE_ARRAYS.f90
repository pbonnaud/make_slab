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

subroutine ALLOCATE_DUPLICATE_ARRAYS(icanal)

!   ************************************************************************************************
!   **                                READ LAMMPS TEMPLATE FILES                                  **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                           **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

    use module_library;

    use module_duplicate;

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ***********************************************************************************************

    integer (kind=4), intent(in) :: icanal;

!   ************************************************************************************************

    integer (kind=4) :: i, j, k;

    integer (kind=4) :: IFOUND;

    character (len=20), dimension(1:100) :: LOCAL_ATOM_LABEL;

!   ************************************************************************************************

    integer (kind=4) :: NFILE_LIBRARY_TMP;    

!   ************************************************************************************************

    integer (kind=4) :: IOSEF_NATOM,     &
                        IOSEF_NBOND,     &
                        IOSEF_NANGLE,    &
                        IOSEF_NDIHEDRAL, &
                        IOSEF_NIMPROPER;

    integer (kind=4) :: IOSEF_NTYPE_NATOM,    &
                        IOSEF_NTYPE_BOND,     &
                        IOSEF_NTYPE_ANGLE,    &
                        IOSEF_NTYPE_DIHEDRAL, &
                        IOSEF_NTYPE_IMPROPER;

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
    CHTITLE = 'Allocate duplicate arrays';                                               !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialization of arrays for cell dimensions (duplicate) ###################################
                                                                                         !
    CELL_AXIS_DUPLICATE(1:3) = 0.0d0;                                                    !
                                                                                         !
    CELL_ANGDEG_DUPLICATE(1:3) = 0.0d0;                                                  !
                                                                                         !
!   ### Initialization duplicate arrays ############################################################
                                                                                         !
    NATOM_DUPLICATE = 0;                                                                 !
                                                                                         !
    NBOND_DUPLICATE = 0;                                                                 !
                                                                                         !
    NANGLE_DUPLICATE = 0;                                                                !
                                                                                         !
    NDIHEDRAL_DUPLICATE = 0;                                                             !
                                                                                         !
    NIMPROPER_DUPLICATE = 0;                                                             !
                                                                                         !
    NTYPE_ATOM_DUPLICATE = 0;                                                            !
                                                                                         ! 
    NTYPE_BOND_DUPLICATE = 0;                                                            !
                                                                                         !
    NTYPE_ANGLE_DUPLICATE = 0;                                                           !
                                                                                         ! 
    NTYPE_DIHEDRAL_DUPLICATE = 0;                                                        !
                                                                                         !
    NTYPE_IMPROPER_DUPLICATE = 0;                                                        !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         ! 
!   ### Initialization of parameter arrays #########################################################
                                                                                         !
    NPARAM_BONDS_DUPLICATE = 0;                                                          !
                                                                                         !
    NPARAM_ANGLES_DUPLICATE = 0;                                                         !
                                                                                         !
    NPARAM_DIHEDRALS_DUPLICATE = 0;                                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Read molecule properties in the library files  with the lammps format ######################
                                                                                         !
    if ( ( TRIM(DUPLICATE_CHFILE_FORMAT) == 'lammps'      ) .OR.  &                      !
         ( TRIM(DUPLICATE_CHFILE_FORMAT) == 'padua'       ) .OR.  &                      !
         ( TRIM(DUPLICATE_CHFILE_FORMAT) == 'lammps-gaff' ) ) then;                      !
                                                                                         !
        if ( NFILE_LIBRARY > 0 ) then;                                                   !
                                                                                         !
            CHEXT = TRIM(DUPLICATE_CHNAME_FILE)//'.template';                            !
                                                                                         !
            IOSEF1 = LEN_TRIM(CHEXT);                                                    !
                                                                                         !
            IOSEF2 = 67 - IOSEF1;                                                        !
                                                                                         !
            if ( IOSEF2 < 0 ) IOSEF2 = 1;                                                !
                                                                                         !
            write(icanal,'(a70)') '| '//TRIM(CHEXT)//REPEAT(' ',IOSEF2)//'|';            !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           ### Initialization of local variables containing properties of the molecule ############
                                                                                         !
            EOF = 0;                                                                     !  
                                                                                         !
!           ### Reading the current library file ###################################################
                                                                                         !
            open(2,file=TRIM(CHEXT));                                                    !
                                                                                         !
            do                                                                           !
                                                                                         !
                read(2,'(a)',iostat=EOF) CHARLINE;                                       !
                                                                                         !
                if ( EOF /= 0 ) EXIT;                                                    !
                                                                                         !
                if ( INDEX(CHARLINE,'atoms') > 0 ) then;                                 !
                                                                                         !
                    read(CHARLINE,*) NATOM_DUPLICATE;                                    !
                                                                                         !
                else if ( INDEX(CHARLINE,'bonds') > 0 ) then;                            !
                                                                                         !
                    read(CHARLINE,*) NBOND_DUPLICATE;                                    !
                                                                                         ! 
                else if ( INDEX(CHARLINE,'angles') > 0 ) then;                           !
                                                                                         !
                    read(CHARLINE,*) NANGLE_DUPLICATE;                                   !
                                                                                         !
                else if ( INDEX(CHARLINE,'dihedrals') > 0 ) then;                        !
                                                                                         !
                    read(CHARLINE,*) NDIHEDRAL_DUPLICATE;                                !
                                                                                         !
                else if ( INDEX(CHARLINE,'impropers') > 0 ) then;                        !
                                                                                         !
                    read(CHARLINE,*) NIMPROPER_DUPLICATE;                                !
                                                                                         !
                else if ( INDEX(CHARLINE,'Types') > 0 ) then;                            !
                                                                                         !
                    read(2,*);

                    do j = 1, NATOM_DUPLICATE;

                        read(2,*) IOSEF1, IOSEF2;

                        if ( IOSEF2 > NTYPE_ATOM_DUPLICATE ) NTYPE_ATOM_DUPLICATE = IOSEF2;

                    end do

                else if ( INDEX(CHARLINE,'Bonds') > 0 ) then;

                    read(2,*);

                    do j = 1, NBOND_DUPLICATE;

                        read(2,*) IOSEF1, IOSEF2;

                        if ( IOSEF2 > NTYPE_BOND_DUPLICATE ) NTYPE_BOND_DUPLICATE = IOSEF2;

                    end do

                else if ( INDEX(CHARLINE,'Angles') > 0 ) then;

                    read(2,*);

                    do j = 1, NANGLE_DUPLICATE;

                        read(2,*) IOSEF1, IOSEF2;

                        if ( IOSEF2 > NTYPE_ANGLE_DUPLICATE ) NTYPE_ANGLE_DUPLICATE = IOSEF2;

                    end do

                else if ( INDEX(CHARLINE,'Dihedrals') > 0 ) then;

                    read(2,*);

                    do j = 1, NDIHEDRAL_DUPLICATE;

                        read(2,*) IOSEF1, IOSEF2;

                        if ( IOSEF2 > NTYPE_DIHEDRAL_DUPLICATE ) NTYPE_DIHEDRAL_DUPLICATE = IOSEF2;

                    end do

                else if ( INDEX(CHARLINE,'Impropers') > 0 ) then;

                    read(2,*);

                    do j = 1, NIMPROPER_DUPLICATE;

                        read(2,*) IOSEF1, IOSEF2;

                        if ( IOSEF2 > NTYPE_IMPROPER_DUPLICATE ) NTYPE_IMPROPER_DUPLICATE = IOSEF2;

                    end do

                end if

            end do

            close(2);

            if ( NATOM_DUPLICATE > 0 ) then;                                             !
                                                                                         !
                write(icanal,'(a25,i8,a37)') '| NATOM_DUPLICATE          : ', &          !
                                             NATOM_DUPLICATE,                 &          !
                                             REPEAT(' ',36)//'|';                    !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
            end if                                                                   !
                                                                                         !
            if ( NBOND_DUPLICATE > 0 ) then;                                             !
                                                                                         !
                write(icanal,'(a25,i8,a37)') '| NBOND_DUPLICATE          : ', &          !
                                             NBOND_DUPLICATE,                 &          !
                                             REPEAT(' ',36)//'|';                    !

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            end if

            if ( NANGLE_DUPLICATE > 0 ) then;

                write(icanal,'(a25,i8,a37)') '| NANGLE_DUPLICATE         : ', &
                                             NANGLE_DUPLICATE,                &
                                             REPEAT(' ',36)//'|';

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            end if

            if ( NDIHEDRAL_DUPLICATE > 0 ) then;

                write(icanal,'(a25,i8,a37)') '| NDIHEDRAL_DUPLICATE      : ', &
                                             NDIHEDRAL_DUPLICATE,             &
                                             REPEAT(' ',36)//'|';

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            end if

                if ( NIMPROPER_DUPLICATE > 0 ) then;

                    write(icanal,'(a25,i8,a37)') '| NIMPROPER_DUPLICATE      : ', &
                                                 NIMPROPER_DUPLICATE,             &
                                                 REPEAT(' ',36)//'|';

                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

                end if

                if ( NTYPE_ATOM_DUPLICATE > 0 ) then;                                    !
                                                                                         !
                    write(icanal,'(a28,i8,a34)') '| NTYPE_ATOM_DUPLICATE    : ', &       !
                                                 NTYPE_ATOM_DUPLICATE,           &       !
                                                 REPEAT(' ',33)//'|';                    !
                                                                                         !
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
                end if                                                                   !

                if ( NTYPE_BOND_DUPLICATE > 0 ) then;

                    write(icanal,'(a25,i8,a37)') '| NTYPE_BOND_DUPLICATE     : ', &
                                                 NTYPE_BOND_DUPLICATE,            &
                                                 REPEAT(' ',36)//'|';

                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

                end if

                if ( NTYPE_ANGLE_DUPLICATE > 0 ) then;

                    write(icanal,'(a25,i8,a37)') '| NTYPE_ANGLE_DUPLICATE    : ', &
                                                 NTYPE_ANGLE_DUPLICATE,           &
                                                 REPEAT(' ',36)//'|';

                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

                end if

                if ( NTYPE_DIHEDRAL_DUPLICATE > 0 ) then;

                    write(icanal,'(a25,i8,a37)') '| NTYPE_DIHEDRAL_DUPLICATE : ', &
                                                 NTYPE_DIHEDRAL_DUPLICATE,        &
                                                 REPEAT(' ',36)//'|';

                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

                end if

                if ( NTYPE_IMPROPER_DUPLICATE > 0 ) then;

                    write(icanal,'(a25,i8,a37)') '| NTYPE_IMPROPER_DUPLICATE : ', &
                                                 NTYPE_IMPROPER_DUPLICATE,        &
                                                 REPEAT(' ',36)//'|';

                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

                end if
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
        end if                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Read molecule properties in the library files with the xyz format ##########################
                                                                                         !
    if ( TRIM(DUPLICATE_CHFILE_FORMAT) == 'xyz' ) then;                                  !
                                                                                         !
        if ( NFILE_LIBRARY > 0 ) then;                                                   !
                                                                                         !
!               ### Build the name of the library file to read #####################################
                                                                                         !
                CHEXT = TRIM(DUPLICATE_CHNAME_FILE)//'.xyz';                             !
                                                                                         !
                IOSEF1 = LEN_TRIM(CHEXT);                                                !
                                                                                         !
                IOSEF2 = 67 - IOSEF1;                                                    !
                                                                                         !
                if ( IOSEF2 < 0 ) IOSEF2 = 1;                                            !
                                                                                         !
                write(icanal,'(a70)') '| '//TRIM(CHEXT)//REPEAT(' ',IOSEF2)//'|';        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
!               ### Initialization of local variables containing properties of the molecule ########
                                                                                         !
                EOF = 0;                                                                 !
                                                                                         !
!               ### Reading the current library file ###############################################
                                                                                         !
                LOCAL_ATOM_LABEL(1:100) = 'XXX';                                         !
                                                                                         !
                open(2,file=TRIM(CHEXT));                                                !
                                                                                         !
                read(2,*) NATOM_DUPLICATE;                                               !
                                                                                         !
                read(2,*) CELL_AXIS_DUPLICATE(1:3), CELL_ANGDEG_DUPLICATE(1:3);          !
                                                                                         !
                do j = 1, NATOM_DUPLICATE;                                               !
                                                                                         !
                    read(2,*) CHOSEF1;                                                   ! 
                                                                                         !
                    if ( NTYPE_ATOM_DUPLICATE == 0 ) then;                               !
                                                                                         !
                        NTYPE_ATOM_DUPLICATE = NTYPE_ATOM_DUPLICATE + 1;                 !
                                                                                         !
                        LOCAL_ATOM_LABEL(NTYPE_ATOM_DUPLICATE) = TRIM(CHOSEF1);          ! 
                                                                                         !
                    else                                                                 !
                                                                                         !
                        IFOUND = 0;                                                      !
                                                                                         !
                        do k = 1, NTYPE_ATOM_DUPLICATE;                                  !
                                                                                         !
                            if ( TRIM(LOCAL_ATOM_LABEL(k)) == TRIM(CHOSEF1) ) then;      !
                                                                                         !
                                IFOUND = 1;                                              !
                                                                                         !
                                EXIT;                                                    !
                                                                                         !
                            end if                                                       ! 
                                                                                         !
                        end do

                        if ( IFOUND == 1 ) CYCLE;

                        NTYPE_ATOM_DUPLICATE = NTYPE_ATOM_DUPLICATE + 1;

                        LOCAL_ATOM_LABEL(NTYPE_ATOM_DUPLICATE) = TRIM(CHOSEF1);

                    end if

                end do

                close(2);                                                                !
                                                                                         !
!               ### Write cell properties ##########################################################
                                                                                         !
                write(icanal,'(a70)') '| Cell properties'//REPEAT(' ',52)//'|';          !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !  
                write(icanal,'(a16,3f12.6,a18)') '| a, b, c [A] : ',       &             !
                                                 CELL_AXIS_DUPLICATE(1:3), &             !
                                                 REPEAT(' ',17)//'|';                    !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                write(icanal,'(a30,3f12.6,a4)') '| alpha, beta, gamma [deg.] : ', &      !
                                                CELL_ANGDEG_DUPLICATE(1:3),       &      !
                                                REPEAT(' ',3)//'|';                      !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               ### Write configuration properties #################################################
                                                                                         !
                write(icanal,'(a29,i8,a33)') '| NATOM_DUPLICATE          : ', &          !
                                             NATOM_DUPLICATE,                 &          !
                                             REPEAT(' ',32)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                write(icanal,'(a27,i8,a35)') '| NTYPE_ATOM_DUPLICATE    : ', &           !
                                             NTYPE_ATOM_DUPLICATE,           &           !
                                             REPEAT(' ',34)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
        end if                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocation of duplicate arrays #############################################################
                                                                                         !
    if ( NATOM_DUPLICATE > 0 ) then;                                                     !
                                                                                         !
        allocate(CONFIG_ATOMID_DUPLICATE(1:NATOM_DUPLICATE));                            !
                                                                                         !
        allocate(CONFIG_ATOM_TYPE_DUPLICATE(1:NATOM_DUPLICATE));                         !
                                                                                         !
!       allocate(CONFIG_QI_LIBRARY(1:IMAX_NATOM,1:NFILE_LIBRARY_TMP));                   !
                                                                                         !
!       allocate(CONFIG_VI_LIBRARY(1:3,1:IMAX_NATOM,1:NFILE_LIBRARY_TMP));               !
                                                                                         !
        allocate(CONFIG_RI_DUPLICATE(1:3,1:NATOM_DUPLICATE));                            !
                                                                                         !
        allocate(CONFIG_NAT_DUPLICATE(1:NATOM_DUPLICATE));                               !
                                                                                         !
        write(icanal,'(a70)') '| Duplicate arrays with the size '// &                    !
                              'of NATOM were allocated'//           &                    !
                              REPEAT(' ',13)//'|';                                       !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_ATOM_DUPLICATE > 0 ) then;                                                !
                                                                                         !
        allocate(ATOM_MASSES_DUPLICATE(1:NTYPE_ATOM_DUPLICATE));                         !
                                                                                         !
        allocate(ATOM_LABEL_DUPLICATE(1:NTYPE_ATOM_DUPLICATE));                          !
                                                                                         !
!       allocate(POTENTIAL_CLASS2_LIBRARY(1:2,1:IMAX_NTYPE_NATOM,1:NFILE_LIBRARY_TMP));  !
                                                                                         !
!       allocate(PAIR_COEFF_CROSS_LIBRARY(1:2,1:1000,1:NFILE_LIBRARY_TMP));              !
                                                                                         !
!       allocate(PAIR_ATOMID_CROSS_LIBRARY(1:2,1:1000,1:NFILE_LIBRARY_TMP));             !
                                                                                         !
        write(icanal,'(a70)') '| Duplicate arrays with the size of '// &                 !
                              'NTYPE_NATOM were allocated'//         &                   !
                              REPEAT(' ',7)//'|';                                        !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        ATOM_LABEL_DUPLICATE(1:NTYPE_ATOM_DUPLICATE) = &                                 !
        LOCAL_ATOM_LABEL(1:NTYPE_ATOM_DUPLICATE);                                        !
                                                                                         !
!       ### Set atom masses ########################################################################
                                                                                         !
        ATOM_MASSES_DUPLICATE(1:NTYPE_ATOM_DUPLICATE) = 0.0d0;                           !
                                                                                         !
        do i = 1, NTYPE_ATOM_DUPLICATE;                                                  !
                                                                                         !
            do j = 1, 200;                                                               !
                                                                                         !
                if ( TRIM(ATOM_LABEL_DUPLICATE(i)) ==     &                              !
                     TRIM(ELEMENT_SYMBOL_NAME(j)) ) then;                                !
                                                                                         !
                    ATOM_MASSES_DUPLICATE(i) = MOLAR_MASS_ELEMENT(j);                    !
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         ! 
        end do                                                                           !
                                                                                         !
        write(icanal,'(a70)') '| Masses were set '//REPEAT(' ',51)//'|';                 !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!       ### Write the list of atom labels and species ##############################################
                                                                                         !
        do i = 1, NTYPE_ATOM_DUPLICATE;                                                  !
                                                                                         !
            write(icanal,'(a6,i4,a5,a10,f12.6,a33)') '| i : ', i, ' --> ',     &         !
                                                     ATOM_LABEL_DUPLICATE(i),  &         !
                                                     ATOM_MASSES_DUPLICATE(i), &         !
                                                     REPEAT(' ',32)//'|';                !
                                                                                         !
        end do                                                                           !
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
end subroutine ALLOCATE_DUPLICATE_ARRAYS
