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

subroutine CHECK_FILE_LIBRARY(icanal) 

!   ************************************************************************************************
!   **                                  CHECK FILE LIBRARY                                        **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                : CANAL ON WHICH OUTPUT DATA ARE WRITTEN                             **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

    use module_library;

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ***********************************************************************************************

    integer (kind=4), intent(in) :: icanal;

!   ************************************************************************************************

    integer (kind=4) :: i, j;

    character (len=250) :: CHNAME_LIBRARY_DIRECTORY, CHNAME_TMP1, CHNAME_TMP2, CHNAME_TMP3;

    logical :: PROBE1, PROBE2, PROBE3, PROBE4;

!   ************************************************************************************************

!   integer (kind=4) :: IOSEF1, IOSEF2;

!   character (len=250) :: CHOSEF1;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the routine #############################################################
                                                                                         ! 
    CHTITLE = 'Check File Library';                                                      !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Check the presence of the library ##########################################################
                                                                                         !
    CHNAME_LIBRARY_DIRECTORY = 'LIB-MOLECULAR-MODELS';                                   !
                                                                                         !
    open(22,file='Tempo');                                                               !
                                                                                         !
    write(22,'(a11)') '#!/bin/bash';                                                     !
                                                                                         !
    write(22,*);                                                                         !
                                                                                         !
    write(22,*) 'if [ -d "'//                    &                                       !
                TRIM(CHNAME_LIBRARY_DIRECTORY)// &                                       !
                '" ]; then touch EXIST; fi';                                             !
                                                                                         !
    close(22);                                                                           !
                                                                                         !
    call system('chmod +x Tempo');                                                       !
                                                                                         !
    call system('./Tempo');                                                              !
                                                                                         !
    inquire(FILE='EXIST',EXIST=PROBE1);                                                  !
                                                                                         !
    write(icanal,'(a11,L1,a58)') '| PROBE1 : ', PROBE1, REPEAT(' ',57)//'|';             !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '| The library directory was not found in the '// &            !
                          'current directory'//REPEAT(' ',7)//'|';                       !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    call system('rm -f EXIST Tempo');                                                    !
                                                                                         !
    if ( PROBE1 .EQV. .FALSE. ) then;                                                    ! IF THE LIBRARY DIRECTORY HAS NOT BEEN FOUND IN THE CURRENT DIRECTORY, THEN
                                                                                         !
        do i = 1, 5;                                                                     ! LOOP OVER THE PARENT DIRECTORIES
                                                                                         !
            CHNAME_LIBRARY_DIRECTORY = '../'//TRIM(CHNAME_LIBRARY_DIRECTORY);            ! BUILD THE PATH TO THE DIRECTORY
                                                                                         !
            open(22,file='Tempo');                                                       ! CREATE A TEMPORARY SCRIPT TO CHECK FILES
                                                                                         !
            write(22,*) '#!/bin/bash';                                                   !
                                                                                         !
            write(22,*);                                                                 !
                                                                                         !
            write(22,*) 'if [ -d "'//                    &                               !
                        TRIM(CHNAME_LIBRARY_DIRECTORY)// &                               !
                        '" ]; then touch EXIST; fi';                                     !
                                                                                         !
            close(22);                                                                   !
                                                                                         !
            call system('chmod +x Tempo');                                               !
                                                                                         !
            call system('./Tempo');                                                      !
                                                                                         !
            inquire(FILE='EXIST',EXIST=PROBE2);                                          !
                                                                                         !
            call system('rm -f EXIST Tempo');                                            !
                                                                                         !
            if ( PROBE2 .EQV. .TRUE. ) EXIT;                                             !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
    write(icanal,'(a70)') '| The library directory was found at the '// &                !
                          'following path : '//REPEAT(' ',11)//'|';                      !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    IOSEF1 = LEN_TRIM(CHNAME_LIBRARY_DIRECTORY);                                         !
                                                                                         !
    IOSEF2 = 67 - IOSEF1;                                                                !
                                                                                         !
    if ( IOSEF2 < 0 ) IOSEF2 = 1;                                                        !
                                                                                         !
    write(icanal,'(a70)') '| '//TRIM(CHNAME_LIBRARY_DIRECTORY)//REPEAT(' ',IOSEF2)//'|'; !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Check directory in the library for molecules to insert #####################################
                                                                                         !
    if ( NFILE_LIBRARY > 0 ) then;                                                       ! IF THERE IS MOLECULES TO INSERT IN THE SIMULATION BOX, THEN
                                                                                         !
        open(22,file='Tempo');                                                           !
                                                                                         !
        write(22,'(a11)') '#!/bin/bash';                                                 !
                                                                                         !
        write(22,*);                                                                     !
                                                                                         !
        CHOSEF1 = 'if [ -d "'//TRIM(CHNAME_LIBRARY_DIRECTORY);                           !
                                                                                         !
        if ( TRIM(CHFILE_FORMAT) == 'xyz' ) then;                                        !
                                                                                         !
            CHOSEF1 = TRIM(CHOSEF1)//'/XYZ';                                             !
                                                                                         !
        else if ( TRIM(CHFILE_FORMAT) == 'lammps' ) then;                                !
                                                                                         !
            CHOSEF1 = TRIM(CHOSEF1)//'/LAMMPS';                                          !
                                                                                         !
        else if ( TRIM(CHFILE_FORMAT) == 'padua' ) then;                                 !
                                                                                         !
            CHOSEF1 = TRIM(CHOSEF1)//'/Padua';                                           !
                                                                                         !
        else if ( TRIM(CHFILE_FORMAT) == 'lammps-gaff' ) then;                           !
                                                                                         !
            CHOSEF1 = TRIM(CHOSEF1)//'/lammps-gaff';                                     !
                                                                                         !
        else                                                                             !
                                                                                         !
            write(icanal,'(a70)') '| /!\ Error /!\'//REPEAT(' ',54)//'|';                !
                                                                                         ! 
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            write(icanal,'(a70)') '| The file format asked by the user is '// &          !
                                  'not recognized by this program |';                    !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                              !
                                                                                         !
            write(icanal,*);                                                             !
                                                                                         !
            write(icanal,'(a14)') 'End of program';                                      !
                                                                                         !
            close(icanal);                                                               !
                                                                                         !
        end if                                                                           !
                                                                                         !
        CHOSEF1 = TRIM(CHOSEF1)//'" ]; then touch EXIST; fi';                            !
                                                                                         !
        write(22,*) TRIM(CHOSEF1);                                                       !
                                                                                         !
        close(22);                                                                       !
                                                                                         !
        call system('chmod +x Tempo');                                                   !
                                                                                         !
        call system('./Tempo');                                                          !
                                                                                         !
        inquire(FILE='EXIST',EXIST=PROBE2);                                              !
                                                                                         !
        call system('rm -f  EXIST Tempo');                                               !
                                                                                         !
        if ( PROBE2 .EQV. .FALSE. ) then;                                                !
                                                                                         !
            stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        else                                                                             !
                                                                                         !
            if ( TRIM(CHFILE_FORMAT) == 'xyz' ) then;                                    !
                                                                                         !
                CHNAME_LIBRARY_DIRECTORY = TRIM(CHNAME_LIBRARY_DIRECTORY)//'/XYZ/';      !
                                                                                         !
            else if ( TRIM(CHFILE_FORMAT) == 'lammps' ) then;                            !
                                                                                         !
                CHNAME_LIBRARY_DIRECTORY = TRIM(CHNAME_LIBRARY_DIRECTORY)//'/LAMMPS/';   !
                                                                                         !
            else if ( TRIM(CHFILE_FORMAT) == 'padua' ) then;                             !
                                                                                         !
                CHNAME_LIBRARY_DIRECTORY = TRIM(CHNAME_LIBRARY_DIRECTORY)//'/Padua/';    !
                                                                                         !
            else if ( TRIM(CHFILE_FORMAT) == 'lammps-gaff' ) then;                       !
                                                                                         !
                CHNAME_LIBRARY_DIRECTORY = TRIM(CHNAME_LIBRARY_DIRECTORY)// &            !
                                           '/lammps-gaff/';                              !
                                                                                         !
            end if                                                                       !
                                                                                         !
        end if                                                                           ! 
                                                                                         !
        IOSEF1 = LEN_TRIM(CHNAME_LIBRARY_DIRECTORY);                                     !
                                                                                         !
        IOSEF2 = 67 - IOSEF1;                                                            !
                                                                                         !
        if ( IOSEF2 < 0 ) IOSEF2 = 1;                                                    ! 
                                                                                         !
        write(icanal,'(a70)') '| '//TRIM(CHNAME_LIBRARY_DIRECTORY)// &                   !
                              REPEAT(' ',IOSEF2)//'|';                                   !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
        do i = 1, NFILE_LIBRARY;                                                         ! LOOP OVER THE NFILES TO READ IN THE LIBRARY
                                                                                         ! 
            IOSEF1 = LEN_TRIM(CHNAME_FILE_LIBRARY(1,i));                                 !
                                                                                         !
            IOSEF2 = 67 - IOSEF1;                                                        !
                                                                                         !
            if ( IOSEF2 < 0 ) IOSEF2 = 1;                                                !
                                                                                         !
            write(icanal,'(a70)') '| '//TRIM(CHNAME_FILE_LIBRARY(1,i))// &               !
                                  REPEAT(' ',IOSEF2)//'|';                               !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            CHNAME_TMP1 = TRIM(CHNAME_LIBRARY_DIRECTORY)// &                             !
                          TRIM(CHNAME_FILE_LIBRARY(1,i));                                !
                                                                                         !
!           ### Set extension of the file for lammps templates #####################################
                                                                                         !
            if ( ( TRIM(CHFILE_FORMAT) == 'lammps'      ) .OR.   &                       !
                 ( TRIM(CHFILE_FORMAT) == 'padua'       ) .OR.   &                       !
                 ( TRIM(CHFILE_FORMAT) == 'lammps-gaff' ) ) then;                        !
                                                                                         !
                CHNAME_TMP1 = TRIM(CHNAME_TMP1)//'.template';                            !
                                                                                         !
            end if                                                                       !
                                                                                         !
!           ### Set extension of the file for the xyz format #######################################
                                                                                         !
            if ( TRIM(CHFILE_FORMAT) == 'xyz' ) then;                                    !
                                                                                         !
                CHNAME_TMP1 = TRIM(CHNAME_TMP1)//'.xyz';                                 !
                                                                                         !
            end if                                                                       !
                                                                                         !
            IOSEF1 = LEN_TRIM(CHNAME_TMP1);                                              !
                                                                                         !
            IOSEF2 = 67 - IOSEF1;                                                        !
                                                                                         !
            if ( IOSEF2 < 0 ) IOSEF2 = 1;                                                !
                                                                                         !
            write(icanal,'(a70)') '| '//TRIM(CHNAME_TMP1)//REPEAT(' ',IOSEF2)//'|';      !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
            inquire(FILE=TRIM(CHNAME_TMP1),EXIST=PROBE3);                                !
                                                                                         !
            if ( PROBE3 .EQV. .TRUE. ) then;                                             !
                                                                                         !
                write(icanal,'(a70)') '| The file has been found in the library'// &     !
                                      REPEAT(' ',29)//'|';                               !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
                if ( ( TRIM(CHFILE_FORMAT) == 'lammps'      ) .OR.  &                    !
                     ( TRIM(CHFILE_FORMAT) == 'padua'       ) .OR.  &                    !
                     ( TRIM(CHFILE_FORMAT) == 'lammps-gaff' )  ) then;                   !
                                                                                         !
                    CHNAME_TMP2 = TRIM(CHNAME_LIBRARY_DIRECTORY)// &                     !
                                  TRIM(CHNAME_FILE_LIBRARY(1,i))// &                     !
                                  '.parameter';                                          !
                                                                                         !
                    inquire(FILE=TRIM(CHNAME_TMP2),EXIST=PROBE4);                        !
                                                                                         !
                    IOSEF1 = LEN_TRIM(CHNAME_TMP2);                                      !
                                                                                         !
                    IOSEF2 = 67 - IOSEF1;                                                !
                                                                                         !
                    if ( IOSEF2 < 0 ) IOSEF2 = 1;                                        !
                                                                                         !
                    write(icanal,'(a70)') '| '//TRIM(CHNAME_TMP2)// &                    !
                                          REPEAT(' ',IOSEF2)//'|';                       !
                                                                                         !
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
!                   stop; !//////////////////////////////////////////////////////////////!
                                                                                         !
                    if ( PROBE4 .EQV. .FALSE. ) then;                                    !
                                                                                         !
                        write(icanal,'(a70)') '| The file has not been found '//     &   !
                                              'in the library'//REPEAT(' ',25)//'|';     !
                                                                                         !
                        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                  !
                                                                                         !
                        write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                  !
                                                                                         !
                        write(icanal,*);                                                 !
                                                                                         !
                        write(icanal,'(a14)') 'End of program';                          !
                                                                                         !
                        close(icanal);                                                   !
                                                                                         !
                        stop; !//////////////////////////////////////////////////////////!
                                                                                         !
                    else                                                                 !
                                                                                         !
                        write(icanal,'(a70)') '| The file has been found '//         &   !
                                              'in the library'//REPEAT(' ',29)//'|';     !
                                                                                         !
                        CHNAME_FILE_LIBRARY(1,i) = TRIM(CHNAME_LIBRARY_DIRECTORY)// &    !
                                                   TRIM(CHNAME_FILE_LIBRARY(1,i));       !
                                                                                         !
                    end if                                                               !
                                                                                         !
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
!                   stop; !//////////////////////////////////////////////////////////////!
                                                                                         !
                else                                                                     !
                                                                                         !
                    CHNAME_FILE_LIBRARY(1,i) = TRIM(CHNAME_LIBRARY_DIRECTORY)// &        !
                                               TRIM(CHNAME_FILE_LIBRARY(1,i));           !
                                                                                         !
                end if                                                                   !
                                                                                         !
            else                                                                         !
                                                                                         !
                write(icanal,'(a70)') '| The file has not been found '//     &           !
                                      'in the library'//REPEAT(' ',25)//'|';             !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                          !
                                                                                         !
                write(icanal,*);                                                         !
                                                                                         !
                write(icanal,'(a14)') 'End of program';                                  !
                                                                                         !
                close(icanal);                                                           !
                                                                                         !
                stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Alocate array containing information about the presence of additional info file ############
                                                                                         !
    allocate(ILIBRARY_FILE_INFO(1:NFILE_LIBRARY));                                       !
                                                                                         !
!   ### Initialization of the array containing information about the presence of info file #########
                                                                                         !
    ILIBRARY_FILE_INFO(1:NFILE_LIBRARY) = 0;                                             !
                                                                                         !
!   ### Check for the additional information file ##################################################
                                                                                         !
    if ( NFILE_LIBRARY > 0 ) then;                                                       !
                                                                                         !
        if ( ( TRIM(CHFILE_FORMAT) == 'lammps'      ) .OR.   &                           !
             ( TRIM(CHFILE_FORMAT) == 'padua'       ) .OR.   &                           !
             ( TRIM(CHFILE_FORMAT) == 'lammps-gaff' ) ) then;                            !
                                                                                         !
            do i = 1, NFILE_LIBRARY;                                                     ! 
                                                                                         !
                CHNAME_TMP3 = TRIM(CHNAME_FILE_LIBRARY(1,i))//   &                       !
                               '.info';                                                  !
                                                                                         !
                IOSEF1 = LEN_TRIM(CHNAME_TMP3);                                          !
                                                                                         !
                IOSEF2 = 67 - IOSEF1;                                                    !
                                                                                         !
                if ( IOSEF2 < 0 ) IOSEF2 = 1;                                            !
                                                                                         !
                write(icanal,'(a70)') '| '//TRIM(CHNAME_TMP3)//REPEAT(' ',IOSEF2)//'|';  !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
                inquire(FILE=TRIM(CHNAME_TMP3),EXIST=PROBE1);                            !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
                if ( PROBE1 .EQV. .FALSE. ) then;                                        !
                                                                                         !
                    write(icanal,'(a70)') '| The file has not been found '//     &       !
                                          'in the library'//REPEAT(' ',25)//'|';         !
                                                                                         !
!                   stop; !//////////////////////////////////////////////////////////////!
                                                                                         !
                else                                                                     !
                                                                                         !
                    write(icanal,'(a70)') '| The file has been found '// &               !
                                          'in the library'//REPEAT(' ',29)//'|';         !
                                                                                         !
                    ILIBRARY_FILE_INFO(i) = 1;                                           !
                                                                                         !
                end if                                                                   !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end do                                                                       !
                                                                                         !
        end if                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Check library directories for molecules to remove ##########################################
                                                                                         !
    if ( NFILE_LIBRARY_REMOVE > 0 ) then;                                                ! IF THE NUMBER OF FILES CORRESPONDING TO THE MOLECULES TO BE REMOVED IS GREATER THAN 0
                                                                                         !
        write(icanal,'(a70)') '| CHECK LIBRARIES FOR REMOVING FILES'// &                 !
                              REPEAT(' ',33)//'|';                                       !
                                                                                         !
        if ( TRIM(CHFILE_FORMAT) == 'xyz' ) then;                                        !
                                                                                         !
            open(22,file='Tempo');                                                       !
            write(22,'(a11)') '#!/bin/bash';                                             !
            write(22,*);                                                                 !
            write(22,*) 'if [ -d "'//                    &                               !
                        TRIM(CHNAME_LIBRARY_DIRECTORY)// &                               !
                        '/XYZ" ]; then touch EXIST; fi';                                 !
            close(22);                                                                   !
                                                                                         !
            call system('chmod +x Tempo');                                               !
                                                                                         !
            call system('./Tempo');                                                      !
                                                                                         !
            inquire(FILE='EXIST',EXIST=PROBE2);                                          !
                                                                                         !
            call system('rm -f  EXIST  Tempo');                                          !
                                                                                         !
            if ( PROBE2 .EQV. .FALSE. ) then;                                            !
                stop;                                                                    !
            else                                                                         !
                CHNAME_LIBRARY_DIRECTORY = &                                             !
                TRIM(CHNAME_LIBRARY_DIRECTORY)//'/XYZ/';                                 !
            end if                                                                       !
                                                                                         !
        else if ( TRIM(CHFILE_FORMAT) == 'lammps' ) then;                                !
                                                                                         !
            open(22,file='Tempo');                                                       !
            write(22,'(a11)') '#!/bin/bash';                                             !
            write(22,*);                                                                 !
            write(22,*) 'if [ -d "'//                    &                               !
                        TRIM(CHNAME_LIBRARY_DIRECTORY)// &                               !
                        '/LAMMPS" ]; then touch EXIST; fi';                              !
            close(22);                                                                   !
                                                                                         !
            call system('chmod +x Tempo');                                               !
                                                                                         !
            call system('./Tempo');
                                                                                         !
            inquire(FILE='EXIST',EXIST=PROBE2);                                          !
                                                                                         !
            call system('rm -f  EXIST Tempo');                                           !
                                                                                         !
            if ( PROBE2 .EQV. .FALSE. ) then;                                            !
                stop;                                                                    !
            else                                                                         !
                CHNAME_LIBRARY_DIRECTORY = &                                             !
                TRIM(CHNAME_LIBRARY_DIRECTORY)//'/LAMMPS/';                              !
            end if                                                                       !
        end if                                                                           !
                                                                                         !
        IOSEF1 = LEN_TRIM(CHNAME_LIBRARY_DIRECTORY);                                     !
                                                                                         !
        IOSEF2 = 67 - IOSEF1;                                                            !
                                                                                         !
        if ( IOSEF2 < 0 ) IOSEF2 = 1;                                                    !
                                                                                         !
        write(icanal,'(a70)') '| '//                           &                         !
                              TRIM(CHNAME_LIBRARY_DIRECTORY)// &                         !
                              REPEAT(' ',IOSEF2)//'|';                                   !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        do i = 1, NFILE_LIBRARY_REMOVE;                                                  !
                                                                                         ! 
            IOSEF1 = LEN_TRIM(CHNAME_FILE_LIBRARY_REMOVE(i));                            !
                                                                                         !
            IOSEF2 = 67 - IOSEF1;                                                        !
                                                                                         !
            if ( IOSEF2 < 0 ) IOSEF2 = 1;                                                !
                                                                                         !
            write(icanal,'(a70)') '| '//                                &                !
                                  TRIM(CHNAME_FILE_LIBRARY_REMOVE(i))// &                !
                                  REPEAT(' ',IOSEF2)//'|';                               !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            CHNAME_TMP1 =                       &                                        !
            TRIM(CHNAME_LIBRARY_DIRECTORY)//    &                                        !
            TRIM(CHNAME_FILE_LIBRARY_REMOVE(i));                                         !
                                                                                         !
            if ( TRIM(CHFILE_FORMAT) == 'lammps' ) CHNAME_TMP1 = &                       !
                                                   TRIM(CHNAME_TMP1)// &
                                                   '.template';

            IOSEF1 = LEN_TRIM(CHNAME_TMP1);

            IOSEF2 = 67 - IOSEF1;

            if ( IOSEF2 < 0 ) IOSEF2 = 1;

            write(icanal,'(a70)') '| '//TRIM(CHNAME_TMP1)//REPEAT(' ',IOSEF2)//'|';

            inquire(FILE=TRIM(CHNAME_TMP1),EXIST=PROBE3);

            if ( PROBE3 .EQV. .TRUE. ) then
                write(icanal,'(a70)') '| THE FILE HAS BEEN FOUND IN THE LIBRARY'// &
                                      REPEAT(' ',29)//'|';
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

                if ( TRIM(CHFILE_FORMAT) == 'lammps' ) then;

                    CHNAME_TMP2 = &
                    TRIM(CHNAME_LIBRARY_DIRECTORY)// &
                    TRIM(CHNAME_FILE_LIBRARY_REMOVE(i))//'.parameter';

                    inquire(FILE=TRIM(CHNAME_TMP2),EXIST=PROBE4);

                    IOSEF1 = LEN_TRIM(CHNAME_TMP2);

                    IOSEF2 = 67 - IOSEF1;

                    if ( IOSEF2 < 0 ) IOSEF2 = 1;

                    write(icanal,'(a70)') '| '// &
                                          TRIM(CHNAME_TMP2)// &
                                          REPEAT(' ',IOSEF2)//'|';

                    if ( PROBE4 .EQV. .FALSE. ) then;

                        write(icanal,'(a70)') '| THE FILE HAS NOT BEEN FOUND IN THE LIBRARY'// &
                                              REPEAT(' ',25)//'|';
                        stop;
                    else
                        write(icanal,'(a70)') '| THE FILE HAS BEEN FOUND IN THE LIBRARY'// &
                                              REPEAT(' ',29)//'|';

                        CHNAME_FILE_LIBRARY_REMOVE(i) =  &
                        TRIM(CHNAME_LIBRARY_DIRECTORY)// &
                        TRIM(CHNAME_FILE_LIBRARY_REMOVE(i));
                    end if

                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

                else
                    CHNAME_FILE_LIBRARY_REMOVE(i) = &
                    TRIM(CHNAME_LIBRARY_DIRECTORY)// &
                    TRIM(CHNAME_FILE_LIBRARY_REMOVE(i));

                end if
            else
                write(icanal,'(a70)') '| THE FILE HAS NOT BEEN FOUND IN THE LIBRARY'// &
                                      REPEAT(' ',25)//'|';
                stop;

            end if

        end do

    end if

!   ### Check the library for the generation of polymer chains #####################################
                                                                                         !
    if ( NGenerate_polymer_species > 0 ) then;                                           !
                                                                                         !
        open(22,file='Tempo');                                                           !
                                                                                         !
        write(22,'(a11)') '#!/bin/bash';                                                 !
                                                                                         !
        write(22,*);                                                                     !
                                                                                         !
        CHOSEF1 = 'if [ -d "'//TRIM(CHNAME_LIBRARY_DIRECTORY);                           !
                                                                                         !
        if ( TRIM(Chfile_format_polymer) == 'xyz-monomers' ) then;                       !
                                                                                         !
            CHOSEF1 = TRIM(CHOSEF1)//'/xyz-monomers';                                    !
                                                                                         !
        end if                                                                           !
                                                                                         !
        CHOSEF1 = TRIM(CHOSEF1)//'" ]; then touch EXIST; fi';                            !
                                                                                         !
        write(22,*) TRIM(CHOSEF1);                                                       !
                                                                                         !
        close(22);                                                                       !
                                                                                         !
        call system('chmod +x Tempo');                                                   !
                                                                                         !
        call system('./Tempo');                                                          !
                                                                                         !
        inquire(FILE='EXIST',EXIST=PROBE2);                                              !
                                                                                         !
        call system('rm -f  EXIST Tempo');                                               !
                                                                                         !
        if ( PROBE2 .EQV. .FALSE. ) then;                                                !
                                                                                         !
            stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        else                                                                             !
                                                                                         !
            if ( TRIM(Chfile_format_polymer) == 'xyz-monomers' ) then;                   !
                                                                                         !
                CHNAME_LIBRARY_DIRECTORY = TRIM(CHNAME_LIBRARY_DIRECTORY)// &            !
                                           '/xyz-monomers/';                             !
                                                                                         !
            end if                                                                       !
                                                                                         !
        end if                                                                           ! 
                                                                                         !
        IOSEF1 = LEN_TRIM(CHNAME_LIBRARY_DIRECTORY);                                     !
                                                                                         !
        IOSEF2 = 67 - IOSEF1;                                                            !
                                                                                         !
        if ( IOSEF2 < 0 ) IOSEF2 = 1;                                                    !
                                                                                         !
        write(icanal,'(a70)') '| '//TRIM(CHNAME_LIBRARY_DIRECTORY)// &                   !
                              REPEAT(' ',IOSEF2)//'|';                                   !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!       do i = 1, NGenerate_polymer_species;                                             !
        do i = 1, NTotal_monomers;                                                       !
                                                                                         !
!           do j = 1, NSPECIES_PER_POLYMER_CHAINS(i);                                    !
                                                                                         !
!               IOSEF1 = LEN_TRIM(Chname_file_library_monomer(i,j));                     !
                IOSEF1 = LEN_TRIM(Chname_file_library_monomer(i));                       !
                                                                                         !
                IOSEF2 = 67 - IOSEF1;                                                    !
                                                                                         !
                write(icanal,'(a70)') '| '//TRIM(Chname_file_library_monomer(i))// &     !
                                      REPEAT(' ',IOSEF2)//'|';                           !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                CHNAME_TMP1 = TRIM(CHNAME_LIBRARY_DIRECTORY)//        &                  !
                              TRIM(Chname_file_library_monomer(i))//'.xyz';              !
                                                                                         !
                IOSEF1 = LEN_TRIM(CHNAME_TMP1);                                          !
                                                                                         !
                IOSEF2 = 67 - IOSEF1;                                                    !
                                                                                         !
                if ( IOSEF2 < 0 ) IOSEF2 = 1;                                            !
                                                                                         !
                write(icanal,'(a70)') '| '//TRIM(CHNAME_TMP1)//REPEAT(' ',IOSEF2)//'|';  !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
                inquire(FILE=TRIM(CHNAME_TMP1),EXIST=PROBE3);                            !
                                                                                         !
                if ( PROBE3 .EQV. .TRUE. ) then;                                         !
                                                                                         !
                    write(icanal,'(a70)') '| The file has been found in the library'// & !
                                          REPEAT(' ',29)//'|';                           !
                                                                                         !
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
!                   ### Add the path to the library directory to the library file ##################
                                                                                         !
                    Chname_file_library_monomer(i)   =      &                            !
                    TRIM(CHNAME_LIBRARY_DIRECTORY)//        &                            !
                    TRIM(Chname_file_library_monomer(i));                                !
                                                                                         !
!                   stop; !//////////////////////////////////////////////////////////////!
                                                                                         !
                else                                                                     !
                                                                                         !
                    write(icanal,'(a70)') '| The file has not been found '//     &       !
                                          'in the library'//REPEAT(' ',25)//'|';         !
                                                                                         !
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
                    stop; !//////////////////////////////////////////////////////////////!
                                                                                         !
                end if                                                                   !
                                                                                         !
!           end do                                                                       !
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
end subroutine CHECK_FILE_LIBRARY











