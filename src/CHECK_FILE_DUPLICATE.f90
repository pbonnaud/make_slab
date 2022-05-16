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

subroutine CHECK_FILE_DUPLICATE(icanal) 

!   ************************************************************************************************
!   **                                  CHECK FILE LIBRARY                                        **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                : CANAL ON WHICH OUTPUT DATA ARE WRITTEN                             **
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

    integer (kind=4) :: i, j;

    character (len=250) :: CHNAME_LIBRARY_DIRECTORY, CHNAME_TMP1, CHNAME_TMP2, CHNAME_TMP3;

    logical :: PROBE1, PROBE2, PROBE3, PROBE4;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the routine #############################################################
                                                                                         ! 
    CHTITLE = 'Check File (duplicate)';                                                  !
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
        if ( TRIM(DUPLICATE_CHFILE_FORMAT) == 'xyz' ) then;                              !
                                                                                         !
            CHOSEF1 = TRIM(CHOSEF1)//'/XYZ';                                             !
                                                                                         !
        else if ( TRIM(DUPLICATE_CHFILE_FORMAT) == 'lammps' ) then;                      !
                                                                                         !
            CHOSEF1 = TRIM(CHOSEF1)//'/LAMMPS';                                          !
                                                                                         !
        else if ( TRIM(DUPLICATE_CHFILE_FORMAT) == 'padua' ) then;                       !
                                                                                         !
            CHOSEF1 = TRIM(CHOSEF1)//'/Padua';                                           !
                                                                                         !
        else if ( TRIM(DUPLICATE_CHFILE_FORMAT) == 'lammps-gaff' ) then;                 !
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
!       stop; !//////////////////////////////////////////////////////////////////////////!
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
            if ( TRIM(DUPLICATE_CHFILE_FORMAT) == 'xyz' ) then;                          !
                                                                                         !
                CHNAME_LIBRARY_DIRECTORY = TRIM(CHNAME_LIBRARY_DIRECTORY)//'/XYZ/';      !
                                                                                         !
            else if ( TRIM(DUPLICATE_CHFILE_FORMAT) == 'lammps' ) then;                  !
                                                                                         !
                CHNAME_LIBRARY_DIRECTORY = TRIM(CHNAME_LIBRARY_DIRECTORY)//'/LAMMPS/';   !
                                                                                         !
            else if ( TRIM(DUPLICATE_CHFILE_FORMAT) == 'padua' ) then;                   !
                                                                                         !
                CHNAME_LIBRARY_DIRECTORY = TRIM(CHNAME_LIBRARY_DIRECTORY)//'/Padua/';    !
                                                                                         !
            else if ( TRIM(DUPLICATE_CHFILE_FORMAT) == 'lammps-gaff' ) then;             !
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
        IOSEF1 = LEN_TRIM(DUPLICATE_CHNAME_FILE);                                        !
                                                                                         !
        IOSEF2 = 67 - IOSEF1;                                                            !
                                                                                         !
        if ( IOSEF2 < 0 ) IOSEF2 = 1;                                                    !
                                                                                         !
        write(icanal,'(a70)') '| '//TRIM(DUPLICATE_CHNAME_FILE)// &                      !
                              REPEAT(' ',IOSEF2)//'|';                                   !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        CHNAME_TMP1 = TRIM(CHNAME_LIBRARY_DIRECTORY)// &                                 !
                      TRIM(DUPLICATE_CHNAME_FILE);                                       !
                                                                                         !
!       ### Set extension of the file for lammps templates #########################################
                                                                                         !
        if ( ( TRIM(DUPLICATE_CHFILE_FORMAT) == 'lammps'      ) .OR.   &                 !
             ( TRIM(DUPLICATE_CHFILE_FORMAT) == 'padua'       ) .OR.   &                 !
             ( TRIM(DUPLICATE_CHFILE_FORMAT) == 'lammps-gaff' ) ) then;                  !
                                                                                         !
            CHNAME_TMP1 = TRIM(CHNAME_TMP1)//'.template';                                !
                                                                                         !
        end if                                                                           !
                                                                                         !
!       ### Set extension of the file for the xyz format ###########################################
                                                                                         !
        if ( TRIM(DUPLICATE_CHFILE_FORMAT) == 'xyz' ) then;                              !
                                                                                         !
            CHNAME_TMP1 = TRIM(CHNAME_TMP1)//'.xyz';                                     !
                                                                                         !
        end if                                                                           !
                                                                                         !
        IOSEF1 = LEN_TRIM(CHNAME_TMP1);                                                  !
                                                                                         !
        IOSEF2 = 67 - IOSEF1;                                                            !
                                                                                         !
        if ( IOSEF2 < 0 ) IOSEF2 = 1;                                                    !
                                                                                         !
        write(icanal,'(a70)') '| '//TRIM(CHNAME_TMP1)//REPEAT(' ',IOSEF2)//'|';          !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
        inquire(FILE=TRIM(CHNAME_TMP1),EXIST=PROBE3);                                    !
                                                                                         !
        if ( PROBE3 .EQV. .TRUE. ) then;                                                 !
                                                                                         !
            write(icanal,'(a70)') '| The file has been found in the library'// &         !
                                  REPEAT(' ',29)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
            if ( ( TRIM(DUPLICATE_CHFILE_FORMAT) == 'lammps'      ) .OR.  &              !
                 ( TRIM(DUPLICATE_CHFILE_FORMAT) == 'padua'       ) .OR.  &              !
                 ( TRIM(DUPLICATE_CHFILE_FORMAT) == 'lammps-gaff' )  ) then;             !
                                                                                         !
                CHNAME_TMP2 = TRIM(CHNAME_LIBRARY_DIRECTORY)//            &              !
                              TRIM(DUPLICATE_CHNAME_FILE)//'.parameter';                 !
                                                                                         !
                inquire(FILE=TRIM(CHNAME_TMP2),EXIST=PROBE4);                            !
                                                                                         !
                IOSEF1 = LEN_TRIM(CHNAME_TMP2);                                          !
                                                                                         !
                IOSEF2 = 67 - IOSEF1;                                                    !
                                                                                         !
                if ( IOSEF2 < 0 ) IOSEF2 = 1;                                            !
                                                                                         !
                write(icanal,'(a70)') '| '//TRIM(CHNAME_TMP2)// &                        !
                                      REPEAT(' ',IOSEF2)//'|';                           !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
                if ( PROBE4 .EQV. .FALSE. ) then;                                        !
                                                                                         !
                    write(icanal,'(a70)') '| The file has not been found '//     &       !
                                          'in the library'//REPEAT(' ',25)//'|';         !
                                                                                         !
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
                    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                      !
                                                                                         !
                    write(icanal,*);                                                     !
                                                                                         !
                    write(icanal,'(a14)') 'End of program';                              !
                                                                                         !
                    close(icanal);                                                       ! 
                                                                                         !
                    stop; !//////////////////////////////////////////////////////////////!
                                                                                         !
                else                                                                     !
                                                                                         !
                    write(icanal,'(a70)') '| The file has been found '//         &       !
                                          'in the library'//REPEAT(' ',29)//'|';         !
                                                                                         !
                    DUPLICATE_CHNAME_FILE = TRIM(CHNAME_LIBRARY_DIRECTORY)// &           !
                                            TRIM(DUPLICATE_CHNAME_FILE);                 !
                                                                                         !
                end if                                                                   !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            else                                                                         !
                                                                                         !
                DUPLICATE_CHNAME_FILE = TRIM(CHNAME_LIBRARY_DIRECTORY)// &               !
                                        TRIM(DUPLICATE_CHNAME_FILE);                     !
                                                                                         !
            end if                                                                       !
                                                                                         !
        else                                                                             !
                                                                                         !
            write(icanal,'(a70)') '| The file has not been found '//     &               !
                                  'in the library'//REPEAT(' ',25)//'|';                 !
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
            stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end if                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine CHECK_FILE_DUPLICATE











