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

subroutine CHECK_INSERT_FILE_LIBRARY(icanal,ILOCAL_INSERTS) 

!   ************************************************************************************************
!   **                                  Check Insert File Library                                 **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                : CANAL ON WHICH OUTPUT DATA ARE WRITTEN                             **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

    use module_library;

    use module_osef;

    use module_inserts;

!   ************************************************************************************************

    implicit none;

!   ***********************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    integer (kind=4), intent(in) :: ILOCAL_INSERTS;

!   ************************************************************************************************

    integer (kind=4) :: i, j;

    character (len=250) :: CHNAME_LIBRARY_DIRECTORY;

    logical :: PROBE1, PROBE2, PROBE3, PROBE4;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the routine #############################################################
                                                                                         ! 
    CHTITLE = 'Check Insert File Library';                                               !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### build the full path to library files #######################################################
                                                                                         !
    if ( NPATH_TO_LIBRARY > 0 ) then;                                                    !
                                                                                         !
        CHNAME_LIBRARY_DIRECTORY = 'XXX';                                                !
                                                                                         !
        do i = 1, NPATH_TO_LIBRARY;                                                      !
                                                                                         !
            if ( TRIM(LAMMPS_INSERT_LIB_NAME(ILOCAL_INSERTS)) == &                       !
                 TRIM(LIB_PATH_NAME(i)) ) then;                                          !
                                                                                         !
                CHNAME_LIBRARY_DIRECTORY = TRIM(PATH_TO_LIBRARY(i));                     !
                                                                                         !
                EXIT;                                                                    !
                                                                                         !
            end if                                                                       !
                                                                                         !
        end do                                                                           !
                                                                                         !
!       write(icanal,*) TRIM(CHNAME_LIBRARY_DIRECTORY);                                  !
                                                                                         !
!       write(icanal,*);                                                                 !
                                                                                         !
    else                                                                                 !
                                                                                         !
        write(icanal,*) 'Not implemented! - Stop';                                       !
                                                                                         !
        stop; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                                         !
    end if                                                                               !          
                                                                                         !
!   stop; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                                         !
!   ### Check the presence of the library directory ################################################
                                                                                         !
!   CHNAME_LIBRARY_DIRECTORY = 'LIB-MOLECULAR-MODELS';                                   !
                                                                                         !
    call FIND_DIRECTORY(icanal,1,CHNAME_LIBRARY_DIRECTORY);                              !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialization of file extensions ##########################################################
                                                                                         !
    NFILE_EXT = 0;                                                                       !
                                                                                         !
    LAMMPS_INSERT_CHFILE_EXT(1:3) = 'XXX';                                               !
                                                                                         !
!   ### Define the sub-directory to consider in the main library directory and file extensions #####
                                                                                         !
    select case(TRIM(LAMMPS_INSERT_CHFILE_FORMAT(ILOCAL_INSERTS)));                      !
                                                                                         !
        case('xyz');                                                                     !
                                                                                         !
!           CHOSEF1 = '/XYZ';                                                            !
                                                                                         !
            NFILE_EXT = 1;                                                               !
                                                                                         !
            LAMMPS_INSERT_CHFILE_EXT(1) = '.xyz';                                        !
                                                                                         !
        case('lammps');                                                                  !
                                                                                         !
!           CHOSEF1 = '/LAMMPS';                                                         !
                                                                                         !
            NFILE_EXT = 3;                                                               !
                                                                                         !
            LAMMPS_INSERT_CHFILE_EXT(1) = '.template';                                   !
                                                                                         !
            LAMMPS_INSERT_CHFILE_EXT(2) = '.parameter';                                  !
                                                                                         !
            LAMMPS_INSERT_CHFILE_EXT(3) = '.info';                                       !
                                                                                         !
        case('padua');                                                                   !
                                                                                         !
!           CHOSEF1 = '/Padua';                                                          !
                                                                                         !
            NFILE_EXT = 3;                                                               !
                                                                                         !
            LAMMPS_INSERT_CHFILE_EXT(1) = '.template';                                   !
                                                                                         !
            LAMMPS_INSERT_CHFILE_EXT(2) = '.parameter';                                  !
                                                                                         !                                                                                         !
            LAMMPS_INSERT_CHFILE_EXT(3) = '.info';                                       !
                                                                                         !
        case('lammps-gaff');                                                             !
                                                                                         !
!           CHOSEF1 = '/lammps-gaff';                                                    !
                                                                                         !
            NFILE_EXT = 3;                                                               !
                                                                                         !
            LAMMPS_INSERT_CHFILE_EXT(1) = '.template';                                   !
                                                                                         !
            LAMMPS_INSERT_CHFILE_EXT(2) = '.parameter';                                  !
                                                                                         !
            LAMMPS_INSERT_CHFILE_EXT(3) = '.info';                                       !
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
!   CHNAME_LIBRARY_DIRECTORY = TRIM(CHNAME_LIBRARY_DIRECTORY)//TRIM(CHOSEF1);            !
                                                                                         !
!   write(icanal,*) TRIM(CHNAME_LIBRARY_DIRECTORY);                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Check that the sub-directory is in the library main directory ##############################
                                                                                         !
    call FIND_DIRECTORY(icanal,1,CHNAME_LIBRARY_DIRECTORY);                              !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the name of the file that need to be found in the library sub-directory ############## 
                                                                                         !
    IOSEF1 = 70 - 2 - 1 - LEN_TRIM(LAMMPS_INSERT_CHNAME_FILE(ILOCAL_INSERTS));           !
                                                                                         !
    if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                        !
                                                                                         !
    write(icanal,'(a70)') '| '//TRIM(LAMMPS_INSERT_CHNAME_FILE(ILOCAL_INSERTS))// &      !
                          REPEAT(' ',IOSEF1)//'|';                                       !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Build file names with their full path and extension ########################################
                                                                                         !
    LAMMPS_INSERT_CHNAME_FILE(ILOCAL_INSERTS) =       &                                  !
    TRIM(CHNAME_LIBRARY_DIRECTORY)//'/'//             &                                  !
    TRIM(LAMMPS_INSERT_CHNAME_FILE(ILOCAL_INSERTS));                                     !
                                                                                         !
!   ### Initialization of the variable related to the info file for the lammps input file ##########
                                                                                         !
    ifile_info = 0;                                                                      !
                                                                                         !
!   ### Loop over the maximum number of files that could be read by the program ####################
                                                                                         ! 
    do i = 1, NFILE_EXT;                                                                 !
                                                                                         !
!       ### Build the file with its extension ######################################################
                                                                                         !
        CHOSEF6 = TRIM(LAMMPS_INSERT_CHNAME_FILE(ILOCAL_INSERTS))// &                    !
                  TRIM(LAMMPS_INSERT_CHFILE_EXT(i));                                     !
                                                                                         !
!       ### Write names of files with paths and extensions #########################################
                                                                                         !
!       if ( TRIM(LAMMPS_INSERT_CHNAME_FILE(ILOCAL_INSERTS)) == 'XXX' ) CYCLE;           !
                                                                                         !
        IOSEF1 = 70 - 2 - 1 - LEN_TRIM(CHOSEF6);                                         !
                                                                                         !
        if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                    !
                                                                                         !
        write(icanal,'(a70)') '| '//TRIM(CHOSEF6)//REPEAT(' ',IOSEF1)//'|';              !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       ### Check that the file is in the sub-directory ############################################
                                                                                         !
        inquire(FILE=TRIM(CHOSEF6),EXIST=PROBE3);                                        !
                                                                                         !
!       ### Write the output message of the previous check #########################################
                                                                                         !
        CHOSEF1 = 'has';                                                                 !
                                                                                         !
        if ( PROBE3 .EQV. .FALSE. ) CHOSEF1 = TRIM(CHOSEF1)//' not';                     !
                                                                                         !
        if ( ( PROBE3 .EQV. .TRUE. ) .AND. ( i == 3 ) ) ifile_info = 1;                  !
                                                                                         !
        CHOSEF1 = '| The file '//TRIM(CHOSEF1)//' been found in the library';            !
                                                                                         !
        IOSEF1 = 70 - 1 - LEN_TRIM(CHOSEF1);                                             !
                                                                                         !
        if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                    !
                                                                                         !
        write(icanal,'(a70)') TRIM(CHOSEF1)//REPEAT(' ',IOSEF1)//'|';                    !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end do                                                                               !
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
end subroutine CHECK_INSERT_FILE_LIBRARY











