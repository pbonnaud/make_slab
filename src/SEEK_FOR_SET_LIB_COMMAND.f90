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

subroutine SEEK_FOR_SET_LIB_COMMAND(icanal)

!   ************************************************************************************************
!   **                            Seek for the set_lib keyword in the file                        **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                : CANAL ON WHICH OUTPUT DATA ARE WRITTEN                             **
!   **                                                                                            **
!   ************************************************************************************************

    use module_osef;

    use module_library;

    use module_slab_keywords;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

!   ************************************************************************************************

    integer (kind=4) :: h, i, j;

!   ************************************************************************************************

    integer (kind=4) :: iblank, icharact;

!   ************************************************************************************************

    character (len=250) :: CHARLINE;

!   ************************************************************************************************

    integer (kind=4) :: EOF;

!   ************************************************************************************************
                                                                                         !
!   ### Initialization of arrays related to regions ################################################
                                                                                         !
    NPATH_TO_LIBRARY = 0;                                                                !
                                                                                         !
!   ### Open the input file ########################################################################
                                                                                         !
    open(1,file='input_100.dat',status='old');                                           !
                                                                                         !
    EOF = 0;                                                                             !
                                                                                         !
    do                                                                                   !
                                                                                         !
        read(1,'(a)',iostat=EOF) CHARLINE;                                               !
                                                                                         !
        if ( EOF /= 0 ) EXIT;                                                            !
                                                                                         !
        if ( INDEX(CHARLINE,'#') > 0 ) CYCLE;                                            !
                                                                                         !
        IOSEF2 = INDEX(CHARLINE,TRIM(SLAB_INPUT_KEYWORDS(7)));                           !
                                                                                         !
        if ( IOSEF2 > 20 ) CYCLE;                                                        !
                                                                                         !
        if ( IOSEF2 > 0 ) then;                                                          !
                                                                                         !
            NPATH_TO_LIBRARY = NPATH_TO_LIBRARY + 1;                                     !
                                                                                         !
        end if                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
    close(1);                                                                            !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the number of regions found in the input file ########################################
                                                                                         !
    write(icanal,'(a36,i4,a30)') '| Number of defined library paths : ', &               !
                                 NPATH_TO_LIBRARY,                       &               !
                                 REPEAT(' ',29)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                                         !
!   ### Allocate arrays containing region properties ###############################################
                                                                                         !
    allocate(LIB_PATH_NAME(1:NPATH_TO_LIBRARY));                                         !
                                                                                         !
    allocate(PATH_TO_LIBRARY(1:NPATH_TO_LIBRARY));                                       !
                                                                                         !
!   ### Initialization of arrays containing path-to-library properties #############################
                                                                                         !
    LIB_PATH_NAME(1:NPATH_TO_LIBRARY) = 'XXX';                                           !
                                                                                         !
    PATH_TO_LIBRARY(1:NPATH_TO_LIBRARY) = 'XXX';                                         !
                                                                                         !
!   stop; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                                         !
!   ### Read path-to-library properties ############################################################
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !
    open(1,file='input_100.dat',status='old');                                           !
                                                                                         !
    EOF = 0;                                                                             !
                                                                                         !
    do                                                                                   !
                                                                                         !
        read(1,'(a)',iostat=EOF) CHARLINE;                                               !
                                                                                         !
        if ( EOF /= 0 ) EXIT;                                                            !
                                                                                         !
        if ( INDEX(CHARLINE,'#') > 0 ) CYCLE;                                            !
                                                                                         !
        IOSEF2 = INDEX(CHARLINE,TRIM(SLAB_INPUT_KEYWORDS(7)));                           !
                                                                                         !
        if ( IOSEF2 > 20 ) CYCLE;                                                        !
                                                                                         !
        if ( IOSEF2 > 0 ) then;                                                          !
                                                                                         !
!           write(icanal,*) TRIM(CHARLINE), IOSEF2;                                      !
                                                                                         !
            IOSEF1 = IOSEF1 + 1;                                                         !
                                                                                         !
!           ### Read first the name of the path to the library #####################################
                                                                                         !
            read(CHARLINE,*) CHOSEF1, LIB_PATH_NAME(IOSEF1);                             !
                                                                                         !
!           ### Get the position and the length of the name in the character string ################
                                                                                         !
            IOSEF3 = INDEX(CHARLINE,TRIM(LIB_PATH_NAME(IOSEF1)));                        !
                                                                                         !
            IOSEF4 = LEN_TRIM(LIB_PATH_NAME(IOSEF1));                                    !
                                                                                         !
!           write(icanal,*) TRIM(LIB_PATH_NAME(IOSEF1)), IOSEF3, IOSEF4;                 !
                                                                                         !
!           ### Read the path to the library in the character string ###############################
                                                                                         !
            IOSEF5 = IOSEF3 + IOSEF4;                                                    !
                                                                                         !
            IOSEF6 = LEN_TRIM(CHARLINE);                                                 !
                                                                                         !
            iblank = 0;                                                                  !
                                                                                         !
            icharact = 0;                                                                !
                                                                                         !
            j = IOSEF5;                                                                  !
                                                                                         !
            do i = IOSEF5, IOSEF6;                                                       !
                                                                                         !
!               write(icanal,*) i, CHARLINE(i:i);                                        !
                                                                                         !
                if ( CHARLINE(i:i) == ' ' ) then;                                        !
                                                                                         !
                    if ( icharact == 1 ) then;                                           !
                                                                                         !
                        j = i;                                                           !
                                                                                         !
                        EXIT;                                                            ! 
                                                                                         !
                    else if ( icharact == 0 ) then;                                      !
                                                                                         !
                        iblank = i;                                                      !
                                                                                         !
                        CYCLE;                                                           !
                                                                                         !
                    end if                                                               !
                                                                                         !
                else                                                                     !
                                                                                         !
                    icharact = 1;                                                        !
                                                                                         !
                end if                                                                   !
                                                                                         !
!               write(icanal,*) 'iblank / icharact : ', iblank, icharact;                !
                                                                                         !
                if ( i == IOSEF6 ) j = IOSEF6;                                           !
                                                                                         !
            end do                                                                       !
                                                                                         !
            IOSEF5 = iblank + 1;                                                         !
                                                                                         !
            PATH_TO_LIBRARY(IOSEF1) = TRIM(CHARLINE(IOSEF5:j));                          !
                                                                                         !
!           write(icanal,*) j, TRIM(PATH_TO_LIBRARY(IOSEF1));                            !
                                                                                         !
        end if                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
    close(1);                                                                            !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the list of library paths found in the input file ####################################
                                                                                         !
    write(icanal,'(a70)') '| List of library paths '//REPEAT('>',45)//'|';               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    do i = 1, NPATH_TO_LIBRARY;                                                          ! 
                                                                                         !
        IOSEF1 = 70 - 5 - 17 - 1 -             &                                         !
                 LEN_TRIM(LIB_PATH_NAME(i)) -  &                                         !
                 LEN_TRIM(PATH_TO_LIBRARY(i));                                           !
                                                                                         !
        if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                    !
                                                                                         !
        write(icanal,'(a70)') '|    '//                  &                               !
                              TRIM(LIB_PATH_NAME(i))//   &                               !
                              ' having the path '//      &                               !
                              TRIM(PATH_TO_LIBRARY(i))// &                               !
                              REPEAT(' ',IOSEF1)//'|';                                   !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   ### Write message for the success in reanding the input file ###################################
                                                                                         !
!   write(icanal,'(a70)') '| Regions have been read '//REPEAT(' ',44)//'|';              !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine SEEK_FOR_SET_LIB_COMMAND











