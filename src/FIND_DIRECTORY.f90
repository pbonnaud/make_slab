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

subroutine FIND_DIRECTORY(icanal,ILOCAL_NLEVELS,CHLOCAL_DIRECTORY) 

!   ************************************************************************************************
!   **                                       Find directory                                       **
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

    integer (kind=4), intent(in) :: ILOCAL_NLEVELS;

!   ************************************************************************************************

    character (len=250),intent(inout) :: CHLOCAL_DIRECTORY;

!   ************************************************************************************************

    integer (kind=4) :: i;

!   ************************************************************************************************

    logical :: PROBE1, PROBE2;

!   ************************************************************************************************
                                                                                         !
!   ### Write the initial name of the directory ####################################################
                                                                                         !
    IOSEF1 = 70 - 33 - 1 - LEN_TRIM(CHLOCAL_DIRECTORY);                                  !
                                                                                         !
    if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                        !
                                                                                         !
    write(icanal,'(a70)') '| Looking for the directory name '// &                        !
                          TRIM(CHLOCAL_DIRECTORY)//             &                        !
                          REPEAT(' ',IOSEF1)//'|';                                       !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the number of levels over which the directory will be looked for #####################
                                                                                         !
    write(icanal,'(a40,i4,a26)') '| The directory will be looked for over ', &           !
                                 ILOCAL_NLEVELS,                             &           !
                                 ' levels'//REPEAT(' ',18)//'|';                         !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    IOSEF10 = ILOCAL_NLEVELS + 1;                                                        !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Check the presence of the library ##########################################################
                                                                                         !
    do i = 1, IOSEF10;                                                                   ! Loop over the parent directories
                                                                                         !
        open(22,file='Tempo');                                                           ! Create a temporary script to check the presence of the directory
                                                                                         !
        write(22,'(a11)') '#!/bin/bash';                                                 !
                                                                                         !
        write(22,*);                                                                     !
                                                                                         !
        write(22,*) 'LIB_PATH="'//TRIM(CHLOCAL_DIRECTORY)//'"';                          !
                                                                                         !
        write(22,*);                                                                     !
                                                                                         !
!       CHOSEF1 = 'if [ -d "'//TRIM(CHLOCAL_DIRECTORY)//'" ]; then touch EXIST; fi';
        CHOSEF1 = 'if [ -d $LIB_PATH ]; then touch EXIST; fi';                           !
                                                                                         !
!       write(icanal,*) TRIM(CHOSEF1);        

!       write(22,*) 'if [ -d "'//TRIM(CHLOCAL_DIRECTORY)//'" ]; then touch EXIST; fi';   !
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
!       write(icanal,*) 'PROBE2 : ', PROBE2;                                             !
                                                                                         !
!       write(icanal,*);                                                                 !
                                                                                         !
!       stop; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                                         !
        call system('rm -f EXIST Tempo');                                                !
                                                                                         !
        if ( PROBE2 .EQV. .TRUE. ) EXIT;                                                 !
                                                                                         !
        CHLOCAL_DIRECTORY = '../'//TRIM(CHLOCAL_DIRECTORY);                              ! Update the path to the directory
                                                                                         !
    end do                                                                               !
                                                                                         !
!   ### Write the location of the directory ########################################################
                                                                                         !
    write(icanal,'(a70)') '| The library directory was found at the '// &                !
                          'following path : '//REPEAT(' ',11)//'|';                      !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    IOSEF1 = LEN_TRIM(CHLOCAL_DIRECTORY);                                                !
                                                                                         !
    IOSEF2 = 67 - IOSEF1;                                                                !
                                                                                         !
    if ( IOSEF2 < 0 ) IOSEF2 = 1;                                                        !
                                                                                         !
    write(icanal,'(a70)') '| '//TRIM(CHLOCAL_DIRECTORY)//REPEAT(' ',IOSEF2)//'|';        !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine FIND_DIRECTORY











