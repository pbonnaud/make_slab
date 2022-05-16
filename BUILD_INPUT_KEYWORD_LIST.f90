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

subroutine BUILD_INPUT_KEYWORD_LIST(icanal)

!   ************************************************************************************************
!   **                            Seek for region keyword in the file                             **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                : CANAL ON WHICH OUTPUT DATA ARE WRITTEN                             **
!   ** PRIMITIVE_CELL_AXIS   : DIMENSIONS OF THE PRIMITIVE CELL                                   **
!   ** PRIMITIVE_CELL_ANGDEG : ANGLES OF THE PRIMITIVE CELL                                       **
!   **                                                                                            **
!   ************************************************************************************************

    use module_osef;

    use module_slab_keywords;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

!   ************************************************************************************************

    integer (kind=4) :: h, i, j;

!   ************************************************************************************************

    character (len=250) :: CHARLINE;

!   ************************************************************************************************

    integer (kind=4) :: EOF;

!   ************************************************************************************************
                                                                                         !
!   ### Initialization of arrays for the list of read input keywords ###############################
                                                                                         !
    NREAD_SLAB_INPUT_KEYWORDS = 0;                                                       !
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
        if ( INDEX(CHARLINE,'classic') > 0 ) EXIT;                                       !
                                                                                         !
        do i = 1, NSLAB_INPUT_KEYWORDS;                                                  !
                                                                                         !
            IOSEF2 = INDEX(CHARLINE,TRIM(SLAB_INPUT_KEYWORDS(i))//' ');                  !
                                                                                         !
            if ( IOSEF2 > 20 ) CYCLE;                                                    !
                                                                                         !
            if ( IOSEF2 == 0 ) CYCLE;                                                    !
                                                                                         !
            NREAD_SLAB_INPUT_KEYWORDS = NREAD_SLAB_INPUT_KEYWORDS + 1;                   !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
    close(1);                                                                            !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the number of input keywords found in the input file #################################
                                                                                         !
    write(icanal,'(a29,i4,a37)') '| Number of input keywords : ', &                      !
                                 NREAD_SLAB_INPUT_KEYWORDS,       &                      !
                                 REPEAT(' ',36)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### Allocate arrays containing the list of input keywords ######################################
                                                                                         !
    allocate(LIST_SLAB_INPUT_KEYWORDS(1:NREAD_SLAB_INPUT_KEYWORDS));                     !
                                                                                         !
!   ### Initialization of arrays containing region properties ######################################
                                                                                         !
    LIST_SLAB_INPUT_KEYWORDS(1:NREAD_SLAB_INPUT_KEYWORDS) = 'XXX';                       !
                                                                                         !
!   ### Read region properties #####################################################################
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
        if ( INDEX(CHARLINE,'classic') > 0 ) EXIT;                                       !
                                                                                         !
        do i = 1, NSLAB_INPUT_KEYWORDS;                                                  !
                                                                                         !
            IOSEF2 = INDEX(CHARLINE,TRIM(SLAB_INPUT_KEYWORDS(i))//' ');                  !
                                                                                         !
            if ( IOSEF2 > 20 ) CYCLE;                                                    !
                                                                                         !
            if ( IOSEF2 == 0 ) CYCLE;                                                    !
                                                                                         !
            IOSEF1 = IOSEF1 + 1;                                                         !
                                                                                         !
            LIST_SLAB_INPUT_KEYWORDS(IOSEF1) = TRIM(SLAB_INPUT_KEYWORDS(i));             !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
    close(1);                                                                            !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the list of input keywords found in the input file ###################################
                                                                                         !
    write(icanal,'(a70)') '| List of input keywords '//REPEAT('>',44)//'|';              !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    do i = 1, NREAD_SLAB_INPUT_KEYWORDS;                                                 ! 
                                                                                         !
        IOSEF1 = 70 - 5 - 1 - LEN_TRIM(LIST_SLAB_INPUT_KEYWORDS(i));                     !
                                                                                         !
        if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                    !
                                                                                         !
        write(icanal,'(a70)') '|    '//                            &                     !
                              TRIM(LIST_SLAB_INPUT_KEYWORDS(i))//  &                     !
                              REPEAT(' ',IOSEF1)//'|';                                   !
                                                                                         !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### Write message for the success in reanding the input file ###################################
                                                                                         !
    write(icanal,'(a70)') '| The list of keywords has been built '//REPEAT(' ',31)//'|'; !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine BUILD_INPUT_KEYWORD_LIST











