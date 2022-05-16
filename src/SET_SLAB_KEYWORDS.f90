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

subroutine SET_SLAB_KEYWORDS(icanal)

!   ************************************************************************************************
!   **                                    READ XYZ FILES                                          **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                                          **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

    use module_slab_keywords;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the current routine #####################################################
                                                                                         !
    CHTITLE = 'Set slab keywords';                                                       !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set the number of keywords considered in the current version of the program ################
                                                                                         !
    NSLAB_INPUT_KEYWORDS = 7;                                                            !
                                                                                         !
!   ### Set slab keywords ##########################################################################
                                                                                         !
!   +++ Dimension 1 --> 20                                                               !
                                                                                         !
    SLAB_INPUT_KEYWORDS(1) = 'region';                                                   !
                                                                                         !
    SLAB_INPUT_KEYWORDS(2) = 'create_box';                                               !
                                                                                         !
    SLAB_INPUT_KEYWORDS(3) = 'insert';                                                   !
                                                                                         !
    SLAB_INPUT_KEYWORDS(4) = 'duplicate';                                                !
                                                                                         !
    SLAB_INPUT_KEYWORDS(5) = 'insert_modify';                                            !
                                                                                         !
    SLAB_INPUT_KEYWORDS(6) = 'lmp_input';                                                !
                                                                                         !  
    SLAB_INPUT_KEYWORDS(7) = 'set_lib';                                                  !
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
end subroutine SET_SLAB_KEYWORDS
