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

subroutine GENERATE_CONFIG_PREFIX(icanal,iconfig,CHCONFIG_PREFIX)

!   ************************************************************************************************
!   **                         Generate the configuration prefix                                  **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal : Canal on which output data are written                                            **
!   **                                                                                            **
!   ************************************************************************************************

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal, iconfig;

!   ************************************************************************************************
 
    character (len=150), intent(out) :: CHCONFIG_PREFIX;

!   ************************************************************************************************

    integer (kind=4) :: h, i, j;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the current routine #####################################################
                                                                                         !
    CHTITLE = 'Generate a configuration prefix';                                         !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Generate the prefix of the current molecular configuration #################################
                                                                                         !
    CHCONFIG_PREFIX = 'XXX';                                                             !
                                                                                         !
    write(CHOSEF1,'(i3.3)') iconfig;                                                     !
                                                                                         !
    CHCONFIG_PREFIX = 'GEN'//TRIM(CHOSEF1);                                              !
                                                                                         ! 
!   ### Write the prefix of the current molecular configuration ####################################
                                                                                         !
    IOSEF1 = 70 - 16 - 1 - LEN_TRIM(CHCONFIG_PREFIX);                                    !
                                                                                         !
    if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                        !
                                                                                         !
    write(icanal,'(a70)') '| The prefix is '//     &                                     !
                          TRIM(CHCONFIG_PREFIX)//  &                                     !
                          REPEAT(' ',IOSEF1)//'|';                                       !  
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine GENERATE_CONFIG_PREFIX











