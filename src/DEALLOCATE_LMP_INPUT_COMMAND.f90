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

subroutine DEALLOCATE_LMP_INPUT_COMMAND(icanal)

!   ************************************************************************************************
!   **                            Seek for region keyword in the file                             **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                : CANAL ON WHICH OUTPUT DATA ARE WRITTEN                             **
!   ** PRIMITIVE_CELL_AXIS   : DIMENSIONS OF THE PRIMITIVE CELL                                   **
!   ** PRIMITIVE_CELL_ANGDEG : ANGLES OF THE PRIMITIVE CELL                                       **
!   **                                                                                            **
!   ************************************************************************************************

    use module_lmp_input;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

!   ************************************************************************************************
                                                                                         !
!   ### Deallocate arrays related to fix properties ################################################
                                                                                         !
    if ( NLMP_INPUT_FIX > 0 ) then;                                                      !
                                                                                         !
        deallocate(LMP_INPUT_RUN_STYLE);                                                 !
                                                                                         !
!       deallocate(LMP_INPUT_RUN_TEMPK_INIT);                                            !
                                                                                         !
!       deallocate(LMP_INPUT_RUN_TEMPK_FINAL);                                           !
                                                                                         !
        deallocate(LMP_INPUT_RUN_TAUNH);                                                 !
                                                                                         !
!       deallocate(LMP_INPUT_RUN_NSTEPS);                                                !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Deallocate arrays related to frozen entities ###############################################
                                                                                         !
                                                                                         !
    if ( NLMP_INPUT_FROZEN > 0 ) then;                                                   !
                                                                                         !
        deallocate(LMP_INPUT_FROZEN_GRPNAME);                                            !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Deallocate arrays related to run properties ################################################
                                                                                         !
    if ( NLMP_INPUT_RUN > 0 ) then;                                                      !
                                                                                         !
        deallocate(LMP_INPUT_RUN_TEMPK_INIT);                                            !
                                                                                         !
        deallocate(LMP_INPUT_RUN_TEMPK_FINAL);                                           !
                                                                                         !
        deallocate(LMP_INPUT_RUN_NSTEPS);                                                !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Deallocate arrays related to group properties ##############################################
                                                                                         !
    if ( NLMP_INPUT_GROUP > 0 ) then;                                                    !
                                                                                         !
        deallocate(LMP_INPUT_GROUP_NAME);                                                !
                                                                                         !
        deallocate(LMP_INPUT_GROUP_NLIST);                                               !
                                                                                         !
        deallocate(LMP_INPUT_GROUP_LIST);                                                !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Write message for the success in reanding the input file ###################################
                                                                                         !
!   write(icanal,'(a70)') '| Insert commands were read '//REPEAT(' ',41)//'|';           !
                                                                                         !
!   write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine DEALLOCATE_LMP_INPUT_COMMAND

