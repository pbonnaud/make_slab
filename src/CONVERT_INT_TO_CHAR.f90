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

subroutine CONVERT_INT_TO_CHAR(icanal,LOCAL_INTEGER,NLOOP,LOCAL_CHARACTER)

!   ************************************************************************************************
!   **          SET MATRIX TO SWITCH FROM THE CARTESIAN AXIS TO TRICLINIC AXIS                    **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal : CANAL ON WHICH OUTPUT DATA ARE WRITTEN                                            **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

!   use module_osef;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    integer (kind=4), intent(in) :: LOCAL_INTEGER, NLOOP;

    character (len=250), intent(out) :: LOCAL_CHARACTER;

!   ************************************************************************************************

    integer (kind=4) :: i, j;

    integer (kind=4) :: IMAX_VAL;

    character (len=250) :: CHAR01;

!   ************************************************************************************************
                                                                                         !
!   ### Initialization of the maximum value ########################################################
                                                                                         !
    IMAX_VAL = 10;                                                                       !
                                                                                         !
    do i = 1, NLOOP-1;                                                                   !
                                                                                         !
        IMAX_VAL = IMAX_VAL * 10;                                                        !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   ### Convert the integer into a character string ################################################
                                                                                         !
    j = NLOOP;                                                                           !
                                                                                         !
    do i = 1, NLOOP;                                                                     !
                                                                                         !
        if ( j < 1000 ) write(CHAR01,'(i3)') j;                                          !
                                                                                         !
        if ( j <  100 ) write(CHAR01,'(i2)') j;                                          !
                                                                                         !
        if ( j <   10 ) write(CHAR01,'(i1)') j;                                          !
                                                                                         !
        CHAR01 = 'i'//TRIM(CHAR01);                                                      !
                                                                                         !
        if ( LOCAL_INTEGER < IMAX_VAL ) then;                                            !
                                                                                         !
            write(LOCAL_CHARACTER,'('//TRIM(CHAR01)//')') LOCAL_INTEGER;                 !
                                                                                         !
        end if                                                                           !
                                                                                         !
        IMAX_VAL = IMAX_VAL / 10;                                                        !
                                                                                         !
        j = j - 1;                                                                       !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   write(icanal,*) TRIM(LOCAL_CHARACTER);                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine CONVERT_INT_TO_CHAR
