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

subroutine CONVERT_REAL_TO_CHAR(icanal,LOCAL_REAL,NLOOP,LOCAL_CHARACTER)

!   ************************************************************************************************
!   **          SET MATRIX TO SWITCH FROM THE CARTESIAN AXIS TO TRICLINIC AXIS                    **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal : CANAL ON WHICH OUTPUT DATA ARE WRITTEN                                            **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    integer (kind=4), intent(in) :: NLOOP;

    real (kind=8), intent(in) :: LOCAL_REAL;

!   ************************************************************************************************

    character (len=250), intent(out) :: LOCAL_CHARACTER;

!   ************************************************************************************************

    integer (kind=4) :: i, j;

    real (kind=8) :: RMAX_VAL;

    character (len=250) :: CHAR01;

!   ************************************************************************************************
                                                                                         !
!   ### Initialization of the maximum value ########################################################
                                                                                         !
    RMAX_VAL = 10.0d0;                                                                   !
                                                                                         !
    do i = 1, NLOOP-1;                                                                   !
                                                                                         !
        RMAX_VAL = RMAX_VAL * 10.0d0;                                                    !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   write(icanal,*) 'IMAX_VAL : ', IMAX_VAL;                                             !
                                                                                         !
!   ### Convert the integer into a character string ################################################
                                                                                         !
    j = NLOOP;                                                                           !
                                                                                         !
    do i = 1, NLOOP;                                                                     !
                                                                                         !
!       write(icanal,*) i, j, RMAX_VAL;                                                  !
                                                                                         !
        if ( j < 1000 ) write(CHAR01,'(i3)') j+3;                                        !
                                                                                         !
        if ( j <  100 ) write(CHAR01,'(i2)') j+3;                                        !
                                                                                         !
        if ( j <   10 ) write(CHAR01,'(i1)') j+3;                                        !
                                                                                         !
        CHAR01 = 'f'//TRIM(CHAR01)//'.2';                                                !
                                                                                         !
        if ( LOCAL_REAL < RMAX_VAL ) then;                                               !
                                                                                         !
!           write(icanal,*) TRIM(CHAR01), LOCAL_REAL;                                    !
                                                                                         !
            write(LOCAL_CHARACTER,'('//TRIM(CHAR01)//')') LOCAL_REAL;                    !
                                                                                         !
        end if                                                                           !
                                                                                         !
        RMAX_VAL = RMAX_VAL / 10.0d0;                                                    !
                                                                                         !
        j = j - 1;                                                                       !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   write(icanal,*) TRIM(LOCAL_CHARACTER), LOCAL_REAL;                                   !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine CONVERT_REAL_TO_CHAR
