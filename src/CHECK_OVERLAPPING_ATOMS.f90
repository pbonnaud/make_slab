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

subroutine CHECK_OVERLAPPING_ATOMS(icanal,       &
                                   ISTART,       &
                                   ISTOP,        &
                                   PASSA,        &
                                   PASSB,        &
                                   RCUTOFF,      &
                                   IOVERLAPPING)

!   ************************************************************************************************
!   **               THIS ROUTINE IS CHECKING IF THERE IS OVERLAPPING ATOMS                       **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                : CANAL ON WHICH OUTPUT DATA ARE WRITTEN                             **
!   ** ISTART : 
!   ** NATOM   : NUMBER OF ATOMS IN THE MOLECULAR CONFIGURATION
!   ** CONFIG_RI : COORDINATES OF ATOMS IN THE MOLECULAR CONFIGURATION 
!   ** PASSA
!   ** PASSB
!   ** RCUTOFF
!   ** IOVERLAPPING : FLAG TO TELL THE PROGRAM IF ATOMS ARE OVERLAPPING 
!   **                                                                                            **
!   ************************************************************************************************

    use module_config;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    integer (kind=4), intent(in) :: ISTART, ISTOP;
 
    real (kind=8), intent(in) :: RCUTOFF;

    real (kind=8), dimension(1:3,1:3), intent(in) :: PASSA, PASSB;

!   ************************************************************************************************

    integer (kind=4), intent(out) :: IOVERLAPPING;

!   ************************************************************************************************

    integer (kind=4) :: i, j;

    real (kind=8) :: RIJ;

    real (kind=8), dimension(1:3) :: DRIJ;

!   ************************************************************************************************
                                                                                         !
    IOVERLAPPING = 0;                                                                    !
                                                                                         !
    do i = ISTART, ISTOP;                                                                !
                                                                                         !
        do j = 1, ISTART-1;                                                              !
                                                                                         !
            DRIJ(1:3) = CONFIG_RI(1:3,j) - CONFIG_RI(1:3,i);                             !
                                                                                         !
            call APPLY_PBC(DRIJ(1:3),DRIJ(1:3),PASSA(1:3,1:3),PASSB(1:3,1:3));           !
                                                                                         !
            RIJ = DSQRT( DOT_PRODUCT(DRIJ(1:3),DRIJ(1:3)) );                             !
                                                                                         !
            if ( RIJ <= RCUTOFF ) then;                                                  !
                                                                                         !
                IOVERLAPPING = IOVERLAPPING + 1;                                         !
                                                                                         !
                EXIT;                                                                    !
                                                                                         !
            end if                                                                       !
                                                                                         !
        end do                                                                           !
                                                                                         !
        if ( IOVERLAPPING > 0 ) EXIT;                                                    !
                                                                                         !
    end do                                                                               !
                                                                                         !
end subroutine CHECK_OVERLAPPING_ATOMS











