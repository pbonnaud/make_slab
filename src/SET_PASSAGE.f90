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

subroutine SET_PASSAGE(icanal,LOCAL_CELL_AXIS,LOCAL_CELL_ANGDEG,LOCAL_PASSA,LOCAL_PASSB)

!   ************************************************************************************************
!   **          SET MATRIX TO SWITCH FROM THE CARTESIAN AXIS TO TRICLINIC AXIS                    **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal : CANAL ON WHICH OUTPUT DATA ARE WRITTEN                                            **
!   ** CELL_AXIS : DIMENSIONS OF THE SIMULATION BOX                                               **
!   ** CELL_ANGDEG : ANGLES OF THE SIMULATION BOX                                                 **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

    use module_physical_constants;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    real (kind=8), dimension(1:3), intent(in) :: LOCAL_CELL_AXIS, LOCAL_CELL_ANGDEG;

    real (kind=8), dimension(1:3,1:3), intent(out) :: LOCAL_PASSA, LOCAL_PASSB;

!   ************************************************************************************************

    real (kind=8), dimension(1:3) :: CELL_ANGRAD, CELL_WORK;

    real (kind=8) :: ISING, ITANG;

    real (kind=8), dimension(1:3) :: CELL_COSANG, CELL_SINANG, CELL_TANANG;

    real (kind=8) :: C1, C2, D1, D2, LOCAL_VOL;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
    CHTITLE = 'Set passage matrices';                                                    !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Convert angles from degrees to radians #####################################################
                                                                                         !
    CELL_ANGRAD(1:3) = LOCAL_CELL_ANGDEG(1:3) * twopi / 360.0d0;                         !
                                                                                         !
    write(icanal,'(a70)') '| Set up transfer matrix : '//REPEAT(' ',42)//'|';            !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a8,3f15.6,a17)') '| ANG : ', CELL_ANGRAD(1:3), REPEAT(' ',16)//'|';   !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Compute sine, cosine, and tangent of cell angles ###########################################
                                                                                         !
    CELL_COSANG(1:3) = DCOS( CELL_ANGRAD(1:3) );                                         !
                                                                                         !
    CELL_SINANG(1:3) = DSIN( CELL_ANGRAD(1:3) );                                         !
                                                                                         !  
    CELL_TANANG(1:3) = DTAN( CELL_ANGRAD(1:3) );                                         !
                                                                                         !
    if ( CELL_ANGRAD(3) ==  90.0d0 ) CELL_TANANG(3) =  1000000000000.0d0;                !
                                                                                         !
    if ( CELL_ANGRAD(3) == -90.0d0 ) CELL_TANANG(3) = -1000000000000.0d0;                !
                                                                                         !
!   ### Write sine, cosine, and tangent of cell angles #############################################
                                                                                         !
    write(icanal,'(a8,3f15.6,a17)') '| COS : ', CELL_COSANG(1:3), REPEAT(' ',16)//'|';   !
                                                                                         !
    write(icanal,'(a8,3f15.6,a17)') '| SIN : ', CELL_SINANG(1:3), REPEAT(' ',16)//'|';   !
                                                                                         !
    write(icanal,'(a8,3f15.6,a17)') '| TAN : ', CELL_TANANG(1:3), REPEAT(' ',16)//'|';   !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Compte constant C1, C2, and the volume of the cell #########################################
                                                                                         !
    C1  = 2.0d0 * CELL_COSANG(1) * CELL_COSANG(2) * CELL_COSANG(3);                      !
                                                                                         !
    C2  = DSQRT( 1.0d0 - CELL_COSANG(1)**2 - CELL_COSANG(2)**2 - CELL_COSANG(3)**2 + C1 );
                                                                                         !
    LOCAL_VOL = LOCAL_CELL_AXIS(1) * LOCAL_CELL_AXIS(2) * LOCAL_CELL_AXIS(3) * C2;       !
                                                                                         !
!   ### Write the volume of the cell ###############################################################
                                                                                         !
    write(icanal,'(a13,f20.6,a37)') '| VOL [A3] : ', LOCAL_VOL, REPEAT(' ',36)//'|';     !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialization of passage matrices #########################################################
                                                                                         !
    LOCAL_PASSA(1:3,1:3) = 0.0d0;                                                        !
                                                                                         !
    LOCAL_PASSB(1:3,1:3) = 0.0d0;                                                        !
                                                                                         !
!   ### Compute the passage matrix A ###############################################################
                                                                                         !
    LOCAL_PASSA(1,1) = LOCAL_CELL_AXIS(1);                                               !
                                                                                         !
    LOCAL_PASSA(2,2) = LOCAL_CELL_AXIS(2) * CELL_SINANG(3);                              !
                                                                                         !
    LOCAL_PASSA(3,3) = LOCAL_VOL / ( LOCAL_CELL_AXIS(1) * LOCAL_CELL_AXIS(2) * CELL_SINANG(3) );
                                                                                         !
    LOCAL_PASSA(2,3) = LOCAL_CELL_AXIS(3) * ( CELL_COSANG(1) - CELL_COSANG(2) * CELL_COSANG(3) ) / CELL_SINANG(3);
                                                                                         !
    LOCAL_PASSA(1,3) = LOCAL_CELL_AXIS(3) * CELL_COSANG(2);                              !
                                                                                         !
    LOCAL_PASSA(1,2) = LOCAL_CELL_AXIS(2) * CELL_COSANG(3);                              !
                                                                                         !
!   ### Write the passage matrix A #################################################################
                                                                                         !
    write(icanal,'(a70)') '| RNO --> RON : '//REPEAT(' ',53)//'|';                       ! 
                                                                                         !
    write(icanal,'(a2,3f15.6,a23)') '| ', LOCAL_PASSA(1,1:3), REPEAT(' ',22)//'|';       !
                                                                                         !
    write(icanal,'(a2,3f15.6,a23)') '| ', LOCAL_PASSA(2,1:3), REPEAT(' ',22)//'|';       !
                                                                                         !
    write(icanal,'(a2,3f15.6,a23)') '| ', LOCAL_PASSA(3,1:3), REPEAT(' ',22)//'|';       !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Compute the passage matrix B ###############################################################

    LOCAL_PASSB(1,1) =   1.0d0 / LOCAL_CELL_AXIS(1);

    LOCAL_PASSB(2,2) =   1.0d0 / ( LOCAL_CELL_AXIS(2) * CELL_SINANG(3) );

    LOCAL_PASSB(3,3) =   LOCAL_CELL_AXIS(1) * LOCAL_CELL_AXIS(2) * CELL_SINANG(3) / LOCAL_VOL;

    LOCAL_PASSB(2,3) = - LOCAL_CELL_AXIS(1) * LOCAL_CELL_AXIS(3) *            &
                       ( CELL_COSANG(1) - CELL_COSANG(2) * CELL_COSANG(3) ) / &
                       ( LOCAL_VOL * CELL_SINANG(3) );

    LOCAL_PASSB(1,3) = LOCAL_CELL_AXIS(2) * LOCAL_CELL_AXIS(3) *                               &
                       ( CELL_COSANG(1) / CELL_TANANG(3) - CELL_COSANG(2) / CELL_SINANG(3) ) / &
                       LOCAL_VOL;

    LOCAL_PASSB(1,2) = - 1.0d0 / ( LOCAL_CELL_AXIS(1) * CELL_TANANG(3) );

!   ### Write the passage matrix B #################################################################

    write(icanal,'(a70)') '| RON --> RNO : '//REPEAT(' ',53)//'|';

    write(icanal,'(a2,3f15.6,a23)') '| ', LOCAL_PASSB(1,1:3), REPEAT(' ',22)//'|';

    write(icanal,'(a2,3f15.6,a23)') '| ', LOCAL_PASSB(2,1:3), REPEAT(' ',22)//'|';

    write(icanal,'(a2,3f15.6,a23)') '| ', LOCAL_PASSB(3,1:3), REPEAT(' ',22)//'|';

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!

!   ### Define the working cell ####################################################################

    CELL_WORK(1:3)  = (/LOCAL_PASSA(1,1),LOCAL_PASSA(2,2),LOCAL_PASSA(3,3)/);            !
                                                                                         !
    write(icanal,'(a14,3f15.6,a11)') '| CELL WORK : ',    &                              !
                                     CELL_WORK(1:3),      &                              !
                                     REPEAT(' ',10)//'|';                                !
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine SET_PASSAGE
