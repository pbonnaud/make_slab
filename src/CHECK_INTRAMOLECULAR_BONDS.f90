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

subroutine CHECK_INTRAMOLECULAR_BONDS(icanal,PASSA,PASSB) 

!   ************************************************************************************************
!   **                                  READ THE INPUT FILE                                       **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal : Canal on which output data are written                                            **
!   **                                                                                            **
!   ** PASSA  :                                                                                   **
!   **                                                                                            **
!   ** PASSB  :                                                                                   **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

    use module_physical_constants;

    use module_size_arrays;

    use module_library;

    use module_config;

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    real (kind=8), dimension(1:3,1:3), intent(in) :: PASSA, PASSB;

!   ************************************************************************************************

    integer (kind=4) :: i, j, k, m, p;

    integer (kind=4) :: LOCAL_NNEIGHBOR;

    real (kind=8) :: RIJ01, RIJ02, DCH_MIN;

    integer (kind=4), dimension(1:20) :: LOCAL_NEIGHBOR_ATOMID;

    real (kind=8), dimension(1:3) :: DRIJ;

    real (kind=8), dimension(1:20) :: LOCAL_NEIGHBOR_RIJ; 

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the routine #############################################################
                                                                                         !
    CHTITLE = 'Check intramolecular bonds';                                              !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Check C-H bonds ############################################################################
                                                                                         !
    do i = 1, NATOM;                                                                     !
                                                                                         !
        if ( TRIM(CONFIG_NAT(i)) /= 'C' ) CYCLE;                                         !

        write(icanal,*) 'XXXXXXX', TRIM(CONFIG_NAT(i)), i;
   
        write(icanal,*);
                                                                                         !
        LOCAL_NNEIGHBOR = 0;                                                             ! 
                                                                                         !
        LOCAL_NEIGHBOR_ATOMID(1:20) = 0;                                                 !
                                                                                         !
        LOCAL_NEIGHBOR_RIJ(1:20) = 100000000.0d0;                                        !
                                                                                         !
        do j = 1, NATOM;                                                                 !
                                                                                         !
            if ( TRIM(CONFIG_NAT(j)) /= 'H' ) CYCLE;                                     !
                                                                                         !
            DRIJ(1:3) = CONFIG_RI(1:3,j) - CONFIG_RI(1:3,i);                             !
                                                                                         !
            call APPLY_PBC(DRIJ(1:3),DRIJ(1:3),PASSA(1:3,1:3),PASSB(1:3,1:3));           !
                                                                                         !
            RIJ02 = DOT_PRODUCT(DRIJ(1:3),DRIJ(1:3));                                    !
                                                                                         !
            if ( RIJ02 > 4.0d0 ) CYCLE;                                                  !
                                                                                         !
            RIJ01 = DSQRT( RIJ02 );                                                      ! 
                                                                                         !
            if ( LOCAL_NNEIGHBOR < 20 ) then;                                            !
                                                                                         !
                LOCAL_NNEIGHBOR = LOCAL_NNEIGHBOR + 1;                                   !
                                                                                         !
                LOCAL_NEIGHBOR_ATOMID(LOCAL_NNEIGHBOR) = j;                              !                
                                                                                         !
                LOCAL_NEIGHBOR_RIJ(LOCAL_NNEIGHBOR) = RIJ01;                             !
                                                                                         !
            end if                                                                       !
                                                                                         !
            if ( LOCAL_NNEIGHBOR == 20 ) then;                                           !
                                                                                         !
                ROSEF1 = MAXVAL( LOCAL_NEIGHBOR_RIJ(1:20) );                             !                
                                                                                         !
                if ( RIJ01 > ROSEF1 ) CYCLE;                                             !
                                                                                         !
                IOSEF1 = 0;                                                              !
                                                                                         !
                do k = 1, LOCAL_NNEIGHBOR;                                               !

                    if ( LOCAL_NEIGHBOR_RIJ(k) /= ROSEF1 ) CYCLE;

                    IOSEF1 = k;

                    EXIT;

                end do
                                                                                         !
                LOCAL_NEIGHBOR_ATOMID(IOSEF1) = j;                                       !
                                                                                         !
                LOCAL_NEIGHBOR_RIJ(IOSEF1) = RIJ01;                                      !
                                                                                         !
            end if                                                                       !
                                                                                         !
        end do                                                                           !
                                                                                         !
        if ( LOCAL_NNEIGHBOR == 0 ) CYCLE;

!       do j = 1, LOCAL_NNEIGHBOR;                                                       !
                                                                                         !
!           ROSEF1 = MINVAL( LOCAL_NEIGHBOR_RIJ(j+1:LOCAL_NNEIGHBOR) );                  !                
                                                                                         !
!           if ( LOCAL_NEIGHBOR_RIJ(j) < ROSEF1 ) CYCLE;                                 !
                                                                                         !
!           IOSEF1 = 0;                                                                  !

!           do k = 1, LOCAL_NNEIGHBOR;

!               if ( LOCAL_NEIGHBOR_RIJ(k) /= ROSEF1 ) CYCLE;

!               IOSEF1 = k;

!               EXIT;

!           end do
                                                                                         !
!           LOCAL_NEIGHBOR_RIJ(IOSEF1) = LOCAL_NEIGHBOR_RIJ(j);                          !
                                                                                         !
!           LOCAL_NEIGHBOR_ATOMID(IOSEF1) = j;

!           LOCAL_NEIGHBOR_RIJ(j) = ROSEF1;                                              !
                                                                                         !
!           LOCAL_NEIGHBOR_ATOMID(j) = IOSEF1;                                           !

!       end do                                                                           !
                                                                                         !
!       DCH_MIN = MINVAL(LOCAL_NEIGHBOR_RIJ(1:LOCAL_NNEIGHBOR));                         !
                                                                                         !
        do j = 1, LOCAL_NNEIGHBOR;

           write(icanal,*) j, LOCAL_NEIGHBOR_RIJ(j), LOCAL_NEIGHBOR_ATOMID(j), CONFIG_NAT(LOCAL_NEIGHBOR_ATOMID(j));

           if ( LOCAL_NEIGHBOR_RIJ(j) > 1.9d0 ) CYCLE;

           ROSEF1 = ABS( LOCAL_NEIGHBOR_RIJ(j) - 1.55019d0 );

           if ( ROSEF1 > 0.1d0 ) CYCLE;

           IOSEF1 = LOCAL_NEIGHBOR_ATOMID(j);                         

           DRIJ(1:3) = CONFIG_RI(1:3,IOSEF1) - CONFIG_RI(1:3,i);                         !
                                                                                         !
           call APPLY_PBC(DRIJ(1:3),DRIJ(1:3),PASSA(1:3,1:3),PASSB(1:3,1:3));            !
                                                                                         !
           RIJ02 = DOT_PRODUCT(DRIJ(1:3),DRIJ(1:3));                                     !
                                                                                         !
           RIJ01 = DSQRT( RIJ02 );                                                       ! 
                                                                                         !
           CONFIG_RI(1:3,IOSEF1) = CONFIG_RI(1:3,i) + DRIJ(1:3) / RIJ01;                 !
                                                                                         !
        end do

        write(icanal,*);

!       stop; !//////////////////////////////////////////////////////////////////////////!

    end do


!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine CHECK_INTRAMOLECULAR_BONDS
