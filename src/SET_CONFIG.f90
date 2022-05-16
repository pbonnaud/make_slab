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

subroutine SET_CONFIG(CHCONF_END_NAME,CELL_AXIS,CELL_ANGDEG,PASSA,PASSB)

!   ************************************************************************************************
!   ** 
!   ************************************************************************************************


    use module_data_in;


!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    character (len=50), intent(in) :: CHCONF_END_NAME;

    real (kind=8), dimension(1:3), intent(in) :: CELL_AXIS, CELL_ANGDEG;

    real (kind=8), dimension(1:3,1:3), intent(in) :: PASSA, PASSB;

!   ************************************************************************************************

    integer (kind=4) :: h, i;

!   ************************************************************************************************

    open(103,file=trim(CHCONF_END_NAME));

    if ( iboxc == 1 ) then
        write(103,*) NWORK+8;
    else
        write(103,*) NWORK;
    end if

    write(103,'(6f15.6)') CELL_AXIS(1:3), CELL_ANGDEG(1:3);

    do h = 1, 2
        do i = 1, NBFINAL(h)
            if ( NATFINAL(i,h) == 'Ow' .OR. NATFINAL(i,h) == 'Hw' ) write(103,'(a3,3f15.6)') NATFINAL(i,h), RFINAL(i,1:3,h);
        end do
    end do

    do h = 1, 2
        do i = 1, NBFINAL(h)
            if ( NATFINAL(i,h) == 'Cw' ) write(103,'(a3,3f15.6)') NATFINAL(i,h), RFINAL(i,1:3,h);
        end do
    end do

    do h = 1, 2
        do i = 1, NBFINAL(h)
            if ( NATFINAL(i,h) == 'Ca' ) write(103,'(a3,3f15.6)') NATFINAL(i,h), RFINAL(i,1:3,h);
        end do
    end do

    do h = 1, 2
        do i = 1, NBFINAL(h)
            if ( NATFINAL(i,h) == 'Os' ) write(103,'(a3,3f15.6)') NATFINAL(i,h), RFINAL(i,1:3,h);
        end do
    end do

    do h = 1, 2
        do i = 1, NBFINAL(h)
            if ( NATFINAL(i,h) == 'Si' ) write(103,'(a3,3f15.6)') NATFINAL(i,h), RFINAL(i,1:3,h);
        end do
    end do

    do h = 1, 2
        do i = 1, NBFINAL(h)
            if ( NATFINAL(i,h) == 'O' ) write(103,'(a3,3f15.6)') NATFINAL(i,h), RFINAL(i,1:3,h);
        end do
    end do

    if ( iboxc == 1 ) then
        write(103,'(a3,3f15.6)') 'Xe', MATMUL(PASSA(1:3,1:3),(/-0.5, -0.5, -0.5/));
        write(103,'(a3,3f15.6)') 'Xe', MATMUL(PASSA(1:3,1:3),(/-0.5, -0.5,  0.5/));
        write(103,'(a3,3f15.6)') 'Xe', MATMUL(PASSA(1:3,1:3),(/-0.5,  0.5, -0.5/));
        write(103,'(a3,3f15.6)') 'Xe', MATMUL(PASSA(1:3,1:3),(/-0.5,  0.5,  0.5/));
        write(103,'(a3,3f15.6)') 'Xe', MATMUL(PASSA(1:3,1:3),(/ 0.5, -0.5, -0.5/));
        write(103,'(a3,3f15.6)') 'Xe', MATMUL(PASSA(1:3,1:3),(/ 0.5, -0.5,  0.5/));
        write(103,'(a3,3f15.6)') 'Xe', MATMUL(PASSA(1:3,1:3),(/ 0.5,  0.5, -0.5/));
        write(103,'(a3,3f15.6)') 'Xe', MATMUL(PASSA(1:3,1:3),(/ 0.5,  0.5,  0.5/));
    end if

    close(103);

end subroutine SET_CONFIG




