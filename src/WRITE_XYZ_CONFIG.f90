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

subroutine WRITE_XYZ_CONFIG(icanal,       &
                            CHEXT,        &
                            CELL_AXIS,    &
                            CELL_ANGDEG) 

!   ************************************************************************************************
!   **                                 WRITE XYZ FILES                                            **
!   ************************************************************************************************
!   **                                                                                            **
!   ** INPUTS:                                                                                    **
!   ** ------                                                                                     **
!   **                                                                                            **
!   ** icanal : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                                          **
!   ** CHEXT  : NAME OF THE LAMMPS CONFIGURATION                                                  **
!   **                                                                                            **
!   ** NATOM       : NUMBER OF ATOMS                                                              **
!   ** CONFIG_NAT  : NATURE OF ATOMS IN THE SIMULATION BOX                                        **
!   ** CONFIG_RIJ  : POSITION OF ATOMS IN THE SIMULATION BOX                                      **
!   ** CELL_AXIS   : CELL DIMENSIONS                                                              **
!   ** CELL_ANGDEG : CELL ANGLES                                                                  **
!   **                                                                                            **
!   ************************************************************************************************

    use module_size_arrays;

    use module_config;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    character (len=150), intent(in) :: CHEXT;

    real (kind=8), dimension(1:3), intent(in) :: CELL_AXIS, CELL_ANGDEG;

!   ************************************************************************************************

    integer (kind=4) :: i;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
    CHTITLE = 'WRITE XYZ CONFIG';                                                        !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
    open(103,file=TRIM(CHEXT)//'.xyz');                                                  !
    write(103,'(i8)') NATOM;                                                             !
    write(103,'(6f15.6)') CELL_AXIS(1:3), CELL_ANGDEG(1:3);                              !
                                                                                         !
    do i = 1, NATOM;                                                                     !
        write(103,'(a3,3f15.6)') CONFIG_NAT(i), CONFIG_RI(1:3,i);                        !
    end do                                                                               !
                                                                                         !
    close(103);                                                                          !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop;

end subroutine WRITE_XYZ_CONFIG

