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

subroutine FIND_ATOM_ELEMENTID(icanal,MASSE_ATOMI,ELEMENTI)

!   ************************************************************************************************
!   **                                   Find the atom element ID                                 **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                                          **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    real (kind=8), intent(in) :: MASSE_ATOMI;

!   ************************************************************************************************

    integer (kind=4), intent(out) :: ELEMENTI;

!   ************************************************************************************************

    integer (kind=4) :: j;

!   real (kind=8) :: ROSEF1;

!   ************************************************************************************************

!   integer (kind=4) :: ILENGTH_TITLE;

!   character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Find the location in the Mendeleiev periodic table of elements #############################
                                                                                         !
    ELEMENTI = 0;                                                                        !
                                                                                         !
    do j = 1, 200;                                                                       !
                                                                                         !
        ROSEF1 = ABS( MASSE_ATOMI - MOLAR_MASS_ELEMENT(j) );                             !
                                                                                         !
        if ( ROSEF1 < 0.1d0 ) then;                                                      !
                                                                                         !
            ELEMENTI = j;                                                                !
                                                                                         !
            EXIT;                                                                        !
                                                                                         !
        end if                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
    if ( ELEMENTI == 0 ) then;                                                           !
                                                                                         !
        write(icanal,'(a70)') '| The atom was not found in the periodic '// &            !
                              'table of elements - stop'//REPEAT(' ',4)//'|';            !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                  !
                                                                                         !
        close(icanal);                                                                   !
                                                                                         !
        stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
    end if                                                                               !
                                                                                         !
end subroutine FIND_ATOM_ELEMENTID
