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

subroutine VECT_PRODUCT(VECTU,VECTV,VECTW)

!   ************************************************************************************************
!   **                                COMPUTE VECTOR PRODUCT                                      **
!   ************************************************************************************************

    implicit none

!   ************************************************************************************************

    real (kind=8), dimension(1:3), intent(in) :: VECTU, VECTV;

    real (kind=8), dimension(1:3), intent(out) :: VECTW;

!   ***********************************************************************************************

    VECTW(1) = VECTU(2) * VECTV(3) - VECTU(3) * VECTV(2);
    VECTW(2) = VECTU(3) * VECTV(1) - VECTU(1) * VECTV(3);
    VECTW(3) = VECTU(1) * VECTV(2) - VECTU(2) * VECTV(1);

end subroutine VECT_PRODUCT




