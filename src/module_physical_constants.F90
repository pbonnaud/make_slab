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



module module_physical_constants

!   ************************************************************************************************

    implicit none;


!   ************************************************************************************************
                                                                                         !
    real (kind=8), parameter :: PI    = 3.14159265358979323846d0;                        !
                                                                                         !
    real (kind=8), parameter :: TWOPI = 2.0d0 * PI;                                      !
                                                                                         !
    real (kind=8), parameter :: RAD_TO_DEG = 360.0d0 / TWOPI;                            ! IN [DEG]
                                                                                         !
    real (kind=8), parameter :: DEG_TO_RAD = TWOPI / 360.0d0;                            ! IN [RAD]
                                                                                         !
    real (kind=8), parameter :: kB    = 1.3806503E-23;                                   ! BOLTZMANN'S CONSTANT IN [J/K]
                                                                                         !
    real (kind=8), parameter :: Na    = 6.02214179E23;                                   ! AVOGADRO'S CONSTANT IN [at/mol]
                                                                                         !
    real (kind=8), parameter :: KJMOL_TO_KCALMOL = 0.239006d0;                           ! CONVERSION CONSTANT IN [KCal/mol] 
                                                                                         !
!   ************************************************************************************************

end module module_physical_constants
