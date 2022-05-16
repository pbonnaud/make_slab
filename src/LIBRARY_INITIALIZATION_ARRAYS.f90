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

subroutine LIBRARY_INITIALIZATION_ARRAYS(icanal,IFILE)

!   ************************************************************************************************
!   **                                READ LAMMPS TEMPLATE FILES                                  **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                           **
!   **                                                                                            **
!   ************************************************************************************************

    use module_size_arrays;

    use module_library;

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ***********************************************************************************************

    integer (kind=4), intent(in) :: icanal, IFILE;

!   ************************************************************************************************
                                                                                         !
!   ### Initialization of library arrays ###########################################################
                                                                                         !
    NATOM_LIBRARY(IFILE)     = 0;                                                        !
                                                                                         !
    NBOND_LIBRARY(IFILE)     = 0;                                                        !
                                                                                         !
    NANGLE_LIBRARY(IFILE)    = 0;                                                        !
                                                                                         !
    NDIHEDRAL_LIBRARY(IFILE) = 0;                                                        !
                                                                                         !
    NIMPROPER_LIBRARY(IFILE) = 0;                                                        !
                                                                                         !
    NTYPE_ATOM_LIBRARY(IFILE)     = 0;                                                   !
                                                                                         !
    NTYPE_BOND_LIBRARY(IFILE)     = 0;                                                   !
                                                                                         !
    NTYPE_ANGLE_LIBRARY(IFILE)    = 0;                                                   !
                                                                                         !
    NTYPE_DIHEDRAL_LIBRARY(IFILE) = 0;                                                   !
                                                                                         !
    NTYPE_IMPROPER_LIBRARY(IFILE) = 0;                                                   !

    if ( IMAX_NATOM > 0 ) then;                                                          !
                                                                                         ! 
        CONFIG_ATOMID_LIBRARY(1:IMAX_NATOM,IFILE)     = 0;                               !
                                                                                         !
        CONFIG_ATOM_TYPE_LIBRARY(1:IMAX_NATOM,IFILE)  = 0;                               !
                                                                                         !
        CONFIG_QI_LIBRARY(1:IMAX_NATOM,IFILE)         = 0.0d0;                           !
                                                                                         !
        CONFIG_VI_LIBRARY(1:3,1:IMAX_NATOM,IFILE)     = 0.0d0;                           !
                                                                                         !
        CONFIG_RI_LIBRARY(1:3,1:IMAX_NATOM,IFILE)     = 0.0d0;                           !
                                                                                         !
        CONFIG_NAT_LIBRARY(1:IMAX_NATOM,IFILE)        = 'XXX';                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( IMAX_NBOND > 0 ) then;                                                          !
                                                                                         ! 
        BOND_TYPE_LIBRARY(1:IMAX_NBOND,IFILE)         = 0;                               !
                                                                                         !
        BOND_ATOMID_LIBRARY(1:2,1:IMAX_NBOND,IFILE)   = 0;                               !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( IMAX_NANGLE > 0 ) then;                                                         !
                                                                                         !
        ANGLE_TYPE_LIBRARY(1:IMAX_NANGLE,IFILE)       = 0;                               !
                                                                                         !
        ANGLE_ATOMID_LIBRARY(1:3,1:IMAX_NANGLE,IFILE) = 0;                               !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( IMAX_NDIHEDRAL > 0 ) then;                                                      !
                                                                                         !
        DIHEDRAL_TYPE_LIBRARY(1:IMAX_NDIHEDRAL,IFILE)       = 0;                         !
                                                                                         !
        DIHEDRAL_ATOMID_LIBRARY(1:4,1:IMAX_NDIHEDRAL,IFILE) = 0;                         !
                                                                                         ! 
    end if                                                                               !
                                                                                         !
    if ( IMAX_NIMPROPER > 0 ) then;                                                      !
                                                                                         !
        IMPROPER_TYPE_LIBRARY(1:IMAX_NIMPROPER,IFILE)       = 0;                         !
                                                                                         !
        IMPROPER_ATOMID_LIBRARY(1:4,1:IMAX_NIMPROPER,IFILE) = 0;                         !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Write an output message for the success of the initialization process ######################
                                                                                         !
    write(icanal,'(a70)') '| Library arrays were initialized'//REPEAT(' ',36)//'|';      !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine LIBRARY_INITIALIZATION_ARRAYS
