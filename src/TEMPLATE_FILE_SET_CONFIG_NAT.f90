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



subroutine TEMPLATE_FILE_SET_CONFIG_NAT(icanal,IFILE)

!   ************************************************************************************************
!   **                                READ LAMMPS TEMPLATE FILES                                  **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                         : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                  **
!   **                                                                                            **
!   ** CHEXT                          : NAME OF THE LAMMPS CONFIGURATION                          **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

    use module_size_arrays;

    use module_library;

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ***********************************************************************************************

    integer (kind=4), intent(in) :: icanal, IFILE;

!   ************************************************************************************************

    character (len=250) :: CHEXT;

!   ************************************************************************************************

    integer (kind=4) :: i, j;

    integer (kind=4) :: EOF, EOF2;

    character (len=250) :: CHARLINE, CHARLINE2;

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the routine #############################################################
                                                                                         !
    CHTITLE = 'Set atom labels to atoms of the template file';                           !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialization of nature of atoms in the configuration of the library file #################
                                                                                         !
    CONFIG_NAT_LIBRARY(1:NATOM_LIBRARY(IFILE),IFILE) = 'XXX';                            !
                                                                                         !
!   ### Set nature of atoms in the configuration of the library file ###############################
                                                                                         !
    do j = 1, NATOM_LIBRARY(IFILE);                                                      !
                                                                                         !
        CONFIG_NAT_LIBRARY(j,IFILE) = &                                                  !
        ATOM_LABEL_LIBRARY(CONFIG_ATOM_TYPE_LIBRARY(j,IFILE),IFILE);                     !
                                                                                         !
!       write(icanal,*) j, CONFIG_ATOM_TYPE_LIBRARY(j,IFILE), &                          !
!                       CONFIG_NAT_LIBRARY(j,IFILE);                                     !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine TEMPLATE_FILE_SET_CONFIG_NAT
