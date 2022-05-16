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

subroutine DEALLOCATE_INSERT_COMMAND(icanal)

!   ************************************************************************************************
!   **                            Seek for region keyword in the file                             **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                : CANAL ON WHICH OUTPUT DATA ARE WRITTEN                             **
!   ** PRIMITIVE_CELL_AXIS   : DIMENSIONS OF THE PRIMITIVE CELL                                   **
!   ** PRIMITIVE_CELL_ANGDEG : ANGLES OF THE PRIMITIVE CELL                                       **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

    use module_osef;

    use module_inserts;

    use module_slab_keywords;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

!   ************************************************************************************************

    integer (kind=4) :: h, i, j;

!   ************************************************************************************************

    character (len=250) :: CHARLINE;

!   ************************************************************************************************

    integer (kind=4) :: EOF;

!   ************************************************************************************************
                                                                                         !
!   ### Allocate arrays for the reading of the insert command ######################################
                                                                                         !
    if ( NLAMMPS_INSERTS > 0 ) then;                                                     !
                                                                                         !
        deallocate(LAMMPS_INSERT_NAME);                                                  !
                                                                                         !
        deallocate(LAMMPS_INSERT_FILE_OPTION);                                           !
                                                                                         !
        deallocate(LAMMPS_INSERT_CHFILE_FORMAT);                                         !
                                                                                         !
        deallocate(LAMMPS_INSERT_CHNAME_FILE);                                           !
                                                                                         !
        deallocate(LAMMPS_INSERT_METHOD);                                                !
                                                                                         !
        deallocate(LAMMPS_INSERT_COM_RG);                                                !
                                                                                         !
        deallocate(LAMMPS_INSERT_FLAG);                                                  !
                                                                                         !
        deallocate(LAMMPS_INSERT_REGION_NAME);                                           !
                                                                                         !
        deallocate(LAMMPS_INSERT_RANDOM_SEED);                                           !
                                                                                         !
        deallocate(LAMMPS_INSERT_REGION_IRANK);                                          !
                                                                                         !
        deallocate(LAMMPS_INSERT_NMOLECULE);                                             !
                                                                                         !
        deallocate(LAMMPS_INSERT_NTYPE_ATOM_MAX);                                        !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Deallocate arrays related to the insert_modify command #####################################
                                                                                         !
    if ( NLAMMPS_INSERTS > 0 ) then;                                                     !
                                                                                         !
        if ( LAMMPS_INSERT_ITETHER == 1 ) then;                                          !
                                                                                         !
            deallocate(LAMMPS_INSERT_MODIFY_STYLE);                                      !
                                                                                         !
            deallocate(LAMMPS_INSERT_TETHER_COEFFS);                                     !
                                                                                         !
        end if                                                                           !
                                                                                         !
        if ( LAMMPS_INSERT_IPOTENTIAL_CLASS2 == 1 ) then;                                !
                                                                                         !
            deallocate(LAMMPS_INSERT_POTENTIAL_CLASS2_CHTYPE);                           !
                                                                                         !
            deallocate(LAMMPS_INSERT_POTENTIAL_CLASS2_NAT);                              !
                                                                                         !
            deallocate(LAMMPS_INSERT_POTENTIAL_CLASS2);                                  !
                                                                                         !
        end if                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Write message for the success in reanding the input file ###################################
                                                                                         !
!   write(icanal,'(a70)') '| Insert commands were read '//REPEAT(' ',41)//'|';           !
                                                                                         !
!   write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine DEALLOCATE_INSERT_COMMAND

