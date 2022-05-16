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

subroutine MAKE_TRANSLATION_MOLECULEIDS(icanal)

!   ************************************************************************************************
!   **                                 Make switch atom IDs                                       **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal        : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                                   **
!   **                                                                                            **
!   ** CHEXT         : NAME OF THE LAMMPS CONFIGURATION                                           **
!   **                                                                                            **
!   ** MATA              : PASSAGE MATRIX                                                         **
!   **                                                                                            **
!   ** FLAG_SORT_MOLEC   : FLAG TO CHECK IF MOLECULES AND ATOMS NEED TO BE SORTED                 ** 
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

    use module_config;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

!   include 'mpif.h'

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

!   ************************************************************************************************

    integer (kind=4) :: ierr;

!   ************************************************************************************************

    integer (kind=4) :: i, j;

    integer (kind=4) :: iswitch;

    integer (kind=4) :: IOSEF1, IOSEF2, IOSEF3, IOSEF4, IOSEF5;

    real (kind=8) :: ROSEF1, ROSEF2, ROSEF3, ROSEF4, ROSEF5;

    character (len=10) :: CHOSEF1, CHOSEF2, CHOSEF3;

    character (len=250) :: CHAIN_LENGTH;

!   ************************************************************************************************

    real (kind=8) :: t_start, t_stop;

    real (kind=8) :: tps1, tps2, tps3;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the routine #############################################################
                                                                                         !
    CHTITLE = 'Spatial translation of molecule IDs';                                     !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         ! 
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Get the time at the beginning of the routine ###############################################
                                                                                         !
    call CPU_TIME(t_start);                                                              !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the number of molecules to translate #################################################
                                                                                         !
    write(icanal,'(a29,i4,a37)') '| NTRANSLATION_MOLECULEIDS : ', &                      !
                                 NTRANSLATION_MOLECULEIDS,        &                      !
                                 REPEAT(' ',36)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Translate the selected molecules ###########################################################
                                                                                         !
    if ( NTRANSLATION_MOLECULEIDS > 0 ) then;                                            !
                                                                                         !
        do j = 1, NTRANSLATION_MOLECULEIDS;                                              !
                                                                                         !
            do i = 1, NATOM;                                                             ! 
                                                                                         !
                if ( CONFIG_MOLECULEID(i) == TRANSLATION_MOLECULEIDS(j) ) then;          !
                                                                                         !
                    CONFIG_RI(1:3,i) = CONFIG_RI(1:3,i) + TRANSLATION_VECTOR(1:3,j);     !
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### WRAP UP THE ROUTINE PROCESS ################################################################

!   if ( my_id == root_process ) then;                                                   !
        call CPU_TIME(t_stop);                                                           !
                                                                                         !
        tps1 = ( t_stop - t_start ) / 3600.0d0;                                          !
                                                                                         !
        tps2 = ( tps1 - INT( tps1 ) ) * 60;                                              !
                                                                                         !
        tps3 = ( tps2 - INT( tps2 ) ) * 60;                                              !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        write(icanal,'(a19,i6,a3,i6,a5,i6,a25)') '| Computing time : ',    &             !
                                                 INT( tps1 ),'[h]',        &             !
                                                 INT( tps2 ),'[min]',      &             !
                                                 INT( tps3 ),'[s]'//REPEAT(' ',21)//'|'; !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                  !

!   end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine MAKE_TRANSLATION_MOLECULEIDS
