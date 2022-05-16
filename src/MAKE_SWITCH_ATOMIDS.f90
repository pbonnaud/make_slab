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

subroutine MAKE_SWITCH_ATOMIDS(icanal) !,             &
!                             CHEXT,              & 
!                             MATA,               &
!                             FLAG_SORT_MOLEC); 

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

!   integer (kind=4) :: ILENGTH_CHEXT, ATOMS_FLAG;

!   integer (kind=4) :: ixy, ixz, iyz;

!   integer (kind=4) :: imcon, ICONNECT, IATOM_TYPE;

!   integer (kind=4) :: IFOUND;

!   real (kind=8) :: XLO, XHI, YLO, YHI, ZLO, ZHI;

!   character (len=150) :: CHARLINE;

!   integer (kind=4) :: EOF, EOF2;

    integer (kind=4) :: IOSEF1, IOSEF2, IOSEF3, IOSEF4, IOSEF5;

    real (kind=8) :: ROSEF1, ROSEF2, ROSEF3, ROSEF4, ROSEF5;

!   real (kind=8), dimension(1:3) :: TAB_ROSEF; 

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
    CHTITLE = 'Switch atom IDs';                                                         !
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
!   ### Write the number of atom IDs to switch #####################################################
                                                                                         !
    write(icanal,*) '| NSWITCH_ATOMIDS : ', NSWITCH_ATOMIDS, REPEAT(' ',1)//'|';

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|'; 

!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Switch atom IDs ############################################################################
                                                                                         !
    if ( NSWITCH_ATOMIDS > 0 ) then;

        do iswitch = 1, NSWITCH_ATOMIDS;                                                     !
                                                                                         !
            write(icanal,*) SWITCH_ATOMIDS(1,iswitch), ' --> ', SWITCH_ATOMIDS(2,iswitch);

            write(icanal,*);

            do i = 1, NATOM;

                if ( CONFIG_ATOM_TYPE(i) == SWITCH_ATOMIDS(1,iswitch) ) then;

                    CONFIG_ATOM_TYPE(i) = SWITCH_ATOMIDS(2,iswitch);

                end if

            end do

            NTYPE_ATOM = NTYPE_ATOM - 1;

        end do

        write(icanal,*) 'Atom types were switched';

        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Switch bond IDs ############################################################################
                                                                                         !
    if ( NSWITCH_BONDIDS > 0 ) then;

        do iswitch = 1, NSWITCH_BONDIDS;                                                     !
                                                                                         !
            write(icanal,*) SWITCH_BONDIDS(1,iswitch), ' --> ', &                            !
                            SWITCH_BONDIDS(2,iswitch);

            write(icanal,*);

            do i = 1, NATOM;

                if ( BOND_TYPE(i) == SWITCH_BONDIDS(1,iswitch) ) then;

                    BOND_TYPE(i) = SWITCH_BONDIDS(2,iswitch);

                end if

            end do

            NTYPE_BOND = NTYPE_BOND - 1;

        end do

        write(icanal,*) 'Bond types were switched';

        write(icanal,*);

    end if

!   stop; !//////////////////////////////////////////////////////////////////////////////!


!   ### Switch angle IDs ############################################################################
                                                                                         !
    if ( NSWITCH_ANGLEIDS > 0 ) then;                                                    !
                                                                                         !
        do iswitch = 1, NSWITCH_ANGLEIDS;                                                    !
                                                                                         !
            write(icanal,*) SWITCH_ANGLEIDS(1,iswitch), ' --> ', &                           !
                            SWITCH_ANGLEIDS(2,iswitch);

            write(icanal,*);

            do i = 1, NATOM;

                if ( ANGLE_TYPE(i) == SWITCH_ANGLEIDS(1,iswitch) ) then;

                    ANGLE_TYPE(i) = SWITCH_ANGLEIDS(2,iswitch);

                end if

            end do

            NTYPE_ANGLE = NTYPE_ANGLE - 1;

        end do

        write(icanal,*) 'Angle types were switched';

        write(icanal,*);

    end if

!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Switch dihedral IDs ########################################################################
                                                                                         !
    if ( NSWITCH_DIHEDRALIDS > 0 ) then;

        do iswitch = 1, NSWITCH_DIHEDRALIDS;                                                 !
                                                                                         !
            write(icanal,*) SWITCH_DIHEDRALIDS(1,iswitch), ' --> ', &                        !
                            SWITCH_DIHEDRALIDS(2,iswitch);                                   !
                                                                                         !
            write(icanal,*);

            do i = 1, NATOM;

                if ( DIHEDRAL_TYPE(i) == SWITCH_DIHEDRALIDS(1,iswitch) ) then;

                    DIHEDRAL_TYPE(i) = SWITCH_DIHEDRALIDS(2,iswitch);

                end if

            end do

            NTYPE_DIHEDRAL = NTYPE_DIHEDRAL - 1;

        end do

        write(icanal,*) 'Dihedral types were switched';

        write(icanal,*);

    end if 

!   stop; !//////////////////////////////////////////////////////////////////////////////!

!   ### WRAP UP THE ROUTINE PROCESS ################################################################

!   if ( my_id == root_process ) then;                                                   !
        call CPU_TIME(t_stop);                                                           !
                                                                                         !
        tps1 = ( t_stop - t_start ) / 3600.0d0;                                          !
        tps2 = ( tps1 - INT( tps1 ) ) * 60;                                              !
        tps3 = ( tps2 - INT( tps2 ) ) * 60;                                              !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        write(icanal,'(a19,i6,a3,i6,a5,i6,a25)') '| COMPUTING TIME : ',    &             !
                                                 INT( tps1 ),'[h]',        &             !
                                                 INT( tps2 ),'[min]',      &             !
                                                 INT( tps3 ),'[s]'//REPEAT(' ',21)//'|'; !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
        write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                  !

!   end if                                                                               !

!   stop;

end subroutine MAKE_SWITCH_ATOMIDS
