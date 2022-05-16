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

subroutine DUPLICATE_ASSIGN_FORCE_FIELD(icanal)

!   ************************************************************************************************
!   **                               DUPLICATE THE ORIGINAL UNIT CELL                             **
!   ************************************************************************************************

    use module_data_in;

    use module_library;

    use module_duplicate;

    use module_config;

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

!   ************************************************************************************************

    integer (kind=4) :: i, j, k, m;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the routine title ####################################################################
                                                                                         !
    CHTITLE = 'Assign force field (duplicate)';                                          !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Assign pair potential parameter style  #####################################################
                                                                                         !
    IPOTENTIAL_CLASS2 = IPOTENTIAL_CLASS2_LIBRARY(1);                                    !
                                                                                         !
    if ( IPOTENTIAL_CLASS2 > 0 ) then;                                                   !
                                                                                         !
        CH_PAIR_STYLE = POTENTIAL_CLASS2_CHTYPE_LIBRARY(1);                              !
                                                                                         !
    end if                                                                               !
                                                                                         !
    write(icanal,'(a70)') '| Pair potential parameter style was assigned'// &            ! 
                          REPEAT(' ',24)//'|';                                           !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### Allocate arrays for pair potentials ########################################################
                                                                                         !
    if ( NTYPE_ATOM > 0 ) then;                                                          !
                                                                                         !
        allocate(POTENTIAL_CLASS2(1:2,1:NTYPE_ATOM));                                    !
                                                                                         ! 
        write(icanal,'(a70)') '| Arrays for pair potentials were allocated'// &          !
                              REPEAT(' ',26)//'|';                                       !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Update arrays for pair potentials ##########################################################
                                                                                         !
    if ( NTYPE_ATOM > 0 ) then;                                                          !
                                                                                         !
        POTENTIAL_CLASS2(1:2,1:NTYPE_ATOM) = &                                           !
        POTENTIAL_CLASS2_LIBRARY(1:2,1:NTYPE_ATOM,1);                                    !
                                                                                         !
        write(icanal,'(a70)') '| Arrays for pair potentials were updated'// &            !
                              REPEAT(' ',28)//'|';                                       !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Assign bond parameter properties ###########################################################
                                                                                         !
    NTYPE_BOND = 0;                                                                      !
                                                                                         ! 
    if ( NTYPE_BOND_LIBRARY(1) > 0 ) then;                                               !
                                                                                         !
        CH_BOND_STYLE = CH_BOND_STYLE_LIBRARY(1);                                        !
                                                                                         !
        NPARAM_BONDS = NPARAM_BONDS_LIBRARY(1);                                          ! 
                                                                                         !
        NTYPE_BOND = NTYPE_BOND_LIBRARY(1);                                              !
                                                                                         !
        write(icanal,'(a70)') '| Bond parameter properties were assigned'// &            !
                              REPEAT(' ',28)//'|';                                       !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Allocate arrays for bond potentials ########################################################
                                                                                         !
    if ( NTYPE_BOND > 0 ) then;                                                          !
                                                                                         !
        allocate(BOND_COEFFS(1:4,1:NTYPE_BOND));                                         !
                                                                                         !
        write(icanal,'(a70)') '| Arrays for bond potentials were allocated'// &          !
                              REPEAT(' ',26)//'|';                                       !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Update arrays for bond potentials ##########################################################
                                                                                         !
    if ( NTYPE_BOND > 0 ) then;                                                          !
                                                                                         !
        BOND_COEFFS(1:NPARAM_BONDS,1:NTYPE_BOND) =          &                            !
        BOND_COEFFS_LIBRARY(1:NPARAM_BONDS,1:NTYPE_BOND,1);                              !
                                                                                         !
        write(icanal,'(a70)') '| Arrays for bond potentials were updated'// &            !
                              REPEAT(' ',28)//'|';                                       !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Assign angle parameter properties ##########################################################
                                                                                         !
    NTYPE_ANGLE = 0;                                                                     !
                                                                                         ! 
    if ( NTYPE_ANGLE_LIBRARY(1) > 0 ) then;                                              !
                                                                                         !
        CH_ANGLE_STYLE = CH_ANGLE_STYLE_LIBRARY(1);                                      !
                                                                                         !
        NPARAM_ANGLES = NPARAM_ANGLES_LIBRARY(1);                                        ! 
                                                                                         !
        NTYPE_ANGLE = NTYPE_ANGLE_LIBRARY(1);                                            !
                                                                                         !
        write(icanal,'(a70)') '| Angle parameter properties were assigned'// &           !
                              REPEAT(' ',27)//'|';                                       !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Allocate arrays for angle potentials #######################################################
                                                                                         !
    if ( NTYPE_ANGLE > 0 ) then;                                                         !
                                                                                         !
        allocate(ANGLE_COEFFS(1:4,1:NTYPE_ANGLE));                                       !
                                                                                         !
        write(icanal,'(a70)') '| Arrays for angle potentials were allocated'// &         !
                              REPEAT(' ',25)//'|';                                       !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Update arrays for angle potentials #########################################################
                                                                                         !
    if ( NTYPE_ANGLE > 0 ) then;                                                         !
                                                                                         !
        ANGLE_COEFFS(1:NPARAM_ANGLES,1:NTYPE_ANGLE) =  &                                 !
        ANGLE_COEFFS_LIBRARY(1:NPARAM_ANGLES,1:NTYPE_ANGLE,1);                           !
                                                                                         !
        write(icanal,'(a70)') '| Arrays for angle potentials were updated'// &           !
                              REPEAT(' ',27)//'|';                                       !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### To be coded ################################################################################

    NTYPE_DIHEDRAL = 0;

    NTYPE_IMPROPER = 0;

!   ### Assign partial atomic charges ##############################################################
                                                                                         !
    do i = 1, NATOM;                                                                     !
                                                                                         !
        IOSEF1 = CONFIG_ATOM_TYPE(i);                                                    !
                                                                                         !
        ROSEF1 = ATOM_MASSE(IOSEF1);                                                     !
                                                                                         !
        IOSEF2 = 0;                                                                      !
                                                                                         !
        ROSEF2 = 0.0d0;                                                                  !
                                                                                         !
        ROSEF4 = 0.0d0;                                                                  !
                                                                                         !
!       write(icanal,*) '--> ', TRIM(CONFIG_NAT(i)), IOSEF1, ROSEF1, '>>>>>>>>>>>>>';    !
!       write(icanal,*);                                                                 !
                                                                                         !
        do j = 1, NATOM_LIBRARY(1);                                                      !
                                                                                         !
            IOSEF3 = CONFIG_ATOM_TYPE_LIBRARY(j,1);                                      !
                                                                                         !
            ROSEF3 = ATOM_MASSES_LIBRARY(IOSEF3,1);                                      !
                                                                                         !
            ROSEF5 = CONFIG_QI_LIBRARY(j,1);                                             !
                                                                                         !
!           write(icanal,*) '    o ', ATOM_LABEL_LIBRARY(IOSEF3,1), CONFIG_ATOM_TYPE_LIBRARY(j,1), ROSEF3, CONFIG_QI_LIBRARY(j,1);
!           write(icanal,*);
                                                                                         !
            ROSEF6 = ABS( ROSEF3 - ROSEF1 );                                             !
                                                                                         !
            if ( ROSEF6 > 0.1d0 ) CYCLE;                                                 !
                                                                                         !
            IOSEF2 = IOSEF3;                                                             !

            ROSEF2 = ROSEF3;

            ROSEF4 = ROSEF5;
            
            EXIT;
 
        end do

        if ( IOSEF2 == 0 ) CYCLE; 

!       write(icanal,*) 'XXX --> ', ROSEF4;
!       write(icanal,*);

        CONFIG_QI(i) = ROSEF4;                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   ### Assign bond types  #########################################################################
                                                                                         !
    if ( NTYPE_BOND == 1 ) then;                                                         !
                                                                                         !
        BOND_TYPE(1:NBOND) = 1;                                                          ! 
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Assign angle types #########################################################################
                                                                                         !
    if ( NTYPE_ANGLE == 1 ) then;                                                        ! 
                                                                                         !
        ANGLE_TYPE(1:NANGLE) = 1;                                                        !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Assign parameters for the input file for lammps ############################################
                                                                                         !
    if ( ILIBRARY_FILE_INFO(1) == 1 ) then;                                              !
                                                                                         !
       CH_UNITS_STYLE = CH_UNITS_STYLE_LIBRARY(1);                                       !
                                                                                         !
       CH_ATOM_STYLE = CH_ATOM_STYLE_LIBRARY(1);                                         !
                                                                                         !
       if ( NTYPE_BOND > 0 ) then;                                                       !
                                                                                         !
           CH_BOND_STYLE = CH_BOND_STYLE_LIBRARY(1);                                     !
                                                                                         !
       end if                                                                            !
                                                                                         !
       if ( NTYPE_ANGLE > 0 ) then;                                                      !
                                                                                         !
           CH_ANGLE_STYLE = CH_ANGLE_STYLE_LIBRARY(1);                                   !
                                                                                         !
       end if                                                                            !
                                                                                         !
       CH_SPECIAL_BONDS = 'XXX';                                                         !
                                                                                         !
       if ( TRIM(CH_SPECIAL_BONDS_LIBRARY(1)) /= 'XXX' ) then;                           !
                                                                                         !
           CH_SPECIAL_BONDS = CH_SPECIAL_BONDS_LIBRARY(1);                               !
                                                                                         !
       end if                                                                            !
                                                                                         !
       POTENTIAL_CLASS2_CHTYPE = 'XXX';                                                  !
                                                                                         !
       if ( TRIM(POTENTIAL_CLASS2_CHTYPE_LIBRARY(1)) /= 'XXX' ) then;                    !
                                                                                         !
           POTENTIAL_CLASS2_CHTYPE = POTENTIAL_CLASS2_CHTYPE_LIBRARY(1);                 !
                                                                                         !
       end if                                                                            !
                                                                                         !
       CH_PAIR_MODIFY = 'XXX';                                                           !
                                                                                         !
       if ( TRIM(CH_PAIR_MODIFY_LIBRARY(1)) /= 'XXX' ) then;                             !
                                                                                         !
           CH_PAIR_MODIFY = CH_PAIR_MODIFY_LIBRARY(1);                                   !
                                                                                         !
       end if                                                                            !
                                                                                         !
       CH_KSPACE_STYLE = 'XXX';                                                          !
                                                                                         !
       if ( TRIM(CH_KSPACE_STYLE_LIBRARY(1)) /= 'XXX' ) then;                            !
                                                                                         !
           CH_KSPACE_STYLE = CH_KSPACE_STYLE_LIBRARY(1);                                 !
                                                                                         !
       end if                                                                            !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine DUPLICATE_ASSIGN_FORCE_FIELD
