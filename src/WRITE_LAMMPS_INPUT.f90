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

subroutine WRITE_LAMMPS_INPUT(icanal,CHEXT)

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

    use module_config;

    use module_osef;

    use module_lmp_input;

    use module_inserts;

!   ************************************************************************************************

    implicit none;

!   ***********************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    character (len=150), intent(in) :: CHEXT;

!   ************************************************************************************************

    integer (kind=4) :: i, j, k;

    integer (kind=4) :: irun; 
 
    integer (kind=4) :: ILOCAL_INPUT_FIX;

    integer (kind=4) :: EOF, EOF2;

    integer (kind=4) :: RANDOM_GENERATOR_SEEDS

    real (kind=8) :: RSEED01, grnd; 

!   ************************************************************************************************

    integer (kind=4), allocatable, dimension(:) :: ILOCAL_TEMPK;

!   ************************************************************************************************

    character (len=250), allocatable, dimension(:) :: CHLOCAL_TEMPK_NAME, &
                                                      CHLOCAL_FIX_NAME;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the routine #############################################################
                                                                                         !
    CHTITLE = 'Write lammps input file';                                                 !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the name of the input file to generate ###############################################
                                                                                         !
    IOSEF1 = 67 - LEN_TRIM(CHEXT);                                                       !
                                                                                         !
    write(icanal,'(a70)') '| '//TRIM(CHEXT)//REPEAT(' ',IOSEF1)//'|';                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the LAMMPS input file ################################################################
                                                                                         !
    open(105,file=TRIM(CHEXT));                                                          !
                                                                                         !
    write(105,'(a100)') '### Input file for lammps simulations '//REPEAT('#',62);        !
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
!   ### Write temperatures of the simulation box ###################################################
                                                                                         !
    write(105,'(a100)') '### Temperature [K] '//REPEAT('#',80);                          !
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
    if ( NLMP_INPUT_RUN > 0 ) then;                                                      !
                                                                                         !
        do i = 1, NLMP_INPUT_RUN;                                                        !
                                                                                         !
            call CONVERT_INT_TO_CHAR(icanal,i,3,CHOSEF1);                                !
                                                                                         !
!           ### Build the lammps command setting the initial temperature ###########################
                                                                                         !
            write(CHOSEF2,'(f6.2)') LMP_INPUT_RUN_TEMPK_INIT(i);                         !
                                                                                         !
            IOSEF1 = 49 - 8 - 12 - 13 - LEN_TRIM(CHOSEF1);                               !
                                                                                         !
            IOSEF2 = 22 - 5 + 3 - LEN_TRIM(CHOSEF2);                                     !
                                                                                         !
            CHOSEF3 = 'variable'//          &                                            !
                      REPEAT(' ',12)//      &                                            !
                      'my_temp_init_'//     &                                            !
                      TRIM(CHOSEF1)//       &                                            !
                      REPEAT(' ',IOSEF1)//  &                                            !
                      'equal'//             &                                            !
                      REPEAT(' ',IOSEF2)//  &                                            !
                      TRIM(CHOSEF2);                                                     !
                                                                                         !
            IOSEF3 = LEN_TRIM(CHOSEF3);                                                  !
                                                                                         !
            call CONVERT_INT_TO_CHAR(icanal,IOSEF3,3,CHOSEF4);                           !
                                                                                         ! 
            CHOSEF5 = 'a'//TRIM(CHOSEF4);                                                !
                                                                                         !
!           ### Write the lammps command setting the initial temperature ###########################
                                                                                         !
            write(105,'('//TRIM(CHOSEF5)//')') TRIM(CHOSEF3);                            ! 
                                                                                         !
            write(105,*);                                                                ! 
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           ### Build the lammps command setting the final temperature #############################
                                                                                         !
            write(CHOSEF2,'(f6.2)') LMP_INPUT_RUN_TEMPK_FINAL(i);                        !
                                                                                         !
            IOSEF1 = 49 - 8 - 12 - 14 - LEN_TRIM(CHOSEF1);                               !
                                                                                         !
            IOSEF2 = 22 - 5 + 3 - LEN_TRIM(CHOSEF2);                                     !
                                                                                         !
            CHOSEF3 = 'variable'//          &                                            !
                      REPEAT(' ',12)//      &                                            !
                      'my_temp_final_'//    &                                            !
                      TRIM(CHOSEF1)//       &                                            !
                      REPEAT(' ',IOSEF1)//  &                                            !
                      'equal'//             &                                            !
                      REPEAT(' ',IOSEF2)//  &                                            !
                      TRIM(CHOSEF2);                                                     !
                                                                                         !
            IOSEF3 = LEN_TRIM(CHOSEF3);                                                  !
                                                                                         !
            call CONVERT_INT_TO_CHAR(icanal,IOSEF3,3,CHOSEF4);                           !
                                                                                         ! 
            CHOSEF5 = 'a'//TRIM(CHOSEF4);                                                !
                                                                                         !
!           ### Write the lammps command setting the final temperature #############################
                                                                                         !
            write(105,'('//TRIM(CHOSEF5)//')') TRIM(CHOSEF3);                            ! 
                                                                                         !
            write(105,*);                                                                ! 
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end do;                                                                          !
                                                                                         !
    else                                                                                 !
                                                                                         !
        write(105,'(a64)') 'variable'//REPEAT(' ',12)//     &                            !
                           'my_temp_init'//REPEAT(' ',10)// &                            !
                           'equal'//REPEAT(' ',16)//'1';                                 !
                                                                                         !
        write(105,*);                                                                    ! 
                                                                                         !
        write(105,'(a64)') 'variable'//REPEAT(' ',12)//     &                            !
                           'my_temp_final'//REPEAT(' ',9)// &                            !
                           'equal'//REPEAT(' ',14)//'400';                               !
                                                                                         !
        write(105,*);                                                                    !
                                                                                         !
    end if                                                                               ! 
                                                                                         !
!   ### Write pressures of the simulation box ######################################################
                                                                                         !
    write(105,'(a100)') '### Pressure in [atm] '//REPEAT('#',78);                        !
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
    write(105,'(a65)') '# variable'//REPEAT(' ',11)//      &                             !
                       'my_pressure_init'//REPEAT(' ',6)// &                             !
                       'equal'//REPEAT(' ',16)//'1';                                     !
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !  
    write(105,'(a65)') '# variable'//REPEAT(' ',11)//       &                            !
                       'my_pressure_final'//REPEAT(' ',5)// &                            !
                       'equal'//REPEAT(' ',16)//'1';                                     !
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the number of steps requested for each run ###########################################
                                                                                         !
    write(105,'(a100)') '### Number of steps '//REPEAT('#',80);                          !
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
    if ( NLMP_INPUT_RUN > 0 ) then;                                                      !
                                                                                         !
        do i = 1, NLMP_INPUT_RUN;                                                        !
                                                                                         !
            call CONVERT_INT_TO_CHAR(icanal,i,3,CHOSEF1);                                !
                                                                                         !
!           ### Build the lammps command setting the number of steps ###############################
                                                                                         !
            call CONVERT_INT_TO_CHAR(icanal,LMP_INPUT_RUN_NSTEPS(i),7,CHOSEF2);          !
                                                                                         !
            IOSEF1 = 49 - 8 - 12 - 7 - LEN_TRIM(CHOSEF1);                                !
                                                                                         !
            IOSEF2 = 22 - 5 - LEN_TRIM(CHOSEF2);                                         !
                                                                                         !
            CHOSEF3 = 'variable'//          &                                            !
                      REPEAT(' ',12)//      &                                            !
                      'Nsteps_'//           &                                            !
                      TRIM(CHOSEF1)//       &                                            !
                      REPEAT(' ',IOSEF1)//  &                                            !
                      'equal'//             &                                            !
                      REPEAT(' ',IOSEF2)//  &                                            !
                      TRIM(CHOSEF2);                                                     !
                                                                                         !
            IOSEF3 = LEN_TRIM(CHOSEF3);                                                  !
                                                                                         !
            call CONVERT_INT_TO_CHAR(icanal,IOSEF3,3,CHOSEF4);                           !
                                                                                         !
            CHOSEF5 = 'a'//TRIM(CHOSEF4);                                                !
                                                                                         !
!           ### Write the lammps command setting the initial temperature ###########################
                                                                                         !
            write(105,'('//TRIM(CHOSEF5)//')') TRIM(CHOSEF3);                            ! 
                                                                                         !
            write(105,*);                                                                ! 
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
       end do                                                                            !
                                                                                         !
    else                                                                                 !
                                                                                         !
        write(105,'(a64)') 'variable'//              &                                   !
                           REPEAT(' ',12)//          &                                   !
                           'Nsteps'//                &                                   !
                           REPEAT(' ',16)//          &                                   !
                           'equal          1000000';                                     !
                                                                                         !
        write(105,*);                                                                    !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Write sampling parameters ##################################################################
                                                                                         !
    write(105,'(a71)') 'variable'//REPEAT(' ',12)// &                                    !
                       'Sampling_frequency_thermo'//REPEAT(' ',4)// &                    !
                       'equal             1000';                                         !
                                                                                         !
    write(105,*);                                                                        !  
                                                                                         !
    write(105,'(a71)') 'variable'//REPEAT(' ',12)// &                                    !
                       'Sampling_frequency'//REPEAT(' ',11)//'equal             1000';   !
                                                                                         !
    write(105,*);                                                                        !  
                                                                                         !
    write(105,'(a71)') 'variable'//REPEAT(' ',12)// &                                    !
                       'Restart_frequency'//REPEAT(' ',12)//'equal          1000000';    !
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
    write(105,'(a100)') '### Simulation properties and parameters '//REPEAT('#',59);     !
                                                                                         ! 
    write(105,*);                                                                        !
                                                                                         !
!   ### Write the unit style in the input file #####################################################
                                                                                         !
    CHOSEF1 = 'units'//REPEAT(' ',15)//TRIM(CH_UNITS_STYLE);                             !
                                                                                         !
    IOSEF1 = LEN_TRIM(CHOSEF1);                                                          !
                                                                                         !
    call CONVERT_INT_TO_CHAR(icanal,IOSEF1,3,CHOSEF2);                                   !
                                                                                         !
    write(105,'(a'//TRIM(CHOSEF2)//')') TRIM(CHOSEF1);                                   !
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
    write(105,'(a21)') 'dimension'//REPEAT(' ',11)//'3';                                 !
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
    write(105,'(a25)') '# newton              off';                                      !
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
    write(105,'(a28)')  '# processors         *  *  *';                                  !
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
    write(105,'(a27)')  'boundary'//REPEAT(' ',12)//'p  p  p';                           !
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
!   ### Write the atom style #######################################################################
                                                                                         !
    CHOSEF1 =  'atom_style'//REPEAT(' ',10)//TRIM(CH_ATOM_STYLE);                        !
                                                                                         !
    IOSEF1 = LEN_TRIM(CHOSEF1);                                                          !
                                                                                         !
    call CONVERT_INT_TO_CHAR(icanal,IOSEF1,3,CHOSEF2);                                   !
                                                                                         !
    write(105,'(a'//TRIM(CHOSEF2)//')') TRIM(CHOSEF1);                                   ! 
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
!   ### Write the bond style #######################################################################
                                                                                         !
    CHOSEF1 = 'bond_style'//REPEAT(' ',10)//TRIM(CH_BOND_STYLE);                         !
                                                                                         !
    IOSEF1 = LEN_TRIM(CHOSEF1);                                                          !
                                                                                         !
    call CONVERT_INT_TO_CHAR(icanal,IOSEF1,3,CHOSEF2);                                   !
                                                                                         !
    write(105,'(a'//TRIM(CHOSEF2)//')') TRIM(CHOSEF1);                                   ! 
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
!   ### Write the angle style ######################################################################
                                                                                         !
    CHOSEF1 = 'angle_style'//REPEAT(' ',9)//TRIM(CH_ANGLE_STYLE);                        !
                                                                                         !
    IOSEF1 = LEN_TRIM(CHOSEF1);                                                          !
                                                                                         !
    call CONVERT_INT_TO_CHAR(icanal,IOSEF1,3,CHOSEF2);                                   !
                                                                                         !
    write(105,'(a'//TRIM(CHOSEF2)//')') TRIM(CHOSEF1);                                   ! 
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
!   ### Write the dihedral style in the lammps input ###############################################
                                                                                         !
    if ( NTYPE_DIHEDRAL > 0 ) then;                                                      !
                                                                                         !
        CHOSEF1 = 'dihedral_style'//REPEAT(' ',6)//TRIM(CH_DIHEDRAL_STYLE);              !
                                                                                         !
        IOSEF1 = LEN_TRIM(CHOSEF1);                                                      !
                                                                                         !
        call CONVERT_INT_TO_CHAR(icanal,IOSEF1,3,CHOSEF2);                               !
                                                                                         !
        write(105,'(a'//TRIM(CHOSEF2)//')') TRIM(CHOSEF1);                               !  
                                                                                         !
        write(105,*);                                                                    !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Write special bonds ########################################################################
                                                                                         !
    if ( TRIM(CH_SPECIAL_BONDS) /= 'XXX' ) then;                                         !
                                                                                         !
        CHOSEF1 = 'special_bonds'//REPEAT(' ',7)//TRIM(CH_SPECIAL_BONDS);                !
                                                                                         !
        IOSEF1 = LEN_TRIM(CHOSEF1);                                                      !
                                                                                         !
        call CONVERT_INT_TO_CHAR(icanal,IOSEF1,3,CHOSEF2);                               !
                                                                                         !
        write(105,'(a'//TRIM(CHOSEF2)//')') TRIM(CHOSEF1);                               !  
                                                                                         !
        write(105,*);                                                                    !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Write pair style ###########################################################################
                                                                                         !
    if ( TRIM(POTENTIAL_CLASS2_CHTYPE) /= 'XXX' ) then;                                  !
                                                                                         !
        CHOSEF1 = 'pair_style'//REPEAT(' ',10)//TRIM(POTENTIAL_CLASS2_CHTYPE)// &        !
                  '  11.0  12.0  12.0';                                                  !
                                                                                         !
        IOSEF1 = LEN_TRIM(CHOSEF1);                                                      !
                                                                                         !
        call CONVERT_INT_TO_CHAR(icanal,IOSEF1,3,CHOSEF2);                               !
                                                                                         !
        write(105,'(a'//TRIM(CHOSEF2)//')') TRIM(CHOSEF1);                               !  
                                                                                         !
        write(105,*);                                                                    !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Write pair modify ##########################################################################
                                                                                         !
    if ( TRIM(CH_PAIR_MODIFY) /= 'XXX' ) then;                                           !
                                                                                         !
        CHOSEF1 = 'pair_modify'//REPEAT(' ',9)//TRIM(CH_PAIR_MODIFY);                    !
                                                                                         !
        IOSEF1 = LEN_TRIM(CHOSEF1);                                                      !
                                                                                         !
        call CONVERT_INT_TO_CHAR(icanal,IOSEF1,3,CHOSEF2);                               !
                                                                                         !
        write(105,'(a'//TRIM(CHOSEF2)//')') TRIM(CHOSEF1);                               !  
                                                                                         !
        write(105,*);                                                                    !
                                                                                         !
    end if                                                                               !  
                                                                                         !
!   ### Set the method for Coulomb interactions ####################################################
                                                                                         !
    write(105,'(a100)') '### Set the method for Coulomb interactions '//REPEAT('#',56);  !
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
    if ( TRIM(CH_KSPACE_STYLE) /= 'XXX' ) then;                                          !
                                                                                         !
        CHOSEF1 = 'kspace_style'//REPEAT(' ',8)//TRIM(CH_KSPACE_STYLE)//' 1.0e-5';       !
                                                                                         !
        IOSEF1 = LEN_TRIM(CHOSEF1);                                                      !
                                                                                         !
        call CONVERT_INT_TO_CHAR(icanal,IOSEF1,3,CHOSEF2);                               !
                                                                                         !
        write(105,'(a'//TRIM(CHOSEF2)//')') TRIM(CHOSEF1);                               ! 
                                                                                         !
        write(105,*);                                                                    !
                                                                                         !
    end if                                                                               !
                                                                                         !
    write(105,'(a100)') '### Read the initial molecular configuration '//REPEAT('#',55); !
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
    write(105,'(a30)') 'read_data'//REPEAT(' ',11)//'atoms.data';                        !
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
!   ### Write pair coefficients for cross interactions #############################################
                                                                                         !
    if ( NPAIR_COEFF_CROSS > 0 ) then;                                                   !
                                                                                         !
        write(105,'(a100)') '### Pair coefficients for cross interactions '// &          !
                            REPEAT('#',55);                                              !
                                                                                         !
        write(105,*);                                                                    !
                                                                                         !
        do i = 1, NPAIR_COEFF_CROSS;                                                     !
                                                                                         !
            if ( ( PAIR_COEFF_CROSS(1,i) == 0.0d0 ) .AND.    &                           ! Write only non-zero cross interactions for van der Waals interactions
                 ( PAIR_COEFF_CROSS(2,i) == 0.0d0 ) ) CYCLE;                             !
                                                                                         !
            write(105,'(a17,2i4,2f15.6)') 'pair_coeff'//REPEAT(' ',7),   &               !
                                          PAIR_ATOMID_CROSS(1:2,i),      &               !
                                          PAIR_COEFF_CROSS(1:2,i);                       !
                                                                                         ! 
        end do                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
!   ### Write minimization command #################################################################
                                                                                         !
    write(105,'(a100)') '### Minimization step '//REPEAT('#',78);                        !
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
    write(105,'(a47)') '# minimize'//REPEAT(' ',12)//'1.0e-4  1.0e-6  100  1000';        !
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
    write(105,'(a23)') '# reset_timestep'//REPEAT(' ',6)//'0';                           !
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
    write(105,'(a28)') 'neighbor'//REPEAT(' ',12)//'2.0  bin';                           !
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
    write(105,'(a52)') '# neigh_modify'//REPEAT(' ',8)// &                               !
                       'delay  0  every  1  check  yes';                                 !
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
!   ### Write commands for setting groups ##########################################################
                                                                                         !
    write(105,'(a100)') '### Setting groups for the simulation '//REPEAT('#',62);        !
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
    do i = 1, NLMP_INPUT_GROUP;                                                          !
                                                                                         !
        CHOSEF1 = 'group'//REPEAT(' ',15)//TRIM(LMP_INPUT_GROUP_NAME(i));                !
                                                                                         !
        IOSEF1 = 35 - LEN_TRIM(CHOSEF1);                                                 !
                                                                                         !
        if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                    !
                                                                                         !
        CHOSEF1 = TRIM(CHOSEF1)//REPEAT(' ',IOSEF1)//' type';                            !
                                                                                         !
        do j = 1, LMP_INPUT_GROUP_NLIST(i);                                              !
                                                                                         !
            CHOSEF2 = TRIM(LMP_INPUT_GROUP_LIST(j,i));                                   ! Get the name of the group to consider in the array
                                                                                         !
            do k = 1, NLAMMPS_INSERTS;                                                   ! Loop over the number of insert commands
                                                                                         !
                if ( TRIM(CHOSEF2) == TRIM(LAMMPS_INSERT_NAME(k)) ) then;                ! If the name of the insert command matches with the group to consider, then get the insert command number
                                                                                         !
                    IOSEF2 = k;                                                          !
                                                                                         !
                    EXIT;                                                                !
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            IOSEF3 = 1;                                                                  !
                                                                                         !
            if ( IOSEF2 > 1 ) IOSEF3 = LAMMPS_INSERT_NTYPE_ATOM_MAX(IOSEF2-1) + 1;       !
                                                                                         !
            IOSEF4 = LAMMPS_INSERT_NTYPE_ATOM_MAX(IOSEF2);                               !
                                                                                         !
            do k = IOSEF3, IOSEF4;                                                       !
                                                                                         !
                call CONVERT_INT_TO_CHAR(icanal,k,3,CHOSEF3);                            !
                                                                                         !
                CHOSEF1 = TRIM(CHOSEF1)//'  '//TRIM(CHOSEF3);                            !
                                                                                         ! 
            end do                                                                       !
                                                                                         !
        end do                                                                           !
                                                                                         !
        IOSEF2 = LEN_TRIM(CHOSEF1);                                                      !
                                                                                         !
        call CONVERT_INT_TO_CHAR(icanal,IOSEF2,3,CHOSEF2);                               !
                                                                                         !
        CHOSEF2 = 'a'//TRIM(CHOSEF2);                                                    !
                                                                                         !
        write(105,'('//TRIM(CHOSEF2)//')') TRIM(CHOSEF1);                                !
                                                                                         !
        write(105,*);                                                                    !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   ### Write commands for group unions ############################################################
                                                                                         !
    do i = 1, NLMP_INPUT_GROUP_UNION;                                                    !
                                                                                         !
        CHOSEF1 = 'group'//REPEAT(' ',15)//TRIM(LMP_INPUT_GROUP_UNION_NAME(i));          !
                                                                                         !
        IOSEF1 = 35 - LEN_TRIM(CHOSEF1);                                                 !
                                                                                         !
        if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                    !
                                                                                         !
        CHOSEF1 = TRIM(CHOSEF1)//REPEAT('  ',IOSEF1)//' union';                          !
                                                                                         !
        do j = 1, LMP_INPUT_GROUP_UNION_NLIST(i);                                        !
                                                                                         !
            CHOSEF2 = TRIM(LMP_INPUT_GROUP_UNION_LIST(j,i));                             ! Get the name of the group to consider in the array
                                                                                         !
            CHOSEF1 = TRIM(CHOSEF1)//'  '//TRIM(CHOSEF2);                                !
                                                                                         !
        end do                                                                           !
                                                                                         !
        IOSEF2 = LEN_TRIM(CHOSEF1);                                                      !
                                                                                         !
        call CONVERT_INT_TO_CHAR(icanal,IOSEF2,3,CHOSEF2);                               !
                                                                                         !
        CHOSEF2 = 'a'//TRIM(CHOSEF2);                                                    !
                                                                                         !
        write(105,'('//TRIM(CHOSEF2)//')') TRIM(CHOSEF1);                                !
                                                                                         !
        write(105,*);                                                                    !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   ### Write simulation properties ################################################################
                                                                                         !
    write(105,'(a100)') '### Simulation properties '//REPEAT('#',74);                    !
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
    write(105,'(a26)') 'run_style'//REPEAT(' ',11)//'verlet';                            ! 
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
!   ### Write the timestep #########################################################################
                                                                                         !
    write(105,'(a23)') 'timestep'//REPEAT(' ',12)//'1.0';                                !
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
!   ### Write LAMMPS command for frozen entities ###################################################
                                                                                         !
    if ( NLMP_INPUT_FROZEN > 0 ) then;                                                   !
                                                                                         !
        do i = 1, NLMP_INPUT_FROZEN;                                                     !
                                                                                         !
            CHOSEF3 = 'velocity'//                              &                        !
                      REPEAT(' ',12)//                          &                        !
                      TRIM(LMP_INPUT_FROZEN_GRPNAME(i))//'  '// &                        !
                      'set  0.0  0.0  0.0';                                              !
                                                                                         !
            IOSEF3 = LEN_TRIM(CHOSEF3);                                                  !
                                                                                         !
            call CONVERT_INT_TO_CHAR(icanal,IOSEF3,3,CHOSEF4);                           !
                                                                                         ! 
            CHOSEF5 = 'a'//TRIM(CHOSEF4);                                                !
                                                                                         !
!           ### Write the lammps command setting the initial temperature ###########################
                                                                                         !
            write(105,'('//TRIM(CHOSEF5)//')') TRIM(CHOSEF3);                            ! 
                                                                                         !
            write(105,*);                                                                ! 
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end do                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NLMP_INPUT_FROZEN > 0 ) then;                                                   !
                                                                                         !
        do i = 1, NLMP_INPUT_FROZEN;                                                     !
                                                                                         !
            call CONVERT_INT_TO_CHAR(icanal,i,3,CHOSEF1);                                !
                                                                                         !
            CHOSEF3 = 'fix'//                                   &                        !
                      REPEAT(' ',17)//                          &                        !
                      'freeze_'//                               &                        !
                      TRIM(CHOSEF1)//'  '//                     &                        !
                      TRIM(LMP_INPUT_FROZEN_GRPNAME(i))//'  '// &                        !
                      'setforce  0.0  0.0  0.0';                                         ! 
                                                                                         !
            IOSEF3 = LEN_TRIM(CHOSEF3);                                                  !
                                                                                         !
            call CONVERT_INT_TO_CHAR(icanal,IOSEF3,3,CHOSEF4);                           !
                                                                                         ! 
            CHOSEF5 = 'a'//TRIM(CHOSEF4);                                                !
                                                                                         !
!           ### Write the lammps command setting the initial temperature ###########################
                                                                                         !
            write(105,'('//TRIM(CHOSEF5)//')') TRIM(CHOSEF3);                            ! 
                                                                                         !
            write(105,*);                                                                ! 
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end do                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Set and write commands for generating initial velocities ###################################
                                                                                         !
    write(105,'(a100)') '### Set initial velocities '//REPEAT('#',73);                   !
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
    if ( NLMP_INPUT_FIX > 0 ) then;                                                      !
                                                                                         !
        do i = 1, NLMP_INPUT_FIX;                                                        !
                                                                                         !
!           ### Generation of a seed for velocity creation #########################################
                                                                                         !
            RANDOM_GENERATOR_SEEDS = 0;                                                  !
                                                                                         !
            do while ( RANDOM_GENERATOR_SEEDS == 0 )                                     !
                                                                                         !
                RSEED01 = grnd() !  rand();                                                        !
                                                                                         !
                ROSEF1 = grnd() !  rand();                                                         !
                                                                                         !
                IOSEF1 = INT( ROSEF1 * 10.0d0 );                                         !
                                                                                         !
                do j = 1, IOSEF1;                                                        !
                                                                                         !
                    RSEED01 = RSEED01 * 10.0d0;                                          !
                                                                                         !
                end do                                                                   !
                                                                                         !
                RANDOM_GENERATOR_SEEDS = INT(RSEED01);                                   !
                                                                                         ! 
            end do                                                                       !
                                                                                         !
!           ### Writing the lammps command for the generation of atomic velocities #################
                                                                                         !
            IOSEF10 = ANINT(LOG10(REAL(IOSEF9)));                                        !
                                                                                         !
            call CONVERT_INT_TO_CHAR(icanal,RANDOM_GENERATOR_SEEDS,10,CHOSEF2);          !
                                                                                         !
            CHOSEF3 = 'velocity'//                             &                         !
                      REPEAT(' ',12)//                         &                         !
                      TRIM(LMP_INPUT_RUN_GRPNAME(i))//'  '//   &                         !
                      'create  ${my_temp_init_1} '//           &                         !
                      TRIM(CHOSEF2);                                                     !  
                                                                                         !
            IOSEF3 = LEN_TRIM(CHOSEF3);                                                  !
                                                                                         !
            call CONVERT_INT_TO_CHAR(icanal,IOSEF3,3,CHOSEF4);                           !
                                                                                         ! 
            CHOSEF5 = 'a'//TRIM(CHOSEF4);                                                !
                                                                                         !
!           ### Write the lammps command setting the initial temperature ###########################
                                                                                         !
            write(105,'('//TRIM(CHOSEF5)//')') TRIM(CHOSEF3);                            ! 
                                                                                         !
            write(105,*);                                                                !    
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end do                                                                           !
                                                                                         !
    else                                                                                 ! 
                                                                                         !
        RANDOM_GENERATOR_SEEDS = 2348934;                                                !
                                                                                         !
        write(105,'(a51,i8)') 'velocity'//                       &                       !
                              REPEAT(' ',12)//                   &                       !
                              'all  create  ${my_temp_init_1} ', &                       !
                              RANDOM_GENERATOR_SEEDS;                                    ! 
                                                                                         !
        write(105,*);                                                                    !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate local arrays for the writing of additional options ################################ 
                                                                                         !
!   if ( NLMP_INPUT_RUN > 0 ) then;                                                      !
    if ( NLMP_INPUT_FIX > 0 ) then;                                                      !
                                                                                         !
!      allocate(ILOCAL_TEMPK(1:NLMP_INPUT_RUN));                                         !
       allocate(ILOCAL_TEMPK(1:NLMP_INPUT_FIX));                                         !
                                                                                         !
!      allocate(CHLOCAL_TEMPK_NAME(1:NLMP_INPUT_RUN));                                   !
       allocate(CHLOCAL_TEMPK_NAME(1:NLMP_INPUT_FIX));                                   !
                                                                                         !
    end if                                                                               !
                                                                                         !
    ILOCAL_INPUT_FIX = NLMP_INPUT_FIX;                                                   !
                                                                                         !
    if ( NLMP_INPUT_FIX == 0 ) ILOCAL_INPUT_FIX =  1;                                    !
                                                                                         !
    allocate(CHLOCAL_FIX_NAME(1:ILOCAL_INPUT_FIX));                                      !
                                                                                         !
!   ### Initialization of local arrays for the writing of additional options #######################
                                                                                         !
    if ( NLMP_INPUT_FIX > 0 ) then;                                                      !
                                                                                         !
       ILOCAL_TEMPK(1:NLMP_INPUT_FIX) = 0;                                               !
                                                                                         !
       CHLOCAL_TEMPK_NAME(1:NLMP_INPUT_FIX) = 'XXX';                                     !
                                                                                         !
    end if                                                                               !
                                                                                         !
    CHLOCAL_FIX_NAME(1:ILOCAL_INPUT_FIX) = 'XXX';                                        !
                                                                                         !
!   ### Loop over the number of runs ###############################################################
                                                                                         !
    if ( NLMP_INPUT_RUN > 0 ) then;                                                      !
                                                                                         !
        do irun = 1, NLMP_INPUT_RUN;                                                     !
                                                                                         !
            call CONVERT_INT_TO_CHAR(icanal,irun,3,CHOSEF1);                             !
                                                                                         !
!           ### Write fix commands for the statistical ensemble ####################################
                                                                                         !
            write(105,'(a100)') '### Set fix commands for '//   &                        !
                                'equilibrium MD simulations '// &                        !
                                REPEAT('#',48);                                          !
                                                                                         !
            write(105,*);                                                                !
                                                                                         !
            if ( NLMP_INPUT_FIX > 0 ) then;                                              !
                                                                                         !
                ILOCAL_TEMPK(1:NLMP_INPUT_FIX) = 0;                                      !
                                                                                         !
                do i = 1, NLMP_INPUT_FIX;                                                !
                                                                                         !
!                   ### Build the fix name for LAMMPS ##############################################
                                                                                         !
                    call CONVERT_INT_TO_CHAR(icanal,i,3,CHOSEF6);                        !
                                                                                         !
                    CHLOCAL_FIX_NAME(i) = 'TSTAT_'//TRIM(CHOSEF6);                       !
                                                                                         !
!                   ### Build the fix command setting the statistical ensemble #####################
                                                                                         !
                    call CONVERT_REAL_TO_CHAR(icanal,LMP_INPUT_RUN_TAUNH(i),3,CHOSEF2);  !
                                                                                         !
                    CHOSEF3 = 'fix'//                                    &               !
                              REPEAT(' ',17)//                           &               !
                              TRIM(CHLOCAL_FIX_NAME(i))//'  '//          &               !
                              TRIM(LMP_INPUT_RUN_GRPNAME(i))//'  '//     &               !
                              TRIM(LMP_INPUT_RUN_STYLE(i))//'  '//       &               !
                              'temp  '//                                 &               !
                              '${my_temp_init_'//TRIM(CHOSEF1)//'}  '//  &               !
                              '${my_temp_final_'//TRIM(CHOSEF1)//'}  '// &               !
                              TRIM(CHOSEF2);                                             !
                                                                                         !
                    IOSEF3 = LEN_TRIM(CHOSEF3);                                          !
                                                                                         !
                    call CONVERT_INT_TO_CHAR(icanal,IOSEF3,3,CHOSEF4);                   !
                                                                                         ! 
                    CHOSEF5 = 'a'//TRIM(CHOSEF4);                                        !
                                                                                         !
!                   ### Write the lammps command setting the initial temperature ###################
                                                                                         !
                    write(105,'('//TRIM(CHOSEF5)//')') TRIM(CHOSEF3);                    ! 
                                                                                         !
                    write(105,*);                                                        ! 
                                                                                         !
!                   stop; !//////////////////////////////////////////////////////////////!
                                                                                         !
!                   ### Additional command related to the computation of the temperature ###########
                                                                                         !
                    if ( TRIM(LMP_INPUT_RUN_GRPNAME(i)) == 'all' ) CYCLE;                !
                                                                                         !
                    ILOCAL_TEMPK(i) = 1;                                                 !
                                                                                         !
                    CHLOCAL_TEMPK_NAME(i) = 'temp_'//TRIM(LMP_INPUT_RUN_GRPNAME(i));     !
                                                                                         !
                    CHOSEF3 = 'compute'//                             &                  !
                              REPEAT(' ',13)//                        &                  !
                              TRIM(CHLOCAL_TEMPK_NAME(i))//'  '//     &                  !
                              TRIM(LMP_INPUT_RUN_GRPNAME(i))//'  '//  &                  !
                              'temp';                                                    !
                                                                                         !
                    IOSEF3 = LEN_TRIM(CHOSEF3);                                          !
                                                                                         !
                    call CONVERT_INT_TO_CHAR(icanal,IOSEF3,3,CHOSEF4);                   !
                                                                                         ! 
                    CHOSEF5 = 'a'//TRIM(CHOSEF4);                                        !
                                                                                         !
!                   ### Write the lammps command setting the initial temperature ###################
                                                                                         !
                    write(105,'('//TRIM(CHOSEF5)//')') TRIM(CHOSEF3);                    ! 
                                                                                         !
                    write(105,*);                                                        ! 
                                                                                         !
!                   stop; !//////////////////////////////////////////////////////////////!
                                                                                         !
                end do                                                                   !
                                                                                         !
            else                                                                         !
                                                                                         !
                CHLOCAL_FIX_NAME(1) = 'TSTAT';                                           !
                                                                                         !
                write(105,'(a79)') 'fix'//REPEAT(' ',17)//                      &        !
                                   'TSTAT  all  nvt  temp  ${my_temp_init}  '// &        !
                                   '${my_temp_final}  1';                                !
                                                                                         !
                write(105,*);                                                            !
                                                                                         !
            end if                                                                       !
                                                                                         !
            write(105,'(a139)') '# fix'//REPEAT(' ',17)//                             &  !
                                'TPSTAT  all  npt  temp  ${my_temp_init}  '//         &  !
                                '${my_temp_final}  1    iso  ${my_pressure_init}  '// &  !
                                '${my_pressure_final}   1000';                           !
                                                                                         !
            write(105,*);                                                                !
                                                                                         !
!           ### Write commands related to thermodynamics properties ################################
                                                                                         !
            write(105,'(a100)') '### Set thermodynamics properties'// &                  !
                                ' to write in the output '//          &                  !
                                REPEAT('#',43);                                          ! 
                                                                                         !
            write(105,*);                                                                !
                                                                                         !
            write(105,'(a100)') 'thermo_style'//REPEAT(' ',8)//                &         !
                                'custom  step  cpu  etotal  ke  pe  evdwl  '// &         !
                                'ecoul  elong  temp press  vol  density';                !
                                                                                         !
            write(105,*);                                                                !
                                                                                         !
            if ( NLMP_INPUT_FIX > 0 ) then;                                              !
                                                                                         !
                do i = 1, NLMP_INPUT_FIX;                                                !
                                                                                         !
                    if ( ILOCAL_TEMPK(i) == 0 ) CYCLE;                                   !
                                                                                         !
                    CHOSEF3 = 'thermo_modify'//            &                             !
                              REPEAT(' ',7)//              &                             !
                              'temp  '//                   &                             !
                              TRIM(CHLOCAL_TEMPK_NAME(i));                               !
                                                                                         !
                    IOSEF3 = LEN_TRIM(CHOSEF3);                                          !
                                                                                         !
                    call CONVERT_INT_TO_CHAR(icanal,IOSEF3,3,CHOSEF4);                   !
                                                                                         ! 
                    CHOSEF5 = 'a'//TRIM(CHOSEF4);                                        !
                                                                                         !
!                   ### Write the lammps command setting the initial temperature ###################
                                                                                         !
                    write(105,'('//TRIM(CHOSEF5)//')') TRIM(CHOSEF3);                    !  
                                                                                         !
                    write(105,*);                                                        ! 
                                                                                         !
!                   stop; !//////////////////////////////////////////////////////////////!
                                                                                         !
                end do                                                                   !
                                                                                         !
            end if                                                                       !
                                                                                         !
            write(105,'(a48)') 'thermo'//REPEAT(' ',14)//      &                         !
                               '${Sampling_frequency_thermo}';                           !
                                                                                         !
            write(105,*);                                                                !
                                                                                         !
            write(105,'(a123)') 'dump'//REPEAT(' ',16)//'TRAJ  all  custom  '// &        !
                                '${Sampling_frequency}  dump.lammpstrj  id  '// &        !
                                'mol  type  element q  x  y  z  vx  vy  vz';             !
                                                                                         !
            write(105,*);                                                                !
                                                                                         !
!           ### Build and write lammps modify ######################################################

!           'dump_modify    TRAJ  element  N  C  C  C  H  H  H  C  C  H  C  C  F  S  N  O'

!           do i = 1, NTYPE_ATOM;

!           ATOM_MASSE(i)

!           end do

            write(105,'(a34)') 'dump_modify'//REPEAT(' ',9) //'TRAJ  sort  id';          !
                                                                                         !
            write(105,*);                                                                !
                                                                                         !
            write(105,'(a49)') 'restart'//REPEAT(' ',13)//      &                        !  
                               '${Restart_frequency}  restart';                          !
                                                                                         !
            write(105,*);                                                                !
                                                                                         !
!           ### Write the run command in the lammps input file #####################################
                                                                                         !
            write(105,'(a100)') '### Run the simulation '//REPEAT('#',77);               !
                                                                                         !
            write(105,*);                                                                !
                                                                                         !
            CHOSEF3 = 'run'//                           &                                !
                      REPEAT(' ',17)//                  &                                !
                      '${Nsteps_'//TRIM(CHOSEF1)//'}';                                   !
                                                                                         !
            IOSEF3 = LEN_TRIM(CHOSEF3);                                                  !
                                                                                         !
            call CONVERT_INT_TO_CHAR(icanal,IOSEF3,3,CHOSEF4);                           !
                                                                                         ! 
            CHOSEF5 = 'a'//TRIM(CHOSEF4);                                                !
                                                                                         !
!           ### Write the lammps command setting the initial temperature ###########################
                                                                                         !
            write(105,'('//TRIM(CHOSEF5)//')') TRIM(CHOSEF3);                            ! 
                                                                                         !
            write(105,*);                                                                ! 
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           ### Write unfix command ################################################################
                                                                                         !
            write(105,'(a100)') '### Set unfix commands '//REPEAT('#',77);               !
                                                                                         !
            write(105,*);                                                                !
                                                                                         !
            do i = 1, ILOCAL_INPUT_FIX;                                                  !
                                                                                         !
                CHOSEF3 = 'unfix'//                  &                                   !
                          REPEAT(' ',15)//           &                                   !
                          TRIM(CHLOCAL_FIX_NAME(i));                                     !
                                                                                         !
                IOSEF3 = LEN_TRIM(CHOSEF3);                                              !
                                                                                         !
                call CONVERT_INT_TO_CHAR(icanal,IOSEF3,3,CHOSEF4);                       !
                                                                                         ! 
                CHOSEF5 = 'a'//TRIM(CHOSEF4);                                            !
                                                                                         !
!               ### Write the lammps command setting the initial temperature #######################
                                                                                         !
                write(105,'('//TRIM(CHOSEF5)//')') TRIM(CHOSEF3);                        ! 
                                                                                         !
                write(105,*);                                                            !
                                                                                         !
            end do                                                                       !
                                                                                         !
        end do                                                                           ! End of the loop over the number of runs 
                                                                                         !
    end if                                                                               ! End if the if command for the loop over the number of runs 
                                                                                         !
!   ### Write the command for generating the final molecular configuration #########################
                                                                                         !
    write(105,'(a100)') '### Write the final molecular configuration '//REPEAT('#',56);  !
                                                                                         !
    write(105,*);                                                                        !
                                                                                         !
    write(105,'(a36)') 'write_data'//REPEAT(' ',10)//'atoms.data_final';                 !
                                                                                         !
    close(105);                                                                          !
                                                                                         !
!   ### Deallocate local arrays for the writing of additional options ############################## 
                                                                                         !
    if ( NLMP_INPUT_FIX > 0 ) then;                                                      !
                                                                                         !
       deallocate(ILOCAL_TEMPK);                                                         !
                                                                                         !
       deallocate(CHLOCAL_TEMPK_NAME);                                                   !
                                                                                         !
    end if                                                                               !
                                                                                         !
    deallocate(CHLOCAL_FIX_NAME);                                                        !
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine WRITE_LAMMPS_INPUT
