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


subroutine SEEK_FOR_LMP_INPUT_COMMAND(icanal)

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

    use module_slab_keywords;

    use module_lmp_input;

    use module_inserts;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

!   ************************************************************************************************

    integer (kind=4) :: h, i, j;

    integer (kind=4) :: NLOCAL_COMMAND,     &
                        ILOCAL_RUN,         &
                        ILOCAL_FIX,         &
                        ILOCAL_FROZEN,      &
                        ILOCAL_GROUP,       &
                        ILOCAL_GROUP_UNION;

!   ************************************************************************************************

    character (len=250) :: CHARLINE;

!   ************************************************************************************************

    integer (kind=4) :: EOF;

!   ************************************************************************************************
                                                                                         !
    NLOCAL_COMMAND = 0;                                                                  !
                                                                                         !
!   ### Count the number of times the same insert command is modified ##############################
                                                                                         !
    open(1,file='input_100.dat',status='old');                                           !
                                                                                         !
    EOF = 0;                                                                             !
                                                                                         !
    do                                                                                   !
                                                                                         !
        read(1,'(a)',iostat=EOF) CHARLINE;                                               !
                                                                                         !
        if ( EOF /= 0 ) EXIT;                                                            !
                                                                                         !
        if ( INDEX(CHARLINE,'#') > 0 ) CYCLE;                                            !
                                                                                         !
        if ( INDEX(CHARLINE,'classic ') > 0 ) EXIT;                                      !
                                                                                         !
        IOSEF2 = INDEX(CHARLINE,TRIM(SLAB_INPUT_KEYWORDS(6))//' ');                      ! The keyword correspond to the lmp_input command 
                                                                                         !
        if ( IOSEF2 > 20 ) CYCLE;                                                        !
                                                                                         !
        if ( IOSEF2 == 0 ) CYCLE;                                                        !
                                                                                         !
        NLOCAL_COMMAND = NLOCAL_COMMAND + 1;                                             !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
    end do                                                                               !
                                                                                         !
    close(1);                                                                            !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the number of lmp_input commands found in the input file #############################
                                                                                         !
    write(icanal,'(a39,i4,a27)') '| The number of lmp_input command is : ', &            !
                                 NLOCAL_COMMAND,                            &            !
                                 REPEAT(' ',26)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialization of the number of lammps run set in the lmp_input command ####################
                                                                                         !
    NLMP_INPUT_FIX = 0;                                                                  !
                                                                                         !
    NLMP_INPUT_RUN = 0;                                                                  !
                                                                                         !
    NLMP_INPUT_FROZEN = 0;                                                               !
                                                                                         !
!   ### Count the number of lammps run set in the lmp_input command ################################
                                                                                         !
    open(1,file='input_100.dat',status='old');                                           !
                                                                                         !
    EOF = 0;                                                                             !
                                                                                         !
    do                                                                                   !
                                                                                         !
        read(1,'(a)',iostat=EOF) CHARLINE;                                               !
                                                                                         !
        if ( EOF /= 0 ) EXIT;                                                            !
                                                                                         !
        if ( INDEX(CHARLINE,'#') > 0 ) CYCLE;                                            !
                                                                                         !
        if ( INDEX(CHARLINE,'classic ') > 0 ) EXIT;                                      !
                                                                                         !
        IOSEF2 = INDEX(CHARLINE,TRIM(SLAB_INPUT_KEYWORDS(6))//' ');                      ! The keyword correspond to the lmp_input command 
                                                                                         !
        if ( IOSEF2 > 20 ) CYCLE;                                                        !
                                                                                         !
        if ( IOSEF2 == 0 ) CYCLE;                                                        !
                                                                                         !
        read(CHARLINE,*) CHOSEF1, CHOSEF2;                                               !
                                                                                         !
        if ( TRIM(CHOSEF2) /= 'dynamics' ) CYCLE;                                        !
                                                                                         !
        read(CHARLINE,*) CHOSEF1, CHOSEF2, CHOSEF3, CHOSEF4;                             !
                                                                                         !
!       ### Check keyword properties in the 3rd rank of the command line ###########################
                                                                                         !
        select case(TRIM(CHOSEF3));                                                      !
                                                                                         !
            case('run');                                                                 !
                                                                                         !
                NLMP_INPUT_RUN = NLMP_INPUT_RUN + 1;                                     !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
!           case default;                                                                !
                                                                                         !
!               write(icanal,'(a70)') '| Not implemented - stop '// &                    !
!                                     REPEAT(' ',44)//'|';                               !
                                                                                         !
!               call CLOSING_PROGRAM(icanal,1);                                          !
                                                                                         !
        end select                                                                       !
                                                                                         !
!       ### Check keyword properties in the 4th rank of the command line ###########################
                                                                                         !
        select case(TRIM(CHOSEF4));                                                      !
                                                                                         !
            case('nvt');                                                                 !
                                                                                         !
                NLMP_INPUT_FIX = NLMP_INPUT_FIX + 1;                                     !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            case('frozen');                                                              !
                                                                                         !
                NLMP_INPUT_FROZEN = NLMP_INPUT_FROZEN + 1;                               !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
!           case default;                                                                !
                                                                                         !
!               write(icanal,'(a70)') '| Not implemented - stop '// &                    !
!                                     REPEAT(' ',44)//'|';                               !
                                                                                         !
!               call CLOSING_PROGRAM(icanal,1);                                          !
                                                                                         !
        end select                                                                       !
                                                                                         !
    end do                                                                               !
                                                                                         !
    close(1);                                                                            !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         ! 
!   ### Write the number of lammps run set in the lmp_input command ################################
                                                                                         !
    if ( NLMP_INPUT_FIX > 0 ) then;                                                      !
                                                                                         !
        write(icanal,'(a44,i4,a22)') '| The number of lammps fixes to be set is : ', &   !
                                     NLMP_INPUT_FIX,                                 &   !
                                     REPEAT(' ',21)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the number of frozen entities in the simulation box ##################################
                                                                                         !
    if ( NLMP_INPUT_FROZEN > 0 ) then;                                                   !
                                                                                         !
        write(icanal,'(a37,i4,a29)') '| The number of frozen entities is : ', &          !
                                     NLMP_INPUT_FROZEN,                       &          !
                                     REPEAT(' ',28)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the number of runs in the simulation box #############################################
                                                                                         !
    if ( NLMP_INPUT_RUN > 0 ) then;                                                      !
                                                                                         !
        write(icanal,'(a26,i4,a40)') '| The number of runs is : ', &                     !
                                     NLMP_INPUT_RUN,               &                     !
                                     REPEAT(' ',39)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set the flag for the generation of the lammps input file ###################################
                                                                                         !
    if ( NLMP_INPUT_FIX > 0 ) then;                                                      !
                                                                                         !
        igenerate_lammps_input = 1;                                                      !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Allocate arrays related to fix properties ##################################################
                                                                                         !
    if ( NLMP_INPUT_FIX > 0 ) then;                                                      !
                                                                                         !
        allocate(LMP_INPUT_RUN_STYLE(1:NLMP_INPUT_FIX));                                 !
                                                                                         !
        allocate(LMP_INPUT_RUN_TAUNH(1:NLMP_INPUT_FIX));                                 !
                                                                                         !
        allocate(LMP_INPUT_RUN_GRPNAME(1:NLMP_INPUT_FIX));                               !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Initialization of arrays related to fix properties #########################################
                                                                                         !
    if ( NLMP_INPUT_FIX > 0 ) then;                                                      !
                                                                                         !
        LMP_INPUT_RUN_STYLE(1:NLMP_INPUT_FIX) = 'XXX';                                   !
                                                                                         !
        LMP_INPUT_RUN_TAUNH(1:NLMP_INPUT_FIX) = 0.0d0;                                   !
                                                                                         !
        LMP_INPUT_RUN_GRPNAME(1:NLMP_INPUT_FIX) = 'XXX';                                 !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Allocate arrays related to frozen entities in the simulation box ###########################
                                                                                         !
    if ( NLMP_INPUT_FROZEN > 0 ) then;                                                   !
                                                                                         !
        allocate(LMP_INPUT_FROZEN_GRPNAME(1:NLMP_INPUT_FROZEN));                         !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Initialization of arrays related to frozen entities in the simulation box ##################
                                                                                         !
    if ( NLMP_INPUT_FROZEN > 0 ) then;                                                   !
                                                                                         !
        LMP_INPUT_FROZEN_GRPNAME(1:NLMP_INPUT_FROZEN) = 'XXX';                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Allocate arrays related to run properties ##################################################
                                                                                         !
    if ( NLMP_INPUT_RUN > 0 ) then;                                                      !
                                                                                         !
        allocate(LMP_INPUT_RUN_TEMPK_INIT(1:NLMP_INPUT_RUN));                            !
                                                                                         !
        allocate(LMP_INPUT_RUN_TEMPK_FINAL(1:NLMP_INPUT_RUN));                           !
                                                                                         !
        allocate(LMP_INPUT_RUN_NSTEPS(1:NLMP_INPUT_RUN));                                !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Initialization of arrays related to run properties #########################################
                                                                                         !
    if ( NLMP_INPUT_RUN > 0 ) then;                                                      !
                                                                                         !
        LMP_INPUT_RUN_TEMPK_INIT(1:NLMP_INPUT_RUN) = 0.0d0;                              !
                                                                                         !
        LMP_INPUT_RUN_TEMPK_FINAL(1:NLMP_INPUT_RUN) = 0.0d0;                             !
                                                                                         !
        LMP_INPUT_RUN_NSTEPS(1:NLMP_INPUT_RUN) = 0;                                      !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Read properties of lammps run set in the lmp_input command #################################
                                                                                         !
    if ( ( NLMP_INPUT_FIX    > 0 ) .OR. &                                                !
         ( NLMP_INPUT_FROZEN > 0 ) .OR. &                                                !
         ( NLMP_INPUT_RUN    > 0 ) ) then;                                               !
                                                                                         !
        ILOCAL_FIX = 0;                                                                  !
                                                                                         !
        ILOCAL_FROZEN = 0;                                                               !
                                                                                         !
        ILOCAL_RUN = 0;                                                                  !
                                                                                         !
        open(1,file='input_100.dat',status='old');                                       !
                                                                                         !
        EOF = 0;                                                                         !
                                                                                         !
        do                                                                               !
                                                                                         !
            read(1,'(a)',iostat=EOF) CHARLINE;                                           !
                                                                                         !
            if ( EOF /= 0 ) EXIT;                                                        !
                                                                                         !
            if ( INDEX(CHARLINE,'#') > 0 ) CYCLE;                                        !
                                                                                         !
            if ( INDEX(CHARLINE,'classic ') > 0 ) EXIT;                                  !
                                                                                         !
            IOSEF2 = INDEX(CHARLINE,TRIM(SLAB_INPUT_KEYWORDS(6))//' ');                  ! The keyword correspond to the lmp_input command 
                                                                                         !
            if ( IOSEF2 > 20 ) CYCLE;                                                    !
                                                                                         !
            if ( IOSEF2 == 0 ) CYCLE;                                                    !
                                                                                         !
            read(CHARLINE,*) CHOSEF1, CHOSEF2;                                           !
                                                                                         !
            if ( TRIM(CHOSEF2) /= 'dynamics ' ) CYCLE;                                   !
                                                                                         !
            IOSEF10 = LEN_TRIM(CHARLINE);                                                !
                                                                                         !
!           ### Read properties after the 3rd keyword ##############################################
                                                                                         !
            read(CHARLINE,*) CHOSEF1, CHOSEF2, CHOSEF3;                                  !
                                                                                         !
            select case(TRIM(CHOSEF3));                                                  !
                                                                                         !
                case('run');                                                             !
                                                                                         !
                    ILOCAL_RUN = ILOCAL_RUN + 1;                                         !
                                                                                         !
                    IOSEF3 = INDEX(CHARLINE,' tempk ');                                  !
                                                                                         !
                    if ( IOSEF3 > 0 ) then;                                              !
                                                                                         !
                        IOSEF3 = IOSEF3 + 6;                                             !
                                                                                         !
                        read(CHARLINE(IOSEF3:IOSEF10),*)                  &              !
                                   LMP_INPUT_RUN_TEMPK_INIT(ILOCAL_RUN),  &              ! 
                                   LMP_INPUT_RUN_TEMPK_FINAL(ILOCAL_RUN);                !
                                                                                         !
                    end if                                                               !
                                                                                         !
                    IOSEF3 = INDEX(CHARLINE,' nsteps ');                                 !
                                                                                         !
                    if ( IOSEF3 > 0 ) then;                                              !
                                                                                         !
                        IOSEF3 = IOSEF3 + 7;                                             !
                                                                                         !
                        read(CHARLINE(IOSEF3:IOSEF10),*)               &                 !
                                     LMP_INPUT_RUN_NSTEPS(ILOCAL_RUN);                   !
                                                                                         !
                    end if                                                               !
                                                                                         !
!                   stop; !//////////////////////////////////////////////////////////////!
                                                                                         !
            end select                                                                   !
                                                                                         !
!           ### Read properties after the 4th keyword ##############################################
                                                                                         !
            read(CHARLINE,*) CHOSEF1, CHOSEF2, CHOSEF3, CHOSEF4;                         !
                                                                                         !
            select case(TRIM(CHOSEF4));                                                  !
                                                                                         !
                case('nvt');                                                             !
                                                                                         !
                    ILOCAL_FIX = ILOCAL_FIX + 1;                                         !
                                                                                         !
                    read(CHARLINE,*) CHOSEF1,                               &            !
                                     CHOSEF2,                               &            !
                                     LMP_INPUT_RUN_GRPNAME(ILOCAL_FIX),     &            !
                                     LMP_INPUT_RUN_STYLE(ILOCAL_FIX),       &            !
                                     LMP_INPUT_RUN_TAUNH(ILOCAL_FIX);                    !
                                                                                         !
!                   stop; !//////////////////////////////////////////////////////////////!
                                                                                         !
                                                                                         !
                case('frozen');                                                          !
                                                                                         !
                    ILOCAL_FROZEN = ILOCAL_FROZEN + 1;                                   !
                                                                                         !
                    read(CHARLINE,*) CHOSEF1,                               &            !
                                     CHOSEF2,                               &            !
                                     LMP_INPUT_FROZEN_GRPNAME(ILOCAL_FROZEN);            ! 
                                                                                         !
                                                                                         !
!               case default;                                                            !
                                                                                         !
!                   write(icanal,'(a70)') '| Not implemented - stop '// &                !
!                                         REPEAT(' ',44)//'|';                           !
                                                                                         !
!                   call CLOSING_PROGRAM(icanal,1);                                      !
                                                                                         !
            end select                                                                   !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end do                                                                           !
                                                                                         !
        close(1);                                                                        !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Write properties of lammps fixes as defined with the lmp_input command #####################
                                                                                         !
    if ( NLMP_INPUT_FIX > 0 ) then;                                                      !
                                                                                         !
        do i = 1, NLMP_INPUT_FIX;                                                        !
                                                                                         !
            IOSEF1 = 70 - 14 - 4 - 11 - 1 - LEN_TRIM(LMP_INPUT_RUN_STYLE(i));            !
                                                                                         !
            write(icanal,'(a14,i4,a52)') '| Lammps fix #',              &                !
                                         i,                             &                !
                                         ' has style '//                &                !
                                         TRIM(LMP_INPUT_RUN_STYLE(i))// &                !
                                         REPEAT(' ',IOSEF1)//'|';                        !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            write(icanal,'(a16,f12.5,a42)') '| Tau_NH [fs] : ',     &                    !
                                            LMP_INPUT_RUN_TAUNH(i), &                    !
                                            REPEAT(' ',41)//'|';                         !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write properties of lammps runs as defined with the lmp_input command ######################
                                                                                         !
    if ( NLMP_INPUT_RUN > 0 ) then;                                                      !
                                                                                         !
        do i = 1, NLMP_INPUT_RUN;                                                        !
                                                                                         !
            IOSEF1 = 70 - 14 - 4 - 1;                                                    !
                                                                                         !
            write(icanal,'(a14,i4,a52)') '| Lammps run #',              &                !
                                         i,                             &                !
                                         REPEAT(' ',IOSEF1)//'|';                        !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            write(icanal,'(a16,f12.5,a42)') '| T_init  [K] : ',          &               !
                                            LMP_INPUT_RUN_TEMPK_INIT(i), &               !
                                            REPEAT(' ',41)//'|';                         !
                                                                                         ! 
            write(icanal,'(a16,f12.5,a42)') '| Tfinal  [K] : ',           &              !
                                            LMP_INPUT_RUN_TEMPK_FINAL(i), &              !
                                            REPEAT(' ',41)//'|';                         !
                                                                                         !
            write(icanal,'(a16,i8,a46)') '| Nsteps      : ',       &                     !
                                         LMP_INPUT_RUN_NSTEPS(i),  &                     !
                                         REPEAT(' ',45)//'|';                            !
                                                                                         ! 
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialization of the variable for counting the number of defined groups ################### 
                                                                                         !
    NLMP_INPUT_GROUP = 0;                                                                !
                                                                                         !
!   ### Count the number of lammps groups set in the lmp_input command #############################
                                                                                         !
    open(1,file='input_100.dat',status='old');                                           !
                                                                                         !
    EOF = 0;                                                                             !
                                                                                         !
    do                                                                                   !
                                                                                         !
        read(1,'(a)',iostat=EOF) CHARLINE;                                               !
                                                                                         !
        if ( EOF /= 0 ) EXIT;                                                            !
                                                                                         !
        if ( INDEX(CHARLINE,'#') > 0 ) CYCLE;                                            !
                                                                                         !
        if ( INDEX(CHARLINE,'classic ') > 0 ) EXIT;                                      !
                                                                                         !
        IOSEF2 = INDEX(CHARLINE,TRIM(SLAB_INPUT_KEYWORDS(6))//' ');                      ! The keyword correspond to the lmp_input command 
                                                                                         !
        if ( IOSEF2 > 20 ) CYCLE;                                                        !
                                                                                         !
        if ( IOSEF2 == 0 ) CYCLE;                                                        !
                                                                                         !
        read(CHARLINE,*) CHOSEF1, CHOSEF2;                                               !
                                                                                         !
        if ( TRIM(CHOSEF2) /= 'group ' ) CYCLE;                                          !
                                                                                         !
        NLMP_INPUT_GROUP = NLMP_INPUT_GROUP + 1;                                         !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
    end do                                                                               !
                                                                                         !
    close(1);                                                                            !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         ! 
!   ### Write the number of lammps run set in the lmp_input command ################################
                                                                                         !
    write(icanal,'(a45,i4,a21)') '| The number of lammps groups to be set is : ', &      !
                                 NLMP_INPUT_GROUP,                                &      !
                                 REPEAT(' ',20)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate arrays for the definition of groups in the lammps simulation ######################
                                                                                         !
    if ( NLMP_INPUT_GROUP > 0 ) then;                                                    !
                                                                                         !
        allocate(LMP_INPUT_GROUP_NAME(1:NLMP_INPUT_GROUP));                              !
                                                                                         !
        allocate(LMP_INPUT_GROUP_NLIST(1:NLMP_INPUT_GROUP));                             !
                                                                                         !
        allocate(LMP_INPUT_GROUP_LIST(1:NLAMMPS_INSERTS,1:NLMP_INPUT_GROUP));            !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Initialization of arrays for the definition of groups in the lammps simulation #############
                                                                                         !
    if ( NLMP_INPUT_GROUP > 0 ) then;                                                    !
                                                                                         !
        LMP_INPUT_GROUP_NAME(1:NLMP_INPUT_GROUP) = 'XXX';                                !
                                                                                         !
        LMP_INPUT_GROUP_NLIST(1:NLMP_INPUT_GROUP) = 0;                                   !
                                                                                         !
        LMP_INPUT_GROUP_LIST(1:NLAMMPS_INSERTS,1:NLMP_INPUT_GROUP) = 'XXX';              ! 
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Read group properties in the input file of the current program ############################# 
                                                                                         !
    if ( NLMP_INPUT_GROUP > 0 ) then;                                                    !
                                                                                         !
        ILOCAL_GROUP = 0;                                                                !
                                                                                         !
        open(1,file='input_100.dat',status='old');                                       !
                                                                                         !
        EOF = 0;                                                                         !
                                                                                         !
        do                                                                               !
                                                                                         !
            read(1,'(a)',iostat=EOF) CHARLINE;                                           !
                                                                                         !
            if ( EOF /= 0 ) EXIT;                                                        !
                                                                                         !
            if ( INDEX(CHARLINE,'#') > 0 ) CYCLE;                                        !
                                                                                         !
            if ( INDEX(CHARLINE,'classic ') > 0 ) EXIT;                                  !
                                                                                         !
            IOSEF2 = INDEX(CHARLINE,TRIM(SLAB_INPUT_KEYWORDS(6))//' ');                  ! The keyword correspond to the lmp_input command 
                                                                                         !
            if ( IOSEF2 > 20 ) CYCLE;                                                    !
                                                                                         !
            if ( IOSEF2 == 0 ) CYCLE;                                                    !
                                                                                         !
            read(CHARLINE,*) CHOSEF1, CHOSEF2;                                           !
                                                                                         !
            if ( TRIM(CHOSEF2) /= 'group ' ) CYCLE;                                      !
                                                                                         !
            ILOCAL_GROUP = ILOCAL_GROUP + 1;                                             !
                                                                                         !
            read(CHARLINE,*) CHOSEF1,                             &                      !
                             CHOSEF2,                             &                      !
                             LMP_INPUT_GROUP_NAME(ILOCAL_GROUP),  &                      !
                             LMP_INPUT_GROUP_NLIST(ILOCAL_GROUP);                        ! 
                                                                                         !
            IOSEF5 = LMP_INPUT_GROUP_NLIST(ILOCAL_GROUP);                                !
                                                                                         !
            if ( IOSEF5 == 0 ) then;                                                     !
                                                                                         !
                write(icanal,'(a70)') '| Not implemented - stop '// &                    !
                                      REPEAT(' ',44)//'|';                               !
                                                                                         !
                call CLOSING_PROGRAM(icanal,1);                                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
            read(CHARLINE,*) CHOSEF1,                                     &              !
                             CHOSEF2,                                     &              !
                             CHOSEF3,                                     &              !
                             IOSEF6,                                      &              ! 
                             LMP_INPUT_GROUP_LIST(1:IOSEF5,ILOCAL_GROUP);                !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end do                                                                           !
                                                                                         !
        close(1);                                                                        !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         ! 
    end if                                                                               !
                                                                                         !
!   ### Write properties of groups for lammps inputs ###############################################
                                                                                         !
    if ( NLMP_INPUT_GROUP > 0 ) then;                                                    !
                                                                                         !
        do i = 1, NLMP_INPUT_GROUP;                                                      !
                                                                                         !
            call CONVERT_INT_TO_CHAR(icanal,i,3,CHOSEF1);                                !
                                                                                         !
            CHOSEF2 = '| Lammps group #'//           &                                   !
                      TRIM(CHOSEF1)//                &                                   !
                      ' with name '//                &                                   !
                      TRIM(LMP_INPUT_GROUP_NAME(i));                                     !
                                                                                         !
            IOSEF1 = 70 - 1 - 1 - LEN_TRIM(CHOSEF2);                                     !
                                                                                         !
            CHOSEF2 = TRIM(CHOSEF2)//' '//REPEAT('^',IOSEF1)//'|';                       ! 
                                                                                         !
            write(icanal,'(a70)') TRIM(CHOSEF2);                                         !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            do j = 1, LMP_INPUT_GROUP_NLIST(i);                                          !
                                                                                         !
                call CONVERT_INT_TO_CHAR(icanal,j,3,CHOSEF1);                            !
                                                                                         !
                CHOSEF2 = '| --> '//TRIM(CHOSEF1)//' '//TRIM(LMP_INPUT_GROUP_LIST(j,i)); !
                                                                                         !
                IOSEF1 = 70 - 1 - LEN_TRIM(CHOSEF2);                                     !
                                                                                         !
                CHOSEF2 = TRIM(CHOSEF2)//REPEAT(' ',IOSEF1)//'|';                        !  
                                                                                         !
                write(icanal,'(a70)') TRIM(CHOSEF2);                                     !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        end do                                                                           !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         ! 
    end if                                                                               !
                                                                                         !
!   ### Initialization of the variable for counting the number of group unions #####################
                                                                                         !
    NLMP_INPUT_GROUP_UNION = 0;                                                          !
                                                                                         !
!   ### Count the number of lammps group unions set in the lmp_input command #######################
                                                                                         !
    open(1,file='input_100.dat',status='old');                                           !
                                                                                         !
    EOF = 0;                                                                             !
                                                                                         !
    do                                                                                   !
                                                                                         !
        read(1,'(a)',iostat=EOF) CHARLINE;                                               !
                                                                                         !
        if ( EOF /= 0 ) EXIT;                                                            !
                                                                                         !
        if ( INDEX(CHARLINE,'#') > 0 ) CYCLE;                                            !
                                                                                         !
        if ( INDEX(CHARLINE,'classic ') > 0 ) EXIT;                                      !
                                                                                         !
        IOSEF2 = INDEX(CHARLINE,TRIM(SLAB_INPUT_KEYWORDS(6))//' ');                      ! The keyword correspond to the lmp_input command 
                                                                                         !
        if ( IOSEF2 > 20 ) CYCLE;                                                        !
                                                                                         !
        if ( IOSEF2 == 0 ) CYCLE;                                                        !
                                                                                         !
        read(CHARLINE,*) CHOSEF1, CHOSEF2;                                               !
                                                                                         !
        if ( TRIM(CHOSEF2) /= 'group_union ' ) CYCLE;                                    !
                                                                                         !
        NLMP_INPUT_GROUP_UNION = NLMP_INPUT_GROUP_UNION + 1;                             !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
    end do                                                                               !
                                                                                         !
    close(1);                                                                            !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         ! 
!   ### Write the number of lammps group unions set in the lmp_input command #######################
                                                                                         !
    write(icanal,'(a51,i4,a15)') '| The number of lammps group '// &                     !
                                 'unions to be set is : ',         &                     !
                                 NLMP_INPUT_GROUP_UNION,           &                     !
                                 REPEAT(' ',14)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate arrays for the definition of group unions in the lammps simulation ################
                                                                                         !
    if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
        allocate(LMP_INPUT_GROUP_UNION_NAME(1:NLMP_INPUT_GROUP_UNION));                  !
                                                                                         !
        allocate(LMP_INPUT_GROUP_UNION_NLIST(1:NLMP_INPUT_GROUP_UNION));                 !
                                                                                         !
        allocate(LMP_INPUT_GROUP_UNION_LIST(1:NLAMMPS_INSERTS,          &                !
                                            1:NLMP_INPUT_GROUP_UNION));                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Initialization of arrays for the definition of group unions in the lammps simulation #######
                                                                                         !
    if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
        LMP_INPUT_GROUP_UNION_NAME(1:NLMP_INPUT_GROUP_UNION) = 'XXX';                    !
                                                                                         !
        LMP_INPUT_GROUP_UNION_NLIST(1:NLMP_INPUT_GROUP_UNION) = 0;                       !
                                                                                         !
        LMP_INPUT_GROUP_UNION_LIST(1:NLAMMPS_INSERTS,1:NLMP_INPUT_GROUP_UNION) = 'XXX';  ! 
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Read group union properties in the input file of the current program #######################
                                                                                         !
    if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
        ILOCAL_GROUP_UNION = 0;                                                          !
                                                                                         !
        open(1,file='input_100.dat',status='old');                                       !
                                                                                         !
        EOF = 0;                                                                         !
                                                                                         !
        do                                                                               !
                                                                                         !
            read(1,'(a)',iostat=EOF) CHARLINE;                                           !
                                                                                         !
            if ( EOF /= 0 ) EXIT;                                                        !
                                                                                         !
            if ( INDEX(CHARLINE,'#') > 0 ) CYCLE;                                        !
                                                                                         !
            if ( INDEX(CHARLINE,'classic ') > 0 ) EXIT;                                  !
                                                                                         !
            IOSEF2 = INDEX(CHARLINE,TRIM(SLAB_INPUT_KEYWORDS(6))//' ');                  ! The keyword correspond to the lmp_input command 
                                                                                         !
            if ( IOSEF2 > 20 ) CYCLE;                                                    !
                                                                                         !
            if ( IOSEF2 == 0 ) CYCLE;                                                    !
                                                                                         !
            read(CHARLINE,*) CHOSEF1, CHOSEF2;                                           !
                                                                                         !
            if ( TRIM(CHOSEF2) /= 'group_union ' ) CYCLE;                                !
                                                                                         !
            ILOCAL_GROUP_UNION = ILOCAL_GROUP_UNION + 1;                                 !
                                                                                         !
            read(CHARLINE,*) CHOSEF1,                                         &          !
                             CHOSEF2,                                         &          !
                             LMP_INPUT_GROUP_UNION_NAME(ILOCAL_GROUP_UNION),  &          !
                             LMP_INPUT_GROUP_UNION_NLIST(ILOCAL_GROUP_UNION);            ! 
                                                                                         !
            IOSEF5 = LMP_INPUT_GROUP_UNION_NLIST(ILOCAL_GROUP_UNION);                    !
                                                                                         !
            if ( IOSEF5 == 0 ) then;                                                     !
                                                                                         !
                write(icanal,'(a70)') '| Not implemented - stop '// &                    !
                                      REPEAT(' ',44)//'|';                               !
                                                                                         !
                call CLOSING_PROGRAM(icanal,1);                                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
            read(CHARLINE,*) CHOSEF1,                                                 &  !
                             CHOSEF2,                                                 &  !
                             CHOSEF3,                                                 &  !
                             IOSEF6,                                                  &  ! 
                             LMP_INPUT_GROUP_UNION_LIST(1:IOSEF5,ILOCAL_GROUP_UNION);    !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end do                                                                           !
                                                                                         !
        close(1);                                                                        !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         ! 
    end if                                                                               !
                                                                                         !
!   ### Write properties of group unions for lammps inputs #########################################
                                                                                         !
    if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
        do i = 1, NLMP_INPUT_GROUP_UNION;                                                !
                                                                                         !
            call CONVERT_INT_TO_CHAR(icanal,i,3,CHOSEF1);                                !
                                                                                         !
            CHOSEF2 = '| Lammps group union #'//           &                             !
                      TRIM(CHOSEF1)//                      &                             !
                      ' with name '//                      &                             !
                      TRIM(LMP_INPUT_GROUP_UNION_NAME(i));                               !
                                                                                         !
            IOSEF1 = 70 - 1 - 1 - LEN_TRIM(CHOSEF2);                                     !
                                                                                         !
            CHOSEF2 = TRIM(CHOSEF2)//' '//REPEAT('^',IOSEF1)//'|';                       ! 
                                                                                         !
            write(icanal,'(a70)') TRIM(CHOSEF2);                                         !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            do j = 1, LMP_INPUT_GROUP_UNION_NLIST(i);                                    !
                                                                                         !
                call CONVERT_INT_TO_CHAR(icanal,j,3,CHOSEF1);                            !
                                                                                         !
                CHOSEF2 = '| --> '//                             &                       !
                          TRIM(CHOSEF1)//' '//                   &                       !
                          TRIM(LMP_INPUT_GROUP_UNION_LIST(j,i));                         !
                                                                                         !
                IOSEF1 = 70 - 1 - LEN_TRIM(CHOSEF2);                                     !
                                                                                         !
                CHOSEF2 = TRIM(CHOSEF2)//REPEAT(' ',IOSEF1)//'|';                        !  
                                                                                         !
                write(icanal,'(a70)') TRIM(CHOSEF2);                                     !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        end do                                                                           !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         ! 
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine SEEK_FOR_LMP_INPUT_COMMAND











