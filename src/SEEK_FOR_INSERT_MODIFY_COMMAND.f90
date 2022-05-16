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


subroutine SEEK_FOR_INSERT_MODIFY_COMMAND(icanal)

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

    use module_regions;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

!   ************************************************************************************************

    integer (kind=4) :: h, i, j;

    integer (kind=4) :: ilocal_pair_style, ilocal_tether;

    integer (kind=4) :: ILOCAL_INSERTS, IMAX_MODIFIED;

!   ************************************************************************************************

    character (len=250) :: CHARLINE;

    character (len=250) :: LOCAL_INSERT_NAME;

!   ************************************************************************************************

    integer (kind=4), allocatable, dimension(:) :: NCOUNT;

!   ************************************************************************************************

    integer (kind=4) :: EOF;

!   ************************************************************************************************
                                                                                         !
!   ### Allocate arrays related to the number of times insert arrays are modified ##################
                                                                                         !
    allocate(LAMMPS_INSERT_NMODIFIED(1:NLAMMPS_INSERTS));                                !
                                                                                         !
!   ### Initialization of arrays related to the number of times inseert arrays are modified ########
                                                                                         !
    LAMMPS_INSERT_NMODIFIED(1:NLAMMPS_INSERTS) = 0;                                      !
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
        if ( INDEX(CHARLINE,'classic') > 0 ) EXIT;                                       !
                                                                                         !
        IOSEF2 = INDEX(CHARLINE,TRIM(SLAB_INPUT_KEYWORDS(5))//' ');                      ! The third keyword correspond to the insert_modify command 
                                                                                         !
        if ( IOSEF2 > 20 ) CYCLE;                                                        !
                                                                                         !
        if ( IOSEF2 == 0 ) CYCLE;                                                        !
                                                                                         !
        read(CHARLINE,*) CHOSEF1, LOCAL_INSERT_NAME;                                     !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!       ### Look for the specified insert command in insert arrays #################################
                                                                                         !
        ILOCAL_INSERTS = 0;                                                              !
                                                                                         !
        do i = 1, NLAMMPS_INSERTS;                                                       ! 
                                                                                         !
            if ( TRIM(LAMMPS_INSERT_NAME(i)) /= TRIM(LOCAL_INSERT_NAME) ) CYCLE;         !
                                                                                         ! 
            ILOCAL_INSERTS = i;                                                          ! 
                                                                                         !
            EXIT;                                                                        !
                                                                                         !
        end do                                                                           !
                                                                                         !
        LAMMPS_INSERT_NMODIFIED(ILOCAL_INSERTS) = &                                      !
        LAMMPS_INSERT_NMODIFIED(ILOCAL_INSERTS) + 1;                                     !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
    end do                                                                               !
                                                                                         !
    close(1);                                                                            !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialization of the variable containing the maximum number of modifications ##############
                                                                                         !
    IMAX_MODIFIED = 0;                                                                   !
                                                                                         !
!   ### Write the number of times each insert command is modified ##################################
                                                                                         !
    do i = 1, NLAMMPS_INSERTS;                                                           ! 
                                                                                         !
        if ( LAMMPS_INSERT_NMODIFIED(i) == 0 ) CYCLE;                                    !
                                                                                         !
        if ( LAMMPS_INSERT_NMODIFIED(i) > IMAX_MODIFIED ) then;                          !
                                                                                         !
            IMAX_MODIFIED = LAMMPS_INSERT_NMODIFIED(i);                                  !
                                                                                         ! 
        end if                                                                           !
                                                                                         !
        CHOSEF1 = '| The insert command name '// &                                       !
                  TRIM(LAMMPS_INSERT_NAME(i))//  &                                       !
                  ' is modified';                                                        !
                                                                                         !
        IOSEF1 = LEN_TRIM(CHOSEF1) + 1;                                                  !
                                                                                         !
        call CONVERT_INT_TO_CHAR(icanal,IOSEF1,3,CHOSEF2);                               !
                                                                                         !
        IOSEF2 = 70 - IOSEF1 - 4;                                                        !
                                                                                         !
        call CONVERT_INT_TO_CHAR(icanal,IOSEF2,3,CHOSEF3);                               !
                                                                                         !
        IOSEF2 = IOSEF2 - 6 - 1;                                                         !
                                                                                         !
        CHOSEF4 = 'a'//TRIM(CHOSEF2)//',i4,a'//TRIM(CHOSEF3);                            !
                                                                                         !
        write(icanal,'('//TRIM(CHOSEF4)//')') TRIM(CHOSEF1)//' ',                &       !
                                              LAMMPS_INSERT_NMODIFIED(i),        &       !
                                              ' times'//REPEAT(' ',IOSEF2)//'|';         !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the maximuum number of modifications #################################################
                                                                                         !
    write(icanal,'(a41,i4,a25)') '| The maximum number of modifications is ', &          !
                                 IMAX_MODIFIED,                               &          !
                                 REPEAT(' ',24)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate arrays related to the style of the insert modify command ##########################
                                                                                         !
    allocate(LAMMPS_INSERT_MODIFY_STYLE(1:IMAX_MODIFIED,1:NLAMMPS_INSERTS));             !
                                                                                         !
!   ### Initialization of arrays related to the style of the insert modify command #################
                                                                                         ! 
    LAMMPS_INSERT_MODIFY_STYLE(1:IMAX_MODIFIED,1:NLAMMPS_INSERTS) = 'XXX';               !
                                                                                         !
!   ### Allocate local array for counting the number of modifications ##############################
                                                                                         !
    allocate(NCOUNT(1:NLAMMPS_INSERTS));                                                 !
                                                                                         ! 
!   ### Initialization of the local array for counting the number of modifications #################
                                                                                         !
    NCOUNT(1:NLAMMPS_INSERTS) = 0;                                                       !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialization of variables and arrays related to the insert_modify command ################
                                                                                         !
    LAMMPS_INSERT_IPOTENTIAL_CLASS2 = 0;                                                 !
                                                                                         !
    LAMMPS_INSERT_ITETHER = 0;                                                           ! 
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Read properties related to the insert command ##############################################
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
        if ( INDEX(CHARLINE,'classic') > 0 ) EXIT;                                       !
                                                                                         !
        IOSEF2 = INDEX(CHARLINE,TRIM(SLAB_INPUT_KEYWORDS(5))//' ');                      ! The third keyword correspond to the insert_modify command 
                                                                                         !
        if ( IOSEF2 > 20 ) CYCLE;                                                        !
                                                                                         !
        if ( IOSEF2 == 0 ) CYCLE;                                                        !
                                                                                         !
        read(CHARLINE,*) CHOSEF1, LOCAL_INSERT_NAME;                                     !
                                                                                         !
!       ### Write the name of the insert command that will be modified #############################
                                                                                         !
        IOSEF1 = 70 - 54 - 1 - 1 - LEN_TRIM(LOCAL_INSERT_NAME);                          !
                                                                                         !
        if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                    !
                                                                                         !
        write(icanal,'(a70)') '| The insert_modify command '// &                         !
                              'will change properties of '//   &                         !
                              TRIM(LOCAL_INSERT_NAME)//        &                         !
                              ' '//REPEAT('>',IOSEF1)//'|';                              !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!       ### Look for the specified insert command in insert arrays #################################
                                                                                         !
        ILOCAL_INSERTS = 0;                                                              !
                                                                                         !
        do i = 1, NLAMMPS_INSERTS;                                                       ! 
                                                                                         !
            if ( TRIM(LAMMPS_INSERT_NAME(i)) /= TRIM(LOCAL_INSERT_NAME) ) CYCLE;         !
                                                                                         ! 
            ILOCAL_INSERTS = i;                                                          ! 
                                                                                         !
            EXIT;                                                                        !
                                                                                         !
        end do                                                                           !
                                                                                         !
        write(icanal,'(a19,i4,a47)') '| ILOCAL_INSERTS : ', &                            !
                                     ILOCAL_INSERTS,        &                            !
                                     REPEAT(' ',46)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!       ### Update the array counting the number of modifications ##################################
                                                                                         !
        NCOUNT(ILOCAL_INSERTS) = NCOUNT(ILOCAL_INSERTS) + 1;                             !  
                                                                                         !
!       ### Get the total length of the current line ###############################################
                                                                                         !
        IOSEF10 = LEN_TRIM(CHARLINE);                                                    !
                                                                                         !
!       ### Read the modification type to perform on the insert command ############################
                                                                                         !
        read(CHARLINE,*) CHOSEF1, CHOSEF2, CHOSEF3;                                      !
                                                                                         !
!       ### Set the style of the current insert_modify command #####################################
                                                                                         !
        IOSEF1 = NCOUNT(ILOCAL_INSERTS);                                                 !
                                                                                         !
        LAMMPS_INSERT_MODIFY_STYLE(IOSEF1,ILOCAL_INSERTS) = TRIM(CHOSEF3);               !
                                                                                         !
!       ### Get modification properties depending on the modification type #########################
                                                                                         !
        select case(TRIM(CHOSEF3));                                                      !
                                                                                         !
            case('pair_coeff');                                                          !
                                                                                         !
!               ### Allocate arrays when the first time the pair_coeff style is encountered ########
                                                                                         !
                if ( LAMMPS_INSERT_IPOTENTIAL_CLASS2 == 0 ) then;                        !
                                                                                         !
!                   ### Allocate arrays for the insert_modify command ##############################
                                                                                         !
                    allocate(LAMMPS_INSERT_POTENTIAL_CLASS2_CHTYPE(1:NLAMMPS_INSERTS));  !
                                                                                         !
                    allocate(LAMMPS_INSERT_POTENTIAL_CLASS2_NAT(1:2,                 &   !
                                                                1:IMAX_MODIFIED,     &   !
                                                                1:NLAMMPS_INSERTS));     !
                                                                                         !
                    allocate(LAMMPS_INSERT_POTENTIAL_CLASS2(1:2,                 &       !
                                                            1:IMAX_MODIFIED,     &       !
                                                            1:NLAMMPS_INSERTS));         !
                                                                                         !
!                   ### Initialization of arrays for the insert_modify command #####################
                                                                                         !
                    LAMMPS_INSERT_POTENTIAL_CLASS2_CHTYPE(1:NLAMMPS_INSERTS) = 'XXX';    !
                                                                                         !
                    LAMMPS_INSERT_POTENTIAL_CLASS2_NAT(1:2,              &               !
                                                       1:IMAX_MODIFIED,  &               !
                                                       1:NLAMMPS_INSERTS) = 'XXX';       !
                                                                                         !
                    LAMMPS_INSERT_POTENTIAL_CLASS2(1:2,                 &                !
                                                   1:IMAX_MODIFIED,     &                !
                                                   1:NLAMMPS_INSERTS) = 0.0d0;           ! 
                                                                                         !
!                   ### Set to 1 the flag saying that pair potential parameters were modified ######
                                                                                         !
                    LAMMPS_INSERT_IPOTENTIAL_CLASS2 = 1;                                 !
                                                                                         !
                end if                                                                   !
                                                                                         !
!               ### Read the style of pair_style parameters to modify ##############################
                                                                                         !
                ilocal_pair_style = 0;                                                   !
                                                                                         !
                if ( INDEX(CHARLINE,' lj/cut/coul/long ') > 0 ) then;                    !
                                                                                         !
                    LAMMPS_INSERT_POTENTIAL_CLASS2_CHTYPE(ILOCAL_INSERTS) = &            !
                    'lj/cut/coul/long';                                                  !
                                                                                         !
                    ilocal_pair_style = 1;                                               !
                                                                                         !
                else if ( INDEX(CHARLINE,' lj/charmmfsw/coul/long ') > 0 ) then;         !
                                                                                         !
                    LAMMPS_INSERT_POTENTIAL_CLASS2_CHTYPE(ILOCAL_INSERTS) = &            !
                    'lj/charmmfsw/coul/long';                                            !     
                                                                                         !
                    ilocal_pair_style = 2;                                               !
                                                                                         !
                else                                                                     !
                                                                                         !
                    write(icanal,'(a70)') '| Not implemented - stop '// &                !
                                          REPEAT(' ',44)//'|';                           !
                                                                                         !
                    call CLOSING_PROGRAM(icanal,1);                                      !
                                                                                         !
                end if                                                                   !
                                                                                         ! 
!               ### Write the style of pair_style parameters to modify #############################
                                                                                         !
                IOSEF1 = 70 - 46 - 1 - &                                                 !
                         LEN_TRIM(LAMMPS_INSERT_POTENTIAL_CLASS2_CHTYPE(ILOCAL_INSERTS));!       
                                                                                         !
                if ( IOSEF1 < 0 ) IOSEF1 = 1;                                            !
                                                                                         !
                write(icanal,'(a70)') '| The pair coeffs to modify have the style is '// &   !
                                      TRIM(LAMMPS_INSERT_POTENTIAL_CLASS2_CHTYPE(ILOCAL_INSERTS))// &
                                      REPEAT(' ',IOSEF1)//'|';                           !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
!               ### Read modified potential parameters #############################################
                                                                                         !
                IOSEF1 = NCOUNT(ILOCAL_INSERTS);                                         !
                                                                                         !
                if ( ilocal_pair_style == 1 ) then;                                      !
                                                                                         !
                    IOSEF3 = INDEX(CHARLINE,' lj/cut/coul/long ') + 19;                  !
                                                                                         !
                else if ( ilocal_pair_style == 2 ) then;                                 !
                                                                                         !
                    IOSEF3 = INDEX(CHARLINE,' lj/charmmfsw/coul/long ') + 25;            !
                                                                                         !
                else                                                                     !
                                                                                         !
                    write(icanal,'(a70)') '| Not implemented - stop '// &                !
                                          REPEAT(' ',44)//'|';                           !
                                                                                         !
                    call CLOSING_PROGRAM(icanal,1);                                      !
                                                                                         !
                end if                                                                   !
                                                                                         !
                select case(TRIM(LAMMPS_INSERT_CHFILE_FORMAT(ILOCAL_INSERTS)));          !
                                                                                         !
                    case('xyz');                                                         !
                                                                                         !
                        read(CHARLINE(IOSEF3:IOSEF10),*)                                  &  !
                           LAMMPS_INSERT_POTENTIAL_CLASS2_NAT(1:2,IOSEF1,ILOCAL_INSERTS), &  !
                           LAMMPS_INSERT_POTENTIAL_CLASS2(1:2,IOSEF1,ILOCAL_INSERTS);        !
                                                                                         !
                    case default;                                                        !
                                                                                         !
                        write(icanal,'(a70)') '| Not implemented - stop '// &            !
                                              REPEAT(' ',44)//'|';                       !
                                                                                         !
                        call CLOSING_PROGRAM(icanal,1);                                  !
                                                                                         !
                end select                                                               !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            case('tether');                                                              !
                                                                                         !
!               ### Allocate arrays the the tether style is encountered for the first time #########
                                                                                         !
                if ( LAMMPS_INSERT_ITETHER == 0 ) then;                                  !
                                                                                         !
!                   ### Allocate arrays related to the tether style ################################
                                                                                         !
                    allocate(LAMMPS_INSERT_TETHER_CHTYPE(1:NLAMMPS_INSERTS));            !
                                                                                         !
                    allocate(LAMMPS_INSERT_TETHER_COEFFS(1:2,1:NLAMMPS_INSERTS));        !
                                                                                         !  
!                   ### Initialization of arrays related to the tether style #######################
                                                                                         !
                    LAMMPS_INSERT_TETHER_CHTYPE(1:NLAMMPS_INSERTS) = 'XXX';              !
                                                                                         !
                    LAMMPS_INSERT_TETHER_COEFFS(1:2,1:NLAMMPS_INSERTS) = 0.0d0;          !
                                                                                         !
!                   ### Set to 1 the flag saying that tethered atoms are employed ##################
                                                                                         !
                    LAMMPS_INSERT_ITETHER = 1;                                           !
                                                                                         !
                end if                                                                   !
                                                                                         !
!               ### Read the style of tether parameters to read ####################################
                                                                                         !
                ilocal_tether = 0;                                                       !
                                                                                         !
                read(CHARLINE,*) CHOSEF1,                                         &      !
                                 CHOSEF2,                                         &      !
                                 CHOSEF3,                                         &      !
                                 LAMMPS_INSERT_TETHER_CHTYPE(ILOCAL_INSERTS),     &      !
                                 LAMMPS_INSERT_TETHER_COEFFS(1:2,ILOCAL_INSERTS);        !
                                                                                         !
            case default;                                                                !
                                                                                         !
                write(icanal,'(a70)') '| Not implemented - stop '// &                    !
                                      REPEAT(' ',44)//'|';                               !
                                                                                         !
                call CLOSING_PROGRAM(icanal,1);                                          !
                                                                                         !
        end select                                                                       !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
    end do                                                                               !
                                                                                         !
    close(1);                                                                            !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write potential modification to insert commands ############################################
                                                                                         !
    do i = 1, NLAMMPS_INSERTS;                                                           !
                                                                                         !
        if ( LAMMPS_INSERT_NMODIFIED(i) == 0 ) CYCLE;                                    !
                                                                                         !
        CHOSEF1 = '| '//TRIM(LAMMPS_INSERT_NAME(i));                                     !
                                                                                         !
        IOSEF1 = 70 - 1 - LEN_TRIM(CHOSEF1);                                             !
                                                                                         !
        if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                    !
                                                                                         !
        write(icanal,'(a70)') TRIM(CHOSEF1)//REPEAT(' ',IOSEF1)//'|';                    !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       ### Write modificactions related to pair potentials ########################################
                                                                                         !
        do j = 1, LAMMPS_INSERT_NMODIFIED(i);                                            !
                                                                                         !
            if ( LAMMPS_INSERT_MODIFY_STYLE(j,i) /= 'pair_coeff' ) CYCLE;                !
                                                                                         !
            CHOSEF1 = '| '//                                            &                !
                      TRIM(LAMMPS_INSERT_MODIFY_STYLE(j,i))//' : '//    &                !
                      TRIM(LAMMPS_INSERT_POTENTIAL_CLASS2_NAT(1,j,i))// &                ! 
                      REPEAT(' ',4)//                                   &                !
                      TRIM(LAMMPS_INSERT_POTENTIAL_CLASS2_NAT(2,j,i));                   !
                                                                                         !
            IOSEF1 = LEN_TRIM(CHOSEF1);                                                  !
                                                                                         !
            call CONVERT_INT_TO_CHAR(icanal,IOSEF1,3,CHOSEF2);                           !
                                                                                         !
            IOSEF2 = 70 - IOSEF1 - 2 * 12;                                               !
                                                                                         !
            call CONVERT_INT_TO_CHAR(icanal,IOSEF2,3,CHOSEF3);                           !
                                                                                         !
            IOSEF2 = IOSEF2 - 1;                                                         !
                                                                                         !
            CHOSEF4 = 'a'//TRIM(CHOSEF2)//',2f12.6,a'//TRIM(CHOSEF3);                    !
                                                                                         !
            write(icanal,'('//TRIM(CHOSEF4)//')')                             &          !
                                     TRIM(CHOSEF1),                           &          !
                                     LAMMPS_INSERT_POTENTIAL_CLASS2(1:2,j,i), &          !
                                     REPEAT(' ',IOSEF2)//'|';                            !
                                                                                         !
        end do                                                                           !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       ### Write modifications related to tethered atoms ##########################################
                                                                                         !
        do j = 1, LAMMPS_INSERT_NMODIFIED(i);                                            !
                                                                                         !
            if ( LAMMPS_INSERT_MODIFY_STYLE(j,i) /= 'tether' ) CYCLE;                    !
                                                                                         !
            CHOSEF1 = '| '//                                   &                         !
                      TRIM(LAMMPS_INSERT_MODIFY_STYLE(j,i))//  &                         !
                      ' with type '//                          &                         !
                      TRIM(LAMMPS_INSERT_TETHER_CHTYPE(i))//   &                         !
                      ' : ';                                                             !
                                                                                         !
            IOSEF1 = LEN_TRIM(CHOSEF1);                                                  !
                                                                                         !
            call CONVERT_INT_TO_CHAR(icanal,IOSEF1,3,CHOSEF2);                           !
                                                                                         !
            IOSEF2 = 70 - IOSEF1 - 2 * 12;                                               !
                                                                                         !
            call CONVERT_INT_TO_CHAR(icanal,IOSEF2,3,CHOSEF3);                           !
                                                                                         !
            IOSEF2 = IOSEF2 - 1;                                                         !
                                                                                         !
            CHOSEF4 = 'a'//TRIM(CHOSEF2)//',2f12.6,a'//TRIM(CHOSEF3);                    !
                                                                                         !
            write(icanal,'('//TRIM(CHOSEF4)//')')                             &          !
                                     TRIM(CHOSEF1),                           &          !
                                     LAMMPS_INSERT_TETHER_COEFFS(1:2,i),      &          !
                                     REPEAT(' ',IOSEF2)//'|';                            !
                                                                                         !
        end do                                                                           !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   ### Deallocate local arrays ####################################################################
                                                                                         !
    deallocate(NCOUNT);                                                                  !
                                                                                         !
!   ### Write message for the success in reanding the input file ###################################
                                                                                         !
    write(icanal,'(a70)') '| Insert commands were read '//REPEAT(' ',41)//'|';           !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine SEEK_FOR_INSERT_MODIFY_COMMAND

