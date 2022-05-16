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

subroutine SEEK_FOR_INSERT_COMMAND(icanal)

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

!   ************************************************************************************************

    character (len=250) :: CHARLINE;

!   ************************************************************************************************

    integer (kind=4) :: EOF;

!   ************************************************************************************************
                                                                                         !
!   ### Initialization of variables and arrays related to the insert command #######################
                                                                                         !
    NLAMMPS_INSERTS = 0;                                                                 ! 
                                                                                         !
!   ### Open the input file ########################################################################
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
        IOSEF2 = INDEX(CHARLINE,TRIM(SLAB_INPUT_KEYWORDS(3))//' ');                      !
                                                                                         !
        if ( IOSEF2 > 20 ) CYCLE;                                                        !
                                                                                         !
        if ( IOSEF2 == 0 ) CYCLE;                                                        !
                                                                                         !
        NLAMMPS_INSERTS = NLAMMPS_INSERTS + 1;                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
    close(1);                                                                            !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the number of regions found in the input file ########################################
                                                                                         !
    write(icanal,'(a30,i4,a36)') '| Number of defined inserts : ', &                     !
                                 NLAMMPS_INSERTS,                  &                     !
                                 REPEAT(' ',35)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### Allocate arrays for the reading of the insert command ######################################
                                                                                         !
    allocate(LAMMPS_INSERT_NAME(1:NLAMMPS_INSERTS));                                     !
                                                                                         !
    allocate(LAMMPS_INSERT_FILE_OPTION(1:NLAMMPS_INSERTS));                              !
                                                                                         !
    allocate(LAMMPS_INSERT_LIB_NAME(1:NLAMMPS_INSERTS));                                 !
                                                                                         !
    allocate(LAMMPS_INSERT_CHFILE_FORMAT(1:NLAMMPS_INSERTS));                            !
                                                                                         !
    allocate(LAMMPS_INSERT_CHNAME_FILE(1:NLAMMPS_INSERTS));                              !
                                                                                         !
    allocate(LAMMPS_INSERT_METHOD(1:NLAMMPS_INSERTS));                                   !
                                                                                         !
    allocate(LAMMPS_INSERT_COM_RG(1:3,1:NLAMMPS_INSERTS));                               !
                                                                                         !
    allocate(LAMMPS_INSERT_FLAG(1:NLAMMPS_INSERTS));                                     !
                                                                                         !
    allocate(LAMMPS_INSERT_REGION_NAME(1:NLAMMPS_INSERTS));                              !
                                                                                         !
    allocate(LAMMPS_INSERT_RANDOM_SEED(1:NLAMMPS_INSERTS));                              !
                                                                                         !
    allocate(LAMMPS_INSERT_REGION_IRANK(1:NLAMMPS_INSERTS));                             !
                                                                                         !
    allocate(LAMMPS_INSERT_NMOLECULE(1:NLAMMPS_INSERTS));                                !
                                                                                         !
    allocate(LAMMPS_INSERT_NTYPE_ATOM_MAX(1:NLAMMPS_INSERTS));                           !
                                                                                         !
!   ### Initialization of arrays for the reading of the insert command #############################
                                                                                         !
    LAMMPS_INSERT_NAME(1:NLAMMPS_INSERTS) = 'XXX';                                       !
                                                                                         !
    LAMMPS_INSERT_FILE_OPTION(1:NLAMMPS_INSERTS) = 'XXX';                                !
                                                                                         !
    LAMMPS_INSERT_LIB_NAME(1:NLAMMPS_INSERTS) = 'XXX';                                   !
                                                                                         !
    LAMMPS_INSERT_CHFILE_FORMAT(1:NLAMMPS_INSERTS) = 'XXX';                              !
                                                                                         !
    LAMMPS_INSERT_CHNAME_FILE(1:NLAMMPS_INSERTS) = 'XXX';                                !
                                                                                         !
    LAMMPS_INSERT_METHOD(1:NLAMMPS_INSERTS) = 'XXX';                                     !
                                                                                         !
    LAMMPS_INSERT_COM_RG(1:3,1:NLAMMPS_INSERTS) = 0.0d0;                                 !
                                                                                         !
    LAMMPS_INSERT_FLAG(1:NLAMMPS_INSERTS) = 0;                                           !
                                                                                         !
    LAMMPS_INSERT_REGION_NAME(1:NLAMMPS_INSERTS) = 'XXX';                                !
                                                                                         !
    LAMMPS_INSERT_RANDOM_SEED(1:NLAMMPS_INSERTS) = 0;                                    !
                                                                                         !
    LAMMPS_INSERT_REGION_IRANK(1:NLAMMPS_INSERTS) = 0;                                   !
                                                                                         !
    LAMMPS_INSERT_NMOLECULE(1:NLAMMPS_INSERTS) = 0;                                      !
                                                                                         !
    LAMMPS_INSERT_NTYPE_ATOM_MAX(1:NLAMMPS_INSERTS) = 0;                                 !
                                                                                         !
!   ### Read properties related to the insert command ##############################################
                                                                                         !
    open(1,file='input_100.dat',status='old');                                           !
                                                                                         !
    IOSEF1 = 0;                                                                          !
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
        IOSEF2 = INDEX(CHARLINE,TRIM(SLAB_INPUT_KEYWORDS(3))//' ');                      ! The third keyword correspond to the insert command 
                                                                                         !
        if ( IOSEF2 > 20 ) CYCLE;                                                        !
                                                                                         !
        if ( IOSEF2 == 0 ) CYCLE;                                                        !
                                                                                         !
        IOSEF1 = IOSEF1 + 1;                                                             !
                                                                                         !
!       write(icanal,*) 'IOSEF1 : ', IOSEF1;                                             !
                                                                                         !
        read(CHARLINE,*) CHOSEF1, LAMMPS_INSERT_NAME(IOSEF1);                            !
                                                                                         !
!       ### Get the total length of the current line ###############################################
                                                                                         !
        IOSEF10 = LEN_TRIM(CHARLINE);                                                    !
                                                                                         !
!       ### Read the option for reading the input molecular configuration ##########################
                                                                                         !
        if ( INDEX(CHARLINE,' lib ') > 0 ) then;                                         !
                                                                                         !
            LAMMPS_INSERT_FILE_OPTION(IOSEF1) = 'lib';                                   !
                                                                                         !
            IBUILD_METHOD = 1;                                                           !
                                                                                         !
            NMOLECULAR_CONFIG = NMOLECULAR_CONFIG + 1;                                   !
                                                                                         !
            LAMMPS_INSERT_FLAG(IOSEF1) = NMOLECULAR_CONFIG;                              ! 
                                                                                         !
!           ### Read the name of the defined library where the file is located #####################
                                                                                         !
            IOSEF3 = INDEX(CHARLINE,' lib ') + 5;                                        !
                                                                                         !
            read(CHARLINE(IOSEF3:IOSEF10),*) LAMMPS_INSERT_LIB_NAME(IOSEF1);             !
                                                                                         !
!           ### Read the format of the library file ################################################
                                                                                         !
!           IOSEF3 = INDEX(CHARLINE,' lib ') + 5;                                        !
            IOSEF3 = INDEX(CHARLINE,TRIM(LAMMPS_INSERT_LIB_NAME(IOSEF1))) + &            !
                     LEN_TRIM(LAMMPS_INSERT_LIB_NAME(IOSEF1)) + 1;                       !
                                                                                         !
            read(CHARLINE(IOSEF3:IOSEF10),*) LAMMPS_INSERT_CHFILE_FORMAT(IOSEF1);        !
                                                                                         !
!           ### Read the name of the library file ##################################################
                                                                                         !
            read(CHARLINE(IOSEF3:IOSEF10),*) CHOSEF1, LAMMPS_INSERT_CHNAME_FILE(IOSEF1); ! 
                                                                                         !
        end if                                                                           !
                                                                                         ! 
!       ### Read the inseertion method #############################################################
                                                                                         !
        if ( INDEX(CHARLINE,' com ') > 0 ) then;                                         !
                                                                                         !
            LAMMPS_INSERT_METHOD(IOSEF1) = 'com';                                        !
                                                                                         !
!           ### Read the location of the center of mass ############################################
                                                                                         !
            IOSEF3 = INDEX(CHARLINE,' com ') + 5;                                        !
                                                                                         !
            read(CHARLINE(IOSEF3:IOSEF10),*) LAMMPS_INSERT_COM_RG(1:3,IOSEF1);           !
                                                                                         !
            LAMMPS_INSERT_NMOLECULE(IOSEF1) = 1;                                         !
                                                                                         !
        else if ( INDEX(CHARLINE,' rand/region ') > 0 ) then;                            !
                                                                                         !
            LAMMPS_INSERT_METHOD(IOSEF1) = 'rand/region';                                !
                                                                                         !
!           ### Read properties of the region to consider for the random insertion #################
                                                                                         !
            IOSEF3 = INDEX(CHARLINE,' rand/region ') + 13;                               !
                                                                                         !
            read(CHARLINE(IOSEF3:IOSEF10),*) LAMMPS_INSERT_REGION_NAME(IOSEF1), &        !
                                             LAMMPS_INSERT_RANDOM_SEED(IOSEF1), &        !
                                             LAMMPS_INSERT_NMOLECULE(IOSEF1);            !        
                                                                                         !
            do i = 1, NLAMMPS_REGIONS;                                                   !
                                                                                         !
                if ( TRIM(LAMMPS_INSERT_REGION_NAME(IOSEF1)) /= &                        !
                          LAMMPS_REGION_NAME(i) ) CYCLE;                                 !
                                                                                         !
                LAMMPS_INSERT_REGION_IRANK(IOSEF1) = i;                                  !
                                                                                         !
            end do                                                                       !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end if                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
    close(1);                                                                            !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the name of insert commands ##########################################################
                                                                                         !
    do i = 1, NLAMMPS_INSERTS;                                                           !
                                                                                         !
        IOSEF1 = 70 - 22 - 4 - 14 - 1 - 1 - LEN_TRIM(LAMMPS_INSERT_NAME(i));             !
                                                                                         !
        if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                    !
                                                                                         !
        write(icanal,'(a22,i4,a44)') '| The insert command #',     &                     !
                                     i,                            &                     !
                                     ' has the name '//            &                     !
                                     TRIM(LAMMPS_INSERT_NAME(i))// &                     !
                                     ' '//REPEAT('>',IOSEF1)//'|';                       !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        if ( TRIM(LAMMPS_INSERT_FILE_OPTION(i)) == 'lib' ) then;                         !
                                                                                         !
            write(icanal,'(a70)') '| The structure is taken from a library'// &          !
                                  REPEAT(' ',30)//'|';                                   ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF2 = 70 - 22 - 1 - LEN_TRIM(LAMMPS_INSERT_LIB_NAME(i));                  !
                                                                                         !
            write(icanal,'(a70)') '| The library name is '//        &                    !
                                  TRIM(LAMMPS_INSERT_LIB_NAME(i))// &                    !
                                  REPEAT(' ',IOSEF2)//'|';                               !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            IOSEF2 = 70 - 23 - 1 - LEN_TRIM(LAMMPS_INSERT_CHFILE_FORMAT(i));             !
                                                                                         !
            if ( IOSEF2 < 0 ) IOSEF2 = 1;                                                !
                                                                                         !
            write(icanal,'(a70)') '| The file format is : '//            &               !
                                  TRIM(LAMMPS_INSERT_CHFILE_FORMAT(i))// &               !
                                  REPEAT(' ',IOSEF2)//'|';                               !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           ### Write the name of the library file #################################################
                                                                                         !
            IOSEF2 = 70 - 21 - 1 - LEN_TRIM(LAMMPS_INSERT_CHNAME_FILE(i));               !
                                                                                         !
            if ( IOSEF2 < 0 ) IOSEF2 = 1;                                                !
                                                                                         !
            write(icanal,'(a70)') '| The file name is : '//            &                 !
                                  TRIM(LAMMPS_INSERT_CHNAME_FILE(i))// &                 !
                                  REPEAT(' ',IOSEF2)//'|';                               !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           ### Write the insertion method #########################################################
                                                                                         !
            if ( TRIM(LAMMPS_INSERT_METHOD(i)) == 'com' ) then;                          ! 
                                                                                         !
                write(icanal,'(a70)') '| The insertion method will be '// &              !
                                      'done through the center of mass'// &              !
                                      REPEAT(' ',7)//'|';                                !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                write(icanal,'(a11,3f12.4,a23)') '| RG [A] : ',               &          !
                                                 LAMMPS_INSERT_COM_RG(1:3,i), &          !
                                                 REPEAT(' ',22)//'|';                    !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            else if ( TRIM(LAMMPS_INSERT_METHOD(i)) == 'rand/region' ) then;             !
                                                                                         !
                write(icanal,'(a70)') '| The insertion method will be '// &              !
                                      'done through random insertion'// &                !
                                      REPEAT(' ',9)//'|';                                !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                IOSEF1 = 70 - 12 - 1 - LEN_TRIM(LAMMPS_INSERT_REGION_NAME(i));           !
                                                                                         !
                write(icanal,'(a70)') '| Region is '//                     &             !
                                      TRIM(LAMMPS_INSERT_REGION_NAME(i))// &             !
                                      REPEAT(' ',IOSEF1)//'|';                           !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                write(icanal,'(a17,i4,a49)') '| Region rank is ',           &            !
                                             LAMMPS_INSERT_REGION_IRANK(i), &            !
                                             REPEAT(' ',48)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                write(icanal,'(a26,i8,a36)') '| Random generator seed : ', &             !
                                             LAMMPS_INSERT_RANDOM_SEED(i), &             !
                                             REPEAT(' ',35)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                write(icanal,'(a39,i8,a23)') '| The number of molecules to insert is ', &!
                                             LAMMPS_INSERT_NMOLECULE(i),                &!
                                             REPEAT(' ',22)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            else                                                                         !
                                                                                         ! 
!               write(icanal,*) 'Not implemented - stop!';                               !
                                                                                         !
!               close(icanal);                                                           !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
                write(icanal,'(a70)') '| Not implemented - stop '// &                    !
                                      REPEAT(' ',44)//'|';                               !
                                                                                         !
                call CLOSING_PROGRAM(icanal,1);                                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
        end if                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   ### Write message for the success in reanding the input file ###################################
                                                                                         !
    write(icanal,'(a70)') '| Insert commands were read '//REPEAT(' ',41)//'|';           !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine SEEK_FOR_INSERT_COMMAND











