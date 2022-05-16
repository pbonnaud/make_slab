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

subroutine READ_INPUT(icanal) !,PRIMITIVE_CELL_AXIS,PRIMITIVE_CELL_ANGDEG)

!   ************************************************************************************************
!   **                                  READ THE INPUT FILE                                       **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                : CANAL ON WHICH OUTPUT DATA ARE WRITTEN                             **
!   ** PRIMITIVE_CELL_AXIS   : DIMENSIONS OF THE PRIMITIVE CELL                                   **
!   ** PRIMITIVE_CELL_ANGDEG : ANGLES OF THE PRIMITIVE CELL                                       **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

    use module_duplicate;

    use module_library;

    use module_physical_constants;

    use module_config;

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

!   ************************************************************************************************

!   real (kind=8), dimension(1:3), intent(out) :: PRIMITIVE_CELL_AXIS, PRIMITIVE_CELL_ANGDEG;

!   ************************************************************************************************

    integer (kind=4) :: h, i, j;

!   ************************************************************************************************

    integer (kind=4) :: ilocal_classic, idefined_centers;

    integer (kind=4) :: ipore, iii, NMOLECULE_INSERTION_MAX;

    logical :: PROBE1;

    real (kind=8) :: CSTE1, CSTE2, CSTE3;

    character (len=200) :: CHCHAIN1;

    character (len=250) :: CHAIN_LENGTH;

    character (len=500) :: CHERROR_MESSAGE;

!   ************************************************************************************************

    character (len=250) :: CHARLINE;

!   ************************************************************************************************

    integer (kind=4) :: EOF, EOF2;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the current routine #####################################################
                                                                                         !
    CHTITLE = 'Read the input file';                                                     !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   ### Check the presence of the input file #######################################################
                                                                                         !
    inquire(FILE='input_100.dat',EXIST=PROBE1);                                          !
                                                                                         !
    if ( PROBE1 .EQV. .FALSE. ) then;                                                    !
                                                                                         !
        CHERROR_MESSAGE = 'No input file - Please check your directory';                 !
                                                                                         !
        call WRITE_ERROR_MESSAGE(icanal,70,CHERROR_MESSAGE);                             !
                                                                                         !
    else                                                                                 !
                                                                                         !
        write(icanal,'(a70)') '| The input file has been found in the directory'// &     !
                              REPEAT(' ',21)//'|';                                       !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialization of variables and arrays for reading the input file ##########################
                                                                                         !
    IBUILD_METHOD = 0;                                                                   !
                                                                                         !
    NMOLECULAR_CONFIG = 0;                                                               !
                                                                                         !
    NREGIONS_INSERTION = 0;                                                              !
                                                                                         !
    NGenerate_polymer_species = 0;                                                       !
                                                                                         !
!   ### Set the keyword list of the input file #####################################################
                                                                                         !
    call BUILD_INPUT_KEYWORD_LIST(icanal);                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Seek for library paths #####################################################################
                                                                                         !
    call SEEK_FOR_SET_LIB_COMMAND(icanal);                                               !
                                                                                         !
!   stop; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                                         !
!   ### Seek for the definitions of regions in the input file ######################################
                                                                                         !
    call SEEK_FOR_REGIONS(icanal);                                                       !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Seek for the create_box command ############################################################
                                                                                         !
    call SEEK_FOR_CREATE_BOX_COMMAND(icanal);                                            !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Seek for the insert command ################################################################ 
                                                                                         !
    call SEEK_FOR_INSERT_COMMAND(icanal);                                                ! 
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Seek for the insert modify command #########################################################
                                                                                         !
    call SEEK_FOR_INSERT_MODIFY_COMMAND(icanal);                                         ! 
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Seek for the lmp_input command #############################################################
                                                                                         !
    call SEEK_FOR_LMP_INPUT_COMMAND(icanal);                                             !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialization of variables related to the input file ######################################
                                                                                         !
    iduplicate = 0;                                                                      !
                                                                                         !
    iregions = 0;                                                                        !
                                                                                         !
    iworking_file_insertion = 0;                                                         !
                                                                                         !
!   ### Open the input file ########################################################################
                                                                                         !
    open(1,file='input_100.dat',status='old');                                           !
                                                                                         !
!   ### Look for keywords ##########################################################################
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
        if ( INDEX(CHARLINE,'duplicate') > 0 ) then;                                     !
                                                                                         !
            iduplicate = 1;                                                              !
                                                                                         !
            iduplicate_ffield = 0;                                                       !
                                                                                         !
            IBUILD_METHOD = 1;                                                           !
                                                                                         !
            write(icanal,'(a70)') '| The program will duplicate a simulation box'// &    !
                                  REPEAT(' ',24)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           ### Get the length of the current line #################################################
                                                                                         !
            IOSEF10 = LEN_TRIM(CHARLINE);                                                !
                                                                                         !
!           ### Write the length of the current character vairable #################################
                                                                                         !
            write(icanal,'(a28,i4,a38)') '| The character length is : ', &               !
                                         IOSEF10,                        &               !
                                         REPEAT(' ',37)//'|';                            !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           ### check if a library was defined for the molecular structure to duplicate ############
                                                                                         !
            if ( INDEX(CHARLINE,' lib ') > 0 ) then;                                     !
                                                                                         !
                write(icanal,'(a70)') '| The structure is taken from a library'// &      !
                                      REPEAT(' ',30)//'|';                               !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               ### Set to 1 the flag related to the building method (library) ##################### 
                                                                                         !
                IBUILD_METHOD = 1;                                                       !
                                                                                         !
!               ### Set to 1 the number of molecular configs to generate ###########################
                                                                                         !
                NMOLECULAR_CONFIG = 1;                                                   !
                                                                                         !
!               ### Read the format of the library file ############################################
                                                                                         !
                IOSEF1 = INDEX(CHARLINE,' lib ') + 5;                                    !
                                                                                         !
                read(CHARLINE(IOSEF1:IOSEF10),*) DUPLICATE_CHFILE_FORMAT;                !
                                                                                         !
                IOSEF2 = 70 - 23 - 1 - LEN_TRIM(DUPLICATE_CHFILE_FORMAT);                !
                                                                                         !
                if ( IOSEF2 < 0 ) IOSEF2 = 1;                                            !
                                                                                         !
                write(icanal,'(a70)') '| The file format is : '//     &                  !
                                      TRIM(DUPLICATE_CHFILE_FORMAT)// &                  !
                                      REPEAT(' ',IOSEF2)//'|';                           !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               ### Read the name of the library file ##############################################
                                                                                         !
                read(CHARLINE(IOSEF1:IOSEF10),*) CHOSEF1, DUPLICATE_CHNAME_FILE;         ! 
                                                                                         !
!               ### Write the name of the library file #############################################
                                                                                         !
                IOSEF2 = 70 - 21 - 1 - LEN_TRIM(DUPLICATE_CHNAME_FILE);                  !
                                                                                         !
                if ( IOSEF2 < 0 ) IOSEF2 = 1;                                            !
                                                                                         !
                write(icanal,'(a70)') '| The file name is : '//     &                    !
                                      TRIM(DUPLICATE_CHNAME_FILE)// &                    !
                                      REPEAT(' ',IOSEF2)//'|';                           !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
!           ### Initialization of duplication numbers in the 3 directions of space #################
                                                                                         !
            DUPLICATE_NX = 1;                                                            !
                                                                                         !
            DUPLICATE_NY = 1;                                                            !
                                                                                         !
            DUPLICATE_NZ = 1;                                                            !
                                                                                         !
!           ### Read duplication over the x-direction ############################################## 
                                                                                         !
            if ( INDEX(CHARLINE,' nx ') > 0 ) then;                                      !
                                                                                         !
                IOSEF1 = INDEX(CHARLINE,' nx ') + 4;                                     !
                                                                                         !
                read(CHARLINE(IOSEF1:IOSEF10),*) DUPLICATE_NX;                           !
                                                                                         !
                write(icanal,'(a35,i4,a31)') '| The structure will be duplicated ', &    !
                                             DUPLICATE_NX,                          &    !
                                             ' times over the x-direction'//        &    !
                                             REPEAT(' ',3)//'|';                         !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
!           ### Read duplication over the y-direction ############################################## 
                                                                                         !
            if ( INDEX(CHARLINE,' ny ') > 0 ) then;                                      !
                                                                                         !
                IOSEF1 = INDEX(CHARLINE,' ny ') + 4;                                     !
                                                                                         !
                read(CHARLINE(IOSEF1:IOSEF10),*) DUPLICATE_NY;                           !
                                                                                         !
                write(icanal,'(a35,i4,a31)') '| The structure will be duplicated ', &    !
                                             DUPLICATE_NY,                          &    !
                                             ' times over the y-direction'//        &    !
                                             REPEAT(' ',3)//'|';                         !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
!           ### Read duplication over the z-direction ############################################## 
                                                                                         !
            if ( INDEX(CHARLINE,' nz ') > 0 ) then;                                      !
                                                                                         !
                IOSEF1 = INDEX(CHARLINE,' nz ') + 4;                                     !
                                                                                         !
                read(CHARLINE(IOSEF1:IOSEF10),*) DUPLICATE_NZ;                           !
                                                                                         !
                write(icanal,'(a35,i4,a31)') '| The structure will be duplicated ', &    !
                                             DUPLICATE_NZ,                          &    !
                                             ' times over the z-direction'//        &    !
                                             REPEAT(' ',3)//'|';                         !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
!           ### Read force-field properties ########################################################
                                                                                         !
            if ( INDEX(CHARLINE,' ffield ') > 0 ) then;                                  !
                                                                                         !
                write(icanal,'(a70)') '| A force field was defined '// &                 !
                                      REPEAT(' ',41)//'|';                               !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               ### Set to 1 the flag for assigning a force field ##################################
                                                                                         !
                iduplicate_ffield = 1;                                                   !
                                                                                         !
!               ### Set to 1 the number of library file to read ####################################
                                                                                         !
                NFILE_LIBRARY = 1;                                                       !
                                                                                         !
!               ### Allocate arrays for insertion files ############################################
                                                                                         !
                if ( NFILE_LIBRARY > 0 ) then;                                           !
                                                                                         !
                    allocate(NMOLECULE_INSERTION(1:NFILE_LIBRARY));                      !
                                                                                         !
                    allocate(CHNAME_FILE_LIBRARY(1:3,1:NFILE_LIBRARY));                  !
                                                                                         !
                    allocate(INSERTION_METHOD_LIBRARY(1:NFILE_LIBRARY));                 !
                                                                                         !
                end if                                                                   !
                                                                                         !
!               ### Read the file format for the force field #######################################
                                                                                         !
                IOSEF1 = INDEX(CHARLINE,' ffield ') + 8;                                 !
                                                                                         !
                read(CHARLINE(IOSEF1:IOSEF10),*) CHFILE_FORMAT;                          !
                                                                                         !
                IOSEF2 = 70 - 43 - 1 - LEN_TRIM(CHFILE_FORMAT);                          !
                                                                                         !
                if ( IOSEF2 < 0 ) IOSEF2 = 1;                                            !
                                                                                         !
                write(icanal,'(a70)') '| The file format for the force field is : '// &  !
                                      TRIM(CHFILE_FORMAT)//                           &  !
                                      REPEAT(' ',IOSEF2)//'|';                           !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               ### Read the name of the file for the force field ##################################
                                                                                         !
                read(CHARLINE(IOSEF1:IOSEF10),*) CHOSEF1, CHNAME_FILE_LIBRARY(1,1);      ! 
                                                                                         !
!               ### Write the name of the file for the force field #################################
                                                                                         !
                IOSEF2 = 70 - 41 - 1 - LEN_TRIM(CHNAME_FILE_LIBRARY(1,1));               !
                                                                                         !
                if ( IOSEF2 < 0 ) IOSEF2 = 1;                                            !
                                                                                         !
                write(icanal,'(a70)') '| The file name for the force field is : '// &    !
                                      TRIM(CHNAME_FILE_LIBRARY(1,1))//              &    !
                                      REPEAT(' ',IOSEF2)//'|';                           !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
        end if                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Read the building method ###################################################################
                                                                                         !
!   ilocal_classic = 0;                                                                  !

!   if ( iduplicate == 0 ) then;
!   if ( ilocal_classic == 1 ) then;                                                     !
                                                                                         !
!   read(1,*) IBUILD_METHOD;                                                             ! SET THE METHOD TO BUILD THE INITIAL MOLECULAR STRUCTURE
                                                                                         !
!   if ( IBUILD_METHOD == 1 ) then;                                                      !
                                                                                         !
!       write(icanal,'(a70)') '| The molecular structure will be built '// &             !
!                             'with molecules coming from a |';                          !
                                                                                         !
!       write(icanal,'(a70)') '| library'//REPEAT(' ',60)//'|';                          !
                                                                                         !
!   else if ( IBUILD_METHOD == 2 ) then;                                                 !
                                                                                         !
!       write(icanal,'(a70)') '| The molecular structure will be built with '// &        !
!                             'molecules coming from   |';                               !
                                                                                         !
!       write(icanal,'(a70)') '| the working directory'//REPEAT(' ',46)//'|';            !
                                                                                         !
!   else if ( IBUILD_METHOD == 3 ) then;                                                 !
                                                                                         !
!       write(icanal,'(a70)') '| The molecular structure will be built from a '// &      !
!                             'configuration file'//REPEAT(' ',4)//'|';                  !
                                                                                         !
!   else                                                                                 !
                                                                                         !
!       write(icanal,'(a70)') '| Not implemented '//REPEAT(' ',51)//'|';                 !
                                                                                         !
!       write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                  !
                                                                                         !
!       write(icanal,*);                                                                 !
                                                                                         !
!       write(icanal,'(a26)') '!!! END OF THE PROGRAM !!!';                              !
                                                                                         !
!       close(icanal);                                                                   !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Read properties of the simulation box ######################################################
                                                                                         !
!   read(1,*) PRIMITIVE_CELL_AXIS(1:3);                                                  ! PROPERTIES OF THE SIMULATION BOX IF CREATED
                                                                                         !
!   read(1,*) PRIMITIVE_CELL_ANGDEG(1:3);                                                !
                                                                                         !
!   ### Write properties of the simulation box #####################################################
                                                                                         !
!   write(icanal,'(a33,3f10.2,a7)') '| PRIMITIVE_CELL_AXIS      [A] : ', &               !
!                                   PRIMITIVE_CELL_AXIS(1:3),            &               !
!                                   REPEAT(' ',6)//'|';                                  !
                                                                                         !
!   write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   write(icanal,'(a33,3f10.2,a7)') '| PRIMITIVE_CELL_ANGDEG [deg.] : ', &               !
!                                   PRIMITIVE_CELL_ANGDEG(1:3),          &               ! 
!                                   REPEAT(' ',6)//'|';                                  !
                                                                                         !
!   write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Read the number of molecular configurations to generate ####################################
                                                                                         !
!   read(1,*) NMOLECULAR_CONFIG;                                                         !
                                                                                         !
!   ### Write the number of molecular configurations to generate ###################################
                                                                                         !
!   write(icanal,'(a50,i8,a12)') '| Number of molecular configuration to generate : ', & !
!                                NMOLECULAR_CONFIG,                                    & !
!                                REPEAT(' ',11)//'|';                                    !
                                                                                         !
!   write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Read the number of files to read for insertion and the library name where they should be read ####
                                                                                         !
!   read(1,*) NFILE_LIBRARY, CHFILE_FORMAT, iuse_moleculeid;                             ! NUMBER OF MOLECULE TYPE TO INSERT AND FILE FORMAT
                                                                                         !
!   ### Allocate arrays for insertion files ########################################################
                                                                                         !
!   if ( NFILE_LIBRARY > 0 ) then;                                                       !
                                                                                         !
!       allocate(NMOLECULE_INSERTION(1:NFILE_LIBRARY));                                  !
                                                                                         !
!       allocate(CHNAME_FILE_LIBRARY(1:3,1:NFILE_LIBRARY));                              !
                                                                                         !
!       allocate(INSERTION_METHOD_LIBRARY(1:NFILE_LIBRARY));                             !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   Read properties of files to insert #############################################################
                                                                                         !
!   idefined_centers = 0;                                                                !
                                                                                         !
!   NMOLECULE_INSERTION_MAX = 0;                                                         !
                                                                                         !
!   if ( NFILE_LIBRARY > 0 ) then;                                                       !
                                                                                         !
!       do i = 1, NFILE_LIBRARY;                                                         ! LOOP OVER THE NUMBER OF FILES TO READ IN THE LIBRARY
                                                                                         !
!           read(1,*) NMOLECULE_INSERTION(i),      &                                     ! READ THE NUMBER OF MOLECULES TO INSERT
!                     CHNAME_FILE_LIBRARY(1,i),    &                                     ! AND THE NAME OF THE CREESPONDING FILE IN THE LIBRARY
!                     INSERTION_METHOD_LIBRARY(i);                                       !
                                                                                         !
!           IOSEF1 = 0;                                                                  !
                                                                                         !
!           if ( INSERTION_METHOD_LIBRARY(i) == 2 ) IOSEF1 = 1;                          ! Centers of mass of molecules to insert were defined
                                                                                         !
!!!!           if ( INSERTION_METHOD_LIBRARY(i) == 3 ) IOSEF1 = 1;                          ! Geometrical centers were defined
                                                                                         !
!           if ( IOSEF1 == 1 ) then;                                                     !
                                                                                         !
!               idefined_centers = 1;                                                    !
                                                                                         !
!               if ( NMOLECULE_INSERTION(i) >  NMOLECULE_INSERTION_MAX ) then;           !
                                                                                         !
!                   NMOLECULE_INSERTION_MAX = NMOLECULE_INSERTION(i);                    ! 
                                                                                         !
!               end if                                                                   !
                                                                                         !
!               do j = 1, NMOLECULE_INSERTION(i);                                        !
                                                                                         !
!                   read(1,*);                                                           !
                                                                                         !
!               end do                                                                   !
                                                                                         !
!           end if                                                                       !
                                                                                         !
!       end do                                                                           !
                                                                                         ! 
!   end if                                                                               !
                                                                                         !
!   ### Write properties of files to insert ########################################################
                                                                                         !
!   write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   IOSEF1 = 13 - LEN_TRIM(CHFILE_FORMAT);                                               !
                                                                                         !
!   write(icanal,'(a48,i4,4x,a14)')                                     &                !
!                   '| NFILE_LIBRARY        | CHFILE_FORMAT        : ', &                !
!                   NFILE_LIBRARY,                                      &                !
!                   TRIM(CHFILE_FORMAT)//REPEAT(' ',IOSEF1)//'|';                        !
                                                                                         !
!   write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Check the existence of a working file for insertion ########################################
                                                                                         !
!   if ( NFILE_LIBRARY > 0 ) then;                                                       ! 
                                                                                         !
!       read(1,'(a)') CHAIN_LENGTH;                                                      !
                                                                                         !
!       read(CHAIN_LENGTH,*) iworking_file_insertion;                                    ! READ IF INSERTION HAS TO BE DONE IN AN EXISTING SIMULATION BOX
                                                                                         !
!       if ( iworking_file_insertion == 1 ) then;                                        !
                                                                                         !
!           read(CHAIN_LENGTH,*) iworking_file_insertion, CHWORKING_FILE_INSERTION;      !
                                                                                         !
!       end if                                                                           ! 
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Write if a working file has been chosen for insertion ######################################
                                                                                         !
!   if ( NFILE_LIBRARY > 0 ) then;                                                       ! 
                                                                                         !
!       if ( iworking_file_insertion == 1 ) then;                                        !
                                                                                         !
!           IOSEF1 = 38 - LEN_TRIM(CHWORKING_FILE_INSERTION);                            !
                                                                                         !
!           write(icanal,'(a70)') '| WORKING FILE FOR INSERTION : '// &                  !
!                                 TRIM(CHWORKING_FILE_INSERTION)//    &                  !
!                                 REPEAT(' ',IOSEF1)//'|';                               !
                                                                                         !  
!           write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!       else                                                                             !
                                                                                         !
!           write(icanal,'(a70)') '| No working file for insertion'// &                  !
!                                 REPEAT(' ',38)//'|';                                   !
                                                                                         !
!           write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!       end if                                                                           !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!

!   ###

!   if ( NFILE_LIBRARY > 0 ) then; 

!       read(1,'(a)') CHAIN_LENGTH;                                                      !
                                                                                         !
!       read(CHAIN_LENGTH,*) iregions;                                                   ! READ IF REGIONS OF INSERTION WHERE DEFINED IN THE INPUT FILE
                                                                                         !
!       if ( iregions == 1 ) then;                                                       ! IF REGIONS WERE DEFINED IN THE INPUT FILE, THEN
                                                                                         !
!           write(icanal,'(a70)') '| REGIONS OF INSERTION WERE DEFINED '// &             !
!                                 REPEAT(' ',33)//'|';                                   !
                                                                                         !
!           write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           read(CHAIN_LENGTH,*) iregions, NREGIONS_INSERTION;                           ! READ THE NUMBER OF REGIONS DEFINED FOR THE INSERTION OF MOLECULES 
                                                                                         !
!           write(icanal,'(a23,i4,a43)') '| NREGIONS_INSERTION : ', &                    !
!                                        NREGIONS_INSERTION,        &                    !
!                                        REPEAT(' ',42)//'|';                            !
!           write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
!                                                                                        !
!           allocate(REGION_AXIS(1:NREGIONS_INSERTION));                                 ! ALLOCATE ARRAY CONTAINING THE CONSTRAINT DIRECTION FOR A GIVEN REGION
                                                                                         !
!           allocate(REGION_BOUNDS(1:2,1:NREGIONS_INSERTION));                           !
                                                                                         !
!           do j = 1, NREGIONS_INSERTION;                                                ! LOOP OVER THE NUMBER OF REGIONS DEFINED FOR THE INSERTION OF MOLECULES
                                                                                         !
!               read(1,'(a)') CHAIN_LENGTH;                                              !
                                                                                         !
!               read(CHAIN_LENGTH,*) CHOSEF1;                                            !
                                                                                         !
!               if ( ( TRIM(CHOSEF1) == 'x' ) .OR. &                                     !
!                    ( TRIM(CHOSEF1) == 'X' ) ) REGION_AXIS(j) = 1;                      !
                                                                                         !
!               if ( ( TRIM(CHOSEF1) == 'y' ) .OR. &                                     !
!                    ( TRIM(CHOSEF1) == 'Y' ) ) REGION_AXIS(j) = 2;                      !
                                                                                         !
!               if ( ( TRIM(CHOSEF1) == 'z' ) .OR. &                                     !
!                    ( TRIM(CHOSEF1) == 'Z' ) ) REGION_AXIS(j) = 3;                      !
                                                                                         !
!               read(CHAIN_LENGTH,*) CHOSEF1, REGION_BOUNDS(1:2,j);                      !
                                                                                         !
!               write(icanal,'(a10,i2,a21,f8.2,a5,f8.2,a16)')             &              !
!                                                '| REGION #',            &              !
!                                                j,                       &              !
!                                                ' WAS DEFINED BETWEEN ', &              !
!                                                REGION_BOUNDS(1,j),      &              !
!                                                ' AND ',                 &              !
!                                                REGION_BOUNDS(2,j),      &              !
!                                                REPEAT(' ',15)//'|';                    !
!               write(icanal,'(a70)') '| ALONG THE '//TRIM(CHOSEF1)//'-AXIS [A]'// &     !
!                                     REPEAT(' ',47)//'|';                               !
                                                                                         !
!               write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!           end do                                                                       !
                                                                                         !
!       end if                                                                           !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Read the number of files to be read for deletion and the library name ######################
                                                                                         !
!   read(1,*) NFILE_LIBRARY_REMOVE, CHFILE_FORMAT_REMOVE;                                !
                                                                                         !
!   if ( NFILE_LIBRARY_REMOVE > 0 ) then;                                                !
                                                                                         !
!       NMOLECULAR_CONFIG = 1;                                                           !
                                                                                         !
!       IOSEF1 = 70 - 48 - 4 - 4 - 1 - LEN_TRIM(CHFILE_FORMAT_REMOVE);                   !
                                                                                         !
!       if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                    !
                                                                                         !
!       write(icanal,'(a48,i4,4x,a14)')                                  &               !
!                    '| NFILE_LIBRARY_REMOVE | CHFILE_FORMAT_REMOVE : ', &               !
!                    NFILE_LIBRARY_REMOVE,                               &               !
!                    TRIM(CHFILE_FORMAT_REMOVE)//REPEAT(' ',IOSEF1)//'|';                !
                                                                                         !
!       write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       allocate(NMOLECULE_DELETION(1:NFILE_LIBRARY_REMOVE));                            !
                                                                                         !
!       allocate(CHNAME_FILE_LIBRARY_REMOVE(1:NFILE_LIBRARY_REMOVE));                    !
                                                                                         !
!       do i = 1, NFILE_LIBRARY_REMOVE;                                                  !
                                                                                         !
!           read(1,*) NMOLECULE_DELETION(i), CHNAME_FILE_LIBRARY_REMOVE(i);              !  
                                                                                         !
!           IOSEF1 = 44 - LEN_TRIM(CHNAME_FILE_LIBRARY_REMOVE(i));                       !
                                                                                         !
!           write(icanal,'(a11,i4,i8,2x,a45)') '|'//REPEAT(' ',10),                  &   !
!                                              i,                                    &   !
!                                              NMOLECULE_DELETION(i),                &   !
!                                              TRIM(CHNAME_FILE_LIBRARY_REMOVE(i))// &   !
!                                              REPEAT(' ',IOSEF1)//'|';                  !
!       end do                                                                           !
                                                                                         !
!       write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         ! 
!       read(1,*) CHWORKING_FILE_REMOVE;                                                 !
                                                                                         !
!       IOSEF1 = 52 - LEN_TRIM(CHWORKING_FILE_REMOVE);                                   !
                                                                                         !
!       write(icanal,'(a70)') '| WORKING FILE : '//         &                            !
!                             TRIM(CHWORKING_FILE_REMOVE)// &                            !
!                             REPEAT(' ',IOSEF1)//'|';                                   !
                                                                                         !
!       write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Read the number of dummy particles to insert in the simulation box #########################
                                                                                         !
!   read(1,*) NTYPE_DUMMY_PARTICLE;                                                      !
                                                                                         !
!   write(icanal,'(a25,i8,a37)') '| NTYPE_DUMMY_PARTICLE : ', &                          !
!                                NTYPE_DUMMY_PARTICLE,        &                          !
!                                REPEAT(' ',36)//'|';                                    !
                                                                                         !
!   write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate arrays for dymmy particles #####

!   if ( NTYPE_DUMMY_PARTICLE > 0 ) then;                                                !
                                                                                         !
!       allocate(NDUMMY_PARTICLE(1:NTYPE_DUMMY_PARTICLE));                               !
                                                                                         !
!       allocate(CHNAME_DUMMY_PARTICLE(1:NTYPE_DUMMY_PARTICLE));                         !
                                                                                         !
!       allocate(RADIUS_DUMMY_PARTICLE(1:NTYPE_DUMMY_PARTICLE));                         !
                                                                                         !
!       allocate(MASSE_DUMMY_PARTICLE(1:NTYPE_DUMMY_PARTICLE));                          !
                                                                                         !
!       do i = 1, NTYPE_DUMMY_PARTICLE;                                                  !
                                                                                         !
!           read(1,*) NDUMMY_PARTICLE(i),       &                                        !
!                     CHNAME_DUMMY_PARTICLE(i), &                                        !
!                     RADIUS_DUMMY_PARTICLE(i), &                                        !
!                     MASSE_DUMMY_PARTICLE(i);                                           !
                                                                                         !
!       end do                                                                           !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Read properties to generate polymer molecules ##############################################
                                                                                         !  
!   read(1,*) NGenerate_polymer_species, &                                               !
!             Chfile_format_polymer;                                                     !
                                                                                         !
!   if ( NGenerate_polymer_species > 0 ) then;                                           !
                                                                                         !
!       write(icanal,'(a13,i4,a53)') '| There is : ',                 &                  !
!                                    NGenerate_polymer_species,       &                  !
!                                    ' polymer species to generate'// &                  !
!                                    REPEAT(' ',24)//'|';                                !
                                                                                         !
!       write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       IOSEF1 = 70 - 45 - LEN_TRIM(Chfile_format_polymer);                              !
                                                                                         !
!       write(icanal,'(a70)') '| The library where monomers are located is '// &         !
!                             TRIM(Chfile_format_polymer)//REPEAT(' ',IOSEF1)//'|';      !
                                                                                         !
!       write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       allocate(NPOLYMER_CHAINS(1:NGenerate_polymer_species));                          !
                                                                                         !
!       allocate(NSPECIES_PER_POLYMER_CHAINS(1:NGenerate_polymer_species));              !
                                                                                         !
!       allocate(NMONOMER_UNITS_PER_CHAINS(1:NGenerate_polymer_species,1:100));          !
                                                                                         !
!       allocate(MONOMER_INSERTION_PARAM(1:NGenerate_polymer_species,1:5));              !
                                                                                         !
!       IOSEF1 = NGenerate_polymer_species * 100;                                        !
                                                                                         !
!       allocate(Chname_file_library_monomer(1:IOSEF1));                                 !
                                                                                         !
!       NTotal_monomers = 0;                                                             !
                                                                                         !
!       do i = 1, NGenerate_polymer_species;                                             !
                                                                                         !
!           write(icanal,'(a18,i4,a48)') '| Polymer species ', i,   &                    !
!                                        ' '//REPEAT('-',46)//'|';                       !
                                                                                         !
!           write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           read(1,*) NPOLYMER_CHAINS(i), NSPECIES_PER_POLYMER_CHAINS(i);                !
                                                                                         !
!           write(icanal,'(a11,i8,a51)') '| There is ', NPOLYMER_CHAINS(i), &            !
!                                        ' polymer chains to insert'//      &            !
!                                        REPEAT(' ',25)//'|';                            !
                                                                                         !
!           write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           write(icanal,'(a11,i4,a55)') '| There is ',                  &               !
!                                        NSPECIES_PER_POLYMER_CHAINS(i), &               !
!                                        ' monomer species per chains'// &               !
!                                        REPEAT(' ',27)//'|';                            !
                                                                                         !
!           write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           do j = 1, NSPECIES_PER_POLYMER_CHAINS(i);                                    !
                                                                                         !
!               NTotal_monomers = NTotal_monomers + 1;                                   !
                                                                                         ! 
!               read(1,*) NMONOMER_UNITS_PER_CHAINS(i,j),   &                            !
!                         Chname_file_library_monomer(NTotal_monomers);                  !
                                                                                         !
!               IOSEF1 = 51 - 13 - 1 - &                                                 !
!                        LEN_TRIM(Chname_file_library_monomer(NTotal_monomers));         !
                                                                                         !
!               write(icanal,'(a11,i8,a51)') '| There is ',                           &  !
!                                            NMONOMER_UNITS_PER_CHAINS(i,j),          &  !
!                                            ' monomers of '//                        &  !
!                                            TRIM(Chname_file_library_monomer(NTotal_monomers))// &  !
!                                            REPEAT(' ',IOSEF1)//'|';                    !
                                                                                         !
!               write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!           end do                                                                       !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           read(1,*) MONOMER_INSERTION_PARAM(i,1:5);                                    ! 
                                                                                         !
!           write(icanal,'(a25,f8.4,a37)') '| MAX_RAND_THETA [DEG] : ',  &               !
!                                          MONOMER_INSERTION_PARAM(i,1), &               !
!                                          REPEAT(' ',36)//'|';                          !
                                                                                         !
!           write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           write(icanal,'(a23,f8.4,a39)') '| MAX_RAND_PHI [DEG] : ',    &               !
!                                          MONOMER_INSERTION_PARAM(i,2), &               !
!                                          REPEAT(' ',38)//'|';                          !
                                                                                         !
!           write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           write(icanal,'(a23,f8.4,a39)') '| MAX_RAND_PSI [DEG] : ',    &               !
!                                          MONOMER_INSERTION_PARAM(i,3), &               !
!                                          REPEAT(' ',38)//'|';                          !
                                                                                         !
!           write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           write(icanal,'(a33,f12.4,a25)') '| Rcut insertion first monomer : ', &       !
!                                           MONOMER_INSERTION_PARAM(i,4),        &       !
!                                           REPEAT(' ',24)//'|';                         !
                                                                                         !
!           write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           write(icanal,'(a33,f12.4,a25)') '| Rcut insertion other monomer : ', &       !
!                                           MONOMER_INSERTION_PARAM(i,5),        &       !
!                                           REPEAT(' ',24)//'|';                         !
                                                                                         !
!           write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!       end do                                                                           !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Read the target temperature for the simulation #############################################
                                                                                         !
!   read(1,*) TARGET_TEMPK;                                                              !
                                                                                         !
!   write(icanal,'(a17,f10.2,a43)') '| TARGET_TEMPK : ', &                               !
!                                   TARGET_TEMPK,        &                               !
!                                   REPEAT(' ',42)//'|';                                 !
                                                                                         !
!   read(1,*) RANDOM_GENERATOR_SEED;                                                     !
                                                                                         !
!   write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### Read the seed for the random generator #####################################################
                                                                                         !
!   write(icanal,'(a26,i16,a28)') '| RANDOM_GENERATOR_SEED : ', &                        !
!                                 RANDOM_GENERATOR_SEED,        &                        !
!                                 REPEAT(' ',27)//'|';                                   !
                                                                                         !
!   write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### Read properties related to the switch of atom IDs ##########################################
                                                                                         !
!   NSWITCH_ATOMIDS = 0;                                                                 !
                                                                                         !
!   read(1,*) NSWITCH_ATOMIDS, SWITCH_FILE_NAME;                                         !
                                                                                         !
!   IOSEF1 = 70 - 11 - 4 - 23 - LEN_TRIM(SWITCH_FILE_NAME);                              !
                                                                                         !
!   if ( NSWITCH_ATOMIDS > 0 ) then;                                                     !
                                                                                         !
!       write(icanal,'(a11,i4,a55)') '| There is ',             &                        !
!                                    NSWITCH_ATOMIDS,           &                        !
!                                    ' atom IDs to swith in '// &                        !
!                                    TRIM(SWITCH_FILE_NAME)//   &                        ! 
!                                    REPEAT(' ',IOSEF1)//'|';                            !
                                                                                         !
!       write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!!!!!!   if ( NSWITCH_ATOMIDS > 0 ) then;                                                     !
                                                                                         !
!       do i = 1, NSWITCH_ATOMIDS;                                                       !
                                                                                         !
!           read(1,*) SWITCH_ATOMIDS(1:2,i);                                             !
                                                                                         !
!       end do                                                                           !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Read properties related to the switch of bond IDs ##########################################
                                                                                         !
!   read(1,*) NSWITCH_BONDIDS;                                                           !
                                                                                         !
!   IOSEF1 = 70 - 11 - 4 - 23 - LEN_TRIM(SWITCH_FILE_NAME);                              !
                                                                                         !
!   if ( NSWITCH_BONDIDS > 0 ) then;                                                     !
                                                                                         !
!       write(icanal,'(a11,i4,a55)') '| There is ',             &                        !
!                                    NSWITCH_BONDIDS,           &                        !
!                                    ' bond IDs to swith in '// &                        !
!                                    TRIM(SWITCH_FILE_NAME)//   &                        !  
!                                    REPEAT(' ',IOSEF1)//'|';                            !
                                                                                         !
!       write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       do i = 1, NSWITCH_BONDIDS;                                                       ! 
                                                                                         !
!           read(1,*) SWITCH_BONDIDS(1:2,i);                                             !
                                                                                         !
!       end do                                                                           !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Read properties related to the switch of angle IDs #########################################
                                                                                         !
!   read(1,*) NSWITCH_ANGLEIDS;                                                          !
                                                                                         !
!   IOSEF1 = 70 - 11 - 4 - 24 - LEN_TRIM(SWITCH_FILE_NAME);                              !
                                                                                         !
!   if ( NSWITCH_ANGLEIDS > 0 ) then;                                                    !
                                                                                         !
!       write(icanal,'(a11,i4,a55)') '| There is ',              &                       !
!                                    NSWITCH_ANGLEIDS,           &                       !
!                                    ' angle IDs to swith in '// &                       !
!                                    TRIM(SWITCH_FILE_NAME)//    &                       !  
!                                    REPEAT(' ',IOSEF1)//'|';                            !
                                                                                         !
!       write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!!!!!   if ( NSWITCH_ANGLEIDS > 0 ) then;                                                    !
                                                                                         !
!       do i = 1, NSWITCH_ANGLEIDS;                                                      !
                                                                                         !
!           read(1,*) SWITCH_ANGLEIDS(1:2,i);                                            !
                                                                                         !
!       end do                                                                           !
                                                                                         !
!   end if                                                                               ! 
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Read properties related to the switch of dihedral IDs ######################################
                                                                                         !
!   read(1,*) NSWITCH_DIHEDRALIDS;                                                       !
                                                                                         !
!   IOSEF1 = 70 - 11 - 4 - 27 - LEN_TRIM(SWITCH_FILE_NAME);                              !
                                                                                         !
!   if ( NSWITCH_DIHEDRALIDS > 0 ) then;                                                 !
                                                                                         !
!       write(icanal,'(a11,i4,a55)') '| There is ',                 &                    !
!                                    NSWITCH_DIHEDRALIDS,           &                    !
!                                    ' dihedral IDs to swith in '// &                    !
!                                    TRIM(SWITCH_FILE_NAME)//       &                    !  
!                                    REPEAT(' ',IOSEF1)//'|';                            !
                                                                                         !
!       write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!!!!   if ( NSWITCH_DIHEDRALIDS > 0 ) then;                                                 !
                                                                                         !
!       do i = 1, NSWITCH_DIHEDRALIDS;                                                   !
                                                                                         !
!           read(1,*) SWITCH_DIHEDRALIDS(1:2,i);                                         !
                                                                                         !
!       end do                                                                           !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Read properties for the spatial translation of designed molecules ##########################
                                                                                         !
!   read(1,*) NTRANSLATION_MOLECULEIDS, TRANSLATION_FILE_NAME;                           !
                                                                                         !
!   IOSEF1 = 70 - 11 - 4 - 38 - LEN_TRIM(TRANSLATION_FILE_NAME);                         !
                                                                                         !
!   if ( NTRANSLATION_MOLECULEIDS > 0 ) then;                                            !
                                                                                         !
!       write(icanal,'(a11,i4,a55)') '| There is ',                            &         !
!                                    NTRANSLATION_MOLECULEIDS,                 &         !
!                                    ' molecules to spatially translate in '// &         !
!                                    TRIM(TRANSLATION_FILE_NAME)//             &         !  
!                                    REPEAT(' ',IOSEF1)//'|';                            !
                                                                                         !
!       write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       do i = 1, NTRANSLATION_MOLECULEIDS;                                              !
                                                                                         !
!           read(1,*) TRANSLATION_MOLECULEIDS(i), TRANSLATION_VECTOR(1:3,i);             ! 
                                                                                         !
!       end do                                                                           !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   close(1);                                                                            !
                                                                                         !
!   ### Allocate arrays for reading location of centers of mass of molecules to insert #############
                                                                                         !
!   if ( idefined_centers == 1 ) then;                                                   !
                                                                                         !
!       allocate(MOLECULE_COM_LIBRARY(1:3,1:NMOLECULE_INSERTION_MAX,1:NFILE_LIBRARY));   !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Initialization of arrays for the reading of centers of mass ################################
                                                                                         !
!   if ( idefined_centers == 1 ) then;                                                   !
                                                                                         !
!       MOLECULE_COM_LIBRARY(1:3,1:NMOLECULE_INSERTION_MAX,1:NFILE_LIBRARY) = 0.0d0;     !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Read centers of mass and geometrical centers for molecules to instert ######################
                                                                                         !
!   if ( idefined_centers == 1 ) then;                                                   !
                                                                                         !
!       open(1,file='input_100.dat',status='old');                                       !
                                                                                         !
!       do i = 1, 10;                                                                    !
                                                                                         !
!           read(1,*);                                                                   !
                                                                                         !
!       end do                                                                           !
                                                                                         !
!       do i = 1, NFILE_LIBRARY;                                                         !
                                                                                         !
!           read(1,*);                                                                   !
                                                                                         !
!           do j = 1, NMOLECULE_INSERTION(i);                                            !
                                                                                         !
!               read(1,*) MOLECULE_COM_LIBRARY(1:3,j,i);                                 !
                                                                                         !
!               write(icanal,*) i, j, MOLECULE_COM_LIBRARY(1:3,j,i);                     !
                                                                                         !
!           end do                                                                       !
                                                                                         !
!       end do                                                                           !
                                                                                         !
!       close(1);                                                                        !
                                                                                         ! 
!   end if                                                                               !
                                                                                         !
!   end if

!   ### Write message for the success in reanding the input file ###################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '| The input file has been read '//REPEAT(' ',38)//'|';        !
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine READ_INPUT











