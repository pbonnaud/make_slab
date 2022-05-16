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


program MAKE_SLAB100

!   ***********************************************************************************************
!   **                        PROGRAM TO CREATE MOLECULAR CONFIGURATIONS                         **
!   ***********************************************************************************************

    use module_data_in;

    use module_library;

    use module_duplicate;

    use module_config;

    use module_slab_keywords;

    use module_inserts;

    use module_simulation_box;

    use module_lmp_input;

    use module_osef;

!   ***********************************************************************************************

    implicit none;

!   ***********************************************************************************************

!   include 'mpif.h'

!   ************************************************************************************************

!   integer (kind=4) :: ierr, num_procs, my_id, status(MPI_STATUS_SIZE);

    integer (kind=4) :: root_process;

!   ************************************************************************************************

    integer (kind=4) :: i, j, k;

    integer (kind=4) :: icanal;

    integer (kind=4) :: IDATA, IFIELD, iprimitive_cell, FLAG_SORT_MOLEC;

    integer (kind=4) :: idum;

    real (kind=8) :: grnd;

    real (kind=8), dimension(1:3) :: PRIMITIVE_CELL_AXIS, PRIMITIVE_CELL_ANGDEG;

    real (kind=8), dimension(1:3) :: CELL_AXIS, CELL_ANGDEG;

    real (kind=8), dimension(1:3,1:3) :: PASSA, PASSB;

    integer (kind=4) :: PRIMITIVE_NATOM, NBRE_REPEAT_PATTERN;

    real (kind=8), dimension(1:3,1:20) :: PRIMITIVE_RI, REPEAT_PATTERN_SI;

    character (len=3), dimension(1:20) :: PRIMITIVE_NAT;

    character (len=150) :: CHEXT, CHCONFIG_PREFIX;

!   ************************************************************************************************

    real (kind=8), dimension(1:6) :: MATA;

!   ### Parametres du substrat #####################################################################

    integer (kind=4) :: ie;

    integer (kind=4) :: ICONF_FINAL;

    character (len=200) :: PATH_DIRECTORY, CH_WORKING_FILE;

!   ### WARNING MESSAGE PARAMETERS #################################################################

    integer (kind=4) :: EOF;

    logical :: PROBE1;

!   ***********************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE, CHPRGM_VERSION, CHDATE;

!   ***********************************************************************************************
                                                                                         !
    root_process = 0;                                                                    ! DEFINE PROCESSOR 0 AS THE ROOT PROCESS
                                                                                         !
!   call MPI_INIT(ierr);                                                                 !
                                                                                         !
!   call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr);                                       !
                                                                                         !
!   call MPI_COMM_SIZE(MPI_COMM_WORLD,num_procs,ierr);                                   !
                                                                                         !
!   write(*,*) 'my_id : ', my_id;                                                        !
                                                                                         !
!   ### Set the canal on which output data of the program will be written ##########################
                                                                                         ! 
    icanal = 99;                                                                         !
                                                                                         !
!   ### Write the title of the program #############################################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        open(icanal,file='OUTPUT_100.dat');                                              !
                                                                                         !
        CHTITLE        = 'Build A Simulation Box';                                       !
                                                                                         !
        CHPRGM_VERSION = 'v1.77';                                                        !   
                                                                                         !
        CHDATE         = '06/17/2021';                                                   !
                                                                                         !
        ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                               !
                                                                                         !
        call MAKE_PROGRAM_TITLE(icanal,          &                                       !
                                70,              &                                       !
                                ILENGTH_TITLE,   &                                       !
                                TRIM(CHTITLE),   &                                       !
                                CHPRGM_VERSION,  &                                       !
                                CHDATE);                                                 !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Set atom properties as found in the Mendeleiev periodic table of elements ##################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        call SET_ATOM_PROPERTIES(icanal);                                                !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Set bond properties among atoms of the Mendeleiev periodic table ###########################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        call SET_BOND_PROPERTIES(icanal);                                                !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Set keyword of the slab program ############################################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        call SET_SLAB_KEYWORDS(icanal);                                                  !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Read the input file of the current program #################################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        call READ_INPUT(icanal);                                                         !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Set the passage matrix of the simulation box to consider ###################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        call SET_PASSAGE(icanal,                  &                                      !
                         LAMMPS_CELL_AXIS(1:3),   &                                      !
                         LAMMPS_CELL_ANGDEG(1:3), &                                      !
                         LAMMPS_PASSA(1:3,1:3),   &                                      !
                         LAMMPS_PASSB(1:3,1:3));                                         !
                                                                                         !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Initialization of flags ####################################################################
                                                                                         !
    igenerate_lammps_input = 0;                                                          !
                                                                                         !
!   ### Initialization of counters #################################################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        ILAMMPS_INSERTS = 0;                                                             !
                                                                                         !
!   end if                                                                               ! 
                                                                                         !  
!   ### Loop over the list of keywords for the generation of the molecular configuration ###########
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        do i = 1, NREAD_SLAB_INPUT_KEYWORDS;                                             !
                                                                                         !
            if ( TRIM(LIST_SLAB_INPUT_KEYWORDS(i)) == &                                  ! keyword region
                 TRIM(SLAB_INPUT_KEYWORDS(1)) ) CYCLE;                                   !
                                                                                         !
            if ( TRIM(LIST_SLAB_INPUT_KEYWORDS(i)) == &                                  ! keyword create_box
                 TRIM(SLAB_INPUT_KEYWORDS(2)) ) CYCLE;                                   !
                                                                                         !
            if ( TRIM(LIST_SLAB_INPUT_KEYWORDS(i)) == &                                  ! keyword insert_modify
                 TRIM(SLAB_INPUT_KEYWORDS(5)) ) CYCLE;                                   !
                                                                                         !
            if ( TRIM(LIST_SLAB_INPUT_KEYWORDS(i)) == &                                  ! keyword lmp_input
                 TRIM(SLAB_INPUT_KEYWORDS(6)) ) CYCLE;                                   !
                                                                                         !
            if ( TRIM(LIST_SLAB_INPUT_KEYWORDS(i)) == &                                  ! keyword set_lib
                 TRIM(SLAB_INPUT_KEYWORDS(7)) ) CYCLE;                                   !
                                                                                         !
            select case(TRIM(LIST_SLAB_INPUT_KEYWORDS(i)));                              !
                                                                                         !
                case('insert');                                                          !
                                                                                         !
                    ILAMMPS_INSERTS = ILAMMPS_INSERTS + 1;                               !
                                                                                         !
                    call CHECK_INSERT_FILE_LIBRARY(icanal,ILAMMPS_INSERTS);              !
                                                                                         !
!                   stop; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                                                         !
                    call ALLOCATE_INSERT_LIBRARY_ARRAYS(icanal,ILAMMPS_INSERTS);         !
                                                                                         !
!                   stop; !//////////////////////////////////////////////////////////////!
                                                                                         !
!                   ### Library file insertion depending on the file format ########################
                                                                                         !
                    select case(iglobal_file);                                           !
                                                                                         !
                    case(1);                                                             !
                                                                                         !
                        call READ_XYZ_TEMPLATE(icanal,1,3);                              ! 
                                                                                         !
!                       stop; !//////////////////////////////////////////////////////////!
                                                                                         !
                        call INSERT_SET_INTERATOMIC_POTENTIALS(icanal,1,1, &             !
                                                               ILAMMPS_INSERTS);         !
                                                                                         !
!                       stop; !//////////////////////////////////////////////////////////!
                                                                                         !
                    case(2);                                                             !
                                                                                         !
                        call READ_LAMMPS_TEMPLATE(icanal,1,1);                           !  
                                                                                         !
                        call READ_INTERATOMIC_POTENTIALS_TEMPLATE(icanal,1,1);           !
                                                                                         !
!                       stop; !//////////////////////////////////////////////////////////!
                                                                                         !
                        if ( ILIBRARY_FILE_INFO(1) == 1 ) then;                          !
                                                                                         !
                            call READ_FILE_INFO_TEMPLATE(icanal,1);                      !
                                                                                         !
                            igenerate_lammps_input = 1;                                  !
                                                                                         !
                        end if                                                           !
                                                                                         !
                        call TEMPLATE_FILE_SET_CONFIG_NAT(icanal,1);                     !
                                                                                         !
                    case default;                                                        !
                                                                                         !
                        write(icanal,'(a70)') '| Not implemented - stop '// &            !
                                              REPEAT(' ',44)//'|';                       !
                                                                                         !
                        call CLOSING_PROGRAM(icanal,1);                                  !
                                                                                         !
                    end select                                                           !
                                                                                         !
!                   stop; !//////////////////////////////////////////////////////////////!
                                                                                         !
!                   ### Initialization of the random generator if needed ###########################
                                                                                         !
                    if ( INSERTION_METHOD_LIBRARY(1) == 3 ) then;                        !
                                                                                         !                                                                               
                        IOSEF1 = LAMMPS_INSERT_RANDOM_SEED(ILAMMPS_INSERTS);             !
                                                                                         !
                        call INITIALIZATION_RANDOM_GENERATOR(icanal,IOSEF1);             !
                                                                                         !
!                       stop; !//////////////////////////////////////////////////////////!
                                                                                         !
                    end if                                                               !
                                                                                         !
                    call INSERT_LIBRARY_MOLECULE(icanal,ILAMMPS_INSERTS);                !
                                                                                         !
!                   stop; !//////////////////////////////////////////////////////////////!
                                                                                         !
                    call DEALLOCATE_INSERT_LIBRARY_ARRAYS(icanal,ILAMMPS_INSERTS);       !
                                                                                         !
!                   stop; !//////////////////////////////////////////////////////////////!
                                                                                         !
                case default;                                                            !
                                                                                         !
                    write(icanal,'(a70)') '| Not implemented - stop '// &                !
                                          REPEAT(' ',44)//'|';                           !
                                                                                         !
                    call CLOSING_PROGRAM(icanal,1);                                      !
                                                                                         !
            end select                                                                   !
                                                                                         ! 
        end do                                                                           !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Modify the configuration for tethered atoms ################################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        if ( NLAMMPS_INSERTS > 0 ) then;                                                 !
                                                                                         !
            do j = 1, NLAMMPS_INSERTS;                                                   !
                                                                                         !
                if ( LAMMPS_INSERT_NMODIFIED(j) == 0 ) CYCLE;                            !
                                                                                         !
                do k = 1, LAMMPS_INSERT_NMODIFIED(j);                                    !
                                                                                         !
                    if ( LAMMPS_INSERT_MODIFY_STYLE(k,j) /= 'tether' ) CYCLE;            !
                                                                                         !
                    call MAKE_TETHERED_ATOMS(icanal,j,k);                                !
                                                                                         !
!                   stop; !//////////////////////////////////////////////////////////////!
                                                                                         !
                end do                                                                   ! 
                                                                                         !
            end do                                                                       !
                                                                                         ! 
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end if                                                                           !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Rearrange group unions #####################################################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                          !
                                                                                         !
            call LMP_INPUT_UPDATE_GROUP_UNION(icanal);                                   !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end if                                                                           !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Generate the prefix for the configuration files ############################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        call GENERATE_CONFIG_PREFIX(icanal,1,CHCONFIG_PREFIX);                           !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Write the XYZ configuration ################################################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        CHEXT = TRIM(CHCONFIG_PREFIX);                                                   !
                                                                                         !
        call WRITE_XYZ_CONFIG(icanal,                  &                                 !
                              CHEXT,                   &                                 !
                              LAMMPS_CELL_AXIS(1:3),   &                                 !
                              LAMMPS_CELL_ANGDEG(1:3));                                  !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Write the generated lammps configuration ###################################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        CHEXT = TRIM(CHCONFIG_PREFIX)//'_atoms.data';                                    !
                                                                                         !
        call WRITE_LAMMPS_CONFIG(icanal,                &                                !
                                 CHEXT,                 &                                !
                                 LAMMPS_CELL_AXIS(1:3), &                                !
                                 LAMMPS_MATA(1:6));                                      !
                                                                                         !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Generate the input file for lammps #########################################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        if ( igenerate_lammps_input == 1 ) then;                                         !
                                                                                         !
            CHEXT = TRIM(CHCONFIG_PREFIX)//'_input';                                     !
                                                                                         !
            call WRITE_LAMMPS_INPUT(icanal,CHEXT);                                       !
                                                                                         !
        end if                                                                           !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Deallocate arrays ##########################################################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        call DEALLOCATE_INPUT_KEYWORD_LIST(icanal);                                      ! 
                                                                                         !
        call DEALLOCATE_REGIONS(icanal);                                                 !
                                                                                         !
        call DEALLOCATE_INSERT_COMMAND(icanal);                                          !
                                                                                         !
        call DEALLOCATE_LMP_INPUT_COMMAND(icanal);                                       !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Closing the program ########################################################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   ! 
                                                                                         !
        write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                  !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',26)// &                                    !
                              ' End of program '//  &                                    !
                              REPEAT(' ',26)//'|';                                       !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                  !
                                                                                         !
        close(icanal);                                                                   !
                                                                                         !
!   end if                                                                               !
                                                                                         !
    stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!SSSSSSSSSSSSSSSSSSSSSS
!SSSSSSSSSSSSSSSSSSSSSS
!SSSSSSSSSSSSSSSSSSSSSS
!SSSSSSSSSSSSSSSSSSSSSS
!SSSSSSSSSSSSSSSSSSSSSS
!SSSSSSSSSSSSSSSSSSSSSS
!SSSSSSSSSSSSSSSSSSSSSS
!SSSSSSSSSSSSSSSSSSSSSS
!SSSSSSSSSSSSSSSSSSSSSS
!SSSSSSSSSSSSSSSSSSSSSS
!SSSSSSSSSSSSSSSSSSSSSS

!   ### Check files in the library #################################################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        if ( IBUILD_METHOD == 1 ) then;                                                  ! If the building method is based on a predefined library, then
                                                                                         !
            call CHECK_FILE_LIBRARY(99);                                                 !
                                                                                         !
        end if                                                                           !
                                                                                         !
        stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Check the file to read for duplication of the molecular structure ##########################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        if ( IBUILD_METHOD == 1 ) then;                                                  !
                                                                                         !
            if ( iduplicate == 1 ) then;                                                 !
                                                                                         !
                call CHECK_FILE_DUPLICATE(99);                                           !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
        end if                                                                           !
                                                                                         !
!   end if                                                                               ! 
                                                                                         !
!   ### Allocate arrays for molecules taken from the library #######################################
                                                                                         ! 
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        if ( IBUILD_METHOD == 1 ) then;                                                  !
                                                                                         !
            call ALLOCATE_LIBRARY_ARRAYS(99);                                            !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end if                                                                           !
                                                                                         !
!   end if                                                                               !  
                                                                                         !
!   ### Allocate arrays for molecular structures to duplicate ######################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        if ( IBUILD_METHOD == 1 ) then;                                                  !
                                                                                         !
            if ( iduplicate == 1 ) then;                                                 !
                                                                                         !
                call ALLOCATE_DUPLICATE_ARRAYS(99);                                      !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
        end if                                                                           !
                                                                                         !
!   end if                                                                               ! 
                                                                                         !
!   ### Read lammps templates ######################################################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        if ( IBUILD_METHOD == 1 ) then;                                                  !
                                                                                         !
            if ( ( TRIM(CHFILE_FORMAT) == 'lammps'      ) .OR. &                         !
                 ( TRIM(CHFILE_FORMAT) == 'padua'       ) .OR. &                         !
                 ( TRIM(CHFILE_FORMAT) == 'lammps-gaff' ) ) then;                        !
                                                                                         !
                if ( NFILE_LIBRARY > 0 ) then;                                           !
                                                                                         !
                    do i = 1, NFILE_LIBRARY;                                             !
                                                                                         !
                        call READ_LAMMPS_TEMPLATE(99,i,1);                               !
                                                                                         !
                    end do                                                               !
                                                                                         !
                end if                                                                   !
                                                                                         !
            end if                                                                       !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end if                                                                           !
                                                                                         !
!   end if                                                                               !  
                                                                                         !
!   ### Read interatomic potential template (lammps) ###############################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        if ( IBUILD_METHOD == 1 ) then;                                                  !
                                                                                         !
            if ( ( TRIM(CHFILE_FORMAT) == 'lammps'      ) .OR. &                         !
                 ( TRIM(CHFILE_FORMAT) == 'padua'       ) .OR. &                         !
                 ( TRIM(CHFILE_FORMAT) == 'lammps-gaff' ) ) then;                        !
                                                                                         !
                if ( NFILE_LIBRARY > 0 ) then;                                           !
                                                                                         !
                    do i = 1, NFILE_LIBRARY;                                             !
                                                                                         !
                        call READ_INTERATOMIC_POTENTIALS_TEMPLATE(99,i,1);               !
                                                                                         !
                    end do                                                               !
                                                                                         !
                end if                                                                   !
                                                                                         !
            end if                                                                       !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end if                                                                           !
                                                                                         !
!   end if                                                                               !  
                                                                                         !
!   ### Initialization of flags ####################################################################
                                                                                         !
    igenerate_lammps_input = 0;                                                          !
                                                                                         !
!   ### Read file info template (lammps) ###########################################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        if ( IBUILD_METHOD == 1 ) then;                                                  !
                                                                                         !
            if ( ( TRIM(CHFILE_FORMAT) == 'lammps'      ) .OR. &                         !
                 ( TRIM(CHFILE_FORMAT) == 'padua'       ) .OR. &                         !
                 ( TRIM(CHFILE_FORMAT) == 'lammps-gaff' ) ) then;                        !
                                                                                         !
                if ( NFILE_LIBRARY > 0 ) then;                                           !
                                                                                         !
                    do i = 1, NFILE_LIBRARY;                                             !
                                                                                         !
                        if ( ILIBRARY_FILE_INFO(i) == 1 ) then;                          !
                                                                                         !
                            call READ_FILE_INFO_TEMPLATE(99,i);                          !
                                                                                         !
                            igenerate_lammps_input = 1;                                  !
                                                                                         !
                        end if                                                           !
                                                                                         !
                    end do                                                               !
                                                                                         !
                end if                                                                   !
                                                                                         !
            end if                                                                       !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end if                                                                           !
                                                                                         !
!   end if                                                                               !  
                                                                                         !
!   ### Assign atom labels to library data #########################################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        if ( IBUILD_METHOD == 1 ) then;                                                  !
                                                                                         !
            if ( ( TRIM(CHFILE_FORMAT) == 'lammps'      ) .OR. &                         !
                 ( TRIM(CHFILE_FORMAT) == 'padua'       ) .OR. &                         !
                 ( TRIM(CHFILE_FORMAT) == 'lammps-gaff' ) ) then;                        !
                                                                                         !
                if ( NFILE_LIBRARY > 0 ) then;                                           !
                                                                                         !
                    do i = 1, NFILE_LIBRARY;                                             !
                                                                                         !
                        CONFIG_NAT_LIBRARY(1:NATOM_LIBRARY(i),i) = 'XXX';                !
                                                                                         !
                        do j = 1, NATOM_LIBRARY(i);                                      !
                                                                                         !
                            CONFIG_NAT_LIBRARY(j,i) = &                                  !
                            ATOM_LABEL_LIBRARY(CONFIG_ATOM_TYPE_LIBRARY(j,i),i);         !
                                                                                         !
                        end do                                                           !
                                                                                         !
                    end do                                                               !
                                                                                         !
                end if                                                                   !
                                                                                         !
            end if                                                                       !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end if                                                                           !
                                                                                         !
!   end if                                                                               !  
                                                                                         !
!   ### Read the file to duplicate in the xyz format ###############################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        if ( IBUILD_METHOD == 1 ) then;                                                  !
                                                                                         !
            if ( iduplicate == 1 ) then;                                                 !
                                                                                         !
                if ( TRIM(DUPLICATE_CHFILE_FORMAT) == 'xyz' ) then;                      !
                                                                                         !
                    call DUPLICATE_READ_XYZ_TEMPLATE(99,1);                              !
                                                                                         !
!                   stop; !//////////////////////////////////////////////////////////////!
                                                                                         !
                end if                                                                   !
                                                                                         !
            end if                                                                       !
                                                                                         !
        end if                                                                           !
                                                                                         !
!   end if                                                                               ! 
                                                                                         !
!   ### Set properties of the primitive simulation box for duplication #############################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        if ( IBUILD_METHOD == 1 ) then;                                                  !
                                                                                         !
            if ( iduplicate == 1 ) then;                                                 !
                                                                                         !
                CELL_AXIS(1) = CELL_AXIS_DUPLICATE(1) * DUPLICATE_NX;                    !
                                                                                         !
                CELL_AXIS(2) = CELL_AXIS_DUPLICATE(2) * DUPLICATE_NY;                    !
                                                                                         ! 
                CELL_AXIS(3) = CELL_AXIS_DUPLICATE(3) * DUPLICATE_NZ;                    !
                                                                                         !
                CELL_ANGDEG(1:3) = CELL_ANGDEG_DUPLICATE(1:3);                           !
                                                                                         !
                call SET_PASSAGE(99,               &                                     !
                                 CELL_AXIS(1:3),   &                                     !
                                 CELL_ANGDEG(1:3), &                                     !
                                 PASSA(1:3,1:3),   &                                     !
                                 PASSB(1:3,1:3));                                        !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
        end if                                                                           !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Buid the supercell #########################################################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        if ( IBUILD_METHOD == 1 ) then;                                                  !
                                                                                         !
            if ( iduplicate == 1 ) then;                                                 !
                                                                                         !
                call DUPLICATE_BUILD_SUPER_CELL(99,PASSA(1:3,1:3),PASSB(1:3,1:3));       !
                                                                                         ! 
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
        end if                                                                           !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Topology of the molecular configuration ####################################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        if ( IBUILD_METHOD == 1 ) then;                                                  !
                                                                                         !
            if ( ( iduplicate == 1 ) .AND. ( iduplicate_ffield == 1 ) ) then;            !
                                                                                         !
                call COMPUTE_CONFIG_TOPOLOGY(99,PASSA(1:3,1:3),PASSB(1:3,1:3));          !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
        end if                                                                           !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Assign the force field on the molecular configuration ######################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        if ( IBUILD_METHOD == 1 ) then;                                                  !
                                                                                         !
            if ( ( iduplicate == 1 ) .AND. ( iduplicate_ffield == 1 ) ) then;            !
                                                                                         !
                call DUPLICATE_ASSIGN_FORCE_FIELD(99);                                   !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
        end if                                                                           !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Loop over the number of configurations to generate #########################################
                                                                                         !
!           do k = 1, NMOLECULAR_CONFIG;                                                 !
                                                                                         !

!   ### Set properties of the simualtion box #######################################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        if ( IBUILD_METHOD == 1 ) then;                                                  !
                                                                                         !
            if ( iduplicate == 1 ) then;                                                 !
                                                                                         !
                MATA(1) = PASSA(1,1);                                                    !
                                                                                         !
                MATA(2) = PASSA(2,2);                                                    !
                                                                                         !
                MATA(3) = PASSA(3,3);                                                    !
                                                                                         !
                MATA(4) = PASSA(2,3);                                                    !
                                                                                         !
                MATA(5) = PASSA(1,3);                                                    !
                                                                                         !
                MATA(6) = PASSA(1,2);                                                    !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
        end if                                                                           !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Generate the prefix for the configuration files ############################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        if ( IBUILD_METHOD == 1 ) then;                                                  !
                                                                                         !
            if ( iduplicate == 1 ) then;                                                 !
                                                                                         !
                call GENERATE_CONFIG_PREFIX(99,1,CHCONFIG_PREFIX);                       !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
        end if                                                                           !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Write the XYZ configuration ################################################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        if ( IBUILD_METHOD == 1 ) then;                                                  !
                                                                                         !
            if ( iduplicate == 1 ) then;                                                 !
                                                                                         !
                CHEXT = TRIM(CHCONFIG_PREFIX);                                           !
                                                                                         !
                call WRITE_XYZ_CONFIG(99,CHEXT,CELL_AXIS(1:3),CELL_ANGDEG(1:3));         !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
        end if                                                                           !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Write the generated lammps configuration ###################################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        if ( IBUILD_METHOD == 1 ) then;                                                  !
                                                                                         !
            if ( iduplicate == 1 ) then;                                                 !
                                                                                         !
                CHEXT = TRIM(CHCONFIG_PREFIX)//'_atoms.data';                            !
                                                                                         !
                call WRITE_LAMMPS_CONFIG(99,             &                               !
                                         CHEXT,          &                               !
                                         CELL_AXIS(1:3), &                               !
                                         MATA(1:6));                                     !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
        end if                                                                           !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   #### Generate the input file for lammps ########################################################
                                                                                         !
!  if ( my_id == root_process ) then;                                                    !
                                                                                         !
        if ( IBUILD_METHOD == 1 ) then;                                                  !
                                                                                         !
            if ( iduplicate == 1 ) then;                                                 !
                                                                                         !
                if ( igenerate_lammps_input == 1 ) then;                                 !
                                                                                         !
                    CHEXT = TRIM(CHCONFIG_PREFIX)//'_input';                             !
                                                                                         !
                    call WRITE_LAMMPS_INPUT(99,CHEXT);                                   !
                                                                                         !
                end if                                                                   !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
        end if                                                                           !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Other modes ################################################################################

    if ( iduplicate == 0 ) then;

!   ### Check files in the library #################################################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        if ( IBUILD_METHOD == 1 ) then;                                                  ! If the building method is based on a predefined library, then
                                                                                         !
            CELL_AXIS(1:3)   = PRIMITIVE_CELL_AXIS(1:3);                                 !
                                                                                         !
            CELL_ANGDEG(1:3) = PRIMITIVE_CELL_ANGDEG(1:3);                               !
                                                                                         !
!           call CHECK_FILE_LIBRARY(99);                                                 !
                                                                                         !
        end if                                                                           !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Build of modify the molecular system #######################################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
!       ### If the chosen method for building or modifying the system is from a library ############
                                                                                         !
        if ( IBUILD_METHOD == 1 ) then;                                                  !
                                                                                         !
            MATA(1:6) = 0.0d0;                                                           !
                                                                                         !
!           ### Allocate arrays for molecules taken from the library ###############################
                                                                                         !
!           call ALLOCATE_LIBRARY_ARRAYS(99);                                            !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           ### Read Library files in the xyz format ###############################################
                                                                                         !
            if ( TRIM(CHFILE_FORMAT) == 'xyz' ) then;                                    !
                                                                                         !
                if ( NFILE_LIBRARY > 0 ) then;                                           ! IF THE NUMBER OF MOLECULES TO INSERT IS GREATER THAN 0,THEN
                                                                                         !
                    do i = 1, NFILE_LIBRARY;                                             ! LOOP OVER THE NUMBER OF FILES TO READ IN THE LIBRARY
                                                                                         !
                        call READ_XYZ_TEMPLATE(99,i,1);                                  ! 
                                                                                         !
                    end do                                                               !
                                                                                         !
                end if                                                                   !
                                                                                         !
            end if                                                                       !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!          ### Read library files in the lammps format #############################################
                                                                                         !
            if ( ( TRIM(CHFILE_FORMAT) == 'lammps'      ) .OR. &                         !
                 ( TRIM(CHFILE_FORMAT) == 'padua'       ) .OR. &                         !
                 ( TRIM(CHFILE_FORMAT) == 'lammps-gaff' ) ) then;                        !
                                                                                         !
                if ( NFILE_LIBRARY > 0 ) then;                                           ! IF THE NUMBER OF MOLECULES TO INSERT IS GREATER THAN 0,THEN
                                                                                         !
!                   do i = 1, NFILE_LIBRARY;                                             ! LOOP OVER THE NUMBER OF FILES TO READ IN THE LIBRARY
                                                                                         !
!                       call READ_LAMMPS_TEMPLATE(99,i,1);                               !
                                                                                         !
!                       stop; !//////////////////////////////////////////////////////////!
                                                                                         !
!                       call READ_INTERATOMIC_POTENTIALS_TEMPLATE(99,i,1);               !
                                                                                         !
!                       stop; !//////////////////////////////////////////////////////////!
                                                                                         !
!                       if ( ILIBRARY_FILE_INFO(i) == 1 ) then;                          !
                                                                                         !
!                           call READ_FILE_INFO_TEMPLATE(99,i);                          !
                                                                                         !
!                           igenerate_lammps_input = 1;                                  !
                                                                                         !
!                       end if                                                           !
                                                                                         !
!                       stop; !//////////////////////////////////////////////////////////!
                                                                                         !
!                       CONFIG_NAT_LIBRARY(1:NATOM_LIBRARY(i),i) = 'XXX';                !
                                                                                         !
!                       do j = 1, NATOM_LIBRARY(i);                                      !
                                                                                         !
!                           CONFIG_NAT_LIBRARY(j,i) = &                                  !
!                           ATOM_LABEL_LIBRARY(CONFIG_ATOM_TYPE_LIBRARY(j,i),i);         !
                                                                                         !
!                       end do                                                           !
                                                                                         !
!                   end do                                                               !
                                                                                         !
                    if ( iworking_file_insertion == 1 ) then;                            !
                                                                                         !
                        call READ_LAMMPS_CONFIG(99,                       &              !
                                                CHWORKING_FILE_INSERTION, &              !
                                                MATA(1:6),                &              !
                                                FLAG_SORT_MOLEC);                        !
                                                                                         !
                        PRIMITIVE_CELL_AXIS(1:3) = MATA(1:3);                            !
                                                                                         ! 
                        if ( MATA(4) == 0.0d0 ) PRIMITIVE_CELL_ANGDEG(1) = 90.0d0;       !
                                                                                         !
                        if ( MATA(5) == 0.0d0 ) PRIMITIVE_CELL_ANGDEG(2) = 90.0d0;       !
                                                                                         !
                        if ( MATA(6) == 0.0d0 ) PRIMITIVE_CELL_ANGDEG(3) = 90.0d0;       ! 
                                                                                         !
                    end if                                                               !
                                                                                         !
!                   stop; !//////////////////////////////////////////////////////////////!
                                                                                         !
                end if                                                                   !
                                                                                         !
                if ( NFILE_LIBRARY_REMOVE > 0 ) then;                                    !
                                                                                         !
                    do i = 1, NFILE_LIBRARY_REMOVE;                                      ! LOOP OVER THE NUMBER OF MOLECULES TO BE REMOVED
                                                                                         !
                        call READ_LAMMPS_TEMPLATE(99,i,0);                               !
                                                                                         !
!                       stop; !//////////////////////////////////////////////////////////!
                                                                                         !
                        call READ_INTERATOMIC_POTENTIALS_TEMPLATE(99,i,0);               !
                                                                                         !
!                       stop; !//////////////////////////////////////////////////////////!
                                                                                         !
                        CONFIG_NAT_LIBRARY(1:NATOM_LIBRARY(i),i) = 'XXX';                !
                                                                                         !
                        do j = 1, NATOM_LIBRARY(i);                                      !
                                                                                         !
                            CONFIG_NAT_LIBRARY(j,i) = &                                  !
                            ATOM_LABEL_LIBRARY(CONFIG_ATOM_TYPE_LIBRARY(j,i),i);         !
                                                                                         ! 
                        end do                                                           ! 
                                                                                         !
                    end do                                                               !
                                                                                         !
!                   stop; !//////////////////////////////////////////////////////////////!
                                                                                         !
                    call READ_LAMMPS_CONFIG(99,                    &                     !
                                            CHWORKING_FILE_REMOVE, &                     !
                                            MATA(1:6),             &                     !
                                            FLAG_SORT_MOLEC);                            !
                                                                                         !
                    PRIMITIVE_CELL_AXIS(1:3) = MATA(1:3);                                !
                                                                                         ! 
                    if ( MATA(4) == 0.0d0 ) PRIMITIVE_CELL_ANGDEG(1) = 90.0d0;           !
                                                                                         !
                    if ( MATA(5) == 0.0d0 ) PRIMITIVE_CELL_ANGDEG(2) = 90.0d0;           !
                                                                                         !
                    if ( MATA(6) == 0.0d0 ) PRIMITIVE_CELL_ANGDEG(3) = 90.0d0;           !
                                                                                         !
                end if                                                                   !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           ### Read library files for the building of polymers ####################################
                                                                                         !
            if ( NGenerate_polymer_species > 0 ) then;                                   !
                                                                                         !
                do i = 1, NTotal_monomers;                                               !
                                                                                         !
                   call READ_XYZ_TEMPLATE(99,i,2);                                       !
                                                                                         !
                end do                                                                   ! 
                                                                                         ! 
            end if                                                                       !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           ### Look for radical location in library monomers ######################################
                                                                                         !
            if ( NGenerate_polymer_species > 0 ) then;                                   !
                                                                                         !
                call FIND_MONOMER_RADICALS(icanal,PASSA,PASSB);                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           ### Set properties of the simulation box ###############################################
                                                                                         !
            call SET_PASSAGE(99,                         &                               !
                             PRIMITIVE_CELL_AXIS(1:3),   &                               !
                             PRIMITIVE_CELL_ANGDEG(1:3), &                               !
                             PASSA(1:3,1:3),             &                               !
                             PASSB(1:3,1:3));                                            !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           ### Allocate arrays for the generation of the molecular configuration ##################
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            if ( ( NFILE_LIBRARY > 0 ) .AND. ( iworking_file_insertion == 0 ) ) then;    !
                                                                                         !
                IOSEF1 = 1;                                                              !
                                                                                         !
            end if                                                                       !
                                                                                         !
            if ( NGenerate_polymer_species > 0 ) IOSEF1 = 1;                             !
                                                                                         !
            if ( IOSEF1 == 1 ) call ALLOCATE_CONFIG_ARRAYS(99);                          !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           ### Initialization of the random generator #############################################
                                                                                         !
!           call INITIALIZATION_RANDOM_GENERATOR(icanal);                                !
            call INITIALIZATION_RANDOM_GENERATOR(icanal,RANDOM_GENERATOR_SEED);          !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           ### Loop over the number of configurations to generate #################################
                                                                                         !
            do k = 1, NMOLECULAR_CONFIG;                                                 !
                                                                                         !
!               ### Insertion of molecules / particles in the simulation box #######################
                                                                                         !
                if ( NFILE_LIBRARY > 0 ) then;                                           !
                                                                                         !
!                   ### Random insertion of molecules in the simulation box ########################
                                                                                         !
                    call RANDOM_INSERTION_MOLECULE(99,PASSA(1:3,1:3),PASSB(1:3,1:3));    !
                                                                                         !
!                   stop; !//////////////////////////////////////////////////////////////!
                                                                                         !
!                   ### Random insetion of particles in the simulation box #########################
                                                                                         !
                    if ( NTYPE_DUMMY_PARTICLE > 0 ) then;                                !
                                                                                         !
                        call RANDOM_INSERTION_DUMMY_PARTICLE(99,              &          !
!                                                            IBUILD_METHOD,   &          ! 
                                                             PASSA(1:3,1:3),  &          !
                                                             PASSB(1:3,1:3));            !                                  
                                                                                         !
                    end if                                                               !
                                                                                         !
!                   stop; !//////////////////////////////////////////////////////////////!
                                                                                         !
                end if                                                                   !
                                                                                         !
!               ### Randomly remove some molecules from the simulation box #########################
                                                                                         !
                if ( NFILE_LIBRARY_REMOVE > 0 ) then;                                    !
                                                                                         !
                    call RANDOM_DELETION_MOLECULE(99,PASSA(1:3,1:3),PASSB(1:3,1:3));     !
                                                                                         !
                end if                                                                   !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
!               ### Generate polymer chains in the simulation box ##################################
                                                                                         !
                if ( NGenerate_polymer_species > 0 ) then;                               !
                                                                                         !
                    call BUILD_POLYMER_CHAINS(99,PASSA(1:3,1:3),PASSB(1:3,1:3));         !
                                                                                         !
                end if                                                                   !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
!               ### Check intramolecular bonds #####################################################
                                                                                         !
                if ( NGenerate_polymer_species > 0 ) then;                               !
                                                                                         !
                    call CHECK_INTRAMOLECULAR_BONDS(icanal,PASSA,PASSB);                 !
                                                                                         !
                end if                                                                   !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!

!               ### Topology ###

!               ### Check the charge of the system (neutral) #######################################
                                                                                         !
                call CHECK_ELECTRONEUTRALITY(99);                                        !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
!               ### Set properties of the simualtion box ###########################################
                                                                                         !
!               MATA(1) = PASSA(1,1);                                                    !
                                                                                         !
!               MATA(2) = PASSA(2,2);                                                    !
                                                                                         !
!               MATA(3) = PASSA(3,3);                                                    !
                                                                                         !
!               MATA(4) = PASSA(2,3);                                                    !
                                                                                         !
!               MATA(5) = PASSA(1,3);                                                    !
                                                                                         !
!               MATA(6) = PASSA(1,2);                                                    !
                                                                                         !
!               ### Generate the prefix for the configuration files ################################
                                                                                         !
                call GENERATE_CONFIG_PREFIX(99,k,CHCONFIG_PREFIX);                       !
                                                                                         !
!               ### Write the XYZ configuration ####################################################
                                                                                         !
                CHEXT = TRIM(CHCONFIG_PREFIX);                                           !
                                                                                         !
                call WRITE_XYZ_CONFIG(99,CHEXT,CELL_AXIS(1:3),CELL_ANGDEG(1:3));         !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
!               ### Write the generated lammps configuration #######################################
                                                                                         !
                CHEXT = TRIM(CHCONFIG_PREFIX)//'_atoms.data';                            !
                                                                                         !
                call WRITE_LAMMPS_CONFIG(99,             &                               !
                                         CHEXT,          &                               !
                                         CELL_AXIS(1:3), &                               !
                                         MATA(1:6));                                     !
                                                                                         !
!               ### Generate the input file for lammps #############################################
                                                                                         !
                if ( igenerate_lammps_input == 1 ) then;                                 !
                                                                                         !
                    CHEXT = TRIM(CHCONFIG_PREFIX)//'_input';                             !
                                                                                         !
                    call WRITE_LAMMPS_INPUT(99,CHEXT);                                   !
                                                                                         !
                end if                                                                   !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end do                                                                       !
                                                                                         !
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
!WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW

!           ### Deallocate library arrays ##########################################################
                                                                                         ! 
            call DEALLOCATE_LIBRARY_ARRAYS(99);                                          !
                                                                                         !
        else if ( IBUILD_METHOD == 2 ) then;                                             !
                                                                                         !
            if ( ( NFILE_LIBRARY > 0 ) .OR. ( NFILE_LIBRARY_REMOVE > 0 ) ) then;         !
                                                                                         !
                call ALLOCATE_CONFIG_LIBRARY_ARRAYS(99);                                 !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
                do i = 1, NFILE_LIBRARY;                                                 ! LOOP OVER THE NUMBER OF FILES TO READ IN THE LIBRARY
                                                                                         !
                    call READ_LAMMPS_CONFIG_LIBRARY(99,i);                               !
                                                                                         !
                    CONFIG_NAT_LIBRARY(1:NATOM_LIBRARY(i),i) = 'XXX';                    !
                                                                                         !
                    do j = 1, NATOM_LIBRARY(i);                                          !
                                                                                         !
                        CONFIG_NAT_LIBRARY(j,i) = &                                      !
                        ATOM_LABEL_LIBRARY(CONFIG_ATOM_TYPE_LIBRARY(j,i),i);             !
                                                                                         !
                    end do                                                               !
                                                                                         !
                end do                                                                   !
                                                                                         !
            end if                                                                       !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         ! 
            call SET_PASSAGE(99,                         &                               !
                             PRIMITIVE_CELL_AXIS(1:3),   &                               !
                             PRIMITIVE_CELL_ANGDEG(1:3), &                               !
                             PASSA(1:3,1:3),             &                               !
                             PASSB(1:3,1:3));                                            !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
            call ALLOCATE_CONFIG_ARRAYS(99);                                             !
                                                                                         !
            call RANDOM_INSERTION_MOLECULE(99,                 &                         !
                                           PASSA(1:3,1:3),     &                         !
                                           PASSB(1:3,1:3));                              !
                                                                                         !
            MATA(1) = PASSA(1,1);                                                        !
            MATA(2) = PASSA(2,2);                                                        !
            MATA(3) = PASSA(3,3);                                                        !
            MATA(4) = PASSA(2,3);                                                        !
            MATA(5) = PASSA(1,3);                                                        !
            MATA(6) = PASSA(1,2);                                                        !
                                                                                         !
            CONFIG_MOLECULEID(1:NATOM) = 0;                                              !
                                                                                         !
            call WRITE_LAMMPS_CONFIG(99,               &                                 !
                                     'atoms.data_new', &                                 !
                                     CELL_AXIS(1:3),   &                                 !
                                     MATA(1:6));                                         !
                                                                                         !
            do j = 1, NATOM;                                                             !
                CONFIG_NAT(j) = &                                                        !
                ATOM_LABEL(CONFIG_ATOM_TYPE(j));                                         !
            end do                                                                       !
                                                                                         !
!           ### Write the xyz molecular configuration ##############################################
                                                                                         !
            CHEXT = 'XXX';                                                               !
                                                                                         !
            call WRITE_XYZ_CONFIG(99,                  &                                 !
                                  CHEXT,               &                                 !
                                  CELL_AXIS(1:3),      &                                 !
                                  CELL_ANGDEG(1:3));                                     !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        else if ( IBUILD_METHOD == 3 ) then;                                             !
                                                                                         !
            CHEXT = '000_atoms.data';                                                    !
                                                                                         !
            if ( NSWITCH_ATOMIDS > 0 ) CHEXT = TRIM(SWITCH_FILE_NAME);                   !
                                                                                         !
            call READ_LAMMPS_CONFIG(99,                    &                             !
                                    CHEXT,                 &                             !
                                    MATA(1:6),             &                             !
                                    FLAG_SORT_MOLEC);                                    !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
            CELL_AXIS(1:3) = MATA(1:3);                                                  !
                                                                                         !
            if ( MATA(4) == 0.0d0 ) CELL_ANGDEG(1) = 90.0d0;                             !
                                                                                         !
            if ( MATA(5) == 0.0d0 ) CELL_ANGDEG(2) = 90.0d0;                             !
                                                                                         !
            if ( MATA(6) == 0.0d0 ) CELL_ANGDEG(3) = 90.0d0;                             !
                                                                                         !
            call SET_PASSAGE(99,                         &                               !
                             CELL_AXIS(1:3),             &                               !
                             CELL_ANGDEG(1:3),           &                               !
                             PASSA(1:3,1:3),             &                               !
                             PASSB(1:3,1:3));                                            !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
            if ( NTYPE_DUMMY_PARTICLE > 0 ) then;                                        !
                                                                                         !
                call RANDOM_INSERTION_DUMMY_PARTICLE(99,              &                  !
!                                                    IBUILD_METHOD,   &                  ! 
                                                     PASSA(1:3,1:3),  &                  !
                                                     PASSB(1:3,1:3));                    !                                  
                                                                                         !
            end if                                                                       !
                                                                                         !
!           ### Switch IDs #########################################################################
                                                                                         !
            call MAKE_SWITCH_ATOMIDS(99);                                                !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           ### Translated molecules ###############################################################
                                                                                         !
            call MAKE_TRANSLATION_MOLECULEIDS(99);                                       !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           ### Write the LAMMPS configuration #####################################################
                                                                                         !
            CHEXT = 'XXX_atoms.data_new';                                                !
                                                                                         !
            call WRITE_LAMMPS_CONFIG(99,CHEXT,CELL_AXIS(1:3),MATA(1:6));                 !
                                                                                         !
!           ### Write the xyz configuration ########################################################
                                                                                         !
            CHEXT = 'XXX';                                                               !
                                                                                         !
            call WRITE_XYZ_CONFIG(99,CHEXT,CELL_AXIS(1:3),CELL_ANGDEG(1:3));             !
                                                                                         !
        else                                                                             !
                                                                                         !
            write(99,'(a70)') '| Not implemented'//REPEAT(' ',52)//'|';                  !
                                                                                         !
            write(99,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
            write(99,'(a70)') '+'//REPEAT('-',68)//'+';                                  !
                                                                                         !
            write(99,*);                                                                 !
                                                                                         !
            write(99,'(a14)') 'End of program';                                          !
                                                                                         !
            close(99);                                                                   !
                                                                                         !
            stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end if                                                                           !
                                                                                         !
!       ### DEALLOCATION OF ARRAYS #################################################################
                                                                                         !
        if ( iregions == 1 ) then;                                                       !
                                                                                         !
            deallocate(REGION_AXIS);                                                     !
                                                                                         !
            deallocate(REGION_BOUNDS);                                                   !
                                                                                         !
        end if                                                                           !
                                                                                         !

    end if



!   end if

!   ### Closing the program ########################################################################
                                                                                         !
!   if ( my_id == root_process ) then;                                                   ! 
                                                                                         !
!       write(99,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!       write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       write(99,'(a70)') '|'//REPEAT(' ',26)//' End of program '//REPEAT(' ',26)//'|';  !
                                                                                         !
!       write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       write(99,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!   XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!   if ( my_id == root_process ) then; 

!   PATH_DIRECTORY = './';                                                               !
                                                                                         !
!   CH_WORKING_FILE = TRIM(CHFSUB);                                                      !
                                                                                         !
!   inquire(FILE=TRIM(TRIM(PATH_DIRECTORY)//TRIM(CH_WORKING_FILE)//'.amc'),EXIST=PROBE1);!
                                                                                         !
!   if ( PROBE1 .EQV. .TRUE. ) then;                                                     !
                                                                                         !
!       write(99,*) 'THE INPUT CONFIGURATION IS AMC FORMAT';                             !
                                                                                         !
!       write(99,*) 'THIS PROGRAM IS NOT DESIGN FOR SUCH FILES';                         !
                                                                                         !
!       write(99,*) 'END OF THE PROGRAM';                                                !
                                                                                         !
!       write(99,*);                                                                     !
                                                                                         !
!       close(99);                                                                       !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   call READ_DATA(99,PATH_DIRECTORY,CH_WORKING_FILE,IDATA);                             !
                                                                                         !
!   call READ_FIELD(99,PATH_DIRECTORY,CH_WORKING_FILE,IFIELD);                           !
                                                                                         !
!   PASSA(1:3,1:3) = 0.0d0;

!   PASSB(1:3,1:3) = 0.0d0;

!   if ( iprimitive_cell == 1 ) then

!       call SET_PASSAGE(99,PRIMITIVE_CELL_AXIS(1:3),   &
!                           PRIMITIVE_CELL_ANGDEG(1:3), &
!                           PASSA(1:3,1:3),             &
!                           PASSB(1:3,1:3));

!       call BUILD_PRIMITIVE_CELL(99,PATH_DIRECTORY, &
!                                 CH_WORKING_FILE,PRIMITIVE_CELL_AXIS(1:3),PRIMITIVE_CELL_ANGDEG(1:3), &
!                                 PASSA(1:3,1:3),PASSB(1:3,1:3),PRIMITIVE_NATOM,PRIMITIVE_NAT(1:20),  &
!                                 PRIMITIVE_RI(1:3,1:20),NBRE_REPEAT_PATTERN,REPEAT_PATTERN_SI(1:3,1:20));

!       call BUILD_SUPER_CELL(99,PATH_DIRECTORY,CH_WORKING_FILE);

!   else
!       if ( IDATA == 1 ) then
!           allocate(DATA(1:espece));

!           do i = 1, espece
!               allocate(DATA(i)%NAT(1:NATSP(i)));
!               allocate(DATA(i)%LABEL(1:NATSP(i)));
!               allocate(DATA(i)%INMOL(1:NATSP(i)));
!               allocate(DATA(i)%RI(1:3,1:NATSP(i)));
!               allocate(DATA(i)%RG(1:3,1:NATSP(i)));
!           end do

!           if ( chslab == 1 ) allocate(SUB(1:NSLAB));

!       end if

!      if ( ( ipsub == 1 ) .AND. & 
!           ( ifsub == 1 ) ) call READ_CONFIG(PATH_DIRECTORY,    &
!                                             CH_WORKING_FILE,   &
!                                             IDATA,             &
!                                             IFIELD,            &
!                                             CELL_AXIS(1:3),    &
!                                             CELL_ANGDEG(1:3));

!      call SET_PASSAGE(99);

!       if ( imodif == 1 ) then
!           call APPLY_MODIFICATION();
!       else
!           if ( ipsub == 1 .AND. ifsub == 1 ) then
!               call MAKE_SUPCELLCONF();                       ! CREATE THE SUPERCELL
!           else if ( ipsub == 1 .AND. ifsub == 0 ) then
!               NWORK = NMULTI * NMOTIF * NMAILLE;
!           end if
!       end if

!       read(CH_WORKING_FILE,*) ICONF_FINAL;

!       ICONF_FINAL = ICONF_FINAL + 1;

!       write(CH_WORKING_FILE,'(i3.3)') ICONF_FINAL;

!       write(99,*) CH_WORKING_FILE;

!       call WRITE_CONFIG(PATH_DIRECTORY,  &
!                         CH_WORKING_FILE, &
!                         CELL_AXIS(1:3),  &
!                         CELL_ANGDEG(1:3));

!       if ( IFIELD == 1 ) call WRITE_FIELD(PATH_DIRECTORY,CH_WORKING_FILE);

!   end if

!!!!   write(99,'(a22)') '!!! END OF PROGRAM !!!';

    close(99);

!   end if

!   ### 

!   call MPI_FINALIZE(ierr);                                                             !

end program MAKE_SLAB100

subroutine sgrnd(seed)

      implicit integer(a-z)

!* Period parameters
      parameter(N     =  624)

      dimension mt(0:N-1)
!                     the array for the state vector
      common /block/mti,mt
      save   /block/
!*
!*      setting initial seeds to mt[N] using
!*      the generator Line 25 of Table 1 in
!*      [KNUTH 1981, The Art of Computer Programming
!*         Vol. 2 (2nd Ed.), pp102]
!*
      mt(0)= iand(seed,-1)
      do 1000 mti=1,N-1
        mt(mti) = iand(69069 * mt(mti-1),-1)
 1000 continue

      return

end subroutine sgrnd

!************************************************************************

double precision function grnd()

      implicit integer(a-z)

!* Period parameters
      parameter(N     =  624)
      parameter(N1    =  N+1)
      parameter(M     =  397)
      parameter(MATA  = -1727483681)
                                 !   constant vector a
!     parameter(UMASK = -2147483648)

      integer (kind=4), parameter :: UMASK = -2147483648;

                                  !  most significant w-r bits
      parameter(LMASK =  2147483647)
                                   ! least significant r bits
! Tempering parameters

      parameter(TMASKB= -1658038656)
      parameter(TMASKC= -272236544)

      dimension mt(0:N-1)
!*                     the array for the state vector
      common /block/mti,mt
      save   /block/
      data   mti/N1/
!*                     mti==N+1 means mt[N] is not initialized
!*
      dimension mag01(0:1)
      data mag01/0, MATA/
      save mag01
!*                        mag01(x) = x * MATA for x=0,1
!*
      TSHFTU(y)=ishft(y,-11)
      TSHFTS(y)=ishft(y,7)
      TSHFTT(y)=ishft(y,15)
      TSHFTL(y)=ishft(y,-18)
!*
      if(mti.ge.N) then
!*                       generate N words at one time
        if(mti.eq.N+1) then
!*                            if sgrnd() has not been called,
          call sgrnd(4357)
!*                              a default initial seed is used
        endif
!*
        do 1000 kk=0,N-M-1
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
 1000   continue
        do 1100 kk=N-M,N-2
            y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
            mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
 1100   continue
        y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
        mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
        mti = 0
      endif
!*
      y=mt(mti)
      mti=mti+1
      y=ieor(y,TSHFTU(y))
      y=ieor(y,iand(TSHFTS(y),TMASKB))
      y=ieor(y,iand(TSHFTT(y),TMASKC))
      y=ieor(y,TSHFTL(y))
!*
      if(y.lt.0) then
        grnd=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)
      else
        grnd=dble(y)/(2.0d0**32-1.0d0)
      endif
!*
      return

end function grnd
