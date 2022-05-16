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

subroutine ALLOCATE_INSERT_LIBRARY_ARRAYS(icanal,ILOCAL_INSERTS)

!   ************************************************************************************************
!   **                      Allocate library arrays for the insert command                        **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                           **
!   **                                                                                            **
!   ************************************************************************************************

    use module_size_arrays;

    use module_library;

    use module_osef;

    use module_inserts;

!   ************************************************************************************************

    implicit none;

!   ***********************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    integer (kind=4), intent(in) :: ILOCAL_INSERTS;

!   ************************************************************************************************

    integer (kind=4) :: i, j, k;

    integer (kind=4) :: IFOUND;

    integer (kind=4) :: ilocal_file

    character (len=20), dimension(1:100) :: LOCAL_ATOM_LABEL;

!   ************************************************************************************************

    integer (kind=4) :: NFILE_LIBRARY_TMP;    

!   ************************************************************************************************

    integer (kind=4) :: IOSEF_NATOM,     &
                        IOSEF_NBOND,     &
                        IOSEF_NANGLE,    &
                        IOSEF_NDIHEDRAL, &
                        IOSEF_NIMPROPER;

    integer (kind=4) :: IOSEF_NTYPE_NATOM,       &
                        IOSEF_NTYPE_BOND,        &
                        IOSEF_NTYPE_ANGLE,       &
                        IOSEF_NTYPE_DIHEDRAL,    &
                        IOSEF_NTYPE_IMPROPER,    &
                        IOSEF_NPAIR_COEFF_CROSS;

    character (len=250) :: CHEXT;

    character (len=250) :: LOCAL_CHNAME_FILE;

!   ************************************************************************************************

    character (len=250) :: CHARLINE;

    integer (kind=4) :: EOF;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the routine #############################################################
                                                                                         !
    CHTITLE = 'Allocate library arrays for the insert command';                          !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the insert command to consider #######################################################
                                                                                         !
    write(icanal,'(a19,i4,a47)') '| ILOCAL_INSERTS : ', &                                !
                                 ILOCAL_INSERTS,        &                                !
                                 REPEAT(' ',46)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set the number of files to read in the library #############################################
                                                                                         !
    NFILE_LIBRARY = 1;                                                                   !
                                                                                         !
    NFILE_LIBRARY_TMP = NFILE_LIBRARY;                                                   !
                                                                                         !
    write(icanal,'(a25,i4,a41)') '| NFILE_LIBRARY_TMP    : ', &                          !
                                 NFILE_LIBRARY_TMP,           &                          !
                                 REPEAT(' ',40)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate arrays for insertion files ########################################################
                                                                                         !
    if ( NFILE_LIBRARY > 0 ) then;                                                       !
                                                                                         !
        allocate(NMOLECULE_INSERTION(1:NFILE_LIBRARY));                                  !
                                                                                         !
        allocate(CHNAME_FILE_LIBRARY(1:3,1:NFILE_LIBRARY));                              !
                                                                                         !
        allocate(INSERTION_METHOD_LIBRARY(1:NFILE_LIBRARY));                             !
                                                                                         !
        allocate(ILIBRARY_FILE_INFO(1:NFILE_LIBRARY));                                   !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Update arrays for insertion files ##########################################################
                                                                                         !
    if ( NFILE_LIBRARY > 0 ) then;                                                       !
                                                                                         !
        NMOLECULE_INSERTION(NFILE_LIBRARY) = LAMMPS_INSERT_NMOLECULE(ILOCAL_INSERTS);    !
                                                                                         !
        ILIBRARY_FILE_INFO(NFILE_LIBRARY) = ifile_info;                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Set properties depending on the insertion method ###########################################
                                                                                         !
    select case(TRIM(LAMMPS_INSERT_METHOD(ILOCAL_INSERTS)));                             !
                                                                                         !
        case('com');                                                                     !
                                                                                         !
            INSERTION_METHOD_LIBRARY(1) = 2;                                             !
                                                                                         !
!           ### Allocate the library array containing the center of mass ###########################
                                                                                         !
            allocate(MOLECULE_COM_LIBRARY(1:3,                      &                    !
                                          1:NMOLECULE_INSERTION(1), &                    !
                                          1:NFILE_LIBRARY));                             !
                                                                                         !
!           ### Update the array containing the defined center of mass #############################
                                                                                         !
            MOLECULE_COM_LIBRARY(1:3,1,1) = LAMMPS_INSERT_COM_RG(1:3,ILOCAL_INSERTS);    !
                                                                                         !
        case('rand/region');                                                             !
                                                                                         !
            INSERTION_METHOD_LIBRARY(1) = 3;                                             !
                                                                                         !
!           ### Allocate the library arrays for random insertion in a region #######################
                                                                                         !
            allocate(LIBRARY_REGION_NAME_INSERTION(1:NFILE_LIBRARY));                    !
                                                                                         !
            allocate(LIBRARY_RANDOM_SEED(1:NFILE_LIBRARY));                              !
                                                                                         !
            allocate(LIBRARY_REGION_IRANK(1:NFILE_LIBRARY));                             !
                                                                                         !
!           ### Update the array containing names of regions for random insertion ##################
                                                                                         !
            LIBRARY_REGION_NAME_INSERTION(1) = &                                         !
            TRIM(LAMMPS_INSERT_REGION_NAME(ILOCAL_INSERTS));                             !
                                                                                         !
            LIBRARY_RANDOM_SEED(1) = LAMMPS_INSERT_RANDOM_SEED(ILOCAL_INSERTS);          !
                                                                                         !
            LIBRARY_REGION_IRANK(1) = LAMMPS_INSERT_REGION_IRANK(ILOCAL_INSERTS);        !
                                                                                         !
        case default;                                                                    !
                                                                                         !
            write(icanal,'(a70)') '| Not implemented - stop '// &                        !
                                  REPEAT(' ',44)//'|';                                   !
                                                                                         !
            call CLOSING_PROGRAM(icanal,1);                                              !
                                                                                         !
    end select                                                                           !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate library arrays for cell dimensions ################################################
                                                                                         !
    allocate(CELL_AXIS_LIBRARY(1:3,1:NFILE_LIBRARY_TMP));                                !
                                                                                         !
!   ### Initialization of arrays for cell dimensions ###############################################
                                                                                         !
    CELL_AXIS_LIBRARY(1:3,1:NFILE_LIBRARY_TMP) = 0.0d0;                                  !
                                                                                         !
!   ### Allocate library arrays ####################################################################
                                                                                         !
    allocate(NATOM_LIBRARY(1:NFILE_LIBRARY_TMP));                                        !
                                                                                         !
    allocate(NBOND_LIBRARY(1:NFILE_LIBRARY_TMP));                                        !
                                                                                         !
    allocate(NANGLE_LIBRARY(1:NFILE_LIBRARY_TMP));                                       !
                                                                                         !
    allocate(NDIHEDRAL_LIBRARY(1:NFILE_LIBRARY_TMP));                                    !
                                                                                         !
    allocate(NIMPROPER_LIBRARY(1:NFILE_LIBRARY_TMP));                                    !
                                                                                         !
    allocate(NTYPE_ATOM_LIBRARY(1:NFILE_LIBRARY_TMP));                                   !
                                                                                         !
    allocate(NTYPE_BOND_LIBRARY(1:NFILE_LIBRARY_TMP));                                   !
                                                                                         !
    allocate(NTYPE_ANGLE_LIBRARY(1:NFILE_LIBRARY_TMP));                                  !
                                                                                         !
    allocate(NTYPE_DIHEDRAL_LIBRARY(1:NFILE_LIBRARY_TMP));                               !
                                                                                         ! 
    allocate(NTYPE_IMPROPER_LIBRARY(1:NFILE_LIBRARY_TMP));                               !
                                                                                         !
    allocate(CH_UNITS_STYLE_LIBRARY(1:NFILE_LIBRARY_TMP));                               !
                                                                                         !
    allocate(CH_SPECIAL_BONDS_LIBRARY(1:NFILE_LIBRARY_TMP));                             !
                                                                                         !
    allocate(CH_PAIR_MODIFY_LIBRARY(1:NFILE_LIBRARY_TMP));                               !
                                                                                         !
    allocate(CH_KSPACE_STYLE_LIBRARY(1:NFILE_LIBRARY_TMP));                              !
                                                                                         !
    allocate(IPOTENTIAL_CLASS2_LIBRARY(1:NFILE_LIBRARY_TMP));                            !
                                                                                         !
    allocate(NPAIR_COEFF_CROSS_LIBRARY(1:NFILE_LIBRARY_TMP));                            ! 
                                                                                         !
    allocate(CH_ATOM_STYLE_LIBRARY(1:NFILE_LIBRARY_TMP));                                !
                                                                                         !
    allocate(CH_BOND_STYLE_LIBRARY(1:NFILE_LIBRARY_TMP));                                ! 
                                                                                         !
    allocate(CH_ANGLE_STYLE_LIBRARY(1:NFILE_LIBRARY_TMP));                               !
                                                                                         !
    allocate(CH_DIHEDRAL_STYLE_LIBRARY(1:NFILE_LIBRARY_TMP));                            !
                                                                                         !
    allocate(POTENTIAL_CLASS2_CHTYPE_LIBRARY(1:NFILE_LIBRARY_TMP));                      !
                                                                                         !
    allocate(NPARAM_BONDS_LIBRARY(1:NFILE_LIBRARY_TMP));                                 !
                                                                                         !
    allocate(NPARAM_ANGLES_LIBRARY(1:NFILE_LIBRARY_TMP));                                !
                                                                                         !
    allocate(NPARAM_DIHEDRALS_LIBRARY(1:NFILE_LIBRARY_TMP));                             !
                                                                                         !
    if ( NGenerate_polymer_species > 0 ) then;                                           !
                                                                                         !
        allocate(ILOCR1_LIBRARY(1:NFILE_LIBRARY_TMP));                                   !
                                                                                         !
        allocate(ILOCR2_LIBRARY(1:NFILE_LIBRARY_TMP));                                   !
                                                                                         !
    end if                                                                               !
                                                                                         !
    write(icanal,'(a70)') '| Library arrays were allocated'//REPEAT(' ',38)//'|';        !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialization of library arrays ###########################################################
                                                                                         !
    NATOM_LIBRARY(1:NFILE_LIBRARY_TMP) = 0;                                              !
                                                                                         !
    NBOND_LIBRARY(1:NFILE_LIBRARY_TMP) = 0;                                              !
                                                                                         !
    NANGLE_LIBRARY(1:NFILE_LIBRARY_TMP) = 0;                                             !
                                                                                         !
    NDIHEDRAL_LIBRARY(1:NFILE_LIBRARY_TMP) = 0;                                          !
                                                                                         !
    NIMPROPER_LIBRARY(1:NFILE_LIBRARY_TMP) = 0;                                          !
                                                                                         !
    NTYPE_ATOM_LIBRARY(1:NFILE_LIBRARY_TMP) = 0;                                         !
                                                                                         !
    NTYPE_BOND_LIBRARY(1:NFILE_LIBRARY_TMP) = 0;                                         !
                                                                                         !
    NTYPE_ANGLE_LIBRARY(1:NFILE_LIBRARY_TMP) = 0;                                        !
                                                                                         !
    NTYPE_DIHEDRAL_LIBRARY(1:NFILE_LIBRARY_TMP) = 0;                                     !
                                                                                         ! 
    NTYPE_IMPROPER_LIBRARY(1:NFILE_LIBRARY_TMP) = 0;                                     !
                                                                                         !
    CH_UNITS_STYLE_LIBRARY(1:NFILE_LIBRARY_TMP) = 'XXX';                                 !
                                                                                         !
    CH_SPECIAL_BONDS_LIBRARY(1:NFILE_LIBRARY_TMP) = 'XXX';                               !
                                                                                         !
    CH_PAIR_MODIFY_LIBRARY(1:NFILE_LIBRARY_TMP) = 'XXX';                                 !
                                                                                         !
    CH_KSPACE_STYLE_LIBRARY(1:NFILE_LIBRARY_TMP) = 'XXX';                                !
                                                                                         !
    IPOTENTIAL_CLASS2_LIBRARY(1:NFILE_LIBRARY_TMP) = 0;                                  !
                                                                                         !
    NPAIR_COEFF_CROSS_LIBRARY(1:NFILE_LIBRARY_TMP) = 0;                                  !  
                                                                                         !
    CH_ATOM_STYLE_LIBRARY(1:NFILE_LIBRARY_TMP) = 'XXX';                                  !
                                                                                         !
    CH_BOND_STYLE_LIBRARY(1:NFILE_LIBRARY_TMP) = 'XXX';                                  !   
                                                                                         !
    CH_ANGLE_STYLE_LIBRARY(1:NFILE_LIBRARY_TMP) = 'XXX';                                 !
                                                                                         !
    CH_DIHEDRAL_STYLE_LIBRARY(1:NFILE_LIBRARY_TMP) = 'XXX';                              !
                                                                                         !
    POTENTIAL_CLASS2_CHTYPE_LIBRARY(1:NFILE_LIBRARY_TMP) = 'XXX';                        !
                                                                                         !
!   ### Initialization of parameter arrays #########################################################
                                                                                         !
    NPARAM_BONDS_LIBRARY(1:NFILE_LIBRARY_TMP) = 0;                                       !
                                                                                         !
    NPARAM_ANGLES_LIBRARY(1:NFILE_LIBRARY_TMP) = 0;                                      !
                                                                                         !
    NPARAM_DIHEDRALS_LIBRARY(1:NFILE_LIBRARY_TMP) = 0;                                   !
                                                                                         !
    write(icanal,'(a70)') '| Library arrays were initialized'//REPEAT(' ',36)//'|';      !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialization of variables ################################################################
                                                                                         !
    if ( ILOCAL_INSERTS == 0 ) then;                                                     !
                                                                                         !
        IMAX_NATOM     = 0;                                                              !
                                                                                         !
        IMAX_NBOND     = 0;                                                              !
                                                                                         !
        IMAX_NANGLE    = 0;                                                              !
                                                                                         !
        IMAX_NDIHEDRAL = 0;                                                              !
                                                                                         !
        IMAX_NIMPROPER = 0;                                                              !
                                                                                         !
        IMAX_NTYPE_NATOM    = 0;                                                         !
                                                                                         !
        IMAX_NTYPE_BOND     = 0;                                                         !
                                                                                         !
        IMAX_NTYPE_ANGLE    = 0;                                                         !
                                                                                         !
        IMAX_NTYPE_DIHEDRAL = 0;                                                         !
                                                                                         !
        IMAX_NTYPE_IMPROPER = 0;                                                         !
                                                                                         !
        IMAX_NPAIR_COEFF_CROSS = 0;                                                      ! 
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Set a flag for the file format #############################################################
                                                                                         !
    ilocal_file = 0;                                                                     !
                                                                                         !
    select case(TRIM(LAMMPS_INSERT_CHFILE_FORMAT(ILOCAL_INSERTS)));                      !
                                                                                         !
        case('xyz');                                                                     !
                                                                                         !
            ilocal_file = 1;                                                             !
                                                                                         !
        case('lammps');                                                                  !
                                                                                         !
            ilocal_file = 2;                                                             !
                                                                                         !
        case('padua');                                                                   !
                                                                                         !
            ilocal_file = 2;                                                             !
                                                                                         !
        case('lammps-gaff');                                                             !
                                                                                         !
            ilocal_file = 2;                                                             !
                                                                                         !
        case default;                                                                    !
                                                                                         !
            write(icanal,'(a70)') '| Not implemented - stop '// &                        !
                                  REPEAT(' ',44)//'|';                                   !
                                                                                         !
            call CLOSING_PROGRAM(icanal,1);                                              !
                                                                                         !
    end select                                                                           !
                                                                                         !
    iglobal_file = ilocal_file;                                                          !
                                                                                         !
!   ### Build the file with its extension ##########################################################
                                                                                         !
    LOCAL_CHNAME_FILE = TRIM(LAMMPS_INSERT_CHNAME_FILE(ILOCAL_INSERTS))// &              !
                        TRIM(LAMMPS_INSERT_CHFILE_EXT(1));                               !
                                                                                         !
!   ### Write names of library files to read #######################################################
                                                                                         !
    IOSEF1 = 70 - 2 - 1 - LEN_TRIM(LOCAL_CHNAME_FILE);                                   !
                                                                                         !
    if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                        ! 
                                                                                         !
    write(icanal,'(a70)') '| '//TRIM(LOCAL_CHNAME_FILE)//REPEAT(' ',IOSEF1)//'|';        !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Put the name of the library file in the common array #######################################
                                                                                         !
    CHNAME_FILE_LIBRARY(1,NFILE_LIBRARY) = TRIM(LOCAL_CHNAME_FILE);                      !
                                                                                         !
    CHNAME_FILE_LIBRARY(2,NFILE_LIBRARY) =            &                                  !
    TRIM(LAMMPS_INSERT_CHNAME_FILE(ILOCAL_INSERTS))// &                                  !
    TRIM(LAMMPS_INSERT_CHFILE_EXT(2));                                                   !
                                                                                         !
    CHNAME_FILE_LIBRARY(3,NFILE_LIBRARY) =            &                                  !
    TRIM(LAMMPS_INSERT_CHNAME_FILE(ILOCAL_INSERTS))// &                                  !
    TRIM(LAMMPS_INSERT_CHFILE_EXT(3));                                                   !
                                                                                         !
!   ### Initialization of local variables containing properties of the molecule ####################
                                                                                         !
    EOF = 0;                                                                             ! 
                                                                                         !
    IOSEF_NATOM     = 0;                                                                 !
                                                                                         ! 
    IOSEF_NBOND     = 0;                                                                 !  
                                                                                         !
    IOSEF_NANGLE    = 0;                                                                 !
                                                                                         !
    IOSEF_NDIHEDRAL = 0;                                                                 !
                                                                                         !
    IOSEF_NIMPROPER = 0;                                                                 !
                                                                                         !
    IOSEF_NTYPE_NATOM    = 0;                                                            !
                                                                                         !
    IOSEF_NTYPE_BOND     = 0;                                                            !
                                                                                         !
    IOSEF_NTYPE_ANGLE    = 0;                                                            !
                                                                                         !
    IOSEF_NTYPE_DIHEDRAL = 0;                                                            !
                                                                                         !
    IOSEF_NTYPE_IMPROPER = 0;                                                            !
                                                                                         !
    IOSEF_NPAIR_COEFF_CROSS = 0;                                                         !
                                                                                         !
    LOCAL_ATOM_LABEL(1:100) = 'XXX';                                                     !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Read library files #########################################################################
                                                                                         !
    select case(ilocal_file);                                                            !
                                                                                         !
        case(1);                                                                         ! Data in the library are in the xyz format
                                                                                         !
!           ### Open the library file ##############################################################
                                                                                         !
            open(2,file=TRIM(LOCAL_CHNAME_FILE));                                        !
                                                                                         !
!           ### Read properties of the library file ################################################
                                                                                         !
            read(2,*) IOSEF_NATOM;                                                       !
                                                                                         !
            read(2,*);                                                                   !
                                                                                         !
!           ### Find atom labels ###################################################################
                                                                                         !
            do j = 1, IOSEF_NATOM;                                                       !
                                                                                         !
                read(2,*) CHOSEF1;                                                       !
                                                                                         !
                if ( IOSEF_NTYPE_NATOM == 0 ) then;                                      !
                                                                                         !
                    IOSEF_NTYPE_NATOM = IOSEF_NTYPE_NATOM + 1;                           !
                                                                                         !
                    LOCAL_ATOM_LABEL(IOSEF_NTYPE_NATOM) = TRIM(CHOSEF1);                 !
                                                                                         !
                else                                                                     !
                                                                                         !
                    IFOUND = 0;                                                          !
                                                                                         !
                    do k = 1, IOSEF_NTYPE_NATOM;                                         !
                                                                                         !
                        if ( TRIM(LOCAL_ATOM_LABEL(k)) /= TRIM(CHOSEF1) ) CYCLE;         !
                                                                                         ! 
                        IFOUND = 1;                                                      !
                                                                                         !
                        EXIT;                                                            !
                                                                                         !
                    end do                                                               ! 
                                                                                         !
                    if ( IFOUND == 1 ) CYCLE;                                            !
                                                                                         !
                    IOSEF_NTYPE_NATOM = IOSEF_NTYPE_NATOM + 1;                           !
                                                                                         !
                    LOCAL_ATOM_LABEL(IOSEF_NTYPE_NATOM) = TRIM(CHOSEF1);                 !
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
!           ### Closing the library file ###########################################################
                                                                                         !
            close(2);                                                                    !
                                                                                         !
        case(2);                                                                         ! Data in the library are in the lammps format
                                                                                         !
!           ### Open the library file ##############################################################
                                                                                         !
            open(2,file=TRIM(LOCAL_CHNAME_FILE));                                        !
                                                                                         !
!           ### Read properties of the molecular configuration from the library  ###################
                                                                                         !
            do                                                                           !
                                                                                         !
                read(2,'(a)',iostat=EOF) CHARLINE;                                       !
                                                                                         !
                if ( EOF /= 0 ) EXIT;                                                    !
                                                                                         !
                if ( INDEX(CHARLINE,'atoms') > 0 ) then;                                 !
                                                                                         !
                    read(CHARLINE,*) IOSEF_NATOM;                                        !
                                                                                         !
                else if ( INDEX(CHARLINE,'bonds') > 0 ) then;                            !
                                                                                         !
                    read(CHARLINE,*) IOSEF_NBOND;                                        !
                                                                                         !
                else if ( INDEX(CHARLINE,'angles') > 0 ) then;                           !
                                                                                         !
                    read(CHARLINE,*) IOSEF_NANGLE;                                       !
                                                                                         !
                else if ( INDEX(CHARLINE,'dihedrals') > 0 ) then;                        !
                                                                                         !
                    read(CHARLINE,*) IOSEF_NDIHEDRAL;                                    !
                                                                                         !
                else if ( INDEX(CHARLINE,'impropers') > 0 ) then;                        !
                                                                                         !
                    read(CHARLINE,*) IOSEF_NIMPROPER;                                    !
                                                                                         !
                else if ( INDEX(CHARLINE,'Types') > 0 ) then;                            !
                                                                                         !
                    read(2,*);                                                           !
                                                                                         !
                    do j = 1, IOSEF_NATOM;                                               !
                                                                                         !
                        read(2,*) IOSEF1, IOSEF2;                                        !
                                                                                         !
                        if ( IOSEF2 > IOSEF_NTYPE_NATOM ) IOSEF_NTYPE_NATOM = IOSEF2;    !
                                                                                         !
                    end do                                                               !
                                                                                         !
                else if ( INDEX(CHARLINE,'Bonds') > 0 ) then;                            !
                                                                                         !
                    read(2,*);                                                           !
                                                                                         !
                    do j = 1, IOSEF_NBOND;                                               !
                                                                                         !
                        read(2,*) IOSEF1, IOSEF2;                                        !
                                                                                         !
                        if ( IOSEF2 > IOSEF_NTYPE_BOND ) IOSEF_NTYPE_BOND = IOSEF2;      !
                                                                                         !
                    end do                                                               !
                                                                                         !
                else if ( INDEX(CHARLINE,'Angles') > 0 ) then;                           !
                                                                                         !
                    read(2,*);                                                           !
                                                                                         !
                    do j = 1, IOSEF_NANGLE;                                              !
                                                                                         !
                        read(2,*) IOSEF1, IOSEF2;                                        !
                                                                                         !
                        if ( IOSEF2 > IOSEF_NTYPE_ANGLE ) IOSEF_NTYPE_ANGLE = IOSEF2;    !
                                                                                         !
                    end do                                                               !
                                                                                         !
                else if ( INDEX(CHARLINE,'Dihedrals') > 0 ) then;                        !
                                                                                         !
                    read(2,*);                                                           !
                                                                                         !
                    do j = 1, IOSEF_NDIHEDRAL;                                           !
                                                                                         !
                        read(2,*) IOSEF1, IOSEF2;                                        !
                                                                                         !
                        if ( IOSEF2 > IOSEF_NTYPE_DIHEDRAL ) then;                       !
                                                                                         !
                            IOSEF_NTYPE_DIHEDRAL = IOSEF2;                               !
                                                                                         !
                        end if                                                           !
                                                                                         !
                    end do                                                               !
                                                                                         !
                else if ( INDEX(CHARLINE,'Impropers') > 0 ) then;                        !
                                                                                         !
                    read(2,*);                                                           !
                                                                                         !
                    do j = 1, IOSEF_NIMPROPER;                                           !
                                                                                         !
                        read(2,*) IOSEF1, IOSEF2;                                        !
                                                                                         !
                        if ( IOSEF2 > IOSEF_NTYPE_IMPROPER ) then;                       !
                                                                                         !
                            IOSEF_NTYPE_IMPROPER = IOSEF2;                               !
                                                                                         !
                        end if                                                           !
                                                                                         !
                    end do                                                               !
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
!           ### Closing the library file ###########################################################
                                                                                         !
            close(2);                                                                    !
                                                                                         !
        case default;                                                                    !
                                                                                         !
            write(icanal,'(a70)') '| Not implemented - stop '// &                        !
                                  REPEAT(' ',44)//'|';                                   !
                                                                                         !
            call CLOSING_PROGRAM(icanal,1);                                              !
                                                                                         !
    end select                                                                           !
                                                                                         !
!   ### Update global arrays #######################################################################
                                                                                         !
    if ( IOSEF_NATOM > IMAX_NATOM ) IMAX_NATOM = IOSEF_NATOM;                            !
                                                                                         !
    if ( IOSEF_NBOND > IMAX_NBOND ) IMAX_NBOND = IOSEF_NBOND;                            !
                                                                                         !
    if ( IOSEF_NANGLE > IMAX_NANGLE ) IMAX_NANGLE = IOSEF_NANGLE;                        !
                                                                                         !
    if ( IOSEF_NDIHEDRAL > IMAX_NDIHEDRAL ) IMAX_NDIHEDRAL = IOSEF_NDIHEDRAL;            !
                                                                                         !
    if ( IOSEF_NIMPROPER > IMAX_NIMPROPER ) IMAX_NIMPROPER = IOSEF_NIMPROPER;            !
                                                                                         !
    if ( IOSEF_NTYPE_NATOM > IMAX_NTYPE_NATOM ) IMAX_NTYPE_NATOM = IOSEF_NTYPE_NATOM;    !
                                                                                         !
    if ( IOSEF_NTYPE_BOND > IMAX_NTYPE_BOND ) IMAX_NTYPE_BOND = IOSEF_NTYPE_BOND;        !
                                                                                         !
    if ( IOSEF_NTYPE_ANGLE > IMAX_NTYPE_ANGLE ) IMAX_NTYPE_ANGLE = IOSEF_NTYPE_ANGLE;    !
                                                                                         !
    if ( IOSEF_NTYPE_DIHEDRAL > IMAX_NTYPE_DIHEDRAL ) then;                              !
                                                                                         !
        IMAX_NTYPE_DIHEDRAL = IOSEF_NTYPE_DIHEDRAL;                                      !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( IOSEF_NTYPE_IMPROPER > IMAX_NTYPE_IMPROPER ) then;                              !
                                                                                         !
        IMAX_NTYPE_IMPROPER = IOSEF_NTYPE_IMPROPER;                                      !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Write properties of the current library file ###############################################
                                                                                         !
    if ( IOSEF_NATOM > 0 ) then;                                                         !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| IOSEF_NATOM          : ', &                      !
                                     IOSEF_NATOM,                 &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( IOSEF_NBOND > 0 ) then;                                                         !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| IOSEF_NBOND          : ', &                      !
                                     IOSEF_NBOND,                 &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( IOSEF_NANGLE > 0 ) then;                                                        !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| IOSEF_NANGLE         : ', &                      !
                                     IOSEF_NANGLE,                &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( IOSEF_NDIHEDRAL > 0 ) then;                                                     !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| IOSEF_NDIHEDRAL      : ', &                      !
                                     IOSEF_NDIHEDRAL,             &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( IOSEF_NIMPROPER > 0 ) then;                                                     !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| IOSEF_NIMPROPER      : ', &                      !
                                     IOSEF_NIMPROPER,             &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( IOSEF_NTYPE_NATOM > 0 ) then;                                                   !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| IOSEF_NTYPE_NATOM    : ', &                      !
                                     IOSEF_NTYPE_NATOM,           &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( IOSEF_NTYPE_BOND > 0 ) then;                                                    !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| IOSEF_NTYPE_BOND     : ', &                      !
                                     IOSEF_NTYPE_BOND,            &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( IOSEF_NTYPE_ANGLE > 0 ) then;                                                   !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| IOSEF_NTYPE_ANGLE    : ', &                      !
                                     IOSEF_NTYPE_ANGLE,           &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( IOSEF_NTYPE_DIHEDRAL > 0 ) then;                                                !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| IOSEF_NTYPE_DIHEDRAL : ', &                      !
                                     IOSEF_NTYPE_DIHEDRAL,        &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( IOSEF_NTYPE_IMPROPER > 0 ) then;                                                !
                                                                                         ! 
        write(icanal,'(a25,i8,a37)') '| IOSEF_NTYPE_IMPROPER : ', &                      !
                                     IOSEF_NTYPE_IMPROPER,        &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Read molecule properties in the library files  with the lammps format ######################
                                                                                         !
!   if ( ( TRIM(CHFILE_FORMAT) == 'lammps'      ) .OR.  &                                !
!        ( TRIM(CHFILE_FORMAT) == 'padua'       ) .OR.  &                                !
!        ( TRIM(CHFILE_FORMAT) == 'lammps-gaff' ) ) then;                                !
                                                                                         !
!       if ( NFILE_LIBRARY > 0 ) then;                                                   !
                                                                                         !
!           do i = 1, NFILE_LIBRARY;                                                     ! LOOP OVER THE NUMBER OF FILES TO READ IN THE LIBRARY
                                                                                         !
!               CHEXT = TRIM(CHNAME_FILE_LIBRARY(i))//'.template';                       !
                                                                                         !
!               IOSEF1 = LEN_TRIM(CHEXT);                                                !
                                                                                         !
!               IOSEF2 = 67 - IOSEF1;                                                    !
                                                                                         !
!               if ( IOSEF2 < 0 ) IOSEF2 = 1;                                            !
                                                                                         !
!               write(icanal,'(a70)') '| '//TRIM(CHEXT)//REPEAT(' ',IOSEF2)//'|';        !
                                                                                         !
!               write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
!           end do                                                                       !
                                                                                         !
!       end if                                                                           !
                                                                                         !
!       ### Read library files of molecules to be removed in the molecular configuration ###############
                                                                                         !
!       if ( NFILE_LIBRARY_REMOVE > 0 ) then;                                            !  
                                                                                         !
!           do i = 1, NFILE_LIBRARY_REMOVE;                                              ! LOOP OVER THE NUMBER OF FILES TO READ IN THE LIBRARY
                                                                                         !
!               CHEXT = TRIM(CHNAME_FILE_LIBRARY_REMOVE(i))//'.template';                !
                                                                                         !
!               IOSEF1 = LEN_TRIM(CHEXT);                                                !
                                                                                         !
!               IOSEF2 = 67 - IOSEF1;                                                    !
                                                                                         !
!               write(icanal,'(a70)') '| '//TRIM(CHEXT)//REPEAT(' ',IOSEF2)//'|';        !
                                                                                         !
!               write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               open(2,file=TRIM(CHEXT));                                                !
                                                                                         !
!               EOF = 0;                                                                 !
                                                                                         !
!               IOSEF_NATOM     = 0;                                                     !
!               IOSEF_NBOND     = 0;                                                     !
!               IOSEF_NANGLE    = 0;                                                     !
!               IOSEF_NDIHEDRAL = 0;                                                     !
!               IOSEF_NIMPROPER = 0;                                                     !

!               IOSEF_NTYPE_NATOM    = 0;
!               IOSEF_NTYPE_BOND     = 0;
!               IOSEF_NTYPE_ANGLE    = 0;
!               IOSEF_NTYPE_DIHEDRAL = 0;
!               IOSEF_NTYPE_IMPROPER = 0;

!               do

!                   read(2,'(a)',iostat=EOF) CHARLINE;

!                   if ( EOF /= 0 ) EXIT;

!                   if ( INDEX(CHARLINE,'atoms') > 0 ) then;

!                       read(CHARLINE,*) IOSEF_NATOM;

!                       if ( IOSEF_NATOM > IMAX_NATOM ) IMAX_NATOM = IOSEF_NATOM;

!                   else if ( INDEX(CHARLINE,'bonds') > 0 ) then;

!                       read(CHARLINE,*) IOSEF_NBOND;

!                       if ( IOSEF_NBOND > IMAX_NBOND ) IMAX_NBOND = IOSEF_NBOND;

!                   else if ( INDEX(CHARLINE,'angles') > 0 ) then;

!                       read(CHARLINE,*) IOSEF_NANGLE;

!                       if ( IOSEF_NANGLE > IMAX_NANGLE ) IMAX_NANGLE = IOSEF_NANGLE;

!                   else if ( INDEX(CHARLINE,'dihedrals') > 0 ) then;

!                       read(CHARLINE,*) IOSEF_NDIHEDRAL;

!                       if ( IOSEF_NDIHEDRAL > IMAX_NDIHEDRAL ) IMAX_NDIHEDRAL = IOSEF_NDIHEDRAL;

!                   else if ( INDEX(CHARLINE,'impropers') > 0 ) then;

!                       read(CHARLINE,*) IOSEF_NIMPROPER;

!                       if ( IOSEF_NIMPROPER > IMAX_NIMPROPER ) IMAX_NIMPROPER = IOSEF_NIMPROPER;

!                   else if ( INDEX(CHARLINE,'Types') > 0 ) then;

!                       read(2,*);

!                       do j = 1, IOSEF_NATOM;

!                           read(2,*) IOSEF1, IOSEF2;

!                           if ( IOSEF2 > IOSEF_NTYPE_NATOM ) IOSEF_NTYPE_NATOM = IOSEF2;

!                       end do

!                       if ( IOSEF_NTYPE_NATOM > IMAX_NTYPE_NATOM ) IMAX_NTYPE_NATOM = IOSEF_NTYPE_NATOM;

!                   else if ( INDEX(CHARLINE,'Bonds') > 0 ) then;

!                       read(2,*);

!                       do j = 1, IOSEF_NBOND;

!                           read(2,*) IOSEF1, IOSEF2;

!                           if ( IOSEF2 > IOSEF_NTYPE_BOND ) IOSEF_NTYPE_BOND = IOSEF2;

!                       end do

!                       if ( IOSEF_NTYPE_BOND > IMAX_NTYPE_BOND ) IMAX_NTYPE_BOND = IOSEF_NTYPE_BOND;

!                   else if ( INDEX(CHARLINE,'Angles') > 0 ) then;

!                       read(2,*);

!                       do j = 1, IOSEF_NANGLE;

!                           read(2,*) IOSEF1, IOSEF2;

!                           if ( IOSEF2 > IOSEF_NTYPE_ANGLE ) IOSEF_NTYPE_ANGLE = IOSEF2;

!                       end do
 
!                       if ( IOSEF_NTYPE_ANGLE > IMAX_NTYPE_ANGLE ) IMAX_NTYPE_ANGLE = IOSEF_NTYPE_ANGLE;

!                   else if ( INDEX(CHARLINE,'Dihedrals') > 0 ) then;

!                       read(2,*);

!                       do j = 1, IOSEF_NDIHEDRAL;

!                           read(2,*) IOSEF1, IOSEF2;

!                           if ( IOSEF2 > IOSEF_NTYPE_DIHEDRAL ) IOSEF_NTYPE_DIHEDRAL = IOSEF2;

!                       end do

!                       if ( IOSEF_NTYPE_DIHEDRAL > IMAX_NTYPE_DIHEDRAL ) IMAX_NTYPE_DIHEDRAL = IOSEF_NTYPE_DIHEDRAL;

!                   else if ( INDEX(CHARLINE,'Impropers') > 0 ) then;

!                       read(2,*);

!                       do j = 1, IOSEF_NIMPROPER;

!                           read(2,*) IOSEF1, IOSEF2;

!                           if ( IOSEF2 > IOSEF_NTYPE_IMPROPER ) IOSEF_NTYPE_IMPROPER = IOSEF2;

!                       end do

!                       if ( IOSEF_NTYPE_IMPROPER > IMAX_NTYPE_IMPROPER ) IMAX_NTYPE_IMPROPER = IOSEF_NTYPE_IMPROPER;

!                   end if

!               end do

!               close(2);
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
!               write(icanal,'(a25,i8,a37)') '| IOSEF_NATOM          : ', &              !
!                                            IOSEF_NATOM,                 &              !
!                                            REPEAT(' ',36)//'|';                        !
                                                                                         !
!               write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

!               write(icanal,'(a25,i8,a37)') '| IOSEF_NBOND          : ', &
!                                            IOSEF_NBOND,                 &
!                                            REPEAT(' ',36)//'|';

!               write(icanal,'(a25,i8,a37)') '| IOSEF_NANGLE         : ', &
!                                            IOSEF_NANGLE,                &
!                                            REPEAT(' ',36)//'|';

!               write(icanal,'(a25,i8,a37)') '| IOSEF_NDIHEDRAL      : ', &
!                                            IOSEF_NDIHEDRAL,             &
!                                            REPEAT(' ',36)//'|';

!               write(icanal,'(a25,i8,a37)') '| IOSEF_NIMPROPER      : ', &
!                                            IOSEF_NIMPROPER,             &
!                                            REPEAT(' ',36)//'|';

!               write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

!               write(icanal,'(a25,i8,a37)') '| IOSEF_NTYPE_NATOM    : ', &
!                                            IOSEF_NTYPE_NATOM,           &
!                                            REPEAT(' ',36)//'|';

!               write(icanal,'(a25,i8,a37)') '| IOSEF_NTYPE_BOND     : ', &
!                                            IOSEF_NTYPE_BOND,            &
!                                            REPEAT(' ',36)//'|';

!               write(icanal,'(a25,i8,a37)') '| IOSEF_NTYPE_ANGLE    : ', &
!                                            IOSEF_NTYPE_ANGLE,           &
!                                            REPEAT(' ',36)//'|';

!               write(icanal,'(a25,i8,a37)') '| IOSEF_NTYPE_DIHEDRAL : ', &
!                                            IOSEF_NTYPE_DIHEDRAL,        &
!                                            REPEAT(' ',36)//'|';

!               write(icanal,'(a25,i8,a37)') '| IOSEF_NTYPE_IMPROPER : ', &
!                                            IOSEF_NTYPE_IMPROPER,        &
!                                            REPEAT(' ',36)//'|';

!               write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

!           end do

!       end if                                                                           !
                                                                                         !


!   else 
                                                                                         !
!      write(icanal,*) 'Not implemented - stop!';                                        !
                                                                                         !
!      stop; !///////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   end if                                                                               !
                                                                                         !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Read monomer properties in the library files ###############################################
                                                                                         !
!   if ( NGenerate_polymer_species > 0 ) then;                                           !
                                                                                         !
!       do i = 1, NTotal_monomers;                                                       !
                                                                                         !
!           CHEXT = TRIM(Chname_file_library_monomer(i))//'.xyz';                        !
                                                                                         !
!           IOSEF1 = LEN_TRIM(CHEXT);                                                    !
                                                                                         !
!           IOSEF2 = 67 - IOSEF1;                                                        !
                                                                                         !
!           if ( IOSEF2 < 0 ) IOSEF2 = 1;                                                !
                                                                                         !
!           write(icanal,'(a70)') '| '//TRIM(CHEXT)//REPEAT(' ',IOSEF2)//'|';            !
                                                                                         !
!           write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           EOF = 0;                                                                     !
                                                                                         !
!           IOSEF_NATOM     = 0;                                                         !
                                                                                         !
!           IOSEF_NBOND     = 0;                                                         !
                                                                                         !
!           IOSEF_NANGLE    = 0;                                                         !
                                                                                         !
!           IOSEF_NDIHEDRAL = 0;                                                         !
                                                                                         !
!           IOSEF_NIMPROPER = 0;                                                         !
                                                                                         !
!           IOSEF_NTYPE_NATOM    = 0;                                                    !
                                                                                         !
!           IOSEF_NTYPE_BOND     = 0;                                                    !
                                                                                         !
!           IOSEF_NTYPE_ANGLE    = 0;                                                    !
                                                                                         !
!           IOSEF_NTYPE_DIHEDRAL = 0;                                                    !
                                                                                         !
!           IOSEF_NTYPE_IMPROPER = 0;                                                    !
                                                                                         !
!           LOCAL_ATOM_LABEL(1:100) = 'XXX';                                             !
                                                                                         !
!           open(2,file=TRIM(CHEXT));                                                    !
                                                                                         !
!           read(2,*) IOSEF_NATOM;                                                       ! 
                                                                                         !
!           read(2,*);                                                                   !
                                                                                         !
!           do j = 1, IOSEF_NATOM;                                                       !
                                                                                         !
!               read(2,*) CHOSEF1;                                                       !
                                                                                         !
!               if ( IOSEF_NTYPE_NATOM == 0 ) then;                                      !

!                   IOSEF_NTYPE_NATOM = IOSEF_NTYPE_NATOM + 1;

!                   LOCAL_ATOM_LABEL(IOSEF_NTYPE_NATOM) = TRIM(CHOSEF1);

!               else

!                   IFOUND = 0;

!                   do k = 1, IOSEF_NTYPE_NATOM;                                         !
                                                                                         !
!                       if ( TRIM(LOCAL_ATOM_LABEL(k)) == TRIM(CHOSEF1) ) then;          !
                                                                                         !
!                           IFOUND = 1;                                                  !
                                                                                         !
!                           EXIT;                                                        !
                                                                                         !
!                       end if                                                           !
                                                                                         !
!                   end do                                                               !
                                                                                         !                            
!                   if ( IFOUND == 1 ) CYCLE;                                            !
                                                                                         !
!                   IOSEF_NTYPE_NATOM = IOSEF_NTYPE_NATOM + 1;                           !
                                                                                         !
!                   LOCAL_ATOM_LABEL(IOSEF_NTYPE_NATOM) = TRIM(CHOSEF1);                 !
                                                                                         !
!               end if                                                                   !
                                                                                         !
!           end do                                                                       !
                                                                                         !
!           close(2);                                                                    !
                                                                                         !
!           if ( IOSEF_NATOM > IMAX_NATOM ) IMAX_NATOM = IOSEF_NATOM;                    !
                                                                                         !
!           write(icanal,'(a25,i8,a37)') '| IOSEF_NATOM          : ', &                  !
!                                        IOSEF_NATOM,                 &                  !
!                                        REPEAT(' ',36)//'|';                            !
                                                                                         !
!           write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           if ( IOSEF_NTYPE_NATOM > IMAX_NTYPE_NATOM ) then;                            !
                                                                                         !
!               IMAX_NTYPE_NATOM = IOSEF_NTYPE_NATOM;                                    !
                                                                                         !
!           end if                                                                       !
                                                                                         !
!           write(icanal,'(a25,i8,a37)') '| IOSEF_NTYPE_NATOM    : ', &                  !
!                                        IOSEF_NTYPE_NATOM,           &                  !
!                                        REPEAT(' ',36)//'|';                            !
                                                                                         !
!           write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!       end do                                                                           !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocation of library arrays ###############################################################
                                                                                         !
    if ( ( IMAX_NATOM > 0 ) .AND. ( NFILE_LIBRARY_TMP > 0 ) ) then;                      !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| IMAX_NATOM           : ', &                      !
                                     IMAX_NATOM,                  &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        allocate(CONFIG_ATOMID_LIBRARY(1:IMAX_NATOM,1:NFILE_LIBRARY_TMP));               !
                                                                                         !
        allocate(CONFIG_ATOM_TYPE_LIBRARY(1:IMAX_NATOM,1:NFILE_LIBRARY_TMP));            !
                                                                                         !
        allocate(CONFIG_QI_LIBRARY(1:IMAX_NATOM,1:NFILE_LIBRARY_TMP));                   !
                                                                                         !
        allocate(CONFIG_VI_LIBRARY(1:3,1:IMAX_NATOM,1:NFILE_LIBRARY_TMP));               !
                                                                                         !
        allocate(CONFIG_RI_LIBRARY(1:3,1:IMAX_NATOM,1:NFILE_LIBRARY_TMP));               !
                                                                                         !
        allocate(CONFIG_NAT_LIBRARY(1:IMAX_NATOM,1:NFILE_LIBRARY_TMP));                  !
                                                                                         !
        write(icanal,'(a70)') '| Library arrays with the size '// &                      !
                              'of NATOM were allocated'//         &                      !
                              REPEAT(' ',15)//'|';                                       !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( ( IMAX_NBOND > 0 ) .AND. ( NFILE_LIBRARY_TMP > 0 ) ) then;                      !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| IMAX_NBOND           : ', &                      !
                                     IMAX_NBOND,                  &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        allocate(BOND_TYPE_LIBRARY(1:IMAX_NBOND,1:NFILE_LIBRARY_TMP));                   !
                                                                                         !
        allocate(BOND_ATOMID_LIBRARY(1:2,1:IMAX_NBOND,1:NFILE_LIBRARY_TMP));             !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    if ( ( IMAX_NANGLE > 0 ) .AND. ( NFILE_LIBRARY_TMP > 0 ) ) then;                     !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| IMAX_NANGLE          : ', &                      !
                                     IMAX_NANGLE,                 &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        allocate(ANGLE_TYPE_LIBRARY(1:IMAX_NANGLE,1:NFILE_LIBRARY_TMP));                 !
                                                                                         !
        allocate(ANGLE_ATOMID_LIBRARY(1:3,1:IMAX_NANGLE,1:NFILE_LIBRARY_TMP));           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    if ( ( IMAX_NDIHEDRAL > 0 ) .AND. ( NFILE_LIBRARY_TMP > 0 ) ) then;                  !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| IMAX_NDIHEDRAL       : ', &                      !
                                     IMAX_NDIHEDRAL,              &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        allocate(DIHEDRAL_TYPE_LIBRARY(1:IMAX_NDIHEDRAL,1:NFILE_LIBRARY_TMP));           !
                                                                                         !
        allocate(DIHEDRAL_ATOMID_LIBRARY(1:4,1:IMAX_NDIHEDRAL,1:NFILE_LIBRARY_TMP));     !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( ( IMAX_NIMPROPER > 0 ) .AND. ( NFILE_LIBRARY_TMP > 0 ) ) then;                  ! 
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| IMAX_NIMPROPER       : ', &                      !
                                     IMAX_NIMPROPER,              &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         ! 
        allocate(IMPROPER_TYPE_LIBRARY(1:IMAX_NIMPROPER,1:NFILE_LIBRARY_TMP));           !
                                                                                         !
        allocate(IMPROPER_ATOMID_LIBRARY(1:4,1:IMAX_NIMPROPER,1:NFILE_LIBRARY_TMP));     !
                                                                                         !
    end if                                                                               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    if ( ( IMAX_NTYPE_NATOM > 0 ) .AND. ( NFILE_LIBRARY_TMP > 0 ) ) then;                !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| IMAX_NTYPE_NATOM     : ', &                      !
                                     IMAX_NTYPE_NATOM,            &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        allocate(ATOM_MASSES_LIBRARY(1:IMAX_NTYPE_NATOM,1:NFILE_LIBRARY_TMP));           !
                                                                                         !
        allocate(ATOM_LABEL_LIBRARY(1:IMAX_NTYPE_NATOM,1:NFILE_LIBRARY_TMP));            !
                                                                                         !
        allocate(POTENTIAL_CLASS2_LIBRARY(1:2,1:IMAX_NTYPE_NATOM,1:NFILE_LIBRARY_TMP));  !
                                                                                         !
        allocate(PAIR_COEFF_CROSS_LIBRARY(1:2,1:1000,1:NFILE_LIBRARY_TMP));              !
                                                                                         !
        allocate(PAIR_ATOMID_CROSS_LIBRARY(1:2,1:1000,1:NFILE_LIBRARY_TMP));             !
                                                                                         !
        write(icanal,'(a70)') '| Library arrays with the size of '// &                   !
                              'NTYPE_NATOM were allocated'//         &                   !
                              REPEAT(' ',9)//'|';                                        !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( ( IMAX_NTYPE_BOND > 0 ) .AND. ( NFILE_LIBRARY_TMP > 0 ) ) then;                 !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| IMAX_NTYPE_BOND      : ', &                      !
                                     IMAX_NTYPE_BOND,             &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        allocate(BOND_COEFFS_LIBRARY(1:4,1:IMAX_NTYPE_BOND,1:NFILE_LIBRARY_TMP));        !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( ( IMAX_NTYPE_ANGLE > 0 ) .AND. ( NFILE_LIBRARY_TMP > 0 ) ) then;                !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| IMAX_NTYPE_ANGLE     : ', &                      !
                                     IMAX_NTYPE_ANGLE,            &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        allocate(ANGLE_COEFFS_LIBRARY(1:4,1:IMAX_NTYPE_ANGLE,1:NFILE_LIBRARY_TMP));      !
                                                                                         !
        allocate(BONDBOND_COEFFS_LIBRARY(1:3,1:IMAX_NTYPE_ANGLE,1:NFILE_LIBRARY_TMP));   !
                                                                                         !
        allocate(BONDANGLE_COEFFS_LIBRARY(1:4,1:IMAX_NTYPE_ANGLE,1:NFILE_LIBRARY_TMP));  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( ( IMAX_NTYPE_DIHEDRAL > 0 ) .AND. ( NFILE_LIBRARY_TMP > 0 ) ) then;             !             
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| IMAX_NTYPE_DIHEDRAL  : ', &                      !
                                     IMAX_NTYPE_DIHEDRAL,         &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        allocate(DIHEDRAL_COEFFS_LIBRARY(1:6,1:IMAX_NTYPE_DIHEDRAL,1:NFILE_LIBRARY_TMP));!
                                                                                         !
        allocate(MIDDLEBONDTORSION_COEFFS_LIBRARY(1:4,1:IMAX_NTYPE_DIHEDRAL, &           !
                                                  1:NFILE_LIBRARY_TMP));                 !
                                                                                         !
        allocate(ENDBONDTORSION_COEFFS_LIBRARY(1:8,1:IMAX_NTYPE_DIHEDRAL, &              !
                                               1:NFILE_LIBRARY_TMP));                    !
                                                                                         !
        allocate(ANGLETORSION_COEFFS_LIBRARY(1:8,1:IMAX_NTYPE_DIHEDRAL, &                !
                                             1:NFILE_LIBRARY_TMP));                      !
                                                                                         !
        allocate(ANGLEANGLETORSION_COEFFS_LIBRARY(1:3,1:IMAX_NTYPE_DIHEDRAL, &           !
                                                  1:NFILE_LIBRARY_TMP));                 !   
                                                                                         !
        allocate(BONDBOND13_COEFFS_LIBRARY(1:3,1:IMAX_NTYPE_DIHEDRAL, &                  !
                                           1:NFILE_LIBRARY_TMP));                        !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( ( IMAX_NTYPE_IMPROPER > 0 ) .AND. ( NFILE_LIBRARY_TMP > 0 ) ) then;             !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| IMAX_NTYPE_IMPROPER  : ', &                      !
                                     IMAX_NTYPE_IMPROPER,         &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        allocate(IMPROPER_COEFFS_LIBRARY(1:2,1:IMAX_NTYPE_IMPROPER, &                    !
                                         1:NFILE_LIBRARY_TMP));                          !
                                                                                         !
        allocate(ANGLEANGLE_COEFFS_LIBRARY(1:6,                   &                      !
                                           1:IMAX_NTYPE_IMPROPER, &                      !
                                           1:NFILE_LIBRARY_TMP));                        !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( IMAX_NPAIR_COEFF_CROSS == 0 ) then;                                             !
                                                                                         !
        NPAIR_COEFF_CROSS_LIBRARY = 0;                                                   !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   if ( ILOCAL_INSERTS == 3 ) stop; !///////////////////////////////////////////////////!
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine ALLOCATE_INSERT_LIBRARY_ARRAYS
