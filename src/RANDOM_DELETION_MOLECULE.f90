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

subroutine RANDOM_DELETION_MOLECULE(icanal, & 
                                    PASSA,  &
                                    PASSB) 

!   ************************************************************************************************
!   **                                  READ THE INPUT FILE                                       **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                    : CANAL ON WHICH OUTPUT DATA ARE WRITTEN                         **
!   ** PASSA                     :  **
!   ** PASSB                     :  **
!   **                                                                                            **
!   ************************************************************************************************

!   use data_in;
    use module_physical_constants;

    use module_size_arrays;

    use module_library;

    use module_config;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    real (kind=8), dimension(1:3,1:3), intent(in) :: PASSA, PASSB;

!   ************************************************************************************************

    integer (kind=4) :: i, j, k, m, n;

    integer (kind=4) :: IFOUND;

    integer (kind=4) :: IINSERTED, IATOM, IOVERLAPPING;

    integer (kind=4) :: IATOM_TYPE_TRANSLATE, IMOLECULEID;

    integer (kind=4) :: IBOND_TYPE_TRANSLATE, IANGLE_TYPE_TRANSLATE, IIMPROPER_TYPE_TRANSLATE;

    integer (kind=4) :: IDIHEDRAL_TYPE_TRANSLATE;

    integer (kind=4) :: MAX_MOLECULEID;

    integer (kind=4), dimension(1:1000) :: NCOUNT_DELETION;

    integer (kind=4), dimension(1:1000,1:1000) :: MOLECULEID_DELETION;

    integer (kind=4) :: CHECK_NTYPE_ATOM, NLIST_ATOM_TYPE;

    integer (kind=4) :: CHECK_NTYPE_BOND, NLIST_BOND_TYPE;

    integer (kind=4), allocatable, dimension(:) :: LIST_ATOM_TYPE,      &
                                                   LIST_BOND_TYPE,      &
                                                   MOLECULE_NTYPE_ATOM, &
                                                   MOLECULE_NTYPE_BOND;

    real (kind=8) :: grnd;

    real (kind=8) :: SUM_MASSES, INV_SUM_MASSES, RIJ;

    real (kind=8) :: CHECK_MOLAR_MASS, DELTA_MOLAR_MASS, MAX_MOLAR_MASS;

    real (kind=8), dimension(1:1000) :: MOLAR_MASS;

    real (kind=8), dimension(1:3) :: RINO, RI, RAND_ANGLES, RG, DRIJ;

    real (kind=8), dimension(1:3,1:3) :: MATROT;

    real (kind=8), allocatable, dimension(:) :: MOLECULE_MOLAR_MASS;


    integer (kind=4) :: NATOM_NEW,     &
                        NBOND_NEW,     &
                        NANGLE_NEW,    &
                        NDIHEDRAL_NEW, &
                        NIMPROPER_NEW; 

    integer (kind=4) :: NTYPE_ATOM_NEW,     &
                        NTYPE_BOND_NEW,     &
                        NTYPE_ANGLE_NEW,    &
                        NTYPE_DIHEDRAL_NEW, &
                        NTYPE_IMPROPER_NEW;

    integer (kind=4) :: MAX_NTYPE_ATOM,     &
                        MAX_NTYPE_BOND,     &
                        MAX_NTYPE_ANGLE,    &
                        MAX_NTYPE_DIHEDRAL, &
                        MAX_NTYPE_IMPROPER;

!   ************************************************************************************************

    integer (kind=4) :: IOSEF1, IOSEF2, IOSEF3, IOSEF4, IOSEF5;

    character (len=100) :: CHOSEF1, CHOSEF2;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the routine #############################################################
                                                                                         !
    CHTITLE = 'Random  Deletion Molecule';                                               !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialization of arrays ###################################################################
                                                                                         !
    MOLAR_MASS(1:1000) = 0.0d0;                                                          !
                                                                                         !
    NCOUNT_DELETION(1:1000) = 0;                                                         !
                                                                                         !
    MOLECULEID_DELETION(1:1000,1:1000) = 0;                                              !
                                                                                         !
    MAX_NTYPE_ATOM     = 0;                                                              !
                                                                                         !
    MAX_NTYPE_BOND     = 0;                                                              !
                                                                                         !
    MAX_NTYPE_ANGLE    = 0;                                                              !
                                                                                         !
    MAX_NTYPE_DIHEDRAL = 0;                                                              !
                                                                                         !
    MAX_NTYPE_IMPROPER = 0;                                                              !
                                                                                         !
    MAX_MOLAR_MASS = 0.0d0;                                                              !
                                                                                         !
!   ### Read properties of library molecules #######################################################
                                                                                         !
    do i = 1, NFILE_LIBRARY_REMOVE;                                                      !
                                                                                         !
        write(icanal,'(a17,i4,a49)') '| LIBRARY FILE # ', i, ' '//REPEAT('/',47)//'|';   !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        IOSEF1 = 67 - LEN_TRIM(CHNAME_FILE_LIBRARY_REMOVE(i));                           !
                                                                                         !           
        write(icanal,'(a70)') '| '//TRIM(CHNAME_FILE_LIBRARY_REMOVE(i))// &              !
                              REPEAT(' ',IOSEF1)//'|';                                   !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        write(icanal,'(a23,i8,a39)') '| NATOM_LIBRARY      : ',     &                    !
                                     NATOM_LIBRARY(i),              &                    !
                                     REPEAT(' ',38)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        write(icanal,'(a23,i8,a39)') '| NTYPE_ATOM_LIBRARY : ',     &                    !
                                     NTYPE_ATOM_LIBRARY(i),         &                    !
                                     REPEAT(' ',38)//'|';                                !

        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
                                                                                         !
        write(icanal,'(a23,i8,a39)') '| NTYPE_BOND_LIBRARY : ',     &                    !  
                                     NTYPE_BOND_LIBRARY(i),         &                    !
                                     REPEAT(' ',38)//'|';                                !

        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
 
        write(icanal,'(a23,i8,a39)') '| NTYPE_ANGLE_LIBRARY : ',    &                    !
                                     NTYPE_ANGLE_LIBRARY(i),        &                    !
                                     REPEAT(' ',38)//'|';                                !

        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

        write(icanal,'(a23,i8,a39)') '| NTYPE_DIHEDRAL_LIBRARY : ', &                    !
                                     NTYPE_DIHEDRAL_LIBRARY(i),     &                    !
                                     REPEAT(' ',38)//'|';                                !

        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

        write(icanal,'(a23,i8,a39)') '| NTYPE_IMPROPER_LIBRARY : ', &                    !
                                     NTYPE_IMPROPER_LIBRARY(i),     &                    !
                                     REPEAT(' ',38)//'|';                                !

        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

        stop;
                                                                                         !
        if ( NTYPE_ATOM_LIBRARY(i) > MAX_NTYPE_ATOM ) MAX_NTYPE_ATOM  = &                !
                                                      NTYPE_ATOM_LIBRARY(i);             !

        if ( NTYPE_BOND_LIBRARY(i) > MAX_NTYPE_BOND ) MAX_NTYPE_BOND  = &                !
                                                      NTYPE_BOND_LIBRARY(i);             !

        if ( NTYPE_ANGLE_LIBRARY(i) > MAX_NTYPE_ANGLE ) MAX_NTYPE_ANGLE = &              !
                                                        NTYPE_ANGLE_LIBRARY(i);          !

        if ( NTYPE_DIHEDRAL_LIBRARY(i) > MAX_NTYPE_DIHEDRAL ) MAX_NTYPE_DIHEDRAL = &     !
                                                              NTYPE_DIHEDRAL_LIBRARY(i); !

        if ( NTYPE_IMPROPER_LIBRARY(i) > MAX_NTYPE_IMPROPER ) MAX_NTYPE_IMPROPER = &     !
                                                              NTYPE_IMPROPER_LIBRARY(i); !    
                                                                                         !
        do j = 1, NATOM_LIBRARY(i);                                                      ! LOOP OVER THE NUMBER OF ATOMS IN THE MOLECULE TYPE I
                                                                                         !
            do k = 1, NTYPE_ATOM_LIBRARY(i);                                             ! LOOP OVER THE NUMBER OF ATOM TYPES IN THE MOLECULE I
                                                                                         !
                if ( TRIM(CONFIG_NAT_LIBRARY(j,i)) == &                                  !
                     TRIM(ATOM_LABEL_LIBRARY(k,i)) ) then;                               !
                    MOLAR_MASS(i) = MOLAR_MASS(i) + ATOM_MASSES_LIBRARY(k,i);            !
                    EXIT;                                                                !
                end if                                                                   !
            end do                                                                       !
                                                                                         !
        end do                                                                           !

        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
                                                                                         !
        write(icanal,'(a23,f12.6,a35)') '| MOLAR_MASS [g/mol] : ', &                     !
                                        MOLAR_MASS(i),             &                     !
                                        REPEAT(' ',34)//'|';                             !
                                                                                         !
        if ( MOLAR_MASS(i) > MAX_MOLAR_MASS ) MAX_MOLAR_MASS = MOLAR_MASS(i);

        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a23,i8,a39)') '| MAX_NTYPE_ATOM     : ', MAX_NTYPE_ATOM, &            !
                                 REPEAT(' ',38)//'|';                                    !
    write(icanal,'(a23,i8,a39)') '| MAX_NTYPE_BOND     : ', MAX_NTYPE_BOND, &            !
                                 REPEAT(' ',38)//'|';                                    !
    write(icanal,'(a23,i8,a39)') '| MAX_NTYPE_ANGLE    : ', MAX_NTYPE_ANGLE, &           !
                                 REPEAT(' ',38)//'|';                                    !
    write(icanal,'(a23,i8,a39)') '| MAX_NTYPE_DIHEDRAL : ', MAX_NTYPE_DIHEDRAL, &        !
                                 REPEAT(' ',38)//'|';                                    !
    write(icanal,'(a23,i8,a39)') '| MAX_NTYPE_IMPROPER : ', MAX_NTYPE_IMPROPER, &        !
                                 REPEAT(' ',38)//'|';                                    !
                                                                                         !  
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a27,f12.6,a31)') '| MAX_MOLAR_MASS [g/mol] : ', &
                                    MAX_MOLAR_MASS,                &
                                    REPEAT(' ',30)//'|';

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate local arrays ######################################################################
                                                                                         !
    call ALLOCATE_PARAMETER_ARRAY_OLD(icanal,         &                                  !
                                      NTYPE_ATOM,     &                                  !
                                      NTYPE_BOND,     &                                  !
                                      NTYPE_ANGLE,    &                                  !
                                      NTYPE_DIHEDRAL, &                                  !
                                      NTYPE_IMPROPER);                                   !
                                                                                         !
    write(icanal,'(a70)') '| LOCAL ARRAYS WERE ALLOCATED'//REPEAT(' ',40)//'|';          !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### ALLOCATE ARRAYS TO CHECK MOLECULES TO REMOVE ###############################################
                                                                                         !
    allocate(LIST_ATOM_TYPE(1:NTYPE_ATOM));                                              !
                                                                                         !
    allocate(LIST_BOND_TYPE(1:NTYPE_BOND));                                              !

!   ### Initialisation of variables and arrays #####################################################
                                                                                         !
    CHECK_NTYPE_ATOM = 0;                                                                !
                                                                                         !
    CHECK_NTYPE_BOND = 0;                                                                ! 
                                                                                         !
    NLIST_ATOM_TYPE = 0;                                                                 !
                                                                                         !
    NLIST_BOND_TYPE = 0;                                                                 !
                                                                                         !
    LIST_ATOM_TYPE(1:NTYPE_ATOM) = 0;                                                    !
                                                                                         !
    LIST_BOND_TYPE(1:NTYPE_BOND) = 0;                                                    !

!   ### FIND THE MAXIMUM NUMBER OF MOLECULES IN THE SIMULATION BOX #################################
                                                                                         !
    MAX_MOLECULEID = MAXVAL( CONFIG_MOLECULEID(1:NATOM) );                               !
                                                                                         !
    write(icanal,'(a19,i8,a43)') '| MAX_MOLECULEID : ', &                                ! 
                                 MAX_MOLECULEID,        &                                !
                                 REPEAT(' ',42)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !

!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### COMPUTE THE MOLAR MASS OF EACH MOLECULES ###################################################
                                                                                         !
    allocate(MOLECULE_MOLAR_MASS(1:MAX_MOLECULEID));                                     !
                                                                                         !
    MOLECULE_MOLAR_MASS(1:MAX_MOLECULEID) = 0.0d0;                                       ! 
                                                                                         !
    do j = 1, MAX_MOLECULEID;                                                            !
                                                                                         !
        do i = 1, NATOM;                                                                 !
                                                                                         !
            if ( CONFIG_MOLECULEID(i) == j ) MOLECULE_MOLAR_MASS(j) = &                  !
                                             MOLECULE_MOLAR_MASS(j) + &                  !
                                             ATOM_MASSE(CONFIG_ATOM_TYPE(i));            !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !

!   ### COMPUTE THE NUMBER OF ATOM TYPE PER MOLECULE ###############################################    

    allocate(MOLECULE_NTYPE_ATOM(1:MAX_MOLECULEID));

    MOLECULE_NTYPE_ATOM(1:MAX_MOLECULEID) = 0;

    do j = 1, MAX_MOLECULEID;                                                            !
                                                                                         !
        if ( MOLECULE_MOLAR_MASS(j) > MAX_MOLAR_MASS ) CYCLE;                            !
                                                                                         !
        CHECK_NTYPE_ATOM = 0;                                                            ! SET TO 0 THE NUMBER OF ATOM TYPES TO CHECK
                                                                                         !
        LIST_ATOM_TYPE(1:NTYPE_ATOM) = 0;                                                ! SET TO 0 THE LIST OF ATOM TYPES
                                                                                         !
        do i = 1, NATOM;                                                                 ! LOOP OVER THE NUMBER OF ATOMS IN THE INPUT MOLECULAR CONFIGURATION
                                                                                         ! TO FIND OTHER ATOMS IN THE CURRENT MOLECULE
                                                                                         !
            if ( CONFIG_MOLECULEID(i) /= j ) CYCLE;                                      !
                                                                                         !
            if ( CHECK_NTYPE_ATOM == 0 ) then;                                           ! 
                                                                                         !
                CHECK_NTYPE_ATOM = CHECK_NTYPE_ATOM + 1;                                 ! UPDATE THE NUMBER OF ATOM TYPES IN THE MOLECULE RELATED TO ATOM I 
                                                                                         !
                LIST_ATOM_TYPE(CHECK_NTYPE_ATOM) = CONFIG_ATOM_TYPE(i);                  ! UPDATE THE LIST OF ATOM TYPES IN THE MOLECULE RELATED TO ATOM I
                                                                                         !
            else                                                                         !
                                                                                         !
                IFOUND = 0;                                                              !
                                                                                         !
                do k = 1, CHECK_NTYPE_ATOM;                                              ! LOOP OVER THE NUMBER OF ATOM TYPES TO CHECK
                                                                                         !
                    if ( LIST_ATOM_TYPE(k) == CONFIG_ATOM_TYPE(i) ) then;                !
                                                                                         !
                        IFOUND = 1;                                                      !
                                                                                         !
                        EXIT;                                                            !

                    end if                                                               !
                                                                                         !
                end do                                                                   !
                                                                                         !
                if ( IFOUND == 0 ) then;                                                 !
                                                                                         !
                    CHECK_NTYPE_ATOM = CHECK_NTYPE_ATOM + 1;                             !

                    if ( CHECK_NTYPE_ATOM > NTYPE_ATOM ) then;

                        write(icanal,*) 'CONFIG_MOLECULEID(i) : ', &
                                        CONFIG_MOLECULEID(i), &
                                        CONFIG_ATOM_TYPE(i);
                        write(icanal,*);
                        write(icanal,*) 'CHECK_NTYPE_ATOM | MAX_NTYPE_ATOM : ', &
                                        CHECK_NTYPE_ATOM, MAX_NTYPE_ATOM;
                        write(icanal,*);
                        write(icanal,*) 'THE PROGRAM IS GOING TO STOP !';
                        stop;

                    end if
                                                                                         !
                    LIST_ATOM_TYPE(CHECK_NTYPE_ATOM) = CONFIG_ATOM_TYPE(i);              !
                                                                                         !
                end if

            end if
                                                                                         !
        end do;                                                                          !

        MOLECULE_NTYPE_ATOM(j) = CHECK_NTYPE_ATOM;

!       write(icanal,*) ' MOLECULE_NTYPE_ATOM ', MOLECULE_NTYPE_ATOM(j), MAX_NTYPE_ATOM;

    end do
                                                                                         !
!   ### COMPUTE THE NUMBER OF BOND TYPE PER MOLECULE ###############################################
                                                                                         !
    allocate(MOLECULE_NTYPE_BOND(1:MAX_MOLECULEID));                                     !
                                                                                         !
    MOLECULE_NTYPE_BOND(1:MAX_MOLECULEID) = 0;                                           !
                                                                                         !
    do j = 1, MAX_MOLECULEID;                                                            !
                                                                                         !
        if ( MOLECULE_MOLAR_MASS(j) > MAX_MOLAR_MASS ) CYCLE;                            !                                                             
                                                                                         !
        CHECK_NTYPE_BOND = 0;                                                            !
                                                                                         !
        LIST_BOND_TYPE(1:NTYPE_BOND) = 0;                                                !
                                                                                         !
        do i = 1, NATOM;                                                                 !
                                                                                         !
            if ( CONFIG_MOLECULEID(i) /= j ) CYCLE;

            do m = 1, NBOND;                                                                 ! LOOP OVER THE NUMBER OF BONDS IN THE INPUT MOLECULAR CONFIGURATION
                                                                                         !
                do n = 1, 2;                                                                 !
                                                                                         !
                    if ( BOND_ATOMID(n,m) /= CONFIG_ATOMID(i) ) CYCLE;

                    IFOUND = 0;

                    if ( CHECK_NTYPE_BOND > 0 ) then;

                        IFOUND = 0;

                        do k = 1, CHECK_NTYPE_BOND;

                            if ( BOND_TYPE(m) /= LIST_BOND_TYPE(k) ) CYCLE;

                            IFOUND = 1;

                            EXIT;

                        end do

                    else

                        IFOUND = 0;

                    end if

                    if ( IFOUND == 0 ) then;

                        CHECK_NTYPE_BOND = CHECK_NTYPE_BOND + 1;

                        LIST_BOND_TYPE(CHECK_NTYPE_BOND) = BOND_TYPE(m);

!                       EXIT;

                    end if

                end do

            end do

        end do

        MOLECULE_NTYPE_BOND(j) = CHECK_NTYPE_BOND;

    end do
                                                                                         !
!   stop;                                                                                !
                                                                                         !
!   ### FIND AND REMOVE MOLECULES ##################################################################
                                                                                         !
    do i = 1, NATOM;                                                                     ! LOOP OVER THE NUMBER OF ATOMS IN THE INPUT CONFIGURATION
                                                                                         !
        if ( MOLECULE_MOLAR_MASS(CONFIG_MOLECULEID(i)) > MAX_MOLAR_MASS ) CYCLE;         ! IF THE MOLAR MASS OF THE MOLECULE TO WHICH THE ATOM BELONGS IS GREATER THAN
                                                                                         ! THE MAXIMUM MOLAR MASS OF MOLECULES THAT HAVE TO BE REMOVED, THEN CYCLE;
                                                                                         !
        if ( CONFIG_MOLECULEID(i) == 0 ) CYCLE;                                          ! IF THE MOLECULE ID IS EQUAL TO 0, THEN CYCLE
                                                                                         !
        CHECK_MOLAR_MASS = ATOM_MASSE(CONFIG_ATOM_TYPE(i));                              ! GET THE MOLAR MASS OF ATOM I BELONGING TO A GIVEN MOLECULE ID
                                                                                         !
        IFOUND = 0;                                                                      !
                                                                                         !
        do j = 1, NFILE_LIBRARY_REMOVE;                                                  ! LOOP OVER THE NUMBER OF FILES RELATED TO MOLECULES THAT HAS TO BE REMOVED
                                                                                         !
            if ( NCOUNT_DELETION(j) == 0 ) CYCLE;                                        ! IF NO MOLECULES J WERE REMOVED YET, THEN CYCLE 
                                                                                         !
            do k = 1, NCOUNT_DELETION(j);                                                ! LOOP OVER THE NUMBER OF ALREADY DELETED MOLECULES
                                                                                         !
                if ( CONFIG_MOLECULEID(i) /= MOLECULEID_DELETION(k,j) ) CYCLE;           ! IF THE ID OF MOLECULE I IS NOT IN THE LIST OF DELETED MOLECULES, THEN CYCLE
                                                                                         !
                IFOUND = 1;                                                              ! IF THE ID OF MOLECULE I IS IN THE LIST OF DELETED MOLECULES, THEN SET IFOUND TO 1
                                                                                         !
                EXIT;                                                                    ! EXIT THE CURRENT LOOP
                                                                                         !
            end do                                                                       !
                                                                                         !
            if ( IFOUND == 1 ) EXIT;                                                     ! IF THE ID OF MOLECULE I WAS FOUND IN THE LIST OF DELETED MOLECULES, THEN EXIT THE LOOP
                                                                                         !
        end do                                                                           !
                                                                                         !
        if ( IFOUND == 1 ) CYCLE;                                                        ! IF THE ID OF MOLECULE I WAS FOUND IN THE LIST, THEN CYCLE
                                                                                         ! ELSE, A NEW MOLECULE HAS BEEN FOUND
                                                                                         !
        do j = 1, NATOM;                                                                 ! LOOP OVER THE NUMBER OF ATOMS IN THE INPUT MOLECULAR CONFIGURATION
                                                                                         ! TO FIND OTHER ATOMS IN THE CURRENT MOLECULE
            if ( j == i ) CYCLE;                                                         ! IF ATOM J IS THE SAME AS ATOM I, THEN LOOP
                                                                                         !
            if ( CONFIG_MOLECULEID(j) /= CONFIG_MOLECULEID(i) ) CYCLE;                   ! IF THE ATOM J DOES NOT BELONG TO THE SAME MOLECULE AS ATOM I, THEN CYCLE
                                                                                         !
            CHECK_MOLAR_MASS = &                                                         ! UPDATE THE MOLAR MASS OF THE CHECKED MOLECULE
            CHECK_MOLAR_MASS + &                                                         !
            ATOM_MASSE(CONFIG_ATOM_TYPE(j));                                             !
                                                                                         !
        end do;                                                                          !

!       stop;

        do j = 1, NFILE_LIBRARY_REMOVE;                                                  ! LOOP OVER THE LIST OF MOLECULES TO REMOVE
                                                                                         !
            DELTA_MOLAR_MASS = ABS( CHECK_MOLAR_MASS - MOLAR_MASS(j) );                  !

            if ( ( DELTA_MOLAR_MASS < 0.1d0 )                   .AND. &                  ! IF THE CHECKED MOLECULE HAS THE SAME MOLAR MASS AS THE MOLECULE J FROM THE LIBRARY
                 ( MOLECULE_NTYPE_ATOM(CONFIG_MOLECULEID(i)) == NTYPE_ATOM_LIBRARY(j) )  .AND. &                  ! AND THE NUMBER OF ATOM TYPES IS THE SAME
                 ( MOLECULE_NTYPE_BOND(CONFIG_MOLECULEID(i)) == NTYPE_BOND_LIBRARY(j) ) ) then;                   ! AND THE NUMBER OF BOND TYPES IS THE SAME, THEN
                 
                                                                                         !
                if ( NCOUNT_DELETION(j) == NMOLECULE_DELETION(j) ) CYCLE;                !
                                                                                         !
                IFOUND = 0;

                do k = 1, NFILE_LIBRARY_REMOVE;

                    do m = 1, NCOUNT_DELETION(k);

                        if ( CONFIG_MOLECULEID(i) == MOLECULEID_DELETION(m,k) ) then;

                            IFOUND = 1;

                            EXIT;

                        end if

                    end do

                    if ( IFOUND == 1 ) EXIT;

                end do

                if ( IFOUND == 1 ) CYCLE;


                write(icanal,'(a70)') '|/// FOUND '//REPEAT('/',58)//'|';                !
                write(icanal,'(a2,i8,2f12.6,a36)') '| ',                 &               !
                                                   CONFIG_MOLECULEID(i), &               !
                                                   CHECK_MOLAR_MASS,     &               !
                                                   MOLAR_MASS(j),        &               !
                                                   REPEAT(' ',35)//'|';                  !
                                                                                         !
                write(icanal,'(a29,f12.6,a29)') '| CHECK_MOLAR_MASS [g/mol] : ', &       !
                                                CHECK_MOLAR_MASS,                &       !
                                                REPEAT(' ',28)//'|';                     !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                write(icanal,'(a21,2i8,a33)') '| CHECK_NTYPE_ATOM : ', &                 ! 
                                              MOLECULE_NTYPE_ATOM(CONFIG_MOLECULEID(i)),        &                 !
                                              NTYPE_ATOM_LIBRARY(j),   &                 !
                                              REPEAT(' ',32)//'|';                       !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                write(icanal,'(a21,2i8,a33)') '| CHECK_NTYPE_BOND : ', &                 ! 
                                              MOLECULE_NTYPE_BOND(CONFIG_MOLECULEID(i)),        &                 !
                                              NTYPE_BOND_LIBRARY(j),   &                 !
                                              REPEAT(' ',32)//'|';                       !  
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                NCOUNT_DELETION(j) = NCOUNT_DELETION(j) + 1;                             ! UPDATE THE NUMBER OF DELETED MOLECULES J 
                                                                                         !
                MOLECULEID_DELETION(NCOUNT_DELETION(j),j) = CONFIG_MOLECULEID(i);        !
                                                                                         !
                write(icanal,'(a20,i8,a42)') '| NCOUNT_DELETION : ', &                   !
                                              NCOUNT_DELETION(j),    &                   !
                                              REPEAT(' ',41)//'|';                       !
                                                                                         !
                                                                                         !
                do k = 1, NATOM;                                                         !
                                                                                         !
                    if ( CONFIG_MOLECULEID(k) /= CONFIG_MOLECULEID(i) ) CYCLE;           !
                                                                                         !
                    write(icanal,'(a2,3i8,a44)') '| ',                 &                 !
                                                 CONFIG_MOLECULEID(k), &                 !
                                                 CONFIG_ATOMID(k),     &                 !
                                                 CONFIG_ATOM_TYPE(k),  &                 !
                                                 REPEAT(' ',43)//'|';                    !
                                                                                         !
                                                                                         !
                end do                                                                   !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end do                                                                               ! 
                                                                                         !
!   ### UPDATE SYSTEM PROPERTIES ###################################################################
                                                                                         !
    NATOM_NEW     = NATOM;                                                               !
    NBOND_NEW     = NBOND;                                                               !
    NANGLE_NEW    = NANGLE;                                                              !
    NDIHEDRAL_NEW = NDIHEDRAL;                                                           !
    NIMPROPER_NEW = NIMPROPER;                                                           !
                                                                                         !
    NTYPE_ATOM_NEW     = NTYPE_ATOM;                                                     !
    NTYPE_BOND_NEW     = NTYPE_BOND;                                                     !
    NTYPE_ANGLE_NEW    = NTYPE_ANGLE;                                                    !
    NTYPE_DIHEDRAL_NEW = NTYPE_DIHEDRAL;                                                 !
    NTYPE_IMPROPER_NEW = NTYPE_IMPROPER;                                                 !
                                                                                         !
!   ### REMOVE SELECTED MOLECULES ##################################################################
                                                                                         !
    do i = 1, NFILE_LIBRARY_REMOVE;                                                      !
                                                                                         !
        IOSEF1 = 1 + LEN_TRIM(CHNAME_FILE_LIBRARY_REMOVE(i));                            !
                                                                                         !
        write(CHOSEF1,'(i8)') IOSEF1;                                                    !
                                                                                         !
        IOSEF2 = 53 - IOSEF1;                                                            !
                                                                                         !
        write(CHOSEF2,'(i8)') IOSEF2+1;                                                  !
                                                                                         !
        write(icanal,'(a'//TRIM(CHOSEF1)//     &                                         !
                     ',2i8,a'//TRIM(CHOSEF2)// &                                         !
                     ')') '|'//TRIM(CHNAME_FILE_LIBRARY_REMOVE(i)), &                    !
                          NCOUNT_DELETION(i),                       &                    !
                          NMOLECULE_DELETION(i),                    &                    !
                          REPEAT(' ',IOSEF2)//'|';                                       !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !

        if ( NMOLECULE_DELETION(i) == NCOUNT_DELETION(i) ) then;

            do j = 1, NCOUNT_DELETION(i);
                write(icanal,'(a12,i4,a54)')  '| DELETING #',           &
                                              MOLECULEID_DELETION(j,i), &
                                              REPEAT(' ',53)//'|';

                do k = 1, NATOM;
                    if ( CONFIG_MOLECULEID(k) /= MOLECULEID_DELETION(j,i) ) CYCLE; 

                    do m = 1, NBOND;
                        do n = 1, 2;
                            if ( BOND_ATOMID(n,m) /= CONFIG_ATOMID(k) ) CYCLE;

                            BOND_ATOMID(n,m) = 0;

                            if ( BOND_TYPE(m) /= 0 ) BOND_COEFFS(1:4,BOND_TYPE(m)) = 0.0d0;

                            BOND_TYPE(m) = 0;

                        end do
                    end do

                    do m = 1, NANGLE;
                        do n = 1, 3;
                            if ( ANGLE_ATOMID(n,m) /= CONFIG_ATOMID(k) ) CYCLE;

                            ANGLE_ATOMID(n,m) = 0;

                            if ( ANGLE_TYPE(m) /= 0 ) then;
                                ANGLE_COEFFS(1:4,ANGLE_TYPE(m)) = -99999.99d0;

                                BONDBOND_COEFFS(1:3,ANGLE_TYPE(m)) = -99999.99d0;

                                BONDANGLE_COEFFS(1:4,ANGLE_TYPE(m)) = -99999.99d0;
                            end if

                            ANGLE_TYPE(m) = 0;
!                           end if
                        end do
                    end do

                    do m = 1, NDIHEDRAL;
                        do n = 1, 4;
                            if ( DIHEDRAL_ATOMID(n,m) /= CONFIG_ATOMID(k) ) CYCLE;

                            DIHEDRAL_ATOMID(n,m) = 0;

                            if ( DIHEDRAL_TYPE(m) /= 0 ) then;

                                DIHEDRAL_COEFFS(1:6,DIHEDRAL_TYPE(m)) = -99999.99d0; 

                                MIDDLEBONDTORSION_COEFFS(1:4,DIHEDRAL_TYPE(m)) = -99999.99d0;

                                ENDBONDTORSION_COEFFS(1:8,DIHEDRAL_TYPE(m)) = -99999.99d0;

                                ANGLETORSION_COEFFS(1:8,DIHEDRAL_TYPE(m)) = -99999.99d0;

                                ANGLEANGLETORSION_COEFFS(1:3,DIHEDRAL_TYPE(m)) = -99999.99d0;

                                BONDBOND13_COEFFS(1:3,DIHEDRAL_TYPE(m)) = -99999.99d0;
 
                            end if

                            DIHEDRAL_TYPE(m) = 0;
!                           end if
                        end do
                    end do

                    do m = 1, NIMPROPER;
                        do n = 1, 4;
                            if ( IMPROPER_ATOMID(n,m) /= CONFIG_ATOMID(k) ) CYCLE;

                            IMPROPER_ATOMID(n,m) = 0;

                            if ( IMPROPER_TYPE(m) /= 0 ) then;
                            
                                IMPROPER_COEFFS(1:2,IMPROPER_TYPE(m)) = -99999.99d0;

                                ANGLEANGLE_COEFFS(1:6,IMPROPER_TYPE(m)) = -99999.99d0;

                            end if

                            IMPROPER_TYPE(m) = 0;
!                           end if
                        end do
                    end do

                    ATOM_LABEL(CONFIG_ATOM_TYPE(k)) = 'XXX';

                    ATOM_MASSE(CONFIG_ATOM_TYPE(k)) = 0.0d0;

                    POTENTIAL_CLASS2(1:2,CONFIG_ATOM_TYPE(k)) = 0.0d0;

                    CONFIG_MOLECULEID(k) = 0;

                    CONFIG_ATOMID(k) = 0;

                    CONFIG_ATOM_TYPE(k) = 0;

                    CONFIG_RI(1:3,k) = 0.0d0;

                    CONFIG_VI(1:3,k) = 0.0d0;

                    CONFIG_QI(k) = 0.0d0;

                    CONFIG_NAT(k) = 'XX';

                end do

                NATOM_NEW     = NATOM_NEW     - NATOM_LIBRARY(i);
                NBOND_NEW     = NBOND_NEW     - NBOND_LIBRARY(i);
                NANGLE_NEW    = NANGLE_NEW    - NANGLE_LIBRARY(i);
                NDIHEDRAL_NEW = NDIHEDRAL_NEW - NDIHEDRAL_LIBRARY(i);
                NIMPROPER_NEW = NIMPROPER_NEW - NIMPROPER_LIBRARY(i);

            end do
                                                                                         !
            NTYPE_ATOM_NEW     = NTYPE_ATOM_NEW     - NTYPE_ATOM_LIBRARY(i);             !
            NTYPE_BOND_NEW     = NTYPE_BOND_NEW     - NTYPE_BOND_LIBRARY(i);             !
            NTYPE_ANGLE_NEW    = NTYPE_ANGLE_NEW    - NTYPE_ANGLE_LIBRARY(i);            !
            NTYPE_DIHEDRAL_NEW = NTYPE_DIHEDRAL_NEW - NTYPE_DIHEDRAL_LIBRARY(i);         !
            NTYPE_IMPROPER_NEW = NTYPE_IMPROPER_NEW - NTYPE_IMPROPER_LIBRARY(i);         !
                                                                                         !
        else                                                                             !
                                                                                         !
            write(icanal,*) 'NUMBER OF EXPECTED DELETIONS DIFFERS FROM THE NUMBER OF DELETED MOLECULES';
            write(icanal,*);
            write(icanal,*) 'EXPECTED : ', NMOLECULE_DELETION(i)
            write(icanal,*) 'DELETED  : ', NCOUNT_DELETION(i)
            write(icanal,*);
            write(icanal,*) 'THE PROGRAM IS GOING TO STOP!';                             !

            stop;                                                                        !
                                                                                         !
        end if                                                                           !
                                                                                         ! 
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
    write(icanal,'(a19,i8,a43)') '| NATOM_NEW      : ', NATOM_NEW,     &                 !
                                 REPEAT(' ',42)//'|';                                    !
    write(icanal,'(a19,i8,a43)') '| NBOND_NEW      : ', NBOND_NEW,     &                 !
                                 REPEAT(' ',42)//'|';                                    !
    write(icanal,'(a19,i8,a43)') '| NANGLE_NEW     : ', NANGLE_NEW,    &                 !
                                 REPEAT(' ',42)//'|';                                    !
    write(icanal,'(a19,i8,a43)') '| NDIHEDRAL_NEW  : ', NDIHEDRAL_NEW, &                 !
                                 REPEAT(' ',42)//'|';                                    !
    write(icanal,'(a19,i8,a43)') '| NIMPROPER_NEW  : ', NIMPROPER_NEW, &                 !
                                 REPEAT(' ',42)//'|';                                    !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
    write(icanal,'(a23,i8,a39)') '| NTYPE_ATOM_NEW     : ', NTYPE_ATOM_NEW,     &        !
                                 REPEAT(' ',38)//'|';                                    ! 
    write(icanal,'(a23,i8,a39)') '| NTYPE_BOND_NEW     : ', NTYPE_BOND_NEW,     &        !
                                 REPEAT(' ',38)//'|';                                    !
    write(icanal,'(a23,i8,a39)') '| NTYPE_ANGLE_NEW    : ', NTYPE_ANGLE_NEW,    &        !
                                 REPEAT(' ',38)//'|';                                    !
    write(icanal,'(a23,i8,a39)') '| NTYPE_DIHEDRAL_NEW : ', NTYPE_DIHEDRAL_NEW, &        !
                                 REPEAT(' ',38)//'|';                                    !
    write(icanal,'(a23,i8,a39)') '| NTYPE_IMPROPER_NEW : ', NTYPE_IMPROPER_NEW, &        !
                                 REPEAT(' ',38)//'|';                                    !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### UPDATE POTENTIAL PARAMETERS ################################################################
                                                                                         !
    ATOM_LABEL_OLD(1:NTYPE_ATOM) = ATOM_LABEL(1:NTYPE_ATOM);                             !
                                                                                         !
    ATOM_MASSE_OLD(1:NTYPE_ATOM) = ATOM_MASSE(1:NTYPE_ATOM);                             !
                                                                                         !
    POTENTIAL_CLASS2_OLD(1:2,1:NTYPE_ATOM) = POTENTIAL_CLASS2(1:2,1:NTYPE_ATOM);         !
                                                                                         !
    write(icanal,'(a70)') '| LABELS AND MOLAR MASSES : '//REPEAT(' ',41)//'|';           !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !
    do i = 1, NTYPE_ATOM;                                                                !
                                                                                         !
        if ( ATOM_MASSE_OLD(i) == 0.0d0 ) CYCLE;                                         !
                                                                                         !
        IOSEF1 = IOSEF1 + 1;                                                             !
                                                                                         !
        ATOM_MASSE(IOSEF1) = ATOM_MASSE_OLD(i);                                          !
                                                                                         !
        ATOM_LABEL(IOSEF1) = ATOM_LABEL_OLD(i);                                          !
                                                                                         !
        write(icanal,'(a2,i8,2x,a3,2x,f12.6,a41)') '| ', &                               !
                                     IOSEF1,             &                               !
                                     ATOM_LABEL(IOSEF1), &                               !
                                     ATOM_MASSE(IOSEF1), &                               !
                                     REPEAT(' ',40)//'|';                                !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
    write(icanal,'(a2,3i8,a44)') '| ',           &                                       !
                                 IOSEF1,         &                                       !
                                 NTYPE_ATOM_NEW, &                                       !
                                 NTYPE_ATOM,     &                                       !
                                 REPEAT(' ',43)//'|';                                    !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
    write(icanal,'(a70)') '| POTENTIAL_CLASS2 : '//REPEAT(' ',48)//'|';                  !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !
    do i = 1, NTYPE_ATOM;                                                                !
                                                                                         !
        if ( ( POTENTIAL_CLASS2_OLD(1,i) == 0.0d0 ) .AND. &                              !
             ( POTENTIAL_CLASS2_OLD(2,i) == 0.0d0 ) ) CYCLE;                             !
                                                                                         !
        IOSEF1 = IOSEF1 + 1;                                                             !
                                                                                         !
        POTENTIAL_CLASS2(1:2,IOSEF1) = POTENTIAL_CLASS2_OLD(1:2,i);                      !
                                                                                         !
        write(icanal,'(a2,i4,2f12.6,a40)') '| ',                         &               !
                                           IOSEF1,                       &               !
                                           POTENTIAL_CLASS2(1:2,IOSEF1), &               !
                                           REPEAT(' ',39)//'|';                          !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      ! 
    write(icanal,'(a2,3i8,a44)') '| ',           &                                       !
                                 IOSEF1,         &                                       !
                                 NTYPE_ATOM_NEW, &                                       !
                                 NTYPE_ATOM,     &                                       !
                                 REPEAT(' ',43)//'|';                                    !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    BOND_COEFFS_OLD(1:4,1:NTYPE_BOND) = BOND_COEFFS(1:4,1:NTYPE_BOND);                   !
                                                                                         !
    write(icanal,'(a70)') '| BOND_COEFFS : '//REPEAT(' ',53)//'|';                       !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !
    do i = 1, NTYPE_BOND;                                                                !
                                                                                         !
        if ( ( BOND_COEFFS_OLD(1,i) == 0.0d0 ) .AND. &                                   !
             ( BOND_COEFFS_OLD(2,i) == 0.0d0 ) .AND. &                                   !
             ( BOND_COEFFS_OLD(3,i) == 0.0d0 ) .AND. &                                   !
             ( BOND_COEFFS_OLD(4,i) == 0.0d0 ) ) CYCLE;                                  !
                                                                                         !
        IOSEF1 = IOSEF1 + 1;                                                             !
                                                                                         !
        BOND_COEFFS(1:4,IOSEF1) = BOND_COEFFS_OLD(1:4,i);                                !
                                                                                         !
        write(icanal,'(a2,i4,4f15.6,a4)') '| ',                    &                     !
                                          IOSEF1,                  &                     !
                                          BOND_COEFFS(1:4,IOSEF1), &                     !
                                          REPEAT(' ',3)//'|';                            !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
    write(icanal,'(a2,3i8,a44)') '| ',           &                                       !
                                 IOSEF1,         &                                       !
                                 NTYPE_BOND_NEW, &                                       !
                                 NTYPE_BOND,     &                                       !
                                 REPEAT(' ',43)//'|';                                    !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    ANGLE_COEFFS_OLD(1:4,1:NTYPE_ANGLE) = ANGLE_COEFFS(1:4,1:NTYPE_ANGLE);               !
                                                                                         !
    BONDBOND_COEFFS_OLD(1:3,1:NTYPE_ANGLE) = BONDBOND_COEFFS(1:3,1:NTYPE_ANGLE);         ! 
                                                                                         !
    BONDANGLE_COEFFS_OLD(1:4,1:NTYPE_ANGLE) = BONDANGLE_COEFFS(1:4,1:NTYPE_ANGLE);       !
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !
    write(icanal,'(a70)') '| ANGLE_COEFFS : '//REPEAT(' ',52)//'|';                      !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    do i = 1, NTYPE_ANGLE;                                                               !
                                                                                         !
        if ( ( ANGLE_COEFFS_OLD(1,i) == -99999.99d0 ) .AND. &                            !
             ( ANGLE_COEFFS_OLD(2,i) == -99999.99d0 ) .AND. &                            !
             ( ANGLE_COEFFS_OLD(3,i) == -99999.99d0 ) .AND. &                            !
             ( ANGLE_COEFFS_OLD(4,i) == -99999.99d0 ) ) CYCLE;                           !
                                                                                         !
        IOSEF1 = IOSEF1 + 1;                                                             !
                                                                                         !
        ANGLE_COEFFS(1:4,IOSEF1) = ANGLE_COEFFS_OLD(1:4,i);                              !
                                                                                         !
        write(icanal,'(a2,i4,4f15.6,a4)') '| ',                     &                    !
                                          IOSEF1,                   &                    !
                                          ANGLE_COEFFS(1:4,IOSEF1), &                    !
                                          REPEAT(' ',3)//'|';                            !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      ! 
    write(icanal,'(a2,3i8,a44)') '| ',            &                                      !
                                 IOSEF1,          &                                      !
                                 NTYPE_ANGLE_NEW, &                                      !
                                 NTYPE_ANGLE,     &                                      !
                                 REPEAT(' ',43)//'|';                                    !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !
    write(icanal,'(a70)') '| BONDBOND_COEFFS : '//REPEAT(' ',49)//'|';                   !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    do i = 1, NTYPE_ANGLE;                                                               !
                                                                                         !
        if ( ( BONDBOND_COEFFS_OLD(1,i) == -99999.99d0 ) .AND. &                         !
             ( BONDBOND_COEFFS_OLD(2,i) == -99999.99d0 ) .AND. &                         !
             ( BONDBOND_COEFFS_OLD(3,i) == -99999.99d0 ) ) CYCLE;                        !
                                                                                         !
        IOSEF1 = IOSEF1 + 1;                                                             !
                                                                                         !
        BONDBOND_COEFFS(1:3,IOSEF1) = BONDBOND_COEFFS_OLD(1:3,i);                        !
                                                                                         !
        write(icanal,'(a2,i4,3f15.6,a19)') '| ',                        &                !
                                           IOSEF1,                      &                !
                                           BONDBOND_COEFFS(1:3,IOSEF1), &                !
                                           REPEAT(' ',18)//'|';                          !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
    write(icanal,'(a2,3i8,a44)') '| ',            &                                      !
                                 IOSEF1,          &                                      !
                                 NTYPE_ANGLE_NEW, &                                      !
                                 NTYPE_ANGLE,     &                                      !
                                 REPEAT(' ',43)//'|';                                    !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '| BONDANGLE_COEFFS : '//REPEAT(' ',48)//'|';                  !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !
    do i = 1, NTYPE_ANGLE;                                                               !
                                                                                         !
        if ( ( BONDANGLE_COEFFS_OLD(1,i) == -99999.99d0 ) .AND. &                        !
             ( BONDANGLE_COEFFS_OLD(2,i) == -99999.99d0 ) .AND. &                        !
             ( BONDANGLE_COEFFS_OLD(3,i) == -99999.99d0 ) .AND. &                        !
             ( BONDANGLE_COEFFS_OLD(4,i) == -99999.99d0 ) ) CYCLE;                       !
                                                                                         !
        IOSEF1 = IOSEF1 + 1;                                                             ! 
                                                                                         !
        BONDANGLE_COEFFS(1:4,IOSEF1) = BONDANGLE_COEFFS_OLD(1:4,i);                      !
                                                                                         !
        write(icanal,'(a2,i4,4f15.6,a4)') '| ',                         &                !
                                          IOSEF1,                       &                !
                                          BONDANGLE_COEFFS(1:4,IOSEF1), &                !
                                          REPEAT(' ',3)//'|';                            !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
    write(icanal,'(a2,3i8,a44)') '| ',            &                                      !
                                 IOSEF1,          &                                      !
                                 NTYPE_ANGLE_NEW, &                                      !
                                 NTYPE_ANGLE,     &                                      !
                                 REPEAT(' ',43)//'|';                                    !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    DIHEDRAL_COEFFS_OLD(1:6,1:NTYPE_DIHEDRAL) = DIHEDRAL_COEFFS(1:6,1:NTYPE_DIHEDRAL);   !
                                                                                         !
    MIDDLEBONDTORSION_COEFFS_OLD(1:4,1:NTYPE_DIHEDRAL) = &                               !
    MIDDLEBONDTORSION_COEFFS(1:4,1:NTYPE_DIHEDRAL);                                      !
                                                                                         !
    ENDBONDTORSION_COEFFS_OLD(1:8,1:NTYPE_DIHEDRAL) = &                                  !
    ENDBONDTORSION_COEFFS(1:8,1:NTYPE_DIHEDRAL);                                         !
                                                                                         !
    ANGLETORSION_COEFFS_OLD(1:8,1:NTYPE_DIHEDRAL) = &                                    !
    ANGLETORSION_COEFFS(1:8,1:NTYPE_DIHEDRAL);                                           !
                                                                                         !
    ANGLEANGLETORSION_COEFFS_OLD(1:3,1:NTYPE_DIHEDRAL) = &                               !
    ANGLEANGLETORSION_COEFFS(1:3,1:NTYPE_DIHEDRAL);                                      !
                                                                                         !
    BONDBOND13_COEFFS_OLD(1:3,1:NTYPE_DIHEDRAL) = &                                      !
    BONDBOND13_COEFFS(1:3,1:NTYPE_DIHEDRAL);                                             !
                                                                                         !
    write(icanal,'(a70)') '| DIHEDRAL_COEFFS :'//REPEAT(' ',50)//'|';                    !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !
    do i = 1, NTYPE_DIHEDRAL;                                                            !
                                                                                         !
        if ( ( DIHEDRAL_COEFFS_OLD(1,i) == -99999.99d0 ) .AND. &                         !
             ( DIHEDRAL_COEFFS_OLD(2,i) == -99999.99d0 ) .AND. &                         !
             ( DIHEDRAL_COEFFS_OLD(3,i) == -99999.99d0 ) .AND. &                         !
             ( DIHEDRAL_COEFFS_OLD(4,i) == -99999.99d0 ) .AND. &                         !
             ( DIHEDRAL_COEFFS_OLD(5,i) == -99999.99d0 ) .AND. &                         !
             ( DIHEDRAL_COEFFS_OLD(6,i) == -99999.99d0 ) ) CYCLE;                        !
                                                                                         !
        IOSEF1 = IOSEF1 + 1;                                                             !
                                                                                         !
        DIHEDRAL_COEFFS(1:6,IOSEF1) = DIHEDRAL_COEFFS_OLD(1:6,i);                        !
                                                                                         !
        write(icanal,'(a2,i4,4f15.6,a4)') '| ',                        &                 !
                                          IOSEF1,                      &                 ! 
                                          DIHEDRAL_COEFFS(1:4,IOSEF1), &                 !
                                          REPEAT(' ',3)//'|';                            !
                                                                                         !
        write(icanal,'(a6,2f15.6,a34)') '|'//REPEAT(' ',5),          &                   !
                                        DIHEDRAL_COEFFS(5:6,IOSEF1), &                   !
                                        REPEAT(' ',33)//'|';                             !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
    write(icanal,'(a2,3i8,a44)') '| ',               &                                   !
                                 IOSEF1,             &                                   !
                                 NTYPE_DIHEDRAL_NEW, &                                   !
                                 NTYPE_DIHEDRAL,     &                                   !
                                 REPEAT(' ',43)//'|';                                    !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '| MIDDLEBONDTORSION_COEFFS : '//REPEAT(' ',40)//'|';          !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !
    do i = 1, NTYPE_DIHEDRAL;                                                            !
                                                                                         !
        if ( ( MIDDLEBONDTORSION_COEFFS_OLD(1,i) == -99999.99d0 ) .AND. &                !
             ( MIDDLEBONDTORSION_COEFFS_OLD(2,i) == -99999.99d0 ) .AND. &                !
             ( MIDDLEBONDTORSION_COEFFS_OLD(3,i) == -99999.99d0 ) .AND. &                !
             ( MIDDLEBONDTORSION_COEFFS_OLD(4,i) == -99999.99d0 ) ) CYCLE;               !
                                                                                         !
        IOSEF1 = IOSEF1 + 1;                                                             ! 
                                                                                         !
        MIDDLEBONDTORSION_COEFFS(1:4,IOSEF1) = MIDDLEBONDTORSION_COEFFS_OLD(1:4,i);      !
                                                                                         !
        write(icanal,'(a2,i4,4f15.6,a4)') '| ',                                 &        !
                                          IOSEF1,                               &        !
                                          MIDDLEBONDTORSION_COEFFS(1:4,IOSEF1), &        !
                                          REPEAT(' ',3)//'|';                            !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
    write(icanal,'(a2,3i8,a44)') '| ',               &                                   !
                                 IOSEF1,             &                                   !
                                 NTYPE_DIHEDRAL_NEW, &                                   !
                                 NTYPE_DIHEDRAL,     &                                   !
                                 REPEAT(' ',43)//'|';                                    !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '| ENDBONDTORSION_COEFFS : '//REPEAT(' ',43)//'|';             !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !
    do i = 1, NTYPE_DIHEDRAL;                                                            !
                                                                                         !
        if ( ( ENDBONDTORSION_COEFFS_OLD(1,i) == -99999.99d0 ) .AND. &                   !
             ( ENDBONDTORSION_COEFFS_OLD(2,i) == -99999.99d0 ) .AND. &                   !
             ( ENDBONDTORSION_COEFFS_OLD(3,i) == -99999.99d0 ) .AND. &                   !
             ( ENDBONDTORSION_COEFFS_OLD(4,i) == -99999.99d0 ) .AND. &                   !
             ( ENDBONDTORSION_COEFFS_OLD(5,i) == -99999.99d0 ) .AND. &                   !
             ( ENDBONDTORSION_COEFFS_OLD(6,i) == -99999.99d0 ) .AND. &                   !
             ( ENDBONDTORSION_COEFFS_OLD(7,i) == -99999.99d0 ) .AND. &                   !
             ( ENDBONDTORSION_COEFFS_OLD(8,i) == -99999.99d0 ) ) CYCLE;                  !
                                                                                         !
        IOSEF1 = IOSEF1 + 1;                                                             !
                                                                                         !
        ENDBONDTORSION_COEFFS(1:8,IOSEF1) = ENDBONDTORSION_COEFFS_OLD(1:8,i);            !
                                                                                         !
                                                                                         !
        write(icanal,'(a2,i4,4f15.6,a4)') '| ',                              &           !
                                          IOSEF1,                            &           !
                                          ENDBONDTORSION_COEFFS(1:4,IOSEF1), &           !
                                          REPEAT(' ',3)//'|';                            !
                                                                                         !
        write(icanal,'(a6,4f15.6,a4)') '|'//REPEAT(' ',5),                &              !
                                       ENDBONDTORSION_COEFFS(5:8,IOSEF1), &              !
                                       REPEAT(' ',3)//'|';                               !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
    write(icanal,'(a2,3i8,a44)') '| ',               &                                   !
                                 IOSEF1,             &                                   !
                                 NTYPE_DIHEDRAL_NEW, &                                   !
                                 NTYPE_DIHEDRAL,     &                                   !
                                 REPEAT(' ',43)//'|';                                    !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '| ANGLETORSION_COEFFS : '//REPEAT(' ',45)//'|';               !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !
    do i = 1, NTYPE_DIHEDRAL;                                                            !
                                                                                         !
        if ( ( ANGLETORSION_COEFFS_OLD(1,i) == -99999.99d0 ) .AND. &                     !
             ( ANGLETORSION_COEFFS_OLD(2,i) == -99999.99d0 ) .AND. &                     !
             ( ANGLETORSION_COEFFS_OLD(3,i) == -99999.99d0 ) .AND. &                     !
             ( ANGLETORSION_COEFFS_OLD(4,i) == -99999.99d0 ) .AND. &                     !
             ( ANGLETORSION_COEFFS_OLD(5,i) == -99999.99d0 ) .AND. &                     !
             ( ANGLETORSION_COEFFS_OLD(6,i) == -99999.99d0 ) .AND. &                     !
             ( ANGLETORSION_COEFFS_OLD(7,i) == -99999.99d0 ) .AND. &                     !
             ( ANGLETORSION_COEFFS_OLD(8,i) == -99999.99d0 ) ) CYCLE;                    !
                                                                                         !
        IOSEF1 = IOSEF1 + 1;                                                             !
                                                                                         !
        ANGLETORSION_COEFFS(1:8,IOSEF1) = ANGLETORSION_COEFFS_OLD(1:8,i);                !
                                                                                         !
        write(icanal,'(a2,i4,4f15.6,a4)') '| ',                            &             !
                                          IOSEF1,                          &             !
                                          ANGLETORSION_COEFFS(1:4,IOSEF1), &             !
                                          REPEAT(' ',3)//'|';                            !
                                                                                         !
        write(icanal,'(a6,4f15.6,a4)') '|'//REPEAT(' ',5),              &                !
                                       ANGLETORSION_COEFFS(5:8,IOSEF1), &                !
                                       REPEAT(' ',3)//'|';                               !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
    write(icanal,'(a2,3i8,a44)') '| ',               &                                   !
                                 IOSEF1,             &                                   !
                                 NTYPE_DIHEDRAL_NEW, &                                   !
                                 NTYPE_DIHEDRAL,     &                                   !
                                 REPEAT(' ',43)//'|';                                    !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '| ANGLEANGLETORSION_COEFFS : '//REPEAT(' ',40)//'|';          !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !
    do i = 1, NTYPE_DIHEDRAL;                                                            !
                                                                                         !
        if ( ( ANGLEANGLETORSION_COEFFS_OLD(1,i) == -99999.99d0 ) .AND. &                !
             ( ANGLEANGLETORSION_COEFFS_OLD(2,i) == -99999.99d0 ) .AND. &                !
             ( ANGLEANGLETORSION_COEFFS_OLD(3,i) == -99999.99d0 ) ) CYCLE;               !
                                                                                         !
        IOSEF1 = IOSEF1 + 1;                                                             !
                                                                                         ! 
        ANGLEANGLETORSION_COEFFS(1:3,IOSEF1) = ANGLEANGLETORSION_COEFFS_OLD(1:3,i);      !
                                                                                         !
        write(icanal,'(a2,i4,3f15.6,a19)') '| ',                                 &       !
                                           IOSEF1,                               &       !
                                           ANGLEANGLETORSION_COEFFS(1:3,IOSEF1), &       !
                                           REPEAT(' ',18)//'|';                          !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
    write(icanal,'(a2,3i8,a44)') '| ',               &                                   !
                                 IOSEF1,             &                                   !
                                 NTYPE_DIHEDRAL_NEW, &                                   !
                                 NTYPE_DIHEDRAL,     &                                   !
                                 REPEAT(' ',43)//'|';                                    !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '| BONDBOND13_COEFFS : '//REPEAT(' ',47)//'|';                 !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !
    do i = 1, NTYPE_DIHEDRAL;                                                            !
                                                                                         !
        if ( ( BONDBOND13_COEFFS_OLD(1,i) == -99999.99d0 ) .AND. &                       !
             ( BONDBOND13_COEFFS_OLD(2,i) == -99999.99d0 ) .AND. &                       !
             ( BONDBOND13_COEFFS_OLD(3,i) == -99999.99d0 ) ) CYCLE;                      !
                                                                                         !
        IOSEF1 = IOSEF1 + 1;                                                             !
                                                                                         !
        BONDBOND13_COEFFS(1:3,IOSEF1) = BONDBOND13_COEFFS_OLD(1:3,i);                    !
                                                                                         !
        write(icanal,'(a2,i4,3f15.6,a19)') '| ',                          &              !
                                           IOSEF1,                        &              !
                                           BONDBOND13_COEFFS(1:3,IOSEF1), &              !
                                           REPEAT(' ',18)//'|';                          !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
    write(icanal,'(a2,3i8,a44)') '| ',               &                                   !
                                 IOSEF1,             &                                   !
                                 NTYPE_DIHEDRAL_NEW, &                                   !
                                 NTYPE_DIHEDRAL,     &                                   !
                                 REPEAT(' ',43)//'|';                                    !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      ! 
                                                                                         !
    IMPROPER_COEFFS_OLD(1:2,1:NTYPE_IMPROPER) = IMPROPER_COEFFS(1:2,1:NTYPE_IMPROPER);   !
                                                                                         ! 
    ANGLEANGLE_COEFFS_OLD(1:6,1:NTYPE_IMPROPER) = &                                      !
    ANGLEANGLE_COEFFS(1:6,1:NTYPE_IMPROPER);                                             !
                                                                                         !
    write(icanal,'(a70)') '| IMPROPER_COEFFS : '//REPEAT(' ',49)//'|';                   !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !  
    do i = 1, NTYPE_IMPROPER;                                                            !
                                                                                         !
        if ( ( IMPROPER_COEFFS_OLD(1,i) == -99999.99d0 ) .AND. &                         !
             ( IMPROPER_COEFFS_OLD(2,i) == -99999.99d0 ) ) CYCLE;                        !
                                                                                         !
        IOSEF1 = IOSEF1 + 1;                                                             !
                                                                                         ! 
        IMPROPER_COEFFS(1:2,IOSEF1) = IMPROPER_COEFFS_OLD(1:2,i);                        !
                                                                                         !
        write(icanal,'(a2,i4,2f15.6,a34)') '| ',                        &                !
                                           IOSEF1,                      &                !
                                           IMPROPER_COEFFS(1:2,IOSEF1), &                !
                                           REPEAT(' ',33)//'|';                          !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
    write(icanal,'(a2,3i8,a44)') '| ',               &                                   !
                                 IOSEF1,             &                                   ! 
                                 NTYPE_IMPROPER_NEW, &                                   !
                                 NTYPE_IMPROPER,     &                                   !
                                 REPEAT(' ',43)//'|';                                    !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '| ANGLEANGLE_COEFFS : '//REPEAT(' ',47)//'|';                 !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !
    do i = 1, NTYPE_IMPROPER;                                                            !
                                                                                         !
        if ( ( ANGLEANGLE_COEFFS_OLD(1,i) == -99999.99d0 ) .AND. &                       !
             ( ANGLEANGLE_COEFFS_OLD(2,i) == -99999.99d0 ) .AND. &                       !
             ( ANGLEANGLE_COEFFS_OLD(3,i) == -99999.99d0 ) .AND. &                       !
             ( ANGLEANGLE_COEFFS_OLD(4,i) == -99999.99d0 ) .AND. &                       !
             ( ANGLEANGLE_COEFFS_OLD(5,i) == -99999.99d0 ) .AND. &                       !
             ( ANGLEANGLE_COEFFS_OLD(6,i) == -99999.99d0 ) ) CYCLE;                      !
                                                                                         !
        IOSEF1 = IOSEF1 + 1;                                                             !
                                                                                         !
        ANGLEANGLE_COEFFS(1:6,IOSEF1) = ANGLEANGLE_COEFFS_OLD(1:6,i);                    !
                                                                                         !
        write(icanal,'(a2,i4,4f15.6,a4)') '| ',                          &               !
                                          IOSEF1,                        &               !
                                          ANGLEANGLE_COEFFS(1:4,IOSEF1), &               !
                                          REPEAT(' ',3)//'|';                            !
                                                                                         !
        write(icanal,'(a6,2f15.6,a34)') '|'//REPEAT(' ',5),           &                  !
                                       ANGLEANGLE_COEFFS(5:6,IOSEF1), &                  !
                                       REPEAT(' ',33)//'|';                              !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
    write(icanal,'(a2,3i8,a44)') '| ',               &                                   !
                                 IOSEF1,             &                                   !
                                 NTYPE_IMPROPER_NEW, &                                   !
                                 NTYPE_IMPROPER,     &                                   !
                                 REPEAT(' ',43)//'|';                                    !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### DEALLOCATE ARRAYS CONTAINING OLD PARAMETERS ################################################
                                                                                         !
    call DEALLOCATE_PARAMETER_ARRAY_OLD(icanal,         &                                !
                                        NTYPE_ATOM,     &                                !
                                        NTYPE_BOND,     &                                !
                                        NTYPE_ANGLE,    &                                !
                                        NTYPE_DIHEDRAL, &                                !
                                        NTYPE_IMPROPER);                                 !
                                                                                         !
!   ### ALLOCATE ARRAYS TO UPDATE THE CONFIGURATION ################################################
                                                                                         !
    call ALLOCATE_CONFIG_ARRAY_OLD(icanal,    &                                          !
                                   NATOM,     &                                          !
                                   NBOND,     &                                          !
                                   NANGLE,    &                                          !  
                                   NDIHEDRAL, &                                          !
                                   NIMPROPER);                                           !
                                                                                         !
!   ### UPDATE #####################################################################################
                                                                                         !
    CONFIG_MOLECULEID_OLD(1:NATOM) = CONFIG_MOLECULEID(1:NATOM);                         !
                                                                                         !
    CONFIG_ATOMID_OLD(1:NATOM) = CONFIG_ATOMID(1:NATOM);                                 !
                                                                                         !
    CONFIG_ATOM_TYPE_OLD(1:NATOM) = CONFIG_ATOM_TYPE(1:NATOM);                           !
                                                                                         !
    CONFIG_RI_OLD(1:3,1:NATOM) = CONFIG_RI(1:3,1:NATOM);                                 !
                                                                                         !
    CONFIG_VI_OLD(1:3,1:NATOM) = CONFIG_VI(1:3,1:NATOM);                                 !
                                                                                         !
    CONFIG_QI_OLD(1:NATOM) = CONFIG_QI(1:NATOM);                                         !
                                                                                         !
    CONFIG_NAT_OLD(1:NATOM) = CONFIG_NAT(1:NATOM);                                       !
                                                                                         !
    write(icanal,'(a70)') '| NATOM : '//REPEAT(' ',59)//'|';                             !
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !
    do i = 1, NATOM;                                                                     !
                                                                                         !
        if ( ( CONFIG_MOLECULEID_OLD(i) == 0     ) .AND. &                               !
             ( CONFIG_ATOMID_OLD(i)     == 0     ) .AND. &                               !
             ( CONFIG_ATOM_TYPE_OLD(i)  == 0     ) .AND. &                               !
             ( CONFIG_RI_OLD(1,i)       == 0.0d0 ) .AND. &                               !
             ( CONFIG_RI_OLD(2,i)       == 0.0d0 ) .AND. &                               !
             ( CONFIG_RI_OLD(3,i)       == 0.0d0 ) .AND. &                               !
             ( CONFIG_VI_OLD(1,i)       == 0.0d0 ) .AND. &                               !
             ( CONFIG_VI_OLD(2,i)       == 0.0d0 ) .AND. &                               !
             ( CONFIG_VI_OLD(3,i)       == 0.0d0 ) .AND. &                               !
             ( CONFIG_QI_OLD(i)         == 0.0d0 ) .AND. &                               !
             ( CONFIG_NAT_OLD(i)        == 'XX'  ) ) CYCLE;                              !

!       if ( CONFIG_ATOM_TYPE_OLD(i) > 14 ) CYCLE;
                                                                                         !
        IOSEF1 = IOSEF1 + 1;                                                             !
                                                                                         !
        CONFIG_MOLECULEID(IOSEF1) = CONFIG_MOLECULEID_OLD(i);                            !
                                                                                         !
        CONFIG_ATOMID(IOSEF1)     = CONFIG_ATOMID_OLD(i);                                !
                                                                                         !
        CONFIG_ATOM_TYPE(IOSEF1)  = CONFIG_ATOM_TYPE_OLD(i);                             !
                                                                                         !
        CONFIG_RI(1:3,IOSEF1)     = CONFIG_RI_OLD(1:3,i);                                !
                                                                                         !
        CONFIG_VI(1:3,IOSEF1)     = CONFIG_VI_OLD(1:3,i);                                !
                                                                                         !
        CONFIG_QI(IOSEF1)         = CONFIG_QI_OLD(i);                                    !
                                                                                         !
        CONFIG_NAT(IOSEF1)        = CONFIG_NAT_OLD(i);                                   !

!       if ( CONFIG_ATOM_TYPE(IOSEF1) > 14 ) write(icanal,*) CONFIG_ATOMID(IOSEF1), CONFIG_ATOM_TYPE(IOSEF1);
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
    write(icanal,'(a2,3i8,a44)') '| ', IOSEF1, NATOM_NEW, NATOM, REPEAT(' ',43)//'|';    !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !

!   write(icanal,*) 'XXXX ', MAXVAL( CONFIG_ATOM_TYPE(1:NATOM_NEW) );
!   write(icanal,*) 'XXXX ', MINVAL( CONFIG_ATOM_TYPE(1:NATOM_NEW) );

!   stop


    write(icanal,'(a70)') '| NBOND : '//REPEAT(' ',59)//'|';                             !
                                                                                         !
    BOND_ATOMID_OLD(1:2,1:NBOND)  = BOND_ATOMID(1:2,1:NBOND);                            !
                                                                                         !
    BOND_TYPE_OLD(1:NBOND) = BOND_TYPE(1:NBOND);                                         !
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !
    do i = 1, NBOND;                                                                     !
                                                                                         !
        if ( ( BOND_ATOMID_OLD(1,i) == 0 ) .AND. &                                       !
             ( BOND_ATOMID_OLD(2,i) == 0 ) .AND. &                                       !
             ( BOND_TYPE_OLD(i) == 0 ) ) CYCLE;                                          !
                                                                                         !
        IOSEF1 = IOSEF1 + 1;                                                             !
                                                                                         !
        BOND_ATOMID(1:2,IOSEF1) = BOND_ATOMID_OLD(1:2,i);                                !
                                                                                         !
        BOND_TYPE(IOSEF1) = BOND_TYPE_OLD(i);                                            !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
    write(icanal,'(a2,3i8,a44)') '| ', IOSEF1, NBOND_NEW, NBOND, REPEAT(' ',43)//'|';    !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '| NANGLE : '//REPEAT(' ',58)//'|';                            !
                                                                                         !
    ANGLE_ATOMID_OLD(1:3,1:NANGLE) = ANGLE_ATOMID(1:3,1:NANGLE);                         !
                                                                                         !
    ANGLE_TYPE_OLD(1:NANGLE) = ANGLE_TYPE(1:NANGLE);                                     !
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !
    do i = 1, NANGLE;                                                                    !
                                                                                         !
        if ( ( ANGLE_ATOMID_OLD(1,i) == 0 ) .AND. &                                      !
             ( ANGLE_ATOMID_OLD(2,i) == 0 ) .AND. &                                      !
             ( ANGLE_ATOMID_OLD(3,i) == 0 ) .AND. &                                      !
             ( ANGLE_TYPE_OLD(i) == 0 ) ) CYCLE;                                         !
                                                                                         !
        IOSEF1 = IOSEF1 + 1;                                                             !
                                                                                         !
        ANGLE_ATOMID(1:3,IOSEF1) = ANGLE_ATOMID_OLD(1:3,i);                              !
                                                                                         !
        ANGLE_TYPE(IOSEF1) = ANGLE_TYPE_OLD(i);                                          !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
    write(icanal,'(a2,3i8,a44)') '| ', IOSEF1, NANGLE_NEW, NANGLE, REPEAT(' ',43)//'|';  !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '| NDIHEDRAL : '//REPEAT(' ',55)//'|';                         !
                                                                                         !
    DIHEDRAL_ATOMID_OLD(1:4,1:NDIHEDRAL) = DIHEDRAL_ATOMID(1:4,1:NDIHEDRAL);             !
                                                                                         !
    DIHEDRAL_TYPE_OLD(1:NDIHEDRAL) = DIHEDRAL_TYPE(1:NDIHEDRAL);                         !
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !
    do i = 1, NDIHEDRAL;                                                                 !
                                                                                         !
        if ( ( DIHEDRAL_ATOMID_OLD(1,i) == 0 ) .AND. &                                   !
             ( DIHEDRAL_ATOMID_OLD(2,i) == 0 ) .AND. &                                   !
             ( DIHEDRAL_ATOMID_OLD(3,i) == 0 ) .AND. &                                   !
             ( DIHEDRAL_ATOMID_OLD(4,i) == 0 ) .AND. &                                   !
             ( DIHEDRAL_TYPE(i) == 0 ) ) CYCLE;                                          !
                                                                                         !
        IOSEF1 = IOSEF1 + 1;                                                             !
                                                                                         !
        DIHEDRAL_ATOMID(1:4,IOSEF1) = DIHEDRAL_ATOMID_OLD(1:4,i);                        !
                                                                                         !
        DIHEDRAL_TYPE(IOSEF1) = DIHEDRAL_TYPE(i);                                        !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      ! 
    write(icanal,'(a2,3i8,a44)') '| ', IOSEF1, NDIHEDRAL_NEW, NDIHEDRAL, &               !
                                 REPEAT(' ',43)//'|';                                    !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      ! 
                                                                                         !
    write(icanal,'(a70)') '| NIMPROPER : '//REPEAT(' ',55)//'|';                         !
                                                                                         !
    IMPROPER_ATOMID_OLD(1:4,1:NIMPROPER) = IMPROPER_ATOMID(1:4,1:NIMPROPER);             !
                                                                                         !
    IMPROPER_TYPE_OLD(1:NIMPROPER) = IMPROPER_TYPE(1:NIMPROPER);                         !
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !
    do i = 1, NIMPROPER;                                                                 !
                                                                                         !
        if ( ( IMPROPER_ATOMID_OLD(1,i) == 0 ) .AND. &                                   !
             ( IMPROPER_ATOMID_OLD(2,i) == 0 ) .AND. &                                   !
             ( IMPROPER_ATOMID_OLD(3,i) == 0 ) .AND. &                                   !
             ( IMPROPER_ATOMID_OLD(4,i) == 0 ) .AND. &                                   !
             ( IMPROPER_TYPE_OLD(i) == 0 ) ) CYCLE;                                      !
                                                                                         !
        IOSEF1 = IOSEF1 + 1;                                                             !
                                                                                         !
        IMPROPER_ATOMID(1:4,IOSEF1) = IMPROPER_ATOMID_OLD(1:4,i);                        !
                                                                                         !
        IMPROPER_TYPE(IOSEF1) = IMPROPER_TYPE_OLD(i);                                    !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
    write(icanal,'(a2,3i8,a44)') '| ', IOSEF1, NIMPROPER_NEW, NIMPROPER, &               !
                                 REPEAT(' ',43)//'|';                                    !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### DEALLOCATE ARRAYS ##########################################################################
                                                                                         !
    call DEALLOCATE_CONFIG_ARRAY_OLD(icanal,    &                                        !
                                     NATOM,     &                                        !
                                     NBOND,     &                                        !
                                     NANGLE,    &                                        ! 
                                     NDIHEDRAL, &                                        !
                                     NIMPROPER);                                         !
                                                                                         !
!   ################################################################################################
                                                                                         !
    NATOM     = NATOM_NEW;                                                               !
    NBOND     = NBOND_NEW;                                                               !
    NANGLE    = NANGLE_NEW;                                                              !
    NDIHEDRAL = NDIHEDRAL_NEW;                                                           !
    NIMPROPER = NIMPROPER_NEW;                                                           !
                                                                                         !
    NTYPE_ATOM     = NTYPE_ATOM_NEW;                                                     !
    NTYPE_BOND     = NTYPE_BOND_NEW;                                                     !
    NTYPE_ANGLE    = NTYPE_ANGLE_NEW;                                                    !
    NTYPE_DIHEDRAL = NTYPE_DIHEDRAL_NEW;                                                 !
    NTYPE_IMPROPER = NTYPE_IMPROPER_NEW;                                                 !
                                                                                         !
!   ### ALLOCATE ARRAYS TO CHECK MOLECULES TO REMOVE ###############################################
                                                                                         !
    deallocate(LIST_ATOM_TYPE);                                                          !
                                                                                         !
    deallocate(MOLECULE_MOLAR_MASS);                                                     !
                                                                                         !
    deallocate(MOLECULE_NTYPE_ATOM);                                                     !
                                                                                         !
    deallocate(MOLECULE_NTYPE_BOND);                                                     !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!

end subroutine RANDOM_DELETION_MOLECULE











