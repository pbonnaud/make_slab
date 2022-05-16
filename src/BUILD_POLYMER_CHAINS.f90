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

subroutine BUILD_POLYMER_CHAINS(icanal,PASSA,PASSB) 

!   ************************************************************************************************
!   **                                  READ THE INPUT FILE                                       **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                    : CANAL ON WHICH OUTPUT DATA ARE WRITTEN                         **
!   **                                                                                            **
!   ** PASSA                     :                                                                **
!   **                                                                                            **
!   ** PASSB                     :                                                                **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

    use module_physical_constants;

    use module_size_arrays;

    use module_library;

    use module_config;

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    real (kind=8), dimension(1:3,1:3), intent(in) :: PASSA, PASSB;

!   ************************************************************************************************

    integer (kind=4) :: i, j, k, m, p;

    integer (kind=4) :: iselect, NRemaining_monomers, ifile_lib, ninsertions;

    integer (kind=4) :: iloc_r1, iloc_r2;

    integer (kind=4) :: IFOUND;

    integer (kind=4) :: IINSERTED, IREGIONS_INSERTION, IATOM, IOVERLAPPING;

    integer (kind=4) :: IATOM_TYPE_TRANSLATE, IMOLECULEID;

    integer (kind=4) :: IBOND_TYPE_TRANSLATE, IANGLE_TYPE_TRANSLATE, IIMPROPER_TYPE_TRANSLATE;

    integer (kind=4) :: IDIHEDRAL_TYPE_TRANSLATE;

    integer (kind=4) :: NMonomer_total;

    integer (kind=4) :: ILOOP_MAX1, NLOOP_MAX1;

    real (kind=8) :: grnd;

    real (kind=8) :: mselect;

    real (kind=8) :: SUM_MASSES, INV_SUM_MASSES, RIJ, RR2C2N, RR1C11;

    real (kind=8) :: TAN_THETA1, TAN_THETA2, COS_THETA1, COS_THETA2, THETA1_RAD, THETA2_RAD;

    real (kind=8) :: COS_PHI1, COS_PHI2, PHI1_RAD, PHI2_RAD;

    real (kind=8) :: NORMU, NORMV;

    real (kind=8) :: NORM_RC2NR2_OLD, NORM_R1C11_NEW;

    real (kind=8) :: DELTA_THETA, DELTA_PHI, DELTA_PSI;

    real (kind=8) :: NEW_THETA, NEW_PHI, NEW_PSI;

    real (kind=8) :: COS_DELTA_THETA, SIN_DELTA_THETA, COS_DELTA_PHI, SIN_DELTA_PHI;

    real (kind=8) :: COS_DELTA_PSI, SIN_DELTA_PSI;

    real (kind=8) :: MAX_RAND_THETA, MAX_RAND_PHI, MAX_RAND_PSI, RCUT1, RCUT2;

    real (kind=8) :: RAND_THETA, RAND_PHI, RAND_PSI;

    real (kind=8), dimension(1:3) :: RINO, RI, RAND_ANGLES, RG, DRIJ;

    real (kind=8), dimension(1:3) :: R1_REF, R2_REF, R2_REF_OLD, RC11, RC2N, RC2N_OLD; 

    real (kind=8), dimension(1:3) :: RC2NR2_OLD, RR1C11_NEW;

    real (kind=8), dimension(1:3,1:3) :: MATROT,           &
                                         MATRIX_ROTATION1, &
                                         MATRIX_ROTATION2, &
                                         MATRIX_ROTATION3;

    integer (kind=4), allocatable, dimension(:) :: NMonomer_units_per_chain_local;

    real (kind=8), allocatable, dimension(:) :: Insertion_probability, Insertion_probability_local;

!   ************************************************************************************************

    integer (kind=4) :: NATOM_NEW, NBOND_NEW, NANGLE_NEW, NDIHEDRAL_NEW, NIMPROPER_NEW;

    integer (kind=4) :: NTYPE_ATOM_NEW,     &
                        NTYPE_BOND_NEW,     &
                        NTYPE_ANGLE_NEW,    &
                        NTYPE_DIHEDRAL_NEW, &
                        NTYPE_IMPROPER_NEW;

!   ************************************************************************************************

!   integer (kind=4) :: IOSEF1, IOSEF2, IOSEF3, IOSEF4, IOSEF5, IOSEF6;

!   real (kind=8) :: ROSEF1;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the routine #############################################################
                                                                                         !
    CHTITLE = 'Build polymer chains';                                                    !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialization of the random generator #####################################################
                                                                                         !
!   call sgrnd(RANDOM_GENERATOR_SEED);                                                   !
                                                                                         !
!   write(icanal,'(a70)') '| The random generator was initialized'//REPEAT(' ',31)//'|'; !
                                                                                         !
!   write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set the maximum number of loop for each insertion try ######################################
                                                                                         !
    NLOOP_MAX1 = 1000;                                                                   !
                                                                                         !
!   ### Generation of polymer chains ###############################################################
                                                                                         !
    NATOM_NEW = 0;                                                                       !
                                                                                         !
    do i = 1, NGenerate_polymer_species;                                                 !
                                                                                         !
        write(icanal,'(a33,i4,a33)') '| Insertion of polymer species # ', i, &           !
                                     ' '//REPEAT('-',31)//'|';                           !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       ### Write the number of monomer species contained by this polymer chain ####################
                                                                                         !
        write(icanal,'(a11,i4,a55)') '| There is ',                                 &    !
                                      NSPECIES_PER_POLYMER_CHAINS(i),               &    !
                                      ' monomer species for that polymer species'// &    !
                                      REPEAT(' ',13)//'|';                               !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       ### Set and write maximum random angles for the generation of polymer chains ###############
                                                                                         !
        MAX_RAND_THETA = MONOMER_INSERTION_PARAM(i,1);                                   ! IN [DEG]
                                                                                         !
        MAX_RAND_PHI = MONOMER_INSERTION_PARAM(i,2);                                     ! IN [DEG]
                                                                                         !
        MAX_RAND_PSI = MONOMER_INSERTION_PARAM(i,3);                                     ! IN [DEG]
                                                                                         ! 
        MAX_RAND_THETA = MAX_RAND_THETA * DEG_TO_RAD;                                    ! IN [RAD]
                                                                                         !
        MAX_RAND_PHI = MAX_RAND_PHI * DEG_TO_RAD;                                        ! IN [RAD]
                                                                                         !
        MAX_RAND_PSI = MAX_RAND_PSI * DEG_TO_RAD;                                        ! IN [RAD]
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
        RCUT1 = MONOMER_INSERTION_PARAM(i,4);                                            ! IN [A]
                                                                                         !
        RCUT2 = MONOMER_INSERTION_PARAM(i,5);                                            ! IN [A]
                                                                                         !
!       ### Allocation and initialization of local arrays containing polymer chain properties ######
                                                                                         !
        IOSEF1 = NSPECIES_PER_POLYMER_CHAINS(i);                                         !
                                                                                         !
        allocate(NMonomer_units_per_chain_local(1:IOSEF1));                              !
                                                                                         !
        NMonomer_units_per_chain_local(1:IOSEF1) = 0;                                    !
                                                                                         !
!       ### Compute and write the initial total number of monomers per polymer chains ##############
                                                                                         !
        IOSEF1 = NSPECIES_PER_POLYMER_CHAINS(i);                                         !
                                                                                         !
        NMonomer_total = SUM(NMONOMER_UNITS_PER_CHAINS(i,1:IOSEF1));                     !
                                                                                         !
        write(icanal,'(a24,i4,a42)') '| This species contains ', &                       !
                                     NMonomer_total,             &                       !
                                     ' monomers to insert'//     &                       !
                                     REPEAT(' ',22)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       ### Allocation and initialization of arrays for insertion probability ######################
                                                                                         !
        allocate(Insertion_probability(1:NSPECIES_PER_POLYMER_CHAINS(i)));               !
                                                                                         !
        allocate(Insertion_probability_local(1:NSPECIES_PER_POLYMER_CHAINS(i)));         !
                                                                                         !
        Insertion_probability(1:NSPECIES_PER_POLYMER_CHAINS(i)) = 0.0d0;                 !
                                                                                         !
        Insertion_probability_local(1:NSPECIES_PER_POLYMER_CHAINS(i)) = 0.0d0;           !
                                                                                         !
!       ### Compute the insertion probability for each monomer species #############################
                                                                                         !
        do j = 1, NSPECIES_PER_POLYMER_CHAINS(i);                                        !
                                                                                         !
            ROSEF1 = REAL( NMONOMER_UNITS_PER_CHAINS(i,j) ) / REAL( NMonomer_total );    !
                                                                                         !
            if ( j > 1 ) Insertion_probability(j) = Insertion_probability(j-1);          !
                                                                                         ! 
            Insertion_probability(j) = Insertion_probability(j) + ROSEF1;                !
                                                                                         !
            write(icanal,'(a36,i4,a3,2f8.4,a11)') '| Insertion probability '// &         !
                                                  'for monomer ',              &         !
                                                  j,                           &         !
                                                  ' : ',                       &         !
                                                  ROSEF1,                      &         !
                                                  Insertion_probability(j),    &         !
                                                  REPEAT(' ',10)//'|';                   ! 
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        end do                                                                           !
                                                                                         !
!       ### Insertion of the polymer chains of the same species ####################################
                                                                                         !
        do j = 1, NPOLYMER_CHAINS(i);                                                    !
                                                                                         !
            write(icanal,'(a26,i4,a40)') '| Insertion of molecule # ', &                 !
                                         j,                            &                 !
                                         ' '//REPEAT('/',38)//'|';                       !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           ### Initialization of the remaining number of monomers to insert #######################
                                                                                         !
            IOSEF1 = NSPECIES_PER_POLYMER_CHAINS(i);                                     !
                                                                                         !
            NMonomer_units_per_chain_local(1:IOSEF1) = &                                 !
            NMONOMER_UNITS_PER_CHAINS(i,1:IOSEF1);                                       !
                                                                                         !
            NRemaining_monomers = SUM(NMonomer_units_per_chain_local(1:IOSEF1));         !
                                                                                         !
            write(icanal,'(a11,i4,a55)') '| There is ',          &                       !
                                         NRemaining_monomers,    &                       !
                                         ' monomers to insert'// &                       !
                                         REPEAT(' ',35)//'|';                            !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           ### Initialization of local insertion probabilities ####################################
                                                                                         !
            IOSEF1 = NSPECIES_PER_POLYMER_CHAINS(i);                                     ! 
                                                                                         !
            Insertion_probability_local(1:IOSEF1) = Insertion_probability(1:IOSEF1);     !
                                                                                         !
!           ### Loop until the given number of monomers is inserted ################################
                                                                                         !
            R2_REF(1:3) = 0.0d0;                                                         !
                                                                                         !
            RC11(1:3) = 0.0d0;                                                           !
                                                                                         !
            RC2N(1:3) = 0.0d0;                                                           !
                                                                                         !
            ninsertions = 0;                                                             !
                                                                                         !
            do while ( NRemaining_monomers > 0 );                                        ! 
                                                                                         !
                write(icanal,'(a24,i4,a42)') '| NRemaining_monomers : ', &               !
                                             NRemaining_monomers,        &               !
                                             ' '//REPEAT('x',40)//'|';                   !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               ### Determine randomly which monomer species to insert #############################
                                                                                         ! 
                mselect = grnd();                                                        !
                                                                                         !
                ifile_lib = 0;                                                           !
                                                                                         !
                do k = 1, NSPECIES_PER_POLYMER_CHAINS(i);                                !
                                                                                         !
                    ifile_lib = ifile_lib + 1;                                           !
                                                                                         ! 
                    if ( mselect < Insertion_probability_local(k) ) then;                !
                                                                                         !
                        iselect = k;                                                     !
                                                                                         !
                        write(icanal,'(a19,2f8.4,a35)') '| Chosen species : ',                & !
                                                        mselect,                              & !
                                                        Insertion_probability_local(iselect), & !
                                                        REPEAT(' ',34)//'|';             !
                                                                                         !
                        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                  !
                                                                                         ! 
                        EXIT;                                                            !
                                                                                         !
                    end if                                                               !
                                                                                         !
                end do                                                                   !
                                                                                         !
!               if ( ninsertions == 1 ) stop; !//////////////////////////////////////////!
                                                                                         !
                write(icanal,'(a28,2i4,a34)') '| Location in library files ', &          !
                                              iselect,                        &          !
                                              ifile_lib,                      &          !
                                              REPEAT(' ',33)//'|';                       ! 
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               if ( ninsertions == 1 ) stop; !//////////////////////////////////////////!
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
!               ### Get the location of R1 and R2 in the current library file ######################
                                                                                         !
                iloc_r1 = ILOCR1_LIBRARY(ifile_lib);                                     !
                                                                                         !
                iloc_r2 = ILOCR2_LIBRARY(ifile_lib);                                     ! 
                                                                                         !
!               ### Look for the location of the first radical in the library array ################
                                                                                         !
!               iloc_r1 = 0;                                                             !
                                                                                         !
!               do m = 1, NATOM_LIBRARY(ifile_lib);                                      !
                                                                                         !
!                   if ( CONFIG_NAT_LIBRARY(m,ifile_lib) /= 'R1' ) CYCLE;                !
                                                                                         !
!                   iloc_r1 = m;                                                         !
                                                                                         !
!               end do                                                                   !
                                                                                         !
!               if ( iloc_r1 == 0 ) then;                                                !
                                                                                         !
!                   write(icanal,'(a70)') '| No radical R1 found '// &                   !
!                                         'in the library file'//    &                   !
!                                         REPEAT(' ',28)//'|';                           !
                                                                                         !
!                   write(icanal,'(a70)') '| The file is maybe not elligible '// &       !
!                                         'for this routine'//REPEAT(' ',19)//'|';       !
                                                                                         !
!                   write(icanal,'(a14)') 'End of Program';                              !
                                                                                         !
!                   close(icanal);                                                       !
                                                                                         !
!                   stop; !//////////////////////////////////////////////////////////////!
                                                                                         !
!               end if                                                                   !
                                                                                         !
!               ### Look for the location of the second radical in the library array ###############
                                                                                         !
!               iloc_r2 = 0;                                                             !
                                                                                         !
!               do m = 1, NATOM_LIBRARY(ifile_lib);                                      !
                                                                                         !
!                   if ( CONFIG_NAT_LIBRARY(m,ifile_lib) /= 'R2' ) CYCLE;                !
                                                                                         !
!                   iloc_r2 = m;                                                         !
                                                                                         !
!               end do                                                                   !
                                                                                         !
!               if ( iloc_r2 == 0 ) then;                                                !
                                                                                         !
!                   write(icanal,'(a70)') '| No radical R2 found in '// &                !
!                                         'the library file'//          &                !
!                                         REPEAT(' ',28)//'|';                           !
                                                                                         !
!                   write(icanal,*) '| The file is maybe not elligible for this routine';
                                                                                         !
!                   write(icanal,'(a14)') 'End of Program';                              !
                                                                                         !
!                   close(icanal);                                                       !
                                                                                         !
!                   stop; !//////////////////////////////////////////////////////////////!
                                                                                         !
!               end if                                                                   !
                                                                                         !
!               ### Set coordinates of radical R1 of the current monomer ###########################
                                                                                         !
                if ( ninsertions == 0 ) then;                                            !
                                                                                         !
                    IFOUND = 1;                                                          !
                                                                                         ! 
                    do while ( IFOUND == 1 )                                             !
                                                                                         !
                        RINO(1:3) = (/grnd() - 0.5d0, grnd() - 0.5d0, grnd() - 0.5d0/);  ! Generate random scaled coordinates 
                                                                                         !
                        RI(1:3) = MATMUL(PASSA(1:3,1:3),RINO(1:3));                      ! Convert scaled coordinates into Cartesian coordinates
                                                                                         !
                        if ( NATOM_NEW > 0 ) then;                                       !
                                                                                         !
                            IOVERLAPPING = 0;                                            !
                                                                                         !
                            do m = 1, NATOM_NEW;                                         !
                                                                                         !
                                DRIJ(1:3) = CONFIG_RI(1:3,m) - RI(1:3);                  !
                                                                                         !
                                call APPLY_PBC(DRIJ(1:3),      &                         !
                                               DRIJ(1:3),      &                         !
                                               PASSA(1:3,1:3), &                         !
                                               PASSB(1:3,1:3));                          !
                                                                                         !
                                RIJ = DSQRT( DOT_PRODUCT(DRIJ(1:3),DRIJ(1:3)) );         !
                                                                                         !
                                if ( RIJ <= RCUT1 ) then;                                !
                                                                                         !
                                    IOVERLAPPING = IOVERLAPPING + 1;                     !
                                                                                         !
                                    EXIT;                                                !
                                                                                         !
                                end if                                                   !
                                                                                         !
                            end do                                                       !
                                                                                         !
                            if ( IOVERLAPPING == 0 ) IFOUND = 0;                         !
                                                                                         !
                        else                                                             !
                                                                                         !
                            IFOUND = 0;                                                  !
                                                                                         !
                        end if                                                           !
                                                                                         !
                    end do                                                               !
                                                                                         !
                    R1_REF(1:3) = RI(1:3);                                               ! Compute coordinates in Cartesian axis
                                                                                         !
                else                                                                     !
                                                                                         !
                    R1_REF(1:3) = RC2N_OLD(1:3);                                         !
                                                                                         !
                end if                                                                   !
                                                                                         !
                write(icanal,'(a21,3f15.6,a4)') '| Set R1_REF (new) : ', &               !
                                                R1_REF(1:3),             &               !
                                                REPEAT(' ',3)//'|';                      !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               if ( ninsertions == 1 ) stop; !//////////////////////////////////////////!
                                                                                         !
!               ### Look for the coordinates of the carbon atom bonded to R1 in the lib file #######
                                                                                         !
                RR1C11 = 10000000;                                                       !
                                                                                         !
                do m = 1, NATOM_LIBRARY(ifile_lib);                                      !
                                                                                         !
                    if ( CONFIG_NAT_LIBRARY(m,ifile_lib) /= 'C' ) CYCLE;                 !
                                                                                         !
                    DRIJ(1:3) = CONFIG_RI_LIBRARY(1:3,m,ifile_lib) - &                   !
                                CONFIG_RI_LIBRARY(1:3,iloc_r1,ifile_lib);                !
                                                                                         !
                    call APPLY_PBC(DRIJ(1:3),DRIJ(1:3),PASSA(1:3,1:3),PASSB(1:3,1:3));   !
                                                                                         !
                    RIJ = DSQRT( DOT_PRODUCT(DRIJ(1:3),DRIJ(1:3)) );                     !
                                                                                         !
                    if ( RIJ < RR1C11 ) then;                                            !
                                                                                         !
                        RR1C11 = RIJ;                                                    !
                                                                                         !
                        RC11(1:3) = DRIJ(1:3);                                           !
                                                                                         !
                    end if                                                               !
                                                                                         !
                end do                                                                   !
                                                                                         !
                RC11(1:3) = RC11(1:3) + R1_REF(1:3);                                     !
                                                                                         !
                write(icanal,'(a19,3f15.6,a6)') '| Set (new) RC11 : ', &                 !
                                                RC11(1:3),             &                 !
                                                REPEAT(' ',5)//'|';                      !
                                                                                         ! 
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               if ( ninsertions == 1 ) stop; !//////////////////////////////////////////!
                                                                                         !
!               ### Look for the coordinates of the carbon atom bonded to R2 in the lib file #######
                                                                                         !
                RR2C2N = 10000000;                                                       !
                                                                                         !
                do p = 1, NATOM_LIBRARY(ifile_lib);                                      !
                                                                                         !
                    if ( CONFIG_NAT_LIBRARY(p,ifile_lib) /= 'C' ) CYCLE;                 !
                                                                                         !
                    DRIJ(1:3) = CONFIG_RI_LIBRARY(1:3,p,ifile_lib) - &                   !
                                CONFIG_RI_LIBRARY(1:3,iloc_r2,ifile_lib);                !
                                                                                         !
                    call APPLY_PBC(DRIJ(1:3),DRIJ(1:3),PASSA(1:3,1:3),PASSB(1:3,1:3));   !
                                                                                         !
                    RIJ = DSQRT( DOT_PRODUCT(DRIJ(1:3),DRIJ(1:3)) );                     !
                                                                                         !
                    if ( RIJ < RR2C2N ) then;                                            !
                                                                                         !
                        RR2C2N = RIJ;                                                    !
                                                                                         !
                        RC2N(1:3) = CONFIG_RI_LIBRARY(1:3,p,ifile_lib) - &               !
                                    CONFIG_RI_LIBRARY(1:3,iloc_r1,ifile_lib);            !
                                                                                         !
                    end if                                                               !
                                                                                         !
                end do                                                                   !
                                                                                         !
                RC2N(1:3) = RC2N(1:3) + R1_REF(1:3);                                     ! 
                                                                                         ! 
                write(icanal,'(a19,3f15.6,a6)') '| Set (new) RC2N : ', &                 !
                                                RC2N(1:3),             &                 !
                                                REPEAT(' ',5)//'|';                      !
                                                                                         ! 
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               if ( ninsertions == 1 ) stop; !//////////////////////////////////////////!
                                                                                         !
!               ### Compute the vector between R1 and C11 ##########################################
                                                                                         !
                if ( ninsertions > 0 ) then;                                             !
                                                                                         !
                    RR1C11_NEW(1:3) = RC11(1:3) - R1_REF(1:3);                           !
                                                                                         !
                    call APPLY_PBC(RR1C11_NEW(1:3),  &                                   ! 
                                   RR1C11_NEW(1:3),  &                                   !
                                   PASSA(1:3,1:3),   &                                   !
                                   PASSB(1:3,1:3));                                      !
                                                                                         !
                end if                                                                   !
                                                                                         !
!               ### Determine the angle between RC2NR2 and RR1C11 ##################################
                                                                                         !           
                if ( ninsertions > 0 ) then;                                             !
                                                                                         !
                    write(icanal,'(a15,3f12.4,a19)') '| RR1C11_NEW : ',  &               !
                                                     RR1C11_NEW(1:3),    &               !
                                                     REPEAT(' ',18)//'|';                !
                                                                                         !
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
                    write(icanal,'(a15,3f12.4,a19)') '| RC2NR2_OLD : ',  &               !
                                                     RC2NR2_OLD(1:3),    &               ! 
                                                     REPEAT(' ',18)//'|';                !
                                                                                         !
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
                    NORM_R1C11_NEW = DOT_PRODUCT(RR1C11_NEW(1:3),RR1C11_NEW(1:3));       !
                                                                                         !
                    NORM_R1C11_NEW = DSQRT( NORM_R1C11_NEW );                            !
                                                                                         !
                    NORM_RC2NR2_OLD = DOT_PRODUCT(RC2NR2_OLD(1:3),RC2NR2_OLD(1:3));      ! 
                                                                                         !
                    NORM_RC2NR2_OLD = DSQRT( NORM_RC2NR2_OLD );                          !
                                                                                         !
                    write(icanal,'(a19,f12.4,a39)') '| NORM_R1C11_NEW : ', &             !
                                                    NORM_R1C11_NEW,        &             !
                                                    REPEAT(' ',38)//'|';                 !
                                                                                         !
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
                    write(icanal,'(a20,f12.4,a38)') '| NORM_RC2NR2_OLD : ', &            !
                                                    NORM_RC2NR2_OLD,        &            !
                                                    REPEAT(' ',37)//'|';                 !
                                                                                         !
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
                    NORMU = DSQRT( RR1C11_NEW(1) * RR1C11_NEW(1) + &                     !
                                   RR1C11_NEW(2) * RR1C11_NEW(2) );                      !
                                                                                         !
                    NORMV = DSQRT( RC2NR2_OLD(1) * RC2NR2_OLD(1) + &                     !
                                   RC2NR2_OLD(2) * RC2NR2_OLD(2) );                      !
                                                                                         !
                    COS_THETA1 = RR1C11_NEW(1) / NORMU;                                  !
                                                                                         !
                    COS_THETA2 = RC2NR2_OLD(1) / NORMV;                                  !
                                                                                         ! 
                    THETA1_RAD = ACOS( COS_THETA1 );                                     !
                                                                                         !
                    THETA2_RAD = ACOS( COS_THETA2 );                                     !
                                                                                         !
                    COS_PHI1 = RR1C11_NEW(3) / NORM_R1C11_NEW;                           !
                                                                                         !
                    COS_PHI2 = RC2NR2_OLD(3) / NORM_RC2NR2_OLD;                          !
                                                                                         !
                    PHI1_RAD = ACOS( COS_PHI1 );                                         !
                                                                                         !
                    PHI2_RAD = ACOS( COS_PHI2 );                                         !
                                                                                         !
                    DELTA_THETA = THETA1_RAD - THETA2_RAD;                               !
                                                                                         !
                    DELTA_PHI = PHI1_RAD - PHI2_RAD;                                     !
                                                                                         !
                    DELTA_PSI = 0.0d0;                                                   !
                                                                                         !
                    write(icanal,'(a16,2f12.4,a30)') '| DELTA_THETA : ',       &         !
                                                     DELTA_THETA,              &         !
                                                     DELTA_THETA * RAD_TO_DEG, &         !
                                                     REPEAT(' ',29)//'|';                !
                                                                                         !
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
                    write(icanal,'(a14,2f12.4,a32)') '| DELTA_PHI : ',       &           !
                                                     DELTA_PHI,              &           !
                                                     DELTA_PHI * RAD_TO_DEG, &           !
                                                     REPEAT(' ',31)//'|';                !
                                                                                         !
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
                    write(icanal,'(a14,2f12.4,a32)') '| DELTA_PSI : ',       &           !
                                                     DELTA_PSI,              &           !
                                                     DELTA_PSI * RAD_TO_DEG, &           !
                                                     REPEAT(' ',31)//'|';                !
                                                                                         !
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
                end if                                                                   !
                                                                                         !
!               if ( ninsertions == 1 ) stop; !//////////////////////////////////////////!
                                                                                         !
!               ### Initialization of variables for the insertion of monomers ######################
                                                                                         ! 
                IATOM = NATOM_NEW;                                                       !
                                                                                         !
                NEW_THETA = 0.0d0;                                                       !
                                                                                         !
                NEW_PHI = 0.0d0;                                                         !
                                                                                         !
                NEW_PSI = 0.0d0;                                                         !
                                                                                         !
                IOVERLAPPING = 1;                                                        !
                                                                                         !
                ILOOP_MAX1 = 0;                                                          !
                                                                                         !
                do while ( IOVERLAPPING == 1 )                                           !
                                                                                         !
                    write(icanal,'(a10,i8,a13,i8,a31)') '| IATOM : ',        &           !
                                                        IATOM,               &           !
                                                        ' NATOM_NEW : ',     &           !
                                                        NATOM_NEW,           &           !  
                                                        REPEAT(' ',30)//'|';             !
                                                                                         !
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
!                   ### Set random angles applied to the current monomer ###############################
                                                                                         !
                    if ( ninsertions > 0 ) then;                                         !
                                                                                         !
                        RAND_THETA = MAX_RAND_THETA * 2.0d0 * ( grnd() - 0.5d0 );        !
                                                                                         !
                        RAND_PHI = MAX_RAND_PHI * 2.0d0 * ( grnd() - 0.5d0 );            !
                                                                                         !
                        RAND_PSI = MAX_RAND_PSI * 2.0d0 * ( grnd() - 0.5d0 );            !
                                                                                         !
                    end if                                                               !
                                                                                         !
!                   ### Update rotation angles to apply to the current monomer #########################
                                                                                         !
                    if ( ninsertions > 0 ) then;                                         !
                                                                                         !
                        NEW_THETA = DELTA_THETA + RAND_THETA;                            !
                                                                                         !
                        NEW_PHI = DELTA_PHI + RAND_PHI;                                  ! 
                                                                                         !
                        NEW_PSI = DELTA_PSI + RAND_PSI;                                  ! 
                                                                                         !
                        write(icanal,'(a14,2f12.4,a32)') '| NEW_THETA : ',       &       !
                                                         NEW_THETA,              &       !
                                                         NEW_THETA * RAD_TO_DEG, &       !
                                                         REPEAT(' ',31)//'|';            !
                                                                                         !
                        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                  ! 
                                                                                         !
                        write(icanal,'(a12,2f12.4,a34)') '| NEW_PHI : ',       &         !
                                                         NEW_PHI,              &         !
                                                         NEW_PHI * RAD_TO_DEG, &         !
                                                         REPEAT(' ',33)//'|';            !
                                                                                         !
                        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                  !
                                                                                         !
                        write(icanal,'(a12,2f12.4,a34)') '| NEW_PSI : ',       &         !
                                                         NEW_PSI,              &         !
                                                         NEW_PSI * RAD_TO_DEG, &         !
                                                         REPEAT(' ',33)//'|';            !
                                                                                         !
                        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                  !
                                                                                         !
                    end if                                                               !
                                                                                         !
!                   if ( ninsertions == 1 ) stop; !//////////////////////////////////////!
                                                                                         !
!                   ### Set the first rotation matrix ##############################################
                                                                                         !
                    COS_DELTA_THETA = DCOS( NEW_THETA );                                 !  
                                                                                         !
                    SIN_DELTA_THETA = DSIN( NEW_THETA );                                 !
                                                                                         !
                    MATRIX_ROTATION1(1,1) =   COS_DELTA_THETA;                           !
                    MATRIX_ROTATION1(2,1) = - SIN_DELTA_THETA;                           !
                    MATRIX_ROTATION1(3,1) =   0.0d0;                                     !
                    MATRIX_ROTATION1(1,2) =   SIN_DELTA_THETA;                           !
                    MATRIX_ROTATION1(2,2) =   COS_DELTA_THETA;                           !
                    MATRIX_ROTATION1(3,2) =   0.0d0;                                     !
                    MATRIX_ROTATION1(1,3) =   0.0d0;                                     !
                    MATRIX_ROTATION1(2,3) =   0.0d0;                                     !
                    MATRIX_ROTATION1(3,3) =   1.0d0;                                     !
                                                                                         !
!                   ### Set the second rotation matrix #############################################
                                                                                         !
                    COS_DELTA_PHI = DCOS( NEW_PHI );                                     !
                                                                                         !
                    SIN_DELTA_PHI = DSIN( NEW_PHI );                                     !
                                                                                         !
                    MATRIX_ROTATION2(1,1) =   COS_DELTA_PHI;                             !
                    MATRIX_ROTATION2(2,1) =   0.0d0;                                     !
                    MATRIX_ROTATION2(3,1) =  -SIN_DELTA_PHI;                             !
                    MATRIX_ROTATION2(1,2) =   0.0d0;                                     !  
                    MATRIX_ROTATION2(2,2) =   1.0d0;                                     !
                    MATRIX_ROTATION2(3,2) =   0.0d0;                                     !
                    MATRIX_ROTATION2(1,3) =   SIN_DELTA_PHI;                             !
                    MATRIX_ROTATION2(2,3) =   0.0d0;                                     !
                    MATRIX_ROTATION2(3,3) =   COS_DELTA_PHI;                             !
                                                                                         !
!                   ### Set the third rotation matrix ##############################################
                                                                                         !
                    COS_DELTA_PSI = DCOS( NEW_PSI );                                     !
                                                                                         !
                    SIN_DELTA_PSI = DSIN( NEW_PSI );                                     !
                                                                                         !
                    MATRIX_ROTATION3(1,1) =   1.0d0;                                     !
                    MATRIX_ROTATION3(2,1) =   0.0d0;                                     !
                    MATRIX_ROTATION3(3,1) =   0.0d0;                                     !
                    MATRIX_ROTATION3(1,2) =   0.0d0;                                     !
                    MATRIX_ROTATION3(2,2) =   COS_DELTA_PSI;                             !
                    MATRIX_ROTATION3(3,2) =  -SIN_DELTA_PSI;                             !
                    MATRIX_ROTATION3(1,3) =   0.0d0;                                     !
                    MATRIX_ROTATION3(2,3) =   SIN_DELTA_PSI;                             !
                    MATRIX_ROTATION3(3,3) =   COS_DELTA_PSI;                             !
                                                                                         !
!                   if ( ninsertions == 1 ) stop; !//////////////////////////////////////!
                                                                                         !
!                   ### Insertion of the monomer ###################################################
                                                                                         !
                    do m = 1, NATOM_LIBRARY(ifile_lib);                                  !
                                                                                         !
                        if ( CONFIG_NAT_LIBRARY(m,ifile_lib) == 'XXX' ) CYCLE;           !
                                                                                         !
                        DRIJ(1:3) = CONFIG_RI_LIBRARY(1:3,m,ifile_lib) - &               !
                                    CONFIG_RI_LIBRARY(1:3,iloc_r1,ifile_lib);            !
                                                                                         !
                        NATOM_NEW = NATOM_NEW + 1;                                       ! Increment the number of atoms in the whole configuration
                                                                                         !
!                       write(icanal,*) 'Before : ', DRIJ(1:3), ninsertions;

!                       ### Application of rotations ###############################################
                                                                                         !
                        if ( ninsertions > 0 ) then;                                     !
                                                                                         !
                            DRIJ(1:3) = MATMUL(MATRIX_ROTATION1(1:3,1:3),DRIJ(1:3));     !
                                                                                         !
                            DRIJ(1:3) = MATMUL(MATRIX_ROTATION2(1:3,1:3),DRIJ(1:3));     !
                                                                                         !
                            DRIJ(1:3) = MATMUL(MATRIX_ROTATION3(1:3,1:3),DRIJ(1:3));     !
                                                                                         !
                        end if                                                           !
                                                                                         !
!                       write(icanal,*) 'After : ', DRIJ(1:3);

                        CONFIG_RI(1:3,NATOM_NEW) = R1_REF(1:3) + DRIJ(1:3);              ! Set coordinate of the first inserted atom
                                                                                         !
                        CONFIG_NAT(NATOM_NEW) = CONFIG_NAT_LIBRARY(m,ifile_lib);         !
                                                                                         !
                        CONFIG_ATOMID(NATOM_NEW) = NATOM_NEW;                            ! 
                                                                                         !
!                       ### Replace radical R1 by a hydrogen atom ##################################
                                                                                         !
                        if ( CONFIG_NAT_LIBRARY(m,ifile_lib) == 'R1' ) then;             !
                                                                                         !
                            if ( ninsertions == 0 ) then;                                !
                                                                                         !
                                CONFIG_NAT(NATOM_NEW) = 'H';                             ! Replace radical R1 by a H atom at the beginning of the polymer chain  
                                                                                         !
                            else                                                         !
                                                                                         !
                                NATOM_NEW = NATOM_NEW - 1;                               ! 
                                                                                         !
                            end if                                                       !
                                                                                         !
                        end if                                                           !
                                                                                         !
!                       ### Replace the radical R2 by a hydrogen atom in the carbon chain ##########
                                                                                         !
                        if ( CONFIG_NAT_LIBRARY(m,ifile_lib) == 'R2' ) then;             !
                                                                                         ! 
                            R2_REF(1:3) = CONFIG_RI(1:3,NATOM_NEW);                      !
                                                                                         !
                            if ( NRemaining_monomers == 1 ) then;                        !
                                                                                         !
                                CONFIG_NAT(NATOM_NEW) = 'H';                             !
                                                                                         !
                            end if                                                       !
                                                                                         !
                        end if                                                           !
                                                                                         !
                    end do                                                               !
                                                                                         !
!                   ### Check overlapping of atoms in the configuration ############################
                                                                                         !
                    call CHECK_OVERLAPPING_ATOMS(icanal,          &                      !
                                                 IATOM+1,         &                      !
                                                 NATOM_NEW,       &                      !
                                                 PASSA(1:3,1:3),  &                      !
                                                 PASSB(1:3,1:3),  &                      !
                                                 RCUT2,          &                       !
                                                 IOVERLAPPING);                          !
                                                                                         !
                    if ( IOVERLAPPING == 1 ) then;                                       ! 
                                                                                         !
                        NATOM_NEW = IATOM;                                               !
                                                                                         !
                        write(icanal,'(a70)') '| There is overlap !!!'// &               !
                                              REPEAT(' ',47)//'|';                       !
                                                                                         !
                        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                  !
                                                                                         !
                    end if                                                               !
                                                                                         !
!                   if ( ninsertions == 1 ) stop; !//////////////////////////////////////!
                                                                                         !
                    ILOOP_MAX1 = ILOOP_MAX1  + 1;                                        !
                                                                                         !
                    if ( ILOOP_MAX1 > NLOOP_MAX1 ) then;                                 !
                                                                                         !
                        write(icanal,*) '| ILOOP_MAX1 | NLOOP_MAX1 : ', ILOOP_MAX1, NLOOP_MAX1;
                                                                                         !
                        stop; !//////////////////////////////////////////////////////////!
                                                                                         !
                    end if                                                               !
                                                                                         !
                end do                                                                   !
                                                                                         !
!               if ( ninsertions == 1 ) stop; !//////////////////////////////////////////!
                                                                                         !
!               ### Find the location of R2 in the config array ####################################
                                                                                         ! 
                IOSEF1 = 0;                                                              !
                                                                                         !
                do m = 1, NATOM_NEW;                                                     !
                                                                                         !
                    if ( CONFIG_NAT(m) /= 'R2' ) CYCLE;                                  !
                                                                                         !
                    IOSEF1 = m;                                                          !
                                                                                         !
                    EXIT;                                                                !
                                                                                         !
                end do                                                                   !
                                                                                         !
                write(icanal,'(a32,i8,a30)') '| Location of R2 in the array : ', &       !
                                             IOSEF1,                             &       !
                                             REPEAT(' ',29)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                if ( IOSEF1 > 0 ) then;                                                  !
                                                                                         !
                    R2_REF_OLD(1:3) = CONFIG_RI(1:3,IOSEF1);                             !
                                                                                         !
                    IOSEF6 = 70 - 6 - 4 - 1 - 1 - LEN_TRIM(CONFIG_NAT(IOSEF1));          !
                                                                                         !
                    write(icanal,'(a6,i4,a60)') '| --> ',                       &        !  
                                                IOSEF1,                         &        !  
                                                ' '//TRIM(CONFIG_NAT(IOSEF1))// &        !
                                                REPEAT(' ',IOSEF6)//'|';                 !
                                                                                         !
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
                    write(icanal,'(a15,3f12.4,a19)') '| R2_REF_OLD : ', &                !
                                                     R2_REF_OLD(1:3),   &                !
                                                     REPEAT(' ',18)//'|';                !
                                                                                         !
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
                end if                                                                   !
                                                                                         !
!               if ( ninsertions == 1 ) stop; !//////////////////////////////////////////!
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
!               ### Find coordinates of the carbon atom bonded to R2 ###############################
                                                                                         !
                if ( NRemaining_monomers > 1 .AND. IOSEF1 > 0 ) then;                    ! 
                                                                                         !
                    RR2C2N = 10000000;                                                   !
                                                                                         !
                    RC2N_OLD(1:3) = 0.0d0;                                               !
                                                                                         !
                    do p = 1, NATOM_NEW;                                                 !
                                                                                         !
                        if ( CONFIG_NAT(p) /= 'C' ) CYCLE;                               !
                                                                                         !
                        DRIJ(1:3) = CONFIG_RI(1:3,p) - CONFIG_RI(1:3,IOSEF1);            !
                                                                                         !
                        call APPLY_PBC(DRIJ(1:3),      &                                 !
                                      DRIJ(1:3),       &                                 !
                                      PASSA(1:3,1:3),  &                                 !
                                      PASSB(1:3,1:3));                                   !
                                                                                         !
                        RIJ = DSQRT( DOT_PRODUCT(DRIJ(1:3),DRIJ(1:3)) );                 !
                                                                                         !
                        if ( RIJ < RR2C2N ) then;                                        !
                                                                                         !
                            IOSEF2 = p;                                                  !
                                                                                         !
                            RR2C2N = RIJ;                                                !
                                                                                         !
                            RC2N_OLD(1:3) = CONFIG_RI(1:3,p);                            !
                                                                                         !
                        end if                                                           !
                                                                                         !
                    end do                                                               !
                                                                                         ! 
                    write(icanal,'(a11,i4)') '| IOSEF2 : ',      &                       !
                                             IOSEF2,             &                       !
                                             REPEAT(' ',54)//'|';                        !
                                                                                         !
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
                    write(icanal,'(a19,3f15.6,a6)') '| Set (old) RC2N : ', &             !
                                                    RC2N_OLD(1:3),         &             !
                                                    REPEAT(' ',5)//'|';                  !
                                                                                         ! 
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
!                   stop; !//////////////////////////////////////////////////////////////!
                                                                                         !
!                   ### Compute the vector between C2N and R2 ######################################
                                                                                         !
                    RC2NR2_OLD(1:3) = R2_REF_OLD(1:3) - RC2N_OLD(1:3);                   !
                                                                                         !
                    call APPLY_PBC(RC2NR2_OLD(1:3),  &                                   ! 
                                   RC2NR2_OLD(1:3),  &                                   !
                                   PASSA(1:3,1:3),   &                                   !
                                   PASSB(1:3,1:3));                                      !
                                                                                         !
                end if                                                                   !
                                                                                         !
!               if ( ninsertions == 1 ) stop; !//////////////////////////////////////////!
                                                                                         !
!               ### Delete the R2 site in the config array (swap coordinates) ######################
                                                                                         ! 
                if ( NRemaining_monomers > 1 .AND. IOSEF1 > 0 ) then;                    !
                                                                                         !
                    if ( IOSEF1 < NATOM_NEW ) then;                                      !
                                                                                         !
                        CONFIG_NAT(IOSEF1) = CONFIG_NAT(NATOM_NEW);                      ! 
                                                                                         !
                        CONFIG_RI(1:3,IOSEF1) = CONFIG_RI(1:3,NATOM_NEW);                !
                                                                                         !
                    end if                                                               !
                                                                                         !
                    NATOM_NEW = NATOM_NEW - 1;                                           !
                                                                                         !
                end if                                                                   !
                                                                                         !
!               ### Update the number of inserted monomers and the remaining monomers to insert ####
                                                                                         !
                ninsertions = ninsertions + 1;                                           !
                                                                                         !
                NRemaining_monomers = NRemaining_monomers - 1;                           !
                                                                                         !
                write(icanal,'(a22,i4,a44)') '| ninsertions (new) : ', &                 !
                                             ninsertions,              &                 !
                                             REPEAT(' ',43)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                write(icanal,'(a30,i4,a36)') '| NRemaining_monomers (new) : ', &         !
                                             NRemaining_monomers,              &         !
                                             REPEAT(' ',35)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          ! 
                                                                                         !
!               if ( ninsertions == 2 ) stop; !//////////////////////////////////////////!
                                                                                         !
!               ### Update local insertion probabilities ###########################################
                                                                                         !
                NMonomer_units_per_chain_local(iselect) = &                              !
                NMonomer_units_per_chain_local(iselect) - 1;                             !
                                                                                         !
                Insertion_probability_local(1:NSPECIES_PER_POLYMER_CHAINS(i)) = 0.0d0;   !
                                                                                         !
                if ( NRemaining_monomers == 0 ) CYCLE;                                   !
                                                                                         !
                do p = 1, NSPECIES_PER_POLYMER_CHAINS(i);                                !
                                                                                         !
                    ROSEF1 = REAL( NMonomer_units_per_chain_local(p) ) / &               !
                             REAL( NRemaining_monomers );                                !
                                                                                         !
                    if ( p > 1 ) Insertion_probability_local(p) = &                      !
                                 Insertion_probability_local(p-1);                       !
                                                                                         ! 
                    Insertion_probability_local(p) = Insertion_probability_local(p) + &  !
                                                     ROSEF1;                             !
                                                                                         !
                                                                                         ! 
                    write(icanal,'(a36,i4,a3,2f8.4,a11)') '| Insertion probability for monomer ', & !
                                                          p,                                      & !
                                                          ' : ',                                  & !
                                                          ROSEF1,                                 & !
                                                          Insertion_probability_local(p),         & !
                                                          REPEAT(' ',10)//'|';           !   
                                                                                         !
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
                end do                                                                   !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end do                                                                       !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end do                                                                           !
                                                                                         !
!       ### Deallocate arrays related to properties of polymer chains ##############################
                                                                                         !
        deallocate(NMonomer_units_per_chain_local);                                      !
                                                                                         !
!       ### Deallocate arrays for insertion probability ############################################
                                                                                         !
        deallocate(Insertion_probability);                                               !
                                                                                         !
        deallocate(Insertion_probability_local);                                         !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         ! 
!   ### Update simulation box properties ###########################################################
                                                                                         !
    NATOM = NATOM_NEW;                                                                   !
                                                                                         !
!   ### Initialization of arrays with the size of the total number of atoms ########################
                                                                                         !
    CONFIG_QI(1:NATOM) = 0.0d0;                                                          !
                                                                                         !
    CONFIG_VI(1:3,1:NATOM) = 0.0d0;                                                      !
                                                                                         !
    CONFIG_ATOM_TYPE(1:NATOM) = 0;                                                       !
                                                                                         !
    CONFIG_MOLECULEID(1:NATOM) = 1;                                                      !
                                                                                         !
    iuse_moleculeid = 1;                                                                 !
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine BUILD_POLYMER_CHAINS

