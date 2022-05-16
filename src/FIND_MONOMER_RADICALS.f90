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

subroutine FIND_MONOMER_RADICALS(icanal,PASSA,PASSB) 

!   ************************************************************************************************
!   **                                  READ THE INPUT FILE                                       **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal : Canal on which output data are written                                            **
!   **                                                                                            **
!   ** PASSA  :                                                                                   **
!   **                                                                                            **
!   ** PASSB  :                                                                                   **
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

!   integer (kind=4) :: iselect, NRemaining_monomers, ifile_lib, ninsertions;

    integer (kind=4) :: iloc_r1, iloc_r2;

    integer (kind=4) :: LOCAL_NFILES;

!   integer (kind=4) :: IFOUND;

!   integer (kind=4) :: IINSERTED, IREGIONS_INSERTION, IATOM, IOVERLAPPING;

!   integer (kind=4) :: IATOM_TYPE_TRANSLATE, IMOLECULEID;

!   integer (kind=4) :: IBOND_TYPE_TRANSLATE, IANGLE_TYPE_TRANSLATE, IIMPROPER_TYPE_TRANSLATE;

!   integer (kind=4) :: IDIHEDRAL_TYPE_TRANSLATE;

!   integer (kind=4) :: NMonomer_total;

!   integer (kind=4) :: ILOOP_MAX1, NLOOP_MAX1;

!   real (kind=8) :: grnd;

!   real (kind=8) :: mselect;

!   real (kind=8) :: SUM_MASSES, INV_SUM_MASSES, RIJ, RR2C2N, RR1C11;

!   real (kind=8) :: TAN_THETA1, TAN_THETA2, COS_THETA1, COS_THETA2, THETA1_RAD, THETA2_RAD;

!   real (kind=8) :: COS_PHI1, COS_PHI2, PHI1_RAD, PHI2_RAD;

!   real (kind=8) :: NORMU, NORMV;

!   real (kind=8) :: NORM_RC2NR2_OLD, NORM_R1C11_NEW;

!   real (kind=8) :: DELTA_THETA, DELTA_PHI, DELTA_PSI;

!   real (kind=8) :: NEW_THETA, NEW_PHI, NEW_PSI;

!   real (kind=8) :: COS_DELTA_THETA, SIN_DELTA_THETA, COS_DELTA_PHI, SIN_DELTA_PHI;

!   real (kind=8) :: COS_DELTA_PSI, SIN_DELTA_PSI;

!   real (kind=8) :: MAX_RAND_THETA, MAX_RAND_PHI, MAX_RAND_PSI, RCUT1, RCUT2;

!   real (kind=8) :: RAND_THETA, RAND_PHI, RAND_PSI;

!   real (kind=8), dimension(1:3) :: RINO, RI, RAND_ANGLES, RG, DRIJ;

!   real (kind=8), dimension(1:3) :: R1_REF, R2_REF, R2_REF_OLD, RC11, RC2N, RC2N_OLD; 

!   real (kind=8), dimension(1:3) :: RC2NR2_OLD, RR1C11_NEW;

!   real (kind=8), dimension(1:3,1:3) :: MATROT,           &
!                                        MATRIX_ROTATION1, &
!                                        MATRIX_ROTATION2, &
!                                        MATRIX_ROTATION3;

!   integer (kind=4), allocatable, dimension(:) :: NMonomer_units_per_chain_local;

!   real (kind=8), allocatable, dimension(:) :: Insertion_probability, Insertion_probability_local;

!   ************************************************************************************************

!   integer (kind=4) :: NATOM_NEW, NBOND_NEW, NANGLE_NEW, NDIHEDRAL_NEW, NIMPROPER_NEW;

!   integer (kind=4) :: NTYPE_ATOM_NEW,     &
!                       NTYPE_BOND_NEW,     &
!                       NTYPE_ANGLE_NEW,    &
!                       NTYPE_DIHEDRAL_NEW, &
!                       NTYPE_IMPROPER_NEW;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the routine #############################################################
                                                                                         !
    CHTITLE = 'Find monomer radicals';                                                   !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set the number of library files to consider ################################################
                                                                                         !
    LOCAL_NFILES = 0;                                                                    !
                                                                                         !
    do i = 1, NGenerate_polymer_species;                                                 !
                                                                                         !
        do j = 1, NSPECIES_PER_POLYMER_CHAINS(i);                                        !
                                                                                         !
            LOCAL_NFILES = LOCAL_NFILES + 1;                                             !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a17,i8,a45)') '| LOCAL_NFILES : ', &                                  !
                                 LOCAL_NFILES,        &                                  !
                                 REPEAT(' ',24)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !  
                                                                                         !
!   ### Initialization of arrays for the location of radicals ######################################
                                                                                         !
    ILOCR1_LIBRARY(1:LOCAL_NFILES) = 0;                                                  !
                                                                                         !
    ILOCR2_LIBRARY(1:LOCAL_NFILES) = 0;                                                  !
                                                                                         !
!   ### Find the location of the first radical R1 ##################################################
                                                                                         !
!   do i = 1, NTotal_monomers;                                                           !
    do i = 1, LOCAL_NFILES;                                                              !
                                                                                         !
        iloc_r1 = 0;                                                                     !
                                                                                         !
        do j = 1, NATOM_LIBRARY(i);                                                      !
                                                                                         !
            if ( CONFIG_NAT_LIBRARY(j,i) /= 'R1' ) CYCLE;                                !
                                                                                         !
!           iloc_r1 = j;                                                                 ! 
                                                                                         !                               
            ILOCR1_LIBRARY(i) = j;                                                       !
                                                                                         !
        end do                                                                           !
                                                                                         !
!       if ( iloc_r1 == 0 ) then;                                                        !
        if ( ILOCR1_LIBRARY(i) == 0 ) then;                                              !
                                                                                         !
            write(icanal,'(a70)') '| No radical R1 found '// &                           !
                                  'in the library file'//    &                           !
                                  REPEAT(' ',28)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '| The file is maybe not elligible '// &               !
                                          'for this routine'//REPEAT(' ',19)//'|';       !
                                                                                         !
            write(icanal,'(a14)') 'End of Program';                                      !
                                                                                         !
            close(icanal);                                                               !
                                                                                         !
            stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end if                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   ### Find the location of the second radical R2 #################################################
                                                                                         !
!   do i = 1, NTotal_monomers;                                                           !
    do i = 1, LOCAL_NFILES;                                                              !
                                                                                         !
!       iloc_r2 = 0;                                                                     !
                                                                                         !
        do j = 1, NATOM_LIBRARY(i);                                                      !
                                                                                         !
            if ( CONFIG_NAT_LIBRARY(j,i) /= 'R2' ) CYCLE;                                !
                                                                                         !
!           iloc_r2 = j;                                                                 !  
                                                                                         !
            ILOCR2_LIBRARY(i) = j;                                                       ! 
                                                                                         !
        end do                                                                           !
                                                                                         !
!       if ( iloc_r2 == 0 ) then;                                                        !
        if ( ILOCR2_LIBRARY(i) == 0 ) then;                                              !
                                                                                         !
            write(icanal,'(a70)') '| No radical R2 found in '// &                        !
                                  'the library file'//          &                        !
                                  REPEAT(' ',28)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '| The file is maybe not elligible '//   &             !   
                                  'for this routine'//REPEAT(' ',19)//'|';               !
                                                                                         !
            write(icanal,'(a14)') 'End of Program';                                      !
                                                                                         !
            close(icanal);                                                               !
                                                                                         !
            stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end if                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine FIND_MONOMER_RADICALS

