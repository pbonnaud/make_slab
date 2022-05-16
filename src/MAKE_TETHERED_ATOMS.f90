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

subroutine MAKE_TETHERED_ATOMS(icanal,ILOCAL_INSERTS,KLOCAL_MODIFIED) 

!   ************************************************************************************************
!   **                                      Make tethered atoms                                   **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal : CANAL ON WHICH OUTPUT DATA ARE WRITTEN                                            **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

    use module_physical_constants;

    use module_size_arrays;

    use module_library;

    use module_config;

    use module_simulation_box;

    use module_regions;

    use module_inserts;

    use module_lmp_input;

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    integer (kind=4), intent(in) :: ILOCAL_INSERTS, KLOCAL_MODIFIED;

!   ************************************************************************************************

    integer (kind=4) :: i, j, k, m;

    integer (kind=4) :: LOCAL_NTYPE_ATOM_MAX, LOCAL_NTYPE_ATOM_MIN, LOCAL_NNTYPE_ATOM;

    integer (kind=4) :: LOCAL_NTETHER;

    integer (kind=4) :: IATOM_GHOST, IBOND_GHOST;

!   integer (kind=4) :: IFOUND;

!   integer (kind=4) :: ILOCAL_TYPE;

!   integer (kind=4) :: IINSERTED, IREGIONS_INSERTION, IATOM, IOVERLAPPING;

!   integer (kind=4) :: IMOLECULEID;

!   integer (kind=4) :: IATOM_TYPE_TRANSLATE,     &
!                       IBOND_TYPE_TRANSLATE,     &
!                       IANGLE_TYPE_TRANSLATE,    &
!                       IDIHEDRAL_TYPE_TRANSLATE, &
!                       IIMPROPER_TYPE_TRANSLATE, &
!                       IPAIR_CROSS_TRANSLATE;

    real (kind=8) :: grnd;

    real (kind=8) :: LOCAL_EPSILONIJ, LOCAL_SIGMAIJ;

!   real (kind=8) :: SUM_MASSES, INV_SUM_MASSES, RIJ;

    real (kind=8), dimension(1:3) :: RINO;

!   ************************************************************************************************

    integer (kind=4) :: NATOM_NEW, NBOND_NEW, NANGLE_NEW, NDIHEDRAL_NEW, NIMPROPER_NEW;

    integer (kind=4) :: NTYPE_ATOM_NEW,        &
                        NTYPE_BOND_NEW,        &
                        NTYPE_ANGLE_NEW,       &
                        NTYPE_DIHEDRAL_NEW,    &
                        NTYPE_IMPROPER_NEW,    &
                        NPAIR_COEFF_CROSS_NEW;

    integer (kind=4) :: NLAMMPS_INSERTS_NEW;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the routine #############################################################
                                                                                         !
    CHTITLE = 'Make tethered atoms';                                                     !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write locations containing insert_modify properties ########################################
                                                                                         !
    write(icanal,'(a45,i4,a21)') '| The LAMMPS insert command to consider is : ', &      !
                                 ILOCAL_INSERTS,                                  &      !
                                 REPEAT(' ',20)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a41,i4,a25)') '| The modification rank to consider is : ', &          !
                                 KLOCAL_MODIFIED,                             &          !
                                 REPEAT(' ',24)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      ! 
                                                                                         !
    LOCAL_NTYPE_ATOM_MIN = 1;                                                            !
                                                                                         !
    LOCAL_NTYPE_ATOM_MAX = LAMMPS_INSERT_NTYPE_ATOM_MAX(ILOCAL_INSERTS);                 !
                                                                                         !
    if ( ILOCAL_INSERTS > 1 ) then;                                                      !
                                                                                         !
        LOCAL_NTYPE_ATOM_MIN = LAMMPS_INSERT_NTYPE_ATOM_MAX(ILOCAL_INSERTS-1) + 1;       !
                                                                                         !
    end if                                                                               !
                                                                                         !
    write(icanal,'(a26,2i4,a36)') '| Atom type to consider : ', &                        !
                                  LOCAL_NTYPE_ATOM_MIN,         &                        !
                                  LOCAL_NTYPE_ATOM_MAX,         &                        !
                                  REPEAT(' ',35)//'|';                                   !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      ! 
                                                                                         !
    LOCAL_NNTYPE_ATOM = LOCAL_NTYPE_ATOM_MAX - LOCAL_NTYPE_ATOM_MIN + 1;                 !
                                                                                         !
    write(icanal,'(a32,i4,a34)') '| Number of types to consider : ', &                   !
                                 LOCAL_NNTYPE_ATOM,                  &                   !
                                 REPEAT(' ',33)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      ! 
                                                                                         ! 
    LOCAL_NTETHER = 0;                                                                   !
                                                                                         !
    do j = LOCAL_NTYPE_ATOM_MIN, LOCAL_NTYPE_ATOM_MAX;                                   !
                                                                                         !
        do i = 1, NATOM;                                                                 ! 
                                                                                         ! 
            if ( CONFIG_ATOM_TYPE(i) /= j ) CYCLE;                                       !
                                                                                         !
            LOCAL_NTETHER = LOCAL_NTETHER + 1;                                           !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a30,i4,a36)') '| Number of atoms to tether : ', &                     !
                                 LOCAL_NTETHER,                    &                     !
                                 REPEAT(' ',35)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      ! 
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate local arrays ######################################################################
                                                                                         !
    call ALLOCATE_CONFIG_ARRAY_OLD(icanal,     &                                         !
                                   NATOM,      &                                         !
                                   NBOND,      &                                         !
                                   NANGLE,     &                                         !
                                   NDIHEDRAL,  &                                         !
                                   NIMPROPER);                                           !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    call ALLOCATE_PARAMETER_ARRAY_OLD(icanal,            &                               ! 
                                      NTYPE_ATOM,        &                               !
                                      NTYPE_BOND,        &                               !
                                      NTYPE_ANGLE,       &                               !
                                      NTYPE_DIHEDRAL,    &                               !
                                      NTYPE_IMPROPER,    &                               !
                                      NPAIR_COEFF_CROSS);                                !
                                                                                         ! 
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write output message for the allocation success ############################################
                                                                                         !
    write(icanal,'(a70)') '| Local arrays were allocated'//REPEAT(' ',40)//'|';          !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Update local arrays ########################################################################
                                                                                         !
    call UPDATE_CONFIG_ARRAY_OLD(icanal,     &                                           !
                                 NATOM,      &                                           !
                                 NBOND,      &                                           !
                                 NANGLE,     &                                           !
                                 NDIHEDRAL,  &                                           !
                                 NIMPROPER);                                             !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    call UPDATE_PARAMETER_ARRAY_OLD(icanal,            &                                 ! 
                                    NTYPE_ATOM,        &                                 !
                                    NTYPE_BOND,        &                                 !
                                    NTYPE_ANGLE,       &                                 !
                                    NTYPE_DIHEDRAL,    &                                 !
                                    NTYPE_IMPROPER,    &                                 !
                                    NPAIR_COEFF_CROSS);                                  !
                                                                                         ! 
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    write(icanal,'(a70)') '| Local arrays were updated'//REPEAT(' ',42)//'|';            !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Deallocate global arrays ###################################################################
                                                                                         !
    call DEALLOCATE_CONFIG_GLOBAL_ARRAYS(icanal,     &                                   !
                                         NATOM,      &                                   !
                                         NBOND,      &                                   !
                                         NANGLE,     &                                   !
                                         NDIHEDRAL,  &                                   !
                                         NIMPROPER);                                     !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    call DEALLOCATE_PARAMETER_GLOBAL_ARRAYS(icanal,            &                         !
                                            NTYPE_ATOM,        &                         !
                                            NTYPE_BOND,        &                         !
                                            NTYPE_ANGLE,       &                         !
                                            NTYPE_DIHEDRAL,    &                         !
                                            NTYPE_IMPROPER,    &                         !
                                            NPAIR_COEFF_CROSS);                          !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    write(icanal,'(a70)') '| Global arrays were deallocated'//REPEAT(' ',37)//'|';       !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialization of the size of new arrays ###################################################
                                                                                         !
    NATOM_NEW = NATOM;                                                                   !
                                                                                         !
    NBOND_NEW = NBOND;                                                                   !
                                                                                         !
    NANGLE_NEW = NANGLE;                                                                 !
                                                                                         !
    NDIHEDRAL_NEW = NDIHEDRAL;                                                           !
                                                                                         !
    NIMPROPER_NEW = NIMPROPER;                                                           !
                                                                                         !
    NTYPE_ATOM_NEW = NTYPE_ATOM;                                                         !
                                                                                         !
    NTYPE_BOND_NEW = NTYPE_BOND;                                                         !
                                                                                         !
    NTYPE_ANGLE_NEW = NTYPE_ANGLE;                                                       !
                                                                                         !
    NTYPE_DIHEDRAL_NEW = NTYPE_DIHEDRAL;                                                 !
                                                                                         !
    NTYPE_IMPROPER_NEW = NTYPE_IMPROPER;                                                 !
                                                                                         !
!   NPAIR_COEFF_CROSS_NEW = NPAIR_COEFF_CROSS;                                           !
    NPAIR_COEFF_CROSS_NEW = NPAIR_COEFF_CROSS + NTYPE_ATOM - LOCAL_NNTYPE_ATOM;          ! Add cross interaction with all other atom types
                                                                                         !
!   ### Modify the size of new arrays for tethered atoms ###########################################
                                                                                         !
    NATOM_NEW = NATOM_NEW + LOCAL_NTETHER;                                               ! Increase the size by the number of ghost atoms added to the system 
                                                                                         !
    NBOND_NEW = NBOND_NEW + LOCAL_NTETHER;                                               ! Increase the size by the number of bonds for tethering atoms
                                                                                         !
    NTYPE_ATOM_NEW = NTYPE_ATOM_NEW + 1;                                                 ! Add a new type of atom for ghost atoms
                                                                                         !
    NTYPE_BOND_NEW = NTYPE_BOND_NEW + 1;                                                 ! Add a new type of bond for tethering atoms 
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Re-allocate global arrays ##################################################################
                                                                                         !
    call ALLOCATE_CONFIG_GLOBAL_ARRAYS(icanal,         &                                 !
                                       NATOM_NEW,      &                                 !
                                       NBOND_NEW,      &                                 !
                                       NANGLE_NEW,     &                                 !
                                       NDIHEDRAL_NEW,  &                                 !
                                       NIMPROPER_NEW);                                   !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    call ALLOCATE_PARAMETER_GLOBAL_ARRAYS(icanal,                &                       !
                                          NTYPE_ATOM_NEW,        &                       !
                                          NTYPE_BOND_NEW,        &                       !
                                          NTYPE_ANGLE_NEW,       &                       !
                                          NTYPE_DIHEDRAL_NEW,    &                       !
                                          NTYPE_IMPROPER_NEW,    &                       !
                                          NPAIR_COEFF_CROSS_NEW);                        !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    write(icanal,'(a70)') '| Global arrays were re-allocated'//REPEAT(' ',36)//'|';      !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Update global arrays #######################################################################
                                                                                         !
    call UPDATE_CONFIG_ARRAY(icanal,     &                                               !
                             NATOM,      &                                               !
                             NBOND,      &                                               !
                             NANGLE,     &                                               !
                             NDIHEDRAL,  &                                               !
                             NIMPROPER);                                                 !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    call UPDATE_PARAMETER_ARRAY(icanal,            &                                     ! 
                                NTYPE_ATOM,        &                                     !
                                NTYPE_BOND,        &                                     !
                                NTYPE_ANGLE,       &                                     !
                                NTYPE_DIHEDRAL,    &                                     !
                                NTYPE_IMPROPER,    &                                     !
                                NPAIR_COEFF_CROSS);                                      !
                                                                                         ! 
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    write(icanal,'(a70)') '| Global arrays were updated'//REPEAT(' ',41)//'|';           !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !  
!   ### Deallocate local arrays ####################################################################
                                                                                         !
    call DEALLOCATE_CONFIG_ARRAY_OLD(icanal,     &                                       !
                                     NATOM,      &                                       !
                                     NBOND,      &                                       !
                                     NANGLE,     &                                       !
                                     NDIHEDRAL,  &                                       !
                                     NIMPROPER);                                         !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    call DEALLOCATE_PARAMETER_ARRAY_OLD(icanal,            &                             ! 
                                        NTYPE_ATOM,        &                             !
                                        NTYPE_BOND,        &                             !
                                        NTYPE_ANGLE,       &                             !
                                        NTYPE_DIHEDRAL,    &                             !
                                        NTYPE_IMPROPER,    &                             !
                                        NPAIR_COEFF_CROSS);                              !
                                                                                         ! 
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    write(icanal,'(a70)') '| Local arrays were deallocated'//REPEAT(' ',38)//'|';        !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Create ghost atoms #########################################################################
                                                                                         !
    IATOM_GHOST = NATOM;                                                                 !
                                                                                         !
    IBOND_GHOST = NBOND;                                                                 !
                                                                                         !
    do j = LOCAL_NTYPE_ATOM_MIN, LOCAL_NTYPE_ATOM_MAX;                                   !
                                                                                         !
        do i = 1, NATOM;                                                                 ! 
                                                                                         ! 
            if ( CONFIG_ATOM_TYPE(i) /= j ) CYCLE;                                       !

!           write(icanal,*) 'IATOM_GHOST : ', IATOM_GHOST
                                                                                         !
            IATOM_GHOST = IATOM_GHOST + 1;                                               !
                                                                                         !
            CONFIG_RI(1:3,IATOM_GHOST) = CONFIG_RI(1:3,i);                               !
                                                                                         !
            CONFIG_NAT(IATOM_GHOST) = 'Gh';                                              ! 
                                                                                         !
            CONFIG_QI(IATOM_GHOST) = 0.0d0;                                              ! Set the partial charge of the ghost atom
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
            CONFIG_VI(1:3,IATOM_GHOST) = 0.0d0;                                          !
                                                                                         !
            CONFIG_ATOMID(IATOM_GHOST) = IATOM_GHOST;                                    ! Set the ID of the ghost atom
                                                                                         ! 
            CONFIG_ATOM_TYPE(IATOM_GHOST) = NTYPE_ATOM_NEW;                              ! Set the type of the ghost atom
                                                                                         !
            CONFIG_MOLECULEID(IATOM_GHOST) = CONFIG_MOLECULEID(i);                       !
                                                                                         !
!           ### Add a random displacement to the real atom in order to avoid overlapping ###########
                                                                                         !
!           RINO(1:3) = (/grnd() - 0.5d0, grnd() - 0.5d0, grnd() - 0.5d0/);              !
                                                                                         !
!           CONFIG_RI(1:3,i) = CONFIG_RI(1:3,i) + RINO(1:3) * 0.1d0;                     !
!           CONFIG_RI(1:3,i) = CONFIG_RI(1:3,i) + RINO(1:3) * 0.1d0;                     !
            CONFIG_RI(1:3,i) = CONFIG_RI(1:3,i) + 0.1d0;                                 !
                                                                                         !
!           ### Update bond id arrays ##############################################################
                                                                                         !
            IBOND_GHOST = IBOND_GHOST + 1;                                               !
                                                                                         !
            BOND_ATOMID(1,IBOND_GHOST) = i;                                              !
                                                                                         !
            BOND_ATOMID(2,IBOND_GHOST) = IATOM_GHOST;                                    !
                                                                                         !
            BOND_TYPE(IBOND_GHOST) = NTYPE_BOND_NEW;                                     !
                                                                                         !
        end do                                                                           !
                                                                                         !
!       ### Add additional pair coeffs for cross interactions ######################################
                                                                                         !
        IOSEF1 = NPAIR_COEFF_CROSS;                                                      !
                                                                                         !
        do i = 1, NTYPE_ATOM;                                                            !
                                                                                         !
!           ### Check that the atom type i does not belong to atoms to be tethered #################
                                                                                         !
            IOSEF2 = 0;                                                                  !
                                                                                         !
            do k = LOCAL_NTYPE_ATOM_MIN, LOCAL_NTYPE_ATOM_MAX;                           !
                                                                                         !
                if ( i /= k ) CYCLE;                                                     !
                                                                                         !
                IOSEF2 = 1;                                                              !
                                                                                         !
            end do                                                                       !
                                                                                         !
            if ( IOSEF2 == 1 ) CYCLE;                                                    !
                                                                                         !
!           ### Compute cross interactions #########################################################
                                                                                         !
            LOCAL_EPSILONIJ = 0.0d0;                                                     !
                                                                                         !
            LOCAL_SIGMAIJ = 0.0d0;                                                       !
                                                                                         !
            if ( ( POTENTIAL_CLASS2(1,i) /= 0.0d0 ) .AND. &                              !
                 ( POTENTIAL_CLASS2(2,i) /= 0.0d0 ) ) then;                              !
                                                                                         !
                LOCAL_EPSILONIJ = POTENTIAL_CLASS2(1,j) * POTENTIAL_CLASS2(1,i);         !
                                                                                         !
                LOCAL_EPSILONIJ = DSQRT( LOCAL_EPSILONIJ );                              !
                                                                                         !
                LOCAL_SIGMAIJ = POTENTIAL_CLASS2(2,j) + POTENTIAL_CLASS2(2,i);           ! 
                                                                                         !
                LOCAL_SIGMAIJ = 0.5d0 * LOCAL_SIGMAIJ;                                   !
                                                                                         !
            end if                                                                       !
                                                                                         !
!           ### Set newly computed cross interactions ##############################################
                                                                                         !
            IOSEF1 = IOSEF1 + 1;                                                         ! 
                                                                                         !
            PAIR_ATOMID_CROSS(1,IOSEF1) = j;                                             !
                                                                                         !
            PAIR_ATOMID_CROSS(2,IOSEF1) = i;                                             !
                                                                                         ! 
            PAIR_COEFF_CROSS(1,IOSEF1) = LOCAL_EPSILONIJ;                                !
                                                                                         !
            PAIR_COEFF_CROSS(2,IOSEF1) = LOCAL_SIGMAIJ;                                  !
                                                                                         !
        end do                                                                           !
                                                                                         !
!       ### Cancel pair interaction between same atoms of the tethered group #######################
                                                                                         !
        POTENTIAL_CLASS2(1:2,j) = 0.0d0;                                                 !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '| Ghost atoms were created'//REPEAT(' ',43)//'|';             !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Update the size of arrays #################################################################
                                                                                         !
    NATOM = NATOM_NEW;                                                                   !
                                                                                         !
    NBOND = NBOND_NEW;                                                                   !
                                                                                         !
    NANGLE = NANGLE_NEW;                                                                 !
                                                                                         !
    NDIHEDRAL = NDIHEDRAL_NEW;                                                           !
                                                                                         !
    NIMPROPER = NIMPROPER_NEW;                                                           !
                                                                                         !
    NTYPE_ATOM = NTYPE_ATOM_NEW;                                                         !
                                                                                         !
    NTYPE_BOND = NTYPE_BOND_NEW;                                                         !
                                                                                         !
    NTYPE_ANGLE = NTYPE_ANGLE_NEW;                                                       !
                                                                                         !
    NTYPE_DIHEDRAL = NTYPE_DIHEDRAL_NEW;                                                 !
                                                                                         !
    NTYPE_IMPROPER = NTYPE_IMPROPER_NEW;                                                 !
                                                                                         !
    NPAIR_COEFF_CROSS = NPAIR_COEFF_CROSS_NEW;                                           !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Update atom labels and masses for the new molecular configuration #########################
                                                                                         !
    ATOM_LABEL(NTYPE_ATOM) = 'Gh';                                                       !
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !
    ATOM_MASSE(NTYPE_ATOM) = 0.1d0;                                                      !
                                                                                         !
    do j = LOCAL_NTYPE_ATOM_MIN, LOCAL_NTYPE_ATOM_MAX;                                   !
                                                                                         !
        if ( ATOM_MASSE(j) < ATOM_MASSE(NTYPE_ATOM) ) CYCLE;                             !
                                                                                         !
        ATOM_MASSE(NTYPE_ATOM) = ATOM_MASSE(j);                                          !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Update pair potential parameters for the new configuration #################################
                                                                                         !
    POTENTIAL_CLASS2(1:2,NTYPE_ATOM) = 0.0d0;                                            ! 
                                                                                         !
    write(icanal,'(a70)') '| Pair potential parameters were updated'// &                 !
                          REPEAT(' ',29)//'|';                                           !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Update bond types and bond coefficients ####################################################
                                                                                         !
    BOND_COEFFS(1:2,NTYPE_BOND) = LAMMPS_INSERT_TETHER_COEFFS(1:2,KLOCAL_MODIFIED);      !
                                                                                         !
    write(icanal,'(a70)') '| Bond types and bond coefficients were updated'// &          !
                          REPEAT(' ',22)//'|';                                           !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write properties of the molecular configuration ############################################
                                                                                         !
    write(icanal,'(a19,i8,a43)') '| NATOM          : ', &                                !
                                 NATOM,                 &                                !
                                 REPEAT(' ',42)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         ! 
    if ( NBOND > 0 ) then;                                                               !
                                                                                         ! 
        write(icanal,'(a19,i8,a43)') '| NBOND          : ', &                            !
                                     NBOND,                 &                            !
                                     REPEAT(' ',42)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NANGLE > 0 ) then;                                                              !
                                                                                         !
        write(icanal,'(a19,i8,a43)') '| NANGLE         : ', &                            !
                                     NANGLE,                &                            !
                                     REPEAT(' ',42)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NDIHEDRAL > 0 ) then;                                                           !
                                                                                         !
        write(icanal,'(a19,i8,a43)') '| NDIHEDRAL      : ', &                            !
                                     NDIHEDRAL,             &                            !
                                     REPEAT(' ',42)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NIMPROPER > 0 ) then;                                                           !
                                                                                         !
        write(icanal,'(a19,i8,a43)') '| NIMPROPER      : ', &                            !
                                     NIMPROPER,             &                            !
                                     REPEAT(' ',42)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_ATOM > 0 ) then;                                                          !
                                                                                         !
        write(icanal,'(a19,i8,a43)') '| NTYPE_ATOM     : ', &                            !
                                     NTYPE_ATOM,            &                            !
                                     REPEAT(' ',42)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_BOND > 0 ) then;                                                          !
                                                                                         !
        write(icanal,'(a19,i8,a43)') '| NTYPE_BOND     : ', &                            !
                                     NTYPE_BOND,            &                            !
                                     REPEAT(' ',42)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_ANGLE > 0 ) then;                                                         !
                                                                                         !
        write(icanal,'(a19,i8,a43)') '| NTYPE_ANGLE    : ', &                            !
                                     NTYPE_ANGLE,           &                            !
                                     REPEAT(' ',42)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_DIHEDRAL > 0 ) then;                                                      !
                                                                                         !
        write(icanal,'(a19,i8,a43)') '| NTYPE_DIHEDRAL : ', &                            !
                                     NTYPE_DIHEDRAL,        &                            !
                                     REPEAT(' ',42)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_IMPROPER > 0 ) then;                                                      !
                                                                                         !
        write(icanal,'(a19,i8,a43)') '| NTYPE_IMPROPER : ', &                            !
                                     NTYPE_IMPROPER,        &                            !
                                     REPEAT(' ',42)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Allocate local arrays for the insert command properties ####################################
                                                                                         !
    allocate(LAMMPS_INSERT_NAME_OLD(1:NLAMMPS_INSERTS));                                 !
                                                                                         !
    allocate(LAMMPS_INSERT_NTYPE_ATOM_MAX_OLD(1:NLAMMPS_INSERTS));                       !
                                                                                         !
!   ### Update local arrays for the insert command properties ######################################
                                                                                         !
    LAMMPS_INSERT_NAME_OLD(1:NLAMMPS_INSERTS) = LAMMPS_INSERT_NAME(1:NLAMMPS_INSERTS);   !
                                                                                         !
    LAMMPS_INSERT_NTYPE_ATOM_MAX_OLD(1:NLAMMPS_INSERTS) = &                              !
    LAMMPS_INSERT_NTYPE_ATOM_MAX(1:NLAMMPS_INSERTS);                                     !
                                                                                         !
!   ### Deallocate global arrays for the insert command properties #################################
                                                                                         !
    deallocate(LAMMPS_INSERT_NAME);                                                      !
                                                                                         !
    deallocate(LAMMPS_INSERT_NTYPE_ATOM_MAX);                                            !
                                                                                         !
!   ### Initialization of the size of new arrays ###################################################
                                                                                         !
    NLAMMPS_INSERTS_NEW = NLAMMPS_INSERTS;                                               !
                                                                                         !
!   ### Modify the size of new arrays for including the insertion of tethered atoms ################
                                                                                         !
    NLAMMPS_INSERTS_NEW = NLAMMPS_INSERTS_NEW + 1;                                       ! 
                                                                                         !
!   ### Re-allocate global arrays ##################################################################
                                                                                         !
    allocate(LAMMPS_INSERT_NAME(1:NLAMMPS_INSERTS_NEW));                                 !
                                                                                         !
    allocate(LAMMPS_INSERT_NTYPE_ATOM_MAX(1:NLAMMPS_INSERTS_NEW));                       !
                                                                                         !
!   ### Update global arrays #######################################################################
                                                                                         !
    LAMMPS_INSERT_NAME(1:NLAMMPS_INSERTS) = LAMMPS_INSERT_NAME_OLD(1:NLAMMPS_INSERTS);   !
                                                                                         !
    LAMMPS_INSERT_NTYPE_ATOM_MAX(1:NLAMMPS_INSERTS) =    &                               !
    LAMMPS_INSERT_NTYPE_ATOM_MAX_OLD(1:NLAMMPS_INSERTS);                                 !
                                                                                         !
!   ### Deallocate local arrays ####################################################################
                                                                                         !
    deallocate(LAMMPS_INSERT_NAME_OLD);                                                  !
                                                                                         !
    deallocate(LAMMPS_INSERT_NTYPE_ATOM_MAX_OLD);                                        !
                                                                                         !
!   ### Create insert command properties for ghost atoms ###########################################
                                                                                         !
    LAMMPS_INSERT_NAME(NLAMMPS_INSERTS_NEW) =           &                                !
    TRIM(LAMMPS_INSERT_NAME(ILOCAL_INSERTS))//'_ghost';                                  !
                                                                                         ! 
    LAMMPS_INSERT_NTYPE_ATOM_MAX(NLAMMPS_INSERTS_NEW) = NTYPE_ATOM;                      !
                                                                                         !
!   ### Update the size of arrays for the insert command properties for ghost atoms ################
                                                                                         !
    NLAMMPS_INSERTS = NLAMMPS_INSERTS_NEW;                                               !
                                                                                         !
!   ### Allocate local arrays for the lmp_input group command properties ###########################
                                                                                         !
    allocate(LMP_INPUT_GROUP_NAME_OLD(1:NLMP_INPUT_GROUP));                              !
                                                                                         ! 
    allocate(LMP_INPUT_GROUP_NLIST_OLD(1:NLMP_INPUT_GROUP));                             !
                                                                                         !
    allocate(LMP_INPUT_GROUP_LIST_OLD(1:NLAMMPS_INSERTS,1:NLMP_INPUT_GROUP));            !
                                                                                         !
!   ### Update local arrays for the lmp_input group command properties #############################
                                                                                         !
    LMP_INPUT_GROUP_NAME_OLD(1:NLMP_INPUT_GROUP) = &                                     !
    LMP_INPUT_GROUP_NAME(1:NLMP_INPUT_GROUP);                                            !
                                                                                         !
    LMP_INPUT_GROUP_NLIST_OLD(1:NLMP_INPUT_GROUP) = &                                    !
    LMP_INPUT_GROUP_NLIST(1:NLMP_INPUT_GROUP);                                           !
                                                                                         !
    LMP_INPUT_GROUP_LIST_OLD(1:NLAMMPS_INSERTS,1:NLMP_INPUT_GROUP) = &                   !
    LMP_INPUT_GROUP_LIST(1:NLAMMPS_INSERTS,1:NLMP_INPUT_GROUP);                          !
                                                                                         !
!   ### Deallocate glocal arrays ###################################################################
                                                                                         !
    deallocate(LMP_INPUT_GROUP_NAME);                                                    !
                                                                                         ! 
    deallocate(LMP_INPUT_GROUP_NLIST);                                                   !
                                                                                         !
    deallocate(LMP_INPUT_GROUP_LIST);                                                    !
                                                                                         !
!   ### Initialization of the size of new arrays ###################################################
                                                                                         !
    NLMP_INPUT_GROUP_NEW = NLMP_INPUT_GROUP;                                             !
                                                                                         !
!   ### Modify the size of new arrays for including the insertion of tethered atoms ################
                                                                                         !
    NLMP_INPUT_GROUP_NEW = NLMP_INPUT_GROUP_NEW + 1;                                     ! 
                                                                                         !
!   ### Re-allocate global arrays ##################################################################
                                                                                         !
    allocate(LMP_INPUT_GROUP_NAME(1:NLMP_INPUT_GROUP_NEW));                              !
                                                                                         !
    allocate(LMP_INPUT_GROUP_NLIST(1:NLMP_INPUT_GROUP_NEW));                             !
                                                                                         !
    allocate(LMP_INPUT_GROUP_LIST(1:NLAMMPS_INSERTS,1:NLMP_INPUT_GROUP_NEW));            !
                                                                                         !
!   ### Update global arrays #######################################################################
                                                                                         !
    LMP_INPUT_GROUP_NAME(1:NLMP_INPUT_GROUP) = &                                         !
    LMP_INPUT_GROUP_NAME_OLD(1:NLMP_INPUT_GROUP);                                        !
                                                                                         !
    LMP_INPUT_GROUP_NLIST(1:NLMP_INPUT_GROUP) = &                                        !
    LMP_INPUT_GROUP_NLIST_OLD(1:NLMP_INPUT_GROUP);                                       !
                                                                                         !
    LMP_INPUT_GROUP_LIST(1:NLAMMPS_INSERTS,1:NLMP_INPUT_GROUP) = &                       !
    LMP_INPUT_GROUP_LIST_OLD(1:NLAMMPS_INSERTS,1:NLMP_INPUT_GROUP);                      !
                                                                                         !
!   ### Deallocate local arrays ####################################################################
                                                                                         !
    deallocate(LMP_INPUT_GROUP_NAME_OLD);                                                !
                                                                                         ! 
    deallocate(LMP_INPUT_GROUP_NLIST_OLD);                                               !
                                                                                         !
    deallocate(LMP_INPUT_GROUP_LIST_OLD);                                                !
                                                                                         !
!   ### Create lmp_input command properties for ghost atoms ########################################
                                                                                         !
    IOSEF1 = 0;                                                                          !
                                                                                         !
    IOSEF2 = 0;                                                                          !
                                                                                         !
    do i = 1, NLMP_INPUT_GROUP;                                                          !
                                                                                         !
        do j = 1, LMP_INPUT_GROUP_NLIST(i);                                              !
                                                                                         !
            CHOSEF1 = TRIM(LMP_INPUT_GROUP_LIST(j,i));                                   !
                                                                                         !
            if ( TRIM(CHOSEF1) /= TRIM(LAMMPS_INSERT_NAME(ILOCAL_INSERTS)) ) CYCLE;      !
                                                                                         !
            IOSEF1 = i;                                                                  !
                                                                                         ! 
            IOSEF2 = j;                                                                  !
                                                                                         !
            EXIT;                                                                        !
                                                                                         !
        end do                                                                           !
                                                                                         !
        if ( IOSEF1 > 0 ) EXIT;                                                          ! 
                                                                                         !
    end do                                                                               !
                                                                                         !
    LMP_INPUT_GROUP_NAME(NLMP_INPUT_GROUP_NEW) = &                                       !
    TRIM(LMP_INPUT_GROUP_NAME(IOSEF1))//'_ghost';                                        !
                                                                                         !
    LMP_INPUT_GROUP_NLIST(NLMP_INPUT_GROUP_NEW) = &                                      !
    LMP_INPUT_GROUP_NLIST(IOSEF1);                                                       !
                                                                                         !
    do i = 1, LMP_INPUT_GROUP_NLIST(IOSEF1);                                             !
                                                                                         !
        LMP_INPUT_GROUP_LIST(i,NLMP_INPUT_GROUP_NEW) = &                                 !
        TRIM(LMP_INPUT_GROUP_LIST(i,IOSEF1));                                            !
                                                                                         !
        if ( i /= IOSEF2 ) CYCLE;                                                        !
                                                                                         !
        LMP_INPUT_GROUP_LIST(i,NLMP_INPUT_GROUP_NEW) = &                                 !
        TRIM(LMP_INPUT_GROUP_LIST(i,IOSEF1))//'_ghost';                                  !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   ### Update the size of arrays for the lmp_input command properties for ghost atoms #############
                                                                                         !
    NLMP_INPUT_GROUP = NLMP_INPUT_GROUP_NEW;                                             !
                                                                                         !
!   ### Update the group of frozen atoms ###########################################################
                                                                                         !
    if ( NLMP_INPUT_FROZEN > 0 ) then;                                                   !
                                                                                         !
       do i = 1, NLMP_INPUT_FROZEN;                                                      !
                                                                                         !
            if ( TRIM(LMP_INPUT_GROUP_NAME(IOSEF1)) /= &                                 !
                 TRIM(LMP_INPUT_FROZEN_GRPNAME(i))  ) CYCLE;                             !
                                                                                         ! 
            LMP_INPUT_FROZEN_GRPNAME(i) = TRIM(LMP_INPUT_FROZEN_GRPNAME(i))//'_ghost';   !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end do                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Allocate local arrays for the lmp_input group union command properties #####################
                                                                                         !
!   if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
!       allocate(LMP_INPUT_GROUP_UNION_NAME_OLD(1:NLMP_INPUT_GROUP_UNION));              !
                                                                                         !
!       allocate(LMP_INPUT_GROUP_UNION_NLIST_OLD(1:NLMP_INPUT_GROUP_UNION));             !
                                                                                         !
!       allocate(LMP_INPUT_GROUP_UNION_LIST_OLD(1:NLAMMPS_INSERTS,          &            !
!                                               1:NLMP_INPUT_GROUP_UNION));              !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Update local arrays for the definition of group unions in the lammps simulation ############
                                                                                         !
!   if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
!       LMP_INPUT_GROUP_UNION_NAME_OLD(1:NLMP_INPUT_GROUP_UNION) = &                     !
!       LMP_INPUT_GROUP_UNION_NAME(1:NLMP_INPUT_GROUP_UNION);                            !
                                                                                         !
!       LMP_INPUT_GROUP_UNION_NLIST_OLD(1:NLMP_INPUT_GROUP_UNION) = &                    !
!       LMP_INPUT_GROUP_UNION_NLIST(1:NLMP_INPUT_GROUP_UNION);                           !
                                                                                         !
!       LMP_INPUT_GROUP_UNION_LIST_OLD(1:NLAMMPS_INSERTS,1:NLMP_INPUT_GROUP_UNION) = &   !
!       LMP_INPUT_GROUP_UNION_LIST(1:NLAMMPS_INSERTS,1:NLMP_INPUT_GROUP_UNION);          !  
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Deallocate glocal arrays ###################################################################
                                                                                         !
!   if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
!       deallocate(LMP_INPUT_GROUP_UNION_NAME);                                          !
                                                                                         !
!       deallocate(LMP_INPUT_GROUP_UNION_NLIST);                                         !
                                                                                         !
!       deallocate(LMP_INPUT_GROUP_UNION_LIST);                                          !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Initialization of the size of new arrays ###################################################
                                                                                         !
!   if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
!       NLMP_INPUT_GROUP_UNION_NEW = NLMP_INPUT_GROUP_UNION;                             !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Modify the size of new arrays for including the insertion of tethered atoms ################
                                                                                         !
!   if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
!       NLMP_INPUT_GROUP_UNION_NEW = NLMP_INPUT_GROUP_UNION_NEW + 1;                     !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Re-allocate global arrays ##################################################################
                                                                                         !
!   if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
!       allocate(LMP_INPUT_GROUP_UNION_NAME(1:NLMP_INPUT_GROUP_UNION_NEW));              !
                                                                                         !
!       allocate(LMP_INPUT_GROUP_UNION_NLIST(1:NLMP_INPUT_GROUP_UNION_NEW));             !
                                                                                         !
!       allocate(LMP_INPUT_GROUP_UNION_LIST(1:NLAMMPS_INSERTS,              &            !
!                                           1:NLMP_INPUT_GROUP_UNION_NEW));              !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Update global arrays #######################################################################
                                                                                         !
!   if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
!       LMP_INPUT_GROUP_UNION_NAME(1:NLMP_INPUT_GROUP_UNION) =    &                      !
!       LMP_INPUT_GROUP_UNION_NAME_OLD(1:NLMP_INPUT_GROUP_UNION);                        !
                                                                                         !
!       LMP_INPUT_GROUP_UNION_NLIST(1:NLMP_INPUT_GROUP_UNION) =    &                     !
!       LMP_INPUT_GROUP_UNION_NLIST_OLD(1:NLMP_INPUT_GROUP_UNION);                       !
                                                                                         !
!       LMP_INPUT_GROUP_UNION_LIST(1:NLAMMPS_INSERTS,1:NLMP_INPUT_GROUP_UNION) =    &    !
!       LMP_INPUT_GROUP_UNION_LIST_OLD(1:NLAMMPS_INSERTS,1:NLMP_INPUT_GROUP_UNION);      !  
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Deallocate local arrays ####################################################################
                                                                                         !
!   if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
!       deallocate(LMP_INPUT_GROUP_UNION_NAME_OLD);                                      !
                                                                                         !
!       deallocate(LMP_INPUT_GROUP_UNION_NLIST_OLD);                                     !
                                                                                         !
!       deallocate(LMP_INPUT_GROUP_UNION_LIST_OLD);                                      !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Create lmp_input command properties for ghost atoms ########################################
                                                                                         !
!   if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
!       IOSEF1 = 0;                                                                      !
                                                                                         !
!       IOSEF2 = 0;                                                                      !
                                                                                         !
!       IOSEF3 = 0;                                                                      !
                                                                                         !
!       do i = 1, NLMP_INPUT_GROUP_UNION;                                                !
                                                                                         !
!           do j = 1, LMP_INPUT_GROUP_UNION_NLIST(i);                                    !
                                                                                         !
!               CHOSEF1 = TRIM(LMP_INPUT_GROUP_UNION_LIST(j,i));                         !
                                                                                         !
!               CHOSEF2 = TRIM(LMP_INPUT_GROUP_NAME(NLMP_INPUT_GROUP));                  !
                                                                                         !
!               if ( TRIM(CHOSEF1) /= TRIM(CHOSEF2) ) CYCLE;                             !
                                                                                         !
!               do k = 1, NLMP_INPUT_GROUP;                                              !
                                                                                         !                       
!                   do m = 1, LMP_INPUT_GROUP_NLIST(k);                                  !
                                                                                         !
!                       CHOSEF2 = TRIM(LMP_INPUT_GROUP_LIST(m,k));                       ! 
                                                 !
!                       if ( TRIM(CHOSEF1) /= TRIM(CHOSEF2) ) CYCLE;


!                   end do              !

!               end do

!               CHOSEF1 = TRIM(  )   //'_ghost';                                       !
                                                                                         !
!               if ( TRIM(CHOSEF1) == 

!               IOSEF1 = i;                                                              !
                                                                                         ! 
!               IOSEF2 = j;                                                              !
                                                                                         !
!               EXIT;                                                                    !
                                                                                         !
!           end do                                                                       !
                                                                                         !
!           if ( IOSEF1 > 0 ) EXIT;                                                      ! 
                                                                                         !
!       end do                                                                           !
                                                                                         !
!       LMP_INPUT_GROUP_UNION_NAME(NLMP_INPUT_GROUP_UNION_NEW) = &                       !
!       TRIM(LMP_INPUT_GROUP_UNION_NAME(IOSEF1))//'_ghost';                              !
                                                                                         !
!       LMP_INPUT_GROUP_UNION_NLIST(NLMP_INPUT_GROUP_UNION_NEW) = &                      !
!       LMP_INPUT_GROUP_UNION_NLIST(IOSEF1);                                             !
                                                                                         !
!       do i = 1, LMP_INPUT_GROUP_UNION_NLIST(IOSEF1);                                   !
                                                                                         !
!           LMP_INPUT_GROUP_UNION_LIST(i,NLMP_INPUT_GROUP_UNION_NEW) = &                 !
!           TRIM(LMP_INPUT_GROUP_UNION_LIST(i,IOSEF1));                                  !
                                                                                         !
!           if ( i /= IOSEF2 ) CYCLE;                                                    !
                                                                                         !
!           LMP_INPUT_GROUP_UNION_LIST(i,NLMP_INPUT_GROUP_UNION_NEW) = &                 !
!           TRIM(LMP_INPUT_GROUP_UNION_LIST(i,IOSEF1))//'_ghost';                        !
                                                                                         !
!       end do                                                                           !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Update the size of arrays for the lmp_input command properties for ghost atoms #############
                                                                                         !
!   if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
!       NLMP_INPUT_GROUP_UNION = NLMP_INPUT_GROUP_UNION_NEW;                             !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine MAKE_TETHERED_ATOMS

