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

subroutine INSERT_LIBRARY_MOLECULE(icanal,ILOCAL_INSERTS) 

!   ************************************************************************************************
!   **                                  READ THE INPUT FILE                                       **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal : CANAL ON WHICH OUTPUT DATA ARE WRITTEN                                            **
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

    use module_simulation_box;

    use module_regions;

    use module_inserts;

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    integer (kind=4), intent(in) :: ILOCAL_INSERTS;

!   ************************************************************************************************

    integer (kind=4) :: i, j, k, m;

    integer (kind=4) :: IFOUND;

    integer (kind=4) :: ILOCAL_TYPE;

    integer (kind=4) :: IINSERTED, IREGIONS_INSERTION, IATOM, IOVERLAPPING;

    integer (kind=4) :: IMOLECULEID;

    integer (kind=4) :: IATOM_TYPE_TRANSLATE,     &
                        IBOND_TYPE_TRANSLATE,     &
                        IANGLE_TYPE_TRANSLATE,    &
                        IDIHEDRAL_TYPE_TRANSLATE, &
                        IIMPROPER_TYPE_TRANSLATE, &
                        IPAIR_CROSS_TRANSLATE;

    real (kind=8) :: grnd;

    real (kind=8) :: SUM_MASSES, INV_SUM_MASSES, RIJ;

    real (kind=8), dimension(1:3) :: RIRAND, RINO, RI, RI_REF, RAND_ANGLES, RG, DRIJ;

    real (kind=8), dimension(1:3,1:3) :: MATROT;

!   ************************************************************************************************

    integer (kind=4) :: NATOM_NEW, NBOND_NEW, NANGLE_NEW, NDIHEDRAL_NEW, NIMPROPER_NEW;

    integer (kind=4) :: NTYPE_ATOM_NEW,        &
                        NTYPE_BOND_NEW,        &
                        NTYPE_ANGLE_NEW,       &
                        NTYPE_DIHEDRAL_NEW,    &
                        NTYPE_IMPROPER_NEW,    &
                        NPAIR_COEFF_CROSS_NEW;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the routine #############################################################
                                                                                         !
    CHTITLE = 'Insert the library molecule';                                             !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialization of LAMMPS configuration properties ##########################################
                                                                                         !
    if ( ILOCAL_INSERTS == 1 ) then;                                                     !
                                                                                         !
        NATOM     = 0;                                                                   !
        NBOND     = 0;                                                                   !
        NANGLE    = 0;                                                                   !
        NDIHEDRAL = 0;                                                                   !
        NIMPROPER = 0;                                                                   !
                                                                                         !
        NTYPE_ATOM     = 0;                                                              !
        NTYPE_BOND     = 0;                                                              !
        NTYPE_ANGLE    = 0;                                                              !
        NTYPE_DIHEDRAL = 0;                                                              !
        NTYPE_IMPROPER = 0;                                                              !
                                                                                         !
        NPAIR_COEFF_CROSS = 0;                                                           !
                                                                                         !
    end if                                                                               !
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
    NATOM_NEW     = NATOM;                                                               !
                                                                                         !
    NBOND_NEW     = NBOND;                                                               !
                                                                                         !
    NANGLE_NEW    = NANGLE;                                                              !
                                                                                         !
    NDIHEDRAL_NEW = NDIHEDRAL;                                                           !
                                                                                         !
    NIMPROPER_NEW = NIMPROPER;                                                           !
                                                                                         !
    NTYPE_ATOM_NEW     = NTYPE_ATOM;                                                     !
                                                                                         !
    NTYPE_BOND_NEW     = NTYPE_BOND;                                                     !
                                                                                         !
    NTYPE_ANGLE_NEW    = NTYPE_ANGLE;                                                    !
                                                                                         !
    NTYPE_DIHEDRAL_NEW = NTYPE_DIHEDRAL;                                                 !
                                                                                         !
    NTYPE_IMPROPER_NEW = NTYPE_IMPROPER;                                                 !
                                                                                         !
    NPAIR_COEFF_CROSS_NEW = NPAIR_COEFF_CROSS;                                           !
                                                                                         !
!   ### Set the size of new arrays #################################################################
                                                                                         !
    call UPDATE_SIZE_OF_ARRAYS(icanal,                &                                  !
                               NATOM_NEW,             &                                  !
                               NBOND_NEW,             &                                  !
                               NANGLE_NEW,            &                                  !
                               NDIHEDRAL_NEW,         &                                  !
                               NIMPROPER_NEW,         &                                  !
                               NTYPE_ATOM_NEW,        &                                  !
                               NTYPE_BOND_NEW,        &                                  !
                               NTYPE_ANGLE_NEW,       &                                  !
                               NTYPE_DIHEDRAL_NEW,    &                                  !
                               NTYPE_IMPROPER_NEW,    &                                  !
                               NPAIR_COEFF_CROSS_NEW);                                   !
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
!   ### Initialization of variables ################################################################
                                                                                         !
    IPOTENTIAL_CLASS2 = 0;                                                               !
                                                                                         !
    IATOM_TYPE_TRANSLATE     = 0;                                                        !
                                                                                         !
    IBOND_TYPE_TRANSLATE     = 0;                                                        !
                                                                                         !
    IANGLE_TYPE_TRANSLATE    = 0;                                                        !
                                                                                         !
    IDIHEDRAL_TYPE_TRANSLATE = 0;                                                        !
                                                                                         !
    IIMPROPER_TYPE_TRANSLATE = 0;                                                        !
                                                                                         !
    IPAIR_CROSS_TRANSLATE = 0;                                                           !
                                                                                         !
!   if ( iworking_file_insertion == 0 ) then;                                            !
                                                                                         !
!       NATOM     = 0;                                                                   !
                                                                                         !
!       NBOND     = 0;                                                                   !
                                                                                         !
!       NANGLE    = 0;                                                                   !
                                                                                         !
!       NDIHEDRAL = 0;                                                                   !
                                                                                         !
!       NIMPROPER = 0;                                                                   ! 
                                                                                         !
!       NTYPE_ATOM     = 0;                                                              !
                                                                                         !
!       NTYPE_BOND     = 0;                                                              !
                                                                                         !
!       NTYPE_ANGLE    = 0;                                                              !
                                                                                         !
!       NTYPE_DIHEDRAL = 0;                                                              !
                                                                                         !
!       NTYPE_IMPROPER = 0;                                                              !
                                                                                         !
!       NPAIR_COEFF_CROSS = 0;                                                           !
                                                                                         !
!   else                                                                                 !
                                                                                         !
        IATOM_TYPE_TRANSLATE     = NTYPE_ATOM;                                           !
                                                                                         !
        IBOND_TYPE_TRANSLATE     = NTYPE_BOND;                                           !
                                                                                         !
        IANGLE_TYPE_TRANSLATE    = NTYPE_ANGLE;                                          !
                                                                                         !
        IDIHEDRAL_TYPE_TRANSLATE = NTYPE_DIHEDRAL;                                       !
                                                                                         !
        IIMPROPER_TYPE_TRANSLATE = NTYPE_IMPROPER;                                       ! 
                                                                                         !
        IPAIR_CROSS_TRANSLATE    = NPAIR_COEFF_CROSS;                                    !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Set the initial molecule ID for insertion ##################################################
                                                                                         !
    IMOLECULEID = 0;                                                                     !
                                                                                         !
!   if ( iworking_file_insertion == 1 ) then;                                            !
                                                                                         !
        do i = 1, NATOM;                                                                 !
                                                                                         !
            if ( CONFIG_MOLECULEID(i) > IMOLECULEID ) then;                              !
                                                                                         !
                IMOLECULEID = CONFIG_MOLECULEID(i);                                      !
                                                                                         !
            end if                                                                       !
                                                                                         !
        end do                                                                           !
                                                                                         !
!   end if                                                                               !
                                                                                         ! 
    write(icanal,'(a16,i8,a46)') '| IMOLECULEID : ', &                                   !
                                 IMOLECULEID,        &                                   !
                                 REPEAT(' ',45)//'|';                                    ! 
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialization of the number of parameters #################################################
                                                                                         !
    IOSEF1 = MAXVAL(NPARAM_BONDS_LIBRARY(1:NFILE_LIBRARY));                              !
                                                                                         !
    IOSEF2 = MAXVAL(NPARAM_ANGLES_LIBRARY(1:NFILE_LIBRARY));                             !
                                                                                         !
    IOSEF3 = MAXVAL(NPARAM_DIHEDRALS_LIBRARY(1:NFILE_LIBRARY));                          !
                                                                                         !
    if ( iworking_file_insertion == 1 ) then;                                            !
                                                                                         !
        if ( IOSEF1 < NPARAM_BONDS ) IOSEF1 = NPARAM_BONDS;                              !
                                                                                         !
        if ( IOSEF2 < NPARAM_ANGLES ) IOSEF2 = NPARAM_ANGLES;                            !
                                                                                         !
        if ( IOSEF3 < NPARAM_DIHEDRALS ) IOSEF3 = NPARAM_DIHEDRALS;                      !
                                                                                         !
    end if                                                                               !
                                                                                         !
    NPARAM_BONDS = IOSEF1;                                                               !
                                                                                         !
    NPARAM_ANGLES = IOSEF2;                                                              !
                                                                                         !
    NPARAM_DIHEDRALS = IOSEF3;                                                           !
                                                                                         !
!   ### Write the number of parameters expected ####################################################
                                                                                         !
    write(icanal,'(a17,i8,a45)') '| NPARAM_BONDS : ', &                                  !
                                 NPARAM_BONDS,        &                                  !
                                 REPEAT(' ',44)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a18,i8,a44)') '| NPARAM_ANGLES : ', &                                 !
                                 NPARAM_ANGLES,        &                                 !
                                 REPEAT(' ',43)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a21,i8,a41)') '| NPARAM_DIHEDRALS : ', &                              !
                                 NPARAM_DIHEDRALS,        &                              !
                                 REPEAT(' ',40)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Insertion of molecules read in the library #################################################
                                                                                         !
    do i = 1, NFILE_LIBRARY;                                                             ! Loop over the molecule read in the library
                                                                                         !
!       ### Write properties of the current library file ###########################################
                                                                                         !
        write(icanal,'(a17,i4,a49)') '| LIBRARY FILE # ', i, ' '//REPEAT('/',47)//'|';   !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        write(icanal,'(a23,i8,a39)') '| NATOM_LIBRARY      : ', &                        !
                                     NATOM_LIBRARY(i),          &                        !
                                     REPEAT(' ',38)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        write(icanal,'(a23,i8,a39)') '| NTYPE_ATOM_LIBRARY : ', &                        !
                                     NTYPE_ATOM_LIBRARY(i),     &                        !
                                     REPEAT(' ',38)//'|';                                ! 
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!       ### Initialization of the mass and coordianates of the center of mass ######################
                                                                                         !
        SUM_MASSES = 0.0d0;                                                              !
                                                                                         !
        RG(1:3) = 0.0d0;                                                                 !
                                                                                         !
!       ### Compute mass and center of mass of the current file ####################################
                                                                                         !
        do j = 1, NATOM_LIBRARY(i);                                                      ! LOOP OVER THE NUMBER OF ATOMS IN THE MOLECULE TYPE I
                                                                                         !
            ILOCAL_TYPE = CONFIG_ATOM_TYPE_LIBRARY(j,i);                                 !
                                                                                         !
!           write(icanal,*) ILOCAL_TYPE, ATOM_MASSES_LIBRARY(ILOCAL_TYPE,i);             !
                                                                                         !
            RG(1:3) = RG(1:3) + ATOM_MASSES_LIBRARY(ILOCAL_TYPE,i) * &                   ! Here, for the computation of the center of mass, we do not take into account periodic
                                CONFIG_RI_LIBRARY(1:3,j,i);                              ! bundary conditions. We assume that library molecules are not cut in the space
                                                                                         !
            SUM_MASSES = SUM_MASSES + ATOM_MASSES_LIBRARY(ILOCAL_TYPE,i);                !
                                                                                         !
!           do k = 1, NTYPE_ATOM_LIBRARY(i);                                             ! LOOP OVER THE NUMBER OF ATOM TYPES IN THE MOLECULE I
                                                                                         !
!               if ( TRIM(CONFIG_NAT_LIBRARY(j,i)) == &                                  !
!                    TRIM(ATOM_LABEL_LIBRARY(k,i)) ) then;                               !
                                                                                         !
!                   RG(1:3) = RG(1:3) + ATOM_MASSES_LIBRARY(k,i) * &                     ! Here, for the computation of the center of mass, we do not take into account periodic
!                                       CONFIG_RI_LIBRARY(1:3,j,i);                      ! bundary conditions. We assume that library molecules are not cut in the space
                                                                                         !
!                   SUM_MASSES = SUM_MASSES + ATOM_MASSES_LIBRARY(k,i);                  !
                                                                                         !
!                   write(icanal,*) j, ATOM_MASSES_LIBRARY(k,i);                         !
                                                                                         !
!                   EXIT;                                                                !
                                                                                         !
!               end if                                                                   !
                                                                                         !
!           end do                                                                       !
                                                                                         !
        end do                                                                           !
                                                                                         !
        INV_SUM_MASSES = 1.0d0 / SUM_MASSES;                                             !
                                                                                         !
!       ### Write the computed mass and coordinates of the center of mass ##########################
                                                                                         !
        write(icanal,'(a15,2f15.6,a25)') '| SUM MASSES : ', &                            !
                                         SUM_MASSES,        &                            !
                                         INV_SUM_MASSES,    &                            !
                                         REPEAT(' ',24)//'|';                            !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        RG(1:3) = RG(1:3) * INV_SUM_MASSES;                                              !
                                                                                         !
        write(icanal,'(a11,3f12.6,a23)') '| RG [A] : ',      &                           !
                                         RG(1:3),            &                           !
                                         REPEAT(' ',22)//'|';                            !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!       ### Create coordinates of atoms ############################################################
                                                                                         !
        IINSERTED = 0;                                                                   !
                                                                                         !
        do while ( IINSERTED < NMOLECULE_INSERTION(i) )                                  ! LOOP AS LONG AS THE NUMBER OF INSERTED MOLECULES IS LOWER THAN THE TRAGETED NUMBER
                                                                                         !
            IREGIONS_INSERTION = 0;                                                      !
                                                                                         !
            RINO(1:3) = 0.0d0;                                                           !
                                                                                         !
            RI(1:3) = 0.0d0;                                                             !
                                                                                         !
!           ### Determination of the center of mass coordinates for the newly created molecule #####
                                                                                         !
            select case(INSERTION_METHOD_LIBRARY(i));                                    !
                                                                                         !
            case(1);                                                                     !
                                                                                         !
                RINO(1:3) = (/grnd() - 0.5d0, grnd() - 0.5d0, grnd() - 0.5d0/);          ! SET RANDOMLY THE COORDINATE OF THE CENTER OF MASS OF THE MOLECULE TO INSERT
                                                                                         !
                RI(1:3) = MATMUL(LAMMPS_PASSA(1:3,1:3),RINO(1:3));                       ! DEFINE COORDINATES IN THE CARTESIAN AXIS
                                                                                         !
            case(2);                                                                     !
                                                                                         !
                IOSEF1 = IINSERTED + 1;                                                  !
                                                                                         !
                RI(1:3) = MOLECULE_COM_LIBRARY(1:3,IOSEF1,i);                            !
                                                                                         !
            case(3);                                                                     !
                                                                                         !
!               RIRAND(1:3) = (/grnd() - 0.5d0, grnd() - 0.5d0, grnd() - 0.5d0/);        !
                RIRAND(1:3) = (/grnd(), grnd(), grnd()/);                                !
                                                                                         !
                IOSEF1 = LIBRARY_REGION_IRANK(i);                                        !
                                                                                         !
                RI(1) = LAMMPS_REGION_ARGS(1,IOSEF1)       + &                           !
                        LAMMPS_REGION_DIMENSIONS(1,IOSEF1) * &                           !
                        RIRAND(1);                                                       !
                                                                                         ! 
                RI(2) = LAMMPS_REGION_ARGS(3,IOSEF1)       + &                           !
                        LAMMPS_REGION_DIMENSIONS(2,IOSEF1) * &                           !
                        RIRAND(2);                                                       ! 
                                                                                         !
                RI(3) = LAMMPS_REGION_ARGS(5,IOSEF1)       + &                           !
                        LAMMPS_REGION_DIMENSIONS(3,IOSEF1) * &                           !
                        RIRAND(3);                                                       !
                                                                                         !
!               write(icanal,*) RI(1:3);                                                 !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            case default;                                                                !
                                                                                         !
                write(icanal,'(a70)') '| Not implemented - stop '// &                    !
                                      REPEAT(' ',44)//'|';                               !
                                                                                         !
                call CLOSING_PROGRAM(icanal,1);                                          !
                                                                                         !
            end select                                                                   !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           ### Write coordinates of the newly created center of mass ##############################
                                                                                         !
            write(icanal,'(a24,3f12.6,a10)') '| Center location [A] : ', &               !
                                             RI(1:3),                    &               !
                                             REPEAT(' ',9)//'|';                         !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           ### Check the location of the new center for defined regions ###########################
                                                                                         !
!           if ( NREGIONS_INSERTION > 0 ) then;                                          ! IF REGIONS HAVE BEEN DEFINED FOR THE INSERTION OF MOLECULES, THEN
                                                                                         !
!               do j = 1, NREGIONS_INSERTION;                                            ! LOOP OVER THE REGIONS THAT WERE DEFINED
                                                                                         !
!                   if ( ( RI(REGION_AXIS(j)) > REGION_BOUNDS(1,j) ) .AND. &             ! TEST THE NEW GENERATED COORDINATES WITH RESPECT TO THE DEFINED REGIONS
!                        ( RI(REGION_AXIS(j)) < REGION_BOUNDS(2,j) ) ) then;             ! IF COORDINATES ARE IN THE DEFINED REGION, THEN
                                                                                         !
!                       IREGIONS_INSERTION = 1;                                          ! SET IREGIONS_INSERTION TO 1
                                                                                         !
!                   else                                                                 !
                                                                                         !
!                       IREGIONS_INSERTION = 0;                                          ! 
                                                                                         !
!                   end if                                                               !
                                                                                         !
!               end do                                                                   !
                                                                                         !
!           else                                                                         ! IF REGIONS WEER NOT DEFINED, THEN 
                                                                                         !
!               IREGIONS_INSERTION = 1;                                                  ! SET IREGIONS_INSERTION TO 1
                                                                                         !
!           end if                                                                       !
                                                                                         !
!           write(icanal,*) 'IREGIONS_INSERTION : ', IREGIONS_INSERTION;
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           if ( IREGIONS_INSERTION == 0 ) CYCLE;                                        ! CYCLE IF COORDINATES ARE NOT LOCATED IN THE DEFINED REGIONS
                                                                                         !
!           ### Check if there is overlapping with the created center of mass ######################
                                                                                         !
            IOVERLAPPING = 0;                                                            !
                                                                                         !
            do j = 1, NATOM;                                                             ! LOOP OVER THE NUMBER OF ATOMS IN THE MOLECULAR CONFIGURATION
                                                                                         !
                DRIJ(1:3) = CONFIG_RI(1:3,j) - RI(1:3);                                  ! 
                                                                                         !
                call APPLY_PBC(DRIJ(1:3),              &                                 !
                               DRIJ(1:3),              &                                 ! 
                               LAMMPS_PASSA(1:3,1:3),  &                                 !
                               LAMMPS_PASSB(1:3,1:3));                                   !
                                                                                         !
                RIJ = DSQRT( DOT_PRODUCT(DRIJ(1:3),DRIJ(1:3)) );                         !
                                                                                         !
                if ( RIJ <= 2.0d0 ) then;                                                !
                                                                                         !
                    IOVERLAPPING = IOVERLAPPING + 1;                                     !
                                                                                         !
                    EXIT;                                                                !
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a17,i4,a49)') '| IOVERLAPPING : ', &                          !
                                         IOVERLAPPING,        &                          !
                                         REPEAT(' ',48)//'|';                            !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            if ( IOVERLAPPING > 0 ) CYCLE;                                               !
                                                                                         ! 
            write(icanal,'(a70)') '| Overlapping (center)   [OK]'//REPEAT(' ',40)//'|';  !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           ### Rotate the molecule around its center of mass ###################################### 
                                                                                         !
            do k = 1, 1000;                                                              ! LOOP OVER THE NUMBER OF TRIES TO ROTATE THE MOLECULE FOR THE INSERTION
                                                                                         !
                if ( INSERTION_METHOD_LIBRARY(i) == 2 ) then;                            !
                                                                                         !
                    RAND_ANGLES(1:3) = 0.0d0;                                            !
                                                                                         !
                else                                                                     !
                                                                                         !
                    RAND_ANGLES(1:3) = (/ TWOPI * grnd(), &                              ! SET RANDOM ANGLES FOR THE ROTATION OF THE MOLECULE AROUND THE CENTER OF MASS 
                                          TWOPI * grnd(), &                              !
                                          TWOPI * grnd() /);                             !
                                                                                         !
                end if                                                                   !
                                                                                         !
                call SET_ROTATION_MATRIX(RAND_ANGLES(1:3), &                             !
                                         MATROT(1:3,1:3));                               !
                                                                                         !
                IATOM = NATOM;                                                           !
                                                                                         !
                do j = 1, NATOM_LIBRARY(i);                                              ! LOOP OVER THE NUMBER OF ATOMS IN THE MOLECULE TYPE I
                                                                                         !
                    IATOM = IATOM + 1;                                                   !
                                                                                         !
                    DRIJ(1:3) = CONFIG_RI_LIBRARY(1:3,j,i) - RG(1:3);                    !
                                                                                         !
                    call APPLY_PBC(DRIJ(1:3),             &                              !
                                   DRIJ(1:3),             &                              !
                                   LAMMPS_PASSA(1:3,1:3), &                              !
                                   LAMMPS_PASSB(1:3,1:3));                               !
                                                                                         !
                    DRIJ(1:3) = MATMUL(MATROT(1:3,1:3),DRIJ(1:3));                       !
                                                                                         !
                    CONFIG_RI(1:3,IATOM) = RI(1:3) + DRIJ(1:3);                          !
                                                                                         !
                    call APPLY_PBC(CONFIG_RI(1:3,IATOM),   &                             !
                                   CONFIG_RI(1:3,IATOM),   &                             ! 
                                   LAMMPS_PASSA(1:3,1:3),  &                             ! 
                                   LAMMPS_PASSB(1:3,1:3));                               !
                                                                                         !
                    CONFIG_NAT(IATOM) = CONFIG_NAT_LIBRARY(j,i);                         ! SET THE NATURE OF THE CONSIDERED ATOM
                                                                                         !
                    CONFIG_QI(IATOM) = CONFIG_QI_LIBRARY(j,i);                           ! SET THE PARTIAL CHARGE OF THE CONSIDERED ATOM
                                                                                         !
!                   stop; !//////////////////////////////////////////////////////////////!
                                                                                         !
                    if ( ( CONFIG_VI_LIBRARY(1,j,i) /= 0.0d0 ) .AND. &                   !
                         ( CONFIG_VI_LIBRARY(2,j,i) /= 0.0d0 ) .AND. &                   !
                         ( CONFIG_VI_LIBRARY(3,j,i) /= 0.0d0 ) ) then;                   !
                                                                                         !
                        CONFIG_VI(1:3,IATOM) = CONFIG_VI_LIBRARY(1:3,j,i);               !
                                                                                         !
                    else                                                                 !
                                                                                         !
                        CONFIG_VI(1:3,IATOM) = 0.0d0;                                    !
                                                                                         !
                    end if                                                               !
                                                                                         !
!                   write(icanal,*) CONFIG_VI(1:3,IATOM);                                !
                                                                                         !
!                   stop; !//////////////////////////////////////////////////////////////!
                                                                                         !
                    CONFIG_ATOMID(IATOM)    = NATOM + CONFIG_ATOMID_LIBRARY(j,i);        ! SET THE ID OF THE ATOM IN THE MOLECULAR CONFIGURATION
                                                                                         !
                    CONFIG_ATOM_TYPE(IATOM) = IATOM_TYPE_TRANSLATE + &                   ! SET THE TYPE OF THE CONSIDERED ATOM
                                              CONFIG_ATOM_TYPE_LIBRARY(j,i);             ! 
                                                                                         !
                end do                                                                   !
                                                                                         !
                if ( IMOLECULEID > 0 ) then;                                             !
                                                                                         !
                    call CHECK_OVERLAPPING_ATOMS(icanal,NATOM+1,               &         !
                                                        IATOM,                 &         !
                                                        LAMMPS_PASSA(1:3,1:3), &         !
                                                        LAMMPS_PASSB(1:3,1:3), &         !
                                                        1.5d0,                 &         !
                                                        IOVERLAPPING);                   !
                end if                                                                   !
                                                                                         !
                if ( IOVERLAPPING == 0 ) then;                                           !
                                                                                         !
                    do j = 1, NATOM_LIBRARY(i);                                          !
                                                                                         !
                        if ( ( CONFIG_VI(1,j) /= 0.0d0 ) .AND. &                         !
                             ( CONFIG_VI(2,j) /= 0.0d0 ) .AND. &                         !
                             ( CONFIG_VI(3,j) /= 0.0d0 ) ) then;                         !
                                                                                         !
                            CONFIG_VI(1:3,j) = MATMUL(MATROT(1:3,1:3),   &               !
                                                      CONFIG_VI(1:3,j));                 !
                                                                                         !
                        end if                                                           !
                                                                                         ! 
                    end do                                                               !
                                                                                         !
                    EXIT;                                                                !
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
!           ### Write output message for overlapping atoms #########################################
                                                                                         !
            write(icanal,'(a17,i4,a49)') '| IOVERLAPPING : ', &                          !
                                         IOVERLAPPING,        &                          !
                                         REPEAT(' ',48)//'|';                            !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            if ( IOVERLAPPING > 0 ) CYCLE;                                               ! IF THERE IS OVERLAP BETWEEN THE INSERTED MOLECULE AND THE MOLECULAR CONFIGURATION, CYCLE
                                                                                         !
            write(icanal,'(a70)') '| Overlapping (rotation) [OK]'//REPEAT(' ',40)//'|';  !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           ### Update the number of atoms and the molecule ID #####################################
                                                                                         !
            IMOLECULEID = IMOLECULEID + 1;                                               ! UPDATE THE MOLECULE ID
                                                                                         !
            IOSEF1 = NATOM + 1;                                                          !
                                                                                         !
            NATOM = IATOM;                                                               ! UPDATE THE NUMBER OF ATOMS IN THE FINAL CONFIGURATION
                                                                                         !
            if ( NATOM > NATOM_NEW ) then;                                               !
                                                                                         !
                write(icanal,'(a70)') '| NATOM is greater than the size '//      &       !
                                      'of allocated arrays'//REPEAT(' ',17)//'|';        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                write(icanal,'(a31,3i8,a15)') '| IOSEF1 | NATOM | NATOM_NEW : ', &       !
                                              IOSEF1, NATOM, NATOM_NEW,          &       !
                                              REPEAT(' ',14)//'|';                       !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                write(icanal,'(a70)') '| --> Stop'//REPEAT(' ',59)//'|';                 !
                                                                                         !
                call CLOSING_PROGRAM(icanal,1);                                          !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
            CONFIG_MOLECULEID(IOSEF1:NATOM) = IMOLECULEID;                               ! APPLY THE NEW MOLECULE ID 
                                                                                         !
            IINSERTED = IINSERTED + 1;                                                   !
                                                                                         !
            IOSEF1 = IOSEF1 - 1;                                                         !  
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           ### Update bond IDs and bond properties ################################################
                                                                                         !
            if ( NBOND_LIBRARY(i) > 0 ) then;                                            ! IF THE NUMBER OF BONDS OF THE INSERTED MOLECULE IS GREATER THAN 0, THEN
                                                                                         !
                IOSEF2 = NBOND + 1;                                                      !
                                                                                         !
                NBOND  = NBOND + NBOND_LIBRARY(i);                                       ! UPDATE THE NUMBER OF BONDS IN THE FINAL CONFIGURATION 
                                                                                         !
                if ( NBOND > NBOND_NEW ) then;                                           !
                                                                                         !
                    write(icanal,'(a70)') '| NBOND is greater than the size '//      &   !
                                          'of allocated arrays'//REPEAT(' ',17)//'|';    !
                                                                                         !
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
                    write(icanal,'(a32,3i8,a16)') '| IOSEF2 | NBOND | NBOND_NEW : ', &   !
                                                  IOSEF2, NBOND, NBOND_NEW,          &   !
                                                  REPEAT(' ',15)//'|';                   !
                                                                                         !
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
                    write(icanal,'(a70)') '| --> Stop'//REPEAT(' ',59)//'|';             !
                                                                                         !
                    call CLOSING_PROGRAM(icanal,1);                                      !
                                                                                         !
                end if                                                                   !
                                                                                         !
                BOND_ATOMID(1:2,IOSEF2:NBOND) = &                                        !
                IOSEF1                        + &                                        !
                BOND_ATOMID_LIBRARY(1:2,1:NBOND_LIBRARY(i),i);                           !
                                                                                         !
                BOND_TYPE(IOSEF2:NBOND) = &                                              !
                IBOND_TYPE_TRANSLATE    + &                                              !
                BOND_TYPE_LIBRARY(1:NBOND_LIBRARY(i),i);                                 !
                                                                                         !
            end if                                                                       !
                                                                                         !
!           write(icanal,*) 'OK BONDS';                                                  !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           ### Update angle IDs and angle properties ##############################################
                                                                                         !
            if ( NANGLE_LIBRARY(i) > 0 ) then;                                           !
                                                                                         !
                IOSEF3 = NANGLE + 1;                                                     !
                                                                                         !
                NANGLE = NANGLE + NANGLE_LIBRARY(i);                                     ! UPDATE THE NUMBER OF ANGLES IN THE FINAL CONFIGURATION
                                                                                         !
                if ( NANGLE > NANGLE_NEW ) then;                                         !
                                                                                         !
                    write(icanal,'(a70)') '| NANGLE is greater than the size '//      &  !
                                          'of allocated arrays'//REPEAT(' ',16)//'|';    !
                                                                                         !
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
                    write(icanal,'(a33,3i8,a13)') '| IOSEF2 | NANGLE | NANGLE_NEW : ', & !
                                                  IOSEF3, NANGLE, NANGLE_NEW,          & !
                                                  REPEAT(' ',12)//'|';                   !
                                                                                         !
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
                    write(icanal,'(a70)') '| --> Stop'//REPEAT(' ',59)//'|';             !
                                                                                         !
                    call CLOSING_PROGRAM(icanal,1);                                      !
                                                                                         !
                end if                                                                   !
                                                                                         !
                ANGLE_ATOMID(1:3,IOSEF3:NANGLE) = &                                      !
                IOSEF1                          + &                                      !
                ANGLE_ATOMID_LIBRARY(1:3,1:NANGLE_LIBRARY(i),i);                         !
                                                                                         ! 
                ANGLE_TYPE(IOSEF3:NANGLE) = &                                            !
                IANGLE_TYPE_TRANSLATE     + &                                            !
                ANGLE_TYPE_LIBRARY(1:NANGLE_LIBRARY(i),i);                               !
                                                                                         !
            end if                                                                       !
                                                                                         !
!           write(icanal,*) 'OK ANGLES';                                                 !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           ### Update dihedral IDs and dihedral properties ########################################
                                                                                         !
            if ( NDIHEDRAL_LIBRARY(i) > 0 ) then;                                        !
                                                                                         !
                IOSEF5 = NDIHEDRAL + 1;                                                  !
                                                                                         !
                NDIHEDRAL = NDIHEDRAL + NDIHEDRAL_LIBRARY(i);                            ! UPDATE THE NUMBER OF DIHEDRALS IN THE FINAL CONFIGURATION
                                                                                         !
                if ( NDIHEDRAL > NDIHEDRAL_NEW ) then;                                   !
                                                                                         !
                    write(icanal,'(a70)') '| NDIHEDRAL is greater than the '// &         !
                                          'size of allocated arrays'//         &         !
                                          REPEAT(' ',13)//'|';                           !
                                                                                         !
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
                    write(icanal,'(a39,3i8,a7)') '| IOSEF5 | NDIHEDRAL '//         &     !
                                                 '| NDIHEDRAL_NEW : ',             &     !
                                                 IOSEF5, NDIHEDRAL, NDIHEDRAL_NEW, &     !
                                                 REPEAT(' ',6)//'|';                     !
                                                                                         !
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
                    write(icanal,'(a70)') '| --> Stop'//REPEAT(' ',59)//'|';             !
                                                                                         !
                    call CLOSING_PROGRAM(icanal,1);                                      !
                                                                                         !
                end if                                                                   !
                                                                                         !
!               write(icanal,*) IOSEF5, NDIHEDRAL, NDIHEDRAL_NEW;                        !
                                                                                         !
                DIHEDRAL_ATOMID(1:4,IOSEF5:NDIHEDRAL) = &                                !
                IOSEF1                                + &                                !
                DIHEDRAL_ATOMID_LIBRARY(1:4,1:NDIHEDRAL_LIBRARY(i),i);                   !
                                                                                         !
                DIHEDRAL_TYPE(IOSEF5:NDIHEDRAL) = &                                      !
                IDIHEDRAL_TYPE_TRANSLATE        + &                                      !
                DIHEDRAL_TYPE_LIBRARY(1:NDIHEDRAL_LIBRARY(i),i);                         !
                                                                                         !
            end if                                                                       !
                                                                                         !
!           write(icanal,*) 'OK DIHEDRALS';                                              !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           ### Update improper IDs and improper properties ########################################
                                                                                         !
            if ( NIMPROPER_LIBRARY(i) > 0 ) then;                                        !
                                                                                         !
                IOSEF4 = NIMPROPER + 1;                                                  !
                                                                                         !
                NIMPROPER = NIMPROPER + NIMPROPER_LIBRARY(i);                            ! UPDATE THE NUMBER OF IMPROPERS IN THE FINAL CONFIGURATION
                                                                                         !
                if ( NIMPROPER > NIMPROPER_NEW ) then;                                   !
                                                                                         !
                    write(icanal,'(a70)') '| NIMPROPER is greater than the '// &         !
                                          'size of allocated arrays'//         &         !
                                          REPEAT(' ',13)//'|';                           !
                                                                                         !
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
                    write(icanal,'(a39,3i8,a7)') '| IOSEF4 | NIMPROPER '//         &     !
                                                 '| NIMPROPER_NEW : ',             &     !
                                                 IOSEF4, NIMPROPER, NIMPROPER_NEW, &     !
                                                 REPEAT(' ',6)//'|';                     !
                                                                                         !
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                      !
                                                                                         !
                    write(icanal,'(a70)') '| --> Stop'//REPEAT(' ',59)//'|';             !
                                                                                         !
                    call CLOSING_PROGRAM(icanal,1);                                      !
                                                                                         !
                end if                                                                   !
                                                                                         !
                IMPROPER_ATOMID(1:4,IOSEF4:NIMPROPER) = &                                !
                IOSEF1                                + &                                !
                IMPROPER_ATOMID_LIBRARY(1:4,1:NIMPROPER_LIBRARY(i),i);                   !
                                                                                         !
                IMPROPER_TYPE(IOSEF4:NIMPROPER) = &                                      !
                IIMPROPER_TYPE_TRANSLATE        + &                                      !
                IMPROPER_TYPE_LIBRARY(1:NIMPROPER_LIBRARY(i),i);                         !   
                                                                                         !
            end if                                                                       !
                                                                                         !
!           write(icanal,*) 'OK IMPROPERS';                                              !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           ### Update the current status of molecule insertion in the simulation box ##############
                                                                                         !
            ROSEF8 = 100.0d0 * REAL(IINSERTED) / REAL(NMOLECULE_INSERTION(i));           !
                                                                                         !
            write(icanal,'(a14,i8,f12.2,a36)') '| IINSERTED : ',    &                    !  
                                               IINSERTED,           &                    !
                                               ROSEF8,              &                    !
                                               REPEAT(' ',35)//'|';                      !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end do                                                                           !
                                                                                         !
!       write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       ### Update the number of atom types ########################################################
                                                                                         !
        IOSEF1 = NTYPE_ATOM + 1;                                                         ! UPDATE THE NUMBER OF ATOMS IN THE FINAL MOLECULAR CONFIGURATION
                                                                                         !
        NTYPE_ATOM = NTYPE_ATOM + NTYPE_ATOM_LIBRARY(i);                                 ! UPDATE THE NUMBER OF ATOM TYPES     
                                                                                         !
        if ( NTYPE_ATOM > NTYPE_ATOM_NEW ) then;                                         !
                                                                                         !
            write(icanal,'(a70)') '| NTYPE_ATOM is greater than the '// &                !
                                  'size of allocated arrays'//          &                !
                                  REPEAT(' ',12)//'|';                                   !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            write(icanal,'(a41,3i8,a5)') '| IOSEF1 | NTYPE_ATOM '// &                    !
                                         '| NTYPE_ATOM_NEW : ',     &                    !
                                         IOSEF1,                    &                    !
                                         NTYPE_ATOM,                &                    !
                                         NTYPE_ATOM_NEW;                                 !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            write(icanal,'(a70)') '| --> Stop'//REPEAT(' ',59)//'|';                     !
                                                                                         !
            call CLOSING_PROGRAM(icanal,1);                                              !
                                                                                         !
        end if                                                                           !
                                                                                         !
!       write(icanal,*) 'OK 2';                                                          !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!       ### Update the maximum number of ntype_atom for each insertion ############################# 
                                                                                         !
        LAMMPS_INSERT_NTYPE_ATOM_MAX(ILOCAL_INSERTS) = NTYPE_ATOM;                       !
                                                                                         !
!AAAAAAAAAAAAAAAAAAA
!AAAAAAAAAAAAAAAAAAA
!AAAAAAAAAAAAAAAAAAA
!AAAAAAAAAAAAAAAAAAA
!AAAAAAAAAAAAAAAAAAA
!AAAAAAAAAAAAAAAAAAA
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!       ### Update atom style ######################################################################
                                                                                         !
        CH_ATOM_STYLE = CH_ATOM_STYLE_LIBRARY(i);                                        !
                                                                                         !
!       ### Update atom labels and masses for the new molecular configuration ######################
                                                                                         !
        ATOM_LABEL(IOSEF1:NTYPE_ATOM) = ATOM_LABEL_LIBRARY(1:NTYPE_ATOM_LIBRARY(i),i);   ! UPDATE THE LIST OF LABEL OF ATOMS 
                                                                                         !
        ATOM_MASSE(IOSEF1:NTYPE_ATOM) = ATOM_MASSES_LIBRARY(1:NTYPE_ATOM_LIBRARY(i),i);  ! UPDATE THE LIST OF MASSES 
                                                                                         !
!       write(icanal,*) 'OK 3';                                                          !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!       ### Update pair potential parameters for the new configuration #############################
                                                                                         !
        if ( IPOTENTIAL_CLASS2_LIBRARY(i) == 1 ) then;                                   !
                                                                                         !
            POTENTIAL_CLASS2(1:2,IOSEF1:NTYPE_ATOM) = &                                  ! UPDATE LENNARD-JONES-POTENTIAL PARAMETERS
            POTENTIAL_CLASS2_LIBRARY(1:2,1:NTYPE_ATOM_LIBRARY(i),i);                     !
                                                                                         ! 
        end if                                                                           !
                                                                                         !
        CH_PAIR_STYLE = TRIM(POTENTIAL_CLASS2_CHTYPE_LIBRARY(i));                        !
                                                                                         !
        write(icanal,'(a70)') '| Pair potential parameters were updated'// &             !
                              REPEAT(' ',29)//'|';                                       !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!       ### Update pair coefficients for cross interactions ########################################
                                                                                         !
        if ( NPAIR_COEFF_CROSS_LIBRARY(i) > 0 ) then;                                    !
                                                                                         !
            IOSEF4 = NPAIR_COEFF_CROSS + 1;                                              !
                                                                                         !
            NPAIR_COEFF_CROSS = NPAIR_COEFF_CROSS + NPAIR_COEFF_CROSS_LIBRARY(i);        ! Update the number of pair coefficients for cross interactions
                                                                                         !
            if ( NPAIR_COEFF_CROSS > NPAIR_COEFF_CROSS_NEW ) then;                       !
                                                                                         !
                write(icanal,'(a70)') '| NPAIR_COEFF_CROSS is greater than '// &         !
                                      'the size of allocated arrays'//         &         !
                                      REPEAT(' ',5)//'|';                                !
                                                                                         ! 
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                write(icanal,'(a55,3i8)') '| IOSEF4 '//                 &                !
                                          '| NPAIR_COEFF_CROSS '//      &                !
                                          '| NPAIR_COEFF_CROSS_NEW : ', &                !
                                          IOSEF4,                       &                !
                                          NPAIR_COEFF_CROSS,            &                !
                                          NPAIR_COEFF_CROSS_NEW;                         !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                write(icanal,'(a70)') '| --> Stop'//REPEAT(' ',59)//'|';                 !
                                                                                         !
                call CLOSING_PROGRAM(icanal,1);                                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
!           ### Update cross interactions parameters for pair potentials ###########################
                                                                                         !
            PAIR_COEFF_CROSS(1:2,IOSEF4:NPAIR_COEFF_CROSS) = &                           !
            PAIR_COEFF_CROSS_LIBRARY(1:2,1:NPAIR_COEFF_CROSS_LIBRARY(i),i);              ! 
                                                                                         !
            PAIR_ATOMID_CROSS(1:2,IOSEF4:NPAIR_COEFF_CROSS) = &                          !
!           IPAIR_CROSS_TRANSLATE                           + &                          !
            IATOM_TYPE_TRANSLATE + &
            PAIR_ATOMID_CROSS_LIBRARY(1:2,1:NPAIR_COEFF_CROSS_LIBRARY(i),i);             !
                                                                                         !
            IPAIR_CROSS_TRANSLATE = &                                                    !
            IPAIR_CROSS_TRANSLATE + &                                                    !
            NTYPE_ATOM_LIBRARY(i);                                                       !
                                                                                         !

        end if                                                                           !
                                                                                         !
        write(icanal,'(a70)') '| Update number of pair interactions '// &                !
                              '(cross) and coefficients'//              &                !
                              REPEAT(' ',8)//'|';                                        !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!       ### Update the flag related to the presence of pair potentials #############################
                                                                                         !
        if ( IPOTENTIAL_CLASS2_LIBRARY(i) == 1 ) IPOTENTIAL_CLASS2 = 1;                  ! 
                                                                                         !
!       ### Update the variable for translating parameters in arrays ###############################
                                                                                         !
        IATOM_TYPE_TRANSLATE = IATOM_TYPE_TRANSLATE + NTYPE_ATOM_LIBRARY(i);             !
                                                                                         !
! QQQQQQQQQQQQQQQQQQQQQQQQQQQQ
! QQQQQQQQQQQQQQQQQQQQQQQQQQQQ
! QQQQQQQQQQQQQQQQQQQQQQQQQQQQ
! QQQQQQQQQQQQQQQQQQQQQQQQQQQQ
! QQQQQQQQQQQQQQQQQQQQQQQQQQQQ
! QQQQQQQQQQQQQQQQQQQQQQQQQQQQ


!       ### Update bond types and bond coefficients ################################################
                                                                                         !
        if ( NTYPE_BOND_LIBRARY(i) > 0 ) then;                                           !
                                                                                         !
            IOSEF2 = NTYPE_BOND + 1;                                                     !
                                                                                         !
            NTYPE_BOND = NTYPE_BOND + NTYPE_BOND_LIBRARY(i);                             ! UPDATE THE NUMBER OF BOND TYPES 
                                                                                         !
            if ( NTYPE_BOND > NTYPE_BOND_NEW ) then;                                     !
                                                                                         !
                write(icanal,'(a70)') '| NTYPE_BOND is greater than the '// &            !
                                      'size of allocated arrays'//          &            !
                                      REPEAT(' ',12)//'|';                               !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                write(icanal,'(a41,3i8,a5)') '| IOSEF2 | NTYPE_BOND '// &                !
                                             '| NTYPE_BOND_NEW : ',     &                !
                                             IOSEF2,                    &                !
                                             NTYPE_BOND,                &                !
                                             NTYPE_BOND_NEW;                             !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                write(icanal,'(a70)') '| --> Stop'//REPEAT(' ',59)//'|';                 !
                                                                                         !
                call CLOSING_PROGRAM(icanal,1);                                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
!           if ( NPARAM_BONDS_LIBRARY(i) > NPARAM_BONDS ) then;                          !
                                                                                         !
!               NPARAM_BONDS = NPARAM_BONDS_LIBRARY(i);                                  !
                                                                                         !
!           end if                                                                       !
                                                                                         !
            CH_BOND_STYLE = CH_BOND_STYLE_LIBRARY(i);                                    !
                                                                                         !
            IOSEF6 = NPARAM_BONDS_LIBRARY(i);                                            !
                                                                                         !
            BOND_COEFFS(1:IOSEF6,IOSEF2:NTYPE_BOND) = &                                  !
            BOND_COEFFS_LIBRARY(1:IOSEF6,1:NTYPE_BOND_LIBRARY(i),i);                     ! 
                                                                                         !
            IBOND_TYPE_TRANSLATE = IBOND_TYPE_TRANSLATE + NTYPE_BOND_LIBRARY(i);         !
                                                                                         !
        end if                                                                           !
                                                                                         ! 
        write(icanal,'(a70)') '| Bond types and bond coefficients were updated'// &      !
                              REPEAT(' ',22)//'|';                                       !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!       ### Update bond types and bond coefficients ################################################
                                                                                         !
        if ( NTYPE_ANGLE_LIBRARY(i) > 0 ) then;                                          !
                                                                                         !
            IOSEF3 = NTYPE_ANGLE + 1;                                                    !
                                                                                         !
            NTYPE_ANGLE = NTYPE_ANGLE + NTYPE_ANGLE_LIBRARY(i);                          ! UPDATE THE NUMBER OF ANGLE TYPES  
                                                                                         !
            if ( NTYPE_ANGLE > NTYPE_ANGLE_NEW ) then;                                   !
                                                                                         !
                write(icanal,'(a70)') '| NTYPE_ANGLE is greater than the '// &           ! 
                                      'size of allocated arrays'//           &           !
                                      REPEAT(' ',11)//'|';                               !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                write(icanal,'(a43,i8,a3)') '| IOSEF3 | NTYPE_ANGLE '// &                !
                                            '| NTYPE_ANGLE_NEW : ',     &                !
                                            IOSEF3,                     &                !
                                            NTYPE_ANGLE,                &                !
                                            NTYPE_ANGLE_NEW;                             !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                write(icanal,'(a70)') '| --> Stop'//REPEAT(' ',59)//'|';                 !
                                                                                         !
                call CLOSING_PROGRAM(icanal,1);                                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
            CH_ANGLE_STYLE = TRIM(CH_ANGLE_STYLE_LIBRARY(i));                            !
                                                                                         !
            IOSEF6 = NPARAM_ANGLES_LIBRARY(i);                                           !
                                                                                         !
            ANGLE_COEFFS(1:IOSEF6,IOSEF3:NTYPE_ANGLE) = &                                !
            ANGLE_COEFFS_LIBRARY(1:IOSEF6,1:NTYPE_ANGLE_LIBRARY(i),i);                   !
                                                                                         !
            BONDBOND_COEFFS(1:3,IOSEF3:NTYPE_ANGLE) = &                                  !
            BONDBOND_COEFFS_LIBRARY(1:3,1:NTYPE_ANGLE_LIBRARY(i),i);                     !
                                                                                         !
            BONDANGLE_COEFFS(1:4,IOSEF3:NTYPE_ANGLE) = &                                 !
            BONDANGLE_COEFFS_LIBRARY(1:4,1:NTYPE_ANGLE_LIBRARY(i),i);                    !
                                                                                         !
            IANGLE_TYPE_TRANSLATE = IANGLE_TYPE_TRANSLATE + NTYPE_ANGLE_LIBRARY(i);      !
                                                                                         !
        end if                                                                           !
                                                                                         !
        write(icanal,'(a70)') '| Angle types and parameters were updated'// &            !
                              REPEAT(' ',28)//'|';                                       !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!       ### Update dihedral types and parameters ###################################################
                                                                                         !
        if ( NTYPE_DIHEDRAL_LIBRARY(i) > 0 ) then;                                       !
                                                                                         !
            IOSEF5 = NTYPE_DIHEDRAL + 1;                                                 !
                                                                                         !
            NTYPE_DIHEDRAL = NTYPE_DIHEDRAL + NTYPE_DIHEDRAL_LIBRARY(i);                 !
                                                                                         !
            if ( NTYPE_DIHEDRAL > NTYPE_DIHEDRAL_NEW ) then;                             !
                                                                                         !
                write(icanal,'(a70)') '| NTYPE_DIHEDRAL is greater than the '// &        !
                                      'size of allocated arrays'//              &        !
                                      REPEAT(' ',7)//'|';                                !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                write(icanal,'(a49,3i8)') '| IOSEF5 | NTYPE_DIHEDRAL '// &               !
                                          '| NTYPE_DIHEDRAL_NEW : ',     &               !
                                          IOSEF5,                        &               !
                                          NTYPE_DIHEDRAL,                &               !
                                          NTYPE_DIHEDRAL_NEW;                            !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                write(icanal,'(a70)') '| --> Stop'//REPEAT(' ',59)//'|';                 !
                                                                                         !
                call CLOSING_PROGRAM(icanal,1);                                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
            CH_DIHEDRAL_STYLE = TRIM(CH_DIHEDRAL_STYLE_LIBRARY(i));                      !
                                                                                         !
            IOSEF6 = NPARAM_DIHEDRALS_LIBRARY(i);                                        !
                                                                                         !
            DIHEDRAL_COEFFS(1:IOSEF6,IOSEF5:NTYPE_DIHEDRAL) = &                          !
            DIHEDRAL_COEFFS_LIBRARY(1:IOSEF6,1:NTYPE_DIHEDRAL_LIBRARY(i),i);             !
                                                                                         !
            MIDDLEBONDTORSION_COEFFS(1:4,IOSEF5:NTYPE_DIHEDRAL) = &                      !
            MIDDLEBONDTORSION_COEFFS_LIBRARY(1:4,1:NTYPE_DIHEDRAL_LIBRARY(i),i);         !
                                                                                         !
            ENDBONDTORSION_COEFFS(1:8,IOSEF5:NTYPE_DIHEDRAL) = &                         !
            ENDBONDTORSION_COEFFS_LIBRARY(1:8,1:NTYPE_DIHEDRAL_LIBRARY(i),i);            !
                                                                                         !
            ANGLETORSION_COEFFS(1:8,IOSEF5:NTYPE_DIHEDRAL) = &                           !
            ANGLETORSION_COEFFS_LIBRARY(1:8,1:NTYPE_DIHEDRAL_LIBRARY(i),i);              !
                                                                                         !
            ANGLEANGLETORSION_COEFFS(1:3,IOSEF5:NTYPE_DIHEDRAL) = &                      !
            ANGLEANGLETORSION_COEFFS_LIBRARY(1:3,1:NTYPE_DIHEDRAL_LIBRARY(i),i);         !
                                                                                         !
            BONDBOND13_COEFFS(1:3,IOSEF5:NTYPE_DIHEDRAL) = &                             !
            BONDBOND13_COEFFS_LIBRARY(1:3,1:NTYPE_DIHEDRAL_LIBRARY(i),i);                !
                                                                                         !
            IDIHEDRAL_TYPE_TRANSLATE = &                                                 !
            IDIHEDRAL_TYPE_TRANSLATE + &                                                 !
            NTYPE_DIHEDRAL_LIBRARY(i);                                                   !
                                                                                         !
        end if                                                                           !
                                                                                         !
        write(icanal,'(a70)') '| Dihedral angle types and parameters were updated'// &   !
                              REPEAT(' ',19)//'|';                                       !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!       ### Update improper types and parameters ###################################################
                                                                                         !
        if ( NTYPE_IMPROPER_LIBRARY(i) > 0 ) then;                                       !
                                                                                         !
            IOSEF4 = NTYPE_IMPROPER + 1;                                                 !
                                                                                         !
            NTYPE_IMPROPER = NTYPE_IMPROPER + NTYPE_IMPROPER_LIBRARY(i);                 ! UPDATE THE NUMBER OF IMPROPER TYPES
                                                                                         !
            if ( NTYPE_IMPROPER > NTYPE_IMPROPER_NEW ) then;                             !
                                                                                         !
                write(icanal,'(a70)') '| NTYPE_IMPROPER is greater than '// &            !
                                      'the size of allocated arrays'//      &            !
                                      REPEAT(' ',9)//'|';                                !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                write(icanal,'(a49,3i8)') '| IOSEF4 '//              &                   !
                                          '| NTYPE_IMPROPER '//      &                   !
                                          '| NTYPE_IMPROPER_NEW : ', &                   !
                                          IOSEF4,                    &                   !
                                          NTYPE_IMPROPER,            &                   ! 
                                          NTYPE_IMPROPER_NEW;                            !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                write(icanal,'(a70)') '| --> Stop'//REPEAT(' ',59)//'|';                 !
                                                                                         !
                call CLOSING_PROGRAM(icanal,1);                                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
            IMPROPER_COEFFS(1:2,IOSEF4:NTYPE_IMPROPER) = &                               !
            IMPROPER_COEFFS_LIBRARY(1:2,1:NTYPE_IMPROPER_LIBRARY(i),i);                  !
                                                                                         !
            ANGLEANGLE_COEFFS(1:6,IOSEF4:NTYPE_IMPROPER) = &                             !
            ANGLEANGLE_COEFFS_LIBRARY(1:6,1:NTYPE_IMPROPER_LIBRARY(i),i);                !
                                                                                         !
            IIMPROPER_TYPE_TRANSLATE = &                                                 !
            IIMPROPER_TYPE_TRANSLATE + &                                                 !
            NTYPE_IMPROPER_LIBRARY(i);                                                   !
                                                                                         !
        end if                                                                           !
                                                                                         !
        write(icanal,'(a70)') '| Update improper types and parameters'// &               !
                              REPEAT(' ',31)//'|';                                       !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!       ### Update pair coefficients for cross interactions ########################################
                                                                                         !
!       if ( NPAIR_COEFF_CROSS_LIBRARY(i) > 0 ) then;                                    !
                                                                                         !
!           IOSEF4 = NPAIR_COEFF_CROSS + 1;                                              !
                                                                                         !
!           NPAIR_COEFF_CROSS = NPAIR_COEFF_CROSS + NPAIR_COEFF_CROSS_LIBRARY(i);        ! Update the number of pair coefficients for cross interactions
                                                                                         !
!           if ( NPAIR_COEFF_CROSS > NPAIR_COEFF_CROSS_NEW ) then;                       !
                                                                                         !
!               write(icanal,'(a70)') '| NPAIR_COEFF_CROSS is greater than '// &         !
!                                     'the size of allocated arrays'//         &         !
!                                     REPEAT(' ',5)//'|';                                !
                                                                                         ! 
!               write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               write(icanal,'(a55,3i8)') '| IOSEF4 '//                 &                !
!                                         '| NPAIR_COEFF_CROSS '//      &                !
!                                         '| NPAIR_COEFF_CROSS_NEW : ', &                !
!                                         IOSEF4,                       &                !
!                                         NPAIR_COEFF_CROSS,            &                !
!                                         NPAIR_COEFF_CROSS_NEW;                         !
                                                                                         !
!               write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               write(icanal,'(a70)') '| --> Stop'//REPEAT(' ',59)//'|';                 !
                                                                                         !
!               call CLOSING_PROGRAM(icanal,1);                                          !
                                                                                         !
!           end if                                                                       !
                                                                                         !
!           ### Update cross interactions parameters for pair potentials ###########################
                                                                                         !
!           PAIR_COEFF_CROSS(1:2,IOSEF4:NPAIR_COEFF_CROSS) = &                           !
!           PAIR_COEFF_CROSS_LIBRARY(1:2,1:NPAIR_COEFF_CROSS_LIBRARY(i),i);              ! 
                                                                                         !
!           PAIR_ATOMID_CROSS(1:2,IOSEF4:NPAIR_COEFF_CROSS) = &                          !
!!          IPAIR_CROSS_TRANSLATE                           + &                          !
!           IATOM_TYPE_TRANSLATE + &
!           PAIR_ATOMID_CROSS_LIBRARY(1:2,1:NPAIR_COEFF_CROSS_LIBRARY(i),i);             !
                                                                                         !
!           IPAIR_CROSS_TRANSLATE = &                                                    !
!           IPAIR_CROSS_TRANSLATE + &                                                    !
!           NTYPE_ATOM_LIBRARY(i);                                                       !
                                                                                         !
!WWWWWWWWWWWWWWWWWWWWWW
!WWWWWWWWWWWWWWWWWWWWWW
!WWWWWWWWWWWWWWWWWWWWWW

!!          write(icanal,'(a70)') '| Update number of pair interactions '// &            !
!!                                '(cross) and coefficients'//              &            !
!!                                REPEAT(' ',8)//'|';                                    !
                                                                                         !
!!          write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!       end if                                                                           !
                                                                                         !
!       write(icanal,'(a70)') '| Update number of pair interactions '// &                !
!                             '(cross) and coefficients'//              &                !
!                             REPEAT(' ',8)//'|';                                        !
                                                                                         !
!       write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!       ### Update the flag related to the presence of pair potentials #############################
                                                                                         !
!       if ( IPOTENTIAL_CLASS2_LIBRARY(i) == 1 ) IPOTENTIAL_CLASS2 = 1;                  ! 
                                                                                         !
!       ### Update input info parameters ###########################################################
                                                                                         !
        CH_UNITS_STYLE = TRIM(CH_UNITS_STYLE_LIBRARY(i));                                !
                                                                                         !
        CH_ATOM_STYLE =  TRIM(CH_ATOM_STYLE_LIBRARY(i));                                 !
                                                                                         !
        POTENTIAL_CLASS2_CHTYPE = TRIM(POTENTIAL_CLASS2_CHTYPE_LIBRARY(i));              !
                                                                                         !
        CH_BOND_STYLE = TRIM(CH_BOND_STYLE_LIBRARY(i));                                  !
                                                                                         !
        if ( NTYPE_ANGLE_LIBRARY(i) > 0 ) then;                                          !
                                                                                         !
            CH_ANGLE_STYLE = TRIM(CH_ANGLE_STYLE_LIBRARY(i));                            !
                                                                                         !
        end if                                                                           ! 
                                                                                         !
        if ( NTYPE_DIHEDRAL_LIBRARY(i) > 0 ) then;                                       !
                                                                                         !
            CH_DIHEDRAL_STYLE = TRIM(CH_DIHEDRAL_STYLE_LIBRARY(i));                      !
                                                                                         !
        end if                                                                           !
                                                                                         !
        CH_SPECIAL_BONDS = TRIM(CH_SPECIAL_BONDS_LIBRARY(i));                            !
                                                                                         !
        CH_PAIR_MODIFY  = TRIM(CH_PAIR_MODIFY_LIBRARY(i));                               !
                                                                                         !
        CH_KSPACE_STYLE = TRIM(CH_KSPACE_STYLE_LIBRARY(i));                              !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end do                                                                               !
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
!   ### Parameter styles ###########################################################################
                                                                                         !
    IOSEF1 = 70 - 18 - 1 - LEN_TRIM(CH_ATOM_STYLE);                                      !
                                                                                         !
    if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                        !
                                                                                         !
    write(icanal,'(a70)') '| CH_ATOM_STYLE : '//    &                                    !
                          TRIM(CH_ATOM_STYLE)//     &                                    !
                          REPEAT(' ',IOSEF1)//'|';                                       !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    if ( IPOTENTIAL_CLASS2  == 1 ) then;                                                 !
                                                                                         !
        IOSEF1 = 70 - 18 - 1 - LEN_TRIM(CH_PAIR_STYLE);                                  !
                                                                                         !
        if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                    !
                                                                                         !
        write(icanal,'(a70)') '| CH_PAIR_STYLE : '//    &                                !
                              TRIM(CH_PAIR_STYLE)//     &                                !
                              REPEAT(' ',IOSEF1)//'|';                                   !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !   
    end if                                                                               !
                                                                                         !
    if ( NTYPE_BOND > 0 ) then;                                                          !
                                                                                         !
        IOSEF1 = 70 - 18 - 1 - LEN_TRIM(CH_BOND_STYLE);                                  !
                                                                                         !
        if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                    !
                                                                                         !
        write(icanal,'(a70)') '| CH_BOND_STYLE : '//   &                                 !
                              TRIM(CH_BOND_STYLE)//    &                                 !
                              REPEAT(' ',IOSEF1)//'|';                                   !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_ANGLE > 0 ) then;                                                         !
                                                                                         !
        IOSEF1 = 70 - 19 - 1 - LEN_TRIM(CH_ANGLE_STYLE);                                 !
                                                                                         !
        if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                    !
                                                                                         !
        write(icanal,'(a70)') '| CH_ANGLE_STYLE : '//  &                                 !
                              TRIM(CH_ANGLE_STYLE)//   &                                 !
                               REPEAT(' ',IOSEF1)//'|';                                  !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_DIHEDRAL > 0 ) then;                                                      !
                                                                                         !
        IOSEF1 = 70 - 22 - 1 - LEN_TRIM(CH_DIHEDRAL_STYLE);                              !
                                                                                         !
        if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                    !
                                                                                         !
        write(icanal,'(a70)') '| CH_DIHEDRAL_STYLE : '// &                               !
                              TRIM(CH_DIHEDRAL_STYLE)//  &                               !
                              REPEAT(' ',IOSEF1)//'|';                                   !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine INSERT_LIBRARY_MOLECULE

