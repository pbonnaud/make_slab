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

subroutine RANDOM_INSERTION_MOLECULE(icanal,PASSA,PASSB) 

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

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    real (kind=8), dimension(1:3,1:3), intent(in) :: PASSA, PASSB;

!   ************************************************************************************************

    integer (kind=4) :: i, j, k, m;

    integer (kind=4) :: IFOUND;

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

    real (kind=8), dimension(1:3) :: RINO, RI, RI_REF, RAND_ANGLES, RG, DRIJ;

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

    integer (kind=4) :: IOSEF1, IOSEF2, IOSEF3, IOSEF4, IOSEF5, IOSEF6;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the routine #############################################################
                                                                                         !
    CHTITLE = 'Random insertion of molecules';                                           !
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
!   ### Allocate local arrays ######################################################################
                                                                                         !
    if ( NATOM > 0 ) then;                                                               !
                                                                                         !
        allocate(CONFIG_MOLECULEID_OLD(1:NATOM));                                        !
                                                                                         !
        allocate(CONFIG_ATOMID_OLD(1:NATOM));                                            !
                                                                                         !
        allocate(CONFIG_RI_OLD(1:3,1:NATOM));                                            !
                                                                                         !
        allocate(CONFIG_VI_OLD(1:3,1:NATOM));                                            !
                                                                                         ! 
        allocate(CONFIG_ATOM_TYPE_OLD(1:NATOM));                                         !
                                                                                         !
        allocate(CONFIG_QI_OLD(1:NATOM));                                                !
                                                                                         ! 
        allocate(CONFIG_NAT_OLD(1:NATOM));                                               !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NBOND > 0 ) then;                                                               !
                                                                                         !
        allocate(BOND_ATOMID_OLD(1:2,1:NBOND));                                          !
                                                                                         !
        allocate(BOND_TYPE_OLD(1:NBOND));                                                !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NANGLE > 0 ) then;                                                              !
                                                                                         !
        allocate(ANGLE_ATOMID_OLD(1:3,1:NANGLE));                                        !
                                                                                         !
        allocate(ANGLE_TYPE_OLD(1:NANGLE));                                              !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NDIHEDRAL > 0 ) then;                                                           !
                                                                                         !
        allocate(DIHEDRAL_ATOMID_OLD(1:4,1:NDIHEDRAL));                                  !
                                                                                         !
        allocate(DIHEDRAL_TYPE_OLD(1:NDIHEDRAL));                                        !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NIMPROPER > 0 ) then;                                                           !
                                                                                         !
        allocate(IMPROPER_ATOMID_OLD(1:4,1:NIMPROPER));                                  !
                                                                                         !
        allocate(IMPROPER_TYPE_OLD(1:NIMPROPER));                                        !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_ATOM > 0 ) then;                                                          !
                                                                                         !
        allocate(ATOM_MASSE_OLD(1:NTYPE_ATOM));                                          !
                                                                                         !
        allocate(ATOM_LABEL_OLD(1:NTYPE_ATOM));                                          !
                                                                                         !
        allocate(POTENTIAL_CLASS2_OLD(1:2,1:NTYPE_ATOM));                                !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_BOND > 0 )  then;                                                         ! 
                                                                                         !
        allocate(BOND_COEFFS_OLD(1:4,NTYPE_BOND));                                       !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_ANGLE > 0 ) then;                                                         !
                                                                                         !
        allocate(ANGLE_COEFFS_OLD(1:4,1:NTYPE_ANGLE));                                   !
                                                                                         !
        allocate(BONDBOND_COEFFS_OLD(1:3,1:NTYPE_ANGLE));                                !
                                                                                         !
        allocate(BONDANGLE_COEFFS_OLD(1:4,1:NTYPE_ANGLE));                               !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_DIHEDRAL > 0 ) then;                                                      !
                                                                                         !
        allocate(DIHEDRAL_COEFFS_OLD(1:6,1:NTYPE_DIHEDRAL));                             !
                                                                                         !
        allocate(MIDDLEBONDTORSION_COEFFS_OLD(1:4,1:NTYPE_DIHEDRAL));                    !
                                                                                         !
        allocate(ENDBONDTORSION_COEFFS_OLD(1:8,1:NTYPE_DIHEDRAL));                       !
                                                                                         !
        allocate(ANGLETORSION_COEFFS_OLD(1:8,1:NTYPE_DIHEDRAL));                         !
                                                                                         !
        allocate(ANGLEANGLETORSION_COEFFS_OLD(1:3,1:NTYPE_DIHEDRAL));                    !
                                                                                         ! 
        allocate(BONDBOND13_COEFFS_OLD(1:3,1:NTYPE_DIHEDRAL));                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_IMPROPER > 0 ) then;                                                      !
                                                                                         !
        allocate(IMPROPER_COEFFS_OLD(1:2,1:NTYPE_IMPROPER));                             !
                                                                                         !
        allocate(ANGLEANGLE_COEFFS_OLD(1:6,1:NTYPE_IMPROPER));                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NPAIR_COEFF_CROSS > 0 ) then;                                                   !
                                                                                         !
        allocate(PAIR_COEFF_CROSS_OLD(1:2,1:NPAIR_COEFF_CROSS));                         !
                                                                                         !
        allocate(PAIR_ATOMID_CROSS_OLD(1:2,1:NPAIR_COEFF_CROSS));                        !
                                                                                         !
    end if                                                                               !
                                                                                         !
    write(icanal,'(a70)') '| Local arrays were allocated'//REPEAT(' ',40)//'|';          !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Update local arrays ########################################################################
                                                                                         !
    if ( NATOM > 0 ) then;                                                               !
                                                                                         !
        CONFIG_MOLECULEID_OLD(1:NATOM) = CONFIG_MOLECULEID(1:NATOM);                     !
                                                                                         !
        CONFIG_ATOMID_OLD(1:NATOM)     = CONFIG_ATOMID(1:NATOM);                         !
                                                                                         !
        CONFIG_ATOM_TYPE_OLD(1:NATOM)  = CONFIG_ATOM_TYPE(1:NATOM);                      !
                                                                                         !
        CONFIG_RI_OLD(1:3,1:NATOM)     = CONFIG_RI(1:3,1:NATOM);                         !
                                                                                         !   
        CONFIG_VI_OLD(1:3,1:NATOM)     = CONFIG_VI(1:3,1:NATOM);                         !
                                                                                         !
        CONFIG_QI_OLD(1:NATOM)         = CONFIG_QI(1:NATOM);                             !
                                                                                         !
        CONFIG_NAT_OLD(1:NATOM)        = CONFIG_NAT(1:NATOM);                            !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NBOND > 0 ) then;                                                               !
                                                                                         !
        BOND_ATOMID_OLD(1:2,1:NBOND)  = BOND_ATOMID(1:2,1:NBOND);                        !
                                                                                         !
        BOND_TYPE_OLD(1:NBOND)        = BOND_TYPE(1:NBOND);                              !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NANGLE > 0 ) then;                                                              !
                                                                                         !
        ANGLE_ATOMID_OLD(1:3,1:NANGLE) = ANGLE_ATOMID(1:3,1:NANGLE);                     !
                                                                                         !  
        ANGLE_TYPE_OLD(1:NANGLE)       = ANGLE_TYPE(1:NANGLE);                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NDIHEDRAL > 0 ) then;                                                           !
                                                                                         !
        DIHEDRAL_ATOMID_OLD(1:4,1:NDIHEDRAL) = DIHEDRAL_ATOMID(1:4,1:NDIHEDRAL);         !
                                                                                         !
        DIHEDRAL_TYPE_OLD(1:NDIHEDRAL)       = DIHEDRAL_TYPE(1:NDIHEDRAL);               !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NIMPROPER > 0 ) then;                                                           !
                                                                                         !
        IMPROPER_ATOMID_OLD(1:4,1:NIMPROPER) = IMPROPER_ATOMID(1:4,1:NIMPROPER);         !
                                                                                         !
        IMPROPER_TYPE_OLD(1:NIMPROPER)       = IMPROPER_TYPE(1:NIMPROPER);               !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_ATOM > 0 ) then;                                                          !
                                                                                         !
        ATOM_LABEL_OLD(1:NTYPE_ATOM)           = ATOM_LABEL(1:NTYPE_ATOM);               !
                                                                                         !
        ATOM_MASSE_OLD(1:NTYPE_ATOM)           = ATOM_MASSE(1:NTYPE_ATOM);               !
                                                                                         !
        POTENTIAL_CLASS2_OLD(1:2,1:NTYPE_ATOM) = POTENTIAL_CLASS2(1:2,1:NTYPE_ATOM);     !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_BOND > 0 ) then;                                                          !
                                                                                         !
        BOND_COEFFS_OLD(1:4,1:NTYPE_BOND) = BOND_COEFFS(1:4,1:NTYPE_BOND);               !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_ANGLE > 0 ) then;                                                         !
                                                                                         !
        ANGLE_COEFFS_OLD(1:4,1:NTYPE_ANGLE)     = ANGLE_COEFFS(1:4,1:NTYPE_ANGLE);       !
                                                                                         !
        BONDBOND_COEFFS_OLD(1:3,1:NTYPE_ANGLE)  = BONDBOND_COEFFS(1:3,1:NTYPE_ANGLE);    !
                                                                                         !
        BONDANGLE_COEFFS_OLD(1:4,1:NTYPE_ANGLE) = BONDANGLE_COEFFS(1:4,1:NTYPE_ANGLE);   !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_DIHEDRAL > 0 ) then;                                                      !
                                                                                         !
        DIHEDRAL_COEFFS_OLD(1:6,1:NTYPE_DIHEDRAL)          = &                           ! 
        DIHEDRAL_COEFFS(1:6,1:NTYPE_DIHEDRAL);                                           !
                                                                                         !
        MIDDLEBONDTORSION_COEFFS_OLD(1:4,1:NTYPE_DIHEDRAL) = &                           !
        MIDDLEBONDTORSION_COEFFS(1:4,1:NTYPE_DIHEDRAL);                                  !
                                                                                         !
        ENDBONDTORSION_COEFFS_OLD(1:8,1:NTYPE_DIHEDRAL)    = &                           !
        ENDBONDTORSION_COEFFS(1:8,1:NTYPE_DIHEDRAL);                                     !
                                                                                         !
        ANGLETORSION_COEFFS_OLD(1:8,1:NTYPE_DIHEDRAL)      = &                           !
        ANGLETORSION_COEFFS(1:8,1:NTYPE_DIHEDRAL);                                       !
                                                                                         !
        ANGLEANGLETORSION_COEFFS_OLD(1:3,1:NTYPE_DIHEDRAL) = &                           !
        ANGLEANGLETORSION_COEFFS(1:3,1:NTYPE_DIHEDRAL);                                  !
                                                                                         !
        BONDBOND13_COEFFS_OLD(1:3,1:NTYPE_DIHEDRAL)        = &                           !
        BONDBOND13_COEFFS(1:3,1:NTYPE_DIHEDRAL);                                         !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_IMPROPER > 0 ) then;                                                      !
                                                                                         !
        IMPROPER_COEFFS_OLD(1:2,1:NTYPE_IMPROPER)   = &                                  !
        IMPROPER_COEFFS(1:2,1:NTYPE_IMPROPER);                                           !
                                                                                         !
        ANGLEANGLE_COEFFS_OLD(1:6,1:NTYPE_IMPROPER) = &                                  !
        ANGLEANGLE_COEFFS(1:6,1:NTYPE_IMPROPER);                                         !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NPAIR_COEFF_CROSS > 0 ) then;                                                   !
                                                                                         !
        PAIR_COEFF_CROSS_OLD(1:2,1:NPAIR_COEFF_CROSS) = &                                !
        PAIR_COEFF_CROSS(1:2,1:NPAIR_COEFF_CROSS);                                       !
                                                                                         !
        PAIR_ATOMID_CROSS_OLD(1:2,1:NPAIR_COEFF_CROSS) = &                               !
        PAIR_ATOMID_CROSS(1:2,1:NPAIR_COEFF_CROSS);                                      !
                                                                                         !
    end if                                                                               ! 
                                                                                         !
    write(icanal,'(a70)') '| Local arrays were updated'//REPEAT(' ',42)//'|';            !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Deallocate global arrays ###################################################################
                                                                                         !
    if ( NATOM > 0 ) then;                                                               !
                                                                                         !
        deallocate(CONFIG_MOLECULEID);                                                   !
                                                                                         !
        deallocate(CONFIG_ATOMID);                                                       !
                                                                                         !
        deallocate(CONFIG_RI);                                                           !
                                                                                         !
        deallocate(CONFIG_VI);                                                           !
                                                                                         !
        deallocate(CONFIG_ATOM_TYPE);                                                    !
                                                                                         !
        deallocate(CONFIG_QI);                                                           !
                                                                                         !
        deallocate(CONFIG_NAT);                                                          !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NBOND > 0 ) then;                                                               !
                                                                                         !
        deallocate(BOND_ATOMID);                                                         !
                                                                                         !
        deallocate(BOND_TYPE);                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NANGLE > 0 ) then;                                                              !
                                                                                         !
        deallocate(ANGLE_ATOMID);                                                        !
                                                                                         !
        deallocate(ANGLE_TYPE);                                                          !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NDIHEDRAL > 0 ) then;                                                           !
                                                                                         !
        deallocate(DIHEDRAL_ATOMID);                                                     !
                                                                                         !
        deallocate(DIHEDRAL_TYPE);

    end if

    if ( NIMPROPER > 0 ) then;

        deallocate(IMPROPER_ATOMID);

        deallocate(IMPROPER_TYPE);

    end if

    if ( NTYPE_ATOM > 0 ) then;

        deallocate(ATOM_MASSE);                                              !

        deallocate(ATOM_LABEL);                                              !

        deallocate(POTENTIAL_CLASS2);                                    !

    end if
    !
    if ( NTYPE_BOND > 0 ) then;

        deallocate(BOND_COEFFS);                                           !

    end if                                                                                         !

    if ( NTYPE_ANGLE > 0 ) then;                                                                                     

        deallocate(ANGLE_COEFFS);                                       !

        deallocate(BONDBOND_COEFFS);                                    !

        deallocate(BONDANGLE_COEFFS);                                   !

    end if
                                                                                         !
    if ( NTYPE_DIHEDRAL > 0 ) then;                                                      !
                                                                                         !
        deallocate(DIHEDRAL_COEFFS);                                                     !
                                                                                         !
        deallocate(MIDDLEBONDTORSION_COEFFS);                                            !
                                                                                         !
        deallocate(ENDBONDTORSION_COEFFS);                                               !
                                                                                         !
        deallocate(ANGLETORSION_COEFFS);                                                 !
                                                                                         !
        deallocate(ANGLEANGLETORSION_COEFFS);                                            !
                                                                                         !  
        deallocate(BONDBOND13_COEFFS);                                                   !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_IMPROPER > 0 ) then;                                                      !
                                                                                         !
        deallocate(IMPROPER_COEFFS);                                                     !
                                                                                         !
        deallocate(ANGLEANGLE_COEFFS);                                                   !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NPAIR_COEFF_CROSS > 0 ) then;                                                   !
                                                                                         !
        deallocate(PAIR_COEFF_CROSS);                                                    !
                                                                                         !
        deallocate(PAIR_ATOMID_CROSS);                                                   !
                                                                                         !
    end if                                                                               !
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
    if ( iworking_file_insertion == 1 ) then;                                            ! If the instertion is made in an existing simulation box
                                                                                         !
        do i = 1, NFILE_LIBRARY;                                                         ! LOOP OVER THE MOLECULE TYPE READ IN THE LIBRARY 
                                                                                         !
            NATOM_NEW = NATOM_NEW              + &                                       !
                        NATOM_LIBRARY(i)       * &                                       !
                        NMOLECULE_INSERTION(i);                                          !
                                                                                         !
            if ( NBOND_LIBRARY(i)     > 0 ) NBOND_NEW              = &                   !
                                            NBOND_NEW              + &                   !
                                            NBOND_LIBRARY(i)       * &                   !
                                            NMOLECULE_INSERTION(i);                      !
                                                                                         !
            if ( NANGLE_LIBRARY(i)    > 0 ) NANGLE_NEW             = &                   !
                                            NANGLE_NEW             + &                   !
                                            NANGLE_LIBRARY(i)      * &                   !
                                            NMOLECULE_INSERTION(i);                      !
                                                                                         !
            if ( NDIHEDRAL_LIBRARY(i) > 0 ) NDIHEDRAL_NEW          = &                   !
                                            NDIHEDRAL_NEW          + &                   !
                                            NDIHEDRAL_LIBRARY(i)   * &                   !
                                            NMOLECULE_INSERTION(i);                      !
                                                                                         !
            if ( NIMPROPER_LIBRARY(i) > 0 ) NIMPROPER_NEW          = &                   !
                                            NIMPROPER_NEW          + &                   !
                                            NIMPROPER_LIBRARY(i)   * &                   !
                                            NMOLECULE_INSERTION(i);                      !
                                                                                         !
            NTYPE_ATOM_NEW = NTYPE_ATOM_NEW     + &                                      !
                             NTYPE_ATOM_LIBRARY(i);                                      !
                                                                                         !
            if ( NTYPE_BOND_LIBRARY(i)     > 0 ) NTYPE_BOND_NEW     = &                  !
                                                 NTYPE_BOND_NEW     + &                  !
                                                 NTYPE_BOND_LIBRARY(i);                  !
                                                                                         !
            if ( NTYPE_ANGLE_LIBRARY(i)    > 0 ) NTYPE_ANGLE_NEW    = &                  !
                                                 NTYPE_ANGLE_NEW    + &                  !
                                                 NTYPE_ANGLE_LIBRARY(i);                 !
                                                                                         !
            if ( NTYPE_DIHEDRAL_LIBRARY(i) > 0 ) NTYPE_DIHEDRAL_NEW = &                  !
                                                 NTYPE_DIHEDRAL_NEW + &                  !
                                                 NTYPE_DIHEDRAL_LIBRARY(i);              !
                                                                                         !
            if ( NTYPE_IMPROPER_LIBRARY(i) > 0 ) NTYPE_IMPROPER_NEW = &                  !
                                                 NTYPE_IMPROPER_NEW + &                  !
                                                 NTYPE_IMPROPER_LIBRARY(i);              !
                                                                                         !
            if ( NPAIR_COEFF_CROSS_LIBRARY(i) > 0 ) then;                                !
                                                                                         !
                NPAIR_COEFF_CROSS_NEW =       &                                          !
                NPAIR_COEFF_CROSS_NEW +       &                                          !
                NPAIR_COEFF_CROSS_LIBRARY(i);                                            !
                                                                                         !
             end if                                                                      !
                                                                                         !
         end do                                                                          !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Write the new size of configuration arrays #################################################
                                                                                         !
    write(icanal,'(a23,i8,a39)') '| NATOM_NEW          : ', &                            !
                                 NATOM_NEW,                 &                            !
                                 REPEAT(' ',38)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    if ( NBOND_NEW > 0 ) then;                                                           !
                                                                                         !
        write(icanal,'(a23,i8,a39)') '| NBOND_NEW          : ', &                        !
                                     NBOND_NEW,                 &                        !
                                     REPEAT(' ',38)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    write(icanal,'(a23,i8,a39)') '| NANGLE_NEW         : ', NANGLE_NEW,         &        !
                                 REPEAT(' ',38)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a23,i8,a39)') '| NDIHEDRAL_NEW      : ', NDIHEDRAL_NEW,      &        !
                                 REPEAT(' ',38)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    if ( NIMPROPER_NEW > 0 ) then;                                                       !
                                                                                         !
        write(icanal,'(a23,i8,a39)') '| NIMPROPER_NEW      : ', &                        !
                                     NIMPROPER_NEW,             &                        !
                                     REPEAT(' ',38)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    write(icanal,'(a23,i8,a39)') '| NTYPE_ATOM_NEW     : ', &                            !
                                 NTYPE_ATOM_NEW,            &                            !
                                 REPEAT(' ',38)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a23,i8,a39)') '| NTYPE_BOND_NEW     : ', NTYPE_BOND_NEW,     &        !
                                 REPEAT(' ',38)//'|';                                    !

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

    write(icanal,'(a23,i8,a39)') '| NTYPE_ANGLE_NEW    : ', NTYPE_ANGLE_NEW,    &        !
                                 REPEAT(' ',38)//'|';                                    !

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

    write(icanal,'(a23,i8,a39)') '| NTYPE_DIHEDRAL_NEW : ', NTYPE_DIHEDRAL_NEW, &        !
                                 REPEAT(' ',38)//'|';                                    !

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

    if ( NTYPE_IMPROPER_NEW > 0 ) then;                                                  !
                                                                                         !
        write(icanal,'(a23,i8,a39)') '| NTYPE_IMPROPER_NEW : ', &                        !
                                     NTYPE_IMPROPER_NEW,        &                        !
                                     REPEAT(' ',38)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NPAIR_COEFF_CROSS_NEW > 0 ) then;                                               !
                                                                                         !
        write(icanal,'(a26,i8,a36)') '| NPAIR_COEFF_CROSS_NEW : ', &                     !
                                     NPAIR_COEFF_CROSS_NEW,        &                     !
                                     REPEAT(' ',35)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Re-allocate global arrays ##################################################################
                                                                                         ! 
    if ( NATOM_NEW > 0 ) then;                                                           !
                                                                                         !
        allocate(CONFIG_MOLECULEID(1:NATOM_NEW));                                        !
                                                                                         !
        allocate(CONFIG_ATOMID(1:NATOM_NEW));                                            !
                                                                                         !
        allocate(CONFIG_RI(1:3,1:NATOM_NEW));                                            !
                                                                                         !
        allocate(CONFIG_VI(1:3,1:NATOM_NEW));                                            !
                                                                                         !
        allocate(CONFIG_ATOM_TYPE(1:NATOM_NEW));                                         !
                                                                                         !
        allocate(CONFIG_QI(1:NATOM_NEW));                                                !
                                                                                         !
        allocate(CONFIG_NAT(1:NATOM_NEW));                                               !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NBOND_NEW > 0 ) then;                                                           !
                                                                                         !
        allocate(BOND_ATOMID(1:2,1:NBOND_NEW));                                          !
                                                                                         !
        allocate(BOND_TYPE(1:NBOND_NEW));                                                !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NANGLE_NEW > 0 ) then;                                                          !
                                                                                         !
        allocate(ANGLE_ATOMID(1:3,1:NANGLE_NEW));                                        !
                                                                                         !
        allocate(ANGLE_TYPE(1:NANGLE_NEW));                                              !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NDIHEDRAL_NEW > 0 ) then;                                                       !
                                                                                         !
        allocate(DIHEDRAL_ATOMID(1:4,1:NDIHEDRAL_NEW));                                  !
                                                                                         !
        allocate(DIHEDRAL_TYPE(1:NDIHEDRAL_NEW));                                        !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NIMPROPER_NEW > 0 ) then;                                                       !
                                                                                         !
        allocate(IMPROPER_ATOMID(1:4,1:NIMPROPER_NEW));                                  !
                                                                                         !
        allocate(IMPROPER_TYPE(1:NIMPROPER_NEW));                                        !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_ATOM_NEW > 0 ) then;                                                      !
                                                                                         !
        allocate(ATOM_MASSE(1:NTYPE_ATOM_NEW));                                          !
                                                                                         !
        allocate(ATOM_LABEL(1:NTYPE_ATOM_NEW));                                          !
                                                                                         !
        allocate(POTENTIAL_CLASS2(1:2,1:NTYPE_ATOM_NEW));                                !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_BOND_NEW > 0 ) then;                                                      !
                                                                                         ! 
        allocate(BOND_COEFFS(1:4,NTYPE_BOND_NEW));                                       !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_ANGLE_NEW > 0 ) then;                                                     !
                                                                                         !
        allocate(ANGLE_COEFFS(1:4,1:NTYPE_ANGLE_NEW));                                   !
                                                                                         !
        allocate(BONDBOND_COEFFS(1:3,1:NTYPE_ANGLE_NEW));                                !
                                                                                         !
        allocate(BONDANGLE_COEFFS(1:4,1:NTYPE_ANGLE_NEW));                               !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_DIHEDRAL_NEW > 0 ) then;                                                  !
                                                                                         ! 
        allocate(DIHEDRAL_COEFFS(1:6,1:NTYPE_DIHEDRAL_NEW));                             !
                                                                                         !
        allocate(MIDDLEBONDTORSION_COEFFS(1:4,1:NTYPE_DIHEDRAL_NEW));                    !
                                                                                         !
        allocate(ENDBONDTORSION_COEFFS(1:8,1:NTYPE_DIHEDRAL_NEW));                       !
                                                                                         !
        allocate(ANGLETORSION_COEFFS(1:8,1:NTYPE_DIHEDRAL_NEW));                         !
                                                                                         !
        allocate(ANGLEANGLETORSION_COEFFS(1:3,1:NTYPE_DIHEDRAL_NEW));                    !
                                                                                         !
        allocate(BONDBOND13_COEFFS(1:3,1:NTYPE_DIHEDRAL_NEW));                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_IMPROPER_NEW > 0 ) then;                                                  !
                                                                                         !
        allocate(IMPROPER_COEFFS(1:2,1:NTYPE_IMPROPER_NEW));                             !
                                                                                         !
        allocate(ANGLEANGLE_COEFFS(1:6,1:NTYPE_IMPROPER_NEW));                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NPAIR_COEFF_CROSS_NEW > 0 ) then;                                               !
                                                                                         !
        allocate(PAIR_COEFF_CROSS(1:2,1:NPAIR_COEFF_CROSS_NEW));                         !
                                                                                         !
        allocate(PAIR_ATOMID_CROSS(1:2,1:NPAIR_COEFF_CROSS_NEW));                        !
                                                                                         !
    end if                                                                               !
                                                                                         !
    write(icanal,'(a70)') '| Global arrays were re-allocated'//REPEAT(' ',36)//'|';      !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Update global arrays #######################################################################
                                                                                         !
    if ( NATOM > 0 ) then;                                                               !
                                                                                         !
        CONFIG_MOLECULEID(1:NATOM) = CONFIG_MOLECULEID_OLD(1:NATOM);                     !
                                                                                         !
        CONFIG_ATOMID(1:NATOM)     = CONFIG_ATOMID_OLD(1:NATOM);                         !
                                                                                         !
        CONFIG_ATOM_TYPE(1:NATOM)  = CONFIG_ATOM_TYPE_OLD(1:NATOM);                      !

        CONFIG_RI(1:3,1:NATOM)     = CONFIG_RI_OLD(1:3,1:NATOM);                         !

        CONFIG_VI(1:3,1:NATOM)     = CONFIG_VI_OLD(1:3,1:NATOM);                         !

        CONFIG_QI(1:NATOM)         = CONFIG_QI_OLD(1:NATOM);                             !

        CONFIG_NAT(1:NATOM)        = CONFIG_NAT_OLD(1:NATOM);                            !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NBOND > 0 ) then;                                                               !
                                                                                         !
        BOND_ATOMID(1:2,1:NBOND)  = BOND_ATOMID_OLD(1:2,1:NBOND);                        !

        BOND_TYPE(1:NBOND)        = BOND_TYPE_OLD(1:NBOND);                              !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NANGLE > 0 ) then;                                                              !
                                                                                         !
        ANGLE_ATOMID(1:3,1:NANGLE) = ANGLE_ATOMID_OLD(1:3,1:NANGLE);                     !

        ANGLE_TYPE(1:NANGLE)       = ANGLE_TYPE_OLD(1:NANGLE);                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NDIHEDRAL > 0 ) then;                                                           !
                                                                                         !
        DIHEDRAL_ATOMID(1:4,1:NDIHEDRAL) = DIHEDRAL_ATOMID_OLD(1:4,1:NDIHEDRAL);         !

        DIHEDRAL_TYPE(1:NDIHEDRAL)       = DIHEDRAL_TYPE_OLD(1:NDIHEDRAL);               !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NIMPROPER > 0 ) then;                                                           !
                                                                                         !
        IMPROPER_ATOMID(1:4,1:NIMPROPER) = IMPROPER_ATOMID_OLD(1:4,1:NIMPROPER);         !

        IMPROPER_TYPE(1:NIMPROPER)       = IMPROPER_TYPE_OLD(1:NIMPROPER);               !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_ATOM > 0 ) then;                                                          !
                                                                                         !
        ATOM_LABEL(1:NTYPE_ATOM)           = ATOM_LABEL_OLD(1:NTYPE_ATOM);               !

        ATOM_MASSE(1:NTYPE_ATOM)           = ATOM_MASSE_OLD(1:NTYPE_ATOM);               !

        POTENTIAL_CLASS2(1:2,1:NTYPE_ATOM) = POTENTIAL_CLASS2_OLD(1:2,1:NTYPE_ATOM);     !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_BOND > 0 ) then;                                                          !
                                                                                         !
        BOND_COEFFS(1:4,1:NTYPE_BOND) = BOND_COEFFS_OLD(1:4,1:NTYPE_BOND);               !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_ANGLE > 0 ) then;                                                         !
                                                                                         !
        ANGLE_COEFFS(1:4,1:NTYPE_ANGLE)     = ANGLE_COEFFS_OLD(1:4,1:NTYPE_ANGLE);       !

        BONDBOND_COEFFS(1:3,1:NTYPE_ANGLE)  = BONDBOND_COEFFS_OLD(1:3,1:NTYPE_ANGLE);    !

        BONDANGLE_COEFFS(1:4,1:NTYPE_ANGLE) = BONDANGLE_COEFFS_OLD(1:4,1:NTYPE_ANGLE);   !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_DIHEDRAL > 0 ) then;                                                      !
                                                                                         !
        DIHEDRAL_COEFFS(1:6,1:NTYPE_DIHEDRAL)          = &                               ! 
        DIHEDRAL_COEFFS_OLD(1:6,1:NTYPE_DIHEDRAL);                                       !

        MIDDLEBONDTORSION_COEFFS(1:4,1:NTYPE_DIHEDRAL) = &                               !
        MIDDLEBONDTORSION_COEFFS_OLD(1:4,1:NTYPE_DIHEDRAL);                              !

        ENDBONDTORSION_COEFFS(1:8,1:NTYPE_DIHEDRAL)    = &                               !
        ENDBONDTORSION_COEFFS_OLD(1:8,1:NTYPE_DIHEDRAL);                                 !

        ANGLETORSION_COEFFS(1:8,1:NTYPE_DIHEDRAL)      = &                               !
        ANGLETORSION_COEFFS_OLD(1:8,1:NTYPE_DIHEDRAL);                                   !

        ANGLEANGLETORSION_COEFFS(1:3,1:NTYPE_DIHEDRAL) = &                               !
        ANGLEANGLETORSION_COEFFS_OLD(1:3,1:NTYPE_DIHEDRAL);                              !

        BONDBOND13_COEFFS(1:3,1:NTYPE_DIHEDRAL)        = &                               !
        BONDBOND13_COEFFS_OLD(1:3,1:NTYPE_DIHEDRAL);                                     !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_IMPROPER > 0 ) then;                                                      !
                                                                                         !
        IMPROPER_COEFFS(1:2,1:NTYPE_IMPROPER)   = &                                      !
        IMPROPER_COEFFS_OLD(1:2,1:NTYPE_IMPROPER);                                       !
                                                                                         !
        ANGLEANGLE_COEFFS(1:6,1:NTYPE_IMPROPER) = &                                      !
        ANGLEANGLE_COEFFS_OLD(1:6,1:NTYPE_IMPROPER);                                     !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NPAIR_COEFF_CROSS > 0 ) then;                                                   !
                                                                                         !
        PAIR_COEFF_CROSS(1:2,1:NPAIR_COEFF_CROSS) = &                                    !
        PAIR_COEFF_CROSS_OLD(1:2,1:NPAIR_COEFF_CROSS);                                   !
                                                                                         !
        PAIR_ATOMID_CROSS(1:2,1:NPAIR_COEFF_CROSS) = &                                   !
        PAIR_ATOMID_CROSS_OLD(1:2,1:NPAIR_COEFF_CROSS);                                  !
                                                                                         !
    end if                                                                               ! 
                                                                                         !
    write(icanal,'(a70)') '| Global arrays were updated'//REPEAT(' ',41)//'|';           !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
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
    if ( iworking_file_insertion == 0 ) then;                                            !
                                                                                         !
        NATOM     = 0;                                                                   !
                                                                                         !
        NBOND     = 0;                                                                   !
                                                                                         !
        NANGLE    = 0;                                                                   !
                                                                                         !
        NDIHEDRAL = 0;                                                                   !
                                                                                         !
        NIMPROPER = 0;                                                                   ! 
                                                                                         !
        NTYPE_ATOM     = 0;                                                              !
                                                                                         !
        NTYPE_BOND     = 0;                                                              !
                                                                                         !
        NTYPE_ANGLE    = 0;                                                              !
                                                                                         !
        NTYPE_DIHEDRAL = 0;                                                              !
                                                                                         !
        NTYPE_IMPROPER = 0;                                                              !
                                                                                         !
        NPAIR_COEFF_CROSS = 0;                                                           !
                                                                                         !
    else                                                                                 !
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
    end if                                                                               !
                                                                                         !
!   ### Set the initial molecule ID for insertion ##################################################
                                                                                         !
    IMOLECULEID = 0;                                                                     !
                                                                                         !
    if ( iworking_file_insertion == 1 ) then;                                            !
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
    end if                                                                               !
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
!   ### Insertion of molecules read in the library #################################################
                                                                                         !
    do i = 1, NFILE_LIBRARY;                                                             ! LOOP OVER THE MOLECULE TYPE READ IN THE LIBRARY
                                                                                         !
        write(icanal,'(a17,i4,a49)') '| LIBRARY FILE # ', i, ' '//REPEAT('/',47)//'|';   !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        SUM_MASSES = 0.0d0;                                                              ! INITIALIZATION OF THE SUM OF MASSES
                                                                                         !
        RG(1:3) = 0.0d0;                                                                 ! INITIALIZATION OF THE CENTER OF MASS
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
        do j = 1, NATOM_LIBRARY(i);                                                      ! LOOP OVER THE NUMBER OF ATOMS IN THE MOLECULE TYPE I
                                                                                         !
            do k = 1, NTYPE_ATOM_LIBRARY(i);                                             ! LOOP OVER THE NUMBER OF ATOM TYPES IN THE MOLECULE I
                                                                                         !
                if ( TRIM(CONFIG_NAT_LIBRARY(j,i)) == &                                  !
                     TRIM(ATOM_LABEL_LIBRARY(k,i)) ) then;                               !
                                                                                         !
                    RG(1:3) = RG(1:3) + ATOM_MASSES_LIBRARY(k,i) * &                     ! Here, for the computation of the center of mass, we do not take into account periodic
                                        CONFIG_RI_LIBRARY(1:3,j,i);                      ! bundary conditions. We assume that library molecules are not cut in the space
                                                                                         !
                    SUM_MASSES = SUM_MASSES + ATOM_MASSES_LIBRARY(k,i);                  !
                                                                                         !
                    EXIT;                                                                !
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
        end do                                                                           !
                                                                                         !
        INV_SUM_MASSES = 1.0d0 / SUM_MASSES;                                             !
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
            if ( INSERTION_METHOD_LIBRARY(i) == 1 ) then;                                !
                                                                                         !
                RINO(1:3) = (/grnd() - 0.5d0, grnd() - 0.5d0, grnd() - 0.5d0/);          ! SET RANDOMLY THE COORDINATE OF THE CENTER OF MASS OF THE MOLECULE TO INSERT
                                                                                         !
                RI(1:3) = MATMUL(PASSA(1:3,1:3),RINO(1:3));                              ! DEFINE COORDINATES IN THE CARTESIAN AXIS
                                                                                         !
            else if ( INSERTION_METHOD_LIBRARY(i) == 2 ) then;                           ! From defined center of mass
                                                                                         !
!               RG(1:3) = 0.5d0 * CELL_AXIS_LIBRARY(1:3,i);                              !
                                                                                         !
!               RI(3) = CELL_AXIS_LIBRARY(3,i);                                          !
                                                                                         !
                IOSEF1 = IINSERTED + 1;                                                  !
                                                                                         !
!               write(icanal,'(a27,f12.6)') '| Geometrical center [A] : ', RG(1:3);

!               write(icanal,*);

!               write(icanal,*) '| Defined center [A] : ', MOLECULE_COM_LIBRARY(1:3,IOSEF1,i);

!               write(icanal,*);

                RI(1:3) = MOLECULE_COM_LIBRARY(1:3,IOSEF1,i);                  !
                                                                                         !
                write(icanal,'(a24,3f12.6,a10)') '| Center location [A] : ', &           !
                                                 RI(1:3),                    &           !
                                                 REPEAT(' ',9)//'|';                     !

!               write(icanal,*);

            end if                                                                       !
                                                                                         !
!           if ( IINSERTED == 1) stop;

!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!           ### Check the location of the new center for fefined regions ###########################
                                                                                         !
            if ( NREGIONS_INSERTION > 0 ) then;                                          ! IF REGIONS HAVE BEEN DEFINED FOR THE INSERTION OF MOLECULES, THEN
                                                                                         !
                do j = 1, NREGIONS_INSERTION;                                            ! LOOP OVER THE REGIONS THAT WERE DEFINED
                                                                                         !
                    if ( ( RI(REGION_AXIS(j)) > REGION_BOUNDS(1,j) ) .AND. &             ! TEST THE NEW GENERATED COORDINATES WITH RESPECT TO THE DEFINED REGIONS
                         ( RI(REGION_AXIS(j)) < REGION_BOUNDS(2,j) ) ) then;             ! IF COORDINATES ARE IN THE DEFINED REGION, THEN
                                                                                         !
                        IREGIONS_INSERTION = 1;                                          ! SET IREGIONS_INSERTION TO 1
                                                                                         !
                    else                                                                 !
                                                                                         !
                        IREGIONS_INSERTION = 0;                                          ! 
                                                                                         !
                    end if                                                               !
                                                                                         !
                end do                                                                   !
                                                                                         !
            else                                                                         ! IF REGIONS WEER NOT DEFINED, THEN 
                                                                                         !
                IREGIONS_INSERTION = 1;                                                  ! SET IREGIONS_INSERTION TO 1
                                                                                         !
            end if                                                                       !
                                                                                         !
!           write(icanal,*) 'IREGIONS_INSERTION : ', IREGIONS_INSERTION;
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
            if ( IREGIONS_INSERTION == 0 ) CYCLE;                                        ! CYCLE IF COORDINATES ARE NOT LOCATED IN THE DEFINED REGIONS
                                                                                         !
            IOVERLAPPING = 0;                                                            !
                                                                                         !
            do j = 1, NATOM;                                                             ! LOOP OVER THE NUMBER OF ATOMS IN THE MOLECULAR CONFIGURATION
                                                                                         !
                DRIJ(1:3) = CONFIG_RI(1:3,j) - RI(1:3);                                  ! 
                                                                                         !
                call APPLY_PBC(DRIJ(1:3),       &                                        !
                               DRIJ(1:3),       &                                        ! 
                               PASSA(1:3,1:3),  &                                        !
                               PASSB(1:3,1:3));                                          !
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
            write(icanal,*) 'IOVERLAPPING : ', IOVERLAPPING;

            write(icanal,'(a70)') '| Overlapping (center)   [OK]'//REPEAT(' ',40)//'|';  !
                                                                                         !
!           write(icanal,*);

!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
            if ( IOVERLAPPING > 0 ) CYCLE;                                               !
                                                                                         ! 
            do k = 1, 1000;                                                              ! LOOP OVER THE NUMBER OF TRIES TO ROTATE THE MOLECULE FOR THE INSERTION
                                                                                         !
                RAND_ANGLES(1:3) = 0.0d0;                                                !
                                                                                         !
                if ( INSERTION_METHOD_LIBRARY(i) == 1 ) then;                            !
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
                    call APPLY_PBC(DRIJ(1:3),      &                                     !
                                   DRIJ(1:3),      &                                     !
                                   PASSA(1:3,1:3), &                                     !
                                   PASSB(1:3,1:3));                                      !
                                                                                         !
                    DRIJ(1:3) = MATMUL(MATROT(1:3,1:3),DRIJ(1:3));                       !
                                                                                         !
                    CONFIG_RI(1:3,IATOM) = RI(1:3) + DRIJ(1:3);                          !
                                                                                         !
                    call APPLY_PBC(CONFIG_RI(1:3,IATOM), &                               !
                                   CONFIG_RI(1:3,IATOM), &                               ! 
                                   PASSA(1:3,1:3),       &                               ! 
                                   PASSB(1:3,1:3));                                      !

!                   write(icanal,*) CONFIG_RI(1:3,IATOM);
                                                                                         ! 
                    CONFIG_NAT(IATOM) = CONFIG_NAT_LIBRARY(j,i);                         ! SET THE NATURE OF THE CONSIDERED ATOM
                                                                                         !
                    CONFIG_QI(IATOM) = CONFIG_QI_LIBRARY(j,i);                           ! SET THE PARTIAL CHARGE OF THE CONSIDERED ATOM
                                                                                         !
!                   write(icanal,*) 'OK QI';

!                   stop;

                    if ( ( CONFIG_VI_LIBRARY(1,j,i) /= 0.0d0 ) .AND. &                   !
                         ( CONFIG_VI_LIBRARY(2,j,i) /= 0.0d0 ) .AND. &                   !
                         ( CONFIG_VI_LIBRARY(3,j,i) /= 0.0d0 ) ) then;                   !
                                                                                         !
                        CONFIG_VI(1:3,IATOM) = CONFIG_VI_LIBRARY(1:3,j,i);               !

                    else

                        CONFIG_VI(1:3,IATOM) = 0.0d0;    
                                                                                         !
                    end if                                                               !
                                                                                         !
!                   write(icanal,*) CONFIG_VI(1:3,IATOM);

!                   stop;                                                                !
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
                    call CHECK_OVERLAPPING_ATOMS(icanal,NATOM+1,                  &      !
                                                        IATOM,                    &      !
                                                        PASSA(1:3,1:3),           &      !
                                                        PASSB(1:3,1:3),           &      !
                                                        1.5d0,                    &      !
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

            write(icanal,*) 'IOVERLAPPING ', IOVERLAPPING;

!           write(icanal,*) 'OK 1 | IOVERLAPPING : ', IOVERLAPPING;
            write(icanal,*) '| Overlapping (rotation) [OK]';

            write(icanal,*);

!           if ( IINSERTED == 1) stop;

                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
            if ( IOVERLAPPING > 0 ) CYCLE;                                               ! IF THERE IS OVERLAP BETWEEN THE INSERTED MOLECULE AND THE MOLECULAR CONFIGURATION, CYCLE
                                                                                         !
            IMOLECULEID = IMOLECULEID + 1;                                               ! UPDATE THE MOLECULE ID
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
            IOSEF1 = NATOM + 1;                                                          !
                                                                                         !
            NATOM = IATOM;                                                               ! UPDATE THE NUMBER OF ATOMS IN THE FINAL CONFIGURATION
                                                                                         !
            if ( NATOM > NATOM_NEW ) then;                                               !
                                                                                         !
                write(icanal,*) '| NATOM IS GREATER THAN THE SIZE OF ALLOCATED ARRAYS';  !
                                                                                         !
                write(icanal,'(a31,3i8)') '| IOSEF1 | NATOM | NATOM_NEW : ', &           !
                                          IOSEF1, NATOM, NATOM_NEW;                      !
                                                                                         !
                write(icanal,*) '--> STOP';                                              !
                                                                                         !
                stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
            CONFIG_MOLECULEID(IOSEF1:NATOM) = IMOLECULEID;                               ! APPLY THE NEW MOLECULE ID 
                                                                                         !
            IINSERTED = IINSERTED + 1;                                                   !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
                IOSEF1 = IOSEF1 - 1;                                                     !  
                                                                                         !
                if ( NBOND_LIBRARY(i) > 0 ) then;                                        ! IF THE NUMBER OF BONDS OF THE INSERTED MOLECULE IS GREATER THAN 0, THEN
                                                                                         !
                    IOSEF2 = NBOND + 1;                                                  !
                                                                                         !
                    NBOND  = NBOND + NBOND_LIBRARY(i);                                   ! UPDATE THE NUMBER OF BONDS IN THE FINAL CONFIGURATION 
                                                                                         !
                    if ( NBOND > NBOND_NEW ) then;

                        write(icanal,*) '| NBOND IS GREATER THAN THE SIZE OF ALLOCATED ARRAYS';

                        write(icanal,*) '| IOSEF2 | NBOND | NBOND_NEW : ', IOSEF2, NBOND, NBOND_NEW;

                        write(icanal,*) '--> STOP';

                        stop; !//////////////////////////////////////////////////////////!
                                                                                         !
                    end if                                                               !
                                                                                         !
                    BOND_ATOMID(1:2,IOSEF2:NBOND) = &                                    !
                    IOSEF1                        + &                                    !
                    BOND_ATOMID_LIBRARY(1:2,1:NBOND_LIBRARY(i),i);                       !
                                                                                         !
                    BOND_TYPE(IOSEF2:NBOND) = &                                          !
                    IBOND_TYPE_TRANSLATE    + &                                          !
                    BOND_TYPE_LIBRARY(1:NBOND_LIBRARY(i),i);                             !
                                                                                         !
                end if                                                                   !
                                                                                         !
!               write(icanal,*) 'OK BONDS';                                              !
                                                                                         !
!               stop;

                if ( NANGLE_LIBRARY(i) > 0 ) then;                                       !
                                                                                         !
                    IOSEF3 = NANGLE + 1;                                                 !
                                                                                         !
                    NANGLE = NANGLE + NANGLE_LIBRARY(i);                                 ! UPDATE THE NUMBER OF ANGLES IN THE FINAL CONFIGURATION
                                                                                         !

                    if ( NANGLE > NANGLE_NEW ) then;  

                        write(icanal,*) '| NANGLE IS GREATER THAN THE SIZE OF ALLOCATED ARRAYS';
                        write(icanal,*) '| IOSEF2 | NANGLE | NANGLE_NEW : ', IOSEF3, NANGLE, NANGLE_NEW;
                        write(icanal,*) '--> STOP';
                        stop;

                    end if        

                    ANGLE_ATOMID(1:3,IOSEF3:NANGLE) = &                                  !
                    IOSEF1                          + &                                  !
                    ANGLE_ATOMID_LIBRARY(1:3,1:NANGLE_LIBRARY(i),i);                     !
                                                                                         ! 
                    ANGLE_TYPE(IOSEF3:NANGLE) = &                                        !
                    IANGLE_TYPE_TRANSLATE     + &                                        !
                    ANGLE_TYPE_LIBRARY(1:NANGLE_LIBRARY(i),i);                           !
                                                                                         !
                end if                                                                   !

!               write(icanal,*) 'OK ANGLES';

!               stop;
                                                                                         !
                if ( NDIHEDRAL_LIBRARY(i) > 0 ) then;                                    !
                                                                                         !
                    IOSEF5 = NDIHEDRAL + 1;                                              !
                                                                                         !
                    NDIHEDRAL = NDIHEDRAL + NDIHEDRAL_LIBRARY(i);                        ! UPDATE THE NUMBER OF DIHEDRALS IN THE FINAL CONFIGURATION
                                                                                         !
                    if ( NDIHEDRAL > NDIHEDRAL_NEW ) then;

                        write(icanal,*) '| NDIHEDRAL IS GREATER THAN THE SIZE OF ALLOCATED ARRAYS';
                        write(icanal,*) '| IOSEF5 | NDIHEDRAL | NDIHEDRAL_NEW : ', IOSEF5, NDIHEDRAL, NDIHEDRAL_NEW;
                        write(icanal,*) '--> STOP';
                        stop;

                    end if

!                   write(icanal,*) IOSEF5, NDIHEDRAL, NDIHEDRAL_NEW;

                    DIHEDRAL_ATOMID(1:4,IOSEF5:NDIHEDRAL) = &                            !
                    IOSEF1                                + &                            !
                    DIHEDRAL_ATOMID_LIBRARY(1:4,1:NDIHEDRAL_LIBRARY(i),i);               !
                                                                                         !
                    DIHEDRAL_TYPE(IOSEF5:NDIHEDRAL) = &                                  !
                    IDIHEDRAL_TYPE_TRANSLATE        + &                                  !
                    DIHEDRAL_TYPE_LIBRARY(1:NDIHEDRAL_LIBRARY(i),i);                     !
                                                                                         !
                end if                                                                   !

!               write(icanal,*) 'OK DIHEDRALS';

!               stop;
                !
                if ( NIMPROPER_LIBRARY(i) > 0 ) then;                                    !
                                                                                         !
                    IOSEF4 = NIMPROPER + 1;                                              !
                                                                                         !
                    NIMPROPER = NIMPROPER + NIMPROPER_LIBRARY(i);                        ! UPDATE THE NUMBER OF IMPROPERS IN THE FINAL CONFIGURATION
                                                                                         !
                    if ( NIMPROPER > NIMPROPER_NEW ) then;                               !
                                                                                         !
                        write(icanal,*) '| NIMPROPER IS GREATER THAN THE SIZE OF ALLOCATED ARRAYS';
                                                                                         !
                        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                  !
                                                                                         !
                        write(icanal,*) '| IOSEF4 | NIMPROPER | NIMPROPER_NEW : ', &     !
                                        IOSEF4, NIMPROPER, NIMPROPER_NEW;                !

                        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                  !
                                                                                         !
                        write(icanal,*) '--> STOP';                                      !
                                                                                         !
                        stop; !//////////////////////////////////////////////////////////!
                                                                                         !
                    end if                                                               !
                                                                                         !
                    IMPROPER_ATOMID(1:4,IOSEF4:NIMPROPER) = &                            !
                    IOSEF1                                + &                            !
                    IMPROPER_ATOMID_LIBRARY(1:4,1:NIMPROPER_LIBRARY(i),i);               !
                                                                                         !
                    IMPROPER_TYPE(IOSEF4:NIMPROPER) = &                                  !
                    IIMPROPER_TYPE_TRANSLATE        + &                                  !
                    IMPROPER_TYPE_LIBRARY(1:NIMPROPER_LIBRARY(i),i);                     ! 
                                                                                         !
                end if                                                                   !

!               write(icanal,*) 'OK IMPROPERS';

!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
!               ### Update the current status of molecule insertion in the simulation box ##########
                                                                                         !
                write(icanal,'(a14,i8,f12.2,a36)') '| IINSERTED : ',             &       !  
                                                   IINSERTED,                    &       !
                                                   100.0d0 * REAL(IINSERTED) /   &       !
                                                   REAL(NMOLECULE_INSERTION(i)), &       !
                                                   REPEAT(' ',35)//'|';                  !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
!           end if                                                                       !
                                                                                         !
        end do                                                                           !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        IOSEF1 = NTYPE_ATOM + 1;                                                         ! UPDATE THE NUMBER OF ATOMS IN THE FINAL MOLECULAR CONFIGURATION
                                                                                         !
!       ### Update the number of atom types ########################################################
                                                                                         !
        NTYPE_ATOM = NTYPE_ATOM + NTYPE_ATOM_LIBRARY(i);                                 ! UPDATE THE NUMBER OF ATOM TYPES     
                                                                                         !
        if ( NTYPE_ATOM > NTYPE_ATOM_NEW ) then;                                         !
                                                                                         !
            write(icanal,*) '| NTYPE_ATOM IS GREATER THAN THE SIZE OF ALLOCATED ARRAYS';
            write(icanal,*) '| IOSEF1 | NTYPE_ATOM | NTYPE_ATOM_NEW : ', IOSEF1, NTYPE_ATOM, NTYPE_ATOM_NEW;
            write(icanal,*) '--> STOP';
                                                                                         !
            stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end if                                                                           !
                                                                                         !
!       write(icanal,*) 'OK 2';

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
!       write(icanal,*) 'OK 3';

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
        IATOM_TYPE_TRANSLATE = IATOM_TYPE_TRANSLATE + NTYPE_ATOM_LIBRARY(i);             !
                                                                                         !
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
                write(icanal,*) '| NTYPE_BOND IS GREATER THAN THE SIZE OF ALLOCATED ARRAYS';

                write(icanal,*) '| IOSEF2 | NTYPE_BOND | NTYPE_BOND_NEW : ', IOSEF2, NTYPE_BOND, NTYPE_BOND_NEW;

                write(icanal,*) '--> STOP';
                                                                                         !
                stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
!           if ( NPARAM_BONDS_LIBRARY(i) > NPARAM_BONDS ) then;                          !
                                                                                         !
!               NPARAM_BONDS = NPARAM_BONDS_LIBRARY(i);                                  !
                                                                                         !
!           end if                                                                       !
                                                                                         !
            CH_BOND_STYLE = CH_BOND_STYLE_LIBRARY(i);

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
        if ( NTYPE_ANGLE_LIBRARY(i) > 0 ) then;                                          !
                                                                                         !
            IOSEF3 = NTYPE_ANGLE + 1;                                                    !
                                                                                         !
            NTYPE_ANGLE = NTYPE_ANGLE + NTYPE_ANGLE_LIBRARY(i);                          ! UPDATE THE NUMBER OF ANGLE TYPES  
                                                                                         !
            if ( NTYPE_ANGLE > NTYPE_ANGLE_NEW ) then;                                   !
                                                                                         !
                write(icanal,*) '| NTYPE_ANGLE IS GREATER THAN THE SIZE OF ALLOCATED ARRAYS';

                write(icanal,*) '| IOSEF3 | NTYPE_ANGLE | NTYPE_ANGLE_NEW : ', IOSEF3, NTYPE_ANGLE, NTYPE_ANGLE_NEW;

                write(icanal,*) '--> STOP';

                stop;

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
                write(icanal,*) '| NTYPE_DIHEDRAL IS GREATER THAN THE SIZE OF ALLOCATED ARRAYS';

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

                write(icanal,*) '| IOSEF5 | NTYPE_DIHEDRAL | NTYPE_DIHEDRAL_NEW : ', IOSEF5, NTYPE_DIHEDRAL, NTYPE_DIHEDRAL_NEW;

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

                write(icanal,*) '--> STOP';

                stop;

            end if

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
                write(icanal,*) '| IOSEF4 | NTYPE_IMPROPER | NTYPE_IMPROPER_NEW : ', &
                                IOSEF4, &
                                NTYPE_IMPROPER, &
                                NTYPE_IMPROPER_NEW;

                write(icanal,*) '--> STOP';

                stop;

            end if

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
                write(icanal,*) '| IOSEF4 | NPAIR_COEFF_CROSS | NPAIR_COEFF_CROSS_NEW : ', &
                                IOSEF4,             &
                                NPAIR_COEFF_CROSS, &
                                NPAIR_COEFF_CROSS_NEW;
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !

                write(icanal,*) '--> STOP';

                close(icanal);

                stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
!           ### Update cross interactions parameters for pair potentials ###########################
                                                                                         !
            PAIR_COEFF_CROSS(1:2,IOSEF4:NPAIR_COEFF_CROSS) = &                           !
            PAIR_COEFF_CROSS_LIBRARY(1:2,1:NPAIR_COEFF_CROSS_LIBRARY(i),i);              ! 
                                                                                         !
            PAIR_ATOMID_CROSS(1:2,IOSEF4:NPAIR_COEFF_CROSS) = &                          !
            IPAIR_CROSS_TRANSLATE                           + &                          !
            PAIR_ATOMID_CROSS_LIBRARY(1:2,1:NPAIR_COEFF_CROSS_LIBRARY(i),i);             !
                                                                                         !
            IPAIR_CROSS_TRANSLATE        = &                                             !
            IPAIR_CROSS_TRANSLATE        + &                                             !
            NTYPE_ATOM_LIBRARY(i);                                                       !
!           NPAIR_COEFF_CROSS_LIBRARY(i);                                                !
                                                                                         !
            write(icanal,'(a70)') '| Update number of pair interactions '// &            !
                                  '(cross) and coefficients'//              &            !
                                  REPEAT(' ',8)//'|';                                    !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        end if                                                                           !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!QQQQQQQQQQQQQQQQQQQQQQQQQQQ
!QQQQQQQQQQQQQQQQQQQQQQQQQQQ
!QQQQQQQQQQQQQQQQQQQQQQQQQQQ
!QQQQQQQQQQQQQQQQQQQQQQQQQQQ
!QQQQQQQQQQQQQQQQQQQQQQQQQQQ
!QQQQQQQQQQQQQQQQQQQQQQQQQQQ
!QQQQQQQQQQQQQQQQQQQQQQQQQQQ
!QQQQQQQQQQQQQQQQQQQQQQQQQQQ
!QQQQQQQQQQQQQQQQQQQQQQQQQQQ



!       ###

        if ( IPOTENTIAL_CLASS2_LIBRARY(i) == 1 ) IPOTENTIAL_CLASS2 = 1;                  ! 
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

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

    write(icanal,'(a19,i8,a43)') '| NBOND          : ', &                                !
                                 NBOND,                 &                                !
                                 REPEAT(' ',42)//'|';                                    !

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

    write(icanal,'(a19,i8,a43)') '| NANGLE         : ', &                                !
                                 NANGLE,                &                                !
                                 REPEAT(' ',42)//'|';                                    !

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

    write(icanal,'(a19,i8,a43)') '| NDIHEDRAL      : ', &                                !
                                 NDIHEDRAL,             &                                !
                                 REPEAT(' ',42)//'|';                                    !

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

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
    write(icanal,'(a19,i8,a43)') '| NTYPE_ATOM     : ', &                                !
                                 NTYPE_ATOM,            &                                !
                                 REPEAT(' ',42)//'|';                                    !

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

    write(icanal,'(a19,i8,a43)') '| NTYPE_BOND     : ', &                                !
                                 NTYPE_BOND,            &                                !
                                 REPEAT(' ',42)//'|';                                    !

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

    write(icanal,'(a19,i8,a43)') '| NTYPE_ANGLE    : ', &                                !
                                 NTYPE_ANGLE,           &                                !
                                 REPEAT(' ',42)//'|';                                    !

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

    write(icanal,'(a19,i8,a43)') '| NTYPE_DIHEDRAL : ', &                                !
                                 NTYPE_DIHEDRAL,        &                                !
                                 REPEAT(' ',42)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
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
    IOSEF1 = 70 - 18 - 1 - LEN_TRIM(CH_PAIR_STYLE);                                      !
                                                                                         !
    if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                        !
                                                                                         !
    write(icanal,'(a70)') '| CH_PAIR_STYLE : '//    &                                    !
                          TRIM(CH_PAIR_STYLE)//     &                                    !
                          REPEAT(' ',IOSEF1)//'|';                                       !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
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
end subroutine RANDOM_INSERTION_MOLECULE

