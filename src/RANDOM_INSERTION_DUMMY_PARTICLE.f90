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

subroutine RANDOM_INSERTION_DUMMY_PARTICLE(icanal,        &
!                                          IBUILD_METHOD, & 
                                           PASSA,         &
                                           PASSB) 

!   ************************************************************************************************
!   **                                  READ THE INPUT FILE                                       **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                    : CANAL ON WHICH OUTPUT DATA ARE WRITTEN                         **
!   ** NFILE_LIBRARY             : NUMBER OF FILES COMING FROM THE LIBRARY                        **
!   ** NMOLECULE_INSERTION       : NUMBER OF INSERTED MOLECULES FOR A GIVE FILE                   **
!   ** CHNAME_FILE_LIBRARY       : FILE NAME CONTAINING THE MOLECULAR CONFIGURATION OF THE        **
!   **                             MOLECULE                                                       **
!   ** NATOM_FILE_LIBRARY        : NUMBER OF ATOMS IN THE FILES FROM THE LIBRARY                  **
!   ** NBOND_LIBRARY             : NUMBER OF BONDS IN THE FILES FROM THE LIBRARY                  **
!   ** NANGLE_LIBRARY            : NUMBER OF ANGLES IN THE FILES FROM THE LIBRARY                 **
!   ** NDIHEDRAL_LIBRARY         : NUMBER OF DIHEDRAL ANGLES IN THE FILES FROM THE LIBRARY        **
!   ** NIMPROPER_LIBRARY         : NUMBER OF IMPROPERS IN THE FILES FROM THE LIBRARY              **
!   ** ATOM_LABEL_LIBRARY        :  **
!   ** CONFIG_NAT_LIBRARY        :  **
!   ** CONFIG_QI_LIBRARY         : CHARGES OF ATOMS IN THE FILES FROM THE LIBRARY  **
!   ** CONFIG_RI_LIBRARY         : COORDINATES OF ATOMS IN THE FILES FROM THE LIBRARY **
!   ** POTENTIAL_CLASS2_LIBRARY  : LENNARD-JONES POTENTIAL PARAMETERS IN THE FILES FROM THE LIBRARY **
!   ** BOND_COEFFS_LIBRARY       :
!   ** ANGLE_COEFFS_LIBRARY      :
!   ** BONDBOND_COEFFS_LIBRARY   :
!   ** BONDANGLE_COEFFS_LIBRARY  : 
!   ** DIHEDRAL_COEFFS_LIBRARY   :
!   ** MIDDLEBONDTORSION_COEFFS_LIBRARY : 
!   ** ENDBONDTORSION_COEFFS_LIBRARY    :
!   ** ANGLETORSION_COEFFS_LIBRARY      :
!   ** ANGLEANGLETORSION_COEFFS_LIBRARY :
!   ** BONDBOND13_COEFFS_LIBRARY        :
!   ** IMPROPER_COEFFS_LIBRARY   :
!   ** ANGLEANGLE_COEFFS_LIBRARY :
!   ** PASSA                     :  **
!   ** PASSB                     :  **
!   ** NTYPE_ATOM_LIBRARY        : NUMBER OF ATOM TYPES IN THE LIBRARY FILE **
!   ** NTYPE_BOND_LIBRARY        : NUMBER OF BOND TYPES IN THE LIBRARY FILE **
!   ** NTYPE_ANGLE_LIBRARY       : NUMBER OF ANGLE TYPES IN THE LIBRARY FILE **
!   ** NTYPE_DIHEDRAL_LIBRARY    : NUMBER OF DIHEDRAL TYPES IN THE LIBRARY FILE **
!   ** NTYPE_IMPROPER_LIBRARY    : NUMBER OF IMPROPER TYPES IN THE LIBRARY FILE **
!   ** ATOM_MASSES_LIBRARY       :
!   ** CONFIG_ATOMID_LIBRARY     :
!   ** CONFIG_ATOM_TYPE_LIBRARY  : 
!   ** BOND_TYPE_LIBRARY         :
!   ** ANGLE_TYPE_LIBRARY        :
!   ** DIHEDRAL_TYPE_LIBRARY     : 
!   ** IMPROPER_TYPE_LIBRARY     :
!   ** BOND_ATOMID_LIBRARY       :
!   ** ANGLE_ATOMID_LIBRARY      :
!   ** DIHEDRAL_ATOMID_LIBRARY   :
!   ** IMPROPER_ATOMID_LIBRARY   :
!   ** IPOTENTIAL_CLASS2_LIBRARY : FLAG
!   **
!   ** NATOM                    : NUMBER OF ATOMS IN THE FINAL CONFIGURATION  ** 
!   ** NBOND                    : NUMBER OF BONDS IN THE FINAL CONFIGURATION   **
!   ** NANGLE                   : NUMBER OF ANGLES IN THE FINAL CONFIGURATION   **
!   ** NDIHEDRAL                : NUMBER OF DIHEDRALS IN THE FINAL CONFIGURATION **                       
!   ** NIMPROPER                : NUMBER OF IMPROPERS IN THE FINAL CONFIGURATION **
!   ** NTYPE_ATOM               : NUMBER OF ATOM TYPES IN THE FINAL MOLECULAR CONFIGURATION          **
!   ** NTYPE_BOND               : NUMBER OF BOND TYPES IN THE FINAL MOLECULAR CONFIGURATION          **
!   ** NTYPE_ANGLE              : NUMBER OF ANGLE TYPES IN THE FINAL MOLECULAR CONFIGURATION **
!   ** NTYPE_DIHEDRAL           : NUMBER OF DIHEDRAL TYPES IN THE FINAL MOLECULAR CONFIGURATION **
!   ** NTYPE_IMPROPER           : NUMBER OF IMPROPER TYPES IN THE FINAL MOLECULAR CONFIGURATION **
!   ** ATOM_MASSE               : MASS OF ATOMS IN THE MOLECULAR CONFIGURATION                       **
!   ** ATOM_LABEL               : LIST OF ATOM LABEL IN THE FINAL MOLECULAR CONFIGURATION            **
!   ** CONFIG_NAT               :
!   ** CONFIG_QI                : CHARGES OFR ATOMS IN THE FINAL MOLECULAR CONFIGURATION
!   ** CONFIG_RI                : COORDINATES OF ATOMS IN THE FINAL MOLECULAR CONFIGURATION
!   ** POTENTIAL_CLASS2         : LENNARD-JONES POTENTIAL PARAMETERS IN THE FINAL MOLECULAR CONFIGURATION **
!   ** BOND_COEFFS              :
!   ** ANGLE_COEFFS             : 
!   ** BONDBOND_COEFFS          :
!   ** BONDANGLE_COEFFS         :
!   ** DIHEDRAL_COEFFS          :
!   ** MIDDLEBONDTORSION_COEFFS :
!   ** ENDBONDTORSION_COEFFS    :
!   ** ANGLETORSION_COEFFS      : 
!   ** ANGLEANGLETORSION_COEFFS :
!   ** BONDBOND13_COEFFS        :
!   ** IMPROPER_COEFFS          :
!   ** ANGLEANGLE_COEFFS        :
!   ** CONFIG_ATOMID            :
!   ** CONFIG_MOLECULEID        :
!   ** CONFIG_ATOM_TYPE         : 
!   ** BOND_TYPE                :
!   ** ANGLE_TYPE               :
!   ** DIHEDRAL_TYPE            :
!   ** IMPROPER_TYPE            :
!   ** BOND_ATOMID              :
!   ** ANGLE_ATOMID             :
!   ** DIHEDRAL_ATOMID          :
!   ** IMPROPER_ATOMID          :
!   ** IPOTENTIAL_CLASS2        : FLAG
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

!   integer (kind=4), intent(in) :: icanal, IBUILD_METHOD;
    integer (kind=4), intent(in) :: icanal;

    real (kind=8), dimension(1:3,1:3), intent(in) :: PASSA, PASSB;

!   ************************************************************************************************

    integer (kind=4) :: i, j, k, m;

    integer (kind=4) :: IFOUND;

    integer (kind=4) :: NLOOP_DUMMY_INSERTION, ILOOP_DUMMY;

    integer (kind=4) :: IINSERTED, IATOM, IOVERLAPPING;

    integer (kind=4) :: IATOM_TYPE_TMP, IMOLECULEID;

    real (kind=8) :: grnd;

    real (kind=8) :: SUM_MASSES, INV_SUM_MASSES, RIJ, EPSILON_LJ96;

    real (kind=8), dimension(1:3) :: RINO, RI, RAND_ANGLES, RG, DRIJ;

    real (kind=8), dimension(1:3,1:3) :: MATROT;

!   ************************************************************************************************

    integer (kind=4) :: NATOM_NEW,      &
                        NTYPE_ATOM_NEW;

    integer (kind=4), allocatable, dimension(:) :: CONFIG_ATOM_TYPE_NEW,   &
                                                   CONFIG_ATOMID_NEW,      &
                                                   CONFIG_MOLECULEID_NEW;

    real (kind=8), allocatable, dimension(:) :: ATOM_MASSE_NEW;

    real (kind=8), allocatable, dimension(:) :: CONFIG_QI_NEW;

    real (kind=8), allocatable, dimension(:,:) :: POTENTIAL_CLASS2_NEW;

    real (kind=8), allocatable, dimension(:,:) :: CONFIG_RI_NEW,  &
                                                  CONFIG_VI_NEW;

    character (len=20), allocatable, dimension(:) :: ATOM_LABEL_NEW;

    character (len=20), allocatable, dimension(:) :: CONFIG_NAT_NEW;

!   ************************************************************************************************

    integer (kind=4) :: IOSEF1, IOSEF2, IOSEF3, IOSEF4, IOSEF5;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
    CHTITLE = 'RANDOM INSERTION DUMMY MOLECULE';                                         !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
    NLOOP_DUMMY_INSERTION = 1000000;                                                     !
                                                                                         !
    IMOLECULEID = 0;                                                                     !
                                                                                         !
    write(icanal,'(a18,i8,a46)') '| IBUILD_METHOD : ', &                                 !
                                 IBUILD_METHOD,        &                                 !
                                 REPEAT(' ',45)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
    write(icanal,'(a25,i8,a37)') '| NTYPE_DUMMY_PARTICLE : ', &                          !
                                NTYPE_DUMMY_PARTICLE,        &                           !
                                REPEAT(' ',36)//'|';                                     !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### INITIALIZATION OF THE RANDOM GENERATOR ##########################################!
                                                                                         !
    call sgrnd(RANDOM_GENERATOR_SEED);                                                   !
                                                                                         ! 
!   ### SET POTENTIAL PARAMETER FOR DUMMY PARTICLES #####################################!
                                                                                         !
    EPSILON_LJ96 = TARGET_TEMPK * kB * Na * 0.001d0;                                     ! IN [kJ/mol] 
    EPSILON_LJ96 = EPSILON_LJ96 * KJMOL_TO_KCALMOL;                                      ! IN [KCal/mol] 
                                                                                         !
!   ### ALLACATE ARRAYS #################################################################!
                                                                                         !
!   if ( IBUILD_METHOD == 3 ) then;

        NATOM_NEW = NATOM;

        NTYPE_ATOM_NEW = NTYPE_ATOM;

        do i = 1, NTYPE_DUMMY_PARTICLE;

            NATOM_NEW = NATOM_NEW + NDUMMY_PARTICLE(i);

            NTYPE_ATOM_NEW = NTYPE_ATOM_NEW + 1;

        end do

        write(icanal,*) '| NATOM_NEW : ', NATOM_NEW;

        write(icanal,*) '| NTYPE_ATOM_NEW : ', NTYPE_ATOM_NEW;

        allocate(CONFIG_RI_NEW(1:3,1:NATOM_NEW));

        allocate(CONFIG_NAT_NEW(1:NATOM_NEW));

        allocate(CONFIG_QI_NEW(1:NATOM_NEW));

        allocate(CONFIG_VI_NEW(1:3,1:NATOM_NEW));

        allocate(CONFIG_ATOMID_NEW(1:NATOM_NEW));

        allocate(CONFIG_ATOM_TYPE_NEW(1:NATOM_NEW));

        allocate(CONFIG_MOLECULEID_NEW(1:NATOM_NEW));

        allocate(ATOM_LABEL_NEW(1:NTYPE_ATOM_NEW));                                      ! UPDATE THE LIST OF LABEL OF ATOMS 
                                                                                         !
        allocate(ATOM_MASSE_NEW(1:NTYPE_ATOM_NEW));                                      ! UPDATE THE LIST OF MASSES 
                                                                                         !
        allocate(POTENTIAL_CLASS2_NEW(1:2,1:NTYPE_ATOM_NEW));                            !
                                                                                         !
        CONFIG_RI_NEW(1:3,1:NATOM)             = CONFIG_RI(1:3,1:NATOM);

        CONFIG_NAT_NEW(1:NATOM)                = CONFIG_NAT(1:NATOM);
 
        CONFIG_QI_NEW(1:NATOM)                 = CONFIG_QI(1:NATOM); 

        CONFIG_VI_NEW(1:3,1:NATOM)             = CONFIG_VI(1:3,1:NATOM);

        CONFIG_ATOMID_NEW(1:NATOM)             = CONFIG_ATOMID(1:NATOM);

        CONFIG_ATOM_TYPE_NEW(1:NATOM)          = CONFIG_ATOM_TYPE(1:NATOM); 

        CONFIG_MOLECULEID_NEW(1:NATOM)         = CONFIG_MOLECULEID(1:NATOM);

        ATOM_LABEL_NEW(1:NTYPE_ATOM)           = ATOM_LABEL(1:NTYPE_ATOM);

        ATOM_MASSE_NEW(1:NTYPE_ATOM)           = ATOM_MASSE(1:NTYPE_ATOM);

        POTENTIAL_CLASS2_NEW(1:2,1:NTYPE_ATOM) = POTENTIAL_CLASS2(1:2,1:NTYPE_ATOM);

!   end if                                                                               !

!   ### SET THE INITIAL MOLECULE ID FOR DUMMY PARTICLE INSERTION ########################!
                                                                                         !
    if ( IBUILD_METHOD == 1 ) then;

    if ( NFILE_LIBRARY > 0 ) then;

        do i = 1, NFILE_LIBRARY;

            IMOLECULEID = IMOLECULEID + NMOLECULE_INSERTION(i);

        end do
                                                                                         !
    end if                                                                               !

    else if ( IBUILD_METHOD == 3 ) then;

        do i = 1, NATOM;

            if ( IMOLECULEID < CONFIG_MOLECULEID(i) ) IMOLECULEID = &
                                                      CONFIG_MOLECULEID(i);

        end do

    end if

    write(icanal,'(a16,i8,a46)') '| IMOLECULEID : ', &                                   !
                                 IMOLECULEID,        &                                   !
                                 REPEAT(' ',45)//'|';                                    !
                                                                                         !
!   ### INSWRTION OF THE DUMMY PARTICLES ################################################!

    IATOM = NATOM;

    IATOM_TYPE_TMP = NTYPE_ATOM;                                                         !
                                                                                         !
    do i = 1, NTYPE_DUMMY_PARTICLE;                                                      ! LOOP OVER THE DUMMY PARTICLE TYPES 
                                                                                         !
        write(icanal,'(a70)') '| CHNAME_DUMMY_PARTICLE : ',   &                          !
                              TRIM(CHNAME_DUMMY_PARTICLE(i)), &                          !
                              REPEAT(' ',10)//'|';                                       ! 
                                                                                         !
        write(icanal,'(a23,i8,a39)') '| NDUMMY_PARTICLE    : ', &                        !
                                     NDUMMY_PARTICLE(i),        &                        !
                                     REPEAT(' ',38)//'|';                                !
                                                                                         !
        write(icanal,'(a26,f12.6,a32)') '| RADIUS_DUMMY_PARTICLE : ', &                  !
                                        RADIUS_DUMMY_PARTICLE(i),     &                  !
                                        REPEAT(' ',31)//'|';                             !
                                                                                         !
        write(icanal,'(a26,f12.6,a32)') '| MASSE_DUMMY_PARTICLE : ', &                   !
                                        MASSE_DUMMY_PARTICLE(i),     &                   !
                                        REPEAT(' ',31)//'|';                             ! 

        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

        IINSERTED = 0;

        ILOOP_DUMMY = 0;

        IATOM_TYPE_TMP = IATOM_TYPE_TMP + 1;                                             ! UPDATE THE NUMBER OF ATOM TYPES     

        do while ( ( IINSERTED   < NDUMMY_PARTICLE(i)    ) .AND. &
                   ( ILOOP_DUMMY < NLOOP_DUMMY_INSERTION ) )

            ILOOP_DUMMY = ILOOP_DUMMY + 1;

            RINO(1) = grnd() - 0.5d0;
            RINO(2) = grnd() - 0.5d0;
            RINO(3) = grnd() - 0.5d0;

            write(icanal,'(a9,3f15.6,a46)') '| RINO : ', RINO(1:3), REPEAT(' ',45)//'|'; !

            RI(1:3) = MATMUL(PASSA(1:3,1:3),RINO(1:3));

            write(icanal,'(a9,3f15.6)') '| RI   : ', RI(1:3);

            IOVERLAPPING = 0;

            do j = 1, IATOM;                                                             ! LOOP OVER ATOMS IN THE MOLECULAR CONFIGURATION 
                                                                                         !
                DRIJ(1:3) = CONFIG_RI(1:3,j) - RI(1:3);                                  !
                                                                                         !
                call APPLY_PBC(DRIJ(1:3),DRIJ(1:3),PASSA(1:3,1:3),PASSB(1:3,1:3));       !
                                                                                         !
                RIJ = DSQRT( DOT_PRODUCT(DRIJ(1:3),DRIJ(1:3)) );                         !
                                                                                         !
                if ( RIJ <= RADIUS_DUMMY_PARTICLE(i) ) then;
                    IOVERLAPPING = IOVERLAPPING + 1;
                    EXIT;
                end if
            end do

            if ( IOVERLAPPING > 0 ) CYCLE;

            IINSERTED = IINSERTED + 1;

            IATOM = IATOM + 1;

            CONFIG_RI_NEW(1:3,IATOM) = RI(1:3);

            CONFIG_NAT_NEW(IATOM)    = TRIM(CHNAME_DUMMY_PARTICLE(i));

            CONFIG_QI_NEW(IATOM)     = 0.0d0;

            CONFIG_VI_NEW(1:3,IATOM) = 0.0d0;

            CONFIG_ATOMID_NEW(IATOM) = IATOM;

            CONFIG_ATOM_TYPE_NEW(IATOM) = IATOM_TYPE_TMP;
 
            IMOLECULEID = IMOLECULEID + 1;                                                   ! UPDATE THE MOLECULE ID
                                                                                             !
            CONFIG_MOLECULEID_NEW(IATOM) = IMOLECULEID;                                      ! APPLY THE NEW MOLECULE ID 

            IOSEF1 = NATOM;

            write(icanal,'(a14,i8,f12.2,a36)') '| IINSERTED : ', &
                                               IINSERTED,        &
                                               100.0d0 * REAL(IINSERTED) / REAL(NDUMMY_PARTICLE(i)), &
                                               REPEAT(' ',35)//'|';
        end do

        IOSEF1 = NTYPE_ATOM     + 1;                                                     ! UPDATE THE NUMBER OF ATOMS IN THE FINAL MOLECULAR CONFIGURATION

        write(icanal,*) 'OK 1';

        write(icanal,*) 'OK 2';

        ATOM_LABEL_NEW(IATOM_TYPE_TMP) = TRIM(CHNAME_DUMMY_PARTICLE(i));                 ! UPDATE THE LIST OF LABEL OF ATOMS 
                                                                                         !
        ATOM_MASSE_NEW(IATOM_TYPE_TMP) = MASSE_DUMMY_PARTICLE(i);                        ! UPDATE THE LIST OF MASSES 
                                                                                         !
        POTENTIAL_CLASS2_NEW(1,IATOM_TYPE_TMP) = EPSILON_LJ96;                           !
                                                                                         !
        POTENTIAL_CLASS2_NEW(2,IATOM_TYPE_TMP) = 2.0d0                    * &            !
                                                 RADIUS_DUMMY_PARTICLE(i) * &            !
                                                 1.5d0**(1.0d0/3.0d0);                   !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end do                                                                               !

!   ####

!   if ( IBUILD_METHOD == 3 ) then;

        deallocate(CONFIG_RI);

        deallocate(CONFIG_NAT);

        deallocate(CONFIG_QI);

        deallocate(CONFIG_VI);

        deallocate(CONFIG_ATOMID);

        deallocate(CONFIG_ATOM_TYPE);

        deallocate(CONFIG_MOLECULEID);

        deallocate(ATOM_LABEL);                                      ! UPDATE THE LIST OF LABEL OF ATOMS 
                                                                                         !
        deallocate(ATOM_MASSE);                                      ! UPDATE THE LIST OF MASSES 
                                                                                         !
        deallocate(POTENTIAL_CLASS2);                                !

        NATOM = NATOM_NEW;

        NTYPE_ATOM = NTYPE_ATOM_NEW;

        allocate(CONFIG_RI(1:3,1:NATOM));

        allocate(CONFIG_NAT(1:NATOM));

        allocate(CONFIG_QI(1:NATOM));

        allocate(CONFIG_VI(1:3,1:NATOM));

        allocate(CONFIG_ATOMID(1:NATOM));

        allocate(CONFIG_ATOM_TYPE(1:NATOM));

        allocate(CONFIG_MOLECULEID(1:NATOM));

        allocate(ATOM_LABEL(1:NTYPE_ATOM));                                      ! UPDATE THE LIST OF LABEL OF ATOMS 
                                                                                         !
        allocate(ATOM_MASSE(1:NTYPE_ATOM));                                      ! UPDATE THE LIST OF MASSES 
                                                                                         !
        allocate(POTENTIAL_CLASS2(1:2,1:NTYPE_ATOM));                            !
                                                                                         !
        CONFIG_RI(1:3,1:NATOM)     = CONFIG_RI_NEW(1:3,1:NATOM);

        CONFIG_NAT(1:NATOM)        = CONFIG_NAT_NEW(1:NATOM);

        CONFIG_QI(1:NATOM)         = CONFIG_QI_NEW(1:NATOM);

        CONFIG_VI(1:3,1:NATOM)     = CONFIG_VI_NEW(1:3,1:NATOM);

        CONFIG_ATOMID(1:NATOM)     = CONFIG_ATOMID_NEW(1:NATOM);

        CONFIG_ATOM_TYPE(1:NATOM)  = CONFIG_ATOM_TYPE_NEW(1:NATOM);

        CONFIG_MOLECULEID(1:NATOM) = CONFIG_MOLECULEID_NEW(1:NATOM);

        ATOM_LABEL(1:NTYPE_ATOM)   = ATOM_LABEL_NEW(1:NTYPE_ATOM);

        ATOM_MASSE(1:NTYPE_ATOM)   = ATOM_MASSE_NEW(1:NTYPE_ATOM);

        POTENTIAL_CLASS2(1:2,1:NTYPE_ATOM) = POTENTIAL_CLASS2_NEW(1:2,1:NTYPE_ATOM);

        deallocate(CONFIG_RI_NEW);

        deallocate(CONFIG_NAT_NEW);

        deallocate(CONFIG_QI_NEW);

        deallocate(CONFIG_VI_NEW);

        deallocate(CONFIG_ATOMID_NEW);

        deallocate(CONFIG_ATOM_TYPE_NEW);

        deallocate(CONFIG_MOLECULEID_NEW);

        deallocate(ATOM_LABEL_NEW);                                      ! UPDATE THE LIST OF LABEL OF ATOMS 
                                                                                         !
        deallocate(ATOM_MASSE_NEW);                                      ! UPDATE THE LIST OF MASSES 
                                                                                         !
        deallocate(POTENTIAL_CLASS2_NEW);                            !

!   end if


    write(icanal,'(a19,i8,a43)') '| NATOM          : ', NATOM,          REPEAT(' ',42)//'|';
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
    write(icanal,'(a19,i8,a43)') '| NTYPE_ATOM     : ', NTYPE_ATOM,     REPEAT(' ',42)//'|';
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';

end subroutine RANDOM_INSERTION_DUMMY_PARTICLE











