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

subroutine ALLOCATE_CONFIG_ARRAYS(icanal)

!   ************************************************************************************************
!   **                                  READ THE INPUT FILE                                       **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                    : CANAL ON WHICH OUTPUT DATA ARE WRITTEN                         **
!   **                                                                                            **
!   ************************************************************************************************

    use module_size_arrays;

    use module_library;

    use module_config;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

!   ************************************************************************************************

    integer (kind=4) :: i, j, k, m;

!   ************************************************************************************************

    integer (kind=4) :: IOSEF1, IOSEF2, IOSEF3, IOSEF4, IOSEF5;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the routine #############################################################
                                                                                         !
    CHTITLE = 'Allocate config arrays';                                                  !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialization of config properties ########################################################
                                                                                         !
    NATOM     = 0;                                                                       !
                                                                                         !
    NBOND     = 0;                                                                       !
                                                                                         !
    NANGLE    = 0;                                                                       !
                                                                                         !
    NDIHEDRAL = 0;                                                                       !
                                                                                         !
    NIMPROPER = 0;                                                                       !
                                                                                         !
    NTYPE_ATOM     = 0;

    NTYPE_BOND     = 0;

    NTYPE_ANGLE    = 0;

    NTYPE_DIHEDRAL = 0;

    NTYPE_IMPROPER = 0;
                                                                                         !
    NPAIR_COEFF_CROSS = 0;                                                               !
                                                                                         !
    write(icanal,'(a70)') '| Config properties were initialized'//REPEAT(' ',33)// '|';  !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Get the number of library files to consider ################################################
                                                                                         !
!   if ( NFILE_LIBRARY > 0 )  IOSEF1 = NFILE_LIBRARY;                                    !
                                                                                         !
!   if ( NTotal_monomers > 0 ) IOSEF1 = NTotal_monomers;                                 !
                                                                                         !
!   ### Build simulation properties, when inserting molecules from a library #######################
                                                                                         !
    if ( NFILE_LIBRARY > 0 ) then;                                                       !
                                                                                         !
        do i = 1, NFILE_LIBRARY;                                                         ! LOOP OVER THE MOLECULE TYPE READ IN THE LIBRARY 
                                                                                         ! 
            NATOM = NATOM               + &                                              !
                    NATOM_LIBRARY(i)    * &                                              !
                    NMOLECULE_INSERTION(i);                                              !
                                                                                         !
            if ( NBOND_LIBRARY(i) > 0 ) then;                                            !
                                                                                         !
                NBOND = NBOND + NBOND_LIBRARY(i) *      &                                !
                                NMOLECULE_INSERTION(i);                                  !
                                                                                         !
            end if                                                                       !
                                                                                         !
            if ( NANGLE_LIBRARY(i)    > 0 ) NANGLE    = NANGLE    + NANGLE_LIBRARY(i)    * &
                                                        NMOLECULE_INSERTION(i);

            if ( NDIHEDRAL_LIBRARY(i) > 0 ) NDIHEDRAL = NDIHEDRAL + NDIHEDRAL_LIBRARY(i) * &
                                                        NMOLECULE_INSERTION(i);

            if ( NIMPROPER_LIBRARY(i) > 0 ) NIMPROPER = NIMPROPER + NIMPROPER_LIBRARY(i) * &
                                                        NMOLECULE_INSERTION(i);

            NTYPE_ATOM     = NTYPE_ATOM     + NTYPE_ATOM_LIBRARY(i);

            if ( NTYPE_BOND_LIBRARY(i)     > 0 ) NTYPE_BOND     = NTYPE_BOND     + &
                                                                  NTYPE_BOND_LIBRARY(i);

            if ( NTYPE_ANGLE_LIBRARY(i)    > 0 ) NTYPE_ANGLE    = NTYPE_ANGLE    + &
                                                                  NTYPE_ANGLE_LIBRARY(i);

            if ( NTYPE_DIHEDRAL_LIBRARY(i) > 0 ) NTYPE_DIHEDRAL = NTYPE_DIHEDRAL + &
                                                                  NTYPE_DIHEDRAL_LIBRARY(i);
                                                                                         !
            if ( NTYPE_IMPROPER_LIBRARY(i) > 0 ) then;                                   !
                                                                                         !
                NTYPE_IMPROPER = NTYPE_IMPROPER +           &                            !
                                 NTYPE_IMPROPER_LIBRARY(i);                              !
                                                                                         !
            end if                                                                       !
                                                                                         !
            if ( NPAIR_COEFF_CROSS_LIBRARY(i) > 0 ) then;                                !
                                                                                         !
                NPAIR_COEFF_CROSS = NPAIR_COEFF_CROSS + NPAIR_COEFF_CROSS_LIBRARY(i);    !
                                                                                         !
            end if                                                                       !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Build configuration properties, when inserting particles ###################################

!   if ( NTYPE_DUMMY_PARTICLE > 0 ) then;

!       do i = 1, NTYPE_DUMMY_PARTICLE;

!           NATOM = NATOM + NDUMMY_PARTICLE(i);

!           NTYPE_ATOM = NTYPE_ATOM + 1;

!       end do

!   end if

!   ### Build configuration properties, when building polymers #####################################
                                                                                         !
    if ( NGenerate_polymer_species > 0 ) then;                                           !
                                                                                         !
!       do i = 1, NTotal_monomers;                                                       ! LOOP OVER THE MOLECULE TYPE READ IN THE LIBRARY 
                                                                                         !
        IOSEF1 = 0;                                                                      !

        do i = 1, NGenerate_polymer_species;

            do j = 1, NSPECIES_PER_POLYMER_CHAINS(i);
 
                IOSEF1 = IOSEF1 + 1;

                NATOM = NATOM + NATOM_LIBRARY(IOSEF1)          * &
                                NMONOMER_UNITS_PER_CHAINS(i,j) * &
                                NPOLYMER_CHAINS(i);                                      !
                                                                                         !
!               if ( NBOND_LIBRARY(i) > 0 ) NBOND = NBOND + NBOND_LIBRARY(i);            !
                                                                                         !
!               if ( NANGLE_LIBRARY(i) > 0 ) NANGLE = NANGLE + NANGLE_LIBRARY(i);            !
                                                                                         !
!               if ( NDIHEDRAL_LIBRARY(i) > 0 ) NDIHEDRAL = NDIHEDRAL + NDIHEDRAL_LIBRARY(i);!

!               if ( NIMPROPER_LIBRARY(i) > 0 ) NIMPROPER = NIMPROPER + NIMPROPER_LIBRARY(i);!

                if ( NTYPE_ATOM_LIBRARY(i) > 0 ) then;

                    NTYPE_ATOM = NTYPE_ATOM + NTYPE_ATOM_LIBRARY(i);

                end if

!               if ( NTYPE_BOND_LIBRARY(i) > 0 ) NTYPE_BOND = NTYPE_BOND + NTYPE_BOND_LIBRARY(i);

!               if ( NTYPE_ANGLE_LIBRARY(i) > 0 ) NTYPE_ANGLE = NTYPE_ANGLE + NTYPE_ANGLE_LIBRARY(i);

!               if ( NTYPE_DIHEDRAL_LIBRARY(i) > 0 ) NTYPE_DIHEDRAL = NTYPE_DIHEDRAL + &
!                                                                     NTYPE_DIHEDRAL_LIBRARY(i);

!               if ( NTYPE_IMPROPER_LIBRARY(i) > 0 ) NTYPE_IMPROPER = NTYPE_IMPROPER + &
!                                                                     NTYPE_IMPROPER_LIBRARY(i);

            end do

        end do

    end if

!   ### Writing configuration properties ###########################################################
                                                                                         !
    write(icanal,'(a19,i8,a43)') '| NATOM          : ', &                                !
                                 NATOM,                 &                                !
                                 REPEAT(' ',42)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    if ( NBOND > 0 ) then;

        write(icanal,'(a19,i8,a43)') '| NBOND          : ', NBOND, REPEAT(' ',42)//'|';

        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|'; 

    end if

    if ( NANGLE > 0 ) then;

        write(icanal,'(a19,i8,a43)') '| NANGLE         : ', NANGLE, REPEAT(' ',42)//'|';

        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|'; 

    end if

    if ( NDIHEDRAL > 0 ) then;

        write(icanal,'(a19,i8,a43)') '| NDIHEDRAL      : ', NDIHEDRAL, REPEAT(' ',42)//'|';

        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|'; 

    end if

    if ( NIMPROPER > 0 ) then;

        write(icanal,'(a19,i8,a43)') '| NIMPROPER      : ', NIMPROPER, REPEAT(' ',42)//'|';

        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|'; 

    end if

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

    write(icanal,'(a19,i8,a43)') '| NTYPE_ATOM     : ', NTYPE_ATOM, REPEAT(' ',42)//'|';

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|'; 

    if ( NTYPE_BOND > 0 ) then;

        write(icanal,'(a19,i8,a43)') '| NTYPE_BOND     : ', NTYPE_BOND, REPEAT(' ',42)//'|';

        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|'; 

    end if

    if ( NTYPE_ANGLE > 0 ) then;

        write(icanal,'(a19,i8,a43)') '| NTYPE_ANGLE    : ', NTYPE_ANGLE, REPEAT(' ',42)//'|';

        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|'; 

    end if

    if ( NTYPE_DIHEDRAL > 0 ) then;

        write(icanal,'(a19,i8,a43)') '| NTYPE_DIHEDRAL : ', NTYPE_DIHEDRAL, REPEAT(' ',42)//'|';

        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|'; 

    end if

    if ( NTYPE_IMPROPER > 0 ) then; 

        write(icanal,'(a19,i8,a43)') '| NTYPE_IMPROPER : ', NTYPE_IMPROPER, REPEAT(' ',42)//'|';

        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|'; 

    end if
                                                                                         !
    if ( NPAIR_COEFF_CROSS > 0 ) then;                                                   !
                                                                                         !
        write(icanal,'(a22,i8,a40)') '| NPAIR_COEFF_CROSS : ', &                         !
                                     NPAIR_COEFF_CROSS,        &                         !
                                     REPEAT(' ',39)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Allocate

    if ( NATOM > 0 ) then;

        allocate(CONFIG_QI(1:NATOM));

        allocate(CONFIG_RI(1:3,1:NATOM));

        allocate(CONFIG_VI(1:3,1:NATOM));

        allocate(CONFIG_NAT(1:NATOM));

        allocate(CONFIG_ATOMID(1:NATOM));

        allocate(CONFIG_MOLECULEID(1:NATOM)); 

        allocate(CONFIG_ATOM_TYPE(1:NATOM));

        CONFIG_QI(1:NATOM)         = 0.0d0;

        CONFIG_RI(1:3,1:NATOM)     = 0.0d0;

        CONFIG_VI(1:3,1:NATOM)     = 0.0d0;

        CONFIG_NAT(1:NATOM)        = 'XXX';

        CONFIG_ATOMID(1:NATOM)     = 0;

        CONFIG_MOLECULEID(1:NATOM) = 0;

        CONFIG_ATOM_TYPE(1:NATOM)  = 0;

    end if 
                                                                                         !
!   ### Allocatation of arrays with the size of NTYPE_ATOM #########################################
                                                                                         !
    if ( NTYPE_ATOM > 0 ) then;                                                          !

        allocate(ATOM_LABEL(1:NTYPE_ATOM));

        allocate(ATOM_MASSE(1:NTYPE_ATOM));

        allocate(POTENTIAL_CLASS2(1:2,1:NTYPE_ATOM));

    end if
                                                                                         !
!   ### Initialization of arrays with the size of NTYPE_ATOM #######################################
                                                                                         !
    if ( NTYPE_ATOM > 0 ) then;

        ATOM_LABEL(1:NTYPE_ATOM)           = 'XXX';

        ATOM_MASSE(1:NTYPE_ATOM)           = 0.0d0;

        POTENTIAL_CLASS2(1:2,1:NTYPE_ATOM) = 0.0d0;

    end if

!   ### Allocation of arrays with the size of the number of pair coeff for cross interactions ######
                                                                                         !
    if ( NPAIR_COEFF_CROSS > 0 ) then;                                                   !
                                                                                         !
        allocate(PAIR_COEFF_CROSS(1:2,1:NPAIR_COEFF_CROSS));                             !
                                                                                         !
        allocate(PAIR_ATOMID_CROSS(1:2,1:NPAIR_COEFF_CROSS));                            !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Initialization of arrays with the size of the number of pair coeff for cross inter. ########
                                                                                         !
    if ( NPAIR_COEFF_CROSS > 0 ) then;                                                   !
                                                                                         ! 
        PAIR_COEFF_CROSS(1:2,1:NPAIR_COEFF_CROSS) = 0.0d0;                               !
                                                                                         !
        PAIR_ATOMID_CROSS(1:2,1:NPAIR_COEFF_CROSS) = 0;                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Allocate arrays with the size of the number of bonds #######################################
                                                                                         !
    if ( NBOND > 0 ) then;                                                               !

        allocate(BOND_ATOMID(1:2,1:NBOND));

        allocate(BOND_TYPE(1:NBOND));

        BOND_ATOMID(1:2,1:NBOND) = 0;

        BOND_TYPE(1:NBOND)       = 0;

    end if

    if ( NTYPE_BOND > 0 ) then;

        allocate(BOND_COEFFS(1:4,1:NTYPE_BOND));

        BOND_COEFFS(1:4,1:NTYPE_BOND) = 0.0d0;

    end if

    if ( NANGLE > 0 ) then;
 
        allocate(ANGLE_ATOMID(1:3,1:NANGLE));

        allocate(ANGLE_TYPE(1:NANGLE));

        ANGLE_ATOMID(1:3,1:NANGLE) = 0;

        ANGLE_TYPE(1:NANGLE)       = 0;

    end if

    if ( NTYPE_ANGLE > 0 ) then;

        allocate(ANGLE_COEFFS(1:4,1:NTYPE_ANGLE));

        allocate(BONDBOND_COEFFS(1:3,1:NTYPE_ANGLE));

        allocate(BONDANGLE_COEFFS(1:4,1:NTYPE_ANGLE));

        ANGLE_COEFFS(1:4,1:NTYPE_ANGLE)     = 0.0d0;

        BONDBOND_COEFFS(1:3,1:NTYPE_ANGLE)  = 0.0d0;

        BONDANGLE_COEFFS(1:4,1:NTYPE_ANGLE) = 0.0d0;

    end if

    if ( NDIHEDRAL > 0 ) then;

        allocate(DIHEDRAL_ATOMID(1:4,1:NDIHEDRAL));

        allocate(DIHEDRAL_TYPE(1:NDIHEDRAL));

        DIHEDRAL_ATOMID(1:4,1:NDIHEDRAL) = 0;

        DIHEDRAL_TYPE(1:NDIHEDRAL)    = 0;

    end if

    if ( NTYPE_DIHEDRAL > 0 ) then;

        allocate(DIHEDRAL_COEFFS(1:6,1:NTYPE_DIHEDRAL));

        allocate(MIDDLEBONDTORSION_COEFFS(1:4,1:NTYPE_DIHEDRAL));

        allocate(ENDBONDTORSION_COEFFS(1:8,1:NTYPE_DIHEDRAL));

        allocate(ANGLETORSION_COEFFS(1:8,1:NTYPE_DIHEDRAL));

        allocate(ANGLEANGLETORSION_COEFFS(1:3,1:NTYPE_DIHEDRAL));

        allocate(BONDBOND13_COEFFS(1:3,1:NTYPE_DIHEDRAL));

        DIHEDRAL_COEFFS(1:6,1:NTYPE_DIHEDRAL)          = 0.0d0;

        MIDDLEBONDTORSION_COEFFS(1:4,1:NTYPE_DIHEDRAL) = 0.0d0;

        ENDBONDTORSION_COEFFS(1:8,1:NTYPE_DIHEDRAL)    = 0.0d0;

        ANGLETORSION_COEFFS(1:8,1:NTYPE_DIHEDRAL)      = 0.0d0;

        ANGLEANGLETORSION_COEFFS(1:3,1:NTYPE_DIHEDRAL) = 0.0d0;

        BONDBOND13_COEFFS(1:3,1:NTYPE_DIHEDRAL)        = 0.0d0;

    end if

    if ( NIMPROPER > 0 ) then;

        allocate(IMPROPER_ATOMID(1:4,1:NIMPROPER));

        allocate(IMPROPER_TYPE(1:NIMPROPER));

        IMPROPER_ATOMID(1:4,1:NIMPROPER) = 0;

        IMPROPER_TYPE(1:NIMPROPER)       = 0;

    end if

    if ( NTYPE_IMPROPER > 0 ) then;                                                      !
                                                                                         !
        allocate(IMPROPER_COEFFS(1:2,1:NTYPE_IMPROPER));                                 !
                                                                                         !
        allocate(ANGLEANGLE_COEFFS(1:6,1:NTYPE_IMPROPER));                               !
                                                                                         !
        IMPROPER_COEFFS(1:2,1:NTYPE_IMPROPER)   = 0.0d0;                                 !
                                                                                         !
        ANGLEANGLE_COEFFS(1:6,1:NTYPE_IMPROPER) = 0.0d0;                                 !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Wrap up of the routine #####################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine ALLOCATE_CONFIG_ARRAYS











