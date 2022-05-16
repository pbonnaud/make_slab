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

module module_data_in

    implicit none;

!   ************************************************************************************************

    type position
        character (len=3) :: nat;
        logical :: wt;
        real (kind=8) :: x, y, z;
        real (kind=8) :: q;
    end type

    type PNTrAZ
       real (kind=8) :: C6, C8, C10, A, b;
       character (len=3) :: spec1, spec2;
    end type

    type molecule
        integer (kind=4) :: MULTI, LABEL;
        real (kind=8) :: sigma, eps, q, M;
        real (kind=8), dimension(1:3) :: RI;
        character (len=3) :: nat, CHDISPREP; 
    end type

    type coord
!       SEQUENCE;
        integer (kind=4), allocatable, dimension(:)  :: LABEL, IDILAT, INMOL, IDBOX;
        character (len=3), allocatable, dimension(:) :: NAT;
        real (kind=8), allocatable, dimension(:,:) :: RI, RG;

        complex (kind=8), allocatable, dimension(:,:,:) :: EIKRI;
    end type

    type substrat
        integer (kind=4) :: LABEL, IDBOX;
        character (len=3) :: NAT;
        real (kind=8), dimension(1:3) :: RI;
    end type

!   ************************************************************************************************

    integer (kind=4) :: RANDOM_GENERATOR_SEED;

    real (kind=8) :: TARGET_TEMPK;

!   ### Set arrays for properties of atoms in the Mendeleiev periodic table of elements ############
                                                                                         !
    real (kind=8), dimension(1:200) :: MOLAR_MASS_ELEMENT, &                             !
                                       ELEMENT_RADIUS;                                   !
                                                                                         !
    character (len=150), dimension(1:200) :: ELEMENT_FULL_NAME,   &                      !
                                             ELEMENT_SYMBOL_NAME;                        !
                                                                                         !
!   ### Set arrays for bond properties #############################################################
                                                                                         !
    real (kind=8), dimension(1:200,1:200) :: ELEMENT_COVALENT_BONDS;                     ! 
                                                                                         !
!   ### Arrays for the building method of the molecular configuration ##############################
                                                                                         !
    integer (kind=4) :: IBUILD_METHOD;                                                   !
                                                                                         !
!   ### Input file parameters ######################################################################

    integer (kind=4) :: NMOLECULAR_CONFIG;

    integer (kind=4) :: iduplicate;

    integer (kind=4) :: iregions, iworking_file_insertion, iuse_moleculeid;

    integer (kind=4) :: igenerate_lammps_input;

    integer (kind=4) :: NREGIONS_INSERTION;

    integer (kind=4), dimension(:), allocatable :: REGION_AXIS;

    real (kind=8), dimension(:,:), allocatable :: REGION_BOUNDS;

    integer (kind=4) :: NSWITCH_ATOMIDS,          &
                        NSWITCH_BONDIDS,          &
                        NSWITCH_ANGLEIDS,         &
                        NSWITCH_DIHEDRALIDS,      &
                        NTRANSLATION_MOLECULEIDS;

    integer (kind=4), dimension(1:1000) :: TRANSLATION_MOLECULEIDS;

    integer (kind=4), dimension(1:2,1:1000) :: SWITCH_ATOMIDS,     &
                                               SWITCH_BONDIDS,     &
                                               SWITCH_ANGLEIDS,    &
                                               SWITCH_DIHEDRALIDS;

    real (kind=8), dimension(1:3,1:1000) :: TRANSLATION_VECTOR;

    character (len=250) :: SWITCH_FILE_NAME,      &
                           TRANSLATION_FILE_NAME;


!   ############    
!   ############    
!   ############    
!   ############    
!   ############    

    integer (kind=4) :: imodif, ipsub, ifsub, idople, iboxc, IROTA, ITRAN;

    integer (kind=4) :: ch_amorph;

    real (kind=8) ::  lmax;

    integer (kind=4) :: nbins;

    integer (kind=4) :: IKEEPGATH, ICENTROSYM, IOXYTETRA, IPLATELET;

    integer (kind=4) :: IUMBRELLA, NLUMBRELLA;

    real (kind=8) :: DELTAZ_KEEPGATH;

    character (len=3) :: CHFSUB;

    integer (kind=4) :: NLIST_SPECIES_MODIF;

    integer (kind=4), dimension(1:10) :: LIST_SPECIES_MODIF;

    real (kind=4) :: DRI_TRANS;

!   ### READ DATA PARAMETERS #######################################################################

    integer (kind=4) :: chslab;

    integer (kind=4) :: espece, NSLAB, NSLAB_FINAL, NSLAB_PLATE;

    integer (kind=4), dimension(1:10) :: nsite, N, NATSP, FIXSITE, NSITE_TOT;

    integer (kind=4), dimension(1:10) :: NATSP_FINAL, NATSP_PLATE;

    type (molecule), dimension(1:100,1:100) :: mol;

    integer (kind=4) :: NLABEL;

    real (kind=8), dimension(1:10) :: ATOM_SIG, ATOM_EPS, ATOM_CHAR, ATOM_MASS;

!   ### READ FIELD PARAMETERS ######################################################################
    integer (kind=4) :: ICOMA;

    real (kind=8), dimension(1:3,1:10) :: COMA_COORD, COMA_COORD_FINAL;

!   ### READ CONFIGURATION PARAMETERS ##############################################################
    integer (kind=4) :: NTOT_ATOMS, NTOT_ATOMS_FINAL;

    type (coord), allocatable, dimension(:) :: DATA, DATA_FINAL, DATA_PLATE;

    type (substrat), allocatable, dimension(:) :: SUB, SUB_FINAL, SUB_PLATE;

    integer (kind=4) :: NATOM_FINAL, NATOM_TMP, NSLAB_TMP;

    integer (kind=4), dimension(1:10) :: NATSP_TMP;

    type (coord), allocatable, dimension(:) :: DATA_TMP;

    type (substrat), allocatable, dimension(:) :: SUB_TMP;

!   ### STARTING CELL PARAMETERS ###################################################################

    integer (kind=4) :: NMAX;

    integer (kind=4) :: NVDWIWS;

    integer (kind=4) :: NMOTIF, NMAILLE, NBCAT;
    integer (kind=4) :: NMULTI, NWORK, NWORK_EXTD;

    integer (kind=4) :: NBNATIN;

    integer (kind=4), dimension(1:2) :: NBINIT, NBFINAL, NBDOPLE; 

    integer (kind=4), dimension(1:20) :: NBATLIST;

    real (kind=8), dimension(1:20) :: QLIST;

    character (len=3), dimension(1:20) :: NATLIST; 

    type (position), dimension(1:30) :: motif;

    integer (kind=4), dimension(1:3) :: MULTI;

    integer (kind=4), dimension(1:100) :: NMOL;

    integer (kind=4) :: nisp, nsusp;

    type(PNTrAZ), dimension(1:100) :: parsub;

    character (len=3), allocatable, dimension(:,:) :: NATINIT, NATFINAL, NATDOPLE;

    integer (kind=4), allocatable, dimension(:,:) :: ADSTATUS;

    real (kind=8), allocatable, dimension(:,:,:) :: RINIT, RFINAL, RDOPLE;

    real (kind=8), dimension(1:3) :: DOPLE_ANGDEG, DOPLE_ANGRAD, DOPLE_TRANS;

!   ### FINAL CELL PARAMETERS ######################################################################

    integer (kind=4) :: NESPECE, NVDWI, NSTOT_WALL;

    integer (kind=4) :: ILJONES, IPNTRAZ;

    integer (kind=4) :: FLAG_BONDS, FLAG_ANGLES;

    integer (kind=4), dimension(1:20) :: ISTRUCT, IBIAS, IBONDS, IANGLES;

    real (kind=8), dimension(1:3,1:10) :: RCDM_MOL, ORIENT_MOL;

    integer (kind=4), dimension(1:10) :: ICDM_MOL;

    real (kind=8), dimension(1:20,1:20) :: MASSMO, BONDS_KR, ANGLES_KR, BONDS_R0, ANGLES_THETA0;

    integer (kind=4), dimension(1:20,1:20) :: MULAT;

    real (kind=8), dimension(1:20) :: MTOT;

    real (kind=8), dimension(1:10,1:10) :: SIGATO, EPSATO, CHARATO;

    real (kind=8), dimension(1:10,1:10) :: AIJATO, BIJATO, C06ATO, C08ATO, C10ATO;

!   ### CELL PARAMETERS  ###########################################################################

    real (kind=8) :: VOL;

    real (kind=8), dimension(1:3) :: HPORE;

end module module_data_in
