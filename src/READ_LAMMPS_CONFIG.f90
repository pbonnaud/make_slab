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


subroutine READ_LAMMPS_CONFIG(icanal,CHEXT,MATA,FLAG_SORT_MOLEC); 

!   ************************************************************************************************
!   **                                 READ LAMMPS FILES                                          **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal          : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                                 **
!   **                                                                                            **
!   ** CHEXT           : NAME OF THE LAMMPS CONFIGURATION                                         **
!   **                                                                                            **
!   ** MATA            : PASSAGE MATRIX                                                           **
!   **                                                                                            **
!   ** FLAG_SORT_MOLEC : FLAG TO CHECK IF MOLECULES AND ATOMS NEED TO BE SORTED                   ** 
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

    use module_config;

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

!   include 'mpif.h'

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    character (len=150), intent(in) :: CHEXT;

!   ************************************************************************************************

    real (kind=8), dimension(1:6), intent(out) :: MATA;

    integer (kind=4), intent(out) :: FLAG_SORT_MOLEC;

!   ************************************************************************************************

    integer (kind=4) :: ierr;

!   ************************************************************************************************

    integer (kind=4) :: i, j;

    integer (kind=4) :: ILENGTH_CHEXT; !, ATOMS_FLAG;

    integer (kind=4) :: ixy, ixz, iyz;

    integer (kind=4) :: imcon, ICONNECT, IATOM_TYPE;

    integer (kind=4) :: IFOUND;

    real (kind=8) :: XLO, XHI, YLO, YHI, ZLO, ZHI;

    character (len=150) :: CHARLINE;

    integer (kind=4) :: EOF, EOF2;

    real (kind=8), dimension(1:3) :: TAB_ROSEF; 

    character (len=250) :: CHAIN_LENGTH;

!   ************************************************************************************************

    real (kind=8) :: t_start, t_stop;

    real (kind=8) :: tps1, tps2, tps3;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the routine #############################################################
                                                                                         !
    CHTITLE = 'Read LAMMPS configuration';                                               !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   ### Get the time at the beginning of the routine ###############################################
                                                                                         !
    call CPU_TIME(t_start);                                                              !
                                                                                         !
!   ### Write the name of the configuration file ###################################################
                                                                                         !
    ILENGTH_CHEXT = LEN_TRIM(CHEXT);                                                     !
                                                                                         ! 
    IOSEF1 = 67 - ILENGTH_CHEXT;                                                         !
                                                                                         !
    write(icanal,'(a70)') '| '//TRIM(CHEXT)//REPEAT(' ',IOSEF1)//'|';                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialisation of variables ################################################################
                                                                                         !
    NATOM     = 0;                                                                       !
                                                                                         !
    NBOND     = 0;

    NANGLE    = 0;

    NDIHEDRAL = 0;

    NIMPROPER = 0;

    NTYPE_ATOM     = 0;

    NTYPE_BOND     = 0;

    NTYPE_ANGLE    = 0; 

    NTYPE_DIHEDRAL = 0;

    NTYPE_IMPROPER = 0; 

    XLO = 0.0d0;

    XHI = 0.0d0;

    YLO = 0.0d0;

    YHI = 0.0d0;

    ZLO = 0.0d0;

    ZHI = 0.0d0;

    MATA(1:6) = 0.0d0;

    FLAG_SORT_MOLEC = 0;

    IPOTENTIAL_CLASS2 = 0;

!   ### Read the configuration file ################################################################
                                                                                         !
    open(2,file=TRIM(CHEXT));                                                            ! 
                                                                                         !
    EOF = 0;                                                                             !
                                                                                         !
    do                                                                                   !
                                                                                         !
        read(2,'(a)',iostat=EOF) CHAIN_LENGTH;                                           !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
        if ( EOF /= 0 ) EXIT;                                                            !
                                                                                         !
        if ( INDEX(CHAIN_LENGTH,'atoms') > 0 ) then;                                     !
                                                                                         !
            read(CHAIN_LENGTH,*) NATOM;                                                  !
                                                                                         !
            write(icanal,'(a19,i8,a43)') '| NATOM          : ', &                        !
                                         NATOM,                 &                        !
                                         REPEAT(' ',42)//'|';                            !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            if ( NATOM > 0 ) then;                                                       !
                                                                                         !
                allocate(CONFIG_ATOMID(1:NATOM));                                        !
                                                                                         !
                allocate(CONFIG_MOLECULEID(1:NATOM));                                    !
                                                                                         !
                allocate(CONFIG_ATOM_TYPE(1:NATOM));                                     !
                                                                                         !
                allocate(CONFIG_QI(1:NATOM));                                            !
                                                                                         !
                allocate(CONFIG_VI(1:3,1:NATOM));                                        !
                                                                                         !
                allocate(CONFIG_RI(1:3,1:NATOM));                                        !
                                                                                         !
                allocate(CONFIG_NAT(1:NATOM));                                           !
                                                                                         !
                CONFIG_ATOMID(1:NATOM)     = 0;                                          !
                                                                                         !
                CONFIG_MOLECULEID(1:NATOM) = 0;                                          !
                                                                                         !
                CONFIG_ATOM_TYPE(1:NATOM)  = 0;                                          !
                                                                                         !
                CONFIG_QI(1:NATOM)         = 0.0d0;                                      !
                                                                                         !
                CONFIG_VI(1:3,1:NATOM)     = 0.0d0;                                      !
                                                                                         ! 
                CONFIG_RI(1:3,1:NATOM)     = 0.0d0;                                      !
                                                                                         !
                CONFIG_NAT(1:NATOM)        = 'XXX';                                      !
                                                                                         !
            end if                                                                       !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        else if ( INDEX(CHAIN_LENGTH,'atom types') > 0 ) then;                           !
                                                                                         !
            read(CHAIN_LENGTH,*) NTYPE_ATOM;                                             !
                                                                                         !
            write(icanal,'(a19,i8,a43)') '| NTYPE_ATOM     : ', &                        !
                                         NTYPE_ATOM,            &                        !
                                         REPEAT(' ',42)//'|';                            !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            if ( NTYPE_ATOM > 0 ) then;                                                  !
                                                                                         !
                allocate(ATOM_MASSE(1:NTYPE_ATOM));                                      !
                                                                                         !
                allocate(ATOM_LABEL(1:NTYPE_ATOM));                                      !
                                                                                         !
                allocate(POTENTIAL_CLASS2(1:2,1:NTYPE_ATOM));                            !
                                                                                         ! 
                ATOM_MASSE(1:NTYPE_ATOM)  = 0.0d0;                                       !
                                                                                         !
                ATOM_LABEL(1:NTYPE_ATOM)  = 'XXX';                                       !
                                                                                         !
                POTENTIAL_CLASS2(1:2,1:NTYPE_ATOM) = 0.0d0;                              !
                                                                                         !
            end if                                                                       !
                                                                                         !
            else if ( INDEX(CHAIN_LENGTH,'bonds') > 0 ) then;                            !
                                                                                         !
                read(CHAIN_LENGTH,*) NBOND;                                              !
                                                                                         !
                write(icanal,'(a19,i8,a43)') '| NBOND          : ', &                    !
                                             NBOND,                 &                    !
                                             REPEAT(' ',42)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                if ( NBOND > 0 ) then;                                                   !
                                                                                         !
                    allocate(BOND_TYPE(1:NBOND));                                        !
                                                                                         !
                    allocate(BOND_ATOMID(1:2,1:NBOND));                                  !
                                                                                         !
                    BOND_TYPE(1:NBOND)       = 0;                                        !
                                                                                         !
                    BOND_ATOMID(1:2,1:NBOND) = 0;                                        !
                                                                                         !
                end if                                                                   !
                                                                                         !
            else if ( INDEX(CHAIN_LENGTH,'bond types') > 0 ) then;                       !
                                                                                         !
                read(CHAIN_LENGTH,*) NTYPE_BOND;                                         !
                                                                                         !
                write(icanal,'(a19,i8,a43)') '| NTYPE_BOND     : ', &                    !
                                             NTYPE_BOND,            &                    !
                                             REPEAT(' ',42)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                if ( NTYPE_BOND > 0 ) then;                                              !
                                                                                         !
                    allocate(BOND_COEFFS(1:4,1:NTYPE_BOND));                             !
                                                                                         ! 
                    BOND_COEFFS(1:4,1:NTYPE_BOND) = 0.0d0;                               !
                                                                                         !
                end if                                                                   !
                                                                                         !
            else if ( INDEX(CHAIN_LENGTH,'angles') > 0 ) then;                           !
                                                                                         !
                read(CHAIN_LENGTH,*) NANGLE;                                             !
                                                                                         !
                write(icanal,'(a19,i8,a43)') '| NANGLE         : ', &                    !
                                             NANGLE,                &                    !
                                             REPEAT(' ',42)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                allocate(ANGLE_TYPE(1:NANGLE));                                          !
                                                                                         !
                allocate(ANGLE_ATOMID(1:3,1:NANGLE));                                    !
                                                                                         !
                ANGLE_TYPE(1:NANGLE)       = 0;                                          !
                                                                                         !
                ANGLE_ATOMID(1:3,1:NANGLE) = 0;                                          !
                                                                                         !
            else if ( INDEX(CHAIN_LENGTH,'angle types') > 0 ) then;                      !
                                                                                         !
                read(CHAIN_LENGTH,*) NTYPE_ANGLE;                                        !
                                                                                         !
                write(icanal,'(a19,i8,a43)') '| NTYPE_ANGLE    : ', &                    !
                                             NTYPE_ANGLE,           &                    !
                                             REPEAT(' ',42)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                if ( NTYPE_ANGLE > 0 ) then;                                             !
                                                                                         !
                    allocate(ANGLE_COEFFS(1:4,1:NTYPE_ANGLE));                           !
                                                                                         !
                    allocate(BONDBOND_COEFFS(1:3,1:NTYPE_ANGLE));                        !
                                                                                         !
                    allocate(BONDANGLE_COEFFS(1:4,1:NTYPE_ANGLE));                       !
                                                                                         !
                    ANGLE_COEFFS(1:4,1:NTYPE_ANGLE)     = 0.0d0;                         !
                                                                                         !
                    BONDBOND_COEFFS(1:3,1:NTYPE_ANGLE)  = 0.0d0;                         !
                                                                                         !
                    BONDANGLE_COEFFS(1:4,1:NTYPE_ANGLE) = 0.0d0;                         !
                                                                                         !
                end if                                                                   !
                                                                                         !
            else if ( INDEX(CHAIN_LENGTH,'dihedrals') > 0 ) then;                        !
                                                                                         !
                read(CHAIN_LENGTH,*) NDIHEDRAL;                                          !
                                                                                         !
                write(icanal,'(a19,i8,a43)') '| NDIHEDRAL      : ', &                    !
                                             NDIHEDRAL,             &                    ! 
                                             REPEAT(' ',42)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                if ( NDIHEDRAL > 0 ) then;                                               !
                                                                                         !
                    allocate(DIHEDRAL_TYPE(1:NDIHEDRAL));                                !
                                                                                         !
                    allocate(DIHEDRAL_ATOMID(1:4,1:NDIHEDRAL));                          !
                                                                                         !
                    DIHEDRAL_TYPE(1:NDIHEDRAL)       = 0;                                !
                                                                                         !
                    DIHEDRAL_ATOMID(1:4,1:NDIHEDRAL) = 0;                                !
                                                                                         !
                end if                                                                   !
                                                                                         !
            else if ( INDEX(CHAIN_LENGTH,'dihedral types') > 0 ) then;                   !
                                                                                         ! 
                read(CHAIN_LENGTH,*) NTYPE_DIHEDRAL;                                     !
                                                                                         !
                write(icanal,'(a19,i8,a43)') '| NTYPE_DIHEDRAL : ', &                    ! 
                                             NTYPE_DIHEDRAL,        &                    !
                                             REPEAT(' ',42)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                if ( NTYPE_DIHEDRAL > 0 ) then;

                    allocate(DIHEDRAL_COEFFS(1:6,1:NTYPE_DIHEDRAL));                         !

                    allocate(MIDDLEBONDTORSION_COEFFS(1:4,1:NTYPE_DIHEDRAL));                !

                    allocate(ENDBONDTORSION_COEFFS(1:8,1:NTYPE_DIHEDRAL));                   !

                    allocate(ANGLETORSION_COEFFS(1:8,1:NTYPE_DIHEDRAL));                     !

                    allocate(ANGLEANGLETORSION_COEFFS(1:3,1:NTYPE_DIHEDRAL));                !

                    allocate(BONDBOND13_COEFFS(1:3,1:NTYPE_DIHEDRAL));                       !
                                                                                         !
                    DIHEDRAL_COEFFS(1:6,1:NTYPE_DIHEDRAL)          = 0.0d0;                  !

                    MIDDLEBONDTORSION_COEFFS(1:4,1:NTYPE_DIHEDRAL) = 0.0d0;                  !

                    ENDBONDTORSION_COEFFS(1:8,1:NTYPE_DIHEDRAL)    = 0.0d0;                  !

                    ANGLETORSION_COEFFS(1:8,1:NTYPE_DIHEDRAL)      = 0.0d0;                  !

                    ANGLEANGLETORSION_COEFFS(1:3,1:NTYPE_DIHEDRAL) = 0.0d0;                  !

                    BONDBOND13_COEFFS(1:3,1:NTYPE_DIHEDRAL)        = 0.0d0;                  !

                end if
                                                                                         !
            else if ( INDEX(CHAIN_LENGTH,'impropers') > 0 ) then;                        !
                                                                                         ! 
                read(CHAIN_LENGTH,*) NIMPROPER;                                          !
                                                                                         !
                write(icanal,'(a19,i8,a43)') '| NIMPROPER      : ', &                    !
                                             NIMPROPER,             &                    !
                                             REPEAT(' ',42)//'|';                        !

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
                                                                                         !
                if ( NIMPROPER > 0 ) then;

                    allocate(IMPROPER_TYPE(1:NIMPROPER));                                    !

                    allocate(IMPROPER_ATOMID(1:4,1:NIMPROPER));                              !
                                                                                         !
                    IMPROPER_TYPE(1:NIMPROPER)       = 0;                                    !

                    IMPROPER_ATOMID(1:4,1:NIMPROPER) = 0;                                    !

                end if
                                                                                         !
            else if ( INDEX(CHAIN_LENGTH,'improper types') > 0 ) then;                   !
                                                                                         !
                read(CHAIN_LENGTH,*) NTYPE_IMPROPER;                                     !
                                                                                         !
                write(icanal,'(a19,i8,a43)') '| NTYPE_IMPROPER : ', &                    !
                                             NTYPE_IMPROPER,        &                    !
                                             REPEAT(' ',42)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

                if ( NTYPE_IMPROPER > 0 ) then;

                    allocate(IMPROPER_COEFFS(1:2,1:NTYPE_IMPROPER));                         !

                    allocate(ANGLEANGLE_COEFFS(1:6,1:NTYPE_IMPROPER));                       !
                                                                                         !
                    IMPROPER_COEFFS(1:2,1:NTYPE_IMPROPER)   = 0.0d0;                         !

                    ANGLEANGLE_COEFFS(1:6,1:NTYPE_IMPROPER) = 0.0d0;                         !

                end if
                                                                                         !
            else if ( INDEX(CHAIN_LENGTH,'xlo xhi') > 0 ) then;                          !
                                                                                         !
                read(CHAIN_LENGTH,*) XLO, XHI;                                           !
                                                                                         !
                MATA(1) = XHI - XLO;                                                     !
                                                                                         !
                write(icanal,'(a16,f15.6,a39)') '| MATA(1) [A] : ', &                    !
                                                MATA(1),            &                    !
                                                REPEAT(' ',38)//'|';                     !

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
                                                                                         !
            else if ( INDEX(CHAIN_LENGTH,'ylo yhi') > 0 ) then;                          !
                                                                                         !
                read(CHAIN_LENGTH,*) YLO, YHI;                                           !
                                                                                         !
                MATA(2) = YHI - YLO;                                                     !
                                                                                         !
                write(icanal,'(a16,f15.6,a39)') '| MATA(2) [A] : ', &                    !
                                                MATA(2),            &                    !
                                                REPEAT(' ',38)//'|';                     !

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
                                                                                         !
            else if ( INDEX(CHAIN_LENGTH,'zlo zhi') > 0 ) then;                          !
                                                                                         !
                read(CHAIN_LENGTH,*) ZLO, ZHI;                                           !
                                                                                         !
                MATA(3) = ZHI - ZLO;                                                     !
                                                                                         !
                write(icanal,'(a16,f15.6,a39)') '| MATA(3) [A] : ', &                    !
                                                MATA(3),            &                    !
                                                REPEAT(' ',38)//'|';                     !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            else if ( INDEX(CHAIN_LENGTH,'xy xz yz') > 0 ) then;                         !
                                                                                         !
                read(CHAIN_LENGTH,*) MATA(6), MATA(5), MATA(4);                          !
                                                                                         !
                write(icanal,'(a16,f15.6,a39)') '| XY TILT [A] : ', &                    !
                                                MATA(6),            &                    !
                                                REPEAT(' ',38)//'|';                     !
 
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

                write(icanal,'(a16,f15.6,a39)') '| XZ TILT [A] : ', &                    !
                                                MATA(5),            &                    !
                                                REPEAT(' ',38)//'|';                     !

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

                write(icanal,'(a16,f15.6,a39)') '| YZ TILT [A] : ', &                    !
                                                MATA(4),            &                    !
                                                REPEAT(' ',38)//'|';                     !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            else if ( INDEX(CHAIN_LENGTH,'Masses') > 0 ) then;

                read(2,*);

                do i = 1, NTYPE_ATOM;

                   read(2,*) IOSEF1, ATOM_MASSE(IOSEF1);

                end do

                do i = 1, NTYPE_ATOM;

                    do j = 1, 200;                                                       !
                                                                                         !
                        ROSEF1 = ABS( ATOM_MASSE(i) - MOLAR_MASS_ELEMENT(j) );           !
                                                                                         !
                        if ( ROSEF1 < 0.1d0 ) ATOM_LABEL(i) = ELEMENT_SYMBOL_NAME(j);    !
                                                                                         !
                    end do                                                               !
                                                                                         !
                end do                                                                   !
                                                                                         !
                write(icanal,'(a70)') '| Masses were read'//REPEAT(' ',51)//'|';         !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            else if ( INDEX(CHAIN_LENGTH,'Pair Coeffs') > 0 ) then;                      !
                                                                                         !
                if ( INDEX(CHAIN_LENGTH,'lj/class2/coul/long') > 0 ) then;               !
                                                                                         !
                    IPOTENTIAL_CLASS2 = 1;                                               !
                                                                                         !
                    IOSEF1 = INDEX(CHAIN_LENGTH,'lj/class2/coul/long');                  !
                                                                                         !
                    IOSEF2 = IOSEF1 + 19;                                                !
                                                                                         !
                else if ( INDEX(CHAIN_LENGTH,'lj/cut/coul/long') > 0 ) then;             !
                                                                                         !
                    IPOTENTIAL_CLASS2 = 2;                                               !
                                                                                         !
                    IOSEF1 = INDEX(CHAIN_LENGTH,'lj/cut/coul/long');                     !
                                                                                         !
                    IOSEF2 = IOSEF1 + 16;                                                !
                                                                                         !
                end if                                                                   !
                                                                                         !
                CH_PAIR_STYLE = CHAIN_LENGTH(IOSEF1:IOSEF2);                             !
                                                                                         !
                read(2,*);                                                               !
                                                                                         !
                do i = 1, NTYPE_ATOM;                                                    !
                                                                                         !
                    read(2,*) IOSEF1, POTENTIAL_CLASS2(1:2,IOSEF1);                      !
                                                                                         !
                end do                                                                   !
                                                                                         !
                write(icanal,'(a70)') '| Pair Coeffs were read'//REPEAT(' ',46)//'|';    !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
!           else if ( INDEX(CHAIN_LENGTH,'Bond Coeffs # class2') > 0 ) then;             !
            else if ( INDEX(CHAIN_LENGTH,'Bond Coeffs') > 0 ) then;                      !
                                                                                         !
                if ( INDEX(CHAIN_LENGTH,'class2') > 0 ) CH_BOND_STYLE = 'class2';        !
                                                                                         !
                if ( INDEX(CHAIN_LENGTH,'harmonic') > 0 ) CH_BOND_STYLE = 'harmonic';    !
                                                                                         !
                NPARAM_BONDS = 1;                                                        !
                                                                                         !
                if ( TRIM(CH_BOND_STYLE) == 'class2' ) NPARAM_BONDS = 4;                 !
                                                                                         !
                if ( TRIM(CH_BOND_STYLE) == 'harmonic' ) NPARAM_BONDS = 2;               !
                                                                                         !
                read(2,*);                                                               !
                                                                                         !
                do i = 1, NTYPE_BOND;                                                    !
                                                                                         !
                    read(2,*) IOSEF1, BOND_COEFFS(1:NPARAM_BONDS,IOSEF1);                !
                                                                                         !
                end do                                                                   !
                                                                                         !
                write(icanal,'(a70)') '| Bond Coeffs were read'//REPEAT(' ',46)//'|';    !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            else if ( INDEX(CHAIN_LENGTH,'Angle Coeffs') > 0 ) then;                     ! 
                                                                                         !
                if ( INDEX(CHAIN_LENGTH,'class2') > 0 ) CH_ANGLE_STYLE = 'class2';       !
                                                                                         !
                if ( INDEX(CHAIN_LENGTH,'harmonic') > 0 ) CH_ANGLE_STYLE = 'harmonic';   !
                                                                                         !
                NPARAM_ANGLES = 1;                                                       !
                                                                                         !
                if ( TRIM(CH_ANGLE_STYLE) == 'class2' ) NPARAM_ANGLES = 4;               !
                                                                                         !
                if ( TRIM(CH_ANGLE_STYLE) == 'harmonic' ) NPARAM_ANGLES = 2;             !
                                                                                         !
                read(2,*);                                                               !
                                                                                         !
                do i = 1, NTYPE_ANGLE;                                                   !
                                                                                         !
                    read(2,*) IOSEF1, ANGLE_COEFFS(1:NPARAM_ANGLES,IOSEF1);              !
                                                                                         !
                end do                                                                   !
                                                                                         !
                write(icanal,'(a70)') '| Angle Coeffs were read'//REPEAT(' ',45)//'|';   !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            else if ( INDEX(CHAIN_LENGTH,'BondBond Coeffs') > 0 ) then;

                read(2,*);

                do i = 1, NTYPE_ANGLE;

                    read(2,*) IOSEF1, BONDBOND_COEFFS(1:3,IOSEF1);

                end do

                write(icanal,'(a70)') '| BondBond Coeffs were read'//REPEAT(' ',42)//'|';

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

!               stop;

            else if ( INDEX(CHAIN_LENGTH,'BondAngle Coeffs') > 0 ) then;

                read(2,*);

                do i = 1, NTYPE_ANGLE;

                    read(2,*) IOSEF1, BONDANGLE_COEFFS(1:4,IOSEF1);

                end do

                write(icanal,'(a70)') '| BondAngle Coeffs were read'//REPEAT(' ',41)//'|';

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

!               stop;

            else if ( INDEX(CHAIN_LENGTH,'Dihedral Coeffs') > 0 ) then;                  !
                                                                                         !
                if ( INDEX(CHAIN_LENGTH,'opls') > 0 ) CH_DIHEDRAL_STYLE = 'opls';        !
                                                                                         !
                if ( INDEX(CHAIN_LENGTH,'class2') > 0 ) CH_DIHEDRAL_STYLE = 'class2';    !
                                                                                         !
                NPARAM_DIHEDRALS = 1;                                                    !
                                                                                         !
                if ( TRIM(CH_DIHEDRAL_STYLE) == 'opls' ) NPARAM_DIHEDRALS = 4;           !
                                                                                         !
                if ( TRIM(CH_DIHEDRAL_STYLE) == 'class2' ) NPARAM_DIHEDRALS = 6;         !
                                                                                         !
                read(2,*);                                                               !
                                                                                         !
                do i = 1, NTYPE_DIHEDRAL;                                                !
                                                                                         !
                    read(2,*) IOSEF1, DIHEDRAL_COEFFS(1:NPARAM_DIHEDRALS,IOSEF1);        !
                                                                                         !
                end do                                                                   !
                                                                                         !
                write(icanal,'(a70)') '| Dihedral Coeffs were read'// &                  !
                                      REPEAT(' ',42)//'|';                               !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            else if ( INDEX(CHAIN_LENGTH,'MiddleBondTorsion Coeffs') > 0 ) then;

                read(2,*);

                do i = 1, NTYPE_DIHEDRAL;

                    read(2,*) IOSEF1, MIDDLEBONDTORSION_COEFFS(1:4,IOSEF1);

                end do

                write(icanal,'(a70)') '| MiddleBondTorsion Coeffs were read'// &
                                      REPEAT(' ',33)//'|';

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|'; 

            else if ( INDEX(CHAIN_LENGTH,'EndBondTorsion Coeffs') > 0 ) then;

                read(2,*);

                do i = 1, NTYPE_DIHEDRAL;

                    read(2,*) IOSEF1, ENDBONDTORSION_COEFFS(1:8,IOSEF1);

                end do

                write(icanal,'(a70)') '| EndBondTorsion Coeffs were read'// &
                                      REPEAT(' ',36)//'|';

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            else if ( INDEX(CHAIN_LENGTH,'AngleTorsion Coeffs # class2') > 0 ) then;

                read(2,*);

                do i = 1, NTYPE_DIHEDRAL;

                    read(2,*) IOSEF1, ANGLETORSION_COEFFS(1:8,IOSEF1);

                end do

                write(icanal,'(a70)') '| AngleTorsion Coeffs were read'// &
                                      REPEAT(' ',38)//'|';

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
 
            else if ( INDEX(CHAIN_LENGTH,'AngleAngleTorsion Coeffs') > 0 ) then;

                read(2,*);

                do i = 1, NTYPE_DIHEDRAL;

                    read(2,*) IOSEF1, ANGLEANGLETORSION_COEFFS(1:3,IOSEF1);

                end do

                write(icanal,'(a70)') '| AngleAngleTorsion Coeffs were read'// &
                                      REPEAT(' ',33)//'|';

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            else if ( INDEX(CHAIN_LENGTH,'BondBond13 Coeffs') > 0 ) then;

                read(2,*);

                do i = 1, NTYPE_DIHEDRAL;

                    read(2,*) IOSEF1, BONDBOND13_COEFFS(1:3,IOSEF1);

                end do

                write(icanal,'(a70)') '| AngleAngleTorsion Coeffs were read'// &
                                      REPEAT(' ',33)//'|';

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            else if ( INDEX(CHAIN_LENGTH,'Improper Coeffs # class2') > 0 ) then;

                read(2,*);

                do i = 1, NTYPE_IMPROPER;

                    read(2,*) IOSEF1, IMPROPER_COEFFS(1:2,IOSEF1);

                end do

                write(icanal,'(a70)') '| Improper Coeffs were read'// &
                                      REPEAT(' ',42)//'|';

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            else if ( INDEX(CHAIN_LENGTH,'AngleAngle Coeffs') > 0 ) then;

                read(2,*);

                do i = 1, NTYPE_IMPROPER;

                    read(2,*) IOSEF1, ANGLEANGLE_COEFFS(1:6,IOSEF1);

                end do

                write(icanal,'(a70)') '| AngleAngle Coeffs were read'// &
                                      REPEAT(' ',40)//'|';

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            else if ( INDEX(CHAIN_LENGTH,'Atoms') > 0 ) then;
  
                ATOMS_FLAG = 0;
  
                if ( INDEX(CHAIN_LENGTH,'charge') > 0 ) then;

                    ATOMS_FLAG = 5;

                    write(icanal,'(a70)') '| The charge option was used to write coordinates'// &
                    REPEAT(' ',20)//'|';

                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

                else if ( INDEX(CHAIN_LENGTH,'full')   > 0 ) then;
 
                    ATOMS_FLAG = 13;

                    write(icanal,'(a70)') '| The full option was used to write coordinates'// &
                                          REPEAT(' ',22)//'|';

                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

                else

                    write(icanal,'(a70)') '| '//REPEAT('/!\',66)//' |';
                    write(icanal,*) '| /!\ WARNING NO ATOM FLAG WAS SET TO '// &
                                    'READ THE CONFIGURATION'//                 &
                                    REPEAT(' ',5)//'/!\ !';
                    write(icanal,'(a70)') '| '//REPEAT('/!\',66)//' |';
                    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
                    write(icanal,'(a70)') '| PLEASE, ADD << # full >> OR << # charge >> '// &
                                          'NEXT TO Atoms IN THE '// &
                                          TRIM(CHEXT)// ' FILE |';
                   write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';
                   write(icanal,*);
                   write(icanal,*) 'END OF PROGRAM';

                   close(icanal);

                   stop;

                end if

                write(icanal,'(a15,i8,a47)') '| ATOMS_FLAG : ', ATOMS_FLAG, REPEAT(' ',46)//'|';

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

                read(2,*);

                do i = 1, NATOM;                                                         !
                                                                                         !
                    read(2,'(a)') CHAIN_LENGTH;                                          !

                    EOF2 = 0;

                    if ( ATOMS_FLAG == 5 ) then;

                        read(CHAIN_LENGTH,*,iostat=EOF2) CONFIG_ATOMID(i),    &
                                                         CONFIG_ATOM_TYPE(i), &
                                                         CONFIG_QI(i),        &
                                                         CONFIG_RI(1:3,i);
                    else if ( ATOMS_FLAG == 13 ) then;

                        read(CHAIN_LENGTH,*,iostat=EOF2) CONFIG_ATOMID(i),     &
                                                         CONFIG_MOLECULEID(i), &
                                                         CONFIG_ATOM_TYPE(i),  &
                                                         CONFIG_QI(i),         &
                                                         CONFIG_RI(1:3,i);

                        if ( CONFIG_ATOMID(i) > FLAG_SORT_MOLEC ) FLAG_SORT_MOLEC = CONFIG_ATOMID(i);

                    else

                        read(CHAIN_LENGTH,*,iostat=EOF2) IOSEF1, IOSEF2, IATOM_TYPE, &
                                                         CONFIG_QI(i), &
                                                         CONFIG_RI(1:3,i);
 
                        if ( EOF2 /= 0 ) then;

                            read(CHAIN_LENGTH,*) IOSEF1, IATOM_TYPE, &
                                                 CONFIG_QI(i), &
                                                 CONFIG_RI(1:3,i);

                        end if

                    end if

                    do j = 1, NTYPE_ATOM;

                        if ( CONFIG_ATOM_TYPE(i) == j ) then;

                            CONFIG_NAT(i) = ATOM_LABEL(j);

                            EXIT;

                        end if

                    end do

                end do

                write(icanal,'(a70)') '| Atoms were read'//REPEAT(' ',52)//'|';

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

!               stop; !//////////////////////////////////////////////////////////////////! 
                                                                                         !
            else if ( INDEX(CHAIN_LENGTH,'Velocities') > 0 ) then;                       !
                                                                                         !
                read(2,*);                                                               !
                                                                                         !
                do i = 1, NATOM;                                                         !
                                                                                         !
                    read(2,*) IOSEF1, CONFIG_VI(1:3,i);                                  !
                                                                                         !
                end do                                                                   !
                                                                                         !
                write(icanal,'(a70)') '| Velocities were read'//REPEAT(' ',47)//'|';     !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            else if ( INDEX(CHAIN_LENGTH,'Bonds') > 0 ) then;                            !

                read(2,*);

                do i = 1, NBOND;

                    read(2,*) IOSEF1, BOND_TYPE(i), BOND_ATOMID(1:2,i);

                end do

                write(icanal,'(a70)') '| Bonds were read'//REPEAT(' ',52)//'|';

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            else if ( INDEX(CHAIN_LENGTH,'Angles') > 0 ) then;                           !
                                                                                         !
                read(2,*);                                                               !
                                                                                         !
                do i = 1, NANGLE;                                                        !
                                                                                         ! 
                    read(2,*) IOSEF1, ANGLE_TYPE(i), ANGLE_ATOMID(1:3,i);                !
                                                                                         !
                end do                                                                   !
                                                                                         !
                write(icanal,'(a70)') '| Angles were read'//REPEAT(' ',51)//'|';         !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            else if ( INDEX(CHAIN_LENGTH,'Dihedrals') > 0 ) then;                        !
                                                                                         !
                read(2,*);

                do i = 1, NDIHEDRAL;

                    read(2,*) IOSEF1, DIHEDRAL_TYPE(i), DIHEDRAL_ATOMID(1:4,i);

                end do

                write(icanal,'(a70)') '| Dihedrals were read'//REPEAT(' ',48)//'|';

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|'; 

!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            else if ( INDEX(CHAIN_LENGTH,'Impropers') > 0 ) then;                        !
                                                                                         !
                read(2,*);

                do i = 1, NIMPROPER;

                    read(2,*) IOSEF1, IMPROPER_TYPE(i), IMPROPER_ATOMID(1:4,i);

                end do

                write(icanal,'(a70)') '| Impropers were read'//REPEAT(' ',48)//'|';

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

!               stop; !//////////////////////////////////////////////////////////////////! 
                                                                                         !
            end if

        end do

        close(2);

!   end if

!   ### SHARE DATA WITH THE SLAVE PROCESSORS #######################################################
!   call MPI_BCAST(NATOM,1,MPI_INTEGER,root_process,MPI_COMM_WORLD,ierr);

!   call MPI_BCAST(NTYPE_ATOM,1,MPI_INTEGER,root_process,MPI_COMM_WORLD,ierr);

!   if ( my_id /= root_process ) then
!       allocate(ATOM_MASSE(1:NTYPE_ATOM));

!       allocate(CONFIG_MOLECULEID(1:NATOM));

!       allocate(CONFIG_ATOM_TYPE(1:NATOM));

!       allocate(CONFIG_RI(1:3,1:NATOM));
!   end if

!   call MPI_BCAST(ATOM_MASSE,     &
!                  NTYPE_ATOM,     &
!                  MPI_INTEGER,    &                                                     !
!                  root_process,   &                                                     !
!                  MPI_COMM_WORLD, &                                                     !
!                  ierr);                                                                !
                                                                                         !
!   call MPI_BCAST(CONFIG_MOLECULEID, &                                                  ! SHARE IDs OF MOLECULES WITH OTHER PROCESSORS
!                  NATOM,             &                                                  !
!                  MPI_INTEGER,       &                                                  !
!                  root_process,      &                                                  !
!                  MPI_COMM_WORLD,    &                                                  !
!                  ierr);                                                                !
                                                                                         !
!   call MPI_BCAST(CONFIG_ATOM_TYPE,NATOM,MPI_INTEGER,root_process,MPI_COMM_WORLD,ierr);

!   do i = 1, 3
!       call MPI_BCAST(CONFIG_RI(i,:),NATOM,MPI_DOUBLE_PRECISION,root_process,MPI_COMM_WORLD,ierr);
!   end do

!   ### WRAP UP THE ROUTINE PROCESS #####################################################!
                                                                                         !
!   if ( my_id == root_process ) then;                                                   !
                                                                                         !
        call CPU_TIME(t_stop);                                                           !
                                                                                         !
        tps1 = ( t_stop - t_start ) / 3600.0d0;                                          !
                                                                                         !
        tps2 = ( tps1 - INT( tps1 ) ) * 60;                                              !
                                                                                         !
        tps3 = ( tps2 - INT( tps2 ) ) * 60;                                              !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        write(icanal,'(a19,i6,a3,i6,a5,i6,a25)') '| COMPUTING TIME : ',    &             !
                                                 INT( tps1 ),'[h]',        &             !
                                                 INT( tps2 ),'[min]',      &             !
                                                 INT( tps3 ),'[s]'//REPEAT(' ',21)//'|'; !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                  !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine READ_LAMMPS_CONFIG
