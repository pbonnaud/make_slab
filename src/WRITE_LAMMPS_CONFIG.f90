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



subroutine WRITE_LAMMPS_CONFIG(icanal,CHEXT,CELL_AXIS,MATA) 

!   ************************************************************************************************
!   **                                 READ LAMMPS FILES                                          **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal     : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                                      **
!   **                                                                                            **
!   ** CHEXT      : NAME OF THE LAMMPS CONFIGURATION                                              **
!   **                                                                                            **
!   ** CELL_AXIS  : CELL DIMENSIONS                                                               **
!   **                                                                                            **
!   ** MATA       : PASSAGE MATRIX                                                                **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

    use module_size_arrays;

    use module_config;

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ***********************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    character (len=150), intent(in) :: CHEXT;

    real (kind=8), dimension(1:3), intent(in) :: CELL_AXIS;

    real (kind=8), dimension(1:6), intent(in) :: MATA;

!   ************************************************************************************************

    integer (kind=4) :: i, j;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the current routine #####################################################
                                                                                         !
    CHTITLE = 'Write LAMMPS configuration';                                              !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the name of the lammps configuration to write ########################################
                                                                                         !
    IOSEF1 = 70 - 2 - 1 - LEN_TRIM(CHEXT);                                               !
                                                                                         !
    if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                        !
                                                                                         !
    write(icanal,'(a70)') '| '//TRIM(CHEXT)//REPEAT(' ',IOSEF1)//'|';                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the LAMMPS configuration #############################################################
                                                                                         !
    open(103,file=TRIM(CHEXT));                                                          !
                                                                                         !
    write(103,'(a52)') '# System description ###############################';           !
                                                                                         !
    write(103,*);                                                                        !
                                                                                         !
    write(103,'(a6,3f15.6)') '#     ', CELL_AXIS(1:3);                                   !
                                                                                         !
    write(103,*);                                                                        !
                                                                                         !
!   ### Write properties of the molecular configuration ############################################
                                                                                         !
    write(103,'(i8,1x,a5)') NATOM, 'atoms';                                              !
                                                                                         !
    if ( NBOND     > 0 ) write(103,'(i8,1x,a5)') NBOND,     'bonds';                     !
                                                                                         !
    if ( NANGLE    > 0 ) write(103,'(i8,1x,a6)') NANGLE,    'angles';                    !
                                                                                         !
    if ( NDIHEDRAL > 0 ) write(103,'(i8,1x,a9)') NDIHEDRAL, 'dihedrals';                 !
                                                                                         !
    if ( NIMPROPER > 0 ) write(103,'(i8,1x,a9)') NIMPROPER, 'impropers';                 !
                                                                                         !
    write(103,*);                                                                        !
                                                                                         !
    if ( NTYPE_ATOM     > 0 ) write(103,'(i8,1x,a10)') NTYPE_ATOM,     'atom types';     !
                                                                                         !
    if ( NTYPE_BOND     > 0 ) write(103,'(i8,1x,a10)') NTYPE_BOND,     'bond types';     !
                                                                                         !
    if ( NTYPE_ANGLE    > 0 ) write(103,'(i8,1x,a11)') NTYPE_ANGLE,    'angle types';    !
                                                                                         !
    if ( NTYPE_DIHEDRAL > 0 ) write(103,'(i8,1x,a14)') NTYPE_DIHEDRAL, 'dihedral types'; !
                                                                                         !
    if ( NTYPE_IMPROPER > 0 ) write(103,'(i8,1x,a14)') NTYPE_IMPROPER, 'improper types'; !
                                                                                         !
    write(103,*);                                                                        !
                                                                                         !
    write(icanal,'(a70)') '| Properties of the molecular configuration were written'// & !
                          REPEAT(' ',13)//'|';                                           !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write properties of the simulation box #####################################################
                                                                                         !
    write(103,'(2f15.6,a11)') -0.5d0 * MATA(1), 0.5d0 * MATA(1), '    xlo xhi';          !
                                                                                         !
    write(103,'(2f15.6,a11)') -0.5d0 * MATA(2), 0.5d0 * MATA(2), '    ylo yhi';          !
                                                                                         !
    write(103,'(2f15.6,a11)') -0.5d0 * MATA(3), 0.5d0 * MATA(3), '    zlo zhi';          ! 
                                                                                         !
    if ( ( ABS(MATA(6)) > 1.0E-4 ) .OR.                                   &              !
         ( ABS(MATA(5)) > 1.0E-4 ) .OR.                                   &              !
         ( ABS(MATA(4)) > 1.0E-4 ) ) write(103,'(3f15.6,a9)') MATA(6),    &              !
                                                              MATA(5),    &              !
                                                              MATA(4),    &              !
                                                              ' xy xz yz';               !
    write(103,*);                                                                        !
                                                                                         !
    write(icanal,'(a70)') '| Properties of the simulation box were written'// &          !
                          REPEAT(' ',22)//'|';                                           !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the list of masses ###################################################################
                                                                                         !
    if ( NTYPE_ATOM > 0 ) then;                                                          !
                                                                                         !
        write(103,'(a6)') 'Masses';                                                      !
                                                                                         !
        write(103,*);                                                                    !
                                                                                         !
        do i = 1, NTYPE_ATOM;                                                            !
                                                                                         !
            write(103,'(i4,f12.6)') i, ATOM_MASSE(i);                                    !
                                                                                         !
        end do                                                                           !
                                                                                         !
        write(103,*);                                                                    !
                                                                                         !
        write(icanal,'(a70)') '| Masses were written'// &                                !
                              REPEAT(' ',48)//'|';                                       !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
    end if                                                                               !
                                                                                         ! 
!   ### Write pair potentials ######################################################################
                                                                                         !
    if ( IPOTENTIAL_CLASS2 > 0 ) then;                                                   !
                                                                                         !
        CHOSEF1 = 'Pair Coeffs # '//TRIM(CH_PAIR_STYLE);                                 !
                                                                                         !
        IOSEF1 = LEN_TRIM(CHOSEF1);                                                      !
                                                                                         !
        if ( IOSEF1 < 1000 ) write(CHOSEF2,'(i3)') IOSEF1;                               !
                                                                                         !
        if ( IOSEF1 <  100 ) write(CHOSEF2,'(i2)') IOSEF1;                               !
                                                                                         !
        if ( IOSEF1 <   10 ) write(CHOSEF2,'(i1)') IOSEF1;                               !
                                                                                         !
        write(103,*);                                                                    !
                                                                                         !
        write(103,'(a'//TRIM(CHOSEF2)//')') TRIM(CHOSEF1);                               !
                                                                                         !
        write(103,*);                                                                    !
                                                                                         !
        do i = 1, NTYPE_ATOM;                                                            !
                                                                                         !
            write(103,'(i4,2f12.4)') i, POTENTIAL_CLASS2(1:2,i);                         !
                                                                                         !
        end do                                                                           !
                                                                                         !
        write(icanal,'(a70)') '| Pair potential parameters were written'// &             !
                              REPEAT(' ',29)//'|';                                       !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Write bond potential parameters ############################################################
                                                                                         !
    if ( NTYPE_BOND > 0 ) then;                                                          !
                                                                                         !
        CHOSEF1 = 'Bond Coeffs # '//TRIM(CH_BOND_STYLE);                                 !
                                                                                         !
        IOSEF1 = LEN_TRIM(CHOSEF1);                                                      !
                                                                                         !
        if ( IOSEF1 < 1000 ) write(CHOSEF2,'(i3)') IOSEF1;                               !
                                                                                         !
        if ( IOSEF1 <  100 ) write(CHOSEF2,'(i2)') IOSEF1;                               !
                                                                                         !
        if ( IOSEF1 <   10 ) write(CHOSEF2,'(i1)') IOSEF1;                               !
                                                                                         !
        write(103,*);                                                                    !
                                                                                         !
        write(103,'(a'//TRIM(CHOSEF2)//')') TRIM(CHOSEF1);                               !
                                                                                         !
        write(103,*);                                                                    !
                                                                                         !
        do i = 1, NTYPE_BOND;                                                            !
                                                                                         !
            write(103,'(i4,4f12.4)') i, BOND_COEFFS(1:NPARAM_BONDS,i);                   !
                                                                                         !
        end do                                                                           !
                                                                                         !
        write(icanal,'(a70)') '| Bond coeffs were written'//REPEAT(' ',43)//'|';         !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Write angle coefficients ###################################################################
                                                                                         !
    if ( NTYPE_ANGLE > 0 ) then;                                                         !
                                                                                         !
        CHOSEF1 = 'Angle Coeffs # '//TRIM(CH_ANGLE_STYLE);                               !
                                                                                         !
        IOSEF1 = LEN_TRIM(CHOSEF1);                                                      !
                                                                                         !
        if ( IOSEF1 < 1000 ) write(CHOSEF2,'(i3)') IOSEF1;                               !
                                                                                         !
        if ( IOSEF1 <  100 ) write(CHOSEF2,'(i2)') IOSEF1;                               !
                                                                                         !
        if ( IOSEF1 <   10 ) write(CHOSEF2,'(i1)') IOSEF1;                               !
                                                                                         !
        write(103,*);                                                                    !
                                                                                         !
        write(103,'(a'//TRIM(CHOSEF2)//')') TRIM(CHOSEF1);                               !
                                                                                         !
        write(103,*);                                                                    !
                                                                                         !
        do i = 1, NTYPE_ANGLE;                                                           !
                                                                                         !
            write(103,'(i4,4f12.4)') i, ANGLE_COEFFS(1:NPARAM_ANGLES,i);                 !
                                                                                         !
        end do                                                                           !
                                                                                         !
        write(icanal,'(a70)') '| Angle coeffs were written'//REPEAT(' ',42)//'|';        !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
!       ### Write additional angle coefficients for class2 potentials ##############################
                                                                                         !
        if ( TRIM(CH_ANGLE_STYLE) == 'class2' ) then;                                    !
                                                                                         !
            write(103,*);                                                                !
                                                                                         !
            write(103,'(a15)') 'BondBond Coeffs';                                        !
                                                                                         !
            write(103,*);                                                                !
                                                                                         !
            do i = 1, NTYPE_ANGLE;                                                       !
                                                                                         !
                write(103,'(i4,3f12.4)') i, BONDBOND_COEFFS(1:3,i);                      !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| BondBond coeffs were written'//REPEAT(' ',42)//'|'; !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
            write(103,*);                                                                !
                                                                                         !
            write(103,'(a16)') 'BondAngle Coeffs';                                       !
                                                                                         !
            write(103,*);                                                                !
                                                                                         !
            do i = 1, NTYPE_ANGLE;                                                       !
                                                                                         !
                write(103,'(i4,4f12.4)') i, BONDANGLE_COEFFS(1:4,i);                     !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| BondAngle coeffs were written'//REPEAT(' ',42)//'|';!
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        end if                                                                           !
                                                                                         !
!       stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Write dihedral coefficients ################################################################
                                                                                         !
    if ( NTYPE_DIHEDRAL > 0 ) then;                                                      !
                                                                                         !
        CHOSEF1 = 'Dihedral Coeffs # '//TRIM(CH_DIHEDRAL_STYLE);                         !
                                                                                         !
        IOSEF1 = LEN_TRIM(CHOSEF1);                                                      !
                                                                                         !
        if ( IOSEF1 < 1000 ) write(CHOSEF2,'(i3)') IOSEF1;                               !
                                                                                         !
        if ( IOSEF1 <  100 ) write(CHOSEF2,'(i2)') IOSEF1;                               !
                                                                                         !
        if ( IOSEF1 <   10 ) write(CHOSEF2,'(i1)') IOSEF1;                               !
                                                                                         !
        write(103,*);                                                                    !
                                                                                         !
        write(103,'(a'//TRIM(CHOSEF2)//')') TRIM(CHOSEF1);                               !
                                                                                         ! 
        write(103,*);                                                                    !
                                                                                         !
        IOSEF2 = 0;                                                                      !
                                                                                         !
        if ( TRIM(CH_DIHEDRAL_STYLE) == 'charmmfsw' ) then;                              !
                                                                                         !
            CHOSEF3 = 'i8,f12.6,2i8,f12.6';                                              !
                                                                                         !
            IOSEF2 = 1;                                                                  ! 
                                                                                         !
        else                                                                             !
                                                                                         !
            write(CHOSEF3,'(i1)') NPARAM_DIHEDRALS;                                      !
                                                                                         !
            CHOSEF3 = 'i8,'//TRIM(CHOSEF3)//'f12.6';                                     !
                                                                                         !
        end if                                                                           !
                                                                                         !
        do i = 1, NTYPE_DIHEDRAL;                                                        !
                                                                                         !
            if ( IOSEF2 == 1 ) then;                                                     !
                                                                                         !
                write(103,'('//TRIM(CHOSEF3)//')') i,                         &          !
                                                   DIHEDRAL_COEFFS(1,i),      &          !
                                                   INT(DIHEDRAL_COEFFS(2,i)), &          !
                                                   INT(DIHEDRAL_COEFFS(3,i)), &          !
                                                   DIHEDRAL_COEFFS(4,i);                 !
                                                                                         !
            else                                                                         !
                                                                                         !
                write(103,'('//TRIM(CHOSEF3)//')') i,  &                                 !
                                                   DIHEDRAL_COEFFS(1:NPARAM_DIHEDRALS,i);!
                                                                                         !
            end if                                                                       !
                                                                                         !
        end do                                                                           !
                                                                                         !
!       ### Write additional parameters related to class2 potential ################################
                                                                                         !
        if ( TRIM(CH_DIHEDRAL_STYLE) == 'class2' ) then;                                 !
                                                                                         !
            write(103,*);                                                                !
                                                                                         !
            write(103,'(a24)') 'MiddleBondTorsion Coeffs';                               !
                                                                                         !
            write(103,*);                                                                !
                                                                                         ! 
            do i = 1, NTYPE_DIHEDRAL;                                                    !
                                                                                         !
                write(103,'(i4,4f12.4)') i, MIDDLEBONDTORSION_COEFFS(1:4,i);             !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(103,*);

            write(103,'(a21)') 'EndBondTorsion Coeffs';

            write(103,*);

            do i = 1, NTYPE_DIHEDRAL;

                write(103,'(i4,8f12.4)') i, ENDBONDTORSION_COEFFS(1:8,i);

            end do

            write(103,*);

            write(103,'(a28)') 'AngleTorsion Coeffs # class2';

            write(103,*);

            do i = 1, NTYPE_DIHEDRAL;

                write(103,'(i4,8f12.4)') i, ANGLETORSION_COEFFS(1:8,i);

            end do

            write(103,*);

            write(103,'(a24)') 'AngleAngleTorsion Coeffs';

            write(103,*);

            do i = 1, NTYPE_DIHEDRAL;

                write(103,'(i4,3f12.4)') i, ANGLEANGLETORSION_COEFFS(1:3,i);

            end do

            write(103,*);

            write(103,'(a17)') 'BondBond13 Coeffs';

            write(103,*);

            do i = 1, NTYPE_DIHEDRAL;

                write(103,'(i4,3f12.4)') i, BONDBOND13_COEFFS(1:3,i);

            end do

        end if

    end if

!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    if ( NTYPE_IMPROPER > 0 ) then;                                                      !
                                                                                         ! 
        write(103,*);                                                                    !
                                                                                         !
        write(103,'(a24)') 'Improper Coeffs # class2';                                   ! 
                                                                                         !
        write(103,*);

        do i = 1, NTYPE_IMPROPER;

            write(103,'(i4,2f12.4)') i, IMPROPER_COEFFS(1:2,i);

        end do

        write(103,*);

        write(103,'(a17)') 'AngleAngle Coeffs';

        write(103,*);

        do i = 1, NTYPE_IMPROPER;

            write(103,'(i4,6f12.4)') i, ANGLEANGLE_COEFFS(1:6,i);

        end do

    end if
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set the format of coordinates ##############################################################
                                                                                         !
    CHOSEF1 = 'Atoms #'//TRIM(CH_ATOM_STYLE);                                            !
                                                                                         !
    if ( INDEX(TRIM(CH_ATOM_STYLE),'charge') > 0 ) ATOMS_FLAG = 5;                       !
                                                                                         !
    if ( INDEX(TRIM(CH_ATOM_STYLE),'full') > 0 ) ATOMS_FLAG = 13;                        !
                                                                                         !
!   if ( ATOMS_FLAG == 5 ) CHOSEF1 = 'Atoms # charge';                                   !
                                                                                         !
!   if ( ATOMS_FLAG == 13 ) CHOSEF1 = 'Atoms # full';                                    !
                                                                                         !
    if ( iuse_moleculeid == 1 ) then;                                                    !
                                                                                         !
        ATOMS_FLAG = 13;                                                                 !
                                                                                         !
        CHOSEF1 = 'Atoms # full';                                                        !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   write(icanal,*) TRIM(CHOSEF1);

!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write properties of coordinates ############################################################
                                                                                         !
    IOSEF1 = LEN_TRIM(CHOSEF1);                                                          !
                                                                                         !
    if ( IOSEF1 < 1000 ) write(CHOSEF2,'(i3)') IOSEF1;                                   !
                                                                                         !
    if ( IOSEF1 <  100 ) write(CHOSEF2,'(i2)') IOSEF1;                                   !
                                                                                         !
    if ( IOSEF1 <   10 ) write(CHOSEF2,'(i1)') IOSEF1;                                   !
                                                                                         !
    write(103,*);                                                                        !
                                                                                         !
    write(103,'(a'//TRIM(CHOSEF2)//')') TRIM(CHOSEF1);                                   !
                                                                                         !
    write(103,*);                                                                        !
                                                                                         !
!   ### Write atomic coordinates of the molecular configuration ####################################
                                                                                         !
    do i = 1, NATOM;                                                                     !
                                                                                         !
        if ( ATOMS_FLAG == 5 ) then;                                                     !
                                                                                         !
            write(103,'(2i8,4f21.15)') CONFIG_ATOMID(i),     &                           !
                                       CONFIG_ATOM_TYPE(i),  &                           !
                                       CONFIG_QI(i),         &                           !
                                       CONFIG_RI(1:3,i);                                 !
                                                                                         !
        else if ( ATOMS_FLAG == 13 ) then;                                               !
                                                                                         ! 
            write(103,'(3i8,4f21.15)') CONFIG_ATOMID(i),     &                           !
                                       CONFIG_MOLECULEID(i), &                           !
                                       CONFIG_ATOM_TYPE(i),  &                           !
                                       CONFIG_QI(i),         &                           !
                                       CONFIG_RI(1:3,i);                                 !
                                                                                         !
        else                                                                             !
                                                                                         !
            stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end if                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   ### Write atomic velocities of the molecular configuration #####################################
                                                                                         !
    write(103,*);                                                                        !
                                                                                         !
    write(103,'(a10)') 'Velocities';                                                     !
                                                                                         !
    write(103,*);                                                                        !
                                                                                         !
    do i = 1, NATOM;                                                                     !
                                                                                         !
        write(103,'(i8,3f21.15)') CONFIG_ATOMID(i), CONFIG_VI(1:3,i);                    !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   ### Write bonds in the molecular configuration #################################################
                                                                                         !
    if ( NBOND > 0 ) then;                                                               !
                                                                                         ! 
        write(103,*);                                                                    !
                                                                                         !
        write(103,'(a5)') 'Bonds';                                                       !
                                                                                         !
        write(103,*);                                                                    !
                                                                                         !
        do i = 1, NBOND;                                                                 !
                                                                                         !
            write(103,'(4i8)') i, BOND_TYPE(i), BOND_ATOMID(1:2,i);                      !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write angles in the molecular configuration ################################################
                                                                                         !
    if ( NANGLE > 0 ) then;                                                              !
                                                                                         !
        write(103,*);                                                                    !
                                                                                         !
        write(103,'(a6)') 'Angles';                                                      !
                                                                                         !
        write(103,*);                                                                    !
                                                                                         !
        do i = 1, NANGLE;                                                                !
                                                                                         !
            write(103,'(5i8)') i, ANGLE_TYPE(i), ANGLE_ATOMID(1:3,i);                    !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write dihedrals in the molecular configuration #############################################
                                                                                         !
    if ( NDIHEDRAL > 0 ) then;                                                           !
                                                                                         !
        write(103,*);                                                                    !
                                                                                         !
        write(103,'(a9)') 'Dihedrals';                                                   !
                                                                                         !
        write(103,*);                                                                    !
                                                                                         !
        do i = 1, NDIHEDRAL;                                                             !
                                                                                         !
            write(103,'(6i8)') i, DIHEDRAL_TYPE(i), DIHEDRAL_ATOMID(1:4,i);              !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Write impropers in the molecular configuration #############################################
                                                                                         !
    if ( NIMPROPER > 0 ) then;                                                           !
                                                                                         !
        write(103,*);

        write(103,'(a9)') 'Impropers';

        write(103,*);

        do i = 1, NIMPROPER;

            write(103,'(6i8)') i, IMPROPER_TYPE(i), IMPROPER_ATOMID(1:4,i);

        end do

    end if

    close(103);
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine WRITE_LAMMPS_CONFIG

