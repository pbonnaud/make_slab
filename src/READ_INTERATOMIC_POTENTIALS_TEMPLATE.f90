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


subroutine READ_INTERATOMIC_POTENTIALS_TEMPLATE(icanal,IFILE,IINSREM)

!   ************************************************************************************************
!   **                                READ LAMMPS TEMPLATE FILES                                  **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                         : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                  **
!   **                                                                                            **
!   ** CHEXT                          : NAME OF THE LAMMPS CONFIGURATION                          **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

    use module_size_arrays;

    use module_library;

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ***********************************************************************************************

    integer (kind=4), intent(in) :: icanal, IFILE, IINSREM;

!   ************************************************************************************************

    character (len=250) :: CHEXT;

!   ************************************************************************************************

    integer (kind=4) :: i, j;

    integer (kind=4) :: EOF, EOF2;

    character (len=250) :: CHARLINE, CHARLINE2;

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the routine #############################################################
                                                                                         !
    CHTITLE = 'Read Interatomic Potentials Template';                                    !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the file number in the library list ##################################################
                                                                                         !
    write(icanal,'(a10,i8,a52)') '| IFILE : ', IFILE, REPEAT(' ',51)//'|'                !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Build and write the name of the file containing parameters in the library ##################
                                                                                         !
    if ( IINSREM == 1 ) then;                                                            !
                                                                                         !
        CHEXT = TRIM(CHNAME_FILE_LIBRARY(2,IFILE));                                      !
                                                                                         !
    else                                                                                 !
                                                                                         ! 
        CHEXT = TRIM(CHNAME_FILE_LIBRARY_REMOVE(IFILE))//'.parameter';                   !
                                                                                         !
    end if                                                                               !
                                                                                         !
    IOSEF1 = 67 - LEN_TRIM(CHEXT);                                                       !
                                                                                         !
    if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                        !
                                                                                         !
    write(icanal,'(a70)') '| '//TRIM(CHEXT)//REPEAT(' ',IOSEF1)//'|';                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialization of arrays ###################################################################
                                                                                         !
    IPOTENTIAL_CLASS2_LIBRARY(IFILE) = 0;                                                !
                                                                                         !
    if ( IMAX_NTYPE_NATOM > 0 ) then;                                                    !
                                                                                         !
        ATOM_MASSES_LIBRARY(1:IMAX_NTYPE_NATOM,IFILE) = 0.0d0                            !
                                                                                         !
        ATOM_LABEL_LIBRARY(1:IMAX_NTYPE_NATOM,IFILE) = 'XXX';                            !
                                                                                         !
        POTENTIAL_CLASS2_LIBRARY(1:2,1:IMAX_NTYPE_NATOM,IFILE) = 0.0d0;                  !
                                                                                         !
    end if                                                                               !   
                                                                                         !
    if ( IMAX_NTYPE_BOND > 0 ) then;                                                     !
                                                                                         !
        BOND_COEFFS_LIBRARY(1:4,1:IMAX_NTYPE_BOND,IFILE) = 0.0d0;                        !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( IMAX_NTYPE_ANGLE > 0 ) then;                                                    !
                                                                                         !
        ANGLE_COEFFS_LIBRARY(1:4,1:IMAX_NTYPE_ANGLE,IFILE) = 0.0d0;                      !
                                                                                         !
        BONDBOND_COEFFS_LIBRARY(1:3,1:IMAX_NTYPE_ANGLE,IFILE) = 0.0d0;                   !
                                                                                         !
        BONDANGLE_COEFFS_LIBRARY(1:4,1:IMAX_NTYPE_ANGLE,IFILE) = 0.0d0;                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( IMAX_NTYPE_DIHEDRAL > 0 ) then;                                                 !
                                                                                         !
        DIHEDRAL_COEFFS_LIBRARY(1:6,1:IMAX_NTYPE_DIHEDRAL,IFILE)          = 0.0d0;       !
                                                                                         !
        MIDDLEBONDTORSION_COEFFS_LIBRARY(1:4,1:IMAX_NTYPE_DIHEDRAL,IFILE) = 0.0d0;       !
                                                                                         ! 
        ENDBONDTORSION_COEFFS_LIBRARY(1:8,1:IMAX_NTYPE_DIHEDRAL,IFILE)    = 0.0d0;       !
                                                                                         !
        ANGLETORSION_COEFFS_LIBRARY(1:8,1:IMAX_NTYPE_DIHEDRAL,IFILE)      = 0.0d0;       !
                                                                                         !
        ANGLEANGLETORSION_COEFFS_LIBRARY(1:3,1:IMAX_NTYPE_DIHEDRAL,IFILE) = 0.0d0;       !
                                                                                         !
        BONDBOND13_COEFFS_LIBRARY(1:3,1:IMAX_NTYPE_DIHEDRAL,IFILE)        = 0.0d0;       !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    if ( IMAX_NTYPE_IMPROPER > 0 ) then;                                                 !
                                                                                         !
        IMPROPER_COEFFS_LIBRARY(1:2,1:IMAX_NTYPE_IMPROPER,IFILE)          = 0.0d0;       !
                                                                                         !
        ANGLEANGLE_COEFFS_LIBRARY(1:6,1:IMAX_NTYPE_IMPROPER,IFILE)        = 0.0d0;       !
                                                                                         !
    end if                                                                               !
                                                                                         !

!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write message confirming the initialization of parameter arrays ############################
                                                                                         !
    write(icanal,'(a70)') '| Parameter arrays were initialized'//REPEAT(' ',34)//'|';    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Read parameters of the library molecule ####################################################
                                                                                         !
    open(5,file=TRIM(CHEXT));                                                            !
                                                                                         ! 
    EOF = 0;                                                                             !
                                                                                         !
    do                                                                                   !
                                                                                         !
        read(5,'(a)',iostat=EOF) CHARLINE;                                               !
                                                                                         !
        if ( EOF /= 0 ) EXIT;                                                            !
                                                                                         !
        if ( INDEX(CHARLINE,'Masses') > 0 ) then;                                        !
                                                                                         !
            read(5,*);                                                                   !
                                                                                         !
            if ( NTYPE_ATOM_LIBRARY(IFILE) == 0 ) CYCLE;                                 !
                                                                                         !
            do i = 1, NTYPE_ATOM_LIBRARY(IFILE);                                         !
                                                                                         !
                read(5,*) IOSEF1, ATOM_MASSES_LIBRARY(IOSEF1,IFILE);                     !
                                                                                         !
            end do                                                                       !
                                                                                         !
            do i = 1, NTYPE_ATOM_LIBRARY(IFILE);                                         !
                                                                                         !
                do j = 1, 200;                                                           !
                                                                                         !
                    ROSEF1 = ABS( ATOM_MASSES_LIBRARY(i,IFILE) - &                       !
                                  MOLAR_MASS_ELEMENT(j) );                               !
                                                                                         !
                    if ( ROSEF1 >= 0.1d0 ) CYCLE;                                        !
                                                                                         !
!                   if ( ROSEF1 < 0.1d0 ) then;                                          !
                                                                                         !
                    ATOM_LABEL_LIBRARY(i,IFILE) = ELEMENT_SYMBOL_NAME(j);                !
                                                                                         !
                    EXIT;                                                                !
                                                                                         !
!                   end if                                                               !
                                                                                         !
                end do                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            write(icanal,'(a70)') '| Template masses were read'// &                      !
                                      REPEAT(' ',42)//'|';                               !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           end if                                                                       !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////! 
                                                                                         !
        else if ( INDEX(CHARLINE,'Pair Coeffs') > 0 ) then;                              !
                                                                                         !
            if ( NTYPE_ATOM_LIBRARY(IFILE) > 0 ) then;                                   !
                                                                                         !
                IPOTENTIAL_CLASS2_LIBRARY(IFILE) = 1;                                    !
                                                                                         !
                if ( INDEX(CHARLINE,'# lj/class2/coul/long') > 0 ) then;                 !
                                                                                         !
                    POTENTIAL_CLASS2_CHTYPE_LIBRARY(IFILE) = 'lj/class2/coul/long';      !
                                                                                         !
                else if ( INDEX(CHARLINE,'# lj/cut/coul/long') > 0 ) then;               !
                                                                                         !
                    POTENTIAL_CLASS2_CHTYPE_LIBRARY(IFILE) = 'lj/cut/coul/long';         !
                                                                                         !
                else if ( INDEX(CHARLINE,'# lj/charmmfsw/coul/long') > 0 ) then;         !
                                                                                         !
                    POTENTIAL_CLASS2_CHTYPE_LIBRARY(IFILE) = 'lj/charmmfsw/coul/long';   !
                                                                                         !
                end if                                                                   !
                                                                                         !
                read(5,*);                                                               !
                                                                                         !
                do i = 1, NTYPE_ATOM_LIBRARY(IFILE);                                     !
                                                                                         !
                    read(5,*) IOSEF1, POTENTIAL_CLASS2_LIBRARY(1:2,IOSEF1,IFILE);        !
                                                                                         !
                end do                                                                   !
                                                                                         !
                write(icanal,'(a70)') '| Pair coefficients were read'// &                !
                                      REPEAT(' ',40) //'|';                              !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        else if ( INDEX(CHARLINE,' Bond Coeffs ') > 0 ) then;                            !
                                                                                         !
            EOF2 = 0;                                                                    !
                                                                                         !
            if ( NTYPE_BOND_LIBRARY(IFILE) > 0 ) then;                                   !
                                                                                         !
                if ( INDEX(CHARLINE,'# class2') > 0 ) then;                              !
                                                                                         !
                    CH_BOND_STYLE_LIBRARY(IFILE) = 'class2';                             !
                                                                                         !
                    NPARAM_BONDS_LIBRARY(IFILE) = 4;                                     ! 
                                                                                         !
                else if ( INDEX(CHARLINE,'# harmonic') > 0 ) then;                       !
                                                                                         !
                    CH_BOND_STYLE_LIBRARY(IFILE) = 'harmonic';                           !
                                                                                         !
                    NPARAM_BONDS_LIBRARY(IFILE) = 2;                                     !
                                                                                         !
                end if                                                                   !
                                                                                         !
                read(5,*);                                                               !
                                                                                         !
                do i = 1, NTYPE_BOND_LIBRARY(IFILE);                                     ! LOOP OVER THE NUMBER OF BOND TYPES IN THE MOLECULAR CONFIGURATION
                                                                                         !
                    IOSEF2 = NPARAM_BONDS_LIBRARY(IFILE);                                !
                                                                                         !
                    read(5,*,iostat=EOF2) IOSEF1,                                     &  !
                                          BOND_COEFFS_LIBRARY(1:IOSEF2,IOSEF1,IFILE);    !
                                                                                         !
                    if ( EOF2 /= 0 ) EXIT;                                               !
                                                                                         ! 
                end do                                                                   !
                                                                                         !
                write(icanal,'(a70)') '| Bond coefficients were read'// &                !
                                      REPEAT(' ',40)//'|';                               !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////! 
                                                                                         !
        else if ( INDEX(CHARLINE,' Angle Coeffs ') > 0 ) then;                           !
                                                                                         !
            EOF2 = 0;                                                                    !
                                                                                         !
            if ( NTYPE_ANGLE_LIBRARY(IFILE) > 0 ) then;                                  !
                                                                                         !
                if ( INDEX(CHARLINE,'# class2') > 0 ) then;                              !
                                                                                         !
                    CH_ANGLE_STYLE_LIBRARY(IFILE) = 'class2';                            !
                                                                                         !
                    NPARAM_ANGLES_LIBRARY(IFILE) = 4;                                    ! 
                                                                                         !
                else if ( INDEX(CHARLINE,'# harmonic') > 0 ) then;                       !
                                                                                         !
                    CH_ANGLE_STYLE_LIBRARY(IFILE) = 'harmonic';                          !
                                                                                         !
                    NPARAM_ANGLES_LIBRARY(IFILE) = 2;                                    !
                                                                                         !
                end if                                                                   !
                                                                                         !
                read(5,*);                                                               !
                                                                                         !
                do i = 1, NTYPE_ANGLE_LIBRARY(IFILE);                                    ! LOOP OVER THE NUMBER OF ANGLE TYPES IN THE MOLECULAR CONFIGURATION
                                                                                         !
                    IOSEF2 = NPARAM_ANGLES_LIBRARY(IFILE);                               !
                                                                                         !
                    read(5,*,iostat=EOF2) IOSEF1,                                 &      !
                                          ANGLE_COEFFS_LIBRARY(1:IOSEF2,IOSEF1,IFILE);   !
                                                                                         !
                    if ( EOF2 /= 0 ) EXIT;                                               !
                                                                                         !
                end do                                                                   !
                                                                                         !     
                write(icanal,'(a70)') '| Angle coefficients were read'// &               !
                                      REPEAT(' ',39)//'|';                               !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        else if ( INDEX(CHARLINE,' BondBond Coeffs') > 0 ) then;                         !
                                                                                         !
            if ( NTYPE_ANGLE_LIBRARY(IFILE) > 0 ) then;                                  !
                                                                                         !
                read(5,*);                                                               !
                                                                                         !
                do i = 1, NTYPE_ANGLE_LIBRARY(IFILE);                                    !
                                                                                         !
                    read(5,*) IOSEF1, BONDBOND_COEFFS_LIBRARY(1:3,IOSEF1,IFILE);         !
                                                                                         !
                end do

                write(icanal,'(a70)') '| BONDBOND COEFFICIENTS WERE READ'// &
                                      REPEAT(' ',36)//'|';

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            end if

        else if ( INDEX(CHARLINE,' BondAngle Coeffs') > 0 ) then;                        !
                                                                                         !
            if ( NTYPE_ANGLE_LIBRARY(IFILE) > 0 ) then;                                  !

                read(5,*);

                do i = 1, NTYPE_ANGLE_LIBRARY(IFILE);

                    read(5,*) IOSEF1, BONDANGLE_COEFFS_LIBRARY(1:4,IOSEF1,IFILE);

                end do

                write(icanal,'(a70)') '| BONDANGLE COEFFICIENTS WERE READ'// &
                                      REPEAT(' ',35)//'|';

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            end if
 
        else if ( INDEX(CHARLINE,' Dihedral Coeffs') > 0 ) then;                         !
                                                                                         !
            if ( NTYPE_DIHEDRAL_LIBRARY(IFILE) > 0 ) then;                               !
                                                                                         !
                if ( INDEX(CHARLINE,'# class2') > 0 ) then;                              !
                                                                                         !
                    CH_DIHEDRAL_STYLE_LIBRARY(IFILE) = 'class2';                         !
                                                                                         !
                    NPARAM_DIHEDRALS_LIBRARY(IFILE) = 6;                                 ! 
                                                                                         !
                else if ( INDEX(CHARLINE,'# opls') > 0 ) then;                           !
                                                                                         !
                    CH_DIHEDRAL_STYLE_LIBRARY(IFILE) = 'opls';                           !
                                                                                         !
                    NPARAM_DIHEDRALS_LIBRARY(IFILE) = 4;                                 !
                                                                                         !
                else if ( INDEX(CHARLINE,'# charmmfsw') > 0 ) then;                      !
                                                                                         !
                    CH_DIHEDRAL_STYLE_LIBRARY(IFILE) = 'charmmfsw';                      !
                                                                                         !
                     NPARAM_DIHEDRALS_LIBRARY(IFILE) = 4;                                !
                                                                                         !
                end if                                                                   !
                                                                                         !
                write(icanal,'(a47,i4,a19)') '| The number of dihedral '//    &          !
                                             'parameters to read is ',        &          !
                                             NPARAM_DIHEDRALS_LIBRARY(IFILE), &          !
                                             REPEAT(' ',18)//'|';                        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
                read(5,*);                                                               !
                                                                                         !
                do i = 1, NTYPE_DIHEDRAL_LIBRARY(IFILE);                                 !
                                                                                         !
                    IOSEF2 = NPARAM_DIHEDRALS_LIBRARY(IFILE);                            !
                                                                                         !
                    read(5,*) IOSEF1, DIHEDRAL_COEFFS_LIBRARY(1:IOSEF2,IOSEF1,IFILE);    !
                                                                                         !
!                   write(icanal,*) IOSEF1, DIHEDRAL_COEFFS_LIBRARY(1:IOSEF2,IOSEF1,IFILE);

                end do                                                                   !
                                                                                         !
                write(icanal,'(a70)') '| Dihedral coefficients were read'// &            !
                                      REPEAT(' ',36)//'|';                               !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
!       ### Read parameters for class2 potential ###################################################
                                                                                         !
        else if ( INDEX(CHARLINE,' MiddleBondTorsion Coeffs') > 0 ) then;                !
                                                                                         !
            if ( NTYPE_DIHEDRAL_LIBRARY(IFILE) > 0 ) then;                               !
                                                                                         !
                read(5,*);

                do i = 1, NTYPE_DIHEDRAL_LIBRARY(IFILE);

                   read(5,*) IOSEF1, MIDDLEBONDTORSION_COEFFS_LIBRARY(1:4,IOSEF1,IFILE);

                end do

                write(icanal,'(a70)') '| MIDDLEBONDTORSION COEFFICIENTS WERE READ'// &
                                      REPEAT(' ',27)//'|';

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            end if

        else if ( INDEX(CHARLINE,' EndBondTorsion Coeffs') > 0 ) then;                   !

            if ( NTYPE_DIHEDRAL_LIBRARY(IFILE) > 0 ) then;

                read(5,*);

                do i = 1, NTYPE_DIHEDRAL_LIBRARY(IFILE);

                    read(5,*) IOSEF1, ENDBONDTORSION_COEFFS_LIBRARY(1:8,IOSEF1,IFILE);

                end do

                write(icanal,'(a70)') '| ENDBONDTORSION COEFFICIENTS WERE READ'// &
                                      REPEAT(' ',30)//'|';

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

            end if

        else if ( INDEX(CHARLINE,' AngleTorsion Coeffs # class2') > 0 ) then

            if ( NTYPE_DIHEDRAL_LIBRARY(IFILE) > 0 ) then;

                read(5,*);

                do i = 1, NTYPE_DIHEDRAL_LIBRARY(IFILE);

                    read(5,*) IOSEF1, ANGLETORSION_COEFFS_LIBRARY(1:8,IOSEF1,IFILE);

                end do
                                                                                         !
                write(icanal,'(a70)') '| ANGLETORSION COEFFICIENTS WERE READ'// &        !
                                      REPEAT(' ',32)//'|';                               !

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            end if                                                                       ! 
                                                                                         !
       else if ( INDEX(CHARLINE,' AngleAngleTorsion Coeffs') > 0 ) then;                 !
                                                                                         !
            if ( NTYPE_DIHEDRAL_LIBRARY(IFILE) > 0 ) then;                               !
                                                                                         !
                read(5,*);                                                               !
                                                                                         !
                do i = 1, NTYPE_DIHEDRAL_LIBRARY(IFILE);                                 !

                    read(5,*) IOSEF1,                                             &      !
                              ANGLEANGLETORSION_COEFFS_LIBRARY(1:3,IOSEF1,IFILE);        !

                end do                                                                   !
                                                                                         !
                write(icanal,'(a70)') '| ANGLEANGLETORSION COEFFICIENTS WERE READ'// &   !
                                      REPEAT(' ',27)//'|';                               !

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
        else if ( INDEX(CHARLINE,' BondBond13 Coeffs') > 0 ) then;                       !
                                                                                         !
            if ( NTYPE_DIHEDRAL_LIBRARY(IFILE) > 0 ) then;                               !
                                                                                         !
                read(5,*);                                                               !
                                                                                         !
                do i = 1, NTYPE_DIHEDRAL_LIBRARY(IFILE);                                 !

                    read(5,*) IOSEF1, BONDBOND13_COEFFS_LIBRARY(1:3,IOSEF1,IFILE);       !

                end do                                                                   !
                                                                                         !
                write(icanal,'(a70)') '| BONDBOND13 COEFFICIENTS WERE READ'// &          !
                                      REPEAT(' ',34)//'|';                               !

                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
        else if ( INDEX(CHARLINE,' Improper Coeffs') > 0 ) then;                         !
                                                                                         !
            EOF2 = 0;                                                                    !
                                                                                         !
            if ( NTYPE_IMPROPER_LIBRARY(IFILE) > 0 ) then;                               !
                                                                                         !
                read(5,*);                                                               !
                                                                                         !
                do i = 1, NTYPE_IMPROPER_LIBRARY(IFILE);                                 !
                                                                                         !
                    read(5,*,iostat=EOF2) IOSEF1,                                    &   !
                                          IMPROPER_COEFFS_LIBRARY(1:2,IOSEF1,IFILE);     !
                                                                                         !
                    if ( EOF2 /= 0 ) EXIT;                                               !
                                                                                         !
                end do                                                                   !
                                                                                         !
                write(icanal,'(a70)') '| IMPROPER COEFFICIENTS WERE READ'// &            !
                                      REPEAT(' ',36)//'|';                               !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
        else if ( INDEX(CHARLINE,' AngleAngle Coeffs') > 0 ) then;                       !
                                                                                         !
            read(5,*);                                                                   !
                                                                                         !
            if ( NTYPE_IMPROPER_LIBRARY(IFILE) > 0 ) then;                               !
                                                                                         !
                do i = 1, NTYPE_IMPROPER_LIBRARY(IFILE);                                 !
                                                                                         !
                    read(5,*) IOSEF1, ANGLEANGLE_COEFFS_LIBRARY(1:6,IOSEF1,IFILE);       !
                                                                                         !
                end do                                                                   !
                                                                                         !
                write(icanal,'(a70)') '| ANGLEANGLE COEFFICIENTS WERE READ'// &          !
                                      REPEAT(' ',34)//'|';                               !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
            end if                                                                       !
                                                                                         !
        end if                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
    close(5);                                                                            !
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine READ_INTERATOMIC_POTENTIALS_TEMPLATE
