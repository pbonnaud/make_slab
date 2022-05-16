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



subroutine ALLOCATE_CONFIG_LIBRARY_ARRAYS(icanal, &
                                          IINSREM)

!   ************************************************************************************************
!   **                                READ LAMMPS TEMPLATE FILES                                  **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                           **
!   ** ILENGTH               : LENGTH OF THE CHARACTER CHAIN CHEXT                                **
!   ** CHEXT                 : NAME OF THE LAMMPS CONFIGURATION                                   **
!   **                                                                                            **
!   ** NATOM_TPLTE           : NUMBER OF ATOMS IN THE TEMPLATE                                    **
!   ** NBOND_TPLTE           : NUMBER OF COVALENT BONDS IN THE MOLECULE DESRCIBED IN THE TEMPLATE **
!   ** NANGLE_TPLTE          : NUMBER OF ANGLES IN THE TEMPLATED MOLECULE                         **
!   ** NDIHEDRAL_TPLTE       : NUMBER OF DIHEDRAL ANGLES IN THE TEMPLATED MOLECULE                **
!   ** NIMPROPER_TPLTE       : NUMBER OF IMPROPER ANGLES IN THE TEMPLATED MOLECULE                **
!   ** TEMPLATE_TYPE         : TYPE OF ATOMS IN THE TEMPLATED MOLECULE                            **
!   ** NTYPE_TEMPLATE        : NUMBER OF ATOM TYPES IN THE TEMPLATED MOLECULE                     **
!   ** NTYPE_BOND_TPLTE      : NUMBER OF BOND TYPES IN THE TEMPLATED MOLECULE                     **
!   ** NTYPE_ANGLE_TPLTE     : NUMBER OF ANGLE TYPES IN THE TEMPLATED MOLECULE                    **
!   ** NTYPE_DIHEDRAL_TPLTE  : NUMBER OF DIHEDRAL TYPES IN THE TEMPLATED MOLECULE                 **
!   ** NTYPE_IMPROPER_TPLTE  : NUMBER OF IMPROPER TYPES IN THE TEMPLATED MOLECULE                 **
!   ** CONFIG_QI_TPLTE       : PARTIAL CHARGES OF ATOMS IN THE TEMPLATED MOLECULE                 **
!   ** CONFIG_RI_TPLTE       : COORDINATES OF ATOMS IN THE TEMPLATED MOLECULE                     **
!   ** CONFIG_ATOMID_TPLTE   : ATOM ID IN THE TEMPLATED MOLECULE                                  **
!   ** BOND_TYPE_TPLTE       : TYPE OF BONDS IN THE MOLECULAR CONFIGURATION                       **
!   ** ANGLE_TYPE_TPLTE      : TYPE OF ANGLES IN THE MOLECULAR CONFIGURATION                      **
!   ** DIHEDRAL_TYPE_TPLTE   : TYPE OF DIHEDRALS IN THE MOLECULAR CONFIGURATION                   **
!   ** IMPROPER_TYPE_TPLTE   : TYPE OF IMPROPER ANGLES IN THE MOLECULAR CONFIGURATION             **
!   ** BOND_ATOMID_TPLTE     : ATOM ID RELATED TO THE BONDS IN THE MOLECULAR CONFIGURATION        **
!   ** ANGLE_ATOMID_TPLTE    : ATOM ID RELATED TO THE ANGLES IN THE MOLECULAR CONFIGURATION       **
!   ** IMPROPER_ATOMID_TPLTE : ATOM ID RELATED TO THE IMPROPER ANGLES IN THE MOLECULAR            **
!   **                         CONFIGURATION                                                      **
!   **                                                                                            **
!   ************************************************************************************************

    use module_size_arrays;

    use module_library;

!   ************************************************************************************************

    implicit none;

!   ***********************************************************************************************

    integer (kind=4), intent(in) :: icanal, IINSREM;

!   ************************************************************************************************

    integer (kind=4) :: i, j;

!   ************************************************************************************************

    integer (kind=4) :: NFILE_LIBRARY_TMP;

!   ************************************************************************************************

    integer (kind=4) :: IOSEF_NATOM, IOSEF_NBOND, IOSEF_NANGLE, IOSEF_NDIHEDRAL, IOSEF_NIMPROPER;

    integer (kind=4) :: IOSEF_NTYPE_NATOM,    &
                        IOSEF_NTYPE_BOND,     &
                        IOSEF_NTYPE_ANGLE,    &
                        IOSEF_NTYPE_DIHEDRAL, &
                        IOSEF_NTYPE_IMPROPER;

    character (len=250) :: CHEXT;

!   ************************************************************************************************

    character (len=250) :: CHARLINE;

    integer (kind=4) :: EOF;

    integer (kind=4) :: IOSEF1, IOSEF2, IOSEF3, IOSEF4;

    character (len=10) :: CHOSEF1, CHOSEF2, CHOSEF3;

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************

    CHTITLE = 'ALLOCATE CONFIG LIBRARY ARRAYS';

    ILENGTH_TITLE = LEN_TRIM(CHTITLE);

    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));


    if ( IINSREM == 1 ) then;
        NFILE_LIBRARY_TMP = NFILE_LIBRARY;
    else
        NFILE_LIBRARY_TMP = NFILE_LIBRARY_REMOVE;
    end if

    write(icanal,*) '| NFILE_LIBRARY        : ', NFILE_LIBRARY;

    write(icanal,*) '| NFILE_LIBRARY_REMOVE : ', NFILE_LIBRARY_REMOVE;


    allocate(NATOM_LIBRARY(1:NFILE_LIBRARY_TMP));
    allocate(NBOND_LIBRARY(1:NFILE_LIBRARY_TMP));
    allocate(NANGLE_LIBRARY(1:NFILE_LIBRARY_TMP));
    allocate(NDIHEDRAL_LIBRARY(1:NFILE_LIBRARY_TMP));
    allocate(NIMPROPER_LIBRARY(1:NFILE_LIBRARY_TMP));
    allocate(NTYPE_ATOM_LIBRARY(1:NFILE_LIBRARY_TMP));
    allocate(NTYPE_BOND_LIBRARY(1:NFILE_LIBRARY_TMP));
    allocate(NTYPE_ANGLE_LIBRARY(1:NFILE_LIBRARY_TMP));
    allocate(NTYPE_DIHEDRAL_LIBRARY(1:NFILE_LIBRARY_TMP));
    allocate(NTYPE_IMPROPER_LIBRARY(1:NFILE_LIBRARY_TMP));

    allocate(IPOTENTIAL_CLASS2_LIBRARY(1:NFILE_LIBRARY_TMP));

    IMAX_NATOM     = 0;
    IMAX_NBOND     = 0;
    IMAX_NANGLE    = 0;
    IMAX_NDIHEDRAL = 0;
    IMAX_NIMPROPER = 0;

    IMAX_NTYPE_NATOM    = 0;
    IMAX_NTYPE_BOND     = 0;
    IMAX_NTYPE_ANGLE    = 0;
    IMAX_NTYPE_DIHEDRAL = 0;
    IMAX_NTYPE_IMPROPER = 0;

    do i = 1, NFILE_LIBRARY_TMP; ! LOOP OVER THE NUMBER OF FILES TO READ IN THE LIBRARY

        if ( IINSREM == 1 ) then; 

            CHEXT = TRIM(CHNAME_FILE_LIBRARY(1,i)); 

        else

            CHEXT = TRIM(CHNAME_FILE_LIBRARY_REMOVE(i)); 

        end if

        IOSEF1 = LEN_TRIM(CHEXT);

        IOSEF2 = 67 - IOSEF1; 

        write(icanal,'(a70)') '| '//TRIM(CHEXT)//REPEAT(' ',IOSEF2)//'|';

        open(2,file=TRIM(CHEXT));

        EOF = 0;

        IOSEF_NATOM     = 0;
        IOSEF_NBOND     = 0;
        IOSEF_NANGLE    = 0;
        IOSEF_NDIHEDRAL = 0;
        IOSEF_NIMPROPER = 0;

        IOSEF_NTYPE_NATOM    = 0;
        IOSEF_NTYPE_BOND     = 0;
        IOSEF_NTYPE_ANGLE    = 0;
        IOSEF_NTYPE_DIHEDRAL = 0;
        IOSEF_NTYPE_IMPROPER = 0;

        do
            read(2,'(a)',iostat=EOF) CHARLINE;

            if ( EOF /= 0 ) EXIT;

            if ( INDEX(CHARLINE,'atoms') > 0 ) then
                read(CHARLINE,*) IOSEF_NATOM;
                if ( IOSEF_NATOM > IMAX_NATOM ) IMAX_NATOM = IOSEF_NATOM;



            else if ( INDEX(CHARLINE,'Atoms') > 0 ) then
                read(2,*);

                do j = 1, IOSEF_NATOM;
                    read(2,*) IOSEF1, IOSEF2;
                    if ( IOSEF2 > IOSEF_NTYPE_NATOM ) IOSEF_NTYPE_NATOM = IOSEF2;
                end do

                if ( IOSEF_NTYPE_NATOM > IMAX_NTYPE_NATOM ) IMAX_NTYPE_NATOM = IOSEF_NTYPE_NATOM;

!           else if ( INDEX(CHARLINE,'bonds') > 0 ) then
!               read(CHARLINE,*) IOSEF_NBOND;
!               if ( IOSEF_NBOND > IMAX_NBOND ) IMAX_NBOND = IOSEF_NBOND;

!           else if ( INDEX(CHARLINE,'angles') > 0 ) then
!               read(CHARLINE,*) IOSEF_NANGLE;
!               if ( IOSEF_NANGLE > IMAX_NANGLE ) IMAX_NANGLE = IOSEF_NANGLE;

!           else if ( INDEX(CHARLINE,'dihedrals') > 0 ) then
!               read(CHARLINE,*) IOSEF_NDIHEDRAL;
!               if ( IOSEF_NDIHEDRAL > IMAX_NDIHEDRAL ) IMAX_NDIHEDRAL = IOSEF_NDIHEDRAL;

!           else if ( INDEX(CHARLINE,'impropers') > 0 ) then
!               read(CHARLINE,*) IOSEF_NIMPROPER;
!               if ( IOSEF_NIMPROPER > IMAX_NIMPROPER ) IMAX_NIMPROPER = IOSEF_NIMPROPER;

!           else if ( INDEX(CHARLINE,'Types') > 0 ) then
!               read(2,*);
!               do j = 1, IOSEF_NATOM;
!                   read(2,*) IOSEF1, IOSEF2;
!                   if ( IOSEF2 > IOSEF_NTYPE_NATOM ) IOSEF_NTYPE_NATOM = IOSEF2;
!               end do

!               if ( IOSEF_NTYPE_NATOM > IMAX_NTYPE_NATOM ) IMAX_NTYPE_NATOM = IOSEF_NTYPE_NATOM;

!           else if ( INDEX(CHARLINE,'Bonds') > 0 ) then
!               read(2,*);
!               do j = 1, IOSEF_NBOND;
!                   read(2,*) IOSEF1, IOSEF2;
!                   if ( IOSEF2 > IOSEF_NTYPE_BOND ) IOSEF_NTYPE_BOND = IOSEF2;
!               end do

!               if ( IOSEF_NTYPE_BOND > IMAX_NTYPE_BOND ) IMAX_NTYPE_BOND = IOSEF_NTYPE_BOND;

!           else if ( INDEX(CHARLINE,'Angles') > 0 ) then
!               read(2,*);
!               do j = 1, IOSEF_NANGLE;
!                   read(2,*) IOSEF1, IOSEF2;
!                   if ( IOSEF2 > IOSEF_NTYPE_ANGLE ) IOSEF_NTYPE_ANGLE = IOSEF2;
!               end do

!               if ( IOSEF_NTYPE_ANGLE > IMAX_NTYPE_ANGLE ) IMAX_NTYPE_ANGLE = IOSEF_NTYPE_ANGLE;

!           else if ( INDEX(CHARLINE,'Dihedrals') > 0 ) then;
!               read(2,*);
!               do j = 1, IOSEF_NDIHEDRAL;
!                   read(2,*) IOSEF1, IOSEF2;
!                   if ( IOSEF2 > IOSEF_NTYPE_DIHEDRAL ) IOSEF_NTYPE_DIHEDRAL = IOSEF2;
!               end do

!               if ( IOSEF_NTYPE_DIHEDRAL > IMAX_NTYPE_DIHEDRAL ) IMAX_NTYPE_DIHEDRAL = IOSEF_NTYPE_DIHEDRAL;

!           else if ( INDEX(CHARLINE,'Impropers') > 0 ) then
!               read(2,*);
!               do j = 1, IOSEF_NIMPROPER;
!                   read(2,*) IOSEF1, IOSEF2;
!                   if ( IOSEF2 > IOSEF_NTYPE_IMPROPER ) IOSEF_NTYPE_IMPROPER = IOSEF2;
!               end do

!               if ( IOSEF_NTYPE_IMPROPER > IMAX_NTYPE_IMPROPER ) IMAX_NTYPE_IMPROPER = IOSEF_NTYPE_IMPROPER;

            end if
        end do

        close(2);

        write(icanal,'(a25,i8,a37)') '| IOSEF_NATOM          : ', IOSEF_NATOM,     REPEAT(' ',36)//'|';
        write(icanal,'(a25,i8,a37)') '| IOSEF_NBOND          : ', IOSEF_NBOND,     REPEAT(' ',36)//'|';
        write(icanal,'(a25,i8,a37)') '| IOSEF_NANGLE         : ', IOSEF_NANGLE,    REPEAT(' ',36)//'|';
        write(icanal,'(a25,i8,a37)') '| IOSEF_NDIHEDRAL      : ', IOSEF_NDIHEDRAL, REPEAT(' ',36)//'|';
        write(icanal,'(a25,i8,a37)') '| IOSEF_NIMPROPER      : ', IOSEF_NIMPROPER, REPEAT(' ',36)//'|';
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
        write(icanal,'(a25,i8,a37)') '| IOSEF_NTYPE_NATOM    : ', IOSEF_NTYPE_NATOM,    REPEAT(' ',36)//'|';
        write(icanal,'(a25,i8,a37)') '| IOSEF_NTYPE_BOND     : ', IOSEF_NTYPE_BOND,     REPEAT(' ',36)//'|';
        write(icanal,'(a25,i8,a37)') '| IOSEF_NTYPE_ANGLE    : ', IOSEF_NTYPE_ANGLE,    REPEAT(' ',36)//'|';
        write(icanal,'(a25,i8,a37)') '| IOSEF_NTYPE_DIHEDRAL : ', IOSEF_NTYPE_DIHEDRAL, REPEAT(' ',36)//'|';
        write(icanal,'(a25,i8,a37)') '| IOSEF_NTYPE_IMPROPER : ', IOSEF_NTYPE_IMPROPER, REPEAT(' ',36)//'|';
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
    end do

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
    write(icanal,'(a25,i8,a37)') '| IMAX_NATOM           : ', IMAX_NATOM,     REPEAT(' ',36)//'|';

    allocate(CONFIG_ATOMID_LIBRARY(1:IMAX_NATOM,1:NFILE_LIBRARY_TMP));
    allocate(CONFIG_ATOM_TYPE_LIBRARY(1:IMAX_NATOM,1:NFILE_LIBRARY_TMP));
    allocate(CONFIG_QI_LIBRARY(1:IMAX_NATOM,1:NFILE_LIBRARY_TMP));
    allocate(CONFIG_VI_LIBRARY(1:3,1:IMAX_NATOM,1:NFILE_LIBRARY_TMP));
    allocate(CONFIG_RI_LIBRARY(1:3,1:IMAX_NATOM,1:NFILE_LIBRARY_TMP));

    allocate(CONFIG_NAT_LIBRARY(1:IMAX_NATOM,1:NFILE_LIBRARY_TMP));

!   write(icanal,'(a25,i8,a37)') '| IMAX_NBOND           : ', IMAX_NBOND,     REPEAT(' ',36)//'|';

!   allocate(BOND_TYPE_LIBRARY(1:IMAX_NBOND,1:NFILE_LIBRARY));
!   allocate(BOND_ATOMID_LIBRARY(1:2,1:IMAX_NBOND,1:NFILE_LIBRARY));

!   write(icanal,'(a25,i8,a37)') '| IMAX_NANGLE          : ', IMAX_NANGLE,    REPEAT(' ',36)//'|';

!   allocate(ANGLE_TYPE_LIBRARY(1:IMAX_NANGLE,1:NFILE_LIBRARY));
!   allocate(ANGLE_ATOMID_LIBRARY(1:3,1:IMAX_NANGLE,1:NFILE_LIBRARY));

!   write(icanal,'(a25,i8,a37)') '| IMAX_NDIHEDRAL       : ', IMAX_NDIHEDRAL, REPEAT(' ',36)//'|';

!   allocate(DIHEDRAL_TYPE_LIBRARY(1:IMAX_NDIHEDRAL,1:NFILE_LIBRARY));
!   allocate(DIHEDRAL_ATOMID_LIBRARY(1:4,1:IMAX_NDIHEDRAL,1:NFILE_LIBRARY));

!   write(icanal,'(a25,i8,a37)') '| IMAX_NIMPROPER       : ', IMAX_NIMPROPER, REPEAT(' ',36)//'|';

!   allocate(IMPROPER_TYPE_LIBRARY(1:IMAX_NIMPROPER,1:NFILE_LIBRARY));
!   allocate(IMPROPER_ATOMID_LIBRARY(1:4,1:IMAX_NIMPROPER,1:NFILE_LIBRARY));

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
    write(icanal,'(a25,i8,a37)') '| IMAX_NTYPE_NATOM     : ', IMAX_NTYPE_NATOM,    REPEAT(' ',36)//'|';

    allocate(ATOM_MASSES_LIBRARY(1:IMAX_NTYPE_NATOM,1:NFILE_LIBRARY_TMP));
    allocate(ATOM_LABEL_LIBRARY(1:IMAX_NTYPE_NATOM,1:NFILE_LIBRARY_TMP));
!   allocate(POTENTIAL_CLASS2_LIBRARY(1:2,1:IMAX_NTYPE_NATOM,1:NFILE_LIBRARY));

!   write(icanal,'(a25,i8,a37)') '| IMAX_NTYPE_BOND      : ', IMAX_NTYPE_BOND,     REPEAT(' ',36)//'|'; 

!   allocate(BOND_COEFFS_LIBRARY(1:4,1:IMAX_NTYPE_BOND,1:NFILE_LIBRARY));

!   write(icanal,'(a25,i8,a37)') '| IMAX_NTYPE_ANGLE     : ', IMAX_NTYPE_ANGLE,    REPEAT(' ',36)//'|';

!   allocate(ANGLE_COEFFS_LIBRARY(1:4,1:IMAX_NTYPE_ANGLE,1:NFILE_LIBRARY));
!   allocate(BONDBOND_COEFFS_LIBRARY(1:3,1:IMAX_NTYPE_ANGLE,1:NFILE_LIBRARY));
!   allocate(BONDANGLE_COEFFS_LIBRARY(1:4,1:IMAX_NTYPE_ANGLE,1:NFILE_LIBRARY));

!   write(icanal,'(a25,i8,a37)') '| IMAX_NTYPE_DIHEDRAL  : ', IMAX_NTYPE_DIHEDRAL, REPEAT(' ',36)//'|';

!   allocate(DIHEDRAL_COEFFS_LIBRARY(1:6,1:IMAX_NTYPE_DIHEDRAL,1:NFILE_LIBRARY));
!   allocate(MIDDLEBONDTORSION_COEFFS_LIBRARY(1:4,1:IMAX_NTYPE_DIHEDRAL,1:NFILE_LIBRARY));
!   allocate(ENDBONDTORSION_COEFFS_LIBRARY(1:8,1:IMAX_NTYPE_DIHEDRAL,1:NFILE_LIBRARY));
!   allocate(ANGLETORSION_COEFFS_LIBRARY(1:8,1:IMAX_NTYPE_DIHEDRAL,1:NFILE_LIBRARY));
!   allocate(ANGLEANGLETORSION_COEFFS_LIBRARY(1:3,1:IMAX_NTYPE_DIHEDRAL,1:NFILE_LIBRARY));
!   allocate(BONDBOND13_COEFFS_LIBRARY(1:3,1:IMAX_NTYPE_DIHEDRAL,1:NFILE_LIBRARY));

!   write(icanal,'(a25,i8,a37)') '| IMAX_NTYPE_IMPROPER  : ', IMAX_NTYPE_IMPROPER, REPEAT(' ',36)//'|';

!   allocate(IMPROPER_COEFFS_LIBRARY(1:2,1:IMAX_NTYPE_IMPROPER,1:NFILE_LIBRARY));
!   allocate(ANGLEANGLE_COEFFS_LIBRARY(1:6,1:IMAX_NTYPE_IMPROPER,1:NFILE_LIBRARY));

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';

end subroutine ALLOCATE_CONFIG_LIBRARY_ARRAYS
