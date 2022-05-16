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

subroutine READ_LAMMPS_CONFIG_LIBRARY(icanal,IFILE)

!   ************************************************************************************************
!   **                                READ LAMMPS TEMPLATE FILES                                  **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                           **
!   **                                                                                            **
!   ************************************************************************************************

    use module_size_arrays;

    use module_library;

!   ************************************************************************************************

    implicit none;

!   ***********************************************************************************************

    integer (kind=4), intent(in) :: icanal, IFILE;

!   ************************************************************************************************

    character (len=250) :: CHEXT;

    integer (kind=4) :: i, j;

    character (len=250) :: CHARLINE;

    integer (kind=4) :: EOF;

    integer (kind=4) :: IOSEF1, IOSEF2, IOSEF3, IOSEF4;

    character (len=10) :: CHOSEF1, CHOSEF2, CHOSEF3;

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
    CHTITLE = 'Read LAMMPS config library';                                              !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Get the name of the file to read ###########################################################

    CHEXT = TRIM(CHNAME_FILE_LIBRARY(1,IFILE));

    IOSEF1 = LEN_TRIM(CHEXT);

    IOSEF2 = 67 - IOSEF1; 

    write(icanal,'(a70)') '| '//TRIM(CHEXT)//REPEAT(' ',IOSEF2)//'|';
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

    NATOM_LIBRARY(IFILE)     = 0;
!   NBOND_LIBRARY(IFILE)     = 0; 
!   NANGLE_LIBRARY(IFILE)    = 0;
!   NDIHEDRAL_LIBRARY(IFILE) = 0;
!   NIMPROPER_LIBRARY(IFILE) = 0;

    NTYPE_ATOM_LIBRARY(IFILE)     = 0;
!   NTYPE_BOND_LIBRARY(IFILE)     = 0;
!   NTYPE_ANGLE_LIBRARY(IFILE)    = 0;
!   NTYPE_DIHEDRAL_LIBRARY(IFILE) = 0;
!   NTYPE_IMPROPER_LIBRARY(IFILE) = 0;

    CONFIG_ATOM_TYPE_LIBRARY(1:IMAX_NATOM,IFILE)  = 0;
!   BOND_TYPE_LIBRARY(1:IMAX_NBOND,IFILE)         = 0;
!   ANGLE_TYPE_LIBRARY(1:IMAX_NANGLE,IFILE)       = 0;
!   DIHEDRAL_TYPE_LIBRARY(1:IMAX_NDIHEDRAL,IFILE) = 0;
!   IMPROPER_TYPE_LIBRARY(1:IMAX_NIMPROPER,IFILE) = 0;

    CONFIG_ATOMID_LIBRARY(1:IMAX_NATOM,IFILE)           = 0;
!   BOND_ATOMID_LIBRARY(1:2,1:IMAX_NBOND,IFILE)         = 0;
!   ANGLE_ATOMID_LIBRARY(1:3,1:IMAX_NANGLE,IFILE)       = 0;
!   DIHEDRAL_ATOMID_LIBRARY(1:4,1:IMAX_NDIHEDRAL,IFILE) = 0;
!   IMPROPER_ATOMID_LIBRARY(1:4,1:IMAX_NIMPROPER,IFILE) = 0;

    IPOTENTIAL_CLASS2_LIBRARY(IFILE) = 0;


    CONFIG_QI_LIBRARY(1:IMAX_NATOM,IFILE) = 0.0d0;

    CONFIG_VI_LIBRARY(1:3,1:IMAX_NATOM,IFILE) = 0.0d0;

    CONFIG_RI_LIBRARY(1:3,1:IMAX_NATOM,IFILE) = 0.0d0;

    open(2,file=TRIM(CHEXT));

    EOF = 0;

    do
        read(2,'(a)',iostat=EOF) CHARLINE;

        if ( EOF /= 0 ) EXIT;

        if ( INDEX(CHARLINE,'atoms') > 0 ) then
            read(CHARLINE,*) NATOM_LIBRARY(IFILE);

            write(icanal,*) 'OK 1';

        else if ( INDEX(CHARLINE,'atom types') > 0 ) then
            read(CHARLINE,*) NTYPE_ATOM_LIBRARY(IFILE);

            write(icanal,*) 'OK 2';

        else if ( INDEX(CHARLINE,'Masses') > 0 ) then
            write(icanal,*) 'OK 3', NTYPE_ATOM_LIBRARY(IFILE);

            read(2,*);
            do i = 1, NTYPE_ATOM_LIBRARY(IFILE);  
                read(2,*) IOSEF1, ATOM_MASSES_LIBRARY(IOSEF1,IFILE);
            end do

            write(icanal,*) 'OK 4';

            do i = 1, NTYPE_ATOM_LIBRARY(IFILE);   !      NTYPE_TEMPLATE
                if ( ABS( ATOM_MASSES_LIBRARY(i,IFILE) -  1.00794d0 ) < 0.1d0 ) ATOM_LABEL_LIBRARY(i,IFILE) = 'H';
                if ( ABS( ATOM_MASSES_LIBRARY(i,IFILE) - 12.01070d0 ) < 0.1d0 ) ATOM_LABEL_LIBRARY(i,IFILE) = 'C';
                if ( ABS( ATOM_MASSES_LIBRARY(i,IFILE) - 14.00670d0 ) < 0.1d0 ) ATOM_LABEL_LIBRARY(i,IFILE) = 'N';
                if ( ABS( ATOM_MASSES_LIBRARY(i,IFILE) - 15.99940d0 ) < 0.1d0 ) ATOM_LABEL_LIBRARY(i,IFILE) = 'O';
                if ( ABS( ATOM_MASSES_LIBRARY(i,IFILE) - 28.08550d0 ) < 0.1d0 ) ATOM_LABEL_LIBRARY(i,IFILE) = 'Si';
                if ( ABS( ATOM_MASSES_LIBRARY(i,IFILE) - 32.06400d0 ) < 0.1d0 ) ATOM_LABEL_LIBRARY(i,IFILE) = 'S';
                if ( ABS( ATOM_MASSES_LIBRARY(i,IFILE) - 40.07800d0 ) < 0.1d0 ) ATOM_LABEL_LIBRARY(i,IFILE) = 'Ca';
                if ( ABS( ATOM_MASSES_LIBRARY(i,IFILE) - 55.84500d0 ) < 0.1d0 ) ATOM_LABEL_LIBRARY(i,IFILE) = 'Fe';
            end do

            write(icanal,*) 'MASSES WERE READ';

!       else if ( INDEX(CHARLINE,'bonds') > 0 ) then
!           read(CHARLINE,*) NBOND_LIBRARY(IFILE);

!       else if ( INDEX(CHARLINE,'angles') > 0 ) then
!           read(CHARLINE,*) NANGLE_LIBRARY(IFILE);

!       else if ( INDEX(CHARLINE,'dihedrals') > 0 ) then
!           read(CHARLINE,*) NDIHEDRAL_LIBRARY(IFILE);

!       else if ( INDEX(CHARLINE,'impropers') > 0 ) then
!           read(CHARLINE,*) NIMPROPER_LIBRARY(IFILE);

!       else if ( INDEX(CHARLINE,'Charges') > 0 ) then
!           read(2,*);
!           do i = 1, NATOM_LIBRARY(IFILE);
!               read(2,*) IOSEF1, CONFIG_QI_LIBRARY(IOSEF1,IFILE);
!           end do

!           write(icanal,'(a70)') '| Charges were read'//REPEAT(' ',50)//'|';

        else if ( INDEX(CHARLINE,'Atoms') > 0 ) then
            read(2,*);
            do i = 1, NATOM_LIBRARY(IFILE);
                read(2,*) CONFIG_ATOMID_LIBRARY(i,IFILE),                                 &
                          CONFIG_ATOM_TYPE_LIBRARY(CONFIG_ATOMID_LIBRARY(i,IFILE),IFILE), &
                          CONFIG_QI_LIBRARY(CONFIG_ATOMID_LIBRARY(i,IFILE),IFILE),        &
                          CONFIG_RI_LIBRARY(1:3,CONFIG_ATOMID_LIBRARY(i,IFILE),IFILE);
            end do

!           write(icanal,'(a70)') '| Coords were read'//REPEAT(' ',51)//'|';

        else if ( INDEX(CHARLINE,'Velocities') > 0 ) then
            read(2,*);
            do i = 1, NATOM_LIBRARY(IFILE);
                read(2,*) IOSEF1, CONFIG_VI_LIBRARY(1:3,IOSEF1,IFILE);
            end do

!       else if ( INDEX(CHARLINE,'Types') > 0 ) then
!           read(2,*);
!           do i = 1, NATOM_LIBRARY(IFILE);
!               read(2,*) IOSEF1, CONFIG_ATOM_TYPE_LIBRARY(IOSEF1,IFILE);

!               if ( CONFIG_ATOM_TYPE_LIBRARY(IOSEF1,IFILE) > NTYPE_ATOM_LIBRARY(IFILE) ) NTYPE_ATOM_LIBRARY(IFILE) = CONFIG_ATOM_TYPE_LIBRARY(IOSEF1,IFILE);
!           end do

!           write(icanal,'(a70)') '| Types were read'//REPEAT(' ',52)//'|';

!       else if ( INDEX(CHARLINE,'Bonds') > 0 ) then
!           read(2,*);
!           do i = 1, NBOND_LIBRARY(IFILE);
!               read(2,*) IOSEF1, BOND_TYPE_LIBRARY(IOSEF1,IFILE), BOND_ATOMID_LIBRARY(1:2,IOSEF1,IFILE)
!               if ( BOND_TYPE_LIBRARY(IOSEF1,IFILE) > NTYPE_BOND_LIBRARY(IFILE) ) NTYPE_BOND_LIBRARY(IFILE) = BOND_TYPE_LIBRARY(IOSEF1,IFILE);
!           end do

!           write(icanal,'(a70)') '| Bonds were read'//REPEAT(' ',52)//'|';

!       else if ( INDEX(CHARLINE,'Angles') > 0 ) then
!           read(2,*);
!           do i = 1, NANGLE_LIBRARY(IFILE);
!               read(2,*) IOSEF1, ANGLE_TYPE_LIBRARY(IOSEF1,IFILE), ANGLE_ATOMID_LIBRARY(1:3,IOSEF1,IFILE)
!               if ( ANGLE_TYPE_LIBRARY(IOSEF1,IFILE) > NTYPE_ANGLE_LIBRARY(IFILE) ) NTYPE_ANGLE_LIBRARY(IFILE) = ANGLE_TYPE_LIBRARY(IOSEF1,IFILE);
!           end do

!           write(icanal,'(a70)') '| Angles were read'//REPEAT(' ',51)//'|';

!        else if ( INDEX(CHARLINE,'Dihedrals') > 0 ) then;
!           read(2,*);
!           do i = 1, NDIHEDRAL_LIBRARY(IFILE);
!               read(2,*) IOSEF1, DIHEDRAL_TYPE_LIBRARY(IOSEF1,IFILE), DIHEDRAL_ATOMID_LIBRARY(1:4,IOSEF1,IFILE);
!               if ( DIHEDRAL_TYPE_LIBRARY(IOSEF1,IFILE) > NTYPE_DIHEDRAL_LIBRARY(IFILE) ) NTYPE_DIHEDRAL_LIBRARY(IFILE) = DIHEDRAL_TYPE_LIBRARY(IOSEF1,IFILE);
!           end do

!           write(icanal,'(a70)') '| Dihedrals were read'//REPEAT(' ',48)//'|';

!        else if ( INDEX(CHARLINE,'Impropers') > 0 ) then
!           read(2,*);
!           do i = 1, NIMPROPER_LIBRARY(IFILE);
!               read(2,*) IOSEF1, IMPROPER_TYPE_LIBRARY(IOSEF1,IFILE), IMPROPER_ATOMID_LIBRARY(1:4,IOSEF1,IFILE);
!               if ( IMPROPER_TYPE_LIBRARY(IOSEF1,IFILE) > NTYPE_IMPROPER_LIBRARY(IFILE) ) NTYPE_IMPROPER_LIBRARY(IFILE) = IMPROPER_TYPE_LIBRARY(IOSEF1,IFILE);
!           end do

!           write(icanal,'(a70)') '| Impropers were read'//REPEAT(' ',48)//'|';
        end if
    end do

    close(2);

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
    write(icanal,'(a25,i8,a37)') '| NATOM_TPLTE          : ', NATOM_LIBRARY(IFILE),     REPEAT(' ',36)//'|';
!   write(icanal,'(a25,i8,a37)') '| NBOND_TPLTE          : ', NBOND_LIBRARY(IFILE),     REPEAT(' ',36)//'|';
!   write(icanal,'(a25,i8,a37)') '| NANGLE_TPLTE         : ', NANGLE_LIBRARY(IFILE),    REPEAT(' ',36)//'|';
!   write(icanal,'(a25,i8,a37)') '| NDIHEDRAL_TPLTE      : ', NDIHEDRAL_LIBRARY(IFILE), REPEAT(' ',36)//'|';
!   write(icanal,'(a25,i8,a37)') '| NIMPROPER_TPLTE      : ', NIMPROPER_LIBRARY(IFILE), REPEAT(' ',36)//'|';
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';

    write(icanal,'(a25,i8,a37)') '| NTYPE_TEMPLATE       : ', NTYPE_ATOM_LIBRARY(IFILE),     REPEAT(' ',36)//'|';
!   write(icanal,'(a25,i8,a37)') '| NTYPE_BOND_TPLTE     : ', NTYPE_BOND_LIBRARY(IFILE),     REPEAT(' ',36)//'|'; 
!   write(icanal,'(a25,i8,a37)') '| NTYPE_ANGLE_TPLTE    : ', NTYPE_ANGLE_LIBRARY(IFILE),    REPEAT(' ',36)//'|';
!   write(icanal,'(a25,i8,a37)') '| NTYPE_DIHEDRAL_TPLTE : ', NTYPE_DIHEDRAL_LIBRARY(IFILE), REPEAT(' ',36)//'|';
!   write(icanal,'(a25,i8,a37)') '| NTYPE_IMPROPER_TPLTE : ', NTYPE_IMPROPER_LIBRARY(IFILE), REPEAT(' ',36)//'|';
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';

end subroutine READ_LAMMPS_CONFIG_LIBRARY
