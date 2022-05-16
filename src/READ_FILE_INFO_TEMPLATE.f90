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

subroutine READ_FILE_INFO_TEMPLATE(icanal,IFILE)

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

    integer (kind=4), intent(in) :: icanal, IFILE;

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
    CHTITLE = 'Read info file template';                                                 !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the position in the array of the library file to read ################################
                                                                                         !
    write(icanal,'(a10,i8,a52)') '| IFILE : ', IFILE, REPEAT(' ',51)//'|'                !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Build and write the name of the file to consider ###########################################
                                                                                         !
    CHEXT = TRIM(CHNAME_FILE_LIBRARY(3,IFILE));                                          !
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
!   ### Initialization of arrays for parameter styles ##############################################
                                                                                         !
    CH_SPECIAL_BONDS_LIBRARY(IFILE) = 'XXX';                                             !
                                                                                         !
    CH_PAIR_MODIFY_LIBRARY(IFILE) = 'XXX';                                               !
                                                                                         !
    CH_KSPACE_STYLE_LIBRARY(IFILE) = 'XXX';                                              !
                                                                                         !
!   ### Read the info file #########################################################################
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
        if ( INDEX(CHARLINE,'units') > 0 ) then;                                         !
                                                                                         !
            read(CHARLINE,*) CHOSEF1, CH_UNITS_STYLE_LIBRARY(IFILE);                     !
                                                                                         !
            IOSEF1 = 70 - 15 - 1 - LEN_TRIM(CH_UNITS_STYLE_LIBRARY(IFILE));              !
                                                                                         !
            if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                !
                                                                                         !
            write(icanal,'(a70)') '| Unit style : '//                   &                !
                                  TRIM(CH_UNITS_STYLE_LIBRARY(IFILE))// &                !
                                  REPEAT(' ',IOSEF1)//'|';                               !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHARLINE,'atom_style') > 0 ) then;                               !
                                                                                         !
            read(CHARLINE,*) CHOSEF1, CH_ATOM_STYLE_LIBRARY(IFILE);                      !
                                                                                         !
            IOSEF1 = 70 - 15 - 1 - LEN_TRIM(CH_ATOM_STYLE_LIBRARY(IFILE));               !
                                                                                         !
            if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                !
                                                                                         !
            write(icanal,'(a70)') '| Atom style : '//                  &                 !
                                  TRIM(CH_ATOM_STYLE_LIBRARY(IFILE))// &                 !
                                  REPEAT(' ',IOSEF1)//'|';                               !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHARLINE,'pair_style') > 0 ) then;                               !  
                                                                                         !
            if ( INDEX(CHARLINE,'lj/charmmfsw/coul/long') > 0 ) then;                    !
                                                                                         !
                POTENTIAL_CLASS2_CHTYPE_LIBRARY(IFILE) = 'lj/charmmfsw/coul/long';       !
                                                                                         !
            end if                                                                       !
                                                                                         !
            IOSEF1 = 70 - 15 - 1 - LEN_TRIM(POTENTIAL_CLASS2_CHTYPE_LIBRARY(IFILE));     !
                                                                                         !
            if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                !
                                                                                         !
            write(icanal,'(a70)') '| Pair style : '//                            &       !
                                  TRIM(POTENTIAL_CLASS2_CHTYPE_LIBRARY(IFILE))// &       !
                                  REPEAT(' ',IOSEF1)//'|';                               !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHARLINE,'bond_style') > 0 ) then;                               !
                                                                                         !
            read(CHARLINE,*) CHOSEF1, CH_BOND_STYLE_LIBRARY(IFILE);                      !
                                                                                         !
            IOSEF1 = 70 - 15 - 1 - LEN_TRIM(CH_BOND_STYLE_LIBRARY(IFILE));               !
                                                                                         !
            if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                !
                                                                                         !
            write(icanal,'(a70)') '| Bond style : '//                  &                 !
                                  TRIM(CH_BOND_STYLE_LIBRARY(IFILE))// &                 !
                                  REPEAT(' ',IOSEF1)//'|';                               !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHARLINE,'angle_style') > 0 ) then;                              !
                                                                                         !
            read(CHARLINE,*) CHOSEF1, CH_ANGLE_STYLE_LIBRARY(IFILE);                     !
                                                                                         !
            IOSEF1 = 70 - 16 - 1 - LEN_TRIM(CH_ANGLE_STYLE_LIBRARY(IFILE));              !
                                                                                         !
            if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                !
                                                                                         !
            write(icanal,'(a70)') '| Angle style : '//                  &                !
                                  TRIM(CH_ANGLE_STYLE_LIBRARY(IFILE))// &                !
                                  REPEAT(' ',IOSEF1)//'|';                               !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHARLINE,'dihedral_style') > 0 ) then;                           !
                                                                                         !
            read(CHARLINE,*) CHOSEF1, CH_DIHEDRAL_STYLE_LIBRARY(IFILE);                  !
                                                                                         !
            IOSEF1 = 70 - 19 - 1 - LEN_TRIM(CH_DIHEDRAL_STYLE_LIBRARY(IFILE));           !
                                                                                         !
            if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                !
                                                                                         !
            write(icanal,'(a70)') '| Dihedral style : '//                  &             !
                                  TRIM(CH_DIHEDRAL_STYLE_LIBRARY(IFILE))// &             !
                                  REPEAT(' ',IOSEF1)//'|';                               !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHARLINE,'special_bonds') > 0 ) then;                            !
                                                                                         !
            read(CHARLINE,*) CHOSEF1, CH_SPECIAL_BONDS_LIBRARY(IFILE);                   !
                                                                                         !
            IOSEF1 = 70 - 18 - 1 - LEN_TRIM(CH_SPECIAL_BONDS_LIBRARY(IFILE));            !
                                                                                         !
            if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                !
                                                                                         !
            write(icanal,'(a70)') '| Special bonds : '//                  &              !
                                  TRIM(CH_SPECIAL_BONDS_LIBRARY(IFILE))// &              !
                                  REPEAT(' ',IOSEF1)//'|';                               !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHARLINE,'pair_modify') > 0 ) then;                              !
                                                                                         !
            if ( INDEX(CHARLINE,'mix') > 0 ) then;                                       !
                                                                                         !
                CH_PAIR_MODIFY_LIBRARY(IFILE) = 'mix';                                   !
                                                                                         !
            end if                                                                       !
                                                                                         !
            if ( INDEX(CHARLINE,'arithmetic') > 0 ) then;                                !
                                                                                         !
                CH_PAIR_MODIFY_LIBRARY(IFILE) = TRIM(CH_PAIR_MODIFY_LIBRARY(IFILE))// &  !
                                                '   arithmetic';                         !
                                                                                         !
            end if                                                                       !
                                                                                         !
            IOSEF1 = 70 - 16 - 1 - LEN_TRIM(CH_PAIR_MODIFY_LIBRARY(IFILE));              !
                                                                                         !
            if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                !
                                                                                         !
            write(icanal,'(a70)') '| Pair modify : '//                  &                !
                                  TRIM(CH_PAIR_MODIFY_LIBRARY(IFILE))// &                !
                                  REPEAT(' ',IOSEF1)//'|';                               !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        else if ( INDEX(CHARLINE,'kspace_style') > 0 ) then;                             !
                                                                                         !
            read(CHARLINE,*) CHOSEF1, CH_KSPACE_STYLE_LIBRARY(IFILE);                    !
                                                                                         !
            IOSEF1 = 70 - 17 - 1 - LEN_TRIM(CH_KSPACE_STYLE_LIBRARY(IFILE));             !
                                                                                         !
            if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                !
                                                                                         !
            write(icanal,'(a70)') '| Kspace style : '//                  &               !
                                  TRIM(CH_KSPACE_STYLE_LIBRARY(IFILE))// &               !
                                  REPEAT(' ',IOSEF1)//'|';                               !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        end if                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
    close(5);                                                                            !
                                                                                         !
!   ### Check the presence of pair coefficients in the info file ###################################
                                                                                         !
    NPAIR_COEFF_CROSS_LIBRARY(IFILE) = 0;                                                !
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
        if ( INDEX(CHARLINE,'pair_coeff') > 0 ) then;                                    !
                                                                                         !
            NPAIR_COEFF_CROSS_LIBRARY(IFILE) = NPAIR_COEFF_CROSS_LIBRARY(IFILE) + 1;     !
                                                                                         !
            CYCLE;                                                                       !
                                                                                         !
        end if                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
    close(5);                                                                            !
                                                                                         !
    write(icanal,'(a11,i4,a55)') '| There is ',                         &                !
                                 NPAIR_COEFF_CROSS_LIBRARY(IFILE),      &                !
                                 ' cross interaction pair potential '// &                !
                                 'parameters to read'//                 &                !
                                 REPEAT(' ',2)//'|';                                     !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    if ( NPAIR_COEFF_CROSS_LIBRARY(IFILE) > 1000 ) then;                                 !
                                                                                         !
        write(icanal,'(a70)') '| The number of pairs found is greater than the '// &     !
                              'maximum number of pairs allowed in the current '//  &     !
                              'program'//REPEAT(' ',1)//'|';                             !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                  !
                                                                                         !
        write(icanal,*);                                                                 !
                                                                                         !
        write(icanal,'(a14)') 'End of program';                                          !
                                                                                         !
        close(icanal);                                                                   !
                                                                                         !
        stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Get atom IDs and pair coefficients #########################################################
                                                                                         !
    NPAIR_COEFF_CROSS_LIBRARY(IFILE) = 0;                                                !
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
        if ( INDEX(CHARLINE,'pair_coeff') > 0 ) then;                                    !
                                                                                         !
            NPAIR_COEFF_CROSS_LIBRARY(IFILE) = NPAIR_COEFF_CROSS_LIBRARY(IFILE) + 1;     !
                                                                                         !
            IOSEF1 = NPAIR_COEFF_CROSS_LIBRARY(IFILE);                                   !
                                                                                         !
            read(CHARLINE,*) CHOSEF1,                                     &              !
                             PAIR_ATOMID_CROSS_LIBRARY(1:2,IOSEF1,IFILE), &              !
                             PAIR_COEFF_CROSS_LIBRARY(1:2,IOSEF1,IFILE);                 !
                                                                                         !
        end if                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
    close(5);                                                                            !
                                                                                         !
    write(icanal,'(a70)') '| Pair coefficients for cross interactions were read'// &     !
                          REPEAT(' ',17)//'|';                                           !
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine READ_FILE_INFO_TEMPLATE
