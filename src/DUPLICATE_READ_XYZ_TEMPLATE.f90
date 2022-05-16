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

subroutine DUPLICATE_READ_XYZ_TEMPLATE(icanal,IFILE)

!   ************************************************************************************************
!   **                                READ XYZ TEMPLATE FILES                                     **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                           **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

    use module_library;

    use module_duplicate;

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ***********************************************************************************************

    integer (kind=4), intent(in) :: icanal, IFILE;

!   ************************************************************************************************

    integer (kind=4) :: i, j;

!   ************************************************************************************************

    integer (kind=4) :: IFOUND;

    integer (kind=4) :: EOF;

    character (len=250) :: CHEXT;

    character (len=250) :: CHARLINE;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the routine title ####################################################################
                                                                                         !
    CHTITLE = 'Read XYZ template (duplicate)';                                           !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write coordinates of templated molecules or monomers in arrays #############################
                                                                                         !
    write(icanal,'(a10,i8,a52)') '| IFILE : ', IFILE, REPEAT(' ',51)//'|'                !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Build the name of the file to read #########################################################
                                                                                         !
    CHEXT = TRIM(DUPLICATE_CHNAME_FILE)//'.xyz';                                         !
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
!   ### Initialization of library arrays ###########################################################
                                                                                         !
!   NATOM_DUPLICATE = 0;                                                                 !
                                                                                         !
!   NBOND_DUPLICATE = 0;                                                                 !
                                                                                         !
!   NANGLE_DUPLICATE = 0;                                                                !
                                                                                         !
!   NDIHEDRAL_DUPLICATE = 0;                                                             !
                                                                                         !
!   NIMPROPER_DUPLICATE = 0;                                                             !
                                                                                         !
!   NTYPE_ATOM_DUPLICATE = 0;                                                            !
                                                                                         !
!   NTYPE_BOND_DUPLICATE = 0;                                                            !
                                                                                         !
!   NTYPE_ANGLE_DUPLICATE = 0;                                                           !
                                                                                         !
!   NTYPE_DIHEDRAL_DUPLICATE = 0;                                                        !
                                                                                         !
!   NTYPE_IMPROPER_DUPLICATE = 0;                                                        !
                                                                                         !
    if ( NATOM_DUPLICATE > 0 ) then;                                                     !
                                                                                         ! 
        CONFIG_ATOMID_DUPLICATE(1:NATOM_DUPLICATE) = 0;                                  !
                                                                                         !
        CONFIG_ATOM_TYPE_DUPLICATE(1:NATOM_DUPLICATE) = 0;                               !
                                                                                         !
!       CONFIG_QI_LIBRARY(1:IMAX_NATOM,IFILE)         = 0.0d0;                           !
                                                                                         !
!       CONFIG_VI_LIBRARY(1:3,1:IMAX_NATOM,IFILE)     = 0.0d0;                           !
                                                                                         !
        CONFIG_RI_DUPLICATE(1:3,1:NATOM_DUPLICATE) = 0.0d0;                              !
                                                                                         !
        CONFIG_NAT_DUPLICATE(1:NATOM_DUPLICATE) = 'XXX';                                 !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   if ( IMAX_NBOND > 0 ) then;                                                          !
                                                                                         ! 
!       BOND_TYPE_LIBRARY(1:IMAX_NBOND,IFILE)         = 0;                               !
                                                                                         !
!       BOND_ATOMID_LIBRARY(1:2,1:IMAX_NBOND,IFILE)   = 0;                               !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   if ( IMAX_NANGLE > 0 ) then;                                                         !
                                                                                         !
!       ANGLE_TYPE_LIBRARY(1:IMAX_NANGLE,IFILE)       = 0;                               !
                                                                                         !
!       ANGLE_ATOMID_LIBRARY(1:3,1:IMAX_NANGLE,IFILE) = 0;                               !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   if ( IMAX_NDIHEDRAL > 0 ) then;                                                      !
                                                                                         !
!       DIHEDRAL_TYPE_LIBRARY(1:IMAX_NDIHEDRAL,IFILE)       = 0;                         !
                                                                                         !
!       DIHEDRAL_ATOMID_LIBRARY(1:4,1:IMAX_NDIHEDRAL,IFILE) = 0;                         !
                                                                                         ! 
!   end if                                                                               !
                                                                                         !
!   if ( IMAX_NIMPROPER > 0 ) then;                                                      !
                                                                                         !
!       IMPROPER_TYPE_LIBRARY(1:IMAX_NIMPROPER,IFILE)       = 0;                         !
                                                                                         !
!       IMPROPER_ATOMID_LIBRARY(1:4,1:IMAX_NIMPROPER,IFILE) = 0;                         !
                                                                                         !
!   end if                                                                               !
                                                                                         !
    write(icanal,'(a70)') '| Library arrays were initialized'//REPEAT(' ',36)//'|';      !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Reading the library file ###################################################################
                                                                                         !
    EOF = 0;                                                                             !
                                                                                         !
    open(2,file=TRIM(CHEXT));                                                            !
                                                                                         !
    read(2,*);                                                                           !
                                                                                         !
    read(2,*);                                                                           !
                                                                                         !
!   ### Read atomic coordinates ####################################################################
                                                                                         !
    do i = 1, NATOM_DUPLICATE;                                                           !
                                                                                         !
        read(2,*,iostat=EOF) CONFIG_NAT_DUPLICATE(i),    &                               !
                             CONFIG_RI_DUPLICATE(1:3,i);                                 ! 
                                                                                         !
        if ( EOF /= 0 ) then;                                                            !
                                                                                         !
            write(icanal,*) i;                                                           !
                                                                                         !
            write(icanal,*);                                                             !
                                                                                         !
            write(icanal,*) 'There was a problem in reading coordinates of the xyz '// & !
                            'file - stop';                                               !
                                                                                         !
            stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end if                                                                           !
                                                                                         !
        CONFIG_ATOMID_DUPLICATE(i) = i;                                                  !
                                                                                         !
    end do                                                                               !
                                                                                         !
    close(2);                                                                            !
                                                                                         !
    write(icanal,'(a70)') '| Coords were read'//REPEAT(' ',51)//'|';                     !                
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set the list of atom labels ################################################################
                                                                                         !
    do i = 1, NATOM_DUPLICATE;                                                           !
                                                                                         !
        IFOUND = 0;                                                                      !
                                                                                         !
        CHOSEF1 = TRIM(CONFIG_NAT_DUPLICATE(i));                                         !
                                                                                         !
        if ( NTYPE_ATOM_DUPLICATE > 0 ) then;                                            !
                                                                                         !
            IFOUND = 0;                                                                  !
                                                                                         !
            do j = 1, NTYPE_ATOM_DUPLICATE;                                              !
                                                                                         !
                if ( TRIM(ATOM_LABEL_DUPLICATE(j)) == TRIM(CHOSEF1) ) then;              !
                                                                                         !
                    CONFIG_ATOM_TYPE_DUPLICATE(i) = j;                                   !
                                                                                         !
                    IFOUND = 1;                                                          !
                                                                                         !
                    EXIT;                                                                !
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
        end if                                                                           !
                                                                                         !
        if ( IFOUND == 1 ) CYCLE;                                                        !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '| The list of atom labels was set'//REPEAT(' ',36)//'|';      !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set atom masses ############################################################################
                                                                                         !
!   ATOM_MASSES_LIBRARY(1:NTYPE_ATOM_LIBRARY(IFILE),IFILE) = 0.0d0;                      !
                                                                                         !
!   do i = 1, NTYPE_ATOM_LIBRARY(IFILE);                                                 !
                                                                                         !
!       do j = 1, 200;                                                                   !
                                                                                         !
!           if ( TRIM(ATOM_LABEL_LIBRARY(i,IFILE)) == &
!                TRIM(ELEMENT_SYMBOL_NAME(j)) ) then;!
                                                                                         !
!               ATOM_MASSES_LIBRARY(i,IFILE) = MOLAR_MASS_ELEMENT(j);                    !
                                                                                         !
!           end if                                                                       !
                                                                                         !
!       end do                                                                           !
                                                                                         ! 
!   end do                                                                               !
                                                                                         !
!   write(icanal,'(a70)') '| Masses were set '//REPEAT(' ',51)//'|';                     !
                                                                                         !
!   write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set units and styles #######################################################################
                                                                                         !
    CH_UNITS_STYLE_DUPLICATE = 'real';                                                   !
                                                                                         !
    CH_ATOM_STYLE_DUPLICATE = 'full';                                                    !
                                                                                         !
!   POTENTIAL_CLASS2_CHTYPE_LIBRARY(IFILE) = 'XXX';                                      !
                                                                                         !
!   IPOTENTIAL_CLASS2_LIBRARY(IFILE) = 0;                                                !
                                                                                         !
!   ### Write properties of the library file #######################################################
                                                                                         !
    if ( NATOM_DUPLICATE > 0 ) then;                                                     !
                                                                                         !
        write(icanal,'(a29,i8,a33)') '| NATOM_DUPLICATE          : ', &                  !
                                     NATOM_DUPLICATE,                 &                  !
                                     REPEAT(' ',32)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   if ( NBOND_LIBRARY(IFILE) > 0 ) then;                                                !
                                                                                         !
!       write(icanal,'(a25,i8,a37)') '| NBOND_TPLTE          : ', &                      !
!                                    NBOND_LIBRARY(IFILE),        &                      !
!                                    REPEAT(' ',36)//'|';                                !
                                                                                         !
!       write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   if ( NANGLE_LIBRARY(IFILE) > 0 ) then;                                               !
                                                                                         !
!       write(icanal,'(a25,i8,a37)') '| NANGLE_TPLTE         : ', &                      !
!                                    NANGLE_LIBRARY(IFILE),       &                      !
!                                    REPEAT(' ',36)//'|';                                !
                                                                                         !
!       write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  ! 
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   if ( NDIHEDRAL_LIBRARY(IFILE) > 0 ) then;                                            !
                                                                                         !
!       write(icanal,'(a25,i8,a37)') '| NDIHEDRAL_TPLTE      : ', &                      !
!                                    NDIHEDRAL_LIBRARY(IFILE),    &                      !
!                                    REPEAT(' ',36)//'|';                                !
!   end if                                                                               !
                                                                                         !
!   if ( NIMPROPER_LIBRARY(IFILE) > 0 ) then;                                            !
                                                                                         !
!       write(icanal,'(a25,i8,a37)') '| NIMPROPER_TPLTE      : ', &                      !
!                                    NIMPROPER_LIBRARY(IFILE),    &                      !
!                                    REPEAT(' ',36)//'|';                                !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    if ( NTYPE_ATOM_DUPLICATE > 0 ) then;                                                !
                                                                                         !
        write(icanal,'(a31,i8,a31)') '| NTYPE_ATOM_DUPLICATE       : ', &                !
                                     NTYPE_ATOM_DUPLICATE,              &                !
                                     REPEAT(' ',30)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   if ( NTYPE_BOND_LIBRARY(IFILE) > 0 ) then;                                           !
                                                                                         !
!       write(icanal,'(a25,i8,a37)') '| NTYPE_BOND_TPLTE     : ', &                      !
!                                    NTYPE_BOND_LIBRARY(IFILE),   &                      !
!                                    REPEAT(' ',36)//'|';                                !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   if ( NTYPE_ANGLE_LIBRARY(IFILE) > 0 ) then;                                          !
                                                                                         !
!       write(icanal,'(a25,i8,a37)') '| NTYPE_ANGLE_TPLTE    : ', &                      !
!                                    NTYPE_ANGLE_LIBRARY(IFILE),  &                      !
!                                    REPEAT(' ',36)//'|';                                !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   if ( NTYPE_DIHEDRAL_LIBRARY(IFILE) > 0 ) then;                                       !
                                                                                         !
!       write(icanal,'(a25,i8,a37)') '| NTYPE_DIHEDRAL_TPLTE : ',   &                    !
!                                    NTYPE_DIHEDRAL_LIBRARY(IFILE), &                    !
!                                    REPEAT(' ',36)//'|';                                !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   if ( NTYPE_IMPROPER_LIBRARY(IFILE) > 0 ) then;                                       !
                                                                                         !
!       write(icanal,'(a25,i8,a37)') '| NTYPE_IMPROPER_TPLTE : ',   &                    !
!                                    NTYPE_IMPROPER_LIBRARY(IFILE), &                    !
!                                    REPEAT(' ',36)//'|';                                !
                                                                                         !
!   end if                                                                               !
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine DUPLICATE_READ_XYZ_TEMPLATE
