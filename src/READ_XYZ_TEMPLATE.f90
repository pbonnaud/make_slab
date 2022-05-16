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

subroutine READ_XYZ_TEMPLATE(icanal,IFILE,IINSREM)

!   ************************************************************************************************
!   **                                READ XYZ TEMPLATE FILES                                     **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                           **
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
    CHTITLE = 'Read XYZ template';                                                       !
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
    if ( IINSREM == 1 ) then;                                                            !
                                                                                         !
        CHEXT = TRIM(CHNAME_FILE_LIBRARY(1,IFILE))//'.xyz';                              !
                                                                                         !
    else if ( IINSREM == 0 ) then;                                                       !
                                                                                         !
        CHEXT = TRIM(CHNAME_FILE_LIBRARY_REMOVE(IFILE))//'.xyz';                         !
                                                                                         !
    else if ( IINSREM == 2 ) then;                                                       !
                                                                                         ! 
        CHEXT = TRIM(Chname_file_library_monomer(IFILE))//'.xyz';                        !
                                                                                         !
    else                                                                                 !
                                                                                         !
        CHEXT = TRIM(CHNAME_FILE_LIBRARY(1,IFILE));                                      !
                                                                                         !
    end if                                                                               !
                                                                                         !
    IOSEF1 = LEN_TRIM(CHEXT);                                                            !
                                                                                         !
    IOSEF2 = 67 - IOSEF1;                                                                !
                                                                                         !
    if ( IOSEF2 < 0 ) IOSEF2 = 1;                                                        !
                                                                                         ! 
    write(icanal,'(a70)') '| '//TRIM(CHEXT)//REPEAT(' ',IOSEF2)//'|';                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    call LIBRARY_INITIALIZATION_ARRAYS(icanal,IFILE);                                    !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialization of library arrays ###########################################################
                                                                                         !
!   NATOM_LIBRARY(IFILE)     = 0;                                                        !
                                                                                         !
!   NBOND_LIBRARY(IFILE)     = 0;                                                        !
                                                                                         !
!   NANGLE_LIBRARY(IFILE)    = 0;                                                        !
                                                                                         !
!   NDIHEDRAL_LIBRARY(IFILE) = 0;                                                        !
                                                                                         !
!   NIMPROPER_LIBRARY(IFILE) = 0;                                                        !
                                                                                         !
!   NTYPE_ATOM_LIBRARY(IFILE)     = 0;                                                   !
                                                                                         !
!   NTYPE_BOND_LIBRARY(IFILE)     = 0;                                                   !
                                                                                         !
!   NTYPE_ANGLE_LIBRARY(IFILE)    = 0;                                                   !
                                                                                         !
!   NTYPE_DIHEDRAL_LIBRARY(IFILE) = 0;                                                   !
                                                                                         !
!   NTYPE_IMPROPER_LIBRARY(IFILE) = 0;                                                   !
                                                                                         !
!   if ( IMAX_NATOM > 0 ) then;                                                          !
                                                                                         ! 
!       CONFIG_ATOMID_LIBRARY(1:IMAX_NATOM,IFILE)     = 0;                               !
                                                                                         !
!       CONFIG_ATOM_TYPE_LIBRARY(1:IMAX_NATOM,IFILE)  = 0;                               !
                                                                                         !
!       CONFIG_QI_LIBRARY(1:IMAX_NATOM,IFILE)         = 0.0d0;                           !
                                                                                         !
!       CONFIG_VI_LIBRARY(1:3,1:IMAX_NATOM,IFILE)     = 0.0d0;                           !
                                                                                         !
!       CONFIG_RI_LIBRARY(1:3,1:IMAX_NATOM,IFILE)     = 0.0d0;                           !
                                                                                         !
!       CONFIG_NAT_LIBRARY(1:IMAX_NATOM,IFILE)        = 'XXX';                           !
                                                                                         !
!   end if                                                                               !
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
!   write(icanal,'(a70)') '| Library arrays were initialized'//REPEAT(' ',36)//'|';      !
                                                                                         !
!   write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Reading the library file ###################################################################
                                                                                         !
    EOF = 0;                                                                             !
                                                                                         !
    open(2,file=TRIM(CHEXT));                                                            !
                                                                                         !
!   ### Read the number of atoms in the configuration to read ######################################
                                                                                         !
    read(2,*,iostat=EOF) NATOM_LIBRARY(IFILE);                                           !
                                                                                         !
    if ( EOF /= 0 ) then;                                                                !
                                                                                         !
        write(icanal,*) 'There was a problem in reading the number of atom of the current file - stop';
                                                                                         !
        stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
    end if                                                                               !
                                                                                         !
    write(icanal,'(a25,i8,a37)') '| NATOM_LIBRARY(IFILE) : ', &                          !
                                 NATOM_LIBRARY(IFILE),        &                          !
                                 REPEAT(' ',36)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Read cell properties of the xyz configuration ##############################################
                                                                                         !
    if ( IINSREM == 2 ) then;                                                            !
                                                                                         !
        read(2,*,iostat=EOF);                                                            !
                                                                                         !
    else                                                                                 !
                                                                                         !
        read(2,*,iostat=EOF) CELL_AXIS_LIBRARY(1:3,IFILE);                               !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( EOF /= 0 ) then;                                                                !
                                                                                         !
        write(icanal,*) 'There was a problem in reading the 2nd line of the xyz file - stop';
                                                                                         !
        stop; !//////////////////////////////////////////////////////////////////////////!
                                                                                         !
    end if                                                                               !
                                                                                         !
    write(icanal,'(a70)') '| 2nd line was read'//REPEAT(' ',50)//'|';                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         ! 
!   ### Read atomic coordinates ####################################################################
                                                                                         !
    do i = 1, NATOM_LIBRARY(IFILE);                                                      !
                                                                                         !
!       write(icanal,*) '--> ', i;

!       read(2,*,iostat=EOF);                                                            !
                                                                                         !  
!       read(2,*,iostat=EOF) CONFIG_NAT_LIBRARY(i,IFILE);                                !

        read(2,*,iostat=EOF) CONFIG_NAT_LIBRARY(i,IFILE),   &                            !
                             CONFIG_RI_LIBRARY(1:3,i,IFILE);                             ! 
                                                                                         !
        if ( EOF /= 0 ) then;                                                            !
                                                                                         !
            write(icanal,*) i;

            write(icanal,*);
 
            write(icanal,*) 'There was a problem in reading coordinates of the xyz file - stop';
                                                                                         !
            stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        end if                                                                           !
                                                                                         !
        CONFIG_ATOMID_LIBRARY(i,IFILE) = i;                                              !
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
!   ### Write cell axis of the library file (if found) #############################################
                                                                                         !
    if ( EOF == 0 ) then;                                                                !
                                                                                         !
        if ( IINSREM /= 2 ) then;                                                        !
                                                                                         !
            write(icanal,'(a26,3f12.6,a8)') '| CELL_AXIS_LIBRARY [A] : ', &              !
                                            CELL_AXIS_LIBRARY(1:3,IFILE), &              !
                                            REPEAT(' ',7)//'|';                          !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        end if                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set the list of atom labels ################################################################
                                                                                         !
    NTYPE_ATOM_LIBRARY(IFILE) = 0;                                                       !
                                                                                         !
    do i = 1, NATOM_LIBRARY(IFILE);                                                      !
                                                                                         !
        IFOUND = 0;                                                                      !
                                                                                         !
        CHOSEF1 = TRIM(CONFIG_NAT_LIBRARY(i,IFILE));                                     !
                                                                                         !
        if ( NTYPE_ATOM_LIBRARY(IFILE) > 0 ) then;                                       !
                                                                                         !
            IFOUND = 0;                                                                  !
                                                                                         !
            do j = 1, NTYPE_ATOM_LIBRARY(IFILE);                                         !
                                                                                         !
                if ( TRIM(ATOM_LABEL_LIBRARY(j,IFILE)) == TRIM(CHOSEF1) ) then;          !
                                                                                         !
                    CONFIG_ATOM_TYPE_LIBRARY(i,IFILE) = j;                               !
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
        NTYPE_ATOM_LIBRARY(IFILE) = NTYPE_ATOM_LIBRARY(IFILE) + 1;                       !
                                                                                         !
        ATOM_LABEL_LIBRARY(NTYPE_ATOM_LIBRARY(IFILE),IFILE) = TRIM(CHOSEF1);             !
                                                                                         !
        CONFIG_ATOM_TYPE_LIBRARY(i,IFILE) = NTYPE_ATOM_LIBRARY(IFILE);                   !
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
    ATOM_MASSES_LIBRARY(1:NTYPE_ATOM_LIBRARY(IFILE),IFILE) = 0.0d0;                      !
                                                                                         !
    do i = 1, NTYPE_ATOM_LIBRARY(IFILE);                                                 !
                                                                                         !
        do j = 1, 200;                                                                   !
                                                                                         !
            if ( TRIM(ATOM_LABEL_LIBRARY(i,IFILE)) == &                                  !
                 TRIM(ELEMENT_SYMBOL_NAME(j)) ) then;                                    !
                                                                                         !
                ATOM_MASSES_LIBRARY(i,IFILE) = MOLAR_MASS_ELEMENT(j);                    !
                                                                                         !
            end if                                                                       !
                                                                                         !
        end do                                                                           !
                                                                                         ! 
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '| Masses were set '//REPEAT(' ',51)//'|';                     !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set units and styles #######################################################################
                                                                                         !
    CH_UNITS_STYLE_LIBRARY(i) = 'real';                                                  !
                                                                                         !
    CH_ATOM_STYLE_LIBRARY(IFILE) = 'full';                                               !
                                                                                         !
    POTENTIAL_CLASS2_CHTYPE_LIBRARY(IFILE) = 'XXX';                                      !
                                                                                         !
    IPOTENTIAL_CLASS2_LIBRARY(IFILE) = 0;                                                !
                                                                                         !
    PAIR_COEFF_CROSS_LIBRARY(1:2,1:1000,IFILE) = 0.0d0;                                  !
                                                                                         !
    PAIR_ATOMID_CROSS_LIBRARY(1:2,1:1000,1:IFILE) = 0;                                   !
                                                                                         !
!   ### Write properties of the library file #######################################################
                                                                                         !
    if ( NATOM_LIBRARY(IFILE) > 0 ) then;                                                !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| NATOM_TPLTE          : ', &                      !
                                     NATOM_LIBRARY(IFILE),        &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NBOND_LIBRARY(IFILE) > 0 ) then;                                                !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| NBOND_TPLTE          : ', &                      !
                                     NBOND_LIBRARY(IFILE),        &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NANGLE_LIBRARY(IFILE) > 0 ) then;                                               !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| NANGLE_TPLTE         : ', &                      !
                                     NANGLE_LIBRARY(IFILE),       &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  ! 
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NDIHEDRAL_LIBRARY(IFILE) > 0 ) then;                                            !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| NDIHEDRAL_TPLTE      : ', &                      !
                                     NDIHEDRAL_LIBRARY(IFILE),    &                      !
                                     REPEAT(' ',36)//'|';                                !
    end if                                                                               !
                                                                                         !
    if ( NIMPROPER_LIBRARY(IFILE) > 0 ) then;                                            !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| NIMPROPER_TPLTE      : ', &                      !
                                     NIMPROPER_LIBRARY(IFILE),    &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
    end if                                                                               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    if ( NTYPE_ATOM_LIBRARY(IFILE) > 0 ) then;                                           !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| NTYPE_TEMPLATE       : ', &                      !
                                     NTYPE_ATOM_LIBRARY(IFILE),   &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_BOND_LIBRARY(IFILE) > 0 ) then;                                           !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| NTYPE_BOND_TPLTE     : ', &                      !
                                     NTYPE_BOND_LIBRARY(IFILE),   &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_ANGLE_LIBRARY(IFILE) > 0 ) then;                                          !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| NTYPE_ANGLE_TPLTE    : ', &                      !
                                     NTYPE_ANGLE_LIBRARY(IFILE),  &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_DIHEDRAL_LIBRARY(IFILE) > 0 ) then;                                       !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| NTYPE_DIHEDRAL_TPLTE : ',   &                    !
                                     NTYPE_DIHEDRAL_LIBRARY(IFILE), &                    !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_IMPROPER_LIBRARY(IFILE) > 0 ) then;                                       !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| NTYPE_IMPROPER_TPLTE : ',   &                    !
                                     NTYPE_IMPROPER_LIBRARY(IFILE), &                    !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine READ_XYZ_TEMPLATE
