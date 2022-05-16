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


subroutine READ_LAMMPS_TEMPLATE(icanal,IFILE,IINSREM)

!   ************************************************************************************************
!   **                                READ LAMMPS TEMPLATE FILES                                  **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                : CANAL ON WHICH THE OUTPUT FILE IS WRITEN                           **
!   **                                                                                            **
!   ************************************************************************************************

    use module_size_arrays;

    use module_library;

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ***********************************************************************************************

    integer (kind=4), intent(in) :: icanal, IFILE, IINSREM;

!   ************************************************************************************************

    character (len=250) :: CHEXT;

    integer (kind=4) :: i, j;

    character (len=250) :: CHARLINE;

    integer (kind=4) :: EOF;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the routine title ####################################################################
                                                                                         !
    CHTITLE = 'Read LAMMPS template';                                                    !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the position of the templated molecule in library arrays #############################
                                                                                         !
    write(icanal,'(a10,i8,a52)') '| IFILE : ', IFILE, REPEAT(' ',51)//'|'                !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### Build the name of the templated file #######################################################
                                                                                         !
    if ( IINSREM == 1 ) then;                                                            !
                                                                                         !                                                                                    
        CHEXT = TRIM(CHNAME_FILE_LIBRARY(1,IFILE));                                      !
                                                                                         !
    else                                                                                 !
                                                                                         !
        CHEXT = TRIM(CHNAME_FILE_LIBRARY_REMOVE(IFILE))//'.template';                    !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Write the name of the considered tempated file #############################################
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
!   ### Initialization of arrays containing properties of the library molecule #####################
                                                                                         !
    NATOM_LIBRARY(IFILE)     = 0;                                                        !
                                                                                         !
    NBOND_LIBRARY(IFILE)     = 0;                                                        !
                                                                                         !
    NANGLE_LIBRARY(IFILE)    = 0;                                                        !
                                                                                         !
    NDIHEDRAL_LIBRARY(IFILE) = 0;                                                        !
                                                                                         !
    NIMPROPER_LIBRARY(IFILE) = 0;                                                        !
                                                                                         !
    NTYPE_ATOM_LIBRARY(IFILE)     = 0;                                                   !
                                                                                         !
    NTYPE_BOND_LIBRARY(IFILE)     = 0;                                                   !
                                                                                         !
    NTYPE_ANGLE_LIBRARY(IFILE)    = 0;                                                   !
                                                                                         !
    NTYPE_DIHEDRAL_LIBRARY(IFILE) = 0;                                                   !
                                                                                         !
    NTYPE_IMPROPER_LIBRARY(IFILE) = 0;                                                   !
                                                                                         !
    if ( IMAX_NATOM > 0 ) then;                                                          !
                                                                                         !
        CONFIG_ATOMID_LIBRARY(1:IMAX_NATOM,IFILE)     = 0;                               !
                                                                                         !
        CONFIG_ATOM_TYPE_LIBRARY(1:IMAX_NATOM,IFILE)  = 0;                               !
                                                                                         !
        CONFIG_QI_LIBRARY(1:IMAX_NATOM,IFILE)         = 0.0d0;                           !
                                                                                         !
        CONFIG_VI_LIBRARY(1:3,1:IMAX_NATOM,IFILE)     = 0.0d0;                           !
                                                                                         !
        CONFIG_RI_LIBRARY(1:3,1:IMAX_NATOM,IFILE)     = 0.0d0;                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( IMAX_NBOND > 0 ) then;                                                          !
                                                                                         !
        BOND_TYPE_LIBRARY(1:IMAX_NBOND,IFILE)         = 0;                               !
                                                                                         !
        BOND_ATOMID_LIBRARY(1:2,1:IMAX_NBOND,IFILE)   = 0;                               !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( IMAX_NANGLE > 0 ) then;                                                         !
                                                                                         !
        ANGLE_TYPE_LIBRARY(1:IMAX_NANGLE,IFILE)       = 0;                               !
                                                                                         !
        ANGLE_ATOMID_LIBRARY(1:3,1:IMAX_NANGLE,IFILE) = 0;                               !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( IMAX_NDIHEDRAL > 0 ) then;                                                      !
                                                                                         !
        DIHEDRAL_TYPE_LIBRARY(1:IMAX_NDIHEDRAL,IFILE)       = 0;                         !
                                                                                         !
        DIHEDRAL_ATOMID_LIBRARY(1:4,1:IMAX_NDIHEDRAL,IFILE) = 0;                         !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( IMAX_NIMPROPER > 0 ) then;                                                      !
                                                                                         !
        IMPROPER_TYPE_LIBRARY(1:IMAX_NIMPROPER,IFILE)       = 0;                         !
                                                                                         !
        IMPROPER_ATOMID_LIBRARY(1:4,1:IMAX_NIMPROPER,IFILE) = 0;                         !
                                                                                         !
    end if                                                                               !
                                                                                         !
    write(icanal,'(a70)') '| Arrays were initialized'//REPEAT(' ',44)//'|';              !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Set the atom style if no info files were found in the library ##############################
                                                                                         !
    if ( ILIBRARY_FILE_INFO(1) == 0 ) then;                                              !
                                                                                         !
        CH_ATOM_STYLE_LIBRARY(IFILE) = 'full';                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Read the library file template #############################################################
                                                                                         !
    open(2,file=TRIM(CHEXT));                                                            !
                                                                                         !
    EOF = 0;                                                                             !
                                                                                         !
    do                                                                                   !
        read(2,'(a)',iostat=EOF) CHARLINE;                                               !
                                                                                         !
        if ( EOF /= 0 ) EXIT;                                                            !
                                                                                         !
        if ( INDEX(CHARLINE,'atoms') > 0 ) then;                                         !
                                                                                         !
            read(CHARLINE,*) NATOM_LIBRARY(IFILE);                                       !
                                                                                         !
            write(icanal,'(a70)') '| Number of atoms was read'//REPEAT(' ',43)//'|';     !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        else if ( INDEX(CHARLINE,'bonds') > 0 ) then;                                    !
                                                                                         !
            read(CHARLINE,*) NBOND_LIBRARY(IFILE);                                       !
                                                                                         !
        else if ( INDEX(CHARLINE,'angles') > 0 ) then;                                   !
                                                                                         !
            read(CHARLINE,*) NANGLE_LIBRARY(IFILE);                                      !
                                                                                         !
        else if ( INDEX(CHARLINE,'dihedrals') > 0 ) then;                                !
                                                                                         !
            read(CHARLINE,*) NDIHEDRAL_LIBRARY(IFILE);                                   !
                                                                                         !
            write(icanal,'(a70)') '| Number of dihedrals was read'//REPEAT(' ',39)//'|'; !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        else if ( INDEX(CHARLINE,'impropers') > 0 ) then;                                !
                                                                                         !
            read(CHARLINE,*) NIMPROPER_LIBRARY(IFILE);                                   !
                                                                                         !
            write(icanal,'(a70)') '| Number of impropers was read'//REPEAT(' ',39)//'|'; !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
!           stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
        else if ( INDEX(CHARLINE,'Charges') > 0 ) then;                                  !
                                                                                         !
            read(2,*);                                                                   !
                                                                                         !
            if ( NATOM_LIBRARY(IFILE) > 0 ) then;                                        !
                                                                                         !
                do i = 1, NATOM_LIBRARY(IFILE);                                          !
                                                                                         !
                    read(2,*) IOSEF1, CONFIG_QI_LIBRARY(IOSEF1,IFILE);                   !
                                                                                         !
                end do                                                                   !
                                                                                         !
                write(icanal,'(a70)') '| Charges were read'//REPEAT(' ',50)//'|';        !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
        else if ( INDEX(CHARLINE,'Coords') > 0 ) then;                                   !
                                                                                         !
            read(2,*);                                                                   !
                                                                                         !
            if ( NATOM_LIBRARY(IFILE) > 0 ) then;                                        !
                                                                                         !
                do i = 1, NATOM_LIBRARY(IFILE);                                          !
                                                                                         !
                    read(2,*) CONFIG_ATOMID_LIBRARY(i,IFILE), &                          !
                              CONFIG_RI_LIBRARY(1:3,          &                          !
                              CONFIG_ATOMID_LIBRARY(i,IFILE), &                          !
                              IFILE);                                                    !
                                                                                         !
                end do                                                                   !
                                                                                         !
                write(icanal,'(a70)') '| Coords were read'//REPEAT(' ',51)//'|';         !                
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
        else if ( INDEX(CHARLINE,'Types') > 0 ) then;                                    !
                                                                                         !
            read(2,*);                                                                   !
                                                                                         !
            if ( NATOM_LIBRARY(IFILE) > 0 ) then;                                        !
                                                                                         !
                do i = 1, NATOM_LIBRARY(IFILE);                                          !
                                                                                         !
                    read(2,*) IOSEF1, CONFIG_ATOM_TYPE_LIBRARY(IOSEF1,IFILE);            !
                                                                                         !
                    if ( CONFIG_ATOM_TYPE_LIBRARY(IOSEF1,IFILE) > &                      !
                         NTYPE_ATOM_LIBRARY(IFILE) )              &                      !
                         NTYPE_ATOM_LIBRARY(IFILE) =              &                      !
                         CONFIG_ATOM_TYPE_LIBRARY(IOSEF1,IFILE);                         !
                                                                                         !
                end do                                                                   !
                                                                                         !
                write(icanal,'(a70)') '| Types were read'//REPEAT(' ',52)//'|';          !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
        else if ( INDEX(CHARLINE,'Bonds') > 0 ) then;                                    !
                                                                                         !
            read(2,*);                                                                   !
                                                                                         !
            if ( NBOND_LIBRARY(IFILE) > 0 ) then;                                        !
                                                                                         !
                do i = 1, NBOND_LIBRARY(IFILE);                                          ! 
                                                                                         !
                    read(2,*) IOSEF1,                               &                    ! 
                              BOND_TYPE_LIBRARY(IOSEF1,IFILE),      &                    !
                              BOND_ATOMID_LIBRARY(1:2,IOSEF1,IFILE);                     !
                                                                                         !
                    if ( BOND_TYPE_LIBRARY(IOSEF1,IFILE) > &                             !
                         NTYPE_BOND_LIBRARY(IFILE) )       &                             !
                         NTYPE_BOND_LIBRARY(IFILE) =       &                             !
                         BOND_TYPE_LIBRARY(IOSEF1,IFILE);                                !
                                                                                         !
                end do                                                                   !
                                                                                         !
                write(icanal,'(a70)') '| Bonds were read'//REPEAT(' ',52)//'|';          !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
        else if ( INDEX(CHARLINE,'Angles') > 0 ) then;                                   !
                                                                                         !
            read(2,*);                                                                   !
                                                                                         !
            if ( NANGLE_LIBRARY(IFILE) > 0 ) then;                                       !
                                                                                         !
                do i = 1, NANGLE_LIBRARY(IFILE);                                         !
                                                                                         !
                    read(2,*) IOSEF1,                                &                   !
                              ANGLE_TYPE_LIBRARY(IOSEF1,IFILE),      &                   !
                              ANGLE_ATOMID_LIBRARY(1:3,IOSEF1,IFILE);                    !
                                                                                         !
                    if ( ANGLE_TYPE_LIBRARY(IOSEF1,IFILE) > &                            !
                         NTYPE_ANGLE_LIBRARY(IFILE) )       &                            !
                         NTYPE_ANGLE_LIBRARY(IFILE) =       &                            !
                         ANGLE_TYPE_LIBRARY(IOSEF1,IFILE);                               !
                                                                                         !
                end do                                                                   !
                                                                                         !
                write(icanal,'(a70)') '| Angles were read'//REPEAT(' ',51)//'|';         !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
         else if ( INDEX(CHARLINE,'Dihedrals') > 0 ) then;                               !
                                                                                         !
            read(2,*);                                                                   !
                                                                                         !
            if ( NDIHEDRAL_LIBRARY(IFILE) > 0 ) then;                                    !
                                                                                         !
                do i = 1, NDIHEDRAL_LIBRARY(IFILE);                                      !
                                                                                         !
                    read(2,*) IOSEF1,                                   &                !
                              DIHEDRAL_TYPE_LIBRARY(IOSEF1,IFILE),      &                !
                              DIHEDRAL_ATOMID_LIBRARY(1:4,IOSEF1,IFILE);                 !
                                                                                         !
                    if ( DIHEDRAL_TYPE_LIBRARY(IOSEF1,IFILE) > &                         !
                         NTYPE_DIHEDRAL_LIBRARY(IFILE) )       &                         !
                         NTYPE_DIHEDRAL_LIBRARY(IFILE) =       &                         !
                         DIHEDRAL_TYPE_LIBRARY(IOSEF1,IFILE);                            !
                                                                                         !
                end do                                                                   !
                                                                                         ! 
                write(icanal,'(a70)') '| Dihedrals were read'//REPEAT(' ',48)//'|';      !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
         else if ( INDEX(CHARLINE,'Impropers') > 0 ) then;                               !
                                                                                         !
            read(2,*);                                                                   !
                                                                                         !
            if ( NIMPROPER_LIBRARY(IFILE) > 0 ) then;                                    ! 
                                                                                         !
                do i = 1, NIMPROPER_LIBRARY(IFILE);                                      !
                                                                                         !
                    read(2,*) IOSEF1,                                   &                !
                              IMPROPER_TYPE_LIBRARY(IOSEF1,IFILE),      &                !
                              IMPROPER_ATOMID_LIBRARY(1:4,IOSEF1,IFILE);                 !
                                                                                         !
                    if ( IMPROPER_TYPE_LIBRARY(IOSEF1,IFILE) > &                         !
                         NTYPE_IMPROPER_LIBRARY(IFILE) )       &                         !
                         NTYPE_IMPROPER_LIBRARY(IFILE) =       &                         !
                         IMPROPER_TYPE_LIBRARY(IOSEF1,IFILE);                            !
                                                                                         !
                end do                                                                   !
                                                                                         !
                write(icanal,'(a70)') '| Impropers were read'//REPEAT(' ',48)//'|';      !
                                                                                         !
                write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                          !
                                                                                         !
!               stop; !//////////////////////////////////////////////////////////////////!
                                                                                         !
            end if                                                                       !
                                                                                         !
        end if                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
    close(2);                                                                            !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a25,i8,a37)') '| NATOM_TPLTE          : ', &                          !
                                 NATOM_LIBRARY(IFILE),        &                          !
                                 REPEAT(' ',36)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a25,i8,a37)') '| NBOND_TPLTE          : ', &                          !
                                 NBOND_LIBRARY(IFILE),        &                          !
                                 REPEAT(' ',36)//'|';                                    !
    write(icanal,'(a25,i8,a37)') '| NANGLE_TPLTE         : ', &                          !
                                 NANGLE_LIBRARY(IFILE),       &                          !
                                 REPEAT(' ',36)//'|';                                    !
                                                                                         !
    if ( NDIHEDRAL_LIBRARY(IFILE) > 0 ) then;                                            !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| NDIHEDRAL_TPLTE      : ', &                      !
                                     NDIHEDRAL_LIBRARY(IFILE),    &                      !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
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
    write(icanal,'(a25,i8,a37)') '| NTYPE_TEMPLATE       : ', &                          !
                                 NTYPE_ATOM_LIBRARY(IFILE),   &                          !
                                 REPEAT(' ',36)//'|';                                    !
    write(icanal,'(a25,i8,a37)') '| NTYPE_BOND_TPLTE     : ', &                          !
                                 NTYPE_BOND_LIBRARY(IFILE),   &                          !
                                 REPEAT(' ',36)//'|';                                    !
    write(icanal,'(a25,i8,a37)') '| NTYPE_ANGLE_TPLTE    : ', &                          !
                                 NTYPE_ANGLE_LIBRARY(IFILE),  &                          !
                                 REPEAT(' ',36)//'|';                                    !
                                                                                         !
    if (  NTYPE_DIHEDRAL_LIBRARY(IFILE) > 0 ) then;                                      !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| NTYPE_DIHEDRAL_TPLTE : ',   &                    !
                                     NTYPE_DIHEDRAL_LIBRARY(IFILE), &                    !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
    if ( NTYPE_IMPROPER_LIBRARY(IFILE) > 0 ) then;                                       !
                                                                                         !
        write(icanal,'(a25,i8,a37)') '| NTYPE_IMPROPER_TPLTE : ',   &                    !
                                     NTYPE_IMPROPER_LIBRARY(IFILE), &                    !
                                     REPEAT(' ',36)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
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
end subroutine READ_LAMMPS_TEMPLATE
