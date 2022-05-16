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

subroutine DEALLOCATE_LIBRARY_ARRAYS(icanal) 

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

    integer (kind=4), intent(in) :: icanal;

!   ************************************************************************************************

    integer (kind=4) :: i, j;

!   ************************************************************************************************

    character (len=250) :: CHARLINE;

    integer (kind=4) :: EOF;

    integer (kind=4) :: IOSEF1, IOSEF2, IOSEF3, IOSEF4;

    character (len=10) :: CHOSEF1, CHOSEF2, CHOSEF3;

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
    CHTITLE = 'DEALLOCATE LIBRARY ARRAYS';                                               !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);

    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));

!   deallocate(CONFIG_ATOMID_LIBRARY,             &
!              CONFIG_ATOM_TYPE_LIBRARY,          &
!              CONFIG_QI_LIBRARY,                 &
!              CONFIG_RI_LIBRARY,                 &
!              CONFIG_NAT_LIBRARY,                &
!              BOND_TYPE_LIBRARY,                 &
!              BOND_ATOMID_LIBRARY,               &
!              ANGLE_TYPE_LIBRARY,                &
!              ANGLE_ATOMID_LIBRARY,              &
!              DIHEDRAL_TYPE_LIBRARY,             &
!              DIHEDRAL_ATOMID_LIBRARY,           &
!              IMPROPER_TYPE_LIBRARY,             &
!              IMPROPER_ATOMID_LIBRARY,           &
!              ATOM_MASSES_LIBRARY,               &
!              ATOM_LABEL_LIBRARY,                &
!              POTENTIAL_CLASS2_LIBRARY,          &
!              BOND_COEFFS_LIBRARY,               &
!              ANGLE_COEFFS_LIBRARY,              &
!              BONDBOND_COEFFS_LIBRARY,           &
!              BONDANGLE_COEFFS_LIBRARY,          &
!              DIHEDRAL_COEFFS_LIBRARY,           &
!              MIDDLEBONDTORSION_COEFFS_LIBRARY,  &
!              ENDBONDTORSION_COEFFS_LIBRARY,     &
!              ANGLETORSION_COEFFS_LIBRARY,       &
!              ANGLEANGLETORSION_COEFFS_LIBRARY,  &
!              BONDBOND13_COEFFS_LIBRARY,         &
!              IMPROPER_COEFFS_LIBRARY,           &
!              ANGLEANGLE_COEFFS_LIBRARY);


   if ( ( IMAX_NATOM > 0 ) .AND. ( NFILE_LIBRARY > 0 ) ) then;                      !

        deallocate(CONFIG_ATOMID_LIBRARY);               !
        deallocate(CONFIG_ATOM_TYPE_LIBRARY);            !
        deallocate(CONFIG_QI_LIBRARY);                   !
        deallocate(CONFIG_VI_LIBRARY);               !
        deallocate(CONFIG_RI_LIBRARY);               !
                                                                                         !
        deallocate(CONFIG_NAT_LIBRARY);                  !

    end if                                                                               !
                                                                                         !
    if ( ( IMAX_NBOND > 0 ) .AND. ( NFILE_LIBRARY > 0 ) ) then;                      !

        deallocate(BOND_TYPE_LIBRARY);                   !
        deallocate(BOND_ATOMID_LIBRARY);             !

    end if 

    if ( ( IMAX_NANGLE > 0 ) .AND. ( NFILE_LIBRARY > 0 ) ) then;                     ! 

        deallocate(ANGLE_TYPE_LIBRARY);                 !
        deallocate(ANGLE_ATOMID_LIBRARY);              !

    end if                                                                               !
                                                                                         !
    if ( ( IMAX_NDIHEDRAL > 0 ) .AND. ( NFILE_LIBRARY > 0 ) ) then;                  !      

        deallocate(DIHEDRAL_TYPE_LIBRARY);           !
        deallocate(DIHEDRAL_ATOMID_LIBRARY);     !

    end if                                                                               !
                                                                                         !
    if ( ( IMAX_NIMPROPER > 0 ) .AND. ( NFILE_LIBRARY > 0 ) ) then;                  ! 

        deallocate(IMPROPER_TYPE_LIBRARY);           !
        deallocate(IMPROPER_ATOMID_LIBRARY);     !

    end if                                                                               !
                                                                                         !
    if ( ( IMAX_NTYPE_NATOM > 0 ) .AND. ( NFILE_LIBRARY > 0 ) ) then;                    !
             
        deallocate(ATOM_MASSES_LIBRARY);           !
        deallocate(ATOM_LABEL_LIBRARY);            !
        deallocate(POTENTIAL_CLASS2_LIBRARY);  !

    end if                                                                               !
                                                                                         !
    if ( ( IMAX_NTYPE_BOND > 0 ) .AND. ( NFILE_LIBRARY > 0 ) ) then;                     !

        deallocate(BOND_COEFFS_LIBRARY);        !

    end if

    if ( ( IMAX_NTYPE_ANGLE > 0 ) .AND. ( NFILE_LIBRARY > 0 ) ) then;

        deallocate(ANGLE_COEFFS_LIBRARY);
        deallocate(BONDBOND_COEFFS_LIBRARY);
        deallocate(BONDANGLE_COEFFS_LIBRARY);

    end if

    if ( ( IMAX_NTYPE_DIHEDRAL > 0 ) .AND. ( NFILE_LIBRARY > 0 ) ) then;

        deallocate(DIHEDRAL_COEFFS_LIBRARY);
        deallocate(MIDDLEBONDTORSION_COEFFS_LIBRARY);
        deallocate(ENDBONDTORSION_COEFFS_LIBRARY);
        deallocate(ANGLETORSION_COEFFS_LIBRARY);
        deallocate(ANGLEANGLETORSION_COEFFS_LIBRARY);
        deallocate(BONDBOND13_COEFFS_LIBRARY);

    end if

    if ( ( IMAX_NTYPE_IMPROPER > 0 ) .AND. ( NFILE_LIBRARY > 0 ) ) then;

        deallocate(IMPROPER_COEFFS_LIBRARY);
        deallocate(ANGLEANGLE_COEFFS_LIBRARY);

    end if

    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';

!   stop;

end subroutine DEALLOCATE_LIBRARY_ARRAYS
