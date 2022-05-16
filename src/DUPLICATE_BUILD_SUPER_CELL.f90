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

subroutine DUPLICATE_BUILD_SUPER_CELL(icanal,PASSA,PASSB)

!   ************************************************************************************************
!   **                               DUPLICATE THE ORIGINAL UNIT CELL                             **
!   ************************************************************************************************

    use module_data_in;

    use module_duplicate;

    use module_config;

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

    real (kind=8), dimension(1:3,1:3), intent(in) :: PASSA, PASSB;

!   ************************************************************************************************

    integer (kind=4) :: i, j, k, m;

    integer (kind=4) :: IATOM;

    real (kind=8), dimension(1:3) :: TRANSR, RINO;

    integer (kind=4) :: INIT_NATOM;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the routine title ####################################################################
                                                                                         !
    CHTITLE = 'Build super cell (duplicate)';                                            !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
    ATOMS_FLAG = 5;   ! Atoms # charge

!   ### Set and write the number of atoms in the supercell #########################################
                                                                                         !
    NATOM = NATOM_DUPLICATE * DUPLICATE_NX * DUPLICATE_NY * DUPLICATE_NZ;                !
                                                                                         !
    write(icanal,'(a10,i8,a52)') '| NATOM : ', NATOM, REPEAT(' ',51)//'|';               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////! 
                                                                                         !
!   ### Allocate arrays of the final molecular configuration ####################################### 
                                                                                         !
    if ( NATOM > 0 ) then;                                                               !
                                                                                         !
        allocate(CONFIG_ATOMID(1:NATOM));                                                !
                                                                                         !
        allocate(CONFIG_MOLECULEID(1:NATOM));                                            !
                                                                                         !
        allocate(CONFIG_ATOM_TYPE(1:NATOM));                                             !
                                                                                         !
        allocate(CONFIG_QI(1:NATOM));                                                    !
                                                                                         !
        allocate(CONFIG_VI(1:3,1:NATOM));                                                !
                                                                                         !
        allocate(CONFIG_RI(1:3,1:NATOM));                                                !
                                                                                         !
        allocate(CONFIG_NAT(1:NATOM));                                                   !
                                                                                         !
    end if                                                                               !
                                                                                         !
    write(icanal,'(a70)') '| Configuration arrays were allocated '//REPEAT(' ',31)//'|'; !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////! 
                                                                                         !
!   ### Initialization of arrays for the final molecular configuration #############################
                                                                                         !
    if ( NATOM > 0 ) then;                                                               !
                                                                                         !
        CONFIG_ATOMID(1:NATOM) = 0;                                                      !
                                                                                         !
        CONFIG_MOLECULEID(1:NATOM) = 0;                                                  !
                                                                                         !
        CONFIG_ATOM_TYPE(1:NATOM) = 0;                                                   !
                                                                                         !
        CONFIG_QI(1:NATOM) = 0.0d0;                                                      !
                                                                                         !
        CONFIG_VI(1:3,1:NATOM) = 0.0d0;                                                  !
                                                                                         ! 
        CONFIG_RI(1:3,1:NATOM) = 0.0d0;                                                  !
                                                                                         !
        CONFIG_NAT(1:NATOM) = 'XXX';                                                     !
                                                                                         !
    end if                                                                               !
                                                                                         !
    write(icanal,'(a70)') '| Configuration arrays were initialized '// &                 !
                          REPEAT(' ',29)//'|';                                           !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////! 
                                                                                         !
!   ### Set constants for the duplication ##########################################################
                                                                                         !
    ROSEF1 = 1.0d0 / REAL(DUPLICATE_NX);                                                 !
                                                                                         !
    ROSEF2 = 1.0d0 / REAL(DUPLICATE_NY);                                                 !
                                                                                         !
    ROSEF3 = 1.0d0 / REAL(DUPLICATE_NZ);                                                 !
                                                                                         !
    ROSEF4 = - 0.5d0 * ROSEF1 * ( REAL(DUPLICATE_NX) - 1.0d0 );                          !
                                                                                         !
    ROSEF5 = - 0.5d0 * ROSEF2 * ( REAL(DUPLICATE_NY) - 1.0d0 );                          !
                                                                                         !
    ROSEF6 = - 0.5d0 * ROSEF3 * ( REAL(DUPLICATE_NZ) - 1.0d0 );                          !
                                                                                         !
!   ### Duplicate coordinates ######################################################################
                                                                                         !
    IATOM = 0;                                                                           !
                                                                                         !
    TRANSR(1) = ROSEF4;                                                                  !
                                                                                         !
    do k = 1, DUPLICATE_NX;                                                              !
                                                                                         !
        TRANSR(2) = ROSEF5;                                                              !
                                                                                         !
        do j = 1, DUPLICATE_NY;                                                          !
                                                                                         !
           TRANSR(3) = ROSEF6;                                                           !
                                                                                         !
            do i = 1, DUPLICATE_NZ;                                                      !
                                                                                         !
                do m = 1, NATOM_DUPLICATE;                                               !
                                                                                         !
                    IATOM = IATOM + 1;                                                   !
                                                                                         !
                    CONFIG_ATOMID(IATOM) = IATOM;                                        !
                                                                                         !
                    CONFIG_NAT(IATOM) = CONFIG_NAT_DUPLICATE(m);                         !
                                                                                         !
                    CONFIG_ATOM_TYPE(IATOM) = CONFIG_ATOM_TYPE_DUPLICATE(m);             !
                                                                                         !
                    CONFIG_RI(1:3,IATOM) = CONFIG_RI_DUPLICATE(1:3,m);                   !
                                                                                         !
                    RINO(1:3) = MATMUL(PASSB(1:3,1:3),CONFIG_RI(1:3,IATOM));             !
                                                                                         !
                    RINO(1:3) = RINO(1:3) + TRANSR(1:3);                                 !
                                                                                         !
                    CONFIG_RI(1:3,IATOM) = MATMUL(PASSA(1:3,1:3),RINO(1:3));             !
                                                                                         !
                end do                                                                   !
                                                                                         ! 
                if ( DUPLICATE_NZ /= 1 ) TRANSR(3) = TRANSR(3) + ROSEF3;                 !
                                                                                         !
            end do                                                                       !
                                                                                         !
            TRANSR(2) = TRANSR(2) + ROSEF2;                                              !
                                                                                         !
        end do                                                                           !
                                                                                         !
        TRANSR(1) = TRANSR(1) + ROSEF1;                                                  !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   ### Update configuration properties for the molecular configuration ############################
                                                                                         !
    NTYPE_ATOM = NTYPE_ATOM_DUPLICATE;                                                   !
                                                                                         !
!   ### Write configuration properties for the molecular configuration #############################
                                                                                         !
    write(icanal,'(a15,i4,a51)') '| NTYPE_ATOM : ', NTYPE_ATOM, REPEAT(' ',50)//'|';     !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### Allocate arrays related to configuration properties ########################################
                                                                                         !
    if ( NTYPE_ATOM > 0 ) then;                                                          !
                                                                                         !
        allocate(ATOM_MASSE(1:NTYPE_ATOM));                                              !
                                                                                         !
        allocate(ATOM_LABEL(1:NTYPE_ATOM));                                              !
                                                                                         !
        write(icanal,'(a70)') '| Lists of atom masses and labels were allocated'// &     !
                              REPEAT(' ',21)//'|';                                       !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Initialization of arrays related to configuration properties ###############################
                                                                                         !
    if ( NTYPE_ATOM > 0 ) then;                                                          !
                                                                                         !
        ATOM_MASSE(1:NTYPE_ATOM)  = 0.0d0;                                               !
                                                                                         !
        ATOM_LABEL(1:NTYPE_ATOM)  = 'XXX';                                               !
                                                                                         !
        write(icanal,'(a70)') '| Lists of atom masses and labels were initialized'// &   !
                              REPEAT(' ',19)//'|';                                       !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Update arrays related to configuration properties ##########################################
                                                                                         !
    if ( NTYPE_ATOM > 0 ) then;                                                          !
                                                                                         !
        ATOM_MASSE(1:NTYPE_ATOM) = ATOM_MASSES_DUPLICATE(1:NTYPE_ATOM);                  !
                                                                                         !
        ATOM_LABEL(1:NTYPE_ATOM) = ATOM_LABEL_DUPLICATE(1:NTYPE_ATOM);                   !
                                                                                         !
        write(icanal,'(a70)') '| Lists of atom masses and labels were updated'// &       !
                              REPEAT(' ',23)//'|';                                       !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Update force field properties ##############################################################
                                                                                         !
    IPOTENTIAL_CLASS2 = 0;                                                               ! 
                                                                                         !
!   ### closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine DUPLICATE_BUILD_SUPER_CELL
