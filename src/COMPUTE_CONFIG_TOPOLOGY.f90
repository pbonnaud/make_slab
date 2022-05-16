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

subroutine COMPUTE_CONFIG_TOPOLOGY(icanal,PASSA,PASSB)

!   ************************************************************************************************
!   **                     Compute the topology of the molecular configuration                    **
!   ************************************************************************************************

    use module_data_in;

    use module_library;

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

!   ************************************************************************************************

    real (kind=8) :: RIJ01, RIJ02;

    real (kind=8), dimension(1:3) :: RI, DRIJ;

!   ************************************************************************************************

    integer (kind=4) :: NBOND_TMP, IBOND_TMP;

    integer (kind=4) :: NANGLE_TMP, IANGLE_TMP;

    integer (kind=4) :: IMOLECULE;

    integer (kind=4) :: ELEMENTI, ELEMENTJ;

    integer (kind=4), allocatable, dimension(:,:) :: BOND_ATOMID_TMP;

    real (kind=8) :: MASSE_ATOMI, MASSE_ATOMJ;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the routine title ####################################################################
                                                                                         !
    CHTITLE = 'Compute the topology of the molecular configuration';                     !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Initialisation of variables ################################################################
                                                                                         !
    NBOND = 0;                                                                           !
                                                                                         !
    NANGLE = 0;                                                                          !
                                                                                         !
    NDIHEDRAL = 0;                                                                       !
                                                                                         !
    NIMPROPER = 0;                                                                       !
                                                                                         !
!   ### Set and allocate temporary arrays ##########################################################
                                                                                         !
    NBOND_TMP = 0;                                                                       !
                                                                                         !
!   ### Find covalent bonds among atoms ############################################################
                                                                                         !
    do i = 1, NATOM;                                                                     !
                                                                                         !
        RI(1:3) = CONFIG_RI(1:3,i);                                                      !
                                                                                         !
        IOSEF1 = CONFIG_ATOM_TYPE(i);                                                    !
                                                                                         !
        MASSE_ATOMI = ATOM_MASSE(IOSEF1);                                                !
                                                                                         !
        call FIND_ATOM_ELEMENTID(icanal,MASSE_ATOMI,ELEMENTI);                           !
                                                                                         !
        do j = 1, NATOM;                                                                 !
                                                                                         !
            if ( i == j ) CYCLE;                                                         !
                                                                                         !
            DRIJ(1:3) = CONFIG_RI(1:3,j) - RI(1:3);                                      ! 
                                                                                         !
            call APPLY_PBC(DRIJ(1:3),DRIJ(1:3),PASSA(1:3,1:3),PASSB(1:3,1:3));           !
                                                                                         !
            RIJ02 = DOT_PRODUCT(DRIJ(1:3),DRIJ(1:3));                                    !
                                                                                         !
            if ( RIJ02 > 4.0d0 ) CYCLE;                                                  !
                                                                                         !
            IOSEF2 = CONFIG_ATOM_TYPE(j);                                                !
                                                                                         !
            MASSE_ATOMJ = ATOM_MASSE(IOSEF2);                                            ! 
                                                                                         !
            call FIND_ATOM_ELEMENTID(icanal,MASSE_ATOMJ,ELEMENTJ);                       !
                                                                                         !
!           ### Compute the distance between atom i and atom j #####################################
                                                                                         !
            RIJ01 = DSQRT(RIJ02);                                                        ! in [A]
                                                                                         !
            ROSEF3 = RIJ01 - ELEMENT_COVALENT_BONDS(ELEMENTI,ELEMENTJ);                  !
                                                                                         !
            if ( ROSEF3 > 0.02d0  ) CYCLE;                                               !
                                                                                         !
            NBOND_TMP = NBOND_TMP + 1;                                                   !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   ### Write the temporary number of bonds ########################################################
                                                                                         !
    write(icanal,'(a21,i4,a45)') '| The algorithm found', &                              !
                                 NBOND_TMP,               &                              !
                                 ' temporary bonds'//     &                              !
                                 REPEAT(' ',28)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### Allocate arrays related to temporary bond atom ids ######################################### 
                                                                                         !
    allocate(BOND_ATOMID_TMP(1:2,1:NBOND_TMP));                                          !
                                                                                         !
!   ### Build temporary array of bond atom ids #####################################################
                                                                                         !
    IBOND_TMP = 0;                                                                       !
                                                                                         !
    do i = 1, NATOM;                                                                     !
                                                                                         !
        RI(1:3) = CONFIG_RI(1:3,i);                                                      !
                                                                                         !
        IOSEF1 = CONFIG_ATOM_TYPE(i);                                                    !
                                                                                         !
        MASSE_ATOMI = ATOM_MASSE(IOSEF1);                                                !
                                                                                         !
        call FIND_ATOM_ELEMENTID(icanal,MASSE_ATOMI,ELEMENTI);                           !
                                                                                         !
        do j = 1, NATOM;                                                                 !
                                                                                         !
            if ( i == j ) CYCLE;                                                         !
                                                                                         !
            DRIJ(1:3) = CONFIG_RI(1:3,j) - RI(1:3);                                      ! 
                                                                                         !
            call APPLY_PBC(DRIJ(1:3),DRIJ(1:3),PASSA(1:3,1:3),PASSB(1:3,1:3));           !
                                                                                         !
            RIJ02 = DOT_PRODUCT(DRIJ(1:3),DRIJ(1:3));                                    !
                                                                                         !
            if ( RIJ02 > 4.0d0 ) CYCLE;                                                  !
                                                                                         !
            IOSEF2 = CONFIG_ATOM_TYPE(j);                                                !
                                                                                         !
            MASSE_ATOMJ = ATOM_MASSE(IOSEF2);                                            ! 
                                                                                         !
            call FIND_ATOM_ELEMENTID(icanal,MASSE_ATOMJ,ELEMENTJ);                       !
                                                                                         !
!           ### Compute the distance between atom i and atom j #####################################
                                                                                         !
            RIJ01 = DSQRT(RIJ02);                                                        ! in [A]
                                                                                         !
            ROSEF3 = RIJ01 - ELEMENT_COVALENT_BONDS(ELEMENTI,ELEMENTJ);                  !
                                                                                         !
            if ( ROSEF3 > 0.02d0  ) CYCLE;                                               !
                                                                                         !
!           ### Check if the bond was already counted ##############################################
                                                                                         !
            if ( IBOND_TMP > 0 ) then;                                                   !
                                                                                         !
                IOSEF1 = 0;                                                              !
                                                                                         !
                IOSEF2 = 0;                                                              !
                                                                                         !
                do k = 1, IBOND_TMP;                                                     !
                                                                                         !
                    if ( BOND_ATOMID_TMP(1,k) == i ) IOSEF1 = 1;                         !
                                                                                         !
                    if ( BOND_ATOMID_TMP(2,k) == i ) IOSEF1 = 2;                         !
                                                                                         !
                    if ( IOSEF1 == 0 ) CYCLE;                                            !
                                                                                         !
                    if ( IOSEF1 == 1 ) then;                                             !
                                                                                         !
                        if ( BOND_ATOMID_TMP(2,k) == j ) IOSEF2 = 2;                     !
                                                                                         !
                    end if                                                               !
                                                                                         !
                    if ( IOSEF1 == 2 ) then;                                             !
                                                                                         !
                        if ( BOND_ATOMID_TMP(1,k) == j ) IOSEF2 = 1;                     !
                                                                                         !
                    end if                                                               !
                                                                                         !
                    if ( ( IOSEF1 > 0 ) .AND. ( IOSEF2 > 0 ) ) EXIT;                     !
                                                                                         !
                end do                                                                   !
                                                                                         !
                if ( ( IOSEF1 > 0 ) .AND. ( IOSEF2 > 0 ) ) CYCLE;                        !
                                                                                         !
            end if                                                                       !
                                                                                         !
!           ### Update bond atom id arrays if the bond was not found ###############################
                                                                                         !
            IBOND_TMP = IBOND_TMP + 1;                                                   !
                                                                                         !
            BOND_ATOMID_TMP(1,IBOND_TMP) = i;                                            ! 
                                                                                         !
            BOND_ATOMID_TMP(2,IBOND_TMP) = j;                                            !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a14,i4,a52)') '| IBOND_TMP : ', IBOND_TMP, REPEAT(' ',51)//'|';       !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### Update and write the number of bonds in the final molecular configuration ##################
                                                                                         !
    NBOND = IBOND_TMP;                                                                   !    
                                                                                         !
    write(icanal,'(a19,i8,a43)') '| NBOND          : ', &                                !
                                 NBOND,                 &                                !
                                 REPEAT(' ',42)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### Allocate arrays related to bond atom ids for the final molecular configuration #############
                                                                                         !
    if ( NBOND > 0 ) then;                                                               !
                                                                                         !
        allocate(BOND_TYPE(1:NBOND));                                                    !
                                                                                         !
        allocate(BOND_ATOMID(1:2,1:NBOND));                                              !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Update arrays related to bond atom ids for the final molecular configuration ###############
                                                                                         !
    if ( NBOND > 0 ) then;                                                               !
                                                                                         !
        BOND_TYPE(1:NBOND) = 0;                                                          !
                                                                                         !
        BOND_ATOMID(1:2,1:NBOND) = BOND_ATOMID_TMP(1:2,1:NBOND);                         !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Initialization of temporary number of angles ###############################################
                                                                                         !
    NANGLE_TMP = 0;                                                                      !
                                                                                         !
!   ### Counting the temporary number of angle in the molecular configuration ######################
                                                                                         !
    do i = 1, NBOND-1;                                                                   !
                                                                                         !
        IOSEF1 = BOND_ATOMID(1,i);                                                       !
                                                                                         !
        IOSEF2 = BOND_ATOMID(2,i);                                                       !
                                                                                         !
        do j = i+1, NBOND;                                                               !
                                                                                         !
            IOSEF3 = 0;                                                                  !
                                                                                         !
            IOSEF4 = 0;                                                                  !
                                                                                         !
            if ( BOND_ATOMID(1,j) == IOSEF1 ) IOSEF3 = 1;                                !  
                                                                                         !
            if ( BOND_ATOMID(2,j) == IOSEF1 ) IOSEF3 = 2;                                !
                                                                                         !
            if ( BOND_ATOMID(1,j) == IOSEF2 ) IOSEF4 = 1;                                !
                                                                                         !
            if ( BOND_ATOMID(2,j) == IOSEF2 ) IOSEF4 = 2;                                !
                                                                                         !
            if ( ( IOSEF3 == 0 ) .AND. ( IOSEF4 == 0 ) ) CYCLE;                          !
                                                                                         !
            if ( ( IOSEF3 > 0 ) .AND. ( IOSEF4 > 0 ) ) CYCLE;                            !
                                                                                         !
            NANGLE_TMP = NANGLE_TMP + 1;                                                 !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a15,i4,a51)') '| NANGLE_TMP : ', NANGLE_TMP, REPEAT(' ',50)//'|';     !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### Update and and write the number of angles in the molecular configuration ###################
                                                                                         !
    NANGLE = NANGLE_TMP;                                                                 !
                                                                                         !
    write(icanal,'(a19,i8,a43)') '| NANGLE         : ', &                                !
                                 NANGLE,                &                                !
                                 REPEAT(' ',42)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### Allocate arrays related to angle atom ids of the molecular configuration ###################
                                                                                         !
    if ( NANGLE > 0 ) then;                                                              !
                                                                                         !
        allocate(ANGLE_TYPE(1:NANGLE));                                                  !
                                                                                         !
        allocate(ANGLE_ATOMID(1:3,1:NANGLE));                                            !
                                                                                         !
    end if                                                                               !
                                                                                         ! 
!   ### Initialization of arrays related to angle atom ids of the molecular configuration ##########
                                                                                         !
    if ( NANGLE > 0 ) then;                                                              !
                                                                                         !
        ANGLE_TYPE(1:NANGLE) = 0;                                                        !
                                                                                         !
        ANGLE_ATOMID(1:3,1:NANGLE) = 0;                                                  !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Build the array of angle atom ids ##########################################################
                                                                                         !
    if ( NANGLE > 0 ) then;                                                              !
                                                                                         !
        IANGLE_TMP = 0;                                                                  !
                                                                                         !
        do i = 1, NBOND-1;                                                               !
                                                                                         !
            IOSEF1 = BOND_ATOMID(1,i);                                                   !
                                                                                         !
            IOSEF2 = BOND_ATOMID(2,i);                                                   !
                                                                                         !
            do j = i+1, NBOND;                                                           !
                                                                                         !
                IOSEF3 = 0;                                                              !
                                                                                         !
                IOSEF4 = 0;                                                              !
                                                                                         !
                if ( BOND_ATOMID(1,j) == IOSEF1 ) IOSEF3 = 1;                            !  
                                                                                         !
                if ( BOND_ATOMID(2,j) == IOSEF1 ) IOSEF3 = 2;                            !
                                                                                         !
                if ( BOND_ATOMID(1,j) == IOSEF2 ) IOSEF4 = 1;                            !
                                                                                         !
                if ( BOND_ATOMID(2,j) == IOSEF2 ) IOSEF4 = 2;                            !
                                                                                         !
                if ( ( IOSEF3 == 0 ) .AND. ( IOSEF4 == 0 ) ) CYCLE;                      !
                                                                                         !
                if ( ( IOSEF3 > 0 ) .AND. ( IOSEF4 > 0 ) ) CYCLE;                        !
                                                                                         !
                IANGLE_TMP = IANGLE_TMP + 1;                                             !
                                                                                         !
                if ( IOSEF3 > 0 ) then;                                                  !
                                                                                         !
                    ANGLE_ATOMID(1,IANGLE_TMP) = IOSEF2;                                 !
                                                                                         !
                    ANGLE_ATOMID(2,IANGLE_TMP) = IOSEF1;                                 !             
                                                                                         !
                    if ( IOSEF3 == 1 ) ANGLE_ATOMID(3,IANGLE_TMP) = BOND_ATOMID(2,j);    !  
                                                                                         !
                    if ( IOSEF3 == 2 ) ANGLE_ATOMID(3,IANGLE_TMP) = BOND_ATOMID(1,j);    !  
                                                                                         !
                end if                                                                   !
                                                                                         !
                if ( IOSEF4 > 0 ) then;                                                  !
                                                                                         !
                    ANGLE_ATOMID(1,IANGLE_TMP) = IOSEF1;                                 !
                                                                                         !
                    ANGLE_ATOMID(2,IANGLE_TMP) = IOSEF2;                                 !             
                                                                                         !
                    if ( IOSEF4 == 1 ) ANGLE_ATOMID(3,IANGLE_TMP) = BOND_ATOMID(2,j);    !  
                                                                                         !
                    if ( IOSEF4 == 2 ) ANGLE_ATOMID(3,IANGLE_TMP) = BOND_ATOMID(1,j);    !  
                                                                                         !
                end if                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### To be coded atom ids for dihedrals and inpropers ######

    NDIHEDRAL = 0;

    NIMPROPER = 0;

!   ### Set molecule ids ###########################################################################
                                                                                         !
    IMOLECULE = 0;                                                                       ! 
                                                                                         !
    if ( NBOND > 0 ) then;                                                               !
                                                                                         !
        do i = 1, NBOND;                                                                 !
                                                                                         !
            IOSEF1 = BOND_ATOMID(1,i);                                                   !
                                                                                         !
            IOSEF2 = BOND_ATOMID(2,i);                                                   !
                                                                                         !
            if ( ( CONFIG_MOLECULEID(IOSEF1) == 0 ) .AND. &                              !
                 ( CONFIG_MOLECULEID(IOSEF2) == 0 ) ) then;                              !
                                                                                         !
                IMOLECULE = IMOLECULE + 1;                                               !
                                                                                         !
                CONFIG_MOLECULEID(IOSEF1) = IMOLECULE;                                   !
                                                                                         !
                CONFIG_MOLECULEID(IOSEF2) = IMOLECULE;                                   !
                                                                                         !
            else if ( ( CONFIG_MOLECULEID(IOSEF1) >  0 ) .AND. &                         !
                      ( CONFIG_MOLECULEID(IOSEF2) == 0 ) ) then;                         !
                                                                                         !
                CONFIG_MOLECULEID(IOSEF2) = CONFIG_MOLECULEID(IOSEF1);                   !
                                                                                         !
            else if ( ( CONFIG_MOLECULEID(IOSEF1) == 0 ) .AND. &                         !
                      ( CONFIG_MOLECULEID(IOSEF2) >  0 ) ) then;                         !
                                                                                         !
                CONFIG_MOLECULEID(IOSEF1) = CONFIG_MOLECULEID(IOSEF2);                   !
                                                                                         !
            else if ( ( CONFIG_MOLECULEID(IOSEF1) > 0 ) .AND. &                          !
                      ( CONFIG_MOLECULEID(IOSEF2) > 0 ) ) then;                          !
                                                                                         !
                if ( CONFIG_MOLECULEID(IOSEF1) < CONFIG_MOLECULEID(IOSEF2) ) then; 

                    CONFIG_MOLECULEID(IOSEF2) = CONFIG_MOLECULEID(IOSEF1);

                else if ( CONFIG_MOLECULEID(IOSEF1) > CONFIG_MOLECULEID(IOSEF2) ) then;

                    CONFIG_MOLECULEID(IOSEF1) = CONFIG_MOLECULEID(IOSEF2);

                end if

            end if

        end do                                                                           !

        IOSEF1 = MAXVAL(CONFIG_MOLECULEID(1:NATOM));

        if ( IOSEF1 > 0 ) ATOMS_FLAG = 13;

        write(icanal,*) 'There is ', IOSEF1, ' molecules in the configuration';

        write(icanal,*);

    end if

!   ### Deallocate arrays related to temporary bond atom ids #######################################
                                                                                         !
    deallocate(BOND_ATOMID_TMP);                                                         !
                                                                                         !
!   ### closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine COMPUTE_CONFIG_TOPOLOGY
