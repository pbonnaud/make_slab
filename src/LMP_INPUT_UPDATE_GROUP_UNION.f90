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

subroutine LMP_INPUT_UPDATE_GROUP_UNION(icanal) 

!   ************************************************************************************************
!   **                                      Make tethered atoms                                   **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal : CANAL ON WHICH OUTPUT DATA ARE WRITTEN                                            **
!   **                                                                                            **
!   ************************************************************************************************

    use module_data_in;

    use module_physical_constants;

    use module_size_arrays;

    use module_library;

    use module_config;

    use module_simulation_box;

    use module_regions;

    use module_inserts;

    use module_lmp_input;

    use module_osef;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

!   integer (kind=4), intent(in) :: ILOCAL_INSERTS, KLOCAL_MODIFIED;

!   ************************************************************************************************

    integer (kind=4) :: i, j, k, m;

    integer (kind=4) :: LOCAL_NCOUNT;

!   integer (kind=4) :: LOCAL_NTYPE_ATOM_MAX, LOCAL_NTYPE_ATOM_MIN, LOCAL_NNTYPE_ATOM;

!   integer (kind=4) :: LOCAL_NTETHER;

!   integer (kind=4) :: IATOM_GHOST, IBOND_GHOST;

    real (kind=8) :: grnd;

!   real (kind=8) :: LOCAL_EPSILONIJ, LOCAL_SIGMAIJ;

!   real (kind=8), dimension(1:3) :: RINO;

!   ************************************************************************************************

!   integer (kind=4) :: NATOM_NEW, NBOND_NEW, NANGLE_NEW, NDIHEDRAL_NEW, NIMPROPER_NEW;

!   integer (kind=4) :: NTYPE_ATOM_NEW,        &
!                       NTYPE_BOND_NEW,        &
!                       NTYPE_ANGLE_NEW,       &
!                       NTYPE_DIHEDRAL_NEW,    &
!                       NTYPE_IMPROPER_NEW,    &
!                       NPAIR_COEFF_CROSS_NEW;

!   integer (kind=4) :: NLAMMPS_INSERTS_NEW;

!   ************************************************************************************************

    integer (kind=4), allocatable, dimension(:) :: NFOUND;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the routine #############################################################
                                                                                         !
    CHTITLE = 'Update group unions of the lmp_input command';                            !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate local arrays for testing if union groups need to be modified ######################
                                                                                         !
    if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
        allocate(NFOUND(1:NLMP_INPUT_GROUP_UNION));                                      !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Initialization of local arrays for testing if union groups need to be modified #############
                                                                                         !
    if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
        NFOUND(1:NLMP_INPUT_GROUP_UNION) = 0;                                            !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Look for group unions where atoms were tethered ############################################
                                                                                         !
    if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
        do i = 1, NLMP_INPUT_GROUP_UNION;                                                !
                                                                                         !
            IOSEF1 = 0;                                                                  !
                                                                                         !
            do j = 1, LMP_INPUT_GROUP_UNION_NLIST(i);                                    !
                                                                                         !
                CHOSEF1 = TRIM(LMP_INPUT_GROUP_UNION_LIST(j,i));                         !
                                                                                         !
                do k = 1, NLMP_INPUT_GROUP;                                              !
                                                                                         !
                    CHOSEF2 = TRIM(LMP_INPUT_GROUP_NAME(k));                             ! Check that the group name belong to the list of name in the group union command
                                                                                         !
                    IOSEF2 = INDEX(CHOSEF2,'_ghost');                                    ! 
                                                                                         !
                    if ( IOSEF2 > 0 ) CYCLE;                                             ! If the group name is a group of ghost atoms, then cycle
                                                                                         !
                    if ( TRIM(CHOSEF1) /= TRIM(CHOSEF2) ) CYCLE;                         ! If the group name does not belong to the list of group union commands, then cycle
                                                                                         !
                    CHOSEF3 = TRIM(CHOSEF2)//'_ghost';                                   !
                                                                                         !
                    do m = 1, NLMP_INPUT_GROUP;                                          !
                                                                                         !
                        if ( TRIM(CHOSEF3) /= TRIM(LMP_INPUT_GROUP_NAME(m)) ) CYCLE;     ! If no group of ghost atoms were created for the current group, then cycle
                                                                                         !
                        IOSEF1 = IOSEF1 + 1;                                             !
                                                                                         !
                        EXIT;                                                            !
                                                                                         !
                    end do                                                               !
                                                                                         !
                    EXIT;                                                                !
                                                                                         !
                end do                                                                   !
                                                                                         !
            end do                                                                       !
                                                                                         !
            if ( IOSEF1 /= LMP_INPUT_GROUP_UNION_NLIST(i) ) CYCLE;                       !
                                                                                         !
            NFOUND(i) = 1;                                                               !
                                                                                         !
        end do                                                                           ! 
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Count and write the number of group union command that need to be modified #################
                                                                                         !
    if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
        LOCAL_NCOUNT = SUM(NFOUND(1:NLMP_INPUT_GROUP_UNION));                            !
                                                                                         !
        write(icanal,'(a17,i4,a49)') '| LOCAL_NCOUNT : ', &                              !
                                     LOCAL_NCOUNT,        &                              !
                                     REPEAT(' ',48)//'|';                                !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  ! 
                                                                                         !
    end if                                                                               !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate local arrays for the lmp_input group union command properties #####################
                                                                                         !
    if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
        allocate(LMP_INPUT_GROUP_UNION_NAME_OLD(1:NLMP_INPUT_GROUP_UNION));              !
                                                                                         !
        allocate(LMP_INPUT_GROUP_UNION_NLIST_OLD(1:NLMP_INPUT_GROUP_UNION));             !
                                                                                         !
        allocate(LMP_INPUT_GROUP_UNION_LIST_OLD(1:NLAMMPS_INSERTS,          &            !
                                                1:NLMP_INPUT_GROUP_UNION));              !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Update local arrays for the definition of group unions in the lammps simulation ############
                                                                                         !
    if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
        LMP_INPUT_GROUP_UNION_NAME_OLD(1:NLMP_INPUT_GROUP_UNION) = &                     !
        LMP_INPUT_GROUP_UNION_NAME(1:NLMP_INPUT_GROUP_UNION);                            !
                                                                                         !
        LMP_INPUT_GROUP_UNION_NLIST_OLD(1:NLMP_INPUT_GROUP_UNION) = &                    !
        LMP_INPUT_GROUP_UNION_NLIST(1:NLMP_INPUT_GROUP_UNION);                           !
                                                                                         !
        LMP_INPUT_GROUP_UNION_LIST_OLD(1:NLAMMPS_INSERTS,1:NLMP_INPUT_GROUP_UNION) = &   !
        LMP_INPUT_GROUP_UNION_LIST(1:NLAMMPS_INSERTS,1:NLMP_INPUT_GROUP_UNION);          !  
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Deallocate glocal arrays ###################################################################
                                                                                         !
    if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
        deallocate(LMP_INPUT_GROUP_UNION_NAME);                                          !
                                                                                         !
        deallocate(LMP_INPUT_GROUP_UNION_NLIST);                                         !
                                                                                         !
        deallocate(LMP_INPUT_GROUP_UNION_LIST);                                          !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Initialization of the size of new arrays ###################################################
                                                                                         !
    if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
        NLMP_INPUT_GROUP_UNION_NEW = NLMP_INPUT_GROUP_UNION;                             !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Modify the size of new arrays for including the insertion of tethered atoms ################
                                                                                         !
    if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
        NLMP_INPUT_GROUP_UNION_NEW = NLMP_INPUT_GROUP_UNION_NEW + LOCAL_NCOUNT;          !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Re-allocate global arrays ##################################################################
                                                                                         !
    if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
        allocate(LMP_INPUT_GROUP_UNION_NAME(1:NLMP_INPUT_GROUP_UNION_NEW));              !
                                                                                         !
        allocate(LMP_INPUT_GROUP_UNION_NLIST(1:NLMP_INPUT_GROUP_UNION_NEW));             !
                                                                                         !
        allocate(LMP_INPUT_GROUP_UNION_LIST(1:NLAMMPS_INSERTS,              &            !
                                            1:NLMP_INPUT_GROUP_UNION_NEW));              !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Update global arrays #######################################################################
                                                                                         !
    if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
        LMP_INPUT_GROUP_UNION_NAME(1:NLMP_INPUT_GROUP_UNION) =    &                      !
        LMP_INPUT_GROUP_UNION_NAME_OLD(1:NLMP_INPUT_GROUP_UNION);                        !
                                                                                         !
        LMP_INPUT_GROUP_UNION_NLIST(1:NLMP_INPUT_GROUP_UNION) =    &                     !
        LMP_INPUT_GROUP_UNION_NLIST_OLD(1:NLMP_INPUT_GROUP_UNION);                       !
                                                                                         !
        LMP_INPUT_GROUP_UNION_LIST(1:NLAMMPS_INSERTS,1:NLMP_INPUT_GROUP_UNION) =    &    !
        LMP_INPUT_GROUP_UNION_LIST_OLD(1:NLAMMPS_INSERTS,1:NLMP_INPUT_GROUP_UNION);      !  
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Deallocate local arrays ####################################################################
                                                                                         !
    if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
        deallocate(LMP_INPUT_GROUP_UNION_NAME_OLD);                                      !
                                                                                         !
        deallocate(LMP_INPUT_GROUP_UNION_NLIST_OLD);                                     !
                                                                                         !
        deallocate(LMP_INPUT_GROUP_UNION_LIST_OLD);                                      !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Create lmp_input command properties for ghost atoms ########################################
                                                                                         !
    if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
        IOSEF1 = NLMP_INPUT_GROUP_UNION;                                                 !
                                                                                         !
        do i = 1, NLMP_INPUT_GROUP_UNION;                                                !
                                                                                         !
            if ( NFOUND(i) == 0 ) CYCLE;                                                 !
                                                                                         !
            IOSEF1 = IOSEF1 + 1;                                                         !
                                                                                         !
            LMP_INPUT_GROUP_UNION_NAME(IOSEF1) =           &                             !
            TRIM(LMP_INPUT_GROUP_UNION_NAME(i))//'_ghost';                               !
                                                                                         !
            LMP_INPUT_GROUP_UNION_NLIST(IOSEF1) = &                                      !
            LMP_INPUT_GROUP_UNION_NLIST(i);                                              !
                                                                                         !
            do j = 1, LMP_INPUT_GROUP_UNION_NLIST(i);                                    !
                                                                                         !
                LMP_INPUT_GROUP_UNION_LIST(j,IOSEF1) =           &                       !
                TRIM(LMP_INPUT_GROUP_UNION_LIST(j,i))//'_ghost';                         !
                                                                                         !
!               write(icanal,*) TRIM(LMP_INPUT_GROUP_UNION_LIST(j,IOSEF1))

            end do                                                                       !
                                                                                         !
        end do                                                                           !
                                                                                         !
    end if                                                                               !
                                                                                         !
!   ### Update the size of arrays for the lmp_input command properties for ghost atoms #############
                                                                                         !
    if ( NLMP_INPUT_GROUP_UNION > 0 ) then;                                              !
                                                                                         !
        NLMP_INPUT_GROUP_UNION = NLMP_INPUT_GROUP_UNION_NEW;                             !
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
end subroutine LMP_INPUT_UPDATE_GROUP_UNION

