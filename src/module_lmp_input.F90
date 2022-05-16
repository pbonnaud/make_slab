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

module module_lmp_input


!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************
                                                                                         !
    integer (kind=4) :: NLMP_INPUT_FIX,         &                                        !
                        NLMP_INPUT_RUN,         &                                        !
                        NLMP_INPUT_FROZEN,      &                                        ! 
                        NLMP_INPUT_GROUP,       &                                        !
                        NLMP_INPUT_GROUP_UNION;                                          !
                                                                                         !
    integer (kind=4) :: NLMP_INPUT_GROUP_NEW,       &                                    !
                        NLMP_INPUT_GROUP_UNION_NEW;                                      !
                                                                                         !
!   ************************************************************************************************
                                                                                         !
    integer (kind=4), allocatable, dimension(:) :: LMP_INPUT_RUN_NSTEPS;                 !
                                                                                         !
    integer (kind=4), allocatable, dimension(:) :: LMP_INPUT_GROUP_NLIST,       &        !
                                                   LMP_INPUT_GROUP_UNION_NLIST;          !
                                                                                         !
    integer (kind=4), allocatable, dimension(:) :: LMP_INPUT_GROUP_NLIST_OLD,       &    !
                                                   LMP_INPUT_GROUP_UNION_NLIST_OLD;      !
                                                                                         !
!   ************************************************************************************************

    real (kind=8), allocatable, dimension(:) :: LMP_INPUT_RUN_TEMPK_INIT,  &
                                                LMP_INPUT_RUN_TEMPK_FINAL, &
                                                LMP_INPUT_RUN_TAUNH;


!   ************************************************************************************************

!   character (len=250), dimension(1:3) :: LAMMPS_INSERT_CHNAME, LAMMPS_INSERT_CHFILE_EXT;

!   ************************************************************************************************
                                                                                         !
    character (len=250), allocatable, dimension(:) :: LMP_INPUT_RUN_STYLE;               !
                                                                                         !
    character (len=250), allocatable, dimension(:) :: LMP_INPUT_GROUP_NAME,       &      !
                                                      LMP_INPUT_GROUP_UNION_NAME;        !
                                                                                         !
    character (len=250), allocatable, dimension(:) :: LMP_INPUT_GROUP_NAME_OLD,       &  !
                                                      LMP_INPUT_GROUP_UNION_NAME_OLD;    !
                                                                                         !
    character (len=250), allocatable, dimension(:) :: LMP_INPUT_RUN_GRPNAME;             !
                                                                                         !
    character (len=250), allocatable, dimension(:) :: LMP_INPUT_FROZEN_GRPNAME;          !
                                                                                         !
    character (len=250), allocatable, dimension(:,:) :: LMP_INPUT_GROUP_LIST,       &    !
                                                        LMP_INPUT_GROUP_UNION_LIST;      !
                                                                                         !
    character (len=250), allocatable, dimension(:,:) :: LMP_INPUT_GROUP_LIST_OLD,      & !
                                                        LMP_INPUT_GROUP_UNION_LIST_OLD;  !
                                                                                         !


end module module_lmp_input
