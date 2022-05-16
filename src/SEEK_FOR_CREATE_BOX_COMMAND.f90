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

subroutine SEEK_FOR_CREATE_BOX_COMMAND(icanal)

!   ************************************************************************************************
!   **                            Seek for region keyword in the file                             **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                : CANAL ON WHICH OUTPUT DATA ARE WRITTEN                             **
!   ** PRIMITIVE_CELL_AXIS   : DIMENSIONS OF THE PRIMITIVE CELL                                   **
!   ** PRIMITIVE_CELL_ANGDEG : ANGLES OF THE PRIMITIVE CELL                                       **
!   **                                                                                            **
!   ************************************************************************************************

    use module_simulation_box;

    use module_osef;

    use module_regions;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

!   ************************************************************************************************

    integer (kind=4) :: h, i, j;

!   ************************************************************************************************

    integer (kind=4) :: ilocal_create_box, ilocal_region;

!   ************************************************************************************************

    character (len=250) :: CHARLINE;

!   ************************************************************************************************

    integer (kind=4) :: EOF;

!   ************************************************************************************************
                                                                                         !
!   ### Initialization of variables and arrays related to the create box command ###################
                                                                                         !
    LAMMPS_CREATE_BOX_RNAME = 'XXX';                                                     ! 
                                                                                         !
    ilocal_create_box = 0;                                                               ! 
                                                                                         !
!   ### Open the input file ########################################################################
                                                                                         !
    open(1,file='input_100.dat',status='old');                                           !
                                                                                         !
    EOF = 0;                                                                             !
                                                                                         !
    do                                                                                   !
                                                                                         !
        read(1,'(a)',iostat=EOF) CHARLINE;                                               !
                                                                                         !
        if ( EOF /= 0 ) EXIT;                                                            !
                                                                                         !
        if ( INDEX(CHARLINE,'#') > 0 ) CYCLE;                                            !
                                                                                         !
        if ( INDEX(CHARLINE,'create_box') > 0 ) then;                                    !
                                                                                         !
            if ( ilocal_create_box == 1 ) then; 

                write(icanal,*) '| [Error] More than one create box command was used in the input file (input_100.dat)'

                write(icanal,*) 'End of program';

                close(icanal);

                stop;

            end if
                                                                                         ! 
            read(CHARLINE,*) CHOSEF1, IOSEF1, LAMMPS_CREATE_BOX_RNAME;                   !
                                                                                         !
            ilocal_create_box = 1;                                                       !
                                                                                         !
        end if                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
    close(1);                                                                            !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the name of the region used for creating the simulation box ##########################
                                                                                         !
    IOSEF1 = 70 - 42 - 1 - LEN_TRIM(LAMMPS_CREATE_BOX_RNAME);                            !
                                                                                         !
    if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                        !
                                                                                         !
    write(icanal,'(a70)') '| The region used for creating the box is '// &               !
                          TRIM(LAMMPS_CREATE_BOX_RNAME)//                &               !
                          REPEAT(' ',IOSEF1)//'|';                                       !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### Find the location of the selected region in the region list ################################
                                                                                         !
    ilocal_region = 0;                                                                   !
                                                                                         !
    do i = 1, NLAMMPS_REGIONS;                                                           !
                                                                                         !
        if ( TRIM(LAMMPS_REGION_NAME(i)) /= TRIM(LAMMPS_CREATE_BOX_RNAME) ) CYCLE;       !
                                                                                         !
        ilocal_region = i;                                                               !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   ### Write the style of the selected region #####################################################
                                                                                         !
    IOSEF1 = 70 - 38 - 1 - LEN_TRIM(LAMMPS_REGION_STYLE(ilocal_region));                 !
                                                                                         !
    if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                        !
                                                                                         !
    write(icanal,'(a70)') '| The style of the selected region is '// &                   !
                          TRIM(LAMMPS_REGION_STYLE(ilocal_region))// &                   !
                          REPEAT(' ',IOSEF1)//'|';                                       !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### Define the simulation box depending on the region style ####################################
                                                                                         !
    select case(TRIM(LAMMPS_REGION_STYLE(ilocal_region)));                               !
                                                                                         !
        case('block');                                                                   !
                                                                                         !
            LAMMPS_MATA(1) = LAMMPS_REGION_ARGS(2,ilocal_region) - &                     ! IN [A]
                             LAMMPS_REGION_ARGS(1,ilocal_region);                        !
                                                                                         !
            LAMMPS_MATA(2) = LAMMPS_REGION_ARGS(4,ilocal_region) - &                     ! IN [A]
                             LAMMPS_REGION_ARGS(3,ilocal_region);                        !
                                                                                         !
            LAMMPS_MATA(3) = LAMMPS_REGION_ARGS(6,ilocal_region) - &                     ! IN [A]
                             LAMMPS_REGION_ARGS(5,ilocal_region);                        !
                                                                                         !
            LAMMPS_MATA(4:6) = 0.0d0;                                                    !
                                                                                         !
            LAMMPS_CELL_AXIS(1:3) = LAMMPS_MATA(1:3);                                    !
                                                                                         !
            LAMMPS_CELL_ANGDEG(1:3) = 90.0d0;                                            !
                                                                                         !
        case default;                                                                    !
                                                                                         !
            write(icanal,*) '| Not implemented - end of program';                        !
                                                                                         !
            close(icanal);                                                               !
                                                                                         !
            stop; !//////////////////////////////////////////////////////////////////////!
                                                                                         !
    end select                                                                           !
                                                                                         ! 
!   ### Write properties of the simulation box #####################################################
                                                                                         !
    write(icanal,'(a16,f15.6,a39)') '| MATA(1) [A] : ',  &                               !
                                    LAMMPS_MATA(1),      &                               !
                                    REPEAT(' ',38)//'|';                                 !
                                                                                         !
    write(icanal,'(a16,f15.6,a39)') '| MATA(2) [A] : ',  &                               !
                                    LAMMPS_MATA(2),      &                               !
                                    REPEAT(' ',38)//'|';                                 !
                                                                                         !
    write(icanal,'(a16,f15.6,a39)') '| MATA(3) [A] : ',  &                               !
                                    LAMMPS_MATA(3),      &                               !
                                    REPEAT(' ',38)//'|';                                 !
                                                                                         !
    write(icanal,'(a16,f15.6,a39)') '| MATA(4) [A] : ',  &                               !
                                    LAMMPS_MATA(4),      &                               !
                                    REPEAT(' ',38)//'|';                                 !
                                                                                         !
    write(icanal,'(a16,f15.6,a39)') '| MATA(5) [A] : ',  &                               !
                                    LAMMPS_MATA(5),      &                               !
                                    REPEAT(' ',38)//'|';                                 !
                                                                                         !
    write(icanal,'(a16,f15.6,a39)') '| MATA(6) [A] : ',  &                               !
                                    LAMMPS_MATA(6),      &                               !
                                    REPEAT(' ',38)//'|';                                 !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !   
                                                                                         !
                                                                                         !
!   ### Write message for the success in reanding the input file ###################################
                                                                                         !
    write(icanal,'(a70)') '| The simulation box was created '//REPEAT(' ',36)//'|';      !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine SEEK_FOR_CREATE_BOX_COMMAND











