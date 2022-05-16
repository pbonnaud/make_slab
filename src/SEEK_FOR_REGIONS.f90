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

subroutine SEEK_FOR_REGIONS(icanal)

!   ************************************************************************************************
!   **                            Seek for region keyword in the file                             **
!   ************************************************************************************************
!   **                                                                                            **
!   ** icanal                : CANAL ON WHICH OUTPUT DATA ARE WRITTEN                             **
!   ** PRIMITIVE_CELL_AXIS   : DIMENSIONS OF THE PRIMITIVE CELL                                   **
!   ** PRIMITIVE_CELL_ANGDEG : ANGLES OF THE PRIMITIVE CELL                                       **
!   **                                                                                            **
!   ************************************************************************************************

    use module_osef;

    use module_regions;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

!   ************************************************************************************************

    integer (kind=4) :: h, i, j;

!   ************************************************************************************************

    character (len=250) :: CHARLINE;

!   ************************************************************************************************

    integer (kind=4) :: EOF;

!   ************************************************************************************************
                                                                                         !
!   ### Initialization of arrays related to regions ################################################
                                                                                         !
    NLAMMPS_REGIONS = 0;                                                                 !
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
        IOSEF2 = INDEX(CHARLINE,'region');                                               !
                                                                                         !
        if ( IOSEF2 > 20 ) CYCLE;                                                        !
                                                                                         !
        if ( IOSEF2 > 0 ) then;                                                          !
                                                                                         !
            NLAMMPS_REGIONS = NLAMMPS_REGIONS + 1;                                       !
                                                                                         !
        end if                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
    close(1);                                                                            !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the number of regions found in the input file ########################################
                                                                                         !
    write(icanal,'(a30,i4,a36)') '| Number of defined regions : ', &                     !
                                 NLAMMPS_REGIONS,                  &                     !
                                 REPEAT(' ',35)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   ### Allocate arrays containing region properties ###############################################
                                                                                         !
    allocate(LAMMPS_REGION_NAME(1:NLAMMPS_REGIONS));                                     !
                                                                                         !
    allocate(LAMMPS_REGION_STYLE(1:NLAMMPS_REGIONS));                                    !
                                                                                         !
    allocate(LAMMPS_REGION_ARGS(1:6,1:NLAMMPS_REGIONS));                                 !
                                                                                         !
    allocate(LAMMPS_REGION_DIMENSIONS(1:3,1:NLAMMPS_REGIONS));                           !
                                                                                         !
!   ### Initialization of arrays containing region properties ######################################
                                                                                         !
    LAMMPS_REGION_NAME(1:NLAMMPS_REGIONS) = 'XXX';                                       !
                                                                                         !
    LAMMPS_REGION_STYLE(1:NLAMMPS_REGIONS) = 'XXX';                                      !
                                                                                         !
    LAMMPS_REGION_ARGS(1:6,1:NLAMMPS_REGIONS) = 0.0d0;                                   !
                                                                                         !
    LAMMPS_REGION_DIMENSIONS(1:3,1:NLAMMPS_REGIONS) = 0.0d0;                             !
                                                                                         !
!   ### Read region properties #####################################################################
                                                                                         !
    IOSEF1 = 0;                                                                          !
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
        IOSEF2 = INDEX(CHARLINE,'region');                                               !
                                                                                         !                                                                                         !
        if ( IOSEF2 > 20 ) CYCLE;                                                        !
                                                                                         !
        if ( IOSEF2 > 0 ) then;                                                          !
                                                                                         !
            IOSEF1 = IOSEF1 + 1;                                                         !
                                                                                         !
            read(CHARLINE,*) CHOSEF1,                     &                              !
                             LAMMPS_REGION_NAME(IOSEF1),  &                              !
                             LAMMPS_REGION_STYLE(IOSEF1);                                !
                                                                                         !
            if ( INDEX(LAMMPS_REGION_STYLE(IOSEF1),'block') > 0 ) then;                  !
                                                                                         !
                read(CHARLINE,*) CHOSEF1,                        &                       !
                                 CHOSEF2,                        &                       !
                                 CHOSEF3,                        &                       !
                                 LAMMPS_REGION_ARGS(1:6,IOSEF1);                         !
                                                                                         !
            end if                                                                       !
                                                                                         !
        end if                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
    close(1);                                                                            !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Write the list of regions found in the input file ##########################################
                                                                                         !
    write(icanal,'(a70)') '| List of regions '//REPEAT('>',51)//'|';                     !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    do i = 1, NLAMMPS_REGIONS;                                                           ! 
                                                                                         !
        IOSEF1 = 70 - 5 - 18 - 1 -                 &                                     !
                 LEN_TRIM(LAMMPS_REGION_NAME(i)) - &                                     !
                 LEN_TRIM(LAMMPS_REGION_STYLE(i));                                       !
                                                                                         !
        if ( IOSEF1 < 0 ) IOSEF1 = 1;                                                    !
                                                                                         !
        write(icanal,'(a70)') '|    '//                      &                           !
                              TRIM(LAMMPS_REGION_NAME(i))//  &                           !
                              ' having the style '//         &                           !
                              TRIM(LAMMPS_REGION_STYLE(i))// &                           !
                              REPEAT(' ',IOSEF1)//'|';                                   !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
        write(icanal,'(a13,2f12.4,a33)') '| xlo, xhi : ',           &                    !
                                         LAMMPS_REGION_ARGS(1:2,i), &                    !
                                         REPEAT(' ',32)//'|';                            !
                                                                                         !
        write(icanal,'(a13,2f12.4,a33)') '| ylo, yhi : ',           &                    !
                                         LAMMPS_REGION_ARGS(3:4,i), &                    !
                                         REPEAT(' ',32)//'|';                            !
                                                                                         !
        write(icanal,'(a13,2f12.4,a33)') '| zlo, zhi : ',           &                    !
                                         LAMMPS_REGION_ARGS(5:6,i), &                    !
                                         REPEAT(' ',32)//'|';                            !
                                                                                         !
        write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                  !
                                                                                         !
!       ### Set the reference point of the region depending on its style ###########################
                                                                                         !
        if ( TRIM(LAMMPS_REGION_STYLE(i)) == 'block' ) then;                             !
                                                                                         !
            LAMMPS_REGION_DIMENSIONS(1,i) = LAMMPS_REGION_ARGS(2,i) - &                  !
                                            LAMMPS_REGION_ARGS(1,i);                     !
                                                                                         !
            LAMMPS_REGION_DIMENSIONS(2,i) = LAMMPS_REGION_ARGS(4,i) - &                  !
                                            LAMMPS_REGION_ARGS(3,i);                     !
                                                                                         !
            LAMMPS_REGION_DIMENSIONS(3,i) = LAMMPS_REGION_ARGS(6,i) - &                  !
                                            LAMMPS_REGION_ARGS(5,i);                     !
                                                                                         !
            write(icanal,'(a19,3f12.5,a15)') '| Lx, Ly, Lz [A] : ',           &          !
                                             LAMMPS_REGION_DIMENSIONS(1:3,i), &          !
                                             REPEAT(' ',14)//'|';                        !
                                                                                         !
            write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                              !
                                                                                         !
        end if                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
!   ### Write message for the success in reanding the input file ###################################
                                                                                         !
    write(icanal,'(a70)') '| Regions have been read '//REPEAT(' ',44)//'|';              !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine SEEK_FOR_REGIONS











