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

subroutine CHECK_ELECTRONEUTRALITY(icanal);

!   ************************************************************************************************
!   **                                                                                            **
!   ** ADD_NLABEL : NUMBER DIFFERENT LABEL TO MODIFY 
!   ** ADD_LABEL  : LABEL TO MODIFY
!   ** ADD_NIONS  : NUMBER OF IONS PER LABEL TO INSERT
!   **
!   ************************************************************************************************

    use module_data_in;

    use module_physical_constants;

    use module_size_arrays;

!   use module_library;

    use module_config;

!   ************************************************************************************************

    implicit none;

!   ************************************************************************************************

    integer (kind=4), intent(in) :: icanal;

!   ************************************************************************************************

    integer (kind=4) :: h, i, j; 

    integer (kind=4) :: IMAX_MOLECULEID, NCHARGED_MOLECULES, NCHARGED_ATOMS;

!   integer (kind=4) :: JLABEL, NNEWSPEC;

!   character (len=3), dimension(1:10) :: ATOM_LABEL;

    real (kind=8) :: SUM_CONFIG, LOCAL_CHARGE_MOLECULE, LOCAL_CORRECTING_CHARGE;

!, SUM_CONFIG_NEW, SUM_ADDED_IONS, ADDED_IONS;

!   real (kind=8) :: QNEWSPEC, QMAX, DIFF;

!   real (kind=8), dimension(1:NLABEL) :: SUM_QLABEL;

    integer (kind=4), allocatable, dimension(:) :: LOCAL_CONFIG_CHARGED;

!   ************************************************************************************************

    integer (kind=4) :: IOSEF1;

!   ************************************************************************************************

    integer (kind=4) :: ILENGTH_TITLE;

    character (len=250) :: CHTITLE;

!   ************************************************************************************************
                                                                                         !
!   ### Write the title of the routine #############################################################
                                                                                         !
    CHTITLE = 'Check electroneutrality';                                                 !
                                                                                         !
    ILENGTH_TITLE = LEN_TRIM(CHTITLE);                                                   !
                                                                                         !
    call MAKE_ROUTINE_TITLE(icanal,70,ILENGTH_TITLE,TRIM(CHTITLE));                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Compute the net charge of the molecular configuration ######################################
                                                                                         !
    SUM_CONFIG = SUM(CONFIG_QI(1:NATOM));                                                !
                                                                                         !
    write(icanal,'(a19,f15.6,a36)') '| SUM_CONFIG [e] : ', &                             !
                                    SUM_CONFIG,            &                             !
                                    REPEAT(' ',35)//'|';                                 !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Get the total number of molecules in the molecular configuration ###########################
                                                                                         !
    IMAX_MOLECULEID = MAXVAL(CONFIG_MOLECULEID(1:NATOM));                                !
                                                                                         !
    write(icanal,'(a11,i8,a51)') '| There is ', IMAX_MOLECULEID,      &                  !
                                 ' molecules in the simulation box'// &                  !
                                 REPEAT(' ',18)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Allocate array for the finding of charged molecules ########################################
                                                                                         !
    allocate(LOCAL_CONFIG_CHARGED(1:IMAX_MOLECULEID));                                   !
                                                                                         !
!   ### Initialization of the array for finding charged molecules ##################################
                                                                                         !
    LOCAL_CONFIG_CHARGED(1:IMAX_MOLECULEID) = 0;                                         !
                                                                                         !
!   ### Compute the number of charged molecules and corresponding number of atoms ##################
                                                                                         !
    NCHARGED_MOLECULES = 0;                                                              !
                                                                                         !
    NCHARGED_ATOMS = 0;                                                                  !
                                                                                         !
    do i = 1, IMAX_MOLECULEID;                                                           !
                                                                                         !
        LOCAL_CHARGE_MOLECULE = 0.0d0;                                                   !
                                                                                         !
        IOSEF1 = 0;                                                                      !
                                                                                         !
        do j = 1, NATOM;                                                                 !
                                                                                         !
            if ( CONFIG_MOLECULEID(j) /= i ) CYCLE;                                      !
                                                                                         !
            LOCAL_CHARGE_MOLECULE = LOCAL_CHARGE_MOLECULE + CONFIG_QI(j);                !
                                                                                         !
            IOSEF1 = IOSEF1 + 1;                                                         !
                                                                                         !
        end do                                                                           !
                                                                                         !
        if ( ABS(LOCAL_CHARGE_MOLECULE) < 1.0E-6 ) CYCLE;                                !
                                                                                         ! 
        NCHARGED_MOLECULES = NCHARGED_MOLECULES + 1;                                     !
                                                                                         !
        NCHARGED_ATOMS = NCHARGED_ATOMS + IOSEF1;                                        ! 
                                                                                         !
        LOCAL_CONFIG_CHARGED(i) = INT ( SIGN(1.0d0,LOCAL_CHARGE_MOLECULE) );             !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a11,i8,a51)') '| There is ',                               &          !
                                 NCHARGED_MOLECULES,                          &          !
                                 ' charged molecules in the simulation box'// &          !
                                 REPEAT(' ',10)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a11,i8,a51)') '| There is ',                            &             !
                                 NCHARGED_ATOMS,                           &             !
                                 ' atoms belonging to charged molecules'// &             !
                                 REPEAT(' ',13)//'|';                                    !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Compute the correcting charge to make the system electroneutral ############################
                                                                                         !
    LOCAL_CORRECTING_CHARGE = 0.0d0;                                                     !
                                                                                         !
    if ( NCHARGED_ATOMS > 0 ) then;                                                      !
                                                                                         !
        LOCAL_CORRECTING_CHARGE = SUM_CONFIG / REAL(NCHARGED_ATOMS);                     ! 
                                                                                         !
    end if                                                                               !
                                                                                         !
    write(icanal,'(a27,f12.6,a31)') '| The correcting charge is ', &                     !
                                    LOCAL_CORRECTING_CHARGE,       &                     !
                                    ' [e]'//REPEAT(' ',26)//'|';                         !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Add the correcting charge to atoms belonging to charged molecules ##########################
                                                                                         !
    do i = 1, IMAX_MOLECULEID;                                                           !
                                                                                         !
        if ( LOCAL_CONFIG_CHARGED(i) == 0 ) CYCLE;                                       !
                                                                                         !
        do j = 1, NATOM;                                                                 !
                                                                                         !
            if ( CONFIG_MOLECULEID(j) /= i ) CYCLE;                                      !
                                                                                         !
            CONFIG_QI(j) = CONFIG_QI(j) - LOCAL_CORRECTING_CHARGE;                       ! 
                                                                                         !
        end do                                                                           !
                                                                                         !
    end do                                                                               !
                                                                                         !
    write(icanal,'(a70)') '| Charges were corrected'//REPEAT(' ',45)//'|';               !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### Compute the net corrected charge of the molecular configuration ############################
                                                                                         !
    SUM_CONFIG = SUM(CONFIG_QI(1:NATOM));                                                !
                                                                                         !
    write(icanal,'(a27,f15.6,a28)') '| SUM_CONFIG (corr.) [e] : ', &                     !
                                    SUM_CONFIG,                    &                     !
                                    REPEAT(' ',27)//'|';                                 !
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
!   ### deallocate local arrays ####################################################################
                                                                                         !
    deallocate(LOCAL_CONFIG_CHARGED);                                                    !
                                                                                         !
!   ### Closing the routine ########################################################################
                                                                                         !
    write(icanal,'(a70)') '|'//REPEAT(' ',68)//'|';                                      !
                                                                                         !
    write(icanal,'(a70)') '+'//REPEAT('-',68)//'+';                                      !
                                                                                         !
!   stop; !//////////////////////////////////////////////////////////////////////////////!
                                                                                         !
end subroutine CHECK_ELECTRONEUTRALITY

