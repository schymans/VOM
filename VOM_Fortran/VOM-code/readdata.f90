!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Read data needed to run the VOM
!   Parameter optimisation algorithm based on a paper by Duan et al.
!   (1993, J. Opt. Theory and Appl., 76, 501--521).
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Currently only includes the reading of command line arguments
! Needs to contain all subroutines that read data in the future
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Copyright (C) 2008 Stan Schymanski
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program. If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



subroutine read_commandline(outputpath_tmp, inputpath_tmp, change_in, change_out )
      use vom_vegwat_mod
      implicit none

     CHARACTER*100, INTENT(out)                    :: outputpath_tmp ! Temporary outputpath 
     CHARACTER*100, INTENT(out)                    :: inputpath_tmp  ! Temporary inputpath
     LOGICAL, INTENT(out)                          :: change_in      ! Change input true/false
     LOGICAL, INTENT(out)                          :: change_out     ! Change output true/false

     CHARACTER(len=100), DIMENSION(:), ALLOCATABLE :: args           ! array for arguments
     INTEGER                                       :: ix             ! Counter
     INTEGER                                       :: num_args       ! Number of arguments


     !count arguments
     num_args = command_argument_count()

     !create array to save arguments
     allocate(args(num_args))  

     do ix = 1, num_args
         !loop over arguments and save them
         call get_command_argument(ix,args(ix))
     end do

     change_in  = .FALSE.
     change_out = .FALSE.

     !loop over saved arguments and check flags
     do ix = 1, num_args

        !check inputpath
        if(args(ix) .eq. "-i") then
           inputpath_tmp = args(ix+1)
           change_in = .True.
        end if

        !check outputpath
        if(args(ix) .eq. "-o") then
           outputpath_tmp = args(ix+1)
           change_out = .True.
        end if

        !check namelist and read again if needed
        if(args(ix) .eq. "-n") then
           sfile_namelist = args(ix+1)
           write(*,*) "Changed vom_namelist:", sfile_namelist
        end if

     end do





end subroutine read_commandline
