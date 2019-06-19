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



subroutine read_commandline()
      use vom_vegwat_mod
      implicit none

     integer :: num_args, ix
     character(len=100), dimension(:), allocatable :: args

     num_args = command_argument_count()
     allocate(args(num_args))  ! I've omitted checking the return status of the allocation 

     do ix = 1, num_args
         call get_command_argument(ix,args(ix))
         ! now parse the argument as you wis
     end do


     do ix = 1, num_args
        if(args(ix) .eq. "-i") then
           i_inputpath = args(ix+1)
           write(*,*) "Changed inputpath to:", i_inputpath
        end if

        if(args(ix) .eq. "-o") then
           i_outputpath = args(ix+1)
           write(*,*) "Changed outputpath to:",i_outputpath
        end if

     end do

end subroutine read_commandline
