!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!random run of VOM
!
!
!
!
!
!
!
!
!
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


subroutine random_samples ()

use vom_sce_mod
use vom_vegwat_mod

      implicit none


      real*8,dimension(6)              :: r        !random values
      integer                          :: n        !number of iterations
      real*8,dimension(6)              :: paramset !random parameterset
      integer                          :: i_loop   !current loop
      real*8                           :: obj      !objective ncp
      CHARACTER(3)  :: str

      !read settings for optimization
      call read_shufflepar()

      !read max and min ranges of parameters
      call read_shufflevar ()

      !initialize random seed 
      call random_seed()
      call transpmodel_init_once(vom_command)

      !open file for output
      open(kfile_random_output, FILE=sfile_random_output)

      !loop for n random samples, needs parallelization
      do i_loop=1, i_iter

write(*,*) "i_loop", i_loop

         !generate random number between 0  and 1
         call random_number(r)

         !make a random parameterset
         paramset=r*(parmax-parmin)+parmin

         !run the model with the random set

         call transpmodel(paramset, vom_npar, obj, 2)




         write(kfile_random_output, *) obj

      end do

      close( kfile_random_output )


end subroutine
