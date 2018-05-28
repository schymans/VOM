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
      
      open(kfile_vd_d , FILE=sfile_vd_d)
      open(kfile_esoil, FILE=sfile_esoil)
      open(kfile_jmax25t, FILE=sfile_jmax25t)
      open(kfile_jmax25g, FILE=sfile_jmax25g)
      open(kfile_vegcov, FILE=sfile_vegcov)
      open(kfile_resp, FILE=sfile_resp)
      open(kfile_lambdat, FILE=sfile_lambdat)
      open(kfile_lambdag, FILE=sfile_lambdag)
      open(kfile_rrt, FILE=sfile_rrt)
      open(kfile_rrg, FILE=sfile_rrg)
      open(kfile_asst, FILE=sfile_asst)
      open(kfile_assg, FILE=sfile_assg)
      open(kfile_su_av, FILE=sfile_su_av)
      open(kfile_zw, FILE=sfile_zw)
      open(kfile_wsnew, FILE=sfile_wsnew)
      open(kfile_spgfcf, FILE=sfile_spgfcf)
      open(kfile_infx, FILE=sfile_infx)
      open(kfile_etmt, FILE=sfile_etmt)
      open(kfile_etmg, FILE=sfile_etmg)
      open(kfile_su1, FILE=sfile_su1)
      open(kfile_topt, FILE=sfile_topt)

      !loop for n random samples, needs parallelization
      do i_loop=1, i_iter

write(*,*) "i_loop", i_loop

         !generate random number between 0  and 1
         call random_number(r)

         !make a random parameterset
         paramset=r*(parmax-parmin)+parmin

         !run the model with the random set

         call transpmodel(paramset, vom_npar, obj, vom_command)




         write(kfile_random_output, *) obj

      end do

      close( kfile_random_output )

      close(kfile_vd_d)
      close(kfile_esoil)
      close(kfile_jmax25t)
      close(kfile_jmax25g)
      close(kfile_vegcov)
      close(kfile_resp)
      close(kfile_lambdat)
      close(kfile_lambdag)
      close(kfile_rrt)
      close(kfile_rrg)
      close(kfile_asst)
      close(kfile_assg)
      close(kfile_su_av)
      close(kfile_zw)
      close(kfile_wsnew)
      close(kfile_spgfcf)
      close(kfile_infx)
      close(kfile_etmt)
      close(kfile_etmg)
      close(kfile_su1)
      close(kfile_topt)









end subroutine
