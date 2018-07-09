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
      INTEGER                          :: iostat   !input-output status

      namelist /outputfiles/ vd_d_out, &
                             esoil_out, &
                             jmax25t_out, &
                             jmax25g_out, &
                             vegcov_out, &
                             resp_out, &
                             lambdat_out, &
                             lambdag_out, &
                             rrt_out, &
                             rrg_out, &
                             asst_out, &
                             assg_out, &
                             su_av_out, &
                             zw_out, &
                             wsnew_out, &
                             spgfcf_out, &
                             infx_out, &
                             etmt_out, &
                             etmg_out, &
                             su1_out, &
                             topt_out

      open(kfile_outputlist, FILE=sfile_outputlist, STATUS='old',          &
     &                     FORM='formatted', IOSTAT=iostat)
      if (iostat .eq. 0) then
        read(kfile_outputlist, outputfiles)
      endif
      close(kfile_outputlist)

      !read settings for optimization
      call read_shufflepar()

      !read max and min ranges of parameters
      call read_shufflevar ()


      !initialize random seed 
      call random_seed()
      call transpmodel_init_once(vom_command)

      !open file for output
      open(kfile_random_output, FILE=trim(adjustl(i_outputpath)) // trim(adjustl(sfile_random_output)))
      open(kfile_random_params, FILE=trim(adjustl(i_outputpath)) // trim(adjustl(sfile_random_params)))

      if( vd_d_out .eqv. .TRUE.) then
         open(kfile_vd_d , FILE=trim(adjustl(i_outputpath)) // trim(adjustl(sfile_vd_d)))
      end if

      if( esoil_out .eqv. .TRUE.) then
         open(kfile_esoil, FILE=trim(adjustl(i_outputpath)) // trim(adjustl(sfile_esoil)))
      end if

      if( jmax25t_out .eqv. .TRUE.) then
         open(kfile_jmax25t, FILE=trim(adjustl(i_outputpath)) // trim(adjustl(sfile_jmax25t)))
      end if
      if( jmax25g_out .eqv. .TRUE.) then
      open(kfile_jmax25g, FILE=trim(adjustl(i_outputpath)) // trim(adjustl(sfile_jmax25g)))
      end if
      if( vegcov_out .eqv. .TRUE.) then
      open(kfile_vegcov, FILE=trim(adjustl(i_outputpath)) // trim(adjustl(sfile_vegcov)))
      end if
      if( resp_out .eqv. .TRUE.) then
         open(kfile_resp, FILE=trim(adjustl(i_outputpath)) // trim(adjustl(sfile_resp)))
      end if
      if( lambdat_out .eqv. .TRUE.) then
         open(kfile_lambdat, FILE=trim(adjustl(i_outputpath)) // trim(adjustl(sfile_lambdat)))
      end if
      if( lambdag_out .eqv. .TRUE.) then
         open(kfile_lambdag, FILE=trim(adjustl(i_outputpath)) // trim(adjustl(sfile_lambdag)))
      end if
      if( rrt_out .eqv. .TRUE.) then
         open(kfile_rrt, FILE=trim(adjustl(i_outputpath)) // trim(adjustl(sfile_rrt)))
      end if
      if( rrg_out .eqv. .TRUE.) then
         open(kfile_rrg, FILE=trim(adjustl(i_outputpath)) // trim(adjustl(sfile_rrg)))
      end if
      if( asst_out .eqv. .TRUE.) then
         open(kfile_asst, FILE=trim(adjustl(i_outputpath)) // trim(adjustl(sfile_asst)))
      end if
      if( assg_out .eqv. .TRUE.) then
         open(kfile_assg, FILE=trim(adjustl(i_outputpath)) // trim(adjustl(sfile_assg)))
      end if
      if( su_av_out .eqv. .TRUE.) then
         open(kfile_su_av, FILE=trim(adjustl(i_outputpath)) // trim(adjustl(sfile_su_av)))
      end if
      if( zw_out .eqv. .TRUE.) then
         open(kfile_zw, FILE=trim(adjustl(i_outputpath)) // trim(adjustl(sfile_zw)))
      end if
      if( wsnew_out .eqv. .TRUE.) then
         open(kfile_wsnew, FILE=trim(adjustl(i_outputpath)) // trim(adjustl(sfile_wsnew)))
      end if
      if( spgfcf_out .eqv. .TRUE.) then
         open(kfile_spgfcf, FILE=trim(adjustl(i_outputpath)) // trim(adjustl(sfile_spgfcf)))
      end if
      if( infx_out .eqv. .TRUE.) then
         open(kfile_infx, FILE=trim(adjustl(i_outputpath)) // trim(adjustl(sfile_infx)))
      end if
      if( etmt_out .eqv. .TRUE.) then
          open(kfile_etmt, FILE=trim(adjustl(i_outputpath)) // trim(adjustl(sfile_etmt)))
      end if
      if( etmg_out .eqv. .TRUE.) then
         open(kfile_etmg, FILE=trim(adjustl(i_outputpath)) // trim(adjustl(sfile_etmg)))
      end if
      if( su1_out .eqv. .TRUE.) then
         open(kfile_su1, FILE=trim(adjustl(i_outputpath)) // trim(adjustl(sfile_su1)))
      end if
      if( topt_out .eqv. .TRUE.) then
         open(kfile_topt, FILE=trim(adjustl(i_outputpath)) // trim(adjustl(sfile_topt)))
      end if

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
         write(kfile_random_params, *) paramset

      end do

      close( kfile_random_output )
      close( kfile_random_params )

      if( vd_d_out .eqv. .TRUE.) then
         close(kfile_vd_d)
      end if
      if( esoil_out .eqv. .TRUE.) then
         close(kfile_esoil)
      end if
      if( jmax25t_out .eqv. .TRUE.) then
         close(kfile_jmax25t)
      end if
      if( jmax25g_out .eqv. .TRUE.) then
         close(kfile_jmax25g)
      end if
      if( vegcov_out .eqv. .TRUE.) then
         close(kfile_vegcov)
      end if
      if( resp_out .eqv. .TRUE.) then
         close(kfile_resp)
      end if
      if( lambdat_out .eqv. .TRUE.) then
         close(kfile_lambdat)
      end if
      if( lambdag_out .eqv. .TRUE.) then
         close(kfile_lambdag)
      end if
      if( rrt_out .eqv. .TRUE.) then
         close(kfile_rrt)
      end if
      if( rrg_out .eqv. .TRUE.) then
         close(kfile_rrg)
      end if
      if( asst_out .eqv. .TRUE.) then
         close(kfile_asst)
      end if
      if( assg_out .eqv. .TRUE.) then
         close(kfile_assg)
      end if
      if( su_av_out .eqv. .TRUE.) then
         close(kfile_su_av)
      end if
      if( zw_out .eqv. .TRUE.) then
         close(kfile_zw)
      end if
      if( wsnew_out .eqv. .TRUE.) then
         close(kfile_wsnew)
      end if
      if( spgfcf_out .eqv. .TRUE.) then
         close(kfile_spgfcf)
      end if
      if( infx_out .eqv. .TRUE.) then
         close(kfile_infx)
      end if
      if( etmt_out .eqv. .TRUE.) then
         close(kfile_etmt)
      end if
      if( etmg_out .eqv. .TRUE.) then
         close(kfile_etmg)
      end if
      if( su1_out .eqv. .TRUE.) then
         close(kfile_su1)
      end if
      if( topt_out .eqv. .TRUE.) then
         close(kfile_topt)
      end if

end subroutine

