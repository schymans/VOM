!***********************************************************************
!        Optimised Vegetation Optimality Model (VOM)
!        Core program to run optimisation (sce) and transpmodel
!-----------------------------------------------------------------------
!        Author: Stan Schymanski, CWR, University of Western Australia
!        05/05/2004
!
!        Now at: MPI for Biogeochemistry, Jena, Germany
!        30/07/2007
!   sschym@bgc-jena.mpg.de
!
!-----------------------------------------------------------------------
!
!  Copyright (C) 2008  Stan Schymanski
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!***********************************************************************

      program vom
      use vom_sce_mod
      implicit none

      REAL*8, ALLOCATABLE :: vom_invar(:)
      REAL*8              :: vom_objfun

      INTEGER             :: beststat
      INTEGER             :: iostat
      LOGICAL             :: exist
      LOGICAL             :: run_sce

      beststat = 0

!     * Parameter definitions

      call read_shufflepar()

      allocate(vom_invar(vom_npar))


      !optimize with sce
      if (vom_command .eq. 1 ) then

         call transpmodel_init_once(vom_command)

         !remove old outputfiles for safety 
         open( kfile_sceout, FILE=trim(adjustl(i_outputpath)) // &
             trim(adjustl(sfile_sceout)), IOSTAT=iostat)
            if (iostat .eq. 0) close(kfile_sceout, status='delete')
         open( kfile_progress, FILE=trim(adjustl(i_outputpath)) // &
             trim(adjustl(sfile_progress)), STATUS='old', IOSTAT=iostat)
            if (iostat .eq. 0) close(kfile_progress, status='delete')
         open( kfile_lastloop, FILE=trim(adjustl(i_outputpath)) // &
             trim(adjustl(sfile_lastloop)), STATUS='old', IOSTAT=iostat)
            if (iostat .eq. 0) close(kfile_lastloop, status='delete')
         open( kfile_lastbest, FILE=trim(adjustl(i_outputpath)) // &
             trim(adjustl(sfile_lastbest)), STATUS='old', IOSTAT=iostat)
            if (iostat .eq. 0) close(kfile_lastbest, status='delete')
         open( kfile_bestpars, FILE=trim(adjustl(i_outputpath)) // &
             trim(adjustl(sfile_bestpars)), STATUS='old', IOSTAT=iostat)
            if (iostat .eq. 0) close(kfile_bestpars, status='delete')
         open( kfile_beststat, FILE=trim(adjustl(i_outputpath)) // &
             trim(adjustl(sfile_beststat)), STATUS='old', IOSTAT=iostat)
            if (iostat .eq. 0) close(kfile_beststat, status='delete')



         !run with initial seed
         write(*,*) "Initializing..."
         call sce_main()

         !run again untill sce_status.txt is written
         run_sce = .TRUE.
         do while (run_sce)
            write(*,*) "Start looping"
            call sce_main()
            open(kfile_beststat, FILE=trim(adjustl(i_outputpath)) // &
             trim(adjustl(sfile_beststat)), STATUS='old', IOSTAT=iostat)
            if (iostat .eq. 0) then
               run_sce = .FALSE.
            end if
         end do

      endif

      !run normally
      if (vom_command .eq. 2 ) then

         write(*,*) "Start single calculation with parameters..."

         call transpmodel_init_once(vom_command)

         open(kfile_pars, FILE=trim(adjustl(i_outputpath)) // &
             trim(adjustl(sfile_pars)),&
              STATUS='old', IOSTAT=iostat)
              if (iostat .ne. 0) then
                write(0,*) "ERROR opening ", sfile_pars
                stop
              endif
              rewind(kfile_pars)
              read(kfile_pars,*) vom_invar(:)
            close(kfile_pars)


        call transpmodel(vom_invar, SIZE(vom_invar), vom_objfun, vom_command)

      endif

      !run normally, no output except for ncp
      if (vom_command .eq. 3 ) then

         write(*,*) "Start calculation of ncp with parameters..."

         open(kfile_pars, FILE=trim(adjustl(i_outputpath)) // &
             trim(adjustl(sfile_pars)), STATUS='old', IOSTAT=iostat)
              if (iostat .ne. 0) then
                write(0,*) "ERROR opening ", sfile_pars
                stop
              endif
              rewind(kfile_pars)
              read(kfile_pars,*) vom_invar(:)
            close(kfile_pars)

        call transpmodel_init_once(vom_command)
        call transpmodel(vom_invar, SIZE(vom_invar), vom_objfun, vom_command)

        write(*,*) "The best carbon profit was: ",vom_objfun
      endif


      if (vom_command .eq. 5 ) then
              write(*,*) "Random sampling of parameters ... "
         call random_samples()

      endif



      deallocate(vom_invar)

      write(*,*) "Program terminated"

      end
