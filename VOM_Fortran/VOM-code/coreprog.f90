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
      use vom_file_mod
      implicit none

      INTEGER       :: vom_command
      REAL*8        :: vom_invar(6)
      REAL*8        :: vom_objfun

      INTEGER       :: npar
      CHARACTER(3)  :: str

      INTEGER :: iostat
      LOGICAL :: exist

!-----------------------------------------------------------------------
! for debug purposes:
! option1='-optimise'
!-----------------------------------------------------------------------

!     * Parameter definitions

      open(kfile_shufflepar, FILE=sfile_shufflepar, STATUS='old')
      read(kfile_shufflepar,*) vom_command

      inquire(FILE=sfile_finalbest, EXIST=exist)
      if (exist .and. vom_command .ne. 3) vom_command = 2

      if (vom_command .eq. 2 .or. vom_command .eq. 3) then
        close(kfile_shufflepar)

!       * single model run with given (optimised) parameters

        if (vom_command .eq. 3) then
          write(*,*) "Start calculation of ncp with parameters..."
        else
          write(*,*) "Calculation of results with optimised parameters..."
        endif

        if (vom_command .eq. 3) then
          open(kfile_pars, FILE=sfile_pars, STATUS='old', IOSTAT=iostat)
          if (iostat .ne. 0) then
            write(*,*) "ERROR reading ", sfile_pars
            stop
          endif
          read(kfile_pars,*) vom_invar(:)
          close(kfile_pars)
          vom_objfun = 0.d0
        else
          if (exist) then
            open(kfile_finalbest, FILE=sfile_finalbest,                &
     &                            STATUS='old', IOSTAT=iostat)
            if (iostat .ne. 0) then
              write(*,*) "ERROR reading ", sfile_finalbest
              stop
            endif
          else
!           * reads input parameters from previous optimisation
            open(kfile_finalbest, FILE=sfile_currentbest,              &
     &                            STATUS='old', IOSTAT=iostat)
            if (iostat .ne. 0) then
              write(*,*) "ERROR reading ", sfile_currentbest
              stop
            endif
          endif
          rewind(kfile_finalbest)
          read(kfile_finalbest,*) vom_invar(:), vom_objfun
          close(kfile_finalbest)
          write(*,'(" The best carbon profit was: ",E13.6)') vom_objfun
        endif

        call transpmodel(vom_invar, SIZE(vom_invar), vom_objfun, vom_command)

        if (vom_command .eq. 3) then
          write(*,*) "Model run COMPLETE"
          write(*,*) ' '
          write(*,'(" The carbon profit achieved is: ",E13.6)') vom_objfun
          write(*,*) "Best ncp is saved in model_output.txt"
        else
          write(*,*) 'Model run COMPLETE'
          write(*,*) ' '
          write(*,'(" The carbon profit achieved is: ",E13.6)') vom_objfun
          write(*,*) "Hourly results are saved in resulthourly.txt"
          write(*,*) "Daily results are saved in resultsdaily.txt"
          write(*,*) "Yearly results are saved in yearly.txt"
          write(*,*) "Soil results are saved in delyudaily.txt, rsurfdaily.txt, ruptkhourly.txt, suvechourly.txt"
        endif

      else
        npar = 0
        do
          read(kfile_shufflepar,*,IOSTAT=iostat) str
          if (iostat .lt. 0) exit
          if (str .eq. 'var') npar = npar + 1
        enddo

        close(kfile_shufflepar)

        if (npar .ne. 6) then
          write(*,*) "ERROR: shuffle.par has to contain 6 parameters (var)"
          stop
        endif

        call sce_main()

      endif

      write(*,*) "Program terminated"

      end
