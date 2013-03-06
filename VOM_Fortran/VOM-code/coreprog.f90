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

      INTEGER       :: command
      REAL*8        :: invar(6)
      REAL*8        :: netass

      INTEGER       :: npar
      CHARACTER(3)  :: str

      INTEGER  :: iostat
      INTEGER  :: nrun

!-----------------------------------------------------------------------
! for debug purposes:
! option1='-optimise'
!-----------------------------------------------------------------------

!     * Parameter definitions

      open(kfile_shufflepar, FILE=sfile_shufflepar, STATUS='old')
      read(kfile_shufflepar,*) command

!     * now with fourth commmand (3 for compute ncp oonly with pars.txt)

      if (command .eq. 3) then
        close(kfile_shufflepar)
        open(kfile_pars, FILE=sfile_pars, STATUS='old', IOSTAT=iostat)

        if (iostat .eq. 0) then
          rewind(kfile_pars)
          read(kfile_pars,*) invar(:)
          close(kfile_pars)

          netass = 0.d0
          nrun = 1

          write(*,*) "Start calculation of ncp with parameters..."

          call transpmodel(invar, size(invar), nrun, netass, command)

          write(*,*) "Model run COMPLETE"
          write(*,*) ' '
          write(*,'(" The carbon profit achieved is: ",E12.6)') netass
          write(*,*) "Best ncp is saved in model_output.txt"
        else
            write(*,*) "ERROR reading ", sfile_pars
            stop
        endif

      else
        open(kfile_finalbest, FILE=sfile_finalbest, STATUS='old', IOSTAT=iostat)

        if (iostat .eq. 0 .or. command .eq. 2) then
          command = 2
          if (iostat .ne. 0) then
            close(kfile_finalbest)
!           * reads input parameters from previous optimisation
            open(kfile_finalbest, FILE=sfile_currentbest)
          endif

!         * model run with optimised parameters

          nrun = 1

          write(*,*) "Calculation of results with optimised parameters..."

          rewind(kfile_finalbest)
          read(kfile_finalbest,*) invar(:), netass
          close(kfile_finalbest)

          write(*,'(" The best carbon profit was: ",E12.6)') netass

        call transpmodel(invar, size(invar), nrun, netass, command)

          write(*,*) 'Model run COMPLETE'
          write(*,*) ' '
          write(*,'(" The carbon profit achieved is: ",E12.6)') netass
          write(*,*) "Hourly results are saved in resulthourly.txt"
          write(*,*) "Daily results are saved in resultsdaily.txt"
          write(*,*) "Yearly results are saved in yearly.txt"
          write(*,*) "Soil results are saved in delyudaily.txt, rsurfdaily.txt, ruptkhourly.txt, suvechourly.txt"
        else

        npar = 0
        do
          read(kfile_shufflepar,*,IOSTAT=iostat) str
          if (iostat .lt. 0) exit
          if (str .eq. 'var') npar = npar + 1
        enddo

        close(kfile_shufflepar)
        close(kfile_finalbest)

        if (npar .ne. 6) then
          write(*,*) "ERROR: shuffle.par has to contain 6 parameters (var)"
          stop
        endif

          call sce()

        endif

      endif

      write(*,*) "Program terminated"

      end
