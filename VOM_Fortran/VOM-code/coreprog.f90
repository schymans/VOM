      program vom
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

      implicit none

!      CHARACTER(60) :: outformat
      INTEGER       :: command
      INTEGER       :: iostatus, stat
      CHARACTER(3)  :: str
      REAL*8        :: invar(6)
      REAL*8        :: netass
      INTEGER       :: nrun
      INTEGER       :: npar

!-----------------------------------------------------------------------
! Parameter definitions for sce 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! for debug purposes:
! option1='-optimise'
!-----------------------------------------------------------------------

      open(1, file='shuffle.par', status='old')
      read(1,*) command
    
!     * now with fourth commmand (3 for compute ncp oonly with pars.txt)
      if (command .eq. 3) then
      close(1)
        open(3, file='pars.txt', status='old', iostat=stat)
      
        if (stat .eq. 0) then
          rewind(3)
          read(3,*) invar(:)
          close(3)
!         * values given, because they are needed for calling transpmodel subroutine 
          netass = 0.d0
          nrun = 1
      
          print *, "Pars.txt read. Start calculation of ncp with parameters..."
      
          call transpmodel(invar, size(invar), nrun, netass, command)

          print *, "Model run COMPLETE"
          write(*,'(" The carbon profit achieved is: ",e12.6)') netass
          print *, "Best ncp is saved in model_output.txt"
!         * programm ends here      
        else
          write(*,*) "ERROR: pars.txt missing."
          stop
!         * programm ends here
        endif
                
      else
        open(2, file='finalbest.txt', status='old', iostat=iostatus)

        if (iostatus .eq. 0 .or. command .eq. 2) then
          command = 2
          if (iostatus .ne. 0) then
            close(2)
            open(2, file='currentbest.txt')    ! reads input parameters from previous optimisation
          endif

!-----------------------------------------------------------------------
! model run with optimised parameters
!-----------------------------------------------------------------------

          nrun = 1
      
          print *,"Calculation of results with optimised parameters..."
          
          rewind(2)
          read(2,*) invar(:), netass
          close(2)

          write(*,'(" The best carbon profit was: ",e12.6)') netass

          call transpmodel(invar, size(invar), nrun, netass, command)

          write(*,'(/" Model run COMPLETE",/)')
          write(*,'(" The carbon profit achieved is: ",e12.6)') netass
          print *, "Hourly results are saved in resulthourly.txt"
          print *, "Daily results are saved in resultsdaily.txt"
          print *, "Yearly results are saved in yearly.txt"
          print *, "Soil results are saved in delyudaily.txt, rsurfdaily.txt, ruptkhourly.txt, suvechourly.txt"
        
    else
          npar = 0
            
          do
            read(1,*,iostat=iostatus) str
            if (iostatus .lt. 0) exit
            if (str .eq. 'var') npar = npar + 1
          enddo
            
          close(1)
      
          if (npar .ne. 6) then
            write(*,*) "ERROR: shuffle.par has to contain 6 parameters (var)"
      close(2)
            stop
          endif

          call sce()
      
        endif
    
      endif 
    
      print *,"Program terminated"

      end