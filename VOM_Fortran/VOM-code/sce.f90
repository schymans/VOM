!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  SHUFFLED COMPLEX EVOLUTION
!    Parameter optimisation algorithm based on a paper by Duan et al. (1993,
!    J. Opt. Theory and Appl., 76, 501--521). 
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  VERSION 2.1        ---  01 February 1999         
!  Written by:        Neil Viney, Centre for Water Research (CWR), The University of WA
!  Modified by Stan Schymanski, CWR, 05 April 2004 (to run with transpmodel)  
!  Extended by Stan Schymanski, SESE, 02 June 2006 to follow Muttil & Liong (2004,
!  Journal of Hydraulic Engineering-Asce 130(12):1202-1205) and Duan et al. (1994, 
!  Journal of Hydrology 158)
!        
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  This implementation MAXIMISES the objective function, which is calculated
!  by the model, not by the optimiser.        The optimiser transfers parameter values
!  to the model subroutine and receives the value of the objective function.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Copyright (C) 2008 Stan Schymanski
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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      module sce_mod

      INTEGER :: success = 0
      INTEGER :: npar
      INTEGER :: ncomp
      INTEGER :: ncompmin
      INTEGER :: ncomp2
      INTEGER :: nopt
      INTEGER :: mopt
      INTEGER :: sopt
      INTEGER :: qopt
      INTEGER :: alpha
      INTEGER :: nrun
      INTEGER :: nloop
      INTEGER :: nsincebest
      INTEGER :: patience

      INTEGER, ALLOCATABLE :: optid(:)

      REAL*8  :: ranscal
      REAL*8  :: focus
      REAL*8  :: worstbest
      REAL*8  :: bestobj
      REAL*8  :: bestincomp
      REAL*8  :: resolution

      REAL*8, ALLOCATABLE :: wgt(:)
      REAL*8, ALLOCATABLE :: cv_(:)
      REAL*8, ALLOCATABLE :: ranarr(:)

      REAL*8, ALLOCATABLE :: shufflevar(:,:)
      REAL*8, ALLOCATABLE :: parval(:)
      REAL*8, ALLOCATABLE :: parmin(:)
      REAL*8, ALLOCATABLE :: parmax(:)
      REAL*8, ALLOCATABLE :: ofvec(:)

      CHARACTER(9), ALLOCATABLE :: parname(:)
      CHARACTER(12)  :: evolution
      CHARACTER(60)  :: outformat, loopformat

!     * allocated variables for optsensitivity()
      REAL*8, ALLOCATABLE :: dataarray(:,:)
      REAL*8, ALLOCATABLE :: shufflevar2(:)

!     * allocated variables for cce()
      INTEGER, ALLOCATABLE :: parentsid(:)
      REAL*8,  ALLOCATABLE :: objfunsub(:)
      REAL*8,  ALLOCATABLE :: invarsub(:,:)
      LOGICAL, ALLOCATABLE :: selected(:)

!     * allocated variables for simplex()
      REAL*8, ALLOCATABLE :: centroid(:)
      REAL*8, ALLOCATABLE :: newpoint(:)

!     * file codes

      INTEGER :: kfile_sceout      = 701
      INTEGER :: kfile_shufflepar  = 702
      INTEGER :: kfile_progress    = 703
      INTEGER :: kfile_lastloop    = 704
      INTEGER :: kfile_currentbest = 705
      INTEGER :: kfile_bestpars    = 706
      INTEGER :: kfile_finalbest   = 707

!     * file names

      CHARACTER(80) :: sfile_sceout      = 'sce.out'
      CHARACTER(80) :: sfile_shufflepar  = 'shuffle.par'
      CHARACTER(80) :: sfile_progress    = 'progress.txt'
      CHARACTER(80) :: sfile_lastloop    = 'lastloop.txt'
      CHARACTER(80) :: sfile_currentbest = 'currentbest.txt'
      CHARACTER(80) :: sfile_bestpars    = 'bestpars.txt'
      CHARACTER(80) :: sfile_finalbest   = 'finalbest.txt'

      end module sce_mod

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine sce ()
      use sce_mod
      implicit none

      INTEGER :: i_, stat, first, numcv
      REAL*8  :: maxcv
      REAL*8, ALLOCATABLE :: sumvar(:)
      CHARACTER(24)  :: logdate
      CHARACTER(300) :: writeformat

!EXTERNAL compar
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initialize the random number generator
! (standard subroutine, based on the date and time)
      if (command .ne. 4) then
        CALL RANDOM_SEED()
      endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! WRITE SCREEN HEADER

      write(*,'(/"SHUFFLED COMPLEX EVOLUTION OPTIMISER")')
      call fdate(logdate)
      write(*,'(/"  Run time:   ",a)') logdate

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! INITIALIZATION

      call sce_init()

      allocate(sumvar(nopt))

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! BEGIN MODEL LOOP (EXECUTE s PASSES)
! CALL model SUBROUTINE        (OPEN run.log FOR CONSOLE OUTPUT);
! CALCULATE OBJECTIVE FUNCTION FOR DEFAULT PARAMETER VALUES
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! INSERTED BY STAN TO ALLOW CONTINUATION OF OPTIMSATION FROM PREVIOUSLY SAVED STEP

      open(kfile_lastloop, file=sfile_lastloop(1:len_trim(sfile_lastloop)), status='old', iostat=stat)
      if (stat .eq. 0) then
        read(kfile_lastloop,*) ncomp2
        read(kfile_lastloop,*) nloop
        read(kfile_lastloop,*) nrun
        read(kfile_lastloop,*) nsincebest
        read(kfile_lastloop,loopformat) ofvec(:)
        do i_ = 1, npar
          read(kfile_lastloop, loopformat) shufflevar(i_,:)
        enddo
        close(kfile_lastloop)
        bestobj = ofvec(1)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! OPEN FILES FOR STORING OBJECTIVE FUNCTION AND PARAMETER VALUES

        open(kfile_sceout, file=sfile_sceout(1:len_trim(sfile_sceout)), status='old', position='append')
        open(kfile_bestpars, file=sfile_bestpars(1:len_trim(sfile_bestpars)), status='old', position='append')
        open(kfile_progress, file=sfile_progress(1:len_trim(sfile_progress)), status='old', position='append')
        write(kfile_progress,'(/"  NEW Run time:   ",a)') logdate
      else

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! OPEN AND EMPTY FILE FOR STORING OBJECTIVE FUNCTION AND PARAMETER VALUES

        open(kfile_sceout, file=sfile_sceout(1:len_trim(sfile_sceout)), status='replace')
        open(kfile_bestpars, file=sfile_bestpars(1:len_trim(sfile_bestpars)), status='replace')
        open(kfile_progress, file=sfile_progress(1:len_trim(sfile_progress)), status='replace')

!       * write file header
        write(kfile_progress,'(/"SHUFFLED COMPLEX EVOLUTION OPTIMISER")')
        write(kfile_progress,'(/"  Run time:   ",a)') logdate

        call initialseed

        close(kfile_sceout)
        close(kfile_bestpars)
        close(kfile_progress)
        return
      endif

! END OF INSERTION
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! BEGIN SCE LOOP

      do while (nrun .lt. 20000 .and. nloop .lt. 500)
        nloop = nloop + 1

!       * Saving the best OF of the worst complex in worstbest for
!       * assessment of gene pool mixing

        first = 1 + (ncomp2 - 1) * mopt
        worstbest = ofvec(first)
        call sortcomp(shufflevar(:,:), shape(shufflevar(:,:)), ofvec(:), size(ofvec(:)))  ! [SORT ENTIRE ARRAYS]

!       * WRITE BEST_PARAMETERS FILE FOR PREVIOUS LOOP

        writeformat = '("Finished ",i4," main loops'
        writeformat(29:66) = ' --- best objective function =",e12.6)'
        write(*,writeformat) nloop, ofvec(1)
        write(kfile_progress,writeformat) nloop, ofvec(1)
        writeformat = '(/"No improvement in OF for",i5," loops")'
        write(*,writeformat) nsincebest
        write(kfile_progress,writeformat) nsincebest
        nsincebest = nsincebest + 1
    
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ASSESS CONVERGENCE
! CHANGED BY STAN TO CALCULATE THE FLUCTUATION RANGE RELATIVELY TO
! THE FEASIBLE RANGE, INSTEAD OF CV:

        numcv = ncomp2 * mopt
        sumvar(:) = sum(shufflevar(optid(:), 1:numcv),2) / numcv  ! mean parameter values
        do i_ = 1, nopt
          cv_(i_) = maxval(abs((shufflevar(optid(i_), 1:numcv) - sumvar(i_))  &
     &          / (parmin(optid(i_)) - parmax(optid(i_)))) * 100.d0)  ! distance from mean in % of feasible range
        enddo
        maxcv = maxval(cv_(:))                ! maximum distance
        writeformat = '("Greatest parameter range: ",f5.2,"%'
        writeformat(38:67) = ' for optimised parameter ",a9)'
        write(*,writeformat) maxcv, parname(optid(maxloc(cv_(:))))
        write(kfile_progress,writeformat) maxcv, parname(optid(maxloc(cv_(:))))
        if (maxcv .ge. resolution) then
          if (nsincebest .le. patience) then
            call writepars
            call run_cce()
            return
          else
            writeformat = '(/"No improvement in OF for",i5," loops"/" '
            writeformat(44:65) = ' About to give up...")'
            write(*,writeformat) nsincebest
            write(kfile_progress,writeformat) nsincebest
            call ck_success()
            if (success .eq. 1) return
          endif
        else
          writeformat = '(/"First Convergence criterion satisfied..."/"'
          writeformat(47:92) = '  parameter ranges are all less than 0.1 %"//)'
          write(*,writeformat)
          write(kfile_progress,writeformat)

!         * STAN'S MODIFICATION TO ASSESS SENSITIVITY OF OBJECTIVE
!         * FUNCTION TO EACH PARAMETER:

          call ck_success()
          if (success .eq. 1) return
        endif
      enddo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! END OF BLOCK THAT CAN BE PARALLELISED!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! TERMINATE PROGRAM

      writeformat = '(/"FAILURE TO CONVERGE..."/"  Number of runs has '
      writeformat(50:96) = 'reached 20000." //"  Program terminated."/,a,a)'
      write(*,writeformat) char(7), char(7)
      write(kfile_progress,writeformat) char(7), char(7)
      call sortcomp(shufflevar(:,:), shape(shufflevar(:,:)), ofvec(:), size(ofvec(:)))  ! [SORT ENTIRE ARRAYS]
      call writepars                              ! PROGRAM STOP
      open(kfile_finalbest, file=sfile_finalbest(1:len_trim(sfile_finalbest)))
      write(kfile_finalbest,outformat) shufflevar(:,1), bestobj
      close(kfile_finalbest)
      close(kfile_sceout)
      close(kfile_bestpars)
      close(kfile_progress)

      deallocate(sumvar)

      return
      end subroutine sce

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine sce_init ()
      use sce_mod
      implicit none

      INTEGER :: i_, j_, ios, command
      INTEGER, ALLOCATABLE :: paropt(:)
      CHARACTER(3)   :: str
      CHARACTER(60)  :: informat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! FIND AND OPEN PARAMETER FILE, shuffle.par
! ESTABLISH NUMBER OF PARAMETERS npar

      open(kfile_shufflepar, file=sfile_shufflepar(1:len_trim(sfile_shufflepar)), status='old')
      read(kfile_shufflepar,*)
      read(kfile_shufflepar,*)
      npar = 0
      do
        read(kfile_shufflepar,*,iostat=ios) str
        if (ios .lt. 0) exit
        if (str .eq. 'var') npar = npar + 1
      enddo
      rewind(kfile_shufflepar)
      allocate(parname(npar))
      allocate(parval(npar))
      allocate(parmin(npar))
      allocate(parmax(npar))
      allocate(paropt(npar))

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  LOAD INITIAL PARAMETER VALUES AND PARAMETER RANGES

      read(kfile_shufflepar,'(i1)') command
      read(kfile_shufflepar,*) ncomp
      read(kfile_shufflepar,*) ncompmin
      read(kfile_shufflepar,*) resolution
      read(kfile_shufflepar,*) patience
      read(kfile_shufflepar,*) alpha
      read(kfile_shufflepar,*) focus
      read(kfile_shufflepar,'(a60)') informat
      do i_ = 1, npar
        read(kfile_shufflepar,informat) parname(i_), parval(i_), parmin(i_), &
     &                                  parmax(i_), paropt(i_)
      enddo
      close(kfile_shufflepar)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  CALCULATE NUMBER OF OPTIMISABLE PARAMETERS

      nopt = sum(paropt(:))
      mopt = 2 * nopt + 1                         ! SCE VARIABLE m
      ncomp = MAX(ncomp, Ceiling(1.d0 + 2.d0 ** nopt / (1.d0 + 2.d0 * nopt)))  ! number of complexes after Muttil(2004) or from shuffle.par
      sopt = mopt * ncomp                         ! SCE VARIABLE s
      qopt = nopt + 1                             ! CCE VARIABLE q
      ncomp2 = ncomp
      allocate(optid(nopt))
      allocate(shufflevar(npar,sopt))
      allocate(ofvec(sopt))
      allocate(wgt(mopt))
      allocate(cv_(nopt))
      allocate(ranarr(nopt))
      allocate(dataarray(nopt*8+1,nopt+1))
      allocate(shufflevar2(npar))
      allocate(parentsid(qopt))
      allocate(objfunsub(qopt))
      allocate(invarsub(npar,qopt))
      allocate(selected(mopt))
      allocate(centroid(npar))
      allocate(newpoint(npar))

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! INITIALISE OUTPUT FORMAT STRING
! INTRODUCED BY STAN TO HAVE ONE COLUMN PER PARAMETER AND ONE FOR OF

      write(str,'(i2)') npar + 1                  ! internal write to convert from integer to string
      outformat = '('//str//'e34.25)'             ! includes a column for each parameter and a column for the value of OF

!     * ADDED BY STAN TO WRITE shufflevar AND ofvec OF LAST LOOP TO FILE
      write(str,'(i3)') sopt                      ! internal write to convert from number to string
      loopformat = '('//str//'e34.25)'            ! includes a column for each set 

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! INITIALISE optid: THE INDEX OF THOSE PARAMETERS THAT ARE OPTIMISABLE

      j_ = 0
      do i_ = 1, npar
        if (paropt(i_) .gt. 0) then
          j_ = j_ + 1
          optid(j_) = i_
        endif
      enddo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! INITIALISE optid: THE INDEX OF THOSE PARAMETERS THAT ARE OPTIMISABLE
! ASSIGN PROBABILITY WEIGHTS

      do i_ = 1, mopt
        wgt(i_) = 2.d0 * (mopt + 1.d0 - i_) / mopt / (mopt + 1.d0)
      enddo

      return
      end subroutine sce_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine initialseed ()
      use sce_mod
      implicit none

      INTEGER :: i_, j_, k_, worstcount                             ! worstcount for counting number of negative objective functions
      INTEGER, ALLOCATABLE :: posarray(:,:)
      REAL*8,  ALLOCATABLE :: initpop(:,:)

      allocate(posarray(2**nopt,nopt))
      allocate(initpop(nopt,5))

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PRINT DIMENSION INFORMATION

      write(*,'("  Number of model parameters:",10x,i3)') npar
      write(*,'("  Number of optimisable parameters:",i7)') nopt
      write(*,'("  Maximum number of complexes:",9x,i3)') ncomp
      write(*,'("  Minimum number of runs per complex:",i5,//)') mopt
      write(kfile_progress,'("  Number of model parameters:",10x,i3)') npar
      write(kfile_progress,'("  Number of optimisable parameters:",i7)') nopt
      write(kfile_progress,'("  Maximum number of complexes:",9x,i3)') ncomp
      write(kfile_progress,'("  Minimum number of runs per complex:",i5,//)') mopt

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! CALCULATING OF USING INITAL GUESS IN SHUFFLE.PAR

      nsincebest = 0
      evolution = 'seed'
      ofvec(:) = -9999.9d0
      shufflevar(:,1) = parval(:)
      nrun = 0
	  worstcount = 0
      call runmodel(shufflevar(:,1), ofvec(1))
      
	  if (ofvec(1) .le. 0) then
	    worstcount = worstcount + 1
	  endif
	  
	  bestobj = ofvec(1)
      bestincomp = bestobj
      open(kfile_currentbest, file=sfile_currentbest(1:len_trim(sfile_currentbest)))
      write(kfile_currentbest,outformat) shufflevar(:,1), bestobj
      close(kfile_currentbest)
      write(*,'("Systematic seed of",i4," parameters for ",i2," complexes. Initial OF= ",e12.6)') nopt, ncomp, ofvec(1)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! GENERATE A SYSTEMATIC ARRAY OF INITIAL PARAMETER VALUES FOLLOWING
! Muttil & Liong (2004, Journal of Hydraulic Engineering-Asce 130(12):1202-1205
! (INSERTED BY STAN)
!        
! NONAXIAL POINTS:

      posarray(1,1) = 1
      posarray(2,1) = 2
      do j_ = 1, nopt
        k_ = optid(j_)
        initpop(j_,1) = 0.125d0 * parmax(k_) + 0.875d0 * parmin(k_)  ! each position j contains the intial perturbation of an optimised parameter
        initpop(j_,2) = 0.125d0 * parmin(k_) + 0.875d0 * parmax(k_)
        initpop(j_,3) = 0.5d0 * (parmin(k_) + parmax(k_))
        initpop(j_,4) = 0.25d0 * parmax(k_) + 0.75d0 * parmin(k_)
        initpop(j_,5) = 0.25d0 * parmin(k_) + 0.75d0 * parmax(k_)
        posarray(2**(j_-1)+1:2**j_, 1:j_-1) = posarray(1:2**(j_-1), 1:j_-1)
        posarray(1:2**(j_-1), j_) = 1
        posarray(2**(j_-1)+1:2**j_, j_) = 2
      enddo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! THIS BLOCK CAN BE PARALLELISED!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      do i_ = 2,sopt
        shufflevar(:,i_) = parval(:)              ! TO SET NON-OPTIMISING PARAMETERS
        if (i_ .le. 4) then
          do j_ = 1, nopt
            k_ = optid(j_)
            shufflevar(k_,i_) = initpop(j_,i_+1)
          enddo
          call runmodel(shufflevar(:,i_), ofvec(i_))
		  
		  if (ofvec(i_) .le. 0) then
	        worstcount = worstcount + 1
	      endif
		  
        elseif (i_ .le. Size(posarray(:,:),1)) then
          do j_ = 1, nopt
            k_ = optid(j_)
            shufflevar(k_,i_) = initpop(j_, posarray(i_-4, j_))
          enddo
          call runmodel(shufflevar(:,i_), ofvec(i_))
		  
		  if (ofvec(i_) .le. 0) then
	        worstcount = worstcount + 1
	      endif
		  
        else

!         * IF MORE POINTS ARE NEEDED, GENERATE RANDOM POINTS

          evolution = 'mutation'
          do while (ofvec(i_) .le. 0.d0)          ! first loop must generate feasible values to start with
            call random_number(ranarr(:))         ! RANDOM ARRAY OF SIZE 1:nopt
            do j_ = 1, nopt
              k_ = optid(j_)

!             * STAN'S MODIFICATION TO GET COMPLETELY RANDOM SEED:

              shufflevar(k_,i_) = parmin(k_) + focus * (parmax(k_) - parmin(k_)) * ranarr(j_)
            enddo
            shufflevar(optid(:),i_) = merge(shufflevar(optid(:),i_), parmin(optid(:)),      &
     &                             shufflevar(optid(:),i_) .gt. parmin(optid(:)))
            shufflevar(optid(:),i_) = merge(shufflevar(optid(:),i_), parmax(optid(:)),      &
     &                             shufflevar(optid(:),i_) .lt. parmax(optid(:)))
            call runmodel(shufflevar(:,i_), ofvec(i_))
			
			if (ofvec(i_) .le. 0) then
	          worstcount = worstcount + 1
	        endif
			
			if (nrun .gt. 100 .and. nrun .eq. worstcount) then       ! program stops after 100 runs without positive objective function
			  success = 2
			  call writepars()
			  exit
		    endif
						
          enddo
        endif
        if (ofvec(i_) .gt. bestobj) then
          bestobj = ofvec(i_)
          open(kfile_currentbest, file=sfile_currentbest(1:len_trim(sfile_currentbest)))
          write(kfile_currentbest,outformat) shufflevar(:,i_), bestobj
          close(kfile_currentbest)
        endif
		
		if (success .eq. 2) exit
		
      enddo
      nloop = -1                                  ! FIRST LOOP IS LOOP ZERO
      call writeloop

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! END OF BLOCK THAT CAN BE PARALLELISED!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      deallocate(posarray)
      deallocate(initpop)

      return
      end subroutine initialseed

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine ck_success ()
      use sce_mod
      implicit none

      call optsensitivity()
      if (success .eq. 1) then
        call sortcomp(shufflevar(:,:), shape(shufflevar(:,:)), ofvec(:), size(ofvec(:)))  ! [SORT ENTIRE ARRAYS]
        call writepars
        write(*,'(/"Optimisation completed successfully."/)')
        write(kfile_progress,'(/"Optimisation completed successfully."/)')
        print *,char(7)
        print *,char(7)
        print *,char(7)
        close(kfile_sceout)
        close(kfile_bestpars)
        close(kfile_progress)
      else
        call sortcomp(shufflevar(:,:), shape(shufflevar(:,:)), ofvec(:), size(ofvec(:)))  ! [SORT ENTIRE ARRAYS]
        call writepars
        nsincebest = 0
        write(kfile_bestpars,outformat) shufflevar(:,1), bestobj
      endif

      return
      end subroutine ck_success

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ADDED BY STAN TO ASSESS SENSITIVITY OF OBJECTIVE FUNCTION TO EACH OPTIMISED
! PARAMETER:
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine optsensitivity ()
      use sce_mod
      implicit none

      INTEGER :: i_, j_, failed10, pos
      REAL*8  :: distmin, distmax, oldpar, newpar, ofvec2, ofchange, parchange
      CHARACTER(300) :: writeformat

      shufflevar2(:) = shufflevar(:,1)
      ofvec2 = ofvec(1)
      failed10 = 0

      print *,"SENSITIVITY ANALYSIS"
      print *,"(changes in % of feasible range)"
      print *
      write(kfile_progress,'("SENSITIVITY ANALYSIS"/"changes in % of feasible range)" )')
      dataarray(1,1:nopt) = shufflevar2(optid(:))
      dataarray(1,nopt+1) = ofvec2
      evolution = 'test'
      pos = 0
      do i_ = 1, nopt
        oldpar = shufflevar2(optid(i_))

!       * TO MAKE SURE THAT PERTURBATIONS DO NOT EXCEED THE FEASIBLE RANGE
!       * AND THE PARAMETER VALUES THEMSELVES

        distmin = oldpar - parmin(optid(i_))
        distmax = parmax(optid(i_)) - oldpar
        writeformat = '(/"change of var: ",a9, "(",e9.3,")")'
        write(*,writeformat) parname(optid(i_)), oldpar
        write(kfile_progress,writeformat) parname(optid(i_)), oldpar
        do j_ = 0, 3
          newpar = oldpar - distmin * 10.d0 ** (-j_)
          shufflevar2(optid(i_)) = newpar
          nrun = nrun + 1
          call transpmodel(shufflevar2(:), size(shufflevar2(:)), nrun, ofvec2, 1)
          if (ofvec2 .gt. bestobj) then
            bestobj = ofvec2
            open(kfile_currentbest, file=sfile_currentbest(1:len_trim(sfile_currentbest)))
            write(kfile_currentbest,outformat) shufflevar2(:), bestobj
            close(kfile_currentbest)
          endif
          write(kfile_sceout,outformat) shufflevar2(optid(:)), ofvec2
          dataarray((i_-1)*8+j_+2,1:nopt) = shufflevar2(optid(:))
          dataarray((i_-1)*8+j_+2,nopt+1) = ofvec2
          ofchange = (ofvec2 - dataarray(1,nopt + 1))               &
     &             / abs(dataarray(1,nopt + 1)) * 100.d0
          parchange = (newpar - oldpar) / (parmax(optid(i_))            &
     &              - parmin(optid(i_))) * 100.d0  ! the change of the parameter in % of feasible range
          writeformat = '(f6.3,"% (",f14.6,")",": change of OF by ",'
          writeformat(44:66) = 'e9.3,"%"," (",e9.3,")")'
          write(*,writeformat) parchange, newpar, ofchange, ofvec2
          write(kfile_progress,writeformat) parchange, newpar, ofchange, ofvec2
          if (ofchange .gt. 1.d-10) then
            shufflevar(:,sopt-pos) = shufflevar2(:)
            ofvec(sopt-pos) = ofvec2
            pos = pos + 1
            if (abs(parchange) .gt. resolution) then
              failed10 = failed10 + 1
            endif
          endif
        enddo

        do j_ = 3, 0, -1
          newpar = oldpar + distmax * 10.d0 ** (-j_)
          shufflevar2(optid(i_)) = newpar
          nrun = nrun + 1
          call transpmodel(shufflevar2(:), size(shufflevar2(:)), nrun, ofvec2, 1)
          if (ofvec2 .gt. bestobj) then
            bestobj = ofvec2
            open(kfile_currentbest, file=sfile_currentbest(1:len_trim(sfile_currentbest)))
            write(kfile_currentbest,outformat) shufflevar2(:), bestobj
            close(kfile_currentbest)
          endif
          write(kfile_sceout,outformat) shufflevar2(optid(:)), ofvec2
          dataarray((i_-1)*8+j_+6,1:nopt) = shufflevar2(optid(:))
          dataarray((i_-1)*8+j_+6,nopt+1) = ofvec2
          ofchange = (ofvec2 - dataarray(1,nopt + 1))               &
     &             / abs(dataarray(1,nopt + 1)) * 100.d0
          parchange = (newpar - oldpar) / (parmax(optid(i_))            &
     &              - parmin(optid(i_))) * 100.d0  ! the change of the parameter in % of feasible range
          write(*,writeformat) parchange, newpar, ofchange, ofvec2
          write(kfile_progress,writeformat) parchange, newpar, ofchange, ofvec2
          if (ofchange .gt. 1.0d-10) then
            shufflevar(:,sopt-pos) = shufflevar2(:)
            ofvec(sopt-pos) = ofvec2
            pos = pos + 1
            if (abs(parchange) .gt. resolution) then
              failed10 = failed10 + 1
            endif
          endif
        enddo
        shufflevar2(optid(i_)) = oldpar
      enddo

!     * ASSESS SECOND CONVERGENCE CRITERIUM: NO PARAMETER CHANGE BY MORE THAN
!     * 10% LEADS TO AN INCREASE IN OBJECTIVE FUNCTION

      if (failed10 .gt. 0) then
        writeformat = '(i2," parameter(s) more than ",f6.3,"% out of '
        writeformat(47:84) = 'optimum."/"Optimisation continued...")'
        write(*,writeformat) failed10, resolution
        write(kfile_progress,writeformat) failed10, resolution
      else
        writeformat = '(//"Second convergence criterion satisfied:"/"'
        writeformat(47:85) = 'no parameter shift by more than ",f6.3,'
        writeformat(86:118) = '"% of max distance leads to"/"an '
        writeformat(119:155) = 'increase in the objective function.")'
        write(*,writeformat) resolution
        write(kfile_progress,writeformat) resolution
        success = 1
      endif

!     * IF CRITERIUM NOT SATISFIED, APPEND NEWLY CREATED DATA TO THE END OF
!     * shufflevar AND ofvec ARRAYS AND RETURN TO MAIN LOOP TO RE-RUN SCE-LOOP
!     * --> deleted

      return
      end subroutine optsensitivity

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine run_cce ()
      use sce_mod
      implicit none

      INTEGER :: m_, first
      CHARACTER(300) :: writeformat

!     * PARTITION THE sopt POINTS INTO ncomp2 COMPLEXES
!     * EXECUTE CCE ALGORITHM

      if (ncomp2 .gt. 1) then
        if (nloop .gt. 0) then
          if (ncomp2 .gt. ncompmin) then
            ncomp2 = ncomp2 - 1                   ! REDUCE NUMBER OF COMPLEXES AS nloop INCREASES
          else
            if (worstbest .le. ofvec(1+(ncomp2-1)*mopt)) then
              writeformat = '("No gene pool mixing ... '
              writeformat(27:64) = 'reducing number of complexes by one.")'
              write(*,writeformat)
              write(kfile_progress,writeformat)
              ncomp2 = ncomp2 - 1
            endif
          endif
        endif
      endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! THIS BLOCK CAN BE PARALLELISED!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      do m_ = 1, ncomp2
        first = 1 + (m_ - 1) * mopt
        writeformat = '("Start of loop",i4,", complex",i2,'
        writeformat(36:55) = '": best OF =",e12.6)'
        write(*,writeformat) nloop + 1, m_, ofvec(first)
        write(kfile_progress,writeformat) nloop + 1, m_ ,ofvec(first)
        if (m_ .eq. 1) then
          bestincomp = -9999.d0                   ! SET LESS THAN bestobj
        else
          bestincomp = ofvec(first)
        endif
        call cce(ofvec(first:m_*mopt), shufflevar(:,first:m_*mopt))
        writeformat(3:7) = '  End'
        write(*,writeformat) nloop + 1, m_, ofvec(first)
        write(kfile_progress,writeformat) nloop + 1, m_, ofvec(first)
      enddo
      worstbest = ofvec(first)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! WRITE shufflevar AND ofvec OF LAST LOOP TO FILE

      call writeloop
      close(kfile_sceout)
      close(kfile_bestpars)
      close(kfile_progress)

      return
      end subroutine run_cce

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine cce (objfun, invar)
      use sce_mod
      implicit none

!     * Declarations
      REAL*8, DIMENSION(mopt),      INTENT(inout) :: objfun
      REAL*8, DIMENSION(npar,mopt), INTENT(inout) :: invar

!     * Definitions
      INTEGER :: i_, j_, nsel, rannum, l

!     * SELECT PARENTS

      do l = 1, mopt
        selected(:) = .false.
        nsel = 0
        do while (nsel .ne. qopt)
          call random_number(ranscal)             ! SCALAR RANDOM NUMBER
          rannum = ceiling((2.d0 * mopt + 1.d0 - sqrt(4.d0 * mopt      &
     &           * (mopt + 1.d0) * (1.d0 - ranscal) + 1.d0)) * 0.5d0)

!         * NOTE: A SIMPLER ALTERNATIVE TO THE ABOVE LINE IS (19.03.2004):

          if (rannum .ge. 1 .and. rannum .le. mopt) then
            if (.not. selected(rannum)) then
              selected(rannum) = .true.
              nsel = nsel + 1
            endif
          endif
        enddo
        nsel = 0
        i_ = 0
        do while (nsel .ne. qopt)
          i_ = i_ + 1
          if (selected(i_)) then
            nsel = nsel + 1
            parentsid(nsel) = i_
          endif
        enddo

!       * GENERATE OFFSPRING AND SORT THE RESULTING COMPLEX

        objfunsub(:) = objfun(parentsid(:))
        invarsub(:,:) = invar(:, parentsid(:))
        do j_ = 1, alpha
          call simplex(invarsub(:,:), objfunsub(:))
        enddo
        objfun(parentsid(:)) = objfunsub(:)
        invar(:, parentsid(:)) = invarsub(:,:)
        call sortcomp(invar, shape(invar), objfun, size(objfun))
      enddo

      return
      end subroutine cce

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine simplex (invar, objfun)
      use sce_mod
      implicit none

!     * Declarations
      REAL*8, DIMENSION(npar,qopt), INTENT(inout) :: invar
      REAL*8, DIMENSION(qopt),      INTENT(inout) :: objfun

!     * Definitions
      INTEGER :: i_, j_
      REAL*8  :: newobjfun, minj, rangej

!     * REFLECTION STEP

      evolution = 'reflection'
      newpoint(:) = invar(:,qopt)
      centroid(optid(:)) = 1.d0 / (qopt - 1) * sum(invar(optid(:), 1:qopt - 1), 2)
      newpoint(optid(:)) = 2.d0 * centroid(optid(:)) - invar(optid(:),qopt)
      if (minval(newpoint(:) - parmin(:)) .lt. 0.d0 .or.                     &
     &    maxval(newpoint(:) - parmax(:)) .gt. 0.d0) then

!       * MUTATION STEP
!       * NB: MUTATION IS BASED ON THE SMALLEST HYPERCUBE OF THE SUBCOMPLEX,
!       * NOT THE SMALLEST HYPERCUBE OF THE COMPLEX AS SUGGESTED BY DUAN ET AL.

        do i_ = 1, nopt
          call random_number(ranscal)             ! SCALAR RANDOM NUMBER
          j_ = optid(i_)
          minj = parmin(j_)
          rangej = parmax(j_) - minj
          newpoint(j_) = minj + rangej * ranscal   ! UNIFORMLY DISTRIBUTED SAMPLE
        enddo
        evolution = 'mutation'

!       * CALCULATE OBJECTIVE FUNCTION

      endif
      call runmodel(newpoint(:), newobjfun)
      if (newobjfun .gt. objfun(qopt)) then
        do j_ = 1, qopt
          if (newobjfun .gt. objfun(j_)) then     ! SORT objfun HERE
            objfun(j_+1:qopt) = objfun(j_:qopt-1) ! IN CASE alpha > 1
            invar(:,j_+1:qopt) = invar(:,j_:qopt-1)
            objfun(j_) = newobjfun
            invar(:,j_) = newpoint(:)
            exit
          endif
        enddo
      else

!       * CONTRACTION STEP

        newpoint(optid(:)) = (centroid(optid(:)) + invar(optid(:),qopt)) * 0.5d0
        evolution = 'contraction'
        call runmodel(newpoint(:), newobjfun)
        if (newobjfun .gt. objfun(qopt)) then
          do j_ = 1, qopt
            if (newobjfun .gt. objfun(j_)) then
              objfun(j_+1:qopt) = objfun(j_:qopt-1)
              invar(:,j_+1:qopt) = invar(:,j_:qopt-1)
              objfun(j_) = newobjfun
              invar(:,j_) = newpoint(:)
              exit
            endif
          enddo
        else

!         * ANOTHER MUTATION STEP

          do i_ = 1, nopt
            call random_number(ranarr(:))         ! RANDOM ARRAY OF SIZE 1:nopt
            j_ = optid(i_)
            minj = minval(invar(j_,:))
            rangej = maxval(invar(j_,:)) - minj
            newpoint(j_) = minj + rangej * ranarr(i_)  ! UNIFORMLY DISTRIBUTED SAMPLE
          enddo
          evolution = 'mutation'
          call runmodel(newpoint(:), newobjfun)
          if (newobjfun .gt. objfun(qopt-1)) then
            do j_ = 1, qopt - 1
              if (newobjfun .gt. objfun(j_)) then
                objfun(j_+1:qopt) = objfun(j_:qopt-1)
                invar(:,j_+1:qopt) = invar(:,j_:qopt-1)
                objfun(j_) = newobjfun
                invar(:,j_) = newpoint(:)
                exit
              endif
            enddo
          elseif (newobjfun .gt. 0.d0) then       ! REPLACE ANYWAY
            objfun(qopt) = newobjfun
            invar(:,qopt) = newpoint(:)
          endif
        endif
      endif

!     * UPDATE 'currentbest' IF APPROPRIATE

      if (newobjfun .gt. bestobj) then
        bestobj = newobjfun
        bestincomp = newobjfun
        open(kfile_currentbest,file=sfile_currentbest(1:len_trim(sfile_currentbest)))
        write(kfile_currentbest,outformat) invar(:,1), bestobj
        close(kfile_currentbest)
        nsincebest = 0
      elseif (newobjfun .gt. bestincomp) then
        bestincomp = newobjfun
      endif

      return
      end subroutine simplex

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine runmodel (invar, objfun)
      use sce_mod
      implicit none

!     * Declarations
      REAL*8, DIMENSION(npar), INTENT(in) :: invar
      REAL*8, INTENT(out) :: objfun

!     * Definitions
      CHARACTER(1)   :: bestmark
      CHARACTER(300) :: writeformat

      nrun = nrun + 1
      call transpmodel(invar(:), npar, nrun, objfun, 1)
      bestmark = ' '
      if (objfun .gt. bestobj) then
        bestmark = '+'
      elseif (objfun .gt. bestincomp) then
        bestmark = '.'
      endif
      if (evolution .ne. 'test') then
        writeformat = '("Run",i6,"  (",a,"):",t28,"OF =",e12.6,1x,a)'
        write(*,writeformat) nrun, trim(evolution), objfun, bestmark
        write(kfile_progress,writeformat) nrun, trim(evolution), objfun, bestmark
        write(kfile_sceout,outformat) invar(optid(:)), objfun
        flush(kfile_progress)
      endif

      return
      end subroutine runmodel

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine sortcomp (invar, dim_invar, objfun, dim_objfun)
      implicit none

!     * Declarations
      INTEGER, DIMENSION(2), INTENT(in) :: dim_invar
      REAL*8, DIMENSION(dim_invar(1),dim_invar(2)), INTENT(inout) :: invar
      INTEGER, INTENT(in) :: dim_objfun
      REAL*8, DIMENSION(dim_objfun), INTENT(inout) :: objfun

!     * Definitions
      INTEGER :: i_, j_(1)

      REAL*8,  DIMENSION(:),   ALLOCATABLE :: objfun2
      REAL*8,  DIMENSION(:,:), ALLOCATABLE :: invar2
      INTEGER, DIMENSION(:),   ALLOCATABLE :: newobjfun

      allocate(objfun2(dim_objfun))
      allocate(invar2(dim_invar(1),dim_invar(2)))
      allocate(newobjfun(dim_objfun))

!     * EXTERNAL compar

      objfun2(:) = objfun(:)
      newobjfun(:) = -99                          ! NEEDED TO SEPARATE EQUAL O.F. VALUES
      invar2(:,:) = invar(:,:)
!     call qsort(objfun(:), dim, 8, compar)       ! USE compar(b,c)=(c-b)/dabs(c-b)
      call qsort(objfun(:), dim_objfun)
      do i_ = 1, dim_objfun
        j_ = minloc(real(dabs(objfun2(:) - objfun(i_))), newobjfun(:) .lt. 0)
        newobjfun(j_(1)) = i_
        invar(:,i_) = invar2(:,j_(1))
      enddo

      deallocate(objfun2)
      deallocate(invar2)
      deallocate(newobjfun)

      return
      end subroutine sortcomp

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!     * WRITE shufflevar AND ofvec OF LAST LOOP TO FILE AND TERMINATE
      subroutine writeloop ()
      use sce_mod
      implicit none

      INTEGER :: i_

      open(kfile_lastloop, file=sfile_lastloop(1:len_trim(sfile_lastloop)))
      write(kfile_lastloop,'(i3)') ncomp2
      write(kfile_lastloop,'(i4)') nloop
      write(kfile_lastloop,'(i10)') nrun
      write(kfile_lastloop,'(i10)') nsincebest
      write(kfile_lastloop,loopformat) ofvec(:)
      do i_ = 1, npar
        write(kfile_lastloop,loopformat) shufflevar(i_,:)
      enddo
      close(kfile_lastloop)

      return
      end subroutine writeloop

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine writepars ()
      use sce_mod
      implicit none

      INTEGER :: i_

      write(*,'(/"PARAMETER|     VALUE   |    MINVAL   |    MAXVAL   |  CV (%)     ")')
      write(kfile_progress,'(/"PARAMETER|     VALUE   |    MINVAL   |    MAXVAL   |  CV (%)     ")')
      do i_ = 1, nopt
        write(*,'(a9,5e14.6)') parname(optid(i_)), shufflevar(optid(i_),1), parmin(optid(i_)), parmax(optid(i_)), cv_(i_)
        write(kfile_progress,'(a9,5e14.6)') parname(optid(i_)), shufflevar(optid(i_),1), parmin(optid(i_)), parmax(optid(i_)), cv_(i_)
      enddo
      write(kfile_progress,'(//)')
      print *
      print *
      if (success .eq. 1) then
        open(kfile_finalbest, file=sfile_finalbest(1:len_trim(sfile_finalbest)))
        write(kfile_finalbest,outformat) shufflevar(:,1), bestobj
        close(kfile_finalbest)
      endif

      if (success .eq. 2) then
        open(kfile_finalbest, file=sfile_finalbest(1:len_trim(sfile_finalbest)))
        write(kfile_finalbest, '("   0.0E+00  0.0E+00  0.0E+00  0.0E+00  0.0E+00  0.0E+00  0.0E+00")')
        close(kfile_finalbest)
      endif

      return
      end subroutine writepars

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!     * SORTS THE VECTOR OBJFUN IN DESCENDING ORDER AND RETURNS IT
      subroutine qsort (objfun, dim_objfun)
      implicit none

!     * Declarations
      INTEGER :: dim_objfun
      REAL*8, DIMENSION(dim_objfun), INTENT(inout) :: objfun

!     * Definitions
      INTEGER :: i_, j_, k_

      do i_ = 2, dim_objfun
        if (objfun(i_) .gt. objfun(i_-1)) then
          k_ = i_ - 2
          do j_ = i_ - 2, 1, -1
            k_ = j_
            if (objfun(i_) .lt. objfun(j_)) exit
            if (j_ .eq. 1) k_ = 0   ! If objfun(i_)>objfun(1), then cycle objfun(i_) to the top
          enddo
          objfun(k_+1:i_) = cshift(objfun(k_+1:i_), -1)
        endif
      enddo

! For debugging:
!!$ do i_ = 1, dim_objfun
!!$  print *, objfun(i_)
!!$ enddo

      return
      end subroutine qsort

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!!$INTEGER function compar(b,c)
!!$  
!!$  implicit none
!!$  
!!$  REAL*8, INTENT(in) :: b,c
!!$  
!!$  if(b.gt.c) then                              ! SORTS b,c IN DECREASING ORDER
!!$     compar = -1
!!$  elseif(b.lt.c) then
!!$     compar = 1
!!$  else
!!$     compar = 0
!!$  endif
!!$  return
!!$end function compar
