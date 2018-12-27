!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! SHUFFLED COMPLEX EVOLUTION
!   Parameter optimisation algorithm based on a paper by Duan et al.
!   (1993, J. Opt. Theory and Appl., 76, 501--521).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! VERSION 2.1 ---  01 February 1999
! Written by Neil Viney, Centre for Water Research (CWR), The University of WA
! Modified by Stan Schymanski, CWR, 05 April 2004 (to run with transpmodel)
! Extended by Stan Schymanski, SESE, 02 June 2006 to follow Muttil & Liong
! (2004, Journal of Hydraulic Engineering-Asce 130(12):1202-1205) and
! Duan et al. (1994, Journal of Hydrology 158)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! This implementation MAXIMISES the objective function, which is
! calculated by the model, not by the optimiser. The optimiser transfers
! parameter values to the model subroutine and receives the value of the
! objective function.
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

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine sce_main ()
      use vom_sce_mod
      implicit none

      INTEGER             :: run_initialseed
      INTEGER             :: ii, first, numcv
      REAL*8              :: maxcv
      CHARACTER(300)      :: writeformat
      CHARACTER(len=135)  :: msg
      INTEGER             :: tmp2(2)
      CHARACTER(len=9)    :: tmp3(1)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! INITIALIZATION

      call sce_init(run_initialseed)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! BEGIN MODEL LOOP (EXECUTE s PASSES)
! CALL model SUBROUTINE        (OPEN run.log FOR CONSOLE OUTPUT);
! CALCULATE OBJECTIVE FUNCTION FOR DEFAULT PARAMETER VALUES
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      if (run_initialseed == 1) then
        call initialseed
        return
      endif

!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! BEGIN SCE LOOP

      do while (nrun .lt. 20000 .and. nloop .lt. 500)

          nloop = nloop + 1

!         * Saving the best OF of the worst complex in worstbest for
!         * assessment of gene pool mixing

          first = 1 + (ncomp2 - 1) * mopt
          worstbest = ofvec(first)
!         * [SORT ENTIRE ARRAYS]
!         * use temporary variable to prevent warning in ifort
          tmp2(:) = SHAPE(shufflevar(:,:))
          call sortcomp(shufflevar(:,:), tmp2(:), ofvec(:), SIZE(ofvec(:)))

!         * WRITE BEST_PARAMETERS FILE FOR PREVIOUS LOOP

          writeformat = '("Finished ",i4," main loops'
          writeformat(29:66) = ' --- best objective function =",e12.6)'
          write(msg,writeformat) nloop, ofvec(1)
          write(kfile_progress,*) TRIM(msg)
          write(kfile_progress,*) " "
          writeformat = '("No improvement in OF for",i5," loops")'
          write(msg,writeformat) nsincebest
          write(kfile_progress,*) TRIM(msg)
          nsincebest = nsincebest + 1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ASSESS CONVERGENCE
! CHANGED BY STAN TO CALCULATE THE FLUCTUATION RANGE RELATIVELY TO
! THE FEASIBLE RANGE, INSTEAD OF CV:

          numcv = ncomp2 * mopt
          sumvar(:) = SUM(shufflevar(optid(:), 1:numcv),2) / numcv  ! mean parameter values
          do ii = 1, nopt
!           * distance from mean in % of feasible range
            cv_(ii) = MAXVAL(ABS((shufflevar(optid(ii), 1:numcv)       &
     &              - sumvar(ii)) / (parmin(optid(ii))                 &
     &              - parmax(optid(ii)))) * 100.d0)
          enddo
          maxcv = MAXVAL(cv_(:))        ! maximum distance
          writeformat = '("Greatest parameter range: ",f5.2,"%'
          writeformat(38:67) = ' for optimised parameter ",a9)'
!         * use temporary variable to prevent warning in ifort
          tmp3(:) = parname(optid(MAXLOC(cv_(:))))
          write(msg,writeformat) maxcv, tmp3
          write(kfile_progress,*) TRIM(msg)
          if (maxcv .ge. i_resolution) then
            if (nsincebest .le. i_patience) then
              call writepars()

              call run_cce()

              return
            else
              write(kfile_progress,*) " "
              writeformat = '("No improvement in OF for",i5," loops")'
              write(msg,writeformat) nsincebest
              write(kfile_progress,*) TRIM(msg)
              write(kfile_progress,*) "  About to give up..."
            endif
          else
            write(kfile_progress,*) " "
            write(kfile_progress,*) "First Convergence criterion satisfied..."
            write(kfile_progress,*) "  parameter ranges are all less than 0.1 %"
            write(kfile_progress,*) " "
            write(kfile_progress,*) " "
          endif

!       * STAN'S MODIFICATION TO ASSESS SENSITIVITY OF OBJECTIVE
!       * FUNCTION TO EACH PARAMETER:

          call ck_success()

          if (success .eq. 1) then
            !close(kfile_sceout)
            close(kfile_bestpars)
            if (kfile_progress .ne. 6) close(kfile_progress)
            exit
          endif
      enddo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! TERMINATE PROGRAM

        if (success .ne. 1) then
          write(kfile_progress,*) " "
          write(kfile_progress,*) "FAILURE TO CONVERGE..."
          write(kfile_progress,*) "  Number of runs has reached 20000."
          write(kfile_progress,*) " "
          write(kfile_progress,*) "  Program terminated."
          write(msg,'(A,A)') CHAR(7), CHAR(7)
          write(kfile_progress,*) TRIM(msg)
!         * [SORT ENTIRE ARRAYS]
!         * use temporary variable to prevent warning in ifort
          tmp2(:) = SHAPE(shufflevar(:,:))
          call sortcomp(shufflevar(:,:), tmp2(:), ofvec(:), SIZE(ofvec(:)))
          call writepars()              ! PROGRAM STOP
          call write_lastbest(shufflevar(:,1), vom_npar, bestobj, 1)
            !close(kfile_sceout)
            close(kfile_bestpars)
          if (kfile_progress .ne. 6) close(kfile_progress)
        endif

      return
      end subroutine sce_main

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine sce_init (run_initialseed)
      use vom_sce_mod
      implicit none

      INTEGER, INTENT(out) :: run_initialseed

      INTEGER             :: ii, jj
      CHARACTER(3)        :: str
      CHARACTER(24)       :: logdate
      CHARACTER(len=135)  :: msg

!EXTERNAL compar
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Initialize the random number generator
! (standard subroutine, based on the date and time)
      if (vom_command .ne. 4) then
        CALL RANDOM_SEED()
      endif
      call fdate(logdate)
!only run the next part the first time for allocation
if( .not. allocated(optid) ) then
      call read_shufflevar()


   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! WRITE SCREEN HEADER

      write(*,*) 'SHUFFLED COMPLEX EVOLUTION OPTIMISER'

      write(*,*) " "
      write(msg,'("  Run time:   ",A)') logdate
      write(*,*) TRIM(msg)

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! READ PARAMETER FILE shufflevar


   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! CALCULATE NUMBER OF OPTIMISABLE PARAMETERS

      nopt = SUM(paropt(:))
      mopt = 2 * nopt + 1               ! SCE VARIABLE m
!     * number of complexes after Muttil(2004) or from shuffle.par
      i_ncomp_ = MAX(i_ncomp_, CEILING(1.d0 + 2.d0 ** nopt             &
     &         / (1.d0 + 2.d0 * nopt)))
      sopt = mopt * i_ncomp_            ! SCE VARIABLE s
      qopt = nopt + 1                   ! CCE VARIABLE q
      ncomp2 = i_ncomp_

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   ! ALLOCATE GLOBAL FIELDS



      allocate(optid(nopt))
      allocate(shufflevar(vom_npar,sopt))
      allocate(ofvec(sopt))
      allocate(wgt(mopt))
      allocate(cv_(nopt))
      allocate(ranarr(nopt))
      !allocate(ranarr_simplex(nopt))
      allocate(dataarray(nopt*8+1,nopt+1))
      allocate(shufflevar2(vom_npar))
      !allocate(parentsid(qopt))
      !allocate(objfunsub(qopt))
      !allocate(invarsub(vom_npar,qopt))
     ! allocate(selected(mopt))
     ! allocate(centroid(vom_npar))
     ! allocate(newpoint(vom_npar))
      allocate(posarray(2**nopt,nopt))
      allocate(initpop(nopt,5))
      allocate(sumvar(nopt))

end if

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! INITIALISE OUTPUT FORMAT STRING
! INTRODUCED BY STAN TO HAVE ONE COLUMN PER PARAMETER AND ONE FOR OF

      write(str,'(i2)') vom_npar + 1    ! internal write to convert from integer to string
      outformat = '('//str//'e24.15)'   ! includes a column for each parameter and a column for the value of OF

!     * ADDED BY STAN TO WRITE shufflevar AND ofvec OF LAST LOOP TO FILE
      write(str,'(i3)') sopt            ! internal write to convert from number to string
      loopformat = '('//str//'e24.15)'  ! includes a column for each set

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! INITIALISE optid: THE INDEX OF THOSE PARAMETERS THAT ARE OPTIMISABLE

      jj = 0
      do ii = 1, vom_npar
        if (paropt(ii) .gt. 0) then
          jj = jj + 1
          optid(jj) = ii
        endif
      enddo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! INITIALISE optid: THE INDEX OF THOSE PARAMETERS THAT ARE OPTIMISABLE
! ASSIGN PROBABILITY WEIGHTS

      do ii = 1, mopt
        wgt(ii) = 2.d0 * (mopt + 1.d0 - ii) / mopt / (mopt + 1.d0)
      enddo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! INSERTED BY STAN TO ALLOW CONTINUATION OF OPTIMSATION FROM PREVIOUSLY SAVED STEP

      call open_output (logdate, run_initialseed)

      return
      end subroutine sce_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine read_shufflepar ()
      use vom_sce_mod
      implicit none

      INTEGER :: iostat

!     * Definition of variable parameters

      namelist /shufflepar/ vom_command, i_ncomp_, i_ncompmin,         &
     &                      i_resolution, i_patience, i_nsimp,         &
     &                      i_focus, i_iter, vom_npar, n_thread

!     * Input of variable parameters from the parameter file

      open(kfile_namelist, FILE=sfile_namelist, STATUS='old',          &
     &                     FORM='formatted', IOSTAT=iostat)
      if (iostat .eq. 0) then
        read(kfile_namelist, shufflepar)
      endif
      close(kfile_namelist)

      if (vom_npar > nparmax) then
        write(0,*) "ERROR: Number of parameters in shufflevar larger as nparmax"
        write(0,*) "HINT: change the parameter nparmax in the module definitions"
        stop
      endif

      return
      end subroutine read_shufflepar

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine read_shufflevar ()
      use vom_sce_mod
      implicit none

      INTEGER :: iostat
      CHARACTER(9) :: parname0 (nparmax)
      REAL*8       :: parval0  (nparmax)
      REAL*8       :: parmin0  (nparmax)
      REAL*8       :: parmax0  (nparmax)
      INTEGER      :: paropt0  (nparmax)

!     * Definition of variable parameters

      namelist /shuffle2par/ parname0, parval0, parmin0, parmax0, paropt0

!     * LOAD INITIAL PARAMETER VALUES AND PARAMETER RANGES
      if( .not. allocated(parname) ) then
         open(kfile_namelist, FILE=sfile_namelist, STATUS='old',     &
     &                     FORM='formatted', IOSTAT=iostat)
         if (iostat .eq. 0) then
           read(kfile_namelist, shuffle2par)
         endif
         close(kfile_namelist)

!     * allocate and set the parameter fields

         allocate(parname(vom_npar))
         allocate(parval(vom_npar))
         allocate(parmin(vom_npar))
         allocate(parmax(vom_npar))
         allocate(paropt(vom_npar))
      end if
      parname(:) = parname0(1:vom_npar)
      parval(:)  = parval0(1:vom_npar)
      parmin(:)  = parmin0(1:vom_npar)
      parmax(:)  = parmax0(1:vom_npar)
      paropt(:)  = paropt0(1:vom_npar)

      return
      end subroutine read_shufflevar

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine open_output (logdate, run_initialseed)
      use vom_sce_mod
      implicit none

      CHARACTER(24), INTENT(in) :: logdate
      INTEGER,      INTENT(out) :: run_initialseed

      INTEGER             :: ii, iostat
      CHARACTER(len=135)  :: msg
      REAL*8, ALLOCATABLE :: tmp_8(:)

      run_initialseed = 0

        open(kfile_lastloop, FILE=trim(adjustl(i_outputpath)) // &
             trim(adjustl(sfile_lastloop)), STATUS='old', IOSTAT=iostat)
        if (iostat .eq. 0) then

!         * read parameter values to continue from last run

          read(kfile_lastloop,*) ncomp2
          read(kfile_lastloop,*) nloop
          read(kfile_lastloop,*) nrun
          read(kfile_lastloop,*) nsincebest
          read(kfile_lastloop,loopformat) ofvec(:)
!         * use temporary variable to prevent warning in ifort
          allocate(tmp_8(sopt))
          do ii = 1, vom_npar
            read(kfile_lastloop, loopformat) tmp_8(:)
            shufflevar(ii,:) = tmp_8(:)
          enddo
          deallocate(tmp_8)
          close(kfile_lastloop)
          bestobj = ofvec(1)

!         * open files for storing objective function and parameter values

            open(kfile_sceout, FILE=trim(adjustl(i_outputpath)) // &
             trim(adjustl(sfile_sceout)), STATUS='old', POSITION='append')
            open(kfile_bestpars, FILE=trim(adjustl(i_outputpath)) // &
             trim(adjustl(sfile_bestpars)), STATUS='old', POSITION='append')
          if (kfile_progress .ne. 6) then
            open(kfile_progress, FILE=trim(adjustl(i_outputpath)) // &
             trim(adjustl(sfile_progress)), STATUS='old', POSITION='append')
          endif
          write(kfile_progress,*) " "
          write(msg,'("  NEW Run time:   ",A)') logdate
          write(kfile_progress,*) TRIM(msg)
        else

!         * open empty files for storing objective function and parameter values

            open(kfile_sceout, FILE=trim(adjustl(i_outputpath)) // &
             trim(adjustl(sfile_sceout)), STATUS='replace')
            open(kfile_bestpars, FILE=trim(adjustl(i_outputpath)) // &
             trim(adjustl(sfile_bestpars)), STATUS='replace')
          if (kfile_progress .ne. 6) then
            open(kfile_progress, FILE=trim(adjustl(i_outputpath)) // &
             trim(adjustl(sfile_progress)), STATUS='replace')
          endif

!         * write output header

          write(kfile_progress,*) "SHUFFLED COMPLEX EVOLUTION OPTIMISER"
          write(kfile_progress,*) " "
          write(msg,'("  Run time:   ",A)') logdate
          write(kfile_progress,*) TRIM(msg)

          bestobj = 0.0
          bestincomp = 0.0
          run_initialseed = 1
        endif

      return
      end subroutine open_output

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine initialseed ()
      use vom_sce_mod
      implicit none

      INTEGER            :: ii
      INTEGER            :: jj, kk
      INTEGER            :: worstcount  ! worstcount for counting number of negative objective functions
      CHARACTER(len=135) :: msg

      do ii = 1,sopt

        if (ii .eq. 1) then

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! PRINT DIMENSION INFORMATION

            write(kfile_progress,*) '  Number of model parameters:         ',vom_npar
            write(kfile_progress,*) '  Number of optimisable parameters:   ',nopt
            write(kfile_progress,*) '  Maximum number of complexes:        ',i_ncomp_
            write(kfile_progress,*) '  Minimum number of runs per complex: ',mopt

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! CALCULATING OF USING INITAL GUESS IN SHUFFLE.PAR

            nsincebest = 0
            evolution = 'seed'
            ofvec(:) = -9999.9d0
            shufflevar(:,1) = parval(:)
            nrun = 0
            worstcount = 0

          call runmodel(shufflevar(:,1), ofvec(1))

            if (ofvec(1) .le. 0.d0) worstcount = worstcount + 1
            bestobj = ofvec(1)
            bestincomp = bestobj
            call write_lastbest(shufflevar(:,1), vom_npar, bestobj, 0)
            write(msg,'("Systematic seed of",i4," parameters for ",i2," complexes. Initial OF= ",e13.6)') nopt, i_ncomp_, ofvec(1)
            write(*,*) TRIM(msg)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! GENERATE A SYSTEMATIC ARRAY OF INITIAL PARAMETER VALUES FOLLOWING
! Muttil & Liong (2004, Journal of Hydraulic Engineering-Asce 130(12):1202-1205
! (INSERTED BY STAN)
!
! NONAXIAL POINTS:

            posarray(1,1) = 1
            posarray(2,1) = 2
            do jj = 1, nopt
              kk = optid(jj)
!             * each position j contains the intial perturbation of an optimised parameter
              initpop(jj,1) = 0.125d0 * parmax(kk) + 0.875d0 * parmin(kk)
              initpop(jj,2) = 0.125d0 * parmin(kk) + 0.875d0 * parmax(kk)
              initpop(jj,3) = 0.5d0 * (parmin(kk) + parmax(kk))
              initpop(jj,4) = 0.25d0 * parmax(kk) + 0.75d0 * parmin(kk)
              initpop(jj,5) = 0.25d0 * parmin(kk) + 0.75d0 * parmax(kk)
              posarray(2**(jj-1)+1:2**jj, 1:jj-1) = posarray(1:2**(jj-1), 1:jj-1)
              posarray(1:2**(jj-1), jj) = 1
              posarray(2**(jj-1)+1:2**jj, jj) = 2
            enddo
        endif

        if (ii .gt. 1) then
            shufflevar(:,ii) = parval(:)  ! TO SET NON-OPTIMISING PARAMETERS
        endif

        if (ii .gt. 1 .and. ii .le. 4) then
            do jj = 1, nopt
              kk = optid(jj)
              shufflevar(kk,ii) = initpop(jj,ii+1)
            enddo

          call runmodel(shufflevar(:,ii), ofvec(ii))

            if (ofvec(ii) .le. 0) worstcount = worstcount + 1
        endif

        if (ii .gt. 4 .and. ii .le. SIZE(posarray(:,:),1)) then
            do jj = 1, nopt
              kk = optid(jj)
              shufflevar(kk,ii) = initpop(jj, posarray(ii-4, jj))
            enddo

          call runmodel(shufflevar(:,ii), ofvec(ii))

            if (ofvec(ii) .le. 0) worstcount = worstcount + 1
        endif

!       * IF MORE POINTS ARE NEEDED, GENERATE RANDOM POINTS

        if (ii .gt. SIZE(posarray(:,:),1)) then
          evolution = 'mutation'

!         * first loop must generate feasible values to start with
          do while (ofvec(ii) .le. 0.d0)

              call random_number(ranarr(:))  ! RANDOM ARRAY OF SIZE 1:nopt
              do jj = 1, nopt
                kk = optid(jj)

!               * STAN'S MODIFICATION TO GET COMPLETELY RANDOM SEED:

                shufflevar(kk,ii) = parmin(kk) + i_focus               &
     &                            * (parmax(kk) - parmin(kk)) * ranarr(jj)
              enddo
              shufflevar(optid(:),ii) = MERGE(shufflevar(optid(:),ii), parmin(optid(:)), &
     &                                  shufflevar(optid(:),ii) .gt. parmin(optid(:)))
              shufflevar(optid(:),ii) = MERGE(shufflevar(optid(:),ii), parmax(optid(:)), &
     &                                  shufflevar(optid(:),ii) .lt. parmax(optid(:)))

            call runmodel(shufflevar(:,ii), ofvec(ii))

              if (ofvec(ii) .le. 0) worstcount = worstcount + 1
!             * program stops after 100 runs without positive objective function
              if (nrun .gt. 100 .and. nrun .eq. worstcount) then
                success = 2
                call writepars()
                exit
              endif

          enddo

        endif

          if (ofvec(ii) .gt. bestobj) then
            bestobj = ofvec(ii)
            call write_lastbest(shufflevar(:,ii), vom_npar, bestobj, 0)
          endif

        if (success .eq. 2) exit
      enddo

        nloop = -1                      ! FIRST LOOP IS LOOP ZERO
        call writeloop()
        !close(kfile_sceout)
        close(kfile_bestpars)
        if (kfile_progress .ne. 6) close(kfile_progress)

      return
      end subroutine initialseed

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine ck_success ()
      use vom_sce_mod
      implicit none

      INTEGER :: tmp2(2)

      call optsensitivity()

        if (success .eq. 1) then
!         * [SORT ENTIRE ARRAYS]
!         * use temporary variable to prevent warning in ifort
          tmp2(:) = SHAPE(shufflevar(:,:))
          call sortcomp(shufflevar(:,:), tmp2(:), ofvec(:), SIZE(ofvec(:)))
          call writepars()
          write(kfile_progress,*) 'Optimisation completed successfully.'
          !close(kfile_sceout)
          close(kfile_bestpars)
          if (kfile_progress .ne. 6) close(kfile_progress)
        else
!         * [SORT ENTIRE ARRAYS]
!         * use temporary variable to prevent warning in ifort
          tmp2(:) = SHAPE(shufflevar(:,:))
          call sortcomp(shufflevar(:,:), tmp2(:), ofvec(:), SIZE(ofvec(:)))
          call writepars()
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
      use vom_sce_mod
      implicit none

      INTEGER              :: ii, jj, kk, failed10, pos
      REAL*8               :: distmin, distmax
      REAL*8               :: oldpar, newpar
      REAL*8               :: ofvec2, ofchange, parchange
      CHARACTER(300)       :: writeformat
      CHARACTER(len=135)   :: msg
      REAL*8, ALLOCATABLE  :: tmp_8(:)

        shufflevar2(:) = shufflevar(:,1)
        ofvec2 = ofvec(1)
        failed10 = 0

        write(kfile_progress,*) 'SENSITIVITY ANALYSIS'
        write(kfile_progress,*) 'changes in % of feasible range)'
        dataarray(1,1:nopt) = shufflevar2(optid(:))
        dataarray(1,nopt+1) = ofvec2
        evolution = 'test'
        pos = 0

      do ii = 1, nopt

          oldpar = shufflevar2(optid(ii))

!         * TO MAKE SURE THAT PERTURBATIONS DO NOT EXCEED THE FEASIBLE RANGE
!         * AND THE PARAMETER VALUES THEMSELVES

          distmin = oldpar - parmin(optid(ii))
          distmax = parmax(optid(ii)) - oldpar
          write(kfile_progress,*) " "
          writeformat = '("change of var: ",A9, "(",E9.3,")")'
          write(msg,writeformat) parname(optid(ii)), oldpar
          write(kfile_progress,*) TRIM(msg)
          writeformat = '(f7.3,"% (",e14.7,")",": change of OF by ",'
          writeformat(44:68) = 'e10.3,"%"," (",e10.3,")")'

          jj = 0                        ! counts from 0 to 3 and than from 3 to 0
          kk = 1                        ! counter: first +1 after 3 steps it becomes -1

        do while (jj .ge. 0)

            if (kk .eq. 1) then
              newpar = oldpar - distmin * 10.d0 ** (-jj)
            else
              newpar = oldpar + distmax * 10.d0 ** (-jj)
            endif
            nrun = nrun + 1
            shufflevar2(optid(ii)) = newpar

          call transpmodel(shufflevar2(:), SIZE(shufflevar2(:)), ofvec2, 1)

            if (ofvec2 .gt. bestobj) then
              bestobj = ofvec2
              call write_lastbest(shufflevar2(:), vom_npar, bestobj, 0)
            endif

!             * use temporary variable to prevent warning in ifort
              allocate(tmp_8(nopt))
              tmp_8(:) = shufflevar2(optid(:))
              write(kfile_sceout,outformat) tmp_8(:), ofvec2
              deallocate(tmp_8)

            if (kk .eq. 1) then
              dataarray((ii-1)*8+jj+2,1:nopt) = shufflevar2(optid(:))
              dataarray((ii-1)*8+jj+2,nopt+1) = ofvec2
            else
              dataarray((ii-1)*8+jj+6,1:nopt) = shufflevar2(optid(:))
              dataarray((ii-1)*8+jj+6,nopt+1) = ofvec2
            endif
            ofchange = (ofvec2 - dataarray(1,nopt + 1))                &
     &               / ABS(dataarray(1,nopt + 1)) * 100.d0
            parchange = (newpar - oldpar) / (parmax(optid(ii))         &
     &                - parmin(optid(ii))) * 100.d0  ! the change of the parameter in % of feasible range
            write(msg,writeformat) parchange, newpar, ofchange, ofvec2
            write(kfile_progress,*) TRIM(msg)
            flush(kfile_progress)

            if (ofchange .gt. 1.0d-10) then
              shufflevar(:,sopt-pos) = shufflevar2(:)
              ofvec(sopt-pos) = ofvec2
              pos = pos + 1
              if (ABS(parchange) .gt. i_resolution) then
                failed10 = failed10 + 1
              endif
            endif

            jj = jj + kk
            if (jj .eq. 4) then
              jj = 3
              kk = -1
            endif
        enddo

          shufflevar2(optid(ii)) = oldpar

      enddo

!     * ASSESS SECOND CONVERGENCE CRITERIUM: NO PARAMETER CHANGE BY MORE THAN
!     * 10% LEADS TO AN INCREASE IN OBJECTIVE FUNCTION

        if (failed10 .gt. 0) then
          writeformat = '(I2," parameter(s) more than ",F6.3,"% out of optimum.")'
          write(msg,writeformat) failed10, i_resolution
          write(kfile_progress,*) TRIM(msg)
          write(kfile_progress,*) "Optimisation continued..."
        else
          write(kfile_progress,*) " "
          write(kfile_progress,*) " "
          write(kfile_progress,*) "Second convergence criterion satisfied:"
          writeformat = '("no parameter shift by more than ",F6.3,"% of max distance leads to")'
          write(msg,writeformat) i_resolution
          write(kfile_progress,*) TRIM(msg)
          write(kfile_progress,*) "an increase in the objective function."
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
      use vom_sce_mod
      use vom_vegwat_mod
  !$ USE omp_lib
      implicit none

      INTEGER              :: m_, first
      CHARACTER(300)       :: writeformat
      CHARACTER(len=135)   :: msg

!     * PARTITION THE sopt POINTS INTO ncomp2 COMPLEXES
!     * EXECUTE CCE ALGORITHM

        if (ncomp2 .gt. 1) then
          if (nloop .gt. 0) then
            if (ncomp2 .gt. i_ncompmin) then
              ncomp2 = ncomp2 - 1       ! REDUCE NUMBER OF COMPLEXES AS nloop INCREASES
            else
              if (worstbest .le. ofvec(1+(ncomp2-1)*mopt)) then
                writeformat = '("No gene pool mixing ... '
                writeformat(27:64) = 'reducing number of complexes by one.")'
                write(msg,writeformat)
                write(kfile_progress,*) TRIM(msg)
                ncomp2 = ncomp2 - 1
              endif
            endif
          endif
        endif

      bestincomp = -9999.d0       ! SET LESS THAN bestobj
      write(kfile_progress,*) "Looping over complexes"

      flush(kfile_progress)
      call OMP_SET_NUM_THREADS(n_thread)
      call vom_dealloc()

      !loop over complexes
      !!$OMP shared( ofvec)
      !$OMP parallel default(shared) &
      !$OMP private( m_, first, msg, writeformat) &
      !$OMP COPYIN( time, error, finish, nyear, nday, nhour, th_, c_testday,   & 
      !$OMP topt_, par_y, srad_y,  vd_d, vd_y, &
      !$OMP rain_y, gammastar, wsnew, wsold, o_pct, pcg_d, c_pcgmin, &
      !$OMP o_wstexp, o_wsgexp, o_lambdatf, o_lambdagf, gstomt, gstomg, &
      !$OMP rlt_h, rlt_d, rlt_y, rlg_h, rlg_d, rlg_y, transpt, transpg, q_tct_d, tct_y, tcg_d, &
      !$OMP tcg_y, jactt, jactg, jmaxt_h, jmaxg_h, jmax25t_d, jmax25g_d, &
      !$OMP asst_h, asst_d, asst_y, assg_h, assg_d, assg_y, &
      !$OMP q_cpcct_d, cpcct_y, cpccg_d, cpccg_y, etmt__, etmt_h, etmt_d, etmt_y, etmg__, etmg_h, &
      !$OMP etmg_d, etmg_y, etm_y, mqt_, mqtnew, mqtold, dmqt, q_mqx, mqsst_, mqsstmin, q_md, &
      !$OMP o_mdstore, o_rtdepth, o_rgdepth, pos_slt, pos_slg, pos_ult, pos_ulg, changef, &
      !$OMP rootlim, posmna, rrt_d, rrt_y, rrg_d, rrg_y, sumruptkt_h,  &
      !$OMP wlayer_, wlayernew, dt_, dtmax, dtsu_count, dtmax_count, esoil__, esoil_h, &
      !$OMP esoil_d, esoil_y, spgfcf__, spgfcf_h, spgfcf_d, inf__, infx__, infx_h, infx_d, &
      !$OMP zw_, zwnew, wc_,  io__, io_h, ioacum, &
      !$OMP ranscal, bestobj, bestincomp, evolution)
      !$OMP do SCHEDULE(DYNAMIC)
      do m_ = 1, ncomp2

          first = 1 + (m_ - 1) * mopt
          writeformat = '("Start of loop",i4,", complex",i2,'
          writeformat(36:55) = '": best OF =",e12.6)'
          write(msg,writeformat) nloop + 1, m_, ofvec(first)
          write(kfile_progress,*) TRIM(msg)
          flush(kfile_progress)

          if (m_ .eq. 1) then
            bestincomp = -9999.d0       ! SET LESS THAN bestobj
          else
            bestincomp = ofvec(first)
          endif

        !start CCE  
        call cce(ofvec(first:m_*mopt), shufflevar(:,first:m_*mopt))

        !Deallocate some variables to avoid problems in parallel
        call vom_dealloc()

          writeformat(3:7) = '  End'
          write(msg,writeformat) nloop + 1, m_, ofvec(first)
          write(kfile_progress,*) TRIM(msg)
          flush(kfile_progress)

      enddo
      !$OMP end do
      !$OMP end parallel

        worstbest = ofvec(first)

!       * WRITE shufflevar AND ofvec OF LAST LOOP TO FILE
        call writeloop()
        !close(kfile_sceout)
        close(kfile_bestpars)
        if (kfile_progress .ne. 6) close(kfile_progress)

      return
      end subroutine run_cce

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine cce (objfun, invar)
      use vom_sce_mod

      implicit none

!     * Declarations
      REAL*8, DIMENSION(mopt),          INTENT(inout) :: objfun
      REAL*8, DIMENSION(vom_npar,mopt), INTENT(inout) :: invar

!     * Definitions
      INTEGER       :: l_
      INTEGER       :: ii, nsel, rannum
      INTEGER       :: tmp2(2)
      REAL*8  :: ranscal2                ! Scalar random number
      LOGICAL, DIMENSION(mopt) :: selected
      INTEGER, DIMENSION(qopt) :: parentsid  ! Pointer to optimisable parameters
      REAL*8, DIMENSION(qopt) :: objfunsub  ! Pointer to optimisable parameters
      REAL*8, DIMENSION(vom_npar,qopt) :: invarsub  ! Pointer to optimisable parameters

!     * SELECT PARENTS

      do l_ = 1, mopt

          selected(:) = .false.
          nsel = 0

          do while (nsel .ne. qopt)
            call random_number(ranscal2)      ! SCALAR RANDOM NUMBER
            rannum = CEILING((2.d0 * mopt + 1.d0 - sqrt(4.d0 * mopt    &
     &             * (mopt + 1.d0) * (1.d0 - ranscal2) + 1.d0)) * 0.5d0)

!           * NOTE: A SIMPLER ALTERNATIVE TO THE ABOVE LINE IS (19.03.2004):

            if (rannum .ge. 1 .and. rannum .le. mopt) then
              if (.not. selected(rannum)) then
                selected(rannum) = .true.
                nsel = nsel + 1
              endif
            endif
          enddo

          nsel = 0
          ii = 0
          do while (nsel .ne. qopt)
            ii = ii + 1
            if (selected(ii)) then
              nsel = nsel + 1
              parentsid(nsel) = ii
            endif
          enddo

!         * GENERATE OFFSPRING AND SORT THE RESULTING COMPLEX

          objfunsub(:) = objfun(parentsid(:))
          invarsub(:,:) = invar(:, parentsid(:))

        call simplex(invarsub(:,:), objfunsub(:))

          objfun(parentsid(:)) = objfunsub(:)
          invar(:, parentsid(:)) = invarsub(:,:)
!         * use temporary variable to prevent warning in ifort
          tmp2(:) = SHAPE(invar)
          call sortcomp(invar, tmp2(:), objfun, SIZE(objfun))

      enddo

      return
      end subroutine cce

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine simplex (invar, objfun)
      use vom_sce_mod

      implicit none

!     * Declarations
      REAL*8, DIMENSION(vom_npar,qopt), INTENT(inout) :: invar
      REAL*8, DIMENSION(qopt),          INTENT(inout) :: objfun

!     * Definitions
      INTEGER       :: l_, kk
      INTEGER       :: ii, jj
      INTEGER       :: max_j_
      REAL*8        :: newobjfun
      REAL*8        :: minj, rangej
      REAL*8, DIMENSION(vom_npar) :: centroid  ! Centroid of parameter sets for simplex procedure
      REAL*8, DIMENSION(vom_npar) :: newpoint  ! New parameter set resulting from simplex procedure
      REAL*8, ALLOCATABLE :: ranarr_simplex(:)  ! Array of random numbers

      allocate(ranarr_simplex(nopt))

      do l_ = 1, i_nsimp

!         * REFLECTION STEP

          evolution = 'reflection'
          newpoint(:) = invar(:,qopt)
          centroid(optid(:)) = 1.d0 / (qopt - 1) * SUM(invar(optid(:), 1:qopt - 1), 2)
          newpoint(optid(:)) = 2.d0 * centroid(optid(:)) - invar(optid(:),qopt)
          if (MINVAL(newpoint(:) - parmin(:)) .lt. 0.d0 .or.           &
     &        MAXVAL(newpoint(:) - parmax(:)) .gt. 0.d0) then

!           * MUTATION STEP
!           * NB: MUTATION IS BASED ON THE SMALLEST HYPERCUBE OF THE SUBCOMPLEX,
!           * NOT THE SMALLEST HYPERCUBE OF THE COMPLEX AS SUGGESTED BY DUAN ET AL.

            do ii = 1, nopt
              call random_number(ranscal)  ! SCALAR RANDOM NUMBER
              jj = optid(ii)
              minj = parmin(jj)
              rangej = parmax(jj) - minj
              newpoint(jj) = minj + rangej * ranscal  ! UNIFORMLY DISTRIBUTED SAMPLE
            enddo
            evolution = 'mutation'
          endif

          kk = 1
        do while (kk .le. 3)

!         * CALCULATE OBJECTIVE FUNCTION
          call runmodel(newpoint(:), newobjfun)

            jj = 1
            if (kk .eq. 1) then
              if (newobjfun .gt. objfun(qopt)) then
                max_j_ = qopt
                kk = 4
              else
!               * CONTRACTION STEP
                newpoint(optid(:)) = (centroid(optid(:)) + invar(optid(:),qopt)) * 0.5d0
                evolution = 'contraction'
              endif
            endif

            if (kk .eq. 2) then
              if (newobjfun .gt. objfun(qopt)) then
                max_j_ = qopt
                kk = 4
              else
!               * ANOTHER MUTATION STEP
                do ii = 1, nopt
                  call random_number(ranarr_simplex(:))  ! RANDOM ARRAY OF SIZE 1:nopt
                  jj = optid(ii)
                  minj = MINVAL(invar(jj,:))
                  rangej = MAXVAL(invar(jj,:)) - minj
                  newpoint(jj) = minj + rangej * ranarr_simplex(ii)  ! UNIFORMLY DISTRIBUTED SAMPLE
                enddo
                evolution = 'mutation'
              endif
            endif

            if (kk .eq. 3) then
              if (newobjfun .gt. objfun(qopt-1)) then
                max_j_ = qopt - 1
                kk = 4
              elseif (newobjfun .gt. 0.d0) then  ! REPLACE ANYWAY
                objfun(qopt) = newobjfun
                invar(:,qopt) = newpoint(:)
              endif
            endif

            if (kk .eq. 4) then
!            * SORT objfun HERE IN CASE alpha > 1
              do while (newobjfun .le. objfun(jj) .and. jj .le. max_j_)
                jj = jj + 1
              enddo
              if (jj .lt. qopt) then
                objfun(jj+1:qopt) = objfun(jj:qopt-1)
                invar(:,jj+1:qopt) = invar(:,jj:qopt-1)
              endif
              if (jj .le. qopt) then
                objfun(jj) = newobjfun
                invar(:,jj) = newpoint(:)
              endif
            endif

            kk = kk + 1

          enddo

!         * UPDATE 'currentbest' IF APPROPRIATE

          if (newobjfun .gt. bestobj) then
            bestobj = newobjfun
            bestincomp = newobjfun
            call write_lastbest(invar(:,1), vom_npar, bestobj, 0)
            nsincebest = 0
          elseif (newobjfun .gt. bestincomp) then
            bestincomp = newobjfun
          endif

      enddo

      deallocate(ranarr_simplex)

      return
      end subroutine simplex

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine runmodel (invar, objfun)
      use vom_sce_mod

      implicit none

!     * Declarations
      REAL*8, DIMENSION(vom_npar), INTENT(in)    :: invar
      REAL*8,                      INTENT(out)   :: objfun

!     * Definitions
      CHARACTER(1)        :: bestmark
      CHARACTER(300)      :: writeformat
      CHARACTER(len=135)  :: msg
      REAL*8, ALLOCATABLE :: tmp_8(:)

        nrun = nrun + 1

      call transpmodel(invar(:), vom_npar, objfun, 1)

        bestmark = ' '
        if (objfun .gt. bestobj) then
          bestmark = '+'
        elseif (objfun .gt. bestincomp) then
          bestmark = '.'
        endif
        if (evolution .ne. 'test') then
          writeformat = '("Run",i6,"  (",a,"):",t28,"OF =",e12.6,1x,a)'
          write(msg,writeformat) nrun, TRIM(evolution), objfun, bestmark
          write(kfile_progress,*) TRIM(msg)

!           * use temporary variable to prevent warning in ifort
            allocate(tmp_8(nopt))
            tmp_8(:) = invar(optid(:))
            write(kfile_sceout,outformat) tmp_8(:), objfun
            flush(kfile_progress)
            deallocate(tmp_8)
        endif

      return
      end subroutine runmodel

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine sortcomp (invar, dim_invar, objfun, dim_objfun)
      implicit none

!     * Declarations
      INTEGER, DIMENSION(2),         INTENT(in)    :: dim_invar
      REAL*8, DIMENSION(dim_invar(1),dim_invar(2)), INTENT(inout) :: invar
      INTEGER,                       INTENT(in)    :: dim_objfun
      REAL*8, DIMENSION(dim_objfun), INTENT(inout) :: objfun

!     * Definitions
      INTEGER :: ii, jj(1)

      REAL*8,  DIMENSION(:),   ALLOCATABLE :: objfun2
      REAL*8,  DIMENSION(:,:), ALLOCATABLE :: invar2
      INTEGER, DIMENSION(:),   ALLOCATABLE :: newobjfun

      allocate(objfun2(dim_objfun))
      allocate(invar2(dim_invar(1),dim_invar(2)))
      allocate(newobjfun(dim_objfun))

!     * EXTERNAL compar

      objfun2(:) = objfun(:)
      newobjfun(:) = -99                ! NEEDED TO SEPARATE EQUAL O.F. VALUES
      invar2(:,:) = invar(:,:)
!     call qsort(objfun(:), dim, 8, compar)  ! USE compar(b,c)=(c-b)/dabs(c-b)
      call qsort(objfun(:), dim_objfun)
      do ii = 1, dim_objfun
        jj = MINLOC(REAL(dabs(objfun2(:) - objfun(ii))), newobjfun(:) .lt. 0)
        newobjfun(jj(1)) = ii
        invar(:,ii) = invar2(:,jj(1))
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
      use vom_sce_mod
      implicit none

      INTEGER             :: ii
      REAL*8, ALLOCATABLE :: tmp_8(:)

        open(kfile_lastloop, FILE=trim(adjustl(i_outputpath)) // &
             trim(adjustl(sfile_lastloop)))
          write(kfile_lastloop,'(i3)')  ncomp2
          write(kfile_lastloop,'(i4)')  nloop
          write(kfile_lastloop,'(i10)') nrun
          write(kfile_lastloop,'(i10)') nsincebest
          write(kfile_lastloop,loopformat) ofvec(:)
          allocate(tmp_8(sopt))
          do ii = 1, vom_npar
!           * use temporary variable to prevent warning in ifort
            tmp_8(:) = shufflevar(ii,:)
            write(kfile_lastloop,loopformat) tmp_8(:)
          enddo
          deallocate(tmp_8)
        close(kfile_lastloop)

      return
      end subroutine writeloop

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine writepars ()
      use vom_sce_mod
      implicit none

      INTEGER            :: ii
      CHARACTER(len=135) :: msg

      write(kfile_progress,*) " "
      write(msg,'("PARAMETER|     VALUE   |    MINVAL   |    MAXVAL   |  CV (%)     ")')
      write(kfile_progress,*) TRIM(msg)
      do ii = 1, nopt
        write(msg,'(a9,5e14.6)') parname(optid(ii)), shufflevar(optid(ii),1), parmin(optid(ii)), parmax(optid(ii)), cv_(ii)
        write(kfile_progress,*) TRIM(msg)
      enddo
      write(kfile_progress,*) ' '

      if (success .eq. 1) then
        call write_lastbest(shufflevar(:,1), vom_npar, bestobj, 1)
      endif

      if (success .eq. 2) then
        call write_lastbest(shufflevar(:,1)*0.d0, vom_npar, 0.d0, 1)
      endif

      return
      end subroutine writepars

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine write_lastbest (var, nvar, obj, final)
      use vom_sce_mod
      implicit none

      INTEGER                              :: nvar
      REAL*8, DIMENSION(nvar) , INTENT(in) :: var
      REAL*8                               :: obj
      INTEGER                              :: final

      REAL         :: stat
      CHARACTER*80 :: sfile

        open(kfile_lastbest, FILE=trim(adjustl(i_outputpath)) // &
             trim(adjustl(sfile_lastbest)))
          write(kfile_lastbest,outformat) var(:), obj
        close(kfile_lastbest)

        if (final == 1) then
          open(kfile_beststat, FILE=trim(adjustl(i_outputpath)) // &
             trim(adjustl(sfile_beststat)))
            write(kfile_beststat,*) 1
          close(kfile_beststat)
        endif

      return
      end subroutine write_lastbest

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     * SORTS THE VECTOR OBJFUN IN DESCENDING ORDER AND RETURNS IT

      subroutine qsort (objfun, dim_objfun)
      implicit none

!     * Declarations
      INTEGER,                       INTENT(in)    :: dim_objfun
      REAL*8, DIMENSION(dim_objfun), INTENT(inout) :: objfun

!     * Definitions
      INTEGER :: ii, jj, kk

      do ii = 2, dim_objfun
        if (objfun(ii) .gt. objfun(ii-1)) then
          kk = ii - 2
          do jj = ii - 2, 1, -1
            kk = jj
            if (objfun(ii) .lt. objfun(jj)) exit
            if (jj .eq. 1) kk = 0       ! If objfun(ii)>objfun(1), then cycle objfun(ii) to the top
          enddo
          objfun(kk+1:ii) = CSHIFT(objfun(kk+1:ii), -1)
        endif
      enddo

! For debugging:
!!$ do ii = 1, dim_objfun
!!$  print *, objfun(ii)
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
!!$  if(b.gt.c) then                     ! SORTS b,c IN DECREASING ORDER
!!$     compar = -1
!!$  elseif(b.lt.c) then
!!$     compar = 1
!!$  else
!!$     compar = 0
!!$  endif
!!$  return
!!$end function compar
