	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	!		
	!  SHUFFLED COMPLEX EVOLUTION
	!    Parameter optimisation algorithm based on a paper by Duan et al. (1993,
	!    J. Opt. Theory and Appl., 76, 501--521). 
	!
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	!  VERSION 2.1	---  01 February 1999 	
	!  Written by:	Neil Viney, Centre for Water Research (CWR), The University of WA
	!  Modified by Stan Schymanski, CWR, 05 April 2004 (to run with transpmodel)  
	!  Extended by Stan Schymanski, SESE, 02 June 2006 to follow Muttil & Liong (2004,
	!  Journal of Hydraulic Engineering-Asce 130(12):1202-1205) and Duan et al. (1994, 
	!  Journal of Hydrology 158)
	!	
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	!
	!  This implementation MAXIMISES the objective function, which is calculated
	!  by the model, not by the optimiser.	The optimiser transfers parameter values
	!  to the model subroutine and receives the value of the objective function.
	!
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	!
	!  Compile with: Compaq Visual Fortran 5.5	
	!
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	!
	!  Copyright (C) 2006 Stan Schymanski
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
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	subroutine sce(success)
		
	implicit none

	INTEGER, PARAMETER :: double=SELECTED_REAL_KIND(13) 	! == REAL*8
	INTEGER :: i,j,k,m,success
	INTEGER :: npar,ncomp,ncompmin,ncomp2,nopt,mopt,sopt,qopt,alpha,&
				beta,nrun,nloop,first
	INTEGER :: numcv,ios,nsincebest,patience
	INTEGER, DIMENSION(:), ALLOCATABLE :: paropt,optid
	INTEGER, DIMENSION(:,:), ALLOCATABLE :: posarray
	INTEGER, DIMENSION(3) :: iarray
	INTEGER, DIMENSION(1) :: rseed
	real(double), DIMENSION(:), ALLOCATABLE :: wgt,cv,ranarr
	real(double) :: maxcv,ranscal,focus
!	REAL(double), DIMENSION(:,:), ALLOCATABLE :: invar
!	REAL(double), DIMENSION(:), ALLOCATABLE :: parval,parmin,parmax,objfun,sumvar
!	REAL(double) :: worstbest,bestobj,bestincomp
	real(double), DIMENSION(:,:), ALLOCATABLE :: invar,initpop
	real(double), DIMENSION(:), ALLOCATABLE :: parval,parmin,parmax,&
				objfun,sumvar
	real(double) :: worstbest,bestobj,bestincomp,resolution
	CHARACTER(9), DIMENSION(:), ALLOCATABLE :: parname
	CHARACTER(15) :: command
	CHARACTER(24) :: logdate
	CHARACTER(16) :: filename
	CHARACTER(12) :: evolution
	CHARACTER(3) :: str
	CHARACTER(60) :: outformat,informat,loopformat
	LOGICAL :: newbest=.true.
	EXTERNAL compar
	integer(4) status

	!
	! WRITE SCREEN HEADER
	!
	write(*,'(/"SHUFFLED COMPLEX EVOLUTION OPTIMISER")')
	call fdate(logdate)
	write(*,'(/"  Run time:   ",a)') logdate
	!
	!
	!  FIND AND OPEN PARAMETER FILE, shuffle.par
	!    ESTABLISH NUMBER OF PARAMETERS npar
	!
	open(1,file='shuffle.par',status='old')
	read(1,*)
	read(1,*)
	npar=0
	do
	  read(1,*,iostat=ios) str
	  if(ios.lt.0) exit
	  if(str.eq.'var') npar=npar+1
	end do
	rewind(1)
	allocate(parname(npar),parval(npar),parmin(npar),parmax(npar),&
		paropt(npar))
	!
	!
	!  LOAD INITIAL PARAMETER VALUES AND PARAMETER RANGES
	!
	read(1,'(a15)') command
	read(1,*) ncomp
	read(1,*) ncompmin
	read(1,*) resolution
	read(1,*) patience
	read(1,*) alpha
	read(1,*) focus
	read(1,'(a60)') informat
	do i=1,npar
	  read(1,informat) parname(i),parval(i),parmin(i), &
		parmax(i),paropt(i)
	end do
	close(1)
	write(*,'(/"  Model command is:",/,4x,a15,//)') command
	!
	!
	!  CALCULATE NUMBER OF OPTIMISABLE PARAMETERS
	!
	nopt=sum(paropt)
	mopt=2*nopt+1 			! SCE VARIABLE m
	ncomp=Max(ncomp,Ceiling(1. + 2.**nopt/(1. + 2.*nopt)))		! number of complexes after Muttil(2004) or from shuffle.par
	sopt=mopt*ncomp 			! SCE VARIABLE s
	qopt=nopt+1 			! CCE VARIABLE q
	beta=mopt
	ncomp2=ncomp
	allocate(optid(nopt),invar(npar,sopt),objfun(sopt),wgt(mopt),&
		cv(nopt),sumvar(nopt),ranarr(nopt),posarray(2**nopt,nopt),&
		initpop(nopt,5))
	!
	!  PRINT DIMENSION INFORMATION
	!
	write(*,'("  Number of model parameters:",10x,i3)') npar
	write(*,'("  Number of optimisable parameters:",i7)') nopt
	write(*,'("  Maximum number of complexes:",9x,i3)') ncomp
	write(*,'("  Minimum number of runs per complex:",i5,//)') mopt
	!
	!
	! OPEN AND EMPTY FILE FOR STORING OBJECTIVE FUNCTION AND PARAMETER VALUES
	!
!	open(8,file='sce.out',access='direct',recl=8*(nopt+1))
	open(8,file='sce.out')
	close(8,status='delete')
	open(3,file='bestpars.txt')
	close(3,status='delete')
	open(2,file='progress.txt')
	close(2,status='delete')
!	open(8,file='sce.out',access='direct',recl=8*(nopt+1))
!	write(8,rec=1) nopt
	open(8,file='sce.out')
	open(3,file='bestpars.txt')
	open(2,file='progress.txt')

	
	write(2,'(/"SHUFFLED COMPLEX EVOLUTION OPTIMISER")')
	write(2,'(/"  Run time:   ",a)') logdate
	write(2,'(/"  Model command is:",/,4x,a15,//)') command
	write(2,'("  Number of model parameters:",10x,i3)') npar
	write(2,'("  Number of optimisable parameters:",i7)') nopt
	write(2,'("  Maximum number of complexes:",9x,i3)') ncomp
	write(2,'("  Minimum number of runs per complex:",i5,//)') mopt
	!
	! INITIALISE OUTPUT FORMAT STRING
	!
	!outformat='(3f12.5,/,2(f9.3,f12.3,/),2f11.5,/,2(f9.3,/),2f7.2,/,2(2f9.5,/),'// &
	!		 'f9.5,f12.5,/,2(2f10.6,/),f9.5,/,2(2f9.5,/),2f9.3,/,f8.3,/,2f9.5,/,'// &
	!		 'f8.5,/,f8.5,f8.2,f8.5,/,3f13.4,/,f13.4,/,3(f8.4,f12.4,/),f9.5,/,'// &
	!		 'f10.4,f10.5,/,f10.6,f13.9,/,2(2f9.5,/),f9.5,f11.7,/,'// & ! END OF SEDS
	!		 'f9.5,/,f7.4,/,f8.6,f9.6,f6.2,/,f4.2,f9.3,f11.7,/,f5.3,f8.3,/,'// &
	!		 'f8.6,/,f5.3,f8.3,/,f10.8,f7.4,/,f7.5,f8.5,/,f7.5,f7.4,/,'// &
	!		 'f8.6,f6.3,/,f8.6,f9.6,/,f4.2,f9.3,f11.7,/,f5.3,f8.3,/,f8.6,/,'// &
	!		 '3(f9.4,/))' 	! ONLY 3 RESERVOIRS!

	!INTRODUCED BY STAN TO HAVE ONE COLUMN PER PARAMETER AND ONE FOR OF
  write(str,'(i2)') npar+1				! internal write to convert from integer to string
  outformat='('//str//'e34.25)'		! includes a column for each parameter and a column for the value of OF
	!	ADDED BY STAN TO WRITE INVAR AND OBJFUN OF LAST LOOP TO FILE
	write(str,'(i3)') sopt		!internal write to convert from number to string
	loopformat='('//str//'e34.6)'		! includes a column for each sublayer 
	!
	!
	! INITIALISE optid: THE INDEX OF THOSE PARAMETERS THAT ARE OPTIMISABLE
	!
	j=0
	do i=1,npar
	  if(paropt(i).gt.0) then
	    j=j+1
	    optid(j)=i
	  endif
	end do
	!
	!
	!  ASSIGN PROBABILITY WEIGHTS
	!
	do i=1,mopt
	  wgt(i)=2.*(mopt+1.-i)/mopt/(mopt+1.)
	end do
	!
	!
	!  INITIALISE RANDOM NUMBER SEED: A NEW SEED EVERY TIME SCE IS INVOKED
	!
	call itime(iarray)
	rseed=(1+iarray(1))*(1+iarray(2))*(1+iarray(3))
	call random_seed(put=rseed)
	!
	!
	!  BEGIN MODEL LOOP (EXECUTE s PASSES)
	!    CALL model SUBROUTINE	(OPEN run.log FOR CONSOLE OUTPUT);
	!    CALCULATE OBJECTIVE FUNCTION FOR DEFAULT PARAMETER VALUES
	!
! INSERTED BY STAN TO ALLOW CONTINUATION OF OPTIMSATION FROM PREVIOUSLY SAVED STEP
	if (command.eq.'-continue      ') then
		open(9,file='lastloop.txt')
		read(9,*) ncomp2
		read(9,loopformat) objfun
		do i=1,npar
			read(9,loopformat) invar(i,:)
		enddo
		close(9)
		bestobj=objfun(1)
		goto 101
	endif
! END OF INSERTION


	nsincebest=0
	evolution='seed'
	objfun=-9999.9
	invar(:,1)=parval
	nrun=0
    call runmodel(invar(:,1),optid,objfun(1),bestobj,bestincomp)
	bestobj=objfun(1)
	bestincomp=bestobj
	open(4,file='currentbest.txt')
	!write(4,outformat) invar(:,1)
	!write(4,outformat) bestobj
	write(4,outformat) invar(:,1),bestobj
!	write(4,*) invar(:,1),bestobj
	close(4)
	write(*,'("Systematic seed of",i4," parameters for ",i2," complexes. Initial OF= ",e12.6)')	&
			nopt,ncomp,objfun(1)

	!	GENERATE A SYSTEMATIC ARRAY OF INITIAL PARAMETER VALUES FOLLOWING
	!	Muttil & Liong (2004, Journal of Hydraulic Engineering-Asce 130(12):1202-1205
	!	(INSERTED BY STAN)
	


	! NONAXIAL POINTS:

	posarray(1,1)=1
	posarray(2,1)=2
	do j=1,nopt
		k=optid(j)
		initpop(j,1)=0.125*parmax(j)+0.875*parmin(j)
		initpop(j,2)=0.125*parmin(j)+0.875*parmax(j)
		initpop(j,3)=0.5*(parmin(j)+parmax(j))
		initpop(j,4)=0.25*parmax(j)+0.75*parmin(j)
		initpop(j,5)=0.25*parmin(j)+0.75*parmax(j)

		posarray(2**(j-1)+1:2**j,1:j-1)=posarray(1:2**(j-1),1:j-1)
		posarray(1:2**(j-1),j)=1
		posarray(2**(j-1)+1:2**j,j)=2
	enddo

	do i=2,sopt
		invar(:,i)=parval		! TO SET NON-OPTIMISING PARAMETERS
		if(i.le.4) then
			do j=1,nopt
				k=optid(j)
				invar(k,i)=initpop(j,i+1)
			enddo
			call runmodel(invar(:,i),optid,objfun(i),bestobj,bestincomp)
			if(objfun(i).gt.bestobj) then
				bestobj=objfun(i)
				open(4,file='currentbest.txt') 
				write(4,outformat) invar(:,i),bestobj
				close(4)
			endif
		elseif(i.le.Size(posarray,1)) then
			do j=1,nopt
				k=optid(j)
				invar(k,i)=initpop(j,posarray(i-4,j))
			enddo
			call runmodel(invar(:,i),optid,objfun(i),bestobj,bestincomp)
			if(objfun(i).gt.bestobj) then
				bestobj=objfun(i)
				open(4,file='currentbest.txt') 
				write(4,outformat) invar(:,i),bestobj
				close(4)
			endif	
		else

!	write(*,'("Start of loop",i4,", complex",i2,": best OF =",e12.6)')	&
!			nrun+1,nrun,objfun(1)		
	!
	!
	!    IF MORE POINTS ARE NEEDED, GENERATE RANDOM POINTS
	!
			evolution='mutation'
			  do while(objfun(i).le.0.0)		! first loop must generate feasible values to start with
			    call random_number(ranarr)			! RANDOM ARRAY OF SIZE 1:nopt
			    do j=1,nopt
				    k=optid(j)
	! STAN'S MODIFICATION TO GET COMPLETELY RANDOM SEED:
					invar(k,i)=parmin(k) + focus*(parmax(k)-parmin(k))*ranarr(j)
			    end do
			    invar(optid,i)=merge(invar(optid,i),parmin(optid), &
										invar(optid,i).gt.parmin(optid))
			    invar(optid,i)=merge(invar(optid,i),parmax(optid), &
										invar(optid,i).lt.parmax(optid))
				call runmodel(invar(:,i),optid,objfun(i),bestobj,bestincomp)
			 end do
		  if(objfun(i).gt.bestobj) then
			bestobj=objfun(i)
			open(4,file='currentbest.txt') 
				write(4,outformat) invar(:,i),bestobj
	!		write(4,*) invar(:,i),bestobj
			close(4)
		  endif
		endif
	end do
	!
	!
	!  BEGIN SCE LOOP
	!
101	nloop=-1					! FIRST LOOP IS LOOP ZERO

	do while(nrun.lt.20000.and.nloop.lt.500)
	  nloop=nloop+1
	  call sortcomp(invar,objfun) 		! [SORT ENTIRE ARRAYS]
	!
	!
	!    WRITE BEST_PARAMETERS FILE FOR PREVIOUS LOOP
	!

	  write(*,'("Finished ",i4," main loops --- best objective function =",'// &
		'e12.6,/)') nloop,objfun(1)
		write(2,'("Finished ",i4," main loops --- best objective function =",'// &
		'e12.6,/)') nloop,objfun(1)
	  if(newbest) then
		write(3,outformat) invar(:,1),bestobj
!		write(3,*) invar(:,1),bestobj
	    newbest=.false.
		nsincebest=0
	  else
		write(*,'(/"No improvement in OF for",i5," loops")') nsincebest
		write(2,'(/"No improvement in OF for",i5," loops")') nsincebest
		nsincebest=nsincebest+1
	  endif
	!write(*,'(/"No improvement in OF for",i5,"loops"
	!
	!    ASSESS CONVERGENCE
	!
!	  numcv=(ncomp2*mopt+1)/2 	! USE HALF OF COMPLEX to CALCULATE CV
!	  sumvar=sum(invar(optid,1:numcv),2)
!	  cv=(numcv*sum(invar(optid,1:numcv)**2,2)/sumvar/sumvar-1.)**.5*100.
!
	!		CHANGED BY STAN TO CALCULATE THE FLUCTUATION RANGE RELATIVELY TO
	!		THE FEASIBLE RANGE, INSTEAD OF CV:
		numcv=ncomp2*mopt
		sumvar=sum(invar(optid,1:numcv),2)/numcv		! mean parameter values
		do i=1,nopt
			cv(i)=maxval(abs((invar(optid(i),1:numcv)-sumvar(i))/&
				(parmin(optid(i))-parmax(optid(i))))*100)		! distance from mean in % of feasible range
		enddo
		maxcv=maxval(cv)												! maximum distance
		print*,'Greatest parameter range: ',maxcv,'% for optimised parameter ',parname(optid(maxloc(cv)))
		write(2,'("Greatest parameter range: ",f5.2,"% for optimised parameter ",a9)') maxcv,parname(optid(maxloc(cv)))
	  if(maxcv.ge.resolution.AND.nsincebest.le.patience) then
!	    print*,cv
!	    print*,'Parameter values:'
!	    print*,invar(optid,1)
!			write(2,'("range values: ")') 
!			write(2,outformat) cv
!			write(2,'("parameter values: ")')
!			write(2,outformat) invar(optid,1)
		call writepars(invar(optid,1),parname(optid),parmin(optid),parmax(optid)) 	! PROGRAM STOP

	  elseif(nsincebest.gt.patience) then
			write(*,'(/"No improvement in OF for",i5," loops"/"  About to give up...")') nsincebest
			write(2,'(/"No improvement in OF for",i5," loops"/"  About to give up...")') nsincebest
		  call writepars(invar(optid,1),parname(optid),parmin(optid),parmax(optid)) 
			call optsensitivity(invar,objfun,optid,parname,parmin,parmax,&
					outformat,nrun,success,bestobj,resolution)
			if (success.eq.1) then
				write(*,'(/"Optimisation completed successfully."/)')
				write(2,'(/"Optimisation completed successfully."/)')
				call sortcomp(invar,objfun) 		! [SORT ENTIRE ARRAYS]
			  call writepars(invar(optid,1),parname(optid),parmin(optid),parmax(optid)) 	! PROGRAM STOP
				print*,char(7)
				print*,char(7)
				print*,char(7)				
				return
			else
				newbest=.true. 
				cycle
			endif
	  else
	    write(*,'(/"First Convergence criterion satisfied..."/"  parameter ranges are '// &
			'all less than 0.1 %"//)')
			write(2,'(/"First Convergence criterion satisfied..."/"  parameter ranges are '// &
			'all less than 0.1 %"//)')
		  call writepars(invar(optid,1),parname(optid),parmin(optid),parmax(optid)) 	! PROGRAM STOP
	!
	!		STAN'S MODIFICATION TO ASSESS SENSITIVITY OF OBJECTIVE
	!		FUNCTION TO EACH PARAMETER:
			call optsensitivity(invar,objfun,optid,parname,parmin,parmax,&
					outformat,nrun,success,bestobj,resolution)
			if (success.eq.1) then
				write(*,'(/"Optimisation completed successfully."/)')
				write(2,'(/"Optimisation completed successfully."/)')
				print*,char(7)
				print*,char(7)
				print*,char(7)
				call sortcomp(invar,objfun) 		! [SORT ENTIRE ARRAYS]
			  call writepars(invar(optid,1),parname(optid),parmin(optid),parmax(optid)) 	! PROGRAM STOP
				return
			else
				newbest=.true. 
				cycle
			endif
	  endif
	!
	!
	!    PARTITION THE sopt POINTS INTO ncomp2 COMPLEXES
	!    EXECUTE CCE ALGORITHM
	!
	  if(ncomp2.gt.1) then
	    if(nloop.gt.0.and.ncomp2.gt.ncompmin) then
		ncomp2=ncomp2-1 		! REDUCE NUMBER OF COMPLEXES AS nloop INCREASES
	    elseif(ncomp2.le.ncomp/2.and.worstbest.le.objfun(1+(ncomp2-1)*mopt)) then
		print*,'No gene pool mixing ... reducing number of complexes by one.'
		write(2,'("No gene pool mixing ... reducing number of complexes by one.")')
		ncomp2=ncomp2-1
	    endif
	  endif
	  do m=1,ncomp2
	    first=1+(m-1)*mopt
	    write(*,'("Start of loop",i4,", complex",i2,": best OF =",e12.6)')	&
			nloop+1,m,objfun(first)
	    write(2,'("Start of loop",i4,", complex",i2,": best OF =",e12.6)')	&
			nloop+1,m,objfun(first)
	    if(m.eq.1) then
		bestincomp=-9999. 			! SET LESS THAN bestobj
	    else
		bestincomp=objfun(first)
	    endif
	    call cce(objfun(first:m*mopt),invar(:,first:m*mopt))
	    write(*,'("  End of loop",i4,", complex",i2,": best OF =",e12.6)')	&
			nloop+1,m,objfun(first) 
			write(2,'("  End of loop",i4,", complex",i2,": best OF =",e12.6)')	&
			nloop+1,m,objfun(first) 
	  end do
	  worstbest=objfun(first)

!	ADDED BY STAN TO WRITE INVAR AND OBJFUN OF LAST LOOP TO FILE
		write(str,'(i3)') sopt		!internal write to convert from number to string
		loopformat='('//str//'e34.6)'		! includes a column for each sublayer 
		open(9,file='lastloop.txt')
		write(9,'(i3)') ncomp2
		write(9,loopformat) objfun
!		write(9,*) objfun
		do i=1,npar
			write(9,loopformat) invar(i,:)
!			write(9,*) invar(i,:)
		enddo
		close(9)

	end do
	!
	!
	!  TERMINATE PROGRAM
	!
	write(*,'(/"FAILURE TO CONVERGE..."/"  Number of runs has reached 20000." '// &
	    '//"  Program terminated."/,a,a)') char(7),char(7)
	write(2,'(/"FAILURE TO CONVERGE..."/"  Number of runs has reached 20000." '// &
	    '//"  Program terminated."/,a,a)') char(7),char(7)
	call writepars(invar(optid,1),parname(optid),parmin(optid),parmax(optid)) 		! PROGRAM STOP
	close(8)
	close(3)
	close(2)
	
	contains


	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	subroutine runmodel(invar,optid,objfun,bestobj,bestincomp)
	implicit none
!	INTEGER :: k,system
	INTEGER, DIMENSION(:), INTENT(in) :: optid
	real(double), DIMENSION(:), INTENT(in) :: invar
!	REAL(double), INTENT(out) :: objfun
	real(double) :: objfun
	real(double), INTENT(inout) :: bestobj,bestincomp
	CHARACTER(1) :: bestmark

	nrun=nrun+1
	call transpmodel(invar,nrun,objfun,'-optimise')
	bestmark=' '
	if(objfun.gt.bestobj) then
		bestmark='+'
	elseif(objfun.gt.bestincomp) then
		bestmark='.'
	endif
	if (evolution.ne.'test') then
		write(*,'("Run",i6,"  (",a,"):",t28,"OF =",e12.6,1x,a)') nrun,trim(evolution), &
			objfun,bestmark
		write(2,'("Run",i6,"  (",a,"):",t28,"OF =",e12.6,1x,a)') nrun,trim(evolution), &
			objfun,bestmark
		write(8,outformat) invar(optid),objfun
!		write(8,*) invar(optid),objfun
	endif
!	write(8,rec=nrun+1) objfun,invar(optid)
	return
	end subroutine runmodel


	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	subroutine cce(objfun,invar)
	! INHERITS:  (optid,qopt,mopt,nrun,m,p)
	implicit none
	INTEGER :: i,nsel,rannum,l
	INTEGER, DIMENSION(:), ALLOCATABLE :: parentsid
	real(double), DIMENSION(:), INTENT(inout) :: objfun
	real(double), DIMENSION(:,:), INTENT(inout) :: invar
	real(double), DIMENSION(:), ALLOCATABLE :: objfunsub
	real(double), DIMENSION(:,:), ALLOCATABLE :: invarsub
	LOGICAL, DIMENSION(size(objfun)) :: selected

	!
	!
	! DIMENSIONALISE
	!
	allocate(parentsid(qopt),objfunsub(qopt),invarsub(npar,qopt))
	!
	!
	!  SELECT PARENTS
	!
	do l=1,beta
		selected=.false.
		nsel=0
		do while(nsel.ne.qopt)
			call random_number(ranscal) 			! SCALAR RANDOM NUMBER
			rannum=ceiling((2.*mopt+1.-sqrt(4.*mopt*(mopt+1.)*(1-ranscal)+1))*0.5)
	! NOTE: A SIMPLER ALTERNATIVE TO THE ABOVE LINE IS (19.03.2004):
	!    rannum=ceiling(mopt*(1-sqrt(ranscal)))
			if(rannum.ge.1.and.rannum.le.mopt) then
		if(.not.selected(rannum)) then
			selected(rannum)=.true.
			nsel=nsel+1
		endif
			endif
		end do
		nsel=0
		i=0
		do while(nsel.ne.qopt)
			i=i+1
			if(selected(i)) then
		nsel=nsel+1
		parentsid(nsel)=i
			endif
		end do
	!
	!
	! GENERATE OFFSPRING AND SORT THE RESULTING COMPLEX
	!
		objfunsub=objfun(parentsid)
		invarsub=invar(:,parentsid)
		do j=1,alpha
			call simplex(invarsub,objfunsub)
		end do
		objfun(parentsid)=objfunsub
		invar(:,parentsid)=invarsub
		call sortcomp(invar,objfun)
	end do
	end subroutine cce	


	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	subroutine simplex(invar,objfun)
	implicit none
	INTEGER :: i,j
	real(double), DIMENSION(:), INTENT(inout) :: objfun
	real(double), DIMENSION(:,:), INTENT(inout) :: invar
	real(double), DIMENSION(size(invar,1)) :: centroid,newpoint
	real(double) :: newobjfun,minj,rangej

	!
	!
	! REFLECTION STEP
	!
	evolution='reflection'
	newpoint=invar(:,qopt)
	centroid(optid)=1./(qopt-1)*sum(invar(optid,1:qopt-1),2)
	newpoint(optid)=2.*centroid(optid)-invar(optid,qopt)
	if(minval(newpoint-parmin).lt.0.or.maxval(newpoint-parmax).gt.0) then
	!
	!
	! MUTATION STEP
	!   NB: MUTATION IS BASED ON THE SMALLEST HYPERCUBE OF THE SUBCOMPLEX, NOT THE
	!		SMALLEST HYPERCUBE OF THE COMPLEX AS SUGGESTED BY DUAN ET AL.
	!
		do i=1,nopt
			j=optid(i)
!			if(newpoint(j).lt.parmin(j).or.newpoint(j).gt.parmax(j)) then				!COMMENTED OUT BY STAN TO GENERATE COMPLETELY RANDOM POINT WITHIN HYPERCUBE
		call random_number(ranscal) 		! SCALAR RANDOM NUMBER
!		minj=minval(invar(j,:))
!		rangej=maxval(invar(j,:))-minj
		minj=parmin(j)
		rangej=parmax(j)-minj
		newpoint(j)=minj+rangej*ranscal 	! UNIFORMLY DISTRIBUTED SAMPLE
		evolution='mutation'
!			endif
		end do
	!
	!
	! CALCULATE OBJECTIVE FUNCTION
	!
	endif
	call runmodel(newpoint,optid,newobjfun,bestobj,bestincomp)
	if(newobjfun.gt.objfun(qopt)) then
		do j=1,qopt
			if(newobjfun.gt.objfun(j)) then 		! SORT objfun HERE	
		objfun(j+1:qopt)=objfun(j:qopt-1) 		! IN CASE alpha > 1
		invar(:,j+1:qopt)=invar(:,j:qopt-1)
		objfun(j)=newobjfun
		invar(:,j)=newpoint
		exit
			endif
		end do
	else
	!
	!
	! CONTRACTION STEP
	!
		newpoint(optid)=(centroid(optid)+invar(optid,qopt))*0.5
		evolution='contraction'
		call runmodel(newpoint,optid,newobjfun,bestobj,bestincomp)
		if(newobjfun.gt.objfun(qopt)) then
			do j=1,qopt
		if(newobjfun.gt.objfun(j)) then
			objfun(j+1:qopt)=objfun(j:qopt-1)
			invar(:,j+1:qopt)=invar(:,j:qopt-1)
			objfun(j)=newobjfun
			invar(:,j)=newpoint
			exit
		endif
			end do
		else
	!
	!
	! ANOTHER MUTATION STEP
	!
			do i=1,nopt
		call random_number(ranarr)		! RANDOM ARRAY OF SIZE 1:nopt
		j=optid(i)
		minj=minval(invar(j,:))
		rangej=maxval(invar(j,:))-minj
!		minj=parmin(j)
!		rangej=parmax(j)-minj
		newpoint(j)=minj+rangej*ranarr(i) 	! UNIFORMLY DISTRIBUTED SAMPLE
			end do
			evolution='mutation'
			call runmodel(newpoint,optid,newobjfun,bestobj,bestincomp)
			if(newobjfun.gt.objfun(qopt-1)) then
		do j=1,qopt-1
			if(newobjfun.gt.objfun(j)) then
				objfun(j+1:qopt)=objfun(j:qopt-1)
				invar(:,j+1:qopt)=invar(:,j:qopt-1)
				objfun(j)=newobjfun
				invar(:,j)=newpoint
				exit
			endif
		end do
			elseif(newobjfun.gt.0.0) then				! REPLACE ANYWAY
		objfun(qopt)=newobjfun
		invar(:,qopt)=newpoint
			endif
		endif
	endif
	!
	!
	! UPDATE 'currentbest' IF APPROPRIATE
	!
	if(newobjfun.gt.bestobj) then
		bestobj=newobjfun
		bestincomp=newobjfun
		open(4,file='currentbest.txt')
		write(4,outformat) invar(:,1),bestobj
!		write(4,*) invar(:,1),bestobj
		close(4)
		newbest=.true.
	elseif(newobjfun.gt.bestincomp) then
		bestincomp=newobjfun
	endif
	return
	end subroutine simplex	



	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	subroutine sortcomp(invar,objfun)
	implicit none

	INTEGER :: i,j(1),dim
	real(double), DIMENSION(:),INTENT(inout) :: objfun
	real(double), DIMENSION(:,:),INTENT(inout) :: invar
	real(double), DIMENSION(size(objfun)) :: objfun2
	real(double), DIMENSION(size(invar,1),size(invar,2)) :: invar2
	INTEGER, DIMENSION(size(objfun)) :: newobjfun
	EXTERNAL compar

	dim=size(objfun)
	objfun2=objfun
	newobjfun=-99 			! NEEDED TO SEPARATE EQUAL O.F. VALUES
	invar2=invar
	call qsort(objfun,dim,8,compar) 			! USE compar(b,c)=(c-b)/dabs(c-b)
	do i=1,dim
		j=minloc(real(dabs(objfun2-objfun(i))),newobjfun.lt.0)
		newobjfun(j)=i
		invar(:,i)=invar2(:,j(1))
	end do
	end subroutine sortcomp



	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	subroutine writepars(invar,parname,parmin,parmax)
	implicit none
	INTEGER :: i
	real(double), DIMENSION(:), INTENT(in) :: invar,parmin,parmax
	CHARACTER(9), DIMENSION(:), INTENT(in) :: parname
	write(*,'(/"PARAMETER|     VALUE   |    MINVAL   |    MAXVAL   |  CV (%)     ")')
	write(2,'(/"PARAMETER|     VALUE   |    MINVAL   |    MAXVAL   |  CV (%)     ")')

	do i=1,nopt
		write(*,'(a9,5e14.6)') parname(i),invar(i),parmin(i),parmax(i),cv(i)
		write(2,'(a9,5e14.6)') parname(i),invar(i),parmin(i),parmax(i),cv(i)
	end do
	write(2,'(//)')
	print*
	print*
!	print*,char(7)
!		stop
	end subroutine writepars

	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	! ADDED BY STAN TO ASSESS SENSITIVITY OF OBJECTIVE FUNCTION TO EACH OPTIMISED 
	! PARAMETER:
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	subroutine optsensitivity(invar,objfun,optid,parname,parmin,parmax,&
						outformat,nrun,success,bestobj,okchange)
	implicit none
	real(double), DIMENSION(:),INTENT(inout) :: objfun
	real(double), DIMENSION(:,:),INTENT(inout) :: invar
	real(double), DIMENSION(size(invar,1)) :: invar2
	real(double), INTENT(inout) :: bestobj
	real(double), intent(in) :: okchange
	
	INTEGER, DIMENSION(:), INTENT(in) :: optid
!	REAL(double), DIMENSION(:), INTENT(inout) :: invar
	INTEGER(4) :: i,j,numvar,nrun,numrec,failed10,success,pos
	real(double), DIMENSION(:), INTENT(in) :: parmin,parmax
	CHARACTER(9), DIMENSION(:), INTENT(in) :: parname
	CHARACTER(60) :: outformat
	
	real(double) :: distmin,distmax,oldpar,newpar,distance,objfun2,&
									ofchange,parchange
	real(double), dimension(:,:), allocatable :: dataarray
	
	numvar=size(optid)
	allocate(dataarray(numvar*8+1,numvar+1))
	invar2=invar(:,1)
	objfun2=objfun(1)
	numrec=size(objfun)
	failed10=0
	
	print *,"SENSITIVITY ANALYSIS"
	print *,"(changes in % of feasible range)"
	print *
	write(2,'("SENSITIVITY ANALYSIS"/"changes in % of feasible range)" )')

!	open(9,file='optsensit.txt')
!	close(9,status='delete')
!	open(9,file='optsensit.txt')
!	write(9,outformat) invar2,objfun2
	dataarray(1,1:numvar)=invar2
	dataarray(1,numvar+1)=objfun2

	evolution='test'
	pos=0
	do i=1,nopt
		oldpar=invar2(optid(i))
! TO MAKE SURE THAT PERTURBATIONS DO NOT EXCEED THE FEASIBLE RANGE
! AND THE PARAMETER VALUES THEMSELVES
		distmin=oldpar-parmin(optid(i))
!		if (distmin.gt.abs(oldpar)) then
!			if (abs(oldpar).ge.1.0E-6) then			! 
!				distmin=oldpar
!			endif
!		endif				
		distmax=parmax(optid(i))-oldpar
!		if (distmax.gt.abs(oldpar)) then
!			if (abs(oldpar).ge.1.0E-6) then
!				distmax=oldpar
!			endif
!		endif

		write(*,'(/"change of var: ",a9, "(",e9.3,")")') parname(optid(i)),&
					oldpar
		write(2,'(/"change of var: ",a9, "(",e9.3,")")') parname(optid(i)),&
					oldpar

		do j=0,3
 			newpar=oldpar-distmin*10.0**(-j)
			invar2(optid(i))=newpar
			nrun=nrun+1
			call transpmodel(invar2,nrun,objfun2,'-optimise')
			if (objfun2.gt.bestobj) then
				bestobj=objfun2
				open(4,file='currentbest.txt')
				write(4,outformat) invar2,bestobj
!				write(4,*) invar2,bestobj
				close(4)
			endif	
			write(8,outformat) invar2(optid),objfun2
!			write(8,*) invar2(optid),objfun2
			dataarray((i-1)*8+j+2,1:numvar)=invar2
			dataarray((i-1)*8+j+2,numvar+1)=objfun2
			ofchange=(objfun2-dataarray(1,numvar+1))/abs(dataarray(1,numvar+&
								1))*100.0
			parchange=(newpar-oldpar)/(parmax(optid(i))-parmin(optid(i)))*100		! the change of the parameter in % of feasible range
			write(*,'(f6.3,"% (",f14.6,")",": change of OF by ",e9.3,"%"," (",&
						e9.3,")")') parchange,newpar,ofchange,objfun2
			write(2,'(f6.3,"% (",f14.6,")",": change of OF by ",e9.3,"%"," (",&
						e9.3,")")') parchange,newpar,ofchange,objfun2

			if (ofchange.gt.1.0E-10) then
				invar(:,numrec-pos)=invar2
				objfun(numrec-pos)=objfun2
				pos=pos+1
				if (abs(parchange).gt.okchange) then
					failed10=failed10+1
				endif
			endif
		enddo

		do j=3,0,-1
			newpar=oldpar+distmax*10.0**(-j)
			invar2(optid(i))=newpar
			nrun=nrun+1
			call transpmodel(invar2,nrun,objfun2,'-optimise')
			if (objfun2.gt.bestobj) then
				bestobj=objfun2
				open(4,file='currentbest.txt')
				write(4,outformat) invar2,bestobj
!				write(4,*) invar2,bestobj
				close(4)
			endif		
			write(8 ,outformat) invar2(optid),objfun2
!			write(8 ,*) invar2(optid),objfun2
			dataarray((i-1)*8+j+6,1:numvar)=invar2
			dataarray((i-1)*8+j+6,numvar+1)=objfun2
			ofchange=(objfun2-dataarray(1,numvar+1))/abs(dataarray(1,numvar+&
								1))*100.0
			parchange=(newpar-oldpar)/(parmax(optid(i))-parmin(optid(i)))*100		! the change of the parameter in % of feasible range
			write(*,'(f6.3,"% (",f14.6,")",": change of OF by ",e9.3,"%"," (",&
						e9.3,")")') parchange,newpar,ofchange,objfun2
			write(2,'(f6.3,"% (",f14.6,")",": change of OF by ",e9.3,"%"," (",&
						e9.3,")")') parchange,newpar,ofchange,objfun2

			if (ofchange.gt.1.0E-10) then
				invar(:,numrec-pos)=invar2
				objfun(numrec-pos)=objfun2
				pos=pos+1
				if (abs(parchange).gt.okchange) then
					failed10=failed10+1
				endif
			endif
		enddo
		invar2(optid(i))=oldpar
	enddo
	close(9)

	! ASSESS SECOND CONVERGENCE CRITERIUM: NO PARAMETER CHANGE BY MORE THAN 10% LEADS
	! TO AN INCREASE IN OBJECTIVE FUNCTION
	if (failed10.gt.0) then
		write(*,'(i2," parameter(s) more than ",f6.3,"% out of optimum."/"Optim&
			isation continued...")') failed10,okchange
		write(2,'(i2," parameter(s) more than ",f6.3,"% out of optimum."/"Optim&
			isation continued...")') failed10,okchange
	else
		write(*,'(//"Second convergence criterion satisfied:"/"no parameter &
			shift by more than ",f6.3,"% of max distance leads to"/"an increase in &
			the objective function.")') okchange
		write(2,'(//"Second convergence criterion satisfied:"/"no parameter &
			shift by more than ",f6.3,"% of max distance leads to"/"an increase in &
			the objective function.")') okchange
		success=1
	endif
	
	! IF CRITERIUM NOT SATISFIED, APPEND NEWLY CREATED DATA TO THE END OF INVAR AND
	! OBJFUN ARRAYS AND RETURN TO MAIN LOOP TO RE-RUN SCE-LOOP
!	if (success.eq.0) then
!		invar(:,numrec-8*numvar:numrec)=transpose(dataarray(2:numvar*8+1,1:numvar))
!		objfun(numrec-8*numvar:numrec)=dataarray(2:numvar*8+1,numvar+1)
!	endif
		
	deallocate(dataarray)
	end subroutine optsensitivity

	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	end



	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	integer function compar(b,c)
	implicit none
	INTEGER, PARAMETER :: double=SELECTED_REAL_KIND(13) 	! == REAL*8
	real(double), INTENT(in) :: b,c
	if(b.gt.c) then 		! SORTS b,c IN DECREASING ORDER
	  compar=-1
	elseif(b.lt.c) then
	  compar=1
	else
	  compar=0
	endif
	return
	end function compar




