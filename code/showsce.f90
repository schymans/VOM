	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	!
	! 	SHOWSCE:
	!   ROUTINE TO SHOW EVOLUTION OF PARAMETER VALUES AND OBJECTIVE FUNCTION
	! 	Written by: Stan Schymanski, Centre for Water Research, UWA
	!
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	!
	!  Compile with: Compaq Visual Fortran 5.5	
	!
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	!
	!  Copyright (C) 2005 Stan Schymanski
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

subroutine showsce(success)
	use avdef
	use avviewer
	use dflib

	implicit none

	INTEGER, PARAMETER :: double=SELECTED_REAL_KIND(13) 	! == REAL*8
	integer, parameter :: XAxis=0,YAxis=1,ZAxis=2
	integer :: nloop,npar,ios,nrec,status,success,hv,onoroff,XComp&
							,YComp,ZComp,WComp,length
	integer(2) :: i,j,nopt
	character(60) :: informat, outformat
	character(3) :: str
	character(1) :: key
	character(9), dimension(:), allocatable :: parname,parname1
	real(double), dimension(:,:), allocatable ::	parval
	real(double), dimension(:), allocatable :: vector,vector1,vector3
	real(double) :: percentile,low,high,blank,blank1,blank2,blank3
	integer(2), dimension(:), allocatable :: paropt,optid

	print *, "Do you want to plot optimisation process? (y/n)"
	key = GETCHARQQ()
	if (key /= 'y') then
		return
	endif

	!==============================================
	! OPEN FILES
	!==============================================
	write(*,'("Opening parameter file...")')
	open(1,file='shuffle.par',status='old')
	open(8, file='sce.out',status='old')
	rewind(8)



	!==============================================
	! READ RECORD LENGTHS AND DEFINE FORMAT SPECS
	!==============================================
	read(1,*)
	read(1,*)
	npar=0
	do
		read(1,*,iostat=ios) str
		if(ios.lt.0) exit
		if(str.eq.'var') npar=npar+1
	end do
		
	allocate(paropt(npar),parname1(npar))
	write(str,'(i3)') npar+1		!internal write to convert from number to string
	outformat='('//str//'f14.6)'		! includes a column for each parameter and a column for the value of OF

	rewind(1)
	read(1,*) 
	read(1,*) 
	read(1,*)
	read(1,*)
	read(1,*)
	read(1,*)
	read(1,*)
	read(1,'(a60)') informat
	do i=1,npar
	  read(1,informat) parname1(i),blank1,blank2,&
		blank3,paropt(i)
	end do
	close(1)
	nopt=sum(paropt)
	allocate(parname(nopt+1),optid(nopt),vector3(npar))

	! FILTER NAMES OF OPTIMISED PARAMETERS ONLY
	
	j=0
	do i=1,npar
	  if(paropt(i).gt.0) then
	    j=j+1
	    optid(j)=i
	  endif
	end do
	parname(1:nopt)=parname1(optid)
	parname(nopt+1)='OF'			! OF=objective function


	! READ NUMBER OF RECORDS IN SCE.OUT
	nrec=0
	do
		read(8,*,iostat=ios) vector3
		if (ios.lt.0) exit
		nrec=nrec+1
	enddo
	rewind(8)


	!==============================================
	! ALLOCATE ARRAY LENGTHS 
	!==============================================
	allocate(parval(nrec,nopt+1),vector(nrec))

	do i=1,nrec
		read(8,*) parval(i,:)
	enddo
	close(8)

	!==============================================
	! LOAD DATA INTO ARRAY VIEWER
	!==============================================
	call faglStartWatch(parval,status)
	!call faglShow(parval,status)
	!call faglEndWatch(parval,status)

	!==============================================
	! DETERMINE VARIABLES FOR DISPLAY
	!==============================================
	onoroff=1
	call favStartViewer (hv, status)
	call favSetArray (hv, parval, status)
	call favSetGraphType (hv, 3, status)		! change integer number to change graph type
	
	zComp=nopt+1
	call favSetAxisLabel (hv, zaxis, parname(zComp), status)
	percentile=90.0
	vector=parval(:,zComp)
	length=size(vector)
	call quantiles(vector,length,percentile,low,high)		
	call favSetZClamp (hv, low, high, status)
	call favSetAxisAutoScale (hv, 0, status)
	call favUpdate (hv, 0, status)

	do
		
		write(*,'("Variable names and their codes:/")')
		do i=1,nopt
			write(*,'("name: ",a," code: ",i3)') parname(i),i
		enddo
		write(*,'(/"Choose variable for x-axis: ")')
		read(*,'(i2)') xComp
		call favUpdate (hv, 0, status)

		call favSetAxisLabel (hv, xaxis, parname(xComp), status)
		write(*,'("Choose variable for y-axis: ")')
		read(*,'(i2)') yComp
		call favUpdate (hv, 0, status)

		call favSetAxisLabel (hv, yaxis, parname(yComp), status)
		write(*,'("Choose variable for colour: ")')
		read(*,'(i2)') wComp
		call favUpdate (hv, 0, status)

		call favSetCompIndex (hv, XComp, YComp, ZComp, WComp, status)
	
		!==============================================
		! SETTING AXES RANGES
		!==============================================
		percentile=90.0
		vector=parval(:,xComp)
		length=size(vector)
		call quantiles(vector,length,percentile,low,high)		
		call favSetXClamp (hv, low, high, status)

		vector=parval(:,yComp)
		length=size(vector)
		call quantiles(vector,length,percentile,low,high)		
		call favSetYClamp (hv, low, high, status)

!		vector=parval(:,zComp)
!		length=size(vector)
!		call quantiles(vector,length,percentile,low,high)		
!		call favSetzClamp (hv, low, high, status)
		
!		call favSetAxisAutoScale (hv, 0, status)
		call favShowWindow (hv, onoroff, status)
		call favUpdate (hv, 0, status)

		print *, "Done!  Plot other variables? (y/n)"

		key = GETCHARQQ()
		if (key /= 'y') then
			print *,""
			print *, "Array Viewer window will be closed."
			print *,""
			print *, "SAVE BEFORE PROCEEDING!"
			print *,""
			print *, "Press any key when ready"
			key = GETCHARQQ()
			call faglEndWatch(parval,status)
			call favEndViewer (hv, status)

			exit  ! break out of the outer do loop
		end if
	
	end do




	return
end
