subroutine plotdata(dataarray,npar,plotstyle,parnames)
	use avdef
	use avviewer
	use dflib
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	!
	! 	ROUTINE TO PLOT DATA WITH ARRAY VISUALISER
	! 	Written by: Stan Schymanski, Centre for Water Research, UWA
	! 	To be used with scestan.f90 
	!
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	implicit none

	integer, parameter :: XAxis=0,YAxis=1,ZAxis=2
	integer :: npar,ios,nrec,status,success,hv,onoroff,XComp&
							,YComp,ZComp,WComp,length
	integer(2) :: i
	character(20) :: outformat
	character(3) :: str
	character(1) :: key
	character(9), dimension(:) :: parname
	real(8), dimension(:,:), allocatable ::	dataarray
	real(8), dimension(:), allocatable :: test,vector1
	real(8) :: percentile,low,high


	!==============================================
	! OPEN FILES
	!==============================================
	write(*,'("Opening parameter file...")')
	open(1,file='shuffle.par',status='old')
	open(8, file='sce.out',status='old')


	!==============================================
	! ALLOCATE ARRAY LENGTHS 
	!==============================================
	allocate(parval(nrec,npar+1),vector1(nrec))
	!FORMAT: ONE COLUMN PER PAREMETER AND ONE FOR OF



	!==============================================
	! READ VALUES
	!==============================================
	read(1,*)
	read(1,*)
	read(1,*)
	do i=1,npar
		read(1,'(7x,a9)') parname(i)
	enddo
	parname(npar+1)='OF'			! OF=objective function
	close(1)
	do i=1,nrec
		read(8,outformat) parval(i,:)
	enddo
	close(8)

	!==============================================
	! LOAD DATA INTO ARRAY VIEWER
	!==============================================
	call faglStartWatch(array,status)
	!call faglShow(parval,status)
	!call faglEndWatch(parval,status)

	!==============================================
	! DETERMINE VARIABLES FOR DISPLAY
	!==============================================
	onoroff=1
	call favStartViewer (hv, status)
	call favSetArray (hv, parval, status)
	call favSetGraphType (hv, 3, status)		! change integer number to change graph type
	call favGetCompIndex (hv, XComp, YComp, ZComp, WComp, status)
	call favSetCompIndex (hv, XComp, YComp, ZComp, WComp, status)
!	call favShowWindow (hv, onoroff, status)
	
	zComp=8
	call favSetAxisLabel (hv, zaxis, parname(zComp), status)
	percentile=90.0
	vector1=parval(:,zComp)
	length=size(vector1)
	call quantiles(vector1,length,percentile,low,high)		
	call favSetZClamp (hv, low, high, status)
	call favSetAxisAutoScale (hv, 0, status)
	call favUpdate (hv, 0, status)

	do
		
		write(*,'("Variable names and their codes:/")')
		do i=1,npar
			write(*,'("name: ",a,", code: ",i3)') parname(i),i
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

		!==============================================
		! SETTING AXES RANGES
		!==============================================
		percentile=90.0
		vector1=parval(:,xComp)
		length=size(vector1)
		call quantiles(vector1,length,percentile,low,high)		
		call favSetXClamp (hv, low, high, status)

		vector1=parval(:,yComp)
		length=size(vector1)
		call quantiles(vector1,length,percentile,low,high)		
		call favSetYClamp (hv, low, high, status)

!		vector1=parval(:,zComp)
!		length=size(vector1)
!		call quantiles(vector1,length,percentile,low,high)		
!		call favSetzClamp (hv, low, high, status)
		
!		call favSetAxisAutoScale (hv, 0, status)

		call favUpdate (hv, 0, status)

		print *, "Done!  Plot other variables? (y/n)"

		key = GETCHARQQ()
		if (key /= 'y') then
			exit  ! break out of the outer do loop
		end if
	
	end do


	return
end
