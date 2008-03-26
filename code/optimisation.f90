!********************************************************************
!     Optimised Vegetation Catchment Feedback Model (ovecafee)
!     Core program to run optimisation (sce) and transpmodel
!--------------------------------------------------------------------
!     Author: Stan Schymanski, CWR, University of Western Australia
!     05/05/2004
!     FOR COMPILATION WITH COMPAQ VISUAL FORTRAN
!  
!  Copyright (C) 2004  Stan Schymanski
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
!
!********************************************************************
program ovecafee
use dfport
use avdef
use avviewer

implicit none

character(60) outformat
character(15) command
character(3) str
INTEGER, PARAMETER :: double=SELECTED_REAL_KIND(13) 	! == REAL*8
real(double), allocatable :: invar(:)
real(double), DIMENSION(1) :: netass
integer nrun,success,npar,ios,i
	
!call GETARG(1, option1)
!--------------------------------------------------------------------
!	for debug purposes:
!	option1='-optimise'
!--------------------------------------------------------------------		
open(1,file='shuffle.par',status='old')
read(1,*) command
	
if (command.ne.'-compute       ') then
	rewind(1)
	close(1)
	call sce(success)
	print *,"Optimisation of vegetation and model run COMPLETE"
else

!--------------------------------------------------------------------
!	model run with optimised parameters
!--------------------------------------------------------------------
	nrun=1
	call showsce(success)
	print *,"Calculation of results with optimised parameters..."
	open(1,file='shuffle.par',status='old')
	read(1,*)
	read(1,*)
	npar=0
	do
		read(1,*,iostat=ios) str
		if(ios.lt.0) exit
		if(str.eq.'var') npar=npar+1
	end do
	close(1)
	allocate(invar(npar))
	write(str,'(i3)') npar+1		!internal write to convert from number to string
	outformat='('//str//'e14.6)'		! includes a column for each parameter and a column for the value of OF
	open(2,file='currentbest.txt')		! reads input parameters from previous optimisation
	rewind(4)
!	read(2,outformat)invar(:),netass
	read(2,*) invar(:),netass
	close(2)
	write(*,'(" The best carbon profit was: ",e12.6)') netass
	call transpmodel(invar,nrun,netass,command)
	write(*,'(/" Model run COMPLETE",/)')
	write(*,'(" The carbon profit achieved is: ",e12.6)') netass
!	nrun=2
!	command='-optimise'
!		call transpmodel(invar,nrun,netass,command)
!	write(*,'(/" Model run COMPLETE",/)')
!	write(*,'(" The carbon profit achieved is: ",e12.6)') netass
	print *, "Hourly results are saved in etasslin.txt"
	print *, "Daily results are saved in etassdaily.txt"
	print *, "Yearly results are saved in yearly.txt"

endif


end
