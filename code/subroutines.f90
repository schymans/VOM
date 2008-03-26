!==============================================
! QUANTILES:
! SUBROUTINE TO CALCULATE THE QUANTILES OF 
! A VECTOR
  !   Written by: Stan Schymanski, Centre for Water Research, UWA
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
!==============================================
! needs the length of the vector to determine allocation
! percentile is the percent of data to be included

subroutine quantiles(vector,length,percentile,low,high)		
	implicit none
	INTEGER, PARAMETER :: double=SELECTED_REAL_KIND(13) 	! == REAL*8
  integer :: length
	real(double), dimension(length),intent(in) :: vector
	real(double) :: percentile,low,high,cut 				
	integer(2), external :: cmp_function

	call qsort(vector,length,8,cmp_function) 
	cut=floor(0.5*length*(1-percentile/100))
	high=vector(cut)
	low=vector(length-cut)
end subroutine

! function to use for qsort
integer function cmp_function(b,c)
	implicit none
	INTEGER, PARAMETER :: double=SELECTED_REAL_KIND(13) 	! == REAL*8
	REAL(double), INTENT(in) :: b,c
	if(b.gt.c) then 		! SORTS b,c IN DECREASING ORDER
	  cmp_function=-1
	elseif(b.lt.c) then
	  cmp_function=1
	else
	  cmp_function=0
	endif
	return
end function cmp_function
