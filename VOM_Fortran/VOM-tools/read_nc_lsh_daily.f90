!********************************************************************
!	First version written by: Kerstin Sickel, MPI for Biogeochemistry, Jena, Germany
!
!       Modified and extended by Andreas Ostrowski, MPI for Biogeochemistry, Jena, Germany
!	
!       First release on 04/03/2009
!
!       aostrow@bgc-jena.mpg.de
!
!--------------------------------------------------------------------
!
!    Copyright (C) 2009 Andreas Ostrowski
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
!********************************************************************
program read_nc_lsh_daily

!=================================
! PURPOSE OF THIS FORTRAN PROGRAM:
!=================================
!
! creating of a dailyweather text file for use in VOM
!
! using lsh dataset (netCDF files), please look at README.txt inside svn folder for more information
!
! reads out variables for a specific region (point), longitude and latidude is to pass
!
! following variables / files are needed:
        ! tmax    - maximum temperature
        ! tmin    - minimum temperature
        ! prec    - precipitation
        ! shumin  - specific humidity
        ! pres    - pressure
        ! dswrf   - downward shortwave radiation
!
! output is a textfile with the following parameters:
        ! doy
        ! day
        ! month
        ! year
        ! tmax
        ! tmin
        ! prec
        ! evap
        ! rg ("=dswrf")
        ! vp (vapour pressure)
        ! pres (air pressure)
        ! ca (ratio of CO2 in air)
!
! for compiling use:
!
! pgf90 -o ncgetlshd read_nc_lsh_daily.f90 -I/usr/local/netcdf-3.5/include -L/usr/local/netcdf-3.5/lib -lnetcdf
!
!
! for extracting:
!
! execute created file (./ncgetlshd) with two parameters (longitude and latidude)
!
! the original dataset is in 1 degree resolution
! therefore the used parameters will be rounded, but you have to use decimal degrees !
! negative coordinates (south or west) are substracted by 1 for extraction of data for the right cell
! example: if you are using "S -0.4" originally the rounding routine would get "0", like also for "N 0.4"
! but "S -0.4" means the cell below, so we have to substract 1 to get the right cell
!
! you have to adjust your path to the data directory (or not) and your required years below !
!
!=================================

      IMPLICIT NONE

!=================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! adjust the path to the data directory !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! if it is the same as this fortran file: clear the path string but keep the inverted commas without space -> ''

      CHARACTER(*), PARAMETER :: datadir = '/Net/Groups/C-Side/BTM/data/lsh/downloadMoreDaily/'

!=================================     

      INTEGER, PARAMETER :: jpi = 360, jpj = 180, jit = 366      
      INTEGER :: ji, jj, jk, x, y
      INTEGER :: time, year, month, day

      CHARACTER(80) :: nameou
      CHARACTER(80) :: filein1,filein2,filein3
      CHARACTER(80) :: filein4,filein5,filein6,filein7
      
      REAL :: ncep1s(jpi,jpj,366)
      REAL :: ncep2s(jpi,jpj,366)
      REAL :: ncep3s(jpi,jpj,366)
      REAL :: ncep4s(jpi,jpj,366)
      REAL :: ncep5s(jpi,jpj,366)
      REAL :: ncep6s(jpi,jpj,366)
      REAL :: ncep7s(jpi,jpj,366)
       
      REAL :: tmax(jpi,jpj,365)
      REAL :: tmin(jpi,jpj,365)
      REAL :: prec(jpi,jpj,365)
      REAL, PARAMETER :: evap = -99.
      REAL, PARAMETER :: ca = 350.
      REAL :: shum(jpi,jpj,365)
      REAL :: pres(jpi,jpj,365)
      REAL :: dswrf(jpi,jpj,365)
      REAL :: dlwrf(jpi,jpj,365)
      REAL :: vp(jpi,jpj,365)
      
      REAL :: cell_lon, cell_lat, u_lon, u_lat
      INTEGER :: icell_lon, icell_lat, ilon, ilat
      REAL :: lat1(jpi), lon1(jpj), time1(jit)

      INCLUDE 'netcdf.inc'

      INTEGER :: ncidin1, ncidin2, ncidin3
      INTEGER :: ncidin4, ncidin5, ncidin6, ncidin7      
      INTEGER :: varlat1, varlon1, vartime1
      INTEGER :: varid1, varid2, varid3
      INTEGER :: varid4, varid5, varid6, varid7
      INTEGER :: status

!=================================
!  parameter query (longitude and latitude)
!=================================

      INTEGER :: argc              ! arg count
      CHARACTER(30) :: argv        ! arg values (make as big as you need)

      write(*,*) "---------------"
      write(*,*) "PROGRAM STARTED"
      write(*,*) "---------------"

      write(*,*) "TAKING PARAMTERS"
      write(*,*) "dataset is in 1 degree resolution"
      write(*,*) "your parameters will be rounded assuming they are given in decimal degrees"
      write(*,*) "negative parameters (south and west coordinates) are subtracted by 1 for calculation"

      argc = COMMAND_ARGUMENT_COUNT() 
      call getarg(0,argv)

      ! number of parameters checking
      if (argc .ne. 2) then
         write(*,*) 'to few arguments'
         stop
      end if
      
      ! taking parameters
      call getarg(1,argv)
      read (argv,*)cell_lon   
      call getarg(2,argv)
      read (argv,*)cell_lat   
      
      ! dimension of parameters checking
      if ((cell_lon.ge.180.).OR.(cell_lon.le.-180.)) then
          write(*,*) 'longitude outside dimension'
          stop 
      end if
      if ((cell_lat.ge.90.).OR.(cell_lat.le.-90.)) then
          write(*,*) 'latitude outside dimension'
          stop 
      end if

      ! parameter border adjusting
      ! longitude
      if (cell_lon.lt.0.) then
         u_lon = cell_lon - 1.
      else
         u_lon = cell_lon
      end if
      
      ! latitude
      if (cell_lat.lt.0.) then
         u_lat = cell_lat - 1.
      else
         u_lat = cell_lat
      end if

      ! just used for creating output file (name)
      icell_lon = int(u_lon)
      icell_lat = int(u_lat)

      write(*,*) "PARAMETERS TAKEN - IF NO MESSAGE APPEARED: ALL FINE"
      write(*,*) "USING:", icell_lon, icell_lat

      ! output for testing purpose
!      write(*,*) "lon=", cell_lon
!      write(*,*) "lat=", cell_lat

!      write(*,*) "int_lon=", icell_lon
!      write(*,*) "int_lat=", icell_lat

!      write(*,*) "u_lon=", u_lon
!      write(*,*) "u_lat=", u_lat

      
!=================================
!  open the netCDF file
!=================================
      
      write(*,*) "CREATING OUTPUT FILE AND OPEN FOR WRITING" 

      ! creating of output file
      ! after creation you have to rename it to dailyweather.prn for use in VOM
      write(nameou,'(a19,SP,i4.3,a1,SP,i3.2,a4)')   &
           & "lsh_output_30y_pres",icell_lon,"_",icell_lat,".txt"
      open(10,file=nameou,status='unknown')
      write(10,'(12a8)')'    Dcum','     Day','   Month','    Year','   T.Max','   T.Min','    Rain','    Evap','    Radn','      VP','    Pres','      Ca'

      !  print*,nameou      
      write(*,*) "OUTPUT FILE CREATED AND OPENED", nameou

!=================================           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! years to proceed, adjust to your needs !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do year=1971,2000
!=================================     
        
         ! The function NF_OPEN opens an existing netCDF dataset for access.

         write(*,*) "-----------------------------------"
         write(*,*) "START OPENING NETCDF FILES FOR YEAR", year
   
         write(filein1,'(a11,i4,a1,i4,a3)')   &
              & "tmax_daily_",year,"-",year,".nc"
         status = nf_open(datadir // filein1, 0, ncidin1)
         if (status.ne.0) write(*,*) "Problem by opening ",filein1,"Message:", status
    
         write(filein2,'(a11,i4,a1,i4,a3)')   &
              & "tmin_daily_",year,"-",year,".nc"
         status = nf_open(datadir // filein2, 0, ncidin2)
         if (status.ne.0) write(*,*) "Problem by opening ",filein2,"Message:", status
       
         write(filein3,'(a11,i4,a1,i4,a3)')   &
              & "prcp_daily_",year,"-",year,".nc"
         status = nf_open(datadir // filein3, 0, ncidin3)
         if (status.ne.0) write(*,*) "Problem by opening ",filein3,"Message:", status
       
         write(filein4,'(a11,i4,a1,i4,a3)')   &
              & "shum_daily_",year,"-",year,".nc"
         status = nf_open(datadir // filein4, 0, ncidin4)
         if (status.ne.0) write(*,*) "Problem by opening ",filein4,"Message:", status
       
         write(filein5,'(a11,i4,a1,i4,a3)')   &
              & "pres_daily_",year,"-",year,".nc"
         status = nf_open(datadir // filein5, 0, ncidin5)
         if (status.ne.0) write(*,*) "Problem by opening ",filein5,"Message:", status
       
         write(filein6,'(a12,i4,a1,i4,a3)')   &
              & "dswrf_daily_",year,"-",year,".nc"
         status = nf_open(datadir // filein6, 0, ncidin6)
         if (status.ne.0) write(*,*) "Problem by opening ",filein6,"Message:", status

         write(filein7,'(a12,i4,a1,i4,a3)')   &
              & "dlwrf_daily_",year,"-",year,".nc"
         status = nf_open(datadir // filein7, 0, ncidin7)
         if (status.ne.0) write(*,*) "Problem by opening ",filein7,"Message:", status    
      
         write(*,*) "OPENING DONE - IF NO MESSAGE APPEARED: ALL FINE"

!=================================
!  Inquire attributes for the grid and variables & checking if correct
!=================================

         ! The function NF_INQ_VARID returns the ID of a netCDF variable, given its name.

         write(*,*) "START INQUIRING NETCDF FILES"

         status = nf_inq_varid(ncidin1, 'longitude', varlon1)
         if (status.ne.0) write(*,*) "Problem in varlon1",status
         status = nf_inq_varid(ncidin1, 'latitude', varlat1)
         if (status.ne.0) write(*,*) "Problem in varlat1",status
         status = nf_inq_varid(ncidin1, 'time', vartime1)
         if (status.ne.0) write(*,*) "Problem in vartime1",status

         status = nf_inq_varid(ncidin1, 'tmax', varid1)
         if (status.ne.0) write(*,*) "Problem in varid1",status
       
         status = nf_inq_varid(ncidin2, 'tmin', varid2)
         if (status.ne.0) write(*,*) "Problem in varid2",status
       
         status = nf_inq_varid(ncidin3, 'prcp', varid3)
         if (status.ne.0) write(*,*) "Problem in varid3",status
       
         status = nf_inq_varid(ncidin4, 'shum', varid4)
         if (status.ne.0) write(*,*) "Problem in varid4",status
       
         status = nf_inq_varid(ncidin5, 'pres', varid5)
         if (status.ne.0) write(*,*) "Problem in varid5",status
       
         status = nf_inq_varid(ncidin6, 'dswrf', varid6)
         if (status.ne.0) write(*,*) "Problem in varid6",status

         status = nf_inq_varid(ncidin7, 'dlwrf', varid7)
         if (status.ne.0) write(*,*) "Problem in varid7",status
 
         write(*,*) "INQUIRING DONE - IF NO MESSAGE APPEARED: ALL FINE"  

!=================================
!  Read the data
!=================================

         ! The members of the NF_GET_VAR_ type family of functions read all the values from a netCDF variable of an open netCDF dataset

         write(*,*) "START READING NETCDF FILES"

         ! not needed ?
         status = nf_get_var_real(ncidin1, varlon1, lon1)
         status = nf_get_var_real(ncidin1, varlat1, lat1)
         status = nf_get_var_real(ncidin1, vartime1, time1)
         
         ! needed !
         status = nf_get_var_real(ncidin1, varid1, ncep1s)              
         status = nf_get_var_real(ncidin2, varid2, ncep2s)              
         status = nf_get_var_real(ncidin3, varid3, ncep3s)              
         status = nf_get_var_real(ncidin4, varid4, ncep4s)              
         status = nf_get_var_real(ncidin5, varid5, ncep5s)              
         status = nf_get_var_real(ncidin6, varid6, ncep6s) 
         status = nf_get_var_real(ncidin7, varid7, ncep7s)  
         
         write(*,*) "READING NETCDF FILES DONE" 
   
!=================================
!  close the netCDF file
!=================================

         ! The function NF_CLOSE closes an open netCDF dataset.
         
         status = nf_close(ncidin1)
         status = nf_close(ncidin2)
         status = nf_close(ncidin3)
         status = nf_close(ncidin4)
         status = nf_close(ncidin5)
         status = nf_close(ncidin6)
         status = nf_close(ncidin7)

         write(*,*) "NETCDF FILES CLOSED - START CALCULATING"
         
!=======================================================================

         ! description of variables and the used conversion

         ! temperature in nc files: K
         ! needed for output: °C 
         ! therefore: -273.15

         ! precipitation in nc files: kg/m2/s
         ! needed for output: mm/d
         ! therefore: *86400 (for time, kg=l=mm)

         ! evap is not available and also not needed for the VOM model
         ! therefore: dummy with -99., already assigned in variable definition above

         ! specific humidity in nc files: kg/kg
         ! needed for calculation of vp
         ! here g/kg is used, therefore *1000 in vp formula, conversion not used anymore, see vp description below

         ! pressure in nc files: Pa
         ! needed for calculation of vp
         ! here hPa is used, therefore /100 in vp formula, conversion not used anymore, see vp description below

         ! downward shortwave radiation in nc files: W(/m2)
         ! needed for output: MJ/d(/m2)
         ! therefore: *86400/1000000

         ! ca is not available and at the moment we use a constant value for the VOM model
         ! therefore: dummy with 350. (ymol), already assigned in variable definition above

	 do ji = 1, 180
              x = 180 + ji	
              do jj = 1, jpj
                 y = 181 - jj
                 do jk = 1, 365
                    tmax(x, jj, jk)  = ncep1s(ji, y, jk) - 273.15
                    tmin(x, jj, jk)  = ncep2s(ji, y, jk) - 273.15
                    prec(x, jj, jk)  = ncep3s(ji, y, jk) * 86400.
                  !  evap = -99.
                    shum(x, jj, jk)  = ncep4s(ji, y, jk)
                    pres(x, jj, jk)  = ncep5s(ji, y, jk)
                    dswrf(x, jj, jk) = (ncep6s(ji, y, jk) * 86400.) / 1000000.
                    dlwrf(x, jj, jk) = (ncep7s(ji, y, jk) * 86400.) / 1000000.	  	          
                 enddo
              enddo
           enddo
           
           do ji = 181, 360
              x = ji - 180	 
              do jj = 1, jpj
                 y = 181 - jj
                 do jk = 1, 365	  	          
                    tmax(x, jj, jk)  = ncep1s(ji ,y, jk) - 273.15
                    tmin(x, jj, jk)  = ncep2s(ji ,y, jk) - 273.15
                    prec(x, jj, jk)  = ncep3s(ji ,y, jk) * 86400.
                    ! evap = -99.
                    shum(x ,jj, jk)  = ncep4s(ji, y, jk)
                    pres(x ,jj, jk)  = ncep5s(ji, y, jk)
                    dswrf(x ,jj, jk) = (ncep6s(ji, y, jk) * 86400.) / 1000000.
                    dlwrf(x ,jj, jk) = (ncep7s(ji, y, jk) * 86400.) / 1000000.
                 enddo
              enddo
           enddo
           
           do ji = 1, jpi
              do jj = 1, jpj
                 do jk = 1, 365
                    
                    ! actually the formula was inteded to use with fitted units inside (as described above), but the outcoming results make no sense ?!
                    	  	          
                    !	       vp(ji,jj,jk) = ((pres(ji,jj,jk)/100.) * (shum(ji,jj,jk)*1000.)) /  &
                    !     & 	   		      (0.378*(shum(ji,jj,jk)*1000.) + 0.622)  
                    
                    ! therefore the following formula is used:
                    ! caclulation is done by using Pa for pressure variables and kg/kg for specific humidity (as the original formula is doing)
                    ! the overall term is divided by 100 to get vp in hPa
                    
                    vp(ji, jj, jk) = (((pres(ji, jj, jk)) * (shum(ji, jj, jk))) /  &
                         & 	   		      (0.378 * (shum(ji, jj, jk)) + 0.622)) / 100.  
                                        
                 enddo
              enddo
           enddo
           
           write(*,*) "CALCULATING DONE - START WRITING OUTPUT" 
           
!=================================
!WRITING THE OUTPUT
!=================================
       
           ilon = 181 + int(u_lon)
           ilat =  90 - int(u_lat)
      
! output for testing purpose
!           write(*,*) "ilon=", ilon
!           write(*,*) "ilat=", ilat
           
           do time=1,365   
              if (time.le. 31 ) then
                 month = 1
                 day   = time
              endif
              if (time .gt. 31 .and. time .le. 59) then
                 month = 2
                 day   = time - 31
              endif
              if (time .gt. 59 .and. time .le. 90) then
                 month = 3
                 day   = time - 59
              endif
              if (time .gt. 90 .and. time .le. 120) then
                 month = 4
                 day   = time - 90
              endif
              if (time .gt. 120 .and. time .le. 151) then
                 month = 5
                 day   = time - 120
              endif
              if (time .gt. 151 .and. time .le. 181) then
                 month = 6
                 day   = time - 151
              endif
              if (time .gt. 181 .and. time .le. 212) then
                 month = 7
                 day   = time - 181
              endif
              if (time .gt. 212 .and. time .le. 243) then
                 month = 8
                 day   = time - 212
              endif
              if (time .gt. 243 .and. time .le. 273) then
                 month = 9
                 day   = time - 243
              endif
              if (time .gt. 273 .and. time .le. 304) then
                 month = 10
                 day   = time - 273
              endif
              if (time .gt. 304 .and. time .le. 334) then
                 month = 11
                 day   = time - 304
              endif
              if (time .gt. 334 ) then
                 month = 12
                 day   = time - 334
              endif
              
              ! maximum temperature in °C
              ! minimium temperature in °C
              ! precipitation in mm / d
              ! evap not needed, therefore as dummy
              ! radiation in  MJ / d (/m2)
              ! vp in hPa
              ! pressure divided by 100 to get hPa  
              ! ca as constant value in ymol

              write(10, '(4i8, 8f8.2)')&
                   &       time,  day,  month, year,  & 
                   &       tmax(ilon, ilat, time), tmin(ilon, ilat, time),   &
                   &       prec(ilon, ilat, time), evap,   &
                   &       dswrf(ilon, ilat, time), vp(ilon, ilat, time), pres(ilon, ilat, time) / 100, ca
              
           enddo
           
           write(*,*) "WRITING OUTPUT DONE FOR YEAR", year 
           
!=================================

        enddo   ! year loop
        
        write(*,*) "----------------"
        write(*,*) "PROGRAM FINISHED"
        write(*,*) "----------------"
        
      END program read_nc_lsh_daily
