!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Write output data from the VOM to netcdf
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Currently only includes the reading of command line arguments
! Needs to contain all subroutines that read data in the future
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
!*-----File opening (saving climate and gstom ass data)-----------------

      subroutine vom_open_output (nc_flag)
      use vom_vegwat_mod
      use netcdf
      implicit none

     LOGICAL, INTENT(in) :: nc_flag
     CHARACTER(len = 100)             :: filename 
     CHARACTER(len = 250)             :: startdate 
     CHARACTER(len = 2)             :: day_tmp 
     CHARACTER(len = 2)             :: month_tmp 
     CHARACTER(len = 4)             :: year_tmp 
     integer :: lon_dimid
     integer :: lon_varid
     integer :: lat_dimid
     integer :: lat_varid
     integer :: dimids(3)

     integer :: time_dimid
     integer :: n_lat = 1
     integer :: n_lon = 1
     integer :: status 

      if(nc_flag .eqv. .True.) then


         ! 
         filename = trim(adjustl(i_outputpath))// &
                  trim(adjustl(nfile_resultsdaily))

         ! Create the file. 
         call check( nf90_create(filename, nf90_clobber, ncid) ) 

         ! Define the dimensions. NetCDF will hand back an ID for each. 
         call check(  nf90_def_dim(ncid, "lat", n_lat, lat_dimid) )
         call check(  nf90_def_dim(ncid, "lon", n_lon, lon_dimid) )
         call check(  nf90_def_dim(ncid, "time", NF90_UNLIMITED, time_dimid)) 

         call check(  nf90_def_var(ncid, "lat", NF90_REAL, lat_dimid, lat_varid) )
         call check(  nf90_def_var(ncid, "lon", NF90_REAL, lon_dimid, lon_varid) )
         call check(  nf90_def_var(ncid, "time", NF90_INT, time_dimid, time_varid) )

         ! Assign units attributes to coordinate variables.
         call check(  nf90_put_att(ncid, lat_varid, "units", "degrees_north") )
         call check(  nf90_put_att(ncid, lon_varid, "units", "degrees_south") )

         ! internal write for datestamp
         write(day_tmp,'(I2)') fday(1)
         write(month_tmp,'(I2)') fmonth(1)
         write(year_tmp,'(I4)') fyear(1)

         startdate = adjustl( "days since ")//trim(adjustl(day_tmp))// & 
                 trim(adjustl("-"))//trim(adjustl( month_tmp))// &
                  trim(adjustl("-"))//trim(adjustl(year_tmp))  

         call check(  nf90_put_att(ncid, time_varid, "units", startdate) )
    
         ! Array with id of dimensions
         dimids = (/ lon_dimid, lat_dimid, time_dimid /)

         ! Defining the variables.
         call check(  nf90_def_var(ncid, 'rain', NF90_REAL, dimids, rain_varid) )
         call check(  nf90_def_var(ncid, 'tairmax', NF90_REAL, dimids, tairmax_varid) )
         call check(  nf90_def_var(ncid, 'tairmin', NF90_REAL, dimids, tairmin_varid) )
         call check(  nf90_def_var(ncid, 'par', NF90_REAL, dimids, par_varid) )
         call check(  nf90_def_var(ncid, 'vd', NF90_REAL, dimids, vd_varid) )
         call check(  nf90_def_var(ncid, 'esoil', NF90_REAL, dimids, esoil_varid) )
         call check(  nf90_def_var(ncid, 'jmax25t', NF90_REAL, dimids, jmax25t_varid) )
         call check(  nf90_def_var(ncid, 'jmax25g', NF90_REAL, dimids, jmax25g_varid) )
         call check(  nf90_def_var(ncid, 'pc', NF90_REAL, dimids, pc_varid) )
         call check(  nf90_def_var(ncid, 'rlt', NF90_REAL, dimids, rlt_varid) )
         call check(  nf90_def_var(ncid, 'rlg', NF90_REAL, dimids, rlg_varid) )
         call check(  nf90_def_var(ncid, 'lambdat', NF90_REAL, dimids, lambdat_varid) )
         call check(  nf90_def_var(ncid, 'lambdag', NF90_REAL, dimids, lambdag_varid) )
         call check(  nf90_def_var(ncid, 'rrt', NF90_REAL, dimids, rrt_varid) )
         call check(  nf90_def_var(ncid, 'rrg', NF90_REAL, dimids, rrg_varid) )
         call check(  nf90_def_var(ncid, 'asst', NF90_REAL, dimids, asst_varid) )
         call check(  nf90_def_var(ncid, 'assg', NF90_REAL, dimids, assg_varid) )
         call check(  nf90_def_var(ncid, 'su_avg', NF90_REAL, dimids, su_avg_varid) )
         call check(  nf90_def_var(ncid, 'zw', NF90_REAL, dimids, zw_varid) )
         call check(  nf90_def_var(ncid, 'ws', NF90_REAL, dimids, ws_varid) )
         call check(  nf90_def_var(ncid, 'spgfcf', NF90_REAL, dimids, spgfcf_varid) )
         call check(  nf90_def_var(ncid, 'infx', NF90_REAL, dimids, infx_varid) )
         call check(  nf90_def_var(ncid, 'etmt', NF90_REAL, dimids, etmt_varid))
         call check(  nf90_def_var(ncid, 'etmg', NF90_REAL, dimids, etmg_varid) )
         call check(  nf90_def_var(ncid, 'su_1', NF90_REAL, dimids, su_1_varid) )
         call check(  nf90_def_var(ncid, 'topt', NF90_REAL, dimids, topt_varid) )
         call check(  nf90_def_var(ncid, 'tcg', NF90_REAL, dimids, tcg_varid) )
         call check(  nf90_def_var(ncid, 'tct', NF90_REAL, dimids, tct_varid) )
         call check(  nf90_def_var(ncid, 'cpccg', NF90_REAL, dimids, cpccg_d_varid) )
         call check(  nf90_def_var(ncid, 'cpcct', NF90_REAL, dimids, cpcct_d_varid) )
         call check(  nf90_def_var(ncid, 'lai_t', NF90_REAL, dimids, lai_t_varid) )
         call check(  nf90_def_var(ncid, 'lai_g', NF90_REAL, dimids, lai_g_varid) )
         call check(  nf90_def_var(ncid, 'ncp_g', NF90_REAL, dimids, ncp_g_varid) )
         call check(  nf90_def_var(ncid, 'ncp_t', NF90_REAL, dimids, ncp_t_varid) )

         ! Add units to the variables.
         call check(  nf90_put_att(ncid, rain_varid, "units", "mm/d") )
         call check(  nf90_put_att(ncid, tairmax_varid, "units", "oC") )
         call check(  nf90_put_att(ncid, tairmin_varid, "units", "oC") )
         call check(  nf90_put_att(ncid, par_varid, "units", "mol m-2 d-1") )
         call check(  nf90_put_att(ncid, vd_varid, "units", "mol mol-1") )
         call check(  nf90_put_att(ncid, esoil_varid, "units", "m/d") )
         call check(  nf90_put_att(ncid, jmax25t_varid, "units", "mol m-2 s-1") )
         call check(  nf90_put_att(ncid, jmax25g_varid, "units", "mol m-2 s-1") )
         call check(  nf90_put_att(ncid, pc_varid, "units", "-") )
         call check(  nf90_put_att(ncid, rlt_varid, "units", "mol m-2 d-1") )
         call check(  nf90_put_att(ncid, rlg_varid, "units", "mol m-2 d-1") )
         call check(  nf90_put_att(ncid, lambdat_varid, "units", "mol mol-1") )
         call check(  nf90_put_att(ncid, lambdag_varid, "units", "mol mol-1") )
         call check(  nf90_put_att(ncid, rrt_varid, "units", "mol m-2 d-1") )
         call check(  nf90_put_att(ncid, rrg_varid, "units", "mol m-2 d-1") )
         call check(  nf90_put_att(ncid, asst_varid, "units", "mol m-2 d-1") )
         call check(  nf90_put_att(ncid, assg_varid, "units", "mol m-2 d-1") )
         call check(  nf90_put_att(ncid, su_avg_varid, "units", "-") )
         call check(  nf90_put_att(ncid, zw_varid, "units", "m") )
         call check(  nf90_put_att(ncid, ws_varid, "units", "m"))
         call check(  nf90_put_att(ncid, spgfcf_varid, "units", "m/d") )
         call check(  nf90_put_att(ncid, infx_varid, "units", "m/d") )
         call check(  nf90_put_att(ncid, etmt_varid, "units", "m/d") )
         call check(  nf90_put_att(ncid, etmg_varid, "units", "m/d") )
         call check(  nf90_put_att(ncid, su_1_varid, "units", "-") )
         call check(  nf90_put_att(ncid, topt_varid, "units", "oC") )
         call check(  nf90_put_att(ncid, tcg_varid, "units", "mol/m2/s") )
         call check(  nf90_put_att(ncid, tct_varid, "units", "mol/m2/s") )
         call check(  nf90_put_att(ncid, cpccg_d_varid, "units", "mol/m2/s")) 
         call check(  nf90_put_att(ncid, cpcct_d_varid, "units", "mol/m2/s")) 
         call check(  nf90_put_att(ncid, lai_t_varid, "units", "-") )
         call check(  nf90_put_att(ncid, lai_g_varid, "units", "-") )
         call check(  nf90_put_att(ncid, ncp_g_varid, "units", "mol m-2") )  
         call check(  nf90_put_att(ncid, ncp_t_varid, "units", "mol m-2") )  

         ! End define mode.
         call check(  nf90_enddef(ncid) )

         ! Write the coordinate variable data. This will put the latitudes
         ! and longitudes of our data grid into the netCDF file.
         call check(  nf90_put_var(ncid, lat_varid, i_lat) )
         !call check(  nf90_put_var(ncid, lon_varid, i_lon) 

         !call check( nf90_put_var(ncid, time_varid, nday  ) )

      else
      !else plain text files instead of netcdf

         open(kfile_resultsdaily, FILE=trim(adjustl(i_outputpath))// &
           trim(adjustl(sfile_resultsdaily)), STATUS='replace')
         write(kfile_resultsdaily,'(A6,A7,A7,A7,A7, 34A15)') 'fyear',      &
         &  'fmonth', 'fday', 'nday', 'nhour', 'rain', 'tairmax', 'tairmin', &
         &  'par', 'vd', 'esoil', 'jmax25t', 'jmax25g', 'pc', 'rlt', 'rlg', &
         &  'lambdat', 'lambdag', 'rrt', 'rrg', 'asst', 'assg', 'su_avg',  &
         &  'zw', 'ws', 'spgfcf', 'infx', 'etmt', 'etmg', 'su_1', 'topt',  &
         &  'tcg', 'tct', 'cpccg_d', 'cpcct_d',  'lai_t', 'lai_g', &
         &'ncp_g', 'ncp_t'


      end if


      open(kfile_resultshourly, FILE=trim(adjustl(i_outputpath))//      &
           trim(adjustl(sfile_resultshourly)), STATUS='replace')
      write(kfile_resultshourly,'(A6,A7,A7,A7,A7,22A15)') 'fyear',     &
     &  'fmonth', 'fday', 'nday', 'nhour', 'rain', 'tair', 'par', 'vd',&
     &  'esoil', 'pc', 'jmax25t', 'jmax25g', 'mqt', 'rl', 'lambdat',   &
     &  'lambdag', 'rr', 'asst', 'assg', 'etmt', 'etmg', 'su_1',       &
     &  'zw', 'ws', 'spgfcf', 'infx'



      open(kfile_resultsyearly, FILE=trim(adjustl(i_outputpath))// &
           trim(adjustl(sfile_resultsyearly)), STATUS='replace')
      write(kfile_resultsyearly,'(A6,18A16)') "nyear", "rain", "par",  &
     &  "srad", "vd", "esoil", "etmt", "etmg", "assg", "rlg", "rrg",   &
     &  "cpccg", "tcg", "etmt", "asst", "rlt", "rrt", "cpcct", "tct"

        open(kfile_rsurfdaily, FILE=trim(adjustl(i_outputpath))// &
             trim(adjustl(sfile_rsurfdaily)), STATUS='replace')
        write(kfile_rsurfdaily,'(2A6,A4,A7,A)') 'fyear', 'fmonth',     &
     &    'fday', 'nday', 'rsurft_sublayer'

        open(kfile_delzhourly, FILE=trim(adjustl(i_outputpath))// &
             trim(adjustl(sfile_delzhourly)), STATUS='replace')
        write(kfile_delzhourly,'(2A6,A4,A7,A5,A)') 'fyear', 'fmonth',  &
     &    'fday', 'nday', 'nhour', 'delz_sublayer'

        open(kfile_ruptkthourly, FILE=trim(adjustl(i_outputpath))// &
             trim(adjustl(sfile_ruptkthourly)), STATUS='replace')
        write(kfile_ruptkthourly,'(2A6,A4,A7,A5,A)') 'fyear', 'fmonth',&
     &    'fday', 'nday', 'nhour', 'ruptkt_sublayer'

        open(kfile_suhourly, FILE=trim(adjustl(i_outputpath))// &
             trim(adjustl(sfile_suhourly)), STATUS='replace')
        write(kfile_suhourly,'(2A6,A4,A7,A5,A)') 'fyear', 'fmonth',    &
     &    'fday', 'nday', 'nhour', 'su_sublayer'



      return
      end subroutine vom_open_output


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     subroutine vom_write_day ( rain, tairmax, tairmin, par,   &
             &  vd, esoil, jmax25t, jmax25g,             &
             &  pc, rlt , rlg, lambdat, lambdag,         &
             &  rrt, rrg , asst, &
             &  assg, su_avg, zw, ws,     &
             &  spgfcf, infx, etmt, etmg, su_1, topt,              &
             & tcg, tct, cpccg, cpcct,                  &
             & lai_t, lai_g, tp_netassg, tp_netasst, rsurft, nc_flag )


      use vom_vegwat_mod
      use netcdf
      implicit none

      REAL*8,  INTENT(in) :: rain
      REAL*8,  INTENT(in) :: tairmax
      REAL*8,  INTENT(in) :: tairmin
      REAL*8,  INTENT(in) :: par
      REAL*8,  INTENT(in) :: vd
      REAL*8,  INTENT(in) :: esoil
      REAL*8,  INTENT(in) :: jmax25t
      REAL*8,  INTENT(in) :: jmax25g
      REAL*8,  INTENT(in) :: pc
      REAL*8,  INTENT(in) :: rlt
      REAL*8,  INTENT(in) :: rlg
      REAL*8,  INTENT(in) :: lambdat
      REAL*8,  INTENT(in) :: lambdag
      REAL*8,  INTENT(in) :: rrt
      REAL*8,  INTENT(in) :: rrg
      REAL*8,  INTENT(in) :: asst
      REAL*8,  INTENT(in) :: assg
      REAL*8,  INTENT(in) :: su_avg
      REAL*8,  INTENT(in) :: zw
      REAL*8,  INTENT(in) :: ws
      REAL*8,  INTENT(in) :: spgfcf
      REAL*8,  INTENT(in) :: infx
      REAL*8,  INTENT(in) :: etmt
      REAL*8,  INTENT(in) :: etmg
      REAL*8,  INTENT(in) :: su_1
      REAL*8,  INTENT(in) :: topt
      REAL*8,  INTENT(in) :: tcg
      REAL*8,  INTENT(in) :: tct
      REAL*8,  INTENT(in) :: cpccg
      REAL*8,  INTENT(in) :: cpcct
      REAL*8,  INTENT(in) :: lai_t
      REAL*8,  INTENT(in) :: lai_g
      REAL*8,  INTENT(in) :: tp_netassg
      REAL*8,  INTENT(in) :: tp_netasst
      REAL*8,  DIMENSION(s_maxlayer), INTENT(in):: rsurft
      LOGICAL, INTENT(in) :: nc_flag

      integer :: start(3)
      integer :: count(3)
      CHARACTER(60) :: dailyformat
      CHARACTER(3)  :: str


     if(nc_flag .eqv. .TRUE.) then
      count = (/ 1, 1, 1 /)
      start = (/ 1, 1, nday /)

      !add timestep
      call check( nf90_put_var(ncid, time_varid, nday, start= (/ nday/)  ) )

      !add variable values
      call check( nf90_put_var(ncid, rain_varid, rain, start = start ) )
      call check( nf90_put_var(ncid, tairmax_varid, tairmax, start = start ) )
      call check( nf90_put_var(ncid, tairmin_varid, tairmin, start = start ) )
      call check( nf90_put_var(ncid, par_varid, par, start = start ) )
      call check( nf90_put_var(ncid, vd_varid, vd, start = start ) )
      call check( nf90_put_var(ncid, esoil_varid, esoil, start = start ) )
      call check( nf90_put_var(ncid, jmax25t_varid, jmax25t, start = start ) )
      call check( nf90_put_var(ncid, jmax25g_varid, jmax25g, start = start ) )
      call check( nf90_put_var(ncid, pc_varid, pc, start = start ) )
      call check( nf90_put_var(ncid, rlt_varid, rlt, start = start ) )
      call check( nf90_put_var(ncid, rlg_varid, rlg, start = start ) )
      call check( nf90_put_var(ncid, lambdat_varid, lambdat, start = start ) )
      call check( nf90_put_var(ncid, lambdag_varid, lambdag, start = start ) )
      call check( nf90_put_var(ncid, rrt_varid, rrt, start = start ) )
      call check( nf90_put_var(ncid, rrg_varid, rrg, start = start ) )
      call check( nf90_put_var(ncid, asst_varid, asst, start = start ) )
      call check( nf90_put_var(ncid, assg_varid, assg, start = start ) )
      call check( nf90_put_var(ncid, su_avg_varid, su_avg, start = start ) )
      call check( nf90_put_var(ncid, zw_varid, zw, start = start ) )
      call check( nf90_put_var(ncid, ws_varid, ws, start = start ) )
      call check( nf90_put_var(ncid, spgfcf_varid, spgfcf, start = start ) )
      call check( nf90_put_var(ncid, infx_varid, infx, start = start ) )
      call check( nf90_put_var(ncid, etmt_varid, etmt, start = start ) )
      call check( nf90_put_var(ncid, etmg_varid, etmg, start = start ) )
      call check( nf90_put_var(ncid, su_1_varid, su_1, start = start ) )
      call check( nf90_put_var(ncid, topt_varid, topt, start = start ) )
      call check( nf90_put_var(ncid, tcg_varid, tcg, start = start ) )
      call check( nf90_put_var(ncid, tct_varid, tct, start = start ) )
      call check( nf90_put_var(ncid, cpccg_d_varid, cpccg, start = start ) )
      call check( nf90_put_var(ncid, cpcct_d_varid, cpcct, start = start ) )
      call check( nf90_put_var(ncid, lai_t_varid, lai_t, start = start ) )
      call check( nf90_put_var(ncid, lai_g_varid, lai_g, start = start ) )
      call check( nf90_put_var(ncid, ncp_g_varid, tp_netassg, start = start ) )
      call check( nf90_put_var(ncid, ncp_t_varid, tp_netasst, start = start ) )

     else

!     * internal write to convert from number to string
      write(str,'(I3)') wlayer_
!     * includes a column for each sublayer
      dailyformat = '(I6,I6,I4,I7,'//str//'E14.6)'

      write(kfile_resultsdaily,'(I6,I7,I7,I7,I7,34E15.5)')             &
     &  fyear(nday), fmonth(nday), fday(nday), nday, nhour-1,          &
     &  rain, tairmax, tairmin, par,   &
     &  vd, esoil, jmax25t, jmax25g,             &
     &  pc, rlt , rlg, lambdat, lambdag,         &
     &  rrt, rrg, asst, &
     &  assg, su_avg, zw, ws,     &
     &  spgfcf, infx, etmt, etmg, su_1, topt,              &
     & tcg, tct, cpccg, cpcct,                  &
     & lai_t, lai_g, tp_netassg, tp_netasst         


     write(kfile_rsurfdaily,dailyformat) fyear(nday), fmonth(nday), &
     &    fday(nday), nday, rsurft(1:wlayer_) 

      end if

      return
      end subroutine vom_write_day

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     subroutine vom_write_year ( n_year, rain_yearly, par_yearly, srad_yearly, vd_yearly, esoil_yearly, &
             & etm_yearly, etmg_yearly, assg_yearly, rlg_yearly, rrg_yearly, cpccg_yearly, tcg_yearly,&
             & etmt_yearly, asst_yearly, rlt_yearly, rrt_yearly, cpcct_yearly, tct_yearly )

      use vom_vegwat_mod
      implicit none

      INTEGER,  INTENT(in) :: n_year
      REAL*8,  INTENT(in) :: rain_yearly
      REAL*8,  INTENT(in) :: par_yearly
      REAL*8,  INTENT(in) :: srad_yearly
      REAL*8,  INTENT(in) :: vd_yearly
      REAL*8,  INTENT(in) :: esoil_yearly
      REAL*8,  INTENT(in) :: etm_yearly
      REAL*8,  INTENT(in) :: etmg_yearly
      REAL*8,  INTENT(in) :: assg_yearly
      REAL*8,  INTENT(in) :: rlg_yearly
      REAL*8,  INTENT(in) :: rrg_yearly
      REAL*8,  INTENT(in) :: cpccg_yearly
      REAL*8,  INTENT(in) :: tcg_yearly
      REAL*8,  INTENT(in) :: etmt_yearly
      REAL*8,  INTENT(in) :: asst_yearly
      REAL*8,  INTENT(in) :: rlt_yearly
      REAL*8,  INTENT(in) :: rrt_yearly
      REAL*8,  INTENT(in) :: cpcct_yearly
      REAL*8,  INTENT(in) :: tct_yearly



      write(kfile_resultsyearly,'(i6,18e16.6)') n_year, rain_yearly,       &
     &    par_yearly, srad_yearly, vd_yearly, esoil_yearly, etm_yearly,     &
     &    etmg_yearly, assg_yearly, rlg_yearly, rrg_yearly, cpccg_yearly, tcg_yearly,                &
     &    etmt_yearly, asst_yearly, rlt_yearly, rrt_yearly, cpcct_yearly, tct_yearly


      return
      end subroutine vom_write_year

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine check(status)
    use netcdf
    integer, intent ( in) :: status


    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  






