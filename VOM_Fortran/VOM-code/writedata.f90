!***********************************************************************
!        Optimised Vegetation Optimality Model (VOM)
!        Write output data from the VOM to netcdf
!        Original code coming from: https://github.com/schymans/VOM
!-----------------------------------------------------------------------
!        Author: Stan Schymanski
!
!        Now at: LIST, Luxembourg Institute of Science and Technology,
!                Belvaux, Luxembourg
!
!        Version: reading of command line arguments
!        Needs to contain all subroutines that read data in the future
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
     CHARACTER(len = 250)             :: refdate 
     CHARACTER(len = 2)             :: day_tmp 
     CHARACTER(len = 2)             :: month_tmp 
     CHARACTER(len = 4)             :: year_tmp 
     CHARACTER(len = 4)             :: firstyear_hourly 

     integer :: lon_dimid
     integer :: lon_varid
     integer :: lon_rsurf_dimid
     integer :: lon_rsurf_varid
     integer :: lon_hourly_dimid
     integer :: lon_hourly_varid
     integer :: lon_yearly_dimid
     integer :: lon_yearly_varid
     integer :: lon_ruptkt_dimid
     integer :: lon_ruptkt_varid
     integer :: lon_suhourly_dimid
     integer :: lon_suhourly_varid

     integer :: lat_dimid
     integer :: lat_varid
     integer :: lat_rsurf_dimid
     integer :: lat_rsurf_varid
     integer :: lat_hourly_dimid
     integer :: lat_hourly_varid
     integer :: lat_yearly_dimid
     integer :: lat_yearly_varid
     integer :: lat_ruptkt_dimid
     integer :: lat_ruptkt_varid
     integer :: lat_suhourly_dimid
     integer :: lat_suhourly_varid

     integer :: z_dimid
     integer :: z_varid
     integer :: z_ruptkt_varid
     integer :: z_suhourly_varid
     integer :: z_ruptkt_dimid
     integer :: z_suhourly_dimid
     integer :: dimids(3)
     integer :: dimids_surf(4)
     integer :: dimids_hourly(3)
     integer :: dimids_yearly(3)
     integer :: dimids_ruptkt(4)
     integer :: dimids_suhourly(4)



     integer :: time_dimid
     integer :: time_rsurf_dimid
     integer :: time_hourly_dimid
     integer :: time_yearly_dimid
     integer :: time_ruptkt_dimid
     integer :: time_suhourly_dimid

     integer :: n_lat = 1
     integer :: n_lon = 1
     integer :: status 
     integer :: i
     integer :: jday_base
     integer :: jday_start
     integer :: jday_startout
     real*8, dimension(s_maxlayer) :: depth 


      if(nc_flag .eqv. .True.) then

         ! **************************************************************
         ! first for dailyresults
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
         call check(  nf90_put_att(ncid, lon_varid, "units", "degrees_east") )


         ! internal write for datestamp
         write(day_tmp,'(I2)') fday(1)
         write(month_tmp,'(I2)') fmonth(1)
         write(year_tmp,'(I4)') fyear(1)

         refdate = adjustl( "days since ")//trim(adjustl("01"))// & 
                 trim(adjustl("-"))//trim(adjustl( "01"))// &
                  trim(adjustl("-"))//trim(adjustl("1900"))//trim(" 00:00:00")

         call get_julday(1900, 1, 1, jday_base)
         call get_julday(fyear(1), fmonth(1), fday(1), jday_start)
         startday = jday_start - jday_base

         call check(  nf90_put_att(ncid, time_varid, "units", refdate) )
    
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
         call check(  nf90_put_var(ncid, lon_varid, i_lon) )


         ! **************************************************************
         ! daily results rsurf
         filename = trim(adjustl(i_outputpath))// &
                  trim(adjustl(nfile_rsurfdaily))

         ! Create the file. 
         call check( nf90_create(filename, nf90_clobber, ncid_rsurf) ) 

         ! Define the dimensions. NetCDF will hand back an ID for each. 
         call check(  nf90_def_dim(ncid_rsurf, "latitude", n_lat, lat_rsurf_dimid) )
         call check(  nf90_def_dim(ncid_rsurf, "longitude", n_lon, lon_rsurf_dimid) )
         call check(  nf90_def_dim(ncid_rsurf, "time", NF90_UNLIMITED, time_rsurf_dimid)) 
         call check(  nf90_def_dim(ncid_rsurf, "level", s_maxlayer, z_dimid) )

         call check(  nf90_def_var(ncid_rsurf, "latitude", NF90_REAL, lat_rsurf_dimid, lat_rsurf_varid) )
         call check(  nf90_def_var(ncid_rsurf, "longitude", NF90_REAL, lon_rsurf_dimid, lon_rsurf_varid) )
         call check(  nf90_def_var(ncid_rsurf, "time", NF90_INT, time_rsurf_dimid, time_rsurf_varid) )
         call check(  nf90_def_var(ncid_rsurf, "level", NF90_REAL, z_dimid, z_varid) )

         ! Assign units attributes to coordinate variables.
         call check(  nf90_put_att(ncid_rsurf, time_rsurf_varid, "units", refdate) )
         call check(  nf90_put_att(ncid_rsurf, lat_rsurf_varid, "units", "degrees_north") )
         call check(  nf90_put_att(ncid_rsurf, lon_rsurf_varid, "units", "degrees_east") )
         call check(  nf90_put_att(ncid_rsurf, z_varid, "units", "meters") )

         dimids_surf = (/  lat_rsurf_dimid, lon_rsurf_dimid,  z_dimid, time_rsurf_dimid /)

         !defining the variable
         call check(  nf90_def_var(ncid_rsurf, 'rsurf', NF90_REAL, dimids_surf, rsurf_varid) )

         ! Add units to the variable.
         call check(  nf90_put_att(ncid_rsurf, rsurf_varid, "units", "m2 m-3 d-1") )  

         ! End define mode.
         call check(  nf90_enddef(ncid_rsurf) )

         call check(  nf90_put_var(ncid_rsurf, lat_rsurf_varid, i_lat ) )
         call check(  nf90_put_var(ncid_rsurf, lon_rsurf_varid, i_lon ))

         depth(1) = 0.5*s_delz(1)
         do i=2, s_maxlayer             
           depth(i) = (sum(s_delz(1:(i-1) )) + 0.5*s_delz(i) ) 
         end do

         call check(  nf90_put_var(ncid_rsurf, z_varid, depth) )

         ! **************************************************************
         ! results hourly
         filename = trim(adjustl(i_outputpath))// &
                  trim(adjustl(nfile_resultshourly))

         ! Create the file. 
         call check( nf90_create(filename, nf90_clobber, ncid_hourly) ) 

         ! Define the dimensions. NetCDF will hand back an ID for each. 
         call check(  nf90_def_dim(ncid_hourly, "lat", n_lat, lat_hourly_dimid) )
         call check(  nf90_def_dim(ncid_hourly, "lon", n_lon, lon_hourly_dimid) )
         call check(  nf90_def_dim(ncid_hourly, "time", NF90_UNLIMITED, time_hourly_dimid)) 

         call check(  nf90_def_var(ncid_hourly, "lat", NF90_REAL, lat_hourly_dimid, lat_hourly_varid) )
         call check(  nf90_def_var(ncid_hourly, "lon", NF90_REAL, lon_hourly_dimid, lon_hourly_varid) )
         call check(  nf90_def_var(ncid_hourly, "time", NF90_INT, time_hourly_dimid, time_hourly_varid) )

         ! Assign units attributes to coordinate variables.
         call check(  nf90_put_att(ncid_hourly, lat_hourly_varid, "units", "degrees_north") )
         call check(  nf90_put_att(ncid_hourly, lon_hourly_varid, "units", "degrees_east") )

         ! internal write for datestamp
         write(firstyear_hourly,'(I4)') i_firstyear
         !refdate = adjustl( "hours since ")//trim(adjustl(day_tmp))// & 
         !        trim(adjustl("-"))//trim(adjustl( month_tmp))// &
         !         trim(adjustl("-"))//trim(adjustl(year_tmp))  
         refdate = adjustl( "hours since ")//trim(adjustl("01"))// & 
                 trim(adjustl("-"))//trim(adjustl( "01"))// &
                  trim(adjustl("-"))//trim(adjustl("1900"))//trim(" 00:00:00")


         call get_julday(i_firstyear, 1, 1, jday_startout)

         starthour = ((jday_startout - jday_base) * 24)  -1

         !write(*,*) starthour    

         call check(  nf90_put_att(ncid_hourly, time_hourly_varid, "units", refdate) )

         ! Array with id of dimensions
         dimids_hourly = (/ lon_hourly_dimid, lat_hourly_dimid, time_hourly_dimid /)

        ! Defining the variables.
         call check(  nf90_def_var(ncid_hourly, 'rain', NF90_REAL, dimids_hourly, rainh_varid) )
         call check(  nf90_def_var(ncid_hourly, 'tair', NF90_REAL, dimids_hourly, tairh_varid) )
         call check(  nf90_def_var(ncid_hourly, 'par', NF90_REAL, dimids_hourly, parh_varid) )
         call check(  nf90_def_var(ncid_hourly, 'vd', NF90_REAL, dimids_hourly, vdh_varid) )
         call check(  nf90_def_var(ncid_hourly, 'esoil', NF90_REAL, dimids_hourly, esoilh_varid) )
         call check(  nf90_def_var(ncid_hourly, 'pc', NF90_REAL, dimids_hourly, pch_varid) )
         call check(  nf90_def_var(ncid_hourly, 'jmax25t', NF90_REAL, dimids_hourly, jmax25th_varid) )
         call check(  nf90_def_var(ncid_hourly, 'jmax25g', NF90_REAL, dimids_hourly, jmax25gh_varid) )
         call check(  nf90_def_var(ncid_hourly, 'mqt', NF90_REAL, dimids_hourly, mqth_varid) )
         call check(  nf90_def_var(ncid_hourly, 'rl', NF90_REAL, dimids_hourly, rlh_varid) )
         call check(  nf90_def_var(ncid_hourly, 'lambdat', NF90_REAL, dimids_hourly, lambdath_varid) )
         call check(  nf90_def_var(ncid_hourly, 'lambdag', NF90_REAL, dimids_hourly, lambdagh_varid) )
         call check(  nf90_def_var(ncid_hourly, 'rr', NF90_REAL, dimids_hourly, rrh_varid) )
         call check(  nf90_def_var(ncid_hourly, 'asst', NF90_REAL, dimids_hourly, assth_varid) )
         call check(  nf90_def_var(ncid_hourly, 'assg', NF90_REAL, dimids_hourly, assgh_varid) )
         call check(  nf90_def_var(ncid_hourly, 'etmt', NF90_REAL, dimids_hourly, etmth_varid))
         call check(  nf90_def_var(ncid_hourly, 'etmg', NF90_REAL, dimids_hourly, etmgh_varid) )
         call check(  nf90_def_var(ncid_hourly, 'su_1', NF90_REAL, dimids_hourly, su1h_varid) )
         call check(  nf90_def_var(ncid_hourly, 'zw', NF90_REAL, dimids_hourly, zwh_varid) )
         call check(  nf90_def_var(ncid_hourly, 'ws', NF90_REAL, dimids_hourly, wsh_varid) )
         call check(  nf90_def_var(ncid_hourly, 'spgfcf', NF90_REAL, dimids_hourly, spgfcfh_varid) )
         call check(  nf90_def_var(ncid_hourly, 'infx', NF90_REAL, dimids_hourly, infxh_varid) )


         ! Add units to the variables.
         call check(  nf90_put_att(ncid_hourly, rainh_varid, "units", "m/s") )
         call check(  nf90_put_att(ncid_hourly, tairh_varid, "units", "oC") )
         call check(  nf90_put_att(ncid_hourly, parh_varid, "units", "mol m-2 s-1") )
         call check(  nf90_put_att(ncid_hourly, vdh_varid, "units", "mol mol-1") )
         call check(  nf90_put_att(ncid_hourly, esoilh_varid, "units", "m/h") )
         call check(  nf90_put_att(ncid_hourly, jmax25th_varid, "units", "mol m-2 s-1") )
         call check(  nf90_put_att(ncid_hourly, jmax25gh_varid, "units", "mol m-2 s-1") )
         call check(  nf90_put_att(ncid_hourly, pch_varid, "units", "-") )
         call check(  nf90_put_att(ncid_hourly, rlh_varid, "units", "mol m-2 d-1") )
         call check(  nf90_put_att(ncid_hourly, mqth_varid, "units", "kg m-2") )
         call check(  nf90_put_att(ncid_hourly, lambdath_varid, "units", "mol mol-1") )
         call check(  nf90_put_att(ncid_hourly, lambdagh_varid, "units", "mol mol-1") )
         call check(  nf90_put_att(ncid_hourly, rrh_varid, "units", "mol m-2 s-1") )
         call check(  nf90_put_att(ncid_hourly, assth_varid, "units", "mol m-2 h-1") )
         call check(  nf90_put_att(ncid_hourly, assgh_varid, "units", "mol m-2 h-1") )
         call check(  nf90_put_att(ncid_hourly, su1h_varid, "units", "-") )
         call check(  nf90_put_att(ncid_hourly, zwh_varid, "units", "m") )
         call check(  nf90_put_att(ncid_hourly, wsh_varid, "units", "m"))
         call check(  nf90_put_att(ncid_hourly, spgfcfh_varid, "units", "m/h") )
         call check(  nf90_put_att(ncid_hourly, infxh_varid, "units", "m/h") )
         call check(  nf90_put_att(ncid_hourly, etmth_varid, "units", "m/h") )
         call check(  nf90_put_att(ncid_hourly, etmgh_varid, "units", "m/h") )

         ! End define mode.
         call check(  nf90_enddef(ncid_hourly) )

         call check(  nf90_put_var(ncid_hourly, lat_hourly_varid, i_lat ) )
         call check(  nf90_put_var(ncid_hourly, lon_hourly_varid, i_lon ))

         ! **************************************************************
         ! results yearly
         filename = trim(adjustl(i_outputpath))// &
                  trim(adjustl(nfile_resultsyearly))

         ! Create the file. 
         call check( nf90_create(filename, nf90_clobber, ncid_yearly) ) 

         ! Define the dimensions. NetCDF will hand back an ID for each. 
         call check(  nf90_def_dim(ncid_yearly, "lat", n_lat, lat_yearly_dimid) )
         call check(  nf90_def_dim(ncid_yearly, "lon", n_lon, lon_yearly_dimid) )
         call check(  nf90_def_dim(ncid_yearly, "time", NF90_UNLIMITED, time_yearly_dimid)) 

         call check(  nf90_def_var(ncid_yearly, "lat", NF90_REAL, lat_yearly_dimid, lat_yearly_varid) )
         call check(  nf90_def_var(ncid_yearly, "lon", NF90_REAL, lon_yearly_dimid, lon_yearly_varid) )
         call check(  nf90_def_var(ncid_yearly, "time", NF90_INT, time_yearly_dimid, time_yearly_varid) )

         ! Assign units attributes to coordinate variables.
         call check(  nf90_put_att(ncid_yearly, lat_yearly_varid, "units", "degrees_north") )
         call check(  nf90_put_att(ncid_yearly, lon_yearly_varid, "units", "degrees_east") )

         refdate = adjustl( "years since ")//trim(adjustl("01"))// & 
                 trim(adjustl("-"))//trim(adjustl( "01"))// &
                  trim(adjustl("-"))//trim(adjustl("1900"))//trim(" 00:00:00")

         call check(  nf90_put_att(ncid_yearly, time_yearly_varid, "units", refdate) )

         ! Array with id of dimensions
         dimids_yearly = (/ lon_yearly_dimid, lat_yearly_dimid, time_yearly_dimid /)


       ! Defining the variables.
         call check(  nf90_def_var(ncid_yearly, 'rain', NF90_REAL, dimids_yearly, rainy_varid) )
         call check(  nf90_def_var(ncid_yearly, 'par', NF90_REAL, dimids_yearly, pary_varid) )
         call check(  nf90_def_var(ncid_yearly, 'vd', NF90_REAL, dimids_yearly, vdy_varid) )
         call check(  nf90_def_var(ncid_yearly, 'srad', NF90_REAL, dimids_yearly, srady_varid) )
         call check(  nf90_def_var(ncid_yearly, 'esoil', NF90_REAL, dimids_yearly, esoily_varid) )
         call check(  nf90_def_var(ncid_yearly, 'etmg', NF90_REAL, dimids_yearly, etmgy_varid) )
         call check(  nf90_def_var(ncid_yearly, 'assg', NF90_REAL, dimids_yearly, assgy_varid) )
         call check(  nf90_def_var(ncid_yearly, 'rlg', NF90_REAL, dimids_yearly, rlgy_varid) )
         call check(  nf90_def_var(ncid_yearly, 'rrg', NF90_REAL, dimids_yearly, rrgy_varid) )
         call check(  nf90_def_var(ncid_yearly, 'cpccg', NF90_REAL, dimids_yearly, cpccgy_varid) )
         call check(  nf90_def_var(ncid_yearly, 'tcg', NF90_REAL, dimids_yearly, tcgy_varid) )
         call check(  nf90_def_var(ncid_yearly, 'etmt', NF90_REAL, dimids_yearly, etmty_varid))
         call check(  nf90_def_var(ncid_yearly, 'asst', NF90_REAL, dimids_yearly, assty_varid) )
         call check(  nf90_def_var(ncid_yearly, 'rlt', NF90_REAL, dimids_yearly, rlty_varid) )
         call check(  nf90_def_var(ncid_yearly, 'rrt', NF90_REAL, dimids_yearly, rrty_varid) )
         call check(  nf90_def_var(ncid_yearly, 'cpcct', NF90_REAL, dimids_yearly, cpccty_varid) )
         call check(  nf90_def_var(ncid_yearly, 'tct', NF90_REAL, dimids_yearly, tcty_varid) )

         ! Add units to the variables.
         call check(  nf90_put_att(ncid_yearly, rainy_varid, "units", "m/year") )
         call check(  nf90_put_att(ncid_yearly, pary_varid, "units", "mol m-2 s-1") )
         call check(  nf90_put_att(ncid_yearly, vdy_varid, "units", "mol mol-1") )
         call check(  nf90_put_att(ncid_yearly, srady_varid, "units", "MJ m-2 year-1") )
         call check(  nf90_put_att(ncid_yearly, esoily_varid, "units", "m/year") )
         call check(  nf90_put_att(ncid_yearly, etmgy_varid, "units", "m/year") )
         call check(  nf90_put_att(ncid_yearly, assgy_varid, "units", "mol m-2 year-1") )
         call check(  nf90_put_att(ncid_yearly, rlgy_varid, "units", "mol m-2 year-1") )
         call check(  nf90_put_att(ncid_yearly, rrgy_varid, "units", "mol m-2 year-1") )
         call check(  nf90_put_att(ncid_yearly, cpccgy_varid, "units", "mol m-2 year-1") )
         call check(  nf90_put_att(ncid_yearly, tcgy_varid, "units", "mol m-2 year-1") )
         call check(  nf90_put_att(ncid_yearly, etmty_varid, "units", "m/year") )
         call check(  nf90_put_att(ncid_yearly, assty_varid, "units", "mol m-2 year-1") )
         call check(  nf90_put_att(ncid_yearly, rlty_varid, "units", "mol m-2 year-1") )
         call check(  nf90_put_att(ncid_yearly, rrty_varid, "units", "mol m-2 year-1") )
         call check(  nf90_put_att(ncid_yearly, cpccty_varid, "units", "mol m-2 year-1") )
         call check(  nf90_put_att(ncid_yearly, tcty_varid, "units", "mol m-2 year-1") )

         ! End define mode.
         call check(  nf90_enddef(ncid_yearly) )

         call check(  nf90_put_var(ncid_yearly, lat_yearly_varid, i_lat ) )
         call check(  nf90_put_var(ncid_yearly, lon_yearly_varid, i_lon ))

         ! **************************************************************
         ! results root water uptake
         filename = trim(adjustl(i_outputpath))// &
                  trim(adjustl(nfile_ruptkthourly))

         ! Create the file. 
         call check( nf90_create(filename, nf90_clobber, ncid_ruptkt) ) 

         ! Define the dimensions. NetCDF will hand back an ID for each. 
         call check(  nf90_def_dim(ncid_ruptkt, "lat", n_lat, lat_ruptkt_dimid) )
         call check(  nf90_def_dim(ncid_ruptkt, "lon", n_lon, lon_ruptkt_dimid) )
         call check(  nf90_def_dim(ncid_ruptkt, "time", NF90_UNLIMITED, time_ruptkt_dimid)) 
         call check(  nf90_def_dim(ncid_ruptkt, "level", s_maxlayer, z_ruptkt_dimid) )

         call check(  nf90_def_var(ncid_ruptkt, "lat", NF90_REAL, lat_ruptkt_dimid, lat_ruptkt_varid) )
         call check(  nf90_def_var(ncid_ruptkt, "lon", NF90_REAL, lon_ruptkt_dimid, lon_ruptkt_varid) )
         call check(  nf90_def_var(ncid_ruptkt, "time", NF90_INT, time_ruptkt_dimid, time_ruptkt_varid) )
         call check(  nf90_def_var(ncid_ruptkt, "level", NF90_REAL, z_ruptkt_dimid, z_ruptkt_varid) )

         ! Assign units attributes to coordinate variables.
         call check(  nf90_put_att(ncid_ruptkt, lat_ruptkt_varid, "units", "degrees_north") )
         call check(  nf90_put_att(ncid_ruptkt, lon_ruptkt_varid, "units", "degrees_east") )
         call check(  nf90_put_att(ncid_ruptkt, z_ruptkt_varid, "units", "meters") )

         !startdate = adjustl( "hours since ")//trim(adjustl(day_tmp))// & 
         !       trim(adjustl("-"))//trim(adjustl( month_tmp))// &
         !         trim(adjustl("-"))//trim(adjustl(year_tmp)) 

         refdate = adjustl( "hours since ")//trim(adjustl("01"))// & 
                 trim(adjustl("-"))//trim(adjustl( "01"))// &
                  trim(adjustl("-"))//trim(adjustl("1900"))//trim(" 00:00:00")

         call check(  nf90_put_att(ncid_ruptkt, time_ruptkt_varid, "units", refdate) )

         ! Array with id of dimensions
         dimids_ruptkt = (/ lat_ruptkt_dimid, lon_ruptkt_dimid,  z_ruptkt_dimid, time_ruptkt_dimid /)

         !defining the variable
         call check(  nf90_def_var(ncid_ruptkt, 'ruptkt', NF90_REAL, dimids_ruptkt, ruptkt_varid) )

         ! Add units to the variable.
         call check(  nf90_put_att(ncid_ruptkt, ruptkt_varid, "units", "m/h") )  

         ! End define mode.
         call check(  nf90_enddef(ncid_ruptkt) )

         call check(  nf90_put_var(ncid_ruptkt, z_ruptkt_varid, depth) )
         call check(  nf90_put_var(ncid_ruptkt, lat_ruptkt_varid, i_lat ) )
         call check(  nf90_put_var(ncid_ruptkt, lon_ruptkt_varid, i_lon ))


         ! **************************************************************
         ! results hourly soil moisture
         filename = trim(adjustl(i_outputpath))// &
                  trim(adjustl(nfile_suhourly))

         ! Create the file. 
         call check( nf90_create(filename, nf90_clobber, ncid_suhourly) ) 

         ! Define the dimensions. NetCDF will hand back an ID for each. 
         call check(  nf90_def_dim(ncid_suhourly, "lat", n_lat, lat_suhourly_dimid) )
         call check(  nf90_def_dim(ncid_suhourly, "lon", n_lon, lon_suhourly_dimid) )
         call check(  nf90_def_dim(ncid_suhourly, "time", NF90_UNLIMITED, time_suhourly_dimid)) 
         call check(  nf90_def_dim(ncid_suhourly, "level", s_maxlayer, z_suhourly_dimid) )

         call check(  nf90_def_var(ncid_suhourly, "lat", NF90_REAL, lat_suhourly_dimid, lat_suhourly_varid) )
         call check(  nf90_def_var(ncid_suhourly, "lon", NF90_REAL, lon_suhourly_dimid, lon_suhourly_varid) )
         call check(  nf90_def_var(ncid_suhourly, "time", NF90_INT, time_suhourly_dimid, time_suhourly_varid) )
         call check(  nf90_def_var(ncid_suhourly, "level", NF90_REAL, z_suhourly_dimid, z_suhourly_varid) )

         ! Assign units attributes to coordinate variables.
         call check(  nf90_put_att(ncid_suhourly, lat_suhourly_varid, "units", "degrees_north") )
         call check(  nf90_put_att(ncid_suhourly, lon_suhourly_varid, "units", "degrees_east") )
         call check(  nf90_put_att(ncid_suhourly, z_suhourly_varid, "units", "meters") )

         !startdate = adjustl( "hours since ")//trim(adjustl(day_tmp))// & 
         !        trim(adjustl("-"))//trim(adjustl( month_tmp))// &
         !         trim(adjustl("-"))//trim(adjustl(year_tmp))   
         refdate = adjustl( "hours since ")//trim(adjustl("01"))// & 
                 trim(adjustl("-"))//trim(adjustl( "01"))// &
                  trim(adjustl("-"))//trim(adjustl("1900"))//trim(" 00:00:00")

         call check(  nf90_put_att(ncid_suhourly, time_suhourly_varid, "units", refdate) )

        ! Array with id of dimensions
         dimids_suhourly = (/ lat_suhourly_dimid, lon_suhourly_dimid,  z_suhourly_dimid, time_suhourly_dimid /)

         !defining the variable
         call check(  nf90_def_var(ncid_suhourly, 'su', NF90_REAL, dimids_suhourly, suhourly_varid) )

         ! Add units to the variable.
         call check(  nf90_put_att(ncid_suhourly, suhourly_varid, "units", "-") )  

         ! End define mode.
         call check(  nf90_enddef(ncid_suhourly) )

         call check(  nf90_put_var(ncid_suhourly, z_suhourly_varid, depth) )
         call check(  nf90_put_var(ncid_suhourly, lat_suhourly_varid, i_lat ) )
         call check(  nf90_put_var(ncid_suhourly, lon_suhourly_varid, i_lon ))

      else
      !else plain text files instead of netcdf

         open(kfile_resultsdaily, FILE=trim(adjustl(i_outputpath))// &
           trim(adjustl(sfile_resultsdaily)), STATUS='replace')
         write(kfile_resultsdaily,'(A6,A7,A7,A7,A7, 42A15)') 'fyear',         &
         &  'fmonth', 'fday', 'nday', 'nhour', 'rain', 'tairmax', 'tairmin',  &
         &  'par', 'vd', 'esoil', 'jmax25t', 'jmax25g', 'jmax25t', 'jmax25g', &
         &  'pc', 'rlt', 'rlg',                                               &
         &  'lambdat', 'lambdag', 'rrt', 'rrg', 'asst', 'assg', 'su_avg',     &
         &  'zw', 'ws', 'spgfcf', 'infx', 'etmt', 'etmg', 'su_1', 'topt',     &
         &  'tcg', 'tct', 'cpccg_d', 'cpcct_d', 'fsunt', 'fshadet',           &
         &  'fsung', 'fshadeg', 'lai_t', 'lai_g', 'lai_tot',                  &
         &   'cai_g', 'ncp_g', 'ncp_t'

        open(kfile_rsurfdaily, FILE=trim(adjustl(i_outputpath))// &
             trim(adjustl(sfile_rsurfdaily)), STATUS='replace')
        write(kfile_rsurfdaily,'(2A6,A4,A7,A)') 'fyear', 'fmonth',     &
     &    'fday', 'nday', 'rsurft_sublayer'


      open(kfile_resultshourly, FILE=trim(adjustl(i_outputpath))//      &
           trim(adjustl(sfile_resultshourly)), STATUS='replace')
      write(kfile_resultshourly,'(A6,A7,A7,A7,A7,24A15)') 'fyear',     &
     &  'fmonth', 'fday', 'nday', 'nhour', 'rain', 'tair', 'par',      &
     &  'gstomt', 'gstomg', 'vd',                                      &
     &  'esoil', 'pc', 'jmax25t', 'jmax25g', 'mqt', 'rl', 'lambdat',   &
     &  'lambdag', 'rr', 'asst', 'assg', 'etmt', 'etmg', 'su_1',       &
     &  'zw', 'ws', 'spgfcf', 'infx'



      open(kfile_resultsyearly, FILE=trim(adjustl(i_outputpath))// &
           trim(adjustl(sfile_resultsyearly)), STATUS='replace')
      write(kfile_resultsyearly,'(A6,18A16)') "nyear", "rain", "par",  &
     &  "srad", "vd", "esoil", "etmt", "etmg", "assg", "rlg", "rrg",   &
     &  "cpccg", "tcg", "etmt", "asst", "rlt", "rrt", "cpcct", "tct"



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

      end if

      return
      end subroutine vom_open_output


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     subroutine vom_write_day ( rain, tairmax, tairmin, par,         &
             &  vd, esoil, jmax25t, jmax25g, jmax25ts, jmax25gs,     &
             &  pc, rlt , rlg, lambdat, lambdag,                     &
             &  rrt, rrg , asst,                                     &
             &  assg, su_avg, zw, ws,                                &
             &  spgfcf, infx, etmt, etmg, su_1, topt,                &
             & tcg, tct, cpccg, cpcct,                               &
             & fsun_t, fshade_t, fsun_g, fshade_g,                   &
             & lai_t, lai_g, lai_tot, caig, tp_netassg, tp_netasst, rsurft, nc_flag )


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
      REAL*8,  INTENT(in) :: jmax25ts
      REAL*8,  INTENT(in) :: jmax25gs      
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
      REAL*8,  INTENT(in) :: fsun_t
      REAL*8,  INTENT(in) :: fshade_t  
      REAL*8,  INTENT(in) :: fsun_g
      REAL*8,  INTENT(in) :: fshade_g            
      REAL*8,  INTENT(in) :: lai_t
      REAL*8,  INTENT(in) :: lai_g
      REAL*8,  INTENT(in) :: lai_tot
      REAL*8,  INTENT(in) :: caig      
      REAL*8,  INTENT(in) :: tp_netassg
      REAL*8,  INTENT(in) :: tp_netasst
      REAL*8,  DIMENSION(s_maxlayer), INTENT(in):: rsurft
      LOGICAL, INTENT(in) :: nc_flag

      integer :: start(3)
      integer :: count(3)

      integer :: start_rsurf(4)
      integer :: count_rsurf(4)

      CHARACTER(60) :: dailyformat
      CHARACTER(3)  :: str


     if(nc_flag .eqv. .TRUE.) then
      count = (/ 1, 1, 1 /)
      start = (/ 1, 1,  nday /)

      !add timestep
      call check( nf90_put_var(ncid, time_varid, ( startday + nday -1 ) , start= (/ nday /)   ) )

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

      count_rsurf = (/ 1, 1, s_maxlayer, 1 /)
      start_rsurf = (/ 1, 1, 1, nday  /)
      !count_rsurf = (/ 1, 1, 1 /)
      !start_rsurf = (/ 1, 1, nday /)

      !add timestep
      call check( nf90_put_var(ncid_rsurf, time_rsurf_varid, ( startday + nday -1 ) , start= (/ nday/)  ) )
      call check( nf90_put_var(ncid_rsurf, rsurf_varid, rsurft, start = start_rsurf, count=count_rsurf ) )

      !call check( nf90_put_var(ncid_rsurf, rsurf_varid, rsurft(1), start = start_rsurf ) )

     else

!     * internal write to convert from number to string
      write(str,'(I3)') wlayer_
!     * includes a column for each sublayer
      dailyformat = '(I6,I6,I4,I7,'//str//'E14.6)'

      write(kfile_resultsdaily,'(I6,I7,I7,I7,I7,42E15.5)')    &
     &  fyear(nday), fmonth(nday), fday(nday), nday, nhour-1, &
     &  rain, tairmax, tairmin, par,                          &
     &  vd, esoil, jmax25t, jmax25g, jmax25ts, jmax25gs,      &
     &  pc, rlt , rlg, lambdat, lambdag,                      &
     &  rrt, rrg, asst,                                       &
     &  assg, su_avg, zw, ws,                                 &
     &  spgfcf, infx, etmt, etmg, su_1, topt,                 &
     &  tcg, tct, cpccg, cpcct,                               &
     &  fsun_t, fshade_t, fsun_g, fshade_g,                   &
     &  lai_t, lai_g, lai_tot, caig, tp_netassg, tp_netasst         


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
             & etmt_yearly, asst_yearly, rlt_yearly, rrt_yearly, cpcct_yearly, tct_yearly, nc_flag )

      use vom_vegwat_mod
      use netcdf
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
      LOGICAL, INTENT(in) :: nc_flag
      integer :: start(3)
      integer :: count(3)


     if(nc_flag .eqv. .TRUE.) then


      count = (/ 1, 1, 1 /)
      start = (/ 1, 1, n_year - fyear(1) + 1 /)


      !add timestep
      call check( nf90_put_var(ncid_yearly, time_yearly_varid, (n_year - 1900 ) , &
         start= (/ n_year - fyear(1) + 1/)  ) ) 



      !add variable values
      call check( nf90_put_var(ncid_yearly, rainy_varid, rain_yearly, start = start ) )
      call check( nf90_put_var(ncid_yearly, pary_varid, par_yearly, start = start ) )
      call check( nf90_put_var(ncid_yearly, vdy_varid, vd_yearly, start = start ) )
      call check( nf90_put_var(ncid_yearly, srady_varid, srad_yearly, start = start ) )
      call check( nf90_put_var(ncid_yearly, esoily_varid, esoil_yearly, start = start ) )
      call check( nf90_put_var(ncid_yearly, etmgy_varid, etmg_yearly, start = start ) )
      call check( nf90_put_var(ncid_yearly, assgy_varid, assg_yearly, start = start ) )
      call check( nf90_put_var(ncid_yearly, rlgy_varid, rlg_yearly, start = start ) )
      call check( nf90_put_var(ncid_yearly, rrgy_varid, rrg_yearly, start = start ) )
      call check( nf90_put_var(ncid_yearly, cpccgy_varid, cpccg_yearly, start = start ) )
      call check( nf90_put_var(ncid_yearly, tcgy_varid, tcg_yearly, start = start ) )
      call check( nf90_put_var(ncid_yearly, etmty_varid, etmt_yearly, start = start ) )
      call check( nf90_put_var(ncid_yearly, assty_varid, asst_yearly, start = start ) )
      call check( nf90_put_var(ncid_yearly, rlty_varid, rlt_yearly, start = start ) )
      call check( nf90_put_var(ncid_yearly, rrty_varid, rrt_yearly, start = start ) )
      call check( nf90_put_var(ncid_yearly, cpccty_varid, cpcct_yearly, start = start ) )
      call check( nf90_put_var(ncid_yearly, tcty_varid, tct_yearly, start = start ) )


     else

        write(kfile_resultsyearly,'(i6,18e16.6)') n_year, rain_yearly,       &
        &    par_yearly, srad_yearly, vd_yearly, esoil_yearly, etm_yearly,     &
        &    etmg_yearly, assg_yearly, rlg_yearly, rrg_yearly, cpccg_yearly, tcg_yearly,                &
        &    etmt_yearly, asst_yearly, rlt_yearly, rrt_yearly, cpcct_yearly, tct_yearly

     end if

      return
      end subroutine vom_write_year

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine vom_write_hourly ( year, month, day, num_day, num_hour, num_hour_tot,  &
          &    rain_hourly, tair_hourly, par_hourly, gstomt_hourly, gstomg_hourly ,vd_hourly, esoil_hourly,    &
          &    pc_hourly, jmax25t_hourly, jmax25g_hourly, mqt_hourly,          &
          &    rl_hourly, lambdat_hourly, lambdag_hourly, rr_hourly,  &
          &    asst_hourly, assg_hourly, etmt_hourly, etmg_hourly, su1_hourly, zw_hourly, ws_hourly, &
          &    spgfcf_hourly, infx_hourly, ruptkt_hourly, su_hourly, nc_flag )
      use vom_vegwat_mod
      use netcdf
      implicit none


      INTEGER,  INTENT(in) :: year
      INTEGER,  INTENT(in) :: month
      INTEGER,  INTENT(in) :: day
      INTEGER,  INTENT(in) :: num_day
      INTEGER,  INTENT(in) :: num_hour
      INTEGER,  INTENT(in) :: num_hour_tot
      REAL*8,  INTENT(in) :: rain_hourly
      REAL*8,  INTENT(in) :: tair_hourly
      REAL*8,  INTENT(in) :: par_hourly
      REAL*8,  INTENT(in) :: gstomt_hourly
      REAL*8,  INTENT(in) :: gstomg_hourly
      REAL*8,  INTENT(in) :: vd_hourly
      REAL*8,  INTENT(in) :: esoil_hourly
      REAL*8,  INTENT(in) :: pc_hourly
      REAL*8,  INTENT(in) :: jmax25t_hourly
      REAL*8,  INTENT(in) :: jmax25g_hourly
      REAL*8,  INTENT(in) :: mqt_hourly
      REAL*8,  INTENT(in) :: rl_hourly
      REAL*8,  INTENT(in) :: lambdat_hourly
      REAL*8,  INTENT(in) :: lambdag_hourly
      REAL*8,  INTENT(in) :: rr_hourly
      REAL*8,  INTENT(in) :: asst_hourly
      REAL*8,  INTENT(in) :: assg_hourly
      REAL*8,  INTENT(in) :: etmt_hourly
      REAL*8,  INTENT(in) :: etmg_hourly
      REAL*8,  INTENT(in) :: su1_hourly
      REAL*8,  INTENT(in) :: zw_hourly
      REAL*8,  INTENT(in) :: ws_hourly
      REAL*8,  INTENT(in) :: spgfcf_hourly
      REAL*8,  INTENT(in) :: infx_hourly
      REAL*8,  DIMENSION(s_maxlayer), INTENT(in):: ruptkt_hourly
      REAL*8,  DIMENSION(s_maxlayer), INTENT(in):: su_hourly
      LOGICAL, INTENT(in) :: nc_flag

      CHARACTER(60) :: hourlyformat
      CHARACTER(3)  :: str
      integer :: start(3)
      integer :: count(3)
      integer :: start_ruptkt(4)
      integer :: count_ruptkt(4)
      integer :: start_suhourly(4)
      integer :: count_suhourly(4)

     if(nc_flag .eqv. .TRUE.) then



      if (year .ge. i_firstyear .and. year .le. i_lastyear) then

      !*****************************
      !general hourly results

      if( (month .eq. 1) .and. (day .eq. 1) .and. (year .eq. i_firstyear) .and. (num_hour .eq. 1)) then
         hourdiff = num_hour_tot - 1
      end if

      count = (/ 1, 1, 1 /)
      start = (/ 1, 1, num_hour_tot-hourdiff /)

      !add timestep
      call check( nf90_put_var(ncid_hourly, time_hourly_varid, (starthour + num_hour_tot - hourdiff ), &
                                start= (/ num_hour_tot-hourdiff /)  ) )




      !add variable values
      call check( nf90_put_var(ncid_hourly, rainh_varid, rain_hourly, start = start ) )
      call check( nf90_put_var(ncid_hourly, tairh_varid, tair_hourly, start = start ) )
      call check( nf90_put_var(ncid_hourly, parh_varid, par_hourly, start = start ) )
      call check( nf90_put_var(ncid_hourly, vdh_varid, vd_hourly, start = start ) )
      call check( nf90_put_var(ncid_hourly, esoilh_varid, esoil_hourly, start = start ) )
      call check( nf90_put_var(ncid_hourly, pch_varid, pc_hourly, start = start ) )
      call check( nf90_put_var(ncid_hourly, jmax25th_varid, jmax25t_hourly, start = start ) )
      call check( nf90_put_var(ncid_hourly, jmax25gh_varid, jmax25g_hourly, start = start ) )
      call check( nf90_put_var(ncid_hourly, mqth_varid, mqt_hourly, start = start ) )
      call check( nf90_put_var(ncid_hourly, rlh_varid, rl_hourly, start = start ) )
      call check( nf90_put_var(ncid_hourly, lambdath_varid, lambdat_hourly, start = start ) )
      call check( nf90_put_var(ncid_hourly, lambdagh_varid, lambdag_hourly, start = start ) )
      call check( nf90_put_var(ncid_hourly, rrh_varid, rr_hourly, start = start ) )
      call check( nf90_put_var(ncid_hourly, assth_varid, asst_hourly, start = start ) )
      call check( nf90_put_var(ncid_hourly, assgh_varid, assg_hourly, start = start ) )
      call check( nf90_put_var(ncid_hourly, etmth_varid, etmt_hourly, start = start ) )
      call check( nf90_put_var(ncid_hourly, etmgh_varid, etmg_hourly, start = start ) )
      call check( nf90_put_var(ncid_hourly, su1h_varid, su1_hourly, start = start ) )
      call check( nf90_put_var(ncid_hourly, zwh_varid, zw_hourly, start = start ) )
      call check( nf90_put_var(ncid_hourly, wsh_varid, ws_hourly, start = start ) )
      call check( nf90_put_var(ncid_hourly, spgfcfh_varid, spgfcf_hourly, start = start ) )
      call check( nf90_put_var(ncid_hourly, infxh_varid, infx_hourly, start = start ) )


      !*****************************
      !root water uptake

      count_ruptkt = (/ 1, 1, s_maxlayer, 1 /)
      start_ruptkt = (/ 1, 1, 1, num_hour_tot-hourdiff /)

      !add time
      call check( nf90_put_var(ncid_ruptkt, time_ruptkt_varid, (starthour + num_hour_tot - hourdiff), &
                    start= (/ num_hour_tot-hourdiff/)  ) )

      !add results root water uptake
      call check( nf90_put_var(ncid_ruptkt, ruptkt_varid, ruptkt_hourly, start = start_ruptkt, count=count_ruptkt ) )



      !*****************************
      !soil moisture 
      count_suhourly = (/ 1, 1, s_maxlayer, 1 /)
      start_suhourly = (/ 1, 1, 1, num_hour_tot-hourdiff /)

     !add time
      call check( nf90_put_var(ncid_suhourly, time_suhourly_varid, (starthour + num_hour_tot - hourdiff), &
            start= (/ num_hour_tot-hourdiff /)  ) )


      !add results soil moisture

      call check( nf90_put_var(ncid_suhourly, suhourly_varid, su_hourly, start = start_suhourly, count=count_suhourly ) )


     end if

     else


         if (year .ge. i_firstyear .and. year .le. i_lastyear) then
!          * internal write to convert from number to string
          write(str,'(i3)') wlayer_
!         * includes a column for each sublayer
          hourlyformat = '(I6,I6,I4,I7,I5,'//str//'E14.6)'

          write(kfile_resultshourly,'(I6,I7,I7,I7,I7,24E15.5)')          &
          &    year, month, day, num_day, num_hour,          &
          &    rain_hourly, tair_hourly, par_hourly, gstomt_hourly, gstomg_hourly, vd_hourly, esoil_hourly,    &
          &    pc_hourly, jmax25t_hourly, jmax25g_hourly, mqt_hourly,          &
          &    rl_hourly, lambdat_hourly, lambdag_hourly, rr_hourly,  &
          &    asst_hourly, assg_hourly, etmt_hourly, etmg_hourly, su1_hourly, zw_hourly, ws_hourly, &
          &    spgfcf_hourly, infx_hourly

          write(kfile_delzhourly,hourlyformat) fyear(nday),            &
          &      fmonth(nday), fday(nday), nday, nhour, s_delz(1:wlayer_)

          write(kfile_ruptkthourly,hourlyformat) fyear(nday),          &
          &      fmonth(nday), fday(nday), nday, nhour, ruptkt_h(1:wlayer_)

          write(kfile_suhourly,hourlyformat) fyear(nday),              &
          &      fmonth(nday), fday(nday), nday, nhour, su__(1:wlayer_)
         endif


      end if

      return
      end subroutine vom_write_hourly



  subroutine check(status)
    use netcdf
    integer, intent ( in) :: status


    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  


  subroutine get_julday(yy, mm, dd, jday)

     implicit none
     integer, intent(in) :: yy
     integer, intent(in) :: mm
     integer, intent(in) :: dd
     integer, intent(out) :: jday

     integer :: part1
     integer :: part2
     integer :: part3    



       part1 = dd-32075+1461*(yy + 4800 + (mm-14)/12)/4 
       part2 = 367*(mm - 2 - (mm-14)/12*12)
       part3 = -3*((yy+4900+(mm-14)/12)/100)/4

       jday = part1 + part2 + part3

  end subroutine get_julday




