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



      subroutine vom_open_output_nc ()
      use vom_vegwat_mod
      use netcdf
      implicit none

     integer :: ncid
     CHARACTER*100                  :: filename 
     integer :: lon_dimid
     integer :: lon_varid
     integer :: lat_dimid
     integer :: lat_varid

     integer :: time_dimid
     integer :: n_lat = 1
     integer :: n_lon = 1

     filename = trim(adjustl(i_outputpath))// &
           trim(adjustl(sfile_resultsdaily))


     ! Create the file. 
     call check( nf90_create(filename, nf90_clobber, ncid) )

     ! Define the dimensions. NetCDF will hand back an ID for each. 
     call check( nf90_def_dim(ncid, "lat", n_lat, lat_dimid) )
     call check( nf90_def_dim(ncid, "lon", n_lon, lon_dimid) )
     call check( nf90_def_dim(ncid, "time", n_lon, NF90_UNLIMITED, time_dimid) )


     call check( nf90_def_var(ncid, "lat", NF90_REAL, lat_dimid, lat_varid) )
     call check( nf90_def_var(ncid, "lon", NF90_REAL, lon_dimid, lon_varid) )

     ! Assign units attributes to coordinate variables.
     call check( nf90_put_att(ncid, lat_varid, "units", "degrees_north") )
     call check( nf90_put_att(ncid, lon_varid, "units", "degrees_south") )

     ! Define the netCDF variables.
     call check( nf90_def_var(ncid, 'rain', NF90_REAL, dimids, rain_varid) )
     call check( nf90_def_var(ncid, 'tairmax', NF90_REAL, dimids, tairmax_varid) )
     call check( nf90_def_var(ncid, 'tairmin', NF90_REAL, dimids, tairmin_varid) )
     call check( nf90_def_var(ncid, 'par', NF90_REAL, dimids, par_varid) )
     call check( nf90_def_var(ncid, 'vd', NF90_REAL, dimids, vd_varid) )
     call check( nf90_def_var(ncid, 'esoil', NF90_REAL, dimids, esoil_varid) )
     call check( nf90_def_var(ncid, 'jmax25t', NF90_REAL, dimids, jmax25t_varid) )
     call check( nf90_def_var(ncid, 'jmax25g', NF90_REAL, dimids, jmax25g_varid) )
     call check( nf90_def_var(ncid, 'pc', NF90_REAL, dimids, pc_varid) )
     call check( nf90_def_var(ncid, 'rlt', NF90_REAL, dimids, rlt_varid) )
     call check( nf90_def_var(ncid, 'rlg', NF90_REAL, dimids, rlg_varid) )
     call check( nf90_def_var(ncid, 'lambdat', NF90_REAL, dimids, lambdat_varid) )
     call check( nf90_def_var(ncid, 'lambdag', NF90_REAL, dimids, lambdag_varid) )
     call check( nf90_def_var(ncid, 'rrt', NF90_REAL, dimids, rrt_varid) )
     call check( nf90_def_var(ncid, 'rrg', NF90_REAL, dimids, rrg_varid) )
     call check( nf90_def_var(ncid, 'asst', NF90_REAL, dimids, asst_varid) )
     call check( nf90_def_var(ncid, 'assg', NF90_REAL, dimids, assg_varid) )
     call check( nf90_def_var(ncid, 'su_avg', NF90_REAL, dimids, su_avg_varid) )
     call check( nf90_def_var(ncid, 'zw', NF90_REAL, dimids, zw_varid) )
     call check( nf90_def_var(ncid, 'ws', NF90_REAL, dimids, ws_varid) )
     call check( nf90_def_var(ncid, 'spgfcf', NF90_REAL, dimids, spgfcf_varid) )
     call check( nf90_def_var(ncid, 'infx', NF90_REAL, dimids, infx_varid) )
     call check( nf90_def_var(ncid, 'etmt', NF90_REAL, dimids, etmt_varid) )
     call check( nf90_def_var(ncid, 'etmg', NF90_REAL, dimids, etmg_varid) )
     call check( nf90_def_var(ncid, 'su_1', NF90_REAL, dimids, su_1_varid) )
     call check( nf90_def_var(ncid, 'topt', NF90_REAL, dimids, topt_varid) )
     call check( nf90_def_var(ncid, 'tcg', NF90_REAL, dimids, tcg_varid) )
     call check( nf90_def_var(ncid, 'tct', NF90_REAL, dimids, tct_varid) )
     call check( nf90_def_var(ncid, 'cpccg_d', NF90_REAL, dimids, cpccg_d_varid) )
     call check( nf90_def_var(ncid, 'cpcct_d', NF90_REAL, dimids, cpcct_d_varid) )
     call check( nf90_def_var(ncid, 'lai_t', NF90_REAL, dimids, lai_t_varid) )
     call check( nf90_def_var(ncid, 'lai_g', NF90_REAL, dimids, lai_g_varid) )
     call check( nf90_def_var(ncid, 'ncp_g', NF90_REAL, dimids, ncp_g_varid) )
     call check( nf90_def_var(ncid, 'ncp_t', NF90_REAL, dimids, ncp_t_varid) )

    ! Assign units attributes to the netCDF variables.
     call check( nf90_put_att(ncid, rain_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, tairmax_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, tairmin_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, par_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, vd_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, esoil_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, jmax25t_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, jmax25g_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, pc_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, rlt_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, rlg_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, lambdat_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, lambdag_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, rrt_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, rrg_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, asst_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, assg_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, su_avg_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, zw_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, ws_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, spgfcf_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, infx_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, etmt_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, etmg_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, su_1_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, topt_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, tcg_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, tct_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, cpccg_d_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, cpcct_d_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, lai_t_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, lai_g_varid, "units", "m/d") )
     call check( nf90_put_att(ncid, ncp_g_varid, "units", "m/d") )  
     call check( nf90_put_att(ncid, ncp_t_varid, "units", "m/d") )  

     ! End define mode.
     call check( nf90_enddef(ncid) )

     ! Write the coordinate variable data. This will put the latitudes
     ! and longitudes of our data grid into the netCDF file.
     call check( nf90_put_var(ncid, lat_varid, i_lat) )
     !call check( nf90_put_var(ncid, lon_varid, i_lon) )


      return
      end subroutine vom_open_output_nc









