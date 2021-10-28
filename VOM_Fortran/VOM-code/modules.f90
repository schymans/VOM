!***********************************************************************
!        Optimised Vegetation Optimality Model (VOM)
!        File and parameter definitions for the VOM
!        Original code coming from: https://github.com/schymans/VOM
!-----------------------------------------------------------------------
!        Author: 
!           Stan Schymanski
!           Now at: LIST, Luxembourg Institute of Science and Technology,
!                Belvaux, Luxembourg
!    
!        Contributors: 
!           Remko Nijzink
!           Now at: LIST, Luxembourg Institute of Science and Technology,
!                Belvaux, Luxembourg
!
!        Version: 
!-----------------------------------------------------------------------
!
!
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


      module vom_file_mod

!     * file codes

      INTEGER :: kfile_dailyweather  = 101
      INTEGER :: kfile_hourlyweather = 102

      INTEGER :: kfile_namelist      = 401
      INTEGER :: kfile_outputlist    = 402

      INTEGER :: kfile_resultshourly = 201
      INTEGER :: kfile_resultsdaily  = 202
      INTEGER :: kfile_resultsyearly = 203
      INTEGER :: kfile_rsurfdaily    = 204
      INTEGER :: kfile_delzhourly    = 205
      INTEGER :: kfile_ruptkthourly  = 206
      INTEGER :: kfile_suhourly      = 207
      INTEGER :: kfile_soilprofile   = 208
      INTEGER :: kfile_model_output  = 209

      INTEGER :: kfile_random_output  = 210

      INTEGER :: kfile_vd_d           = 311 
      INTEGER :: kfile_esoil          = 312 
      INTEGER :: kfile_jmax25t        = 313
      INTEGER :: kfile_jmax25g        = 314
      INTEGER :: kfile_vegcov         = 315
      INTEGER :: kfile_resp           = 316
      INTEGER :: kfile_lambdat        = 317
      INTEGER :: kfile_lambdag        = 318
      INTEGER :: kfile_rrt            = 319
      INTEGER :: kfile_rrg            = 320
      INTEGER :: kfile_asst           = 321
      INTEGER :: kfile_assg           = 322
      INTEGER :: kfile_su_av          = 323
      INTEGER :: kfile_zw             = 324
      INTEGER :: kfile_wsnew          = 325
      INTEGER :: kfile_spgfcf         = 326
      INTEGER :: kfile_infx           = 327
      INTEGER :: kfile_etmt           = 328
      INTEGER :: kfile_etmg           = 329
      INTEGER :: kfile_su1            = 330
      INTEGER :: kfile_topt           = 331
      INTEGER :: kfile_random_params  = 332
      INTEGER :: kfile_perc_cov       = 333



      INTEGER :: kfile_sceout        = 701
      INTEGER :: kfile_progress      = 702
      INTEGER :: kfile_lastloop      = 703
      INTEGER :: kfile_lastbest      = 704
      INTEGER :: kfile_bestpars      = 705
      INTEGER :: kfile_beststat      = 706
      INTEGER :: kfile_pars          = 707

!     * file names

      CHARACTER(len=*),parameter :: sfile_dailyweather  = 'dailyweather.prn'
      CHARACTER(len=*),parameter :: sfile_hourlyweather = 'hourlyweather.prn'

      CHARACTER*100              :: sfile_namelist      = 'vom_namelist'
      CHARACTER(len=*),parameter :: sfile_outputlist    = 'output_namelist'
      CHARACTER(len=*),parameter :: sfile_resultshourly = 'results_hourly.txt'
      CHARACTER(len=*),parameter :: sfile_resultsdaily  = 'results_daily.txt'
      CHARACTER(len=*),parameter :: sfile_resultsyearly = 'results_yearly.txt'
      CHARACTER(len=*),parameter :: sfile_rsurfdaily    = 'rsurf_daily.txt'
      CHARACTER(len=*),parameter :: sfile_delzhourly    = 'delz_hourly.txt'
      CHARACTER(len=*),parameter :: sfile_ruptkthourly  = 'ruptkt_hourly.txt'
      CHARACTER(len=*),parameter :: sfile_suhourly      = 'su_hourly.txt'
      CHARACTER(len=*),parameter :: sfile_soilprofile   = 'soilprofile.par'
      CHARACTER(len=*),parameter :: sfile_model_output  = 'model_output.txt'

      CHARACTER(len=*),parameter :: nfile_resultshourly = 'results_hourly.nc'
      CHARACTER(len=*),parameter :: nfile_resultsdaily  = 'results_daily.nc'
      CHARACTER(len=*),parameter :: nfile_resultsyearly = 'results_yearly.nc'
      CHARACTER(len=*),parameter :: nfile_rsurfdaily    = 'rsurf_daily.nc'
      CHARACTER(len=*),parameter :: nfile_delzhourly    = 'delz_hourly.nc'
      CHARACTER(len=*),parameter :: nfile_ruptkthourly  = 'ruptkt_hourly.nc'
      CHARACTER(len=*),parameter :: nfile_suhourly      = 'su_hourly.nc'

      CHARACTER(len=*),parameter :: sfile_random_output  = 'random_ncp.txt'
      CHARACTER(len=*),parameter :: sfile_random_params  = 'random_params.txt'

      CHARACTER(len=*),parameter :: sfile_vd_d           = 'vpd.txt'
      CHARACTER(len=*),parameter :: sfile_esoil          = 'esoil.txt'
      CHARACTER(len=*),parameter :: sfile_jmax25t        = 'jmax25t.txt'
      CHARACTER(len=*),parameter :: sfile_jmax25g        = 'jmax25g.txt'
      CHARACTER(len=*),parameter :: sfile_vegcov         = 'veg_cover.txt'
      CHARACTER(len=*),parameter :: sfile_resp           = 'leaf_resp.txt'
      CHARACTER(len=*),parameter :: sfile_lambdat        = 'lambdat.txt'
      CHARACTER(len=*),parameter :: sfile_lambdag        = 'lambdag.txt'
      CHARACTER(len=*),parameter :: sfile_rrt            = 'root_resp_t.txt'
      CHARACTER(len=*),parameter :: sfile_rrg            = 'root_resp_g.txt'
      CHARACTER(len=*),parameter :: sfile_asst           = 'asst.txt'
      CHARACTER(len=*),parameter :: sfile_assg           = 'assg.txt'
      CHARACTER(len=*),parameter :: sfile_su_av          = 'su_av.txt'
      CHARACTER(len=*),parameter :: sfile_zw             = 'wat_table.txt'
      CHARACTER(len=*),parameter :: sfile_wsnew          = 'soil_wat_storage.txt'
      CHARACTER(len=*),parameter :: sfile_spgfcf         = 'seepage.txt'
      CHARACTER(len=*),parameter :: sfile_infx           = 'infilt.txt'
      CHARACTER(len=*),parameter :: sfile_etmt           = 'etmt.txt'
      CHARACTER(len=*),parameter :: sfile_etmg           = 'etmg.txt'
      CHARACTER(len=*),parameter :: sfile_su1            = 'su1.txt'
      CHARACTER(len=*),parameter :: sfile_topt           = 'temp_opt.txt'
      CHARACTER(len=*),parameter :: sfile_perc_cov       = 'perc_cov.txt'

      CHARACTER(len=*),parameter :: sfile_sceout        = 'sce_out.txt'
      CHARACTER(len=*),parameter :: sfile_progress      = 'sce_progress.txt'
      CHARACTER(len=*),parameter :: sfile_lastloop      = 'sce_lastloop.txt'
      CHARACTER(len=*),parameter :: sfile_lastbest      = 'sce_lastbest.txt'
      CHARACTER(len=*),parameter :: sfile_bestpars      = 'sce_bestpars.txt'
      CHARACTER(len=*),parameter :: sfile_beststat      = 'sce_status.txt'
      CHARACTER(len=*),parameter :: sfile_pars          = 'pars.txt'

      CHARACTER*100  :: i_outputpath     ! Path for model outputs
      CHARACTER*100  :: i_inputpath      ! Path to model inputs



     integer :: ncid
     integer :: ncid_rsurf
     integer :: ncid_hourly
     integer :: ncid_yearly
     integer :: ncid_ruptkt
     integer :: ncid_suhourly


     integer :: rain_varid
     integer :: tairmax_varid
     integer :: tairmin_varid
     integer :: par_varid
     integer :: vd_varid
     integer :: esoil_varid
     integer :: jmax25t_varid
     integer :: jmax25g_varid
     integer :: pc_varid
     integer :: rlt_varid
     integer :: rlg_varid
     integer :: lambdat_varid
     integer :: lambdag_varid
     integer :: rrt_varid
     integer :: rrg_varid
     integer :: asst_varid
     integer :: assg_varid
     integer :: su_avg_varid
     integer :: zw_varid
     integer :: ws_varid
     integer :: spgfcf_varid
     integer :: infx_varid
     integer :: etmt_varid
     integer :: etmg_varid
     integer :: su_1_varid
     integer :: topt_varid
     integer :: tcg_varid
     integer :: tct_varid
     integer :: cpccg_d_varid
     integer :: cpcct_d_varid
     integer :: lai_t_varid
     integer :: lai_g_varid
     integer :: ncp_g_varid
     integer :: ncp_t_varid
     integer :: rsurf_varid

     integer :: rainh_varid
     integer :: tairh_varid
     integer :: parh_varid
     integer :: vdh_varid
     integer :: esoilh_varid
     integer :: jmax25th_varid
     integer :: jmax25gh_varid
     integer :: pch_varid
     integer :: rlh_varid
     integer :: mqth_varid
     integer :: lambdath_varid
     integer :: lambdagh_varid
     integer :: rrh_varid
     integer :: assth_varid
     integer :: assgh_varid
     integer :: su1h_varid
     integer :: zwh_varid
     integer :: wsh_varid
     integer :: spgfcfh_varid
     integer :: infxh_varid
     integer :: etmth_varid
     integer :: etmgh_varid
 
     integer :: rainy_varid
     integer :: pary_varid
     integer :: vdy_varid
     integer :: srady_varid
     integer :: esoily_varid
     integer :: etmgy_varid
     integer :: assgy_varid
     integer :: rlgy_varid
     integer :: rrgy_varid
     integer :: cpccgy_varid
     integer :: tcgy_varid
     integer :: etmty_varid
     integer :: assty_varid
     integer :: rlty_varid
     integer :: rrty_varid
     integer :: cpccty_varid
     integer :: tcty_varid

     integer :: ruptkt_varid
     integer :: suhourly_varid

     integer :: time_varid
     integer :: time_rsurf_varid
     integer :: time_hourly_varid
     integer :: time_yearly_varid
     integer :: time_ruptkt_varid
     integer :: time_suhourly_varid

     integer :: startday
     integer :: starthour
     integer :: startyear
     integer :: hourdiff

      end module vom_file_mod

!     ******************************************************************
!     * Module defining variables and parameters for the vegetation
!     * model (transpmodel).
!     ******************************************************************

      module vegmod
      implicit none

      INTEGER :: optmode                ! Indicator of optimisation mode
      REAL*8  :: time                   ! Seconds of hour
      REAL*8  :: error                  ! Cumulative error in water balance
      INTEGER :: finish                 ! flag to finish all loops

      REAL*8, ALLOCATABLE :: output_mat(:,:)

      REAL*8, PARAMETER :: p_a     = 1.6d0        ! Ratio of diffusivities of water vapour to CO2 in air
      REAL*8, PARAMETER :: p_pi    = 3.14159d0    ! Pi-constant
      REAL*8, PARAMETER :: p_mpbar = 10.2d0       ! Conversion factor from MPa to bar
      REAL*8, PARAMETER :: p_E     = 2.7182818d0  ! Eurler's number
      REAL*8, PARAMETER :: p_R_    = 8.314d0      ! Molar gas konstant
      REAL*8, PARAMETER :: l_E_    = 2.45d0       ! Latent heat of vaporization (MJ/kg)
      REAL*8, PARAMETER :: srad2par_h = 2.0699d0  ! Conversion from srad to par hourly (mol/MJ)
      REAL*8, PARAMETER :: srad2par_d = 2.0804d0  ! Conversion from srad to par daily (mol/MJ)
      REAL*8, PARAMETER :: rho_wat = 1000.0d0     ! Density of water (kg/m3)
      REAL*8, PARAMETER :: Gsc     = 0.0820d0     ! Solar constant (MJ m-2 day-1)
      REAL*8, PARAMETER :: X0      = 0.26d0       ! parameter to split PAR          
      REAL*8, PARAMETER :: Y0      = 0.96d0       ! parameter to split PAR           
      REAL*8, PARAMETER :: Y1      =  0.05d0      ! parameter to split PAR      
      
      INTEGER :: nyear                  ! Year
      INTEGER :: nday                   ! Day since start of run
      INTEGER :: nhour                  ! Hour of day
      INTEGER :: th_                    ! Hour since start of run
      INTEGER :: c_testday              ! Number of days for initial check if netass>0

      INTEGER, ALLOCATABLE :: fyear(:)    ! Year for each day
      INTEGER, ALLOCATABLE :: fmonth(:)   ! Month for each day
      INTEGER, ALLOCATABLE :: fday(:)     ! Day of month
      INTEGER, ALLOCATABLE :: dayyear(:)  ! Day of year

!     * climate

      REAL*8, ALLOCATABLE :: phi_zenith(:)  ! zenith angle
      
      REAL*8, ALLOCATABLE :: tair_h(:)    ! Hourly air temperature (K)
      REAL*8, ALLOCATABLE :: tairmin_d(:) ! Daily minimum temperature (K)
      REAL*8, ALLOCATABLE :: tairmax_d(:) ! Daily maximum temperature (K)

      REAL*8              :: topt_      ! Optimal temperature in temperature response curve (K)

      REAL*8, ALLOCATABLE :: press_d(:) ! Daily air pressure (Pa)
      REAL*8, ALLOCATABLE :: press_h(:) ! Hourly air pressure (Pa)

      REAL*8, ALLOCATABLE :: par_h(:)   ! Hourly photosynthetically active radiation (mol/m2/s)
      REAL*8, ALLOCATABLE :: par_d(:)   ! Daily photosynthetically active radiation (mol/m2/d)
      REAL*8              :: par_y      ! Annual photosynthetically active radiation (mol/m2/y)
      
      REAL*8, ALLOCATABLE :: pardiff_h(:) ! Hourly diffuse photosynthetically active radiation (mol/m2/s)
      REAL*8, ALLOCATABLE :: pardir_h(:)  ! Hourly direct photosynthetically active radiation (mol/m2/s)
      
      REAL*8, ALLOCATABLE :: par_et_h(:)  ! Hourly extraterrestrial radiation (mol/m2/s)
      
      REAL*8, ALLOCATABLE :: srad_d(:)  ! Daily shortwave radiation  (MJ/m2/d)
      REAL*8              :: srad_y     ! Annual shortwave radiation (MJ/m2/y)

      REAL*8, ALLOCATABLE :: ca_h(:)    ! Hourly atmospheric CO2 mole fraction
      REAL*8, ALLOCATABLE :: ca_d(:)    ! Daily atmospheric CO2 mole fraction

      REAL*8, ALLOCATABLE :: vp_d(:)    ! Daily absolute vapour pressure (hPa)
      REAL*8, ALLOCATABLE :: vp_h(:)    ! Hourly absolute vapour pressure (hPa)

      REAL*8, ALLOCATABLE :: vd_h(:)    ! Hourly atmospheric vapour deficit (VPD/air pressure) (Pa/Pa)
      REAL*8              :: vd_d       ! Mean daily atmospheric vapour deficit (mol/mol)
      REAL*8              :: vd_y       ! Mean annual atmospheric vapour deficit (mol/mol)

      REAL*8, ALLOCATABLE :: rain_h(:)  ! Hourly rainfall rate (m/s)
      REAL*8, ALLOCATABLE :: rain_d(:)  ! Daily rainfall (mm/d)
      REAL*8              :: rain_y     ! Annual rainfall (mm/y)

!     * soil

      REAL*8, ALLOCATABLE :: c_hhydrst(:)  ! Hydrostatic head in each layer relative to soil surface (m)

      REAL*8  :: gammastar              ! CO2 compensation point (mol/mol)

      REAL*8  :: wsnew                  ! Total soil water store at next time step (m)
      REAL*8  :: wsold                  ! Previous total soil water storage (m)

      REAL*8  :: o_cait                 ! Crown area index perennial vegetation (0-1)
      REAL*8  :: caig_d(3)              ! Crown area index seasonal vegetation (caig_d(2) is actual value)
      REAL*8  :: c_caigmin              ! Minimum grass crown area index; initial point for growth (-)

      REAL*8 :: fpar_lg(3)                !local fraction of absorbed radiation grasses (-)
      REAL*8 :: fpar_lt(3)                !local fraction of absorbed radiation trees (-)
     
      REAL*8 :: fpard_lg                !mean local fraction of absorbed radiation grasses per day (-)
      REAL*8 :: fpard_lt                !mean local fraction of absorbed radiation trees pea day (-)
      
      REAL*8 :: frac_sung(3)            !fraction sunlit leaves grasses (-)
      REAL*8 :: frac_shadeg(3)          !fraction shaded leaves grasses (-)

      REAL*8 :: frac_sunt(3)            !fraction sunlit leaves trees (-)
      REAL*8 :: frac_shadet(3)          !fraction shaded leaves trees (-)
      
!     * leaf

      REAL*8  :: o_wstexp               ! Exponent for calculating lambdat_d (-)
      REAL*8  :: o_wsgexp               ! Exponent for calculating lambdag (-)
      REAL*8  :: o_lambdatf             ! Factor for calculating lambdat_d (mol/mol/m)
      REAL*8  :: o_lambdagf             ! Factor for calculating lambdag_d (mol/mol/m)
      REAL*8  :: lambdat_d              ! Target dE/dA for calculating gstomt (mol/mol)
      REAL*8  :: lambdag_d              ! Target dE/dA for calculating gstomg (mol/mol)
      REAL*8  :: gstomt                 ! Tree stomatal conductance (mol/m2/s)
      REAL*8  :: gstomg(3,3)            ! Grass stomatal conductance (mol/m2/s)

      REAL*8  :: rlt_h(3,3)             ! Tree leaf respiration for different values of Jmax (rlt_h(2) is actual value) (mol/h)
      REAL*8  :: rlt_d                  ! Daily tree leaf respiration (mol/d)
      REAL*8  :: rlt_y                  ! Annual tree leaf respiration (mol/y)
      REAL*8  :: rlg_h(3,3)           ! Grass leaf respiration (mol/h)
      REAL*8  :: rlg_d                  ! Daily grass leaf respiration (mol/d)
      REAL*8  :: rlg_y                  ! Annual grass leaf respiration (mol/y)

      REAL*8  :: transpt                ! Tree transpiration rate (mol/m2/s)
      REAL*8  :: transpg(3,3)         ! Grass transpiration rate (mol/m2/s)

      REAL*8  :: q_tct_d(3)             ! Tree foliage turnover costs (mol/m2/s)
      REAL*8  :: tct_y                  ! Annual tree foliage turnover costs (mol/m2/y)
      REAL*8  :: tcg_d(3, 3)            ! Grass foliage turnover costs (mol/m2/s)
      REAL*8  :: tcg_y                  ! Annual grass foliage turnover costs (mol/m2/y)

      REAL*8  :: jactt(3,3)             ! Electron transport rates for different values of Jmax (jactt(2) is actual value) (mol/m2/s)
      REAL*8  :: jactg(3,3)           ! Grass electron transport rate (mol/m2/s)

      REAL*8  :: jmaxt_h(3)             ! Tree photosynthetic electron transport capacity (mol/m2/s)
      REAL*8  :: jmaxg_h(3)             ! Grass electron transport capacity (mol/m2/s)

      REAL*8  :: jmax25t_d(3)           ! Tree photosynthetic electron transport capacity at 25oC (mol/m2/s)
      REAL*8  :: jmax25g_d(3)           ! Grass photosynthetic electron transport capacity at 25oC (mol/m2/s)
 
      REAL*8  :: lai_lt(3)              ! Local leaf area index trees (-)
      REAL*8  :: lai_lg(3)              ! Local leaf area index grasses (-)

!     * plant water

      REAL*8  :: asst_h(3,3)            ! Tree hourly assimilation rate for different values of Jmax (asst_h(2) is actual value) (mol/m2/h)
      REAL*8  :: asst_d(3,3)            ! Daily tree assimilation (mol/m2/d)
      REAL*8  :: asst_y                 ! Annual tree assimilation (mol/m2/y)
      REAL*8  :: assg_h(3,3,3)          ! Hourly grass assimilation (mol/m2/h)
      REAL*8  :: assg_d(3,3,3)          ! Daily grass assimilation (mol/m2/d)
      REAL*8  :: assg_y                 ! Annual grass assimilation (mol/m2/y)

      REAL*8  :: q_cpcct_d              ! Tree water transport costs as a function of projected cover and rooting depth (mol/m2/s)
      REAL*8  :: cpcct_y                ! Annual tree water transport costs (mol/m2/y)
      REAL*8  :: cpccg_d(3)             ! Grass water transport costs (mol/m2/s)
      REAL*8  :: cpccg_y                ! Annual grass water transport costs (mol/m2/y)

      REAL*8  :: etmt__                 ! Transpiration rate (m/s)
      REAL*8  :: etmt_h                 ! Hourly transpiration (m/h)
      REAL*8  :: etmt_d                 ! Daily transpiration rate (m/d)
      REAL*8  :: etmt_y                 ! Annual tree transpiration (mm/y)
      REAL*8  :: etmg__(3,3)          ! Grass transpiration rate (m/s)
      REAL*8  :: etmg_h                 ! Hourly grass transpiration (m/h)
      REAL*8  :: etmg_d                 ! Daily grass transpiration (m/d)
      REAL*8  :: etmg_y                 ! Annual grass transpiration (mm/y)
      REAL*8  :: etm_y                  ! Annual total transpiration (mm/y)

      REAL*8  :: mqt_                   ! Tree water content (kg/m2)
      REAL*8  :: mqtnew                 ! Tree water content in next time step (kg/m2)
      REAL*8  :: mqtold                 ! Previous tree water content (kg/m2)
      REAL*8  :: dmqt                   ! Rate of change in tree water content (kg/m2/s)
      REAL*8  :: q_mqx                  ! Tree maximum water content per ground area (kg/m2)
      REAL*8  :: mqsst_                 ! Tree water content at steady state (kg/m2)
      REAL*8  :: mqsstmin               ! Tree water content at turgor loss point (kg/m2)

      REAL*8  :: q_md                     ! Tree dry mass per unit ground area (kg/m2)
      REAL*8  :: o_mdstore                ! Wood water storage parameter of trees  (kg/m2)

!     * roots

      REAL*8  :: o_rtdepth               ! Tree rooting depth (m)
      REAL*8  :: o_rgdepth               ! Grass rooting depth (m)

      INTEGER             :: pos_slt    ! Lowest soil layer containing tree roots 
      INTEGER             :: pos_slg    ! Lowest soil layer containing grass roots
      INTEGER             :: pos_ult    ! Lowest soil layer containing tree roots within unsaturated zone
      INTEGER             :: pos_ulg    ! Lowest soil layer containing grass roots within unsaturated zone

      REAL*8              :: changef    ! Change factor for adjusting root surface area (-)

      REAL*8, ALLOCATABLE :: rsurft_(:)    ! Root surface area of trees in each layer (m2/m3)
      REAL*8, ALLOCATABLE :: rsurftnew(:)  ! Adjusted root surface area of trees in each layer for next day (m2/m3)
      REAL*8, ALLOCATABLE :: rsurfg_(:)    ! Root surface area of grasses in each layer (m2/m3)
      REAL*8, ALLOCATABLE :: rsurfgnew(:)  ! Adjusted root surface area of grasses in each layer for next day (m2/m3)

      REAL*8              :: rootlim(3,3,3)  ! Indicator whether root surface are was limiting root water uptake (-)

      REAL*8, ALLOCATABLE :: rsoil(:)   ! Resistance to water flow towards roots in each soil layer (s)

      REAL*8, ALLOCATABLE :: refft(:)   ! Relative root water uptake efficiency for trees in each layer (-)
      REAL*8, ALLOCATABLE :: reffg(:)   ! Relative root water uptake efficiency for grasses in each layer (-)
      INTEGER             :: posmna(3)  ! Pointer to variable values that achieved maximum net assimilation (-)
 
      REAL*8              :: rrt_d      ! Tree root respiration rate (mol/m2/s)
      REAL*8              :: rrt_y      ! Annual tree root respiration (mol/m2/y)
      REAL*8              :: rrg_d      ! Grass root respiration (mol/m2/s)
      REAL*8              :: rrg_y      ! Annual grass root respiration (mol/m2/y)

      REAL*8, ALLOCATABLE :: prootm(:)  ! Root hydraulic head in each layer (m)

      REAL*8              :: sumruptkt_h  ! Hourly total tree root water uptake  (m/h)
      REAL*8, ALLOCATABLE :: ruptkt__(:)  ! Root water uptake rate perennial veg (m/s)
      REAL*8, ALLOCATABLE :: ruptkt_h(:)  ! Hourly root water uptake by trees in each layer (m/h)
      REAL*8, ALLOCATABLE :: ruptkt_d(:)  ! Daily root water uptake by trees in each layer (m/d)
      REAL*8, ALLOCATABLE :: ruptkg__(:)  ! Root water uptake rate seasonal veg (m/s)
      REAL*8, ALLOCATABLE :: ruptkg_h(:)  ! Hourly root water uptake by grasses in each layer (m/h)
      REAL*8, ALLOCATABLE :: ruptkg_d(:)  ! Daily root water uptake by grasses in each layer (m/d)
      REAL*8, ALLOCATABLE :: perc_cov_veg(:)  ! Daily coverage of vegetation (-)


!     ****************************
!     * input parameters input.par
!     ****************************

      REAL*8  :: i_alpha                ! Initial slope of electron transport curve (-)
      REAL*8  :: i_cpccf                ! Water transport costs per m root depth and m^2 cover (mol/m^3/s)
      REAL*8  :: i_tcfg                 ! Turnover cost factor for foliage grasses (tc=i_tcf*LAI) (mol/m^2/s)
      REAL*8  :: i_tcft                 ! Turnover cost factor for foliage trees (tc=i_tcf*LAI) (mol/m^2/s)
      INTEGER :: i_maxyear              ! Number of years to process
      INTEGER :: i_testyear             ! Number of years after which to perform initial test of netass
      REAL*8  :: i_ha                   ! Temperature response parameter (J/mol)
      REAL*8  :: i_hd                   ! Temperature response parameter (J/mol)
      REAL*8  :: i_toptf                ! Parameter to calculate adaptation of topt (range 0-1 for no to full adaptation) (-)
      REAL*8  :: i_toptstart            ! Start temperature parameter for topt to calculate jmax (K)
      REAL*8  :: i_rlratio              ! Ratio of leaf respiration to photosynthetic capacity (-)

!     * Catchment parameters

      REAL*8  :: i_lat                 ! geogr. latitude
      REAL*8  :: i_lon                 ! geogr. longitude

!     * Soil parameters

!     * Vertical Resolution

!     * Vegetation Parameters


      REAL*8  :: i_mdtf                 ! Total dry mass of living tissues of trees per unit pc (kg/m^2)
      REAL*8  :: i_mqxtf                ! Total water storage capacity in living tissues of trees per unit md (kg/m^2)
      REAL*8  :: i_rrootm               ! Root water uptake resistivity in soil (s)
      REAL*8  :: i_rsurfmin             ! Minimum root area per m^3 to be maintained (m2/m3)
      REAL*8  :: i_rsurf_               ! Initial root surface area per m^3 (m2/m3)
      REAL*8  :: i_rootrad              ! Average fine root radius (m)
      REAL*8  :: i_prootmg              ! Constant root balance pressure of 1.5 MPa in grasses (m)
      REAL*8  :: i_growthmax            ! Parameter determining maximum daily growth increment of root surface area (m2/m3/d)
      REAL*8  :: i_incrcovg             ! parameter determining maximum increment percentage of grass cover (-)
      REAL*8  :: i_incrjmax             ! parameter determining maximum increment percentage of jmax25 (-)
      REAL*8  :: i_jmax_ini             ! parameter determining the start value of jmax25 (mol/m2/s)
      REAL*8  :: i_incrlait             ! parameter determining maximum increment percentage of lai trees (-)
      REAL*8  :: i_incrlaig             ! parameter determining maximum increment percentage of lai grasses
      REAL*8  :: i_chi_g                ! ratio projected areas of canopy elements on horizontal and vertical surfaces (-)
      REAL*8  :: i_chi_t                ! ratio projected areas of canopy elements on horizontal and vertical surfaces (-)
      REAL*8  :: i_alpha_abs                ! ratio projected areas of canopy elements on horizontal and vertical surfaces (-)
      REAL*8  :: i_trans_vegcov         ! fraction of radiative energy reaching soil under full cover (0-1) (-)

      INTEGER :: i_firstyear            ! First year for the generation of hourly output in computation mode
      INTEGER :: i_lastyear             ! Last year for the generation of hourly output in computation mode

      INTEGER :: i_write_h              ! Flag to write out hourly input values after conversation from daily values
      INTEGER :: i_read_pc              ! Flag to read vegetation cover as input
      LOGICAL :: i_write_nc             ! Flag to write to netcdf instead of plain text
      INTEGER :: i_lai_function         ! Switch to use 1) linear or 2) exponential LAI estimate, as function of cover
      INTEGER :: i_no_veg               ! Flag to switch vegetation off (1=no vegetation)


!     * Derived parameters

      REAL*8  :: c_epsln                ! Soil porosity (-)

      INTEGER :: c_maxhour              ! Number of hours to process (h)
      INTEGER :: c_maxday               ! Number of days to process (d)

      !$OMP threadprivate( time, error, finish, nyear, nday, nhour, th_, c_testday,   & 
      !$OMP topt_, par_y, srad_y,   &
      !$OMP vd_d, vd_y, rain_y, gammastar, wsnew, wsold, o_cait, caig_d, c_caigmin, &
      !$OMP o_wstexp, o_wsgexp, o_lambdatf, o_lambdagf, lambdat_d, lambdag_d, gstomt, gstomg, &
      !$OMP rlt_h, rlt_d, rlt_y, rlg_h, rlg_d, rlg_y, transpt, transpg, q_tct_d, tct_y, tcg_d, &
      !$OMP tcg_y, jactt, jactg, jmaxt_h, jmaxg_h, jmax25t_d, jmax25g_d, &
      !$OMP asst_h, asst_d, asst_y, assg_h, assg_d, assg_y, &
      !$OMP q_cpcct_d, cpcct_y, cpccg_d, cpccg_y, etmt__, etmt_h, etmt_d, etmt_y, etmg__, etmg_h, &
      !$OMP etmg_d, etmg_y, etm_y, mqt_, mqtnew, mqtold, dmqt, q_mqx, mqsst_, mqsstmin, q_md, &
      !$OMP o_mdstore, o_rtdepth, o_rgdepth, pos_slt, pos_slg, pos_ult, pos_ulg, changef, &
      !$OMP rootlim, posmna, &
      !$OMP ruptkt__, rsurft_, rsurftnew, prootm, ruptkt_d, ruptkt_h, ruptkg_h, ruptkg_d, &
      !$OMP refft, reffg, ruptkg__, rsurfg_, rsurfgnew, rsoil,      &  
      !$OMP rrt_d, rrt_y, rrg_d, rrg_y, sumruptkt_h, output_mat)


      end module vegmod

!     ******************************************************************
!     * Module defining variables and parameters for the water balance 
!     * model (watbal).
!     ******************************************************************

      module watmod
      implicit none

      INTEGER :: wlayer_                ! Number of soil layers in unsaturated zone
      INTEGER :: wlayernew              ! Number of layers in unsaturated zone for next time step

      REAL*8  :: dt_                     ! Length of time step (s)
      REAL*8  :: dtmax                   ! Maximum time step (s)
      INTEGER :: dtsu_count, dtmax_count ! For speed analysis, count of how often dtsu or dtmax is time limiting

      REAL*8  :: esoil__                ! Bare soil evaporation rate (m/s)
      REAL*8  :: esoil_h                ! Hourly soil evaporation (m/h)
      REAL*8  :: esoil_d                ! Daily soil evaporation rate (m/d)
      REAL*8  :: esoil_y                ! Annual soil evaporation rate (mm/y)

      REAL*8  :: spgfcf__               ! Seepage face flow rate (m/s)
      REAL*8  :: spgfcf_h               ! Hourly seepage face flow (m/h)
      REAL*8  :: spgfcf_d               ! Daily seepage face flow (m/d)

      REAL*8  :: inf__                  ! Infiltration rate (m/s)
      REAL*8  :: infx__                 ! Infiltration excess runoff rate (m/s)
      REAL*8  :: infx_h                 ! Hourly infiltration excess runoff (m/h)
      REAL*8  :: infx_d                 ! Daily infiltration excess runoff (m/d)

      REAL*8  :: zw_                    ! Elevation of water table (m)
      REAL*8  :: zwnew                  ! Elevation of water table at next time step (m)

      REAL*8  :: wc_                    ! Total soil water content (m)

      REAL*8, ALLOCATABLE :: cH2Ol_s(:) ! Soil water content in each layer (could also be called theta(:)) (m)

      REAL*8, ALLOCATABLE :: qbl(:)     ! Water flux across bottom boundary of each layer (positive upwards) (m/s)

      REAL*8, ALLOCATABLE :: pcap_(:)   ! Soil matric head (m)
      REAL*8, ALLOCATABLE :: pcapnew(:) ! Matric pressure head in each layer at next time step (m)

      REAL*8, ALLOCATABLE :: iovec(:)   ! Mass balance for each soil layer (m/s)
      REAL*8              :: io__       ! Input minus output of water in soil domain (m/s)
      REAL*8              :: io_h       ! Hourly input minus output of water in soil domain (m/d)
      REAL*8              :: ioacum     ! Cumulative mass balance of soil water (m)

      REAL*8, ALLOCATABLE :: su__(:)    ! Soil saturation degree in each layer (-)
      REAL*8, ALLOCATABLE :: sunew(:)   ! Soil saturation degree in each layer at next time step (-)
      REAL*8, ALLOCATABLE :: sueq(:)    ! Soil saturation degree above water table in hydrostatic equilibrium (-)
      REAL*8, ALLOCATABLE :: dsu(:)     ! Rate of change in soil saturation degree in each layer (1/s)

      REAL*8, ALLOCATABLE :: kunsat_(:)    ! Unsaturated hydraulic conductivity in each soil layer (m/s)
      REAL*8, ALLOCATABLE :: kunsatnew(:)  ! Unsaturated hydraulic conductivity in each layer at next time step (m/s)

      REAL*8, ALLOCATABLE :: s_ksat(:)    ! Saturated hydraulic conductivity (m/s)
      REAL*8, ALLOCATABLE :: s_thetas(:)  ! Saturated soil water content (-)
      REAL*8, ALLOCATABLE :: s_thetar(:)  ! Residual soil water content (-)
      REAL*8, ALLOCATABLE :: s_nvg(:)     ! Van Genuchten soil parameter n (-)
      REAL*8, ALLOCATABLE :: s_avg(:)     ! Van Genuchten soil parameter a (1/m)
      REAL*8, ALLOCATABLE :: c_mvg(:)     ! Van Genuchten soil parameter m (-)

!     ****************************
!     * input parameters input.par
!     ****************************

!     * Catchment parameters

      REAL*8  :: i_cgs                  ! Capital Gamma S (length scale for seepage outflow) (m)
      REAL*8  :: i_zr                   ! Average channel elevation (m)
      REAL*8  :: i_go                   ! Slope close to channel (radians)

!     * Soil parameters

      REAL*8  :: i_ksat                 ! Saturated hydraulic conductivity (m/s)
      REAL*8  :: i_thetar               ! Residual soil water content (-)
      REAL*8  :: i_thetas               ! Saturated soil water content (-)
      REAL*8  :: i_nvg                  ! Van Genuchten soil parameter n (-)
      REAL*8  :: i_avg                  ! Van Genuchten soil parameter a (1/m)
      REAL*8  :: i_mvg                  ! Van Genuchten soil parameter m (-)

      !$OMP threadprivate(wlayer_, wlayernew, dt_, dtmax, dtsu_count, dtmax_count, esoil__, esoil_h, &
      !$OMP i_cgs, i_zr, i_go, i_ksat, i_thetar, i_thetas, i_nvg, i_avg, &
      !$OMP esoil_d, esoil_y, spgfcf__, spgfcf_h, spgfcf_d, inf__, infx__, infx_h, infx_d, &
      !$OMP pcap_, su__, sunew, kunsat_, qbl, dsu, pcapnew, kunsatnew, sueq, cH2Ol_s, iovec,   &
      !$OMP zw_, zwnew, wc_, io__, io_h, ioacum)

      end module watmod

!     ******************************************************************
!     * Module defining variables and parameters common to both the
!     * vegetation model (transpmodel) and water balance model (watbal).
!     ******************************************************************

      module vom_vegwat_mod
      use vom_file_mod
      use vegmod
      use watmod
      implicit none

      INTEGER :: s_maxlayer             ! Number of soil layers

      REAL*8  :: i_cz                   ! Average soil elevation (m)

      REAL*8  :: i_delz                 ! Thickness of each soil layer (m)
      REAL*8, ALLOCATABLE :: s_delz(:)  ! Thickness of each soil layer (m)

      !$OMP threadprivate( i_cz )


      end module vom_vegwat_mod

!     ******************************************************************
!     * Module defining parameters and variables for
!     * SHUFFLED COMPLEX EVOLUTION (sce)
!     * Parameter optimisation algorithm based on a paper by Duan et al.
!     * (1993, J. Opt. Theory and Appl., 76, 501--521).
!     ******************************************************************

      module vom_sce_mod
      use vom_file_mod
      implicit none

      INTEGER :: n_thread               ! Number of threads to be used

      INTEGER :: success                ! Indicator wheter optimisation ended successfully

      INTEGER :: ncomp2                 ! Number of complexes
      INTEGER :: nopt                   ! Number of optimised parameters
      INTEGER :: mopt                   ! Minimum number of runs per complex
      INTEGER :: sopt                   ! SCE variable s
      INTEGER :: qopt                   ! SCE variable q
      INTEGER :: nrun                   ! Number of runs performed so far
      INTEGER :: nloop                  ! Number of loops performed so far
      INTEGER :: nsincebest             ! Number of runs since last improvement in objective function

      INTEGER, ALLOCATABLE :: optid(:)  ! Pointer to optimisable parameters

      REAL*8  :: ranscal                ! Scalar random number
      REAL*8  :: worstbest              ! Best OF of the worst complex in worstbest for assessment of gene pool mixing
      REAL*8  :: bestobj                ! Best OF value
      REAL*8  :: bestincomp             ! Best OF value in complex

      REAL*8, ALLOCATABLE :: wgt(:)     ! Probability weights for each optimisable parameter
      REAL*8, ALLOCATABLE :: cv_(:)     ! Coefficient of variation for each optimisable parameter in last loop
      REAL*8, ALLOCATABLE :: ranarr(:)  ! Array of random numbers
      !REAL*8, ALLOCATABLE :: ranarr_simplex(:)  ! Array of random numbers

      REAL*8, ALLOCATABLE :: shufflevar(:,:)  ! Population of parameter sets
      REAL*8, ALLOCATABLE :: ofvec(:)   ! Population of objective function values related to shufflevar

      CHARACTER(12)  :: evolution       ! Character string explaining how parameter set was generated
      CHARACTER(60)  :: outformat, loopformat

!     * allocated variables for sce()
      REAL*8, ALLOCATABLE :: sumvar(:)

!     * allocated variables for initialseed()
      INTEGER, ALLOCATABLE :: posarray(:,:)
      REAL*8,  ALLOCATABLE :: initpop(:,:)

!     * allocated variables for optsensitivity()
      REAL*8, ALLOCATABLE :: dataarray(:,:)  ! Parameter sets and objective functions used for sensitivity analysis
      REAL*8, ALLOCATABLE :: shufflevar2(:)  ! Parameter sets used for sensitivity analysis

!     * allocated variables for cce()
      !INTEGER, ALLOCATABLE :: parentsid(:)  ! Pointer to optimisable parameters
      !REAL*8,  ALLOCATABLE :: objfunsub(:)  ! Subset of objective function values for members of population
      !REAL*8,  ALLOCATABLE :: invarsub(:,:) ! Subset of parameter sets for members of population
      !LOGICAL, ALLOCATABLE :: selected(:)

!     * allocated variables for simplex()
!      REAL*8, ALLOCATABLE :: centroid(:)  ! Centroid of parameter sets for simplex procedure
 !     REAL*8, ALLOCATABLE :: newpoint(:)  ! New parameter set resulting from simplex procedure

      INTEGER, parameter  :: nparmax = 16

!     ************************************
!     * namelist parameters for shufflepar
!     ************************************

      INTEGER :: vom_command            ! Indicator of optimisation mode (0 for -optimise, 1 for -continue, 2 for compute, 3 for compute ncp only with pars.txt, 4 for optimise without random_seed)

      INTEGER :: i_ncomp_               ! Initial number of complexes
      INTEGER :: i_ncompmin             ! Minimum number of complexes
      REAL*8  :: i_resolution           ! Convergence criterion (fraction of max variation when optimisation stops)
      INTEGER :: i_patience             ! Number of runs without improvement until optimisation is aborted
      INTEGER :: i_nsimp                ! Number of simplex runs per complex
      REAL*8  :: i_focus                ! Spread of the random seed around the initial values (if <1, then limited)
      INTEGER :: i_iter                 ! Maximum iterations in case of random runs
      INTEGER :: vom_npar               ! Number of model parameters carried through
      LOGICAL :: sce_restart            ! Restarting SCE or starting from scratch
      REAL*8  :: runtime_limit          ! Maximum runtime for sce; 1day (= 1440 minutes)

!     ************************************
!     * namelist parameters for shufflevar
!     ************************************

      CHARACTER(9), ALLOCATABLE :: parname(:) ! Parameter names
      REAL*8,       ALLOCATABLE :: parval(:)  ! Initial parameter values read from shuffle.par
      REAL*8,       ALLOCATABLE :: parmin(:)  ! Minimum parameter values defining search domain
      REAL*8,       ALLOCATABLE :: parmax(:)  ! Maximum parameter values defining search domain
      INTEGER,      ALLOCATABLE :: paropt(:)


      LOGICAL                   :: vd_d_out   !flag for ouput file vd_d
      LOGICAL                   :: esoil_out  !flag for ouput file esoil
      LOGICAL                   :: jmax25t_out!flag for ouput file jmax25t
      LOGICAL                   :: jmax25g_out!flag for ouput file jmax25g
      LOGICAL                   :: vegcov_out !flag for ouput file vegcov
      LOGICAL                   :: resp_out   !flag for ouput file resp_out
      LOGICAL                   :: lambdat_out!flag for ouput file lambdat
      LOGICAL                   :: lambdag_out!flag for ouput file lambdag
      LOGICAL                   :: rrt_out    !flag for ouput file rrt
      LOGICAL                   :: rrg_out    !flag for ouput file rrg
      LOGICAL                   :: asst_out   !flag for ouput file asst
      LOGICAL                   :: assg_out   !flag for ouput file assg
      LOGICAL                   :: su_av_out  !flag for ouput file su_av
      LOGICAL                   :: zw_out     !flag for ouput file zw
      LOGICAL                   :: wsnew_out  !flag for ouput file wsnew
      LOGICAL                   :: spgfcf_out !flag for ouput file spgfcf
      LOGICAL                   :: infx_out   !flag for ouput file infx
      LOGICAL                   :: etmt_out   !flag for ouput file etmt
      LOGICAL                   :: etmg_out   !flag for ouput file etmg
      LOGICAL                   :: su1_out    !flag for ouput file su1
      LOGICAL                   :: topt_out   !flag for ouput file topt
        
      !$OMP threadprivate(ranscal, bestobj, bestincomp, evolution)


      end module vom_sce_mod
