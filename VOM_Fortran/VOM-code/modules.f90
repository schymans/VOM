!     ******************************************************************
!     * File definitions for VOM
!     ******************************************************************

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

      CHARACTER(len=*),parameter :: sfile_namelist      = 'vom_namelist'
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

      CHARACTER(len=*),parameter :: sfile_sceout        = 'sce_out.txt'
      CHARACTER(len=*),parameter :: sfile_progress      = 'sce_progress.txt'
      CHARACTER(len=*),parameter :: sfile_lastloop      = 'sce_lastloop.txt'
      CHARACTER(len=*),parameter :: sfile_lastbest      = 'sce_lastbest.txt'
      CHARACTER(len=*),parameter :: sfile_bestpars      = 'sce_bestpars.txt'
      CHARACTER(len=*),parameter :: sfile_beststat      = 'sce_status.txt'
      CHARACTER(len=*),parameter :: sfile_pars          = 'pars.txt'

      CHARACTER*100  :: i_outputpath     ! Constant root balance pressure of 1.5 MPa in grasses
      CHARACTER*100  :: i_inputpath      ! Constant root balance pressure of 1.5 MPa in grasses

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

      REAL*8, ALLOCATABLE :: tair_h(:)    ! Hourly air temperature (K)
      REAL*8, ALLOCATABLE :: tairmin_d(:) ! Daily minimum temperature (K)
      REAL*8, ALLOCATABLE :: tairmax_d(:) ! Daily maximum temperature (K)

      REAL*8              :: topt_      ! Optimal temperature in temperature response curve

      REAL*8, ALLOCATABLE :: press_d(:) ! Daily air pressure (Pa)

      REAL*8, ALLOCATABLE :: par_h(:)   ! Hourly photosynthetically active radiation (mol/m2/s)
      REAL*8, ALLOCATABLE :: par_d(:)   ! Daily photosynthetically active radiation
      REAL*8              :: par_y      ! Annual photosynthetically active radiation

      REAL*8, ALLOCATABLE :: srad_d(:)  ! Daily shortwave radiation
      REAL*8              :: srad_y     ! Annual shortwave radiation

      REAL*8, ALLOCATABLE :: ca_h(:)    ! Hourly atmospheric CO2 mole fraction
      REAL*8, ALLOCATABLE :: ca_d(:)    ! Daily atmospheric CO2 mole fraction

      REAL*8, ALLOCATABLE :: vp_d(:)    ! Daily absolute vapour pressure (Pa)

      REAL*8, ALLOCATABLE :: vd_h(:)    ! Hourly atmospheric vapour deficit (VPD/air pressure)
      REAL*8              :: vd_d       ! Mean daily atmospheric vapour deficit
      REAL*8              :: vd_y       ! Mean annual atmospheric vapour deficit

      REAL*8, ALLOCATABLE :: rain_h(:)  ! Hourly rainfall rate (m/s)
      REAL*8, ALLOCATABLE :: rain_d(:)  ! Daily rainfall
      REAL*8              :: rain_y     ! Annual rainfall

!     * soil

      REAL*8, ALLOCATABLE :: c_hhydrst(:)  ! Hydrostatic head in each layer relative to soil surface

      REAL*8  :: gammastar              ! CO2 compensation point

      REAL*8  :: wsnew                  ! Total soil water store at next time step
      REAL*8  :: wsold                  ! Previous total soil water storage

      REAL*8  :: o_pct = 0.300000d+00   ! Projected cover perennial vegetation (0-1)
      REAL*8  :: pcg_d(3)               ! Projected cover seasonal vegetation (pcg_d(2) is actual value)
      REAL*8  :: c_pcgmin               ! Minimum grass pc; initial point for growth

!     * leaf

      REAL*8  :: o_wstexp   = -0.564496d+00  ! Exponent for calculating lambdat_d
      REAL*8  :: o_wsgexp   = -0.132889d+01  ! Exponent for calculating lambdag
      REAL*8  :: o_lambdatf =  0.160181d+04  ! Factor for calculating lambdat_d
      REAL*8  :: o_lambdagf =  0.779827d+03  ! Factor for calculating lambdag_d
      REAL*8  :: lambdat_d              ! Target dE/dA for calculating gstomt
      REAL*8  :: lambdag_d              ! Target dE/dA for calculating gstomg
      REAL*8  :: gstomt                 ! Tree stomatal conductance
      REAL*8  :: gstomg(3,3)            ! Grass stomatal conductance

      REAL*8  :: rlt_h(3)               ! Tree leaf respiration for different values of Jmax (rlt_h(2) is actual value)
      REAL*8  :: rlt_d                  ! Daily tree leaf respiration
      REAL*8  :: rlt_y                  ! Annual tree leaf respiration
      REAL*8  :: rlg_h(3,3)             ! Grass leaf respiration
      REAL*8  :: rlg_d                  ! Daily grass leaf respiration
      REAL*8  :: rlg_y                  ! Annual grass leaf respiration

      REAL*8  :: transpt                ! Tree transpiration rate
      REAL*8  :: transpg(3,3)           ! Grass transpiration rate (mol/m2/s)

      REAL*8  :: q_tct_d                ! Tree foliage turnover costs
      REAL*8  :: tct_y                  ! Annual tree foliage turnover costs
      REAL*8  :: tcg_d(3)               ! Grass foliage turnover costs
      REAL*8  :: tcg_y                  ! Annual grass foliage turnover costs

      REAL*8  :: jactt(3)               ! Electron transport rates for different values of Jmax (jactt(2) is actual value)
      REAL*8  :: jactg(3,3)             ! Grass electron transport rate

      REAL*8  :: jmaxt_h(3)             ! Tree photosynthetic electron transport capacity
      REAL*8  :: jmaxg_h(3)             ! Grass electron transport capacity

      REAL*8  :: jmax25t_d(3)           ! Tree photosynthetic electron transport capacity at 25oC
      REAL*8  :: jmax25g_d(3)           ! Grass photosynthetic electron transport capacity at 25oC

!     * plant water

      REAL*8  :: asst_h(3)              ! Tree hourly assimilation rate for different values of Jmax (asst_h(2) is actual value)
      REAL*8  :: asst_d(3)              ! Daily tree assimilation
      REAL*8  :: asst_y                 ! Annual tree assimilation
      REAL*8  :: assg_h(3,3)            ! Hourly grass assimilation
      REAL*8  :: assg_d(3,3)            ! Daily grass assimilation
      REAL*8  :: assg_y                 ! Annual grass assimilation

      REAL*8  :: q_cpcct_d              ! Tree water transport costs as a function of projected cover and rooting depth (mol/m2/s)
      REAL*8  :: cpcct_y                ! Annual tree water transport costs
      REAL*8  :: cpccg_d(3)             ! Grass water transport costs
      REAL*8  :: cpccg_y                ! Annual grass water transport costs

      REAL*8  :: etmt__                 ! Transpiration rate (m/s)
      REAL*8  :: etmt_h                 ! Hourly transpiration
      REAL*8  :: etmt_d                 ! Daily transpiration rate
      REAL*8  :: etmt_y                 ! Annual tree transpiration
      REAL*8  :: etmg__(3,3)            ! Grass transpiration rate (m/s)
      REAL*8  :: etmg_h                 ! Hourly grass transpiration
      REAL*8  :: etmg_d                 ! Daily grass transpiration
      REAL*8  :: etmg_y                 ! Annual grass transpiration
      REAL*8  :: etm_y                  ! Annual total transpiration

      REAL*8  :: mqt_                   ! Tree water content
      REAL*8  :: mqtnew                 ! Tree water content in next time step
      REAL*8  :: mqtold                 ! Previous tree water content
      REAL*8  :: dmqt                   ! Rate of change in tree water content
      REAL*8  :: q_mqx                  ! Tree maximum water content per ground area
      REAL*8  :: mqsst_                 ! Tree water content at steady state
      REAL*8  :: mqsstmin               ! Tree water content at turgor loss point

      REAL*8  :: q_md                     ! Tree dry mass per unit ground area
      REAL*8  :: o_mdstore = 1.000000d+02 ! Wood water storage parameter of trees

!     * roots

      REAL*8  :: o_rtdepth = 0.300000d+01 ! Tree rooting depth (m)
      REAL*8  :: o_rgdepth = 0.100000d+01 ! Grass rooting depth

      INTEGER             :: pos_slt    ! Lowest soil layer containing tree roots
      INTEGER             :: pos_slg    ! Lowest soil layer containing grass roots
      INTEGER             :: pos_ult    ! Lowest soil layer containing tree roots within unsaturated zone
      INTEGER             :: pos_ulg    ! Lowest soil layer containing grass roots within unsaturated zone

      REAL*8              :: changef    ! Change factor for adjusting root surface area

      REAL*8, ALLOCATABLE :: rsurft_(:)    ! Root surface area of trees in each layer
      REAL*8, ALLOCATABLE :: rsurftnew(:)  ! Adjusted root surface area of trees in each layer for next day
      REAL*8, ALLOCATABLE :: rsurfg_(:)    ! Root surface area of grasses in each layer
      REAL*8, ALLOCATABLE :: rsurfgnew(:)  ! Adjusted root surface area of grasses in each layer for next day

      REAL*8              :: rootlim(3,3)  ! Indicator whether root surface are was limiting root water uptake

      REAL*8, ALLOCATABLE :: rsoil(:)   ! Resistance to water flow towards roots in each soil layer

      REAL*8, ALLOCATABLE :: refft(:)   ! Relative root water uptake efficiency for trees in each layer
      REAL*8, ALLOCATABLE :: reffg(:)   ! Relative root water uptake efficiency for grasses in each layer
      INTEGER             :: posmna(2)  ! Pointer to variable values that achieved maximum net assimilation

      REAL*8              :: rrt_d      ! Tree root respiration rate (mol/m2/s)
      REAL*8              :: rrt_y      ! Annual tree root respiration
      REAL*8              :: rrg_d      ! Grass root respiration
      REAL*8              :: rrg_y      ! Annual grass root respiration

      REAL*8, ALLOCATABLE :: prootm(:)  ! Root hydraulic head in each layer

      REAL*8              :: sumruptkt_h  ! Hourly total tree root water uptake
      REAL*8, ALLOCATABLE :: ruptkt__(:)  ! Root water uptake rate perennial veg (m/s)
      REAL*8, ALLOCATABLE :: ruptkt_h(:)  ! Hourly root water uptake by trees in each layer
      REAL*8, ALLOCATABLE :: ruptkt_d(:)  ! Daily root water uptake by trees in each layer
      REAL*8, ALLOCATABLE :: ruptkg__(:)  ! Root water uptake rate seasonal veg (m/s)
      REAL*8, ALLOCATABLE :: ruptkg_h(:)  ! Hourly root water uptake by grasses in each layer
      REAL*8, ALLOCATABLE :: ruptkg_d(:)  ! Daily root water uptake by grasses in each layer


!     ****************************
!     * input parameters input.par
!     ****************************

      REAL*8  :: i_alpha   = 0.3d0      ! Initial slope of electron transport curve
      REAL*8  :: i_cpccf   = 1.2d-6     ! Water transport costs per m root depth and m^2 cover
      REAL*8  :: i_tcf     = 2.2d-7     ! Turnover cost factor for foliage (tc=i_tcf*LAI)
      INTEGER :: i_maxyear = 30         ! Number of years to process
      INTEGER :: i_testyear = 1         ! Number of years after which to perform initial test of netass
      REAL*8  :: i_ha      = 43790.0d0  ! Temperature response parameter
      REAL*8  :: i_hd      = 2.0d5      ! Temperature response parameter
      REAL*8  :: i_toptf   = 0.0d0      ! Parameter to calculate adaptation of topt (range 0-1 for no to full adaptation)
      REAL*8  :: i_toptstart = 305.0d0  ! Start parameter for topt to calculate jmax(temp in K)
      REAL*8  :: i_rlratio = 0.07d0     ! Ratio of leaf respiration to photosynthetic capacity

!     * Catchment parameters

      REAL*8  :: i_lat     = 12.5d0     ! geogr. latitude

!     * Soil parameters

!     * Vertical Resolution

!     * Vegetation Parameters

      REAL*8  :: i_mdtf    = 10000.0d0  ! Total dry mass of living tissues of trees per unit pc (g/m^2)
      REAL*8  :: i_mqxtf   = 1.0d0      ! Total water storage capacity in living tissues of trees per unit md
      REAL*8  :: i_rrootm  = 1.02d8     ! Root water uptake resistivity in soil
      REAL*8  :: i_rsurfmin = 0.03d0    ! Minimum root area per m^3 to be maintained
      REAL*8  :: i_rsurf_  = 0.3d0      ! Initial root surface area per m^3
      REAL*8  :: i_rootrad = 0.3d-3     ! Average fine root radius
      REAL*8  :: i_prootmg = 150.0d0    ! Constant root balance pressure of 1.5 MPa in grasses
      REAL*8  :: i_growthmax = 0.1d0    ! Parameter determining maximum daily growth increment of root surface area
      REAL*8  :: i_incrcovg = 0.02d0    ! parameter determining maximum increment percentage of grass cover
      REAL*8  :: i_incrjmax = 0.01d0    ! parameter determining maximum increment percentage of jmax25

      INTEGER :: i_firstyear = 2000     ! First year for the generation of hourly output in computation mode
      INTEGER :: i_lastyear = 2000      ! Last year for the generation of hourly output in computation mode

      INTEGER :: i_write_h = 0          ! Flag to write out hourly input values after conversation from daily values




!     * Derived parameters

      REAL*8  :: c_epsln                ! Soil porosity

      INTEGER :: c_maxhour              ! Number of hours to process
      INTEGER :: c_maxday               ! Number of days to process

      !$OMP threadprivate( time, error, finish, nyear, nday, nhour, th_, c_testday,   & 
      !$OMP topt_, par_y, srad_y,   &
      !$OMP vd_d, vd_y, rain_y, gammastar, wsnew, wsold, o_pct, pcg_d, c_pcgmin, &
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

      REAL*8  :: esoil__                ! Bare soil evaporation rate
      REAL*8  :: esoil_h                ! Hourly soil evaporation
      REAL*8  :: esoil_d                ! Daily soil evaporation rate
      REAL*8  :: esoil_y                ! Annual soil evaporation rate

      REAL*8  :: spgfcf__               ! Seepage face flow rate
      REAL*8  :: spgfcf_h               ! Hourly seepage face flow
      REAL*8  :: spgfcf_d               ! Daily seepage face flow

      REAL*8  :: inf__                  ! Infiltration rate
      REAL*8  :: infx__                 ! Infiltration excess runoff rate
      REAL*8  :: infx_h                 ! Hourly infiltration excess runoff
      REAL*8  :: infx_d                 ! Daily infiltration excess runoff

      REAL*8  :: zw_                    ! Elevation of water table
      REAL*8  :: zwnew                  ! Elevation of water table at next time step

      REAL*8  :: wc_

      REAL*8, ALLOCATABLE :: cH2Ol_s(:) ! Soil water content in each layer (could also be called theta(:))

      REAL*8, ALLOCATABLE :: qbl(:)     ! Water flux across bottom boundary of each layer (positive upwards)

      REAL*8, ALLOCATABLE :: pcap_(:)   ! Soil matric head (m)
      REAL*8, ALLOCATABLE :: pcapnew(:) ! Matric pressure head in each layer at next time step

      REAL*8, ALLOCATABLE :: iovec(:)   ! Mass balance for each soil layer
      REAL*8              :: io__       ! Input minus output of water in soil domain
      REAL*8              :: io_h       ! Hourly input minus output of water in soil domain
      REAL*8              :: ioacum     ! Cumulative mass balance of soil water

      REAL*8, ALLOCATABLE :: su__(:)    ! Soil saturation degree in each layer
      REAL*8, ALLOCATABLE :: sunew(:)   ! Soil saturation degree in each layer at next time step
      REAL*8, ALLOCATABLE :: sueq(:)    ! Soil saturation degree above water table in hydrostatic equilibrium
      REAL*8, ALLOCATABLE :: dsu(:)     ! Rate of change in soil saturation degree in each layer

      REAL*8, ALLOCATABLE :: kunsat_(:)    ! Unsaturated hydraulic conductivity in each soil layer
      REAL*8, ALLOCATABLE :: kunsatnew(:)  ! Unsaturated hydraulic conductivity in each layer at next time step

      REAL*8, ALLOCATABLE :: s_ksat(:)    ! Saturated hydraulic conductivity
      REAL*8, ALLOCATABLE :: s_thetas(:)  ! Saturated soil water content
      REAL*8, ALLOCATABLE :: s_thetar(:)  ! Residual soil water content
      REAL*8, ALLOCATABLE :: s_nvg(:)     ! Van Genuchten soil parameter n
      REAL*8, ALLOCATABLE :: s_avg(:)     ! Van Genuchten soil parameter a
      REAL*8, ALLOCATABLE :: c_mvg(:)     ! Van Genuchten soil parameter m

!     ****************************
!     * input parameters input.par
!     ****************************

!     * Catchment parameters

      REAL*8  :: i_cgs     = 10.0d0     ! Capital Gamma S (length scale for seepage outflow) (m)
      REAL*8  :: i_zr      = 10.0d0     ! Average channel elevation
      REAL*8  :: i_go      = 0.033d0    ! Slope close to channel in radians

!     * Soil parameters

      REAL*8  :: i_ksat    = 1.23d-5    ! Saturated hydraulic conductivity
      REAL*8  :: i_thetar  = 0.065d0    ! Residual soil water content
      REAL*8  :: i_thetas  = 0.41d0     ! Saturated soil water content
      REAL*8  :: i_nvg     = 1.89d0     ! Van Genuchten soil parameter n
      REAL*8  :: i_avg     = 7.5d0      ! Van Genuchten soil parameter a
      REAL*8  :: i_mvg                  ! Van Genuchten soil parameter m

      !$OMP threadprivate(wlayer_, wlayernew, dt_, dtmax, dtsu_count, dtmax_count, esoil__, esoil_h, &
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

      REAL*8  :: i_cz      = 15.0d0     ! Average soil elevation in m

      REAL*8  :: i_delz    = 0.5d0      ! Thickness of each soil layer (m)
      REAL*8, ALLOCATABLE :: s_delz(:)  ! Thickness of each soil layer


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

      INTEGER :: n_thread = 1           ! Number of threads to be used

      INTEGER :: success = 0            ! Indicator wheter optimisation ended successfully

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

      INTEGER, parameter  :: nparmax = 9

!     ************************************
!     * namelist parameters for shufflepar
!     ************************************

      INTEGER :: vom_command            ! Indicator of optimisation mode (0 for -optimise, 1 for -continue, 2 for compute, 3 for compute ncp only with pars.txt, 4 for optimise without random_seed)

      INTEGER :: i_ncomp_     = 2       ! Initial number of complexes
      INTEGER :: i_ncompmin   = 2       ! Minimum number of complexes
      REAL*8  :: i_resolution = 1.0     ! Convergence criterion (fraction of max variation when optimisation stops)
      INTEGER :: i_patience   = 10      ! Number of runs without improvement until optimisation is aborted
      INTEGER :: i_nsimp      = 3       ! Number of simplex runs per complex
      REAL*8  :: i_focus      = 1.0     ! Spread of the random seed around the initial values (if <1, then limited)
      INTEGER :: i_iter       = 10      ! Maximum iterations in case of random runs
      INTEGER :: vom_npar     = 8       ! Number of model parameters carried through

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
