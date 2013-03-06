!     ******************************************************************
!     * File definitions for VOM
!     ******************************************************************

      module vom_file_mod

!     * file codes

      INTEGER :: kfile_dailyweather  = 101
      INTEGER :: kfile_hourlyweather = 102

      INTEGER :: kfile_inputpar      = 201
      INTEGER :: kfile_resultshourly = 201
      INTEGER :: kfile_resultsdaily  = 202
      INTEGER :: kfile_yearly        = 203
      INTEGER :: kfile_rsurfdaily    = 204
      INTEGER :: kfile_delyuhourly   = 205
      INTEGER :: kfile_ruptkhourly   = 206
      INTEGER :: kfile_suvechourly   = 207
      INTEGER :: kfile_soilprofile   = 208
      INTEGER :: kfile_model_output  = 209

      INTEGER :: kfile_sceout        = 701
      INTEGER :: kfile_progress      = 702
      INTEGER :: kfile_lastloop      = 703
      INTEGER :: kfile_currentbest   = 704
      INTEGER :: kfile_bestpars      = 705
      INTEGER :: kfile_finalbest     = 706
      INTEGER :: kfile_pars          = 707
      INTEGER :: kfile_shufflepar    = 708

!     * file names

      CHARACTER(len=*),parameter :: sfile_dailyweather  = 'dailyweather.prn'
      CHARACTER(len=*),parameter :: sfile_hourlyweather = 'hourlyweather.prn'

      CHARACTER(len=*),parameter :: sfile_inputpar      = 'input.par'
      CHARACTER(len=*),parameter :: sfile_resultshourly = 'resultshourly.txt'
      CHARACTER(len=*),parameter :: sfile_resultsdaily  = 'resultsdaily.txt'
      CHARACTER(len=*),parameter :: sfile_yearly        = 'yearly.txt'
      CHARACTER(len=*),parameter :: sfile_rsurfdaily    = 'rsurfdaily.txt'
      CHARACTER(len=*),parameter :: sfile_delyuhourly   = 'delyuhourly.txt'
      CHARACTER(len=*),parameter :: sfile_ruptkhourly   = 'ruptkhourly.txt'
      CHARACTER(len=*),parameter :: sfile_suvechourly   = 'suvechourly.txt'
      CHARACTER(len=*),parameter :: sfile_soilprofile   = 'soilprofile.par'
      CHARACTER(len=*),parameter :: sfile_model_output  = 'model_output.txt'

      CHARACTER(len=*),parameter :: sfile_sceout        = 'sce.out'
      CHARACTER(len=*),parameter :: sfile_progress      = 'progress.txt'
      CHARACTER(len=*),parameter :: sfile_lastloop      = 'lastloop.txt'
      CHARACTER(len=*),parameter :: sfile_currentbest   = 'currentbest.txt'
      CHARACTER(len=*),parameter :: sfile_bestpars      = 'bestpars.txt'
      CHARACTER(len=*),parameter :: sfile_finalbest     = 'finalbest.txt'
      CHARACTER(len=*),parameter :: sfile_pars          = 'pars.txt'
      CHARACTER(len=*),parameter :: sfile_shufflepar    = 'shuffle.par'

      end module vom_file_mod

!***********************************************************************
!*  Module defining variables and parameters for the vegetation model
!*  (transpmodel).
!*----------------------------------------------------------------------

      module vegmod
      implicit none

      INTEGER :: optmode                ! Indicator of optimisation mode
      REAL*8  :: time                   ! Seconds of hour
      REAL*8  :: error                  ! Cumulative error in water balance

      REAL*8, PARAMETER :: p_a     = 1.6d0        ! Ratio of diffusivities of water vapour to CO2 in air
      REAL*8, PARAMETER :: p_pi    = 3.14159d0    ! Pi-constant
      REAL*8, PARAMETER :: p_mpbar = 10.2d0       ! Conversion factor from MPa to bar
      REAL*8, PARAMETER :: p_E     = 2.7182818d0  ! Eurler's number
      REAL*8, PARAMETER :: p_R_    = 8.314d0      ! Molar gas konstant

      INTEGER :: nyear                  ! Year
      INTEGER :: nday                   ! Day of year
      INTEGER :: nhour                  ! Hour of day
      INTEGER :: testday                ! Number of days for initial check if netass>0

      INTEGER, ALLOCATABLE :: fyear(:)   ! Year for each day
      INTEGER, ALLOCATABLE :: fmonth(:)  ! Month for each day
      INTEGER, ALLOCATABLE :: fday(:)    ! Day of month
      INTEGER, ALLOCATABLE :: dayyear(:) ! Day of year

!     * climate

      REAL*8              :: tair__     ! Air temperature (K)
      REAL*8, ALLOCATABLE :: tair_h(:)  ! Hourly air temperature (K)
      REAL*8, ALLOCATABLE :: tairmin(:) ! Daily minimum temperature (K)
      REAL*8, ALLOCATABLE :: tairmax(:) ! Daily maximum temperature (K)

      REAL*8              :: topt_      ! Optimal temperature in temperature response curve

      REAL*8, ALLOCATABLE :: press(:)   ! Daily air pressure (Pa)

      REAL*8              :: par__      ! Photosynthetically active radiation (mol/m2/s)
      REAL*8, ALLOCATABLE :: par_h(:)   ! Hourly photosynthetically active radiation
      REAL*8, ALLOCATABLE :: parvec(:)  ! Daily photosynthetically active radiation
      REAL*8              :: par_y      ! Annual photosynthetically active radiation

      REAL*8, ALLOCATABLE :: srad__(:)  ! Daily shortwave radiation
      REAL*8              :: srad_y     ! Annual shortwave radiation

      REAL*8              :: ca__       ! CO2 mole fraction
      REAL*8, ALLOCATABLE :: ca_h(:)    ! Hourly atmospheric CO2 mole fraction
      REAL*8, ALLOCATABLE :: cavec(:)   ! Daily atmospheric CO2 mole fraction

      REAL*8              :: vp_        ! Atmospheric vapour deficit (VPD/air pressure)
      REAL*8, ALLOCATABLE :: vpvec(:)   ! Daily vapour pressure (Pa)

      REAL*8              :: vd__       ! Atmospheric vapour deficit (VPD/air pressure)
      REAL*8, ALLOCATABLE :: vd_h(:)    ! Hourly vapour deficit (VPD/air pressure)
      REAL*8              :: vd_d       ! Mean daily atmospheric vapour deficit
      REAL*8              :: vd_y       ! Mean annual atmospheric vapour deficit

      REAL*8              :: rain__     ! Rainfall rate (m/s)
      REAL*8, ALLOCATABLE :: rain_h(:)  ! Hourly rainfall
      REAL*8, ALLOCATABLE :: rainvec(:) ! Daily rainfall
      REAL*8              :: rain_y     ! Annual rainfall

!     * soil

      REAL*8, ALLOCATABLE :: phydrostaticvec(:) ! Hydrostatic head in each layer relative to soil surface

      REAL*8  :: gammastar_             ! CO2 compensation point

      REAL*8  :: wsnew                  ! Total soil water store at next time step
      REAL*8  :: wsold                  ! Previous total soil water storage

      REAL*8  :: pc_                    ! Projected cover perennial vegetation (0-1)
      REAL*8  :: pcg_(3)                ! Projected cover seasonal vegetation (pcg_(2) is actual value)
      REAL*8  :: pcgmin                 ! Minimum grass pc; initial point for growth

!     * leaf

      REAL*8  :: wsexp                  ! Exponent for calculating lambda_
      REAL*8  :: wsgexp                 ! Exponent for calculating lambdag
      REAL*8  :: lambdafac              ! Factor for calculating lambda_
      REAL*8  :: lambdagfac             ! Factor for calculating lambdag
      REAL*8  :: lambda_                ! Target dE/dA for calculating gstom__
      REAL*8  :: lambdag_               ! Target dE/dA for calculating gstomg__
      REAL*8  :: gstom__                ! Tree stomatal conductance
      REAL*8  :: gstomg__(3,3)          ! Grass stomatal conductance

      REAL*8  :: rl__(3)                ! Tree leaf respiration for different values of Jmax (rl__(2) is actual value)
      REAL*8  :: rl_d                   ! Daily tree leaf respiration
      REAL*8  :: rlt_y                  ! Annual tree leaf respiration
      REAL*8  :: rlg__(3,3)             ! Grass leaf respiration
      REAL*8  :: rlg_d                  ! Daily grass leaf respiration
      REAL*8  :: rlg_y                  ! Annual grass leaf respiration

      REAL*8  :: transp_                ! Tree transpiration rate
      REAL*8  :: transpg(3,3)           ! Grass transpiration rate (mol/m2/s)

      REAL*8  :: tc_                    ! Tree foliage turnover costs
      REAL*8  :: tct_y                  ! Annual tree foliage turnover costs
      REAL*8  :: tcg(3)                 ! Grass foliage turnover costs
      REAL*8  :: tcg_y                  ! Annual grass foliage turnover costs

      REAL*8  :: jact_(3)               ! Electron transport rates for different values of Jmax (jact_(2) is actual value)
      REAL*8  :: jactg(3,3)             ! Grass electron transport rate

      REAL*8  :: jmax__(3)              ! Tree photosynthetic electron transport capacity
      REAL*8  :: jmaxg__(3)             ! Grass electron transport capacity

      REAL*8  :: jmax25_(3)             ! Tree photosynthetic electron transport capacity at 25oC
      REAL*8  :: jmax25g(3)             ! Grass photosynthetic electron transport capacity at 25oC

!     * plant water

      REAL*8  :: ass_h(3)               ! Tree hourly assimilation rate for different values of Jmax (ass_h(2) is actual value)
      REAL*8  :: ass_d(3)               ! Daily tree assimilation
      REAL*8  :: asst_y                 ! Annual tree assimilation
      REAL*8  :: assg_h(3,3)            ! Hourly grass assimilation
      REAL*8  :: assg_d(3,3)            ! Daily grass assimilation
      REAL*8  :: assg_y                 ! Annual grass assimilation

      REAL*8  :: cpcc_                  ! Tree water transport costs as a function of projected cover and rooting depth (mol/m2/s)
      REAL*8  :: cpcct_y                ! Annual tree water transport costs
      REAL*8  :: cpccg(3)               ! Grass water transport costs
      REAL*8  :: cpccg_y                ! Annual grass water transport costs

      REAL*8  :: etm__                  ! Transpiration rate (m/s)
      REAL*8  :: etm_h                  ! Hourly transpiration
      REAL*8  :: etm_d                  ! Daily transpiration rate
      REAL*8  :: etmt_y                 ! Annual tree transpiration
      REAL*8  :: etmg__(3,3)            ! Grass transpiration rate (m/s)
      REAL*8  :: etmg_h                 ! Hourly grass transpiration
      REAL*8  :: etmg_d                 ! Daily grass transpiration
      REAL*8  :: etmg_y                 ! Annual grass transpiration
      REAL*8  :: etm_y                  ! Annual total transpiration

      REAL*8  :: mq_                    ! Tree water content
      REAL*8  :: mqnew                  ! Tree water content in next time step
      REAL*8  :: mqold                  ! Previous tree water content
      REAL*8  :: dmq                    ! Rate of change in tree water content
      REAL*8  :: mqx_                   ! Tree maximum water content per ground area
      REAL*8  :: mqss_                  ! Tree water content at steady state
      REAL*8  :: mqssmin                ! Tree water content at turgor loss point

      REAL*8  :: md_                    ! Tree dry mass per unit ground area
      REAL*8  :: mdstore                ! Wood water storage parameter of trees

!     * roots

      REAL*8              :: rootdepth  ! Tree rooting depth (m)
      REAL*8              :: rgdepth    ! Grass rooting depth

      INTEGER             :: pos_       ! Lowest soil layer containing tree roots
      INTEGER             :: posg       ! Lowest soil layer containing grass roots
      INTEGER             :: postemp_   ! Lowest soil layer containing tree roots within unsaturated zone
      INTEGER             :: postempg   ! Lowest soil layer containing grass roots within unsaturated zone

      REAL*8              :: changef    ! Change factor for adjusting root surface area

      REAL*8, ALLOCATABLE :: rsurfvec(:)    ! Root surface area of trees in each layer
      REAL*8, ALLOCATABLE :: rsurfnewvec(:) ! Adjusted root surface area of trees in each layer for next day
      REAL*8, ALLOCATABLE :: rsurfg_(:)     ! Root surface area of grasses in each layer
      REAL*8, ALLOCATABLE :: rsurfgnew(:)   ! Adjusted root surface area of grasses in each layer for next day

      REAL*8              :: rootlim(3,3)  ! Indicator whether root surface are was limiting root water uptake

      REAL*8, ALLOCATABLE :: rsoilvec(:)  ! Resistance to water flow towards roots in each soil layer

      REAL*8, ALLOCATABLE :: reff(:)    ! Relative root water uptake efficiency for trees in each layer
      REAL*8, ALLOCATABLE :: reffg(:)   ! Relative root water uptake efficiency for grasses in each layer
      INTEGER             :: pos1(1)    ! Pointer to variable values that achieved maximum assimilation
      INTEGER             :: pos2(2)    ! Pointer to variable values that achieved maximum net assimilation

      REAL*8              :: rr_        ! Tree root respiration rate (mol/m2/s)
      REAL*8              :: rrt_y      ! Annual tree root respiration
      REAL*8              :: rrg        ! Grass root respiration
      REAL*8              :: rrg_y      ! Annual grass root respiration

      REAL*8, ALLOCATABLE :: prootmvec(:)  ! Root hydraulic head in each layer

      REAL*8              :: ruptk_h       ! Hourly tree root water uptake
      REAL*8, ALLOCATABLE :: ruptkvec(:)   ! Root water uptake rate perennial veg (m/s)
      REAL*8, ALLOCATABLE :: ruptkvec_h(:) ! Hourly root water uptake by trees in each layer
      REAL*8, ALLOCATABLE :: ruptkvec_d(:) ! Daily root water uptake by trees in each layer
      REAL*8, ALLOCATABLE :: ruptkg(:)     ! Root water uptake rate seasonal veg (m/s)
      REAL*8, ALLOCATABLE :: hruptkg(:)    ! Hourly root water uptake by grasses in each layer
      REAL*8, ALLOCATABLE :: ruptkg_d(:)   ! Daily root water uptake by grasses in each layer

!     ****************************
!     * input parameters input.par
!     ****************************

      REAL*8  :: alpha                  ! Initial slope of electron transport curve
      REAL*8  :: cpccf                  ! Water transport costs per m root depth and m^2 cover
      REAL*8  :: tcf                    ! Turnover cost factor for foliage (tc=tcf*LAI)
      INTEGER :: maxyear                ! Number of years to process
      INTEGER :: testyear               ! Number of years after which to perform initial test of netass 
      REAL*8  :: ha_                    ! Temperature response parameter
      REAL*8  :: hd_                    ! Temperature response parameter
      REAL*8  :: toptfac                ! Parameter to calculate adaptation of topt (range 0-1 for no to full adaptation)
      REAL*8  :: toptstart              ! Start parameter for topt to calculate jmax(temp in K)
      REAL*8  :: rlratio                ! Ratio of leaf respiration to photosynthetic capacity

!     * Catchment parameters

      REAL*8  :: lat_                   ! geogr. latitude

!     * Soil parameters

!     * Vertical Resolution

!     * Vegetation Parameters

      REAL*8  :: mdf                    ! Total dry mass of living tissues of trees per unit pc (g/m^2)
      REAL*8  :: mqxf                   ! Total water storage capacity in living tissues of trees per unit md
      REAL*8  :: rrootm                 ! Root water uptake resistivity in soil
      REAL*8  :: rsurfmin               ! Minimum root area per m^3 to be maintained
      REAL*8  :: rsurfinit              ! Initial root surface area per m^3
      REAL*8  :: rootrad                ! Average fine root radius
      REAL*8  :: prootmg                ! Constant root balance pressure of 1.5 MPa in grasses
      REAL*8  :: growthmax              ! Parameter determining maximum daily growth increment of root surface area

      INTEGER :: firstyear              ! First year for the generation of hourly output in computation mode
      INTEGER :: lastyear               ! Last year for the generation of hourly output in computation mode

!     * Derived parameters

      REAL*8  :: epsln_                 ! Soil porosity

      INTEGER :: maxhour                ! Number of hours to process
      INTEGER :: maxday                 ! Number of days to process

      end module vegmod

!***********************************************************************
!*  Module defining variables and parameters for the water balance 
!*  model (watbal).
!*----------------------------------------------------------------------
!*  Author: Stan Schymanski, Max Planck Institute for Biogeochemistry
!*  Email: sschym@bgc-jena.mpg.de
!*  06/2008
!*----------------------------------------------------------------------

!     ************************************
!     * Parameters for Water balance model
!     ************************************

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

      REAL*8  :: ys_                    ! Elevation of water table
      REAL*8  :: ysnew                  ! Elevation of water table at next time step

      REAL*8  :: wc_

      REAL*8, ALLOCATABLE :: cH2Ol_s(:) ! Soil water content in each layer (could also be called thetavec(:))

      REAL*8, ALLOCATABLE :: qblvec(:)  ! Water flux across bottom boundary of each layer (positive upwards)

      REAL*8, ALLOCATABLE :: pcapvec(:)    ! Soil matric head (m)
      REAL*8, ALLOCATABLE :: pcapnewvec(:) ! Matric pressure head in each layer at next time step

      REAL*8              :: io__       ! Mass balance of soil water
      REAL*8, ALLOCATABLE :: iovec(:)   ! Mass balance for each soil layer
      REAL*8              :: io_h       ! Hourly mass balance of soil water
      REAL*8              :: ioacum     ! Cumulative mass balance of soil water

      REAL*8, ALLOCATABLE :: suvec_(:)   ! Soil saturation degree in each layer
      REAL*8, ALLOCATABLE :: sunewvec(:) ! Soil saturation degree in each layer at next time step
      REAL*8, ALLOCATABLE :: sueqvec(:)  ! Soil saturation degree above water table in hydrostatic equilibrium
      REAL*8, ALLOCATABLE :: dsuvec(:)   ! Rate of change in soil saturation degree in each layer

      REAL*8, ALLOCATABLE :: kunsatvec(:)    ! Unsaturated hydraulic conductivity in each soil layer
      REAL*8, ALLOCATABLE :: kunsatnewvec(:) ! Unsaturated hydraulic conductivity in each layer at next time step

      REAL*8, ALLOCATABLE :: ksatvec(:)    ! Saturated hydraulic conductivity
      REAL*8, ALLOCATABLE :: thetasvec(:)  ! Saturated soil water content
      REAL*8, ALLOCATABLE :: thetarvec(:)  ! Residual soil water content
      REAL*8, ALLOCATABLE :: nvgvec(:)     ! Van Genuchten soil parameter n
      REAL*8, ALLOCATABLE :: avgvec(:)     ! Van Genuchten soil parameter a
      REAL*8, ALLOCATABLE :: mvgvec(:)     ! Van Genuchten soil parameter m

      REAL*8  :: daylength              ! Day length (hours)
      REAL*8  :: netassg_d(3,3)         ! Daily grass net carbon profit
      REAL*8  :: dtmq                   ! Maximum timestep allowed by tree water content change

!     ****************************
!     * input parameters input.par
!     ****************************

!     * Catchment parameters

      REAL*8  :: cgs                    ! Capital Gamma S (length scale for seepage outflow) (m)
      REAL*8  :: zr_                    ! Average channel elevation
      REAL*8  :: go_                    ! Slope close to channel in radians 

!     * Soil parameters

      REAL*8  :: ksat_                  ! Saturated hydraulic conductivity
      REAL*8  :: thetar_                ! Residual soil water content
      REAL*8  :: thetas_                ! Saturated soil water content
      REAL*8  :: nvg_                   ! Van Genuchten soil parameter n
      REAL*8  :: avg_                   ! Van Genuchten soil parameter a
      REAL*8  :: mvg_                   ! Van Genuchten soil parameter m

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

!     * unused
      REAL*8  :: oa_                    ! O2 mole fraction
      REAL*8  :: k25                    ! Temperature response parameter
      REAL*8  :: kopt                   ! Temperature response parameter
      REAL*8  :: zs_                    ! Reference elevation point (set to 0)
      REAL*8  :: inf_h                  ! Hourly infiltration rate
      REAL*8  :: ruptk_d                ! Daily tree root water uptake
      REAL*8  :: epan_y                 ! Annual pan evaporation
      REAL*8  :: mqssnew                ! Tree water content in next time step
      REAL*8  :: jmax_d                 ! Mean daily tree electron transport capacity
      REAL*8  :: jmaxg_d                ! Daily grass electron transport capacity
      REAL*8  :: gstom_d                ! Mean daily tree stomatal conductance
      REAL*8  :: gstomg_d               ! Daily average stomatal conductance
      REAL*8, PARAMETER :: g___   = 9.81d0       ! Gravitational acceleration
      REAL*8, PARAMETER :: rho    = 1000.d0      ! Density of water
      REAL*8, PARAMETER :: degree = 0.0174533d0  ! Conversion from degree to radians
      REAL*8, ALLOCATABLE :: epslnvec(:)     ! Soil porosity
      REAL*8, ALLOCATABLE :: gammastarvec(:) ! Hourly value for CO2 compensation point
      REAL*8, ALLOCATABLE :: netassvec_(:)   ! Daily net assimilation by trees
      REAL*8, ALLOCATABLE :: netassvecg(:)   ! Daily net assimilation by grasses
      REAL*8, ALLOCATABLE :: epan__(:)       ! Daily pan evaporation
!     REAL*8, ALLOCATABLE :: wsnewvec(:)     ! Soil water content in each layer at next time step
!     REAL*8, ALLOCATABLE :: delznewvec(:)
!     REAL*8  :: yu__
!     REAL*8  :: yunew_
!     REAL*8  :: omgo_
!     REAL*8  :: omgonew
!     REAL*8  :: omgu_
!     REAL*8  :: omgunew
!d    REAL*8  :: netass_y
!d    REAL*8  :: gpp_y

      INTEGER :: maxlayer               ! Number of soil layers

      REAL*8  :: cz                     ! Average soil elevation in m
      INTEGER :: parsaved = 0           ! Indicator whether input parameter values have been read and saved

      REAL*8              :: delz_      ! Thickness of each soil layer (m)
      REAL*8, ALLOCATABLE :: delzvec(:) ! Thickness of each soil layer

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

      REAL*8, ALLOCATABLE :: shufflevar(:,:)  ! Population of parameter sets
      REAL*8, ALLOCATABLE :: ofvec(:)   ! Population of objective function values related to shufflevar

      CHARACTER(12)  :: evolution       ! Character string explaining how parameter set was generated
      CHARACTER(60)  :: outformat, loopformat

!     * allocated variables for optsensitivity()
      REAL*8, ALLOCATABLE :: dataarray(:,:)  ! Parameter sets and objective functions used for sensitivity analysis
      REAL*8, ALLOCATABLE :: shufflevar2(:)  ! Parameter sets used for sensitivity analysis

!     * allocated variables for cce()
      INTEGER, ALLOCATABLE :: parentsid(:)  ! Pointer to optimisable parameters
      REAL*8,  ALLOCATABLE :: objfunsub(:)  ! Subset of objective function values for members of population
      REAL*8,  ALLOCATABLE :: invarsub(:,:) ! Subset of parameter sets for members of population
      LOGICAL, ALLOCATABLE :: selected(:)

!     * allocated variables for simplex()
      REAL*8, ALLOCATABLE :: centroid(:)  ! Centroid of parameter sets for simplex procedure
      REAL*8, ALLOCATABLE :: newpoint(:)  ! New parameter set resulting from simplex procedure

!     ************************************
!     * input parameters for shuffle.par
!     ************************************

      INTEGER :: command                ! Indicator of optimisation mode (0 for -optimise, 1 for -continue, 2 for compute, 3 for compute ncp only with pars.txt, 4 for optimise without random_seed)

      INTEGER :: ncomp                  ! Initial number of complexes
      INTEGER :: ncompmin               ! Minimum number of complexes
      REAL*8  :: resolution             ! Convergence criterion (fraction of max variation when optimisation stops)
      INTEGER :: patience               ! Number of runs without improvement until optimisation is aborted
      INTEGER :: nsimp                  ! Number of simplex runs per complex
      REAL*8  :: focus                  ! Spread of the random seed around the initial values (if <1, then limited)

      INTEGER :: npar                   ! Number of model parameters carried through
      CHARACTER(9), ALLOCATABLE :: parname(:) ! Parameter names
      REAL*8,       ALLOCATABLE :: parval(:)  ! Initial parameter values read from shuffle.par
      REAL*8,       ALLOCATABLE :: parmin(:)  ! Minimum parameter values defining search domain
      REAL*8,       ALLOCATABLE :: parmax(:)  ! Maximum parameter values defining search domain
      INTEGER,      ALLOCATABLE :: paropt(:)

      end module vom_sce_mod
