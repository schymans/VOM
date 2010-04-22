!***********************************************************************
!*  Module defining variables and parameters for the vegetation model
!*  (transpmodel) and water balance model (watbal).
!*----------------------------------------------------------------------
!*  Author: Stan Schymanski, Max Planck Institute for Biogeochemistry
!*  Email: sschym@bgc-jena.mpg.de
!*  06/2008
!*----------------------------------------------------------------------

      module vegwatbal
      implicit none

      INTEGER :: optmode

!*=====<Parameters>=====================================================
!* set O2 conc [mol/mol], CO2 conc and diffusion conversion factor from
!* CO2 to H2O as parameters. 


!TODO: replace the following variable names
! M___ by layer
!

      REAL*8  :: oa
      REAL*8  :: ca_
      REAL*8  :: alpha
      REAL*8  :: rlratio
      REAL*8  :: k25
      REAL*8  :: kopt
      REAL*8  :: ha_
      REAL*8  :: hd
      REAL*8  :: tstart
      REAL*8  :: toptfac
      REAL*8  :: cpccf                       ! water transport costs per m root depth and m^2 cover
      REAL*8  :: tcf

      INTEGER :: parsaved = 0
      INTEGER :: ny_
      INTEGER :: nyt
      INTEGER :: firstyear
      INTEGER :: lastyear
      INTEGER :: dtsu_count, dtmax_count   ! for speed analysis, count of how often dtsu or dtmax is time limiting

      REAL*8, PARAMETER :: a__    = 1.6d0
      REAL*8, PARAMETER :: pi     = 3.14159d0
      REAL*8, PARAMETER :: mpbar  = 10.2d0   ! conversion factor from MPa to bar
      REAL*8, PARAMETER :: g___   = 9.81d0
      REAL*8, PARAMETER :: rho    = 1000.d0
      REAL*8, PARAMETER :: degree = 0.0174533d0
      REAL*8, PARAMETER :: E__    = 2.7182818d0
      REAL*8, PARAMETER :: R__    = 8.314d0

      INTEGER :: N__
      INTEGER :: Nh_
      INTEGER :: M___
      INTEGER :: d___
      INTEGER :: h__
      INTEGER :: yr_
      INTEGER :: dtest
      INTEGER :: pos_
      INTEGER :: posg
      INTEGER :: postemp_
      INTEGER :: postempg
      INTEGER :: pos1(1)
      INTEGER :: pos2(2)

      REAL*8, ALLOCATABLE  :: srad(:)
      REAL*8, ALLOCATABLE  :: tmin(:)
      REAL*8, ALLOCATABLE  :: tmax(:)
      REAL*8, ALLOCATABLE  :: rainvec(:)
      REAL*8, ALLOCATABLE  :: netassvec(:)
      REAL*8, ALLOCATABLE  :: epan(:)
      REAL*8, ALLOCATABLE  :: vpvec(:)
      REAL*8, ALLOCATABLE  :: netassvecg(:)
      REAL*8, ALLOCATABLE  :: press(:)
      REAL*8, ALLOCATABLE  :: cavec(:)

      INTEGER, ALLOCATABLE :: year(:)
      INTEGER, ALLOCATABLE :: month(:)
      INTEGER, ALLOCATABLE :: day(:)
      INTEGER, ALLOCATABLE :: dayyear(:)

      REAL*8, ALLOCATABLE  :: parvec(:)
      REAL*8, ALLOCATABLE  :: parh(:)
      REAL*8, ALLOCATABLE  :: vdh(:)
      REAL*8, ALLOCATABLE  :: tairh(:)
      REAL*8, ALLOCATABLE  :: gammastarvec(:)
      REAL*8, ALLOCATABLE  :: rainh(:)
      REAL*8, ALLOCATABLE  :: cah(:)

      REAL*8  :: jact_(3)
      REAL*8  :: lambda_
      REAL*8  :: gstom__
      REAL*8  :: rl__(3)
      REAL*8  :: vd__
      REAL*8  :: par_
      REAL*8  :: rain_
      REAL*8  :: transp_
      REAL*8  :: hass_(3)
      REAL*8  :: pc_
      REAL*8  :: cpcc_
      REAL*8  :: rr_
      REAL*8  :: lambdafac
      REAL*8  :: gammastar
      REAL*8  :: jmax__(3)
      REAL*8  :: jmax25_(3)
!d      REAL*8  :: netassyr
      REAL*8  :: rainyr
      REAL*8  :: epanyr
      REAL*8  :: paryr
      REAL*8  :: radyr
      REAL*8  :: vdyr
      REAL*8  :: etyr
      REAL*8  :: evapyr
!d      REAL*8  :: gppyr
      REAL*8  :: assg_y
      REAL*8  :: rlg_y
      REAL*8  :: rrg_y
      REAL*8  :: cpccg_y
      REAL*8  :: tcg_y
      REAL*8  :: asst_y
      REAL*8  :: rlt_y
      REAL*8  :: rrt_y
      REAL*8  :: cpcct_y
      REAL*8  :: tct_y
      REAL*8  :: daylength
      REAL*8  :: tair
      REAL*8  :: vp_
      REAL*8  :: topt_
      REAL*8  :: rootdepth
      REAL*8  :: mq_
      REAL*8  :: mqssmin
      REAL*8  :: dmq
      REAL*8  :: mqnew
      REAL*8  :: mqss_
      REAL*8  :: mqold
      REAL*8  :: hruptk_
      REAL*8  :: tc_
      REAL*8  :: mqssnew
      REAL*8  :: pcg_(3)
      REAL*8  :: lambdagfac
      REAL*8  :: wsgexp
      REAL*8  :: cpccg(3)
      REAL*8  :: tcg(3)
      REAL*8  :: jmax25g(3)
      REAL*8  :: wsexp
      REAL*8  :: rgdepth
      REAL*8  :: lambdag
      REAL*8  :: rrg
      REAL*8  :: jmaxg__(3)
      REAL*8  :: rlg__(3,3)
      REAL*8  :: jactg(3,3)
      REAL*8  :: gstomg__(3,3)
      REAL*8  :: etmg__(3,3)
      REAL*8  :: transpg(3,3)
      REAL*8  :: hassg(3,3)
      REAL*8  :: hetmg
      REAL*8  :: assg_d(3,3)
      REAL*8  :: jmaxg_d
      REAL*8  :: gstomg_d
      REAL*8  :: etmg_d
      REAL*8  :: rlg_d
      REAL*8  :: netassg_d(3,3)
      REAL*8  :: rootlim(3,3)
      REAL*8  :: changef

      REAL*8  :: md_
      REAL*8  :: mqx_
      REAL*8  :: rsurfmin
      REAL*8  :: rsurfinit
      REAL*8  :: rootrad
      REAL*8  :: rrootm
      REAL*8  :: pcgmin
      REAL*8  :: prootmg
      REAL*8  :: mdf
      REAL*8  :: mqxf
      REAL*8  :: mdstore
      REAL*8  :: growthmax

      REAL*8  :: ruptk_d
      REAL*8  :: jmax_d
      REAL*8  :: rl_d
      REAL*8  :: gstom_d
      REAL*8  :: ass_d(3)
      REAL*8  :: spgfcf_d
      REAL*8  :: infx_d
      REAL*8  :: etm_d
      REAL*8  :: esoil_d
      REAL*8  :: vd_d
      REAL*8  :: time
      REAL*8  :: dtmq
      REAL*8  :: dtmax

!***********************************************************************
!*    Parameters for Water balance model
!***********************************************************************

      INTEGER :: nlayers_
      INTEGER :: nlayersnew

      REAL*8  :: dt_
      REAL*8  :: cz
      REAL*8  :: cgs
      REAL*8  :: lat                         ! geogr. latitude

      REAL*8  :: ksat_
      REAL*8  :: epsln_
      REAL*8  :: nvg_
      REAL*8  :: mvg_
      REAL*8  :: avg_
      REAL*8  :: zs_
      REAL*8  :: zr_
      REAL*8  :: go_
      REAL*8  :: thetas_
      REAL*8  :: thetar_
      REAL*8  :: wsold
      REAL*8  :: delz_

      REAL*8, ALLOCATABLE :: ksatvec(:)
      REAL*8, ALLOCATABLE :: epslnvec(:)
      REAL*8, ALLOCATABLE :: nvgvec(:)
      REAL*8, ALLOCATABLE :: mvgvec(:)
      REAL*8, ALLOCATABLE :: avgvec(:)
      REAL*8, ALLOCATABLE :: thetasvec(:)
      REAL*8, ALLOCATABLE :: thetarvec(:)
      REAL*8, ALLOCATABLE :: sueqvec(:)
      REAL*8, ALLOCATABLE :: cH2Ol_s(:)
      REAL*8, ALLOCATABLE :: iovec(:)

      REAL*8, ALLOCATABLE :: pcapvec(:)
      REAL*8, ALLOCATABLE :: suvec_(:)
      REAL*8, ALLOCATABLE :: ruptkvec(:)
      REAL*8, ALLOCATABLE :: sunewvec(:)
      REAL*8, ALLOCATABLE :: kunsatvec(:)
      REAL*8, ALLOCATABLE :: delzvec(:)
      REAL*8, ALLOCATABLE :: rsurfvec(:)
      REAL*8, ALLOCATABLE :: rsurfnewvec(:)
      REAL*8, ALLOCATABLE :: qblvec(:)
      REAL*8, ALLOCATABLE :: dsuvec(:)
      REAL*8, ALLOCATABLE :: phydrostaticvec(:)
      REAL*8, ALLOCATABLE :: delznewvec(:)
      REAL*8, ALLOCATABLE :: wsnewvec(:)
      REAL*8, ALLOCATABLE :: prootmvec(:)
      REAL*8, ALLOCATABLE :: pcapnewvec(:)
      REAL*8, ALLOCATABLE :: hruptkvec(:)
      REAL*8, ALLOCATABLE :: ruptkvec_d(:)
      REAL*8, ALLOCATABLE :: reff(:)
      REAL*8, ALLOCATABLE :: reffg(:)
      REAL*8, ALLOCATABLE :: ruptkg(:)
      REAL*8, ALLOCATABLE :: hruptkg(:)
      REAL*8, ALLOCATABLE :: ruptkg_d(:)
      REAL*8, ALLOCATABLE :: rsurfg_(:)
      REAL*8, ALLOCATABLE :: rsurfgnew(:)
      REAL*8, ALLOCATABLE :: rsoilvec(:)
      REAL*8, ALLOCATABLE :: kunsatnewvec(:)

      REAL*8  :: ys_
      REAL*8  :: yu__
      REAL*8  :: omgu_
      REAL*8  :: omgo_
      REAL*8  :: etm__
      REAL*8  :: esoil__
      REAL*8  :: ysnew
      REAL*8  :: yunew_
      REAL*8  :: wsnew_
      REAL*8  :: io_
      REAL*8  :: iocum
      REAL*8  :: spgfcf__
      REAL*8  :: hspgfcf
      REAL*8  :: infx__
      REAL*8  :: hinfx
      REAL*8  :: hio
      REAL*8  :: hesoil
      REAL*8  :: hetm_
      REAL*8  :: inf_
      REAL*8  :: omgunew
      REAL*8  :: omgonew
      REAL*8  :: error
      REAL*8  :: hinf_

!     * file coes

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

!     * file names

      CHARACTER(80) :: sfile_dailyweather  = 'dailyweather.prn'
      CHARACTER(80) :: sfile_hourlyweather = 'hourlyweather.prn'
      CHARACTER(80) :: sfile_inputpar      = 'input.par'
      CHARACTER(80) :: sfile_resultshourly = 'resultshourly.txt'
      CHARACTER(80) :: sfile_resultsdaily  = 'resultsdaily.txt'
      CHARACTER(80) :: sfile_yearly        = 'yearly.txt'
      CHARACTER(80) :: sfile_rsurfdaily    = 'rsurfdaily.txt'
      CHARACTER(80) :: sfile_delyuhourly   = 'delyuhourly.txt'
      CHARACTER(80) :: sfile_ruptkhourly   = 'ruptkhourly.txt'
      CHARACTER(80) :: sfile_suvechourly   = 'suvechourly.txt'
      CHARACTER(80) :: sfile_soilprofile   = 'soilprofile.par'

      end module vegwatbal
