!***********************************************************************
!*  Transpiration model and layered water balance
!*----------------------------------------------------------------------
!*  Author: Stan Schymanski, CWR, University of Western Australia
!*  03/2006
!*  now at: Max Planck Institute for Biogeochemistry
!*  Email: sschym@bgc-jena.mpg.de
!*  02/2008
!*  Version: big leaf, trees and grass, layered unsaturated zone
!*  optimised root profile, pcg and Jmax25
!*----------------------------------------------------------------------
!*
!* Numbers in the commented parentheses refer to the equation numeration
!* in Schymanski (2007): PhD thesis, University of W.A.
!* and in the document 'equations.pdf' that comes with the documentation.
!*
!*----------------------------------------------------------------------
!*  Copyright (C) 2008  Stan Schymanski
!*
!*    This program is free software: you can redistribute it and/or modify
!*    it under the terms of the GNU General Public License as published by
!*    the Free Software Foundation, either version 3 of the License, or
!*    (at your option) any later version.
!*
!*    This program is distributed in the hope that it will be useful,
!*    but WITHOUT ANY WARRANTY; without even the implied warranty of
!*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!*    GNU General Public License for more details.
!*
!*    You should have received a copy of the GNU General Public License
!*    along with this program.  If not, see http://www.gnu.org/licenses.
!*
!***********************************************************************

      subroutine transpmodel(invar, dim_invar, tp_netass, option1)
      use vom_vegwat_mod
      implicit none

      INTEGER, INTENT(in)    :: dim_invar
      REAL*8,  INTENT(inout) :: tp_netass
      INTEGER, INTENT(in)    :: option1
      REAL*8, DIMENSION(dim_invar), INTENT(in) :: invar

      tp_netass = 0.d0

      call transpmodel_init(invar, dim_invar, option1)

!     * DAILY LOOPS

      do while (nday  .lt. testday)
        nday = nday + 1

        call vom_daily_init()

!     * HOURLY LOOPS (loops through each hour of daily dataset)
      do nhour = 1, 24

      call vom_hourly_init()

!     * calculate gstom, et and ass

      call vom_gstom()

!     * SUB-HOURLY LOOPS

      do while (time .lt. 3600.d0)

!       * setting variables from previous loop

        call vom_subhourly_init()

!       * root water uptake

        call vom_rootuptake()

        if (md_ .gt. 0.d0) then

!         * steady-state tissue water (mqss)

          if (wlayer_ .ge. 1) then
            call vom_mqss(mqss_)
          else
            mqss_ = 0.9d0 * mqx_
          endif
          mqssmin = MIN(mqssmin,mqss_)

!         * transpiration, gstom and tissue water

          call vom_tissue_water_et(tp_netass)
          if (finish .eq. 1) return

        endif

!       * water balance and conditions at next time step

        call vom_subhourly()

        time = time + dt_
        mqnew = mq_ + dmq * dt_

!       * adding up hourly fluxes

        call vom_add_hourly()

!       * END OF HOUR

      enddo

!     * rl does not need to be included here as ass=-rl if j=0 (at night)
      tp_netass = tp_netass + ass_h(2) - 3600.d0 * (cpcc_ + rr_ + tc_) &
     &          + assg_h(2,2) - 3600.d0 * (cpccg(2) + rrg + tcg(2))
      ass_d(:)    = ass_d(:)    + ass_h(:)
      assg_d(:,:) = assg_d(:,:) + assg_h(:,:)
      ruptkvec_d(:) = ruptkvec_d(:) + ruptkvec_h(:)
      ruptkg_d(:)   = ruptkg_d(:)   + ruptkg_h(:)

      if (optmode .eq. 0) then

        call vom_add_daily()
        call vom_write_hourly()

!       * check water balance

        call vom_check_water()
        if (finish .eq. 1) return

      endif
        enddo

!       * END OF DAY

        call transpmodel_daily_step(tp_netass)

      enddo

!     * END OF DAILY LOOPS

      call transpmodel_yearly_step(tp_netass)

      return
      end subroutine transpmodel

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine transpmodel_daily_step (tp_netass)
      use vom_vegwat_mod
      implicit none

      REAL*8,  INTENT(inout) :: tp_netass

      if (optmode .eq. 0) then
        call vom_write_dayyear()
        call vom_add_yearly()
      endif

!     * ADJUSTMENT OF JMAX25 and PC

      call vom_adapt_foliage()

!     * ADJUSTMENT OF ROOT SURFACE

      call vom_adapt_roots()

      if ((nday .eq. testday) .and. (nday .lt. maxday)) then
        if (tp_netass .le. 0.d0) then
!         * estimates how bad the carbon loss would be instead of
!           running through the whole set
          tp_netass = tp_netass / testyear * maxyear
        else
          testday = maxday
        endif
      endif

      return
      end subroutine transpmodel_daily_step

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine transpmodel_yearly_step (tp_netass)
      use vom_vegwat_mod
      implicit none

      REAL*8, INTENT(in) :: tp_netass

      if (optmode .eq. 0) then
        print *,'Cumulative error in water balance (initial Ws+Input-Output-final Ws, in m): ',error
        print *,'Number of times dtsu was limiting: ',dtsu_count
        print *,'Number of times dtmax was limiting: ',dtmax_count

        close(kfile_resultshourly)
        close(kfile_resultsdaily)
        close(kfile_yearly)
        close(kfile_rsurfdaily)
        close(kfile_delyuhourly)
        close(kfile_ruptkhourly)
        close(kfile_suvechourly)
      endif

      if (optmode .eq. 2) then
        open(kfile_model_output, FILE=sfile_model_output, STATUS='replace')
        write(kfile_model_output,'(E13.6)') tp_netass
        close(kfile_model_output)
      endif

      return
      end subroutine transpmodel_yearly_step

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine transpmodel_init_once ()
      use vom_vegwat_mod
      implicit none

!     * Reading input parameter

      call vom_read_input()

!     * allocate vector sizes

      call vom_alloc()

!     * File opening (saving climate and gstom ass data)

      if (optmode .eq. 0) call vom_open_output()

!     * PARAMETER READING FROM SOILPROFILE.PAR

      call vom_get_soilprofile()

!     * Climate and Calendar data reading

      call vom_get_hourly_clim()

      return
      end subroutine transpmodel_init_once

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine transpmodel_init (vom_invar, vom_insize, vom_command)
      use vom_vegwat_mod
      implicit none

      INTEGER, INTENT(in)    :: vom_insize
      INTEGER, INTENT(in)    :: vom_command
      REAL*8, DIMENSION(vom_insize), INTENT(in) :: vom_invar

!dd   if (vom_insize .lt. 8) then
      if (vom_insize .lt. 6) then
        write(0,*) "ERROR: Number of input parameters less than 6."
        stop
      endif

      if (vom_command .eq. 2) then
        optmode = 0
      elseif (vom_command .eq. 3) then
        optmode = 2
      else
        optmode = 1
      endif

!*----------------------------------------------------------------------
!*     Optimised parameters reading from vom_invar
!*----------------------------------------------------------------------

      lambdagfac = vom_invar(1)
      wsgexp     = vom_invar(2)
      lambdafac  = vom_invar(3)
      wsexp      = vom_invar(4)
      pc_        = vom_invar(5)
      rootdepth  = vom_invar(6)
!dd   mdstore    = vom_invar(7)
!dd   rgdepth    = vom_invar(8)

      if (parsaved .ne. 1) then
        call transpmodel_init_once ()
        parsaved = 1
      endif

!***********************************************************************
!*  Calculation of vegetation parameters
!*
!*----------------------------------------------------------------------
!*  Equations in equations.pdf
!*  (numeration in the commented parentheses)
!***********************************************************************

!*-----Initial values---------------------------------------------------

      call vom_init_vegpar()

      nday = 0
      testday = testyear * 365
      if (optmode .eq. 0 .or. optmode .eq. 2) testday = maxday
      finish = 0

      return
      end subroutine transpmodel_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----PARAMETER READING FROM INPUT.PAR---------------------------------

      subroutine vom_read_input ()
      use vom_vegwat_mod
      implicit none

      INTEGER :: iostat

!     * Input of variable parameters from the parameter file

      open(kfile_inputpar, FILE=sfile_inputpar, STATUS='old')

      read(kfile_inputpar,*) alpha
      read(kfile_inputpar,*) cpccf
      read(kfile_inputpar,*) tcf
      read(kfile_inputpar,*) maxyear
      read(kfile_inputpar,*) testyear
      read(kfile_inputpar,*) ha_
      read(kfile_inputpar,*) hd_
      read(kfile_inputpar,*) toptfac
      read(kfile_inputpar,*) toptstart
      read(kfile_inputpar,*) rlratio

!     * Catchment parameters

      read(kfile_inputpar,*) lat_
      read(kfile_inputpar,*) cz
      read(kfile_inputpar,*) cgs
      read(kfile_inputpar,*) zr_
      read(kfile_inputpar,*) go_

!     * Soil parameters

      read(kfile_inputpar,*) ksat_
      read(kfile_inputpar,*) thetar_
      read(kfile_inputpar,*) thetas_
      read(kfile_inputpar,*) nvg_
      read(kfile_inputpar,*) avg_

!     * Vertical Resolution

      read(kfile_inputpar,*) delz_

!     * Vegetation Parameters

      read(kfile_inputpar,*) mdf
      read(kfile_inputpar,*) mqxf
      read(kfile_inputpar,*) rrootm
      read(kfile_inputpar,*) rsurfmin
      read(kfile_inputpar,*) rsurfinit
      read(kfile_inputpar,*) rootrad
      read(kfile_inputpar,*) prootmg
      read(kfile_inputpar,*) growthmax
      read(kfile_inputpar,*) mdstore
      read(kfile_inputpar,*) rgdepth
      read(kfile_inputpar,*) firstyear
      read(kfile_inputpar,*) lastyear
!     read(kfile_inputpar,*) write_h
      close(kfile_inputpar)

      epsln_ = thetas_ - thetar_        ! epsilon, porosity see Reggiani (2000)
      mvg_ = 1.d0 - (1.d0 / nvg_)       ! van Genuchten soil parameter m

      maxday = maxyear * 365
      maxhour = maxday * 24

!     * The file soilprofile.par contain information about thickness and
!       soil properties in each soil layer, with the layer number in the
!       first column.

      maxlayer = 0

      open(kfile_soilprofile, FILE=sfile_soilprofile,                  &
     &                        STATUS='old', IOSTAT=iostat)
      if (iostat .eq. 0) then
        read(kfile_soilprofile,*) maxlayer
      endif
      close(kfile_soilprofile)

!     * number of soil layers maxlayer assuming same thickness everywhere
      if (maxlayer .eq. 0) maxlayer = ceiling(cz / delz_)

      return
      end subroutine vom_read_input

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----allocate vector sizes--------------------------------------------

      subroutine vom_alloc ()
      use vom_vegwat_mod
      implicit none

      allocate(dayyear(maxday))
      allocate(fday(maxday))
      allocate(fmonth(maxday))
      allocate(fyear(maxday))

      allocate(tairmax(maxday))
      allocate(tairmin(maxday))
      allocate(rainvec(maxday))
      allocate(srad__(maxday))
      allocate(vpvec(maxday))
      allocate(press(maxday))
      allocate(cavec(maxday))

      allocate(par_h(maxhour))
      allocate(vd_h(maxhour))
      allocate(tair_h(maxhour))
      allocate(rain_h(maxhour))
      allocate(ca_h(maxhour))

      allocate(parvec(maxday))

      allocate(pcapvec(maxlayer))
      allocate(suvec_(maxlayer))
      allocate(ruptkvec(maxlayer))
      allocate(sunewvec(maxlayer))
      allocate(kunsatvec(maxlayer))
      allocate(delzvec(maxlayer))
      allocate(rsurfvec(maxlayer))
      allocate(rsurfnewvec(maxlayer))
      allocate(qblvec(maxlayer))
      allocate(dsuvec(maxlayer))
      allocate(phydrostaticvec(maxlayer))
      allocate(prootmvec(maxlayer))
      allocate(pcapnewvec(maxlayer))
      allocate(ruptkvec_d(maxlayer))
      allocate(ruptkvec_h(maxlayer))
      allocate(ruptkg_h(maxlayer))
      allocate(ruptkg_d(maxlayer))
      allocate(reff(maxlayer))
      allocate(reffg(maxlayer))
      allocate(ruptkg(maxlayer))
      allocate(rsurfg_(maxlayer))
      allocate(rsurfgnew(maxlayer))
      allocate(rsoilvec(maxlayer))
      allocate(kunsatnewvec(maxlayer))

      allocate(ksatvec(maxlayer))
      allocate(thetasvec(maxlayer))
      allocate(thetarvec(maxlayer))
      allocate(sueqvec(maxlayer))
      allocate(cH2Ol_s(maxlayer))
      allocate(iovec(maxlayer))
      allocate(avgvec(maxlayer))
      allocate(nvgvec(maxlayer))
      allocate(mvgvec(maxlayer))

      return
      end subroutine vom_alloc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----File opening (saving climate and gstom ass data)-----------------

      subroutine vom_open_output ()
      use vom_vegwat_mod
      implicit none

      open(kfile_resultshourly, FILE=sfile_resultshourly, STATUS='replace')
      write(kfile_resultshourly,'(a6,a7,a7,a7,a7,22a15)') 'year',      &
     &  'month', 'day', 'dcum', 'hour', 'rain', 'tair', 'par', 'vd',   &
     &  'esoil', 'pc', 'jmax25_t', 'jmax25_g', 'mq', 'rl', 'lambda_t', &
     &  'lambda_g', 'rr', 'ass_t', 'ass_g', 'het_t', 'het_g', 'su_1',  &
     &  'ys', 'Ws', 'spgfcf', 'infx'

      open(kfile_resultsdaily, FILE=sfile_resultsdaily, STATUS='replace')
      write(kfile_resultsdaily,'(a6,a7,a7,a7,a7,25a15)') 'year',       &
     &  'month', 'day', 'dcum', 'hour', 'rain', 'tairmax', 'tairmin',  &
     &  'par', 'vd', 'esoil', 'jmax25_t', 'jmax25_g', 'pc', 'rlt+rlg', &
     &  'lambda_t', 'lambda_g', 'rr_t', 'rr_g', 'ass_t', 'ass_g',      &
     &  'su_avg', 'ys', 'ws', 'spgfcf', 'infx', 'etm_t', 'etm_g',      &
     &  'su_1', 'topt'

      open(kfile_yearly, FILE=sfile_yearly, STATUS='replace')
      write(kfile_yearly,'(a6,18a16)') "year", "rain_y",               &
     &  "par_y", "srad_y", "vd_y", "esoil_y", "etm_y", "etmg_y",       &
     &  "assg_y", "rlg_y", "rrg_y", "cpccg_y", "tcg_y",                &
     &  "etmt_y", "asst_y", "rlt_y", "rrt_y", "cpcct_y", "tc_y"

      open(kfile_rsurfdaily, FILE=sfile_rsurfdaily, STATUS='replace')
      write(kfile_rsurfdaily,*) ' year', ' month', ' day', '   dcum',  &
     &  '  rsurfsublayer'

      open(kfile_delyuhourly, FILE=sfile_delyuhourly, STATUS='replace')
      write(kfile_delyuhourly,*) ' year', ' month', ' day', '   dcum',  &
     &  ' hour', '  delyusublayer'

      open(kfile_ruptkhourly, FILE=sfile_ruptkhourly, STATUS='replace')
      write(kfile_ruptkhourly,*) ' year', ' month', ' day', '   dcum', &
     &  ' hour', '  delyusublayer'

      open(kfile_suvechourly, FILE=sfile_suvechourly, STATUS='replace')
      write(kfile_suvechourly,*) ' year', ' month', ' day', '   dcum', &
     &  ' hour', '  susublayer'

      return
      end subroutine vom_open_output

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----PARAMETER READING FROM SOILPROFILE.PAR---------------------------

      subroutine vom_get_soilprofile ()
      use vom_vegwat_mod
      implicit none

      INTEGER :: iostat, j

!     * The file soilprofile.par can contain information about thickness
!       and soil properties in each soil layer, with the layer number in
!       the first column.

      open(kfile_soilprofile, FILE=sfile_soilprofile,                  &
     &                        STATUS='old', IOSTAT=iostat)
      if (iostat .eq. 0) then
        do j = 1, maxlayer
          read(kfile_soilprofile,*) maxlayer, delzvec(j), ksatvec(j),  &
     &      nvgvec(j), avgvec(j), thetasvec(j), thetarvec(j)
          mvgvec(j) = 1.d0 - (1.d0 / nvgvec(j))  ! van Genuchten soil parameter m
        enddo
      else
        delzvec(:)   = delz_
        ksatvec(:)   = ksat_
        nvgvec(:)    = nvg_
        avgvec(:)    = avg_
        thetasvec(:) = thetas_
        thetarvec(:) = thetar_
        mvgvec(:)    = mvg_
      endif
      close(kfile_soilprofile)

      return
      end subroutine vom_get_soilprofile

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----Climate and Calendar data reading--------------------------------

      subroutine vom_get_hourly_clim ()
      use vom_vegwat_mod
      implicit none

      INTEGER :: ii, i, h, oldh, stat
      INTEGER :: dummyint1, dummyint2, dummyint3, dummyint4
      LOGICAL :: exist

      inquire(FILE=sfile_hourlyweather, EXIST=exist)

!     * Creating hourly climate data from daily data

      if (.not. exist) then
        open(kfile_dailyweather, FILE=sfile_dailyweather,              &
     &                           STATUS='old', IOSTAT=stat)
        read(kfile_dailyweather,*)
        do i = 1, maxday
          read(kfile_dailyweather,'(4i8,7f8.2)') dayyear(i), fday(i),  &
     &      fmonth(i), fyear(i), tairmax(i), tairmin(i), rainvec(i),   &
     &      srad__(i), vpvec(i), press(i), cavec(i)
        enddo
        close(kfile_dailyweather)

!       * Calculation of derived parameters

        call vom_calc_derived()

        if (write_h == 1) exist=.TRUE.
      endif

!     * Reading hourly climate data if available

      if (exist) then
        open(kfile_hourlyweather, FILE=sfile_hourlyweather,            &
     &                            STATUS='old', IOSTAT=stat)
        read(kfile_hourlyweather,*)
        ii = 1
        oldh = 99
        do i = 1, maxhour
          read(kfile_hourlyweather,'(5i8,5e11.3)') h, dummyint1,       &
     &      dummyint2, dummyint3, dummyint4, tair_h(i), vd_h(i),       &
     &      par_h(i), rain_h(i), ca_h(i)
          if (h .lt. oldh) then
            dayyear(ii) = dummyint1
            fday(ii)    = dummyint2
            fmonth(ii)  = dummyint3
            fyear(ii)   = dummyint4
            ii = ii + 1
          endif
          oldh = h
        enddo
        close(kfile_hourlyweather)
      endif

      return
      end subroutine vom_get_hourly_clim

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----Calculation of derived parameters--------------------------------

      subroutine vom_calc_derived ()
      use vom_vegwat_mod
      implicit none

      INTEGER :: in, ik, ii
      REAL*8  :: sunr, suns
      REAL*8  :: tairmean
      REAL*8  :: dtair
      REAL*8  :: daylength              ! Day length (hours)

      if (write_h == 1) then
        open(kfile_hourlyweather, FILE=sfile_hourlyweather, STATUS='new')
        write(kfile_hourlyweather,'(5a8,5a11)') 'hour', 'dayyear', 'day', &
     &    'month', 'year', 'tair_h', 'vd_h', 'par_h', 'rain_h', 'ca_h'
      endif

      parvec(:) = 2.0804d0 * srad__(:)       ! (Out[17]), par in mol/m2 if srad was MJ/m2
      do in = 1, maxday
        daylength = 12.d0 - 7.639437d0 * ASIN((0.397949d0              &
     &            * COS(0.172142d0 + 0.017214d0 * dayyear(in))         &
     &            * TAN(0.017453d0 * lat_))                            &
     &            / SQRT(0.920818d0 - 0.079182d0                       &
     &            * COS(0.344284d0 + 0.034428d0 * dayyear(in))))  ! (Out[22]), in hours
!       * sets time of sunrise and sunset
        sunr = 12d0 - 0.5d0 * daylength
        suns = 12d0 + 0.5d0 * daylength
        tairmean = (tairmax(in) + tairmin(in)) / 2.d0
        dtair = tairmax(in) - tairmin(in)
        vp_ = vpvec(in) * 100.d0             ! vp in Pa

!       * Loop through every hour of day, where ik=hour
        do ik = 1, 24
          ii = in * 24 + ik - 24
!         * (derived from 3.52+3.53) (Out[38], accounts for diurnal
!           variation in air temperature
          tair__ = tairmean + dtair * (0.0138d0                        &
     &           * COS(3.513d0 - ((-1.d0 + ik) * p_pi) / 3.d0) + 0.0168d0 &
     &           * COS(0.822d0 - ((-1.d0 + ik) * p_pi) / 4.d0) + 0.0984d0 &
     &           * COS(0.360d0 - ((-1.d0 + ik) * p_pi) / 6.d0) + 0.4632d0 &
     &           * COS(3.805d0 - ((-1.d0 + ik) * p_pi) / 12.d0))
          tair_h(ii) = tair__

          ca_h(ii) = cavec(in)

!         vd__ = 0.006028127d0 * 2.718282d0 ** ((17.27d0 * tair__)     &
!    &       / (237.3d0 + tair__)) - 9.869233d-6 * vp_  ! (derived from 3.54+3.55) (Out[52]), accounts for diurnal variation in vapour deficit
          vd__ = (((0.6108d0 * p_E ** (17.27d0 * tair__ / (tair__      &
     &         + 237.3d0))) * 1000) - vp_) / (press(in) * 100.d0)
          if (vd__ .le. 0.d0) vd__ = 0.d0
          vd_h(ii) = vd__

!         * average rainfall in hour ii (m/s)
          rain_h(ii) = rainvec(in) / (24.d0 * 3600.d0 * 1000.d0)

          if (sunr .le. ik .and. ik + 1 .le. suns) then
!           * ([Out30]), in mol/m2/s (derived from 3.51) accounts for
!             diurnal variation in global irradiance
            par_h(ii) = (-0.000873d0 * parvec(in) * COS(0.017453d0 * lat_) &
     &                * SQRT(0.920818d0 - 0.079182d0                   &
     &                * COS(0.034428d0 * (10.d0 + dayyear(in))))       &
     &                * COS(0.2618d0 * ik) - 0.000347d0 * parvec(in)   &
     &                * COS(0.017214d0 * (10.d0 + dayyear(in)))        &
     &                * SIN(0.017453d0 * lat_)) / (-1.250192d0 * daylength &
     &                * COS(0.017214d0 * (10.d0 + dayyear(in)))        &
     &                * SIN(0.017453d0 * lat_) + 24.d0                 &
     &                * COS(0.017453d0 * lat_)                         &
     &                * SQRT(0.920818d0 - 0.079182d0                   &
     &                * COS(0.034428d0 * (10.d0 + dayyear(in))))       &
     &                * (1.d0 - (0.158363d0                            &
     &                * COS(0.017214d0 * (10.d0 + dayyear(in))) ** 2.d0 &
     &                * TAN(0.017453d0 * lat_) ** 2.d0)                &
     &                / (0.920818d0 - 0.079182d0 * COS(0.034428d0      &
     &                * (10.d0 + dayyear(in))))) ** 0.5d0)
          else
            par_h(ii) = 0.d0
          endif

          if (write_h == 1) then
            write(kfile_hourlyweather,'(5i8,5e11.3)') ik, dayyear(in), &
     &        fday(in), fmonth(in), fyear(in), tair_h(ii), vd_h(ii),   &
     &        par_h(ii), rain_h(ii), ca_h(ii)
          endif

        enddo
      enddo

      if (write_h == 1) close(kfile_hourlyweather)

      return
      end subroutine vom_calc_derived

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----Initial values---------------------------------------------------

      subroutine vom_init_vegpar ()
      use vom_vegwat_mod
      implicit none

      REAL*8  :: dummy

      topt_ = toptstart

!     * Set soil moisture and vegetation parameters to initial conditions

      call waterbalance_init()

      if (rootdepth .gt. cz) then
        write(*,*) 'Root depth greater than soil depth'
        rootdepth = cz
      endif

      wsold  = SUM(cH2Ol_s(:))               ! initial soil water storage
      wsnew = wsold

!     * Set vegetation parameters

      md_   = pc_ * mdf + mdstore
      mqx_  = md_ * mqxf
      mqnew = 0.95d0 * mqx_                  ! initial wood water storage
      mqold = mqnew
      rsurfnewvec(:) = 0.d0

!     * Determining the position of the bottom of the tree root zone

      pos_ = 0
      dummy = 0
      do while (rootdepth .gt. dummy)
        pos_ = pos_ + 1
        dummy = dummy + delzvec(pos_)
      enddo

!     * Determining the position of the bottom of the tree root zone

      posg = 0
      dummy = 0
      do while (rgdepth .gt. dummy)
        posg = posg + 1
        dummy = dummy + delzvec(posg)
      enddo

      rsurfgnew(1:posg) = rsurfinit * delzvec(1:posg)
      rsurfgnew(posg+1:maxlayer) = 0.d0
      if (posg .gt. wlayernew) then
        rsurfgnew(wlayernew+1:posg) = rsurfmin * delzvec(wlayernew+1:pos_)
      endif

!     * root surface density (root surface area/soil volume) in each sublayer

      rsurfnewvec(1:pos_) = rsurfinit * delzvec(1:pos_)
      if (pos_ .gt. wlayernew) then
        rsurfnewvec(wlayernew+1:pos_) = rsurfmin * delzvec(wlayernew+1:pos_)
      endif
      jmax25_(2)   = 0.0003d0
      jmax25g(2)   = 0.0003d0
      pcgmin       = 0.02d0                  ! minimum grass pc; initial point for growth
      pcg_(2)      = MIN(1.d0 - pc_, pcgmin)
      pcg_(:)      = pcg_(2) + (/-0.02,0.0,0.02/)  ! vector with values varying by 1%
      pcg_(3)      = MIN(MAX(pcgmin, pcg_(3)), 1.d0 - pc_)
      rootlim(:,:) = 0.d0

!     * Direct costs

!     * (3.38)  foliage tunrover costs, assuming crown LAI of 2.5
      tc_ = tcf * pc_ * 2.5d0

!     * Setting yearly, daily and hourly parameters

      nyear          = fyear(1)
      rain_y         = 0.d0
      par_y          = 0.d0
      srad_y         = 0.d0
      vd_y           = 0.d0
      etm_y          = 0.d0
      esoil_y        = 0.d0             ! = yearly esoil
      ruptkvec_d(:)  = 0.d0
      ruptkg_d(:)    = 0.d0
      ass_d(:)       = 0.d0
      assg_d(:,:)    = 0.d0
      ioacum         = 0.d0
!     * for grasses
      etmg_y         = 0.d0
      assg_y         = 0.d0
      rlg_y          = 0.d0
      rrg_y          = 0.d0
      cpccg_y        = 0.d0
      tcg_y          = 0.d0
!     * for trees
      etmt_y         = 0.d0
      asst_y         = 0.d0
      rlt_y          = 0.d0
      rrt_y          = 0.d0
      cpcct_y        = 0.d0
      tct_y          = 0.d0

      return
      end subroutine vom_init_vegpar

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine vom_daily_init ()
      use vom_vegwat_mod
      implicit none

      rsurfvec(:) = rsurfnewvec(:)
      rsurfg_(:)  = rsurfgnew(:)
      lambda_     = lambdafac * (SUM(pcapnewvec(1:pos_)) / pos_) ** wsexp  ! (3.45)
      lambdag_    = lambdagfac * pcapnewvec(1) ** wsgexp  ! (3.44)
!     * vector with values varying by 1%
      jmax25_(:)  = jmax25_(2) * (/0.99,1.0,1.01/)
!     * making sure that the values don't become too low, otherwise
!       they could never pick up again
      jmax25_(:)  = MAX(jmax25_(:), 50.0d-6)
      jmax25g(:)  = jmax25g(2) * (/0.99,1.0,1.01/)
      jmax25g(:)  = MAX(jmax25g(:), 50.0d-6)
      pcg_(:)     = pcg_(2) + (/-0.02,0.0,0.02/)  ! vector with values varying by 1%
      pcg_(:)     = MAX(pcg_(:), 0.d0)
      pcg_(3)     = MIN(MAX(pcgmin, pcg_(3)), 1.d0 - pc_)
!     * (3.38) foliage turnover costs, assuming LAI/pc of 2.5
      tcg(:)      = tcf * pcg_(:) * 2.5d0
!     * (3.40), (Out[190])  root respiration [mol/s]
      rr_         = 2.55d-7 * SUM(rsurfvec(1:pos_))

      if (pos_ .gt. wlayernew) then
!       * (3.42, 2.45e-10 from (Out[165])) costs of water distribution and storage
        cpcc_ = cpccf * pc_ * rootdepth + mdstore * 2.45d-10
      else
        cpcc_ = cpccf * pc_ * SUM(delzvec(1:pos_)) + mdstore * 2.45d-10
      endif

      if (wlayernew .lt. posg) then
        cpccg(:) = cpccf * pcg_(:) * rgdepth  ! (3.42) water transport costs
      else
        cpccg(:) = cpccf * pcg_(:) * SUM(delzvec(1:posg))
      endif

!     * (3.40), (Out[190]) root respiration grasses [mol/s]
      rrg = 2.55d-7 * SUM(rsurfg_(1:posg))
!     * resetting the minimum steady-state tissue water content to
!       its maximum value
      mqssmin = mqx_

!     * used for daily recalculation
      if (optmode .eq. 0) then
        tairmax(nday) = -9999.d0
        tairmin(nday) =  9999.d0
        rainvec(nday) =     0.d0
        parvec(nday)  =     0.d0
        srad__(nday)  =     0.d0
      endif

      if (optmode .eq. 0) then
        vd_d     = 0.d0
        etm_d    = 0.d0
        etmg_d   = 0.d0
        esoil_d  = 0.d0
        spgfcf_d = 0.d0
        infx_d   = 0.d0
        rl_d     = 0.d0
        rlg_d    = 0.d0
      endif

      return
      end subroutine vom_daily_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine vom_hourly_init ()
      use vom_vegwat_mod
      implicit none

      INTEGER :: ii

      ii         = nday * 24 + nhour - 24
      rain__     = rain_h(ii)
      tair__     = tair_h(ii)
      vd__       = vd_h(ii)
      par__      = par_h(ii)
      ca__       = ca_h(ii) / 1.0d6

!     * (Out[274], derived from (3.25))
      gammastar_ = 0.00004275d0                                        &
     &           * p_E ** ((18915.d0 * (-25.d0 + tair_h(ii)))          &
     &           / (149.d0 * p_R_ * (273.d0 + tair_h(ii))))

!     * (Out[310], derived from (3.26)) Temperature dependence of Jmax
      jmax__(:) = (p_E ** ((ha_ * (-25.d0 + tair__) * (-273.d0 + topt_ &
     &          + 273.d0 * p_R_ * topt_)) / ((25.d0 + 273.d0 * p_R_    &
     &          * topt_) * (tair__ + 273.d0 * p_R_ * topt_))) * ((-1.d0 &
     &          + p_E ** (-(hd_ * (-298.d0 + topt_)) / (25.d0 + 273.d0 &
     &          * p_R_ * topt_))) * ha_ + hd_) * jmax25_(:)) / ((-1.d0 &
     &          + p_E ** ((hd_ * (273.d0 + tair__ - topt_)) / (tair__  &
     &          + 273.d0 * p_R_ * topt_))) * ha_ + hd_)
!     * (3.24), (Out[312])
      rl__(:) = ((ca__ - gammastar_) * pc_ * jmax__(:) * rlratio)      &
     &        / (4.d0 * (ca__ + 2.d0 * gammastar_) * (1.d0 + rlratio))
!     * (Out[310], derived from (3.26)) Temperature dependence of Jmax
      jmaxg__(:) = (p_E ** ((ha_ * (-25.d0 + tair__) * (-273.d0 + topt_ &
     &           + 273.d0 * p_R_ * topt_)) / ((25.d0 + 273.d0 * p_R_   &
     &           * topt_) * (tair__ + 273.d0 * p_R_ * topt_))) * ((-1.d0 &
     &           + p_E ** (-(hd_ * (-298.d0 + topt_)) / (25.d0         &
     &           + 273.d0 * p_R_ * topt_))) * ha_ + hd_) * jmax25g(:)) &
     &           / ((-1.d0 + p_E ** ((hd_ * (273.d0 + tair__ - topt_)) &
     &           / (tair__ + 273.d0 * p_R_ * topt_))) * ha_ + hd_)
      rlg__(1,:) = ((ca__ - gammastar_) * pcg_(1) * jmaxg__(:)         &
     &           * rlratio) / (4.d0 * (ca__ + 2.d0 * gammastar_)       &
     &           * (1.d0 + rlratio))  ! (3.24), (Out[312])
      rlg__(2,:) = ((ca__ - gammastar_) * pcg_(2) * jmaxg__(:)         &
     &           * rlratio) / (4.d0 * (ca__ + 2.d0 * gammastar_)       &
     &           * (1.d0 + rlratio))  ! (3.24), (Out[312])
      rlg__(3,:) = ((ca__ - gammastar_) * pcg_(3) * jmaxg__(:)         &
     &           * rlratio) / (4.d0 * (ca__ + 2.d0 * gammastar_)       &
     &           * (1.d0 + rlratio))  ! (3.24), (Out[312])

!     * daily recalculation for resultsdaily
      if (optmode .eq. 0) then
        rainvec(nday) = rainvec(nday) + rain__ * 3600.d0 * 1000.d0  ! mm/d
        parvec(nday)  = parvec(nday) + par__ * 3600.d0  ! in mol/m2/d
        srad__(nday)  = parvec(nday) / 2.0804d0  ! MJ/m2/d

        if (tair__ .gt. tairmax(nday)) tairmax(nday) = tair__
        if (tair__ .lt. tairmin(nday)) tairmin(nday) = tair__
      endif

      if (optmode .eq. 0) then
        mqold    = mqnew
        spgfcf_h = 0.d0
        infx_h   = 0.d0
        io_h     = 0.d0
        esoil_h  = 0.d0
        etm_h    = 0.d0
        etmg_h   = 0.d0
        ruptk_h  = 0.d0
      endif

      time          = 0.d0
      ass_h(:)      = 0.d0                   ! hourly assimilation
      assg_h(:,:)   = 0.d0
      ruptkvec_h(:) = 0.d0
      ruptkg_h(:)   = 0.d0

      return
      end subroutine vom_hourly_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----calculate gstom, et and ass -------------------------------------

      subroutine vom_gstom ()
      use vom_vegwat_mod
      implicit none

      REAL*8 :: cond1, cond2
      REAL*8 :: cond3(3,3)
      REAL*8 :: part1, part2, part3, part4, part5
      REAL*8 :: part6, part7, part8, part9

      if (par__ .gt. 0.d0) then
!       * adaptation of topt to air temperature during sunlight
        topt_ = topt_ + toptfac * (tair__ + 273.d0 - topt_)
        jact_(:)   = (1.d0 - p_E ** (-(alpha * par__) / jmax__(:)))    &
     &             * jmax__(:) * pc_         ! (3.23), (Out[311])
        jactg(1,:) = (1.d0 - p_E ** (-(alpha * par__) / jmaxg__(:)))   &
     &             * jmaxg__(:) * pcg_(1)    ! (3.23), (Out[311])
        jactg(2,:) = (1.d0 - p_E ** (-(alpha * par__) / jmaxg__(:)))   &
     &             * jmaxg__(:) * pcg_(2)    ! (3.23), (Out[311])
        jactg(3,:) = (1.d0 - p_E ** (-(alpha * par__) / jmaxg__(:)))   &
     &             * jmaxg__(:) * pcg_(3)    ! (3.23), (Out[311])

        cond1      = (2.d0 * p_a * vd__) / (ca__ + 2.d0 * gammastar_)
        cond2      = (4.d0 * ca__ * rl__(2) + 8.d0 * gammastar_        &
     &             * rl__(2)) / (ca__ - gammastar_)
        cond3(:,:) = (4.d0 * ca__ * rlg__(:,:) + 8.d0 * gammastar_     &
     &             * rlg__(:,:)) / (ca__ - gammastar_)

        if (vd__ .gt. 0.d0 .and. lambda_ .gt. cond1 .and. jact_(2) .gt. cond2) then

!          gstom__ = MAX(0.d0, (0.25d0 * (p_a * (ca__ * (jact_(2) - 4.d0 &
!     &            * rl__(2)) - 4.d0 * gammastar_ * (jact_(2) + 2.d0    &
!     &            * rl__(2))) * vd__ * (ca__ * lambda_ + 2.d0          &
!     &            * gammastar_ * lambda_ - p_a * vd__)                 &
!     &            + 1.7320508075688772d0 * SQRT(p_a * gammastar_       &
!     &            * jact_(2) * (ca__ * (jact_(2) - 4.d0 * rl__(2))     &
!     &            - gammastar_ * (jact_(2) + 8.d0 * rl__(2))) * vd__   &
!     &            * (ca__ * lambda_ + 2.d0 * gammastar_ * lambda_      &
!     &            - 2.d0 * p_a * vd__) ** 2.d0 * (ca__ * lambda_       &
!     &            + 2.d0 * gammastar_ * lambda_ - p_a * vd__))))       &
!     &            / (p_a * (ca__ + 2.d0 * gammastar_) ** 2.d0 * vd__   &
!     &            * (ca__ * lambda_ + 2.d0 * gammastar_ * lambda_      &
!     &            - p_a * vd__)))

          part1 = ca__ + 2.d0 * gammastar_
          part2 = part1 * lambda_ - p_a * vd__
          part3 = p_a * vd__ * part2

          part4 = ca__ * (jact_(2) - 4.d0 * rl__(2))
          part5 = gammastar_ * jact_(2)
          part6 = gammastar_ * 8.d0 * rl__(2)
          part7 = part4 - part5 - part6

          part8 = Sqrt(part5 * part7 * (part2 - p_a * vd__) ** 2.d0 * part3)
          part9 = part7 - 3.d0 * part5 + 1.7320508075688772d0 * part8 / part3

          gstom__ = 0.25d0 * part9 / part1**2.d0
          gstom__ = Max(0.d0, gstom__)       ! (Out[314])

        else
          gstom__ = 0.d0
        endif
        transp_ = p_a * vd__ * gstom__       ! (3.28) transpiration rate in mol/s
        etm__ = (transp_ * 18.d0) / (10.d0 ** 6.d0)  ! transpiration rate in m/s

        where (vd__ .gt. 0.d0 .and. lambdag_ .gt. cond1 .and. jactg(:,:) .gt. cond3(:,:))
          gstomg__(:,:) = MAX(0.d0,(0.25d0 * (p_a * (ca__ * (jactg(:,:) &
     &                  - 4.d0 * rlg__(:,:)) - 4.d0 * gammastar_       &
     &                  * (jactg(:,:) + 2.d0 * rlg__(:,:))) * vd__     &
     &                  * (ca__ * lambdag_ + 2.d0 * gammastar_         &
     &                  * lambdag_ - p_a * vd__) + 1.7320508075688772d0 &
     &                  * SQRT(p_a * gammastar_ * jactg(:,:) * (ca__   &
     &                  * (jactg(:,:) - 4.d0 * rlg__(:,:))             &
     &                  - gammastar_ * (jactg(:,:) + 8.d0              &
     &                  * rlg__(:,:))) * vd__ * (ca__ * lambdag_       &
     &                  + 2.d0 * gammastar_ * lambdag_ - 2.d0 * p_a    &
     &                  * vd__) ** 2.d0 * (ca__ * lambdag_ + 2.d0      &
     &                  * gammastar_ * lambdag_ - p_a * vd__))))       &
     &                  / (p_a * (ca__ + 2.d0 * gammastar_) ** 2.d0    &
     &                  * vd__ * (ca__ * lambdag_ + 2.d0 * gammastar_  &
     &                  * lambdag_ - p_a * vd__)))  ! (Out[314])
        elsewhere
          gstomg__(:,:) = 0.d0
        endwhere
        transpg(:,:) = p_a * vd__ * gstomg__(:,:)  ! (3.28) transpiration rate in mol/s
        etmg__(:,:) = (transpg(:,:) * 18.d0) / (10.d0 ** 6.d0)  ! transpiration rate in m/s
      else
        par__         = 0.d0
        jact_(:)      = 0.d0
        gstom__       = 0.d0
        etm__         = 0.d0
        jactg(:,:)    = 0.d0
        gstomg__(:,:) = 0.d0
        etmg__(:,:)   = 0.d0
      endif

      return
      end subroutine vom_gstom

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*----- setting variables from previous loop----------------------------

      subroutine vom_subhourly_init ()
      use vom_vegwat_mod
      implicit none

      if (wlayernew .lt. pos_) then
        rsurfvec(wlayernew+1:pos_) = rsurfmin * delz_
      endif
      if (wlayernew .lt. posg) then
        rsurfg_(wlayernew+1:posg) = rsurfmin * delz_
      endif

      mq_          = mqnew
      ys_          = ysnew
      wlayer_      = wlayernew
      suvec_(:)    = sunewvec(:)
      pcapvec(:)   = pcapnewvec(:)
      kunsatvec(:) = kunsatnewvec(:)

      return
      end subroutine vom_subhourly_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----root water uptake------------------------------------------------

      subroutine vom_rootuptake ()
      use vom_vegwat_mod
      implicit none

      INTEGER :: i

      if (wlayernew .ge. 1) then
        postemp_ = MIN(pos_, wlayer_)
        do i = 1, postemp_
          phydrostaticvec(i) = (i - 0.5d0) * delz_  ! (Out[238]) hydrostatic head for (3.34)
        enddo
        if (md_ .gt. 0.d0) then
          prootmvec(1:postemp_) = (p_mpbar * (-mq_ + mqx_) * (750.d0   &
     &                          - (750.d0 * mqx_) / (md_ + mqx_)       &
     &                          + (md_ + mqx_) / mqx_)) / (md_ + mqx_) &
     &                          - phydrostaticvec(1:postemp_)  ! (Out[239])
        else
!         * set tissue suction to the same as in grasses, if no storage capacity
          prootmvec(1:postemp_) = prootmg
        endif

!       * soil resistance, (Out[ 241] with svolume=delzvec(1:postemp_)); derived from (3.32)
        rsoilvec(1:postemp_) = SQRT(p_pi / 2.d0) * SQRT((rootrad       &
     &                       * delzvec(1:postemp_))                    &
     &                       / rsurfvec(1:postemp_))                   &
     &                       / kunsatvec(1:postemp_)

!       * root water uptake, Chapter 3.3.3.3 (Out[242])
        if (md_ .gt. 0.d0) then
          ruptkvec(1:postemp_) = (-pcapvec(1:postemp_)                 &
     &                         + prootmvec(1:postemp_))                &
     &                         * rsurfvec(1:postemp_) / (rrootm        &
     &                         + rsoilvec(1:postemp_))
          ruptkvec(postemp_+1:maxlayer) = 0.d0
        else  ! if no storage, uptake happens only when etm__>0

          if (etm__ .gt. 0.d0) then
            ruptkvec(1:postemp_) = MAX(0.d0,(-pcapvec(1:postemp_)      &
     &                           + prootmvec(1:postemp_))              &
     &                           * rsurfvec(1:postemp_) / (rrootm      &
     &                           + rsoilvec(1:postemp_)))
            ruptkvec(postemp_+1:maxlayer) = 0.d0

            if (SUM(ruptkvec(:)) .gt. 0.d0) then
              if (etm__ .gt. SUM(ruptkvec(:))) then
                changef = 1.d0
                etm__   = SUM(ruptkvec(:))
                transp_ = etm__ * 55555.555555555555d0  ! (Out[249]) mol/s=m/s*10^6 g/m/(18g/mol)
                gstom__ = transp_ / (p_a * vd__)
              endif
!             * Setting SUM(ruptkvec)=etm__ and distributing according to relative uptake:
              ruptkvec(:) = etm__ * (ruptkvec(:) / (SUM(ruptkvec(:))))
            else
              ruptkvec(:) = 0.d0
              changef     = 1.d0
              etm__       = 0.d0
              transp_     = 0.d0
              gstom__     = 0.d0
            endif

          else
            ruptkvec(:) = 0.d0
          endif

        endif

        postempg = MIN(posg, wlayer_)
        if (MAXVAL(etmg__(:,:)) .gt. 0.d0) then
!         * root uptake by grasses can not be negative, as storage negligible
          ruptkg(1:posg) = MAX(0.d0,((-pcapvec(1:postempg)             &
     &                   + (prootmg - phydrostaticvec(1:postempg)))    &
     &                   * rsurfg_(:)) / (rrootm + (SQRT(p_pi / 2.d0)  &
     &                   * SQRT(rootrad * delzvec(1:postempg)          &
     &                   / rsurfg_(:))) / kunsatvec(1:postempg)))
          ruptkg(postempg+1:maxlayer) = 0.d0
          if (SUM(ruptkg(:)) .gt. 0.d0) then
            where (etmg__(:,:) .gt. SUM(ruptkg(:)))
              rootlim(:,:)  = 1.d0
              etmg__(:,:)   = SUM(ruptkg(:))
              transpg(:,:)  = etmg__(:,:) * 55555.555555555555d0  ! (Out[249]) mol/s=m/s*10^6 g/m/(18g/mol)
              gstomg__(:,:) = transpg(:,:) / (p_a * vd__)
            end where
            ruptkg(1:postempg) = etmg__(2,2) * (ruptkg(1:postempg)     &
     &                         / (SUM(ruptkg(:))))
          else
            ruptkg(:)      = 0.d0
            etmg__(:,:)    = 0.d0
            transpg(:,:)   = 0.d0
            gstomg__(:,:)  = 0.d0
          endif
        else
          ruptkg(:) = 0.d0
        endif
      else
        ruptkg(:)     = 0.d0
        ruptkvec(:)   = 0.d0
        etmg__(:,:)   = 0.d0
        transpg(:,:)  = 0.d0
        gstomg__(:,:) = 0.d0
      endif

      return
      end subroutine vom_rootuptake

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----steady-state tissue water (mqss) --------------------------------

      subroutine vom_mqss (mqss_out)
      use vom_vegwat_mod
      implicit none

      REAL*8, INTENT(out) :: mqss_out

      REAL*8 :: sum1, sum2, mul1, mul2

!     * (Out[257]) steady-state Mq

!      mqss_out = MAX(0.9d0 * mqx_,(mqx_ * (p_mpbar * (md_ * md_ + 752.d0 &
!     &         * md_ * mqx_ + mqx_ * mqx_) * SUM((rsurfvec(1:postemp_) &
!     &         / (rrootm + rsoilvec(1:postemp_)))) - (md_ + mqx_) * (md_ &
!     &         + mqx_) * (etm__ - SUM(((-phydrostaticvec(1:postemp_)   &
!     &         - pcapvec(1:postemp_)) * rsurfvec(1:postemp_))          &
!     &         / (rrootm + rsoilvec(1:postemp_)))))) / (p_mpbar * (md_ &
!     &         * md_ + 752.d0 * md_ * mqx_ + mqx_ * mqx_)              &
!     &         * SUM((rsurfvec(1:postemp_) / (rrootm + rsoilvec(1:postemp_))))))

      sum1 = Sum(rsurfvec(1:postemp_) / (rrootm + rsoilvec(1:postemp_)))
      mul1 = p_mpbar * (md_ * md_ + 752.d0 * md_ * mqx_ + mqx_ * mqx_) * sum1

      sum2 = Sum(((-phydrostaticvec(1:postemp_) - pcapvec(1:postemp_)) &
     &     * rsurfvec(1:postemp_)) / (rrootm + rsoilvec(1:postemp_)))
      mul2 = (md_ + mqx_) * (md_ + mqx_) * (etm__ - sum2)

      mqss_out = mqx_ * (mul1 - mul2) / mul1
      mqss_out = Max(0.9d0 * mqx_, mqss_out)

      return
      end subroutine vom_mqss

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----transpiration, gstom and tissue water ---------------------------

      subroutine vom_tissue_water_et (netass)
      use vom_vegwat_mod
      implicit none

      REAL*8,  INTENT(inout) :: netass
      character(len=135) :: msg

!     * makes sure that tissue water does not get below 0.9mqx
      if (mq_ .le. 0.9d0 * mqx_) then
        if (wlayer_ .ge. 1) then
          if (etm__ .gt. 0.9d0 * SUM(ruptkvec(:))) then
            if (SUM(ruptkvec(:)) .ge. 0.d0) then
              etm__ = SUM(ruptkvec(:))
              transp_ = etm__ * 55555.555555555555d0  ! (Out[249]) mol/s=m/s*10^6 g/m/(18g/mol)
              gstom__ = transp_ / (p_a * vd__)
            else
              write(msg,'(a20,i2,a1,i2,a1,i4)') 'vegetation dies on: ', &
     &          fday(nday), '/', fmonth(nday), '/', fyear(nday)
              write(*,*) TRIM(msg)
              netass = 0.d0
!             * if tissues water depleted, but still loosing water -> death
              finish = 1
              return
            endif
            call vom_mqss(mqss_)
            mqssmin = MIN(mqssmin, mqss_)
          endif
        else
          etm__   = 0.d0
          transp_ = 0.d0
          gstom__ = 0.d0
        endif
      endif
      if (wlayer_ .ge. 0) then
!       * (3.35), 1.e6 to convert from m (=1000kg/m2) to g/m2; (Out[250])
        dmq = (SUM(ruptkvec(:)) - etm__) * 1.d6
      else
        dmq = -etm__ * 1.d6
      endif

      return
      end subroutine vom_tissue_water_et

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*----water balance and conditions at next time step--------------------

      subroutine vom_subhourly ()
      use vom_vegwat_mod
      implicit none

      REAL*8  :: dtss
      REAL*8  :: dtmq         ! Maximum timestep allowed by tree water content change

      dtmq = 99999.d0

      if (md_ .gt. 0.d0) then
!       * avoids mq from becoming larger than mqx or smaller than 0.9mqx
        if (dmq .gt. 0.d0) then
          dtmq = (mqx_ - mq_) / dmq
        elseif (dmq .lt. 0.d0) then
          dtmq = (0.9d0 * mqx_ - mq_) / dmq
        endif

        if (ABS(mq_ - mqss_) .gt. mqx_ / 1.d6) then
          dtss = (mq_ - mqss_) / (1.d6 * (etm__ - SUM(ruptkvec(:))))
          if (dtss .le. 0.d0) dtss = 99999.d0
        else
          dtss = 99999.d0
        endif
      else
        dtss = 99999.d0
      endif

      dtmax = MIN(dtss, dtmq, 3600.d0 - time)

!     * waterbalance uses dtmax for the determination of dt
      call waterbalance()

      return
      end subroutine vom_subhourly

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine vom_add_hourly ()
      use vom_vegwat_mod
      implicit none

      REAL*8 :: ass__(3)
      REAL*8 :: assg__(3,3)

      ass__(:) = (4.d0 * ca__ * gstom__ + 8.d0 * gammastar_ * gstom__  &
     &         + jact_(:) - 4.d0 * rl__(:) - SQRT((-4.d0 * ca__        &
     &         * gstom__ + 8.d0 * gammastar_ * gstom__ + jact_(:)      &
     &         - 4.d0 * rl__(:)) ** 2.d0 + 16.d0 * gammastar_          &
     &         * gstom__ * (8.d0 * ca__ * gstom__ + jact_(:) + 8.d0    &
     &         * rl__(:)))) / 8.d0  ! (3.22) ; (Out[319])
      ass_h(:) = ass_h(:) + ass__(:) * dt_
      assg__(:,:) = (4.d0 * ca__ * gstomg__(:,:) + 8.d0 * gammastar_   &
     &            * gstomg__(:,:) + jactg(:,:) - 4.d0 * rlg__(:,:)     &
     &            - SQRT((-4.d0 * ca__ * gstomg__(:,:) + 8.d0          &
     &            * gammastar_ * gstomg__(:,:) + jactg(:,:) - 4.d0     &
     &            * rlg__(:,:)) ** 2.d0 + 16.d0 * gammastar_           &
     &            * gstomg__(:,:) * (8.d0 * ca__ * gstomg__(:,:)       &
     &            + jactg(:,:) + 8.d0 * rlg__(:,:)))) / 8.d0  ! (3.22); (Out[319])
      assg_h(:,:) = assg_h(:,:) + assg__(:,:) * dt_
      ruptkvec_h(:) = ruptkvec_h(:) + ruptkvec(:) * dt_
      ruptkg_h(:)   = ruptkg_h(:)   + ruptkg(:)   * dt_
      if (optmode .eq. 0) then
        spgfcf_h = spgfcf_h + dt_ * spgfcf__
        infx_h   = infx_h   + dt_ * infx__
        io_h     = io_h     + dt_ * io__
        esoil_h  = esoil_h  + dt_ * esoil__
        etm_h    = etm_h    + dt_ * etm__
        etmg_h   = etmg_h   + dt_ * etmg__(2,2)
        ruptk_h  = ruptk_h  + dt_ * SUM(ruptkvec(:))
      endif

      return
      end subroutine vom_add_hourly

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine vom_add_daily ()
      use vom_vegwat_mod
      implicit none

      vd_d     = vd_d     + vd__
      etm_d    = etm_d    + etm_h
      etmg_d   = etmg_d   + etmg_h
      esoil_d  = esoil_d  + esoil_h
      spgfcf_d = spgfcf_d + spgfcf_h
      infx_d   = infx_d   + infx_h
      rl_d     = rl_d     + rl__(2) * 3600.d0  ! rl_d in mol/day
      rlg_d    = rlg_d    + rlg__(2,2) * 3600.d0

      return
      end subroutine vom_add_daily

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine vom_write_hourly ()
      use vom_vegwat_mod
      implicit none

      CHARACTER(60) :: hourlyformat
      CHARACTER(3)  :: str

      if (fyear(nday) .ge. firstyear .and. fyear(nday) .le. lastyear) then
!       * internal write to convert from number to string
        write(str,'(i3)') wlayer_
!       * includes a column for each sublayer
        hourlyformat = '(i6,i6,i4,i7,i5,'//str//'e14.6)'
        write(kfile_resultshourly,'(i6,i7,i7,i7,i7,22e15.5)')          &
     &    fyear(nday), fmonth(nday), fday(nday), nday, nhour, rain__,  &
     &    tair__, par__, vd__, esoil_h, pc_ + pcg_(2), jmax25_(2),     &
     &    jmax25g(2), mq_, rl__(2) + rlg__(2,2), lambda_, lambdag_,    &
     &    rr_ + rrg, ass_h(2), assg_h(2,2), etm_h, etmg_h, suvec_(1),  &
     &    ys_, wsnew, spgfcf_h, infx_h
        write(kfile_delyuhourly,hourlyformat) fyear(nday),             &
     &    fmonth(nday), fday(nday), nday, nhour, delzvec(1:wlayer_)
        write(kfile_ruptkhourly,hourlyformat) fyear(nday),             &
     &    fmonth(nday), fday(nday), nday, nhour, ruptkvec_h(1:wlayer_)
        write(kfile_suvechourly,hourlyformat) fyear(nday),             &
     &    fmonth(nday), fday(nday), nday, nhour, suvec_(1:wlayer_)
      endif

      return
      end subroutine vom_write_hourly

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*----- check water balance --------------------------------------------

      subroutine vom_check_water ()
      use vom_vegwat_mod
      implicit none

      REAL*8  :: error1
      character(len=135) :: msg

      ioacum = ioacum + io_h
      wsnew  = SUM(cH2Ol_s(:))
      error  = wsold + ioacum - wsnew

!     * gives an error message if accumulated error exceeds 1 mm

      if (abs(error) .gt. 1.d-3) then
        write(msg,*) 'Error in water balance [mm]:', error, 'io=', io__, &
     &    'wsold=', wsold, 'wsnew=', wsnew
        write(*,*) TRIM(msg)
        finish = 1
      elseif (md_ .gt. 0.d0) then
        error1 = mqold + (ruptk_h - etm_h) * 1.d6 - mqnew
        if (abs(error1 / mqnew) .gt. 1.d-6) then
          write(msg,*) 'Error in tree water balance [%]:', error1 * 100.d0, &
     &      'mqold=', mqold, 'mqnew=', mqnew, 'hruptk=', ruptk_h, 'hetm=', etm_h
          write(*,*) TRIM(msg)
          finish = 1
        endif
      endif

      return
      end subroutine vom_check_water

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine vom_write_dayyear ()
      use vom_vegwat_mod
      implicit none

      CHARACTER(60) :: dailyformat
      CHARACTER(3)  :: str

!     * internal write to convert from number to string
      write(str,'(i3)') wlayer_
!     * includes a column for each sublayer
      dailyformat = '(i6,i6,i4,i7,'//str//'e14.6)'

      write(kfile_resultsdaily,'(i6,i7,i7,i7,i7,25e15.5)')             &
     &  fyear(nday), fmonth(nday), fday(nday), nday, nhour,            &
     &  rainvec(nday), tairmax(nday), tairmin(nday), parvec(nday),     &
     &  vd_d / 24.d0, esoil_d, jmax25_(2), jmax25g(2), pc_ + pcg_(2),  &
     &  rl_d + rlg_d, lambda_, lambdag_, rr_ * 3600.d0 * 24.d0,        &
     &  rrg * 3600.d0 * 24.d0, ass_d(2), assg_d(2,2),                  &
     &  SUM(suvec_(1:wlayer_)) / wlayer_, ys_, wsnew, spgfcf_d,        &
     &  infx_d, etm_d, etmg_d, suvec_(1), topt_
      write(kfile_rsurfdaily,dailyformat) fyear(nday), fmonth(nday),   &
     &  fday(nday), nday, rsurfvec(1:wlayer_)

      if (fyear(nday) .ne. nyear) then
!       * for calculation of vd_y a -1 is added to nday for using dayyear of correct year
        write(kfile_yearly,'(i6,18e16.6)') nyear, rain_y,              &
     &    par_y, srad_y, vd_y / (dayyear(nday-1)), esoil_y, etm_y,     &
     &    etmg_y, assg_y, rlg_y, rrg_y, cpccg_y, tcg_y,                &
     &    etmt_y, asst_y, rlt_y, rrt_y, cpcct_y, tct_y
      endif

!     * WRITING THE ACCUMULATED DATA FROM THE LAST YEAR TO FILE:

      if (nday .eq. maxday) then
!       * call subroutine there to get yearly data for the output
        call vom_add_yearly()
        write(kfile_yearly,'(i6,18e16.6)') nyear, rain_y,              &
     &    par_y, srad_y, vd_y / (dayyear(nday)), esoil_y, etm_y,       &
     &    etmg_y, assg_y, rlg_y, rrg_y, cpccg_y, tcg_y,                &
     &    etmt_y, asst_y, rlt_y, rrt_y, cpcct_y, tct_y
      endif

      return
      end subroutine vom_write_dayyear

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine vom_add_yearly ()
      use vom_vegwat_mod
      implicit none

      if (fyear(nday) .eq. nyear) then
        rain_y   = rain_y + rainvec(nday)    ! in [mm]
        par_y    = par_y + parvec(nday)
        srad_y   = srad_y + srad__(nday)     ! srad originally in MJ/day
        vd_y     = vd_y + vd_d / 24.d0
        etm_y    = etm_y + (etm_d + etmg_d) * 1000.d0  ! in[mm]
        esoil_y  = esoil_y + esoil_d * 1000.d0  ! in [mm]
!       * for grasses
        etmg_y   = etmg_y + etmg_d * 1000.d0 ! in [mm]
        assg_y   = assg_y + assg_d(2,2)
        rlg_y    = rlg_y + rlg_d
        rrg_y    = rrg_y + rrg * 3600.d0 * 24.d0
        cpccg_y  = cpccg_y + cpccg(2) * 3600.d0 * 24.d0
        tcg_y    = tcg_y + tcg(2) * 3600.d0 * 24.d0
!       * for trees
        etmt_y   = etmt_y + etm_d * 1000.d0   ! in [mm]
        asst_y   = asst_y + ass_d(2)
        rlt_y    = rlt_y + rl_d
        rrt_y    = rrt_y + rr_ * 3600.d0 * 24.d0
        cpcct_y  = cpcct_y + cpcc_ * 3600.d0 * 24.d0
        tct_y    = tct_y + tc_ * 3600.d0 * 24.d0
      else
        nyear    = fyear(nday)
        rain_y   = rainvec(nday)
        par_y    = parvec(nday)
        srad_y   = srad__(nday)              ! srad originally in MJ/day
        vd_y     = vd_d / 24.d0
        etm_y    = (etmg_d + etm_d) * 1000.d0
        esoil_y  = esoil_d * 1000.d0
!       * for grasses
        etmg_y   = etmg_d * 1000.d0
        assg_y   = assg_d(2,2)
        rlg_y    = rlg_d
        rrg_y    = rrg * 3600.d0 * 24.d0
        cpccg_y  = cpccg(2) * 3600.d0 * 24.d0
        tcg_y    = tcg(2) * 3600.d0 * 24.d0
!       * for trees
        etmt_y = etm_d * 1000.d0
        asst_y   = ass_d(2)
        rlt_y    = rl_d
        rrt_y    = rr_ * 3600.d0 * 24.d0
        cpcct_y  = cpcc_ * 3600.d0 * 24.d0
        tct_y    = tc_ * 3600.d0 * 24.d0
      endif

      return
      end subroutine vom_add_yearly

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*------ADJUSTMENT OF JMAX25 and PC-------------------------------------

      subroutine vom_adapt_foliage ()
      use vom_vegwat_mod
      implicit none

      REAL*8  :: netassg_d(3,3)         ! Daily grass net carbon profit
      INTEGER :: pos1(1)                ! Pointer to variable values that achieved maximum assimilation

      pos1(:)        = MAXLOC(ass_d(:))
      jmax25_(2)     = jmax25_(pos1(1))
      ass_d(:)       = 0.d0
      netassg_d(1,:) = assg_d(1,:) - 3600.d0 * 24.d0 * (cpccg(1) + rrg + tcg(1))
      netassg_d(2,:) = assg_d(2,:) - 3600.d0 * 24.d0 * (cpccg(2) + rrg + tcg(2))
      netassg_d(3,:) = assg_d(3,:) - 3600.d0 * 24.d0 * (cpccg(3) + rrg + tcg(3))
      pos2(:)        = MAXLOC(netassg_d(:,:))
      pcg_(2)        = MIN(1.d0 - pc_, pcg_(pos2(1)))
      jmax25g(2)     = jmax25g(pos2(2))
      assg_d(:,:)    = 0.d0

      return
      end subroutine vom_adapt_foliage

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*------ADJUSTMENT OF ROOT SURFACE--------------------------------------

      subroutine vom_adapt_roots ()
      use vom_vegwat_mod
      implicit none

      REAL*8 :: maxval_tmp

!     *-----PERENNIAL VEGETATION---------------

      reff(:) = 0.d0
      if (md_ .gt. 0.d0) then  ! if md_=0, then changef is calculated elsewhere
        changef = (0.95d0 * mqx_ - mqssmin) / (0.05d0 * mqx_)  ! (3.47)
      else
!       * changef for md_=0 is either 0 or 1. Change it to be either -1 or + 1:
        changef = 2.d0 * changef - 1.d0
      endif
      reff(1:pos_) = 0.5d0 * ruptkvec_d(1:pos_) / rsurfvec(1:pos_)     &
     &             / (MAXVAL(ruptkvec_d(1:pos_) / rsurfvec(1:pos_)))  ! (3.48)
      where (ruptkvec_d(1:pos_) .lt. 0.d0)
        reff(:) = 0.d0
      end where
      if (changef .lt. 0.d0) then
        reff(:) = 1.d0 - reff(:)
      endif

!     * rsurf=(2*epsln_/rootrad) if all pores filled by roots

      rsurfnewvec(1:pos_) = MIN(2.d0 * epsln_ / rootrad * delzvec(1:pos_), &
     &                    MAX(rsurfmin * delzvec(1:pos_), rsurfvec(1:pos_) &
     &                    + rsurfvec(1:pos_) * growthmax * changef     &
     &                    * reff(1:pos_) * delzvec(1:pos_)))

!     *-----SEASONAL VEGETATION---------------

!     * rootlim is either 0 or 1. Change it to be either -1 or + 1:
      rootlim(pos2(1),pos2(2)) = 2.d0 * rootlim(pos2(1),pos2(2)) - 1.d0

      reffg(:) = 0.d0
      maxval_tmp = MAXVAL(ruptkg_d(1:posg) / rsurfg_(1:posg))
      if (maxval_tmp .ne. 0.d0) then
        reffg(1:posg) = 0.5d0 * ruptkg_d(1:posg) / rsurfg_(1:posg) / maxval_tmp  ! (3.48)
      endif

!     * if roots are going to be reduced, reverse effectivity vector

      if (rootlim(pos2(1),pos2(2)) .lt. 0.d0) then
        reffg(:) = 1.d0 - reffg(:)
      endif

!     * maximum rsurfg depends on rsurf of trees in same layer.

      rsurfgnew(1:posg) = MIN(2.d0 * epsln_ / rootrad                  &
     &                  * delzvec(1:posg) - rsurfvec(1:posg),          &
     &                  MAX(rsurfmin * delzvec(1:posg),                &
     &                  rsurfg_(1:posg) + rsurfg_(1:posg) * growthmax  &
     &                  * rootlim(pos2(1),pos2(2)) * reffg(1:posg)))
      rsurfgnew(posg+1:maxlayer) = 0.d0

      rootlim(:,:)  = 0.d0
      ruptkvec_d(:) = 0.d0
      changef       = 0.d0

      return
      end subroutine vom_adapt_roots
