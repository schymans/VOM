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
!*    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!*
!***********************************************************************


      subroutine transpmodel(invar, dim_invar, nrun, netass, option1)
      use vegwatbal
      implicit none

      INTEGER, INTENT(in)    :: dim_invar
      INTEGER, INTENT(in)    :: nrun
      REAL*8,  INTENT(inout) :: netass
      INTEGER, INTENT(in)    :: option1
      REAL*8, DIMENSION(dim_invar), INTENT(in) :: invar

      INTEGER :: finish

      finish = 0

!dd      if (dim_invar .lt. 8) then
      if (dim_invar .lt. 6) then
        write(*,*) "ERROR: Number of input parameters less than 6."
        stop
      endif

      netass = 0.d0

      if (option1 .eq. 2) then
        optmode = 0
      else
        optmode = 1
      endif
	  
      if (option1 .eq. 3) optmode = 2

!*----------------------------------------------------------------------
!*     Optimised parameters reading from invar
!*----------------------------------------------------------------------

      lambdagfac = invar(1)
      wsgexp     = invar(2)
      lambdafac  = invar(3)
      wsexp      = invar(4)
      pc_        = invar(5)
      rootdepth  = invar(6)
!dd      mdstore    = invar(7)
!dd      rgdepth    = invar(8)

      if (parsaved .ne. 1) then

!*-----Reading input parameter------------------------------------------

        call vom_read_input()

!*-----allocate vector sizes--------------------------------------------

        call vom_alloc()

!*-----File opening (saving climate and gstom ass data)-----------------

        if (optmode .eq. 0) call vom_open_output()
		
        if (optmode .eq. 2) call vom_open_ncp_output()

!*-----PARAMETER READING FROM SOILPROFILE.PAR---------------------------

        call vom_get_soilprofile()

!*-----Climate and Calendar data reading--------------------------------

        call vom_get_hourly_clim()

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

!*-----DAILY LOOPS------------------------------------------------------

      d___ = 0
      dtest = nyt * 365
      if (optmode .eq. 0 .or. optmode .eq. 2) dtest = N__
      do while (d___  .lt. dtest)
        d___ = d___ + 1

        call vom_daily_init()

!*-----HOURLY LOOPS-----------------------------------------------------

!       * loops through each hour of daily dataset
        do h__ = 1, 24

          call vom_hourly_init()

!*-----calculate gstom, et and ass--------------------------------------

          call vom_gstom()

!*-----SUB-HOURLY LOOPS-------------------------------------------------

          do while (time .lt. 3600.d0)

!*-----setting variables from previous loop-----------------------------

            call vom_subhourly_init()

!*-----root water uptake------------------------------------------------

            call vom_rootuptake()

            if (md_ .gt. 0.d0) then

!*-----steady-state tissue water (mqss)---------------------------------

              if (nlayers_ .ge. 1) then
                call vom_mqss(mqss_)
              else
                mqss_ = 0.9d0 * mqx_
              endif
              mqssmin = MIN(mqssmin,mqss_)
            
!*-----transpiration, gstom and tissue water----------------------------

              call vom_tissue_water_et(finish, netass)
              if (finish .eq. 1) return
            endif

!*-----water balance and conditions at next time step-------------------

            call vom_subhourly()

            time = time + dt_
            mqnew = mq_ + dmq * dt_

!*-----adding up hourly fluxes------------------------------------------

            call vom_add_hourly()

!*-----END OF HOUR------------------------------------------------------

          enddo

!         * rl does not need to be included here as ass=-rl if j=0 (at night)
          netass = netass + hass_(2) - 3600.d0 * (cpcc_ + rr_ + tc_)   &
     &           + hassg(2,2) - 3600.d0 * (cpccg(2) + rrg + tcg(2))
          ass_d(:) = ass_d(:) + hass_(:)
          assg_d(:,:) = assg_d(:,:) + hassg(:,:)
          ruptkvec_d(:) = ruptkvec_d(:) + hruptkvec(:)
          ruptkg_d(:) = ruptkg_d(:) + hruptkg(:)

          if (optmode .eq. 0) then

            call vom_add_daily()
            call vom_write_hourly()

!*-----check water balance----------------------------------------------

            call vom_check_water(finish)
            if (finish .eq. 1) return

          endif
        enddo

!*-----END OF DAY-------------------------------------------------------

        if (optmode .eq. 0) then
          call vom_write_dayyear()
		  call vom_add_yearly()
		endif

!*-----ADJUSTMENT OF JMAX25 and PC--------------------------------------

        call vom_adapt_foliage()

!*-----ADJUSTMENT OF ROOT SURFACE---------------------------------------

        call vom_adapt_roots()

        if ((d___ .eq. dtest) .and. (d___ .lt. N__)) then
          if (netass .le. 0.d0) then
!           * estimates how bad the carbon loss would be instead of
!             running through the whole set
            netass = netass / nyt * ny_
          else
            dtest = N__
          endif
        endif

      enddo

!*-----END OF DAILY LOOPS-----------------------------------------------

      if (optmode .eq. 0) then
        print *,"Cumulative error in water balance (initial Ws+Input-Output-final Ws, in m): ",error
        print *,"Number of times dtsu was limiting: ",dtsu_count
        print *,"Number of times dtmax was limiting: ",dtmax_count
        close(kfile_resultshourly)
        close(kfile_resultsdaily)
        close(kfile_yearly)
        close(kfile_rsurfdaily)
        close(kfile_delyuhourly)
        close(kfile_ruptkhourly)
        close(kfile_suvechourly)
      endif
	  
      if (optmode .eq. 2) then
        call vom_write_model_output(netass)
        close(kfile_model_output)
      endif

      return
      end subroutine transpmodel

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----PARAMETER READING FROM INPUT.PAR---------------------------------

      subroutine vom_read_input ()
      use vegwatbal
      implicit none

      INTEGER :: stat

      open(kfile_inputpar, status='old',                               &
     &     file=sfile_inputpar(1:len_trim(sfile_inputpar)))
      read(kfile_inputpar,*) oa
      read(kfile_inputpar,*) alpha
      read(kfile_inputpar,*) cpccf
      read(kfile_inputpar,*) tcf
      read(kfile_inputpar,*) ny_
      read(kfile_inputpar,*) nyt
      read(kfile_inputpar,*) k25
      read(kfile_inputpar,*) kopt
      read(kfile_inputpar,*) ha_
      read(kfile_inputpar,*) hd
      read(kfile_inputpar,*) toptfac
      read(kfile_inputpar,*) tstart
      read(kfile_inputpar,*) rlratio

!*-----Catchment parameters --------------------------------------------

      read(kfile_inputpar,*) lat
      read(kfile_inputpar,*) cz
      read(kfile_inputpar,*) zs_
      read(kfile_inputpar,*) cgs
      read(kfile_inputpar,*) zr_
      read(kfile_inputpar,*) go_    

!*-----Soil parameters -------------------------------------------------

      read(kfile_inputpar,*) ksat_
      read(kfile_inputpar,*) thetar_
      read(kfile_inputpar,*) thetas_
      read(kfile_inputpar,*) nvg_
      read(kfile_inputpar,*) avg_
      epsln_ = thetas_ - thetar_             ! epsilon, porosity see Reggiani (2000)
      mvg_ = 1.d0 - (1.d0 / nvg_)            ! van Genuchten soil parameter m

!*-----Vertical Resolution ---------------------------------------------

      read(kfile_inputpar,*) delz_

!*-----Vegetation Parameters -------------------------------------------

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
      close(kfile_inputpar)

!     * The file soilprofile.par contain information about thickness and
!       soil properties in each soil layer, with the layer number in the
!       first column.

      M___ = 0

      open(kfile_soilprofile, status='old', iostat=stat,               &
     &     file=sfile_soilprofile(1:len_trim(sfile_soilprofile)))
      if (stat .eq. 0) then
        read(kfile_soilprofile,*) M___
      endif
      close(kfile_soilprofile)

!     * number of soil layers M___ assuming same thickness everywhere
      if (M___ .eq. 0) M___ = ceiling(cz / delz_)
      N__ = ceiling(ny_ * 365.d0)
      Nh_ = N__ * 24

      return
      end subroutine vom_read_input

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----allocate vector sizes--------------------------------------------

      subroutine vom_alloc ()
      use vegwatbal
      implicit none

      allocate(srad(N__))
      allocate(tmin(N__))
      allocate(tmax(N__))
      allocate(rainvec(N__))
      allocate(press(N__))
      allocate(year(N__))
      allocate(month(N__))
      allocate(day(N__))
      allocate(dayyear(N__))
      allocate(vpvec(N__))
      allocate(epan(N__))
      allocate(cavec(N__))
      allocate(parvec(N__))
      allocate(netassvec(N__))
      allocate(netassvecg(N__))
      allocate(parh(Nh_))
      allocate(vdh(Nh_))
      allocate(tairh(Nh_))
      allocate(gammastarvec(Nh_))
      allocate(rainh(Nh_))
      allocate(cah(Nh_))
      allocate(pcapvec(M___))
      allocate(suvec_(M___))
      allocate(ruptkvec(M___))
      allocate(sunewvec(M___))
      allocate(kunsatvec(M___))
      allocate(delzvec(M___))
      allocate(rsurfvec(M___))
      allocate(rsurfnewvec(M___))
      allocate(qblvec(M___))
      allocate(dsuvec(M___))
      allocate(phydrostaticvec(M___))
      allocate(delznewvec(M___))
      allocate(wsnewvec(M___))
      allocate(prootmvec(M___))
      allocate(pcapnewvec(M___))
      allocate(ruptkvec_d(M___))
      allocate(hruptkvec(M___))
      allocate(hruptkg(M___))
      allocate(ruptkg_d(M___))
      allocate(reff(M___))
      allocate(reffg(M___))
      allocate(ruptkg(M___))
      allocate(rsurfg_(M___))
      allocate(rsurfgnew(M___))
      allocate(rsoilvec(M___))
      allocate(kunsatnewvec(M___))

      allocate(ksatvec(M___))
      allocate(epslnvec(M___))
      allocate(thetasvec(M___))
      allocate(thetarvec(M___))
      allocate(sueqvec(M___))
      allocate(cH2Ol_s(M___))
      allocate(iovec(M___))
      allocate(avgvec(M___))
      allocate(nvgvec(M___))
      allocate(mvgvec(M___))

      return
      end subroutine vom_alloc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----File opening (saving climate and gstom ass data)-----------------

      subroutine vom_open_output ()
      use vegwatbal
      implicit none

      open(kfile_resultshourly, status='replace',                      &
     &     file=sfile_resultshourly(1:len_trim(sfile_resultshourly)))
      write(kfile_resultshourly,'(a6,a7,a7,a7,a7,24a15)') 'year',      &
     &  'month', 'day', 'dcum', 'hour', 'rain', 'tair', 'par', 'vd',   &
     &  'esoil', 'pc', 'jmax25_t', 'jmax25_g', 'mq', 'rl', 'lambda_t', &
     &  'lambda_g', 'rr', 'ass_t', 'ass_g', 'het_t', 'het_g', 'su_1',  &
     &  'ys', 'Ws', 'omgo', 'spgfcf', 'infx'

      open(kfile_resultsdaily, status='replace',                       &
     &     file=sfile_resultsdaily(1:len_trim(sfile_resultsdaily)))
      write(kfile_resultsdaily,'(a6,a7,a7,a7,a7,25a15)') 'year',       &
     &  'month', 'day', 'dcum', 'hour', 'rain', 'tmax', 'tmin', 'par', &
     &  'vd', 'esoil', 'jmax25_t', 'jmax25_g', 'pc', 'rlt+rlg',        &
     &  'lambda_t', 'lambda_g', 'rr_t', 'rr_g', 'ass_t', 'ass_g',      &
     &  'su_avg', 'ys', 'ws', 'spgfcf', 'infx', 'etm_t', 'etm_g',      &
     &  'su_1', 'topt'

      open(kfile_yearly, status='replace',                             &
     &     file=sfile_yearly(1:len_trim(sfile_yearly)))
      write(kfile_yearly,'(a6,19a16)') "year", "rainyr", "epanyr",     &
     &  "paryr", "radyr", "vdyr", "esoilyr", "etyr", "etmgyr",         &
	 &  "assgyr", "rlgyr", "rrgyr", "cpccgyr", "tcgyr",                &
	 &  "etmtyr", "asstyr", "rltyr", "rrtyr", "cpcctyr", "tcyr"

      open(kfile_rsurfdaily, status='replace',                         &
     &     file=sfile_rsurfdaily(1:len_trim(sfile_rsurfdaily)))
      write(kfile_rsurfdaily,*) ' year', ' month', ' day', '   dcum',  &
     &  '  rsurfsublayer'

      open(kfile_delyuhourly, status='replace',                         &
     &     file=sfile_delyuhourly(1:len_trim(sfile_delyuhourly)))
      write(kfile_delyuhourly,*) ' year', ' month', ' day', '   dcum',  &
     &  ' hour', '  delyusublayer'

      open(kfile_ruptkhourly, status='replace',                        &
     &     file=sfile_ruptkhourly(1:len_trim(sfile_ruptkhourly)))
      write(kfile_ruptkhourly,*) ' year', ' month', ' day', '   dcum', &
     &  ' hour', '  delyusublayer'

      open(kfile_suvechourly, status='replace',                        &
     &     file=sfile_suvechourly(1:len_trim(sfile_suvechourly)))
      write(kfile_suvechourly,*) ' year', ' month', ' day', '   dcum', &
     &  ' hour', '  susublayer'

      return
      end subroutine vom_open_output

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----File opening (saving ncp)----------------------------------------

      subroutine vom_open_ncp_output ()
      use vegwatbal
      implicit none

      open(kfile_model_output, file=sfile_model_output(1:len_trim(sfile_model_output)))

      return
      end subroutine vom_open_ncp_output

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----PARAMETER READING FROM SOILPROFILE.PAR---------------------------

      subroutine vom_get_soilprofile ()
      use vegwatbal
      implicit none

      INTEGER :: stat, j

!     * The file soilprofile.par can contain information about thickness
!       and soil properties in each soil layer, with the layer number in
!       the first column.

      open(kfile_soilprofile, status='old', iostat=stat,               &
     &     file=sfile_soilprofile(1:len_trim(sfile_soilprofile)))
      if (stat .eq. 0) then
        do j = 1, M___
          read(kfile_soilprofile,*) M___, delzvec(j), ksatvec(j),      &
     &      nvgvec(j), avgvec(j), thetasvec(j), thetarvec(j)
          epslnvec(j) = thetasvec(j) - thetarvec(j)  ! porosity
          mvgvec(j) = 1.d0 - (1.d0 / nvgvec(j))  ! van Genuchten soil parameter m
        enddo
      else
        delzvec(:)   = delz_
        ksatvec(:)   = ksat_
        nvgvec(:)    = nvg_
        avgvec(:)    = avg_
        thetasvec(:) = thetas_
        thetarvec(:) = thetar_
        epslnvec(:)  = epsln_
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
      use vegwatbal
      implicit none

      INTEGER :: ii, i, h, oldh, stat
      INTEGER :: dummyint1, dummyint2, dummyint3, dummyint4

      open(kfile_hourlyweather, status='old', iostat=stat,             &
     &     file=sfile_hourlyweather(1:len_trim(sfile_hourlyweather)))
      if (stat .ne. 0) then
        close(kfile_hourlyweather)

!       * Creating hourly climate data from daily data

        open(kfile_dailyweather, status='old', iostat=stat,            &
     &       file=sfile_dailyweather(1:len_trim(sfile_dailyweather)))
        read(kfile_dailyweather,*)
        do i = 1, N__
          read(kfile_dailyweather,'(4i8,8f8.2)') dayyear(i), day(i),   &
     &      month(i), year(i), tmax(i), tmin(i), rainvec(i), epan(i),  &
     &      srad(i), vpvec(i), press(i), cavec(i)
        enddo
        close(kfile_dailyweather)

        open(kfile_hourlyweather, iostat=stat,                         &
     &       file=sfile_hourlyweather(1:len_trim(sfile_hourlyweather)))
        write(kfile_hourlyweather,'(5a8,5a11)') '   hour', ' dayyear', &
     &    '     day', '   month', '    year', '       tair',            &
     &    '         vd', '       parh', '      rainh', '      cah'

!       * Calculation of derived parameters

        call vom_calc_derived()

        close(kfile_hourlyweather)

      endif
     

!       * Reading hourly climate data if available

			open(kfile_hourlyweather, status='old', iostat=stat,             &
     &     file=sfile_hourlyweather(1:len_trim(sfile_hourlyweather)))
			read(kfile_hourlyweather,*)
			ii = 1
			oldh = 99
			do i = 1, Nh_
				read(kfile_hourlyweather,'(5i8,5e11.3)') h, dummyint1,       &
	 &      dummyint2, dummyint3, dummyint4, tairh(i), vdh(i),         &
	 &      parh(i), rainh(i), cah(i)
				if (h .lt. oldh) then
					dayyear(ii) = dummyint1
					day(ii)     = dummyint2
					month(ii)   = dummyint3
					year(ii)    = dummyint4
					ii = ii + 1
				endif
				oldh = h
				gammastarvec(i) = 0.00004275d0 * E__ ** ((18915.d0 * (-25.d0 &
	 &                    + tairh(i))) / (149.d0 * R__ * (273.d0       &
	 &                    + tairh(i))))      ! (Out[274], derived from (3.25))
			enddo
			close(kfile_hourlyweather)
      

      return
      end subroutine vom_get_hourly_clim

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----Calculation of derived parameters--------------------------------

      subroutine vom_calc_derived ()
      use vegwatbal
      implicit none

      INTEGER :: in, ik, ii
      REAL*8  :: sunr, suns
      REAL*8  :: tamean
      REAL*8  :: dtr

      parvec(:) = 2.0804d0 * srad(:)         ! (Out[17]), par in mol/m2 if srad was MJ/m2
      do in = 1, N__
        daylength = 12.d0 - 7.639437d0 * ASIN((0.397949d0              &
     &            * COS(0.172142d0 + 0.017214d0 * dayyear(in))         &
     &            * TAN(0.017453d0 * lat))                             &
     &            / SQRT(0.920818d0 - 0.079182d0                       &
     &            * COS(0.344284d0 + 0.034428d0 * dayyear(in))))  ! (Out[22]), in hours
!       * sets time of sunrise and sunset
        sunr = 12d0 - 0.5d0 * daylength
        suns = 12d0 + 0.5d0 * daylength 
        tamean = (tmax(in) + tmin(in)) / 2.d0
        dtr = tmax(in) - tmin(in)
        vp_ = vpvec(in) * 100.d0             ! vp in Pa

!       * Loop through every hour of day, where ik=hour
        do ik = 1, 24
          ii = in * 24 + ik - 24
!         * (derived from 3.52+3.53) (Out[38], accounts for diurnal
!           variation in air temperature
          tair = tamean + dtr * (0.0138d0                              &
     &         * COS(3.513d0 - ((-1.d0 + ik) * Pi) / 3.d0) + 0.0168d0  &
     &         * COS(0.822d0 - ((-1.d0 + ik) * Pi) / 4.d0) + 0.0984d0  &
     &         * COS(0.360d0 - ((-1.d0 + ik) * Pi) / 6.d0) + 0.4632d0  &
     &         * COS(3.805d0 - ((-1.d0 + ik) * Pi) / 12.d0))
          tairh(ii) = tair
          cah(ii) = cavec(in)
!         vd__ = 0.006028127d0 * 2.718282d0 ** ((17.27d0 * tair)       &
!    &       / (237.3d0 + tair)) - 9.869233d-6 * vp_  ! (derived from 3.54+3.55) (Out[52]), accounts for diurnal variation in vapour deficit

          vd__ = (((0.6108d0 * E__ ** (17.27d0 * tair / (tair + 237.3d0))) &
     &         * 1000) - vp_) / (press(in) * 100.d0)

          if (vd__ .le. 0.d0) vd__ = 0.d0
          vdh(ii) = vd__
!         * average rainfall in hour ii (m/s)
          rainh(ii) = rainvec(in) / (24.d0 * 3600.d0 * 1000.d0)
          if (sunr .le. ik .and. ik + 1 .le. suns) then
!           * ([Out30]), in mol/m2/s (derived from 3.51) accounts for
!             diurnal variation in global irradiance
            parh(ii) = (-0.000873d0 * parvec(in) * COS(0.017453d0 * lat) &
     &               * SQRT(0.920818d0 - 0.079182d0                    &
     &               * COS(0.034428d0 * (10.d0 + dayyear(in))))        &
     &               * COS(0.2618d0 * ik) - 0.000347d0 * parvec(in)    &
     &               * COS(0.017214d0 * (10.d0 + dayyear(in)))         &
     &               * SIN(0.017453d0 * lat)) / (-1.250192d0 * daylength &
     &               * COS(0.017214d0 * (10.d0 + dayyear(in)))         &
     &               * SIN(0.017453d0 * lat) + 24.d0                   &
     &               * COS(0.017453d0 * lat)                           &
     &               * SQRT(0.920818d0 - 0.079182d0                    &
     &               * COS(0.034428d0 * (10.d0 + dayyear(in))))        &
     &               * (1.d0 - (0.158363d0                             &
     &               * COS(0.017214d0 * (10.d0 + dayyear(in))) ** 2.d0 &
     &               * TAN(0.017453d0 * lat) ** 2.d0)                  &
     &               / (0.920818d0 - 0.079182d0 * COS(0.034428d0       &
     &               * (10.d0 + dayyear(in))))) ** 0.5d0)
          else 
            parh(ii) = 0.d0
          endif

          write(kfile_hourlyweather,'(5i8,5e11.3)') ik, dayyear(in),   &
     &      day(in), month(in), year(in), tairh(ii), vdh(ii),          &
     &      parh(ii), rainh(ii), cah(ii)

!         * (Out[274], derived from (3.25))
          gammastarvec(ii) = 0.00004275d0 * E__ ** ((18915.d0          &
     &                     * (-25.d0 + tairh(ii))) / (149.d0 * R__     &
     &                     * (273.d0 + tairh(ii))))
        enddo
      enddo

      return
      end subroutine vom_calc_derived

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----Initial values---------------------------------------------------

      subroutine vom_init_vegpar ()
        use vegwatbal
        implicit none

        INTEGER :: init
        REAL*8  :: dummy

        topt_ = tstart

        !     * Set soil moisture and vegetation parameters to initial conditions 

        !     * Set init=1 to signalise to waterbalance that this is the inital step
        init = 1
        call waterbalance(init)

        if (rootdepth .gt. cz) then
           write(*,*) 'Root depth greater than soil depth'
           rootdepth = cz
        endif

        wsold  = SUM(cH2Ol_s(:))               ! initial soil water storage
        wsnew_ = wsold

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
        rsurfgnew(posg+1:M___) = 0.d0
        if (posg .gt. nlayersnew) then
           rsurfgnew(nlayersnew+1:posg) = rsurfmin                        &
                &                               * delzvec(nlayersnew+1:pos_)      &
                &                               * omgunew
        endif
        !     * root surface density (root surface area/soil volume) in each sublayer
        rsurfnewvec(1:pos_) = rsurfinit * delzvec(1:pos_)
        if (pos_ .gt. nlayersnew) then
           rsurfnewvec(nlayersnew+1:pos_) = rsurfmin                      &
                &                                 * delzvec(nlayersnew+1:pos_)    &
                &                                 * omgunew
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

        yr_            = year(1)
        !d      netassyr       = 0.d0
        !d      gppyr          = 0.d0
        rainyr         = 0.d0
        paryr          = 0.d0
        radyr          = 0.d0
        vdyr           = 0.d0
        etyr           = 0.d0
        epanyr         = 0.d0
        evapyr         = 0.d0        ! = yearly esoil
        ruptkvec_d(:)  = 0.d0
        ruptkg_d(:)    = 0.d0
        ass_d(:)       = 0.d0
        assg_d(:,:)    = 0.d0
        netassvec(:)   = 0.d0
        netassvecg(:)  = 0.d0
        iocum          = 0.d0
        ! for grasses
        etmg_y         = 0.d0
        assg_y         = 0.d0
        rlg_y          = 0.d0
        rrg_y          = 0.d0
        cpccg_y        = 0.d0
        tcg_y          = 0.d0
        ! for trees
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
      use vegwatbal
      implicit none

      rsurfvec(:) = rsurfnewvec(:)
      rsurfg_(:)  = rsurfgnew(:)
      lambda_     = lambdafac * (SUM(pcapnewvec(1:pos_)) / pos_) ** wsexp  ! (3.45) 
      lambdag     = lambdagfac * pcapnewvec(1) ** wsgexp  ! (3.44)
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

      if (pos_ .gt. nlayersnew) then
!       * (3.42, 2.45e-10 from (Out[165])) costs of water distribution and storage
        cpcc_ = cpccf * pc_ * rootdepth + mdstore * 2.45d-10
      else
        cpcc_ = cpccf * pc_ * SUM(delzvec(1:pos_)) + mdstore * 2.45d-10
      endif

      if (nlayersnew .lt. posg) then
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
        tmax(d___)    = -9999.d0
        tmin(d___)    =  9999.d0
        rainvec(d___) =     0.d0
        parvec(d___)  =     0.d0
        srad(d___)    =     0.d0
      endif

      if (optmode .eq. 0) then
        ruptk_d  = 0.d0
        vd_d     = 0.d0
        jmax_d   = 0.d0
        jmaxg_d  = 0.d0
        gstom_d  = 0.d0
        gstomg_d = 0.d0
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
      use vegwatbal
      implicit none

      INTEGER :: ii

      ii        = d___ * 24 + h__ - 24
      rain_     = rainh(ii)
      gammastar = gammastarvec(ii)
      tair      = tairh(ii)
      vd__      = vdh(ii)
      par_      = parh(ii)
      ca_       = cah(ii) / 1.0d6

!     * (Out[310], derived from (3.26)) Temperature dependence of Jmax
      jmax__(:) = (E__ ** ((ha_ * (-25.d0 + tair) * (-273.d0 + topt_   &
     &          + 273.d0 * R__ * topt_)) / ((25.d0 + 273.d0 * R__ * topt_) &
     &          * (tair + 273.d0 * R__ * topt_))) * ((-1.d0 + E__ ** (-(hd &
     &          * (-298.d0 + topt_)) / (25.d0 + 273.d0 * R__ * topt_)))  &
     &          * ha_ + hd) * jmax25_(:)) / ((-1.d0 + E__ ** ((hd * (273.d0 &
     &          + tair - topt_)) / (tair + 273.d0 * R__ * topt_))) * ha_ + hd)
      rl__(:) = ((ca_ - gammastar) * pc_ * jmax__(:) * rlratio) / (4.d0 &
     &        * (ca_ + 2.d0 * gammastar) * (1.d0 + rlratio))  ! (3.24), (Out[312])
!     * (Out[310], derived from (3.26)) Temperature dependence of Jmax
      jmaxg__(:) = (E__ ** ((ha_ * (-25.d0 + tair) * (-273.d0 + topt_  &
     &           + 273.d0 * R__ * topt_)) / ((25.d0 + 273.d0 * R__ * topt_) &
     &           * (tair + 273.d0 * R__ * topt_))) * ((-1.d0 + E__ ** (-(hd &
     &           * (-298.d0 + topt_)) / (25.d0 + 273.d0 * R__ * topt_))) &
     &           * ha_ + hd) * jmax25g(:)) / ((-1.d0 + E__ ** ((hd * (273.d0 &
     &           + tair - topt_)) / (tair + 273.d0 * R__ * topt_))) * ha_ + hd)
      rlg__(1,:) = ((ca_ - gammastar) * pcg_(1) * jmaxg__(:) * rlratio) &
     &           / (4.d0 * (ca_ + 2.d0 * gammastar) * (1.d0 + rlratio))  ! (3.24), (Out[312])
      rlg__(2,:) = ((ca_ - gammastar) * pcg_(2) * jmaxg__(:) * rlratio) &
     &           / (4.d0 * (ca_ + 2.d0 * gammastar) * (1.d0 + rlratio))  ! (3.24), (Out[312])
      rlg__(3,:) = ((ca_ - gammastar) * pcg_(3) * jmaxg__(:) * rlratio) &
     &           / (4.d0 * (ca_ + 2.d0 * gammastar) * (1.d0 + rlratio))  ! (3.24), (Out[312])

!     * daily recalculation for resultsdaily
      if (optmode .eq. 0) then
        rainvec(d___) = rainvec(d___) + rain_ * 3600.d0 * 1000.d0  ! mm/d
        parvec(d___)  = parvec(d___) + par_ * 3600.d0  ! in mol/m2/d
        srad(d___)    = parvec(d___) / 2.0804d0  ! MJ/m2/d

        if (tair .gt. tmax(d___)) tmax(d___) = tair
        if (tair .lt. tmin(d___)) tmin(d___) = tair
      endif

      if (optmode .eq. 0) then
        mqold   = mqnew
        hspgfcf = 0.d0
        hinf_   = 0.d0
        hinfx   = 0.d0
        hio     = 0.d0
        hesoil  = 0.d0
        hetm_   = 0.d0
        hetmg   = 0.d0
        hruptk_ = 0.d0
      endif

      time         = 0.d0
      hass_(:)     = 0.d0                    ! hourly assimilation
      hassg(:,:)   = 0.d0
      hruptkvec(:) = 0.d0
      hruptkg(:)   = 0.d0

      return
      end subroutine vom_hourly_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----calculate gstom, et and ass -------------------------------------

      subroutine vom_gstom ()
      use vegwatbal
      implicit none

      REAL*8 :: cond1, cond2
      REAL*8 :: cond3(3,3)
      REAL*8 :: part1, part2, part3, part4, part5
      REAL*8 :: part6, part7, part8, part9

      if (par_ .gt. 0.d0) then
!       * adaptation of topt to air temperature during sunlight
        topt_ = topt_ + toptfac * (tair + 273.d0 - topt_)
        jact_(:)   = (1.d0 - E__ ** (-(alpha * par_) / jmax__(:)))     &
     &             * jmax__(:) * pc_         ! (3.23), (Out[311])
        jactg(1,:) = (1.d0 - E__ ** (-(alpha * par_) / jmaxg__(:)))    &
     &             * jmaxg__(:) * pcg_(1)    ! (3.23), (Out[311])
        jactg(2,:) = (1.d0 - E__ ** (-(alpha * par_) / jmaxg__(:)))    &
     &             * jmaxg__(:) * pcg_(2)    ! (3.23), (Out[311])
        jactg(3,:) = (1.d0 - E__ ** (-(alpha * par_) / jmaxg__(:)))    &
     &             * jmaxg__(:) * pcg_(3)    ! (3.23), (Out[311])

        cond1      = (2.d0 * a__ * vd__) / (ca_ + 2.d0 * gammastar)
        cond2      = (4.d0 * ca_ * rl__(2) + 8.d0 * gammastar * rl__(2)) &
     &             / (ca_ - gammastar)
        cond3(:,:) = (4.d0 * ca_ * rlg__(:,:) + 8.d0 * gammastar * rlg__(:,:)) &
     &             / (ca_ - gammastar)

        if (vd__ .gt. 0.d0 .and. lambda_ .gt. cond1 .and. jact_(2) .gt. cond2) then

!          gstom__ = MAX(0.d0, (0.25d0 * (a__ * (ca_ * (jact_(2) - 4.d0 &
!     &            * rl__(2)) - 4.d0 * gammastar * (jact_(2) + 2.d0     &
!     &            * rl__(2))) * vd__ * (ca_ * lambda_ + 2.d0           &
!     &            * gammastar * lambda_ - a__ * vd__)                  &
!     &            + 1.7320508075688772d0 * SQRT(a__ * gammastar        &
!     &            * jact_(2) * (ca_ * (jact_(2) - 4.d0 * rl__(2))      &
!     &            - gammastar * (jact_(2) + 8.d0 * rl__(2))) * vd__    &
!     &            * (ca_ * lambda_ + 2.d0 * gammastar * lambda_        &
!     &            - 2.d0 * a__ * vd__) ** 2.d0 * (ca_ * lambda_        &
!     &            + 2.d0 * gammastar * lambda_ - a__ * vd__))))        &
!     &            / (a__ * (ca_ + 2.d0 * gammastar) ** 2.d0 * vd__     &
!     &            * (ca_ * lambda_ + 2.d0 * gammastar * lambda_ - a__ * vd__)))

          part1 = ca_ + 2.d0 * gammastar
          part2 = part1 * lambda_ - a__ * vd__
          part3 = a__ * vd__ * part2

          part4 = ca_ * (jact_(2) - 4.d0 * rl__(2))
          part5 = gammastar * jact_(2)
          part6 = gammastar * 8.d0 * rl__(2)
          part7 = part4 - part5 - part6

          part8 = Sqrt(part5 * part7 * (part2 - a__ * vd__) ** 2.d0 * part3)
          part9 = part7 - 3.d0 * part5 + 1.7320508075688772d0 * part8 / part3

          gstom__ = 0.25d0 * part9 / part1**2.d0
          gstom__ = Max(0.d0, gstom__)       ! (Out[314])

        else
          gstom__ = 0.d0
        endif
        transp_ = a__ * vd__ * gstom__       ! (3.28) transpiration rate in mol/s
        etm__ = (transp_ * 18.d0) / (10.d0 ** 6.d0)  ! transpiration rate in m/s

        where (vd__ .gt. 0.d0 .and. lambdag .gt. cond1 .and. jactg(:,:) .gt. cond3(:,:))
          gstomg__(:,:) = MAX(0.d0,(0.25d0 * (a__ * (ca_ * (jactg(:,:) &
     &                  - 4.d0 * rlg__(:,:)) - 4.d0 * gammastar        &
     &                  * (jactg(:,:) + 2.d0 * rlg__(:,:))) * vd__     &
     &                  * (ca_ * lambdag + 2.d0 * gammastar * lambdag  &
     &                  - a__ * vd__) + 1.7320508075688772d0           &
     &                  * SQRT(a__ * gammastar * jactg(:,:) * (ca_     &
     &                  * (jactg(:,:) - 4.d0 * rlg__(:,:)) - gammastar &
     &                  * (jactg(:,:) + 8.d0 * rlg__(:,:))) * vd__     &
     &                  * (ca_ * lambdag + 2.d0 * gammastar * lambdag  &
     &                  - 2.d0 * a__ * vd__) ** 2.d0 * (ca_ * lambdag  &
     &                  + 2.d0 * gammastar * lambdag - a__ * vd__))))  &
     &                  / (a__ * (ca_ + 2.d0 * gammastar) ** 2.d0      &
     &                  * vd__ * (ca_ * lambdag + 2.d0 * gammastar     &
     &                  * lambdag - a__ * vd__)))  ! (Out[314])
        elsewhere
          gstomg__(:,:) = 0.d0
        endwhere
        transpg(:,:) = a__ * vd__ * gstomg__(:,:)  ! (3.28) transpiration rate in mol/s
        etmg__(:,:) = (transpg(:,:) * 18.d0) / (10.d0 ** 6.d0)  ! transpiration rate in m/s
      else
        par_          = 0.d0
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
      use vegwatbal
      implicit none

!!$      if (nlayersnew .ge. 1) then
!!$        if (nlayers_ .ge. 1) then         ! Need to account for the change in rsurf due to a change in unsaturated soil volume.
!!$          rsurfvec(1:nlayersnew) = rsurfvec(1:nlayersnew)            &
!!$     &                           / (omgu_ * delzvec(1:nlayersnew))   &
!!$     &                           * omgunew * delznewvec(1:nlayersnew)
!!$          rsurfg_(1:nlayersnew) = rsurfg_(1:nlayersnew) / (omgu_     &
!!$     &                          * delzvec(1:nlayersnew)) * omgunew   &
!!$     &                          * delznewvec(1:nlayersnew)
!!$        else                              ! If nlayers at the previous time step was 0, rsurf was set to rsurfmin.
!!$          rsurfvec(1:nlayersnew) = rsurfvec(1:nlayersnew) * omgunew  &
!!$     &                           * delzvec(1:nlayersnew)
!!$          rsurfg_(1:nlayersnew) = rsurfg_(1:nlayersnew) * omgunew    &
!!$     &                          * delzvec(1:nlayersnew)
!!$        endif
!!$      endif

      if (nlayersnew .lt. pos_) then
        rsurfvec(nlayersnew+1:pos_) = rsurfmin * delz_
      endif
      if (nlayersnew .lt. posg) then
        rsurfg_(nlayersnew+1:posg) = rsurfmin * delz_
      endif

      mq_          = mqnew
      ys_          = ysnew
      suvec_(:)    = sunewvec(:)
      yu__         = yunew_
      nlayers_     = nlayersnew
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
      use vegwatbal
      implicit none

      INTEGER :: i

      if (nlayersnew .ge. 1) then
        postemp_ = MIN(pos_, nlayers_)
        do i = 1, postemp_
          phydrostaticvec(i) = (i - 0.5d0) * delz_  ! (Out[238]) hydrostatic head for (3.34)
        enddo
        if (md_.gt.0.d0) then
          prootmvec(1:postemp_) = (mpbar * (-mq_ + mqx_) * (750.d0       &
     &                        - (750.d0 * mqx_) / (md_ + mqx_) + (md_  &
     &                        + mqx_) / mqx_)) / (md_ + mqx_)          &
     &                        - phydrostaticvec(1:postemp_)  ! (Out[239])
        else
         prootmvec(1:postemp_) = prootmg      ! set tissue suction to the same as in grasses, if no storage capacity
        endif
!       * soil resistance, (Out[ 241] with svolume=omgu*delzvec(1:postemp_)); derived from (3.32)
        rsoilvec(1:postemp_) = SQRT(Pi / 2.d0) * SQRT((rootrad         &
     &                       * omgu_ * delzvec(1:postemp_))            &
     &                       / rsurfvec(1:postemp_))                   &
     &                       / kunsatvec(1:postemp_)

!       * root water uptake, Chapter 3.3.3.3 (Out[242])
        if (md_ .gt. 0.d0) then
          ruptkvec(1:postemp_) = (-pcapvec(1:postemp_)                   &
     &                         + prootmvec(1:postemp_))                  &
     &                         * rsurfvec(1:postemp_) / (rrootm          &
     &                         + rsoilvec(1:postemp_))
          ruptkvec(postemp_+1:M___) = 0.d0
        else      ! if no storage, uptake happens only when etm__>0
          if (etm__ .gt. 0.d0) then
            ruptkvec(1:postemp_) = MAX(0.d0,(-pcapvec(1:postemp_)        &
     &                         + prootmvec(1:postemp_))                  &
     &                         * rsurfvec(1:postemp_) / (rrootm          &
     &                         + rsoilvec(1:postemp_)))
            ruptkvec(postemp_+1:M___) = 0.d0
            
            if (SUM(ruptkvec(:)) .gt. 0.d0) then
              if (etm__ .gt. SUM(ruptkvec(:))) then
                changef = 1.d0
                etm__   = SUM(ruptkvec(:))
                transp_ = etm__ * 55555.555555555555d0  ! (Out[249]) mol/s=m/s*10^6 g/m/(18g/mol)
                gstom__ = transp_ / (a__ * vd__)
              endif
              ! Setting SUM(ruptkvec)=etm__ and distributing according to relative uptake:
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

        postempg = MIN(posg, nlayers_)
        if (MAXVAL(etmg__(:,:)) .gt. 0.d0) then
!         * root uptake by grasses can not be negative, as storage negligible
          ruptkg(1:posg) = MAX(0.d0,((-pcapvec(1:postempg)             &
     &                   + (prootmg - phydrostaticvec(1:postempg)))    &
     &                   * rsurfg_(:)) / (rrootm + (SQRT(Pi / 2.d0)    &
     &                   * SQRT(rootrad * omgu_ * delzvec(1:postempg)  &
     &                   / rsurfg_(:))) / kunsatvec(1:postempg)))
          ruptkg(postempg+1:M___) = 0.d0
          if (SUM(ruptkg(:)) .gt. 0.d0) then
            where (etmg__(:,:) .gt. SUM(ruptkg(:)))
              rootlim(:,:)  = 1.d0
              etmg__(:,:)   = SUM(ruptkg(:))
              transpg(:,:)  = etmg__(:,:) * 55555.555555555555d0  ! (Out[249]) mol/s=m/s*10^6 g/m/(18g/mol)
              gstomg__(:,:) = transpg(:,:) / (a__ * vd__)
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
      use vegwatbal
      implicit none

      REAL*8, INTENT(out) :: mqss_out

      REAL*8 :: sum1, sum2, mul1, mul2

!     * (Out[257]) steady-state Mq

!      mqss_out = MAX(0.9d0 * mqx_,(mqx_ * (mpbar * (md_ * md_ + 752.d0 &
!     &         * md_ * mqx_ + mqx_ * mqx_) * SUM((rsurfvec(1:postemp_) &
!     &         / (rrootm + rsoilvec(1:postemp_)))) - (md_ + mqx_) * (md_ &
!     &         + mqx_) * (etm__ - SUM(((-phydrostaticvec(1:postemp_)   &
!     &         - pcapvec(1:postemp_)) * rsurfvec(1:postemp_))          &
!     &         / (rrootm + rsoilvec(1:postemp_)))))) / (mpbar * (md_   &
!     &         * md_ + 752.d0 * md_ * mqx_ + mqx_ * mqx_)              &
!     &         * SUM((rsurfvec(1:postemp_) / (rrootm + rsoilvec(1:postemp_))))))

      sum1 = Sum(rsurfvec(1:postemp_) / (rrootm + rsoilvec(1:postemp_)))
      mul1 = mpbar * (md_ * md_ + 752.d0 * md_ * mqx_ + mqx_ * mqx_) * sum1

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

      subroutine vom_tissue_water_et (finish, netass)
      use vegwatbal
      implicit none

      INTEGER, INTENT(inout) :: finish
      REAL*8,  INTENT(inout) :: netass

!     * makes sure that tissue water does not get below 0.9mqx
      if (mq_ .le. 0.9d0 * mqx_) then
        if (nlayers_ .ge. 1) then
          if (etm__ .gt. 0.9d0 * SUM(ruptkvec(:))) then
            if (SUM(ruptkvec(:)) .ge. 0.d0) then
              etm__ = SUM(ruptkvec(:))
              transp_ = etm__ * 55555.555555555555d0  ! (Out[249]) mol/s=m/s*10^6 g/m/(18g/mol)
              gstom__ = transp_ / (a__ * vd__)
            else
              write(*,'(a20,i2,a1,i2,a1,i4)') 'vegetation dies on: ',  &
     &          day(d___), '/', month(d___), '/', year(d___)
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
      if (nlayers_ .ge. 0) then
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
      use vegwatbal
      implicit none

      INTEGER :: init
      REAL*8  :: dtss

      init = 0

      dtmq = 99999.d0
      if (md_ .gt. 0.d0) then
!     * avoids mq from becoming larger than mqx or smaller than 0.9mqx
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
      call waterbalance(init)

      return
      end subroutine vom_subhourly

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine vom_add_hourly ()
      use vegwatbal
      implicit none

      REAL*8 :: ass__(3)
      REAL*8 :: assg__(3,3)

      ass__(:) = (4.d0 * ca_ * gstom__ + 8.d0 * gammastar * gstom__    &
     &         + jact_(:) - 4.d0 * rl__(:) - SQRT((-4.d0 * ca_ * gstom__ &
     &         + 8.d0 * gammastar * gstom__ + jact_(:) - 4.d0          &
     &         * rl__(:)) ** 2.d0 + 16.d0 * gammastar * gstom__        &
     &         * (8.d0 * ca_ * gstom__ + jact_(:) + 8.d0 * rl__(:)))) / 8.d0  ! (3.22) ; (Out[319])
      hass_(:) = hass_(:) + ass__(:) * dt_
      assg__(:,:) = (4.d0 * ca_ * gstomg__(:,:) + 8.d0 * gammastar     &
     &            * gstomg__(:,:) + jactg(:,:) - 4.d0 * rlg__(:,:)     &
     &            - SQRT((-4.d0 * ca_ * gstomg__(:,:) + 8.d0 * gammastar &
     &            * gstomg__(:,:) + jactg(:,:) - 4.d0 * rlg__(:,:))    &
     &            ** 2.d0 + 16.d0 * gammastar * gstomg__(:,:) * (8.d0  &
     &            * ca_ * gstomg__(:,:) + jactg(:,:) + 8.d0 * rlg__(:,:)))) / 8.d0  ! (3.22); (Out[319])
      hassg(:,:) = hassg(:,:) + assg__(:,:) * dt_
      hruptkvec(:) = hruptkvec(:) + ruptkvec(:) * dt_
      hruptkg(:) = hruptkg(:) + ruptkg(:) * dt_
      if (optmode .eq. 0) then
        hspgfcf = hspgfcf + dt_ * spgfcf__
        hinfx   = hinfx   + dt_ * infx__
        hio     = hio     + dt_ * io_
        hesoil  = hesoil  + dt_ * esoil__
        hetm_   = hetm_   + dt_ * etm__
        hetmg   = hetmg   + dt_ * etmg__(2,2) 
        hruptk_ = hruptk_ + dt_ * SUM(ruptkvec(:))
        hinf_   = hinf_   + dt_ * inf_
      endif

      return
      end subroutine vom_add_hourly

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine vom_add_daily ()
      use vegwatbal
      implicit none

      netassvec(d___)  = netassvec(d___) + hass_(2) - 3600.d0 * (cpcc_ + rr_ + tc_)
!     * rl does not need to be included here as ass=-rl if j=0 (at night)
      netassvecg(d___) = netassvecg(d___) + hassg(2,2) - 3600.d0       &
     &                 * (cpccg(2) + rrg + tcg(2))

      ruptk_d  = ruptk_d  + SUM(ruptkvec(:)) * 3600.d0
      vd_d     = vd_d     + vd__
      jmax_d   = jmax_d   + jmax__(2)
      jmaxg_d  = jmaxg_d  + jmaxg__(2)
      gstom_d  = gstom_d  + gstom__
      gstomg_d = gstomg_d + gstomg__(2,2)
      etm_d    = etm_d    + hetm_
      etmg_d   = etmg_d   + hetmg
      esoil_d  = esoil_d  + hesoil
      spgfcf_d = spgfcf_d + hspgfcf
      infx_d   = infx_d   + hinfx
      rl_d     = rl_d     + rl__(2) * 3600.d0  ! rl_d in mol/day
      rlg_d    = rlg_d    + rlg__(2,2) * 3600.d0

      return
      end subroutine vom_add_daily

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine vom_write_hourly ()
      use vegwatbal
      implicit none

      CHARACTER(60) :: hourlyformat
      CHARACTER(3)  :: str

      if (year(d___) .ge. firstyear .and. year(d___) .le. lastyear) then
!       * internal write to convert from number to string
        write(str,'(i3)') nlayers_
!       * includes a column for each sublayer
        hourlyformat = '(i6,i6,i4,i7,i5,'//str//'e14.6)'
        write(kfile_resultshourly,'(i6,i7,i7,i7,i7,24e15.5)') year(d___), &
     &    month(d___), day(d___), d___, h__, rain_, tair, par_, vd__,  &
     &    hesoil, pc_ + pcg_(2), jmax25_(2), jmax25g(2), mq_,          &
     &    rl__(2) + rlg__(2,2), lambda_, lambdag, rr_ + rrg, hass_(2), &
     &    hassg(2,2), hetm_, hetmg, suvec_(1), ys_, wsnew_, omgo_,     &
     &    hspgfcf, hinfx
        write(kfile_delyuhourly,hourlyformat) year(d___), month(d___),  &
     &    day(d___), d___, h__, delzvec(1:nlayers_)
        write(kfile_ruptkhourly,hourlyformat) year(d___), month(d___), &
     &    day(d___), d___, h__, hruptkvec(1:nlayers_)
        write(kfile_suvechourly,hourlyformat) year(d___), month(d___), &
     &    day(d___), d___, h__, suvec_(1:nlayers_)
      endif

      return
      end subroutine vom_write_hourly

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*----- check water balance --------------------------------------------

      subroutine vom_check_water (finish)
      use vegwatbal
      implicit none

      INTEGER, INTENT(inout) :: finish

      REAL*8  :: error1

      iocum = iocum + hio
      wsnew_ = SUM(cH2Ol_s(:)) 
      error = wsold + iocum - wsnew_
!     * gives an error message if accumulated error exceeds 1 mm
      if (abs(error) .gt. 1.d-3) then
        write(*,*) 'Error in water balance [mm]:', error, 'io=', io_,  &
     &    'wsold=', wsold, 'wsnew=', wsnew_
        finish = 1
      elseif (md_ .gt. 0.d0) then
        error1 = mqold + (hruptk_ - hetm_) * 1.d6 - mqnew
        if (abs(error1 / mqnew) .gt. 1.d-6) then
          write(*,*) 'Error in tree water balance [%]:', error1 * 100.d0, &
     &      'mqold=', mqold, 'mqnew=', mqnew, 'hruptk=', hruptk_, 'hetm=', hetm_
          finish = 1
        endif
      endif

      return
      end subroutine vom_check_water

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine vom_write_dayyear ()
      use vegwatbal
      implicit none

      CHARACTER(60) :: dailyformat
      CHARACTER(3)  :: str

!     * internal write to convert from number to string
      write(str,'(i3)') nlayers_
!     * includes a column for each sublayer
      dailyformat = '(i6,i6,i4,i7,'//str//'e14.6)'

      write(kfile_resultsdaily,'(i6,i7,i7,i7,i7,25e15.5)') year(d___), &
     &  month(d___), day(d___), d___, h__, rainvec(d___), tmax(d___),  &
     &  tmin(d___), parvec(d___), vd_d / 24.d0, esoil_d, jmax25_(2),   &
     &  jmax25g(2), pc_ + pcg_(2), rl_d + rlg_d, lambda_, lambdag,     &
     &  rr_ * 3600.d0 * 24.d0, rrg * 3600.d0 * 24.d0, ass_d(2),        &
     &  assg_d(2,2), SUM(suvec_(1:nlayers_)) / nlayers_, ys_, wsnew_,  &
     &  spgfcf_d, infx_d, etm_d, etmg_d, suvec_(1), topt_
      write(kfile_rsurfdaily,dailyformat) year(d___), month(d___),     &
     &  day(d___), d___, rsurfvec(1:nlayers_)

      if (year(d___) .ne. yr_) then
!     * for calculation of vdyr a -1 is added to d___ for using dayyear of correct year
        write(kfile_yearly,'(i6,19e16.6)') yr_, rainyr, epanyr, paryr, &
     &  radyr, vdyr / (dayyear(d___ - 1)), evapyr, etyr, etmg_y,    &
	 &  assg_y, rlg_y, rrg_y, cpccg_y, tcg_y,                          &
	 &  etmt_y, asst_y, rlt_y, rrt_y, cpcct_y, tct_y
      endif

! WRITING THE ACCUMULATED DATA FROM THE LAST YEAR TO FILE:

      if (d___ .eq. N__) then
!	  * one more call of subroutine to write last day to yearly
		call vom_add_yearly()
	    write(kfile_yearly,'(i6,19e16.6)') yr_, rainyr, epanyr, paryr, &
     &  radyr, vdyr / (dayyear(d___)), evapyr, etyr, etmg_y,               &
	 &  assg_y, rlg_y, rrg_y, cpccg_y, tcg_y,                          &
	 &  etmt_y, asst_y, rlt_y, rrt_y, cpcct_y, tct_y
      endif

      return
      end subroutine vom_write_dayyear

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine vom_add_yearly ()
        use vegwatbal
        implicit none

        if (year(d___) .eq. yr_) then
           rainyr   = rainyr + rainvec(d___)    ! in [mm]
           epanyr   = epanyr + epan(d___)       ! epan originally in [mm]/day
           paryr    = paryr + parvec(d___)
           radyr    = radyr + srad(d___)        ! srad originally in MJ/day
           vdyr     = vdyr + vd_d / 24.d0
           etyr     = etyr + (etm_d + etmg_d) * 1000.d0  ! in[mm]
           evapyr   = evapyr + esoil_d * 1000.d0  ! in [mm]
           !d        netassyr = netassyr + ass_d(2) - (cpcc_ + rr_) * 3600.d0       &
           !d     &           * 24.d0 + assg_d(2,2) - (cpccg(2) + rrg) * 3600.d0 * 24.d0
           !d        gppyr    = gppyr + (ass_d(2) + rl_d) + assg_d(2,2) + rlg_d
           ! for grasses
           etmg_y   = etmg_y + etmg_d * 1000.d0 ! in [mm]
           assg_y   = assg_y + assg_d(2,2)
           rlg_y    = rlg_y + rlg_d
           rrg_y    = rrg_y + rrg * 3600.d0 * 24.d0
           cpccg_y  = cpccg_y + cpccg(2) * 3600.d0 * 24.d0
           tcg_y    = tcg_y + tcg(2) * 3600.d0 * 24.d0
           ! for trees
           etmt_y   = etmt_y + etm_d * 1000.d0   ! in [mm]
           asst_y   = asst_y + ass_d(2)
           rlt_y    = rlt_y + rl_d
           rrt_y    = rrt_y + rr_ * 3600.d0 * 24.d0
           cpcct_y  = cpcct_y + cpcc_ * 3600.d0 * 24.d0
           tct_y    = tct_y + tc_ * 3600.d0 * 24.d0
        else
           yr_      = year(d___)
           rainyr   = rainvec(d___)
           epanyr   = epan(d___)                ! epan originally in [mm]/day
           paryr    = parvec(d___)
           radyr    = srad(d___)                ! srad originally in MJ/day
           vdyr     = vd_d / 24.d0
           etyr     = (etmg_d + etm_d) * 1000.d0
           evapyr   = esoil_d * 1000.d0
           !d        netassyr = ass_d(2) - (cpcc_ + rr_) * 3600.d0 * 24.d0          &
           !d    &           + assg_d(2,2) - (cpccg(2) + rrg) * 3600.d0 * 24.d0
           !d        gppyr    = (ass_d(2) + rl_d) + assg_d(2,2) + rlg_d
           ! for grasses
           etmg_y   = etmg_d * 1000.d0
           assg_y   = assg_d(2,2)
           rlg_y    = rlg_d
           rrg_y    = rrg * 3600.d0 * 24.d0
           cpccg_y  = cpccg(2) * 3600.d0 * 24.d0
           tcg_y    = tcg(2) * 3600.d0 * 24.d0
           ! for trees
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

      subroutine vom_write_model_output (netass)
      use vegwatbal
      implicit none

      REAL*8, INTENT(in) :: netass

      write(kfile_model_output,'(e12.6)') netass

      return
      end subroutine vom_write_model_output

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*------ADJUSTMENT OF JMAX25 and PC-------------------------------------

      subroutine vom_adapt_foliage ()
      use vegwatbal
      implicit none

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
      use vegwatbal
      implicit none


!*-----PERENNIAL VEGETATION---------------
      reff(:) = 0.d0
      if (md_ .gt. 0.d0) then   !if md_=0, then changef is calculated elsewhere
        changef = (0.95d0 * mqx_ - mqssmin) / (0.05d0 * mqx_)  ! (3.47)
      else 
        !     * changef for md_=0 is either 0 or 1. Change it to be either -1 or + 1:
        changef = 2.d0*changef - 1.d0
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
      rsurfnewvec(1:pos_) = MIN(2.d0 * epsln_ / rootrad * delzvec(1:pos_),    &
     &                 MAX(rsurfmin * delzvec(1:pos_), rsurfvec(1:pos_)       &
     &                 + rsurfvec(1:pos_) * growthmax * changef               &
     &                 * reff(1:pos_) * delzvec(1:pos_)))
     
!*-----SEASONAL VEGETATION---------------
!     * rootlim is either 0 or 1. Change it to be either -1 or + 1:
      rootlim(pos2(1), pos2(2)) = 2.d0*rootlim(pos2(1), pos2(2)) - 1.d0
      reffg(:) = 0.d0
      reffg(1:posg) = 0.5d0 * ruptkg_d(1:posg) / rsurfg_(1:posg)     &
     &             / (MAXVAL(ruptkg_d(1:posg) / rsurfg_(1:posg)))  ! (3.48)
!     * if roots are going to be reduced, reverse effectivity vector     
      if (rootlim(pos2(1), pos2(2)).lt.0.d0) then
        reffg(:) = 1.d0 - reffg(:)
      endif
!     * maximum rsurfg depends on rsurf of trees in same layer.
      rsurfgnew(1:posg) = MIN(2.d0 * epsln_ / rootrad * delzvec(1:posg) &
     &                  - rsurfvec(1:posg), MAX(rsurfmin                &
     &                  * delzvec(1:posg), rsurfg_(1:posg)              &
     &                  + rsurfg_(1:posg) * growthmax * rootlim(pos2(1),&
     &                   pos2(2)) * reffg(1:posg)))
      rsurfgnew(posg+1:M___) = 0.d0
      rootlim(:,:)           = 0.d0
      ruptkvec_d(:)          = 0.d0
      changef = 0.d0

      return
      end subroutine vom_adapt_roots
