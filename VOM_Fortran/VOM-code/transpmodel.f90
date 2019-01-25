!***********************************************************************
!*  Transpiration model and layered water balance
!*----------------------------------------------------------------------
!*  Author: Stan Schymanski, CWR, University of Western Australia
!*  03/2006
!*  now at: Max Planck Institute for Biogeochemistry
!*  Email: sschym@bgc-jena.mpg.de
!*  02/2008
!*  Version: big leaf, trees and grass, layered unsaturated zone
!*  optimised root profile, pcg_d and Jmax25
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
  !$ USE omp_lib
      implicit none

      INTEGER, INTENT(in)    :: dim_invar
      REAL*8,  INTENT(inout) :: tp_netass
      INTEGER, INTENT(in)    :: option1
      REAL*8, DIMENSION(dim_invar), INTENT(in) :: invar
      !REAL*8, ALLOCATABLE, DIMENSION(:,:) :: output_mat

      tp_netass = 0.d0

      call transpmodel_init(invar, dim_invar, option1)

!     * DAILY LOOPS

      do while (nday  .lt. c_testday)
        nday = nday + 1

        call vom_daily_init()

!     * HOURLY LOOPS (loops through each hour of daily dataset)
      do nhour = 1, 24
        th_ = nday * 24 + nhour - 24

      call vom_hourly_init()

!     * calculate gstom, et and ass

      call vom_gstom()

!     * SUB-HOURLY LOOPS

      do while (time .lt. 3600.d0)

!       * setting variables from previous loop

        call vom_subhourly_init()

!       * root water uptake

        call vom_rootuptake()

        if (q_md .gt. 0.d0) then

!         * steady-state tissue water (mqss)

          if (wlayer_ .ge. 1) then
            call vom_mqss(mqsst_)
          else
            mqsst_ = 0.9d0 * q_mqx
          endif
          mqsstmin = MIN(mqsstmin,mqsst_)

!         * transpiration, gstom and tissue water

          call vom_tissue_water_et(tp_netass)
          if (finish .eq. 1) return

        endif

!       * water balance and conditions at next time step

        call vom_subhourly()

        time = time + dt_
        mqtnew = mqt_ + dmqt * dt_

!       * adding up hourly fluxes

        call vom_add_hourly()

!       * END OF HOUR

      enddo

!     * rl does not need to be included here as ass=-rl if j=0 (at night)
      tp_netass = tp_netass + asst_h(2) - 3600.d0 * (q_cpcct_d + rrt_d &
     &          + q_tct_d) + assg_h(2,2) - 3600.d0 * (cpccg_d(2)       &
     &          + rrg_d + tcg_d(2))

      asst_d(:)   = asst_d(:)   + asst_h(:)
      assg_d(:,:) = assg_d(:,:) + assg_h(:,:)
      ruptkt_d(:) = ruptkt_d(:) + ruptkt_h(:)
      ruptkg_d(:) = ruptkg_d(:) + ruptkg_h(:)

      !if (optmode .eq. 0) then

       !formatted output for single model run
       if (option1 .eq. 2) then
        call vom_add_daily()
        call vom_write_hourly()

!       * check water balance

        call vom_check_water()
        if (finish .eq. 1) return

      endif

       !formatted output for multiple runs
       if (option1 .eq. 5) then
        call vom_add_daily()
        !call vom_write_hourly() !replace with a new subroutine

!       * check water balance

        call vom_check_water()
        if (finish .eq. 1) return

      endif



        enddo

!       * END OF DAY


           !if(nday .eq. 1) then
              !allocate matrix for 21 variables with lengt of timeseries
           !   allocate( output_mat (21, c_testday ) )
           !end if

           !last daily step
           call transpmodel_daily_step(tp_netass, option1)

      enddo

!     * END OF DAILY LOOPS

      call transpmodel_last_step(tp_netass, option1)

      return
      end subroutine transpmodel

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++ Post-daily step

      subroutine transpmodel_daily_step (tp_netass, option1)
      use vom_vegwat_mod
      implicit none

      REAL*8,  INTENT(inout) :: tp_netass
      INTEGER, INTENT(in)    :: option1
      !REAL*8, DIMENSION(21, c_maxday ), INTENT(inout) :: output_mat

      !if (optmode .eq. 0) then
       !formatted output for single model run
       if (option1 .eq. 2) then
        call vom_write_dayyear( tp_netass )
        call vom_add_yearly()
      endif

       !formatted output for multiple runs
      !if (optmode .eq. 5) then
      if (option1 .eq. 5) then
        call vom_save_dayyear() !replace with new routine
        call vom_add_yearly()
      endif

!     * ADJUSTMENT OF JMAX25 and PC

      call vom_adapt_foliage()

!     * ADJUSTMENT OF ROOT SURFACE

      call vom_adapt_roots()

      if ((nday .eq. c_testday) .and. (nday .lt. c_maxday)) then
        if (tp_netass .le. 0.d0) then
!         * estimates how bad the carbon loss would be instead of
!           running through the whole set
          tp_netass = tp_netass / i_testyear * i_maxyear
        else
          c_testday = c_maxday
        endif
      endif

      return
      end subroutine transpmodel_daily_step

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++ Post-yearly step

      subroutine transpmodel_last_step (tp_netass, option)
      use vom_vegwat_mod
      implicit none

      REAL*8, INTENT(in) :: tp_netass
      integer, INTENT(in) :: option

      if (option .eq. 2) then
        print *,'Cumulative error in water balance (initial Ws+Input-Output-final Ws, in m): ',error
        print *,'Number of times dtsu was limiting: ',dtsu_count
        print *,'Number of times dtmax was limiting: ',dtmax_count

        close(kfile_resultshourly)
        close(kfile_resultsdaily)
        close(kfile_resultsyearly)
          close(kfile_rsurfdaily)
          close(kfile_delzhourly)
          close(kfile_ruptkthourly)
          close(kfile_suhourly)

        write(*,*) "Model run COMPLETE"
        write(*,*) " "
        write(*,*) "The carbon profit achieved is: ",tp_netass
        write(*,*) "Hourly results are saved in resulthourly.txt"
        write(*,*) "Daily results are saved in resultsdaily.txt"
        write(*,*) "Yearly results are saved in yearly.txt"
        write(*,*) "Soil results are saved in delyudaily.txt, rsurfdaily.txt, ruptkhourly.txt, suvechourly.txt"
      endif

      if (option .eq. 3) then
        open(kfile_model_output, FILE=trim(adjustl(i_outputpath))// &
             trim(adjustl(sfile_model_output)), STATUS='replace')
        write(kfile_model_output,'(E13.6)') tp_netass
        close(kfile_model_output)

        write(*,*) "Model run COMPLETE"
        write(*,*) " "
        write(*,*) "The carbon profit achieved is: ",tp_netass
        write(*,*) "Best ncp is saved in model_output.txt"
      endif

      return
      end subroutine transpmodel_last_step

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine transpmodel_init_once (vom_command)
      use vom_vegwat_mod
      implicit none

      INTEGER, INTENT(in) :: vom_command

!     * Reading input parameter

      call vom_read_input()

!     * allocate vector sizes

      call vom_alloc()

!     * File opening (saving climate and gstom ass data)

      if (vom_command .eq. 2) call vom_open_output()

!     * PARAMETER READING FROM SOILPROFILE.PAR

      call vom_get_soilprofile()

!     * Climate and Calendar data reading

      call vom_get_hourly_clim()

!     * get timeseries of vegetation cover
      if(i_read_pc == 1) then
         call vom_get_perc_cov()
      end if

      return
      end subroutine transpmodel_init_once

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine transpmodel_init (vom_invar, vom_npar, vom_command)
      use vom_vegwat_mod
  !$ USE omp_lib
      implicit none

      INTEGER, INTENT(in)    :: vom_npar
      INTEGER, INTENT(in)    :: vom_command
      REAL*8, DIMENSION(vom_npar), INTENT(in) :: vom_invar

      if( .not. allocated(pcap_) ) then

          allocate(pcap_(s_maxlayer))
          allocate(su__(s_maxlayer))
          allocate(ruptkt__(s_maxlayer))
          allocate(sunew(s_maxlayer))
          allocate(kunsat_(s_maxlayer))

          allocate(rsurft_(s_maxlayer))
          allocate(rsurftnew(s_maxlayer))
          allocate(qbl(s_maxlayer))
          allocate(dsu(s_maxlayer))
          allocate(prootm(s_maxlayer))

          allocate(pcapnew(s_maxlayer))
          allocate(ruptkt_d(s_maxlayer))
          allocate(ruptkt_h(s_maxlayer))
          allocate(ruptkg_h(s_maxlayer))
          allocate(ruptkg_d(s_maxlayer))
          allocate(refft(s_maxlayer))
          allocate(reffg(s_maxlayer))
          allocate(ruptkg__(s_maxlayer))
          allocate(rsurfg_(s_maxlayer))
          allocate(rsurfgnew(s_maxlayer))
          allocate(rsoil(s_maxlayer))
          allocate(kunsatnew(s_maxlayer))
          allocate(sueq(s_maxlayer))
          allocate(cH2Ol_s(s_maxlayer))
          allocate(iovec(s_maxlayer))

          allocate( output_mat (21, c_testday ) )

      end if

      if (vom_command .eq. 2) then
        optmode = 0
      elseif (vom_command .eq. 5) then
        optmode = 0
      elseif (vom_command .eq. 3) then
        optmode = 2
      else
        optmode = 1
      endif

!*----------------------------------------------------------------------
!*     Optimised parameters reading from vom_invar
!*----------------------------------------------------------------------


      if (vom_npar .ge. 1) o_lambdagf = vom_invar(1)
      if (vom_npar .ge. 2) o_wsgexp   = vom_invar(2)
      if (vom_npar .ge. 3) o_lambdatf = vom_invar(3)
      if (vom_npar .ge. 4) o_wstexp   = vom_invar(4)
      if (vom_npar .ge. 5) o_pct      = vom_invar(5)
      if (vom_npar .ge. 6) o_rtdepth  = vom_invar(6)
      if (vom_npar .ge. 7) o_mdstore  = vom_invar(7)
      if (vom_npar .ge. 8) o_rgdepth  = vom_invar(8)
      if (vom_npar .ge. 9)  i_cgs     = vom_invar(9)
      if (vom_npar .ge. 10) i_zr      = vom_invar(10)
      if (vom_npar .ge. 11) i_go      = vom_invar(11)
      if (vom_npar .ge. 12) i_ksat    = vom_invar(12)
      if (vom_npar .ge. 13) i_thetar  = vom_invar(13)
      if (vom_npar .ge. 14) i_thetas  = vom_invar(14)
      if (vom_npar .ge. 15) i_nvg     = vom_invar(15)
      if (vom_npar .ge. 16) i_avg     = vom_invar(16)


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
      c_testday = i_testyear * 365
      if (optmode .eq. 0 .or. optmode .eq. 2) c_testday = c_maxday
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

!     * Definition of variable parameters

      namelist /inputpar/ i_alpha, i_cpccf, i_tcf, i_maxyear,          &
     &                    i_testyear, i_ha, i_hd, i_toptf,             &
     &                    i_toptstart, i_rlratio, i_mdtf, i_mqxtf,     &
     &                    i_rrootm, i_rsurfmin, i_rsurf_, i_rootrad,   &
     &                    i_prootmg, i_growthmax, i_incrcovg,          &
     &                    i_incrjmax,                                  &
     &                    i_firstyear,i_lastyear, i_write_h, i_read_pc,&
     &                    i_inputpath, i_outputpath,                   &
     &                    o_lambdagf, o_wsgexp, o_lambdatf, o_wstexp,  &
     &                    o_pct, o_rtdepth, o_mdstore, o_rgdepth

      namelist /input2par/ i_lat, i_cz, i_cgs, i_zr, i_go, i_ksat,     &
     &                     i_thetar, i_thetas, i_nvg, i_avg, i_delz

!     * Input of variable parameters from the parameter file

      open(kfile_namelist, FILE=sfile_namelist, STATUS='old',          &
     &                     FORM='formatted', IOSTAT=iostat)
      if (iostat .eq. 0) then
        read(kfile_namelist, inputpar)
        read(kfile_namelist, input2par)
      endif
      close(kfile_namelist)

      c_epsln = i_thetas - i_thetar     ! epsilon, porosity see Reggiani (2000)
      i_mvg = 1.d0 - (1.d0 / i_nvg)     ! van Genuchten soil parameter m

      c_maxday = i_maxyear * 365
      c_maxhour = c_maxday * 24

!     * The file soilprofile.par contain information about thickness and
!       soil properties in each soil layer, with the layer number in the
!       first column.

      s_maxlayer = 0

      open(kfile_soilprofile, FILE=trim(adjustl(i_inputpath))// &
           trim(adjustl(sfile_soilprofile)),                  &
     &                        STATUS='old', IOSTAT=iostat)
      if (iostat .eq. 0) then
        read(kfile_soilprofile,*) s_maxlayer
      endif
      close(kfile_soilprofile)

!     * Check if i_cz is a multiple of i_delz
!       Raise a warning and correct if this is not the case
      if ( MOD(i_cz, i_delz) .gt. 0.0) then

        write(*,*) "WARNING: i_cz not a multiple of i_delz"
        write(*,*) "Setting i_cz", i_cz, "to", i_cz + MOD(i_cz, i_delz)  
        i_cz = i_cz + MOD(i_cz, i_delz) 

      end if

!     * Check if i_cz - i_zr is a multiple of i_delz
!       Raise a warning and correct if this is not the case
      if ( MOD(i_cz - i_zr, i_delz) .gt. 0.0) then

        write(*,*) "WARNING: i_cz-i_zr not a multiple of i_delz"
        write(*,*) "Correcting i_zr", i_zr, "to", i_zr + MOD(i_cz - i_zr, i_delz)  
        i_zr = i_zr + MOD(i_cz - i_zr, i_delz) 

      end if


!     * number of soil layers s_maxlayer assuming same thickness everywhere
      if (s_maxlayer .eq. 0) s_maxlayer = ceiling(i_cz / i_delz)

      return
      end subroutine vom_read_input

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----allocate vector sizes--------------------------------------------

      subroutine vom_alloc ()
      use vom_vegwat_mod
      implicit none

      allocate(dayyear(c_maxday))
      allocate(fday(c_maxday))
      allocate(fmonth(c_maxday))
      allocate(fyear(c_maxday))

      allocate(tairmax_d(c_maxday))
      allocate(tairmin_d(c_maxday))
      allocate(rain_d(c_maxday))
      allocate(srad_d(c_maxday))
      allocate(vp_d(c_maxday))
      allocate(press_d(c_maxday))
      allocate(ca_d(c_maxday))

      allocate(par_h(c_maxhour))
      allocate(vd_h(c_maxhour))
      allocate(tair_h(c_maxhour))
      allocate(rain_h(c_maxhour))
      allocate(ca_h(c_maxhour))

      allocate(par_d(c_maxday))

      allocate(s_delz(s_maxlayer))
      allocate(s_ksat(s_maxlayer))
      allocate(s_avg(s_maxlayer))
      allocate(s_nvg(s_maxlayer))
      allocate(s_thetas(s_maxlayer))
      allocate(s_thetar(s_maxlayer))
      allocate(c_mvg(s_maxlayer))
      allocate(c_hhydrst(s_maxlayer))

      return
      end subroutine vom_alloc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine vom_dealloc ()
      use vom_vegwat_mod
      implicit none

         if(allocated(pcap_)) then

              deallocate(pcap_)
              deallocate(su__)
              deallocate(ruptkt__)
              deallocate(sunew)
              deallocate(kunsat_)

              deallocate(rsurft_)
              deallocate(rsurftnew)
              deallocate(qbl)
              deallocate(dsu)
              deallocate(prootm)

              deallocate(pcapnew)
              deallocate(ruptkt_d)
              deallocate(ruptkt_h)
              deallocate(ruptkg_h)
              deallocate(ruptkg_d)
              deallocate(refft)
              deallocate(reffg)
              deallocate(ruptkg__)
              deallocate(rsurfg_)
              deallocate(rsurfgnew)
              deallocate(rsoil)
              deallocate(kunsatnew)
              deallocate(sueq)
              deallocate(cH2Ol_s)
              deallocate(iovec)

              deallocate(output_mat)

         end if

      return
      end subroutine vom_dealloc



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----File opening (saving climate and gstom ass data)-----------------

      subroutine vom_open_output ()
      use vom_vegwat_mod
      implicit none


      open(kfile_resultshourly, FILE=trim(adjustl(i_outputpath))//      &
           trim(adjustl(sfile_resultshourly)), STATUS='replace')
      write(kfile_resultshourly,'(A6,A7,A7,A7,A7,22A15)') 'fyear',     &
     &  'fmonth', 'fday', 'nday', 'nhour', 'rain', 'tair', 'par', 'vd',&
     &  'esoil', 'pc', 'jmax25t', 'jmax25g', 'mqt', 'rl', 'lambdat',   &
     &  'lambdag', 'rr', 'asst', 'assg', 'etmt', 'etmg', 'su_1',       &
     &  'zw', 'ws', 'spgfcf', 'infx'

      open(kfile_resultsdaily, FILE=trim(adjustl(i_outputpath))// &
           trim(adjustl(sfile_resultsdaily)), STATUS='replace')
      write(kfile_resultsdaily,'(A6,A7,A7,A7,A7,26A15)') 'fyear',      &
     &  'fmonth', 'fday', 'nday', 'nhour', 'rain', 'tairmax', 'tairmin', &
     &  'par', 'vd', 'esoil', 'jmax25t', 'jmax25g', 'pc', 'rl',        &
     &  'lambdat', 'lambdag', 'rrt', 'rrg', 'asst', 'assg', 'su_avg',  &
     &  'zw', 'ws', 'spgfcf', 'infx', 'etmt', 'etmg', 'su_1', 'topt', 'ncp'

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
!*-----PARAMETER READING FROM SOILPROFILE.PAR---------------------------

      subroutine vom_get_soilprofile ()
      use vom_vegwat_mod
      implicit none

      INTEGER :: iostat, i, j

!     * The file soilprofile.par can contain information about thickness
!       and soil properties in each soil layer, with the layer number in
!       the first column.

      open(kfile_soilprofile, FILE=trim(adjustl(i_inputpath))// &
           trim(adjustl(sfile_soilprofile)),                  &
     &                        STATUS='old', IOSTAT=iostat)
      if (iostat .eq. 0) then
        do j = 1, s_maxlayer
          read(kfile_soilprofile,*) s_maxlayer, s_delz(j), s_ksat(j),  &
     &      s_nvg(j), s_avg(j), s_thetas(j), s_thetar(j)
        enddo
        c_mvg(:) = 1.d0 - (1.d0 / s_nvg(:))  ! van Genuchten soil parameter m
      else
        s_delz(:)   = i_delz
        s_ksat(:)   = i_ksat
        s_nvg(:)    = i_nvg
        s_avg(:)    = i_avg
        s_thetas(:) = i_thetas
        s_thetar(:) = i_thetar
        c_mvg(:)    = i_mvg
      endif
      close(kfile_soilprofile)

      do i = 1, s_maxlayer
        c_hhydrst(i) = (i - 0.5d0) * s_delz(i)  ! (Out[238]) hydrostatic head for (3.34)
      enddo

      return
      end subroutine vom_get_soilprofile

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----READING TIMESERIES OF SEASONAL VEGETATION COVER---------------------------

      subroutine vom_get_perc_cov ()
      use vom_vegwat_mod
      implicit none

      INTEGER :: iostat, i, j
      INTEGER, ALLOCATABLE :: fyear_pc(:)    ! Year for each day
      INTEGER, ALLOCATABLE :: fmonth_pc(:)   ! Month for each day
      INTEGER, ALLOCATABLE :: fday_pc(:)     ! Day of month
      INTEGER, ALLOCATABLE :: dayyear_pc(:)  ! Day of year



     allocate(dayyear_pc( c_maxday    )   )
     allocate(fday_pc( c_maxday    )   )
     allocate(fmonth_pc( c_maxday    )   )
     allocate(fyear_pc( c_maxday    )   )     
     allocate(perc_cov_veg( c_maxday    )   )


      open(kfile_perc_cov, FILE=trim(adjustl(i_inputpath))// &
           trim(adjustl(sfile_perc_cov)),                  &
     &                        STATUS='old', IOSTAT=iostat)
      if (iostat .eq. 0) then
         do i = 1, c_maxday
          read(kfile_perc_cov,'(4i8,f8.2)') dayyear_pc(i), fday_pc(i), &
     &      fmonth_pc(i), fyear_pc(i), perc_cov_veg(i)
        enddo
      end if 
      close(kfile_perc_cov)

      !check if timeseries match

     if( (fyear_pc(1) .ne. fyear(1)) .or. &
         (fday_pc(1) .ne. fday(1)) .or. &
         (fmonth_pc(1) .ne. fmonth(1)) ) then
     stop 'startdate of perc_cov doesnot match with dailyweather'
     end if

     if( (fyear_pc(c_maxday) .ne. fyear(c_maxday)) .or. &
         (fday_pc(c_maxday) .ne. fday(c_maxday)) .or. &
         (fmonth_pc(c_maxday) .ne. fmonth(c_maxday)) ) then
     stop 'enddate of perc_cov doesnot match with dailyweather'
     end if


      return
      end subroutine vom_get_perc_cov


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
      LOGICAL :: exist_daily
      CHARACTER(len=99) :: str

      inquire(FILE=trim(adjustl(i_inputpath))// &
              trim(adjustl(sfile_hourlyweather)), EXIST=exist)

!     * Creating hourly climate data from daily data

      inquire(FILE=trim(adjustl(i_inputpath))// &
              trim(adjustl(sfile_dailyweather)), EXIST=exist_daily)
      if (.not. exist_daily) then
       stop "dailyweather.prn does not exist"
      end if


      if (.not. exist) then
        open(kfile_dailyweather, FILE=trim(adjustl(i_inputpath))// &
             trim(adjustl(sfile_dailyweather)),              &
     &                           STATUS='old', IOSTAT=stat)
        read(kfile_dailyweather,'(A)') str
        if (LEN(TRIM(str)) .ne. 88) then
          write(0,*) "ERROR: ", TRIM(sfile_dailyweather), ": Wrong file format"
          stop
        endif
        do i = 1, c_maxday
          read(kfile_dailyweather,'(4i8,7f8.2)') dayyear(i), fday(i),  &
     &      fmonth(i), fyear(i), tairmax_d(i), tairmin_d(i), rain_d(i),&
     &      srad_d(i), vp_d(i), press_d(i), ca_d(i)
        enddo
        close(kfile_dailyweather)

!       * Calculation of derived parameters

        call vom_calc_derived()

        if (i_write_h == 1) exist=.TRUE.
      endif

!     * Reading hourly climate data if available

      if (exist) then
        open(kfile_hourlyweather, FILE=trim(adjustl(i_inputpath))// &
             trim(adjustl(sfile_hourlyweather)),            &
     &                            STATUS='old', IOSTAT=stat)
        read(kfile_hourlyweather,*)
        ii = 1
        oldh = 99
        do i = 1, c_maxhour
          read(kfile_hourlyweather,'(5i8,5e11.3)') h, dummyint1,       &
     &      dummyint2, dummyint3, dummyint4, tair_h(i), vd_h(i),       &
     &      par_h(i), rain_h(i), ca_h(i)
          if (par_h(i) .lt. 0.d0) par_h(i) = 0.d0
          ca_h(i) = ca_h(i) / 1.0d6
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
      INTEGER :: in1, in2
      REAL*8  :: sunr, suns
      REAL*8  :: tairmean
      REAL*8  :: dtair
      REAL*8  :: daylength              ! Day length (hours)
      REAL*8  :: vp__                   ! Absolute vapour pressure in the air (Pa)

        in1 = 1
        in2 = c_maxday

      if (i_write_h == 1) then
        open(kfile_hourlyweather, FILE=trim(adjustl(i_inputpath))// &
             trim(adjustl(sfile_hourlyweather)), STATUS='new')
        write(kfile_hourlyweather,'(5a8,5a11)') 'hour', 'dayyear', 'fday', &
     &    'fmonth', 'fyear', 'tair_h', 'vd_h', 'par_h', 'rain_h', 'ca_h'
      endif

      do in = in1, in2
        par_d(in) = 2.0804d0 * srad_d(in)  ! (Out[17]), par in mol/m2 if srad was MJ/m2
        daylength = 12.d0 - 7.639437d0 * ASIN((0.397949d0              &
     &            * COS(0.172142d0 + 0.017214d0 * dayyear(in))         &
     &            * TAN(0.017453d0 * i_lat))                           &
     &            / SQRT(0.920818d0 - 0.079182d0                       &
     &            * COS(0.344284d0 + 0.034428d0 * dayyear(in))))  ! (Out[22]), in hours
!       * sets time of sunrise and sunset
        sunr = 12d0 - 0.5d0 * daylength
        suns = 12d0 + 0.5d0 * daylength
        tairmean = (tairmax_d(in) + tairmin_d(in)) / 2.d0
        dtair = tairmax_d(in) - tairmin_d(in)
        vp__ = vp_d(in) * 100.d0             ! vp in Pa

!       * Loop through every hour of day, where ik=hour
        do ik = 1, 24
          ii = in * 24 + ik - 24
!         * (derived from 3.52+3.53) (Out[38], accounts for diurnal
!           variation in air temperature
          tair_h(ii) = tairmean + dtair * (0.0138d0                           &
     &               * COS(3.513d0 - ((-1.d0 + ik) * p_pi) / 3.d0) + 0.0168d0 &
     &               * COS(0.822d0 - ((-1.d0 + ik) * p_pi) / 4.d0) + 0.0984d0 &
     &               * COS(0.360d0 - ((-1.d0 + ik) * p_pi) / 6.d0) + 0.4632d0 &
     &               * COS(3.805d0 - ((-1.d0 + ik) * p_pi) / 12.d0))

          ca_h(ii) = ca_d(in) / 1.0d6

!         vd_h(ii) = 0.006028127d0 * 2.718282d0 ** ((17.27d0 * tair_h(ii)) &
!    &             / (237.3d0 + tair_h(ii))) - 9.869233d-6 * vp__
!         * (derived from 3.54+3.55) (Out[52]), accounts for diurnal variation in vapour deficit
          vd_h(ii) = (((0.6108d0 * p_E ** (17.27d0 * tair_h(ii)        &
     &             / (tair_h(ii) + 237.3d0))) * 1000) - vp__)          &
     &             / (press_d(in) * 100.d0)
          if (vd_h(ii) .le. 0.d0) vd_h(ii) = 0.d0

!         * average rainfall in hour ii (m/s)
          rain_h(ii) = rain_d(in) / (24.d0 * 3600.d0 * 1000.d0)

          if (sunr .le. ik .and. ik + 1 .le. suns) then
!           * ([Out30]), in mol/m2/s (derived from 3.51) accounts for
!             diurnal variation in global irradiance
            par_h(ii) = (-0.000873d0 * par_d(in) * COS(0.017453d0 * i_lat) &
     &                * SQRT(0.920818d0 - 0.079182d0                   &
     &                * COS(0.034428d0 * (10.d0 + dayyear(in))))       &
     &                * COS(0.2618d0 * ik) - 0.000347d0 * par_d(in)    &
     &                * COS(0.017214d0 * (10.d0 + dayyear(in)))        &
     &                * SIN(0.017453d0 * i_lat)) / (-1.250192d0 * daylength &
     &                * COS(0.017214d0 * (10.d0 + dayyear(in)))        &
     &                * SIN(0.017453d0 * i_lat) + 24.d0                &
     &                * COS(0.017453d0 * i_lat)                        &
     &                * SQRT(0.920818d0 - 0.079182d0                   &
     &                * COS(0.034428d0 * (10.d0 + dayyear(in))))       &
     &                * (1.d0 - (0.158363d0                            &
     &                * COS(0.017214d0 * (10.d0 + dayyear(in))) ** 2.d0 &
     &                * TAN(0.017453d0 * i_lat) ** 2.d0)               &
     &                / (0.920818d0 - 0.079182d0 * COS(0.034428d0      &
     &                * (10.d0 + dayyear(in))))) ** 0.5d0)
            if (par_h(ii) .lt. 0.d0) par_h(ii) = 0.d0
          else
            par_h(ii) = 0.d0
          endif

          if (i_write_h == 1) then
            write(kfile_hourlyweather,'(5i8,5e11.3)') ik, dayyear(in), &
     &        fday(in), fmonth(in), fyear(in), tair_h(ii), vd_h(ii),   &
     &        par_h(ii), rain_h(ii), ca_h(ii)
          endif

        enddo
      enddo

      if (i_write_h == 1) close(kfile_hourlyweather)

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

      topt_ = i_toptstart

!     * Set soil moisture and vegetation parameters to initial conditions

      call waterbalance_init()

      if (o_rtdepth .gt. i_cz) then
        write(*,*) 'Root depth greater than soil depth'
        o_rtdepth = i_cz
      endif

      wsold  = SUM(cH2Ol_s(:))               ! initial soil water storage
      wsnew = wsold

!     * Set vegetation parameters

      q_md   = o_pct * i_mdtf + o_mdstore
      q_mqx  = q_md * i_mqxtf
      mqtnew = 0.95d0 * q_mqx                  ! initial wood water storage
      mqtold = mqtnew
      rsurftnew(:) = 0.d0

!     * Determining the position of the bottom of the tree root zone

      pos_slt = 0
      dummy = 0
      do while (o_rtdepth .gt. dummy)
        pos_slt = pos_slt + 1
        dummy = dummy + s_delz(pos_slt)
      enddo

!     * Determining the position of the bottom of the tree root zone

      pos_slg = 0
      dummy = 0
      do while (o_rgdepth .gt. dummy)
        pos_slg = pos_slg + 1
        dummy = dummy + s_delz(pos_slg)
      enddo

      rsurfgnew(1:pos_slg) = i_rsurf_ * s_delz(1:pos_slg)
      rsurfgnew(pos_slg+1:s_maxlayer) = 0.d0
      if (pos_slg .gt. wlayernew) then
        rsurfgnew(wlayernew+1:pos_slg) = i_rsurfmin * s_delz(wlayernew+1:pos_slt)
      endif

!     * root surface density (root surface area/soil volume) in each sublayer

      rsurftnew(1:pos_slt) = i_rsurf_ * s_delz(1:pos_slt)
      if (pos_slt .gt. wlayernew) then
        rsurftnew(wlayernew+1:pos_slt) = i_rsurfmin * s_delz(wlayernew+1:pos_slt)
      endif
      jmax25t_d(2) = 0.0003d0
      jmax25g_d(2) = 0.0003d0
      c_pcgmin     = 0.02d0             ! minimum grass pc; initial point for growth

      if(i_read_pc == 1) then
         pcg_d(:) = perc_cov_veg( 1 )
         !adjust value if perennial + seasonal > 1
         if( (pcg_d(1) + o_pct) .gt. 1.0) then
            pcg_d(:) = 1.d0 - o_pct
         end if

      else       
         pcg_d(2)     = MIN(1.d0 - o_pct, c_pcgmin)
         pcg_d(:)     = pcg_d(2) + (/-i_incrcovg,0.0d0,i_incrcovg/)  ! vector with values varying by 1%
         pcg_d(3)     = MIN(MAX(c_pcgmin, pcg_d(3)), 1.d0 - o_pct)
      end if

      rootlim(:,:) = 0.d0

!     * Direct costs

!     * (3.38)  foliage tunrover costs, assuming crown LAI of 2.5
      q_tct_d = i_tcf * o_pct * 2.5d0

!     * Setting yearly, daily and hourly parameters

      nyear       = fyear(1)
      rain_y      = 0.d0
      par_y       = 0.d0
      srad_y      = 0.d0
      vd_y        = 0.d0
      etm_y       = 0.d0
      esoil_y     = 0.d0
      ruptkt_d(:) = 0.d0
      ruptkg_d(:) = 0.d0
      asst_d(:)   = 0.d0
      assg_d(:,:) = 0.d0
      ioacum      = 0.d0
!     * for grasses
      etmg_y      = 0.d0
      assg_y      = 0.d0
      rlg_y       = 0.d0
      rrg_y       = 0.d0
      cpccg_y     = 0.d0
      tcg_y       = 0.d0
!     * for trees
      etmt_y      = 0.d0
      asst_y      = 0.d0
      rlt_y       = 0.d0
      rrt_y       = 0.d0
      cpcct_y     = 0.d0
      tct_y       = 0.d0

      return
      end subroutine vom_init_vegpar

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine vom_daily_init ()
      use vom_vegwat_mod
      implicit none

      rsurft_(:)   = rsurftnew(:)
      rsurfg_(:)   = rsurfgnew(:)
      lambdat_d    = o_lambdatf * (SUM(pcapnew(1:pos_slt)) / pos_slt) ** o_wstexp  ! (3.45)
      lambdag_d    = o_lambdagf * pcapnew(1) ** o_wsgexp  ! (3.44)
!     * vector with values varying by 1%
      jmax25t_d(:) = jmax25t_d(2) * (/1.0d0-i_incrjmax,1.0d0,1.0d0+i_incrjmax/)
!     * making sure that the values don't become too low, otherwise
!       they could never pick up again
      jmax25t_d(:) = MAX(jmax25t_d(:), 50.0d-6)
      jmax25g_d(:) = jmax25g_d(2) * (/1.0d0-i_incrjmax,1.0d0,1.0d0+i_incrjmax/)
      jmax25g_d(:) = MAX(jmax25g_d(:), 50.0d-6)


      if( i_read_pc == 1) then   
         pcg_d(:) = perc_cov_veg(nday)

         !adjust value if perennial + seasonal > 1
         if( (pcg_d(1) + o_pct) .gt. 1.0) then
            pcg_d(:) = 1.d0 - o_pct
         end if

      else
         pcg_d(:)     = pcg_d(2) + (/-i_incrcovg,0.0d0,i_incrcovg/)  ! perc. change grass cover
         pcg_d(:)     = MAX(pcg_d(:), 0.d0)
         pcg_d(3)     = MIN(MAX(c_pcgmin, pcg_d(3)), 1.d0 - o_pct)
      end if



!     * (3.38) foliage turnover costs, assuming LAI/pc of 2.5
      tcg_d(:)     = i_tcf * pcg_d(:) * 2.5d0
!     * (3.40), (Out[190])  root respiration [mol/s]
      rrt_d        = 2.55d-7 * SUM(rsurft_(1:pos_slt))

      if (pos_slt .gt. wlayernew) then
!       * (3.42, 2.45e-10 from (Out[165])) costs of water distribution and storage
        q_cpcct_d = i_cpccf * o_pct * o_rtdepth + o_mdstore * 2.45d-10
      else
        q_cpcct_d = i_cpccf * o_pct * SUM(s_delz(1:pos_slt)) + o_mdstore * 2.45d-10
      endif

      if (wlayernew .lt. pos_slg) then
        cpccg_d(:) = i_cpccf * pcg_d(:) * o_rgdepth  ! (3.42) water transport costs
      else
        cpccg_d(:) = i_cpccf * pcg_d(:) * SUM(s_delz(1:pos_slg))
      endif

!     * (3.40), (Out[190]) root respiration grasses [mol/s]
      rrg_d = 2.55d-7 * SUM(rsurfg_(1:pos_slg))
!     * resetting the minimum steady-state tissue water content to
!       its maximum value
      mqsstmin = q_mqx

!     * used for daily recalculation
      if (optmode .eq. 0) then
        tairmax_d(nday) = -9999.d0
        tairmin_d(nday) =  9999.d0
        rain_d(nday)    =     0.d0
        par_d(nday)     =     0.d0
        srad_d(nday)    =     0.d0
      endif

      if (optmode .eq. 0) then
        vd_d     = 0.d0
        etmt_d   = 0.d0
        etmg_d   = 0.d0
        esoil_d  = 0.d0
        spgfcf_d = 0.d0
        infx_d   = 0.d0
        rlt_d    = 0.d0
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

!     * (Out[274], derived from (3.25))
      gammastar = 0.00004275d0                                         &
     &          * p_E ** ((18915.d0 * (-25.d0 + tair_h(th_)))          &
     &          / (149.d0 * p_R_ * (273.d0 + tair_h(th_))))

!     * (Out[310], derived from (3.26)) Temperature dependence of Jmax
      jmaxt_h(:) = (p_E ** ((i_ha * (-25.d0 + tair_h(th_)) * (-273.d0  &
     &           + topt_ + 273.d0 * p_R_ * topt_)) / ((25.d0 + 273.d0  &
     &           * p_R_ * topt_) * (tair_h(th_) + 273.d0 * p_R_        &
     &           * topt_))) * ((-1.d0 + p_E ** (-(i_hd * (-298.d0      &
     &           + topt_)) / (25.d0 + 273.d0 * p_R_ * topt_))) * i_ha  &
     &           + i_hd) * jmax25t_d(:)) / ((-1.d0 + p_E ** ((i_hd     &
     &           * (273.d0 + tair_h(th_) - topt_)) / (tair_h(th_)      &
     &           + 273.d0 * p_R_ * topt_))) * i_ha + i_hd)
!     * (3.24), (Out[312])
      rlt_h(:) = ((ca_h(th_) - gammastar) * o_pct * jmaxt_h(:)         &
     &         * i_rlratio) / (4.d0 * (ca_h(th_) + 2.d0 * gammastar)   &
     &         * (1.d0 + i_rlratio))
!     * (Out[310], derived from (3.26)) Temperature dependence of Jmax
      jmaxg_h(:) = (p_E ** ((i_ha * (-25.d0 + tair_h(th_)) * (-273.d0  &
     &           + topt_ + 273.d0 * p_R_ * topt_)) / ((25.d0 + 273.d0  &
     &           * p_R_ * topt_) * (tair_h(th_) + 273.d0 * p_R_        &
     &           * topt_))) * ((-1.d0 + p_E ** (-(i_hd * (-298.d0      &
     &           + topt_)) / (25.d0 + 273.d0 * p_R_ * topt_))) * i_ha  &
     &           + i_hd) * jmax25g_d(:)) / ((-1.d0 + p_E ** ((i_hd     &
     &           * (273.d0 + tair_h(th_) - topt_)) / (tair_h(th_)      &
     &           + 273.d0 * p_R_ * topt_))) * i_ha + i_hd)
      rlg_h(1,:) = ((ca_h(th_) - gammastar) * pcg_d(1) * jmaxg_h(:)    &
     &           * i_rlratio) / (4.d0 * (ca_h(th_) + 2.d0 * gammastar) &
     &           * (1.d0 + i_rlratio))  ! (3.24), (Out[312])
      rlg_h(2,:) = ((ca_h(th_) - gammastar) * pcg_d(2) * jmaxg_h(:)    &
     &           * i_rlratio) / (4.d0 * (ca_h(th_) + 2.d0 * gammastar) &
     &           * (1.d0 + i_rlratio))  ! (3.24), (Out[312])
      rlg_h(3,:) = ((ca_h(th_) - gammastar) * pcg_d(3) * jmaxg_h(:)    &
     &           * i_rlratio) / (4.d0 * (ca_h(th_) + 2.d0 * gammastar) &
     &           * (1.d0 + i_rlratio))  ! (3.24), (Out[312])

!     * daily recalculation for resultsdaily
      if (optmode .eq. 0) then
        rain_d(nday) = rain_d(nday) + rain_h(th_) * 3600.d0 * 1000.d0  ! mm/d
        par_d(nday)  = par_d(nday) + par_h(th_) * 3600.d0  ! in mol/m2/d
        srad_d(nday) = par_d(nday) / 2.0804d0  ! MJ/m2/d

        if (tair_h(th_) .gt. tairmax_d(nday)) tairmax_d(nday) = tair_h(th_)
        if (tair_h(th_) .lt. tairmin_d(nday)) tairmin_d(nday) = tair_h(th_)
      endif

      if (optmode .eq. 0) then
        mqtold      = mqtnew
        spgfcf_h    = 0.d0
        infx_h      = 0.d0
        io_h        = 0.d0
        esoil_h     = 0.d0
        etmt_h      = 0.d0
        etmg_h      = 0.d0
        sumruptkt_h = 0.d0
      endif

      time        = 0.d0
      asst_h(:)   = 0.d0
      assg_h(:,:) = 0.d0
      ruptkt_h(:) = 0.d0
      ruptkg_h(:) = 0.d0

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

      if (par_h(th_) .gt. 0.d0) then
!       * adaptation of topt to air temperature during sunlight
        topt_ = topt_ + i_toptf * (tair_h(th_) + 273.d0 - topt_)
        jactt(:)   = (1.d0 - p_E ** (-(i_alpha * par_h(th_))           &
     &             / jmaxt_h(:))) * jmaxt_h(:) * o_pct  ! (3.23), (Out[311])
        jactg(1,:) = (1.d0 - p_E ** (-(i_alpha * par_h(th_))           &
     &             / jmaxg_h(:))) * jmaxg_h(:) * pcg_d(1)  ! (3.23), (Out[311])
        jactg(2,:) = (1.d0 - p_E ** (-(i_alpha * par_h(th_))           &
     &             / jmaxg_h(:))) * jmaxg_h(:) * pcg_d(2)  ! (3.23), (Out[311])
        jactg(3,:) = (1.d0 - p_E ** (-(i_alpha * par_h(th_))           &
     &             / jmaxg_h(:))) * jmaxg_h(:) * pcg_d(3)  ! (3.23), (Out[311])

        cond1      = (2.d0 * p_a * vd_h(th_)) / (ca_h(th_) + 2.d0 * gammastar)
        cond2      = (4.d0 * ca_h(th_) * rlt_h(2) + 8.d0 * gammastar   &
     &             * rlt_h(2)) / (ca_h(th_) - gammastar)
        cond3(:,:) = (4.d0 * ca_h(th_) * rlg_h(:,:) + 8.d0 * gammastar &
     &             * rlg_h(:,:)) / (ca_h(th_) - gammastar)

        if (vd_h(th_) .gt. 0.d0 .and. lambdat_d .gt. cond1 .and. jactt(2) .gt. cond2) then

!         gstomt = MAX(0.d0, (0.25d0 * (p_a * (ca_h(th_) * (jactt(2)   &
!    &           - 4.d0 * rlt_h(2)) - 4.d0 * gammastar * (jactt(2)     &
!    &           + 2.d0 * rlt_h(2))) * vd_h(th_) * (ca_h(th_)          &
!    &           * lambdat_d + 2.d0 * gammastar * lambdat_d - p_a      &
!    &           * vd_h(th_)) + 1.7320508075688772d0 * SQRT(p_a        &
!    &           * gammastar * jactt(2) * (ca_h(th_) * (jactt(2)       &
!    &           - 4.d0 * rlt_h(2)) - gammastar * (jactt(2) + 8.d0     &
!    &           * rlt_h(2))) * vd_h(th_) * (ca_h(th_) * lambdat_d     &
!    &           + 2.d0 * gammastar * lambdat_d - 2.d0 * p_a           &
!    &           * vd_h(th_)) ** 2.d0 * (ca_h(th_) * lambdat_d + 2.d0  &
!    &           * gammastar * lambdat_d - p_a * vd_h(th_)))))         &
!    &           / (p_a * (ca_h(th_) + 2.d0 * gammastar) ** 2.d0       &
!    &           * vd_h(th_) * (ca_h(th_) * lambdat_d + 2.d0           &
!    &           * gammastar * lambdat_d - p_a * vd_h(th_))))

          part1 = ca_h(th_) + 2.d0 * gammastar
          part2 = part1 * lambdat_d - p_a * vd_h(th_)
          part3 = p_a * vd_h(th_) * part2

          part4 = ca_h(th_) * (jactt(2) - 4.d0 * rlt_h(2))
          part5 = gammastar * jactt(2)
          part6 = gammastar * 8.d0 * rlt_h(2)
          part7 = part4 - part5 - part6

          part8 = SQRT(part5 * part7 * (part2 - p_a * vd_h(th_)) ** 2.d0 * part3)
          part9 = part7 - 3.d0 * part5 + 1.7320508075688772d0 * part8 / part3

          gstomt = 0.25d0 * part9 / part1**2.d0
          gstomt = MAX(0.d0, gstomt)    ! (Out[314])

        else
          gstomt = 0.d0
        endif
        transpt = p_a * vd_h(th_) * gstomt  ! (3.28) transpiration rate in mol/s
        etmt__ = (transpt * 18.d0) / (10.d0 ** 6.d0)  ! transpiration rate in m/s

        where (vd_h(th_) .gt. 0.d0 .and. lambdag_d .gt. cond1 .and. jactg(:,:) .gt. cond3(:,:))
          gstomg(:,:) = MAX(0.d0,(0.25d0 * (p_a * (ca_h(th_)           &
     &                * (jactg(:,:) - 4.d0 * rlg_h(:,:)) - 4.d0        &
     &                * gammastar * (jactg(:,:) + 2.d0 * rlg_h(:,:)))  &
     &                * vd_h(th_) * (ca_h(th_) * lambdag_d + 2.d0      &
     &                * gammastar * lambdag_d - p_a * vd_h(th_))       &
     &                + 1.7320508075688772d0 * SQRT(p_a * gammastar    &
     &                * jactg(:,:) * (ca_h(th_) * (jactg(:,:) - 4.d0   &
     &                * rlg_h(:,:)) - gammastar * (jactg(:,:) + 8.d0   &
     &                * rlg_h(:,:))) * vd_h(th_) * (ca_h(th_)          &
     &                * lambdag_d + 2.d0 * gammastar * lambdag_d       &
     &                - 2.d0 * p_a * vd_h(th_)) ** 2.d0 * (ca_h(th_)   &
     &                * lambdag_d + 2.d0 * gammastar * lambdag_d - p_a &
     &                * vd_h(th_))))) / (p_a * (ca_h(th_) + 2.d0       &
     &                * gammastar) ** 2.d0 * vd_h(th_) * (ca_h(th_)    &
     &                * lambdag_d + 2.d0 * gammastar * lambdag_d       &
     &                - p_a * vd_h(th_))))  ! (Out[314])
        elsewhere
          gstomg(:,:) = 0.d0
        endwhere
        transpg(:,:) = p_a * vd_h(th_) * gstomg(:,:)  ! (3.28) transpiration rate in mol/s
        etmg__(:,:) = (transpg(:,:) * 18.d0) / (10.d0 ** 6.d0)  ! transpiration rate in m/s
      else
        jactt(:)    = 0.d0
        gstomt      = 0.d0
        etmt__      = 0.d0
        jactg(:,:)  = 0.d0
        gstomg(:,:) = 0.d0
        etmg__(:,:) = 0.d0
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

      if (wlayernew .lt. pos_slt) then
        rsurft_(wlayernew+1:pos_slt) = i_rsurfmin * s_delz(wlayernew+1:pos_slt)
      endif
      if (wlayernew .lt. pos_slg) then
        rsurfg_(wlayernew+1:pos_slg) = i_rsurfmin * s_delz(wlayernew+1:pos_slt)
      endif

      mqt_       = mqtnew
      zw_        = zwnew
      wlayer_    = wlayernew
      su__(:)    = sunew(:)
      pcap_(:)   = pcapnew(:)
      kunsat_(:) = kunsatnew(:)

      return
      end subroutine vom_subhourly_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----root water uptake------------------------------------------------

      subroutine vom_rootuptake ()
      use vom_vegwat_mod
      implicit none

      if (wlayernew .ge. 1) then
        pos_ult = MIN(pos_slt, wlayer_)
        if (q_md .gt. 0.d0) then
          prootm(1:pos_ult) = (p_mpbar * (-mqt_ + q_mqx) * (750.d0     &
     &                      - (750.d0 * q_mqx) / (q_md + q_mqx)        &
     &                      + (q_md + q_mqx) / q_mqx)) / (q_md + q_mqx)&
     &                      - c_hhydrst(1:pos_ult)  ! (Out[239])
        else
!         * set tissue suction to the same as in grasses, if no storage capacity
          prootm(1:pos_ult) = i_prootmg
        endif

!       * soil resistance, (Out[ 241] with svolume=s_delz(1:pos_ult)); derived from (3.32)
        rsoil(1:pos_ult) = SQRT(p_pi / 2.d0) * SQRT((i_rootrad         &
     &                   * s_delz(1:pos_ult)) / rsurft_(1:pos_ult))    &
     &                   / kunsat_(1:pos_ult)

!       * root water uptake, Chapter 3.3.3.3 (Out[242])
        if (q_md .gt. 0.d0) then
          ruptkt__(1:pos_ult) = (-pcap_(1:pos_ult) + prootm(1:pos_ult))&
     &                        * rsurft_(1:pos_ult) / (i_rrootm         &
     &                        + rsoil(1:pos_ult))
          ruptkt__(pos_ult+1:s_maxlayer) = 0.d0
        else  ! if no storage, uptake happens only when etmt__>0

          if (etmt__ .gt. 0.d0) then
            ruptkt__(1:pos_ult) = MAX(0.d0,(-pcap_(1:pos_ult)          &
     &                          + prootm(1:pos_ult)) * rsurft_(1:pos_ult) &
     &                          / (i_rrootm + rsoil(1:pos_ult)))
            ruptkt__(pos_ult+1:s_maxlayer) = 0.d0

            if (SUM(ruptkt__(:)) .gt. 0.d0) then
              if (etmt__ .gt. SUM(ruptkt__(:))) then
                changef = 1.d0
                etmt__   = SUM(ruptkt__(:))
                transpt = etmt__ * 55555.555555555555d0  ! (Out[249]) mol/s=m/s*10^6 g/m/(18g/mol)
                gstomt = transpt / (p_a * vd_h(th_))
              endif
!             * Setting SUM(ruptkt__)=etmt__ and distributing according to relative uptake:
              ruptkt__(:) = etmt__ * (ruptkt__(:) / (SUM(ruptkt__(:))))
            else
              ruptkt__(:) = 0.d0
              changef     = 1.d0
              etmt__      = 0.d0
              transpt     = 0.d0
              gstomt      = 0.d0
            endif

          else
            ruptkt__(:) = 0.d0
          endif

        endif

        pos_ulg = MIN(pos_slg, wlayer_)
        if (MAXVAL(etmg__(:,:)) .gt. 0.d0) then
!         * root uptake by grasses can not be negative, as storage negligible
          ruptkg__(1:pos_slg) = MAX(0.d0,((-pcap_(1:pos_ulg)           &
     &                        + (i_prootmg - c_hhydrst(1:pos_ulg)))    &
     &                        * rsurfg_(:)) / (i_rrootm + (SQRT(p_pi / 2.d0)  &
     &                        * SQRT(i_rootrad * s_delz(1:pos_ulg)     &
     &                        / rsurfg_(:))) / kunsat_(1:pos_ulg)))
          ruptkg__(pos_ulg+1:s_maxlayer) = 0.d0
          if (SUM(ruptkg__(:)) .gt. 0.d0) then
            where (etmg__(:,:) .gt. SUM(ruptkg__(:)))
              rootlim(:,:)  = 1.d0
              etmg__(:,:)   = SUM(ruptkg__(:))
              transpg(:,:)  = etmg__(:,:) * 55555.555555555555d0  ! (Out[249]) mol/s=m/s*10^6 g/m/(18g/mol)
              gstomg(:,:)   = transpg(:,:) / (p_a * vd_h(th_))
            end where
            ruptkg__(1:pos_ulg) = etmg__(2,2) * (ruptkg__(1:pos_ulg)   &
     &                          / (SUM(ruptkg__(:))))
          else
            ruptkg__(:)  = 0.d0
            etmg__(:,:)  = 0.d0
            transpg(:,:) = 0.d0
            gstomg(:,:)  = 0.d0
          endif
        else
          ruptkg__(:) = 0.d0
        endif
      else
        ruptkg__(:)  = 0.d0
        ruptkt__(:)  = 0.d0
        etmg__(:,:)  = 0.d0
        transpg(:,:) = 0.d0
        gstomg(:,:)  = 0.d0
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

!     mqss_out = MAX(0.9d0 * q_mqx,(q_mqx * (p_mpbar * (q_md * q_md    &
!    &         + 752.d0 * q_md * q_mqx + q_mqx * q_mqx)                &
!    &         * SUM((rsurft_(1:pos_ult) / (i_rrootm + rsoil(1:pos_ult)))) &
!    &         - (q_md + q_mqx) * (q_md + q_mqx) * (etmt__             &
!    &         - SUM(((-c_hhydrst(1:pos_ult) - pcap_(1:pos_ult))       &
!    &         * rsurft_(1:pos_ult)) / (i_rrootm + rsoil(1:pos_ult)))))) &
!    &         / (p_mpbar * (q_md * q_md + 752.d0 * q_md * q_mqx       &
!    &         + q_mqx * q_mqx) * SUM((rsurft_(1:pos_ult) / (i_rrootm  &
!    &         + rsoil(1:pos_ult))))))

      sum1 = SUM(rsurft_(1:pos_ult) / (i_rrootm + rsoil(1:pos_ult)))
      mul1 = p_mpbar * (q_md * q_md + 752.d0 * q_md * q_mqx + q_mqx * q_mqx) * sum1

      sum2 = SUM(((-c_hhydrst(1:pos_ult) - pcap_(1:pos_ult))           &
     &     * rsurft_(1:pos_ult)) / (i_rrootm + rsoil(1:pos_ult)))
      mul2 = (q_md + q_mqx) * (q_md + q_mqx) * (etmt__ - sum2)

      mqss_out = q_mqx * (mul1 - mul2) / mul1
      mqss_out = MAX(0.9d0 * q_mqx, mqss_out)

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
      CHARACTER(len=135) :: msg

!     * makes sure that tissue water does not get below 0.9mqx
      if (mqt_ .le. 0.9d0 * q_mqx) then
        if (wlayer_ .ge. 1) then
          if (etmt__ .gt. 0.9d0 * SUM(ruptkt__(:))) then
            if (SUM(ruptkt__(:)) .ge. 0.d0) then
              etmt__ = SUM(ruptkt__(:))
              transpt = etmt__ * 55555.555555555555d0  ! (Out[249]) mol/s=m/s*10^6 g/m/(18g/mol)
              gstomt = transpt / (p_a * vd_h(th_))
            else
              write(msg,'(A20,I2,A1,I2,A1,I4)') 'vegetation dies on: ', &
     &          fday(nday), '/', fmonth(nday), '/', fyear(nday)
              write(*,*) TRIM(msg)
              netass = 0.d0
!             * if tissues water depleted, but still loosing water -> death
              finish = 1
              return
            endif
            call vom_mqss(mqsst_)
            mqsstmin = MIN(mqsstmin, mqsst_)
          endif
        else
          etmt__  = 0.d0
          transpt = 0.d0
          gstomt  = 0.d0
        endif
      endif
      if (wlayer_ .ge. 0) then
!       * (3.35), 1.e6 to convert from m (=1000kg/m2) to g/m2; (Out[250])
        dmqt = (SUM(ruptkt__(:)) - etmt__) * 1.d6
      else
        dmqt = -etmt__ * 1.d6
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

      if (q_md .gt. 0.d0) then
!       * avoids mq from becoming larger than mqx or smaller than 0.9mqx
        if (dmqt .gt. 0.d0) then
          dtmq = (q_mqx - mqt_) / dmqt
        elseif (dmqt .lt. 0.d0) then
          dtmq = (0.9d0 * q_mqx - mqt_) / dmqt
        endif

        if (ABS(mqt_ - mqsst_) .gt. q_mqx / 1.d6) then
          dtss = (mqt_ - mqsst_) / (1.d6 * (etmt__ - SUM(ruptkt__(:))))
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

      REAL*8 :: asst__(3)
      REAL*8 :: assg__(3,3)

      asst__(:) = (4.d0 * ca_h(th_) * gstomt + 8.d0 * gammastar        &
     &          * gstomt + jactt(:) - 4.d0 * rlt_h(:) - SQRT((-4.d0    &
     &          * ca_h(th_) * gstomt + 8.d0 * gammastar * gstomt       &
     &          + jactt(:) - 4.d0 * rlt_h(:)) ** 2.d0 + 16.d0          &
     &          * gammastar * gstomt * (8.d0 * ca_h(th_) * gstomt      &
     &          + jactt(:) + 8.d0 * rlt_h(:)))) / 8.d0  ! (3.22) ; (Out[319])
      asst_h(:) = asst_h(:) + asst__(:) * dt_
      assg__(:,:) = (4.d0 * ca_h(th_) * gstomg(:,:) + 8.d0 * gammastar &
     &            * gstomg(:,:) + jactg(:,:) - 4.d0 * rlg_h(:,:)       &
     &            - SQRT((-4.d0 * ca_h(th_) * gstomg(:,:) + 8.d0       &
     &            * gammastar * gstomg(:,:) + jactg(:,:) - 4.d0        &
     &            * rlg_h(:,:)) ** 2.d0 + 16.d0 * gammastar            &
     &            * gstomg(:,:) * (8.d0 * ca_h(th_) * gstomg(:,:)      &
     &            + jactg(:,:) + 8.d0 * rlg_h(:,:)))) / 8.d0  ! (3.22); (Out[319])
      assg_h(:,:) = assg_h(:,:) + assg__(:,:) * dt_
      ruptkt_h(:) = ruptkt_h(:) + ruptkt__(:) * dt_
      ruptkg_h(:) = ruptkg_h(:) + ruptkg__(:) * dt_
      if (optmode .eq. 0) then
        spgfcf_h    = spgfcf_h    + dt_ * spgfcf__
        infx_h      = infx_h      + dt_ * infx__
        io_h        = io_h        + dt_ * io__
        esoil_h     = esoil_h     + dt_ * esoil__
        etmt_h      = etmt_h      + dt_ * etmt__
        etmg_h      = etmg_h      + dt_ * etmg__(2,2)
        sumruptkt_h = sumruptkt_h + dt_ * SUM(ruptkt__(:))
      endif

      return
      end subroutine vom_add_hourly

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine vom_add_daily ()
      use vom_vegwat_mod
      implicit none

      vd_d     = vd_d     + vd_h(th_)
      etmt_d   = etmt_d   + etmt_h
      etmg_d   = etmg_d   + etmg_h
      esoil_d  = esoil_d  + esoil_h
      spgfcf_d = spgfcf_d + spgfcf_h
      infx_d   = infx_d   + infx_h
      rlt_d    = rlt_d    + rlt_h(2)   * 3600.d0  ! rlt_d in mol/day
      rlg_d    = rlg_d    + rlg_h(2,2) * 3600.d0

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

      if (fyear(nday) .ge. i_firstyear .and. fyear(nday) .le. i_lastyear) then
!       * internal write to convert from number to string
        write(str,'(i3)') wlayer_
!       * includes a column for each sublayer
        hourlyformat = '(I6,I6,I4,I7,I5,'//str//'E14.6)'

        write(kfile_resultshourly,'(I6,I7,I7,I7,I7,22E15.5)')          &
     &    fyear(nday), fmonth(nday), fday(nday), nday, nhour,          &
     &    rain_h(th_), tair_h(th_), par_h(th_), vd_h(th_), esoil_h,    &
     &    o_pct + pcg_d(2), jmax25t_d(2), jmax25g_d(2), mqt_,          &
     &    rlt_h(2) + rlg_h(2,2), lambdat_d, lambdag_d, rrt_d + rrg_d,  &
     &    asst_h(2), assg_h(2,2), etmt_h, etmg_h, su__(1), zw_, wsnew, &
     &    spgfcf_h, infx_h

          write(kfile_delzhourly,hourlyformat) fyear(nday),            &
     &      fmonth(nday), fday(nday), nday, nhour, s_delz(1:wlayer_)

          write(kfile_ruptkthourly,hourlyformat) fyear(nday),          &
     &      fmonth(nday), fday(nday), nday, nhour, ruptkt_h(1:wlayer_)

          write(kfile_suhourly,hourlyformat) fyear(nday),              &
     &      fmonth(nday), fday(nday), nday, nhour, su__(1:wlayer_)
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

      REAL*8             :: error1
      CHARACTER(len=135) :: msg

      ioacum = ioacum + io_h
      wsnew  = SUM(cH2Ol_s(:))
      error  = wsold + ioacum - wsnew

!     * gives an error message if accumulated error exceeds 1 mm

      if (abs(error) .gt. 1.d-3) then
        write(msg,*) 'Error in water balance [mm]:', error, 'io=',     &
     &    io__, 'wsold=', wsold, 'wsnew=', wsnew
        write(*,*) TRIM(msg)
        finish = 1
      elseif (q_md .gt. 0.d0) then
        error1 = mqtold + (sumruptkt_h - etmt_h) * 1.d6 - mqtnew
        if (abs(error1 / mqtnew) .gt. 1.d-6) then
          write(msg,*) 'Error in tree water balance [%]:',             &
     &      error1 * 100.d0, 'mqtold=', mqtold, 'mqtnew=', mqtnew,     &
     &      'hruptk=', sumruptkt_h, 'hetm=', etmt_h
          write(*,*) TRIM(msg)
          finish = 1
        endif
      endif

      return
      end subroutine vom_check_water

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


      subroutine vom_write_dayyear ( tp_netass )
      use vom_vegwat_mod
      implicit none

      REAL*8,  INTENT(in) :: tp_netass
      CHARACTER(60) :: dailyformat
      CHARACTER(3)  :: str

!     * internal write to convert from number to string
      write(str,'(I3)') wlayer_
!     * includes a column for each sublayer
      dailyformat = '(I6,I6,I4,I7,'//str//'E14.6)'

      write(kfile_resultsdaily,'(I6,I7,I7,I7,I7,26E15.5)')             &
     &  fyear(nday), fmonth(nday), fday(nday), nday, nhour-1,          &
     &  rain_d(nday), tairmax_d(nday), tairmin_d(nday), par_d(nday),   &
     &  vd_d / 24.d0, esoil_d, jmax25t_d(2), jmax25g_d(2),             &
     &  o_pct + pcg_d(2), rlt_d + rlg_d, lambdat_d, lambdag_d,         &
     &  rrt_d * 3600.d0 * 24.d0, rrg_d * 3600.d0 * 24.d0, asst_d(2),   &
     &  assg_d(2,2), SUM(su__(1:wlayer_)) / wlayer_, zw_, wsnew,       &
     &  spgfcf_d, infx_d, etmt_d, etmg_d, su__(1), topt_, tp_netass

        write(kfile_rsurfdaily,dailyformat) fyear(nday), fmonth(nday), &
     &    fday(nday), nday, rsurft_(1:wlayer_)

      if (fyear(nday) .ne. nyear) then
!       * for calculation of vd_y a -1 is added to nday for using dayyear of correct year
        write(kfile_resultsyearly,'(i6,18e16.6)') nyear, rain_y,       &
     &    par_y, srad_y, vd_y / (dayyear(nday-1)), esoil_y, etm_y,     &
     &    etmg_y, assg_y, rlg_y, rrg_y, cpccg_y, tcg_y,                &
     &    etmt_y, asst_y, rlt_y, rrt_y, cpcct_y, tct_y
      endif

!     * WRITING THE ACCUMULATED DATA FROM THE LAST YEAR TO FILE:

      if (nday .eq. c_maxday) then
!       * call subroutine there to get yearly data for the output
        call vom_add_yearly()
        write(kfile_resultsyearly,'(i6,18e16.6)') nyear, rain_y,       &
     &    par_y, srad_y, vd_y / (dayyear(nday)), esoil_y, etm_y,       &
     &    etmg_y, assg_y, rlg_y, rrg_y, cpccg_y, tcg_y,                &
     &    etmt_y, asst_y, rlt_y, rrt_y, cpcct_y, tct_y
      endif

      return
      end subroutine vom_write_dayyear

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*----- saving fluxes to a single matrix -------------------------------

     subroutine vom_save_dayyear ()
     !subroutine to create one single matrix with all daily output variables
     !The saved matrix will be needed to write txt-files containing results from 
     !all runs. 

      use vom_vegwat_mod
      use vom_sce_mod
      implicit none

      CHARACTER(60) :: dailyformat
      CHARACTER(3)  :: str
      REAL*8, dimension(:), allocatable        :: tmp
      !REAL*8, DIMENSION(21, c_maxday ), intent(inout)        :: output

!     * internal write to convert from number to string
      write(str,'(I3)') wlayer_
!     * includes a column for each sublayer
      dailyformat = '(I6,I6,I4,I7,'//str//'E14.6)'

       !transpiration
       output_mat(1, nday) = etmt_d

       !daily atmospheric vapour deficit
       output_mat(1, nday) = vd_d / 24.d0

       !daily soil evaporation rate
       output_mat(2, nday) = esoil_d

       ! Tree photosynthetic electron transport capacity at 25oC
       output_mat(3, nday) = jmax25t_d(2)

       ! Grass photosynthetic electron transport capacity at 25oC
       output_mat(4, nday) = jmax25g_d(2)

       ! Projected cover perennial vegetation plus actual cover seasonal vegetation 
       output_mat(5, nday) = o_pct + pcg_d(2)

       ! Daily tree plus grass leaf respiration
       output_mat(6, nday) = rlt_d + rlg_d

       ! Target dE/dA for calculating gstomt (tree)
       output_mat(7, nday) = lambdat_d

       ! Target dE/dA for calculating gstomt (grass)
       output_mat(8, nday) = lambdag_d

       ! Tree root respiration rate
       output_mat(9, nday) = rrt_d * 3600.d0 * 24.d0

       ! Grass root respiration rate
       output_mat(10, nday) = rrg_d * 3600.d0 * 24.d0

       ! Daily tree assimilation
       output_mat(11, nday) = asst_d(2)

       ! Daily grass assimilation
       output_mat(12, nday) = assg_d(2,2)

       ! Average soil moisture
       output_mat(13, nday) = SUM(su__(1:wlayer_)) / wlayer_

       ! Elevation of water table
       output_mat(14, nday) = zw_

       ! Total soil water storage
       output_mat(15, nday) = wsnew

       ! Daily seepage face flow
       output_mat(16, nday) = spgfcf_d

       ! Daily infiltration excess runoff
       output_mat(17, nday) = infx_d

       ! Daily transpiration rate (tree)
       output_mat(18, nday) = etmt_d

       ! Daily transpiration rate (grass)
       output_mat(19, nday) = etmg_d

       ! Soil saturation degree in first layer
       output_mat(20, nday) = su__(1)

       ! Optimal temperature in temperature response curve
       output_mat(21, nday) = topt_
 

    !write to file
    if(nday .eq. c_maxday) then


      if( vd_d_out .eqv. .TRUE.) then
       write(kfile_vd_d,*) output_mat(1,:)
      end if
      if( esoil_out .eqv. .TRUE.) then
       write(kfile_esoil,*) output_mat(2,:)
      end if
      if( jmax25t_out .eqv. .TRUE.) then
       write(kfile_jmax25t,*) output_mat(3,:)
      end if
      if( jmax25g_out .eqv. .TRUE.) then
       write(kfile_jmax25g,*) output_mat(4,:)
      end if
      if( vegcov_out .eqv. .TRUE.) then
       write(kfile_vegcov,*) output_mat(5,:)
      end if
      if( resp_out .eqv. .TRUE.) then
       write(kfile_resp,*) output_mat(6,:)
      end if
      if( lambdat_out .eqv. .TRUE.) then
       write(kfile_lambdat,*) output_mat(7,:)
      end if
      if( lambdag_out .eqv. .TRUE.) then
       write(kfile_lambdag,*) output_mat(8,:)
      end if
      if( rrt_out .eqv. .TRUE.) then
       write(kfile_rrt,*) output_mat(9,:)
      end if
      if( rrg_out .eqv. .TRUE.) then
       write(kfile_rrg,*) output_mat(10,:)
      end if
      if( asst_out .eqv. .TRUE.) then
       write(kfile_asst,*) output_mat(11,:)
      end if
      if( assg_out .eqv. .TRUE.) then
       write(kfile_assg,*) output_mat(12,:)
      end if
      if( su_av_out .eqv. .TRUE.) then
       write(kfile_su_av,*) output_mat(13,:)
      end if
      if( zw_out .eqv. .TRUE.) then
       write(kfile_zw,*) output_mat(14,:)
      end if
      if( wsnew_out .eqv. .TRUE.) then
       write(kfile_wsnew,*) output_mat(15,:)
      end if
      if( spgfcf_out .eqv. .TRUE.) then
       write(kfile_spgfcf,*) output_mat(16,:)
      end if
      if( infx_out .eqv. .TRUE.) then
       write(kfile_infx,*) output_mat(17,:)
      end if
      if( etmt_out .eqv. .TRUE.) then
       write(kfile_etmt,*) output_mat(18,:)
      end if
      if( etmg_out .eqv. .TRUE.) then
       write(kfile_etmg,*) output_mat(19,:)
      end if
      if( su1_out .eqv. .TRUE.) then
       write(kfile_su1,*) output_mat(20,:)
      end if
      if( topt_out .eqv. .TRUE.) then
       write(kfile_topt,*) output_mat(21,:)
      end if

    end if 

      return
      end subroutine vom_save_dayyear

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine vom_add_yearly ()
      use vom_vegwat_mod
      implicit none

      if (fyear(nday) .eq. nyear) then
        rain_y   = rain_y  + rain_d(nday)  ! in [mm]
        par_y    = par_y   + par_d(nday)
        srad_y   = srad_y  + srad_d(nday)  ! srad originally in MJ/day
        vd_y     = vd_y    + vd_d / 24.d0
        etm_y    = etm_y   + (etmt_d + etmg_d) * 1000.d0  ! in[mm]
        esoil_y  = esoil_y + esoil_d           * 1000.d0  ! in [mm]
!       * for grasses
        etmg_y   = etmg_y  + etmg_d * 1000.d0 ! in [mm]
        assg_y   = assg_y  + assg_d(2,2)
        rlg_y    = rlg_y   + rlg_d
        rrg_y    = rrg_y   + rrg_d      * 3600.d0 * 24.d0
        cpccg_y  = cpccg_y + cpccg_d(2) * 3600.d0 * 24.d0
        tcg_y    = tcg_y   + tcg_d(2)   * 3600.d0 * 24.d0
!       * for trees
        etmt_y   = etmt_y  + etmt_d * 1000.d0  ! in [mm]
        asst_y   = asst_y  + asst_d(2)
        rlt_y    = rlt_y   + rlt_d
        rrt_y    = rrt_y   + rrt_d     * 3600.d0 * 24.d0
        cpcct_y  = cpcct_y + q_cpcct_d * 3600.d0 * 24.d0
        tct_y    = tct_y   + q_tct_d   * 3600.d0 * 24.d0
      else
        nyear    = fyear(nday)
        rain_y   = rain_d(nday)
        par_y    = par_d(nday)
        srad_y   = srad_d(nday)              ! srad originally in MJ/day
        vd_y     = vd_d / 24.d0
        etm_y    = (etmg_d + etmt_d) * 1000.d0
        esoil_y  = esoil_d * 1000.d0
!       * for grasses
        etmg_y   = etmg_d * 1000.d0
        assg_y   = assg_d(2,2)
        rlg_y    = rlg_d
        rrg_y    = rrg_d      * 3600.d0 * 24.d0
        cpccg_y  = cpccg_d(2) * 3600.d0 * 24.d0
        tcg_y    = tcg_d(2)   * 3600.d0 * 24.d0
!       * for trees
        etmt_y = etmt_d * 1000.d0
        asst_y   = asst_d(2)
        rlt_y    = rlt_d
        rrt_y    = rrt_d     * 3600.d0 * 24.d0
        cpcct_y  = q_cpcct_d * 3600.d0 * 24.d0
        tct_y    = q_tct_d   * 3600.d0 * 24.d0
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

      REAL*8  :: netcg_d(3,3)           ! Daily grass net carbon profit
      INTEGER :: posma(1)               ! Pointer to variable values that achieved maximum assimilation

      posma(:)     = MAXLOC(asst_d(:))
      jmax25t_d(2) = jmax25t_d(posma(1))
      asst_d(:)    = 0.d0
      netcg_d(1,:) = assg_d(1,:) - 3600.d0 * 24.d0 * (cpccg_d(1) + rrg_d + tcg_d(1))
      netcg_d(2,:) = assg_d(2,:) - 3600.d0 * 24.d0 * (cpccg_d(2) + rrg_d + tcg_d(2))
      netcg_d(3,:) = assg_d(3,:) - 3600.d0 * 24.d0 * (cpccg_d(3) + rrg_d + tcg_d(3))
      posmna(:)    = MAXLOC(netcg_d(:,:))
      pcg_d(2)     = MIN(1.d0 - o_pct, pcg_d(posmna(1)))
      jmax25g_d(2) = jmax25g_d(posmna(2))
      assg_d(:,:)  = 0.d0

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

      refft(:) = 0.d0
      if (q_md .gt. 0.d0) then          ! if q_md=0, then changef is calculated elsewhere
        changef = (0.95d0 * q_mqx - mqsstmin) / (0.05d0 * q_mqx)  ! (3.47)
      else
!       * changef for q_md=0 is either 0 or 1. Change it to be either -1 or + 1:
        changef = 2.d0 * changef - 1.d0
      endif
      refft(:) = 0.d0
      maxval_tmp = MAXVAL(ruptkt_d(1:pos_slt) / rsurft_(1:pos_slt))
      if (maxval_tmp .ne. 0.d0) then
        refft(1:pos_slt) = 0.5d0 * ruptkt_d(1:pos_slt) / rsurft_(1:pos_slt) / maxval_tmp  ! (3.48)
      endif
      where (ruptkt_d(1:pos_slt) .lt. 0.d0)
        refft(:) = 0.d0
      end where
      if (changef .lt. 0.d0) then
        refft(:) = 1.d0 - refft(:)
      endif

!     * rsurf=(2*c_epsln/i_rootrad) if all pores filled by roots

      rsurftnew(1:pos_slt) = MIN(2.d0 * c_epsln / i_rootrad * s_delz(1:pos_slt), &
     &                     MAX(i_rsurfmin * s_delz(1:pos_slt), rsurft_(1:pos_slt) &
     &                     + rsurft_(1:pos_slt) * i_growthmax * changef &
     &                     * refft(1:pos_slt) * s_delz(1:pos_slt)))

!     *-----SEASONAL VEGETATION---------------

!     * rootlim is either 0 or 1. Change it to be either -1 or + 1:
      rootlim(posmna(1),posmna(2)) = 2.d0 * rootlim(posmna(1),posmna(2)) - 1.d0

      reffg(:) = 0.d0
      maxval_tmp = MAXVAL(ruptkg_d(1:pos_slg) / rsurfg_(1:pos_slg))
      if (maxval_tmp .ne. 0.d0) then
        reffg(1:pos_slg) = 0.5d0 * ruptkg_d(1:pos_slg) / rsurfg_(1:pos_slg) / maxval_tmp  ! (3.48)
      endif

!     * if roots are going to be reduced, reverse effectivity vector

      if (rootlim(posmna(1),posmna(2)) .lt. 0.d0) then
        reffg(:) = 1.d0 - reffg(:)
      endif

!     * maximum rsurfg depends on rsurf of trees in same layer.

      rsurfgnew(1:pos_slg) = MIN(2.d0 * c_epsln / i_rootrad            &
     &                     * s_delz(1:pos_slg) - rsurft_(1:pos_slg),   &
     &                     MAX(i_rsurfmin * s_delz(1:pos_slg),         &
     &                     rsurfg_(1:pos_slg) + rsurfg_(1:pos_slg)     &
     &                     * i_growthmax * rootlim(posmna(1),posmna(2))&
     &                     * reffg(1:pos_slg)))
      rsurfgnew(pos_slg+1:s_maxlayer) = 0.d0

      rootlim(:,:) = 0.d0
      ruptkt_d(:)  = 0.d0
      changef      = 0.d0

      return
      end subroutine vom_adapt_roots
