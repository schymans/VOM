!***********************************************************************
!        Optimised Vegetation Optimality Model (VOM)
!        Layered water balance
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
!        Version: 06/2008  REW with layered unsaturated zone, no routing
!-----------------------------------------------------------------------
!
!        Numbers in the commented parentheses refer to the equation numeration
!        in Schymanski (2007): PhD thesis, University of W.A.
!        and in the document 'equations.pdf' that comes with the documentation.
!
!        This subroutine ('waterbalance') uses the variables defined in
!        modules.f90, some of which have to be allocated prior to calling
!       'waterbalance'.
!        When calling the 'waterbalance', the 'init' variable must be given
!        (1 to generate initial conditions, otherwise 0).
!        This subroutine calculates the water balance for a time step <=dtmax.
!        The calling program must transfer the new variables ('new' in their names)
!        to the old variables for the next time step, before calling 'waterbalance'
!        again. Example: call waterbalance(init) , time=time+dt, zw=zwnew...
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

      subroutine waterbalance()
      use vom_vegwat_mod
      implicit none

      INTEGER :: jj

!     * CURRENT SOIL WATER CONTENT

      cH2Ol_s(:) = (-su__(:) * s_thetar(:) + su__(:) * s_thetas(:)     &
     &           + s_thetar(:)) * s_delz(:)
      wc_ = SUM(cH2Ol_s)

!     * FLUXES (inf, infx, qbl, esoil__, spgfcf__)

      call waterbalance_fluxes()

!     * changes in water storage (waterbalance)

      if( i_no_veg .eq. 0) then
         io__     = inf__ - esoil__ - spgfcf__ - SUM(ruptkt__(:)) - SUM(ruptkg__(:))  ! (3.19)
      else
         io__     = inf__ - esoil__ - spgfcf__  ! (no vegetation, WB still needs to close)
      end if

      iovec(:) = 0.d0
      iovec(1) = qbl(1) + inf__ - esoil__ - ruptkt__(1) - ruptkg__(1)
      if (wlayer_ .eq. 1) then
        iovec(1) = iovec(1) - spgfcf__
      endif
      if (wlayer_ .gt. 2) then
        do jj = 2, wlayer_ - 1
          iovec(jj) = qbl(jj) - qbl(jj-1) - ruptkt__(jj) - ruptkg__(jj)
        enddo
      endif
      if (wlayer_ .gt. 1) then
        iovec(wlayer_) = qbl(wlayer_) - qbl(wlayer_-1)                 &
     &                 - ruptkt__(wlayer_) - ruptkg__(wlayer_) - spgfcf__
      endif

!     * change in saturation degree

      dsu(:) = -iovec(:) / ((s_thetar(:) - s_thetas(:)) * s_delz(:))

!     * calculation of maximal time step size

      call waterbalance_timestep()

!     Calculating state variables at next time step
!     (sunew, wlayernew, pcapnew, kunsatnew, zwnew)

      call waterbalance_update_state()

      call waterbalance_diag()

      return
      end subroutine waterbalance

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----Run initial run--------------------------------------------------

      subroutine waterbalance_init ()
      use vom_vegwat_mod
      implicit none

      INTEGER :: jj

      dtsu_count  = 0
      dtmax_count = 0
      zw_         = i_cz
      wlayer_     = 0

!     * Adapt the number of unsaturated layers such that ys = zr

      do while (zw_ .gt. i_zr)
        wlayer_ = wlayer_ + 1
        zw_     = zw_ - s_delz(wlayer_)
      enddo
      pcap_(wlayer_+1:s_maxlayer) = 0.d0

!     * Equilibrium pressure head

      do jj = wlayer_, 1, -1
        pcap_(jj) = 0.5d0 * s_delz(jj+1) + 0.5d0 * s_delz(jj) + pcap_(jj+1)
      enddo

!     * equilibrium su above a saturated layer (after eq_sueq1 in Watbal3)

      do jj = 1, s_maxlayer - 1
        sueq(jj) = (1.d0 / ((0.5d0 * (s_delz(jj+1) +  s_delz(jj))      &
     &           * s_avg(jj)) ** s_nvg(jj) + 1.d0)) ** c_mvg(jj)
      enddo
      sueq(s_maxlayer) = (1.d0 / ((0.5d0 * s_delz(jj) * s_avg(jj))     &
     &                 ** s_nvg(jj) + 1.d0)) ** c_mvg(jj)
      su__(:) = (1.d0 / ((pcap_(:) * s_avg(:)) ** s_nvg(:) + 1.d0)) ** c_mvg(:)

!     * unsat. hydrol. cond. as a function of su

      kunsat_(:) = ((-su__(:) ** (1.d0 / c_mvg(:)) + 1.d0) ** c_mvg(:) &
     &           - 1.d0) ** 2.d0 * s_ksat(:) * SQRT(su__(:))
      cH2Ol_s(:) = (-su__(:) * s_thetar(:) + su__(:) * s_thetas(:)     &
     &           + s_thetar(:)) * s_delz(:)

      zwnew     = zw_
      wlayernew = wlayer_
      pcapnew   = pcap_
      sunew     = su__
      kunsatnew = kunsat_

      return
      end subroutine waterbalance_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----FLUXES (inf, infx, qbl, esoil__, spgfcf__)--------------------

      subroutine waterbalance_fluxes ()
      use vom_vegwat_mod
      implicit none

      REAL*8  :: dummy
      INTEGER :: ii, jj

!     * infiltration

      if (rain_h(th_) .gt. 0.d0) then
        if (wlayer_ .ge. 1) then
          inf__ = MIN((s_ksat(1) + kunsat_(1)) / 2.d0 *(1.d0           &
     &          + (2.d0 * pcap_(1)) / s_delz(1)), rain_h(th_))  ! (3.6), (Out[60])
        else
          inf__ = 0.d0
        endif
        infx__ = rain_h(th_) - inf__
      else
        inf__  = 0.d0
        infx__ = 0.d0
      endif

!     * unsaturated flow

      qbl(:) = 0.d0
      if (wlayer_ .gt. 1) then
!       * Runoff occurs only from the layer 'wlayer_',
!         therefore no downward flow into the layers below is allowed.
        do jj = 1, wlayer_ - 1
          qbl(jj) = -0.5d0 * (2.d0 * (pcap_(jj+1) - pcap_(jj)) / (s_delz(jj+1) &
     &            + s_delz(jj)) + 1.d0) * (kunsat_(jj+1) + kunsat_(jj))
        enddo
      endif

!     * soil evaporation

      esoil__  = (par_h(th_)/(srad2par_h * l_E_* rho_wat) ) * &
                  (1.d0 - (1.d0-i_trans_vegcov) * (o_cait + caig_d(2))) * su__(1)

!     * Seepage face flow as a function of zw_ following eq_spgfcf in Watbal3.

      spgfcf__ = 0.d0
      if (zw_ .gt. i_zr) then
        spgfcf__ = MAX(0.d0, 0.5d0 * (SQRT(i_cz - i_zr)                &
     &           - SQRT(i_cz - zw_)) * (zw_ - i_zr) * s_ksat(wlayer_)  &
     &           / (SQRT(i_cz - i_zr) * i_cgs * COS(i_go)))
      endif

!     * MAKING SURE THAT NO SUBLAYER 'OVERFLOWS'
!     * 1.d-16 makes sure that 0 does not get transformed to tiny positive

      if (MAXVAL(su__(1:wlayer_)) .ge. 1.d0) then

        if (wlayer_ .gt. 1) then
          if (su__(1) .ge. 0.99d0) then
            dummy = esoil__ - inf__ + ruptkt__(1) + ruptkg__(1)
            if (qbl(1) - dummy .gt. 0.d0) then
              qbl(1) = dummy - 1.d-16   ! (Out[156])+ruptkg__(1)
            endif
          endif
        endif

        if (wlayer_ .gt. 2) then
          do ii = 2, wlayer_ - 1
            if (su__(ii) .ge. 0.99d0) then
              dummy = qbl(ii-1) + ruptkt__(ii) + ruptkg__(ii)
              if (qbl(ii) - dummy .gt. 0.d0) then
                qbl(ii) = dummy - 1.d-16 ! (Out[158])+ruptkg__(ii)
              endif
            endif
          enddo
        endif

        if (wlayer_ .gt. 1) then
          if (su__(wlayer_) .ge. 1.d0) then
            dummy = -qbl(wlayer_-1) - ruptkt__(wlayer_) - ruptkg__(wlayer_)
!           * make sure that any surplus water runs off
            spgfcf__ = MAX(spgfcf__, dummy + 1.d-16)
          endif
        else
          if (su__(wlayer_) .ge. 1.d0) then
            dummy = inf__ - esoil__ - ruptkt__(1) - ruptkg__(1)
            spgfcf__ = MAX(spgfcf__, dummy + 1.d-16)
          endif
        endif

      endif

      return
      end subroutine waterbalance_fluxes

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----CALCULATION OF MAXIMAL TIME STEP SIZE----------------------------

      subroutine waterbalance_timestep ()
      use vom_vegwat_mod
      implicit none

      REAL*8  :: dtsu
      INTEGER :: jj

!!$    INTEGER :: ii
!!$    LOGICAL :: isnand, isinfd
!!$
!!$!*  for debugging in pgf90:
!!$    do ii = 1, wlayer_
!!$      if (isnand(dsu(ii))) then
!!$        print *, "Its a NaN"
!!$      elseif (isinfd(dsu(ii))) then
!!$        print *, "Its a Inf"
!!$      endif
!!$    enddo

      dtsu = 999999.d0
      do jj = 1, wlayer_
        if (dsu(jj) .lt. 0.d0) then
          dtsu = MIN(dtsu, -0.1d0 * su__(jj) / dsu(jj))
        endif
        if (dsu(jj) .gt. 0.d0) then
          dtsu = MIN(dtsu, 0.1d0 * su__(jj) / dsu(jj),                 &
     &               (1.d0 - su__(jj)) / dsu(jj))
        endif
        
        !dsu equal to zero
        if( abs(dsu(jj)) .lt. epsilon(dsu(jj)) ) then 
           dtsu = MIN(dtsu, 999999.d0)
        end if
        
      enddo

!     * LENGTH OF TIME STEP

      dt_ = MAX(0.d0, MIN(dtsu, dtmax))

      if( dt_ .le. epsilon( MIN(dtsu, dtmax)  )  ) then
         stop "Error dt=0.0s"      
      end if


      if (dt_ .eq. dtsu) then
        dtsu_count = dtsu_count + 1
      elseif (dt_ .eq. dtmax) then
        dtmax_count = dtmax_count + 1
      endif

      return
      end subroutine waterbalance_timestep

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----Calculating state variables at next time step--------------------

      subroutine waterbalance_update_state ()
      use vom_vegwat_mod
      implicit none

      INTEGER :: jj

      sunew(:) = su__(:) + dt_ * dsu(:)
      jj = s_maxlayer
      wlayernew = s_maxlayer

!     * Find the lowest unsaturated layer and set the one below it as wlayernew

      do while (sunew(jj) .gt. 0.999999d0 .and. jj .gt. 1)
        jj = jj - 1
        wlayernew = jj
      enddo

!     * If there is not enough moisture in the layer above the saturated
!       layers, the water table is in the top saturated layer.

      if (sunew(wlayernew) .lt. sueq(wlayernew)) then
        wlayernew = wlayernew + 1
      endif

      sunew(wlayernew+1:s_maxlayer) = 1.d0

      pcapnew(:) = 1.d0 / s_avg(:) * (sunew(:) ** (-1.d0 / c_mvg(:))   &
     &           - 1.d0) ** (1.d0 / s_nvg(:))

      kunsatnew(:) = s_ksat(:) * ((-sunew(:) ** (1.d0                  &
     &             / c_mvg(:)) + 1.d0) ** c_mvg(:)                     &
     &             - 1.d0) ** 2.d0 * SQRT(sunew(:))

!     * Position of the water table after eq_wt1 in watbal3:

      if (wlayernew .lt. s_maxlayer) then
        zwnew = SUM(s_delz(wlayernew:s_maxlayer)) - 2.d0               &
     &        * s_delz(wlayernew) * pcapnew(wlayernew)                 &
     &        / (s_delz(wlayernew+1) + s_delz(wlayernew))
      else
        zwnew = s_delz(wlayernew) - 2.d0 * s_delz(wlayernew)           &
     &        * pcapnew(wlayernew) / s_delz(wlayernew)
      endif

      return
      end subroutine waterbalance_update_state

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      subroutine waterbalance_diag ()
      use vom_vegwat_mod
      implicit none

      character(len=135) :: msg
      REAL*8  :: wcnew

!!$    INTEGER :: ii
!!$    LOGICAL :: isnand, isinfd
!!$
!!$!*  For debugging in pgf90:
!!$    do ii = 1, nlayers
!!$      if (isnand(sunew(ii))) then
!!$        print *, "Its a NaN"
!!$      elseif (isinfd(sunew(ii))) then
!!$        print *, "Its a Inf"
!!$      endif
!!$    enddo

!     * CHECK WATER BALANCE

      cH2Ol_s(:) = (-sunew(:) * s_thetar(:) + sunew(:) * s_thetas(:)   &
     &           + s_thetar(:)) * s_delz(:)

      wcnew = SUM(cH2Ol_s(:))

!!$    print*,"errorstep=",(wc+dt*io-wcnew)

      if (ABS(wc_ + dt_ * io__ - wcnew) .gt. 1.d-6) then
       ! write(msg,*) "error=", (wc_ + dt_ * io__ - wcnew), " ys=", zwnew
       ! write(*,*) TRIM(msg)
       ! write(msg,*) "sum(iovec) = ", SUM(iovec(:)), "; io = ", io__
       ! write(*,*) TRIM(msg)
       ! write(msg,*) "day = ", nday, "; hour = ", nhour
       ! write(*,*) TRIM(msg)



      endif

      return
      end subroutine waterbalance_diag
