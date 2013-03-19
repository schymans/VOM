!***********************************************************************
!*  Layered water balance
!*----------------------------------------------------------------------
!*  Author: Stan Schymanski, Max Planck Institute for Biogeochemistry
!*  Email: sschym@bgc-jena.mpg.de
!*  06/2008
!*  Version: REW with layered unsaturated zone, no routing
!*----------------------------------------------------------------------
!*
!* Numbers in the commented parentheses refer to the equation numeration
!* in the document 'equations.pdf' that comes with the documentation.
!*
!*----------------------------------------------------------------------
!*
!* This subroutine ('waterbalance') uses the variables defined in
!* modules.f90, some of which have to be allocated prior to calling
!* 'waterbalance'.
!* When calling the 'waterbalance', the 'init' variable must be given
!* (1 to generate initial conditions, otherwise 0).
!* This subroutine calculates the water balance for a time step <=dtmax.
!* The calling program must transfer the new variables ('new' in their names)
!* to the old variables for the next time step, before calling 'waterbalance'
!* again. Example: call waterbalance(init) , time=time+dt, ys=ysnew...
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
!*    along with this program. If not, see http://www.gnu.org/licenses.
!*
!***********************************************************************

      subroutine waterbalance()
      use vom_vegwat_mod
      implicit none

      INTEGER :: j

!     * CURRENT SOIL WATER CONTENT

      cH2Ol_s(:) = (-suvec_(:) * thetarvec(:) + suvec_(:)              &
     &           * thetasvec(:) + thetarvec(:)) * delzvec(:)
      wc_ = SUM(cH2Ol_s)

!     * FLUXES (inf, infx, qblvec, esoil__, spgfcf__)

      call waterbalance_fluxes()

!     * changes in water storage (waterbalance)

      io__     = inf__ - esoil__ - spgfcf__ - SUM(ruptkvec(:)) - SUM(ruptkg(:))  ! (3.19)
      iovec(:) = 0.d0
      iovec(1) = qblvec(1) + inf__ - esoil__ - ruptkvec(1) - ruptkg(1)
      if (wlayer_ .eq. 1) then
        iovec(1) = iovec(1) - spgfcf__
      endif
      if (wlayer_ .gt. 2) then
        do j = 2, wlayer_ - 1
          iovec(j) = qblvec(j) - qblvec(j-1) - ruptkvec(j) - ruptkg(j)
        enddo
      endif
      if (wlayer_ .gt. 1) then
        iovec(wlayer_) = qblvec(wlayer_) - qblvec(wlayer_-1)           &
     &                 - ruptkvec(wlayer_) - ruptkg(wlayer_) - spgfcf__
      endif

!     * change in saturation degree

      dsuvec(:) = -iovec(:) / ((thetarvec(:) - thetasvec(:)) * delzvec(:))

!     * calculation of maximal time step size

      call waterbalance_timestep()

!     Calculating state variables at next time step
!     (sunewvec, wlayernew, pcapnewvec, kunsatnewvec, ysnew)

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

      INTEGER :: j

      dtsu_count = 0
      dtmax_count = 0
      ys_ = cz
      wlayer_ = 0

!     * Adapt the number of unsaturated layers such that ys = zr

      do while (ys_ .gt. zr_)
        wlayer_ = wlayer_ + 1
        ys_ = ys_ - delzvec(wlayer_)
      enddo
      pcapvec(wlayer_+1:maxlayer) = 0.d0

!     * Equilibrium pressure head

      do j = wlayer_, 1, -1
        pcapvec(j) = 0.5d0 * delzvec(j+1) + 0.5d0 * delzvec(j) + pcapvec(j+1)
      enddo

!     * equilibrium su above a saturated layer (after eq_sueq1 in Watbal3)

      do j = 1, maxlayer - 1
        sueqvec(j) = (1.d0 / ((0.5d0 * (delzvec(j+1) +  delzvec(j))    &
     &             * avgvec(j)) ** nvgvec(j) + 1.d0)) ** mvgvec(j)
      enddo
      sueqvec(maxlayer) = (1.d0 / ((0.5d0 * delzvec(j) * avgvec(j))    &
     &              ** nvgvec(j) + 1.d0)) ** mvgvec(j)
      suvec_(:) = (1.d0 / ((pcapvec(:) * avgvec(:)) ** nvgvec(:)       &
     &          + 1.d0)) ** mvgvec(:)

!     * unsat. hydrol. cond. as a function of su

      kunsatvec(:) = ((-suvec_(:) ** (1.d0 / mvgvec(:)) + 1.d0)        &
     &             ** mvgvec(:) - 1.d0) ** 2.d0 * ksatvec(:) * SQRT(suvec_(:))
      cH2Ol_s(:) = (-suvec_(:) * thetarvec(:) + suvec_(:)              &
     &           * thetasvec(:) + thetarvec(:)) * delzvec(:)

      ysnew        = ys_
      wlayernew    = wlayer_
      pcapnewvec   = pcapvec
      sunewvec     = suvec_
      kunsatnewvec = kunsatvec

      return
      end subroutine waterbalance_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!*-----FLUXES (inf, infx, qblvec, esoil__, spgfcf__)--------------------

      subroutine waterbalance_fluxes ()
      use vom_vegwat_mod
      implicit none

      REAL*8  :: dummy
      INTEGER :: i, j

!     * infiltration

      if (rain__ .gt. 0.d0) then
        if (wlayer_ .ge. 1) then
          inf__ = MIN((ksatvec(1) + kunsatvec(1)) / 2.d0 *(1.d0        &
     &          + (2.d0 * pcapvec(1)) / delzvec(1)), rain__)  ! (3.6), (Out[60])
        else
          inf__ = 0.d0
        endif
        infx__ = rain__ - inf__
      else
        inf__  = 0.d0
        infx__ = 0.d0
      endif

!     * unsaturated flow

      qblvec(:) = 0.d0
      if (wlayer_ .gt. 1) then
!       * Runoff occurs only from the layer 'wlayer_',
!         therefore no downward flow into the layers below is allowed.
        do j = 1, wlayer_ - 1
          qblvec(j) = -0.5d0 * (2.d0 * (pcapvec(j+1) - pcapvec(j))     &
     &              / (delzvec(j+1) + delzvec(j)) + 1.d0)              &
     &              * (kunsatvec(j+1) + kunsatvec(j))
        enddo
      endif

!     * soil evaporation

      esoil__  = 0.0002d0 * (1.d0 - 0.8d0 * (pc_ + pcg_(2))) * par__ * suvec_(1)

!     * Seepage face flow as a function of ys_ following eq_spgfcf in Watbal3.

      spgfcf__ = 0.d0
      if (ys_ .gt. zr_) then
        spgfcf__ = MAX(0.d0, 0.5d0 * (SQRT(cz - zr_) - SQRT(cz - ys_)) &
     &           * (ys_ - zr_) * ksatvec(wlayer_) / (SQRT(cz - zr_)    &
     &           * cgs * COS(go_)))
      endif

!     * MAKING SURE THAT NO SUBLAYER 'OVERFLOWS'
!     * 1.d-16 makes sure that 0 does not get transformed to tiny positive

      if (MAXVAL(suvec_(1:wlayer_)) .ge. 1.d0) then

        if (wlayer_ .gt. 1) then
          if (suvec_(1) .ge. 0.99d0) then
            dummy = esoil__ - inf__ + ruptkvec(1) + ruptkg(1)
            if (qblvec(1) - dummy .gt. 0.d0) then
              qblvec(1) = dummy - 1.d-16     ! (Out[156])+ruptkg(1)
            endif
          endif
        endif

        if (wlayer_ .gt. 2) then
          do i = 2, wlayer_ - 1
            if (suvec_(i) .ge. 0.99d0) then
              dummy = qblvec(i-1) + ruptkvec(i) + ruptkg(i)
              if (qblvec(i) - dummy .gt. 0.d0) then
                qblvec(i) = dummy - 1.d-16   ! (Out[158])+ruptkg(i)
              endif
            endif
          enddo
        endif

        if (wlayer_ .gt. 1) then
          if (suvec_(wlayer_) .ge. 1.d0) then
            dummy = -qblvec(wlayer_-1) - ruptkvec(wlayer_) - ruptkg(wlayer_)
!           * make sure that any surplus water runs off
            spgfcf__ = MAX(spgfcf__, dummy + 1.d-16)
          endif
        else
          if (suvec_(wlayer_) .ge. 1.d0) then
            dummy = inf__ - esoil__ - ruptkvec(1) - ruptkg(1)
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
      INTEGER :: i, j

!$    LOGICAL :: isnand, isinfd
!$
!$!*  for debugging in pgf90:
!$    do i = 1, wlayer_
!$      if (isnand(dsuvec(i))) then
!$        print *, "Its a NaN"
!$      elseif (isinfd(dsuvec(i))) then
!$        print *, "Its a Inf"
!$      endif
!$    enddo

      dtsu = 999999.d0
      do j = 1, wlayer_
        if (dsuvec(j) .lt. 0.d0) then
          dtsu = MIN(dtsu, -0.1d0 * suvec_(j) / dsuvec(j))
        endif
        if (dsuvec(j) .gt. 0.d0) then
          dtsu = MIN(dtsu, 0.1d0 * suvec_(j) / dsuvec(j),              &
     &               (1.d0 - suvec_(j)) / dsuvec(j))
        endif
      enddo

!     * LENGTH OF TIME STEP

      dt_ = MAX(0.d0, MIN(dtsu, dtmax))

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

      INTEGER :: j

      sunewvec(:) = suvec_(:) + dt_ * dsuvec(:)
      j = maxlayer
      wlayernew = maxlayer

!     * Find the lowest unsaturated layer and set the one below it as wlayernew

      do while (sunewvec(j) .gt. 0.999999d0 .and. j .gt. 1)
        j = j - 1
        wlayernew = j
      enddo

!     * If there is not enough moisture in the layer above the saturated
!       layers, the water table is in the top saturated layer.

      if (sunewvec(wlayernew) .lt. sueqvec(wlayernew)) then
        wlayernew = wlayernew + 1
      endif

      sunewvec(wlayernew+1:maxlayer) = 1.d0

      pcapnewvec(:) = 1.d0 / avgvec(:) * (sunewvec(:) ** (-1.d0        &
     &              / mvgvec(:)) - 1.d0) ** (1.d0 / nvgvec(:))

      kunsatnewvec(:) = ksatvec(:) * ((-sunewvec(:) ** (1.d0           &
     &                / mvgvec(:)) + 1.d0) ** mvgvec(:) - 1.d0)        &
     &                ** 2.d0 * SQRT(sunewvec(:))

!     * Position of the water table after eq_wt1 in watbal3:

      if (wlayernew .lt. maxlayer) then
        ysnew = SUM(delzvec(wlayernew:maxlayer)) - 2.d0                &
     &        * delzvec(wlayernew) * pcapnewvec(wlayernew)             &
     &        / (delzvec(wlayernew+1) + delzvec(wlayernew))
      else
        ysnew = delzvec(wlayernew) - 2.d0 * delzvec(wlayernew)         &
     &        * pcapnewvec(wlayernew) / delzvec(wlayernew)
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

!$    INTEGER :: i
!$    LOGICAL :: isnand, isinfd
!$
!$!*  For debugging in pgf90:
!$    do i = 1, nlayers
!$      if (isnand(sunewvec(i))) then
!$        print *, "Its a NaN"
!$      elseif (isinfd(sunewvec(i))) then
!$        print *, "Its a Inf"
!$      endif
!$    enddo

!     * CHECK WATER BALANCE

      cH2Ol_s(:) = (-sunewvec(:) * thetarvec(:) + sunewvec(:)          &
     &           * thetasvec(:) + thetarvec(:)) * delzvec(:)

      wcnew = SUM(cH2Ol_s(:))

!$    print*,"errorstep=",(wc+dt*io-wcnew)

      if (ABS(wc_ + dt_ * io__ - wcnew) .gt. 1.d-6) then
        write(msg,*) "error=", (wc_ + dt_ * io__ - wcnew), " ys=", ysnew
        write(*,*) TRIM(msg)
        write(msg,*) "sum(iovec) = ", SUM(iovec(:)), "; io = ", io__
        write(*,*) TRIM(msg)
        write(msg,*) "day = ", nday, "; hour = ", nhour
        write(*,*) TRIM(msg)
      endif

      return
      end subroutine waterbalance_diag
