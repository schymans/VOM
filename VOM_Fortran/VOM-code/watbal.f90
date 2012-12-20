      subroutine waterbalance(init)
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

      use watmod
      use vegwatmod
      implicit none

      INTEGER, INTENT(in) :: init

      REAL*8  :: dummy
      LOGICAL :: isnand
      LOGICAL :: isinfd
      REAL*8  :: dtsu
      REAL*8  :: wc_
      REAL*8  :: wcnew
      INTEGER :: i, j

!     * INITIALISATION

!     * Run initial run
      if (init .eq. 1) then
        dtsu_count = 0
        dtmax_count = 0
        ys_ = cz
        nlayers_ = 0
!       * Adapt the number of unsaturated layers such that ys = zr
        do while (ys_ .gt. zr_)
          nlayers_ = nlayers_ + 1
          ys_ = ys_ - delzvec(nlayers_)
        enddo
        pcapvec(nlayers_+1:M___) = 0.d0
!       * Equilibrium pressure head
        do j = nlayers_, 1, -1
          pcapvec(j) = 0.5d0 * delzvec(j+1) + 0.5d0 * delzvec(j) + pcapvec(j+1)
        enddo
!       * equilibrium su above a saturated layer (after eq_sueq1 in Watbal3)
        do j = 1, M___ - 1
          sueqvec(j) = (1.d0 / ((0.5d0 * (delzvec(j+1) +  delzvec(j))  &
     &               * avgvec(j)) ** nvgvec(j) + 1.d0)) ** mvgvec(j)
        enddo
        sueqvec(M___) = (1.d0 / ((0.5d0 * delzvec(j) * avgvec(j))      &
     &                ** nvgvec(j) + 1.d0)) ** mvgvec(j)
        suvec_(:) = (1.d0 / ((pcapvec(:) * avgvec(:)) ** nvgvec(:)     &
     &            + 1.d0)) ** mvgvec(:)
!       * unsat. hydrol. cond. as a function of su
        kunsatvec(:) = ((-suvec_(:) ** (1.d0 / mvgvec(:)) + 1.d0)      &
     &               ** mvgvec(:) - 1.d0) ** 2.d0 * ksatvec(:)         &
     &               * SQRT(suvec_(:))
        cH2Ol_s(:) = (-suvec_(:) * thetarvec(:) + suvec_(:)            &
     &             * thetasvec(:) + thetarvec(:)) * delzvec(:)

        ysnew        = ys_
        nlayersnew   = nlayers_
        pcapnewvec   = pcapvec
        sunewvec     = suvec_
        kunsatnewvec = kunsatvec
!        omgu_        = 1.d0                  ! In this version, the hill slope is vertical, thus there is no saturated surface fraction
!        omgunew      = omgu_

        return
      endif

!     * CURRENT SOIL WATER CONTENT

      cH2Ol_s(:) = (-suvec_(:) * thetarvec(:) + suvec_(:)             &
     &           * thetasvec(:) + thetarvec(:)) * delzvec(:)
      wc_ = SUM(cH2Ol_s)

!     * FLUXES

      if (rain_ .gt. 0.d0) then
        if (nlayers_ .ge. 1) then
          inf_ = MIN((ksatvec(1) + kunsatvec(1)) / 2.d0 *(1.d0        &
     &         + (2.d0 * pcapvec(1)) / delzvec(1)), rain_)  ! (3.6), (Out[60])
        else
          inf_ = 0.d0
        endif
        infx__ = rain_ - inf_
      else 
        inf_ = 0.d0
        infx__ = 0.d0
      endif

      qblvec(:) = 0.d0
      if (nlayers_ .gt. 1) then
!       * Runoff occurs only from the layer 'nlayers_',
!         therefore no downward flow into the layers below is allowed.
        do j = 1, nlayers_ - 1
          qblvec(j) = -0.5d0 * (2.d0 * (pcapvec(j+1) - pcapvec(j))     &
     &              / (delzvec(j+1) + delzvec(j)) + 1.d0)              &
     &              * (kunsatvec(j+1) + kunsatvec(j))
        enddo
      endif

      esoil__  = 0.0002d0 * (1.d0 - 0.8d0 * (pc_ + pcg_(2))) * par_ * suvec_(1)
!     * Seepage face flow as a function of ys_ following eq_spgfcf in Watbal3.
      spgfcf__ = 0.d0
      if (ys_ .gt. zr_) then
        spgfcf__ = MAX(0.d0, 0.5d0 * (sqrt(cz - zr_) - sqrt(cz - ys_))  &
     &           * (ys_ - zr_) * ksatvec(nlayers_) / (sqrt(cz - zr_)    &
     &           * cgs * Cos(go_)))
      endif

!     * MAKING SURE THAT NO SUBLAYER 'OVERFLOWS'

      if (MAXVAL(suvec_(1:nlayers_)) .ge. 1.d0) then

        if (nlayers_ .gt. 1) then
          if (suvec_(1) .ge. 0.99d0) then
            dummy = esoil__ - inf_ + ruptkvec(1) + ruptkg(1)
            if (qblvec(1) - dummy .gt. 0.d0) then
              qblvec(1) = dummy - 1.d-16              ! (Out[156])+ruptkg(1)
            endif
          endif
        endif

        if (nlayers_ .gt. 2) then
          do i = 2, nlayers_ - 1
            if (suvec_(i) .ge. 0.99d0) then
              dummy = qblvec(i-1) + ruptkvec(i) + ruptkg(i)     ! 1.d-16 makes sure that 0 does not get transformed to tiny positive
              if (qblvec(i) - dummy .gt. 0.d0) then
                qblvec(i) = dummy - 1.d-16          ! (Out[158])+ruptkg(i)
              endif
            endif
          enddo
        endif

        if (nlayers_ .gt. 1) then
          if (suvec_(nlayers_) .ge. 1.d0) then
            dummy = -qblvec(nlayers_-1) - ruptkvec(nlayers_) - ruptkg(nlayers_) 
!           * make sure that any surplus water runs off
            spgfcf__ = max(spgfcf__, dummy + 1.d-16)
          endif
        else
          if (suvec_(nlayers_) .ge. 1.d0) then
            dummy = inf_ - esoil__ - ruptkvec(1) - ruptkg(1)
            spgfcf__ = max(spgfcf__, dummy + 1.d-16)
          endif
        endif

      endif

!     * CHANGES IN SOIL MOISTURE IN THE UNSATURATED ZONE

      io_      = inf_ - esoil__ - spgfcf__ - SUM(ruptkvec(:)) - SUM(ruptkg(:))  ! (3.19)
      iovec(:) = 0.d0
      iovec(1) = qblvec(1) + inf_ - esoil__ - ruptkvec(1) - ruptkg(1)
      if (nlayers_ .eq. 1) then
        iovec(1) = iovec(1) - spgfcf__
      endif
      if (nlayers_ .gt. 2) then
        do j = 2, nlayers_ - 1     
          iovec(j) = qblvec(j) - qblvec(j-1) - ruptkvec(j) - ruptkg(j)
        enddo
      endif
      if (nlayers_ .gt. 1) then
        iovec(nlayers_) = qblvec(nlayers_) - qblvec(nlayers_-1)        &
     &                  - ruptkvec(nlayers_) - ruptkg(nlayers_) - spgfcf__
      endif

      dsuvec(:) = -iovec(:) / ((thetarvec(:) - thetasvec(:)) * delzvec(:))

!     * CALCULATION OF MAXIMAL TIME STEP SIZE

      dtsu = 999999.d0
      do j = 1, nlayers_
        if (dsuvec(j) .lt. 0.d0) then
          dtsu = MIN(dtsu, -0.1d0 * suvec_(j) / dsuvec(j))
        endif
        if (dsuvec(j) .gt. 0.d0) then
          dtsu = MIN(dtsu, 0.1d0 * suvec_(j) / dsuvec(j),              &
     &              (1.d0 - suvec_(j)) / dsuvec(j))
        endif
      enddo

!!$!* for debugging in pgf90:
!!$    if (isnand(dsuvec(i))) then
!!$     print *, "Its a NaN"
!!$    elseif (isinfd(dsuvec(i))) then
!!$     print *, "Its a Inf" 
!!$    endif

!     * LENGTH OF TIME STEP

      dt_ = MAX(0.d0, MIN(dtsu, dtmax))
      if(dt_.eq.dtsu) then
        dtsu_count = dtsu_count + 1
      elseif(dt_.eq.dtmax) then
        dtmax_count = dtmax_count + 1
      endif

!*----- Calculating state variables at next time step-------------------

      sunewvec(:) = suvec_(:) + dt_ * dsuvec(:)
      j = M___
      nlayersnew = M___
!     * Find the lowest unsaturated layer and set the one below it as nlayersnew
      do while (sunewvec(j) .gt. 0.999999d0 .and. j .gt. 1)
        j = j - 1
        nlayersnew = j  
      enddo
!     * If there is not enough moisture in the layer above the saturated
!       layers, the water table is in the top saturated layer.
      if (sunewvec(nlayersnew) .lt. sueqvec(nlayersnew)) then
        nlayersnew = nlayersnew + 1
      endif
      sunewvec(nlayersnew + 1:M___) = 1.d0
      pcapnewvec(:) = 1.d0 / avgvec(:) * (sunewvec(:) ** (-1.d0        &
     &              / mvgvec(:)) - 1.d0) ** (1.d0 / nvgvec(:))

      kunsatnewvec(:) = ksatvec(:) * ((-sunewvec(:) ** (1.d0           &
     &                / mvgvec(:)) + 1.d0) ** mvgvec(:) - 1.d0)        &
     &                ** 2.d0 * SQRT(sunewvec(:))
!    Position of the water table after eq_wt1 in watbal3:
      if (nlayersnew .lt. M___) then
        ysnew = SUM(delzvec(nlayersnew:M___)) - 2.d0                   &
     &        * delzvec(nlayersnew) * pcapnewvec(nlayersnew)           &
     &        / (delzvec(nlayersnew+1) + delzvec(nlayersnew))
      else
        ysnew = delzvec(nlayersnew) - 2.d0 * delzvec(nlayersnew)       &
     &        * pcapnewvec(nlayersnew) / delzvec(nlayersnew)
      endif

      cH2Ol_s(:) = (-sunewvec(:) * thetarvec(:) + sunewvec(:)          &
     &           * thetasvec(:) + thetarvec(:)) * delzvec(:)

!!$  ! For debugging in pgf90: 
!!$   do i=1,nlayers
!!$    if (isnand(sunewvec(i))) then
!!$     print *, "Its a NaN"
!!$    elseif (isinfd(sunewvec(i))) then
!!$     print *, "Its a Inf" 
!!$    endif
!!$    enddo

!     * CHECK WATER BALANCE

      wcnew = SUM(cH2Ol_s(:))

!!$ print*,"errorstep=",(wc+dt*io-wcnew)

      if (ABS(wc_ + dt_ * io_ - wcnew) .gt. 1.d-6) then
        print *, "error=", (wc_ + dt_ * io_ - wcnew), " ys=", ysnew
        print *, "sum(iovec) = ", SUM(iovec(:)), "; io = ", io_
        print *, "day = ", d___, "; hour = ", h__
      endif

      return
      end subroutine waterbalance
