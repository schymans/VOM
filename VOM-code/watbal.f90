!********************************************************************
!*  Layered water balance
!*--------------------------------------------------------------------
!*  Author: Stan Schymanski, Max Planck Institute for Biogeochemistry
!*  Email: sschym@bgc-jena.mpg.de
!*  06/2008
!*  Version: REW with layered unsaturated zone, no routing
!*--------------------------------------------------------------------
!*
!* Numbers in the commented parentheses refer to the equation numeration
!* in the document 'equations.pdf' that comes with the documentation.
!*
!*--------------------------------------------------------------------
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
!*--------------------------------------------------------------------
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
!********************************************************************

!subroutine waterbalance(time,dttarget,optmode)
subroutine waterbalance(init) 
 use vegwatbal
 implicit none

 REAL*8 dttarget,yutarget,dummy,timenew
 INTEGER(1), INTENT(in) :: init
 logical isnand
 logical isinfd 


 !* INITIALISATION
 if(init.eq.1) then   ! Run initial run 
  ysnew=zr              ! catchment filled up to channel bottom
  omgunew=1.d0              ! (3.1)
  omgonew=1.d0-omgunew
  yunew=(cz - ysnew)/omgunew          ! (3.2)
  nlayersnew=Ceiling(yunew/delyu)         ! (3.3), but the bottom layer <=nlayers
  delyunewvec(:)=delyu
  delyunewvec(nlayersnew)=yunew - (nlayersnew - 1.d0)*delyu
  depth=ysnew+delyunewvec(nlayersnew)/2.d0
  sunewvec(nlayersnew)=(1.d0 + (alphavg*(depth - ysnew))**nvg)**(-mvg) ! (Out[58])
  do i=nlayersnew-1,1,-1           
   depth=depth+(delyunewvec(i)+delyunewvec(i+1))/2.d0
   sunewvec(i)=(1.d0 + (alphavg*(depth - ysnew))**nvg)**(-mvg)  ! (Out[58])
  enddo
  sunewvec(nlayersnew+1:M)= 1.d0
  pcapnewvec(1:nlayersnew)= (-1.d0 + sunewvec(1:nlayersnew)**(-1.d0/mvg))**& ! (Out[54])
   (1.d0/nvg)/alphavg
  pcapnewvec(nlayersnew+1:M)=0.d0
  return
 endif


 !* CURRENT SOIL WATER CONTENT

 if(nlayers.ge.1) then
  wc=epsln*(ys+omgu*Sum(delyuvec(1:nlayers)*suvec(1:nlayers))) 
 else
  wc=epsln*ys
 endif

 if(nlayers.gt.1) then
  sumsutop=Sum(suvec(1:nlayers-1))
 else
  sumsutop=0.d0
 endif

 !* FLUXES

 if(rain.gt.0.d0) then
  if(nlayers.ge.1) then
   inf=Min((ksat+kunsatvec(1))/2.d0*omgu*(1.d0+(2.d0*&   ! (3.6), (Out[60]) 
    pcapvec(1))/delyuvec(1)),omgu*rain)
  else
   inf=0.d0
  endif
  infx=rain-inf
 else 
  inf=0.d0
  infx=0.d0
 endif
 if(nlayers.ge.1) then
  do i=1,nlayers-1
   qblvec(i)= -(omgu*(1.d0 + (-pcapvec(i) + pcapvec(i+1))/& ! (3.4), (Out[62]) 
    (0.5d0*delyuvec(i) + 0.5d0*delyuvec(i+1)))*&
    0.5d0*(kunsatvec(i) + kunsatvec(i+1)))
  enddo
  qblvec(nlayers)= -(omgu*(1.d0 + (-pcapvec(nlayers)/&   ! (3.5), (Out[62]) 
   (0.5d0*delyuvec(nlayers))))*0.5d0*&
   (kunsatvec(nlayers)+ksat))
  if(yu.eq.cz.and.qblvec(nlayers).gt.0.d0) then
   qblvec(nlayers)=0.d0    ! If the bottom boundary is bedrock, then no water can flow upwards.
  endif
  qblvec(nlayers+1:M)=0.d0
 else
  qblvec=0.d0
 endif
 esoil= 0.0002d0*(1.d0-0.8d0*(pc+pcg(2)))*par*suvec(1)*omgu  ! (3.9), (Out[141])
 esoils= 0.0002d0*(1.d0-0.8d0*(pc+pcg(2)))*par*omgo    ! (3.10), (Out[142], but suvec(1) instead of par removed)
 spgfcf= (0.5d0*ksat*omgo*(ys - zr))/(cgs*Cos(go))    ! (3.7), (Out[61])

 io=inf-esoil-esoils-spgfcf-Sum(ruptkvec)-Sum(ruptkg)  ! (3.19)


 !* MAKING SURE THAT NO SUBLAYER 'OVERFLOWS'
 !
 if(nlayers.ge.1) then
  if(Maxval(suvec(1:nlayers)).gt.0.999d0) then
   if(suvec(1).gt.0.999d0) then
    if(-esoil+inf+qblvec(1)-ruptkvec(1)-&
     ruptkg(1).gt.0.d0) then
     qblvec(1)=esoil - inf + ruptkvec(1)+ruptkg(1)  ! (Out[156]) +ruptkg(1)
    endif
   endif
   do i=2,nlayers
    if(suvec(i).gt.0.999d0) then
     if(-qblvec(i-1) + qblvec(i) - ruptkvec(i)-&
      ruptkg(i).gt.0.d0) then
      qblvec(i)=qblvec(i-1)+ruptkvec(i)+ruptkg(i)  ! (Out[158])+ruptkg(i)
     endif
    endif
   enddo
  endif
 endif

 !* CHANGES IN SOIL MOISTURE IN THE UNSATURATED ZONE
 !
 dtsu=99999.d0
 if(nlayers.ge.1) then
  if(nlayers.ge.2) then
   dsuvec(nlayers)=(qblvec(nlayers) - &
    qblvec(nlayers-1) -&   ! (3.17); (Out[151]) with added grass root uptake (ruptkg) to the equation
    ruptkvec(nlayers)-ruptkg(nlayers))/&
    (delyuvec(nlayers)*epsln*omgu)

   if (dsuvec(nlayers).gt.0.d0) then
    dtsu=Min(dtsu,0.9d0*(1.d0-suvec(nlayers))/dsuvec(nlayers),&
     suvec(nlayers)/(dsuvec(nlayers)*10.d0))
   elseif (dsuvec(nlayers).lt.0.d0) then
    dtsu=Min(dtsu,0.1d0*(-suvec(nlayers)/dsuvec(nlayers)))
   endif
   do i=2,nlayers-1            ! (3.17); (Out[150]) with added ruptkg(i) 
    dsuvec(i)= (qblvec(i) - qblvec(i-1) - ruptkvec(i)-&
     ruptkg(i))/(delyuvec(i)*epsln*omgu)
    if (dsuvec(i).gt.0.d0) then
     dtsu=Min(dtsu,0.9d0*(1.d0-suvec(i))/dsuvec(i),suvec(i)/&
      (dsuvec(i)*10.d0))
    elseif (dsuvec(i).lt.0.d0) then
     dtsu=Min(dtsu,0.1d0*(-suvec(i)/dsuvec(i)))
    endif


!!$!* for debugging in pgf90:
!!$    if (isnand(dsuvec(i))) then
!!$     print *, "Its a NaN"
!!$    elseif (isinfd(dsuvec(i))) then
!!$     print *, "Its a Inf" 
!!$    endif



   enddo
  endif
  dsuvec(1)=(-esoil + inf + qblvec(1) - ruptkvec(1)-&    ! (3.16), (Out[149]) with added ruptkg(1)
   ruptkg(1))/(delyuvec(1)*epsln*omgu)
  if (dsuvec(1).gt.0.d0) then
   dtsu=Min(dtsu,0.9d0*(1.d0-suvec(1))/dsuvec(1),&
    suvec(1)/(dsuvec(1)*10.d0))
  elseif (dsuvec(1).lt.0.d0) then
   dtsu=Min(dtsu,0.1d0*(-suvec(1)/dsuvec(1)))
  endif
  dsuvec(nlayers+1:M)=0.d0
 else
  dsuvec=0.d0
 endif

 if(nlayers.ge.2) then
  sumdsutop=Sum(dsuvec(1:nlayers-1))
  sumsutop=Sum(suvec(1:nlayers-1))
 else
  sumdsutop=0.d0
  sumsutop=0.d0
 endif


 !* LENGTH OF TIME STEP

 dt=Min(1.d0,dtsu,dtmax) !dt must be small enough

 if(ys.le.zr) then
  !* Testing for increase or decrease of yu:
  yutarget=(-(cz*epsln) + dt*io + wc - delyu*epsln*(dt*sumdsutop + sumsutop - &
   dt*(-1 + nlayers)*dsuvec(nlayers) + suvec(nlayers) - &
   nlayers*suvec(nlayers)))/(epsln*(-1 + dt*dsuvec(nlayers) + &
   suvec(nlayers)))
 else
  !* Making sure that Sqrt term is not negative  
  dummy=-1.d0 
  do while (dummy.lt.0.d0)
   dummy=(epsln*(4*(-(cz*epsln) + dt*io + wc)*(cz - zr)*(-1 + &
    dt*dsuvec(nlayers) + suvec(nlayers)) + delyu**2*epsln*(dt*sumdsutop + &
    sumsutop - dt*(-1 + nlayers)*dsuvec(nlayers) + suvec(nlayers) - &
    nlayers*suvec(nlayers))**2))
   dt=0.5d0*dt
  enddo
  yutarget=(delyu*epsln*(-(dt*sumdsutop) - sumsutop + dt*(-1 + &
   nlayers)*dsuvec(nlayers) - suvec(nlayers) + nlayers*suvec(nlayers)) - &
   Sqrt(epsln*(4*(-(cz*epsln) + dt*io + wc)*(cz - zr)*(-1 + &
   dt*dsuvec(nlayers) + suvec(nlayers)) + delyu**2*epsln*(dt*sumdsutop + &
   sumsutop - dt*(-1 + nlayers)*dsuvec(nlayers) + suvec(nlayers) - &
   nlayers*suvec(nlayers))**2)))/(2.*epsln*(-1 + dt*dsuvec(nlayers) + &
   suvec(nlayers)))
 endif


!!$ if(d.eq.62.and.h.eq.22.and.time.lt.50.d0) then
!!$  print*,"wait"
!!$ endif


 if(yutarget.gt.yu) then
  if(yu.eq.Ceiling(yu/delyu)*delyu) then
   yutarget=yu+delyu
  else
   yutarget=Ceiling(yu/delyu)*delyu
  endif
  if(yu.lt.cz-zr.and.yutarget.gt.cz-zr) then
   yutarget=cz-zr
  endif

 elseif(yutarget.lt.yu) then
  if(yu.eq.Floor(yu/delyu)*delyu) then
   yutarget=yu-delyu
  else
   yutarget=Floor(yu/delyu)*delyu
  endif
  if(yu.gt.cz-zr.and.yutarget.lt.cz-zr) then
   yutarget=cz-zr
  endif
 endif


 if(yutarget.gt.yu) then 
  if(ys.lt.zr) then
   if(yutarget.gt.cz) then
    yutarget=cz                  !* PREVENTING NEGATIVE YS
   endif
   dttarget=(-wc + epsln*(cz + yutarget*(-1.d0 + suvec(nlayers)) + delyu*(sumsutop + &
    suvec(nlayers) - nlayers*suvec(nlayers))))/(io - &
    epsln*yutarget*dsuvec(nlayers) + delyu*epsln*(-sumdsutop + (-1.d0 + &
    nlayers)*dsuvec(nlayers)))
  else
   dttarget=(cz**2.d0*epsln + wc*zr - cz*(wc + epsln*zr) + &
    epsln*yutarget*(yutarget*(-1.d0 + suvec(nlayers)) + delyu*(sumsutop + &
    suvec(nlayers) - nlayers*suvec(nlayers))))/(cz*io - io*zr - &
    epsln*yutarget*(yutarget*dsuvec(nlayers) + delyu*(sumdsutop + &
    dsuvec(nlayers) - nlayers*dsuvec(nlayers))))
  endif
 elseif(yutarget.lt.yu) then
  if(ys.gt.zr) then
   dttarget=(cz**2.d0*epsln + wc*zr - cz*(wc + epsln*zr) + &
    epsln*yutarget*(yutarget*(-1.d0 + suvec(nlayers)) + delyu*(sumsutop + &
    suvec(nlayers) - nlayers*suvec(nlayers))))/(cz*io - io*zr - &
    epsln*yutarget*(yutarget*dsuvec(nlayers) + delyu*(sumdsutop + &
    dsuvec(nlayers) - nlayers*dsuvec(nlayers))))
  else
   dttarget=(-wc + epsln*(cz + yutarget*(-1.d0 + suvec(nlayers)) + delyu*(sumsutop + &
    suvec(nlayers) - nlayers*suvec(nlayers))))/(io - &
    epsln*yutarget*dsuvec(nlayers) + delyu*epsln*(-sumdsutop + (-1.d0 + &
    nlayers)*dsuvec(nlayers)))
  endif

 else
  dttarget=99999.d0
 endif

 if(dttarget.le.0.d0) then
  dttarget=99999.d0
 endif

 dt=Max(0.d0,Min(dtsu,dttarget,dtmax))

 !*----- Calculating state variables at next time step----------------------- 
 !   

 if(ys.le.zr.and.yutarget.ge.cz-zr) then
  yunew=(-(cz*epsln) + dt*io + wc - delyu*epsln*(dt*sumdsutop + sumsutop - &
   dt*(-1 + nlayers)*dsuvec(nlayers) + suvec(nlayers) - &
   nlayers*suvec(nlayers)))/(epsln*(-1 + dt*dsuvec(nlayers) + &
   suvec(nlayers))) 
 else
  yunew=(delyu*epsln*(-(dt*sumdsutop) - sumsutop + dt*(-1 + &
   nlayers)*dsuvec(nlayers) - suvec(nlayers) + nlayers*suvec(nlayers)) - &
   Sqrt(epsln*(4*(-(cz*epsln) + dt*io + wc)*(cz - zr)*(-1 + &
   dt*dsuvec(nlayers) + suvec(nlayers)) + delyu**2*epsln*(dt*sumdsutop + &
   sumsutop - dt*(-1 + nlayers)*dsuvec(nlayers) + suvec(nlayers) - &
   nlayers*suvec(nlayers))**2)))/(2.*epsln*(-1 + dt*dsuvec(nlayers) + &
   suvec(nlayers)))
 endif




!!$if(yunew.gt.zr.and.yu.lt.zr) then
!!$Print*,"yunew>zr>yu"
!!$endif
!!$
!!$if(yunew.lt.zr.and.yu.gt.zr) then
!!$Print*,"yunew<zr<yu"
!!$endif
!!$
!!$if(Ceiling(yunew/delyu).ne.Ceiling(yu/delyu)) then
!!$Print*,"nlayers=",Ceiling(yu/delyu),"nlayersnew=",Ceiling(yunew/delyu)
!!$endif



 !*Rounding layers of less than 0.001 mm away from delyu
 if(Abs(yunew - nlayers*delyu).lt.1.e-9) then
  yunew=nlayers*delyu
 endif
 if(Abs(yunew - (nlayers-1)*delyu).lt.1.e-9) then
  yunew=(nlayers-1)*delyu
 endif

!* Calculating other state variables at next time step
 if(yunew.ge.cz-zr) then
  omgunew=1.d0
 else
  omgunew=yunew/(cz - zr)
 endif

 omgonew=1.d0-omgunew
 ysnew=cz-omgunew*yunew
 nlayersnew=Ceiling(yunew/delyu) ! The bottom layer is smaller than or equal to delyu
 delyunewvec(:)=delyu
 if(nlayersnew.ge.1) then 
  delyunewvec(nlayersnew)=yunew - (nlayersnew - 1)*delyu
 endif

 if(nlayers.ge.1) then
  sunewvec(1:nlayers)=suvec(1:nlayers)+dt*dsuvec(1:nlayers)
  pcapnewvec(1:nlayers)= (-1.d0 + sunewvec(1:nlayers)**(-1.d0/mvg))**& ! (Out[54])
   (1.d0/nvg)/alphavg
  sunewvec(nlayersnew+1:M)= 1.d0
  pcapvec(nlayersnew+1:M)=0.d0

 else
  sunewvec=1.d0
  pcapnewvec=0.d0
 endif

 !*ADAPTING SOIL MOISTURE TO NEW NUMBER OF LAYERS
 if(nlayersnew.gt.nlayers) then
  sunewvec(nlayersnew)=sunewvec(nlayers)
  pcapnewvec(nlayersnew)= (-1.d0 + sunewvec(nlayersnew)**(-1.d0/mvg))**& ! (Out[54])
   (1.d0/nvg)/alphavg
 endif



 kunsatnewvec=ksat*Sqrt(sunewvec)*(-1.d0+(1.d0-sunewvec**&    ! (3.14), (Out[55])
  (1.d0/mvg))**mvg)**2.d0



!!$  ! For debugging in pgf90: 
!!$   do i=1,nlayers
!!$    if (isnand(sunewvec(i))) then
!!$     print *, "Its a NaN"
!!$    elseif (isinfd(sunewvec(i))) then
!!$     print *, "Its a Inf" 
!!$    endif
!!$    enddo
!!$


 !* CHECK WATER BALANCE
 if(nlayersnew.ge.1) then
  wcnew=epsln*(ysnew+omgunew*Sum(delyunewvec(1:nlayersnew)*&
   sunewvec(1:nlayersnew))) 
 else
  wcnew=epsln*ysnew
 endif

!!$ print*,"errorstep=",(wc+dt*io-wcnew)
 if(Abs(wc+dt*io-wcnew).gt.1.d-6) then
  print*,"error=",(wc+dt*io-wcnew)," yu=",yu,"; yunew=",yunew
 endif


end subroutine waterbalance
