!********************************************************************
!*     Transpiration model for given stomatal conductance
!*     and canopy efficiency Jmax and layered water balance
!*--------------------------------------------------------------------
!*     Author: Stan Schymanski, CWR, University of Western Australia
!*     03/2006
!*     Version: big leaf, trees and grass, layered unsaturated zone
!*		optimised root profile
!*     
!*     FOR COMPILATION WITH COMPAQ VISUAL FORTRAN
!*  
!*  Copyright (C) 2006  Stan Schymanski!*!*    This program is free software: you can redistribute it and/or modify!*    it under the terms of the GNU General Public License as published by!*    the Free Software Foundation, either version 3 of the License, or!*    (at your option) any later version.!*!*    This program is distributed in the hope that it will be useful,!*    but WITHOUT ANY WARRANTY; without even the implied warranty of!*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the!*    GNU General Public License for more details.!*!*    You should have received a copy of the GNU General Public License!*    along with this program.  If not, see <http://www.gnu.org/licenses/>.!*
!*
!********************************************************************
!
  subroutine transpmodel(invar,nrun,netass,option1)
	use DFLIB
	implicit none

	character(9) option1
	character(60) dailyformat,hourlyformat
	character(3) str
!*
!*=====<Parameters>===================================================
!* set O2 conc [mol/mol], CO2 conc and diffusion conversion factor from
!* CO2 to H2O as parameters. 
  
	INTEGER, PARAMETER :: double=SELECTED_REAL_KIND(13) 	! == REAL*8
	real(double) a,pi,E,r,g,rho,degree,mpbar

	real(double), save:: oa,ca,alpha,rlratio,k25,kopt,ha,hd,topt,&
						cpccf,tcf

	integer, save:: ny,nyt
	parameter(a=1.6,pi=3.14159,mpbar=10.2,&		!mpbar=conversion factor from MPa to bar
					g=9.81,rho=1000.,degree=0.0174533)		!cpccf=water transport costs per m root depth and m^2 cover; 
	parameter(E=2.7182818,R=8.314)
	integer, save:: N,Nh,M,optmode
	integer nrun,d,h,i,in,ii,ik,yr,dtest,pos,posg,postemp,pos1(1),&
				pos2(2)
	real(double), allocatable, save:: srad(:),rhmax(:),rhmin(:),tmin(:),&
				tmax(:),rainvec(:),netassvec(:),epan(:),vpvec(:),&
				netassvecg(:)
	integer, allocatable, save:: year(:),month(:),day(:),dayyear(:)
	real(double), allocatable, save:: avparvec(:),avvdvec(:),vdvec(:),&
		parvec(:),parh(:),vdh(:),tairh(:),gammastarvec(:)
	real(double) jact(3),lambda,gstom,rl(3),vd,par,rain,transp,&
     	ass(3),hass(3),pc,cpcc,sunr,suns,rr,&		!p is a parameter for calculating daylength
     	lambdafac,gammastar,jmax(3),jmax25(3),&
     	netassyr,rainyr,epanyr,paryr,radyr,vdyr,etyr,evapyr,&
		gppyr,daylength,tamean,dtr,tair,vp,&
		rootdepth,mq,mqssmin,time,dmq,mqnew,mqss,changef,&
		mqold,hruptk,error1,tc,mqssnew,pcg(3),lambdagfac,wsgexp,&
		cpccg(3),tcg(3),jmax25g(3),wsexp,rgdepth,&
		lambdag,rrg,jmaxg(3),rlg(3,3),jactg(3,3),gstomg(3,3),&
		etmg(3,3),transpg(3,3),assg(3,3),hassg(3,3),hetmg,&
		assg_d(3,3),jmaxg_d,gstomg_d,etmg_d,rlg_d,netassg_d(3,3),rootlim(3,3)
	
	real(double), save:: md,mqx,rcond,rsurfmin,rsurfinit,rootrad,	&
						rrootm,pcgmin,prootmg,mdf,mqxf,mdstore,growthmax
	real(double) ruptk_d,&
     		jmax_d,rl_d,gstom_d,ass_d(3),spgfcf_d,infx_d,etm_d,&
			esoil_d,vd_d
			
	real(double), INTENT(inout) :: netass
	real(double),  INTENT(in) :: invar(10)
	character(100) :: informat
!*********************************************************************
!*     Parameters for Water balance model
!*********************************************************************
	integer nlayers,nlayersnew
	real(double), save:: dt,cz,cgs,lat		!lat=latitude
	real(double), save:: ksat,epsln,nvg,mvg,alphavg,zs,zr,go,thetas,&
		thetar,wsold,delyu,ets
	real(double), allocatable, save:: pcapvec(:),suvec(:),ruptkvec(:),&
		sunewvec(:),wsvec(:),kunsatvec(:),delyuvec(:),rsurfvec(:),&
		rsurfnewvec(:),qblvec(:),dsuvec(:),phydrostaticvec(:),&
		delyunewvec(:),wsnewvec(:),prootmvec(:),pcapnewvec(:),&
		hruptkvec(:),ruptkvec_d(:),reff(:),ruptkg(:),rsurfg(:),&
		rsurfgnew(:),rsoilvec(:)
	real(double) ys,yu,omgu,omgo,depth,etm,esoil,esoils,dys,dtys,&
			domgu,dtsu,ysnew,yunew,wsnew,io,iocum,spgfcf,hspgfcf,infx,&
			hinfx,hio,hesoil,hrain,hetm,inf,dyu,omgunew,&
			omgonew,error,dtmq,hinf,dtss

!*********************************************************************
!*     Program options
!*********************************************************************
	if (option1.lt.'-optimise') then
		optmode=0
!		write(6,*) option1
	else
		optmode=1
!		write(6,*) option1
	endif
	
	ets=0.		!we assume that no water uptake in the saturated zone is happening
!*--------------------------------------------------------------------
!*
!*    PARAMETER READING FROM INPUT.PAR
!*
!*--------------------------------------------------------------------
	if (nrun.eq.1) then
		open(1,file='input6.par',status='old')
		read(1,*) oa
		read(1,*) ca
		read(1,*) alpha
		read(1,*) cpccf	
		read(1,*) tcf
		read(1,*) ny
		read(1,*) nyt

		read(1,*) k25
		read(1,*) kopt	
		read(1,*) ha
		read(1,*) hd
		read(1,*) topt
		read(1,*) rlratio

!*-----Catchment parameters --------------------------------------------

		read(1,*) lat
		read(1,*) cz
		read(1,*) zs
		read(1,*) cgs	
		read(1,*) zr
		read(1,*) go


!		lat=12.5							!geogr. latitude in degrees
!		cz=30.0										 !Capital z
!		zero=0.0
!		cgs=10.0      ! Capital Gamma S (length scale for seepage outflow REG) (m)
!		zr=cz-5.0 !/10000.0
!		go=0.033		! slope close to channel in radians (assumed 1/30 slope)

!*-----Soil parameters --------------------------------------------

		read(1,*) ksat	
		read(1,*) thetar
		read(1,*) thetas
		read(1,*) nvg
		read(1,*) alphavg


		epsln=thetas-thetar		! epsilon, porosity see Reggiani
		mvg=1-(1/nvg)		! van Genuchten soil parameter m
!*
!*-----Vertical Resolution--- -------------------------------------
		read(1,*) delyu

!*-----Vegetation Parameters--- -------------------------------------
		read(1,*) mdf
		read(1,*) mqxf
		read(1,*) rrootm
		read(1,*) rsurfmin
		read(1,*) rsurfinit
		read(1,*) rootrad
		read(1,*) prootmg
		read(1,*) growthmax
		rcond=1./rrootm				!root conductivity in s^-1 (ruptk=rcond(proot-psoil)rsurf
		close(1)
!*
!*
!*	End of Runoff model initialisation
!*********************************************************************
!*-----calculate vector sizes-----------------------------------------
		N=ceiling(ny*365.25)
		Nh=N*24
		M=ceiling(cz/delyu)				! maximum number of soil sublayers

!*-----allocate vector sizes------------------------------------------
		allocate(srad(N),rhmax(N),rhmin(N),tmin(N),tmax(N),rainvec(N))
		allocate(year(N),month(N),day(N),dayyear(N),vpvec(N),epan(N))
		allocate(avparvec(N),avvdvec(N))
		allocate(vdvec(N),parvec(N))
		allocate(netassvec(N),netassvecg(N))
		allocate(parh(Nh),vdh(Nh),tairh(Nh),gammastarvec(Nh))
		allocate(pcapvec(M),suvec(M),ruptkvec(M),sunewvec(M),wsvec(M),&
				kunsatvec(M),delyuvec(M),rsurfvec(M),rsurfnewvec(M),&
				qblvec(M),dsuvec(M),phydrostaticvec(M),delyunewvec(M),&
				wsnewvec(M),prootmvec(M),pcapnewvec(M),ruptkvec_d(M),&
				hruptkvec(M),reff(M),ruptkg(M),rsurfg(M),rsurfgnew(M),&
				rsoilvec(M))

!*--------------------------------------------------------------------
!*
!*     File opening
!*
!*--------------------------------------------------------------------
!-----saving climate and gstom ass data -----------------------------
		if (optmode.eq.0) then
			open(201,file='resultshourly.txt')
			close(201,status='delete')
			open(201,file='resultshourly.txt')
			write(201,101)'year','month','day','dcum','hour','rain','tmax',&
     			'tmin','epan','par','vd','esoil','inf','jmax25','mq'&
     			,'rl','lambda','lambdag','rr','ass','su','ys',&
     			'Ws','omgo','spgfcf','infx','het','yu','pc'

			open(202,file='etassdaily.txt')
			close(202,status='delete')
			open(202,file='etassdaily.txt')
			write(202,101)'year','month','day','dcum','hour','rain','tmax',&
     			'tmin','epan','esoil','par','vd','ruptk_d','jmax25','jmax'&
     			,'rl','lambda','gstom','rr','ass','su','ys',&
     			'Ws','omgo','spgfcf','infx','etm','sutop','pc'

			open(203,file='yearly.txt')
			close(203,status='delete')
			open(203,file='yearly.txt')
			write(203,'(a6,9a15)') "year","rainyr","epanyr","paryr","radyr","vdyr","etyr",&
				"esoilyr","netassyr","gppyr"

			open(204,file='rsurfdaily.txt')
			close(204,status='delete')
			open(204,file='rsurfdaily.txt')
			write(204,*)' year',' month',' day','   dcum','  rsurfsublayer'

			open(205,file='delyudaily.txt')
			close(205,status='delete')
			open(205,file='delyudaily.txt')
			write(205,*)' year',' month',' day','   dcum',' hour','  delyusublayer'

			open(206,file='ruptkhourly.txt')
			close(206,status='delete')
			open(206,file='ruptkhourly.txt')
			write(206,*)' year',' month',' day','   dcum',' hour','  delyusublayer'

			open(207,file='suvechourly.txt')
			close(207,status='delete')
			open(207,file='suvechourly.txt')
			write(207,*)' year',' month',' day','   dcum',' hour','  susublayer'

		endif		! only creates file if not in 'optimise' mode
!*
!*--------------------------------------------------------------------
!*     Climate and Calendar data reading
!*--------------------------------------------------------------------
		informat='(4i8,8f8.2)'
		open(1,file='datadrill.prn') 
		read(1,*)
		read(1,*)
		do i=1,N
			read(1,informat) dayyear(i),day(i),month(i),year(i),tmax(i),&
							tmin(i),rainvec(i),epan(i),srad(i),vpvec(i)
		enddo		    			
	
		close(  1)

!*********************************************************************
!*     Calculation of vegetation parameters
!*     
!*--------------------------------------------------------------------
!*     Equations in Fortrform1.nb
!*********************************************************************
!*--------------------------------------------------------------------
!*
!*	Calculation of derived parameters
!*
!*--------------------------------------------------------------------
!		rainvec(:)=0.6*rainvec(:)
		
		vdvec=0.006028127*2.718282**((17.27*tmax)/(237.3 +&
					tmax)) - 9.869233e-6*vpvec*100.	![mol/mol]
		parvec=2.0804*srad		! par in mol/m2 if srad was MJ/m2

		do in=1,N
			daylength=12. - 7.639437*ASin((0.397949*Cos(0.172142 + &
				0.017214*dayyear(in))*Tan(0.017453*lat))/&
				Sqrt(0.920818 - 0.079182*Cos(0.344284 + 0.034428*&
				dayyear(in))))
			sunr=12-0.5*daylength 	! sets time of sunrise and sunset
			suns=12+0.5*daylength	
			tamean=(tmax(in)+tmin(in))/2.
			dtr=tmax(in)-tmin(in)
			vp=vpvec(in)*100.	! vp in Pa

			do ik=1,24	!Loop through every hour of day, where ik=hour
				ii=in*24+ik-24
				tair=tamean + dtr*(0.0138*Cos(3.513 - ((-1 + ik)*Pi)/&
					3.) + 0.0168*Cos(0.822 - ((-1 + ik)*Pi)/4.) + &
					0.0984*Cos(0.36 - ((-1 + ik)*Pi)/6.) + &
					0.4632*Cos(3.805 - ((-1 + ik)*Pi)/12.))
				tairh(ii)=tair
				vd=0.006028127*2.718282**((17.27*tair)/(237.3 +&
					tair)) - 9.869233e-6*vp
				if(vd.le.0.)then
					vd=0.0
				endif
				gammastarvec(ii)=0.00004275*E ** ((18915.*(-25. + tair))/&
					(149.*R*(273. + tair)))
				vdh(ii)=vd
				
				if (sunr.le.ik .and. ik+1.le.suns) then
				parh(ii)=(-0.000873*parvec(in)*Cos(0.017453*lat)*&		! in mol/m2/s
					Sqrt(0.920818 - 0.079182*Cos(0.034428*&
					(10. + dayyear(in))))*Cos(0.2618*ik) - 0.000347*&
					parvec(in)*Cos(0.017214*&
					(10. + dayyear(in)))*Sin(0.017453*lat))/&
					(-1.250192*daylength*Cos(0.017214*(10. +&
					dayyear(in)))*Sin(0.017453*lat) + 24.*&
					Cos(0.017453*lat)*Sqrt(0.920818 - 0.079182*&
					Cos(0.034428*(10. + dayyear(in))))*&
					(1. - (0.158363*Cos(0.017214*(10. + dayyear(in)))**2*&
					Tan(0.017453*lat)**2)/(0.920818 - 0.079182*&
					Cos(0.034428*(10. + dayyear(in)))))**0.5)
				else 
					parh(ii)=0.0
				endif
			enddo
		enddo
	endif

!*--------------------------------------------------------------------
!*
!*     Initial values 
!*
!*--------------------------------------------------------------------
	NETASS=0.0


	
!*--------------------------------------------------------------------
!*
!*     Set soil moisture and vegetation parameters to initial conditions 
!*
!*--------------------------------------------------------------------
	
	ysnew=zr				! catchment filled up to channel bottom
	omgunew=1.0
	omgonew=1.0-omgunew
	yunew=(cz - ysnew)/omgunew
	
	nlayersnew=Floor(yunew/delyu)
	delyunewvec(:)=delyu
	delyunewvec(nlayersnew)=yunew - (nlayersnew - 1)*delyu
	
	depth=ysnew+delyunewvec(nlayersnew)/2.
	sunewvec(nlayersnew)=(1 + (alphavg*(depth - ysnew))**nvg)**(-mvg)
	do i=nlayersnew-1,1,-1
		depth=depth+(delyunewvec(i)+delyunewvec(i+1))/2
		sunewvec(i)=(1 + (alphavg*(depth - ysnew))**nvg)**(-mvg)
	enddo
	sunewvec(nlayersnew+1:M)= 1.

	wsnewvec(1:nlayersnew)=sunewvec(1:nlayersnew)*epsln*omgunew*&
		delyunewvec(1:nlayersnew)
	wsold=ysnew*epsln + Sum(wsnewvec(1:nlayersnew))		!initial soil water storage
	wsnew=wsold

	ys=ysnew
	omgu=omgunew
	omgo=omgonew
	yu=yunew
	nlayers=nlayersnew
	delyuvec=delyunewvec
	suvec=sunewvec
	wsvec=wsnewvec

	pcapvec(1:nlayersnew)= (-1 + sunewvec(1:nlayersnew)**(-1/mvg))**&
				(1/nvg)/alphavg
	pcapvec(nlayersnew+1:M)=0.

	write(str,'(i3)') nlayers		!internal write to convert from number to string
	dailyformat='(i6,i6,i4,i7,'//str//'e14.6)'		! includes a column for each sublayer 
	hourlyformat='(i6,i6,i4,i7,i5,'//str//'e14.6)'		! includes a column for each sublayer

!*--------------------------------------------------------------------
!*     Optimised parameters reading from invar
!*--------------------------------------------------------------------

	lambdagfac=invar(1)
	wsgexp=invar(2)	
	lambdafac=invar(3)
	wsexp=invar(4)
	pc=invar(5)
	rootdepth=invar(6)
	mdstore=invar(7)
	rgdepth=invar(8)

!*--------------------------------------------------------------------
!*	Setting yearly, daily and hourly parameters
!*--------------------------------------------------------------------
	yr=year(1)
	netassyr=0.0	
	gppyr=0.0
	rainyr=0.0
	paryr=0.0
	radyr=0.0
	vdyr=0.0
	etyr=0.0
	epanyr=0.0
	evapyr=0.0

	ruptk_d=0.
	ruptkvec_d=0.
	vd_d=0.
	jmax_d=0.
	gstom_d=0.
	etm_d=0.
	esoil_d=0.
	ass_d=0.
	assg_d=0.
	spgfcf_d=0.
	infx_d=0.
	rl_d=0.

	hass=0.
	hassg=0.
	hspgfcf=0.
	hinfx=0.
	hio=0.
	hesoil=0.
	hrain=0.
	hetm=0.
	hetmg=0.
	hruptk=0.
	hruptkvec=0.
	hinf=0.

	netassvec=0.
	netassvecg=0.
	iocum=0.
!*--------------------------------------------------------------------
!*     Set vegetation parameters 
!*--------------------------------------------------------------------
	md=pc*mdf+mdstore
	mqx=md*mqxf
	mqnew=0.95*mqx				!initial wood water storage 
	mqold=mqnew
	rsurfnewvec(:)=0.
	pos=Ceiling(rootdepth/delyu)
	posg=Ceiling(rgdepth/delyu)
	rsurfgnew(1:posg)=rsurfinit*omgunew*delyunewvec(1:posg)
	rsurfgnew(posg+1:M)=0.
	if(posg.gt.nlayersnew) then
		rsurfgnew(nlayersnew+1:posg)=rsurfmin*&
			delyunewvec(nlayersnew+1:pos)*omgunew
	endif
	rsurfnewvec(1:pos)=rsurfinit*delyunewvec(1:pos)*omgunew			!root surface density (root surface area/soil volume) in each sublayer
	if(pos.gt.nlayersnew) then
		rsurfnewvec(nlayersnew+1:pos)=rsurfmin*&
			delyunewvec(nlayersnew+1:pos)*omgunew
	endif
	

	jmax25(2)=0.0003
	jmax25g(2)=0.0003

	pcgmin=0.02			! minimum grass pc; initial point for growth
	pcg(2)=Min(1.-pc,pcgmin)
	pcg= pcg(2)+(/-0.02,0.0,0.02/)		! vector with values varying by 1%
	pcg(3)=Min(Max(pcgmin,pcg(3)),1.0-pc)
	rootlim=0.

!*--------------------------------------------------------------------
!*     Direct costs
!*--------------------------------------------------------------------
!	cpcc=cpccf*pc*rootdepth+mdstore*2.45e-10		! costs of water distribution and storage
	tc=tcf*pc*2.5			! foliage tunrover costs, assuming crown LAI of 2.5
!*--------------------------------------------------------------------
!*    Daily loops
!*--------------------------------------------------------------------
	d=0
	dtest=nyt*365
902	do while (d.lt.dtest)
		d=d+1
		rsurfvec=rsurfnewvec
		rsurfg=rsurfgnew
!		avpar=avparvec(d)
!		avvd=avvdvec(d)
!		jmax25=Max(0.,jmaxfac*avpar**parexp*avvd**vdexp*&
!     					wsnew**wsfac)

		lambda=lambdafac*(Sum(pcapvec(1:pos))/pos)**wsexp
		lambdag=lambdagfac*pcapvec(1)**wsgexp

		jmax25= jmax25(2)*(/0.99,1.0,1.01/)		! vector with values varying by 1%
		jmax25=Max(jmax25,50.0e-6)	! 14.7 micromol/m2/s is the minimum observed by Niinemets et al. (1999)
		jmax25g= jmax25g(2)*(/0.99,1.0,1.01/)
		jmax25g=Max(jmax25g,50.0e-6)

		pcg= pcg(2)+(/-0.02,0.0,0.02/)		! vector with values varying by 1%
		pcg= Max(pcg,0.)
		pcg(3)=Min(Max(pcgmin,pcg(3)),1.0-pc)
!		pcg=Max(pcg,pcgmin)
!		mdg=pcg(2)*3.0*100.0		! dry mass of grasses, assuming leaf area of 3m2 per m2 pcg 100g/m^2 as for deciduous trees in Eamus98
!		mqxg=mdg
!		cpccg=cpccf*pcg*delyu	! costs of water distribution in grasses
		tcg=tcf*pcg*2.5		! foliage turnover costs, assuming LAI/pc of 2.5
		
		rr=2.55e-7*Sum(rsurfvec(1:pos))			! root respiration [mol/s]

		if(pos.gt.nlayersnew) then
			cpcc=cpccf*pc*rootdepth+mdstore*2.45e-10		! costs of water distribution and storage
		else
			cpcc=cpccf*pc*Sum(delyunewvec(1:pos))+mdstore*2.45e-10
		endif

		if(nlayersnew.lt.posg) then
			cpccg=cpccf*pcg*rgdepth
		else
			cpccg=cpccf*pcg*Sum(delyunewvec(1:posg))
		endif

		rrg=2.55e-7*Sum(rsurfg(1:posg))	! root respiration grasses [mol/s] 
		rain=rainvec(d)/(24.0*3600.0*1000.0) 	! in m/s
		
		mqssmin=mqx			! resetting the minimum steady-state tissue water content to its maximum value
!*--------------------------------------------------------------------
!*    Hourly loops
!*--------------------------------------------------------------------	
	do h=1,24		! loops through each hour of daily dataset

		ii=d*24+h-24
		gammastar=gammastarvec(ii)
		tair=tairh(ii)
		vd=vdh(ii)
		par=parh(ii)
		jmax=(E**((ha*(-25. + tair)*(-273. + topt + 273.*r*topt))/&
            ((25. + 273.*r*topt)*(tair + 273.*r*topt)))*&
         ((-1. + E**(-(hd*(-298. + topt))/(25. + 273.*r*topt)))*ha +& 
           hd)*jmax25)/&
		 ((-1. + E**((hd*(273. + tair - topt))/(tair + 273.*r*topt)))*&
          ha + hd)
		rl=((ca - gammastar)*pc*jmax*rlratio)/&
			(4.*(ca + 2.*gammastar)*(1. + rlratio))

		jmaxg=(E**((ha*(-25. + tair)*(-273. + topt + 273.*r*topt))/&
            ((25. + 273.*r*topt)*(tair + 273.*r*topt)))*&
         ((-1. + E**(-(hd*(-298. + topt))/(25. + 273.*r*topt)))*ha +& 
           hd)*jmax25g)/&
		 ((-1. + E**((hd*(273. + tair - topt))/(tair + 273.*r*topt)))*&
          ha + hd)
		rlg(1,:)=((ca - gammastar)*pcg(1)*jmaxg*rlratio)/&
			(4.*(ca + 2.*gammastar)*(1. + rlratio))
		rlg(2,:)=((ca - gammastar)*pcg(2)*jmaxg*rlratio)/&
			(4.*(ca + 2.*gammastar)*(1. + rlratio))
		rlg(3,:)=((ca - gammastar)*pcg(3)*jmaxg*rlratio)/&
			(4.*(ca + 2.*gammastar)*(1. + rlratio))



!*-----calculate gstom, et and ass --------------------------------------- 
		if (par.gt.0.0)then
			jact=(1. - E**(-(alpha*par)/jmax))*jmax*pc
			jactg(1,:)=(1. - E**(-(alpha*par)/jmaxg))*jmaxg*pcg(1)
			jactg(2,:)=(1. - E**(-(alpha*par)/jmaxg))*jmaxg*pcg(2)
			jactg(3,:)=(1. - E**(-(alpha*par)/jmaxg))*jmaxg*pcg(3)

			if (vd.gt.0.0.and.lambda.gt.(2.*a*vd)/(ca + 2.*gammastar)&
				.and.jact(2).gt.(4*ca*rl(2) + 8*gammastar*rl(2))/&
				(ca - gammastar)) then

				gstom=Max(0.,(0.25*(a*(ca*(jact(2) - 4.*rl(2)) - 4.*&
					gammastar*(jact(2) + 2.*rl(2)))*vd*(ca*lambda + 2.*&
					gammastar*lambda - a*vd) + 1.7320508075688772*&
					Sqrt(a*gammastar*jact(2)*(ca*(jact(2) - 4.*rl(2)) - &
					gammastar*(jact(2) + 8.*rl(2)))*&
				   vd*(ca*lambda + 2.*gammastar*lambda - 2.*a*vd) ** &
				   2.*(ca*lambda + 2.*gammastar*lambda - a*vd))))/&
				   (a*(ca +	2.*gammastar) ** 2*vd*(ca*lambda + 2.*&
				   gammastar*lambda - a*vd)))
						
			else
				gstom=0.0
			endif

			transp=a*vd*gstom				!transpiration rate in mol/s
			etm=(transp*18.)/(10.0**6.0)	!transpiration rate in m/s

			where (vd.gt.0.0.and.lambdag.gt.(2.*a*vd)/(ca + 2.*gammastar)&
				.and.jactg.gt.(4*ca*rlg + &
				8*gammastar*rlg)/(ca - gammastar)) 

				gstomg=Max(0.,(0.25*(a*(ca*(jactg-4.*rlg)-4.*&
				  gammastar*(jactg+2.*rlg))*vd*(ca*lambdag+2.*&
				  gammastar*lambdag - a*vd) + 1.7320508075688772*&
				  Sqrt(a*gammastar*jactg*(ca*(jactg-4.*&
				  rlg)-gammastar*(jactg + 8.*rlg))*&
				  vd*(ca*lambdag + 2.*gammastar*lambdag - 2.*a*vd) ** &
				  2.*(ca*lambdag + 2.*gammastar*lambdag - a*vd))))/&
				  (a*(ca +	2.*gammastar) ** 2*vd*(ca*lambdag + 2.*&
				  gammastar*lambdag - a*vd)))
			
			elsewhere
				gstomg=0.0
			endwhere

			transpg=a*vd*gstomg				!transpiration rate in mol/s
			etmg=(transpg*18.)/(10.0**6.0)	!transpiration rate in m/s

		else
			par=0.0
			jact=0.0
			gstom=0.0
			etm=0.0
			jactg=0.0
			gstomg=0.0
			etmg=0.0
		endif
		
!*-----SUB-HOURLY LOOPS --------------------------------------- 
		time=0.
		hass=0.		!hourly assimilation
		hassg=0.
		
		do while (time.lt.3600.)
	!*********************************************************************
	!*	Integrated Multi-layer Soil and Vegetation Water Balance Model
	!*--------------------------------------------------------------------
	!*	Model derived from waterbalanceHS7.nb
	!*********************************************************************
	!*----- setting variables from previous loop---------------------------
			rsurfvec(1:nlayersnew)=rsurfvec(1:nlayersnew)/(omgu*&
				delyuvec(1:nlayersnew))*omgunew*&
				delyunewvec(1:nlayersnew)
			if(nlayersnew.lt.pos) then
				rsurfvec(nlayersnew+1:pos)=rsurfmin*delyu
			endif
			rsurfg(1:nlayersnew)=rsurfg(1:nlayersnew)/(omgu*&
				delyuvec(1:nlayersnew))*omgunew*&
				delyunewvec(1:nlayersnew)
			if(nlayersnew.lt.posg) then
				rsurfvec(nlayersnew+1:pos)=rsurfmin*delyu
			endif
			mq=mqnew
			ys=ysnew
			suvec=sunewvec
			omgu=omgunew
			omgo=omgonew
			yu=yunew
			nlayers=nlayersnew
			delyuvec=delyunewvec
			wsvec=wsnewvec
			suvec=sunewvec

	!*-----soil capillary pressure, infiltration and runoff--------------- 

!			rsurfvec(1:pos)=rsurf*delyunewvec(1:pos)*omgunew
!			if (pos.le.nlayersnew) then
!				rsurfvec(pos)=rsurf*(rootdepth-delyu*(pos-1))*omgu
!			endif

			pcapvec(1:nlayers)= (-1 + suvec(1:nlayers)**(-1/mvg))**&
				(1/nvg)/alphavg
			pcapvec(nlayers+1:M)=0.
			kunsatvec=ksat*Sqrt(suvec)*(-1.+(1.-suvec**&
					(1./mvg))**mvg)**2.

			postemp=Min(pos,nlayers)
			phydrostaticvec(1:postemp)= ((/1:postemp:1/)-0.5)*delyu			
			prootmvec(1:postemp)= (mpbar*(-mq+mqx)*(750 - (750*mqx)/&
				(md + mqx) + (md + mqx)/mqx))/(md + mqx) - &
				phydrostaticvec(1:postemp)
			rsoilvec(1:postemp)=(Sqrt(Pi/2.)*Sqrt((rootrad*omgu*&
				delyuvec(1:postemp))/rsurfvec(1:postemp)))/&
				kunsatvec(1:postemp)

			ruptkvec(1:postemp)=((-pcapvec(1:postemp) + &
				prootmvec(1:postemp))*rsurfvec(1:postemp))/(rrootm +&
				rsoilvec(1:postemp))

			ruptkvec(postemp+1:M)=0.							

			if(Maxval(etmg).gt.0.) then
				ruptkg(1:posg)=Max(0.0,((-pcapvec(1:posg) + &			! root uptake by grasses can not be negative, as storage negligible
					(prootmg-phydrostaticvec(1:posg)))*rsurfg)/&
					(rrootm +(Sqrt(Pi/2.)*Sqrt(rootrad*omgu*&
					delyuvec(1:posg)/rsurfg))/kunsatvec(1:posg)))
				ruptkg(posg+1:M)=0.
				if(Sum(ruptkg).gt.0.) then
					where(etmg.gt.Sum(ruptkg)) 
						rootlim=1.0
						etmg=Sum(ruptkg)
						transpg=etmg*55555.555555555555	!mol/s=m/s*10^6 g/m/(18g/mol)
						gstomg=transpg/(a*vd)
					endwhere

					ruptkg(1:posg)=etmg(2,2)*(ruptkg(1:posg)/(Sum(ruptkg)))
				else
					ruptkg=0.0
					rootlim=0.0
					etmg=0.0
					transpg=0.0
					gstomg=0.0
				endif

				
			else
				ruptkg=0.0
			endif
	
			if (rain.gt.0.) then
				inf=Min((ksat+kunsatvec(1))/2.*omgu*(1+(2.*&
					pcapvec(1))/delyuvec(1)),omgu*rain)
				infx=rain-inf
			else
				inf=0.0
				infx=0.0
			endif

			do i=1,nlayers-1
!				qblvec(i)= -(omgu*(1. + (-pcapvec(i) + pcapvec(i+1))/&
!						(0.5*delyuvec(i) + 0.5*delyuvec(i+1)))*&
!						Max(kunsatvec(i), kunsatvec(i+1)))
				qblvec(i)= -(omgu*(1. + (-pcapvec(i) + pcapvec(i+1))/&
						(0.5*delyuvec(i) + 0.5*delyuvec(i+1)))*&
						0.5*(kunsatvec(i) + kunsatvec(i+1)))
			enddo
			qblvec(nlayers)= -(omgu*(1. + (-pcapvec(nlayers)/&
						(0.5*delyuvec(nlayers))))*0.5*&
						(kunsatvec(nlayers)+ksat))
					
			esoil= 0.0002*(1-0.8*(pc+pcg(2)))*par*suvec(1)*omgu
			esoils= 0.0002*(1-0.8*(pc+pcg(2)))*par*omgo
			spgfcf= (0.5*ksat*omgo*(ys - zr))/(cgs*Cos(go))
			
903			dys=(esoils + spgfcf + qblvec(nlayers))/(epsln*(-1. +&
				Sum(suvec(1:nlayers))/nlayers))

			dtys=99999
			if (dys.lt.0) then					!preventing ys from becoming negative or time step of being too large
				if (ys.le.0.1*delyu) then
					qblvec(nlayers)=0.0
					dys=0.
					dtys=99999.
				else
					dtys=Min(0.1*(-ys/dys),-0.01*delyu/dys)
				endif
			elseif (dys.gt.0) then
				dtys=0.01*delyu/dys
			endif
						
			if(ys.gt.zr) then
				domgu= -dys/(2.*Sqrt((cz - ys)*(cz - zr)))
				dyu= (dys*(-cz + zr))/(2.*Sqrt((cz - ys)*(cz - zr)))
				if(dys.le.0.) then
					dtys=Min(dtys,-(ys-zr)/dys)
				endif
			elseif(ys.lt.zr) then
				domgu=0
				dyu=-dys	
				if(dys.gt.0.) then
					dtys=Min(dtys,(zr-ys)/dys)
				endif
			else
				domgu=0
				dyu=-dys
			endif

!* MAKING SURE THAT NO SUBLAYER 'OVERFLOWS'
			if(Maxval(suvec(1:nlayers)).gt.0.999) then
				if(suvec(1).gt.0.999) then
					if(-esoil+inf+qblvec(1)-ruptkvec(1)-&
						ruptkg(1).gt.0.) then
						qblvec(1)=esoil - inf + ruptkvec(1)+ruptkg(1)
					endif
				endif
				do i=2,nlayers-1
					if(suvec(i).gt.0.999) then
						if(-qblvec(i-1) + qblvec(i) - ruptkvec(i)-&
								ruptkg(i).gt.0.) then
							qblvec(i)=qblvec(i-1)+ruptkvec(i)+ruptkg(i)
						endif
					endif
				enddo
!				if(suvec(nlayers).gt.0.999) then
!					if(qblvec(nlayers) - qblvec(nlayers-1) -&
!							ruptkvec(nlayers)-ruptkg(nlayers).gt.0.) then
!						qblvec(nlayers)=qblvec(nlayers-1) +&
!							ruptkvec(nlayers)+ruptkg(nlayers)
!						go to 903
!					endif
!				endif
			endif
	!*-----steady-state tissue water (mqss) ---------------------------------------------- 

			mqss= Max(0.9*mqx,(mqx*(mpbar*(md*md+752.*md*mqx+mqx*mqx)*&
				Sum((rsurfvec(1:postemp)/&
				(rrootm + rsoilvec(1:postemp)))) - (md + mqx)*&
				(md+mqx)*(etm - Sum(((-phydrostaticvec(1:postemp) -&
				pcapvec(1:postemp))*rsurfvec(1:postemp))/&
				(rrootm + rsoilvec(1:postemp))))))/&
				(mpbar*(md*md + 752.*md*mqx + mqx*mqx)*&
				Sum((rsurfvec(1:postemp)/(rrootm + &
				rsoilvec(1:postemp))))))

			mqssmin=Min(mqssmin,mqss)


	!*-----transpiration, gstom and tissue water ---------------------------------------------- 
			if (mq.le.0.9*mqx) then 			! makes sure that tissue water does not get below 0.9mqx
				if (etm.gt.0.9*Sum(ruptkvec(1:nlayers))) then
					if (Sum(ruptkvec(1:nlayers)).ge.0.) then
						etm=Sum(ruptkvec(1:nlayers))
						transp=etm*55555.555555555555	!mol/s=m/s*10^6 g/m/(18g/mol)
						gstom=transp/(a*vd)
					else
!						where(ruptkvec.lt.0) rsurfvec=rsurfmin*delyuvec*omgu
!						prootmvec(1:nlayers)= (mpbar*(-mq+mqx)*(750 - &
!							(750*mqx)/(md + mqx) + (md + mqx)/mqx))/&
!							(md + mqx) - phydrostaticvec(1:nlayers)
!						ruptkvec(1:nlayers)= (-pcapvec(1:nlayers) + &
!							prootmvec(1:nlayers))*rcond*rsurfvec(1:nlayers)
!						etm=0.
!						gstom=0.
!						if (Sum(ruptkvec).le.0.) then
							write(6,'(a20,i2,a1,i2,a1,i4)') &
								'vegetation	dies on: ',&
								day(d),'/',month(d),'/',year(d)
							netass=0.
							RETURN								!if tissues water depleted, but still loosing water -> death
!						endif
							
					endif
					mqss= Max(0.9*mqx,(mqx*(mpbar*(md*md+752.*md*mqx+&
						mqx*mqx)*Sum((rsurfvec(1:postemp)/&
						(rrootm + rsoilvec(1:postemp)))) - (md + mqx)*&
						(md+mqx)*(etm-Sum(((-phydrostaticvec(1:postemp)&
						-pcapvec(1:postemp))*rsurfvec(1:postemp))/&
						(rrootm + rsoilvec(1:postemp))))))/&
						(mpbar*(md*md + 752.*md*mqx + mqx*mqx)*&
						Sum((rsurfvec(1:postemp)/(rrootm + &
						rsoilvec(1:postemp))))))

					mqssmin=Min(mqssmin,mqss)
				endif
			endif				
			dmq=(Sum(ruptkvec(1:nlayers))-etm)*1.e6		!rate of change in tissue moisture content in g/
			dtmq=99999.
			if (dmq.gt.0.) then			!avoids mq from becoming larger than mqx or smaller than 0.9mqx
				dtmq=(mqx-mq)/dmq
			elseif (dmq.lt.0.) then
				dtmq=(0.9*mqx-mq)/dmq
			endif
	!*-----change in soil moisture ---------------------------------------------- 
			dtsu=99999.

			dsuvec(nlayers)=(qblvec(nlayers) - qblvec(nlayers-1) -&
				ruptkvec(nlayers)-ruptkg(nlayers))/&
				(delyuvec(nlayers)*epsln*omgu)
			if (dsuvec(nlayers).gt.0.) then
				dtsu=Min(dtsu,0.9*(1.-suvec(nlayers))/dsuvec(nlayers),&
					suvec(nlayers)/(dsuvec(nlayers)*10.))
			elseif (dsuvec(nlayers).lt.0.) then
				dtsu=Min(dtsu,0.1*(-suvec(nlayers)/dsuvec(nlayers)))
			endif
		
			do i=2,nlayers-1
				dsuvec(i)= (qblvec(i) - qblvec(i-1) - ruptkvec(i)-&
					ruptkg(i))/(delyuvec(i)*epsln*omgu)
				if (dsuvec(i).gt.0.) then
					dtsu=Min(dtsu,0.9*(1.-suvec(i))/dsuvec(i),suvec(i)/&
						(dsuvec(i)*10.))
				elseif (dsuvec(i).lt.0.) then
					dtsu=Min(dtsu,0.1*(-suvec(i)/dsuvec(i)))
				endif
			enddo

			dsuvec(1)=(-esoil + inf + qblvec(1) - ruptkvec(1)-&
				ruptkg(1))/(delyuvec(1)*epsln*omgu)
			if (dsuvec(1).gt.0.) then
				dtsu=Min(dtsu,0.9*(1.-suvec(1))/dsuvec(1),&
					suvec(1)/(dsuvec(1)*10.))
			elseif (dsuvec(1).lt.0) then
				dtsu=Min(dtsu,0.1*(-suvec(1)/dsuvec(1)))
			endif
			
			
	!*----- Calculating maximum time step -------------------------------------	
901			if(Abs(mq-mqss).gt.mqx/1.e6) then
				dtss=(mq-mqss)/(1.e6*(etm-Sum(ruptkvec(1:nlayers))))
				if(dtss.le.0.) then
					dtss=9999.
				endif
			else
				dtss=9999.
			endif
			dt=Min(dtss,dtmq,dtsu,dtys,3600.-time)
			if(dt.eq.dtss) then			
				pcapnewvec(1:nlayers)= (-1 + (suvec(1:nlayers)+&
					dsuvec(1:nlayers)*dt)**(-1/mvg))**(1/nvg)/alphavg
				pcapnewvec(nlayers+1:M)=0.
				mqssnew=Max(0.9*mqx,(mqx*(mpbar*(md*md+752.*md*mqx+&
					mqx*mqx)*Sum((rsurfvec(1:postemp)/&
					(rrootm + rsoilvec(1:postemp)))) - (md + mqx)*&
					(md+mqx)*(etm-Sum(((-phydrostaticvec(1:postemp)&
					-pcapvec(1:postemp))*rsurfvec(1:postemp))/&
					(rrootm + rsoilvec(1:postemp))))))/&
					(mpbar*(md*md + 752.*md*mqx + mqx*mqx)*&
					Sum((rsurfvec(1:postemp)/(rrootm + &
					rsoilvec(1:postemp))))))
				if(Abs(mqssnew-mqss).gt.mqx/1.e4) then
	!				mqss=0.9*mqx			
					mqss=mqss+0.5*(mqssnew-mqss)
					go to 901
				endif
			endif

	!*----- Calculating state variables at next time step-----------------------	
			
			time=time+dt
			mqnew=mq+dmq*dt
			ysnew=ys+dys*dt
			sunewvec(1:nlayers)=suvec(1:nlayers)+dt*dsuvec(1:nlayers)

			if(ysnew.gt.zr) then
				omgunew=(cz - ysnew)/Sqrt((cz - ysnew)*(cz - zr))
			else 
				omgunew=1.
			endif
			omgonew=1.-omgunew
			yunew=(cz-ysnew)/omgunew
			
			io=dt*(inf-esoil-esoils-spgfcf-Sum(ruptkvec(1:nlayers))-&
				Sum(ruptkg(1:nlayers)))

			nlayersnew=Floor(yunew/delyu)
			delyunewvec(:)=delyu
			delyunewvec(nlayersnew)=yunew - (nlayersnew - 1)*delyu

			wsnewvec(1:nlayersnew-1)= sunewvec(1:nlayersnew-1)*epsln*&
					omgunew*delyunewvec(1:nlayersnew-1)
			wsnewvec(nlayersnew)= Sum(wsvec(1:nlayers))+ys*epsln+io -&
					Sum(wsnewvec(1:nlayersnew-1))-ysnew*epsln
			
			sunewvec(nlayersnew)= wsnewvec(nlayersnew)/&
					(epsln*delyunewvec(nlayersnew)*omgunew)
			sunewvec(nlayersnew+1:M)= 1.

	!*----- adding up hourly fluxes------------------------------------
			ass= (4.*ca*gstom + 8.*gammastar*gstom + jact - 4.*rl - &
				Sqrt((-4.*ca*gstom + 8.*gammastar*gstom + jact - &
                  4.*rl)**2. + &
               16.*gammastar*gstom*(8.*ca*gstom + jact + 8.*rl)))/8.
			hass=hass+ass*dt

			assg= (4.*ca*gstomg + 8.*gammastar*gstomg+jactg - 4.*rlg-&
				Sqrt((-4.*ca*gstomg + 8.*gammastar*gstomg + jactg - &
                  4.*rlg)**2. + &
               16.*gammastar*gstomg*(8.*ca*gstomg + jactg + 8.*rlg)))/8.
			hassg=hassg+assg*dt

			hruptkvec=hruptkvec+ruptkvec*dt

			if (optmode.eq.0) then
				hspgfcf=hspgfcf+dt*spgfcf
				hinfx=hinfx+dt*infx
				hio=hio+io
				hesoil=hesoil+dt*(esoil+esoils)
				hrain=hrain+dt*rain
				hetm=hetm+dt*etm
				hetmg=hetmg+dt*etmg(2,2)	
				hruptk=hruptk+dt*Sum(ruptkvec(1:nlayers))
				hinf=hinf+inf*dt
			endif		
		
!*------END OF HOUR----------------------------------------------------
		enddo
		netass=netass+hass(2)-3600.0*(cpcc+rr+tc)+hassg(2,2)-3600.0*&
				(cpccg(2)+rrg+tcg(2))		! rl does not need to be included here as ass=-rl if j=0 (at night)
		ass_d=ass_d+hass	
		assg_d=assg_d+hassg	
		ruptkvec_d=ruptkvec_d+hruptkvec

		if (optmode.eq.0) then
			netassvec(d)=netassvec(d)+hass(2)-3600.0*(cpcc+rr+tc)
			netassvecg(d)=netassvecg(d)+hassg(2,2)-3600.0*&
				(cpccg(2)+rrg+tcg(2))		! rl does not need to be included here as ass=-rl if j=0 (at night)

!*----- summary------------------------------------------------------
			ruptk_d=ruptk_d+Sum(ruptkvec(1:nlayers))*3600.
			vd_d=vd_d+vd
			jmax_d=jmax_d+jmax(2)
			jmaxg_d=jmaxg_d+jmaxg(2)
			gstom_d=gstom_d+gstom
			gstomg_d=gstomg_d+gstomg(2,2)
			etm_d=etm_d+hetm
			etmg_d=etmg_d+hetmg
			esoil_d=esoil_d+hesoil
			spgfcf_d=spgfcf_d+hspgfcf
			infx_d=infx_d+hinfx
			rl_d=rl_d+rl(2)*3600.		! rl_d in mol/day
			rlg_d=rlg_d+rlg(2,2)*3600.

			write(str,'(i3)') nlayers		!internal write to convert from number to string
			dailyformat='(i6,i6,i4,i7,'//str//'e14.6)'		! includes a column for each sublayer 
			hourlyformat='(i6,i6,i4,i7,i5,'//str//'e14.6)'		! includes a column for each sublayer
			if(year(d).gt.2000.and.year(d).le.2005)then
				write(201,102)year(d),month(d),day(d),d,h,hrain,&
     				tmax(d),tmin(d),epan(d),par,vd,hesoil,hinf,&
					jmax25(2),mq,rl(2)+rlg(2,2),lambda,lambdag,rr+rrg,&
					hass(2)+hassg(2,2),suvec(1),ys,&
					wsnew,omgo,hspgfcf,hinfx,hetm+hetmg,yu,pc+pcg(2)
				write(205,hourlyformat) year(d),month(d), day(d),d,h,delyuvec(1:nlayers)
				write(206,hourlyformat) year(d),month(d), day(d),d,h,hruptkvec(1:nlayers)
				write(207,hourlyformat) year(d),month(d), day(d),d,h,suvec(1:nlayers)


			endif


!*----- check water balance -----------------------------------------
			iocum=iocum+hio
			wsnew=ysnew*epsln + Sum(wsnewvec(1:nlayersnew))	
			error=wsold+iocum-wsnew
			if(abs(error/wsold).gt.1.e-6) then ! gives an error message if accumulated error exceeds 10^-4 of ws
				write(6,*)'Error in water balance [%]:',error*100.,'in=',in,&
     					'io=',io,'wsold=',wsold,'wsnew=',wsnew
			endif
			error1=mqold+(hruptk-hetm)*1.e6-mqnew
			if(abs(error1/mqnew).gt. 1.e-6) then
				write(6,*)'Error in tree water balance [%]:',error1*100.,'mqold=',mqold,&
     					'mqnew=',mqnew,'hruptk=',hruptk,'hetm=',hetm
			endif
			mqold=mqnew
			hspgfcf=0.
			hinfx=0.
			hio=0.
			hesoil=0.
			hrain=0.
			hetm=0.
			hetmg=0.
			hruptk=0.
			hinf=0.
		endif
	
		hass=0.
		hassg=0.
		hruptkvec=0.

	enddo
!*------END OF DAY----------------------------------------------------
	if (optmode.eq.0) then
		write(202,102)year(d),month(d),day(d),d,h,rainvec(d),&
     		tmax(d),tmin(d),epan(d),esoil_d,parvec(d),vd_d/24.,&
			ruptk_d,jmax25(2),jmax25g(2),rl_d,lambda,lambdag,&
			rr*3600.*24.+rrg*3600*24.,ass_d(2)+assg_d(2,2),&
			Sum(suvec(1:nlayers))/nlayers,ys,wsnew,omgo,spgfcf_d,&
			infx_d,etm_d+etmg_d,suvec(1),pcg(2)

		write(204,dailyformat) year(d),month(d), day(d),d,rsurfvec(1:nlayers)

		if(year(d).eq.yr) then
			rainyr=rainyr+rainvec(d)	! in [mm]
			epanyr=epanyr+epan(d)	! epan originally in [mm]/day
			paryr=paryr+parvec(d)
			radyr=radyr+srad(d)		! srad originally in MJ/day
			vdyr=vdyr+vd_d/24.
			etyr=etyr+(etm_d+etmg_d)*1000.	! in[mm]
			evapyr=evapyr+esoil_d*1000.	! in [mm]
!			netassyr=netassyr+ass_d(2)-cpcc*3600.0*24.-rr*3600.0*24.
			netassyr=netassyr+ass_d(2)-(cpcc+rr)*3600.0*24.+&
				assg_d(2,2)	-(cpccg(2)+rrg)*3600.0*24.0
!			gppyr=gppyr+(ass_d(2)+rl_d)
			gppyr=gppyr+(ass_d(2)+rl_d)+assg_d(2,2)+rlg_d
		else
			write(203,'(i6,9e15.6)') yr,rainyr,epanyr,paryr,radyr,vdyr/&
							(dayyear(d)),etyr,evapyr,netassyr,gppyr
			yr=year(d)
			rainyr=rainvec(d)
			epanyr=epan(d)	! epan originally in [mm]/day
			paryr=parvec(d)
			radyr=srad(d)		! srad originally in MJ/day
			vdyr=vd_d/24.
			etyr=(etmg_d+etm_d)*1000.
			evapyr=esoil_d*1000.
			netassyr=ass_d(2)-(cpcc+rr)*3600.0*24.+assg_d(2,2)-&
				(cpccg(2)+rrg)*3600.0*24.0	
			gppyr=(ass_d(2)+rl_d)+assg_d(2,2)+rlg_d
		endif	
		ruptk_d=0.
		vd_d=0.
		jmax_d=0.
		jmaxg_d=0.
		gstom_d=0.
		gstomg_d=0.
		etm_d=0.
		etmg_d=0.
		esoil_d=0.
		spgfcf_d=0.
		infx_d=0.
		rl_d=0.		
		rlg_d=0.
	endif
!*------ADJUSTMENT OF JMAX25 and PC------------------------------------
		pos1=Maxloc(ass_d)		
		jmax25(2)=jmax25(pos1(1))
		ass_d=0.

		netassg_d(1,:)=assg_d(1,:)-3600.0*24.*(cpccg(1)+rrg+tcg(1))
		netassg_d(2,:)=assg_d(2,:)-3600.0*24.*(cpccg(2)+rrg+tcg(2))
		netassg_d(3,:)=assg_d(3,:)-3600.0*24.*(cpccg(3)+rrg+tcg(3))
		pos2=Maxloc(netassg_d)
		pcg(2)=Min(1.0-pc,pcg(pos2(1)))
		jmax25g(2)=jmax25g(pos2(2))
		assg_d=0.

!*------ADJUSTMENT OF ROOT SURFACE------------------------------------
!		rsurfnew=0.
		reff=0.
		changef= (0.95*mqx - mqssmin)/(0.05*mqx)
		reff(1:pos) = 0.5*ruptkvec_d(1:pos)/&
			rsurfvec(1:pos)/(Maxval(ruptkvec_d(1:pos)/rsurfvec(1:pos)))
		where(ruptkvec_d(1:pos).lt.0.) 
			reff=0.
		endwhere


		if(changef.lt.0.) then
			reff=1.-reff
		endif
		rsurfnewvec(1:pos)=Min(2*epsln/rootrad*delyu*omgu,&			!rsurf=(2*epsln/rootrad) if all pores filled by roots
			Max(rsurfmin,rsurfvec(1:pos)+growthmax*changef*&
			reff(1:pos)*omgu*delyu))
		where(rsurfvec.gt.1.) 
			rsurfnewvec=Min(2*epsln/rootrad*delyu*omgu,&
				Max(rsurfmin*delyu*omgu,rsurfvec(1:pos)+&
				rsurfvec*growthmax*changef*reff(1:pos)*omgu*delyuvec))
		endwhere

		rsurfgnew(1:posg)=Min(2*epsln/rootrad*delyuvec(1:posg)*&
			omgu-rsurfvec(1:posg),&
			Max(rsurfmin*delyuvec(1:posg)*omgu,rsurfg(1:posg)*&
			(0.9+0.2*rootlim(pos2(1),pos2(2)))))			! maximum rsurfg depends on rsurf of trees in same layer.

!		rsurfgnew(1:posg)=Minval(2*epsln/rootrad*delyuvec(1:posg)*&
!			omgu-rsurfvec(1:posg),&
!			Max(rsurfmin*delyuvec(1:posg)*omgu,rsurfg(1:posg)*&
!			(0.9+0.2*rootlim(pos2(1),pos2(2)))))			! maximum rsurfg depends on rsurf of trees in same layer.
		rsurfgnew(posg+1:M)=0.
		rootlim=0.
		ruptkvec_d=0.
	enddo
!*------END OF DAILY LOOPS----------------------------------------------
!	netass=SUM(netassvec)
!	netass=SUM(netassvecg)
!	netass=SUM(netassvec)+SUM(netassvecg)
	if (d.lt.N) then
		if (netass.le.0.) then
			netass=netass/nyt*ny		!estimates how bad the carbon loss 
										!would be instead of running through
										!the whole set
		else			
			dtest=N
			goto 902
		endif
	endif

101	format(a6,a7,a7,a7,a7,24a15)
102	format(i6,i7,i7,i7,i7,24e15.5)
103	format(e15.5)
	close (201)
	close (202)
	close (203)
	close (204)
	close (205)
	close (206)
	close (207)
	return
	end










