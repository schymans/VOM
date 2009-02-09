!********************************************************************
!*  Module defining variables and parameters for the vegetation model
!*  (transpmodel) and water balance model (watbal).
!*--------------------------------------------------------------------
!*  Author: Stan Schymanski, Max Planck Institute for Biogeochemistry
!*  Email: sschym@bgc-jena.mpg.de
!*  06/2008
!*--------------------------------------------------------------------
module vegwatbal
implicit none
save
 CHARACTER(60) dailyformat,hourlyformat
 CHARACTER(3) str
 !*
 !*=====<Parameters>===================================================
 !* set O2 conc [mol/mol], CO2 conc and diffusion conversion factor from
 !* CO2 to H2O as parameters. 
 !  
 REAL*8 a,pi,E,r,g,rho,degree,mpbar
 REAL*8, SAVE :: oa,ca,alpha,rlratio,k25,kopt,ha,hd,&
  tstart,toptfac,cpccf,tcf     ! cpccf=water transport costs per m root depth and m^2 cover; 
 INTEGER, SAVE :: parsaved,ny,nyt,firstyear,lastyear
 PARAMETER(a=1.6d0,pi=3.14159d0,mpbar=10.2d0,&      ! mpbar=conversion factor from MPa to bar
  g=9.81d0,rho=1000.d0,degree=0.0174533d0)   
 PARAMETER(E=2.7182818d0,R=8.314d0)
 INTEGER, SAVE :: N,Nh,M
 INTEGER d,h,oldh,i,in,ii,ik,yr,dtest,pos,posg,postemp,postempg,pos1(1),&
  pos2(2),stat,dummyint1,dummyint2,dummyint3,dummyint4
 REAL*8, ALLOCATABLE, SAVE :: srad(:),rhmax(:),rhmin(:),tmin(:),&
  tmax(:),rainvec(:),netassvec(:),epan(:),vpvec(:),&
  netassvecg(:),press(:)
 INTEGER, ALLOCATABLE, SAVE :: year(:),month(:),day(:),dayyear(:)
 REAL*8, ALLOCATABLE, SAVE :: avparvec(:),&
  parvec(:),parh(:),vdh(:),tairh(:),gammastarvec(:),rainh(:)
 REAL*8 jact(3),lambda,gstom,rl(3),vd,par,rain,transp,&
  ass(3),hass(3),pc,cpcc,sunr,suns,rr,sumwsu,cumerror,&    
  lambdafac,gammastar,jmax(3),jmax25(3),&
  netassyr,rainyr,epanyr,paryr,radyr,vdyr,etyr,evapyr,&
  gppyr,daylength,tamean,dtr,tair,vp,topt,&
  rootdepth,mq,mqssmin,dmq,mqnew,mqss,changef,&
  mqold,hruptk,error1,tc,mqssnew,pcg(3),lambdagfac,wsgexp,&
  cpccg(3),tcg(3),jmax25g(3),wsexp,rgdepth,&
  lambdag,rrg,jmaxg(3),rlg(3,3),jactg(3,3),gstomg(3,3),&
  etmg(3,3),transpg(3,3),assg(3,3),hassg(3,3),hetmg,&
  assg_d(3,3),jmaxg_d,gstomg_d,etmg_d,rlg_d,netassg_d(3,3),rootlim(3,3)
 REAL*8, SAVE :: md,mqx,rcond,rsurfmin,rsurfinit,rootrad, &
  rrootm,pcgmin,prootmg,mdf,mqxf,mdstore,growthmax
 REAL*8 ruptk_d, jmax_d,rl_d,gstom_d,ass_d(3),spgfcf_d,infx_d,etm_d,&
  esoil_d,vd_d,time,dtmq,dtmax

 CHARACTER(100) :: informat
 !*********************************************************************
 !*     Parameters for Water balance model
 !*********************************************************************
 !
 INTEGER nlayers,nlayersnew
 REAL*8, SAVE :: dt,cz,cgs,lat        ! lat=geogr. latitude
 REAL*8, SAVE :: ksat,epsln,nvg,mvg,alphavg,zs,zr,go,thetas,&
  thetar,wsold,delyu,ets
 REAL*8, ALLOCATABLE, SAVE :: pcapvec(:),suvec(:),ruptkvec(:),&
  sunewvec(:),wsvec(:),kunsatvec(:),delyuvec(:),rsurfvec(:),&
  rsurfnewvec(:),qblvec(:),dsuvec(:),phydrostaticvec(:),&
  delyunewvec(:),wsnewvec(:),prootmvec(:),pcapnewvec(:),&
  hruptkvec(:),ruptkvec_d(:),reff(:),ruptkg(:),rsurfg(:),&
  rsurfgnew(:),rsoilvec(:),kunsatnewvec(:)
 REAL*8 ys,yu,omgu,omgo,depth,etm,esoil,esoils,dys,dtys,dtomgu,dtyu,&
  domgu,dyu,dtsu,ysnew,yunew,wsnew,io,iocum,spgfcf,hspgfcf,infx,&
  hinfx,hio,hesoil,hetm,inf,omgunew,wc,wcnew,sumsutop,sumdsutop,&
  omgonew,error,hinf,dtss

end module vegwatbal
