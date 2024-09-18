 SUBROUTINE THERMAL (itypwat ,lb      ,nl_soil,dtime  ,trsmx0 ,&
                     zlnd    ,zsno    ,csoilc ,dewmx  ,capr   ,&
                     cnfac   ,csol    ,porsl  ,phi0   ,bsw    ,&
                     dkmg    ,dkdry   ,dksatu ,lai    ,sai    ,&
                     z0m     ,displa  ,sqrtdi ,rootfr ,effcon ,&
                     vmax25  ,slti    ,hlti   ,shti   ,hhti   ,&
                     trda    ,trdm    ,trop   ,gradm  ,binter ,&
                     extkn   ,hu      ,ht     ,hq     ,us     ,&
                     vs      ,tm      ,qm     ,rhoair ,psrf   ,&
                     pco2m   ,po2m    ,coszen ,parsun ,parsha ,&
                     sabvsun ,sabvsha ,sabg   ,frl    ,extkb  ,&
                     extkd   ,thermk  ,fsno   ,sigf   ,dz     ,&
                     z       ,zi      ,tlsun  ,tlsha  ,tss    ,&
                     wice    ,wliq    ,ldew   ,scv    ,snowdp ,&
                     imelt   ,taux    ,tauy   ,fsena  ,fevpa  ,&
                     lfevpa  ,fsenl   ,fevpl  ,etr    ,fseng  ,&
                     fevpg   ,olrg    ,fgrnd  ,rootr  ,qseva  ,&
                     qsdew   ,qsubl   ,qfros  ,sm     ,tref   ,&
                     qref    ,trad    ,rst    ,assim  ,respc  ,&
                     errore  ,emis    ,z0ma   ,zol    ,rib    ,&
                     ustar   ,qstar   ,tstar  ,u10m   ,v10m   ,&
                     f10m    ,fm      ,fh     ,fq     ,rstfac)
!=======================================================================
! this is the main subroutine to execute the calculation
! of thermal processes and surface fluxes
!
! Original author : Yongjiu Dai, 09/15/1999; 08/30/2002
!
! FLOW DIAGRAM FOR THERMAL.F90
!
! THERMAL ===> qsadv
!              groundfluxes
!              eroot                      |dewfraction
!              leaftemone |               |qsadv
!              leaftemtwo |  ---------->  |moninobukini
!                                         !moninobuk
!                                         |stomata
!
!              groundTem     ---------->   meltf
!
!=======================================================================
  use precision
  use phycon_module, only : denh2o,roverg,hvap,hsub,rgas,cpair,&
                            stefnc,denice,tfrz,vonkar,grav 
  implicit none
!---------------------Argument------------------------------------------
  integer, INTENT(in) :: &
        lb,          &! lower bound of array
        nl_soil,     &! upper bound of array
        itypwat       ! land water type (0=soil, 1=urban or built-up, 2=wetland,
! 3=land ice, 4=deep lake, 5=shallow lake)
  real(r8), INTENT(in) :: &
        dtime,       &! model time step [second]
        trsmx0,      &!max transpiration for moist soil+100% veg.  [mm/s]
        zlnd,        &!roughness length for soil [m]
        zsno,        &!roughness length for snow [m]
        csoilc,      &!drag coefficient for soil under canopy [-]
        dewmx,       &!maximum dew
        capr,        &!tuning factor to turn first layer T into surface T
        cnfac,       &!Crank Nicholson factor between 0 and 1
! soil physical parameters
        csol(1:nl_soil), &! heat capacity of soil solids [J/(m3 K)]
        porsl(1:nl_soil),&! soil porosity [-]
        phi0(1:nl_soil), &! soil water suction, negative potential [m]
        bsw(1:nl_soil),  &! clapp and hornbereger "b" parameter [-]
        dkmg(1:nl_soil), &! thermal conductivity of soil minerals [W/m-K]
        dkdry(1:nl_soil),&! thermal conductivity of dry soil [W/m-K]
        dksatu(1:nl_soil), &! thermal conductivity of saturated soil [W/m-K]
! vegetation parameters
        lai,         &! adjusted leaf area index for seasonal variation [-]
        sai,         &! stem area index  [-]
        z0m,         &! roughness length, momentum [m]
        displa,      &! displacement height [m]
        sqrtdi,      &! inverse sqrt of leaf dimension [m**-0.5]
        rootfr(1:nl_soil),&! root fraction
        effcon,      &! quantum efficiency of RuBP regeneration (mol CO2/mol quanta)
        vmax25,      &! maximum carboxylation rate at 25 C at canopy top
        slti,        &! slope of low temperature inhibition function      [s3]
        hlti,        &! 1/2 point of low temperature inhibition function  [s4]
        shti,        &! slope of high temperature inhibition function     [s1]
        hhti,        &! 1/2 point of high temperature inhibition function [s2]
        trda,        &! temperature coefficient in gs-a model             [s5]
        trdm,        &! temperature coefficient in gs-a model             [s6]
        trop,        &! temperature coefficient in gs-a model
        gradm,       &! conductance-photosynthesis slope parameter
        binter,      &! conductance-photosynthesis intercept
        extkn,       &! coefficient of leaf nitrogen allocation
! atmospherical variables and observational height
        hu,          &! observational height of wind [m]
        ht,          &! observational height of temperature [m]
        hq,          &! observational height of humidity [m]
        us,          &! wind component in eastward direction [m/s]
        vs,          &! wind component in northward direction [m/s]
        tm,          &! temperature at agcm reference height [kelvin]
        qm,          &! specific humidity at agcm reference height [kg/kg]
        rhoair,      &! density air [kg/m3]
        psrf,        &! atmosphere pressure at the surface [pa]
        pco2m,       &! CO2 concentration in atmos. (35 pa)
        po2m,        &! O2 concentration in atmos. (20900 pa)
! radiative fluxes
        coszen,      &! cosine of the solar zenith angle
        parsun,      &! photosynthetic active radiation by sunlit leaves (W m-2)
        parsha,      &! photosynthetic active radiation by shaded leaves (W m-2)
        sabvsun,     &! solar radiation absorbed by vegetation [W/m2]
        sabvsha,     &! solar radiation absorbed by vegetation [W/m2]
        sabg,        &! solar radiation absorbed by ground [W/m2]
        frl,         &! atmospheric infrared (longwave) radiation [W/m2]
        extkb,       &! (k, g(mu)/mu) direct solar extinction coefficient
        extkd,       &! diffuse and scattered diffuse PAR extinction coefficient
        thermk,      &! canopy gap fraction for tir radiation
! state variable (1)
        fsno,        &! fraction of ground covered by snow
        sigf,        &! fraction of veg cover, excluding snow-covered veg [-]
        dz(lb:nl_soil),  &! layer thickiness [m]
        z (lb:nl_soil),  &! node depth [m]
        zi(lb-1:nl_soil)  ! interface depth [m]
! state variables (2)
  real(r8), INTENT(inout) :: &
        tlsun,       &! sunlit leaf temperature [K]
        tlsha,       &! shaded leaf temperature [K]
        tss (lb:nl_soil),&! soil temperature [K]
        wice(lb:nl_soil),&! ice lens [kg/m2]
        wliq(lb:nl_soil),&! liqui water [kg/m2]
        ldew,        &! depth of water on foliage [kg/(m2 s)]
        scv,         &! snow cover, water equivalent [mm, kg/m2]
        snowdp        ! snow depth [m]
  integer, INTENT(out) :: & 
       imelt(lb:nl_soil)  ! flag for melting or freezing [-]
! Output fluxes
  real(r8), INTENT(out) :: &
        taux,        &! wind stress: E-W [kg/m/s**2]
        tauy,        &! wind stress: N-S [kg/m/s**2]
        fsena,       &! sensible heat from canopy height to atmosphere [W/m2]
        fevpa,       &! evapotranspiration from canopy height to atmosphere [mm/s]
        lfevpa,      &! latent heat flux from canopy height to atmosphere [W/m2]
        fsenl,       &! ensible heat from leaves [W/m2]
        fevpl,       &! evaporation+transpiration from leaves [mm/s]
        etr,         &! transpiration rate [mm/s]
        fseng,       &! sensible heat flux from ground [W/m2]
        fevpg,       &! evaporation heat flux from ground [mm/s]
        olrg,        &! outgoing long-wave radiation from ground+canopy
        fgrnd,       &! ground heat flux [W/m2]
	rootr(1:nl_soil),&! root resistance of a layer, all layers add to 1
        qseva,       &! ground surface evaporation rate (mm h2o/s)
        qsdew,       &! ground surface dew formation (mm h2o /s) [+]
        qsubl,       &! sublimation rate from snow pack (mm h2o /s) [+]
        qfros,       &! surface dew added to snow pack (mm h2o /s) [+]
        sm,          &! rate of snowmelt [kg/(m2 s)]
        tref,        &! 2 m height air temperature [kelvin]
        qref,        &! 2 m height air specific humidity
        trad,        &! radiative temperature [K]
        rst,         &! stomatal resistance (s m-1)
       assim,       &! assimilation
        respc,       &! respiration
! additional variables required by coupling with WRF or RSM model
        emis,        &! averaged bulk surface emissivity
        z0ma,        &! effective roughness [m]
        zol,         &! dimensionless height (z/L) used in Monin-Obukhov theory
        rib,         &! bulk Richardson number in surface layer
        ustar,       &! u* in similarity theory [m/s]
        qstar,       &! q* in similarity theory [kg/kg]
        tstar,       &! t* in similarity theory [K]
        u10m,        &! 10m u-velocity
        v10m,        &! 10m v-velocity
        f10m,        &! integral of profile function for momentum at 10m
        fm,          &! integral of profile function for momentum
        fh,          &! integral of profile function for heat
        fq,          &! integral of profile function for moisture
        rstfac       !factor of soil water stress
!---------------------Local Variables-----------------------------------
  integer i,j
  real(r8) :: &
       cgrnd,        &! deriv. of soil energy flux wrt to soil temp [w/m2/k]
       cgrndl,       &! deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
       cgrnds,       &! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
       degdT,        &! d(eg)/dT
       dqgdT,        &! d(qg)/dT
       dlrad,        &! downward longwave radiation blow the canopy [W/m2]
       eg,           &! water vapor pressure at temperature T [pa]
       egsmax,       &! max. evaporation which soil can provide at one time step
       egidif,       &! the excess of evaporation over "egsmax"
       emg,          &! ground emissivity (0.97 for snow,
! glaciers and water surface; 0.96 for soil and wetland)
       errore,       &! energy balnce error [w/m2]
       etrc,         &! maximum possible transpiration rate [mm/s]
       fac,          &! soil wetness of surface layer
       fact(lb:nl_soil), &! used in computing tridiagonal matrix
       fsun,         &! fraction of sunlit canopy
       hr,           &! relative humidity
       htvp,         &! latent heat of vapor of water (or sublimation) [j/kg]
       olru,         &! olrg excluding dwonwelling reflection [W/m2]
       olrb,         &! olrg assuming blackbody emission [W/m2]
       psit,         &! negative potential of soil
       par,          &! PAR absorbed by canopy [W/m2]
       qg,           &! ground specific humidity [kg/kg]
       qsatg,        &! saturated humidity [kg/kg]
       qsatgdT,      &! d(qsatg)/dT
       qred,         &! soil surface relative humidity
!       rstfac,       &! factor of soil water stress
       sabv,         &! solar absorbed by canopy [W/m2]
       thm,          &! intermediate variable (tm+0.0098*ht)
       th,           &! potential temperature (kelvin)
       thv,          &! virtual potential temperature (kelvin)
       tl,           &! leaf temperature
       tg,           &! ground surface temperature [K]
       tssbef(lb:nl_soil), &! soil/snow temperature before update
       tinc,         &! temperature difference of two time step
       ur,           &! wind speed at reference height [m/s]
       ulrad,        &! upward longwave radiation above the canopy [W/m2]
       wice0(lb:nl_soil),&! ice mass from previous time-step
       wliq0(lb:nl_soil),&! liquid mass from previous time-step
       wx,           &! patitial volume of ice and water of surface layer
       xmf            ! total latent heat of phase change of ground water
  real(r8) :: z0ma_g,zol_g,rib_g,ustar_g,qstar_g,tstar_g
  real(r8) :: f10m_g, fm_g,fh_g,fq_g,temp1,temp2,temp12m,temp22m,um,obu
!=======================================================================
! [1] Initial set and propositional variables
!=======================================================================
! fluxes
      taux   = 0.;  tauy   = 0.    
      fsena  = 0.;  fevpa  = 0.  
      lfevpa = 0.;  fsenl  = 0.    
      fevpl  = 0.;  etr    = 0.  
      fseng  = 0.;  fevpg  = 0.    
      dlrad  = 0.;  ulrad  = 0. 
      cgrnds = 0.;  cgrndl = 0.    
      cgrnd  = 0.;  tref   = 0. 
      qref   = 0.;  rst    = 2.0e4
      assim  = 0.;  respc  = 0. 
      emis   = 0.;  z0ma   = 0.
      zol    = 0.;  rib    = 0.
      ustar  = 0.;  qstar  = 0.
      tstar  = 0.;  rootr  = 0.
! temperature and water mass from previous time step
      tg = tss(lb)
      tssbef(lb:) = tss(lb:)
      wice0(lb:) = wice(lb:)
      wliq0(lb:) = wliq(lb:)
! emissivity
      emg = 0.96
      if(scv>0. .OR. itypwat==3) emg = 0.97
! latent heat, assumed that the sublimation occured only as wliq=0
      htvp = hvap
      if(wliq(lb)<=0. .AND. wice(lb)>0.) htvp = hsub
! potential temperatur at the reference height
      thm = tm + 0.0098*ht              ! intermediate variable equivalent to
! tm*(pgcm/psrf)**(rgas/cpair)
      th = tm*(100000./psrf)**(rgas/cpair) ! potential T
      thv = th*(1.+0.61*qm)             ! virtual potential T
      ur = max(0.1,sqrt(us*us+vs*vs))   ! limit set to 0.1
!=======================================================================
! [2] specific humidity and its derivative at ground surface
!=======================================================================
      qred = 1.
      call qsadv(tg,psrf,eg,degdT,qsatg,qsatgdT)
      if(itypwat<=1)then            ! soil ground
         wx   = (wliq(1)/denh2o + wice(1)/denice)/dz(1)
         if(porsl(1)<1.e-6)then     ! bed rock
            fac  = 0.001
         else 
            fac  = min(1.,wx/porsl(1))
            fac  = max( fac, 0.001 )
         endif
         psit = -phi0(1) * fac ** (- bsw(1) )   ! psit = max(smpmin, psit)
         hr   = exp(psit/roverg/tg)
         qred = (1.-fsno)*hr + fsno
      endif
      qg = qred*qsatg  
      dqgdT = qred*qsatgdT
      if(qsatg > qm .AND. qm > qred*qsatg)then
        qg = qm; dqgdT = 0.
      endif
!=======================================================================
! [3] Compute sensible and latent fluxes and their derivatives with respect
!     to ground temperature using ground temperatures from previous time step.
!=======================================================================
      if(sigf <= 0.999) then
         call groundfluxes (zlnd,zsno,hu,ht,hq, &
                            us,vs,tm,qm,rhoair,psrf, &
                            ur,thm,th,thv,tg,qg,dqgdT,htvp, &
                            fsno,sigf,cgrnd,cgrndl,cgrnds, &
                            taux,tauy,fsena,fevpa,fseng,fevpg,tref,qref, &
              z0ma_g,zol_g,rib_g,ustar_g,qstar_g,tstar_g,f10m_g,fm_g,fh_g,fq_g)
      end if
!=======================================================================
! [4] Canopy temperature, fluxes from the canopy
!=======================================================================
      par = parsun + parsha
      sabv = sabvsun + sabvsha
      if(sigf >= 0.001) then
! soil water strees factor on stomatal resistance
         call eroot (nl_soil,trsmx0,porsl,bsw,phi0,rootfr,&
                             dz,tss,wliq,rootr,etrc,rstfac)
!         print*,'THERMAL.F90-test-rstfac',rstfac
! fraction of sunlit and shaded leaves of canopy
         fsun = ( 1. - exp(-extkb*lai) ) / max( extkb*lai, 1.e-6 )
         if(coszen<=0.0 .OR. sabv<1.) fsun = 0.
         if(fsun.le.0.1)then
            fsun = 0.
            tl = tlsha
            call leaftemone (dtime ,csoilc ,dewmx  ,htvp   ,lai    ,&
                 sai     ,displa   ,sqrtdi ,z0m    ,effcon ,vmax25 ,&
                 slti    ,hlti     ,shti   ,hhti   ,trda   ,trdm   ,&
                 trop    ,gradm    ,binter ,extkn  ,extkb  ,extkd  ,&
                 hu      ,ht       ,hq     ,us     ,vs     ,thm    ,&
                 th      ,thv      ,qm     ,psrf   ,rhoair ,par    ,&
                 sabv    ,frl      ,thermk ,rstfac ,po2m   ,pco2m  ,&
                 sigf    ,etrc     ,tg     ,qg     ,dqgdT  ,emg    ,&
                 tl      ,ldew     ,taux   ,tauy   ,fseng  ,fevpg  ,&
                 cgrnd   ,cgrndl   ,cgrnds ,tref   ,qref   ,rst    ,&
                 assim   ,respc    ,fsenl  ,fevpl  ,etr    ,dlrad  ,&
                 ulrad   ,z0ma     ,zol    ,rib    ,ustar  ,qstar  ,&
                 tstar   ,f10m     ,fm     ,fh     ,fq)
                 tlsun = tl
		 tlsha = tl
         else
            call leaftemtwo (dtime ,csoilc  ,dewmx  ,htvp   ,lai    ,&
                 sai     ,displa   ,sqrtdi  ,z0m    ,effcon ,vmax25 ,&
                 slti    ,hlti     ,shti    ,hhti   ,trda   ,trdm   ,&
                 trop    ,gradm    ,binter  ,extkn  ,extkb  ,extkd  ,&
                 hu      ,ht       ,hq      ,us     ,vs     ,thm    ,&
                 th      ,thv      ,qm      ,psrf   ,rhoair ,parsun ,&
                 parsha  ,sabvsun  ,sabvsha ,frl    ,fsun   ,thermk ,&
                 rstfac  ,po2m     ,pco2m   ,sigf   ,etrc   ,tg     ,&
                 qg      ,dqgdT    ,emg     ,tlsun  ,tlsha  ,ldew   ,&
                 taux    ,tauy     ,fseng   ,fevpg  ,cgrnd  ,cgrndl ,&
                 cgrnds  ,tref     ,qref    ,rst    ,assim  ,respc  ,&
                 fsenl   ,fevpl    ,etr     ,dlrad  ,ulrad  ,z0ma   ,&
                 zol     ,rib      ,ustar   ,qstar  ,tstar  ,f10m   ,&
                 fm      ,fh       ,fq)
         endif
      endif
! equate canopy temperature to air over bareland.
! required as sigf=0 carried over to next time step
      if(sigf < 0.001)then
         tlsun = tm
         tlsha = tm
         ldew = 0.
      endif
!=======================================================================
! [5] Gound temperature
!=======================================================================
      call groundtem (itypwat,lb,nl_soil,dtime, &
                      capr,cnfac,csol,porsl,dkmg,dkdry,dksatu, &
                      sigf,dz,z,zi,tss,wice,wliq,scv,snowdp, &
                      frl,dlrad,sabg,fseng,fevpg,cgrnd,htvp,emg, &
                      imelt,sm,xmf,fact)
!=======================================================================
! [6] Correct fluxes to present soil temperature
!=======================================================================
      tg = tss(lb)
      tinc = tss(lb) - tssbef(lb)
      fseng = fseng + tinc*cgrnds 
      fevpg = fevpg + tinc*cgrndl
! calculation of evaporative potential; flux in kg m-2 s-1.
! egidif holds the excess energy if all water is evaporated
! during the timestep.  this energy is later added to the sensible heat flux.
      egsmax = (wice(lb)+wliq(lb)) / dtime
      egidif = max( 0., fevpg - egsmax )
      fevpg = min ( fevpg, egsmax )
      fseng = fseng + htvp*egidif
! total fluxes to atmosphere
      fsena = fsenl + fseng
      fevpa = fevpl + fevpg
      lfevpa= hvap*fevpl + htvp*fevpg   ! w/m2 (accouting for sublimation)
!      print*,'THERMAL.F90-test-lfevpa',lfevpa
!      print*,'THERMAL.F90-test-fvap',hvap,'THERMAL.F90-test-ftvp',htvp
      qseva = 0.
      qsubl = 0.
      qfros = 0.
      qsdew = 0.
      if(fevpg >= 0.)then
! not allow for sublimation in melting (melting ==> evap. ==> sublimation)
         qseva = min(wliq(lb)/dtime, fevpg)
         qsubl = fevpg - qseva
      else
         if(tg < tfrz)then
            qfros = abs(fevpg)
         else
            qsdew = abs(fevpg)
         endif
      endif
! ground heat flux
      fgrnd = sabg + dlrad + (1.-sigf)*emg*frl &
            - emg*stefnc*tssbef(lb)**3*(tssbef(lb) + 4.*tinc) &
            - (fseng+fevpg*htvp)
! outgoing long-wave radiation from canopy + ground
      olrg = ulrad &
           + (1.-sigf)*(1.-emg)*frl &
           + (1.-sigf)*emg*stefnc * tssbef(lb)**4 &
! for conservation we put the increase of ground longwave to outgoing
           + 4.*emg*stefnc*tssbef(lb)**3*tinc
! averaged bulk surface emissivity
      olrb = stefnc*tssbef(lb)**3*((1.-sigf)*tssbef(lb) + 4.*tinc)
      olru = ulrad + emg*olrb
      olrb = ulrad + olrb
      emis = olru / olrb
! radiative temperature
      trad = (olrg/stefnc)**0.25
! additonal variables required by WRF and RSM model
      if(sigf < 0.001)then
         ustar = ustar_g
         tstar = tstar_g
         qstar = qstar_g
         rib   = rib_g
         zol   = zol_g
	 z0ma  = z0ma_g
         f10m  = f10m_g
         fm    = fm_g
         fh    = fh_g
         fq    = fq_g
      else if(sigf <= 0.7)then
	 z0ma  = sigf*z0ma  + (1.-sigf)*z0ma_g
! assumed um ~= ur here
         um = ur
         ustar =   sqrt(max(1.e-6,sqrt(taux*taux+tauy*tauy))/rhoair)
         tstar = - fsena/(cpair*ustar*rhoair)
         qstar = - fevpa/(ustar*rhoair)
         zol = (hu-displa)*vonkar*grav*(tstar+0.61*th*qstar)/(ustar**2*thv)
         if(zol .ge. 0.)then   !stable
            zol = min(2.,max(zol,1.e-6))
         else                  !unstable
            zol = max(-100.,min(zol,-1.e-6))
         endif
         obu = (hu-displa)/zol
         call moninobuk(hu,ht,hq,displa,z0ma,z0ma,z0ma,obu,um,&
              ustar,temp1,temp2,temp12m,temp22m,f10m,fm,fh,fq)
         rib = min(5.,zol*ustar**2/(vonkar*temp1*um**2))
      else
         ustar = ustar
         tstar = tstar
         qstar = qstar
         rib   = rib
         zol   = zol
	 z0ma  = z0ma
         f10m  = f10m
         fm    = fm
         fh    = fh
         fq    = fq
      endif
      u10m  = us/ur * ustar/vonkar * f10m
      v10m  = vs/ur * ustar/vonkar * f10m
!=======================================================================
! [7] energy balance error
!=======================================================================
      errore = sabv + sabg + frl - olrg - fsena - lfevpa - xmf
      do j = lb, nl_soil
         errore = errore - (tss(j)-tssbef(j))/fact(j)
      enddo
! if(abs(errore)>.2)then
! write(6,*) 'THERMAL.F90 : energy  balance violation'
! write(6,100) errore,sabv,sabg,frl,olrg,fsenl,fseng,hvap*fevpl,htvp*fevpg,xmf
! endif
100   format(10f10.3)
 END SUBROUTINE THERMAL
