 subroutine CLMMAIN (dtime, doalb, dolai, dosst, &
            nl_soil, maxsnl, dlon, dlat, itypwat, ivt, oro, &
! soil information
            albsol, csol, porsl, phi0, bsw, &
            dkmg, dksatu, dkdry, hksati, &
! vegetation information
            z0m, displa, sqrtdi, &
            effcon, vmax25, slti, hlti, shti, hhti, &
            trda, trdm, trop, gradm, binter, extkn, &       
            chil, ref, tran, rootfr, &
! atmospheric forcing
            frl, sols, soll, solsd, solld, &
            pco2m, po2m, us, vs, tm, qm, &
            prc, prl, psrf, rhoair, &
            hu, ht, hq, &
! land surface variables required for restart
            year, jday, msec, &
            z, dz, tss, wliq, wice, &
            tg, tlsun, tlsha, ldew, sag, scv, snowdp, &
            fveg, fsno, sigf, green, lai, sai, &
            coszen, albg, albv, alb, ssun, ssha, thermk, extkb, extkd, &
! fluxes
            taux,  tauy, &
            fsena, fevpa, lfevpa, fsenl, fevpl, etr, &
            fseng, fevpg, olrg, fgrnd, trad, tref, qref, &
            rsur, rnof, rst, assim, respc, &
            parsun, parsha, sabvsun, sabvsha, sabg, sabvg, xerr, zerr, &
! TUNABLE modle constants
            zlnd,   zsno,   csoilc, dewmx,  wtfact, & 
            capr,   cnfac,  ssi,    wimp,   pondmx, &  
            smpmax, smpmin, trsmx0, tcrit, & 
! time-varying vegetation parameters from read-in file
! additional variables required by coupling with WRF model
            emis, z0ma, zol, rib, ustar, qstar, tstar, &
            u10m, v10m, f10m, fm, fh, fq, rstfac )
!=======================================================================
!
! Main subroutine, advance time information
!
! Original author : Yongjiu Dai, 09/30/1999; 08/30/2002
!
!    FLOW DIAGRAM FOR CLMMAIN
!
!    CLMMAIN ===> netsolar                 |> all surface
!
!                 leafinterception         |]
!                 newsnow                  |] itypwat = 0 (soil ground)
!                 THERMAL                  |]         = 1 (urban & built-up)
!                 WATER                    |]         = 2 (wetland)
!                 snowcompaction           |]         = 3 (land ice)
!                 snowlayerscombine        |]
!                 snowlayersdivide         |]
!                 snowage                  |]
!
!                 LAKE                     |> lake
!                 SOCEAN                   |> ocean and sea ice
!
!                 orb_coszen               |> all surface
!                 EcoModel (lai_empirical) |> land
!                 snowfraction             |> land
!                 albland                  |> land
!                 albocean                 |> ocean & sea ice
!
!=======================================================================
  use precision
  use phycon_module, only : tfrz 
  implicit none
! ------------------------ Dummy Argument ------------------------------
  integer, INTENT(in) :: & 
        nl_soil    , &! number of soil layers
        maxsnl     , &! maximum number of snow layers
        ivt        , &! land cover type of USGS classification or others
        itypwat       ! land water type (0=soil, 1=urban and built-up,
! 2=wetland, 3=land ice, 4=deep lake, 5=shallow lake, 99 = ocean)
  logical, INTENT(in) :: doalb   !true if time for surface albedo calculation
  logical, INTENT(in) :: dolai   !true if time for leaf area index calculation
  logical, INTENT(in) :: dosst   !true to update sst/ice/snow before calculation
! Parameters
! ----------------------
  real(r8), INTENT(in) :: &
        dtime      , &! model time step [second]
        zlnd       , &!roughness length for soil [m]
        zsno       , &!roughness length for snow [m]
        csoilc     , &!drag coefficient for soil under canopy [-]
        dewmx      , &!maximum dew
        wtfact     , &!fraction of model area with high water table
        capr       , &!tuning factor to turn first layer T into surface T
        cnfac      , &!Crank Nicholson factor between 0 and 1
        ssi        , &!irreducible water saturation of snow
        wimp       , &!water impremeable if porosity less than wimp
        pondmx     , &!ponding depth (mm)
        smpmax     , &!wilting point potential in mm
        smpmin     , &!restriction for min of soil poten.  (mm)
        trsmx0     , &!max transpiration for moist soil+100% veg.  [mm/s]
        tcrit      , &!critical temp. to determine rain or snow
        dlon       , &! logitude in radians
        dlat       , &! latitude in radians
! soil physical parameters
	albsol         , &!soil albedo for different coloured soils [-]
        csol(nl_soil)  , &! heat capacity of soil solids [J/(m3 K)]
        porsl(nl_soil) , &! fraction of soil that is voids [-]
        phi0(nl_soil)  , &! minimum soil suction [mm]
        bsw(nl_soil)   , &! clapp and hornbereger "b" parameter [-]
        dkmg(nl_soil)  , &! thermal conductivity of soil minerals [W/m-K]
        dksatu(nl_soil), &! thermal conductivity of saturated soil [W/m-K]
        dkdry(nl_soil) , &! thermal conductivity for dry soil  [J/(K s m)]
        hksati(nl_soil), &! hydraulic conductivity at saturation [mm h2o/s]
! vegetation static, dynamic, derived parameters
        z0m        , &! aerodynamic roughness length [m]
        displa     , &! displacement height [m]
        sqrtdi     , &! inverse sqrt of leaf dimension [m**-0.5]
        effcon     , &! quantum efficiency of RuBP regeneration (mol CO2/mol quanta)
        vmax25     , &! maximum carboxylation rate at 25 C at canopy top
        slti       , &! slope of low temperature inhibition function      [s3]
        hlti       , &! 1/2 point of low temperature inhibition function  [s4]
        shti       , &! slope of high temperature inhibition function     [s1]
        hhti       , &! 1/2 point of high temperature inhibition function [s2]
        trda       , &! temperature coefficient in gs-a model             [s5]
        trdm       , &! temperature coefficient in gs-a model             [s6]
        trop       , &! temperature coefficient in gs-a model
        gradm      , &! conductance-photosynthesis slope parameter
        binter     , &! conductance-photosynthesis intercep
        extkn      , &! coefficient of leaf nitrogen allocation
        chil       , &! leaf angle distribution factor
        ref(2,2)   , &! leaf reflectance (iw=iband, il=life and dead)
        tran(2,2)  , &! leaf transmittance (iw=iband, il=life and dead)
        rootfr(nl_soil) ! fraction of roots in each soil layer
! Forcing
! ----------------------
  real(r8), INTENT(in) :: &
        frl        , &! atmospheric infrared (longwave) radiation [W/m2]
        sols       , &! atm vis direct beam solar rad onto srf [W/m2]
        soll       , &! atm nir direct beam solar rad onto srf [W/m2]
        solsd      , &! atm vis diffuse solar rad onto srf [W/m2]
        solld      , &! atm nir diffuse solar rad onto srf [W/m2]
        pco2m      , &! partial pressure of CO2 at observational height [pa]
        po2m       , &! partial pressure of O2 at observational height [pa]
        us         , &! wind speed in eastward direction [m/s]
        vs         , &! wind speed in northward direction [m/s]
        tm         , &! temperature at agcm reference height [kelvin]
        qm         , &! specific humidity at agcm reference height [kg/kg]
        prc        , &! convective precipitation [mm/s]
        prl        , &! large scale precipitation [mm/s]
        psrf       , &! atmosphere pressure at the surface [pa]
        rhoair     , &! density air [kg/m3]
        hu         , &! observational height of wind [m]
        ht         , &! observational height of temperature [m]
        hq            ! observational height of humidity [m]
! Variables required for restart run
! ----------------------------------------------------------------------
  integer, INTENT(in) :: &
	year,jday,msec ! next time-step /year/julian day/second in a day/
  real(r8), INTENT(inout) :: oro  ! ocean(0)/seaice(2)/ flag
  real(r8), INTENT(inout) :: &
        z(maxsnl+1:nl_soil)   , &! layer depth (m)
        dz(maxsnl+1:nl_soil)  , &! layer thickness (m)
        tss(maxsnl+1:nl_soil) , &! soil + snow layer temperature [K]
        wliq(maxsnl+1:nl_soil), &! liquid water (kg/m2)
        wice(maxsnl+1:nl_soil), &! ice lens (kg/m2)
        tg         , &! ground surface temperature [k]
        tlsun      , &! sunlit leaf temperature [K]
        tlsha      , &! shaded leaf temperature [K]
        ldew       , &! depth of water on foliage [kg/m2/s]
        sag        , &! non dimensional snow age [-]
        scv        , &! snow mass (kg/m2)
        snowdp     , &! snow depth (m)
        fveg       , &! fraction of vegetation cover
        fsno       , &! fractional snow cover
        sigf       , &! fraction of veg cover, excluding snow-covered veg [-]
        green      , &! greenness
        lai        , &! leaf area index
        sai        , &! stem area index
        coszen     , &! cosine of solar zenith angle
        albg(2,2)  , &! albedo, ground [-]
        albv(2,2)  , &! albedo, vegetation [-]
        alb(2,2)   , &! averaged albedo [-]
        ssun(2,2)  , &! sunlit canopy absorption for solar radiation
        ssha(2,2)  , &! shaded canopy absorption for solar radiation
        thermk     , &! canopy gap fraction for tir radiation
        extkb      , &! (k, g(mu)/mu) direct solar extinction coefficient
        extkd         ! diffuse and scattered diffuse PAR extinction coefficient
! Fluxes
! ----------------------------------------------------------------------
  real(r8), INTENT(out) :: &
        taux       , &! wind stress: E-W [kg/m/s**2]
        tauy       , &! wind stress: N-S [kg/m/s**2]
        fsena      , &! sensible heat from canopy height to atmosphere [W/m2]
        fevpa      , &! evapotranspiration from canopy height to atmosphere [mm/s]
        lfevpa     , &! latent heat flux from canopy height to atmosphere [W/2]
        fsenl      , &! ensible heat from leaves [W/m2]
        fevpl      , &! evaporation+transpiration from leaves [mm/s]
        etr        , &! transpiration rate [mm/s]
        fseng      , &! sensible heat flux from ground [W/m2]
        fevpg      , &! evaporation heat flux from ground [mm/s]
        olrg       , &! outgoing long-wave radiation from ground+canopy
        fgrnd      , &! ground heat flux [W/m2]
        xerr       , &! water balance error at current time-step [mm/s]
        zerr       , &! energy balnce errore at current time-step [W/m2]
        tref       , &! 2 m height air temperature [K]
        qref       , &! 2 m height air specific humidity
        trad       , &! radiative temperature [K]
        rsur       , &! surface runoff (mm h2o/s)
        rnof       , &! total runoff (mm h2o/s)
        rst        , &! canopy stomatal resistance
        assim      , &! canopy assimilation
        respc      , &! canopy respiration
        parsun     , &! PAR by sunlit leaves [W/m2]
        parsha     , &! PAR by shaded leaves [W/m2]
        sabvsun    , &! solar absorbed by sunlit vegetation [W/m2]
        sabvsha    , &! solar absorbed by shaded vegetation [W/m2]
        sabg       , &! solar absorbed by ground  [W/m2]
        sabvg      , &! solar absorbed by ground + vegetation [W/m2]
        emis       , &! averaged bulk surface emissivity
        z0ma       , &! effective roughness [m]
        zol        , &! dimensionless height (z/L) used in Monin-Obukhov theory
        rib        , &! bulk Richardson number in surface layer
        ustar      , &! u* in similarity theory [m/s]
        qstar      , &! q* in similarity theory [kg/kg]
        tstar      , &! t* in similarity theory [K]
        u10m       , &! 10m u-velocity
        v10m       , &! 10m v-velocity
        f10m       , &! integral of profile function for momentum at 10m
        fm         , &! integral of profile function for momentum
        fh         , &! integral of profile function for heat
        fq         , &! integral of profile function for moisture
        rstfac        !factor of soil water stress
! ----------------------- Local  Variables -----------------------------
   real(r8) :: &
        calday     , &! Julian cal day (1.xx to 365.xx)
	endwb      , &! water mass at the end of time step
	errore     , &! energy balnce errore (Wm-2)
	errorw     , &! water balnce errore (mm)
        fiold(maxsnl+1:nl_soil), &! fraction of ice relative to the total water
        orb_coszen , &! cosine of the solar zenith angle
        pg         , &! water onto ground including canopy runoff [kg/(m2 s)]
        pg_rain    , &! liquid water onto ground [kg/(m2 s)]
        pg_snow    , &! ice onto ground [kg/(m2 s)]
        qseva      , &! ground surface evaporation rate (mm h2o/s)
        qsdew      , &! ground surface dew formation (mm h2o /s) [+]
        qsubl      , &! sublimation rate from snow pack (mm h2o /s) [+]
        qfros      , &! surface dew added to snow pack (mm h2o /s) [+]
        rootr(1:nl_soil)      , &! root resistance of a layer, all layers add to 1.0
	scvold     , &! snow cover for previous time step [mm]
        sm         , &! rate of snowmelt [kg/(m2 s)]
        ssw        , &! water volumetric content of soil surface layer [m3/m3]
        tssub(7)   , &! surface/sub-surface temperatures [K]
        tssea      , &! sea surface temperature [K]
	totwb      , &! water mass at the begining of time step
        wt         , &! fraction of vegetation buried (covered) by snow [-]
        zi(maxsnl:nl_soil) ! interface level below a "z" level (m)
  integer snl      , &! number of snow layers
        imelt(maxsnl+1:nl_soil), &! flag for: melting=1, freezing=2, Nothing happended=0
        lb         , &! lower bound of arrays
        j             ! do looping index
!======================================================================
!  [1] Solar absorbed by vegetation and ground
!======================================================================
      call netsolar (itypwat,sigf,albg,albv,alb,ssun,ssha,&
                     sols,soll,solsd,solld,&
                     parsun,parsha,sabvsun,sabvsha,sabg,sabvg)
!======================================================================
if(itypwat<=3)then     ! <=== not lake and ocean (itypwat = 0, 1, 2, 3)
!======================================================================
!initial set
      scvold = scv        !snow mass at previous time step
      snl = 0
      do j=maxsnl+1,0
         if(wliq(j)+wice(j)>0.) snl=snl-1
      enddo
      zi(0)=0.
      if(snl<0)then
      do j = -1, snl, -1
         zi(j)=zi(j+1)-dz(j+1)
      enddo
      endif
      do j = 1,nl_soil
         zi(j)=zi(j-1)+dz(j)
      enddo
      totwb = ldew + scv + sum(wice(1:)+wliq(1:))
      if(snl<0) fiold(snl+1:0)=wice(snl+1:0)/(wliq(snl+1:0)+wice(snl+1:0))
!----------------------------------------------------------------------
! [2] Canopy interception and precipitation onto ground surface
!----------------------------------------------------------------------
      call leafinterception (dtime,dewmx,chil, &
                             prc,prl,tm,scv,sigf,lai,sai,ldew,pg)
!----------------------------------------------------------------------
! [3] Initilize new snow nodes for snowfall / sleet
!----------------------------------------------------------------------
      call newsnow (itypwat,maxsnl,dtime,tm,tg,pg,tcrit, &
                    zi(:0),z(:0),dz(:0),tss(:0),wliq(:0),wice(:0),fiold(:0),&
                    snl,sag,scv,snowdp,pg_rain,pg_snow)
!----------------------------------------------------------------------
! [4] Energy AND Water balance
!----------------------------------------------------------------------
      lb  = snl + 1           ! lower bound of array
      CALL THERMAL  &
	   ( itypwat    ,lb       ,nl_soil  ,dtime   ,trsmx0   ,&
             zlnd       ,zsno     ,csoilc   ,dewmx   ,capr     ,&
             cnfac      ,csol     ,porsl    ,phi0    ,bsw      ,&
             dkmg       ,dkdry    ,dksatu   ,lai     ,sai      ,&
             z0m        ,displa   ,sqrtdi   ,rootfr  ,effcon   ,&
             vmax25     ,slti     ,hlti     ,shti    ,hhti     ,&
             trda       ,trdm     ,trop     ,gradm   ,binter   ,&
             extkn      ,hu       ,ht       ,hq      ,us       ,&
             vs         ,tm       ,qm       ,rhoair  ,psrf     ,&
             pco2m      ,po2m     ,coszen   ,parsun  ,parsha   ,&
             sabvsun    ,sabvsha  ,sabg     ,frl     ,extkb    ,&
             extkd      ,thermk   ,fsno     ,sigf    ,dz(lb:)  ,&
             z(lb:)     ,zi(lb-1:),tlsun    ,tlsha   ,tss(lb:) ,&
             wice(lb:)  ,wliq(lb:),ldew     ,scv     ,snowdp   ,&
             imelt(lb:) ,taux     ,tauy     ,fsena   ,fevpa    ,&
             lfevpa     ,fsenl    ,fevpl    ,etr     ,fseng    ,&
             fevpg      ,olrg     ,fgrnd    ,rootr   ,qseva    ,&
             qsdew      ,qsubl    ,qfros    ,sm      ,tref     ,&
             qref       ,trad     ,rst      ,assim   ,respc    ,&
             errore     ,emis     ,z0ma     ,zol     ,rib      ,&
             ustar      ,qstar    ,tstar    ,u10m    ,v10m     ,&
             f10m       ,fm       ,fh       ,fq      ,rstfac  )
      CALL WATER    ( itypwat     ,lb       ,nl_soil ,dtime    ,&
             z(lb:)     ,dz(lb:)  ,zi(lb-1:),bsw     ,porsl    ,&
             phi0       ,hksati   ,rootr    ,tss(lb:),wliq(lb:),&
             wice(lb:)  ,pg_rain  ,sm       ,etr     ,qseva    ,&
             qsdew      ,qsubl    ,qfros    ,rsur    ,rnof     ,&
             wtfact     ,pondmx   ,ssi      ,wimp    ,smpmin   )
      if(snl<0)then
! Compaction rate for snow
! Natural compaction and metamorphosis. The compaction rate
! is recalculated for every new timestep
         lb  = snl + 1   ! lower bound of array
         call snowcompaction (lb,dtime,&
                         imelt(lb:0),fiold(lb:0),tss(lb:0),&
                         wliq(lb:0),wice(lb:0),dz(lb:0))
! Combine thin snow elements
         lb = maxsnl + 1
         call snowlayerscombine (lb,snl,&
                         z(lb:1),dz(lb:1),zi(lb-1:1),&
                         wliq(lb:1),wice(lb:1),tss(lb:1),scv,snowdp)
! Divide thick snow elements
         if(snl<0) &
         call snowlayersdivide (lb,snl,&
                         z(lb:0),dz(lb:0),zi(lb-1:0),&
                         wliq(lb:0),wice(lb:0),tss(lb:0))
      endif
! Set zero to the empty node
      if(snl==0) sag=0.
      wice(maxsnl+1:snl)=0.
      wliq(maxsnl+1:snl)=0.
      tss (maxsnl+1:snl)=0.
      z   (maxsnl+1:snl)=0.
      dz  (maxsnl+1:snl)=0.
      lb = snl + 1
      tg = tss(lb)
      ssw = min(1.,1.e-3*wliq(1)/dz(1))
! ----------------------------------------
! Update the snow age
! ----------------------------------------
      call snowage (dtime, tg, scv, scvold, sag)
! ----------------------------------------
! energy balance
! ----------------------------------------
      zerr=errore
! if(abs(errore)>.2) write(6,*) 'Warning: energy balance violation ',errore,ivt
! ----------------------------------------
! water balance
! ----------------------------------------
      endwb=sum(wice(1:)+wliq(1:))+ldew+scv
      errorw=(endwb-totwb)-(prc+prl-fevpa-rnof)*dtime
      if(itypwat>1) errorw=0.      !wetland, glacier
      xerr=errorw/dtime
! if(abs(errorw)>1.e-3) write(6,*) 'Warning: water balance violation', errorw,ivt
!======================================================================
else if(itypwat <= 5)then   ! <=== is lake (itypwat = 4 or 5)
!======================================================================
      scvold = scv          ! snow mass at previous time step
      zi(0)=0.
      do j = 1,nl_soil
      zi(j)=zi(j-1)+dz(j)
      enddo
      CALL LAKE (nl_soil,itypwat,dlat,dtime,&
           z(1:),dz(1:),zi(0:),hu,ht,hq,us,vs,tm,qm,prc,prl,rhoair,psrf,&
           sabg,frl,tg,tss(1:),wliq(1:),wice(1:),scv,snowdp,trad,tref,qref,&
           taux,tauy,fsena,fevpa,lfevpa,fseng,fevpg,olrg,fgrnd,tcrit,&
           emis,z0ma,zol,rib,ustar,qstar,tstar,u10m,v10m,f10m,fm,fh,fq)
! null data for lake component
           snl = 0
           z   (:0) = 0.0
           dz  (:0) = 0.0
           tss (:0) = 0.0
           wliq(:) = 0.0
           wice(:) = 0.0
           tlsun = tm
           tlsha = tm
           ldew  = 0.0
           fsenl = 0.0
           fevpl = 0.0
           etr   = 0.0
           rsur  = 0.0
           rnof  = 0.0
           rst   = -9999.
           assim = 0.0
           respc = 0.0
	   zerr=0.
	   xerr=0.
! Update the snow age
           call snowage (dtime, tg, scv, scvold, sag)
           ssw   = 0.0
!======================================================================
else                     ! <=== is ocean (itypwat >= 99)
!======================================================================
! simple ocean-sea ice model
    tssea = tg
    tssub (1:7) = tss (1:7) 
    CALL SOCEAN (dosst,dtime,oro,hu,ht,hq,&
                 us,vs,tm,qm,rhoair,psrf,sabg,frl,tssea,tssub(1:7),scv,&
                 taux,tauy,fsena,fevpa,lfevpa,fseng,fevpg,tref,qref,&
                 z0ma,zol,rib,ustar,qstar,tstar,u10m,v10m,f10m,fm,fh,fq,emis,olrg)
! null data for sea component
                 z   (:) = 0.0
                 dz  (:) = 0.0
                 tss (:) = 0.0; tss(1:7) = tssub(1:7)
                 wliq(:) = 0.0
                 wice(:) = 0.0
                 tg      = tssea
                 tlsun   = tm
                 tlsha   = tm
                 ldew    = 0.0
                 sag     = 0.0
                 snowdp  = scv/1000.*20. 
                 trad    = tssea
                 fsenl   = 0.0
                 fevpl   = 0.0  
                 etr     = 0.0  
                 fgrnd   = 0.0
                 rsur    = 0.0
                 rnof    = 0.0
                 rst     = -9999.
                 assim   = 0.0
                 respc   = 0.0
                 xerr    = 0.0
                 zerr    = 0.0
!======================================================================
endif
!======================================================================
! Preparation for the next time step
! 1) time-varying parameters for vegatation
! 2) fraction of snow cover
! 3) solar zenith angle and
! 4) albedos
!======================================================================
! cosine of solar zenith angle
    calday = float(jday) + float(msec)/86400.
    coszen = orb_coszen(calday,dlon,dlat)
    if(itypwat <= 5)then   ! land grid
! need to update lai and sai, fveg, green, they are done once in a day only
    if(dolai)then
! call EcoModel(ivt)
       call lai_empirical(ivt,nl_soil,rootfr,tss(1:),lai,sai,fveg,green)
    endif
! fraction of snow cover.
    call snowfraction (fveg,z0m,snowdp,wt,sigf,fsno)
! albedos
! we supposed call it every time-step, just because
! other vegeation related parameters are needed to create
!*if(doalb)then
    call albland (itypwat,albsol,chil,ref,tran,&
                  fveg,green,lai,sai,coszen,wt,fsno,scv,sag,ssw,tg,&
                  alb,albg,albv,ssun,ssha,thermk,extkb,extkd)
!*endif
    else                   ! ocean grid
!*if(doalb)then
    call albocean (oro,scv,coszen,alb)
!*endif
! null data for sea component
    lai = 0.0
    sai = 0.0
    green = 0.0
    fveg = 0.0
    sigf = 0.0
    fsno = 0.0
    albg(:,:) = alb(:,:)
    albv(:,:) = 0.0
    ssun(:,:) = 0.0
    ssha(:,:) = 0.0
    thermk = 0.0
    extkb = 0.0
    extkd = 0.0
    endif
!----------------------------------------------------------------------
 end subroutine CLMMAIN
