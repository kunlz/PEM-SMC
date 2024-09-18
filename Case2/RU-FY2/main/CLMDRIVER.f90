 subroutine CLMDRIVER (nl_soil,maxsnl,numpatch,idate,deltim,&
                       nftune,nfcon,nforc,nfvar,nfldv,&
                       ftune,fcon,forc,fvar,fldv,rstfac,&
                       dolai,doalb,dosst,oro) 
!=======================================================================
!
! CLM MODEL DRIVER
!
! Original author : Yongjiu Dai, 09/30/1999; 08/30/2002
!
!=======================================================================
 use precision
 use phycon_module, only : tfrz, rgas, vonkar
 implicit none
! ------------------- arguments  ----------------------------------
  integer, INTENT(in) :: nl_soil  ! number of soil layers
  integer, INTENT(in) :: maxsnl   ! max number of snow layers
  integer, INTENT(in) :: numpatch ! number of clm grid points
  integer, INTENT(in) :: idate(3) ! model calendar for next time step (year, julian day, seconds)
  real(r8), INTENT(in) :: deltim  ! time-step
  integer, INTENT(in) :: nftune   ! number of clm tunable constants
  integer, INTENT(in) :: nfcon    ! number of time constant variables
  integer, INTENT(in) :: nforc    ! number of forcing variables
  integer, INTENT(in) :: nfvar    ! number of time varying variables
  integer, INTENT(in) :: nfldv    ! number of output variables
  real(r8), INTENT(in)    :: ftune(nftune) ! clm tunable constants
  real(r8), INTENT(in)    :: fcon (numpatch,nfcon) ! time constant variables (clm restatrt required)
  real(r8), INTENT(in)    :: forc (numpatch,nforc) ! forcing variables
  real(r8), INTENT(inout) :: fvar (numpatch,nfvar) ! time varying variables (clm restatrt required)
  real(r8), INTENT(out)   :: fldv (numpatch,nfldv) ! output variables
  real(r8), INTENT(out)   :: rstfac(numpatch) !factor of water stress
  logical, INTENT(in) :: dolai    ! true if time for time-varying vegetation paramter
  logical, INTENT(in) :: doalb    ! true if time for surface albedo calculation
  logical, INTENT(in) :: dosst    ! true if time for update sst/ice/snow
  real(r8), INTENT(inout) :: oro(numpatch)         ! ocean(0)/seaice(2)/ flag
! ----------------------------------------------------------------
! I. Time invariant model variables
! ----------------------------------------------------------------
  real(r8) dlat     (numpatch) ! latitude in radians
  real(r8) dlon     (numpatch) ! longitude in radians
  integer  itypwat  (numpatch) ! land water type
  integer  ivt      (numpatch) ! land cover type of classification of USGS etc.
! Soil physical parameters
  real(r8) albsol        (numpatch) ! soil albedo for different coloured soils [-]
  real(r8) csol  (nl_soil,numpatch) ! heat capacity of soil solids [J/(m3 K)]
  real(r8) porsl (nl_soil,numpatch) ! fraction of soil that is voids [-]
  real(r8) phi0  (nl_soil,numpatch) ! minimum soil suction [mm]
  real(r8) bsw   (nl_soil,numpatch) ! clapp and hornbereger "b" parameter [-]
  real(r8) dkmg  (nl_soil,numpatch) ! thermal conductivity of soil minerals [W/m-K]
  real(r8) dksatu(nl_soil,numpatch) ! thermal conductivity of saturated soil [W/m-K]
  real(r8) dkdry (nl_soil,numpatch) ! thermal conductivity for dry soil  [W/(m-K)]
  real(r8) hksati(nl_soil,numpatch) ! hydraulic conductivity at saturation [mm h2o/s]
! Vegetation static parameters
  real(r8) z0m      (numpatch) ! aerodynamic roughness length [m]
  real(r8) displa   (numpatch) ! displacement height [m]
  real(r8) sqrtdi   (numpatch) ! inverse sqrt of leaf dimension [m**-0.5]
  real(r8) effcon   (numpatch) ! quantum efficiency of RuBP regeneration (molCO2/molquanta)
  real(r8) vmax25   (numpatch) ! maximum carboxylation rate at 25 C at canopy top
  real(r8) slti     (numpatch) ! s3: slope of low temperature inhibition function
  real(r8) hlti     (numpatch) ! s4: 1/2 point of low temperature inhibition function
  real(r8) shti     (numpatch) ! s1: slope of high temperature inhibition function
  real(r8) hhti     (numpatch) ! s2: 1/2 point of high temperature inhibition function
  real(r8) trda     (numpatch) ! s5: temperature coefficient in gs-a model
  real(r8) trdm     (numpatch) ! s6: temperature coefficient in gs-a model
  real(r8) trop     (numpatch) ! temperature coefficient in gs-a model
  real(r8) gradm    (numpatch) ! conductance-photosynthesis slope parameter
  real(r8) binter   (numpatch) ! conductance-photosynthesis intercep
  real(r8) extkn    (numpatch) ! coefficient of leaf nitrogen allocation
  real(r8) chil     (numpatch) ! leaf angle distribution factor
  real(r8) ref  (2,2,numpatch) ! leaf reflectance (iw=iband, il=life and dead)
  real(r8) tran (2,2,numpatch) ! leaf transmittance (iw=iband, il=life and dead)
  real(r8) rootfr(nl_soil,numpatch)  ! fraction of roots in each soil layer
! CLM time step and TUNABLE constants
  real(r8) dtime               ! CLM time step [seconds]
  real(r8) zlnd                ! roughness length for soil [m]
  real(r8) zsno                ! roughness length for snow [m]
  real(r8) csoilc              ! drag coefficient for soil under canopy [-]
  real(r8) dewmx               ! maximum dew
  real(r8) wtfact              ! fraction of model area with high water table
  real(r8) capr                ! tuning factor to turn first layer T into surface T
  real(r8) cnfac               ! Crank Nicholson factor between 0 and 1
  real(r8) ssi                 ! irreducible water saturation of snow
  real(r8) wimp                ! water impremeable if porosity less than wimp
  real(r8) pondmx              ! ponding depth (mm)
  real(r8) smpmax              ! wilting point potential in mm
  real(r8) smpmin              ! restriction for min of soil poten. (mm)
  real(r8) trsmx0              ! max transpiration for moist soil+100% veg.  [mm/s]
  real(r8) tcrit               ! critical temp. to determine rain or snow
! -----------------------------------------------------------------
! II. Time-varying state variables which reaquired by restart run
! -----------------------------------------------------------------
! Currrnt calendar
  integer year                 ! current year of model run
  integer jday                 ! current julian day of model run
  integer msec                 ! current seconds of model run (0 - 86400)
! Main land surface variables
  real(r8) z   (maxsnl+1:nl_soil,numpatch) ! node depth [m]
  real(r8) dz  (maxsnl+1:nl_soil,numpatch) ! interface depth [m]
  real(r8) tss (maxsnl+1:nl_soil,numpatch) ! soil temperature [K]
  real(r8) wliq(maxsnl+1:nl_soil,numpatch) ! liquid water in layers [kg/m2]
  real(r8) wice(maxsnl+1:nl_soil,numpatch) ! ice lens in layers [kg/m2]
  real(r8) tg       (numpatch) ! ground surface temperature [K]
  real(r8) tlsun    (numpatch) ! sunlit leaf temperature [K]
  real(r8) tlsha    (numpatch) ! shaded leaf temperature [K]
  real(r8) ldew     (numpatch) ! depth of water on foliage [mm]
  real(r8) sag      (numpatch) ! non dimensional snow age [-]
  real(r8) scv      (numpatch) ! snow cover, water equivalent [mm]
  real(r8) snowdp   (numpatch) ! snow depth [meter]
! Vegetation dynamic parameters
  real(r8) fveg     (numpatch) ! fraction of vegetation cover
  real(r8) fsno     (numpatch) ! fraction of snow cover on ground
  real(r8) sigf     (numpatch) ! fraction of veg cover, excluding snow-covered veg [-]
  real(r8) green    (numpatch) ! leaf greenness
  real(r8) lai      (numpatch) ! leaf area index
  real(r8) sai      (numpatch) ! stem area index
! Radiation  related (albedoes)
  real(r8) coszen   (numpatch) ! cosine of solar zenith angle
  real(r8) albg (2,2,numpatch) ! albedo, ground [-]
  real(r8) albv (2,2,numpatch) ! albedo, vegetation [-]
  real(r8) alb  (2,2,numpatch) ! averaged albedo [-]
  real(r8) ssun (2,2,numpatch) ! sunlit canopy absorption for solar radiation (0-1)
  real(r8) ssha (2,2,numpatch) ! shaded canopy absorption for solar radiation (0-1)
  real(r8) thermk   (numpatch) ! canopy gap fraction for tir radiation
  real(r8) extkb    (numpatch) ! (k, g(mu)/mu) direct solar extinction coefficient
  real(r8) extkd    (numpatch) ! diffuse and scattered diffuse PAR extinction coefficient
! Additional variables required by reginal model (WRF & RSM)
  real(r8) trad     (numpatch) ! radiative temperature of surface [K]
  real(r8) tref     (numpatch) ! 2 m height air temperature [kelvin]
  real(r8) qref     (numpatch) ! 2 m height air specific humidity
  real(r8) rst      (numpatch) ! canopy stomatal resistance (s/m)
  real(r8) emis     (numpatch) ! averaged bulk surface emissivity
  real(r8) z0ma     (numpatch) ! effective roughness [m]
  real(r8) zol      (numpatch) ! dimensionless height (z/L) used in Monin-Obukhov theory
  real(r8) rib      (numpatch) ! bulk Richardson number in surface layer
  real(r8) ustar    (numpatch) ! u* in similarity theory [m/s]
  real(r8) qstar    (numpatch) ! q* in similarity theory [kg/kg]
  real(r8) tstar    (numpatch) ! t* in similarity theory [K]
  real(r8) u10m     (numpatch) ! 10m u-velocity
  real(r8) v10m     (numpatch) ! 10m v-velocity
  real(r8) f10m     (numpatch) ! integral of profile function for momentum at 10m
  real(r8) fm       (numpatch) ! integral of profile function for momentum
  real(r8) fh       (numpatch) ! integral of profile function for heat
  real(r8) fq       (numpatch) ! integral of profile function for moisture
! -----------------------------------------------------------------
! III. Forcing
! -----------------------------------------------------------------
! Forcing
  real(r8) pco2m    (numpatch) ! CO2 concentration in atmos. (35 pa)
  real(r8) po2m     (numpatch) ! O2 concentration in atmos. (20900 pa)
  real(r8) us       (numpatch) ! wind in eastward direction [m/s]
  real(r8) vs       (numpatch) ! wind in northward direction [m/s]
  real(r8) tm       (numpatch) ! temperature at reference height [kelvin]
  real(r8) qm       (numpatch) ! specific humidity at reference height [kg/kg]
  real(r8) prc      (numpatch) ! convective precipitation [mm/s]
  real(r8) prl      (numpatch) ! large scale precipitation [mm/s]
  real(r8) psrf     (numpatch) ! atmospheric pressure at the surface [pa]
  real(r8) pbot     (numpatch) ! atm bottom level pressure (or reference height) (pa)
  real(r8) sols     (numpatch) ! atm vis direct beam solar rad onto srf [W/m2]
  real(r8) soll     (numpatch) ! atm nir direct beam solar rad onto srf [W/m2]
  real(r8) solsd    (numpatch) ! atm vis diffuse solar rad onto srf [W/m2]
  real(r8) solld    (numpatch) ! atm nir diffuse solar rad onto srf [W/m2]
  real(r8) frl      (numpatch) ! atmospheric infrared (longwave) radiation [W/m2]
  real(r8) hu       (numpatch) ! observational height of wind [m]
  real(r8) ht       (numpatch) ! observational height of temperature [m]
  real(r8) hq       (numpatch) ! observational height of humidity [m]
!  real(r8) rstfac   (numpatch) ! factor of water stress
  real(r8) rhoair   (numpatch) ! air density [kg/m3]
! -----------------------------------------------------------------
! IV. Fluxes
! -----------------------------------------------------------------
  real(r8) taux     (numpatch) ! wind stress: E-W [kg/m/s2]
  real(r8) tauy     (numpatch) ! wind stress: N-S [kg/m/s2]
  real(r8) fsena    (numpatch) ! sensible heat from canopy height to atmosphere [W/m2]
  real(r8) lfevpa   (numpatch) ! latent heat flux from canopy height to atmosphere [W/m2]
  real(r8) fevpa    (numpatch) ! evapotranspiration from canopy to atmosphere [mm/s]
  real(r8) fsenl    (numpatch) ! sensible heat from leaves [W/m2]
  real(r8) fevpl    (numpatch) ! evaporation+transpiration from leaves [mm/s]
  real(r8) etr      (numpatch) ! transpiration rate [mm/s]
  real(r8) fseng    (numpatch) ! sensible heat flux from ground [W/m2]
  real(r8) fevpg    (numpatch) ! evaporation heat flux from ground [mm/s]
  real(r8) fgrnd    (numpatch) ! ground heat flux [W/m2]
  real(r8) sabvsun  (numpatch) ! solar absorbed by sunlit vegetation [W/m2]
  real(r8) sabvsha  (numpatch) ! solar absorbed by shaded vegetation [W/m2]
  real(r8) sabg     (numpatch) ! solar absorbed by ground  [W/m2]
  real(r8) olrg     (numpatch) ! outgoing long-wave radiation from ground+canopy [W/m2]
  real(r8) rnet     (numpatch) ! net radiation by surface [W/m2]
  real(r8) xerr     (numpatch) ! the error of water banace [mm/s]
  real(r8) zerr     (numpatch) ! the error of energy balance [W/m2]
  real(r8) rsur     (numpatch) ! surface runoff (mm h2o/s)
  real(r8) rnof     (numpatch) ! total runoff (mm h2o/s)
  real(r8) assim    (numpatch) ! canopy assimilation rate (mol m-2 s-1)
  real(r8) respc    (numpatch) ! canopy respiration (mol m-2 s-1)
! -----------------------------------------------------------------
! V. Local declaration
! -----------------------------------------------------------------
  integer i,j,lb,ub,jm         ! loop/array indices
  real(r8) parsun   (numpatch) ! PAR by sunlit leaves [W/m2]
  real(r8) parsha   (numpatch) ! PAR by shaded leaves [W/m2]
  real(r8) sabvg    (numpatch) ! solar absorbed by ground + vegetation [W/m2]
  real(r8) work(numpatch)
! ======================================================================
!  [1] Transfer the time invariant and time-varying variables
! ======================================================================
! Time invariant model variables
      if(nfcon /= 9*nl_soil+29) call abort
!     print*,'-------- clm fcon -------'
!     do ub=1,nfcon
!     work(:)= fcon(1:numpatch,ub)
!     print 100, ub, maxval(work), minval(work)
!     enddo
100   format(1h (i6, 2e16.3), '        nfcon')
      do i = 1, numpatch
         ub = 1
         dlat   (i) = fcon(i,ub)              ; ub = ub + 1                     !1
         dlon   (i) = fcon(i,ub)              ; ub = ub + 1                     !2
         itypwat(i) = nint(fcon(i,ub))        ; ub = ub + 1                     !3
         ivt    (i) = nint(fcon(i,ub))        ; ub = ub + 1                     !4
         albsol (i) = fcon(i,ub)              ; lb = ub + 1; ub = ub + nl_soil  !5
         csol   (1:nl_soil,i) = fcon(i,lb:ub) ; lb = ub + 1; ub = ub + nl_soil  !1_
         porsl  (1:nl_soil,i) = fcon(i,lb:ub) ; lb = ub + 1; ub = ub + nl_soil  !2_
         phi0   (1:nl_soil,i) = fcon(i,lb:ub) ; lb = ub + 1; ub = ub + nl_soil  !3_
         bsw    (1:nl_soil,i) = fcon(i,lb:ub) ; lb = ub + 1; ub = ub + nl_soil  !4_
         dkmg   (1:nl_soil,i) = fcon(i,lb:ub) ; lb = ub + 1; ub = ub + nl_soil  !5_
         dksatu (1:nl_soil,i) = fcon(i,lb:ub) ; lb = ub + 1; ub = ub + nl_soil  !6_
         dkdry  (1:nl_soil,i) = fcon(i,lb:ub) ; lb = ub + 1; ub = ub + nl_soil  !7_
         hksati (1:nl_soil,i) = fcon(i,lb:ub) ; ub = ub + 1                     !8_
         z0m     (i) = fcon(i,ub)             ; ub = ub + 1                     !6
         displa  (i) = fcon(i,ub)             ; ub = ub + 1                     !7
         sqrtdi  (i) = fcon(i,ub)             ; ub = ub + 1                     !8
         effcon  (i) = fcon(i,ub)             ; ub = ub + 1                     !9
         vmax25  (i) = fcon(i,ub)             ; ub = ub + 1                     !10
         slti    (i) = fcon(i,ub)             ; ub = ub + 1                     !11
         hlti    (i) = fcon(i,ub)             ; ub = ub + 1                     !12
         shti    (i) = fcon(i,ub)             ; ub = ub + 1                     !13
         hhti    (i) = fcon(i,ub)             ; ub = ub + 1                     !14
         trda    (i) = fcon(i,ub)             ; ub = ub + 1                     !15
         trdm    (i) = fcon(i,ub)             ; ub = ub + 1                     !16
         trop    (i) = fcon(i,ub)             ; ub = ub + 1                     !17
         gradm   (i) = fcon(i,ub)             ; ub = ub + 1                     !18
         binter  (i) = fcon(i,ub)             ; ub = ub + 1                     !19
         extkn   (i) = fcon(i,ub)             ; ub = ub + 1                     !20
         chil    (i) = fcon(i,ub)             ; ub = ub + 1                     !21
         ref (1,1,i) = fcon(i,ub)             ; ub = ub + 1                     !22
         ref (1,2,i) = fcon(i,ub)             ; ub = ub + 1                     !23
         ref (2,1,i) = fcon(i,ub)             ; ub = ub + 1                     !24
         ref (2,2,i) = fcon(i,ub)             ; ub = ub + 1                     !25
         tran(1,1,i) = fcon(i,ub)             ; ub = ub + 1                     !26
         tran(1,2,i) = fcon(i,ub)             ; ub = ub + 1                     !27
         tran(2,1,i) = fcon(i,ub)             ; ub = ub + 1                     !28
         tran(2,2,i) = fcon(i,ub)             ; lb = ub + 1; ub = ub + nl_soil  !29
         rootfr(1:nl_soil,i) = fcon(i,lb:ub)                                    !9_
      enddo
!      print*,'CLMDRIVER:test:itypwat',itypwat
! CLM time step and TUNABLE constants
      if(nftune /= 14) call abort
      dtime  = deltim
      zlnd   = ftune(1)
      zsno   = ftune(2)
      csoilc = ftune(3)
      dewmx  = ftune(4)
      wtfact = ftune(5)
      capr   = ftune(6)
      cnfac  = ftune(7)
      ssi    = ftune(8)
      wimp   = ftune(9)
      pondmx = ftune(10) 
      smpmax = ftune(11) 
      smpmin = ftune(12) 
      trsmx0 = ftune(13) 
      tcrit  = ftune(14)  
!     print*,'-------- clm ftune -------'
!     print 101, maxval(ftune), minval(ftune)
101   format(1h (2e16.3), '        nftune ')
! Time-varying variables
      year = idate(1)      
      jday = idate(2)      
      msec = idate(3)      
      jm = nl_soil+abs(maxsnl) 
      if(nfvar /= 5*jm+51) call abort
      do i = 1, numpatch
         lb = 1
         ub = jm
         z   (maxsnl+1:nl_soil,i) = fvar(i,lb:ub) ; lb = ub+1; ub = ub+jm !1_
         dz  (maxsnl+1:nl_soil,i) = fvar(i,lb:ub) ; lb = ub+1; ub = ub+jm !2_
         tss (maxsnl+1:nl_soil,i) = fvar(i,lb:ub) ; lb = ub+1; ub = ub+jm !3_
         wliq(maxsnl+1:nl_soil,i) = fvar(i,lb:ub) ; lb = ub+1; ub = ub+jm !4_
         wice(maxsnl+1:nl_soil,i) = fvar(i,lb:ub) ; ub = ub+1             !5_
         tg      (i) = fvar(i,ub)                 ; ub = ub+1             !1
         tlsun   (i) = fvar(i,ub)                 ; ub = ub+1             !2
         tlsha   (i) = fvar(i,ub)                 ; ub = ub+1             !3
         ldew    (i) = fvar(i,ub)                 ; ub = ub+1             !4
         sag     (i) = fvar(i,ub)                 ; ub = ub+1             !5
         scv     (i) = fvar(i,ub)                 ; ub = ub+1             !6
         snowdp  (i) = fvar(i,ub)                 ; ub = ub+1             !7
         fveg    (i) = fvar(i,ub)                 ; ub = ub+1             !8
         fsno    (i) = fvar(i,ub)                 ; ub = ub+1             !9
         sigf    (i) = fvar(i,ub)                 ; ub = ub+1             !10
         green   (i) = fvar(i,ub)                 ; ub = ub+1             !11
         lai     (i) = fvar(i,ub)                 ; ub = ub+1             !12
         sai     (i) = fvar(i,ub)                 ; ub = ub+1             !13
         coszen  (i) = fvar(i,ub)                 ; ub = ub+1             !14
         albg(1,1,i) = fvar(i,ub)                 ; ub = ub+1             !15
         albg(1,2,i) = fvar(i,ub)                 ; ub = ub+1             !16
         albg(2,1,i) = fvar(i,ub)                 ; ub = ub+1             !17
         albg(2,2,i) = fvar(i,ub)                 ; ub = ub+1             !18
         albv(1,1,i) = fvar(i,ub)                 ; ub = ub+1             !19
         albv(1,2,i) = fvar(i,ub)                 ; ub = ub+1             !20
         albv(2,1,i) = fvar(i,ub)                 ; ub = ub+1             !21
         albv(2,2,i) = fvar(i,ub)                 ; ub = ub+1             !22
         alb (1,1,i) = fvar(i,ub)                 ; ub = ub+1             !23
         alb (1,2,i) = fvar(i,ub)                 ; ub = ub+1             !24
         alb (2,1,i) = fvar(i,ub)                 ; ub = ub+1             !25
         alb (2,2,i) = fvar(i,ub)                 ; ub = ub+1             !26
         ssun(1,1,i) = fvar(i,ub)                 ; ub = ub+1             !27
         ssun(1,2,i) = fvar(i,ub)                 ; ub = ub+1             !28
         ssun(2,1,i) = fvar(i,ub)                 ; ub = ub+1             !29
         ssun(2,2,i) = fvar(i,ub)                 ; ub = ub+1             !30
         ssha(1,1,i) = fvar(i,ub)                 ; ub = ub+1             !31
         ssha(1,2,i) = fvar(i,ub)                 ; ub = ub+1             !32
         ssha(2,1,i) = fvar(i,ub)                 ; ub = ub+1             !33
         ssha(2,2,i) = fvar(i,ub)                 ; ub = ub+1             !34
         thermk  (i) = fvar(i,ub)                 ; ub = ub+1             !35
         extkb   (i) = fvar(i,ub)                 ; ub = ub+1             !36
         extkd   (i) = fvar(i,ub)                                         !37
      enddo
!     print*,'-------- clm fvar in -------'
!     do ub=1,nfvar
!     work(:)= fvar(1:numpatch,ub)
!     print 200, ub, maxval(work), minval(work)
!     enddo
200   format(1h (i6, 2e16.3), '        nfvar before')
! ======================================================================
!  [2] atmospheric fields to force clm
! ======================================================================
      if(nforc /= 18) call abort
      do i = 1, numpatch         !clm vector index
         pco2m(i)  = forc(i,1)   !CO2 concentration in atmos. (35 pa)
         po2m(i)   = forc(i,2)   !O2 concentration in atmos. (20900 pa)
         us(i)     = forc(i,3)   !wind in eastward direction [m/s]
         vs(i)     = forc(i,4)   !wind in northward direction [m/s]
         tm(i)     = forc(i,5)   !temperature at reference height [kelvin]
         qm(i)     = forc(i,6)   !specific humidity at reference height [kg/kg]
         prc(i)    = forc(i,7)   !convective precipitation [mm/s]
         prl(i)    = forc(i,8)   !large scale precipitation [mm/s]
         pbot(i)   = forc(i,9)   !atm bottom level pressure (or reference height) (pa)
         psrf(i)   = forc(i,10)  !atmospheric pressure at the surface [pa]
         sols(i)   = forc(i,11)  !atm vis direct beam solar rad onto srf [W/m2]
         soll(i)   = forc(i,12)  !atm nir direct beam solar rad onto srf [W/m2]
         solsd(i)  = forc(i,13)  !atm vis diffuse solar rad onto srf [W/m2]
         solld(i)  = forc(i,14)  !atm nir diffuse solar rad onto srf [W/m2]
         frl(i)    = forc(i,15)  !atmospheric infrared (longwave) radiation [W/m2]
         hu(i)     = forc(i,16)  !observational height of wind [m]
         ht(i)     = forc(i,17)  !observational height of temperature [m]
         hq(i)     = forc(i,18)  !observational height of humidity [m]
         rhoair(i) = (pbot(i)-0.378*qm(i)*pbot(i)/(0.622+0.378*qm(i)))/(rgas*tm(i))
      end do
!     print*,'-------- clm forc -------'
!     do ub=1,nforc
!     work(:)= forc(1:numpatch,ub)
!     print 300, ub, maxval(work), minval(work)
!     enddo
300   format(1h (i6, 2e16.3), '        nforc')
! ======================================================================
! [2] Main driver for CLM
! ======================================================================
      do i = 1, numpatch
      CALL CLMMAIN (dtime, doalb, dolai, dosst, &
      nl_soil, maxsnl, dlon(i), dlat(i), itypwat(i), ivt(i), oro(i), &
! soil information
      albsol(i), csol(1:,i), porsl(1:,i), phi0(1:,i), &
      bsw(1:,i), dkmg(1:,i), dksatu(1:,i), dkdry(1:,i), hksati(1:,i), &
! vegetation information
      z0m(i), displa(i), sqrtdi(i), &
      effcon(i), vmax25(i), slti(i), hlti(i), shti(i), hhti(i), &
      trda(i), trdm(i), trop(i), gradm(i), binter(i), extkn(i), &
      chil(i), ref(1:,1:,i),tran(1:,1:,i), rootfr(1:,i), &
! atmospheric forcing
      frl(i), sols(i), soll(i), solsd(i), solld(i), &
      pco2m(i), po2m(i), us(i), vs(i), tm(i), qm(i), &
      prc(i), prl(i), psrf(i), rhoair(i), &
      hu(i), ht(i), hq(i), &
! model variables needed by restart run
      year, jday, msec, &
      z(maxsnl+1:,i), dz(maxsnl+1:,i), tss(maxsnl+1:,i), &
      wliq(maxsnl+1:,i), wice(maxsnl+1:,i), &
      tg(i), tlsun(i), tlsha(i), ldew(i), &
      sag(i), scv(i), snowdp(i), &
      fveg(i), fsno(i), sigf(i), green(i), lai(i), sai(i), &
      coszen(i), albg(1:,1:,i), albv(1:,1:,i), alb(1:,1:,i), &
      ssun(1:,1:,i), ssha(1:,1:,i), thermk(i), extkb(i), extkd(i), &
! fluxes
      taux(i), tauy(i), &
      fsena(i), fevpa(i), lfevpa(i), fsenl(i), fevpl(i), etr(i), &
      fseng(i), fevpg(i), olrg(i), fgrnd(i), trad(i), tref(i), qref(i), &
      rsur(i), rnof(i), rst(i), assim(i), respc(i), &
      parsun(i),parsha(i),sabvsun(i),sabvsha(i),sabg(i),sabvg(i), &
      xerr(i), zerr(i), &
! TUNABLE modle constants
      zlnd, zsno, csoilc, dewmx, wtfact, &
      capr, cnfac, ssi, wimp, pondmx, &
      smpmax, smpmin, trsmx0, tcrit, &
! time-varying vegetation from read-in file
! additional variables required by coupling with regional model (WRF & RSM)
      emis(i), z0ma(i), zol(i), rib(i), ustar(i), qstar(i), tstar(i), &
      u10m(i), v10m(i), f10m(i), fm(i), fh(i), fq(i), rstfac(i) &
                    )
!      print*,'CLMDRIVER.F90-test-rstfac',rstfac(i)
      enddo
! ======================================================================
! [3] the model variables for restart run
! ======================================================================
      jm = nl_soil+abs(maxsnl) 
      if(nfvar /= 5*jm+51) call abort
      do i = 1, numpatch
         lb = 1
         ub = jm
         fvar(i,lb:ub) = z   (maxsnl+1:nl_soil,i); lb = ub+1; ub = ub+jm !1_
         fvar(i,lb:ub) = dz  (maxsnl+1:nl_soil,i); lb = ub+1; ub = ub+jm !2_
         fvar(i,lb:ub) = tss (maxsnl+1:nl_soil,i); lb = ub+1; ub = ub+jm !3_
         fvar(i,lb:ub) = wliq(maxsnl+1:nl_soil,i); lb = ub+1; ub = ub+jm !4_
         fvar(i,lb:ub) = wice(maxsnl+1:nl_soil,i); ub = ub+1             !5_
         fvar(i,ub)    = tg      (i)             ; ub = ub+1             !1
         fvar(i,ub)    = tlsun   (i)             ; ub = ub+1             !2
         fvar(i,ub)    = tlsha   (i)             ; ub = ub+1             !3
         fvar(i,ub)    = ldew    (i)             ; ub = ub+1             !4
         fvar(i,ub)    = sag     (i)             ; ub = ub+1             !5
         fvar(i,ub)    = scv     (i)             ; ub = ub+1             !6
         fvar(i,ub)    = snowdp  (i)             ; ub = ub+1             !7
         fvar(i,ub)    = fveg    (i)             ; ub = ub+1             !8
         fvar(i,ub)    = fsno    (i)             ; ub = ub+1             !9
         fvar(i,ub)    = sigf    (i)             ; ub = ub+1             !10
         fvar(i,ub)    = green   (i)             ; ub = ub+1             !11
         fvar(i,ub)    = lai     (i)             ; ub = ub+1             !12
         fvar(i,ub)    = sai     (i)             ; ub = ub+1             !13
         fvar(i,ub)    = coszen  (i)             ; ub = ub+1             !14
         fvar(i,ub)    = albg(1,1,i)             ; ub = ub+1             !15
         fvar(i,ub)    = albg(1,2,i)             ; ub = ub+1             !16
         fvar(i,ub)    = albg(2,1,i)             ; ub = ub+1             !17
         fvar(i,ub)    = albg(2,2,i)             ; ub = ub+1             !18
         fvar(i,ub)    = albv(1,1,i)             ; ub = ub+1             !19
         fvar(i,ub)    = albv(1,2,i)             ; ub = ub+1             !20
         fvar(i,ub)    = albv(2,1,i)             ; ub = ub+1             !21
         fvar(i,ub)    = albv(2,2,i)             ; ub = ub+1             !22
         fvar(i,ub)    = alb (1,1,i)             ; ub = ub+1             !23
         fvar(i,ub)    = alb (1,2,i)             ; ub = ub+1             !24
         fvar(i,ub)    = alb (2,1,i)             ; ub = ub+1             !25
         fvar(i,ub)    = alb (2,2,i)             ; ub = ub+1             !26
         fvar(i,ub)    = ssun(1,1,i)             ; ub = ub+1             !27
         fvar(i,ub)    = ssun(1,2,i)             ; ub = ub+1             !28
         fvar(i,ub)    = ssun(2,1,i)             ; ub = ub+1             !29
         fvar(i,ub)    = ssun(2,2,i)             ; ub = ub+1             !30
         fvar(i,ub)    = ssha(1,1,i)             ; ub = ub+1             !31
         fvar(i,ub)    = ssha(1,2,i)             ; ub = ub+1             !32
         fvar(i,ub)    = ssha(2,1,i)             ; ub = ub+1             !33
         fvar(i,ub)    = ssha(2,2,i)             ; ub = ub+1             !34
         fvar(i,ub)    = thermk  (i)             ; ub = ub+1             !35
         fvar(i,ub)    = extkb   (i)             ; ub = ub+1             !36
         fvar(i,ub)    = extkd   (i)             ; ub = ub+1             !37
! Additional variables required by reginal model (WRF & RSM)
         fvar(i,ub)    = trad    (i)             ; ub = ub+1             !38
         fvar(i,ub)    = tref    (i)             ; ub = ub+1             !39
         fvar(i,ub)    = qref    (i)             ; ub = ub+1             !40
         fvar(i,ub)    = rst     (i)             ; ub = ub+1             !41
         fvar(i,ub)    = emis    (i)             ; ub = ub+1             !42
         fvar(i,ub)    = z0ma    (i)             ; ub = ub+1             !43
         fvar(i,ub)    = zol     (i)             ; ub = ub+1             !44
         fvar(i,ub)    = rib     (i)             ; ub = ub+1             !45
         fvar(i,ub)    = ustar   (i)             ; ub = ub+1             !46
         fvar(i,ub)    = qstar   (i)             ; ub = ub+1             !47
         fvar(i,ub)    = tstar   (i)             ; ub = ub+1             !48
         fvar(i,ub)    = fm      (i)             ; ub = ub+1             !49
         fvar(i,ub)    = fh      (i)             ; ub = ub+1             !50
         fvar(i,ub)    = fq      (i)             ;                       !51
      enddo
!     print*,'-------- clm fvar out -------'
!     do ub=1,nfvar
!     work(:)= fvar(1:numpatch,ub)
!     print 400, ub, maxval(work), minval(work)
!     enddo
400   format(1h (i6, 2e16.3), '        nfvar after')
! ======================================================================
! [4] Return required surface fields to atmospheric model:
! ======================================================================
      if(nfldv /= 92) call abort
      do i=1,numpatch
         fldv(i, 1) = taux   (i)  !wind stress: E-W [kg/m/s2]
         fldv(i, 2) = tauy   (i)  !wind stress: N-S [kg/m/s2]
         fldv(i, 3) = fsena  (i)  !sensible heat from canopy height to atmosphere [W/m2]
!         print*,'clmdriver-test-H:',fsena(i)
         fldv(i, 4) = lfevpa (i)  !latent heat flux from canopy height to atmosphere [W/m2]
!         print*,'clmdriver-test-LE:',lfevpa(i)
         fldv(i, 5) = fevpa  (i)  !evapotranspiration from canopy to atmosphere [mm/s]
         fldv(i, 6) = fsenl  (i)  !sensible heat from leaves [W/m2]
         fldv(i, 7) = fevpl  (i)  !evaporation+transpiration from leaves [mm/s]
         fldv(i, 8) = etr    (i)  !transpiration rate [mm/s]
         fldv(i, 9) = fseng  (i)  !sensible heat flux from ground [W/m2]
         fldv(i,10) = fevpg  (i)  !evaporation heat flux from ground [mm/s]
         fldv(i,11) = fgrnd  (i)  !ground heat flux [W/m2]
         fldv(i,12) = sabvsun(i)  !solar absorbed by sunlit canopy [W/m2]
         fldv(i,13) = sabvsha(i)  !solar absorbed by shaded [W/m2]
         fldv(i,14) = sabg   (i)  !solar absorbed by ground  [W/m2]
         fldv(i,15) = olrg   (i)  !outgoing long-wave radiation from ground+canopy [W/m2]
         fldv(i,16) = sabvg  (i) + frl(i) - olrg(i) !net radiation [W/m2]
         fldv(i,17) = xerr   (i)  !the error of water banace [mm/s]
         fldv(i,18) = zerr   (i)  !the error of energy balance [W/m2]
         fldv(i,19) = rsur   (i)  !surface runoff [mm/s]
         fldv(i,20) = rnof   (i)  !total runoff [mm/s]
         fldv(i,21) = assim  (i)  !canopy assimilation rate [mol m-2 s-1]
         fldv(i,22) = respc  (i)  !respiration (plant+soil) [mol m-2 s-1]
                                             lb = 22+1; ub = 22+nl_soil !
         fldv(i,lb:ub) = tss(1:nl_soil,i)  ; lb = ub+1; ub = ub+nl_soil !23_32
         fldv(i,lb:ub) = wliq(1:nl_soil,i) ; lb = ub+1; ub = ub+nl_soil !33_42
         fldv(i,lb:ub) = wice(1:nl_soil,i) ;            ub = ub+1       !43_52
         fldv(i,ub) = tg     (i)  ; ub = ub+1           !53
         fldv(i,ub) = tlsun  (i)  ; ub = ub+1           !54
         fldv(i,ub) = tlsha  (i)  ; ub = ub+1           !55
         fldv(i,ub) = ldew   (i)  ; ub = ub+1           !56
         fldv(i,ub) = scv    (i)  ; ub = ub+1           !57
         fldv(i,ub) = snowdp (i)  ; ub = ub+1           !58
         fldv(i,ub) = fsno   (i)  ; ub = ub+1           !59
         fldv(i,ub) = sigf   (i)  ; ub = ub+1           !60
         fldv(i,ub) = green  (i)  ; ub = ub+1           !61
         fldv(i,ub) = lai    (i)  ; ub = ub+1           !62
         fldv(i,ub) = sai    (i)  ; ub = ub+1           !63
         fldv(i,ub) = alb(1,1,i)  ; ub = ub+1           !64
         fldv(i,ub) = alb(1,2,i)  ; ub = ub+1           !65
         fldv(i,ub) = alb(2,1,i)  ; ub = ub+1           !66
         fldv(i,ub) = alb(2,2,i)  ; ub = ub+1           !67
         fldv(i,ub) = emis   (i)  ; ub = ub+1           !68
         fldv(i,ub) = z0ma   (i)  ; ub = ub+1           !69
         fldv(i,ub) = trad   (i)  ; ub = ub+1           !70
         fldv(i,ub) = ustar  (i)  ; ub = ub+1           !71
         fldv(i,ub) = tstar  (i)  ; ub = ub+1           !72
         fldv(i,ub) = qstar  (i)  ; ub = ub+1           !73
         fldv(i,ub) = zol    (i)  ; ub = ub+1           !74
         fldv(i,ub) = rib    (i)  ; ub = ub+1           !75
         fldv(i,ub) = fm     (i)  ; ub = ub+1           !76
         fldv(i,ub) = fh     (i)  ; ub = ub+1           !77
         fldv(i,ub) = fq     (i)  ; ub = ub+1           !78
! diagnostic variables
         fldv(i,ub) = tref   (i)  ; ub = ub+1           !79
         fldv(i,ub) = qref   (i)  ; ub = ub+1           !80
         fldv(i,ub) = u10m   (i)  ; ub = ub+1           !81
         fldv(i,ub) = v10m   (i)  ; ub = ub+1           !82
         fldv(i,ub) = f10m   (i)  ; ub = ub+1           !83
! forcing
         fldv(i,ub) = us     (i)  ; ub = ub+1           !84
         fldv(i,ub) = vs     (i)  ; ub = ub+1           !85
         fldv(i,ub) = tm     (i)  ; ub = ub+1           !86
         fldv(i,ub) = qm     (i)  ; ub = ub+1           !87
         fldv(i,ub) = prc    (i)  ; ub = ub+1           !88
         fldv(i,ub) = prl    (i)  ; ub = ub+1           !89
         fldv(i,ub) = pbot   (i)  ; ub = ub+1           !90
         fldv(i,ub) = frl    (i)  ; ub = ub+1           !91
         fldv(i,ub) = sols(i)+soll(i)+solsd(i)+solld(i) !92
      enddo
!     print*,'-------- clm fldv out -------'
!     do ub=1,nfldv
!     work(:)= fldv(1:numpatch,ub)
!     print 500, ub, maxval(work), minval(work)
!     enddo
500   format(1h (i6, 2e16.3), '        nfldv -')
 end subroutine CLMDRIVER
