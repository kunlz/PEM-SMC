 subroutine WATER (itypwat ,lb     ,nl_soil ,dtime  ,&
                   z       ,dz     ,zi      ,bsw    ,porsl  ,&
                   phi0    ,hksati ,rootr   ,tss    ,wliq   ,&
                   wice    ,pg_rain,sm      ,etr    ,qseva  ,&
                   qsdew   ,qsubl  ,qfros   ,rsur   ,rnof   ,&
                   wtfact  ,pondmx ,ssi     ,wimp   ,smpmin  )
	           
!=======================================================================
! this is the main subroutine to execute the calculation of
! hydrological processes
!
! Original author : Yongjiu Dai, 09/15/1999; 08/30/2002
!
! FLOW DIAGRAM FOR WATER.F90
!
! WATER ===> snowwater
!            surfacerunoff
!            soilwater
!            subsurfacerunoff
!
!=======================================================================
  use precision
  use phycon_module, only : denice, denh2o, tfrz
  implicit none
!-----------------------Argument---------- ------------------------------
  integer, INTENT(in) :: &
        lb               , &! lower bound of array
        nl_soil          , &! upper bound of array
        itypwat             ! land water type (0=soil, 1=urban or built-up, 2=wetland,
! 3=land ice, 4=deep lake, 5=shallow lake)
  real(r8), INTENT(in) :: &
        dtime            , &! time step (s)
        wtfact           , &! fraction of model area with high water table
        pondmx           , &! ponding depth (mm)
        ssi              , &! irreducible water saturation of snow
        wimp             , &! water impremeable if porosity less than wimp
        smpmin           , &! restriction for min of soil poten. (mm)
        z (lb:nl_soil)   , &! layer depth (m)
        dz(lb:nl_soil)   , &! layer thickness (m)
        zi(lb-1:nl_soil) , &! interface level below a "z" level (m)
        bsw(1:nl_soil)   , &! Clapp-Hornberger "B"
        porsl(1:nl_soil) , &! saturated volumetric soil water content(porosity)
        phi0(1:nl_soil)  , &! saturated soil suction (mm)
        hksati(1:nl_soil), &! hydraulic conductivity at saturation (mm h2o/s)
        rootr(1:nl_soil) , &! root resistance of a layer, all layers add to 1.0
        tss(lb:nl_soil)  , &! soil/snow skin temperature (K)
        pg_rain          , &! rainfall after removal of interception (mm h2o/s)
        sm               , &! snow melt (mm h2o/s)
        etr              , &! actual transpiration (mm h2o/s)
        qseva            , &!ground surface evaporation rate (mm h2o/s)
        qsdew            , &!ground surface dew formation (mm h2o /s) [+]
        qsubl            , &!sublimation rate from snow pack (mm h2o /s) [+]
        qfros               !surface dew added to snow pack (mm h2o /s) [+]
  real(r8), INTENT(inout) :: &
        wice(lb:nl_soil) , &! ice lens (kg/m2)
        wliq(lb:nl_soil)    ! liquid water (kg/m2)
  real(r8), INTENT(out) :: &
        rsur             , &! surface runoff (mm h2o/s)
        rnof                ! total runoff (mm h2o/s)
!
!-----------------------Local Variables------------------------------
!
  integer i                 ! loop counter
  real(r8) :: & 
  eff_porosity(1:nl_soil), &! effective porosity = porosity - vol_ice
       hk(1:nl_soil)     , &! hydraulic conductivity (mm h2o/s)
       dhkdw(1:nl_soil)  , &! d(hk)/d(vol_liq)
       dwat(1:nl_soil)   , &! change in soil water
       gwat              , &! net water input from top
       qinfl             , &! infiltration rate (mm h2o/s)
       rsubst            , &! subsurface runoff (mm h2o/s)
       vol_liq(1:nl_soil), &! partitial volume of liquid water in layer
       vol_ice(1:nl_soil), &! partitial volume of ice lens in layer
       zmm (1:nl_soil)   , &! layer depth (mm)
       dzmm(1:nl_soil)   , &! layer thickness (mm)
       zimm(0:nl_soil)      ! interface level below a "z" level (mm)
!=======================================================================
! [1] update the liquid water within snow layer and the water onto soil
!=======================================================================
      if (lb>=1)then
         gwat = pg_rain + sm - qseva
      else
         call snowwater (lb,dtime,ssi,wimp,&
         pg_rain,qseva,qsdew,qsubl,qfros,dz(lb:0),wice(lb:0),wliq(lb:0),gwat)
      endif
!=======================================================================
! [2] surface runoff and infiltration
!=======================================================================
  if(itypwat<=1)then   ! soil ground only
! porosity of soil, partitial volume of ice and liquid
      do i = 1, nl_soil
         vol_ice(i) = min(porsl(i), wice(i)/(dz(i)*denice))
         eff_porosity(i) = porsl(i)-vol_ice(i)
         vol_liq(i) = min(eff_porosity(i), wliq(i)/(dz(i)*denh2o))
      enddo
! surface runoff including water table
      if (gwat > 0.) then
      call surfacerunoff (nl_soil,wtfact,wimp,bsw,porsl,phi0,hksati,&
           z(1:),dz(1:),zi(0:),vol_liq,vol_ice,eff_porosity,gwat,rsur)
      else
           rsur = 0.
      endif
! infiltration into surface soil layer
      qinfl = gwat - rsur 
!=======================================================================
! [3] determine the change of soil water
!=======================================================================
! convert length units from m to mm
      zmm(1:) = z(1:)*1000.
      dzmm(1:) = dz(1:)*1000.
      zimm(0:) = zi(0:)*1000.
      call soilwater (nl_soil,dtime,wimp,smpmin,porsl,phi0,bsw,hksati, &
                      zmm,dzmm,zimm,tss(1:),vol_liq,vol_ice,eff_porosity, &
                      qinfl,etr,rootr,dwat,hk,dhkdw)
! update the mass of liquid water
      do i= 1, nl_soil
         wliq(i) = max(0.,wliq(i)+dwat(i)*dzmm(i))
      enddo
!=======================================================================
! [4] subsurface runoff and the corrections
!=======================================================================
      call  subsurfacerunoff (nl_soil,dtime,pondmx,dzmm(1:),&
            wliq(1:),eff_porosity(1:),hk(1:),dhkdw(1:),dwat(1:),rsubst)
! total runoff
      rnof = rsubst + rsur                
! renew the ice and liquid mass due to condensation
      if(lb >= 1)then
         wliq(1) = wliq(1) + qsdew*dtime
         wice(1) = wice(1) + (qfros-qsubl)*dtime
      endif
!=======================================================================
! [6] arbitrarily hydrological processes in wetland and glacier
!=======================================================================
  else                          
      if(itypwat==2)then        ! WETLAND
         rsur=0.
         qinfl=gwat
         rsubst=0.
         rnof=0.  
         do i = 1, nl_soil
            if(tss(i)>=tfrz)then
               wice(i) = 0.0
               wliq(i) = 1000.*dz(i)
            else
               wice(i) = 1000.*dz(i)
               wliq(i) = 0.0
            endif
         enddo
      endif
      if(itypwat==3)then        ! LAND ICE
         rsur=gwat
         qinfl=0.
         rsubst=0.
         rnof=rsur
         wice(1:nl_soil) = 1000.*dz(1:nl_soil) 
         wliq(1:nl_soil) = 0.0
      endif
  endif
!-----------------------------------------------------------------------
 end subroutine WATER
