 subroutine iniTimeConst (nl_soil, ivt   , isc   , sand  , clay  , rockdep &
                         ,itypwat, zsoi  , dzsoi , albsol, csol  , porsl   &
                         ,phi0   , bsw   , dkmg  , dksatu, dkdry , hksati  &
                         ,z0m    , displa, sqrtdi, effcon, vmax25, slti    &
                         ,hlti   , shti  , hhti  , trda  , trdm  , trop    &
                         ,gradm  , binter, extkn , chil  , ref   , tran    &
                         ,rootfr , zlnd  , zsno  , csoilc, dewmx , wtfact  &
                         ,capr   , cnfac , ssi   , wimp  , pondmx, smpmax  &
                         ,smpmin , trsmx0, tcrit )  
!===========================================================================
! Initialize time invariant model variables
! Original author: Yongjiu Dai, 09/15/1999; 08/30/2002
!===========================================================================
  use precision
  implicit none
!--------------------------- Input
  integer, INTENT(in) ::   &!
        nl_soil          , &!number of model soil layers
        ivt              , &!index for vegetation type [-]
        isc                 !index for soil color type [-]
  real(r8), INTENT(in) ::  &!
	rockdep             !depth to bed rock
  real(r8), INTENT(inout) ::  &!
        sand(1:nl_soil)  , &!percent of snad
        clay(1:nl_soil)     !percent of clay
!--------------------------- Output
  real(r8), INTENT(out) :: &!--- Soil layer thickness, depths
        zsoi(1:nl_soil)  , &!soil layer depth [m]
       dzsoi(1:nl_soil)     !soil node thickness [m]
  integer, INTENT(out) ::  &!
	itypwat             !land water type (0=soil, 1=urban, 2=wetland,
!3=land ice, 4=deep lake, 5=shallow lake)
  real(r8), INTENT(out) :: &!--- Soil parameters
        albsol,            &!soil albedo for different coloured soils [-]
        csol  (1:nl_soil), &!heat capacity of soil solids [J/(m3 K)]
        porsl (1:nl_soil), &!fraction of soil that is voids [-]
        phi0  (1:nl_soil), &!minimum soil suction [mm]
        bsw   (1:nl_soil), &!clapp and hornbereger "b" parameter [-]
        dkmg  (1:nl_soil), &!thermal conductivity of soil minerals [W/m-K]
        dksatu(1:nl_soil), &!thermal conductivity of saturated soil [W/m-K]
        dkdry (1:nl_soil), &!thermal conductivity for dry soil  [W/(m-K)]
        hksati(1:nl_soil)   !hydraulic conductivity at saturation [mm h2o/s]
  real(r8), INTENT(out) :: &!--- Vegetation static parameters
        z0m              , &!aerodynamic roughness length [m]
        displa           , &!displacement height [m]
        sqrtdi           , &!inverse sqrt of leaf dimension [m**-0.5]
        effcon           , &!quantum efficiency of RuBP regeneration (mol CO2/mol quanta)
        vmax25           , &!maximum carboxylation rate at 25 C at canopy top (mol CO2/m2s)
        shti             , &!slope of high temperature inhibition function     (s1)
        hhti             , &!1/2 point of high temperature inhibition function (s2)
        slti             , &!slope of low temperature inhibition function      (s3)
        hlti             , &!1/2 point of low temperature inhibition function  (s4)
        trda             , &!temperature coefficient in gs-a model             (s5)
        trdm             , &!temperature coefficient in gs-a model             (s6)
        trop             , &!temperature coefficient in gs-a model         (273+25)
        gradm            , &!conductance-photosynthesis slope parameter
        binter           , &!conductance-photosynthesis intercep
        extkn            , &!coefficient of leaf nitrogen allocation
        chil             , &!leaf angle distribution factor
        ref(2,2)         , &!leaf reflectance (iw=iband, il=life and dead)
        tran(2,2)        , &!leaf transmittance (iw=iband, il=life and dead)
        rootfr(1:nl_soil)   !fraction of roots in each soil layer
  real(r8), INTENT(out) :: &!--- Initialize TUNABLE constants
        zlnd             , &!Roughness length for soil [m]
        zsno             , &!Roughness length for snow [m]
        csoilc           , &!Drag coefficient for soil under canopy [-]
        dewmx            , &!maximum dew
        wtfact           , &!Fraction of model area with high water table
        capr             , &!Tuning factor to turn first layer T into surface T
        cnfac            , &!Crank Nicholson factor between 0 and 1
        ssi              , &!Irreducible water saturation of snow
        wimp             , &!Water impremeable if porosity less than wimp
        pondmx           , &!Ponding depth (mm)
        smpmax           , &!Wilting point potential in mm
        smpmin           , &!Restriction for min of soil poten. (mm)
        trsmx0           , &!Max transpiration for moist soil+100% veg. [mm/s]
	tcrit               !critical temp. to determine rain or snow
!--------------------------- Local variables
  integer i, j              !indices
  integer idlak             !index (1=deep lake, 0=shallow)
  integer lsa               !my sensitivity parameters input file number
  character(LEN=256) :: fsa !my sensitivity parameters input file name
  real(r8),allocatable :: senparameters(:)
  real(r8) bd(1:nl_soil)    !bulk density of dry soil material [kg/m^3]
  real(r8) dkm(1:nl_soil)   !
  real(r8) dzlak(1:nl_soil) !
  real(r8) zlak(1:nl_soil)  !
  real(r8) zsoih(0:nl_soil) !interface level below a zsoi level [m]
!--------------------------- Data block
! Soil albedo for different colored soils (saturated soil, visible beam) [-]
  real(r8), dimension(8)  :: &
  solour = (/0.12, 0.11, 0.10, 0.09, 0.08, 0.07, 0.06, 0.05/)
!---------------------------
!  USGS Land Use/Land Cover System Legend
! 0  Ocean
! 1  Urban and Built-Up Land
! 2  Dryland Cropland and Pasture
! 3  Irrigated Cropland and Pasture
! 4  Mixed Dryland/Irrigated Cropland and Pasture
! 5  Cropland/Grassland Mosaic
! 6  Cropland/Woodland Mosaic
! 7  Grassland
! 8  Shrubland
! 9  Mixed Shrubland/Grassland
!10  Savanna
!11  Deciduous Broadleaf Forest
!12  Deciduous Needleleaf Forest
!13  Evergreen Broadleaf Forest
!14  Evergreen Needleleaf Forest
!15  Mixed Forest
!16  Inland Water
!17  Herbaceous Wetland
!18  Wooded Wetland
!19  Barren or Sparsely Vegetated
!20  Herbaceous Tundra
!21  Wooded Tundra
!22  Mixed Tundra
!23  Bare Ground Tundra
!24  Snow or Ice
  integer, parameter :: numUSGS = 24
  real(r8), dimension(numUSGS) :: &!
  z0m_usgs        , &!roughness
  displa_usgs     , &!zero-plane-distance
  sqrtdi_usgs     , &!inverse sqrt of leaf dimension [m**-0.5]
  chil_usgs       , &!leaf angle distribution factor
  ref_s_usgs      , &!leaf reflectance
  ref_sd_usgs     , &!leaf reflectance
  ref_l_usgs      , &!leaf reflectance
  ref_ld_usgs     , &!leaf reflectance
  tran_s_usgs     , &!leaf transmittance
  tran_sd_usgs    , &!leaf transmittance
  tran_l_usgs     , &!leaf transmittance
  tran_ld_usgs    , &!leaf transmittance
  vmax0_usgs      , &!maximum carboxylation rate at 25 C at canopy top
  effcon_usgs     , &!quantum efficiency
  gradm_usgs      , &!conductance-photosynthesis slope parameter
  binter_usgs     , &!conductance-photosynthesis intercept
  respcp_usgs     , &!respiration fraction
  shti_usgs       , &!slope of high temperature inhibition function (s1)
  slti_usgs       , &!slope of low temperature inhibition function (s3)
  trda_usgs       , &!temperature coefficient in gs-a model (s5)
  trdm_usgs       , &!temperature coefficient in gs-a model (s6)
  trop_usgs       , &!temperature coefficient in gs-a model (273.16+25)
  hhti_usgs       , &!1/2 point of high temperature inhibition function (s2)
  hlti_usgs       , &!1/2 point of low temperature inhibition function (s4)
  extkn_usgs      , &!coefficient of leaf nitrogen allocation
  d50_usgs        , &!depth at 50% roots
  d95_usgs        , &!depth at 95% roots
  beta_usgs          !coefficient of root profile
! fenpei senparameters
  allocate(senparameters(6))
!----------------------------------------------------------------------
!my sensitivity parameters
!define my sensitivity parameters input file number and name
  lsa=21
  fsa='/group_homes/lzu_public/home/u120220909911/Summer/Last/LE/input_step.txt'
!open my sensitivity parameters file
  open(unit=lsa,file=fsa,status='old',form='formatted',action='read')
!read in the parameters value to senparameters matrix
  CALL SAparametersread(lsa,senparameters)
!test the senparameters
!  print*,'iniTimeConst.F90-test-senparameters-1',senparameters(1)
!  print*,'iniTimeConst.F90-test-senparameters-2',senparameters(2)
  z0m_usgs    =(/0.100,  0.100,  0.100,  0.100,  0.100,  0.100,  0.100,  0.050,&
                 0.050,  0.100,  2.000,  1.700,  3.500,  1.700,  2.000,  0.100,&
                 0.100,  3.500,  0.050,  0.100,  0.100,  0.100,  0.100,  0.100/)
  displa_usgs =(/0.667,  0.667,  0.667,  0.667,  0.667,  0.667,  0.667,  0.333,&
                 0.333,  0.667, 13.333, 11.333, 23.333, 11.333, 13.333,  0.667,&
                 0.667, 23.333,  0.333,  0.667,  0.667,  0.667,  0.667,  0.667/)
  sqrtdi_usgs(:)=5.0
  chil_usgs  =(/-0.300, -0.300, -0.300, -0.300, -0.300, -0.300, -0.300,  0.010,&
                 0.010, -0.300,  0.250,  0.010,  0.100,  0.010,  0.125, -0.300,&
                -0.300,  0.100,  0.010, -0.300, -0.300, -0.300, -0.300, -0.300/)
  ref_s_usgs  =(/0.105,  0.105,  0.105,  0.105,  0.105,  0.105,  0.105,  0.100,&
                 0.100,  0.105,  0.100,  0.070,  0.100,  0.070,  0.070,  0.105,&
                 0.105,  0.100,  0.100,  0.105,  0.105,  0.105,  0.105,  0.105/)
  ref_sd_usgs =(/0.360,  0.360,  0.360,  0.360,  0.360,  0.360,  0.360,  0.160,&
                 0.160,  0.360,  0.160,  0.160,  0.160,  0.160,  0.160,  0.360,&
                 0.360,  0.160,  0.160,  0.360,  0.360,  0.360,  0.360,  0.360/)
  ref_l_usgs  =(/0.580,  0.580,  0.580,  0.580,  0.580,  0.580,  0.580,  0.450,&
                 0.450,  0.580,  0.450,  0.350,  0.450,  0.350,  0.400,  0.580,&
                 0.580,  0.450,  0.450,  0.580,  0.580,  0.580,  0.580,  0.580/)
  ref_ld_usgs =(/0.580,  0.580,  0.580,  0.580,  0.580,  0.580,  0.580,  0.390,&
                 0.390,  0.580,  0.390,  0.390,  0.390,  0.390,  0.390,  0.580,&
                 0.580,  0.390,  0.390,  0.580,  0.580,  0.580,  0.580,  0.580/)
  tran_s_usgs =(/0.070,  0.070,  0.070,  0.070,  0.070,  0.070,  0.070,  0.070,&
                 0.070,  0.070,  0.050,  0.050,  0.050,  0.050,  0.050,  0.070,&
                 0.070,  0.050,  0.070,  0.070,  0.070,  0.070,  0.070,  0.070/)
  tran_sd_usgs=(/0.220,  0.220,  0.220,  0.220,  0.220,  0.220,  0.220,  0.001,&
                 0.001,  0.220,  0.001,  0.001,  0.001,  0.001,  0.001,  0.220,&
                 0.220,  0.001,  0.001,  0.220,  0.220,  0.220,  0.220,  0.220/)
  tran_l_usgs =(/0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,&
                 0.250,  0.250,  0.250,  0.100,  0.250,  0.100,  0.150,  0.250,&
                 0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250,  0.250/)
  tran_ld_usgs=(/0.380,  0.380,  0.380,  0.380,  0.380,  0.380,  0.380,  0.001,&
                 0.001,  0.380,  0.001,  0.001,  0.001,  0.001,  0.001,  0.380,&
                 0.380,  0.001,  0.001,  0.380,  0.380,  0.380,  0.380,  0.380/)
  vmax0_usgs( 1: 7)=100.0;   vmax0_usgs( 8: 9)=60.0  
  vmax0_usgs(10:13)=100.0;   vmax0_usgs   (14)=60.0 
  vmax0_usgs   (15)=80.0;    vmax0_usgs(16:18)=100.0 
  vmax0_usgs   (19)=60.0;    vmax0_usgs(20:24)=30.0  
  effcon_usgs(1:19)=0.08;  effcon_usgs(20:24)=0.05
  gradm_usgs (1:19)=9.0;   gradm_usgs (20:24)=4.0
  binter_usgs(1:19)=0.01;  binter_usgs(20:24)=0.04
  respcp_usgs(1:19)=0.015; respcp_usgs(20:24)=0.025
  shti_usgs(:)=0.3
  slti_usgs(:)=0.2
  trda_usgs(:)=1.3
  trdm_usgs(:)=328.0
  trop_usgs(:)=298.0
  hhti_usgs=(/308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 313.0,&
              313.0, 308.0, 311.0, 303.0, 313.0, 303.0, 307.0, 308.0,&
              308.0, 313.0, 313.0, 313.0, 313.0, 313.0, 313.0, 308.0/)
  hlti_usgs=(/281.0, 281.0, 281.0, 281.0, 281.0, 281.0, 281.0, 283.0,&
              283.0, 281.0, 283.0, 278.0, 288.0, 278.0, 281.0, 281.0,&
              281.0, 288.0, 283.0, 288.0, 288.0, 288.0, 288.0, 281.0/)
  extkn_usgs(:)=0.5
  d50_usgs  =(/23.0,  21.0,  23.0,  22.0,  15.7,  19.0,   9.3,  47.0,&
	       28.2,  21.7,  16.0,  16.0,  15.0,  15.0,  15.5,   1.0,&
	        9.3,  15.5,  27.0,   9.0,   9.0,   9.0,   9.0,   1.0/)
  d95_usgs =(/121.0, 104.0, 121.0, 112.5,  80.8, 103.8,  49.0, 302.0,&
	      175.5,  99.3,  95.0,  95.0,  91.0,  91.0,  93.0,   1.0,&
               49.0,  93.0, 112.0,  29.0,  29.0,  29.0,  29.0,   1.0/)
  beta_usgs=(/-1.757, -1.835, -1.757, -1.796, -1.577, -1.738,&
              -1.359, -3.245, -2.302, -1.654, -1.681, -1.681,&
              -1.632, -1.632, -1.656, -1.000, -1.359, -1.656,&
              -2.051, -2.621, -2.621, -2.621, -2.621, -1.000 /)
!-----------------------------------------------------------------------
! land water type for USGS classification
!-----------------------------------------------------------------------
         i=ivt
                         itypwat=0  ! soil
      if(i==1)           itypwat=1  ! urban and built-up
      if(i==17.or.i==18) itypwat=2  ! wetland
      if(i==24)          itypwat=3  ! land ice
      if(i==16)          itypwat=4  ! deep lake
      if(i==0)           itypwat=99 ! ocean
!-----------------------------------------------------------------------
! soil layer thickness, depths (m)
!-----------------------------------------------------------------------
! ------ Non Lake ------ !
      if(itypwat<4)then
         do j = 1, nl_soil
           zsoi(j) = 0.025*(exp(0.5*(j-0.5))-1.)  !node depths
         end do
         dzsoi(1)  = 0.5*(zsoi(1)+zsoi(2))        !=zsoih(1)
         dzsoi(nl_soil)= zsoi(nl_soil)-zsoi(nl_soil-1)
         do j = 2,nl_soil-1
           dzsoi(j)= 0.5*(zsoi(j+1)-zsoi(j-1))    !thickness b/n two interfaces
         end do
!         print*,'iniTimeConst.F90-test-dzsoi',dzsoi
         zsoih(0)   = 0.
         zsoih(nl_soil) = zsoi(nl_soil) + 0.5*dzsoi(nl_soil)
         do j = 1, nl_soil-1
            zsoih(j)= 0.5*(zsoi(j)+zsoi(j+1))     !interface depths
         enddo
!         print*,'iniTimeConst.F90-test-zsoih',zsoih
! ------ Lake ------ !
      else if (itypwat<=5)then
         idlak = 1                                !assumed all lakes are deep lake
         if(idlak == 1) then                                                 
            dzlak = (/1., 2., 3., 4., 5., 7., 7., 7., 7., 7./)
             zlak = (/0.5, 1.5, 4.5, 8.0, 12.5, 18.5, 25.5, 32.5, 39.5, 46.5/)
         else
            dzlak = (/.25, 0.5, 0.75, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/)   
             zlak = (/ 0.125,  0.5,  1.125,  2.,  3.,  4.,  5.,  6.,  7.,  8./)
         end if
         zsoih(0) = 0.
         do j = 1, nl_soil
            dzsoi(j) = dzlak(j)
             zsoi(j) =  zlak(j)
             zsoih(j) = zsoih(j-1)+dzsoi(j)
         enddo
! ------ Ocean ------ !
      else
            dzsoi(:) = -999.9
             zsoi(:) = -999.9
            zsoih(:) = -999.9
      endif
!-----------------------------------------------------------------------
! soil thermal and hydraulic properties
!-----------------------------------------------------------------------
! saturated soil albedo for visible beam
            albsol = solour(isc)
! soil thermal and hydraulic properties
      if (itypwat<2)then  ! not wetland, glacier and lake
         do j = 1, nl_soil
	    if(zsoi(j)<=rockdep)then  !NON ROCK
! for all organic case or data missing, assigned to "loam"
            if(sand(j)*clay(j).lt.0.01)then
               sand(j) = 43.
               clay(j) = 18.
            endif
!          porsl(j) = 0.489 - 0.00126*sand(j)
           if(j<=5)then
                porsl(j)=0.489 - 0.00126*sand(j)
           else
                porsl(j)=senparameters(6)
           endif
!           print*,'iniTimeConst.F90-test-porsl-3',porsl(j)
!           phi0(j) = 10. * ( 10.**(1.88-0.0131*sand(j)) )
           if(j<=5)then
               phi0(j)=senparameters(4)
           else
               phi0(j) = 10. * ( 10.**(1.88-0.0131*sand(j)) )
           endif
!           phi0(j)=senparameters(5)
!           print*,'iniTimeConst.F90-test-phi0-4',phi0(j)
            bsw(j) = 2.91 + 0.159*clay(j)
!          bsw(j)=senparameters(3)
!          print*,'iniTimeConst.F90-test-bsw-6',bsw(j)
         hksati(j) = 0.0070556 * ( 10.**(-0.884+0.0153*sand(j)) ) ! mm/s
!          hksati(j)=senparameters(4)
!          print*,'iniTimeConst.F90-test-hksati-2',hksati(j)
             bd(j) = (1.- porsl(j))*2.7e3
           csol(j) = (2.128*sand(j)+2.385*clay(j))/(sand(j)+clay(j))*1.e6  ! J/(m3 K)
           dkm(j)  = (8.80*sand(j)+2.92*clay(j)) / (sand(j)+clay(j))       ! W/(m K)
           dkmg(j) = dkm(j)**(1.-porsl(j))
         dksatu(j) = dkmg(j)*0.57**porsl(j)
          dkdry(j) = (.135*bd(j) + 64.7) / (2.7e3 - 0.947*bd(j))
              if(zsoih(j)>rockdep)then
                 porsl(j) = porsl(j)*(rockdep-zsoih(j-1))/dzsoi(j)
              endif
	    else                      !BEDROCK
          porsl(j:) = 0.
           phi0(j:) = 1.5e-5
            bsw(j:) = 0
         hksati(j:) = 0.
             bd(j:) = 2.7e3
           csol(j:) = 2700.*750. !J/(m3 K)
           dkmg(j:) = 1.0        ! not used
         dksatu(j:) = 2.9        ! W/(m K)
          dkdry(j:) = 2.9        ! W/(m K) [J.R. Garratt 1992, pp291]
            exit
	    endif
         end do
      else                ! wetland, glacier, lake and ocean
          porsl(:) = 1.0
           phi0(:) = 0.0
            bsw(:) = 0.0
           csol(:) = 4.186e06
           dkmg(:) = 1.0
         dksatu(:) = 0.6
          dkdry(:) = 0.6
         hksati(:) = 0.0
      endif
!-----------------------------------------------------------------------
! vegetation static parameters
! the values for glacier and lake are assigned arbitrarily (not used)
!-----------------------------------------------------------------------
           i = ivt
    if( i > 0 )then              ! land grids
         z0m = z0m_usgs   (i)
!      z0m=senparameters(14)
!      print*,'iniTimeConst.F90-test-z0m-36',z0m
     displa = displa_usgs(i)
      sqrtdi = sqrtdi_usgs(i)  
!      sqrtdi=senparameters(28)
!      print*,'iniTimeConst.F90-test-sqrtdi-17',sqrtdi
!      vmax25 = vmax0_usgs(i)*1.e-6
       vmax25=senparameters(1)*1.e-6
!      print*,'iniTimeConst.F90-test-vmax25-19',vmax25
!      effcon = effcon_usgs(i)
       effcon=senparameters(2)
!       print*,'iniTimeConst.F90-test-effcon-18',effcon
      slti   =   slti_usgs(i) 
!      slti=senparameters(24)
!      print*,'iniTimeConst.F90-test-slti-14',slti
      hlti   =   hlti_usgs(i)
!      hlti=senparameters(25)
!      print*,'iniTimeConst.F90-test-hlti-15',hlti
      shti   =   shti_usgs(i)
!      shti=senparameters(26)
!      print*,'iniTimeConst.F90-test-shti-16',shti
      hhti   =   hhti_usgs(i)       
!      hhti=senparameters(27)
!      print*,'iniTimeConst.F90-test-hhti-20',hhti
     trda   =   trda_usgs(i)  
!      trda=senparameters(31)
!      print*,'iniTimeConst.F90-test-trda-21',trda
      trdm   =   trdm_usgs(i)       
!      trdm=senparameters(32)
!      print*,'iniTimeConst.F90-test-trdm-22',trdm
     trop   =   trop_usgs(i)  
!      trop=senparameters(33)
!      print*,'iniTimeConst.F90-test-trop-23',trop
!      gradm  =  gradm_usgs(i)
      gradm=senparameters(3)
!      print*,'iniTimeConst.F90-test-gradm-24',gradm
!      binter = binter_usgs(i)
      binter=senparameters(5)
!      print*,'iniTimeConst.F90-test-binter-25',binter
      extkn  =  extkn_usgs(i)
!      extkn=senparameters(33)
!      print*,'iniTimeConst.F90-test-extkn-26',extkn
      chil = chil_usgs(i) 
!      chil=senparameters(15)
!     print*,'iniTimeConst.F90-test-chil-27',chil
      ref(1,1) = ref_s_usgs(i) 
!     ref(1,1)=senparameters(16)
!     print*,'iniTimeConst.F90-test-ref(1,1)-28',ref(1,1)
      ref(2,1) = ref_l_usgs(i)
!     ref(2,1)=senparameters(18)
!     print*,'iniTimeConst.F90-test-ref(2,1)-30',ref(2,1)
      ref(1,2) = ref_sd_usgs(i)       
!     ref(1,2)=senparameters(17)
!     print*,'iniTimeConst.F90-test-ref(1,2)-29',ref(1,2)
      ref(2,2) = ref_ld_usgs(i)  
!     ref(2,2)=senparameters(19)
!     print*,'iniTimeConst.F90-test-ref(2,2)-31',ref(2,2)
      tran(1,1) = tran_s_usgs(i)      
!     tran(1,1)=senparameters(20)
!     print*,'iniTimeConst.F90-test-tran(1,1)-32',tran(1,1)
      tran(2,1) = tran_l_usgs(i)      
!     tran(2,1)=senparameters(22)
!     print*,'iniTimeConst.F90-test-tran(2,1)-34',tran(2,1)
      tran(1,2) = tran_sd_usgs(i) 
!     tran(1,2)=senparameters(21)
!     print*,'iniTimeConst.F90-test-tran(1,2)-33',tran(1,2)
      tran(2,2) = tran_ld_usgs(i) 
!     tran(2,2)=senparameters(23)
!     print*,'iniTimeConst.F90-test-tran(2,2)-35',tran(2,2)
! The definition of global root distribution is based on
! Schenk and Jackson, 2002: The Global Biogeography of Roots.
! Ecological Monagraph 72(3): 311-328.
      if(itypwat>=3)then  !glacier or lake
      rootfr(:)=0.
      else
      rootfr(1)=1./(1.+(zsoih(  1)*100./d50_usgs(i))**beta_usgs(i)) 
      rootfr(nl_soil)=1.-1./(1.+(zsoih(nl_soil-1)*100./d50_usgs(i))**beta_usgs(i)) 
      do j=2,nl_soil-1
      rootfr(j)=1./(1.+(zsoih(j  )*100./d50_usgs(i))**beta_usgs(i)) &
               -1./(1.+(zsoih(j-1)*100./d50_usgs(i))**beta_usgs(i))
      enddo
      endif
    else                         ! ocean grids
      z0m       = -999.9
      displa    = -999.9 
      sqrtdi    = -999.9
      vmax25    = -999.9
      effcon    = -999.9
      slti      = -999.9
      hlti      = -999.9
      shti      = -999.9
      hhti      = -999.9
      trda      = -999.9
      trdm      = -999.9
      trop      = -999.9
      gradm     = -999.9
      binter    = -999.9
      extkn     = -999.9
      chil      = -999.9
      ref(:,:)  = -999.9
      tran(:,:) = -999.9
      rootfr(:) = -999.9
    endif
!-----------------------------------------------------------------------
! Initialize TUNABLE constants
!-----------------------------------------------------------------------
      zlnd   = 0.01    !Roughness length for soil [m]
!      zlnd=senparameters(1)
!      print*,'iniTimeConst.F90-test-zlnd-8',zlnd
      zsno   = 0.0024  !Roughness length for snow [m]
!      zsno=senparameters(10)
!      print*,'iniTimeConst.F90-test-zsno-11',zsno
      csoilc = 0.004   !Drag coefficient for soil under canopy [-]
!      csoilc=senparameters(9)
!      print*,'iniTimeConst.F90-test-csoilc-10',csoilc
      dewmx  = 0.1     !maximum dew
!      dewmx=senparameters(13)
!      print*,'iniTimeConst.F90-test-dewmx-1',dewmx
      wtfact = 0.3     !Fraction of model area with high water table
!     wtfact=senparameters(6)
!     print*,'iniTimeConst.F90-test-wtfact-5',wtfact
     capr   = 0.34    !Tuning factor to turn first layer T into surface T
!      capr=senparameters(11)
!     print*,'iniTimeConst.F90-test-capr-12',capr
      cnfac  = 0.5     !Crank Nicholson factor between 0 and 1
!     cnfac=senparameters(12)
!     print*,'iniTimeConst.F90-test-cnfac-13',cnfac
     ssi    = 0.033   !Irreducible water saturation of snow
!     ssi=senparameters(40)
!     print*,'iniTimeConst.F90-test-ssi-37',ssi
      wimp   = 0.05    !Water impremeable if porosity less than wimp
!      wimp=senparameters(7)
!      print*,'iniTimeConst.F90-test-wimp-7',wimp
      pondmx = 10.0    !Ponding depth (mm)
!      pondmx=senparameters(8)
!      print*,'iniTimeConst.F90-test-pondmx-9',pondmx
      smpmax = -1.5e5  !Wilting point potential in mm
!      smpmax=senparameters(34)
!      print*,'iniTimeConst.F90-test-smpmax-38',smpmax
      smpmin = -1.e8   !Restriction for min of soil poten. (mm)
!      smpmin=senparameters(35)
!      print*,'iniTimeConst.F90-test-smpmin-39',smpmin
      trsmx0 = 2.e-4   !Max transpiration for moist soil+100% veg. [mm/s]
!      trsmx0=senparameters(36)
!      print*,'iniTimeConst.F90-test-trsmx0-40',trsmx0
      tcrit  = 0.      !critical temp. to determine rain or snow
 end subroutine iniTimeConst
