 subroutine initialize (start_yr,start_jday,start_sec,greenwich,&
                        lon_points,lat_points,maxpatch,nl_soil,maxsnl,&
                        lusrf,&
                        lhistTimeConst,lhistTimeConst2,lhistTimeVar,&
                        lhistTimeVar2,nftune,nfcon,nfvar,nfldv,numpatch,&
                        fldxy)
! ======================================================================
! initialization routine for land surface model.
!
! original author : Yongjiu Dai, 09/15/1999; 08/30/2002
! ======================================================================
   use precision
   use phycon_module
   implicit none
! ----------------------------------------------------------------------
   integer, INTENT(in) :: start_yr       ! starting date for run in year
   integer, INTENT(in) :: start_jday     ! starting date for run in julian day
   integer, INTENT(in) :: start_sec      ! starting time of day for run in seconds
   logical, INTENT(in) :: greenwich      ! true: greenwich time, false: local time
   integer, INTENT(in) :: lon_points     ! number of longitude points on model grid
   integer, INTENT(in) :: lat_points     ! number of latitude points on model grid
   integer, INTENT(in) :: maxpatch       ! max number of patches in a grid
   integer, INTENT(in) :: nl_soil        ! number of soil layers
   integer, INTENT(in) :: maxsnl         ! max number of snow layers
   integer, INTENT(in) :: lusrf          ! logical unit number of surface data
   integer, INTENT(in) :: lhistTimeConst ! logical unit number of restart time-invariant file
   integer, INTENT(in) :: lhistTimeConst2 ! logical unit number of mine
   integer, INTENT(in) :: lhistTimeVar   ! logical unit number of restart time-varying file
   integer, INTENT(in) :: lhistTimeVar2
   integer, INTENT(in) :: nftune         ! number of clm tunable constants
   integer, INTENT(in) :: nfcon          ! number of time constant variables
   integer, INTENT(in) :: nfvar          ! number of time varying variables
   integer, INTENT(in) :: nfldv          ! number of time varying variables
   integer, INTENT(out) :: numpatch      ! total number of patches of grids
! ------------------------ local variables -----------------------------
! surface classification and soil information
  real(r8) latixy  (lon_points,lat_points)          ! latitude in radians
  real(r8) longxy  (lon_points,lat_points)          ! longitude in radians
  real(r8) sand2d  (lon_points,lat_points,nl_soil)  ! percentage of sand
  real(r8) clay2d  (lon_points,lat_points,nl_soil)  ! percentage of clay
  real(r8) rock2d  (lon_points,lat_points)          ! depth to bed rock
  integer  soic2d  (lon_points,lat_points)          ! soil color
  integer  surf2d  (lon_points,lat_points,maxpatch) ! land cover type
  real(r8) fpatch2d(lon_points,lat_points,maxpatch) !subgrid weights
  real(r8) zlnd    ! roughness length for soil [m]
  real(r8) zsno    ! roughness length for snow [m]
  real(r8) csoilc  ! drag coefficient for soil under canopy [-]
  real(r8) dewmx   ! maximum dew
  real(r8) wtfact  ! fraction of model area with high water table
  real(r8) capr    ! tuning factor to turn first layer T into surface T
  real(r8) cnfac   ! Crank Nicholson factor between 0 and 1
  real(r8) ssi     ! irreducible water saturation of snow
  real(r8) wimp    ! water impremeable if porosity less than wimp
  real(r8) pondmx  ! ponding depth (mm)
  real(r8) smpmax  ! wilting point potential in mm
  real(r8) smpmin  ! restriction for min of soil poten. (mm)
  real(r8) trsmx0  ! max transpiration for moist soil+100% veg.  [mm/s]
  real(r8) tcrit   ! critical temp. to determine rain or snow
  real(r8), allocatable :: rockdep (:) ! depth to bedrock
  real(r8), allocatable :: sand  (:,:) ! percent sand
  real(r8), allocatable :: clay  (:,:) ! percent clay
  integer,  allocatable :: isc     (:) ! color classes for soil albedos
  real(r8), allocatable :: work    (:) !
  integer,  allocatable :: ixy_patch(:)  ! patch longitude index
  integer,  allocatable :: jxy_patch(:)  ! patch latitude index
  integer,  allocatable :: mxy_patch(:)  ! patch subgrid index of lnd point
  real(r8), allocatable :: wtxy_patch(:) ! patch weight
  real(r8), allocatable :: dlat    (:) ! latitude in radians
  real(r8), allocatable :: dlon    (:) ! longitude in radians
  integer , allocatable :: itypwat (:) ! land water type
  integer , allocatable :: ivt     (:) ! land cover type of classification of USGS etc.
  real(r8), allocatable :: albsol  (:) ! soil albedo for different coloured soils [-]
  real(r8), allocatable :: csol  (:,:) ! heat capacity of soil solids [J/(m3 K)]
  real(r8), allocatable :: porsl (:,:) ! fraction of soil that is voids [-]
  real(r8), allocatable :: phi0  (:,:) ! minimum soil suction [mm]
  real(r8), allocatable :: bsw   (:,:) ! clapp and hornbereger "b" parameter [-]
  real(r8), allocatable :: dkmg  (:,:) ! thermal conductivity of soil minerals [W/m-K]
  real(r8), allocatable :: dksatu(:,:) ! thermal conductivity of saturated soil [W/m-K]
  real(r8), allocatable :: dkdry (:,:) ! thermal conductivity for dry soil  [W/(m-K)]
  real(r8), allocatable :: hksati(:,:) ! hydraulic conductivity at saturation [mm h2o/s]
  real(r8), allocatable :: z0m     (:) ! aerodynamic roughness length [m]
  real(r8), allocatable :: displa  (:) ! displacement height [m]
  real(r8), allocatable :: sqrtdi  (:) ! inverse sqrt of leaf dimension [m**-0.5]
  real(r8), allocatable :: effcon  (:) ! quantum efficiency of RuBP regeneration
  real(r8), allocatable :: vmax25  (:) ! maximum carboxylation rate at 25 C at canopy top
  real(r8), allocatable :: slti    (:) ! s3: slope of low temperature inhibition function
  real(r8), allocatable :: hlti    (:) ! s4: 1/2 point of low temperature inhibition function
  real(r8), allocatable :: shti    (:) ! s1: slope of high temperature inhibition function
  real(r8), allocatable :: hhti    (:) ! s2: 1/2 point of high temperature inhibition function
  real(r8), allocatable :: trda    (:) ! s5: temperature coefficient in gs-a model
  real(r8), allocatable :: trdm    (:) ! s6: temperature coefficient in gs-a model
  real(r8), allocatable :: trop    (:) ! temperature coefficient in gs-a model
  real(r8), allocatable :: gradm   (:) ! conductance-photosynthesis slope parameter
  real(r8), allocatable :: binter  (:) ! conductance-photosynthesis intercep
  real(r8), allocatable :: extkn   (:) ! coefficient of leaf nitrogen allocation
  real(r8), allocatable :: chil    (:) ! leaf angle distribution factor
  real(r8), allocatable :: ref (:,:,:) ! leaf reflectance (iw=iband, il=life and dead)
  real(r8), allocatable :: tran(:,:,:) ! leaf transmittance (iw=iband, il=life and dead)
  real(r8), allocatable :: rootfr(:,:) ! fraction of roots in each soil layer
! Time-varying state variables which reaquired by restart run
  real(r8), allocatable :: z      (:,:) ! node depth [m]
  real(r8), allocatable :: dz     (:,:) ! interface depth [m]
  real(r8), allocatable :: tss    (:,:) ! soil temperature [K]
  real(r8), allocatable :: wliq   (:,:) ! liquid water in layers [kg/m2]
  real(r8), allocatable :: wice   (:,:) ! ice lens in layers [kg/m2]
  real(r8), allocatable :: tg       (:) ! ground surface temperature [K]
  real(r8), allocatable :: tlsun    (:) ! sunlit leaf temperature [K]
  real(r8), allocatable :: tlsha    (:) ! shaded leaf temperature [K]
  real(r8), allocatable :: ldew     (:) ! depth of water on foliage [mm]
  real(r8), allocatable :: sag      (:) ! non dimensional snow age [-]
  real(r8), allocatable :: scv      (:) ! snow cover, water equivalent [mm]
  real(r8), allocatable :: snowdp   (:) ! snow depth [meter]
  real(r8), allocatable :: fveg     (:) ! fraction of vegetation cover
  real(r8), allocatable :: fsno     (:) ! fraction of snow cover on ground
  real(r8), allocatable :: sigf     (:) ! fraction of veg cover, excluding snow-covered veg [-]
  real(r8), allocatable :: green    (:) ! leaf greenness
  real(r8), allocatable :: lai      (:) ! leaf area index
  real(r8), allocatable :: sai      (:) ! stem area index
  real(r8), allocatable :: coszen   (:) ! cosine of solar zenith angle
  real(r8), allocatable :: albg (:,:,:) ! albedo, ground [-]
  real(r8), allocatable :: albv (:,:,:) ! albedo, vegetation [-]
  real(r8), allocatable :: alb  (:,:,:) ! averaged albedo [-]
  real(r8), allocatable :: ssun (:,:,:) ! sunlit canopy absorption for solar radiation (0-1)
  real(r8), allocatable :: ssha (:,:,:) ! shaded canopy absorption for solar radiation (0-1)
  real(r8), allocatable :: thermk   (:) ! canopy gap fraction for tir radiation
  real(r8), allocatable :: extkb    (:) ! (k, g(mu)/mu) direct solar extinction coefficient
  real(r8), allocatable :: extkd    (:) ! diffuse and scattered diffuse PAR extinction coefficient
  real(r8), allocatable :: trad     (:) ! radiative temperature of surface [K]
  real(r8), allocatable :: tref     (:) ! 2 m height air temperature [kelvin]
  real(r8), allocatable :: qref     (:) ! 2 m height air specific humidity
  real(r8), allocatable :: rst      (:) ! canopy stomatal resistance (s/m)
  real(r8), allocatable :: emis     (:) ! averaged bulk surface emissivity
  real(r8), allocatable :: z0ma     (:) ! effective roughness [m]
  real(r8), allocatable :: zol      (:) ! dimensionless height (z/L) used in Monin-Obukhov theory
  real(r8), allocatable :: rib      (:) ! bulk Richardson number in surface layer
  real(r8), allocatable :: ustar    (:) ! u* in similarity theory [m/s]
  real(r8), allocatable :: qstar    (:) ! q* in similarity theory [kg/kg]
  real(r8), allocatable :: tstar    (:) ! t* in similarity theory [K]
  real(r8), allocatable :: fm       (:) ! integral of profile function for momentum
  real(r8), allocatable :: fh       (:) ! integral of profile function for heat
  real(r8), allocatable :: fq       (:) ! integral of profile function for moisture
  integer   idate(3)                    ! calendar (year, julian day, seconds)
  integer   numpatch_lat(lat_points)    ! number of patches of grids at lon. strip
  real(r8)  ftune(nftune)               ! clm tunable constants
  real(r8), allocatable :: fcon(:,:)    ! time constant variables
  real(r8), allocatable :: fvar(:,:)    ! time varying variables
  real(r8), intent(out) :: fldxy(lon_points,lat_points,nfldv)
  integer   year                        ! current year of model run
  integer   jday                        ! current julian day of model run
  integer   msec                        ! current seconds of model run (0-86400)
  real(r8)  pi                          ! pie
  real(r8)  a                           !
  real(r8)  calday                      ! Julian cal day (1.xx to 365.xx)
  real(r8)  orb_coszen                  ! cosine of the solar zenith angle
  integer   i,j,k,l,m,npatch            ! indices
  integer   jm,lb,ub                    ! indices
! ----------------------------------------------------------------------
! [1] READ IN LAND INFORMATION
! read time-invariant boundary data on [lon_points] x [lat_points] grid.
!    o first [lat_points] values: number of longitude points for each latitude.
!      this allows for variable longitudinal resolution for each latitude
! remaining data is for each grid cell:
!    o 1st : latitude at center of grid cell (degrees)
!    o 2th : longitude at center of grid cell (degrees)
!    o 3th : soil color (1 to 8) for use with soil albedos
!    o 4th : depth to bed rock
!    o 5th : soil texture, %sand, for thermal and hydraulic properties
!    o 6th : soil texture, %clay, for thermal and hydraulic properties
!    o 7th : surface type, for use as multiple subgrid point
!    o 8th : subgrid weight
!
! ------------------ note ------------------
! the model required the same longitudinal resolution for
! all latitude strip. For the variable longitudinal resolution
! cases, please assign the surface types and soil character:
! 0 to ocean grids and -999 to the land grids
! which are not included in the calcultion.
! ----------------------------------------------------------------------
      do j = 1, lat_points
         do i = 1, lon_points
            read(lusrf)&
            latixy(i,j),longxy(i,j),&
            soic2d(i,j),rock2d(i,j),&
            sand2d(i,j,1:nl_soil), clay2d(i,j,1:nl_soil),&
            surf2d(i,j,1:maxpatch),fpatch2d(i,j,1:maxpatch)
         end do
      end do
      close (lusrf)
! convert latitudes and longitudes from degress to radians
      pi = 4.*atan(1.)        
      latixy(:,:) = latixy(:,:)*pi/180. 
      longxy(:,:) = longxy(:,:)*pi/180. 
! ----------------------------------------------------------------------
! [2] MAPPING and ALLOCATE
! Build 1d subgrid patch <-> 2d grid mapping indices and weights
!
! Build mapping indices and weights: [lon_points]x[lat_points] 2d grid <->
! <-> [numpatch] vector of subgrid patches.
! The land surface model works by gathering all the land points on a
! [lon_points]x[lat_points] grid into a vector, and then expanded into
! a vector of [numpatch] subgrid patches, allowing
! for up to [maxpatch] subgrid patches per land point.
! [ixy], [jxy], [patch], and [land] are indices for the mapping:
! [lon_points]x[lat_points] grid <-> [numpatch] vector of subgrid points.
!
!-----------------------------------------------------------------------
! Find total number of patches [numpatch] allowing for multiple subgrid
! patches in a grid cell.
! --------------------------------------------------------------------
      npatch = 0
      numpatch_lat(:) = 0
      do j = 1, lat_points
         do i = 1, lon_points
            do m = 2, maxpatch                
               if(fpatch2d(i,j,m)> 0.)then
                  npatch = npatch+1 !subgrid patch number
                  numpatch_lat(j) = numpatch_lat(j) + 1
               endif
            end do
         end do
      end do
      numpatch = npatch
      if(numpatch.ne.sum(numpatch_lat))then
         write(6,*) 'Total number of patches NOT as the summation of numpatch_lat'
         call abort
      endif
!      write(6,*) 'Total land patches = ', numpatch
! --------------------------------------------------------------------
! Allocates memory for CLM 1d [numpatch] variables
! --------------------------------------------------------------------
      allocate (isc                  (numpatch))      
      allocate (rockdep              (numpatch))  
      allocate (sand         (nl_soil,numpatch))    
      allocate (clay         (nl_soil,numpatch))    
      allocate (work                 (numpatch))
      allocate (ixy_patch            (numpatch))
      allocate (jxy_patch            (numpatch))
      allocate (mxy_patch            (numpatch))
      allocate (wtxy_patch           (numpatch))
      allocate (dlat                 (numpatch))
      allocate (dlon                 (numpatch))
      allocate (itypwat              (numpatch))
      allocate (ivt                  (numpatch))
      allocate (albsol               (numpatch))
      allocate (csol         (nl_soil,numpatch))
      allocate (porsl        (nl_soil,numpatch))
      allocate (phi0         (nl_soil,numpatch))
      allocate (bsw          (nl_soil,numpatch))
      allocate (dkmg         (nl_soil,numpatch))
      allocate (dksatu       (nl_soil,numpatch))
      allocate (dkdry        (nl_soil,numpatch))
      allocate (hksati       (nl_soil,numpatch))
      allocate (z0m                  (numpatch))
      allocate (displa               (numpatch))
      allocate (sqrtdi               (numpatch))
      allocate (effcon               (numpatch))
      allocate (vmax25               (numpatch))
      allocate (slti                 (numpatch))
      allocate (hlti                 (numpatch))
      allocate (shti                 (numpatch))
      allocate (hhti                 (numpatch))
      allocate (trda                 (numpatch))
      allocate (trdm                 (numpatch))
      allocate (trop                 (numpatch))
      allocate (gradm                (numpatch))
      allocate (binter               (numpatch))
      allocate (extkn                (numpatch))
      allocate (chil                 (numpatch))
      allocate (ref              (2,2,numpatch))
      allocate (tran             (2,2,numpatch))
      allocate (rootfr       (nl_soil,numpatch))
      allocate (z   (maxsnl+1:nl_soil,numpatch))
      allocate (dz  (maxsnl+1:nl_soil,numpatch))
      allocate (tss (maxsnl+1:nl_soil,numpatch))
      allocate (wliq(maxsnl+1:nl_soil,numpatch))
      allocate (wice(maxsnl+1:nl_soil,numpatch))
      allocate (tg                   (numpatch))
      allocate (tlsun                (numpatch))
      allocate (tlsha                (numpatch))
      allocate (ldew                 (numpatch))
      allocate (sag                  (numpatch))
      allocate (scv                  (numpatch))
      allocate (snowdp               (numpatch))
      allocate (fveg                 (numpatch))
      allocate (fsno                 (numpatch))
      allocate (sigf                 (numpatch))
      allocate (green                (numpatch))
      allocate (lai                  (numpatch))
      allocate (sai                  (numpatch))
      allocate (coszen               (numpatch))
      allocate (albg             (2,2,numpatch))
      allocate (albv             (2,2,numpatch))
      allocate (alb              (2,2,numpatch))
      allocate (ssun             (2,2,numpatch))
      allocate (ssha             (2,2,numpatch))
      allocate (thermk               (numpatch))
      allocate (extkb                (numpatch))
      allocate (extkd                (numpatch))
      allocate (trad                 (numpatch))
      allocate (tref                 (numpatch))
      allocate (qref                 (numpatch))
      allocate (rst                  (numpatch))
      allocate (emis                 (numpatch))
      allocate (z0ma                 (numpatch))
      allocate (zol                  (numpatch))
      allocate (rib                  (numpatch))
      allocate (ustar                (numpatch))
      allocate (qstar                (numpatch))
      allocate (tstar                (numpatch))
      allocate (fm                   (numpatch))
      allocate (fh                   (numpatch))
      allocate (fq                   (numpatch))
      allocate (fcon(numpatch,nfcon))
      allocate (fvar(numpatch,nfvar))
! --------------------------------------------------------------------
! Build 1d land vector and 1d patch vector mapping components
! --------------------------------------------------------------------
! Determine land vector and patch vector mapping components
      wtxy_patch(:)  = 0.
      npatch = 0
      do j = 1, lat_points
         do i = 1, lon_points
            do m = 2, maxpatch                           
               if(fpatch2d(i,j,m)>0.)then                
               npatch = npatch+1                      
               ixy_patch(npatch)  = i             !patch longitude index
               jxy_patch(npatch)  = j             !patch latitude index
               mxy_patch(npatch)  = m             !patch subgrid index of lnd point
              wtxy_patch(npatch)  = fpatch2d(i,j,m)  !patch weight
                    dlat(npatch)  = latixy(i,j)   !latitude in radians
                    dlon(npatch)  = longxy(i,j)   !longitude in radians
                     ivt(npatch)  = surf2d(i,j,m) !land cover type
                     isc(npatch)  = soic2d(i,j)   !soil color index
                 rockdep(npatch)  = rock2d(i,j)   !depth to bed rock
                   sand(:,npatch) = sand2d(i,j,:) !percent of sand
                   clay(:,npatch) = clay2d(i,j,:) !percent of clay
               end if
            end do
         end do
      end do
      if(numpatch.ne.npatch)then
         write(6,*) 'the number of patches is not identical ', numpatch, npatch
         call abort
      endif
! --------------------------------------------------------------------
! [3]
! INITIALIZE TIME INVARIANT VARIABLES
! ----------------------------------------------------------------------
      do i = 1, numpatch
      CALL iniTimeConst(nl_soil,ivt(i),isc(i),sand(1:,i),clay(1:,i),rockdep(i)&
          ,itypwat(i),z(1:,i),dz(1:,i),albsol(i),csol(1:,i),porsl(1:,i)&
          ,phi0(1:,i),bsw(1:,i),dkmg(1:,i),dksatu(1:,i),dkdry(1:,i),hksati(1:,i)&
          ,z0m(i),displa(i),sqrtdi(i),effcon(i),vmax25(i),slti(i),hlti(i)&
          ,shti(i),hhti(i),trda(i),trdm(i),trop(i), gradm(i),binter(i),extkn(i)&
          ,chil(i),ref(1:,1:,i),tran(1:,1:,i),rootfr(1:,i)&
          ,zlnd,zsno,csoilc,dewmx,wtfact,capr,cnfac&
          ,ssi,wimp,pondmx,smpmax,smpmin,trsmx0,tcrit)  
      enddo
!      print*,'initialize.F90-test-vmax25',vmax25
! ----------------------------------------------------------------------
! [4]
! INITIALIZE TIME-VARYING VARIABLES, as subgrid vectors of length [numpatch]
! initial run: create the time-varying variables based on :
!              i) observation (NOT CODING CURRENTLY), or
!             ii) some already-known information (NO CODING CURRENTLY), or
!            iii) arbitrarily
! continuation run: time-varying data read in from restart file
! ----------------------------------------------------------------------
!4.1 current time of model run
      year = start_yr
      jday = start_jday
      msec = start_sec
      if(.not. greenwich)then
! convert local time to GMT
! off-line cases: input local time was assumed the time at the first grid point
!                 if not this case, please make your change
      a = longxy(1,1)/(15.*pi/180.)*3600.
      msec = msec - int(a)
      if(msec<0)then
         jday = jday -1
         msec = 86400 + msec
      endif
      if(jday<1)then
         year=year-1
         if((mod(year,4)==0 .AND. mod(year,100)/=0) .OR. mod(year,400)==0)then
             jday = 366
          else
             jday = 365
          endif
      endif
      endif
!4.2 cosine of solar zenith angle
      calday = float(jday)+float(msec)/86400.
      do i = 1, numpatch
      coszen(i) = orb_coszen(calday,dlon(i),dlat(i))
      enddo
!4.4 LEAF area index
! CREAT fraction of vegetation cover, greenness, leaf area index, stem index
      lai(:)=0.0; sai(:)=0.0; green(:)=0.0; fveg(:)=0.0
      do i = 1, numpatch
         tss(1:,i) = 283.
! call EcoModel()
         if(ivt(i)>0)then
            call lai_empirical(ivt(i),nl_soil,rootfr(1:,i)&
                              ,tss(1:,i),lai(i),sai(i),fveg(i),green(i))
         endif
      enddo
!4.5 initialize time-varying variables, as subgrid vectors of length [numpatch]
      do i = 1, numpatch
      CALL iniTimeVar(nl_soil,maxsnl,itypwat(i)&
          ,porsl(1:,i),albsol(i),z0m(i),chil(i),ref(1:,1:,i),tran(1:,1:,i)&
	  ,z(maxsnl+1:,i),dz(maxsnl+1:,i)&
          ,tss(maxsnl+1:,i),wliq(maxsnl+1:,i),wice(maxsnl+1:,i)&
          ,tg(i),tlsun(i),tlsha(i),ldew(i),sag(i),scv(i)&
          ,snowdp(i),fveg(i),fsno(i),sigf(i),green(i),lai(i),sai(i),coszen(i)&
          ,albg(1:,1:,i),albv(1:,1:,i),alb(1:,1:,i),ssun(1:,1:,i),ssha(1:,1:,i)&
          ,thermk(i),extkb(i),extkd(i)&
          ,trad(i),tref(i),qref(i),rst(i),emis(i),z0ma(i),zol(i),rib(i)&
          ,ustar(i),qstar(i),tstar(i),fm(i),fh(i),fq(i)&
          )
      enddo
!------------------------------------------------------------------------------
! [5]
! Transfer the time invariant and time-varying variables
!------------------------------------------------------------------------------
! Time invariant model variables
      if(nfcon /= 9*nl_soil+29) call abort
      do i = 1, numpatch
         ub = 1
         fcon(i,ub)    = dlat   (i)           ; ub = ub + 1                    !1
         fcon(i,ub)    = dlon   (i)           ; ub = ub + 1                    !2
         fcon(i,ub)    = itypwat(i)           ; ub = ub + 1                    !3
         fcon(i,ub)    = ivt    (i)           ; ub = ub + 1                    !4
         fcon(i,ub)    = albsol (i)           ; lb = ub + 1; ub = ub + nl_soil !5
         fcon(i,lb:ub) = csol   (1:nl_soil,i) ; lb = ub + 1; ub = ub + nl_soil !1_
         fcon(i,lb:ub) = porsl  (1:nl_soil,i) ; lb = ub + 1; ub = ub + nl_soil !2_
         fcon(i,lb:ub) = phi0   (1:nl_soil,i) ; lb = ub + 1; ub = ub + nl_soil !3_
         fcon(i,lb:ub) = bsw    (1:nl_soil,i) ; lb = ub + 1; ub = ub + nl_soil !4_
         fcon(i,lb:ub) = dkmg   (1:nl_soil,i) ; lb = ub + 1; ub = ub + nl_soil !5_
         fcon(i,lb:ub) = dksatu (1:nl_soil,i) ; lb = ub + 1; ub = ub + nl_soil !6_
         fcon(i,lb:ub) = dkdry  (1:nl_soil,i) ; lb = ub + 1; ub = ub + nl_soil !7_
         fcon(i,lb:ub) = hksati (1:nl_soil,i) ; ub = ub + 1                    !8_
         fcon(i,ub)    = z0m     (i)          ; ub = ub + 1                    !6
         fcon(i,ub)    = displa  (i)          ; ub = ub + 1                    !7
         fcon(i,ub)    = sqrtdi  (i)          ; ub = ub + 1                    !8
         fcon(i,ub)    = effcon  (i)          ; ub = ub + 1                    !9
         fcon(i,ub)    = vmax25  (i)          ; ub = ub + 1                    !10
         fcon(i,ub)    = slti    (i)          ; ub = ub + 1                    !11
         fcon(i,ub)    = hlti    (i)          ; ub = ub + 1                    !12
         fcon(i,ub)    = shti    (i)          ; ub = ub + 1                    !13
         fcon(i,ub)    = hhti    (i)          ; ub = ub + 1                    !14
         fcon(i,ub)    = trda    (i)          ; ub = ub + 1                    !15
         fcon(i,ub)    = trdm    (i)          ; ub = ub + 1                    !16
         fcon(i,ub)    = trop    (i)          ; ub = ub + 1                    !17
         fcon(i,ub)    = gradm   (i)          ; ub = ub + 1                    !18
         fcon(i,ub)    = binter  (i)          ; ub = ub + 1                    !19
         fcon(i,ub)    = extkn   (i)          ; ub = ub + 1                    !20
         fcon(i,ub)    = chil    (i)          ; ub = ub + 1                    !21
         fcon(i,ub)    = ref (1,1,i)          ; ub = ub + 1                    !22
         fcon(i,ub)    = ref (1,2,i)          ; ub = ub + 1                    !23
         fcon(i,ub)    = ref (2,1,i)          ; ub = ub + 1                    !24
         fcon(i,ub)    = ref (2,2,i)          ; ub = ub + 1                    !25
         fcon(i,ub)    = tran(1,1,i)          ; ub = ub + 1                    !26
         fcon(i,ub)    = tran(1,2,i)          ; ub = ub + 1                    !27
         fcon(i,ub)    = tran(2,1,i)          ; ub = ub + 1                    !28
         fcon(i,ub)    = tran(2,2,i)          ; lb = ub + 1; ub = ub + nl_soil !29
         fcon(i,lb:ub) = rootfr(1:nl_soil,i)                                   !9_
      enddo
!      print*,'initialize.F90-test-fcon90-vmax',fcon(1,90)
!      print*,'initialize.F90-test-fcon90-vmax',fcon(2,90)
!      print*,'-------- clm fcon -------'
!      do ub=1,nfcon
!      work(:)= fcon(1:numpatch,ub)
!      print 100, ub, maxval(work), minval(work)
!      enddo
!100   format(1h (i6, 2e16.3), '        nfcon')
! CLM time step and TUNABLE constants
      if(nftune /= 14) call abort
      ftune(1)  = zlnd   
      ftune(2)  = zsno   
      ftune(3)  = csoilc 
      ftune(4)  = dewmx  
      ftune(5)  = wtfact 
      ftune(6)  = capr   
      ftune(7)  = cnfac  
      ftune(8)  = ssi    
      ftune(9)  = wimp   
      ftune(10) = pondmx 
      ftune(11) = smpmax 
      ftune(12) = smpmin 
      ftune(13) = trsmx0 
      ftune(14) = tcrit  
!      print*,'-------- clm ftune -------'
!      print 101, maxval(ftune), minval(ftune)
!101   format(1h (2e16.3), '        nftune ')
! --------------------------------------------
! write out as a restart file [histTimeConst]
! --------------------------------------------
      CALL rstTimeConstWrite (lat_points,lhistTimeConst,nfcon,nftune&
                             ,numpatch,numpatch_lat,ixy_patch,jxy_patch,mxy_patch,wtxy_patch&
                             ,fcon,ftune)
!---------------------------------------------
! write out as a restart file[histTimeConst2]
!---------------------------------------------
!      CALL rstTimeConstWrite2 (lat_points,lhistTimeConst2,nfcon,nftune&
!                              ,numpatch,numpatch_lat,ixy_patch,jxy_patch,mxy_patch,wtxy_patch&
!                            ,fcon,ftune)
!      write (6,*)
!      write (6,*) ('successfully to initialize the land Time invariant ')
! the model variables for restart run
      idate(1) = year     
      idate(2) = jday     
      idate(3) = msec     
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
!      print*,'-------- clm fvar out -------'
!      do ub=1,nfvar
!      work(:)= fvar(1:numpatch,ub)
!      print 102, ub, maxval(work), minval(work)
!      enddo
!102   format(1h (i6, 2e16.3), '        nfvar after')
! average subgrid albedos, srf temperature, etc. for atmospheric model
      fldxy(:,:,:) = 0.0
      do k = 1, numpatch
         i = ixy_patch(k)
         j = jxy_patch(k)
         fldxy(i,j,53) = fldxy(i,j,53) + wtxy_patch(k)*tg     (k)
         fldxy(i,j,64) = fldxy(i,j,64) + wtxy_patch(k)*alb(1,1,k)
         fldxy(i,j,65) = fldxy(i,j,65) + wtxy_patch(k)*alb(1,2,k)
         fldxy(i,j,66) = fldxy(i,j,66) + wtxy_patch(k)*alb(2,1,k)
         fldxy(i,j,67) = fldxy(i,j,67) + wtxy_patch(k)*alb(2,2,k)
         fldxy(i,j,70) = fldxy(i,j,70) + wtxy_patch(k)*tg     (k)
      enddo
      fldxy(:,:,75) = -0.1 
      fldxy(:,:,76) = alog(30.) 
      fldxy(:,:,77) = alog(30.) 
      fldxy(:,:,78) = alog(30.) 
! ------------------------------------------------------------
! write out the model variables for restart run [histTimeVar]
! ------------------------------------------------------------
      CALL rstTimeVarWrite (lhistTimeVar,nfvar,numpatch,idate,fvar)
      CALL rstTimeVarWrite2(lhistTimeVar2,nfvar,numpatch,idate,fvar)
!      write (6,*)
!      write (6,*) ('successfully to initialize the land Time-vraying variables')
! --------------------------------------------------
! Deallocates memory for CLM 1d [numpatch] variables
! --------------------------------------------------
      deallocate (isc    )      
      deallocate (rockdep)  
      deallocate (sand   )    
      deallocate (clay   )    
      deallocate (work   )
      deallocate (ixy_patch )
      deallocate (jxy_patch )
      deallocate (mxy_patch )
      deallocate (wtxy_patch)
      deallocate (dlat   )
      deallocate (dlon   )
      deallocate (itypwat)
      deallocate (ivt    )
      deallocate (albsol )
      deallocate (csol   )
      deallocate (porsl  )
      deallocate (phi0   )
      deallocate (bsw    )
      deallocate (dkmg   )
      deallocate (dksatu )
      deallocate (dkdry  )
      deallocate (hksati )
      deallocate (z0m    )
      deallocate (displa )
      deallocate (sqrtdi )
      deallocate (effcon )
      deallocate (vmax25 )
      deallocate (slti   )
      deallocate (hlti   )
      deallocate (shti   )
      deallocate (hhti   )
      deallocate (trda   )
      deallocate (trdm   )
      deallocate (trop   )
      deallocate (gradm  )
      deallocate (binter )
      deallocate (extkn  )
      deallocate (chil   )
      deallocate (ref    )
      deallocate (tran   )
      deallocate (rootfr )
      deallocate (z      )
      deallocate (dz     )
      deallocate (tss    )
      deallocate (wliq   )
      deallocate (wice   )
      deallocate (tg     )
      deallocate (tlsun  )
      deallocate (tlsha  )
      deallocate (ldew   )
      deallocate (sag    )
      deallocate (scv    )
      deallocate (snowdp )
      deallocate (fveg   )
      deallocate (fsno   )
      deallocate (sigf   )
      deallocate (green  )
      deallocate (lai    )
      deallocate (sai    )
      deallocate (coszen )
      deallocate (albg   )
      deallocate (albv   )
      deallocate (alb    )
      deallocate (ssun   )
      deallocate (ssha   )
      deallocate (thermk )
      deallocate (extkb  )
      deallocate (extkd  )
      deallocate (trad   )
      deallocate (tref   )
      deallocate (qref   )
      deallocate (rst    )
      deallocate (emis   )
      deallocate (z0ma   )
      deallocate (zol    )
      deallocate (rib    )
      deallocate (ustar  )
      deallocate (qstar  )
      deallocate (tstar  )
      deallocate (fm     )
      deallocate (fh     )
      deallocate (fq     )
      deallocate (fcon   )
      deallocate (fvar   )
 end subroutine initialize
