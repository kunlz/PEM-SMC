 subroutine newsnow (itypwat,maxsnl,dtime,tm,tg,pg,tcrit, & 
		     zi,z,dz,tss,wliq,wice,fiold,snl,sag,scv,snowdp,pg_rain,pg_snow)
!=======================================================================
! add new snow nodes.
! Original author : Yongjiu Dai, 09/15/1999; 08/31/2002
!=======================================================================
!
  use precision
  use phycon_module, only : tfrz
  implicit none
! ------------------------ Dummy Argument ------------------------------
  integer, INTENT(in) :: maxsnl  ! maximum number of snow layers
  integer, INTENT(in) :: itypwat ! land water type (0=soil, 1=urban and built-up,
! 2=wetland, 3=land ice, 4=deep lake, 5=shallow lake)
  real(r8), INTENT(in) :: dtime  ! model time step [second]
  real(r8), INTENT(in) :: tm     ! temperature at agcm reference height [kelvin]
  real(r8), INTENT(in) :: tg     ! ground surface temperature [k]
  real(r8), INTENT(in) :: pg     ! water onto ground including canopy runoff [kg/(m2 s)]
  real(r8), INTENT(in) :: tcrit  ! critical temp. to determine rain or snow
  real(r8), INTENT(inout) ::    zi(maxsnl:0)   ! interface level below a "z" level (m)
  real(r8), INTENT(inout) ::     z(maxsnl+1:0) ! layer depth (m)
  real(r8), INTENT(inout) ::    dz(maxsnl+1:0) ! layer thickness (m)
  real(r8), INTENT(inout) ::   tss(maxsnl+1:0) ! soil + snow layer temperature [K]
  real(r8), INTENT(inout) ::  wliq(maxsnl+1:0) ! liquid water (kg/m2)
  real(r8), INTENT(inout) ::  wice(maxsnl+1:0) ! ice lens (kg/m2)
  real(r8), INTENT(inout) :: fiold(maxsnl+1:0) ! fraction of ice relative to the total water
   integer, INTENT(inout) :: snl               ! number of snow layers
  real(r8), INTENT(inout) :: sag               ! non dimensional snow age [-]
  real(r8), INTENT(inout) :: scv               ! snow mass (kg/m2)
  real(r8), INTENT(inout) :: snowdp            ! snow depth (m)
  real(r8), INTENT(out) :: pg_rain  ! liquid water onto ground [kg/(m2 s)]
  real(r8), INTENT(out) :: pg_snow  ! ice onto ground [kg/(m2 s)]
! ----------------------- Local  Variables -----------------------------
  real(r8) bifall     ! bulk density of newly fallen dry snow [kg/m3]
  real(r8) flfall     ! fraction of liquid water within falling precip.
  real(r8) dz_snowf   ! layer thickness rate change due to precipitation [mm/s]
  integer newnode     ! signification when new snow node is set, (1=yes, 0=non)
!-----------------------------------------------------------------------
! the upper limit of air temperature is set for snowfall, this cut-off
! was selected based on Fig. 1, Plate 3-1, of Snow Hydrology (1956).
! the percentage of liquid water by mass, which is arbitrarily set to
! vary linearly with air temp, from 0% at 273.16 to 40% max at 275.16.
      newnode = 0
      if(tm>tfrz+tcrit)then
        flfall = 1.        ! fraction of liquid water within falling precip.
        pg_snow = 0.       ! ice onto ground (mm/s)
        pg_rain = pg       ! liquid water onto ground (mm/s)
        dz_snowf = 0.      ! rate of snowfall, snow depth/s (m/s)
      else
        if(tm<=tfrz)then
          flfall=0.
        else if(tm<=tfrz+2.)then
          flfall=-54.632+0.2*tm
        else
          flfall=0.4
        endif
! use Alta relationship, Anderson(1976); LaChapelle(1961),
! U.S.Department of Agriculture Forest Service, Project F,
! Progress Rep. 1, Alta Avalanche Study Center:Snow Layer Densification.
        if(tm>tfrz+2.)then
          bifall=169.
        else if(tm>tfrz-15.)then
          bifall=50.+1.7*(tm-tfrz+15.)**1.5
        else
          bifall=50.
        endif
        pg_snow = pg*(1.-flfall)                 
        pg_rain = pg*flfall
        dz_snowf = pg_snow/bifall                
        snowdp = snowdp + dz_snowf*dtime         
        scv = scv + pg_snow*dtime      ! snow water equivalent (mm)
        if(itypwat==2 .AND. tg>tfrz)then  ! snowfall on warmer wetland
           scv=0.; snowdp=0.; sag=0.
        endif
      endif
! when the snow accumulation exceeds 10 mm, initialize a snow layer
      if(snl==0 .AND. pg_snow>0.0 .AND. snowdp>=0.01)then  
         snl = -1
         newnode = 1
         dz(0)  = snowdp             ! meter
         z (0)  = -0.5*dz(0)
         zi(-1) = -dz(0)
! currently, the water temperature for the precipitation is simply set
! as the surface air temperature
         sag = 0.                    ! snow age
         tss (0) = min(tfrz, tm)     ! K
         wice(0) = scv               ! kg/m2
         wliq(0) = 0.                ! kg/m2
         fiold(0) = 1.
!**      write(6,*) 'snow layer is built'
      endif
! the change of ice partial density of surface node due to precipitation
! only ice part of snowfall is added here, the liquid part will be added latter
      if(snl<0 .AND. newnode==0)then
         wice(snl+1) = wice(snl+1)+dtime*pg_snow
         dz(snl+1) = dz(snl+1)+dz_snowf*dtime
         z(snl+1) = zi(snl+1) - 0.5*dz(snl+1)
         zi(snl) = zi(snl+1) - dz(snl+1)
      endif
 end subroutine newsnow
