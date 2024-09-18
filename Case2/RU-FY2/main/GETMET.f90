 subroutine GETMET (lumet,lon_points,lat_points,nforc,forcxy,jday,msec)
! ======================================================================
!
! Access meteorological data
!
! Original author : Yongjiu Dai
! ======================================================================
   use precision
   implicit none
   integer, INTENT(in) :: lumet      ! logical unit number of atmospheric forcing data
   integer, INTENT(in) :: lon_points ! number of longitude points on model grid
   integer, INTENT(in) :: lat_points ! number of latitude points on model grid
   integer, INTENT(in) :: nforc      ! number of atmospheric forcing variables
   integer, INTENT(in) :: jday       ! julian day
   integer, INTENT(in) :: msec       ! second in a day
   real(r8), INTENT(out) :: forcxy(lon_points,lat_points,nforc)
! local
   real(r8) solar                    ! incident solar radiation [W/m2]
   real(r8) frl                      ! atmospheric infrared (longwave) radiation [W/m2]
   real(r8) prcp                     ! precipitation [mm/s]
   real(r8) tm                       ! temperature at reference height [kelvin]
   real(r8) us                       ! wind component in eastward direction [m/s]
   real(r8) vs                       ! wind component in northward direction [m/s]
   real(r8) pres                     ! atmospheric pressure [pa]
   real(r8) qm                       ! specific humidity at reference height [kg/kg]
   real(r8) dlon                     ! Centered longitude (radians)
   real(r8) dlat                     ! Centered latitude (radians)
   real(r8) sunang                   ! cosine of solar zenith angle
   real(r8) cloud                    ! cloud fraction
   real(r8) difrat                   ! diffuse fraction
   real(r8) vnrat                    ! band fraction
   real(r8) calday                   ! Julian cal day (1.xx to 365.xx)
   real(r8) orb_coszen               ! cosine of solar zenith angle
   integer i, j                      ! looping index
   integer IYEAR,IMONTH,IDAY,IHOUR
   real(r8) rnet,twet,esat,esatdT,qsat,qsatdT,us_d 
   real(r8) rrref, ddd, ggg
   character(len=2) site
   real(r8) time, VPA
   integer flag, Rflag
! ------------------ note ------------------
! the model required the same longitudinal resolution for
! all latitude strip. For the variable longitudinal resolution
! cases, please assign the meteorological values to
! -999 to the land grids which are not included in the calcultion.
! ----------------------------------------------------------------------
![1] PILPS's Valdai (obs height: wind 10 m, tem & hum 2 m)
! & CABAUW (obs height = 20 m for all) & HAPEX (Obs heigh = 2 m for all)
      read (lumet,*) solar, frl, prcp, tm, us, vs, pres, qm
!10    format (2f7.1, e14.3, 3f10.3, f10.1, e12.3)
!      print*,'solar',solar
!      print*,'frl',frl
!      print*,'prcp',prcp
!      print*,'tm',tm
!      print*,'us',us
!      print*,'vs',vs
!      print*,'pres',pres
!      print*,'qm',qm
     dlon=56.447
     dlat=32.9019
!--------------------------------
!  ABRACOS site Reserva Jaru & Fazenda N.S.
! (obs height forest: zwind=52.5)
! (obs height pasture: zwind=5.4, ztem=5.1)
!     read(lumet,30) site,iyear,iday,ihour,&
!                      solar,rrref,rnet,twet,tm,us,ddd,ggg,prcp
30    format(1x,a2,1x,i4,1x,i3,1x,i4,1x,9(f9.2))
!
! for missing data
!     if (prcp.lt.0.0) prcp=0.0
!     if (solar.lt.0.0) solar=0.0
!     if (us.lt.0.0) us = 1.0
!     vs=0.0
!
!     if(rrref.eq.-99.99) rrref=solar*0.122
!     if(rrref.lt.0.0) rrref=0.0
!
!     if(rnet.eq.-99.99) then
!       if(site.eq.'RD' .or. site.eq.'rd') rnet = 0.79*solar - 26.08
!       if(site.eq.'FD' .or. site.eq.'fd') rnet = 0.79*solar - 16.26
!     endif
!
!     pres =1000.0*100.0
! specific humidity
!     tm = tm + 273.16 - 50.
!     if(tm.lt.274.) stop 'meterological data error 2'
!     twet  = twet + 273.16
!     call qsadv(twet,pres,esat,esatdT,qsat,qsatdT)
!     qm = qsat - 3.5*2.8704E2*(tm-twet)/2.50036e6
!
!     !for Arme, Abracos, Tucson, FIFE 88 & 89:
!     !if(istep.gt.11760)then  ! only used for FIFE Forcing data
!     frl = rnet-(solar-rrref)+5.67e-8*tm**4
!     !endif
!
!     prcp=prcp/3600
! ----------------------------------------------------------------------
      do j = 1, lat_points
         do i = 1, lon_points
           forcxy(i,j,1) = pres*355.e-06
           forcxy(i,j,2) = pres*0.209
           forcxy(i,j,3) = us
           forcxy(i,j,4) = vs
           forcxy(i,j,5) = tm
           forcxy(i,j,6) = qm
           forcxy(i,j,7) = 0.
           forcxy(i,j,8) = prcp
           forcxy(i,j,9) = pres
           forcxy(i,j,10) = pres
! -----
! as the downward solar is in full band, an empirical expression
! will be used to divide fractions of band and incident
! Julian calday (1.xx to 365.xx)
! lon: logitude, lat: latitude
           calday = float(jday) + float(msec)/86400.
           sunang = orb_coszen(calday,dlon,dlat)
           cloud = (1160.*sunang-solar)/(963.*sunang)
           cloud = max(cloud,0.)
           cloud = min(cloud,1.)
           cloud = max(0.58,cloud)
           difrat = 0.0604/(sunang-0.0223)+0.0683
           if(difrat.lt.0.) difrat = 0.
           if(difrat.gt.1.) difrat = 1.
           difrat = difrat+(1.0-difrat)*cloud
           vnrat = (580.-cloud*464.)/((580.-cloud*499.)+(580.-cloud*464.))
           forcxy(i,j,11) = (1.0-difrat)*vnrat*solar
!          print*,'11:',forcxy(i,j,11)
           forcxy(i,j,12) = (1.0-difrat)*(1.0-vnrat)*solar
!           print*,'12',forcxy(i,j,12)
           forcxy(i,j,13) = difrat*vnrat*solar
!           print*,'13',forcxy(i,j,13)
           forcxy(i,j,14) = difrat*(1.0-vnrat)*solar
!           print*,'14',forcxy(i,j,14)
! -----
           forcxy(i,j,15) = frl
!           print*,'15',forcxy(i,j,15)
           forcxy(i,j,16) = 25.
           forcxy(i,j,17) = 25.
           forcxy(i,j,18) = 25.
         enddo 
      enddo 
 end subroutine GETMET
