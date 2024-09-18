  subroutine vec2xy(lat_points,lon_points,numpatch,&
		    ixy_patch,jxy_patch,wtxy_patch,itypwat,nfcon,nforc,nfldv,&
                    fcon,forcxy,fldv,fldxy_r)
! ----------------------------------------------------------------------
! perfrom grid-average from subgrid 1d vector
! subgrid to grid average mapping: average a subgrid input vector [fldv]
! of length to a 2-d [lon_points] x [lat_points] output array [fldxy_r]
!
! Created by Yongjiu Dai
!--------------!-------------------------------------------------------
! 01: taux     ! wind stress: E-W [kg/m/s2]
! 02: tauy     ! wind stress: N-S [kg/m/s2]
! 03: fsena    ! sensible heat from canopy height to atmosphere [W/m2]
! 04: lfevpa   ! latent heat flux from canopy height to atmosphere [W/m2]
! 05: fevpa    ! evapotranspiration from canopy to atmosphere [mm/s]
! 06: fsenl    ! sensible heat from leaves [W/m2]
! 07: fevpl    ! evaporation+transpiration from leaves [mm/s]
! 08: etr      ! transpiration rate [mm/s]
! 09: fseng    ! sensible heat flux from ground [W/m2]
! 10: fevpg    ! evaporation heat flux from ground [mm/s]
! 11: fgrnd    ! ground heat flux [W/m2]
! 12: sabvsun  ! solar absorbed by sunlit canopy [W/m2]
! 13: sabvsha  ! solar absorbed by shaded [W/m2]
! 14: sabg     ! solar absorbed by ground  [W/m2]
! 15: olrg     ! outgoing long-wave radiation from ground+canopy [W/m2]
! 16: rnet     ! net radiation [W/m2]
! 17: xerr     ! the error of water banace [mm/s]
! 18: zerr     ! the error of energy balance [W/m2]
! 19: rsur     ! surface runoff [mm/s]
! 20: rnof     ! total runoff [mm/s]
! 21: assim    ! canopy assimilation rate [mol m-2 s-1]
! 22: respc    ! respiration (plant+soil) [mol m-2 s-1]
!--------------!-------------------------------------------------------
! 23:32: tss   ! soil temperature [K]
! 33:42: wliq  ! liquid water in soil layers [kg/m2]
! 43:52: wice  ! ice lens in soil layers [kg/m2]
!--------------!-------------------------------------------------------
! 53: tg       ! ground surface temperature [K]
! 54: tlsun    ! sunlit leaf temperature [K]
! 55: tlsha    ! shaded leaf temperature [K]
! 56: ldew     ! depth of water on foliage [mm]
! 57: scv      ! snow cover, water equivalent [mm]
! 58: snowdp   ! snow depth [meter]
! 59: fsno     ! fraction of snow cover on ground
! 60: sigf     ! fraction of veg cover, excluding snow-covered veg [-]
! 61: green    ! leaf greenness
! 62: lai      ! leaf area index
! 63: sai      ! stem area index
! 64: alb(1,1) ! averaged albedo [visible, direct]
! 65: alb(1,2) ! averaged albedo [visible, diffuse]
! 66: alb(2,1) ! averaged albedo [near-infrared, direct]
! 67: alb(2,2) ! averaged albedo [near-infrared,diffuse]
! 68: emis     ! averaged bulk surface emissivity
! 69: z0ma     ! effective roughness [m]
!--------------!-------------------------------------------------------
! 70: trad     ! radiative temperature of surface [K]
! 71: ustar    ! u* in similarity theory [m/s]
! 72: tstar    ! t* in similarity theory [kg/kg]
! 73: qstar    ! q* in similarity theory [kg/kg]
! 74: zol      ! dimensionless height (z/L) used in Monin-Obukhov theory
! 75: rib      ! bulk Richardson number in surface layer
! 76: fm       ! integral of profile function for momentum
! 77: fh       ! integral of profile function for heat
! 78: fq       ! integral of profile function for moisture
!--------------!-------------------------------------------------------
! 79: tref     ! 2 m height air temperature [kelvin]
! 80: qref     ! 2 m height air specific humidity [kg/kg]
! 81: u10m     ! 10m u-velocity [m/s]
! 82: v10m     ! 10m v-velocity [m/s]
! 83: f10m     ! integral of profile function for momentum at 10m [-]
!--------------!-------------------------------------------------------
! 84: us       ! wind in eastward direction [m/s]
! 85: vs       ! wind in northward direction [m/s]
! 86: tm       ! temperature at reference height [kelvin]
! 87: qm       ! specific humidity at reference height [kg/kg]
! 88: prc      ! convective precipitation [mm/s]
! 89: prl      ! large scale precipitation [mm/s]
! 90: pbot     ! atmospheric pressure at the surface [pa]
! 91: frl      ! atmospheric infrared (longwave) radiation [W/m2]
! 92: solar    ! downward solar radiation at surface [W/m2]
!--------------!-------------------------------------------------------
      use precision
      use phycon_module, only: vonkar, stefnc, cpair, rgas, grav
      implicit none
! arguments:
      integer, intent(in) :: lon_points            ! number of longitude points on model grid
      integer, intent(in) :: lat_points            ! number of latitude points on model grid
		
      integer, intent(in) :: numpatch              ! total number of patches of grids
      integer, intent(in) :: ixy_patch(numpatch)   ! patch longitude index
      integer, intent(in) :: jxy_patch(numpatch)   ! patch latitude index
      integer, intent(in) :: itypwat(numpatch)     ! land water type
      real(r8), intent(in) :: wtxy_patch(numpatch) ! patch weight
      integer, INTENT(in) :: nfcon                 ! number of time constant variables
      integer, INTENT(in) :: nforc                 ! number of forcing variables
      integer, intent(in) :: nfldv                 ! number of output variables
      real(r8), INTENT(in) :: fcon(numpatch,nfcon) ! time constant variables
      real(r8), intent(in) :: fldv(numpatch,nfldv) ! output variables
      real(r8), intent(in) :: forcxy(lon_points,lat_points,nforc) ! xy gridded forcing
      real(r8), intent(out) :: fldxy_r(lon_points,lat_points,nfldv) ! xy gridded output
! local variables
      integer  :: i,j,k,l                            ! indices
      real(r8) :: a                                  !
      real(r8) :: sumwt(lon_points,lat_points,nfldv) ! sum of wt
      real(r8) :: maxwt(lon_points,lat_points)       ! maximum of patches
      integer :: maxfp(lon_points,lat_points)        ! index of maximum patches
      real(r8) rhoair,thm,th,thv,ur,displa,zldis,hu,ht,hq
      real(r8) z0m,z0h,z0q,us,vs,tm,qm,pbot,psrf
      real(r8) obu,temp1,temp2,temp12m,temp22m
      real(r8) um,thvstar,beta,zii,wc,wc2
      integer ivt
!----------------------------------------------------------------------
      fldxy_r(:,:,:) = 0.0
      sumwt(:,:,:) = 0.0
      maxwt(:,:) = 0.0
      do k = 1, numpatch
         i = ixy_patch(k)
         j = jxy_patch(k)
         maxfp(i,j) = k 
      enddo
! Get the number of patch with the largest fraction of grid within 1d clm array
      do k = 1, numpatch
         i = ixy_patch(k)
         j = jxy_patch(k)
         a = maxwt(i,j)
         maxwt(i,j) = max(a,wtxy_patch(k))
         if((maxwt(i,j)-a).gt.0.0) maxfp(i,j) = k 
      enddo
! Mapping the 1d [numpatch] to grid [lon_points x lat_points]
      do L = 1, 69
         do k = 1, numpatch
            i = ixy_patch(k)
            j = jxy_patch(k)
            if(L.le.22)then             ! fluxes
![1-22] Grid averages by area-weight over grid patches
               fldxy_r(i,j,L) = fldxy_r(i,j,L) + wtxy_patch(k)*fldv(k,L)
               sumwt(i,j,L) = sumwt(i,j,L) + wtxy_patch(k)
            else if(L.le.52)then        ! soil temperature and water
![23-32 33-42 43-52] Area-weight over grid patches but excluding lake and ocean patches
               if(itypwat(k).le.3)then  ! lake and ocean excluded
               fldxy_r(i,j,L) = fldxy_r(i,j,L) + wtxy_patch(k)*fldv(k,L)
               sumwt(i,j,L) = sumwt(i,j,L) + wtxy_patch(k)
               endif
            else                        ! clm state variables
![53-69] Grid averages by area-weight over grid patches
               fldxy_r(i,j,L) = fldxy_r(i,j,L) + wtxy_patch(k)*fldv(k,L)
               sumwt(i,j,L) = sumwt(i,j,L) + wtxy_patch(k)
            endif
         enddo
      enddo
      do L = 1, 69
         do j = 1, lat_points
            do i = 1, lon_points
!              if(sumwt(i,j,L).gt.1.0 .or. sumwt(i,j,L).lt.0.0)then
!                 write(6,*) 'summation of fraction patches = ', sumwt(i,j,l), i,j,j
!                 call abort
!              endif
               if(sumwt(i,j,L).gt.0.)then
                  fldxy_r(i,j,L) = fldxy_r(i,j,L)/sumwt(i,j,L)
               else
                  fldxy_r(i,j,L) = -9999.
               endif
            enddo
         enddo
      enddo
![70-78] Retrieve through averaged fluxes
!     do L = 70, 78
!        do k = 1, numpatch
!           i = ixy_patch(k)
!           j = jxy_patch(k)
!           if(k.eq.maxfp(i,j))then  ! take values as that at the largest patches
!                 fldxy_r(i,j,L) = fldv(k,L)
!                 sumwt(i,j,L) = 1.0
!           endif
!        enddo
!     enddo
      do L = 70, nfldv
         do k = 1, numpatch
            i = ixy_patch(k)
            j = jxy_patch(k)
            sumwt(i,j,L) = sumwt(i,j,L) + wtxy_patch(k)
         enddo
      enddo
      do j = 1, lat_points
         do i = 1, lon_points
            if(sumwt(i,j,70).gt.0.)then         !For land only defined
               z0m = fldxy_r(i,j,69)
               z0h = fldxy_r(i,j,69)
               z0q = fldxy_r(i,j,69)
               displa = 2./3.*z0m/0.07
               hu = max(forcxy(i,j,16),5.+displa)
               ht = max(forcxy(i,j,17),5.+displa)
               hq = max(forcxy(i,j,18),5.+displa)
               zldis = hu-displa
               us = forcxy(i,j,3)
               vs = forcxy(i,j,4)
               tm = forcxy(i,j,5)
               qm = forcxy(i,j,6)
               pbot = forcxy(i,j,9)
               psrf = forcxy(i,j,10)
               rhoair = (pbot-0.378*qm*pbot/(0.622+0.378*qm))/(rgas*tm)
               fldxy_r(i,j,70) = (fldxy_r(i,j,15)/stefnc)**0.25 
               fldxy_r(i,j,71) = sqrt(max(1.e-6,sqrt(fldxy_r(i,j,1)**2+fldxy_r(i,j,2)**2))/rhoair) 
               fldxy_r(i,j,72) = -fldxy_r(i,j,3)/(rhoair*fldxy_r(i,j,71))/cpair
               fldxy_r(i,j,73) = -fldxy_r(i,j,5)/(rhoair*fldxy_r(i,j,71))
               thm = tm + 0.0098*ht
               th = tm*(100000./psrf)**(rgas/cpair)
               thv = th*(1.+0.61*qm)       
               fldxy_r(i,j,74) = zldis*vonkar*grav&
                   * (fldxy_r(i,j,72)+0.61*th*fldxy_r(i,j,73))&
                   / (fldxy_r(i,j,71)**2*thv)
               if(fldxy_r(i,j,74) .ge. 0.)then   !stable
                  fldxy_r(i,j,74) = min(2.,max(fldxy_r(i,j,74),1.e-6))
               else                              !unstable
                  fldxy_r(i,j,74) = max(-100.,min(fldxy_r(i,j,74),-1.e-6))
               endif
               beta = 1.
               zii = 1000.
               thvstar=fldxy_r(i,j,72)+0.61*th*fldxy_r(i,j,73)
               ur = sqrt(us*us+vs*vs)
               if(fldxy_r(i,j,74) .ge. 0.)then
                  um = max(ur,0.1)
               else
                  wc = (-grav*fldxy_r(i,j,71)*thvstar*zii/thv)**(1./3.)
                 wc2 = beta*beta*(wc*wc)
                  um = max(0.1,sqrt(ur*ur+wc2))
               endif
               obu = zldis/fldxy_r(i,j,74)
               call moninobuk(hu,ht,hq,displa,z0m,z0h,z0q,&
                    obu,um,fldxy_r(i,j,71),temp1,temp2,temp12m,temp22m,&
                    fldxy_r(i,j,83),fldxy_r(i,j,76),fldxy_r(i,j,77),fldxy_r(i,j,78))
               fldxy_r(i,j,75) = fldxy_r(i,j,74)*vonkar**3*fldxy_r(i,j,71)**2/(temp1*um**2)
               fldxy_r(i,j,75) = min(5.,fldxy_r(i,j,75)) 
            else
               fldxy_r(i,j,70:78) = -9999.
            endif
         enddo
      enddo
![79-83] Grass land only (for matching the routine meteorological obs.)
      do L = 79, 83
         do k = 1, numpatch
            i = ixy_patch(k)
            j = jxy_patch(k)
            ivt = nint(fcon(k,4)) 
            if(ivt.eq.7 .or. ivt.eq.0)then
               fldxy_r(i,j,L) = fldv(k,L)
            endif
         enddo
      enddo
![84-92] Meteorological forcing
      do j = 1, lat_points
         do i = 1, lon_points
            if(sumwt(i,j,84).gt.0.)then      !For land only defined
               fldxy_r(i,j,84) = forcxy(i,j,3)
               fldxy_r(i,j,85) = forcxy(i,j,4)
               fldxy_r(i,j,86) = forcxy(i,j,5)
               fldxy_r(i,j,87) = forcxy(i,j,6)
               fldxy_r(i,j,88) = forcxy(i,j,7)
               fldxy_r(i,j,89) = forcxy(i,j,8)
               fldxy_r(i,j,90) = forcxy(i,j,9)
               fldxy_r(i,j,91) = forcxy(i,j,15)
               fldxy_r(i,j,92) = forcxy(i,j,11)+forcxy(i,j,12)+forcxy(i,j,13)+forcxy(i,j,14)
            else
               fldxy_r(i,j,84:92) = -9999.
            endif
         enddo
      enddo
  end subroutine vec2xy
! ----------------------------------------------------------------------
! EOP
