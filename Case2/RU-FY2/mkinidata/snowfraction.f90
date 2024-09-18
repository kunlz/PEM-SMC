 subroutine snowfraction (fveg,z0m,snowdp,wt,sigf,fsno)
!=======================================================================
! Original author : Yongjiu Dai, September 15, 1999
! Provide snow cover fraction
!=======================================================================
  use precision
  implicit none
! dummy arguments
  real(r8), INTENT(in) :: snowdp ! snow depth [m]
  real(r8), INTENT(in) :: z0m    ! aerodynamic roughness length [m]
  real(r8), INTENT(in) :: fveg   ! fractional vegetation cover [-]
  real(r8), INTENT(out) :: wt    ! fraction of vegetation covered with snow [-]
  real(r8), INTENT(out) :: sigf  ! fraction of veg cover, excluding snow-covered veg [-]
  real(r8), INTENT(out) :: fsno  ! fraction of soil covered by snow [-]
!-----------------------------------------------------------------------
      if(fveg > 0.001) then
! Fraction of vegetation buried (covered) by snow
         wt = 0.1*snowdp/z0m
         wt = wt/(1.+wt)
! Fraction of vegetation cover free of snow
         sigf = (1.-wt)*fveg
! Fraction of soil covered by snow
         fsno = snowdp/(0.1+snowdp)
      else
         wt = 0.
         sigf = 0.
         fsno = snowdp/(0.1+snowdp)
      endif
      if(sigf < 0.001) sigf = 0.
      if(sigf > 0.999) sigf = 1.
 end subroutine snowfraction
