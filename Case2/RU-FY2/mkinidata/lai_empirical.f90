 subroutine lai_empirical(ivt,nl_soil,rootfr,t,lai,sai,fveg,green) 
!-----------------------------------------------------------------------
! provides leaf and stem area parameters
! Original author : Yongjiu Dai, 08/31/2002
!-----------------------------------------------------------------------
 use precision
 implicit none
 integer, intent(in)  :: ivt      !land cover type
 integer, intent(in)  :: nl_soil  !number of soil layers
 real(r8), intent(in)  :: rootfr(1:nl_soil)  !root fraction
 real(r8), intent(in)  :: t(1:nl_soil)  !soil temperature
 real(r8), intent(out) :: lai     !leaf area index
 real(r8), intent(out) :: sai     !Stem area index
 real(r8), intent(out) :: fveg    !fractional cover of vegetation
 real(r8), intent(out) :: green   !greenness
!local variable
 real(r8) f      !
 real(r8) roota  !accumulates root fraction
 integer jrt     !number of soil layers with 90% root fraction
 integer j       !number of soil layers
!-----------------------------------------------------------------------
! Maximum fractional cover of vegetation [-]
  real(r8), dimension(24), parameter :: &
  vegc=(/1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, &
	 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, &
	 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0 /)
! Maximum leaf area index, the numbers are based on the data of
! "worldwide histrorical estimates of leaf area index, 1932-2000" :
! http://www.daac.ornl.gov/global_vegetation/HistoricalLai/data"
  real(r8), dimension(24), parameter :: &
   xla=(/1.50, 3.29, 4.18, 3.50, 2.50, 3.60, 2.02, 1.53, &
         2.00, 0.85, 4.43, 4.42, 4.56, 4.50, 4.50, 0.00, &
         4.00, 3.63, 0.00, 0.64, 1.60, 1.00, 0.00, 0.00 /) 
! Minimum leaf area index
  real(r8), dimension(24), parameter :: &
! xla0=(/1.00, 0.50, 0.50, 0.50, 1.00, 0.50, 0.50, 0.50, &
!        0.50, 0.50, 0.50, 0.50, 4.00, 4.00, 4.00, 0.00, &
!        3.00, 3.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /)
  xla0=(/1.00, 0.50, 0.50, 0.50, 1.00, 0.50, 0.50, 0.50, &
         0.50, 0.30, 0.50, 0.50, 4.00, 0.50, 4.00, 0.00, &
         3.00, 3.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00 /)
! Stem area index [-]
  real(r8), dimension(24), parameter :: & 
! sai0=(/0.50, 0.50, 0.50, 0.50, 1.00, 1.00, 2.00, 1.00, &
!        1.00, 1.00, 2.00, 2.00, 2.00, 2.00, 2.00, 0.00, &
!        2.00, 2.00, 0.00, 0.50, 0.50, 0.50, 0.00, 0.00 /)
  sai0=(/0.20, 0.20, 0.30, 0.30, 0.50, 0.50, 1.00, 0.50, &
         1.00, 0.50, 2.00, 2.00, 2.00, 2.00, 2.00, 0.00, &
         2.00, 2.00, 0.00, 0.10, 0.10, 0.10, 0.00, 0.00 /)
!-----------------------------------------------------------------------
  roota = 0.
  jrt = 1
  do j = 1, nl_soil
     roota = roota + rootfr(j)
     if(roota>0.9)then
        jrt = j
        exit
     endif
  enddo
! Adjust leaf area index for seasonal variation
  f = max(0.0,1.-0.0016*max(298.-t(jrt),0.0)**2)
  lai = xla(ivt) + (xla0(ivt)-xla(ivt))*(1.-f)
! Sum leaf area index and stem area index
  sai = sai0(ivt)
! Fractional vegetation cover
  fveg = vegc(ivt)
  green = 0.0
  if(fveg > 0.) green = 1.0
 end subroutine lai_empirical
