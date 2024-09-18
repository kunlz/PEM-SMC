 subroutine ticktime (dtime, idate)
!=======================================================================
! Original author: Yongjiu Dai, September 15, 1999
!
!      Notes: Greenwich time
!
!=======================================================================
  implicit none
  real, INTENT(in) :: dtime          ! model time step (second)
  integer, INTENT(inout) :: idate(3) ! year and julian day, sec, respectively
  integer maxday                     ! days of one year (365 or 366)
!-----------------------------------------------------------------------
  idate(3) = idate(3) + nint(dtime)
  if(idate(3)>86400)then
     idate(3) = idate(3) - 86400
     if((mod(idate(1),4)==0 .AND. mod(idate(1),100)/=0) .OR. &
                                  mod(idate(1),400)==0)then
         maxday = 366
     else
         maxday = 365
     endif
     idate(2) = idate(2) + 1
     if(idate(2)>maxday) then
        idate(1) = idate(1) + 1
        idate(2) = 1
     endif
  endif
 end subroutine ticktime
