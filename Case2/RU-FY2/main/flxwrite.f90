 subroutine flxwrite (luout,lon_points,lat_points,nfldv,fldxy)
!=======================================================================
! Original version: Yongjiu Dai, September 15, 1999
!=======================================================================
 use precision
 implicit none
 integer, INTENT(in) :: luout      ! logical unit number of output file
 integer, INTENT(in) :: lon_points ! number of longitude points on model grid
 integer, INTENT(in) :: lat_points ! number of latitude points on model grid
 integer, INTENT(in) ::  nfldv     ! number of output variables
 real(r8), INTENT(in) :: fldxy(lon_points,lat_points,nfldv)
 real(r4) :: a(lon_points,lat_points) ! convert to single precision for writing out
 integer i, j, l
! ----------------------------------------------------------------------
     do l = 1, nfldv
        do j = 1, lat_points
	   do i = 1, lon_points
              a(i,j) = fldxy(i,j,l)
           enddo
        enddo
        write(luout) a
     enddo
     close (luout)
  end subroutine flxwrite
! ----------------------------------------------------------------------
! EOP
subroutine flxwrite2(luout2,lon_points,lat_points,nfldv,fldxy)
!-------------------------------------------------------------
! my subroutine to output fluxs
    use precision
    implicit none
    integer,INTENT(in) :: luout2
    integer,INTENT(in) :: lon_points
    integer,INTENT(in) :: lat_points
    integer,INTENT(in) ::nfldv
    real(r8),INTENT(in) ::fldxy(lon_points,lat_points,nfldv)
    real(r4) :: a(lon_points,lat_points)
    real(r4) :: b(1,nfldv)
    integer i,j,l
!--------------------------------
    print*,'flxwrite.f90-test(NAN):write_before:H',fldxy(1,1,3)
    do l=1,nfldv
       do j=1,lat_points
          do i=1,lon_points
             if(l==21 .or. l==22) then
                a(i,j)=fldxy(i,j,l)*1000000
             else
                a(i,j)=fldxy(i,j,l)
            endif
          enddo
       enddo
       b(1,l)=a(1,1)
    enddo
    print*,'flxwrite.f90-test(NAN):write_after:H:',b(1,3)
    write(luout2,100) b
100 format(1x,4f11.4,1x,f15.9,1x,f11.4,1x,f15.9,&
           1x,f15.9,1x,f11.4,1x,f15.9,1x,8f11.4,&
           1x,2f15.9,1x,69f11.4,1x,f11.1,1x,2f11.4)    
    close(luout2)
end subroutine flxwrite2
subroutine SAoutputwrite(lsao,mstep,SA_output)
!---------------------------------------------------------------------------
!output my sensitivity parameters
!original author:xucong
    use precision
    implicit none
    integer,INTENT(in):: lsao
    integer,INTENT(in):: mstep
    integer i
    real(r8):: SA_output(mstep)
    do i=1,mstep
       write(lsao,*) SA_output(i)
    end do
    close(lsao)
end subroutine SAoutputwrite
