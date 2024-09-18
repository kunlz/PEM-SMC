  subroutine celledge (lat_points,lon_points,longxy,latixy,&
                       edgen,edgee,edges,edgew,center_at_dateline,latn,lonw)
! -----------------------------------------------------------------------
! southern and western edges of grid cells - regional grid
! (can become global as special case)
!
! Latitudes -- southern/northern edges for each latitude strip.
! Grids oriented North to South: the southern
! and northern edges of latitude strip [j] are:
!        northern = latn(j  )
!        southern = latn(j+1)
! [latn] must be dimensioned latn(lat_points+1)
!
! Longitudes -- western edges. Longitudes for the western edge of the
! cells must increase continuously and span 360 degrees. Assume that
! grid starts at Dateline with western edge on Dateline Western edges
! correspond to [longxy] (longitude at center of cell) and range from
! -180 to 180 with negative longitudes west of Greenwich.
! Partial grids that do not span 360 degrees are allowed so long as they
! have the convention of Grid 1 with
!      western edge of grid: >= -180 and < 180
!      eastern edge of grid: > western edge  and <= 180
! [lonw] must be dimensioned lonw(lon+1,lat) because each latitude
! strip can have variable longitudinal resolution
!
! -----------------------------------------------------------------------
      use precision
      implicit none
! arguments:
      integer , intent(in) :: lat_points        ! dimension: number of latitude points
      integer , intent(in) :: lon_points        ! dimension: number of longitude points
      real(r8), intent(in) :: longxy(lon_points,lat_points) ! longitude at center of grid cell
      real(r8), intent(in) :: latixy(lon_points,lat_points) ! latitude at center of grid cell
      real(r8), intent(in) :: edgen             ! northern edge of grid (degrees)
      real(r8), intent(in) :: edgee             ! eastern edge of grid (degrees)
      real(r8), intent(in) :: edges             ! southern edge of grid (degrees)
      real(r8), intent(in) :: edgew             ! western edge of grid (degrees)
      real(r8), intent(out):: latn(lat_points+1)! grid cell latitude, northern edge (degrees)
      real(r8), intent(out):: lonw(lon_points+1)! grid cell longitude, western edge (degrees)
      logical, intent(in) :: center_at_dateline ! logical for 1st grid center at dateline
! local variables
      integer i,j                               ! indices
      real(r8) dx2                              ! 1/2 grid width
!------------------------------------------------------------------------
! latitudes
      latn(1) = edgen 
      latn(lat_points+1) = edges 
      do j = 2, lat_points
         latn(j) = (latixy(1,j-1) + latixy(1,j)) / 2.
      end do
! longitudes
! western edge of first grid cell -- since grid starts with western
! edge on Dateline, lonw(1,j)=-180.
! On a global grid lonw(numlon+1,j)=lonw(1,j)+180.
      lonw(1) = edgew
      if(center_at_dateline)then
         lonw(2) = (longxy(1,1)+longxy(2,1))/2.
      else
         dx2 = longxy(1,1) - lonw(1)
         lonw(2) = longxy(1,1) + dx2
      endif
      do i = 3, lon_points+1
         dx2 = longxy(i-1,1) - lonw(i-1)
         lonw(i) = longxy(i-1,1) + dx2
      end do
    return
  end subroutine celledge
! -----------------------------------------------------------------------
! EOP
