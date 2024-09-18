  subroutine rdgrid(lugrid,edgen,edgee,edges,edgew,lon_points,lat_points,&
                           latixy,longxy,center_at_dateline,latn,lonw,area)
! ----------------------------------------------------------------------
! read land model grid
! region (global) latitude grid goes from:
!                          NORTHERN edge (POLE) to SOUTHERN edge (POLE)
! region (global) longitude grid starts at:
!                          WESTERN edge
!                         (DATELINE with western edge)
! surface grid edges -- grids do not have to be global.
! to allow this, grids must define the north, east, south, west edges:
!    edgen: northern edge of grid : > -90 and <= 90 (degrees)
!    edgee: eastern edge of grid  : > western edge and <= 180
!    edges: southern edge of grid : >= -90  and <  90
!    edgew: western edge of grid  : >= -180 and < 180
!
! ----------------------------------------------------------------------
      use precision
      implicit none
! arguments:
      integer, intent(in) :: lugrid     ! logical unit number of grid file
      integer, intent(in) :: lon_points ! number of input data longitudes
      integer, intent(in) :: lat_points ! number of input data latitudes
      real(r8), intent(in) :: edgen     ! northern edge of grid (degrees)
      real(r8), intent(in) :: edges     ! southern edge of grid (degrees)
      real(r8), intent(in) :: edgee     ! eastern edge of grid (degrees)
      real(r8), intent(in) :: edgew     ! western edge of grid (degrees)
      real(r8), intent(out) :: latixy(lon_points,lat_points) ! input grid: latitude (deg)
      real(r8), intent(out) :: longxy(lon_points,lat_points) ! input grid: longitude (deg)
      real(r8), intent(out) :: latn(lat_points+1) ! grid cell latitude, northern edge (deg)
      real(r8), intent(out) :: lonw(lon_points+1) ! grid cell longitude, western edge (deg)
      real(r8), intent(out) :: area(lon_points,lat_points)   ! input grid: cell area
      logical, intent(out) :: center_at_dateline  ! input grid: cell area
! local variables:
      integer  :: i,j                   ! indices
      real(r8) :: lon(lon_points)       ! read-in longitude array (full grid)
      real(r8) :: lat(lat_points)       ! read-in latitude array (full grid)
      real(r8) :: dx
! ----------------------------------------------------------------------
      write (6,*) 'Attempting to read land grid data .....'
      write (6,'(72a1)') ("-",i=1,60)
      read(lugrid,*) lat
      read(lugrid,*) lon
      print*, lat
      print*, lon
! determine grid longitude and latitudes
      do j = 1, lat_points
         longxy(:,j) = lon(:)
      end do
      center_at_dateline = .false.
! for grids starting at greenwich and centered on Greenwich'
      if(abs(longxy(1,1)+180.).lt.1.e-6)then
         center_at_dateline = .true.
      endif
      do i = 1, lon_points
         latixy(i,:) = lat(:)
      end do
! determine if grid has pole points - if so, make sure that north pole
! is non-land and south pole is land
      if(abs((latixy(1,lat_points)-90.)) < 1.e-6) then
         write(6,*) 'GRID: model has pole_points' 
      else
         write(6,*) 'GRID: model does not have pole_points' 
      endif
! define land grid edges and grid cell areas
      call celledge (lat_points,lon_points,longxy,latixy,&
                     edgen,edgee,edges,edgew,center_at_dateline,latn,lonw)
      call cellarea (lat_points,lon_points,latn,lonw,&
                                edgen,edgee,edges,edgew,area)
  end subroutine rdgrid
! ----------------------------------------------------------------------
! EOP
