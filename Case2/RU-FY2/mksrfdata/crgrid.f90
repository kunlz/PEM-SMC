  subroutine crgrid(edgen,edgee,edges,edgew,lon_points,lat_points,&
                    latixy,longxy,center_at_dateline,latn,lonw,area)
! ----------------------------------------------------------------------
! generate land model grid when mode is offline.
! surface grid edges -- grids do not have to be global.
! to allow this, grids must define the north, east, south, west edges:
!    edgen: northern edge of grid : > -90 and <= 90 (degrees)
!    edgee: eastern edge of grid  : > western edge and <= 180
!    edges: southern edge of grid : >= -90  and <  90
!    edgew: western edge of grid  : >= -180 and < 180
!
! region (global) latitude grid goes from:
!                          NORTHERN edge (POLE) to SOUTHERN edge (POLE)
! region (global) longitude grid starts at:
!                          WESTERN edge
!                         (DATELINE with western edge)
!
! west of Greenwich defined negative
! for global grids, the western edge of the longitude grid starts
! at the dateline
!
! ----------------------------------------------------------------------
      use precision
      implicit none
! arguments:
      integer,  intent(in) :: lon_points         ! input number of longitudes
      integer,  intent(in) :: lat_points         ! input number of latitudes
      real(r8), intent(in) :: edgen              ! northern edge of grid (degrees)
      real(r8), intent(in) :: edgee              ! eastern edge of grid (degrees)
      real(r8), intent(in) :: edges              ! southern edge of grid (degrees)
      real(r8), intent(in) :: edgew              ! western edge of grid (degrees)
      real(r8), intent(out) :: latixy(lon_points,lat_points) ! input grid: latitude (degrees)
      real(r8), intent(out) :: longxy(lon_points,lat_points) ! input grid: longitude (degrees)
      real(r8), intent(out) :: latn(lat_points+1) ! grid cell latitude, northern edge (degrees)
      real(r8), intent(out) :: lonw(lon_points+1) ! grid cell longitude, western edge (degrees)
      real(r8), intent(out) :: area(lon_points,lat_points)   ! input grid: cell area
      logical, intent(out) :: center_at_dateline  ! logical for 1st grid center at dateline
! local variables:
      real(r8) :: dx                             ! land model cell width
      real(r8) :: dy                             ! land model cell length
      integer  :: i,j                            ! indices
! ----------------------------------------------------------------------
! determine grid longitudes and latitudes in increments of dx and dy
      dx = (edgee-edgew)/lon_points
      dy = (edgen-edges)/lat_points
      do j = 1, lat_points
         do i = 1, lon_points
            latixy(i,j) = edgen - (2*j-1)*dy/2.
            longxy(i,j) = edgew + (2*i-1)*dx/2.
         end do
      end do
      center_at_dateline = .false.
! define edges and area of land model grid cells
      call celledge (lat_points,lon_points,longxy,latixy,&
                     edgen,edgee,edges,edgew,center_at_dateline,latn,lonw)
      call cellarea (lat_points,lon_points,latn,lonw,&
                                edgen,edgee,edges,edgew,area)
  end subroutine crgrid
! ----------------------------------------------------------------------
! EOP
