  subroutine cellarea(lat_points,lon_points,latn,lonw,&
                                 edgen,edgee,edges,edgew,area)
! ----------------------------------------------------------------------
! Area of grid cells (square kilometers) - regional grid
! (can become global as special case)
      use precision
      implicit none
! arguments
      integer , intent(in) :: lat_points        ! number of latitude points
      integer , intent(in) :: lon_points        ! number of longitude points
      real(r8), intent(in) :: latn(lat_points+1)! grid cell latitude, northern edge (deg)
      real(r8), intent(in) :: lonw(lon_points+1)! grid cell longitude, western edge (deg)
      real(r8), intent(in) :: edgen             ! northern edge of grid (deg)
      real(r8), intent(in) :: edges             ! southern edge of grid (deg)
      real(r8), intent(in) :: edgew             ! western edge of grid (deg)
      real(r8), intent(in) :: edgee             ! eastern edge of grid (deg)
      real(r8), intent(out):: area(lon_points,lat_points) ! cell area (km**2)
! local variables
      integer i,j                               ! indices
      real(r8) deg2rad                          ! pi/180
      real(r8) global                           ! summed area
      real(r8) dx                               ! cell width: E-W
      real(r8) dy                               ! cell width: N-S
      real(r8) pi                               ! 3.14159265358979323846
      real(r8) re                               ! radius of earth (km)
      real(r8) error                            ! true area for error check
! -----------------------------------------------------------------------
      re = 6.37122e6 * 0.001                    ! kilometer
      pi = 4.*atan(1.)
      deg2rad = pi/180.
      global = 0.
      do j = 1, lat_points
         do i = 1, lon_points
            dx = (lonw(i+1)-lonw(i))*deg2rad
            if(latn(j)>latn(j+1)) then          ! north to south grid
               dy = sin(latn(j)*deg2rad) - sin(latn(j+1)*deg2rad)
            else                                ! south to north grid
               dy = sin(latn(j+1)*deg2rad) - sin(latn(j)*deg2rad)
            end if
            area(i,j) = dx*dy*re*re
            global = global + area(i,j)
         end do
      end do
! make sure total area from grid cells is same as area of grid
! as defined by its edges
      dx = (edgee - edgew) * deg2rad
      dy = sin(edgen*deg2rad) - sin(edges*deg2rad)
      error = dx*dy*re*re
      if(abs(global-error)/error > 0.00001) then
         write (6,*) 'CELLAREA error: correct area is ',error, &
              ' but summed area of grid cells is ', global
         call abort
      end if
      return
  end subroutine cellarea
! -----------------------------------------------------------------------
! EOP
