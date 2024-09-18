PROGRAM mksrfdata
! ======================================================================
! Surface grid edges:
! The model domain was defined with the north, east, south, west edges:
!          edgen: northern edge of grid : > -90 and <= 90 (degrees)
!          edgee: eastern edge of grid  : > western edge and <= 180
!          edges: southern edge of grid : >= -90  and <  90
!          edgew: western edge of grid  : >= -180 and < 180
!
! Region (global) latitude grid goes from:
!                 NORTHERN edge (POLE) to SOUTHERN edge (POLE)
! Region (global) longitude grid starts at:
!                 WESTERN edge (DATELINE with western edge)
!                 West of Greenwich defined negative for global grids,
!                 the western edge of the longitude grid starts at the dateline
!
! Land surface properties was mapped from USGS "raw" data with
!                 30 arc seconds resolution:
!              -  elevation height
!              -  land water mask
!              -  land cover type
!              -  soil texture of top soil layer (0-30cm) (FAO+STATSGO)
!              -  soil texture of bottom soil layer (30-100cm) (FAO+STATSGO)
!
! Following surface dataset will be created:
!             (1) model grid (longitude, latitude)
!             (2) soil color
!             (3) depth to bedrock
!             (4) soil texture (with vertical profile)
!             (5) land cover type (with patches in grid)
!             (6) fraction of patches
!
! Created by Yongjiu Dai
! ======================================================================
      use precision
      implicit none
!----------------------------------------------------------------------
! Define the dimension of model array
!----------------------------------------------------------------------
      integer nl_soil_       ! number of soil layers
      integer maxsnl_        ! max number of snow layers
      integer nfcon_         ! number of time constant variables
      integer nftune_        ! number of clm tunable constants
      integer nfvar_         ! number of time varying variables
      integer nforc_         ! number of forcing variables
      integer nfldv_         ! number of output fluxes
      integer nflai_         ! number of leaf time varying variables
      integer maxpatch_      ! number of clm grid points
      integer nlandcateg_    ! number of land cover categories
      integer nsoilcateg_    ! number of soil texture categories
      parameter(nl_soil_    = 10)
      parameter(maxsnl_     = -5)
      parameter(nfcon_      = 9*nl_soil_+29)
      parameter(nftune_     = 14)
      parameter(nfvar_      = 5*(nl_soil_-maxsnl_)+51)
      parameter(nforc_      = 18)
      parameter(nfldv_      = 92)
      parameter(nflai_      = 4)
      parameter(maxpatch_   = 25)
      parameter(nlandcateg_ = 25)
      parameter(nsoilcateg_ = 17)
! ----------------------------------------------------------------------
      integer, parameter :: nlandcateg = nlandcateg_ ! number of land cover category
      integer, parameter :: nsoilcateg = nsoilcateg_ ! number of soil texture category
      integer, parameter :: nl_soil = nl_soil_   ! number of soil layers
      integer, parameter :: maxpatch = maxpatch_ ! number of soil texture category
! local variables:
      real(r8), allocatable :: latn(:)        ! grid cell latitude, northern edge (degrees)
      real(r8), allocatable :: lonw(:)        ! grid cell longitude, western edge (degrees)
      real(r8), allocatable :: area(:,:)      ! input grid: cell area
      integer,  allocatable :: numtype(:,:,:) ! numbers of land categories in grid
      integer,  allocatable :: numsola(:,:,:) ! numbers of srf soil categories in grid
      integer,  allocatable :: numsolb(:,:,:) ! numbers of deep soil categories in grid
      real(r8), allocatable :: latixy(:,:)    ! grid latitude at center of points (degrees)
      real(r8), allocatable :: longxy(:,:)    ! grid longitude at center of points (degrees)
      real(r8), allocatable :: sand2d(:,:,:)  ! percentage of sand
      real(r8), allocatable :: clay2d(:,:,:)  ! percentage of clay
      real(r8), allocatable :: rock2d(:,:)    ! depth to bed rock
      integer , allocatable :: soic2d(:,:)    ! soil color
      integer , allocatable :: surf2d(:,:,:)  ! land cover type
      real(r8), allocatable :: fpatch2d(:,:,:)!subgrid weights
      character(LEN=72) :: fgridname  ! file names of grid longitude/latitude
      character(LEN=72) :: fdemname   ! file names of USGS 30s elevation
      character(LEN=72) :: fmaskname  ! file names of USGS 30s land water mask
      character(LEN=72) :: flandname  ! file names of USGS 30s land cover type
      character(LEN=72) :: fsolaname  ! file names of USGS 30s soil texture srf 30 cm
      character(LEN=72) :: fsolbname  ! file names of USGS 30s soil texture srf 30-100 cm
      character(LEN=72) :: fsurdat    ! file names of CLM surface information
      integer           :: lon_points ! number of input data longitudes
      integer           :: lat_points ! number of input data latitudes
      real(r8)          :: edgen      ! northern edge of grid (degrees)
      real(r8)          :: edgee      ! eastern edge of grid (degrees)
      real(r8)          :: edges      ! southern edge of grid (degrees)
      real(r8)          :: edgew      ! western edge of grid (degrees)
      logical           :: center_at_dateline  ! logical for 1st grid
      character(LEN=72) :: lndname(5) ! file names of USGS "raw" data
      integer           :: iunit(5)   ! logical number of USGS "raw" data
      integer           :: lugrid     ! logical unit number of grid file
      integer           :: lusrf      ! logical unit number of CLM surface info
      real(r8)          :: zsoi(nl_soil)    ! soil layer depth [m]
      real(r8)          :: dzsoi(nl_soil)   ! soil node thickness [m]
      real(r8)          :: zsoih(0:nl_soil) ! interface level below a zsoi level [m]
      real(r8)          :: a1, a2, a3, a4, f1, f2, ferr, ferr1, ferr2
      integer           :: i, j, k, l, ll, np, np1, np2, loca(1)
! relative amounts of sand (s), and clay (c) in the < 2 mm fraction of
! the component layer was then estimated using table:
      real(r8), dimension(17) ::  s = (/92.,82.,58.,17.,10.,43.,58.,10.,32.,52.,&
                                         6.,22., 0., 0., 0., 0., 0./)
      real(r8), dimension(17) ::  c = (/ 3., 6.,10.,13., 5.,18.,27.,34.,34.,42.,&
                                        47.,58., 0., 0., 0., 0., 0./)
      namelist /mksrfexp/ fgridname, fdemname, fmaskname, &
                          flandname, fsolaname, fsolbname, fsurdat, &
                          lon_points, lat_points, edgen, edgee, edges, edgew
! ----------------------------------------------------------------------
      read(5,mksrfexp)
      allocate (latn(lat_points+1))
      allocate (lonw(lon_points+1))
      allocate (area(lon_points,lat_points))
      allocate (numtype(lon_points,lat_points,nlandcateg)) 
      allocate (numsola(lon_points,lat_points,nsoilcateg)) 
      allocate (numsolb(lon_points,lat_points,nsoilcateg)) 
      allocate (latixy(lon_points,lat_points))
      allocate (longxy(lon_points,lat_points))
      allocate (soic2d(lon_points,lat_points))     
      allocate (rock2d(lon_points,lat_points))
      allocate (sand2d(lon_points,lat_points,nl_soil))      
      allocate (clay2d(lon_points,lat_points,nl_soil))
      allocate (surf2d(lon_points,lat_points,maxpatch))     
      allocate (fpatch2d(lon_points,lat_points,maxpatch))
! Creat model grid
      CALL crgrid (edgen,edgee,edges,edgew,lon_points,lat_points,&
                         latixy,longxy,center_at_dateline,latn,lonw,area)
! Read in the 30 arc second
! Mapping land cover type and soil type to model grid
      iunit(1) = 11
      iunit(2) = 12
      iunit(3) = 13
      iunit(4) = 14
      iunit(5) = 15
      lndname(1) = fdemname
      lndname(2) = fmaskname
      lndname(3) = flandname
      lndname(4) = fsolaname
      lndname(5) = fsolbname
      CALL rdlanddata (lon_points,lat_points,nlandcateg,nsoilcateg,&
                       edgen,edgee,edges,edgew,center_at_dateline,&
                       latn,lonw,iunit,lndname,numtype,numsola,numsolb)
! CLM default soil depth, thickness, depth of layer interface
      do k = 1, nl_soil
         zsoi(k) = 0.025*(exp(0.5*(k-0.5))-1.)  ! node depths
      end do
      dzsoi(1) = 0.5*(zsoi(1)+zsoi(2))          ! =zsoih(1)
      dzsoi(nl_soil) = zsoi(nl_soil)-zsoi(nl_soil-1)
      do k = 2, nl_soil-1
         dzsoi(k) = 0.5*(zsoi(k+1)-zsoi(k-1))   ! thickness b/n two interfaces
      end do
      zsoih(0) = 0.
      zsoih(nl_soil) = zsoi(nl_soil) + 0.5*dzsoi(nl_soil)
      do k = 1, nl_soil-1
         zsoih(k)= 0.5*(zsoi(k)+zsoi(k+1))      ! interface depths
      enddo
! Land cover type, pactes fraction, and soil properties
      do j = 1, lat_points
         do i = 1, lon_points
! LAND cover type with fractions of grid patches
            np = sum(numtype(i,j,:))
!            print*,'mksrfdata.F90-test-veg_numtype',numtype(i,j,:)
            do L = 1, nlandcateg
               surf2d(i,j,L) = L-1     
               fpatch2d(i,j,L) = float(numtype(i,j,L))/float(np)
               if(fpatch2d(i,j,L).lt.0. .or. fpatch2d(i,j,L).gt.1.)then
                  print*, i, j, L, np, numtype(i,j,L), fpatch2d(i,j,L)
                  print*, 'MKSRFDATA error 1'
                  call abort
               endif
            enddo
! ocean coast grid merging: ocean>50% => ocean full grid
! ocean<50%: discard ocean patch and carry to the largest land patch of grid
            if(fpatch2d(i,j,1).gt.0.5)then
               fpatch2d(i,j,1) = 1.0
               fpatch2d(i,j,2:) = 0.
            else
               Loca = maxloc(fpatch2d(i,j,2:))
               L=loca(1) + 1
               fpatch2d(i,j,L) = fpatch2d(i,j,L) + fpatch2d(i,j,1)
               fpatch2d(i,j,1) = 0.
            endif
! all covers smaller than 1% are discarded and carried to
! the biggest patch of grid, with the exception of grass that
! is assumed to always be at least 1% (the deficits are compensated by
! deducting value from the largest patch) for comparing the routine
! meteorological observation (grass_id+1=8 for USGS classification)
            do L = 1, nlandcateg
               if(fpatch2d(i,j,L).lt.1./100. .and. L.ne.8)then
                  Loca = maxloc(fpatch2d(i,j,:))
                  LL=Loca(1) 
                  if(LL.eq.L) stop 'evenly distribution in grid'
                  fpatch2d(i,j,LL) = fpatch2d(i,j,LL) + fpatch2d(i,j,L)
                  fpatch2d(i,j,L) = 0.
               endif
            enddo
            if(fpatch2d(i,j,1).lt.1.0e-6)then
               if(fpatch2d(i,j,8).lt.1./100.)then 
                  Loca = maxloc(fpatch2d(i,j,:))
                  LL=Loca(1) 
                  fpatch2d(i,j,LL) = fpatch2d(i,j,LL) + fpatch2d(i,j,8) - 1./100.
                  fpatch2d(i,j,8) = 1./100.
               endif
            endif
            ferr = sum(fpatch2d(i,j,:)) - 1.
            if(abs(ferr).gt.1.0e-6)then
               print*, 'MKSRFDATA 2: fractional not eq 1', ferr
               call abort
            endif
! SOIL sand and clay percentages
            a1 = 0.
            a2 = 0.
            a3 = 0.
            a4 = 0.
            ferr1 = -1.
            ferr2 = -1.
            np1 = sum(numsola(i,j,:))
!            print*,'mksrfdata.F90:test-soil-numsola',numsola(i,j,:)
            np2 = sum(numsolb(i,j,:))
!            print*,'mksrfdata.F90:test-soil-numsolb',numsolb(i,j,:)
            do L = 1, nsoilcateg
               f1 = float(numsola(i,j,L))/float(np1)
               f2 = float(numsolb(i,j,L))/float(np2)
               a1 = a1 + f1*s(L)     
               a2 = a2 + f1*c(L)
               a3 = a3 + f2*s(L)     
               a4 = a4 + f2*c(L)
               ferr1 = ferr1 + f1
               ferr2 = ferr2 + f2
               if(f1.lt.0. .or. f1.gt.1.)then
                  print*, i, j, L, np1, numsola(i,j,L), f1
                  print*, 'MKSRFDATA error 3'
                  call abort
               endif
               if(f2.lt.0. .or. f2.gt.1.)then
                  print*, i, j, L, np2, numsolb(i,j,L), f1
                  print*, 'MKSRFDATA error 4'
                  call abort
               endif
            enddo
            if(abs(ferr1).gt.1.0e-6)then
               print*, 'MKSRFDATA 5: fractional not eq 1', ferr1
               call abort
            endif
            if(abs(ferr2).gt.1.0e-6)then
               print*, 'MKSRFDATA 6: fractional not eq 1', ferr2
               call abort
            endif
            do k = 1, nl_soil            
               if(zsoih(k).lt.0.3)then
                  sand2d(i,j,k) = a1 
                  clay2d(i,j,k) = a2
               else
                  sand2d(i,j,k) = a3 
                  clay2d(i,j,k) = a4
               endif
               if(sand2d(i,j,k).gt.100. .or. clay2d(i,j,k).gt.100.)then
                  print*, '%sand = ', sand2d(i,j,k), '%clay = ', clay2d(i,j,k)
                  call abort
               endif
            enddo
! SOIL color index: we define the color index depended on the
! dominant land cover type (based on USGS land classification).
! when you have the data, don't forget to replace it with the real values
! the index are from light to dark (1-light <-> 8-dark)
                        soic2d(i,j) = 5  
            Loca = maxloc(fpatch2d(i,j,:))
            L=loca(1) 
            if(L.eq. 1) soic2d(i,j) = 1  ! 0  Ocean (not used)
            if(L.eq. 2) soic2d(i,j) = 1  ! 1  Urban and Built-Up Land
            if(L.eq. 3) soic2d(i,j) = 3  ! 2  Dryland Cropland and Pasture
            if(L.eq. 4) soic2d(i,j) = 5  ! 3  Irrigated Cropland and Pasture
            if(L.eq. 5) soic2d(i,j) = 6  ! 4  Mixed Dryland/Irrigated Cropland and Pasture
            if(L.eq. 6) soic2d(i,j) = 7  ! 5  Cropland/Grassland Mosaic
            if(L.eq. 7) soic2d(i,j) = 7  ! 6  Cropland/Woodland Mosaic
            if(L.eq. 8) soic2d(i,j) = 8  ! 7  Grassland
            if(L.eq. 9) soic2d(i,j) = 3  ! 8  Shrubland
            if(L.eq.10) soic2d(i,j) = 4  ! 9  Mixed Shrubland/Grassland
            if(L.eq.11) soic2d(i,j) = 2  !10  Savanna
            if(L.eq.12) soic2d(i,j) = 8  !11  Deciduous Broadleaf Forest
            if(L.eq.13) soic2d(i,j) = 8  !12  Deciduous Needleleaf Forest
            if(L.eq.14) soic2d(i,j) = 8  !13  Evergreen Broadleaf Forest
            if(L.eq.15) soic2d(i,j) = 8  !14  Evergreen Needleleaf Forest
            if(L.eq.16) soic2d(i,j) = 8  !15  Mixed Forest
            if(L.eq.17) soic2d(i,j) = 1  !16  Water Bodies (not used)
            if(L.eq.18) soic2d(i,j) = 8  !17  Herbaceous Wetland
            if(L.eq.19) soic2d(i,j) = 8  !18  Wooded Wetland
            if(L.eq.20) soic2d(i,j) = 1  !19  Barren or Sparsely Vegetated
            if(L.eq.21) soic2d(i,j) = 7  !20  Herbaceous Tundra
            if(L.eq.22) soic2d(i,j) = 7  !21  Wooded Tundra
            if(L.eq.23) soic2d(i,j) = 7  !22  Mixed Tundra
            if(L.eq.24) soic2d(i,j) = 7  !23  Bare Ground Tundra
            if(L.eq.25) soic2d(i,j) = 1  !24  Snow or Ice (not used)
! SOIL depth to bedrock => when you have the data,
! don't forget to change it with a real value
            rock2d(i,j) = 6.0
         enddo
      enddo
! Write out the surface type data
! PLEASE PAY ATTENTION: when you read the following file, the variables
! have to be defined with the same precision as in this program, .i.e., (r8)
! I struggled for that kind of confusion for several days.
      lusrf = 16
      OPEN(lusrf,file=fsurdat,form='unformatted',status='unknown')
      do j = 1, lat_points
         do i = 1, lon_points
         write(lusrf) latixy(i,j),&
                      longxy(i,j),&
                      soic2d(i,j),&
                      rock2d(i,j),&
                     (sand2d(i,j,l),l=1,nl_soil),&
                     (clay2d(i,j,l),l=1,nl_soil),&
                     (surf2d(i,j,l),l=1,maxpatch),&
                   (fpatch2d(i,j,l),l=1,maxpatch)
         end do
      end do
!      print*, 'Successful in surface data making'
END PROGRAM mksrfdata
! ----------------------------------------------------------------------
! EOP
