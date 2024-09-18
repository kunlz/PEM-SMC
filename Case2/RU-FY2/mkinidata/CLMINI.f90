PROGRAM CLMINI
! ======================================================================
! The Common Land Model was developed in cooperation with
!     Beijing Normal University and IRI         (Dai)
!     Georgia Institute of Technology           (Dickinson)
!     National Center for Atmospheric Research  (Bonan, Oleson)
!     University of Arizona                     (Zeng)
!     University of Texas at Austin             (Yang)
!     GSFC/NASA                                 (Houser, Bosilovich)
!     COLA                                      (Dirmeyer, Schlosser)
!     Colorado State University at Fort Collins (Denning, Baker)
!
! Reference:
!     [1] Dai et al., 2003: The Common Land Model (CLM).
!         Bull. of Amer. Meter. Soc., 84: 1013-1023
!     [2] Dai et al., 2004: A two-big-leaf model for canopy temperature,
!         photosynthesis and stomatal conductance. Journal of Climate
!
!     Created by Yongjiu Dai Februay 2004
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
! ----------------local variables ---------------------------------
      integer, parameter :: nl_soil  = nl_soil_  ! number of soil layers
      integer, parameter :: maxsnl   = maxsnl_   ! max number of snow layers
      integer, parameter :: nfcon    = nfcon_    ! number of time constant variables
      integer, parameter :: nftune   = nftune_   ! number of clm tunable constants
      integer, parameter :: nfvar    = nfvar_    ! number of time varying variables
      integer, parameter :: nfldv    = nfldv_    ! number of output fluxes
      integer, parameter :: maxpatch = maxpatch_ ! max number of patches in a grid
      integer :: start_yr       ! starting date for run in year
      integer :: start_jday     ! starting date for run in julian day
      integer :: start_sec      ! starting time of day for run in seconds
      logical :: greenwich      ! true: greenwich time, false: local time
      integer :: lusrf          ! logical unit number of surface data
      integer :: lulai          ! logical unit number of LAI data
      integer :: lusoil         ! logical unit number of soil initial
      integer :: lumet          ! logical unit number of meteorological forcing
      integer :: lhistTimeConst ! logical unit number of restart time-invariant file
      integer :: lhistTimeConst2 !logical unit number of restart time-invariant
      integer :: lhistTimeVar   ! logical unit number of restart time-varying file
      integer :: lhistTimeVar2
      integer :: luout          ! logical unit number of output
      integer :: lon_points     ! number of longitude points on model grid
      integer :: lat_points     ! number of latitude points on model grid
      integer :: numpatch       ! total number of patches of grids
      real    :: deltim         ! time step of model (second)
      integer :: mstep          ! model step for simulation [-]
      character(LEN=256) :: site           ! site name
      character(LEN=256) :: fsurdat        ! file name of surface data
      character(LEN=256) :: flaidat        ! file name of time-varying vegetation data
      character(LEN=256) :: fsoildat       ! file name of soil state initial data
      character(LEN=256) :: fmetdat        ! file name meteorological data
      character(LEN=256) :: fhistTimeConst ! file name of time-invariant file
      character(LEN=256) :: fhistTimeVar   ! file name of time-varying file
      character(LEN=256) :: foutdat        ! file name of output file
      character(LEN=256) :: finfolist      ! file name of run information
      character(LEN=256) :: cdate          ! character for date
      character(LEN=256) :: fhistTimeVar_name
      character(LEN=256) :: fhistTimeVar_name2
      character(LEN=256) :: fsoildat_name
      real(r8), allocatable :: fldxy(:,:,:)
      namelist /clminiexp/ site,&
                           greenwich,start_yr,start_jday,start_sec,&
                           fsurdat,flaidat,fsoildat,fmetdat,&
                           fhistTimeConst,fhistTimeVar,&
                           foutdat,finfolist,&
                           lon_points,lat_points,deltim,mstep
! ----------------------------------------------------------------------
      read (5,clminiexp)
      lusrf          = 110
      lulai          = 120
      lusoil         = 130
      lumet          = 140
      lhistTimeConst = 150
      lhistTimeConst2 =15
      lhistTimeVar   = 160
      lhistTimeVar2  =18
      luout          = 170
      allocate (fldxy(lat_points,lon_points,nfldv))
      OPEN(unit=lusrf,file=fsurdat,status='old',form='unformatted',action='read')
      write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') start_yr,start_jday,start_sec
      OPEN(unit=lhistTimeConst,file=fhistTimeConst,status='unknown',&
                               form='unformatted',action='write')
!      OPEN(unit=lhistTimeConst2,file='/home/u120220909911/RU-FY2/Third/psuade/NEE/LH/RU-FY2/output/rstTimeConst.txt',status='new',&
!                               form='formatted',action='write')
!      fhistTimeVar_name = trim(fhistTimeVar)//'-'//trim(cdate)
       fhistTimeVar_name = fhistTimeVar
!      fhistTimeVar_name2= trim(fhistTimeVar)//'-'//trim(cdate)//'.txt'
      OPEN(unit=lhistTimeVar,file=fhistTimeVar_name,status='unknown',&
                               form='unformatted',action='write')
!      OPEN(unit=lhistTimeVar2,file=fhistTimeVar_name2,status='new',&
!                               form='formatted',action='write')
      CALL initialize (start_yr,start_jday,start_sec,greenwich,&
                       lon_points,lat_points,maxpatch,nl_soil,maxsnl,&
                       lusrf,&
                       lhistTimeConst,lhistTimeConst2,lhistTimeVar,&
                       lhistTimeVar2,nftune,nfcon,nfvar,nfldv,numpatch,&
                       fldxy)
      OPEN(180,file=finfolist,form='formatted')
      write(180,*) '&clmexp'
      write(180,*) 'site           = ', "'"//trim(site)//"'"             !1
      write(180,*) 'flaidat        = ', "'"//trim(flaidat)//"'"          !2
      write(180,*) 'fmetdat        = ', "'"//trim(fmetdat)//"'"          !3
      write(180,*) 'fhistTimeConst = ', "'"//trim(fhistTimeConst)//"'"   !4
      write(180,*) 'fhistTimeVar   = ', "'"//trim(fhistTimeVar_name)//"'"!5
      write(180,*) 'foutdat        = ', "'"//trim(foutdat)//"'"          !6
      write(180,*) 'lhistTimeConst = ', lhistTimeConst                   !7
      write(180,*) 'lhistTimeVar   = ', lhistTimeVar                     !8
      write(180,*) 'lulai          = ', lulai                            !9
      write(180,*) 'lumet          = ', lumet                            !10
      write(180,*) 'luout          = ', luout                            !11
      write(180,*) 'lon_points     = ', lon_points                       !12
      write(180,*) 'lat_points     = ', lat_points                       !13
      write(180,*) 'numpatch       = ', numpatch                         !14
      write(180,*) 'deltim         = ', deltim                           !15
      write(180,*) 'mstep          = ', mstep                            !16
      write(180,*) '/'
      deallocate (fldxy)
!      write(6,*) 'CLM Initialization Execution Completed'
END PROGRAM CLMINI
! ----------------------------------------------------------------------
! EOP
