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
