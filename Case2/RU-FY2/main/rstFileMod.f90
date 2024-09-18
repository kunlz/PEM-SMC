subroutine rstTimeConstRead(lat_points,lhistTimeConst,nfcon,nftune,&
           numpatch,numpatch_lat,ixy_patch,jxy_patch,mxy_patch,wtxy_patch,&
           fcon,ftune)
!-----------------------------------------------------------------------
! Read in time constant variables reaquired by restart run
! Original author : Yongjiu Dai
!-----------------------------------------------------------------------
      use precision
      implicit none
      integer  lat_points           ! number of latitude points on model grid
      integer  lhistTimeConst       ! logical unit of file
      integer  nftune               ! number of clm tunable constants
      integer  nfcon                ! number of time constant variables
      integer  numpatch             ! total number of clm land points
      integer  numpatch_lat(lat_points) ! number of clm land points at lon. strip
      integer  ixy_patch (numpatch) ! longitude index for each patch point
      integer  jxy_patch (numpatch) ! latitude index for each patch point
      integer  mxy_patch (numpatch) ! patch subgrid index of lnd point
      real(r8) wtxy_patch(numpatch) ! subgrid weight for each patch point
      real(r8) fcon(numpatch,nfcon) ! time constant variables of CLM
      real(r8) ftune(nftune)        ! clm tunable constants
      read(lhistTimeConst) &
        numpatch_lat, &!number of clm land points at lon. strip
           ixy_patch, &!longitude index for each patch point
           jxy_patch, &!latitude index for each patch point
           mxy_patch, &!patch subgrid index of lnd point
          wtxy_patch   !subgrid weight for each patch point
      read(lhistTimeConst) &
           fcon,      &!time constant variables of CLM
           ftune       !clm tunable constants
      close(lhistTimeConst)
end subroutine rstTimeConstRead
subroutine rstTimeConstWrite(lat_points,lhistTimeConst,nfcon,nftune,&
           numpatch,numpatch_lat,ixy_patch,jxy_patch,mxy_patch,wtxy_patch,&
           fcon,ftune)
!-----------------------------------------------------------------------
! Write out time constant variables reaquired by restart run
! Original author : Yongjiu Dai
!-----------------------------------------------------------------------
      use precision
      implicit none
      integer  lat_points           ! number of latitude points on model grid
      integer  lhistTimeConst       ! logical unit of file
      integer  nftune               ! number of clm tunable constants
      integer  nfcon                ! number of time constant variables
      integer  numpatch             ! total number of clm land points
      integer  numpatch_lat(lat_points) ! number of clm land points at lon. strip
      integer  ixy_patch (numpatch) ! longitude index for each patch point
      integer  jxy_patch (numpatch) ! latitude index for each patch point
      integer  mxy_patch (numpatch) ! patch subgrid index of lnd point
      real(r8) wtxy_patch(numpatch) ! subgrid weight for each patch point
      real(r8) fcon(numpatch,nfcon) ! time constant variables of CLM
      real(r8) ftune(nftune)        ! clm tunable constants
      write(lhistTimeConst) &
         numpatch_lat, &!number of clm land points at lon. strip
            ixy_patch, &!longitude index for each patch point
            jxy_patch, &!latitude index for each patch point
            mxy_patch, &!patch subgrid index of lnd point
           wtxy_patch   !subgrid weight for each patch point
      write(lhistTimeConst) &
            fcon,      &!time constant variables of CLM
            ftune       !clm tunable constants
      close(lhistTimeConst)
end subroutine rstTimeConstWrite
subroutine rstTimeVarRead(lhistTimeVar,nfvar,numpatch,idate,fvar)
!-----------------------------------------------------------------------
! Read in time-varying variables reaquired by restart run
! Original author : Yongjiu Dai
!-----------------------------------------------------------------------
      use precision
      implicit none
      integer lhistTimeVar
      integer nfvar                 ! number of time varying variables
      integer numpatch              ! total number of clm land points
      integer idate(3)              ! calendar (year, julian day, seconds)
      real(r8) fvar(numpatch,nfvar) ! time varying variablesl
      read(lhistTimeVar) idate,fvar            
      close(lhistTimeVar)
end subroutine rstTimeVarRead
subroutine rstTimeVarWrite2(lhistTimeVar2,nfvar,numpatch,idate,fvar)
!-----------------------------------------------------------------------
! write out time-varing variables reaquired by restart run
! original author:Cong Xu
!-----------------------------------------------------------------------
   use precision
   implicit none
   integer lhistTimeVar2
   integer nfvar
   integer numpatch
   integer idate(3)
   real(r8) fvar(numpatch,nfvar)
   integer i
   do i=1,numpatch
      write(lhistTimeVar2,100) fvar(i,:)
100   format(1x,126f14.6)
   end do
   close(lhistTimeVar2)
end subroutine rstTimeVarWrite2
subroutine rstTimeVarWrite(lhistTimeVar,nfvar,numpatch,idate,fvar)
!-----------------------------------------------------------------------
! write out time-varying variables reaquired by restart run
! Original author : Yongjiu Dai
!-----------------------------------------------------------------------
      use precision
      implicit none
      integer lhistTimeVar
      integer nfvar                 ! number of time varying variables
      integer numpatch              ! total number of clm land points
      integer idate(3)              ! calendar (year, julian day, seconds)
      real(r8) fvar(numpatch,nfvar) ! time varying variablesl
      write(lhistTimeVar) idate,fvar
      close(lhistTimeVar)
end subroutine rstTimeVarWrite
!subroutine myVarianceWrite(lmyfile,rstfac)
!-------------------------------------------------------------------------
! write out my interesting variance
! Original author: Cong Xu
!-------------------------------------------------------------------------
!     use precision
!     implicit none
!     integer lmyfile !index of myfile
!    integer rstfac !facotr of water stress
!     write(lmyfile)
