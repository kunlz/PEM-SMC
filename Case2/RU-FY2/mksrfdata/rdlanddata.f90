  subroutine rdlanddata(ix,jx,nlandcateg,nsoilcateg,&
                        edgen,edgee,edges,edgew,center_at_dateline,&
                        latn,lonw,iunit,lndname,numtype,numsola,numsolb)
! ----------------------------------------------------------------------
! Creates land model surface dataset from original "raw" data files -
!     USGS data with 30 arc seconds resolution:
!  -  elevation height
!  -  land cover type
!  -  soil texture (FAO+STATSGO)
!
! created by Yongjiu Dai
! ----------------------------------------------------------------------
      use precision
      implicit none
! arguments:
      integer, intent(in) :: ix           ! number of model longitude grid points
      integer, intent(in) :: jx           ! model  of model latitude grid points
      integer, intent(in) :: nlandcateg   ! number of land cover category
      integer, intent(in) :: nsoilcateg   ! number of soil texture category
      real(r8), intent(in) :: edgen       ! northern edge of grid (degrees)
      real(r8), intent(in) :: edgee       ! eastern edge of grid (degrees)
      real(r8), intent(in) :: edges       ! southern edge of grid (degrees)
      real(r8), intent(in) :: edgew       ! western edge of grid (degrees)
      real(r8),intent(in) :: latn(jx+1)   ! grid cell latitude, northern edge (degrees)
      real(r8),intent(in) :: lonw(ix+1)   ! grid cell longitude, western edge (degrees)
      integer, intent(in) :: iunit(5)     ! logical number of USGS "raw" data
      character(LEN=72), intent(in) :: lndname(5) ! file names of USGS "raw" data
      logical, intent(in) :: center_at_dateline   ! logical for 1st grid center at dateline
      integer, intent(out) :: numtype(ix,jx,nlandcateg) ! numbers of land categories in grid
      integer, intent(out) :: numsola(ix,jx,nsoilcateg) ! numbers of srf soil categories in grid
      integer, intent(out) :: numsolb(ix,jx,nsoilcateg) ! numbers of deep soil categories in grid
! local variables:
      integer, parameter :: nlat=21600    ! 180*(60*2)
      integer, parameter :: nlon=43200    ! 360*(60*2)
      character(LEN=2) inter1_chr(nlon)   ! terrain height
      character(LEN=1) landmask_chr(nlon) ! land water mask
      character(LEN=1) landtype_chr(nlon) ! land cover type
      character(LEN=1) landsola_chr(nlon) ! soil texture type in upper 30cm
      character(LEN=1) landsolb_chr(nlon) ! soil texture type in 30-100cm
      integer inter1(nlon)
      integer landmask(nlon)
      integer landtype(nlon)
      integer landsola(nlon)
      integer landsolb(nlon)
      real(r8) dx
      real(r8) dy
      real(r8) lat_i(nlat)
      real(r8) lon_i(nlon) 
      integer ia
      integer i, j, ii, jj, l
      integer nrow_start, nrow_end, ncol_start, ncol_end
      integer nrow, ncol, length
! ----------------------------------------------------------------------
! .. (1) elevation:
      inquire(iolength=length) inter1_chr
!      print*, iunit(1), lndname(1)
      open(iunit(1),file=lndname(1),access='direct',recl=length,& 
                    form='unformatted',status='old') 
!      print 15, iunit(1),lndname(1),length
! .. (2) land-water mask (not used);
! .. (3) vegetation;
! .. (4) soil 0 - 30 cm;
! .. (5) soil 30 - 100 cm;
      inquire(iolength=length) landtype_chr
      do i = 3, 5
        open(iunit(i),file=lndname(i),access='direct',recl=length,&
                      form='unformatted',status='old') 
      enddo
!        print 15, iunit(i),lndname(i),length
!      enddo
!15    format('==> open direct-access: fort.',i2.2,'  file=',&
!                  a50,' record length=',i6)
!   ---------------------------------------------------------------
!   buffer in lanuse and terrain data
!   ---------------------------------------------------------------
      dx = 1./120.
      dy = 1./120.
      do j = 1, nlat
         lat_i(j) = 90. - (2*j-1)*dy/2.
      enddo
      do i = 1, nlon
	 lon_i(i) = -180. + (2*i-1)*dx/2.
      enddo
      numtype(:,:,:) = 0
      numsola(:,:,:) = 0
      numsolb(:,:,:) = 0
      nrow_start = nint((90.-edgen)*120.) + 1
      nrow_end=nrow_start
!      nrow_end   = nint((90.-edges)*120.)
!      nrow_start=nrow_end
      if(nrow_start.lt.1) nrow_start = 1
      if(nrow_end.gt.nlat) nrow_end = nlat
!      if(lat_i(nrow_start).le.latn(1) .and. lat_i(nrow_start).gt.latn(2))then
!         print*, ' nrow_start point is between latn(1) and latn(2) '
!      else
!	 call abort
!      endif
!      if(lat_i(nrow_end).le.latn(jx) .and. lat_i(nrow_end).gt.latn(jx+1))then
!         print*, ' nrow_end point is between latn(jx) and latn(jx+1) '
!      else
!         call abort
!      endif
      ncol_start = nint((180.+edgew)*120.) + 1
      ncol_end   = ncol_start
!      ncol_end   = nint((180.+edgee)*120.)
      if(ncol_start.lt.1) ncol_start = 1
      if(ncol_end.gt.nlon) ncol_end = nlon
!     if(lon_i(ncol_start).ge.lonw(1) .and. lon_i(ncol_start).lt.lonw(2))then
!        print*, '  ncol_start point is between lonw(1) and lonw(2) '
!      else
!         call abort
!      endif
!      if(lon_i(ncol_end).ge.lonw(ix) .and. lon_i(ncol_end).lt.lonw(ix+1))then
!         print*, ' ncol_end point is between lonw(ix) and lonw(ix+1) '
!      else
!         call abort
!     endif
!      print 16 , nrow_start, nrow_end
!16    format(3x,'start_record=',i6,'  end-record=',i6)
      do nrow = nrow_start,nrow_end 
!         print*,'rdlanddata.F90-test-nrow',nrow,'nrow_start',nrow_start,&
!                'nrow_end',nrow_end
         read(iunit(1),rec=nrow,err=100) inter1_chr
!* read(iunit(2),rec=nrow,err=100) landmask_chr
         read(iunit(3),rec=nrow,err=100) landtype_chr
         read(iunit(4),rec=nrow,err=100) landsola_chr
         read(iunit(5),rec=nrow,err=100) landsolb_chr
! determine model grid - number of latitude
         do j = 1, jx
            if(lat_i(nrow).le.latn(j) .and. lat_i(nrow).gt.latn(j+1))then 
               jj = j
               exit
            endif
         enddo
         do ncol = ncol_start, ncol_end
!            print*,'rdlanddata.F90-test-ncol',ncol,'ncol_start',ncol_start,&
!                        'ncol_end',ncol_end
            inter1(ncol) = ia(inter1_chr(ncol),2,-9999)
!*landmask(ncol) = ichar(landmask_chr(ncol))
            landtype(ncol) = ichar(landtype_chr(ncol))
            landsola(ncol) = ichar(landsola_chr(ncol))                   
            landsolb(ncol) = ichar(landsolb_chr(ncol))
! convert ocean ID to 0 (Ross ice shelf north)
            if((inter1(ncol).eq.-9999) .and.&
                            (float(nrow)/120..lt.175.))then
!*    landmask(ncol) = 0
                landtype(ncol) = 0
                landsola(ncol) = 14
                landsolb(ncol) = 14
            endif
! determine model grid - number of longitude
            do i = 1, ix
               if(lon_i(ncol).ge.lonw(i) .and. lon_i(ncol).lt.lonw(i+1))then
                  ii = i
                  exit
               endif
            enddo
! mapping land cover type to model grid
!            L = landtype(ncol) + 1
            L=15
!            print*,'rdlanddata.F90-test-landcover-L',L
            if(L.gt.25) stop 'land categ gt 25'
            numtype(ii,jj,L) = numtype(ii,jj,L)+1 
! mapping upper level soil texture type (0 - 30 cm)
!            L = landsola(ncol)
            L=4
!            print*,'rdlandata.F90-test:landsola-L',L
            if(L.gt.17) stop 'surface soil categ gt 17'
            numsola(ii,jj,L) = numsola(ii,jj,L)+1 
! mapping bottom layer soil texture (30 - 100 cm)
            L = landsolb(ncol)
!            L=2
!            print*,'rdlandata.F90-test:landsolb-L',L
            if(L.gt.17) stop 'bottom soil categ gt 17'
            numsolb(ii,jj,L) = numsolb(ii,jj,L)+1 
         enddo
      enddo
      goto 400
100   print 101,nrow,iunit(i)
101   format(' record =',i8,',  error occured on unit: ',i3)
400   continue
  end subroutine rdlanddata
!-----------------------------------------------------------------------
!EOP
  integer function ia(chr,n,ispval)
!  purpose: to convert a n-bytes character (chr) to integer ia.
!        ** the integer data file is saved as a n-byte character
!           data file. this function is used to recover the
!           character data to the integer data.
!
!  n      --- the number of bytes in chr
!  ispval --- default value for the negative integer.
      character*(*) chr
      integer bit_1, bit_2
      bit_1 = '200'O     ! BINARY '10000000'
      bit_2 = '377'O     ! BINARY '11111111'
      ia    = 0
      ii1 = ichar(chr(1:1))
! .. get the sign -- isn=0 positive, isn=1 negative:
      jj  = iand(ii1,bit_1)
      isn = ishft(jj,-7)
! .. for negative number:
!    because the negative integers are represented by the supplementary
!    binary code inside machine.
        if (isn.eq.1) then
          do m = n+1,4
             nbit = (m-1)*8
             jj = ishft(bit_2,nbit)
             ia = ieor(jj,ia)
          end do
        endif
!   .. get the byte from chr:
         do m = 1,n
           ii2 = ichar(chr(m:m))
! new IBM xlf 8.1 compiler fix: thanks to Jim Edwards
           if (ii2.lt.0) ii2 = ii2 + 256
           mshft = (n-m)*8
           ia2   = ishft(ii2,mshft)
!   .. the abs(integer):
           ia = ieor(ia,ia2)
         end do
      if (ia.lt.0) ia = ispval
      return
  end function ia
!-----------------------------------------------------------------------
!EOP
