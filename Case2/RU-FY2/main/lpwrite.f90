  subroutine lpwrite(idate_p,idate,lhistTimeVar,lhistTimeVar2,luout,luout2,foutdat,lwrite)
  implicit none
  integer, intent(in) :: idate_p(3)
  integer, intent(in) :: idate(3)
  integer, intent(in) :: lhistTimeVar
  integer, intent(in) :: lhistTimeVar2
  integer, intent(in) :: luout
  integer, intent(in) :: luout2
  character(LEN=256), intent(in) :: foutdat 
  logical, intent(inout) :: lwrite
  integer :: months(12)
  integer :: j_month_e
  logical :: leapyear
  integer i, l
  character(LEN=256) fhistTimeVar_name
  character(LEN=256) fhistTimeVar_name2
  character(LEN=256) fout
  character(LEN=256) fout2
  character(LEN=256) cdate
     lwrite = .true.
! Open for model time varying data (model state variables) and history filed
  if(lwrite)then
     write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') idate_p(1),idate_p(2),idate_p(3)
     fhistTimeVar_name = trim(foutdat)//'-rstTimeVar'//'-'//trim(cdate)
     fhistTimeVar_name2= '/home/u120220909911/RU-FY2/Third/psuade/NEE/LH/RU-FY2/output/rst/Valdai'//&
                           '-rstTimeVar'//'-'//trim(cdate)//'.txt'
     fout=trim(foutdat)//'-'//trim(cdate)
     fout2=trim(foutdat)//'-'//trim(cdate)//'.txt'
     open(unit=lhistTimeVar,file=fhistTimeVar_name,form='unformatted',&
                            status='new',action='write')
     open(unit=lhistTimeVar2,file=fhistTimeVar_name2,form='formatted',&
                            status='new',action='write')
     open(unit=luout,file=fout,access='sequential',form='unformatted',&
                     status='new',action='write')
     open(unit=luout2,file=fout2,access='sequential',form='formatted',&
                     status='new',action='write')
  endif
  end subroutine lpwrite
! ------------------------------------------------------------------------
! EOP
