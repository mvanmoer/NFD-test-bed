! -------------------------------------------------------
! PUTFIELD - write a data field to the run history file.
! -------------------------------------------------------
! Arguments:
!
!   name       input   char *    name of the field (only 4 chars used)
!   datatime   input   real      model integration time
!   field      input   real      array of data
!   nx,ny,nz   input   integers  dimensions of field()
!
! If name = "*", the output file is closed (nothing else done)
!
! ==NOTE== make sure you pass a character string, e.g. "NAME"
!          rather than just a character, even if that string
!          is 1-char long.  A string is assumed as input.
!
! A header integer written to the file is set to 0 if the
! file is created from Fortran, and 1 from C.
!
!
    subroutine putfield(name,datatime,field,nx,ny,nz)
    implicit none

    integer                   :: nx,ny,nz
    character*(*)                name
    real                      :: datatime
    real,dimension (nx,ny,nz) :: field

    integer,parameter :: iounit=11,sourcetype=0
    integer,save      :: count=0
    integer           :: ilen
    character         :: tname*4
!
! ... If name = '*' close file and return
!
    if (name.eq.'*') then
      print*,'Output history file closed.'
      close(iounit)
      count = -1
      return
    endif
    if (count.lt.0) then
      print*,'putfield: error: file has already been closed.'
      stop
    endif
!
! ... First call: open file, write file data type as Fortran.
!
    if (count.eq.0) then
      print*,'Writing unformatted data file RunHistory.dat'
      open(iounit,access='stream',file='RunHistory.raw',status='unknown',form='unformatted')
      rewind(iounit)
      !write(iounit) sourcetype
    endif
    count = count+1

!    write(tname,'(a)') name(1:min(4,len(name)))
!    write(6,10) count,tname,nx,ny,nz,datatime
!10	format(' Writing field ',i4,': ',a,'(',i3,',',i3,',',i3,') for T=',f6.1)

! ... Write header
!    write(iounit) count,tname,datatime,nx,ny,nz

! ... Write array
    write(iounit) field

    return
    end
