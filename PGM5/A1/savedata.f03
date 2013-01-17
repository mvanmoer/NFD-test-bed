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
! mvm 20130117 -- changing from B. Jewett's
!  -- moving to Fortran 2003 in order to use streaming access
!  -- writing each timestep as a separate file, instead of one file containing
!     all variables and all timesteps.
!  -- changing to module
!  -- changed subroutine name
!  -- changed to take simulation step instead of simulation time
module savedata
    implicit none
    public
contains
    subroutine savefield(name,step,field,nx,ny,nz)

    integer, intent(in)                     :: nx, ny, nz
    character(len=*), intent(in)            :: name
    integer, intent(in)                     :: step 
    real,dimension (nx,ny,nz), intent(in)   :: field

    integer,parameter :: iounit=11
    character(len=256) :: outfile
    
    ! awaiting gfortran 4.7.0 for this
    ! character(:), allocatable :: outfile

    outfile = name//"."//char(step)//".raw"
    print*,'Writing unformatted data file: '//outfile
    open(unit=iounit,access='stream',file=outfile,status='new')
    write(iounit) field   
    close(iounit) 

    end subroutine savefield
end module savedata
