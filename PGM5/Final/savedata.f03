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
!  -- added interface to handle writing single scalar OR combining separate
!     scalars as interleaved vectorcomponents
module savedata
    use mesh_type
    implicit none
    public
    interface savefield
        module procedure saveScalar
        module procedure save3DVectorComponents
    end interface savefield
contains
    subroutine saveScalar(name,step,field,m)

        type (mesh), intent(in) :: m
        character(len=*), intent(in)            :: name
        integer, intent(in)                     :: step 
        real,dimension (m%nx,m%ny,m%nz), intent(in)   :: field

        integer,parameter :: iounit=11
        character(len=256) :: outfile
        character(len=256) :: int2str

        ! awaiting gfortran 4.7.0 for this
        ! character(:), allocatable :: outfile
    
        ! sigh...
        write(int2str,*) step
        outfile = name//"."//trim(adjustl(int2str))//".raw"
 
        print*,'Writing unformatted scalar file: '//outfile
        open(unit=iounit,access='stream',file=outfile,status='replace')
        write(iounit) field   
        close(iounit) 

    end subroutine saveScalar
    
    subroutine save3DVectorComponents(name,step,x,y,z,m)
        type (mesh), intent(in) :: m 
        character(len=*), intent(in) :: name
        integer, intent(in) :: step
        real,dimension(m%nx,m%ny,m%nz), intent(in) :: x, y, z

        integer, parameter :: iounit=11
        character(len=256) :: outfile
        character(len=256) :: int2str
        integer :: i, j, k

        write(int2str,*) step
        outfile = name//"."//trim(adjustl(int2str))//".vec"

        print*,'Writing unformated, interleaved 3-component file: '//outfile       
        open(unit=iounit,access='stream',file=outfile,status='replace')

        do k = 1, m%nz
            do j = 1, m%ny
                do i = 1, m%nx
                    write(iounit) x(i,j,k)
                    write(iounit) y(i,j,k)
                    write(iounit) z(i,j,k)
                end do
            end do
        end do

        close(iounit)       

    end subroutine save3DVectorComponents
end module savedata
