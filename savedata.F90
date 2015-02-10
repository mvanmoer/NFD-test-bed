! 20130117  
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
#ifdef __GFORTRAN__
        open(unit=iounit,access='stream',file=outfile,status='replace')
#else
        open(unit=iounit,access='sequential',form='binary',file=outfile,status='replace')
#endif
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
        outfile = name//"."//trim(adjustl(int2str))//".vec.raw"

        print*,'Writing unformated, interleaved 3-component file: '//outfile       
#ifdef __GFORTRAN__
        open(unit=iounit,access='stream',file=outfile,status='replace')
#else
        open(unit=iounit,access='sequential',form='binary',file=outfile,status='replace')
#endif

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
