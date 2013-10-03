! 20130916 - file to read various simulation parameters from file
! to avoid recompiling for ever parameter change.
module read_parameters_module 
    use mesh_type
    implicit none
    public
contains
    subroutine read_parameters(infile,m,dt,nstep,nplot)
        type (mesh), intent(inout)              :: m
        character(len=*), intent(in)            :: infile 
        integer, intent(inout)                  :: nstep, nplot
        real, intent(inout)                     :: dt

        integer,parameter :: iounit=11
 
        print *,'Loading parameters from: '//infile
        open(unit=iounit,file=infile,status='old')
        rewind(iounit)
        read(iounit,*) m%nx, m%ny, m%nz, m%dx, m%dy, m%dz, dt, nstep, nplot
        close(iounit) 

    end subroutine read_parameters 
end module read_parameters_module 
