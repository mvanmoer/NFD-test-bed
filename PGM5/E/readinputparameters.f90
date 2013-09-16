! 20130916 - file to read various simulation parameters from file
! to avoid recompiling for ever parameter change.
module read_parameters_module 
    implicit none
    public
contains
    subroutine read_parameters(infile,nx,ny,nz,dx,dy,dz,dt,nstep,nplot)

        character(len=*), intent(in)            :: infile 
        integer, intent(inout)                  :: nx, ny, nz, nstep, nplot
        real, intent(inout)                     :: dx, dy, dz, dt

        integer,parameter :: iounit=11
 
        print *,'Loading parameters from: '//infile
        open(unit=iounit,file=infile,status='old')
        rewind(iounit)
        read(iounit,*) nx, ny, nz, dx, dy, dz, dt, nstep, nplot
        close(iounit) 

    end subroutine read_parameters 
end module read_parameters_module 
