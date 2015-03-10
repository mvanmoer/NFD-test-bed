! advection.f90 
! module contains subroutine for performing advection 
! along 1D array. Uses a Piecewise-Linear method.
module advection
  implicit none
  public
contains
  subroutine advect1D(s1, s2, vel, dx, dt)
    ! advect1D - performs advection along 1D array
    ! s1  -- (-1:nx+2) input vector of current timestep scalars
    ! s2  -- (-1:nx+2) output vector of next timestep scalars
    ! vel -- (0:nx+1)  input vector of staggered C-grid velocities
    ! dx  -- physical grid spacing
    ! dt  -- time delta
     
    real, dimension(-1:), intent(inout) :: s1
    real, dimension(-1:), intent(inout) :: s2
    real, dimension(0:), intent(in) :: vel
    real, intent(in) :: dx, dt
    real :: dtOverdx
    integer :: i
    integer :: length
   
    length = size(s1) - 4 ! num ghost cells on s1/s2
    dtOverdx = dt / dx
  
    call piecewise_linear()
   
    return
  contains
    subroutine piecewise_linear()
      real :: r1, r2, dsi, dsminus1, dsplus1, F1, F2, a, b, c, d

      do i = 1, length
         r1 = abs(dtOverdx * vel(i))
         r2 = abs(dtOverdx * vel(i + 1))
         dsi = (s1(i + 1) - s1(i - 1)) * 0.5
         dsminus1 = (s1(i) - s1(i - 2)) * 0.5
         dsplus1 = (s1(i + 2) - s1(i)) * 0.5

         if (vel(i) >= 0) then
            F1 = r1 * (s1(i - 1) + (1 - r1) * 0.5 * dsminus1)
         else
            F1 = r1 * (-s1(i) + (1 - r1) * 0.5 * dsi)
         end if
         
         if (vel(i + 1) >= 0) then
            F2 = r2 * (s1(i) + (1 - r2) * 0.5 * dsi)
         else
            F2 = r2 * (-s1(i + 1) + (1 - r2) * 0.5 * dsplus1)
         end if
         
         s2(i) = s1(i) - (F2 - F1) + dtOverdx * s1(i) * (vel(i + 1) - vel(i))
         
      end do

      ! extremely important to NOT copy ghost points!!!!
      s1(1:length) = s2(1:length)
    end subroutine piecewise_linear
     
  end subroutine advect1D
end module advection
