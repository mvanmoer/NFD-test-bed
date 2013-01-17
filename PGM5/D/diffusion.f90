! diffusion.f90
! ATMS502 Fall 2011 Program 5
! Mark Van Moer

module diffusion
  implicit none
  public
contains
  subroutine diff(t, u3, v3, w3, u1, v1, w1, dx, dy, dz, nx, ny, nz, CTdt)
    integer, intent(in) :: nx, ny, nz
    real, intent(in) :: dx, dy, dz
    real, dimension(-1:nx+2,-1:ny+2,-1:nz+2), intent(inout) :: t
    real, dimension(-1:nx+2,-1:ny+2,-1:nx+2) :: tprime
    real, dimension(0:nx+1,0:ny+1,0:nz+1), intent(inout) :: u3,v3,w3,u1,v1,w1
    real, intent(in) :: CTdt

    ! low res
    real :: K_t = 5.0
    real :: K_m = 50.0

    real :: del_xx = 0.0, del_yy = 0.0, del_zz = 0.0, tBar = 300.0
    integer :: i, j, k
    
    tprime = t - tBar

!$OMP   PARALLEL DO PRIVATE (i,j,k)
    ! Theta
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             del_xx = (t(i+1,j,k) - 2.0 * t(i,j,k) + t(i-1,j,k)) / dx**2
             del_yy = (t(i,j+1,k) - 2.0 * t(i,j,k) + t(i,j-1,k)) / dy**2
             del_zz = (tprime(i,j,k+1) - 2.0 * tprime(i,j,k) &
                  + tprime(i,j,k-1)) / dz**2
             t(i,j,k) = t(i,j,k) + K_t * (del_xx + del_yy + del_zz)
          end do
       end do
    end do
!$OMP   END PARALLEL DO
    del_xx = 0.0
    del_yy = 0.0
    del_zz = 0.0

!$OMP PARALLEL DO PRIVATE (i,j,k)    
    ! U
    do k = 1, nz
       do j = 1, ny
          do i = 2, nx
             del_xx = (u1(i+1,j,k) - 2.0 * u1(i,j,k) + u1(i-1,j,k)) / dx**2
             del_yy = (u1(i,j+1,k) - 2.0 * u1(i,j,k) + u1(i,j-1,k)) / dy**2
             del_zz = (u1(i,j,k+1) - 2.0 * u1(i,j,k) + u1(i,j,k-1)) / dz**2
             u3(i,j,k) = u3(i,j,k) + CTdt * K_m * (del_xx + del_yy + del_zz)
          end do
       end do
    end do
!$OMP END PARALLEL DO
    del_xx = 0.0
    del_yy = 0.0
    del_zz = 0.0

!$OMP PARALLEL DO PRIVATE (i,j,k)
    ! V
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             del_xx = (v1(i+1,j,k) - 2.0 * v1(i,j,k) + v1(i-1,j,k)) / dx**2
             del_yy = (v1(i,j+1,k) - 2.0 * v1(i,j,k) + v1(i,j-1,k)) / dy**2
             del_zz = (v1(i,j,k+1) - 2.0 * v1(i,j,k) + v1(i,j,k-1)) / dz**2
             v3(i,j,k) = v3(i,j,k) + CTdt * K_m * (del_xx + del_yy + del_zz)
          end do
       end do
    end do

!$OMP END PARALLEL DO
    del_xx = 0.0
    del_yy = 0.0
    del_zz = 0.0

!$OMP PARALLEL DO PRIVATE (i,j,k)

    ! W
    do k = 2, nz
       do j = 1, ny
          do i = 1, nx
             del_xx = (w1(i+1,j,k) - 2.0 * w1(i,j,k) + w1(i-1,j,k)) / dx**2
             del_yy = (w1(i,j+1,k) - 2.0 * w1(i,j,k) + w1(i,j-1,k)) / dy**2
             del_zz = (w1(i,j,k+1) - 2.0 * w1(i,j,k) + w1(i,j,k-1)) / dz**2
             w3(i,j,k) = w3(i,j,k) + CTdt * K_m * (del_xx + del_yy + del_zz)
          end do
       end do
    end do
!$OMP END PARALLEL DO

   end subroutine diff
end module diffusion
