! diffusion.f90
! ATMS502 Fall 2011 Program 5
! Mark Van Moer
!
! Test E
!
module diffusion
  use mesh_type
  use fson
  implicit none
  public
contains
  subroutine diff(t, u3, v3, w3, u1, v1, w1, m, CTdt, fsonFile)
    type (mesh), intent(in) :: m
    real, dimension(-1:m%nx+2,-1:m%ny+2,-1:m%nz+2), intent(inout) :: t
    real, dimension(-1:m%nx+2,-1:m%ny+2,-1:m%nx+2) :: tprime
    real, dimension(0:m%nx+1,0:m%ny+1,0:m%nz+1), intent(inout) :: u3,v3,w3,u1,v1,w1
    real, intent(in) :: CTdt
    character(len=*), intent(in) :: fsonFile

    real :: K_t
    real :: K_m

    real :: del_xx = 0.0, del_yy = 0.0, del_zz = 0.0, tBar = 300.0
    integer :: i, j, k
    type(fson_value), pointer :: config
  
    config => fson_parse(fsonFile)
    call fson_get(config, "K_t", K_t)
    call fson_get(config, "K_m", K_m)
    call fson_destroy(config)
 
    tprime = t - tBar

    ! Theta
!$omp parallel do private (i,j,k)
    do k = 1, m%nz
       do j = 1, m%ny
          do i = 1, m%nx
             del_xx = (t(i+1,j,k) - 2.0 * t(i,j,k) + t(i-1,j,k)) / m%dx**2
             del_yy = (t(i,j+1,k) - 2.0 * t(i,j,k) + t(i,j-1,k)) / m%dy**2
             del_zz = (tprime(i,j,k+1) - 2.0 * tprime(i,j,k) &
                  + tprime(i,j,k-1)) / m%dz**2
             t(i,j,k) = t(i,j,k) + K_t * (del_xx + del_yy + del_zz)
          end do
       end do
    end do
!$omp end parallel do
    del_xx = 0.0
    del_yy = 0.0
    del_zz = 0.0

    ! U
!$omp parallel do private (i,j,k)
    do k = 1, m%nz
       do j = 1, m%ny
          do i = 2, m%nx
             del_xx = (u1(i+1,j,k) - 2.0 * u1(i,j,k) + u1(i-1,j,k)) / m%dx**2
             del_yy = (u1(i,j+1,k) - 2.0 * u1(i,j,k) + u1(i,j-1,k)) / m%dy**2
             del_zz = (u1(i,j,k+1) - 2.0 * u1(i,j,k) + u1(i,j,k-1)) / m%dz**2
             u3(i,j,k) = u3(i,j,k) + CTdt * K_m * (del_xx + del_yy + del_zz)
          end do
       end do
    end do
!$omp end parallel do
    del_xx = 0.0
    del_yy = 0.0
    del_zz = 0.0

    ! V
!$omp parallel do private (i,j,k)
    do k = 1, m%nz
       do j = 1, m%ny
          do i = 1, m%nx
             del_xx = (v1(i+1,j,k) - 2.0 * v1(i,j,k) + v1(i-1,j,k)) / m%dx**2
             del_yy = (v1(i,j+1,k) - 2.0 * v1(i,j,k) + v1(i,j-1,k)) / m%dy**2
             del_zz = (v1(i,j,k+1) - 2.0 * v1(i,j,k) + v1(i,j,k-1)) / m%dz**2
             v3(i,j,k) = v3(i,j,k) + CTdt * K_m * (del_xx + del_yy + del_zz)
          end do
       end do
    end do
!$omp end parallel do
    del_xx = 0.0
    del_yy = 0.0
    del_zz = 0.0


    ! W
!$omp parallel do private (i,j,k)
    do k = 2, m%nz
       do j = 1, m%ny
          do i = 1, m%nx
             del_xx = (w1(i+1,j,k) - 2.0 * w1(i,j,k) + w1(i-1,j,k)) / m%dx**2
             del_yy = (w1(i,j+1,k) - 2.0 * w1(i,j,k) + w1(i,j-1,k)) / m%dy**2
             del_zz = (w1(i,j,k+1) - 2.0 * w1(i,j,k) + w1(i,j,k-1)) / m%dz**2
             w3(i,j,k) = w3(i,j,k) + CTdt * K_m * (del_xx + del_yy + del_zz)
          end do
       end do
    end do
!$omp end parallel do
   end subroutine diff
end module diffusion
