! diffusion.f90
module diffusion
  use mesh_type
  use fson
  implicit none
  public
contains
  subroutine diff(t, u3, v3, w3, u1, v1, w1, m, CTdt, fsonFile)
    type (mesh), intent(in) :: m
    real, dimension(-1:m%nx+2,-1:m%ny+2,-1:m%nz+2), intent(inout) :: t
    real, dimension(-1:m%nx+2,-1:m%ny+2,-1:m%nz+2) :: tprime
    real, dimension(0:m%nx+1,0:m%ny+1,0:m%nz+1), intent(inout) :: u3,v3,w3,u1,v1,w1
    real, intent(in) :: CTdt
    character(len=*), intent(in) :: fsonFile

    real :: K_t
    real :: K_m

    real :: tBar = 300.0 
    ! Redundant to cut out loop dependencies
    real :: del_xx = 0.0, del_yy = 0.0, del_zz = 0.0
    real :: del_xxu = 0.0, del_yyu = 0.0, del_zzu = 0.0
    real :: del_xxv = 0.0, del_yyv = 0.0, del_zzv = 0.0
    real :: del_xxw = 0.0, del_yyw = 0.0, del_zzw = 0.0
    integer :: i, j, k
    type(fson_value), pointer :: config
  
    config => fson_parse(fsonFile)
    call fson_get(config, "K_t", K_t)
    call fson_get(config, "K_m", K_m)
    call fson_destroy(config)
 
    tprime = t - tBar

    ! Theta
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

    ! U
    do k = 1, m%nz
       do j = 1, m%ny
          do i = 2, m%nx
             del_xxu = (u1(i+1,j,k) - 2.0 * u1(i,j,k) + u1(i-1,j,k)) / m%dx**2
             del_yyu = (u1(i,j+1,k) - 2.0 * u1(i,j,k) + u1(i,j-1,k)) / m%dy**2
             del_zzu = (u1(i,j,k+1) - 2.0 * u1(i,j,k) + u1(i,j,k-1)) / m%dz**2
             u3(i,j,k) = u3(i,j,k) + CTdt * K_m * (del_xx + del_yy + del_zz)
          end do
       end do
    end do

    ! V
    do k = 1, m%nz
       do j = 1, m%ny
          do i = 1, m%nx
             del_xxv = (v1(i+1,j,k) - 2.0 * v1(i,j,k) + v1(i-1,j,k)) / m%dx**2
             del_yyv = (v1(i,j+1,k) - 2.0 * v1(i,j,k) + v1(i,j-1,k)) / m%dy**2
             del_zzv = (v1(i,j,k+1) - 2.0 * v1(i,j,k) + v1(i,j,k-1)) / m%dz**2
             v3(i,j,k) = v3(i,j,k) + CTdt * K_m * (del_xx + del_yy + del_zz)
          end do
       end do
    end do

    ! W
    do k = 2, m%nz
       do j = 1, m%ny
          do i = 1, m%nx
             del_xxw = (w1(i+1,j,k) - 2.0 * w1(i,j,k) + w1(i-1,j,k)) / m%dx**2
             del_yyw = (w1(i,j+1,k) - 2.0 * w1(i,j,k) + w1(i,j-1,k)) / m%dy**2
             del_zzw = (w1(i,j,k+1) - 2.0 * w1(i,j,k) + w1(i,j,k-1)) / m%dz**2
             w3(i,j,k) = w3(i,j,k) + CTdt * K_m * (del_xx + del_yy + del_zz)
          end do
       end do
    end do
   end subroutine diff
end module diffusion
