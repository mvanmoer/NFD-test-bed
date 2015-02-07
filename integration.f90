! integration.f90
! Preps 3D arrays for calls to 1D advection routine
! using Strang splitting.
module integration
  use mesh_type
  implicit none
  public
contains
  subroutine integrate(t1, t2, u1, v1, w1, u2, v2, w2, u3, v3, w3, & 
       m, dt, CTdt)
    ! t1 -- physical domain scalars at current time step
    ! t2 -- physical domain scalars at next time step
    ! u  -- u-component of wind field
    ! v  -- v-component of wind field
    ! w  -- w-component of wind field
    ! dx, dy, dz -- grid spacing
    ! dt -- time delta   
    ! nx, ny, nz -- dims
    use advection
    type (mesh), intent(in) :: m
    real, dimension(-1:m%nx+2,-1:m%ny+2,-1:m%nz+2), intent(inout) :: t1, t2
    real, dimension(0:m%nx+1,0:m%ny+1,0:m%nz+1), intent(in) :: u1, v1, w1, u2, v2, w2
    real, dimension(0:m%nx+1,0:m%ny+1,0:m%nz+1), intent(inout) :: u3, v3, w3
    real, intent(in) :: dt, CTdt
    integer :: i, j, k
    real :: uu_x, vu_y, wu_z, uv_x, vv_y, wv_z, uw_x, vw_y, ww_z

    uu_x = 0.0
    vu_y = 0.0
    wu_z = 0.0
    uv_x = 0.0
    vv_y = 0.0
    wv_z = 0.0
    uw_x = 0.0
    vw_y = 0.0
    ww_z = 0.0

    call uadvect()
    call vadvect()
    call wadvect()
    call strang()

  contains
    subroutine uadvect()
      ! U advection  
      !$omp parallel do private (i, j, k)
      do k = 1, m%nz
         do j = 1, m%ny
            do i = 2, m%nx
               uu_x = (u2(i+1,j,k)**2 - u2(i-1,j,k)**2)/(4.0*m%dx)
       
               vu_y = ((v2(i,j+1,k) + v2(i-1,j+1,k)) * &
                    (u2(i,j+1,k) - u2(i,j,k)) + &
                    (v2(i,j,k) + v2(i-1,j,k)) * &
                    (u2(i,j,k) - u2(i,j-1,k))) / (4.0 * m%dy)
               
               wu_z = ((w2(i,j,k+1) + w2(i-1,j,k+1)) * &
                    (u2(i,j,k+1) - u2(i,j,k)) + &
                    (w2(i,j,k) + w2(i-1,j,k)) * &
                    (u2(i,j,k) - u2(i,j,k-1))) / (4.0 * m%dz)
               
               u3(i,j,k) = u3(i,j,k) - CTdt * (uu_x + vu_y + wu_z)
            end do
         end do
      end do
      !$omp end parallel do
    end subroutine uadvect

    subroutine vadvect()
      ! V advection
      !$omp parallel do private (i, j, k)
      do k = 1, m%nz
         do j = 1, m%ny
            do i = 1, m%nx
               uv_x = ((u2(i+1,j,k) + u2(i+1,j-1,k)) * &
                    (v2(i+1,j,k) - v2(i,j,k)) + &
                       (u2(i,j,k) + u2(i,j-1,k)) * &
                       (v2(i,j,k) - v2(i-1,j,k))) / (4.0 * m%dx)
              
               vv_y = (v2(i,j+1,k)**2 - v2(i,j-1,k)**2)/(4.0*m%dy)
               
               wv_z = ((w2(i,j,k+1) + w2(i,j-1,k+1)) * &
                    (v2(i,j,k+1) - v2(i,j,k)) + &
                    (w2(i,j,k) + w2(i,j-1,k)) * &
                    (v2(i,j,k) - v2(i,j,k-1))) / (4.0 * m%dz)

               v3(i,j,k) = v3(i,j,k) - CTdt * (uv_x + vv_y + wv_z)
            end do
         end do
      end do
      !$omp end parallel do
    end subroutine vadvect

    subroutine wadvect
      ! W advection
      !$omp parallel do private (i, j, k)
      do k = 2, m%nz
         do j = 1, m%ny
            do i = 1, m%nx
               uw_x = ((u2(i+1,j,k) + u2(i+1,j,k-1)) * &
                    (w2(i+1,j,k) - w2(i,j,k)) + &
                    (u2(i,j,k) + u2(i,j,k-1)) * &
                    (w2(i,j,k) - w2(i-1,j,k))) / (4.0 * m%dx)
              
               vw_y = ((v2(i,j+1,k) + v2(i,j+1,k-1)) * &
                    (w2(i,j+1,k) - w2(i,j,k)) + &
                    (v2(i,j,k) + v2(i,j,k-1)) * &
                    (w2(i,j,k) - w2(i,j-1,k))) / (4.0 * m%dy)
               
               ww_z = (w2(i,j,k+1)**2 - w2(i,j,k-1)**2)/(4.0*m%dz)
               
               w3(i,j,k) = w3(i,j,k) - CTdt * (uw_x + vw_y + ww_z)
            end do
         end do
      end do
      !$omp end parallel do

    end subroutine wadvect


    subroutine strang()
      ! Strang splitting - 
      ! Fx(dt/2),Fy(dt/2),Fz(dt),Fy(dt/2),Fx(dt/2)

      ! do x dir
      do j = 1, m%ny
         do k = 1, m%nz
            call advect1D(t1(:,j,k), t2(:,j,k), u1(:,j,k), m%dx, dt/2.0)
         end do
      end do

      ! Do y dir
      do i = 1, m%nx
         do k = 1, m%nz
            call advect1D(t1(i,:,k),t2(i,:,k), v1(i,:,k), m%dy, dt/2.0)
         end do
      end do

      ! do z dir
      do i = 1, m%nx
         do j = 1, m%ny
            call advect1D(t1(i,j,:), t2(i,j,:), w1(i,j,:), m%dz, dt)
         end do
      end do

      ! do y dir
      do i = 1, m%nx
         do k = 1, m%nz
            call advect1D(t1(i,:,k),t2(i,:,k), v1(i,:,k), m%dy, dt/2.0)
         end do
      end do

      ! do x dir
      do j = 1, m%ny
         do k = 1, m%nz
            call advect1D(t1(:,j,k),t2(:,j,k), u1(:,j,k), m%dx, dt/2.0)
         end do
      end do
      return
    end subroutine strang
  end subroutine integrate
end module integration
