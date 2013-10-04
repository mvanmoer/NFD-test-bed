! Pressure Gradient Forces and Bouyancy
! Program 5
! ATMS 502 Fall 2011
! Mark Van Moer
! Computes the pgf and bouyancy components for u,v,w
module pgf_bouyancy
  use mesh_type
  implicit none
  public
contains
  subroutine pgf(u3, v3, w3, p3, p1, t2, rho_t, rho_w, &
       m, CTdt)
    type (mesh), intent(in) :: m
    real, dimension(0:m%nx+1,0:m%ny+1,0:m%nz+1), intent(inout) :: u3
    real, dimension(0:m%nx+1,0:m%ny+1,0:m%nz+1), intent(inout) :: v3
    real, dimension(0:m%nx+1,0:m%ny+1,0:m%nz+1), intent(inout) :: w3
    real, dimension(-1:m%nx+2,-1:m%ny+2,-1:m%nz+2), intent(inout) :: t2
    real, dimension(0:m%nx+1,0:m%ny+1,1:m%nz), intent(inout) :: p3, p1
    real, dimension(1:m%nz), intent(in) :: rho_t
    real, dimension(1:m%nz+1), intent(in) :: rho_w
    real, intent(in) :: CTdt
    real, dimension(1:m%nx,1:m%ny,1:m%nz) :: tprime
 
    real :: tBar = 300.0, g = 9.81
 
    ! low res
    real :: c_s
    integer :: i, j, k

    integer, parameter :: iounit = 11
    open(unit=iounit,file='pgf.input',status='old')
    rewind(iounit)
    read(iounit,*) c_s
    close (iounit)

    tprime = t2(1:m%nx,1:m%ny,1:m%nz) - tBar

! mvm - is it better to have 
! a single PARALLEL that some how goes over all three of these loops
! since they are independent ?
! or is it already smart enough to do that?
! that is, is there an implied sequence point after each directive?

    !$OMP PARALLEL DO PRIVATE (i,j,k)
    ! u3 terms
    do k = 1, m%nz
       do j = 1, m%ny
          do i = 2, m%nx 
             u3(i,j,k) = u3(i,j,k) - &
                  CTdt * (1.0/rho_t(k) * (p1(i,j,k) - p1(i-1,j,k))/m%dx)
          end do
       end do
    end do
    !$OMP END PARALLEL DO


    !$OMP PARALLEL DO PRIVATE (i,j,k)
    !       ! v3 terms
    do k = 1, m%nz
       do j = 1, m%ny 
          do i = 1, m%nx
             v3(i,j,k) = v3(i,j,k) - &
                  CTdt * (1.0/rho_t(k) * (p1(i,j,k) - p1(i,j-1,k))/m%dy)
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE (i,j,k)
    ! w3 terms
    do k = 2, m%nz
       do j = 1, m%ny
          do i = 1, m%nx
             w3(i,j,k) = w3(i,j,k) - &
                  (CTdt * 1.0/rho_w(k) * (p1(i,j,k) - p1(i,j,k-1))/m%dz) + &
                  (g * CTdt * ((tprime(i,j,k) + tprime(i,j,k-1))/(2.0 * tBar)))
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    ! set BCs here
    call uvw_bc()

    !$OMP PARALLEL DO PRIVATE (i,j,k)
    ! p3 terms
    do k = 1, m%nz
       do j = 1, m%ny
          do i = 1, m%nx
             p3(i,j,k) = p1(i,j,k) - (CTdt * c_s**2 * &
                  (rho_t(k) * ((u3(i+1,j,k) - u3(i,j,k))/m%dx) + &
                  rho_t(k) * ((v3(i,j+1,k) - v3(i,j,k))/m%dy) + &
                  (rho_w(k+1) * w3(i,j,k+1) - rho_w(k) * w3(i,j,k))/m%dz))                   
          end do
       end do
    end do
    !$OMP END PARALLEL DO            

  contains
    subroutine uvw_bc()
      ! reflect x
      u3(0,:,:) = -u3(3,:,:)
      u3(1,:,:) = -u3(2,:,:)
      u3(m%nx+1,:,:) = -u3(m%nx,:,:)
      ! periodic y
      u3(:,0,:) = u3(:,m%ny,:)
      u3(:,m%ny+1,:) = u3(:,1,:)
      ! zero gradient z
      u3(:,:,0) = u3(:,:,1)
      u3(:,:,m%nz+1) = u3(:,:,m%nz)

      ! reflect x
      v3(0,:,:) = v3(2,:,:)
      v3(m%nx+1,:,:) = v3(m%nx-1,:,:)
      ! periodic y
      v3(:,0,:) = v3(:,m%ny,:)
      v3(:,m%ny+1,:) = v3(:,1,:)
      ! zero gradient z
      v3(:,:,0) = v3(:,:,1)
      v3(:,:,m%nz+1) = v3(:,:,m%nz)

      ! reflect x
      w3(0,:,:) = w3(2,:,:)
      w3(m%nx+1,:,:) = w3(m%nx-1,:,:)
      ! periodic y
      w3(:,0,:) = w3(:,m%ny,:)
      w3(:,m%ny+1,:) = w3(:,1,:)
      ! rigid z
      w3(:,:,0) = 0.0
      w3(:,:,1) = 0.0
      w3(:,:,m%nz+1) = 0.0
    end subroutine uvw_bc
  end subroutine pgf
end module pgf_bouyancy
