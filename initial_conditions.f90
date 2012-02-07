!  initial_conditions.f90
!  Mark Van Moer, ATMS502 Fall 2011
!  Sets the initial conditions.
!  First set the "base state vertical profiles"

module initial_conditions
  implicit none
  public    
contains  
  subroutine ic(t, u, v, rho_t, rho_w, nx, ny, nz, dx, dy, dz)

    ! t -- "theta", 3D array of scalars
    ! u -- 3D array of momentum u-component
    ! v -- 3D array of momentum v-component
    ! w -- 3D array of momentum w-component
    ! p -- 3D array of pressure
    ! rho_t -- 1D array of density at Theta, u, pprime levels
    ! rho_w -- 1D array of density at w levels

    real, dimension(:,:,:), intent(inout) :: t, u, v
    real, dimension(:), intent(inout) :: rho_t, rho_w
    real :: pi
    real, intent(in) :: dx, dy, dz
    integer, intent(in) :: nx, ny, nz
 
    ! Final
    real, parameter :: x1_0 =   25.0, y1_0 = 7775.0, z1_0 = 1500.0
    real, parameter :: x2_0 = 9525.0, y2_0 = 4525.0, z2_0 = 4000.0
    real, parameter :: xr1 = 3500.0, yr1 = 999999.0, zr1 = 3500.0
    real, parameter :: xr2 = 2500.0, yr2 = 2500.0, zr2 = 2500.0
    real, parameter :: dtPrime1 = -10.0
    real, parameter :: dtPrime2 = -8.5
    real, parameter :: Thetabar = 300.0
    real, parameter :: deltaV1 = -30.0
    real, parameter :: deltaV2 = 30.0
    real, parameter :: g = 9.81
    real, parameter :: c_p = 1004.0
    real, parameter :: R_d = 287.0
    real, parameter :: P_0 = 100000.0
    real, parameter :: uperturb = 1.0
    real, dimension(nz) :: zed, Tbar, Pbar

    real :: ranval, rand
    integer :: i, j, k
    real :: r1, r2, x, y, z

    pi = 4.0 * atan(1.0)

    ! base state vertical profiles
    ! each z-slice is the same
   
    do k = 1, nz
       ! height in meters at Theta,u,pPrime levels
       zed(k) = dz / 2.0 + dz * real(k - 1)

       ! temperature in Kelvins 
       Tbar(k) = Thetabar - g/c_p * zed(k)

       ! pressure in Pascals
       Pbar(k) = P_0 * (Tbar(k)/Thetabar)**(c_p/R_d)

       ! density
       rho_t(k) = Pbar(k)/(R_d * Tbar(k))
    end do

    ! rho_w is staggered from rho_t
    do k = 2, nz
       rho_w(k) = 0.5 * (rho_t(k) + rho_t(k-1))
    end do
    rho_w(1) = 9999.0       ! flag -- unused
    rho_w(nz+1) = 9999.0    ! flag -- unused
       
    ! Loop over physical domain.
    t = Thetabar
    v = 0.0
    do k = 1, nz    
       do j = 1, ny
          do i = 1, nx

             ! physical coordinate generation
             x = real(i - 1) * dx + dx / 2.0
             y = real(j - 1) * dy + dy / 2.0
             z = real(k - 1) * dz + dz / 2.0
             
             ! radius only one perturbation for Tests A,B,C
             r1 = sqrt(((x - x1_0)/xr1)**2 + ((y - y1_0)/yr1)**2 + &
                  ((z - z1_0)/zr1)**2)
             r2 = sqrt(((x - x2_0)/xr2)**2 + ((y - y2_0)/yr2)**2 + &
                  ((z - z2_0)/zr2)**2)
                          
             
             if (r1 <= 1.0) then 
                t(i,j,k) = t(i,j,k) + (dtPrime1 * (cos(r1 * pi) + 1.0)/2.0)
                v(i,j,k) = v(i,j,k) + (deltaV1 * (cos(r1 * pi) + 1.0)/2.0)
              else if (r2 <= 1.0) then
                 t(i,j,k) = t(i,j,k) + (dtPrime2 * (cos(r2 * pi) + 1.0)/2.0)
                 v(i,j,k) = v(i,j,k) + (deltaV2 * (cos(r2 * pi) + 1.0)/2.0)
             end if
             
          enddo
       enddo
    enddo


    ! Test F U perturbations
    call srand(0)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx + 1
             ranval = rand(0)
             u(i,j,k) = (ranval - 0.5) * uperturb
          end do
       end do
    end do
 
  end subroutine ic
end module initial_conditions
