!  initial_conditions.f90
!  Mark Van Moer, ATMS502 Fall 2011
!  Sets the initial conditions.
!  First set the "base state vertical profiles"

module initial_conditions
  use mesh_type
  implicit none
  public    
  type thermal
    real :: x, y, z
    real :: rx, ry, rz
    real :: dtPrime
    real :: dU, dV, dW
  end type thermal
contains  
  subroutine ic(t, u, v, rho_t, rho_w, m)

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
    type (mesh), intent(in) :: m
    integer :: numthermals
 
    real, parameter :: Thetabar = 300.0
    real, parameter :: g = 9.81
    real, parameter :: c_p = 1004.0
    real, parameter :: R_d = 287.0
    real, parameter :: P_0 = 100000.0
    real, parameter :: uperturb = 1.0
    real, dimension(m%nz) :: zed, Tbar, Pbar
    real, parameter :: flag = 9999.0

    real :: ranval, rand
    integer :: i, j, k, ii
    real :: d, x, y, z

    type (thermal), allocatable, dimension(:) :: therms
    integer, parameter :: iounit = 11

    ! read thermal data from input.. feel like this should be wrapped, somehow
    open(unit=iounit,file='thermals.input',status='old')
    rewind(iounit)
    read(iounit,*) numthermals
    allocate(therms(numthermals))    

    do i = 1, numthermals
      read(iounit,*) therms(i)%x, therms(i)%y, therms(i)%z, &
                     therms(i)%rx, therms(i)%ry, therms(i)%rz, &
                     therms(i)%dtPrime, therms(i)%dU, therms(i)%dV, therms(i)%dW
    end do
    close (iounit)
    
    pi = 4.0 * atan(1.0)

    ! base state vertical profiles
    ! each z-slice is the same
   
    do k = 1, m%nz
       ! height in meters at Theta,u,pPrime levels
       zed(k) = m%dz * 0.5 + m%dz * real(k - 1)

       ! temperature in Kelvins 
       Tbar(k) = Thetabar - g/c_p * zed(k)

       ! pressure in Pascals
       Pbar(k) = P_0 * (Tbar(k)/Thetabar)**(c_p/R_d)

       ! density
       rho_t(k) = Pbar(k)/(R_d * Tbar(k))
    end do

    ! rho_w is staggered from rho_t
    do k = 2, m%nz
       rho_w(k) = 0.5 * (rho_t(k) + rho_t(k-1))
    end do
    rho_w(1) = flag
    rho_w(m%nz+1) = flag
       
    ! Loop over physical domain.
    t = Thetabar
    v = 0.0
    
    do ii = 1, numthermals 
    do k = 1, m%nz    
       do j = 1, m%ny
          do i = 1, m%nx

             ! physical coordinate generation
             x = real(i - 1) * m%dx + m%dx * 0.5
             y = real(j - 1) * m%dy + m%dy * 0.5
             z = real(k - 1) * m%dz + m%dz * 0.5
                        
             d = sqrt(((x - therms(ii)%x)/therms(ii)%rx)**2 + &
                       ((y - therms(ii)%y)/therms(ii)%ry)**2 + &
                       ((z - therms(ii)%z)/therms(ii)%rz)**2)
             
             if (d <= 1.0) then 
                t(i,j,k) = t(i,j,k) + (therms(ii)%dtPrime * (cos(d * pi) + 1.0)*0.5)
                v(i,j,k) = v(i,j,k) + (therms(ii)%dV * (cos(d * pi) + 1.0)*0.5)
             end if
             
          enddo
       enddo
    enddo
    enddo 
  end subroutine ic
end module initial_conditions
