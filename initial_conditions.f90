!  initial_conditions.f90
!  Sets the initial conditions.
!  First sets the "base state vertical profiles"
module initial_conditions
  use mesh_type
  use fson
  use fson_path_m
  use fson_value_m
  implicit none
  public    
  type thermal
    real :: x, y, z
    real :: rx, ry, rz
    real :: dtPrime
    real :: dU, dV, dW
  end type thermal
  ! needs to be global so callback can see it.
  type (thermal), allocatable, dimension(:) :: therms
contains  
  subroutine ic(t, u, v, rho_t, rho_w, m, fsonFile)
    ! t -- "theta", 3D array of scalars
    ! u -- 3D array of momentum u-component
    ! v -- 3D array of momentum v-component
    ! rho_t -- 1D array of density at Theta, u, pprime levels
    ! rho_w -- 1D array of density at w levels
    ! fsonFile -- path to json config file


    real, dimension(:,:,:), intent(inout) :: t, u, v
    real, dimension(:), intent(inout) :: rho_t, rho_w
    type (mesh), intent(in) :: m
    character(len=*), intent(in) :: fsonFile

    real :: pi
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

    type(fson_value), pointer :: config, therm_array
    character(len=80) :: therm_json_path
    pi = 4.0 * atan(1.0)

    config => fson_parse(fsonFile)
    ! can't seem to pass a literal string
    therm_json_path = "thermals"
    call fson_path_get(config, therm_json_path, therm_array)
    numthermals = fson_value_count(therm_array)
    allocate(therms(numthermals))    
    call fson_get(config, therm_json_path, therm_callback)
    call fson_destroy(config)
    call fson_destroy(therm_array) 

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
! U perturbations
    call random_seed()
    do k = 1, m%nz
       do j = 1, m%ny
          do i = 1, m%nx + 1
             call random_number(ranval)
             u(i,j,k) = (ranval - 0.5) * therms(1)%dU 
          end do
       end do
    end do
  end subroutine ic
  subroutine therm_callback(element, i, c)
    type(fson_value), pointer :: element
    integer :: i, c
    call fson_get(element, "x", therms(i)%x)
    call fson_get(element, "y", therms(i)%y)
    call fson_get(element, "z", therms(i)%z)
    call fson_get(element, "rx", therms(i)%rx)
    call fson_get(element, "ry", therms(i)%ry)
    call fson_get(element, "rz", therms(i)%rz)
    call fson_get(element, "dtPrime", therms(i)%dtPrime)
    call fson_get(element, "dU", therms(i)%dU)
    call fson_get(element, "dV", therms(i)%dV)
    call fson_get(element, "dW", therms(i)%dW)
  end subroutine therm_callback
end module initial_conditions
