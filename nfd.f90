! 3D nonlinear quasi-compressible flow simulation.
! 
! Mark Van Moer
! Based on work done for ATMS 502, Fall 2011
! University of Illinois at Urbana-Champaign
! taught by Brian Jewett.
!  
program nfd 
  use initial_conditions
  use boundary_conditions
  use integration
  use pgf_bouyancy
  use diffusion
  use savedata
  use mesh_type
  use fson
  implicit none
  type(mesh) m
  real :: dt

  real :: CTdt ! dt for centered time
 
  real, dimension(:,:,:), allocatable :: u1, u2, u3, v1, v2, v3, w1, w2, w3
  real, dimension(:,:,:), allocatable :: t1, t2, tprime, t1a
  real, dimension(:,:,:), allocatable :: p1, p2, p3
  real, dimension(:), allocatable :: rho_t
  real, dimension(:), allocatable :: rho_w

  real :: Thetabar = 300.0

  ! strings
  character(len=*), parameter :: fmt='(I5,1X,F4.0,3X,F7.3,1X,F7.3,2X,F7.3,1X,F7.3,2X,F7.3,1X,F6.3,2X,F8.4,2X,F6.4,2X,F8.2,1X,F7.2)'

  integer :: n, nstep, nplot
  integer :: i, j, k
  character(len=80) :: fsonFile
  type(fson_value), pointer :: config
  
  call get_command_argument(1, fsonFile)

  config => fson_parse(fsonFile)
  call fson_get(config, "nx", m%nx)
  call fson_get(config, "ny", m%ny)
  call fson_get(config, "nz", m%nz)
  call fson_get(config, "dx", m%dx)
  call fson_get(config, "dy", m%dy)
  call fson_get(config, "dz", m%dz)
  call fson_get(config, "dt", dt)
  call fson_get(config, "nstep", nstep)
  call fson_get(config, "nplot", nplot)

  allocate(u1(0:m%nx+1,0:m%ny+1,0:m%nz+1), &
           u2(0:m%nx+1,0:m%ny+1,0:m%nz+1), &
           u3(0:m%nx+1,0:m%ny+1,0:m%nz+1), &
           v1(0:m%nx+1,0:m%ny+1,0:m%nz+1), &
           v2(0:m%nx+1,0:m%ny+1,0:m%nz+1), &
           v3(0:m%nx+1,0:m%ny+1,0:m%nz+1), &
           w1(0:m%nx+1,0:m%ny+1,0:m%nz+1), &
           w2(0:m%nx+1,0:m%ny+1,0:m%nz+1), &
           w3(0:m%nx+1,0:m%ny+1,0:m%nz+1))

  allocate(t1(-1:m%nx+2,-1:m%ny+2,-1:m%nz+2), &
           t2(-1:m%nx+2,-1:m%ny+2,-1:m%nz+2), &
           tprime(-1:m%nx+2,-1:m%ny+2,-1:m%nz+2), &
           t1a(-1:m%nx+2,-1:m%ny+2,-1:m%nz+2))

  allocate(p1(0:m%nx+1,0:m%ny+1,m%nz), &
           p2(0:m%nx+1,0:m%ny+1,m%nz), &
           p3(0:m%nx+1,0:m%ny+1,m%nz))

  allocate(rho_t(1:m%nz), rho_w(1:m%nz+1))

  ! set initial conditions
  call ic(t1(1:m%nx,1:m%ny,1:m%nz), u1(1:m%nx+1,1:m%ny,1:m%nz), v1(1:m%nx,1:m%ny+1,1:m%nz), &
       rho_t, rho_w, m, fsonFile)

  ! everything else is zeros
  p1 = 0.0
  p2 = 0.0
  p3 = 0.0
  u2 = 0.0
  u3 = 0.0
  v2 = 0.0
  v3 = 0.0
  w1 = 0.0
  w2 = 0.0
  w3 = 0.0

  ! plot initial condition
  tprime = t1 - Thetabar
  t2 = t1
  call plotunstaggered(0)


  print*, "Step  Time   umin    umax    vmin     vmax     wmin    wmax   tmin      tmax        pmin     pmax"
  print*, "-------------------------------------------------------------------------------------------------"
  write(*,fmt) 0, 0.0, minval(u1), maxval(u1), minval(v1), maxval(v1), &
       minval(w1), maxval(w1), minval(tprime(1:m%nx,1:m%ny,1:m%nz)), &
       maxval(tprime), minval(p1), maxval(p1)

  call bc(u1, u2, u3, v1, v2, v3, w1, w2, w3, &
       t1, t2, p1, p2, m)

  CTdt = dt

  !integration loop
  do n = 1, nstep
  
     ! start the leapfrogging
     t1a = t1
     t2 = t1
     u3 = u1
     v3 = v1
     w3 = w1
     p3 = p1

     call integrate(t1, t2, u1, v1, w1, u2, v2, w2, &
          u3, v3, w3, m, dt, CTdt)

     ! t1a is unadvected t1, t2 is advected t1
     call diff(t1, u3, v3, w3, u1, v1, w1, &
          m, CTdt, fsonFile) 

     ! call PGF
     call pgf(u3, v3, w3, p3, p1, t2, rho_t, rho_w, &
          m, CTdt, fsonFile)

     ! array update
     ! -- currently doing theta update in advect1D
     if (n .ne. 1) then
        u1 = u2
        v1 = v2
        w1 = w2
        p1 = p2
     end if

     u2 = u3
     v2 = v3
     w2 = w3
     p2 = p3

     ! for CT, first iteration just works like FT
     if (n == 1) CTdt = 2.0 * dt

     ! call bc for next iteration
     call bc(u1, u2, u3, v1, v2, v3, w1, w2, w3, &
          t1, t2, p1, p2, m)

     tprime = t1 - Thetabar   

     if (mod(n,nplot) == 0) then
        call plotunstaggered(n)
     end if

     !print intermediate values every 0..9, then by tens
     if (n < 10 .or. mod(n, 10) == 0 .or. n == nstep) then     
        write(*,fmt) n, dt * real(n), &
             minval(u3(1:m%nx+1,1:m%ny,1:m%nz)), maxval(u3(1:m%nx+1,1:m%ny,1:m%nz)), &
             minval(v3(1:m%nx,1:m%ny+1,1:m%nz)), maxval(v3(1:m%nx,1:m%ny+1,1:m%nz)), &
             minval(w3(1:m%nx,1:m%ny,1:m%nz+1)), maxval(w3(1:m%nx,1:m%ny,1:m%nz+1)), &
             minval(tprime(1:m%nx,1:m%ny,1:m%nz)), maxval(tprime(1:m%nx,1:m%ny,1:m%nz)), &
             minval(p3(1:m%nx,1:m%ny,1:m%nz)), maxval(p3(1:m%nx,1:m%ny,1:m%nz))
     end if

  end do ! time loop finished     
  deallocate(u1, u2, u3, v1, v2, v3, w1, w2, w3, t1, t2, tprime, t1a, p1, p2, p3, rho_t, rho_w)

contains
  ! unstaggers u,v,w and plots along with theta prime and p
  ! Tried OMP triple-loops here, but this was causing crashes
  ! Going with regular loops.
  subroutine plotunstaggered(step)
    integer, intent(in) :: step
    real, dimension(1:m%nx,1:m%ny,1:m%nz) :: plotU,plotV,plotW

    if (step < 1.0) then
       ! for initial conditions, plot u1,w1,v1
       do k = 1, m%nz
          do j = 1, m%ny
             do i = 1,m%nx  
                plotW(i,j,k) = 0.5 * (w1(i,j,k+1) + w1(i,j,k))
             end do
          end do
       end do

       do k = 1,m%nz
          do j = 1,m%ny
             do i = 1,m%nx      
                plotU(i,j,k) = 0.5 * (u1(i+1,j,k) + u1(i,j,k))
             end do
          end do
       end do

       do k = 1,m%nz
          do j = 1,m%ny
             do i = 1,m%nx
                plotV(i,j,k) = 0.5 * (v1(i,j+1,k) + v1(i,j,k))
             end do
          end do
       end do
    else
       ! all other steps plot u3,v3,w3
       do k = 1, m%nz
          do j = 1, m%ny
             do i = 1,m%nx  
                plotW(i,j,k) = 0.5 * (w3(i,j,k+1) + w3(i,j,k))
             end do
          end do
       end do

       do k = 1,m%nz
          do j = 1,m%ny
             do i = 1,m%nx      
                plotU(i,j,k) = 0.5 * (u3(i+1,j,k) + u3(i,j,k))
             end do
          end do
       end do

       do k = 1,m%nz
          do j = 1,m%ny
             do i = 1,m%nx
                plotV(i,j,k) = 0.5 * (v3(i,j+1,k) + v3(i,j,k))
             end do
          end do
       end do
    end if

    call savefield("test", step, plotU, plotV, plotW, m)
    call savefield("T", step, tprime(1:m%nx,1:m%ny,1:m%nz), m) 
    call savefield("P", step, p3(1:m%nx,1:m%ny,1:m%nz), m)
  end subroutine plotunstaggered
end program nfd 
