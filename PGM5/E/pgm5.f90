!     Program 5 - 3D nonlinear quasi-compressible flow
!     ATMS502 Fall 2011
!     Mark Van Moer
!
!     Make with supplied makefile
!
!  Test E
!
program pgm5
  use initial_conditions
  use boundary_conditions
  use integration
  use pgf_bouyancy
  use diffusion
  use savedata
  use read_parameters_module 
  implicit none

  !integer, parameter :: nx = 33, ny = 33, nz = 16
  integer :: nx, ny, nz
  real, parameter :: dx = 500.0, dy = 500.0, dz = 500.0
  real :: dt = 1.0

  real :: CTdt ! dt for centered time

  ! Ah, these now have to become heap, not stack.
  ! momentum for n-1,n,n+1
  !real, dimension(0:nx+1,0:ny+1,0:nz+1) :: u1, u2, u3, v1, v2, v3, w1, w2, w3
  ! scalar
  !real, dimension(-1:nx+2,-1:ny+2,-1:nz+2) :: t1, t2, tprime, t1a
  ! pressure
  !real, dimension(0:nx+1,0:ny+1,nz) :: p1, p2, p3
  ! density
  !real, dimension(1:nz) :: rho_t
  !real, dimension(1:nz+1) :: rho_w
 
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

  call read_parameters("inputfile", nx, ny, nz) 
  allocate(u1(0:nx+1,0:ny+1,0:nz+1), &
           u2(0:nx+1,0:ny+1,0:nz+1), &
           u3(0:nx+1,0:ny+1,0:nz+1), &
           v1(0:nx+1,0:ny+1,0:nz+1), &
           v2(0:nx+1,0:ny+1,0:nz+1), &
           v3(0:nx+1,0:ny+1,0:nz+1), &
           w1(0:nx+1,0:ny+1,0:nz+1), &
           w2(0:nx+1,0:ny+1,0:nz+1), &
           w3(0:nx+1,0:ny+1,0:nz+1))

  allocate(t1(-1:nx+2,-1:ny+2,-1:nz+2), &
           t2(-1:nx+2,-1:ny+2,-1:nz+2), &
           tprime(-1:nx+2,-1:ny+2,-1:nz+2), &
           t1a(-1:nx+2,-1:ny+2,-1:nz+2))

  allocate(p1(0:nx+1,0:ny+1,nz), &
           p2(0:nx+1,0:ny+1,nz), &
           p3(0:nx+1,0:ny+1,nz))

  allocate(rho_t(1:nz), rho_w(1:nz+1))

   
  nstep = 450
  nplot = 50

  ! set initial conditions
  call ic(t1(1:nx,1:ny,1:nz), u1(1:nx+1,1:ny,1:nz), v1(1:nx,1:ny+1,1:nz), &
       rho_t, rho_w, nx, ny, nz, dx, dy, dz)

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
  call plotunstaggered(0)


  print*, "Step  Time   umin    umax    vmin     vmax     wmin    wmax   tmin      tmax        pmin     pmax"
  print*, "-------------------------------------------------------------------------------------------------"
  write(*,fmt) 0, 0.0, minval(u1), maxval(u1), minval(v1), maxval(v1), &
       minval(w1), maxval(w1), minval(tprime(1:nx,1:ny,1:nz)), &
       maxval(tprime), minval(p1), maxval(p1)

  call bc(u1, u2, u3, v1, v2, v3, w1, w2, w3, &
       t1, t2, p1, p2, nx, ny, nz)

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

     ! call advect
     ! Passes Tests A1-3 for theta advection
     ! Passes Test D for Theta only advection
     call integrate(t1, t2, u1, v1, w1, u2, v2, w2, &
          u3, v3, w3, dx, dy, dz, dt, nx, ny, nz, CTdt)
    
     ! call diffusion
     ! Passes Test C
     ! t1a is unadvected t1, t2 is advected t1
     call diff(t1a, u3, v3, w3, u1, v1, w1, &
          dx, dy, dz, nx, ny, nz, CTdt) 

     ! call PGF
     ! Passes Test B
     call pgf(u3, v3, w3, p3, p1, t2, rho_t, rho_w, &
          dx, dy, dz, nx, ny, nz, CTdt)

     

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
          t1, t2, p1, p2, nx, ny, nz)

     tprime = t1 - Thetabar   

     if (mod(n,nplot) == 0) then
        call plotunstaggered(n)
     end if

     !print intermediate values every 0..9, then by tens
     if (n < 10 .or. mod(n, 10) == 0 .or. n == nstep) then     
        write(*,fmt) n, dt * real(n), &
             minval(u3(1:nx+1,1:ny,1:nz)), maxval(u3(1:nx+1,1:ny,1:nz)), &
             minval(v3(1:nx,1:ny+1,1:nz)), maxval(v3(1:nx,1:ny+1,1:nz)), &
             minval(w3(1:nx,1:ny,1:nz+1)), maxval(w3(1:nx,1:ny,1:nz+1)), &
             minval(tprime(1:nx,1:ny,1:nz)), maxval(tprime(1:nx,1:ny,1:nz)), &
             minval(p3(1:nx,1:ny,1:nz)), maxval(p3(1:nx,1:ny,1:nz))
     end if

  end do ! time loop finished     
  deallocate(u1, u2, u3, v1, v2, v3, w1, w2, w3, t1, t2, tprime, t1a, p1, p2, p3, rho_t, rho_w)

contains
  ! unstaggers u,v,w and plots along with theta prime and p
  ! Tried OMP triple-loops here, but this was causing the crashes
  ! Going with regular loops.
  subroutine plotunstaggered(step)
    integer, intent(in) :: step
    real, dimension(1:nx,1:ny,1:nz) :: plotU,plotV,plotW

    if (step < 1.0) then
       ! for initial conditions, plot u1,w1,v1
       do k = 1, nz
          do j = 1, ny
             do i = 1,nx  
                plotW(i,j,k) = 0.5 * (w1(i,j,k+1) + w1(i,j,k))
             end do
          end do
       end do

       do k = 1,nz
          do j = 1,ny
             do i = 1,nx      
                plotU(i,j,k) = 0.5 * (u1(i+1,j,k) + u1(i,j,k))
             end do
          end do
       end do

       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
                plotV(i,j,k) = 0.5 * (v1(i,j+1,k) + v1(i,j,k))
             end do
          end do
       end do
    else
       ! all other steps plot u3,v3,w3
       do k = 1, nz
          do j = 1, ny
             do i = 1,nx  
                plotW(i,j,k) = 0.5 * (w3(i,j,k+1) + w3(i,j,k))
             end do
          end do
       end do

       do k = 1,nz
          do j = 1,ny
             do i = 1,nx      
                plotU(i,j,k) = 0.5 * (u3(i+1,j,k) + u3(i,j,k))
             end do
          end do
       end do

       do k = 1,nz
          do j = 1,ny
             do i = 1,nx
                plotV(i,j,k) = 0.5 * (v3(i,j+1,k) + v3(i,j,k))
             end do
          end do
       end do
    end if

!   call savefield("U", step, plotU, nx, ny, nz)
!   call savefield("V", step, plotV, nx, ny, nz)
!   call savefield("W", step, plotW, nx, ny, nz)
    call savefield("test", step, plotU, plotV, plotW, nx, ny, nz)
    call savefield("T", step, tprime(1:nx,1:ny,1:nz), nx, ny, nz) 
    call savefield("P", step, p3(1:nx,1:ny,1:nz), nx, ny, nz)

  end subroutine plotunstaggered
end program pgm5
