! boundary_conditions.f90
! Mark Van Moer, ATMS502 Fall 2011
! Sets the boundary conditions. Slightly altered from 
! bjewett's website version and modularized.
! removed q1, q2, t0, bctype, noslip
module boundary_conditions
  use mesh_type
  implicit none
  public
contains
  subroutine bc(u1, u2, u3, v1, v2, v3, w1, w2, w3, t1, t2, &
       p1, p2, m)
    
    type (mesh), intent(in) :: m
    ! mvm - can I make these assumed shape?
    real, dimension(0:m%nx+1,0:m%ny+1,0:m%nz+1), intent(inout) :: u1, u2, u3, &
         v1, v2, v3, w1, w2, w3
    real, dimension(-1:m%nx+2,-1:m%ny+2,-1:m%nz+2), intent(inout) :: t1, t2
    real, dimension(0:m%nx+1,0:m%ny+1,1:m%nz), intent(inout) :: p1, p2
    !integer, intent(in) :: nx, ny, nz
    
    integer :: i, j, k

!$OMP     PARALLEL DO PRIVATE (i,j,k)
      do k = 1, m%nz
        do j = 1, m%ny
          u1(1   ,j,k) = -u1(2 ,j,k)
          u2(1   ,j,k) = -u2(2 ,j,k)
          u3(1   ,j,k) = -u3(2 ,j,k)
          u1(0   ,j,k) = -u1(3 ,j,k)
          u2(0   ,j,k) = -u2(3 ,j,k)
          u3(0   ,j,k) = -u3(3 ,j,k)
          u1(m%nx+1,j,k) = -u1(m%nx,j,k)
          u2(m%nx+1,j,k) = -u2(m%nx,j,k)
          u3(m%nx+1,j,k) = -u3(m%nx,j,k)
        enddo
        do i = 1,m%nx
          u1(i,   0,k) = u1(i,  m%ny,k)
          u2(i,   0,k) = u2(i,  m%ny,k)
          u3(i,   0,k) = u3(i,  m%ny,k)
          u1(i,m%ny+1,k) = u1(i, 1  ,k)
          u2(i,m%ny+1,k) = u2(i, 1  ,k)
          u3(i,m%ny+1,k) = u3(i, 1  ,k)
        enddo
      enddo
!$OMP     END PARALLEL DO
!$OMP     PARALLEL DO PRIVATE (i,j,k)
      do k = 1,m%nz
        do i = 1,m%nx
          v1(i,m%ny+1,k) =  v1(i  ,1 ,k)
          v2(i,m%ny+1,k) =  v2(i  ,1 ,k)
          v3(i,m%ny+1,k) =  v3(i  ,1 ,k)
          v1(i   ,0,k) =  v1(i  ,m%ny,k)
          v2(i   ,0,k) =  v2(i  ,m%ny,k)
          v3(i   ,0,k) =  v3(i  ,m%ny,k)
        enddo
        do j = 1,m%ny+1
          v1(0   ,j,k) = v1(2   ,j,k)
          v2(0   ,j,k) = v2(2   ,j,k)
          v3(0   ,j,k) = v3(2   ,j,k)
          v1(m%nx+1,j,k) = v1(m%nx-1,j,k)
          v2(m%nx+1,j,k) = v2(m%nx-1,j,k)
          v3(m%nx+1,j,k) = v3(m%nx-1,j,k)
        enddo
      enddo
!$OMP     END PARALLEL DO
!$OMP     PARALLEL DO PRIVATE (i,j,k)
      do k = 2, m%nz
        do j = 1, m%ny
          w1(0   ,j,k) = w1(2   ,j,k)
          w2(0   ,j,k) = w2(2   ,j,k)
          w3(0   ,j,k) = w3(2   ,j,k)
          w1(m%nx+1,j,k) = w1(m%nx-1,j,k)
          w2(m%nx+1,j,k) = w2(m%nx-1,j,k)
          w3(m%nx+1,j,k) = w3(m%nx-1,j,k)
        enddo
        do i = 1, m%nx
          w1(i,   0,k) = w1(i,  m%ny,k)
          w2(i,   0,k) = w2(i,  m%ny,k)
          w3(i,   0,k) = w3(i,  m%ny,k)
          w1(i,m%ny+1,k) = w1(i,  1 ,k)
          w2(i,m%ny+1,k) = w2(i,  1 ,k)
          w3(i,m%ny+1,k) = w3(i,  1 ,k)
        enddo
      enddo
!$OMP     END PARALLEL DO
!$OMP     PARALLEL DO PRIVATE (i,j,k)
      do k = 1, m%nz
        do j = 1, m%ny
          t1( 0  ,j,k) = t1(2   ,j,k)
          t1(-1  ,j,k) = t1(3   ,j,k)
          t1(m%nx+1,j,k) = t1(m%nx-1,j,k)
          t1(m%nx+2,j,k) = t1(m%nx-2,j,k)
          t2( 0  ,j,k) = t2(2   ,j,k)
          t2(-1  ,j,k) = t2(3   ,j,k)
          t2(m%nx+1,j,k) = t2(m%nx-1,j,k)
          t2(m%nx+2,j,k) = t2(m%nx-2,j,k)
       enddo
        do i = 1,m%nx
          t1(i,   0,k) = t1(i,m%ny  ,k)
          t2(i,   0,k) = t2(i,m%ny  ,k)
          t1(i,  -1,k) = t1(i,m%ny-1,k)
          t2(i,  -1,k) = t2(i,m%ny-1,k)
          t1(i,m%ny+1,k) = t1(i,  1 ,k)
          t2(i,m%ny+1,k) = t2(i,  1 ,k)
          t1(i,m%ny+2,k) = t1(i,  2 ,k)
          t2(i,m%ny+2,k) = t2(i,  2 ,k)
       enddo
      enddo
!$OMP     END PARALLEL DO
!$OMP     PARALLEL DO PRIVATE (i,j,k)
      do k = 1,m%nz
        do j = 1,m%ny
          p1(0   ,j,k) = p1(2   ,j,k)
          p2(0   ,j,k) = p2(2   ,j,k)
          p1(m%nx+1,j,k) = p1(m%nx-1,j,k)
          p2(m%nx+1,j,k) = p2(m%nx-1,j,k)
        enddo
        do i = 1,m%nx
          p1(i,   0,k) = p1(i,m%ny  ,k)
          p2(i,   0,k) = p2(i,m%ny  ,k)
          p1(i,m%ny+1,k) = p1(i,  1 ,k)
          p2(i,m%ny+1,k) = p2(i,  1 ,k)
        enddo
      enddo
!$OMP     END PARALLEL DO

!
! ... 0-gradient top, bottom
!

!$OMP   PARALLEL DO PRIVATE (i,j)
    do j = 0,m%ny+1
      do i = 0,m%nx+1
        u1(i,j,m%nz+1) = u1(i,j,m%nz)
        u2(i,j,m%nz+1) = u2(i,j,m%nz)
        u3(i,j,m%nz+1) = u3(i,j,m%nz)
        u1(i,j,0   ) = u1(i,j,1 )
        u2(i,j,0   ) = u2(i,j,1 )
        u3(i,j,0   ) = u3(i,j,1 )
      enddo
      do i = 0,m%nx+1
        v1(i,j,m%nz+1) = v1(i,j,m%nz)
        v2(i,j,m%nz+1) = v2(i,j,m%nz)
        v3(i,j,m%nz+1) = v3(i,j,m%nz)
        v1(i,j,0   ) = v1(i,j,1 )
        v2(i,j,0   ) = v2(i,j,1 )
        v3(i,j,0   ) = v3(i,j,1 )
      enddo
      do i = 0,m%nx+1
        w1(i,j,0   ) = 0.
        w1(i,j,1   ) = 0.
        w1(i,j,m%nz+1) = 0.
        w2(i,j,0   ) = 0.
        w2(i,j,1   ) = 0.
        w2(i,j,m%nz+1) = 0.
        w3(i,j,0   ) = 0.
        w3(i,j,1   ) = 0.
        w3(i,j,m%nz+1) = 0.
      enddo
    enddo
!$OMP   END PARALLEL DO
!$OMP   PARALLEL DO PRIVATE (i,j)
    do j = -1,m%ny+2
          do i = -1,m%nx+2
            t1(i,j,0   ) = t1(i,j,1 )
            t1(i,j,-1  ) = t1(i,j,1 )
            t1(i,j,m%nz+1) = t1(i,j,m%nz)
            t1(i,j,m%nz+2) = t1(i,j,m%nz)
            t2(i,j,0   ) = t2(i,j,1 )
            t2(i,j,-1  ) = t2(i,j,1 )
            t2(i,j,m%nz+1) = t2(i,j,m%nz)
            t2(i,j,m%nz+2) = t2(i,j,m%nz)
         enddo
    enddo
!$OMP   END PARALLEL DO
  end subroutine bc
end module boundary_conditions

