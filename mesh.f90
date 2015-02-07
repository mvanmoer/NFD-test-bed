! mesh.f90
! Simple type representing the grid mesh.
module mesh_type
  implicit none
  public
  type mesh
    integer :: nx, ny, nz
    real :: dx, dy, dz
  end type mesh
end module mesh_type
