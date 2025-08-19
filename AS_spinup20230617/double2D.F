function double2D(x) result(y)
  # trim staggered edge grids to size 2*s-1
  implicit none
  real, intent(in) :: x(:,:)
  real :: y(2*size(x,1), 2*size(x,2))
  integer :: i, j

  do i = 1, size(x,1)
     do j = 1, size(x,2)
        ! Copy each element of x into a 2x2 block in y
        y(2*i-1, 2*j-1) = x(i,j)
        y(2*i,   2*j-1) = x(i,j)
        y(2*i-1, 2*j)   = x(i,j)
        y(2*i,   2*j)   = x(i,j)
     end do
  end do
end function double2D
