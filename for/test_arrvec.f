program test_arrvec
  use arrvec_m
  class(arrvec_t), allocatable :: psi
  integer, dimension(2) :: i, j
  integer :: c, nx = 10, ny = 10

  i = (/ 0, nx - 1 /)
  j = (/ 0, ny - 1 /)

  allocate(psi)
  call psi%ctor(2)
  call psi%init(0, i, j)
  call psi%init(1, i, j)

  print*, psi%at(0)%p%a
  print*, psi%at(0)%p%a(1,1)
  psi%at(0)%p%a(1,1) = 10
  print*, psi%at(0)%p%a(1,1)
  print*, psi%at(-2)%p%a(1,1)
end
