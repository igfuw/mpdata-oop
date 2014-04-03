program test
  use solver_mpdata_m
  use cyclic_m
  implicit none

  if (command_argument_count() /= (9)) then
    print*, "expecting 9 arguments (nx, ny, Cx, Cy, nt, it, f.in, f.out, dec)"
    stop 1
  end if

  block
!listing23
    type(mpdata_t) :: slv
    type(cyclic_t), target :: bcx, bcy
    integer :: nx, ny, nt, it
    real(real_t) :: Cx, Cy
    real(real_t), pointer :: ptr(:,:)
!listing24
    character (len=666) :: fname
    integer :: dec

    nx = arg2int(1)
    ny = arg2int(2)
    Cx = arg2real(3)
    Cy = arg2real(4)
    nt = arg2int(5)
    it = arg2int(6)
    dec = arg2int(9)
    call get_command_argument(7, fname)

!listing25
    call slv%ctor(it, bcx, bcy, nx, ny)

    ptr => slv%state() 
    call read_file(fname, ptr)

    ptr => slv%courant(0) 
    ptr = Cx

    ptr => slv%courant(1) 
    ptr = Cy

    call slv%solve(nt)
!listing26

    block
      logical :: error = .FALSE.
      real(real_t), pointer :: tmp(:,:)
      real :: diff
      character (len=666) :: fname
      allocate(tmp(nx,ny))
      call get_command_argument(8, fname)
      call read_file(fname, tmp)
      diff = maxval(abs(slv%state() - tmp))
      print*, "diff=", diff
      if (diff >= .5 * 10.**(-dec)) error = .TRUE.
      deallocate(tmp)
      if (error) then
        !print*, slv%state()
        stop 1
      end if
    end block
  end block

  contains

  function arg2int(argno) result(return)
    implicit none
    integer :: argno, return
    character (len=666) :: tmp
    call get_command_argument(argno, tmp)
    read(tmp,*)return
  end function

  function arg2real(argno) result(return)
    implicit none
    integer :: argno
    real(real_t) :: return
    character (len=666) :: tmp
    call get_command_argument(argno, tmp)
    read(tmp,*)return
  end function

  subroutine read_file(fname, array)
    character (len=*), intent (in) :: fname
    real(real_t), pointer :: array(:,:)
    integer :: un
    open(newunit=un, file=fname, status='old', action='read' )
    block    
      integer :: i, j
      do i=1, size(array, 1)
        read(un, *) (array (i, j), j=1, size(array, 2))
      end do
    end block
    close (un)
   end subroutine read_file
end program
