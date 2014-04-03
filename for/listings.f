!listing00
! code licensed under the terms of GNU GPL v3
! copyright holder: University of Warsaw
!listing01
module real_m
  implicit none
  integer, parameter :: real_t = kind(0.d0) 
end module
!listing02
module arrvec_m
  use real_m
  implicit none

  type :: arr_t
    real(real_t), allocatable :: a(:,:)
  end type

  type :: arrptr_t
    class(arr_t), pointer :: p
  end type

  type :: arrvec_t
    class(arr_t), allocatable :: arrs(:)
    class(arrptr_t), allocatable :: at(:)
    integer :: length
    contains
    procedure :: ctor => arrvec_ctor
    procedure :: init => arrvec_init
  end type

  contains

  subroutine arrvec_ctor(this, n)
    class(arrvec_t) :: this
    integer, intent(in) :: n

    this%length = n
    allocate(this%at( -n : n-1 ))
    allocate(this%arrs( 0 : n-1 ))
  end subroutine

  subroutine arrvec_init(this, n, i, j)
    class(arrvec_t), target :: this
    integer, intent(in) :: n
    integer, intent(in) :: i(2), j(2)

    allocate(this%arrs(n)%a( i(1) : i(2), j(1) : j(2) ))
    this%at(n)%p => this%arrs(n)
    this%at(n - this%length)%p => this%arrs(n)
  end subroutine
end module
!listing03
module arakawa_c_m
  implicit none

  type :: half_t
  end type

  type(half_t) :: h

  interface operator (+)
    module procedure ph
  end interface

  interface operator (-)
    module procedure mh
  end interface

  contains

  elemental function ph(i, h) result (return)
    integer, intent(in) :: i
    type(half_t), intent(in) :: h
    integer :: return
    return = i 
  end function

  elemental function mh(i, h) result (return)
    integer, intent(in) :: i 
    type(half_t), intent(in) :: h
    integer :: return
    return = i - 1
  end function
end module
!listing04
module pi_m
  use real_m
  implicit none
  contains
  function pi(d, arr, i, j) result(return)
    integer, intent(in) :: d
    real(real_t), allocatable, target :: arr(:,:)
    real(real_t), pointer :: return(:,:)
    integer, intent(in) :: i(2), j(2)
    select case (d) 
      case (0) 
        return => arr( i(1) : i(2), j(1) : j(2) )   
      case (1) 
        return => arr( j(1) : j(2), i(1) : i(2) )   
    end select
  end function

  pure function span(d, i, j) result(return)
    integer, intent(in) :: i(2), j(2)
    integer, intent(in) :: d
    integer :: return
    select case (d)
      case (0)
        return = i(2) - i(1) + 1
      case (1)
        return = j(2) - j(1) + 1
    end select
  end function
end module
!listing05
module donorcell_m
  use real_m
  use arakawa_c_m
  use pi_m
  use arrvec_m
  implicit none
  contains 
!listing06
  elemental function F(psi_l, psi_r, C) result (return)
    real(real_t) :: return
    real(real_t), intent(in) :: psi_l, psi_r, C
    return = (                                         &
      (C + abs(C)) * psi_l +                           &
      (C - abs(C)) * psi_r                             &
    ) / 2
  end function
!listing07
  function donorcell(d, psi, C, i, j) result (return)
    integer :: d
    integer, intent(in) :: i(2), j(2) 
    real(real_t) :: return(span(d, i, j), span(d, j, i))
    real(real_t), allocatable, intent(in) :: psi(:,:), C(:,:)           
    return = (                                         &
      F(                                               &
        pi(d, psi, i,   j),                            &
        pi(d, psi, i+1, j),                            &
        pi(d,   C, i+h, j)                             &
      ) -                                              &
      F(                                               &
        pi(d, psi, i-1, j),                            &
        pi(d, psi, i,   j),                            &
        pi(d, C,   i-h, j)                             &
      )                                                &
    )
  end function
!listing08
  subroutine donorcell_op(psi, n, C, i, j)  
    class(arrvec_t), allocatable :: psi
    class(arrvec_t), pointer :: C
    integer, intent(in) :: n
    integer, intent(in) :: i(2), j(2) 
    
    real(real_t), pointer :: ptr(:,:)
    ptr => pi(0, psi%at(n+1)%p%a, i, j)
    ptr = pi(0, psi%at(n)%p%a, i, j) - (               &
      donorcell(0, psi%at(n)%p%a, C%at(0)%p%a, i, j) + &
      donorcell(1, psi%at(n)%p%a, C%at(1)%p%a, j, i)   &
    )
  end subroutine
!listing09
end module
!listing10
module mpdata_m
  use arrvec_m
  use arakawa_c_m
  use pi_m
  implicit none
  contains 
!listing11
  function mpdata_frac(nom, den) result (return)
    real(real_t), intent(in) :: nom(:,:), den(:,:)
    real(real_t) :: return(size(nom, 1), size(nom, 2))
    where (den > 0)
      return = nom / den
    elsewhere
      return = 0
    end where
  end function
!listing12
  function mpdata_A(d, psi, i, j) result (return)
    integer :: d
    real(real_t), allocatable, intent(in) :: psi(:,:)
    integer, intent(in) :: i(2), j(2)
    real(real_t) :: return(span(d, i, j), span(d, j, i))
    return = mpdata_frac(                              &
      pi(d, psi, i+1, j) - pi(d, psi, i, j),           &
      pi(d, psi, i+1, j) + pi(d, psi, i, j)            &
    )  
  end function
!listing13
  function mpdata_B(d, psi, i, j) result (return)
    integer :: d
    real(real_t), allocatable, intent(in) :: psi(:,:) 
    integer, intent(in) :: i(2), j(2)
    real(real_t) :: return(span(d, i, j), span(d, j, i))
    return = mpdata_frac(                              &
      pi(d, psi, i+1, j+1) + pi(d, psi, i,   j+1)      &
    - pi(d, psi, i+1, j-1) - pi(d, psi, i,   j-1),     &
      pi(d, psi, i+1, j+1) + pi(d, psi, i,   j+1)      &
    + pi(d, psi, i+1, j-1) + pi(d, psi, i,   j-1)      &
    ) / 2
  end function
!listing14
  function mpdata_C_bar(d, C, i, j) result (return)
    integer :: d
    real(real_t), allocatable, intent(in) :: C(:,:) 
    integer, intent(in) :: i(2), j(2)
    real(real_t) :: return(span(d, i, j), span(d, j, i))

    return = (                                         &
      pi(d, C, i+1, j+h) + pi(d, C, i,   j+h) +        &
      pi(d, C, i+1, j-h) + pi(d, C, i,   j-h)          &
    ) / 4               
  end function
!listing15
  function mpdata_C_adf(d, psi, i, j, C) result (return)
    integer :: d
    integer, intent(in) :: i(2), j(2)
    real(real_t) :: return(span(d, i, j), span(d, j, i))
    real(real_t), allocatable, intent(in) :: psi(:,:) 
    class(arrvec_t), pointer :: C
    return =                                           &
      abs(pi(d, C%at(d)%p%a, i+h, j))                  &
      * (1 - abs(pi(d, C%at(d)%p%a, i+h, j)))          &
      * mpdata_A(d, psi, i, j)                         &
      - pi(d, C%at(d)%p%a, i+h, j)                     &
      * mpdata_C_bar(d, C%at(d-1)%p%a, i, j)           &
      * mpdata_B(d, psi, i, j)
  end function
!listing16
end module
!listing17
module halo_m
  use arakawa_c_m
  implicit none

  interface ext
    module procedure ext_n
    module procedure ext_h
  end interface 

  contains

  function ext_n(r, n) result (return)
    integer, intent(in) :: r(2)
    integer, intent(in) :: n
    integer :: return(2)
    
    return = (/ r(1) - n, r(2) + n /)
  end function

  function ext_h(r, h) result (return)
    integer, intent(in) :: r(2)
    type(half_t), intent(in) :: h
    integer :: return(2)
    
    return = (/ r(1) - h, r(2) + h /)
  end function
end module
!listing18
module bcd_m
  use arrvec_m
  implicit none

  type, abstract :: bcd_t
    contains
    procedure(bcd_fill_halos), deferred :: fill_halos
    procedure(bcd_init), deferred :: init
  end type
 
  abstract interface 
    subroutine bcd_fill_halos(this, a, j)
      import :: bcd_t, real_t
      class(bcd_t ) :: this
      real(real_t), allocatable :: a(:,:) 
      integer :: j(2)
    end subroutine

    subroutine bcd_init(this, d, n, hlo)
      import :: bcd_t
      class(bcd_t) :: this
      integer :: d, n, hlo
    end subroutine
  end interface
end module
!listing19
module solver_m
  use arrvec_m
  use bcd_m
  use arakawa_c_m
  use halo_m
  implicit none

  type, abstract :: solver_t
    class(arrvec_t), allocatable :: psi, C
    integer :: n, hlo
    integer :: i(2), j(2) 
    class(bcd_t), pointer :: bcx, bcy
    contains
    procedure :: solve   => solver_solve
    procedure :: state   => solver_state
    procedure :: courant => solver_courant
    procedure :: cycle   => solver_cycle
    procedure(solver_advop), deferred :: advop
  end type 

  abstract interface
    subroutine solver_advop(this)
      import solver_t
      class(solver_t), target :: this
    end subroutine
  end interface

  contains

  subroutine solver_ctor(this, bcx, bcy, nx, ny, hlo)
    use arakawa_c_m
    use halo_m
    class(solver_t) :: this
    class(bcd_t), intent(in), target :: bcx, bcy
    integer, intent(in) :: nx, ny, hlo

    this%n = 0
    this%hlo = hlo
    this%bcx => bcx
    this%bcy => bcy

    this%i = (/ 0, nx - 1 /)
    this%j = (/ 0, ny - 1 /)

    call bcx%init(0, nx, hlo)
    call bcy%init(1, ny, hlo)

    allocate(this%psi)
    call this%psi%ctor(2)
    block
      integer :: n
      do n=0, 1
        call this%psi%init(                            &
          n, ext(this%i, hlo), ext(this%j, hlo)        &
        )
      end do
    end block

    allocate(this%C)
    call this%C%ctor(2)
    call this%C%init(                                  &
      0, ext(this%i, h), ext(this%j, hlo)              &
    )
    call this%C%init(                                  &
      1, ext(this%i, hlo), ext(this%j, h)              &
    )
  end subroutine

  function solver_state(this) result (return)
    class(solver_t) :: this
    real(real_t), pointer :: return(:,:)
    return => this%psi%at(this%n)%p%a(                 &
      this%i(1) : this%i(2),                           &
      this%j(1) : this%j(2)                            &
    )
  end function

  function solver_courant(this, d) result (return)
    class(solver_t) :: this
    integer :: d
    real(real_t), pointer :: return(:,:)
    return => this%C%at(d)%p%a
  end function

  subroutine solver_cycle(this)
    class(solver_t) :: this
    this%n = mod(this%n + 1 + 2, 2) - 2
  end subroutine

  subroutine solver_solve(this, nt)
    class(solver_t) :: this
    integer, intent(in) :: nt
    integer :: t

    do t = 0, nt-1 
      call this%bcx%fill_halos(                        &
        this%psi%at(this%n)%p%a, ext(this%j, this%hlo) &
      )
      call this%bcy%fill_halos(                        &
        this%psi%at(this%n)%p%a, ext(this%i, this%hlo) &
      )
      call this%advop()
      call this%cycle()
    end do
  end subroutine
end module
!listing20
module cyclic_m
  use bcd_m
  use pi_m
  implicit none
  
  type, extends(bcd_t) :: cyclic_t
    integer :: d
    integer :: left_halo(2), rght_halo(2) 
    integer :: left_edge(2), rght_edge(2) 
    contains
    procedure :: init => cyclic_init
    procedure :: fill_halos => cyclic_fill_halos
  end type

  contains

  subroutine cyclic_init(this, d, n, hlo)
    class(cyclic_t) :: this
    integer :: d, n, hlo

    this%d = d
    this%left_halo = (/ -hlo, -1 /) 
    this%rght_halo = (/ n, n-1+hlo /) 
    this%left_edge = (/ 0, hlo-1 /)
    this%rght_edge = (/ n-hlo, n-1 /)
  end subroutine

  subroutine cyclic_fill_halos(this, a, j)
    class(cyclic_t) :: this
    real(real_t), pointer :: ptr(:,:)
    real(real_t), allocatable :: a(:,:)
    integer :: j(2)
    ptr => pi(this%d, a, this%left_halo, j) 
    ptr =  pi(this%d, a, this%rght_edge, j)
    ptr => pi(this%d, a, this%rght_halo, j) 
    ptr =  pi(this%d, a, this%left_edge, j)
  end subroutine
end module
!listing21
module solver_donorcell_m
  use donorcell_m
  use solver_m
  implicit none
  
  type, extends(solver_t) :: donorcell_t
    contains
    procedure :: ctor => donorcell_ctor
    procedure :: advop => donorcell_advop
  end type

  contains
  
  subroutine donorcell_ctor(this, bcx, bcy, nx, ny)
    class(donorcell_t) :: this
    class(bcd_t), intent(in), target :: bcx, bcy
    integer, intent(in) :: nx, ny
    call solver_ctor(this, bcx,bcy, nx,ny, 1)
  end subroutine

  subroutine donorcell_advop(this)
    class(donorcell_t), target :: this
    class(arrvec_t), pointer :: C
    C => this%C
    call donorcell_op(                                 &
      this%psi, this%n, C, this%i, this%j              &
    )
  end subroutine
end module
!listing22
module solver_mpdata_m
  use solver_m
  use mpdata_m
  use donorcell_m
  use halo_m
  implicit none
  
  type, extends(solver_t) :: mpdata_t
    integer :: n_iters, n_tmp
    integer :: im(2), jm(2)
    class(arrvec_t), pointer :: tmp(:) 
    contains
    procedure :: ctor => mpdata_ctor
    procedure :: advop => mpdata_advop
  end type

  contains

  subroutine mpdata_ctor(this, n_iters, bcx, bcy, nx, ny)
    class(mpdata_t) :: this
    class(bcd_t), target :: bcx, bcy
    integer, intent(in) :: n_iters, nx, ny
    integer :: c

    call solver_ctor(this, bcx, bcy, nx, ny, 1)

    this%n_iters = n_iters
    this%n_tmp = min(n_iters - 1, 2)
    if (n_iters > 0) allocate(this%tmp(0:this%n_tmp)) 

    associate (i => this%i, j => this%j, hlo => this%hlo)
      do c=0, this%n_tmp - 1
        call this%tmp(c)%ctor(2)
        call this%tmp(c)%init(0, ext(i, h), ext(j, hlo))
        call this%tmp(c)%init(1, ext(i, hlo), ext(j, h))
      end do

      this%im = (/ i(1) - 1, i(2) /)
      this%jm = (/ j(1) - 1, j(2) /)
    end associate
  end subroutine

  subroutine mpdata_advop(this)
    class(mpdata_t), target :: this
    integer :: step

    associate (i => this%i, j => this%j, im => this%im,&
      jm => this%jm, psi => this%psi, n => this%n,     &
      hlo => this%hlo, bcx => this%bcx, bcy => this%bcy&
    )
      do step=0, this%n_iters-1
        if (step == 0) then
          block
            class(arrvec_t), pointer :: C
            C => this%C
            call donorcell_op(psi, n, C, i, j)
          end block
        else
          call this%cycle()
          call bcx%fill_halos(                         &
            psi%at( n )%p%a, ext(j, hlo)               &
          )
          call bcy%fill_halos(                         &
            psi%at( n )%p%a, ext(i, hlo)               &
          )

          block
            class(arrvec_t), pointer :: C_corr, C_unco
            real(real_t), pointer :: ptr(:,:)

            ! chosing input/output for antidiff. C
            if (step == 1) then
              C_unco => this%C
              C_corr => this%tmp(0)
            else if (mod(step, 2) == 1) then
              C_unco => this%tmp(1) ! odd step
              C_corr => this%tmp(0) ! even step
            else
              C_unco => this%tmp(0) ! odd step
              C_corr => this%tmp(1) ! even step
            end if

            ! calculating the antidiffusive velo
            ptr => pi(0, C_corr%at( 0 )%p%a, im+h, j)
            ptr = mpdata_C_adf(                        &
              0, psi%at( n )%p%a, im, j, C_unco        &
            )      
            call bcy%fill_halos(                       &
              C_corr%at(0)%p%a, ext(i, h)              &
            )

            ptr => pi(0, C_corr%at( 1 )%p%a, i, jm+h)
            ptr = mpdata_C_adf(                        &
              1, psi%at( n )%p%a, jm, i, C_unco        &
            )
            call bcx%fill_halos(                       &
              C_corr%at(1)%p%a, ext(j, h)              &
            )

            ! donor-cell step
            call donorcell_op(psi, n, C_corr, i, j) 
          end block
        end if
      end do
    end associate
  end subroutine
end module
!listing23
