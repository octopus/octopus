subroutine specie_init_dim(nspecies, str, s)
  integer, intent(in) :: nspecies
  character(len=*), intent(in) :: str
  type(specie_type), pointer :: s(:)

  integer :: i, j

  do i = 1, nspecies
    call oct_parse_block_str(str, i-1, 0, s(i)%label)
    call oct_parse_block_double(str, i-1, 1, s(i)%weight)
    s(i)%weight =  units_inp%mass%factor * s(i)%weight ! units conversion

    s(i)%local = .true. ! In 1D or 2D, potential has to be local.
    call oct_parse_block_double(str, i-1, 2, s(i)%Z_val)
    call oct_parse_block_str   (str, i-1, 3, s(i)%user_def)
    ! convert to C string
    j = len(trim(s(i)%user_def))
    s(i)%user_def(j+1:j+1) = achar(0) 
  end do

end subroutine specie_init_dim

subroutine specie_end(ns, s)
  integer, intent(in) :: ns
  type(specie_type), pointer :: s(:)

  integer :: i

  sub_name = 'specie_end'; call push_sub()

  if(associated(s)) then ! sanity check
    deallocate(s); nullify(s)
  end if

  call pop_sub()
  return
end subroutine specie_end

real(r8) function specie_get_local(s, x) result(l)
  type(specie_type), intent(IN) :: s
  real(r8), intent(in) :: x(conf%dim)

  real(r8) :: xx(3)

  xx = 0._r8
  xx(1:conf%dim) = x(:)
  l = oct_parse_potential(xx(1), xx(2), xx(3), sqrt(sum(xx(:)**2)), s%user_def)

end function specie_get_local

! returns the gradient of the local potential
! thsi version uses numerical differences
subroutine specie_get_dlocal(s, x, dl)
  type(specie_type), intent(IN) :: s
  real(r8), intent(in) :: x(conf%dim)
  real(r8), intent(out) :: dl(conf%dim)

  real(r8), parameter :: Delta = 1e-4_r8

  integer :: i
  real(r8) :: l1, l2, xx(3)

  xx = 0._r8
  do i = 1, conf%dim
    xx(1:conf%dim) = x(1:conf%dim)
    xx(i) = xx(i) - Delta
    l1 = oct_parse_potential(xx(1), xx(2), xx(3), sqrt(sum(xx(1:conf%dim)**2)), s%user_def)
    xx(i) = xx(i) + M_TWO*Delta
    l2 = oct_parse_potential(xx(1), xx(2), xx(3), sqrt(sum(xx(1:conf%dim)**2)), s%user_def)
    dl(i) = (l2 - l1)/(M_TWO*Delta)
  end do

end subroutine specie_get_dlocal
