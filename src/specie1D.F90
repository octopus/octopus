subroutine specie1D_init(nspecies, str, s)
  integer, intent(in) :: nspecies
  character(len=*), intent(in) :: str
  type(specie_type), pointer :: s(:)

  integer :: i, j

  do i = 1, nspecies
    s(i)%index = i
    s(i)%local = .true.
    s(i)%nlcc  = .false.

    call loct_parse_block_string(str, i-1, 0, s(i)%label)
    call loct_parse_block_float (str, i-1, 1, s(i)%weight)
    s(i)%weight =  units_inp%mass%factor * s(i)%weight ! units conversion

    s(i)%local = .true. ! In 1D or 2D, potential has to be local.
    call loct_parse_block_float (str, i-1, 2, s(i)%Z_val)
    call loct_parse_block_string(str, i-1, 3, s(i)%user_def)
    ! convert to C string
    j = len(trim(s(i)%user_def))
    s(i)%user_def(j+1:j+1) = achar(0) 
  end do

end subroutine specie1D_init
