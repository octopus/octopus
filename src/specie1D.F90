subroutine specie1D_init(s, location, line)
  type(specie_type), intent(inout) :: s
  integer,           intent(in)    :: location, line

  integer :: j

  if(location.ne.1) then
    message(1) = "In 1 or 2 dimensions all species must be defined in the %Species block"
    call write_fatal(1)
  end if

  call loct_parse_block_float ("Species", line-1, 1, s%weight)
  s%weight =  units_inp%mass%factor * s%weight ! units conversion

  call loct_parse_block_float ("Species", line-1, 2, s%Z_val)
  call loct_parse_block_string("Species", line-1, 3, s%user_def)

  ! convert to C string
  j = len(trim(s%user_def))
  s%user_def(j+1:j+1) = achar(0) 

end subroutine specie1D_init
