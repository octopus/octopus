module pulpo
  use global
  use liboct

contains
  subroutine pulpo_print()
    character(len=2) :: lang
    
    ! some white space
    message(1) = ''; message(2) = ''
    call write_info(2)

    call oct_parse_str('RecipeLang', 'en', lang)
    call lowcase(lang)
    select case(lang)
    case('en')
      call oct_printRecipe(1)
    case('es')
      call oct_printRecipe(2)
    case default
      write(message(1), '(a,a,a)') "Language '", trim(lang), "' is not recognized"
      message(2) = "  RecipeLang = en | es"
      call write_fatal(2)
    end select

    ! DISCLAIMER
    call oct_printRecipe(0)

  end subroutine pulpo_print
end module pulpo
