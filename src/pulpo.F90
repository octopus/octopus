module pulpo
  use global
  use fdf

contains
  subroutine pulpo_print()
    character(len=2) :: lang
    
    ! some white space
    message(1) = ''; message(2) = ''
    call write_info(2)

    lang = fdf_string('RecipeLang', 'en')
    call lowcase(lang)
    select case(lang)
    case('en')
      call printRecipe(1)
    case('es')
      call printRecipe(2)
    case default
      write(message(1), '(a,a,a)') "Language '", trim(lang), "' is not recognized"
      message(2) = "  RecipeLang = en | es"
      call write_fatal(2)
    end select

    ! DISCLAIMER
    call printRecipe(0)

  end subroutine pulpo_print
end module pulpo
