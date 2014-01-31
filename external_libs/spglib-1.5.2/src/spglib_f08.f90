
   module spglib_f08
      
   implicit none
   
   private

   interface 
   
   function spglib_get_symmetry( rotation, translation, max_size, lattice, &
                              & position, types, num_atom, symprec)
      integer, intent(inout) :: rotation
      real(8), intent(inout) :: translation
      integer, intent(in) :: max_size
      real(8), intent(in) :: lattice, position
      integer, intent(in) :: types
      integer, intent(in) :: num_atom
      real(8), intent(in) :: symprec
      integer :: spglib_get_symmetry
   end function spglib_get_symmetry
   
   function spglib_get_multiplicity( lattice, position, types, num_atom, symprec)
      real(8), intent(in) :: lattice, position
      integer, intent(in) :: types
      integer, intent(in) :: num_atom
      real(8), intent(in) :: symprec
      integer :: spglib_get_multiplicity
   end function spglib_get_multiplicity
   
   function spglib_get_international( symbol, lattice, position, types, num_atom, symprec)
      character(len=11), intent(out) :: symbol
      real(8), intent(in) :: lattice, position
      integer, intent(in) :: types
      integer, intent(in) :: num_atom
      real(8), intent(in) :: symprec
      integer :: spglib_get_international ! the number corresponding to 'symbol'. 0 on failure
   end function spglib_get_international
   
   function spglib_get_schoenflies( symbol, lattice, position, types, num_atom, symprec)
      character(len=10), intent(out) :: symbol
      real(8), intent(in) :: lattice, position
      integer, intent(in) :: types
      integer, intent(in) :: num_atom
      real(8), intent(in) :: symprec
      integer :: spglib_get_schoenflies ! the number corresponding to 'symbol'. 0 on failure
   end function spglib_get_schoenflies

   end interface


   public :: &
      & spglib_get_symmetry, spglib_get_multiplicity, &
      & spglib_get_international, spglib_get_schoenflies

   end module spglib_f08
