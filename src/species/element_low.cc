#include "element.hpp"

#include "string_f.h" /* fortran <-> c string compatibility issues */

#include "fortran_types.h"

extern "C" void FC_FUNC_(element_init, ELEMENT_INIT)(pseudopotential::element ** el, STR_F_TYPE const symbol_f STR_ARG1){
  char * symbol_c;
  TO_C_STR1(symbol_f, symbol_c);

  *el = new pseudopotential::element(symbol_c);
  
  free(symbol_c);
}

extern "C" void FC_FUNC_(element_end, ELEMENT_END)(pseudopotential::element ** el){
  delete *el;
}

extern "C" double FC_FUNC_(element_mass, ELEMENT_MASS)(pseudopotential::element ** el){
  return (*el)->mass();
}

extern "C" double FC_FUNC_(element_vdw_radius, ELEMENT_VDW_RADIUS)(pseudopotential::element ** el){
  return (*el)->vdw_radius();
}

extern "C" fint FC_FUNC_(element_atomic_number, ELEMENT_ATOMIC_NUMBER)(pseudopotential::element ** el){
  return (*el)->atomic_number();
}

extern "C" fint FC_FUNC_(element_valid_low, ELEMENT_VALID_LOW)(pseudopotential::element ** el){
  return fint((*el)->valid());
}
