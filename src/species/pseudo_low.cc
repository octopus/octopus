#include "string_f.h" /* fortran <-> c string compatibility issues */

#include "base.hpp"
#include "qso.hpp"
#include "upf.hpp"

extern "C" void pseudo_init(pseudopotential::base * pseudo, STR_F_TYPE const filename_f STR_ARG1){
  char * filename;
  TO_C_STR1(filename_f, filename);


  free(filename);
  
}

