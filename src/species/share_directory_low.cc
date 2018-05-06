#include "share_directory.hpp"

#include "string_f.h" /* fortran <-> c string compatibility issues */

#include "fortran_types.h"

extern "C" void FC_FUNC_(share_directory_set, SHARE_DIRECTORY_SET)(STR_F_TYPE dir_f STR_ARG1){
  char * dir_c;
  TO_C_STR1(dir_f, dir_c);
  pseudopotential::share_directory::set(dir_c);
  free(dir_c);
}
