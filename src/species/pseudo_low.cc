#include <sys/stat.h>
#include <algorithm>

#include "string_f.h" /* fortran <-> c string compatibility issues */

#include "base.hpp"
#include "qso.hpp"
#include "upf.hpp"


extern "C" void FC_FUNC_(pseudo_init, PSEUDO_INIT)(pseudopotential::base * pseudo, STR_F_TYPE const filename_f, int * ierr STR_ARG1){
  char * filename_c;
  TO_C_STR1(filename_f, filename_c);
  std::string filename(filename_c);
  free(filename_c);

  *ierr = 0;
  
  struct stat statbuf;
  bool found_file = !stat(filename.c_str(), &statbuf);

  if(!found_file) *ierr = 1;
  
  std::string extension = filename.substr(filename.find_last_of(".") + 1);
  std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
  
  std::cout << "  Opening file " << filename << std::endl;
  
  pseudo = NULL;
  
  if(extension == "xml"){
    pseudo = new pseudopotential::qso(filename);
  } else if(extension == "upf") {
    pseudo = new pseudopotential::upf(filename);
  } else {
    std::cerr << "Unknown pseudopotential type" << std::endl;
    exit(1);
  }

  assert(pseudo);
  
}

extern "C" void FC_FUNC_(pseudo_end, PSEUDO_END)(pseudopotential::base * pseudo){

  delete pseudo;
}

