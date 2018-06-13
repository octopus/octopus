/*
 Copyright (C) 2018 Xavier Andrade

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.
  
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
  
 You should have received a copy of the GNU Lesser General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <sys/stat.h>
#include <algorithm>

#include "string_f.h" /* fortran <-> c string compatibility issues */

#include "fortran_types.h"

#include "base.hpp"
#include "qso.hpp"
#include "upf1.hpp"
#include "upf2.hpp"
#include "psml.hpp"
#include "detect_format.hpp"

extern "C" fint FC_FUNC_(pseudo_detect_format, PSEUDO_DETECT_FORMAT)(STR_F_TYPE const filename_f STR_ARG1){

  char * filename_c;
  TO_C_STR1(filename_f, filename_c);
  
  pseudopotential::format ft = pseudopotential::detect_format(filename_c);

  free(filename_c);
  return fint(ft);
}

extern "C" void FC_FUNC_(pseudo_init, PSEUDO_INIT)(pseudopotential::base ** pseudo, STR_F_TYPE const filename_f, fint * format, fint * ierr STR_ARG1){
  char * filename_c;
  TO_C_STR1(filename_f, filename_c);
  std::string filename(filename_c);
  free(filename_c);

  *ierr = 0;
  
  struct stat statbuf;
  bool found_file = !stat(filename.c_str(), &statbuf);

  if(!found_file){
    *ierr = fint(pseudopotential::status::FILE_NOT_FOUND);
    return;
  }
  
  std::string extension = filename.substr(filename.find_last_of(".") + 1);
  std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
  
  *pseudo = NULL;

  try{
    switch(pseudopotential::format(*format)){
    case pseudopotential::format::QSO:
      *pseudo = new pseudopotential::qso(filename);
      break;
    case pseudopotential::format::UPF1:
      *pseudo = new pseudopotential::upf1(filename);
      break;
    case pseudopotential::format::UPF2:
      *pseudo = new pseudopotential::upf2(filename);
      break;
    case pseudopotential::format::PSML:
      *pseudo = new pseudopotential::psml(filename);
      break;
    default:
      *ierr = fint(pseudopotential::status::UNKNOWN_FORMAT);
      return;
    }
  } catch(pseudopotential::status stat) {
    *ierr = fint(stat);
    return;
  }

  assert(*pseudo);
  
}

extern "C" void FC_FUNC_(pseudo_end, PSEUDO_END)(pseudopotential::base ** pseudo){
  delete *pseudo;
}

extern "C" fint FC_FUNC_(pseudo_type, PSEUDO_TYPE)(const pseudopotential::base ** pseudo){
  return fint((*pseudo)->type());
}

extern "C" fint FC_FUNC_(pseudo_format, PSEUDO_FORMAT)(pseudopotential::base ** pseudo){
  return fint((*pseudo)->format());
}

extern "C" double FC_FUNC_(pseudo_valence_charge, PSEUDO_VALENCE_CHARGE)(const pseudopotential::base ** pseudo){
  return (*pseudo)->valence_charge();
}

extern "C" double FC_FUNC_(pseudo_mesh_spacing, PSEUDO_MESH_SPACING)(const pseudopotential::base ** pseudo){
  return (*pseudo)->mesh_spacing();
}

extern "C" fint FC_FUNC_(pseudo_mesh_size, PSEUDO_MESH_SIZE)(const pseudopotential::base ** pseudo){
  return (*pseudo)->mesh_size();
}

extern "C" double FC_FUNC_(pseudo_mass, PSEUDO_MASS)(const pseudopotential::base ** pseudo){
  return (*pseudo)->mass();
}

extern "C" fint FC_FUNC_(pseudo_lmax, PSEUDO_LMAX)(const pseudopotential::base ** pseudo){
  return (*pseudo)->lmax();
}

extern "C" fint FC_FUNC_(pseudo_llocal, PSEUDO_LLOCAL)(const pseudopotential::base ** pseudo){
  return (*pseudo)->llocal();
}

extern "C" fint FC_FUNC_(pseudo_nchannels, PSEUDO_NCHANNELS)(const pseudopotential::base ** pseudo){
  return (*pseudo)->nchannels();
}

extern "C" fint FC_FUNC_(pseudo_nprojectors, PSEUDO_NPROJECTORS)(const pseudopotential::base ** pseudo){
  return (*pseudo)->nprojectors();
}

extern "C" void FC_FUNC_(pseudo_grid, PSEUDO_GRID)(const pseudopotential::base ** pseudo, double * grid){
  std::vector<double> val;
  (*pseudo)->grid(val);
  for(unsigned i = 0; i < val.size(); i++) grid[i] = val[i];
}

extern "C" void FC_FUNC_(pseudo_grid_weights, PSEUDO_GRID_WEIGHT)(const pseudopotential::base ** pseudo, double * weight){
  std::vector<double> val;
  (*pseudo)->grid_weights(val);
  for(unsigned i = 0; i < val.size(); i++) weight[i] = val[i];
}

extern "C" void FC_FUNC_(pseudo_local_potential, PSEUDO_LOCAL_POTENTIAL)(const pseudopotential::base ** pseudo, double * local_potential){
  std::vector<double> val;
  (*pseudo)->local_potential(val);
  for(unsigned i = 0; i < val.size(); i++) local_potential[i] = val[i];
}

extern "C" void FC_FUNC_(pseudo_projector, PSEUDO_PROJECTOR)(const pseudopotential::base ** pseudo, fint * l, fint * ic, double * projector){
  std::vector<double> val;
  (*pseudo)->projector(*l, *ic - 1, val);
  for(unsigned i = 0; i < val.size(); i++) projector[i] = val[i];
}

extern "C" double FC_FUNC_(pseudo_dij, PSEUDO_DIJ)(const pseudopotential::base ** pseudo, fint * l, fint * ic, fint * jc){
  return (*pseudo)->d_ij(*l, *ic - 1, *jc - 1);
}

extern "C" void FC_FUNC_(pseudo_radial_potential, PSEUDO_RADIAL_POTENTIAL)(const pseudopotential::base ** pseudo, fint * l,double * radial_potential){
  std::vector<double> val;
  (*pseudo)->radial_potential(*l, val);
  for(unsigned i = 0; i < val.size(); i++) radial_potential[i] = val[i];
}

extern "C" fint FC_FUNC_(pseudo_has_radial_function_low, PSEUDO_HAS_RADIAL_FUNCTION_LOW)(const pseudopotential::base ** pseudo, fint * l){
  return fint((*pseudo)->has_radial_function(*l));
}

extern "C" void FC_FUNC_(pseudo_radial_function, PSEUDO_RADIAL_FUNCTION)(const pseudopotential::base ** pseudo, fint * l, double * radial_function){
  std::vector<double> val;
  (*pseudo)->radial_function(*l, val);
  for(unsigned i = 0; i < val.size(); i++) radial_function[i] = val[i];
}

extern "C" fint FC_FUNC_(pseudo_has_nlcc_low, PSEUDO_HAS_NLCC_LOW)(const pseudopotential::base ** pseudo){
  return fint((*pseudo)->has_nlcc());
}

extern "C" void FC_FUNC_(pseudo_nlcc_density, PSEUDO_NLCC_DENSITY)(const pseudopotential::base ** pseudo, double * nlcc_density){
  std::vector<double> val;
  (*pseudo)->nlcc_density(val);
  for(unsigned i = 0; i < val.size(); i++) nlcc_density[i] = val[i];
}

extern "C" fint FC_FUNC_(pseudo_has_density_low, PSEUDO_HAS_DENSITY_LOW)(const pseudopotential::base ** pseudo){
  return fint((*pseudo)->has_density());
}

extern "C" void FC_FUNC_(pseudo_density, PSEUDO_DENSITY)(const pseudopotential::base ** pseudo, double * density){
  std::vector<double> val;
  (*pseudo)->density(val);
  for(unsigned i = 0; i < val.size(); i++) density[i] = val[i];
}

extern "C" fint FC_FUNC_(pseudo_nwavefunctions, PSEUDO_NWAVEFUNCTIONS)(const pseudopotential::base ** pseudo){
  return (*pseudo)->nwavefunctions();
}

extern "C" void FC_FUNC_(pseudo_wavefunction, PSEUDO_WAVEFUNCTION)
  (const pseudopotential::base ** pseudo, const fint * index, fint * n, fint * l, double * occ, double * wavefunction){

  std::vector<double> val;
  (*pseudo)->wavefunction(*index - 1, *n, *l, *occ, val);
  for(unsigned i = 0; i < val.size(); i++) wavefunction[i] = val[i];
}

extern "C" fint FC_FUNC_(pseudo_exchange, PSEUDO_EXCHANGE)(const pseudopotential::base ** pseudo){
  return fint((*pseudo)->exchange());
}

extern "C" fint FC_FUNC_(pseudo_correlation, PSEUDO_CORRELATION)(const pseudopotential::base ** pseudo){
  return fint((*pseudo)->correlation());
}

extern "C" fint FC_FUNC_(pseudo_has_total_angular_momentum_low, PSEUDO_HAS_TOTAL_ANGULAR_MOMENTUM_LOW)(const pseudopotential::base ** pseudo){
  return fint((*pseudo)->has_total_angular_momentum());
}

extern "C" fint FC_FUNC_(pseudo_projector_2j, PSEUDO_PROJECTOR_2J)(const pseudopotential::base ** pseudo, fint * l, fint * ic){
  return (*pseudo)->projector_2j(*l, *ic - 1);
}

extern "C" fint FC_FUNC_(pseudo_wavefunction_2j, PSEUDO_WAVEFUNCTION_2J)(const pseudopotential::base ** pseudo, fint * ii){
  return (*pseudo)->wavefunction_2j(*ii);
}


