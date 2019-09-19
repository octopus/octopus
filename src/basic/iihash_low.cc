/*
 Copyright (C) 2019 X. Andrade

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.

*/

#include <config.h>

#include <unordered_map>

#include <fortran_types.h>

typedef std::unordered_map<fint, fint> map_type;

extern "C" void FC_FUNC_(iihash_map_init, IIHASH_MAP_INIT)(map_type ** map){
  *map = new map_type;
}


extern "C" void FC_FUNC_(iihash_map_end, IIHASH_MAP_END)(map_type ** map){
  delete *map;
}


extern "C" void FC_FUNC_(iihash_map_insert, IIHASH_MAP_INSERT)(map_type ** map, const fint * key, const fint * val){
  (**map)[*key] = *val;
}

extern "C" void FC_FUNC_(iihash_map_lookup, IIHASH_MAP_LOOKUP)(const map_type ** map, const fint * key, fint * found, fint * val){

  auto it = (*map)->find(*key);

  if(it == (*map)->end()){
    *found = 0;
  } else {
    *found = 1;
    *val = it->second;
  }
  
}

