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

#include <fortran_types.h>
#include <unordered_map>
#include <cassert>
#include <iostream>

struct alloc_cache {
  typedef std::unordered_multimap<fint8, void *> map;
  map list;
  fint8 max_size;
  fint8 current_size;
  fint8 hits;
  fint8 misses;
  double vol_hits;
  double vol_misses;

  alloc_cache(fint8 arg_max_size){
    max_size = arg_max_size;
    current_size = 0;
    hits = 0;
    misses = 0;
    vol_hits = 0.0;
    vol_misses = 0.0;
  }

};

extern "C" void FC_FUNC(alloc_cache_init, ALLOC_CACHE_INIT)(alloc_cache ** cache, fint8 max_size){
  *cache = new alloc_cache(max_size);
}

extern "C" void FC_FUNC(alloc_cache_end, ALLOC_CACHE_END)( alloc_cache ** cache, fint8 * hits, fint8 * misses,
							   double * vol_hits, double * vol_misses){

  assert((*cache)->list.empty());

  delete *cache;
  *hits = (*cache)->hits;
  *misses = (*cache)->misses;
  *vol_hits = (*cache)->vol_hits++;
  *vol_misses = (*cache)->vol_misses++;
}

extern "C" void FC_FUNC(alloc_cache_put_low, ALLOC_CACHE_PUT_LOW)(alloc_cache ** cache, const fint8 * size, void ** loc, fint * put){
  if( (*cache)->current_size + *size <=  (*cache)->max_size ){
    (*cache)->current_size += *size;
    (*cache)->list.insert(alloc_cache::map::value_type(*size, *loc));
    *put = 1;
  } else {
    *put = 0;
  }

}

#define CACHE_ALLOC_ANY_SIZE -1

extern "C" void FC_FUNC(alloc_cache_get_low, ALLOC_CACHE_GET_LOW)(alloc_cache ** cache, const fint8 * size, fint * found, void ** loc){
  if(*size == CACHE_ALLOC_ANY_SIZE){
    auto pos = (*cache)->list.begin();
    *found = pos != (*cache)->list.end();

    if(*found) {
      *loc = pos->second;
      (*cache)->list.erase(pos);
    }
    
  } else {
    auto pos = (*cache)->list.find(*size);
    *found = (pos != (*cache)->list.end());
    
    if(*found){
      (*cache)->current_size -= *size;
      
      assert(pos->first == *size);
      *loc = pos->second;
      (*cache)->list.erase(pos);
      (*cache)->hits++;
      (*cache)->vol_hits += *size;
    } else {
      *loc = NULL;
      (*cache)->misses++;
      (*cache)->vol_misses += *size;
    }

  }
  
}

