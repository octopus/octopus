/*
 Copyright (C) 2016 X. Andrade

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

 $Id$
*/

#include <fortran_types.h>
#include <algorithm>
#include <new>

//Functor object that compares indices based on an array
template <typename TT>
class compare{

public:

  compare(const TT * array_val){
    array = array_val;
  }
  
  TT operator()(const fint ii, const fint jj){
    return array[ii] < array[jj];
  }
  
private:
  
  const TT * array;

};

//Worker function that sorts an array and returns the order
template <typename TT>
void sort2(const fint size, TT * array, fint * indices){

  // first sort the indices
  
  for(fint ii = 0; ii < size; ii++){
    indices[ii] = ii;
  }

  std::sort(indices, indices + size, compare<TT>(array));  


  // now sort the array

  TT * array_copy = new TT[size];
  
  std::copy(array, array + size, array_copy);

  for(fint ii = 0; ii < size; ii++){
    array[ii] = array_copy[indices[ii]];
    indices[ii]++; //convert indices to fortran convention
  }

  delete [] array_copy;
}

//Fortran interfaces

extern "C" void FC_FUNC(isort1, ISORT1)(const fint * size, fint * array){
  std::sort(array, array + *size);
}

extern "C" void FC_FUNC(isort2, ISORT2)(const fint * size, fint * array, fint * indices){
  sort2<fint>(*size, array, indices);
}

extern "C" void FC_FUNC(dsort1, DSORT1)(const fint * size, double * array){
  std::sort(array, array + *size);
}

extern "C" void FC_FUNC(dsort2, DSORT2)(const fint * size, double * array, fint * indices){
  sort2<double>(*size, array, indices);
}
