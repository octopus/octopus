/*
  Copyright (C) 2006 octopus team

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
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.

  $Id: io_csv.c 2146 2006-05-23 17:36:00Z xavier $
*/

/* 
The functions in this file read an array from an ascii matrix (csv) file.
Format with values "valueXYZ" as follows

File values.csv:
--------
value111 value112 value113
value121 value122 value123
value131 value132 value133

value211 value212 value213
value221 value222 value223
value231 value232 value233

value311 value312 value313
value321 value322 value323
value331 value332 value333
--------

That is, every XY-plane as a table of values and all XY-planes separated by an
empty row. 

The given matrix is interpolated/stretched to fit the calculation
box defined in input file. 

Calculation box shape must be "parallelepiped".

The delimiter can be a tab, a comma or a space.
*/

#include <config.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#ifdef HAVE_ERRNO_H
#include <errno.h>
#else
static int errno = -1;
#endif

#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif

#include "string_f.h"

#ifndef HAVE_UINT32_T
#if SIZEOF_UNSIGNED_INT == 4
typedef unsigned int uint32_t;
#elif SIZEOF_UNSIGNED_LONG == 4
typedef unsigned long uint32_t;
#else
#error no suitable 32-bit integer type found
#endif
#endif

#ifndef HAVE_UINT64_T
#if SIZEOF_UNSIGNED_LONG_LONG == 8
typedef unsigned long long uint64_t;
#elif SIZEOF_UNSIGNED_LONG == 8
typedef unsigned long uint64_t;
#else
#error no suitable 64-bit integer type found
#endif
#endif

#include "io_csv.h"

typedef char byte;

static const int size_of[6]      = {4, 8, 8, 16, 4, 8};

void FC_FUNC_(read_csv,READ_CSV)
     (unsigned long * np, byte * f, int * output_type, int * ierr, STR_F_TYPE fname STR_ARG1)
{
  char * filename;
  int i;
  FILE * fd;
  char * buf;
  char * c;
  const char sep [] = "\t\n ,";
  
  int buf_size = 65536;
  buf = (char *) malloc(buf_size*sizeof(char));
  assert(buf != NULL);
    
  TO_C_STR1(fname, filename);
  fd = fopen(filename, "r");
  if(fd == NULL) {
    *ierr = 2;
    return;
  }
  
  free(filename);
  
  if ( (*output_type) == TYPE_FLOAT ) {
    i = 0;
    while(fgets(buf,buf_size*sizeof(char),fd) != NULL){
      float d;
      c = strtok(buf, sep);
      while (c != NULL) {
        assert(i/8 < *np);
        d = strtof(c, (char **) NULL);
        c = (char *) strtok((char *) NULL, sep);
        memcpy(f + i, &d, size_of[(*output_type)]);
        i += size_of[(*output_type)];
      }
    }
  }
  else if ( (*output_type) == TYPE_DOUBLE ) {
    double d;
    i = 0;
    while(fgets(buf,buf_size*sizeof(char),fd) != NULL)
    {
      c = strtok(buf, sep);
      while(c != NULL) {
        assert(i/8 < *np);
        d = strtod(c, (char **) NULL);
        memcpy(f + i, &d, size_of[(*output_type)]);
        c = (char *) strtok((char *) NULL, sep);
        i += size_of[(*output_type)];
      }
    }
  }
  
  free(buf);
  fclose(fd);
}

void FC_FUNC_(get_info_csv,GET_INFO_CSV)
     (unsigned long * dims, int * ierr, STR_F_TYPE fname STR_ARG1)
{
  char * filename;
  char * buf;
  char * c;
  FILE * fd;
  int buf_size = 65536;
  const char sep [] = "\n\t ,";
  
  unsigned long curr_dims [3] = {0, 0, 0};
  unsigned long prev_dims [2] = {0, 0};
  
  TO_C_STR1(fname, filename);
  fd = fopen(filename, "r");
  if (fd == NULL) {
    *ierr = 2;
    return;
  }
  free(filename);
  
  buf = (char *) malloc(buf_size*sizeof(char));
  assert(buf != NULL);
  while(fgets(buf,buf_size*sizeof(char),fd) != NULL) {
    c = strtok(buf, sep);
    
    prev_dims[0] = curr_dims[0];
    
    curr_dims[0] = 0;
    
    /** count the number of columns i.e. the size in x-direction**/
    while(c != NULL) {
      curr_dims[0]++;
      c = (char *) strtok((char *) NULL, sep);
    }

    /** The number of columns must be the same on all non-empty rows **/
    /** This only checks that the number of columns is correct
        within each z-block **/
    if(prev_dims[0] > 0 && curr_dims[0] > 0)
      assert(curr_dims[0] == prev_dims[0]); 
    
    /** If the previous line was empty and the current one is non-empty, then it signifies
        the start of a new z block i.e. a new xy-plane **/
    if (prev_dims[0] == 0 && curr_dims[0] != 0) {
      prev_dims[1] = curr_dims[1];
      curr_dims[1] = 0;
      curr_dims[2]++;

    /** The number of rows i.e. the size in y-direction must be the same within all z-blocks **/
      if (prev_dims[1] > 0 && curr_dims[1] > 0)
        assert(prev_dims[1] == curr_dims[1]);
    }
    
    /** If the current row is non-empty, increase the y-dimension **/

    if (curr_dims[0] > 0)
      curr_dims[1]++;
  }
  
  dims[0] = curr_dims[0];
  dims[1] = curr_dims[1];
  dims[2] = curr_dims[2];
  
  free(buf);
  fclose(fd);
}

