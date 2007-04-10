/*
  Copyright (C) 2006 X. Andrade

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

  $Id: recipes.c 2146 2006-05-23 17:36:00Z xavier $
*/

/* 

The functions in this file write and read an array to binary file.

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
#error no suitable 32 bits integer type found
#endif
#endif

#ifndef HAVE_UINT64_T
#if SIZEOF_UNSIGNED_LONG_LONG == 8
typedef unsigned long long uint64_t;
#elif SIZEOF_UNSIGNED_LONG == 8
typedef unsigned long uint64_t;
#else
#error no suitable 64 bits integer type found
#endif
#endif


typedef char byte;

#define type_float 0
#define type_double 1
#define type_float_complex 2
#define type_double_complex 3

static const int size_of[4]={4, 8, 8, 16};

inline void inf_error(const char * msg, int * ierr){
#ifdef HAVE_PERROR
  perror(msg);
#else
  printf(msg);
  printf(": I/O Error.\n");
#endif
}

typedef union {
  float f[2];
  double d[2];
} multi;


static void convert ( multi * in, multi * out, int t_in, int t_out);

/* how to convert a complex to a real */ 
#define c2r hypot /* take the modulus */


/* THE HEADER OF THE FILE */
typedef struct {
  /* text to identify the file */
  char text[7];

  /* version of the format */
  uint8_t version;

  /* value of 1 in different formats, to recognize endianness */
  uint32_t one_32;
  float one_f;
  uint64_t one_64;
  double one_d;

  /* the size of the array stored */
  uint64_t np;

  /* the type of the wfs */
  uint32_t type;

  /* extra values for future versions*/
  uint32_t extra[5];

} header_t;

inline void init_header(header_t * h){
  int i;

  strcpy(h -> text, "pulpo");
  h -> version = 0;
  h -> one_32  = 1;
  h -> one_f   = 1.0;
  h -> one_64  = 1;
  h -> one_d   = 1.0;
  for(i=0;i<5;i++) h -> extra[i] = 0;
}

inline int check_header(header_t * h, int * correct_endianness){
  if( strcmp("pulpo", h -> text) != 0 ) return 5;
  if( h -> version != 0 ) return 5;
  return 0;
}


void FC_FUNC_(write_binary, WRITE_BINARY)
     (int * np, void * f, int * type, int * ierr, STR_F_TYPE fname STR_ARG1)
{
  header_t h;
  char * filename;
  int fd;
  ssize_t moved;

  *ierr = 0;

  filename = TO_C_STR1(fname);
  fd = open (filename, O_WRONLY | O_CREAT | O_TRUNC, 
	     S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH );
  if( fd < 0 ) {
    inf_error("octopus.write_binary", ierr);
    *ierr = 2;
    return;
  }
  free(filename);

  /* create header */
  init_header(&h);
  h.np = *np;
  h.type = *type;

  /* write header */
  moved = write(fd, &h, sizeof(header_t));

  if(moved < sizeof(header_t)){
    inf_error("octopus.write_binary", ierr);
    close(fd);
    return;
  }

  /* now write the values */
  moved = write(fd, f, (*np)*size_of[(*type)]);

  if(moved < (*np)*size_of[(*type)]){
    inf_error("octopus.write_binary", ierr);
  }

  /* close the file */
  close(fd);

}

void FC_FUNC_(read_binary,READ_BINARY)
     (int * np, byte * f, int * output_type, int * ierr, STR_F_TYPE fname STR_ARG1)
{
  header_t h;
  char * filename;
  int fd, i;
  ssize_t moved;
  int correct_endianness;
  byte * read_f;

  filename = TO_C_STR1(fname);
  fd = open(filename, O_RDONLY);
  if(fd < 0){
    inf_error("octopus.read_binary", ierr);
    *ierr = 2;
    return;
  }

  free(filename);

  /* read header */
  moved = read(fd, &h, sizeof(header_t));
  if ( moved != sizeof(header_t) ) { 
    /* we couldn't read the complete header */
    *ierr = 3;
    return;
  }

  *ierr = check_header(&h, &correct_endianness);
  if( *ierr != 0 ) return;

  /* check whether the sizes match */ 
  if( h.np != *np ){ *ierr = 4; return; }

  if( h.type == *output_type){
    /* format is the same, we just read */
    read_f = f;
  } else {
    /*format is not the same, we store into a temporary array */
    read_f =(byte *) malloc((*np)*size_of[h.type]);
  }

  /* now read the values and close the file */
  moved = read(fd, read_f, (*np)*size_of[h.type]);

  if ( moved != (*np)*size_of[h.type]) { 
    /* we couldn't read the whole dataset */
    *ierr = 3;
    return;
  }
    
  close(fd);

  /* convert values if it is necessary */
  if( h.type != *output_type ){
    for(i=0; i < *np ; i++) 
      convert( (multi *) (read_f + i*size_of[h.type]), 
	       (multi *) (f + i*size_of[*output_type]), 
	       h.type, *output_type);
    free(read_f);

    /* set the error code according to the conversion done (see src/out_inc.F90 ) */
    if ( h.type == type_float )          *ierr = -1;
    if ( h.type == type_float_complex )  *ierr = -2;
    if ( h.type == type_double )         *ierr = -3;
    if ( h.type == type_double_complex ) *ierr = -4;

  }
  
}

/*
  This function converts between types.
*/

static void convert ( multi * in, multi * out, int t_in, int t_out){

  /* real types */
  if(t_in == type_float && t_out == type_double ) {
    out->d[0] = in->f[0]; return;
  }

  if(t_in == type_double && t_out == type_float ){
    out->f[0] = in->d[0]; return;
  }

  /* complex types */
  if(t_in == type_float_complex && t_out == type_double_complex ){
    out->d[0] = in->f[0]; 
    out->d[1] = in->f[1]; 
    return;
  }
  if(t_in == type_double_complex && t_out == type_float_complex ){
    out->f[0] = in->d[0]; 
    out->f[1] = in->d[1]; 
    return;
  }

  /* real to complex */
  if(t_in == type_float && t_out == type_float_complex ){
    out->f[0] = in->f[0]; 
    out->f[1] = (float) 0.0;
    return;
  }

  if(t_in == type_double && t_out == type_float_complex ){
    out->f[0] = in->d[0]; 
    out->f[1] = (float) 0.0;
    return;
  }

  if(t_in == type_float && t_out == type_double_complex ){
    out->d[0] = in->f[0]; 
    out->d[1] = (double) 0.0;
    return;
  }

  if(t_in == type_double && t_out == type_double_complex ){
    out->d[0] = in->d[0]; 
    out->d[1] = (double) 0.0;
    return;
  }

  /* complex to real */
  if(t_in == type_float_complex && t_out == type_float ){
    out->f[0] = c2r(in->f[0], in->f[1]);
    return;
  }

  if(t_in == type_double_complex && t_out == type_float ){
    out->f[0] = c2r(in->d[0], in->d[1]);
    return;
  }

  if(t_in == type_float_complex && t_out == type_double ){
    out->d[0] = c2r(in->f[0], in->f[1]);
    return;
  }

  if(t_in == type_double_complex && t_out == type_double ){
    out->d[0] = c2r(in->d[0], in->d[1]);
    return;
  }

}
