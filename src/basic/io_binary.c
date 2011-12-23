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

#include "io_binary.h"

typedef char byte;

static const int size_of[6]      = {4, 8, 8, 16, 4, 8};
static const int base_size_of[6] = {4, 8, 4, 8, 4, 8};
static const int is_integer[6]   = {0, 0, 0, 0, 1, 1};

static inline void inf_error(const char * msg){
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

/* A very basic endian conversion routine. This can be improved a lot,
   but it is not necessary, as restart files are converted at most
   once per run */

static inline void endian_convert (const int size, char * a){
  char tmp[8];
  int ii;

  for(ii = 0; ii < size; ii++) tmp[ii] = a[ii];
  for(ii = 0; ii < size; ii++) a[ii] = tmp[size - 1 - ii];
}

static void convert ( multi * in, multi * out, int t_in, int t_out);

/* how to convert a complex to a real */ 
#define c2r hypot /* take the modulus */


/* THE HEADER OF THE FILE */
typedef struct {
  /* text to identify the file */
  char text[7]; /*7 bytes*/

  /* version of the format */
  uint8_t version; /* 8 bytes*/

  /* value of 1 in different formats, to recognize endianness */
  uint32_t one_32; /*12 bytes */
  float one_f; /* 16 bytes*/
  uint64_t one_64; /* 24 bytes */
  double one_d; /* 32 bytes */

  /* the size of the array stored */
  uint64_t np; /* 40 bytes */ 

  /* the type of the wfs */
  uint32_t type; /* 44 bytes */

  /* extra values for future versions*/
  uint32_t extra[5]; /* 64 bytes */

} header_t;

static inline void init_header(header_t * h){
  int i;
  strcpy(h -> text, "pulpo");
  h -> text[6] = 0;
  h -> version = 0;
  h -> one_32  = 1;
  h -> one_f   = 1.0;
  h -> one_64  = 1;
  h -> one_d   = 1.0;
  for(i=0;i<5;i++) h -> extra[i] = 0;
}

static inline int check_header(header_t * h, int * correct_endianness){
  if( strcmp("pulpo", h -> text) != 0 ) return 5;
  if( h -> version != 0 ) return 5;

  /* Check the endianness of integer values and fix header
     components */

  if(h -> one_32 != 1) {
    endian_convert(4, (char *) &(h -> one_32));
    if( h -> one_32 != 1 ) return 5;
    endian_convert(4, (char *) &(h -> type));
  }

  if(h -> one_64 != 1) {
    endian_convert(8, (char *) &(h -> one_64));
    if( h -> one_64 != 1 ) return 5;
    endian_convert(8, (char *) &(h -> np));
  }

  /* Check the endianness of floating point values  */
  *correct_endianness = 0;
  if ( base_size_of[h -> type] == 4 ) {
    if(h -> one_f != 1.0) {
      endian_convert(4, (char *) &(h -> one_f));
      if(h -> one_f != 1.0) return 5;
      *correct_endianness = 1;
    }
  } else {
    if(h -> one_d != 1.0) {
      endian_convert(8, (char *) &(h -> one_d));
      if(h -> one_d != 1.0) return 5;
      *correct_endianness = 1;
    }
  }

  return 0;
}


void FC_FUNC_(write_binary,WRITE_BINARY)
     (const int * np, void * f, int * type, int * ierr, STR_F_TYPE fname STR_ARG1)
{
  header_t * h;
  char * filename;
  int fd;
  ssize_t moved;

  h = (header_t *) malloc(sizeof(header_t));
  assert(h != NULL);

  *ierr = 0;

  TO_C_STR1(fname, filename);
  fd = open (filename, O_WRONLY | O_CREAT | O_TRUNC, 
	     S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH );
  free(filename);
  if( fd < 0 ) {
    inf_error("octopus.write_binary");
    *ierr = 2;
    free(h);
    return;
  }

  /* create header */
  init_header(h);
  h->np = *np;
  h->type = *type;

  /* write header */
  moved = write(fd, h, sizeof(header_t));

  if(moved < sizeof(header_t)){
    inf_error("octopus.write_binary");
    close(fd);
    free(h);
    return;
  }

  /* now write the values */
  moved = write(fd, f, (*np)*size_of[(*type)]);

  if(moved < (*np)*size_of[(*type)]){
    inf_error("octopus.write_binary");
  }

  /* close the file */
  close(fd);
  free(h);

}
 
void FC_FUNC_(read_binary,READ_BINARY)
     (const int * np, const int * offset, byte * f, int * output_type, int * ierr, STR_F_TYPE fname STR_ARG1)
{
  header_t * h;
  char * filename;
  int fd, i;
  ssize_t moved;
  int correct_endianness;
  byte * read_f;

  TO_C_STR1(fname, filename);
  fd = open(filename, O_RDONLY);
  free(filename);
  if(fd < 0){
    *ierr = 2;
    return;
  }

  h = (header_t *) malloc(sizeof(header_t));
  assert(h != NULL);

  /* read header */
  moved = read(fd, h, sizeof(header_t));
  if ( moved != sizeof(header_t) ) { 
    /* we couldn't read the complete header */
    *ierr = 3;
    return;
  }

  *ierr = check_header(h, &correct_endianness);
  if( *ierr != 0 ) return;

  /* check whether the sizes match */ 
  if( h->np < *np + *offset ){ *ierr = 4; return; }

  if( h->type == *output_type){
    /* format is the same, we just read */
    read_f = f;
  } else {
    /*format is not the same, we store into a temporary array */
    read_f =(byte *) malloc((*np)*size_of[h->type]);
  }

  /* set the start point */
  if(*offset != 0) lseek(fd, (*offset)*size_of[h->type], SEEK_CUR);

  /* now read the values and close the file */
  moved = read(fd, read_f, (*np)*size_of[h->type]);

  if ( moved != (*np)*size_of[h->type]) { 
    /* we couldn't read the whole dataset */
    *ierr = 3;
    return;
  }
    
  close(fd);
  
  /* convert endianness */
  
  if(correct_endianness) {
    for(i=0; i < (*np)*size_of[h->type] ; i+=base_size_of[h->type]) 
      endian_convert(base_size_of[h->type], (char *) (read_f + i));
  }

  /* convert values if it is necessary */
  if( h->type != *output_type ){
    
    if(is_integer[h->type] || is_integer[*output_type]){
      *ierr = 5;
    } else {

      for(i=0; i < *np ; i++) 
	convert( (multi *) (read_f + i*size_of[h->type]), 
		 (multi *) (f + i*size_of[*output_type]), 
		 h->type, *output_type);
      free(read_f);

      /* set the error code according to the conversion done (see src/out_inc.F90 ) */
      if ( h->type == TYPE_FLOAT )          *ierr = -1;
      if ( h->type == TYPE_FLOAT_COMPLEX )  *ierr = -2;
      if ( h->type == TYPE_DOUBLE )         *ierr = -3;
      if ( h->type == TYPE_DOUBLE_COMPLEX ) *ierr = -4;
    }
  }
  
  free(h);
}

/*
  This function converts between types.
*/

static void convert ( multi * in, multi * out, int t_in, int t_out){

  /* real types */
  if(t_in == TYPE_FLOAT && t_out == TYPE_DOUBLE ) {
    out->d[0] = in->f[0]; return;
  }

  if(t_in == TYPE_DOUBLE && t_out == TYPE_FLOAT ){
    out->f[0] = in->d[0]; return;
  }

  /* complex types */
  if(t_in == TYPE_FLOAT_COMPLEX && t_out == TYPE_DOUBLE_COMPLEX ){
    out->d[0] = in->f[0]; 
    out->d[1] = in->f[1]; 
    return;
  }
  if(t_in == TYPE_DOUBLE_COMPLEX && t_out == TYPE_FLOAT_COMPLEX ){
    out->f[0] = in->d[0]; 
    out->f[1] = in->d[1]; 
    return;
  }

  /* real to complex */
  if(t_in == TYPE_FLOAT && t_out == TYPE_FLOAT_COMPLEX ){
    out->f[0] = in->f[0]; 
    out->f[1] = (float) 0.0;
    return;
  }

  if(t_in == TYPE_DOUBLE && t_out == TYPE_FLOAT_COMPLEX ){
    out->f[0] = in->d[0]; 
    out->f[1] = (float) 0.0;
    return;
  }

  if(t_in == TYPE_FLOAT && t_out == TYPE_DOUBLE_COMPLEX ){
    out->d[0] = in->f[0]; 
    out->d[1] = (double) 0.0;
    return;
  }

  if(t_in == TYPE_DOUBLE && t_out == TYPE_DOUBLE_COMPLEX ){
    out->d[0] = in->d[0]; 
    out->d[1] = (double) 0.0;
    return;
  }

  /* complex to real */
  if(t_in == TYPE_FLOAT_COMPLEX && t_out == TYPE_FLOAT ){
    out->f[0] = c2r(in->f[0], in->f[1]);
    return;
  }

  if(t_in == TYPE_DOUBLE_COMPLEX && t_out == TYPE_FLOAT ){
    out->f[0] = c2r(in->d[0], in->d[1]);
    return;
  }

  if(t_in == TYPE_FLOAT_COMPLEX && t_out == TYPE_DOUBLE ){
    out->d[0] = c2r(in->f[0], in->f[1]);
    return;
  }

  if(t_in == TYPE_DOUBLE_COMPLEX && t_out == TYPE_DOUBLE ){
    out->d[0] = c2r(in->d[0], in->d[1]);
    return;
  }

}


void FC_FUNC_(get_info_binary,GET_INFO_BINARY)
     (int * np, int * type, int * ierr, STR_F_TYPE fname STR_ARG1)
{
  header_t * h;
  char * filename;
  int fd;
  ssize_t moved;
  int correct_endianness;

  TO_C_STR1(fname, filename);
  fd = open(filename, O_RDONLY);
  if(fd < 0){
    *ierr = 2;
    return;
  }

  free(filename);

  h = (header_t *) malloc(sizeof(header_t));
  assert(h != NULL);

  /* read header */
  moved = read(fd, h, sizeof(header_t));

  close(fd);

  if ( moved != sizeof(header_t) ) { 
    /* we couldn't read the complete header */
    *ierr = 3;
    return;
  }

  *ierr = check_header(h, &correct_endianness);
  if( *ierr != 0 ) return;

  *np  = h->np;
  *type = (int) h->type;

  free(h);
}
