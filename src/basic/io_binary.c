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
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
  02110-1301, USA.

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
#include <fortran_types.h>

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

static inline void endian_convert (const int size, char * aa){
  char tmp[8];
  int ii;

  for(ii = 0; ii < size; ii++) tmp[ii] = aa[ii];
  for(ii = 0; ii < size; ii++) aa[ii] = tmp[size - 1 - ii];
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

static inline void init_header(header_t * hp){
  int ii;

  strcpy(hp -> text, "pulpo");

  hp -> text[6] = 0;
  hp -> version = 0;
  hp -> one_32  = 1;
  hp -> one_f   = 1.0;
  hp -> one_64  = 1;
  hp -> one_d   = 1.0;

  for(ii=0;ii<5;ii++) hp -> extra[ii] = 0;
}

static inline int check_header(header_t * hp, int * correct_endianness){
  if( strcmp("pulpo", hp -> text) != 0 ) return 5;
  if( hp -> version != 0 ) return 5;

  /* Check the endianness of integer values and fix header
     components */

  if(hp -> one_32 != 1) {
    endian_convert(4, (char *) &(hp -> one_32));
    if( hp -> one_32 != 1 ) return 5;
    endian_convert(4, (char *) &(hp -> type));
  }

  if(hp -> one_64 != 1) {
    endian_convert(8, (char *) &(hp -> one_64));
    if( hp -> one_64 != 1 ) return 5;
    endian_convert(8, (char *) &(hp -> np));
  }

  /* Check the endianness of floating point values  */
  *correct_endianness = 0;
  if ( base_size_of[hp -> type] == 4 ) {
    if(hp -> one_f != 1.0) {
      endian_convert(4, (char *) &(hp -> one_f));
      if(hp -> one_f != 1.0) return 5;
      *correct_endianness = 1;
    }
  } else {
    if(hp -> one_d != 1.0) {
      endian_convert(8, (char *) &(hp -> one_d));
      if(hp -> one_d != 1.0) return 5;
      *correct_endianness = 1;
    }
  }

  return 0;
}

void io_write_header(const fint * np, fint * type, fint * ierr, fint * iio, STR_F_TYPE fname STR_ARG1)
{
  char * filename;
  header_t * hp;
  int fd;
  ssize_t moved;

  hp = (header_t *) malloc(sizeof(header_t));
  assert(hp != NULL);
  assert(np > 0);

  *ierr = 0;
  TO_C_STR1(fname, filename);
  fd = open (filename, O_WRONLY | O_CREAT | O_TRUNC, 
	     S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH );
  *iio+=100;
  if( fd < 0 ) {
    printf("Filename is %s\n",filename);
    inf_error("octopus.write_header in creating the header");
    *ierr = 2;
    free(hp);
    return;
  }

  /* create header */
  init_header(hp);
  hp->np = *np;
  hp->type = *type;

  /* write header */
  moved = write(fd, hp, sizeof(header_t));

  if(moved < sizeof(header_t)){
    /* we couldn't write the complete header */
    inf_error("octopus.write_header in writing the header");
    *ierr = 3;
  }

  free(hp);
  free(filename);
  close(fd);
  *iio += 1;
}

void FC_FUNC_(write_header,WRITE_HEADER)(const fint * np, fint * type, fint * ierr, fint * iio, STR_F_TYPE fname STR_ARG1)
{ 
  unsigned long fname_len;
  fname_len = l1;
  io_write_header(np, type, ierr, iio, fname, fname_len);
}

void FC_FUNC_(write_binary,WRITE_BINARY)
     (const fint * np, void * ff, fint * type, fint * ierr, fint * iio, fint * nhd, fint * flpe, STR_F_TYPE fname STR_ARG1)
{
  char * filename;
  int fd, ii;
  ssize_t moved;
  unsigned long fname_len;

  assert(np > 0);
  *ierr = 0;
  
  fname_len = l1;
  TO_C_STR1(fname, filename);

  if(*nhd != 1){
    io_write_header(np, type, ierr, iio, fname, fname_len);
  }
  
  fd = open (filename, O_WRONLY, 
	     S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH );    
  iio += 100; 
  free(filename);
  if( fd < 0 ) {
    inf_error("octopus.write_binary in opening the file");
    *ierr = 2;
    return;
  }

  /* skip the header and go until the end */
  lseek(fd, 0, SEEK_END);
  
  /* flip endianness*/
  if (*flpe == 1){
    for(ii=0; ii < (*np)*size_of[(*type)] ; ii+=base_size_of[(*type)]) 
      endian_convert(base_size_of[(*type)], (char *) (ff + ii));
  }
  
  /* now write the values */
  moved = write(fd, ff, (*np)*size_of[(*type)]);

  if(moved < (*np)*size_of[(*type)]){
    /* we couldn't write the whole dataset */
    inf_error("octopus.write_binary in actual writing");
    *ierr = 3;
  }

  /* close the file */
  close(fd);
  iio++;
  return;
}

/* this function neither allocates nor deallocates 'hp' */
void io_read_header(header_t * hp, int * correct_endianness, fint * ierr, fint * iio, STR_F_TYPE fname STR_ARG1)
{
  char * filename;
  int fd;
  ssize_t moved;

  *ierr = 0;
  TO_C_STR1(fname, filename);
  fd = open(filename, O_RDONLY);
  *iio += 100;
  free(filename);
  if(fd < 0){
    *ierr = 2;
    return;
  }

  assert(hp != NULL);

  /* read header */
  moved = read(fd, hp, sizeof(header_t));
  if ( moved != sizeof(header_t) ) { 
    /* we couldn't read the complete header */
    *ierr = 3;
    return;
  }

  *ierr = check_header(hp, correct_endianness);
  if( *ierr != 0 ){
    return;
  }

  close(fd);
  *iio++;
}

void FC_FUNC_(read_binary,READ_BINARY)
     (const fint * np, const fint * offset, byte * ff, fint * output_type, fint * ierr, fint * iio, STR_F_TYPE fname STR_ARG1)
{
  header_t * hp;
  char * filename;
  unsigned long fname_len;
  int fd, ii;
  ssize_t moved;
  int correct_endianness;
  byte * read_f;
  
  assert(np > 0);

  /* read the header */
  fname_len = l1;
  hp = (header_t *) malloc(sizeof(header_t));
  assert(hp != NULL);
  io_read_header(hp, &correct_endianness, ierr, iio, fname, fname_len);
  if (*ierr != 0) {
     free(hp);
     return;
  }
  
  /* check whether the sizes match */ 
  if( hp->np < *np + *offset ){ 
    *ierr = 4;
    free(hp);
    return; 
  }

  TO_C_STR1(fname, filename);
  fd = open(filename, O_RDONLY);
  *iio += 100;
  free(filename);
  
  if(fd < 0){
    *ierr = 2;
    free(hp);
    return;
  }

  if( hp->type == *output_type){
    /* format is the same, we just read */
    read_f = ff;
  } else {
    /*format is not the same, we store into a temporary array */
    read_f =(byte *) malloc((*np)*size_of[hp->type]);
  }

  /* set the start point */
  if(*offset != 0)
    lseek(fd, (*offset)*size_of[hp->type]+sizeof(header_t), SEEK_SET);
  else
    lseek(fd, sizeof(header_t), SEEK_SET);
  
  /* now read the values and close the file */
  moved = read(fd, read_f, (*np)*size_of[hp->type]);

  close(fd);
  *iio++;
  
  if ( moved != (*np)*size_of[hp->type]) { 
    /* we couldn't read the whole dataset */
    *ierr = 3;
    free(hp);
    if(hp->type != *output_type) {
       free(read_f);
    }
    return;
  }
    
  /* convert endianness */
  
  if(correct_endianness) {
    for(ii=0; ii < (*np)*size_of[hp->type] ; ii+=base_size_of[hp->type]) 
      endian_convert(base_size_of[hp->type], (char *) (read_f + ii));
  }

  /* convert values if it is necessary */
  if( hp->type != *output_type ){
    
    if(is_integer[hp->type] || is_integer[*output_type]){
      *ierr = 5;
    } else {

      for(ii=0; ii < *np ; ii++) 
	convert( (multi *) (read_f + ii*size_of[hp->type]), 
		 (multi *) (ff + ii*size_of[*output_type]), 
		 hp->type, *output_type);

      /* set the error code according to the conversion done (see src/basic/io_binary.h) */
      if ( hp->type == TYPE_FLOAT )          *ierr = -1;
      if ( hp->type == TYPE_FLOAT_COMPLEX )  *ierr = -2;
      if ( hp->type == TYPE_DOUBLE )         *ierr = -3;
      if ( hp->type == TYPE_DOUBLE_COMPLEX ) *ierr = -4;
    }
    free(read_f);
  }
  
  free(hp);
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
     (fint * np, fint * type, fint * file_size, fint * ierr, fint * iio, STR_F_TYPE fname STR_ARG1)
{
  header_t * hp;
  int correct_endianness;
  unsigned long fname_len;
  char * filename;
  struct stat st;

  hp = (header_t *) malloc(sizeof(header_t));
  assert(hp != NULL);

  /* read header */
  fname_len = l1;
  io_read_header(hp, &correct_endianness, ierr, iio, fname, fname_len);

  *np  = hp->np;
  *type = (int) hp->type;
  free(hp);

  TO_C_STR1(fname, filename);
  stat(filename, &st);
  *file_size = (int) st.st_size;
  free(filename);
}
