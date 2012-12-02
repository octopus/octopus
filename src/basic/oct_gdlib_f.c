/*
 Copyright (C) 2002-2006 the octopus team

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

 $Id: oct_gsl_f.c 1935 2006-03-20 21:02:08Z appel $
*/

#include <config.h>

#ifdef HAVE_GDLIB

#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <gd.h>

#include "string_f.h" /* Fortran <-> C string compatibility issues */

/* ---------------------- Interface to GD functions ------------------------ */
gdImagePtr FC_FUNC_(oct_gdimage_create_from, OCT_GDIMAGE_CREATE_FROM)
  (STR_F_TYPE name STR_ARG1)
{
  char *name_c, *ext;
  FILE *in;
  gdImagePtr im;

  TO_C_STR1(name, name_c);

  if((in = fopen(name_c, "rb")) == NULL) 
    return NULL; /* could not open file */

  /* get extension of filename */
  for(ext=name_c+strlen(name_c); *ext!='.' && ext>=name_c; ext--){
    *ext = tolower(*ext);
  }
  if(ext < name_c || ext == name_c+strlen(name_c))
    return NULL; /* could not find file type */

  /* get rid of . in extension */
  ext++;

  /* load image file */
  im = NULL;
#ifdef HAVE_GD_JPEG
  if((strcmp(ext, "jpg") == 0) || (strcmp(ext, "JPG") == 0) ||
     (strcmp(ext, "jpeg") == 0) || (strcmp(ext, "JPEG") == 0))
    im = gdImageCreateFromJpeg(in);
#endif

#ifdef HAVE_GD_PNG
  if ((strcmp(ext, "png") == 0) || (strcmp(ext, "PNG") == 0))
    im = gdImageCreateFromPng(in);
#endif

#ifdef HAVE_GD_GIF
  if ((strcmp(ext, "gif") == 0) || (strcmp(ext, "GIF") == 0))
    im = gdImageCreateFromGif(in);
#endif

  fclose(in);

  free(name_c);
  return im;
}

int FC_FUNC_(oct_gdimage_sx, OCT_GDIMAGE_SX)
  (const gdImagePtr *im)
{
  assert(*im != NULL);

  return gdImageSX(*im);
}

int FC_FUNC_(oct_gdimage_sy, OCT_GDIMAGE_SY)
  (const gdImagePtr *im)
{
  assert(*im != NULL);

  return gdImageSY(*im);
}

void FC_FUNC_(oct_gdimage_get_pixel_rgb, OCT_GDIMAGE_GET_PIXEL_RGB)
  (const gdImagePtr *im, const int *x, const int *y, int *r, int *g, int *b)
{
  int color;

  assert(*im != NULL);

  if(gdImageBoundsSafe(*im, *x, *y)){
    color = gdImageGetPixel(*im, *x, *y);
    *r = gdImageRed  (*im, color);
    *g = gdImageGreen(*im, color);
    *b = gdImageBlue (*im, color);
  }else{
    *r = 0; *g = 0; *b = 0;
  }
}

#else
/* this is to avoid an empty source file (not allowed by ANSI C)*/
void useless(){}
#endif
/* defined HAVEGDLIB */
