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
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.

*/

#include <config.h>

#ifdef HAVE_GDLIB

#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <gd.h>

/* ---------------------- Interface to GD functions ------------------------ */
gdImagePtr gdlib_image_create_from(char * name)
{
  char *ext;
  FILE *in;
  gdImagePtr im;

  if((in = fopen(name, "rb")) == NULL) {
    return NULL; /* could not open file */
  }

  /* get extension of filename */
  for(ext=name+strlen(name); *ext!='.' && ext>=name; ext--){
    *ext = tolower(*ext);
  }
  if(ext < name || ext == name+strlen(name)) {
    fclose(in);
    return NULL; /* could not find file type */
  }

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

  return im;
}

int gdlib_image_sx(const gdImagePtr *im)
{
  assert(*im != NULL);

  return gdImageSX(*im);
}

int gdlib_image_sy(const gdImagePtr *im)
{
  assert(*im != NULL);

  return gdImageSY(*im);
}

void gdlib_image_get_pixel_rgb(const gdImagePtr *im, const int *x, const int *y, int *r, int *g, int *b)
{
  int color;

  assert(*im != NULL);

  if(gdImageBoundsSafe(*im, *x, *y)){
    color = gdImageGetPixel(*im, *x, *y);
    *r = gdImageRed  (*im, color);
    *g = gdImageGreen(*im, color);
    *b = gdImageBlue (*im, color);
  }else{
    /* this will happen for boundary points */
    //    fprintf(stderr, "Illegal pixel coordinate %d %d\n", *x, *y);
    *r = 0; *g = 0; *b = 0;
  }
}

#else
/* this is to avoid an empty source file (not allowed by ANSI C)*/
void useless(){}
#endif
/* defined HAVEGDLIB */
