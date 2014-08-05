/*
 Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch

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

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <locale.h>
#include <dirent.h>
#include <sys/time.h>
#include <gsl/gsl_rng.h>

#include "string_f.h"

unsigned long int random_seed()
{
 unsigned long int seed;
 FILE *devrandom;

 if ((devrandom = fopen("/dev/random","r")) == NULL) {
#ifdef HAVE_GETTIMEOFDAY
   struct timeval tv;
   gettimeofday(&tv, 0);
   seed = tv.tv_sec + tv.tv_usec;
#else
   seed = 0;
#endif
 } else {
   fread(&seed, sizeof(seed), 1, devrandom);
   fclose(devrandom);
 }

 return seed;
}

void FC_FUNC_(oct_printrecipe, OCT_PRINTRECIPE)
  (STR_F_TYPE _dir, STR_F_TYPE filename STR_ARG2)
{

#if HAVE_SCANDIR && HAVE_ALPHASORT
  char *lang, *tmp, dir[512];
  struct dirent **namelist;
  int ii, nn;
  gsl_rng *rng;

  /* get language */
  lang = getenv("LANG");
  if(lang == NULL) lang = "en";

  /* convert directory from Fortran to C string */
  TO_C_STR1(_dir, tmp);
  strcpy(dir, tmp);
  free(tmp);

  strcat(dir, "/recipes");

  /* check out if lang dir exists */
  nn = scandir(dir, &namelist, 0, alphasort);
  if (nn < 0){
    printf("Directory does not exist: %s", dir);
    return;
  }

  for(ii=0; ii<nn; ii++)
    if(strncmp(lang, namelist[ii]->d_name, 2) == 0){
      strcat(dir, "/");
      strcat(dir, namelist[ii]->d_name);
      break;
    }

  if(ii == nn)
    strcat(dir, "/en"); /* default */

  /* clean up */
  for(ii=0; ii<nn; ii++)
    free(namelist[ii]);
  free(namelist);

  /* now we read the recipes */
  nn = scandir(dir, &namelist, 0, alphasort);
	
  /* initialize random numbers */
  gsl_rng_env_setup();
  rng = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(rng, random_seed());
  ii = gsl_rng_uniform_int(rng, nn - 2);
  gsl_rng_free(rng);

  strcat(dir, "/");
  strcat(dir, namelist[ii+2]->d_name); /* skip ./ and ../ */

  /* clean up again */
  for(ii=0; ii<nn; ii++)
    free(namelist[ii]);
  free(namelist);

  TO_F_STR2(dir, filename);

#else
  printf("Sorry, recipes cannot be printed unless scandir and alphasort are available with your C compiler.\n");
#endif
}
