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
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 02111-1307, USA.

 $Id$
*/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <locale.h>
#include <dirent.h>
#include <time.h>

#include "string_f.h"

void FC_FUNC_(oct_printrecipe, OCT_PRINTRECIPE)
		 (STR_F_TYPE _dir, STR_F_TYPE filename STR_ARG2)
{

#if HAVE_SCANDIR && HAVE_ALPHASORT
  char *lang, *tmp, dir[512];
	struct dirent **namelist;
	int i, n;

	/* get language */
	lang = getenv("LANG");
	if(lang == NULL) lang = "en";

	/* convert directory from Fortran to C string */
	tmp = TO_C_STR1(_dir);
	strcpy(dir, tmp);
	free(tmp);

	strcat(dir, "/recipes");

	/* check out if lang dir exists */
	n = scandir(dir, &namelist, 0, alphasort);
	if (n < 0){
		printf("Directory does not exist: %s", dir);
		return;
	}

	for(i=0; i<n; i++)
		if(strncmp(lang, namelist[i]->d_name, 2) == 0){
			strcat(dir, "/");
			strcat(dir, namelist[i]->d_name);
			break;
		}

	if(i == n)
		strcat(dir, "/en"); /* default */

	/* clean up */
	for(i=0; i<n; i++)
		free(namelist[i]);
	free(namelist);

	/* now we read the recipes */
	n = scandir(dir, &namelist, 0, alphasort);
	
	/* initialize random numbers */
	srand((unsigned int)time(NULL));
	i = (int) ((double) (n-3 + 1.0) * (rand()/(RAND_MAX + 1.0)));

	strcat(dir, "/");
	strcat(dir, namelist[i+2]->d_name); /* skip ./ and ../ */

	/* clean up again */
	for(i=0; i<n; i++)
		free(namelist[i]);
	free(namelist);

	TO_F_STR2(dir, filename);
#endif
}
