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
*/

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <locale.h>
#include <dirent.h>
#include <time.h>

int F90_FUNC_(print_file, PRINT_FILE)
     (char *name)
{
  FILE *pf;
  int c;

  pf = fopen(name, "r");
  if (pf == NULL) return 1;

	while((c = fgetc(pf)) != EOF)
		fputc(c, stdout);
	fclose(pf); 

	fflush(stdout); 
	return 0;
}

void F90_FUNC_(oct_printrecipe, OCT_PRINTRECIPE)
		 (char *_dir)
{
	char *lang, dir[512];
	struct dirent **namelist;
	int i, n;
	FILE *f;

	/* get language */
	lang = getenv("LANG");
	if(lang == NULL) lang = "en";

	/* check out if lang dir exists */
	strcpy(dir, _dir);
	strcat(dir, "/recipes");

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
	
	/* output selected file */
	F90_FUNC_(print_file, PRINT_FILE) (dir);

	/* print disclaimer */
	strcpy(dir, _dir);
	strcat(dir, "/recipes/disclaimer.txt");
	printf("\n\n");
	F90_FUNC_(print_file, PRINT_FILE) (dir);
	printf("\n");
}
