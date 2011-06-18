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
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifdef HAVE_DIRENT_H
#include <dirent.h>
#endif

#include "string_f.h" /* Fortran <-> c string compatibility issues */

#ifndef F2003_C_PTR
/* Pointer-handling functions */

void FC_FUNC_(set_null, SET_NULL)
  (void ** ptr)
{
  *ptr = NULL;
}

int FC_FUNC_(is_null_int, is_null_int)
  (void ** ptr)
{
  return (*ptr == NULL);
}
#endif

/* *********************** interface functions ********************** */

void FC_FUNC_(oct_mkdir, OCT_MKDIR)
		 (STR_F_TYPE name STR_ARG1)
{
  struct stat buf;
  char *name_c;

  TO_C_STR1(name, name_c);
  if(!*name_c) return;

  if(stat(name_c, &buf) == 0){
    free(name_c);
    return;
  }

#ifndef _WIN32
  mkdir(name_c, 0775);
#else
  mkdir(name_c);
#endif

  free(name_c);
}

void FC_FUNC_(oct_stat, OCT_STAT)
	(int *ierr, STR_F_TYPE name STR_ARG1)
{
  char *name_c;
  struct stat statbuf;

  TO_C_STR1(name, name_c);
  /* for the moment we are only interested in the return value
     could be extended if required. */
  *ierr = stat(name_c, &statbuf);
}

void FC_FUNC_(oct_rm, OCT_RM)
     (STR_F_TYPE name STR_ARG1)
{
  char *name_c;

  TO_C_STR1(name, name_c);
  unlink(name_c);
  free(name_c);
}

void FC_FUNC_(oct_getcwd, OCT_GETCWD)
  (STR_F_TYPE name STR_ARG1)
{
  char s[256];
  getcwd(s, 256);
  TO_F_STR1(s, name);
}

void FC_FUNC_(oct_getenv, OCT_GETENV)
  (STR_F_TYPE var, STR_F_TYPE value STR_ARG2)
{
  char *name_c, *var_c;

  TO_C_STR1(var, name_c);
  var_c = getenv(name_c);
  free(name_c);

  if(var_c != NULL){
    TO_F_STR2(var_c, value);
  }else{
    TO_F_STR2("", value);
  }
}

/* this function gets a string of the form '1-12, 34' and fills
	 array l with the 1 if the number is in the list, or 0 otherwise */
void FC_FUNC_(oct_wfs_list, OCT_WFS_LIST)
		 (STR_F_TYPE str, int l[16384] STR_ARG1)
{
  int i, i1, i2;
  char c[20], *c1, *str_c, *s;

  TO_C_STR1(str, str_c);
  s = str_c;
  
  /* clear list */
  for(i=0; i<16384; i++)
    l[i] = 0;
  
  while(*s){
    /* get integer */
    for(c1 = c; isdigit(*s) || isspace(*s); s++)
      if(isdigit(*s)) *c1++ = *s;
    *c1 = '\0';
    i1 = atoi(c) - 1;
    
    if(*s == '-'){ /* range */
      s++;
      for(c1 = c; isdigit(*s) || isspace(*s); s++)
	if(isdigit(*s)) *c1++ = *s;
      *c1 = '\0';
      i2 = atoi(c) - 1;
    } else /* single value */
      i2 = i1;
    
    for(i=i1; i<=i2; i++)
      if(i>=0 && i<16384)
				l[i] = 1;

    if(*s) s++;
  }

  free(str_c);
}

/* ------------------------------ from varia.c ------------------------------- */
#include "varia.h"

void FC_FUNC_(oct_fft_optimize, OCT_FFT_OPTIMIZE)
  (int *n, int *p, int *par)
{
  fft_optimize(n, *p, *par);
}

void FC_FUNC_(oct_progress_bar, OCT_PROGRESS_BAR)
  (int *a, int *max)
{
  progress_bar(*a, *max);
}

/* ------------------------------ some stuff  -------------------------------- */
void FC_FUNC_(oct_gettimeofday, OCT_GETTIMEOFDAY)
  (int *sec, int *usec)
{
#ifdef HAVE_GETTIMEOFDAY
  struct timeval tv;
  struct timezone tz;

  gettimeofday(&tv, &tz);

  /* The typecast below should use long. However, this causes incompatibilities
     with Fortran integers. 
     Using int will cause wrong results when tv.tv_sec exceeds INT_MAX=2147483647 */
  *sec  = (int) tv.tv_sec;
  *usec = (int) tv.tv_usec;
/*
  char str[sizeof("HH:MM:SS")];
  time_t local;
  local = tv.tv_sec;
  strftime(str, sizeof(str), "%T", localtime(&local));
  printf("%s.%06ld \n", str, (long) tv.tv_usec);
  printf("%ld.%06ld \n", (long) tv.tv_sec, (long) tv.tv_usec);
*/
#else
  *sec  = 0;
  *usec = 0; 
#endif
}

double FC_FUNC_(oct_clock, OCT_CLOCK)
  ()
{
#ifdef HAVE_GETTIMEOFDAY
  int sec, usec;
  FC_FUNC_(oct_gettimeofday, OCT_GETTIMEOFDAY) (&sec, &usec);
  return sec + 1.0e-6*usec;
#else
  return (double) clock()/CLOCKS_PER_SEC;
#endif
}

void FC_FUNC_(oct_nanosleep, OCT_NANOSLEEP)
	(int *sec, int *usec)
{
#ifdef linux
  /* Datatypes should be long instead of int (see comment in gettimeofday) */
  struct timespec req, rem;
  req.tv_sec  = (time_t) *sec;
  req.tv_nsec = (long)   *usec;
  nanosleep(&req, &rem);
#endif
}


void FC_FUNC_(oct_sysname, OCT_SYSNAME)
		 (STR_F_TYPE name STR_ARG1)
{
  char *name_c;
  
  sysname(&name_c);
  TO_F_STR1(name_c, name);
  free(name_c);
}


int FC_FUNC_(number_of_lines, NUMBER_OF_LINES)
  (STR_F_TYPE name STR_ARG1)
{

  FILE *pf;
  int c, i;
  char *name_c;

  TO_C_STR1(name, name_c);
  pf = fopen(name_c, "r");
  free(name_c);
  
  if (pf != NULL) {
    i = 0;
    while ((c = getc(pf)) != EOF) {
      if (c == '\n') i++;
    }
    fclose(pf); 
    return i;
  }else{
    return -1;
  }
}


/* Given a string in C, it breaks it line by line and returns each 
   as a Fortran string. Returns 0 if string does not have more lines. 
*/
void FC_FUNC_(break_c_string, BREAK_C_STRING)
  (char **str, char **s, STR_F_TYPE line_f STR_ARG1)
{
  char *c, line[256]; /* hopefully no line is longer than 256 characters ;) */

  if(*s == NULL) *s = *str;

  if(*s == NULL || **s == '\0'){
    *s = (char *)(0);
    return;
  }

  for(c=line; **s!='\0' && **s!='\n'; (*s)++, c++)
    *c = **s;
  *c = '\0';
  if(**s=='\n') (*s)++;

  TO_F_STR1(line, line_f);
}

/* 

This function searches in directory given by dirname all files that have the following name:

*_<real_number>_<integer>*

It returns the value of <real_number> found that is closest to freq (or abs(freq))
and for which the value of <integer> matches with the tag argument.

The value found is returned in the freq argument.

ierr results:
0 : value found
1 : no matching file found
2 : cannot open the directory
 
*/

void FC_FUNC_(oct_search_file_lr, OCT_SEARCH_FILE_LR)
     (double * freq, int * tag, int * ierr, STR_F_TYPE dirname STR_ARG1)
{
#if HAVE_DIRENT_H && HAVE_CLOSEDIR && HAVE_READDIR && HAVE_STRCHR && HAVE_STRTOD

  DIR * dir;
  struct dirent * ent;
  char * name_c;
  char * num_start, * num_end;
  double read_value, min;
  int found_something, read_tag;

  TO_C_STR1(dirname, name_c);
  dir = opendir(name_c);

  if(dir == NULL){
    *ierr = 2;
    return;
  }
  free(name_c);

  ent = NULL;
  found_something = 0;

  while(1) {
    ent = readdir(dir);
    if( ent == NULL ) break;

    num_start = strchr(ent -> d_name, '_');

    if(num_start != NULL) {
      num_start++; /* now this points to the beginning of the number */

      /* take the numerical value from the string */
      read_value = strtod(num_start, &num_end);
      
      if ( num_end == num_start ) continue; /* no number found */

      /* check that we have the correct tag */
      if(num_end[0] == '_') {

	num_start = num_end + 1;
	read_tag = (int) strtol(num_start, &num_end, 10);
	if ( num_end == num_start ) continue; /* no tag found */
	if ( read_tag != *tag ) continue; /* tag does not match */
 
      } else continue;


      /* if this is the first number we found or the value is closer than previous */
      if ( !found_something || (fabs(fabs(min)-fabs(*freq)) > fabs(fabs(read_value)-fabs(*freq))) ) {
	min = read_value;
	found_something = 1;
      }
    }

  }

  closedir(dir);

  if(found_something) {
    *ierr = 0;
    *freq = min ;
  } else {
    *ierr = 1;
  }

#else
#warning directory search not compiled
  printf("Warning: Directory search not implemented.\n");

#endif

}

double FC_FUNC_(oct_hypotd, OCT_HYPOTD)
  (double *x, double *y)
{
  return hypot(*x, *y);
}

float FC_FUNC_(oct_hypotf, OCT_HYPOTF)
  (float *x, float *y)
{
  return hypotf(*x, *y);
}

void * FC_FUNC_(get_memory_usage, GET_MEMORY_USAGE)()
{
#ifdef linux
  static size_t pagesize = 0;
  FILE *f;
  int pid;
  unsigned long mem;
  char s[256];
  
  if(pagesize == 0)
    pagesize = sysconf(_SC_PAGESIZE);
  
  pid = getpid();
  sprintf(s, "%s%d%s", "/proc/", pid, "/statm");
  if((f = fopen(s, "r")) == (FILE *)NULL) return (void *)(-1);
  fscanf(f, "%lu", &mem);
  fclose(f);
  
  return (void *) (mem*pagesize);
#elif HAVE_SBRK
  return sbrk(0);
#else
  return 0;
#endif
}

void FC_FUNC_(loct_exit_failure, LOCT_EXIT_FAILURE)(const int * exit_code){
  exit(EXIT_FAILURE);
}
