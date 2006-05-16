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

#include "string_f.h" /* fortran <-> c string compatibility issues */

/* *********************** interface functions ********************** */

void FC_FUNC_(oct_mkdir, OCT_MKDIR)
		 (STR_F_TYPE name STR_ARG1)
{
  struct stat buf;
  char *name_c;

  name_c = TO_C_STR1(name);
  if(!*name_c || stat(name_c, &buf) == 0) return;
  mkdir(name_c, 0775);
  free(name_c);
}

void FC_FUNC_(oct_stat, OCT_STAT)
	(int *ierr, STR_F_TYPE name STR_ARG1)
{
  char *name_c;
	struct stat statbuf;

  name_c = TO_C_STR1(name);
	/* for the moment we are only interested in the return value
		 could be extended if required. */
  *ierr = stat(name_c, &statbuf);
}

void FC_FUNC_(oct_rm, OCT_RM)
     (STR_F_TYPE name STR_ARG1)
{
  char *name_c;

  name_c = TO_C_STR1(name);
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
  char *name_c;
  char *var_c;

  var_c = (char *) malloc(256);
  var_c[0]='\0';
  name_c = TO_C_STR1(var);
  if(getenv(name_c) != NULL) var_c = getenv(name_c);
  TO_F_STR2(var_c, value);
}

/* this function gets a string of the form '1-12, 34' and fills
	 array l with the 1 if the number is in the list, or 0 otherwise */
void FC_FUNC_(oct_wfs_list, OCT_WFS_LIST)
		 (STR_F_TYPE str, int l[32] STR_ARG1)
{
  int i, i1, i2;
  char c[20], *c1, *str_c, *s;

  str_c = TO_C_STR1(str);
  s = str_c;
  
  /* clear list */
  for(i=0; i<32; i++)
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
    }else /* single value */
      i2 = i1;
    
    for(i=i1; i<=i2; i++)
      if(i>=0 && i<1024)
	l[i/32] |= 1 << (i % 32);
    
    if(*s) s++;
  }

  free(str_c);
}


/* from ylm.c */
double ylm(double x, double y, double z, int l, int m);

double FC_FUNC_(oct_ylm, OCT_YLM)
  (double *x, double *y, double *z, int *l, int *m)
{
  return ylm(*x, *y, *z, *l, *m);
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

/* -------------------------- interface to METIS ----------------------------- */
#if defined(HAVE_METIS)
#include <metis.h>

void FC_FUNC_(oct_metis_part_mesh_nodal, OCT_METIS_PART_MESH_NODAL)
  (int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag, int *nparts, 
   int *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshNodal(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart);
}

void FC_FUNC_(oct_metis_part_mesh_dual, OCT_METIS_PART_MESH_DUAL)
  (int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag, int *nparts, 
   int *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshDual(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart);
}

void FC_FUNC_(oct_metis_part_graph_recursive, OCT_METIS_PART_GRAPH_RECURSIVE)
  (int *n, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
   idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts,
   int *options, int *edgecut, idxtype *part)
{
  METIS_PartGraphRecursive(n, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}

void FC_FUNC_(oct_metis_part_graph_kway, OCT_METIS_PART_GRAPH_KWAY)
  (int *n, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
   idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts,
   int *options, int *edgecut, idxtype *part)
{
  METIS_PartGraphKway(n, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}

void FC_FUNC_(oct_metis_part_graph_vkway, OCT_METIS_PART_GRAPH_VKWAY)
  (int *n, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
   idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts,
   int *options, int *edgecut, idxtype *part)
{
  METIS_PartGraphVKway(n, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
#endif

/* ------------------------------ some stuff  -------------------------------- */
double FC_FUNC_(oct_clock, OCT_CLOCK)
  ()
{
  return (double) clock();
}

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

void FC_FUNC_(oct_nanosleep, OCT_CLOCK)
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


/* this function is *not* portable. Should get rid of this! */
int FC_FUNC_(oct_getmem, OCT_GETMEM)
     ()
{
#ifdef linux
  static size_t pagesize = 0;
  FILE *f;
  int pid;
  long mem;
  char s[256];
  
  if(pagesize == 0)
    pagesize = sysconf(_SC_PAGESIZE);
  
  pid = getpid();
  sprintf(s, "%s%d%s", "/proc/", pid, "/statm");
  if((f = fopen(s, "r")) == NULL) return -1;
  fscanf(f, "%lu", &mem);
  fclose(f);
  
  return (mem*pagesize) >> 10;
#else
  return -1;
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

  name_c = TO_C_STR1(name);
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
