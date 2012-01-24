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

#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>

#ifdef HAVE_UNAME
#include <sys/utsname.h>
#endif

#ifdef HAVE_IOCTL
#include <sys/ioctl.h>
#endif

#if defined(HAVE_TCGETPGRP)
#include <termios.h>
#endif

#include "varia.h"


/**
 * Gets the name of the machine 
 */
void sysname(char **c)
{
#ifdef HAVE_UNAME
  struct utsname name;
  uname(&name);
  *c = (char *)malloc(sizeof(name.nodename) + sizeof(name.sysname) + 4);
  strcpy(*c, name.nodename);
  strcat(*c, " (");
  strcat(*c, name.sysname);
  strcat(*c, ")");
#else
  *c = (char *)malloc(8);
  strcpy(*c, "unknown");
#endif
}


/**
 * Optimizes the order of the FFT grid.
 * The best FFT grid dimensions are given by 2^a*3^b*5^c*7^d*11^e*13^f
 * where a,b,c,d are arbitrary and e,f are 0 or 1.
 * (http://www.fftw.org/doc/Complex-DFTs.html)
 * par is the parity: the result must satisfy n % 2 == par, provided par >= 0.
 */
void fft_optimize(int *n, int par)
{
  if(*n <= 2) return;

  for(;; (*n)++){
    int i, n2;

    if((par >= 0) && (*n % 2 != par)) continue;
    
    /* For debugging:                 */
    /* printf("%i has factors ", *n); */

    n2 = *n;
    for(i = 2; i <= n2; i++){
      if(n2 % i == 0){
        /* For debugging:    */
	/* printf("%i ", i); */
	if(i > 13) break;
	n2 = n2 / i;
	if(i != 11 && i != 13) i--;
      }
    }
    /* For debugging: */
    /* printf("\n");  */
    if(n2 == 1) return;
  }
}

/**
 * returns true if process is in the foreground
 * copied from openssh scp source 
 */
static int foreground_proc(void)
{
#if defined(HAVE_TCGETPGRP) || defined(HAVE_IOCTL)
  static pid_t pgrp = -1;
  int ctty_pgrp;
  
  if (pgrp == -1) pgrp = getpgrp();
  
#if defined(HAVE_TCGETPGRP)
  return ((ctty_pgrp = tcgetpgrp(STDOUT_FILENO)) != -1 && ctty_pgrp == pgrp);
#else defined(HAVE_IOCTL)
  return ((ioctl(STDOUT_FILENO, TIOCGPGRP, &ctty_pgrp) != -1 && ctty_pgrp == pgrp));
#endif

#else
  return 0;
#endif
}

int getttywidth(void)
{
#ifdef HAVE_IOCTL  
  struct winsize winsize;
  
  if (ioctl(fileno(stdout), TIOCGWINSZ, &winsize) != -1)
    return (winsize.ws_col ? winsize.ws_col : 80);
  else
#endif
    return (80);
}

/**
 * displays progress bar with a percentage 
 */
void progress_bar(int actual, int max)
{
  static struct timeval start;
  static int old_pos, next_print;
  struct timeval now;
  char buf[512], fmt[64];
  int i, j, ratio, barlength, remaining;
  double elapsed;
  
  if(actual < 0) {
    (void) gettimeofday(&start, (struct timezone *) 0);
    actual  = 0;
    old_pos = 0;
    next_print = 10;
  }

  if(max > 0){
    ratio = 100 * actual / max;
    if(ratio < 0)   ratio = 0;
    if(ratio > 100) ratio = 100;
  }else
    ratio = 100;

  if(foreground_proc() == 0){
    if(old_pos == 0){
      printf("ETA: ");
    }

    barlength = getttywidth() - 6;

    j = actual*(barlength - 1)/max;
    if(j > barlength || actual == max) j = barlength;
    if(j < 1) j = 1;

    if(j > old_pos){
      for(i=old_pos+1; i<=j; i++)
	if(i*100/barlength >= next_print){
	  printf("%1d", next_print/10 % 10);
	  next_print += 10;
	}else
	  printf(".");
      old_pos = j;
      if(j == barlength) printf("\n");
    }

  }else{

    sprintf(buf, "%d", max);
    i = strlen(buf);
    if(i<3) i=3;
    sprintf(fmt, "\r[%%%dd/%%%dd]", i, i);
    sprintf(buf, fmt, actual, max);
    sprintf(buf + strlen(buf), " %3d%%" , ratio);
    
    barlength = getttywidth() - strlen(buf) - 16;
    if (barlength > 0) {
      i = barlength * ratio / 100;
      sprintf(buf + strlen(buf),
	      "|%.*s%*s|", i,
	      "*******************************************************"
	      "*******************************************************"
	      "*******************************************************"
	      "*******************************************************"
	      "*******************************************************"
	      "*******************************************************"
	      "*******************************************************",
	      barlength - i, "");
    }
    
    /* time information now */
    (void) gettimeofday(&now, (struct timezone *) 0);
    elapsed = now.tv_sec - start.tv_sec;
  
    if(elapsed <= 0.0 || actual <= 0) {
      sprintf(buf + strlen(buf),
	      "     --:-- ETA");
    }else{
      remaining = (int)(max / (actual / elapsed) - elapsed);
      if(remaining < 0) remaining = 0;
    
      i = remaining / 3600;
      if(i)
	sprintf(buf + strlen(buf), "%4d:", i);
      else
	sprintf(buf + strlen(buf), "     ");
      i = remaining % 3600;
      sprintf(buf + strlen(buf), "%02d:%02d%s", i / 60, i % 60, " ETA");
    }
    printf("%s", buf);

  }
  
  fflush(stdout);
}

#undef timersub
