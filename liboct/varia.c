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

#include <math.h>
#include <sys/time.h>
#include <sys/ioctl.h>
#include <stdio.h>
#include <termios.h>
#include <unistd.h>

#include "config.h"
#include "varia.h"

/* optimizes the order of the fft
	 p is the maximum prime allowed in n */
void fft_optimize(int *n, int p, int par)
{
	if(*n <= 2) return;

  for(;; (*n)++){
		int i, n2;

		if((par > 0) && (*n % 2 != par)) continue;

		n2 = *n;
		for(i = 2; i<=n2; i++){
			if(n2 % i == 0){
				if(i > p) break;
				n2 = n2 / i; 
				i--; 
			}
		}
		if(i > n2) return;
	}
}

/* returns true if process is in the foreground
	 copied from openssh scp source */
static int foreground_proc(void)
{
	static __pid_t pgrp = -1; /* WARNING should be __pid_t */
	int ctty_pgrp;
	
	if (pgrp == -1)
		pgrp = getpgrp();
	
#ifdef HAVE_TCGETPGRP
	return ((ctty_pgrp = tcgetpgrp(STDOUT_FILENO)) != -1 &&
					ctty_pgrp == pgrp);
#else
	return ((ioctl(STDOUT_FILENO, TIOCGPGRP, &ctty_pgrp) != -1 &&
					 ctty_pgrp == pgrp));
#endif
}

int getttywidth(void)
{
	struct winsize winsize;
	
	if (ioctl(fileno(stdout), TIOCGWINSZ, &winsize) != -1)
		return (winsize.ws_col ? winsize.ws_col : 80);
	else
		return (80);
}

/* displays progress bar with a percentage */
void progress_bar(int actual, int max)
{
	static struct timeval start;
	struct timeval now, td, wait;
	char buf[512], fmt[64];
	int i, ratio, barlength, remaining;
	double elapsed;

	if(actual < 0) {
		(void) gettimeofday(&start, (struct timezone *) 0);
		actual = 0;
	}

	if (foreground_proc() == 0)
		return;
	
	if(max > 0){
		ratio = 100 * actual / max;
		if(ratio < 0)   ratio = 0;
		if(ratio > 100) ratio = 100;
	}else
		ratio = 100;

	sprintf(buf, "%d", max);
	i = strlen(buf);
	if(i<3) i=3;
	sprintf(fmt, "\r[%%%dd/%%%dd]", i, i);
	snprintf(buf, sizeof(buf), fmt, actual, max);
	snprintf(buf + strlen(buf), sizeof(buf) - strlen(buf), " %3d%%" , ratio);

	barlength = getttywidth() - strlen(buf) - 16;
	if (barlength > 0) {
		i = barlength * ratio / 100;
		snprintf(buf + strlen(buf), sizeof(buf) - strlen(buf),
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
	
	if(elapsed < 0.0 || actual <= 0) {
		snprintf(buf + strlen(buf), sizeof(buf) - strlen(buf),
						 "     --:-- ETA");
	}else{
		remaining = (int)(max / (actual / elapsed) - elapsed);
		if(remaining < 0) remaining = 0;

		i = remaining / 3600;
		if(i)
			snprintf(buf + strlen(buf), sizeof(buf) - strlen(buf),
							 "%4d:", i);
		else
			snprintf(buf + strlen(buf), sizeof(buf) - strlen(buf),
							 "     ");
		i = remaining % 3600;
		snprintf(buf + strlen(buf), sizeof(buf) - strlen(buf),
						 "%02d:%02d%s", i / 60, i % 60,
						 " ETA");
	}

	printf("%s", buf);
	fflush(stdout);
}

#undef timersub
