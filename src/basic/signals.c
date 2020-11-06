/*
 Copyright (C) 2016 X. Andrade

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
#ifdef HAVE_SIGNAL_H
#include <signal.h>
#endif
#include <stdlib.h>
#include <unistd.h>
#include <fortran_types.h>
#include "string_f.h" /* fortran <-> c string compatibility issues */
#include <string.h>

void FC_FUNC_(block_signals, BLOCK_SIGNALS)(){
#if defined(HAVE_SIGACTION) && defined(HAVE_SIGNAL_H)
  struct sigaction act;

  act.sa_handler = SIG_IGN;
  sigemptyset(&act.sa_mask);
  act.sa_flags = 0;

  sigaction(SIGINT, &act, 0);
  sigaction(SIGTERM, &act, 0);
#endif    
}

void FC_FUNC_(unblock_signals, UNBLOCK_SIGNALS)(){
#if defined(HAVE_SIGACTION) && defined(HAVE_SIGNAL_H)
  struct sigaction act;

  act.sa_handler = SIG_DFL;
  sigemptyset(&act.sa_mask);
  act.sa_flags = 0;

  sigaction(SIGINT, &act, 0);
  sigaction(SIGTERM, &act, 0);

#endif
}

void handle_segv(int *);

#if defined(HAVE_SIGACTION) && defined(HAVE_SIGNAL_H)
					       
void segv_handler(int signum, siginfo_t * si, void * vd){
  handle_segv(&signum);
  signal(signum, SIG_DFL);
  kill(getpid(), signum);
}

#endif
						
void FC_FUNC_(trap_segfault, TRAP_SEGFAULT)(){
#if defined(HAVE_SIGACTION) && defined(HAVE_SIGNAL_H)
  struct sigaction act;

  sigemptyset(&act.sa_mask);
  act.sa_sigaction = segv_handler;
  act.sa_flags = SA_SIGINFO;

  sigaction(SIGTERM, &act, 0);
  sigaction(SIGKILL, &act, 0);
  sigaction(SIGSEGV, &act, 0);
  sigaction(SIGABRT, &act, 0);
  sigaction(SIGINT,  &act, 0);
  sigaction(SIGBUS,  &act, 0);
  sigaction(SIGILL,  &act, 0);
  sigaction(SIGTSTP, &act, 0);
  sigaction(SIGQUIT, &act, 0);
  sigaction(SIGFPE,  &act, 0);
  sigaction(SIGHUP,  &act, 0);
  
#endif
}

void FC_FUNC_(get_signal_description, GET_SIGNAL_DESCRIPTION)(fint * signum, STR_F_TYPE const signame STR_ARG1){
#if defined(HAVE_STRSIGNAL) && defined(HAVE_STRING_H)
  TO_F_STR1(strsignal(*signum), signame);
#else
  TO_F_STR1("(description not available)", signame);
#endif
}
