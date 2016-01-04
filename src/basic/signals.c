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

 $Id$
*/

#include <config.h>
#ifdef HAVE_SIGNAL_H
#include <signal.h>
#endif
#include <stdlib.h>
#include <unistd.h>

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

void FC_FUNC_(dump_call_stack, DUMP_CALL_STACK)(void);

#if defined(HAVE_SIGACTION) && defined(HAVE_SIGNAL_H)
					       
void segv_handler(int signum, siginfo_t * si, void * vd){
  FC_FUNC_(dump_call_stack, DUMP_CALL_STACK)();
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

  sigaction(SIGSEGV, &act, 0);
  sigaction(SIGABRT, &act, 0);
  
#endif
}
