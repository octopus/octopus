/*
 Copyright (C) 2006 X. Andrade

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

 $Id: recipes.c 2146 2006-05-23 17:36:00Z xavier $
*/

#include <config.h>
#ifdef HAVE_SIGNAL_H
#include <signal.h>
#endif
#include <stdlib.h>


void FC_FUNC_(block_signals, BLOCK_SIGNALS)(){
#if HAVE_SIGACTION && HAVE_SIGNAL_H
  struct sigaction act;

  act.sa_handler = SIG_IGN;
  sigemptyset(&act.sa_mask);
  act.sa_flags = 0;

  sigaction(SIGINT, &act, 0);
  sigaction(SIGTERM, &act, 0);
#endif    
}

void FC_FUNC_(unblock_signals, UNBLOCK_SIGNALS)(){
#if HAVE_SIGACTION && HAVE_SIGNAL_H
  struct sigaction act;

  act.sa_handler = SIG_DFL;
  sigemptyset(&act.sa_mask);
  act.sa_flags = 0;

  sigaction(SIGINT, &act, 0);
  sigaction(SIGTERM, &act, 0);

#endif
}

