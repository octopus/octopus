#include <config.h>
#ifdef HAVE_SIGNAL_H
#include <signal.h>
#endif
#include <stdlib.h>


void FC_FUNC(block_signals, BLOCK_SIGNALS)(){
#if HAVE_SIGACTION && HAVE_SIGNAL_H
  struct sigaction act;

  act.sa_handler = SIG_IGN;
  sigemptyset(&act.sa_mask);
  act.sa_flags = 0;

  sigaction(SIGINT, &act, 0);
  sigaction(SIGTERM, &act, 0);
#endif    
}

void FC_FUNC(unblock_signals, UNBLOCK_SIGNALS)(){
#if HAVE_SIGACTION && HAVE_SIGNAL_H
  struct sigaction act;

  act.sa_handler = SIG_DFL;
  sigemptyset(&act.sa_mask);
  act.sa_flags = 0;

  sigaction(SIGINT, &act, 0);
  sigaction(SIGTERM, &act, 0);

#endif
}

