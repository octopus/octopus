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

	$Id: oscillator_strength_clarg.c 2516 2006-10-24 21:31:59Z acastro $
*/

#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "string_f.h"

/* GENERAL FUNCTIONS AND VARIABLES */

char **argv;
int argc;

void FC_FUNC_(set_number_clarg, SET_NUMBER_CLP)(int *nargc)
{
  argc = *nargc+1;
  argv = (char **)malloc(argc*sizeof(char *));
}

void FC_FUNC_(set_clarg, SET_CLARG)(int *i, STR_F_TYPE arg STR_ARG1)
{ argv[*i] = TO_C_STR1(arg);}


void oscillator_strength_help(){
  printf("Usage: oct-oscillator_strength [OPTIONS] [w]\n");
  printf("\n");
  printf("Options:\n");
  printf("  -h              Print this help and exits.\n");
  printf("  -m <mode>       Select the run mode: .\n");
  printf("  -O <operator>   Selects the operator to be analyzed:\n");
  printf("                    o If <operatot> is a pair of integers in the form '(l,m)', then\n");
  printf("                      the operator will be the (l,m) multipole.\n");
  printf("                    o If <operatot> is x, y, or z, then the response operator to\n");
  printf("                      be analyzed will be the dipole in the given direction.\n");
  printf("                    o If the -O option is not given in the command line, then the\n");
  printf("                      observation operator O will be the same as the perturbation\n");
  printf("                      operator that defines the initial kick.\n");
  printf("  -s <dw>         Limits of the search interval: [w-dw,w+dw]\n");
  printf("  -r <r>          Number of resonances to search for.\n");
  printf("  -n <N>          Number of frequencies in which the search interval\n");
  printf("                    is discretized (default 1000)\n");
  printf("  -k <k>          Process, or generate, the k-th order response.\n");
  printf("  -t <time>       The signal analysis will be done by integrating in the \n");
  printf("                    time interval [0, <time>]. If this argument is absent,\n");
  printf("                    it makes use of all the time-window present in the multipoles\n");
  printf("                    files.\n");
  exit(-1);
}

/***************************************************************/




/* FUNCTIONS TO BE USED BY THE PROGRAM oct-oscillator-strength */

void FC_FUNC_(getopt_oscillator_strength, GETOPT_OSCILLATOR_STRENGTH)
  (int *mode, double *omega, double *searchinterval, int *order, 
   int *nresonances, int *nfrequencies, double *time, int *l, int *m, 
   STR_F_TYPE ffile STR_ARG1)
{
  int c;

  /* This line would be present if we wanted to make the omega a 
     mandatory argument. But for the moment I think it should not be mandatory.
     if(argc==1) oscillator_strength_help(); */

  while (1) {
    c = getopt(argc, argv, "hm:s:k:O:r:n:t:f:");
    if (c == -1) break;
    switch (c) {

    case 'h':
      oscillator_strength_help();
      break;

    case 'm':
      *mode = (int)atoi(optarg);
      break;

    case 's':
      *searchinterval = (double)atof(optarg);
      break;

    case 'O':
      c = sscanf(optarg, "(%d,%d)", l, m);
      if(c != 2){
        switch (optarg[0]){
        case 'x':
          *l = 0;
          *m = 1;
          break;
        case 'y':
          *l = 0;
          *m = 2;
          break;
        case 'z':
          *l = 0;
          *m = 3;
          break;
        default:
          printf("Problem reading the -O option value.\n\n");
          oscillator_strength_help();
        }
      }
      break;

    case 'k':
      *order = (int)atoi(optarg);
      break;

    case 'r':
      *nresonances = (int)atoi(optarg);
      break;

    case 'n':
      *nfrequencies = (int)atoi(optarg);
      break;

    case 't':
      *time = (double)atof(optarg);
      break;

    case 'f':
      TO_F_STR1(optarg, ffile);
      break;

    case '?':
      oscillator_strength_help();
      break;

    }
  }
  if (optind < argc) {
    while (optind < argc) *omega = (double)atof(argv[optind++]);
  }

}
/***************************************************************/
