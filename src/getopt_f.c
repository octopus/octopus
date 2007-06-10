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
  printf("Usage: oct-oscillator_strength [OPTIONS] w\n");
  printf("Options:\n");
  printf("  -h              Print this help and exits.\n");
  printf("  -f <file>       Name of the 'multipoles' file to process.\n");
  printf("  -s <dw>         Limits of the search interval: [w-dw,w+dw]\n");
  exit(-1);
}

/***************************************************************/




/* FUNCTIONS TO BE USED BY THE PROGRAM oct-oscillator-strength */

void FC_FUNC_(getopt_oscillator_strength, GETOPT_OSCILLATOR_STRENGTH)
  (double *omega, STR_F_TYPE filename, double *searchinterval STR_ARG1)
{
  int c;

  if(argc==1) oscillator_strength_help();

  while (1) {
    c = getopt(argc, argv, "hf:s:");
    if (c == -1) break;
    switch (c) {

    case 'h':
      oscillator_strength_help();
      break;

    case 'f':
      TO_F_STR1(optarg, filename);
      break;

    case 's':
      *searchinterval = (double)atof(optarg);
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
