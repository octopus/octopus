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
#include <unistd.h>
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
{
  char *c;
  TO_C_STR1(arg, c)
  argv[*i] = c;
}


/* FUNCTIONS TO BE USED BY THE PROGRAM oct-oscillator-strength */

void oscillator_strength_help(){
  printf("Usage: oct-oscillator_strength [OPTIONS] [w]\n");
  printf("\n");
  printf("Options:\n");
  printf("  -h              Print this help and exits.\n");
  printf("  -m <mode>       Select the run mode:\n");
  printf("                    1 (default) analyzes the signal present in an 'ot' file.\n");
  printf("                      This should have been generated by this same utility\n");
  printf("                      (run mode 2).\n");
  printf("                    2 Reads a number of 'multipoles' files, which should be\n");
  printf("                      present in the working directory, and be called\n");
  printf("                      'multipoles.1', 'multipoles.2', ..., and generate an 'ot'\n");
  printf("                      file with the k-th order response of a given operator O.\n");
  printf("                      The order k is decided by the '-k' option. The operator\n");
  printf("                      is decided by the '-O' option.\n");
  printf("                    3 Peforms an analysis of the second-order response of an\n");
  printf("                      operator O, present in the working directory, and\n");
  printf("                      previously generated with run mode 2. It also reads a\n");
  printf("                      file with a list of frequecies around which the search\n");
  printf("                      for resonances is performed.\n");
  printf("                    4 Reads an 'ot' file, and generates an 'omega' file with\n");
  printf("                      either the sine or cosine Fourier transform of the\n");
  printf("                      signal present in 'ot'.\n");
  printf("  -O <operator>   Selects the operator to be analyzed:\n");
  printf("                    o If <operatot> is a pair of integers in the form '(l,m)'\n");
  printf("                      then the operator will be the (l,m) multipole.\n");
  printf("                    o If <operatot> is x, y, or z, then the response operator\n");
  printf("                      to be analyzed will be the dipole in the given direction.\n");
  printf("                    o If the -O option is not given in the command line, then\n");
  printf("                      the observation operator O will be the same as the\n");
  printf("                      perturbation operator that defines the initial kick.\n");
  printf("  -f <file>       This is the file where the frequencies needed in run mode\n");
  printf("                  3 are stored.\n");
  printf("  -d <gamma>      gamma is the damping factor used in the SOS formulae that");
  printf("                  produce (hyper)-polarizabilities.");
  printf("  -s <dw>         Limits of the search interval: [w-dw,w+dw]\n");
  printf("  -r <r>          Number of resonances to search for.\n");
  printf("  -n <N>          Number of frequencies in which the search interval\n");
  printf("                    is discretized (default 1000)\n");
  printf("  -k <k>          Process, or generate, the k-th order response.\n");
  printf("  -t <time>       The signal analysis will be done by integrating in the \n");
  printf("                    time interval [0, <time>]. If this argument is absent,\n");
  printf("                    it makes use of all the time-window present in the\n");
  printf("                    multipoles files.\n");
  exit(-1);
}

/***************************************************************/

void FC_FUNC_(getopt_oscillator_strength, GETOPT_OSCILLATOR_STRENGTH)
  (int *mode, double *omega, double *searchinterval, int *order, 
   int *nresonances, int *nfrequencies, double *time, int *l, int *m, double *damping,
   STR_F_TYPE ffile STR_ARG1)
{
  int c;

  /* This line would be present if we wanted to make the omega a 
     mandatory argument. But for the moment I think it should not be mandatory.
     if(argc==1) oscillator_strength_help(); */

  while (1) {
    c = getopt(argc, argv, "hm:s:k:O:r:n:t:d:f:");
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

    case 'd':
      *damping = (double)atof(optarg);
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


/* FUNCTIONS TO BE USED BY THE PROGRAM oct-harmonic-spectrum */


void harmonic_spectrum_help(){
  printf("Usage: oct-harmonic-spectrum [OPTIONS] \n");
  printf("\n");
  printf("Options:\n");
  printf("  -h              Print this help and exits.\n");
  printf("  -w <freq>       Specifies the fundamental frequency.\n");
  printf("  -p <pol>        Specifies the direction of the light polarization.\n");
  printf("                  The oct-harmonic-spectrum utility program needs to know\n");
  printf("                  the direction along which the emission radiation is\n");
  printf("                  considered to be polarized. It may be linearly polarized\n");
  printf("                  or circularly polarized. The valid options are:\n");
  printf("                     'x' : Linearly polarized field in the x direction.\n");
  printf("                     'y' : Linearly polarized field in the x direction.\n");
  printf("                     'z' : Linearly polarized field in the x direction.\n");
  printf("                     '+' : Circularly polarized field, counterclockwise.\n");
  printf("                     '-' : Circularly polarized field, clockwise.\n");
  printf("                  The default is 'x'\n");
  printf("  -m <mode>       Whether the harmonic spectrum is computed by taking the\n");
  printf("                  second derivative of the dipole moment numerically, or by\n");
  printf("                  making use of the acceleration operator, stored in the\n:");
  printf("                  'acceleration' file. The options are:\n");
  printf("                     '1' : use the dipole, take second derivative numerically.\n");
  printf("                     '2' : use the acceleration file.\n");
  printf("                  The default is '1'\n");
  exit(-1);
}

void FC_FUNC_(getopt_harmonic_spectrum, GETOPT_HARMONIC_SPECTRUM)
     (double *w0, int *m, STR_F_TYPE pol STR_ARG1)
{
  int c;

  while (1) {
    c = getopt(argc, argv, "hw:p:m:");
    if (c == -1) break;
    switch (c) {

    case 'h':
      harmonic_spectrum_help();
      break;

    case 'w':
      *w0 = (double)atof(optarg);
      break;

    case 'p':
      TO_F_STR1(optarg, pol);
      break;

    case 'm':
      *m = (int)atoi(optarg);
      break;

    case '?':
      harmonic_spectrum_help();
    }
  }

}


/***************************************************************/
