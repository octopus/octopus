
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

#include <gsl/gsl_complex_math.h>


/* ------------------------------------------------------ */
gsl_complex gsl_complex_step_real (double a)
{        
  gsl_complex z;
	
  if (a < 0)
	{
		GSL_SET_COMPLEX (&z, 0, 0);
	}
  else
	{
		GSL_SET_COMPLEX (&z, 1, 0);
	}

  return z;
}


/* ------------------------------------------------------ */
gsl_complex gsl_complex_min_real (gsl_complex a, gsl_complex b)
{
  gsl_complex z;
  double min;
	
	/* just consider real parts */
  min = GSL_REAL(a) < GSL_REAL(b) ? GSL_REAL(a) : GSL_REAL(b);
  GSL_SET_COMPLEX (&z, min, 0);
	
  return z;
}


/* ------------------------------------------------------ */
gsl_complex gsl_complex_max_real (gsl_complex a, gsl_complex b)
{
  gsl_complex z;
  double max;
	
	/* just consider real parts */
  max = GSL_REAL(a) > GSL_REAL(b) ? GSL_REAL(a) : GSL_REAL(b);
  GSL_SET_COMPLEX (&z, max, 0);
	
  return z;
}

