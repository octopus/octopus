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

 $Id: operate.c 2146 2006-05-23 17:36:00Z xavier $
*/

#include <stdio.h>

#include <config.h>

#include "beak.h"

#define _XOPEN_SOURCE 600
#include <stdlib.h>

int FC_FUNC_(op_is_available, OP_IS_AVAILABLE)
  (const int * opid, const int * type){
  int result = 1;
  
#ifndef HAVE_VEC
  if( *opid == OP_VEC ) result = 0;
#endif

#ifndef HAVE_AS
  if( *opid == OP_AS ) result = 0;
#endif

#ifdef SINGLE_PRECISION
  if( *opid == OP_AS ) result = 0;
  if( *opid == OP_VEC ) result = 0
#endif
  return result;
}
