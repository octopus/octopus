/*
  Copyright (C) 2009 X. Andrade

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

#ifdef HAVE_PAPI
#include <papi.h>

#define NUM_EVENTS 1
static int papi_available;

void FC_FUNC_(papi_init, PAPI_INIT)(void){
  int events[NUM_EVENTS] = {PAPI_FP_OPS};

  if(PAPI_start_counters(events, NUM_EVENTS) == PAPI_OK){
    papi_available = 1;
  } else {
    papi_available = 0;
  }
}

void FC_FUNC_(papi_end, PAPI_END)(void){
  long_long values[NUM_EVENTS];

  if(papi_available) PAPI_stop_counters(values, NUM_EVENTS);
}

void FC_FUNC_(papi_get_count_and_reset, PAPI_GET_COUNT_AND_RESET)(double * fp){
  long_long values[NUM_EVENTS];

  if(papi_available) {
    PAPI_read_counters(values, NUM_EVENTS);
    fp[0] = (double) values[0];
  } else {
    fp[0] = 0.0;
  }
}

#endif
