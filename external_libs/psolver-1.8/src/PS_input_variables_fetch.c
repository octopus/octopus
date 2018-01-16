#include <config.h>

#define _GNU_SOURCE

#include <fcntl.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char ps_input_vars_def[] = 
#include "PS_input_variables_definition-inc.h"
  ;

// function to recuperate the size of the definition of the input variables
void FC_FUNC(getpsinputdefsize, GETPSINPUTDEFSIZE)(int *pt_len) {
  *pt_len = strlen(ps_input_vars_def);
}

void FC_FUNC(getpsinputdef, GETPSINPUTDEF)(char *to)
{
  memcpy(to,ps_input_vars_def, sizeof(char) * strlen(ps_input_vars_def));
}

/*example of the symbol to export for the database*/
void FC_FUNC_(get_ps_inputvars, GET_PS_INPUTVARS)(char* db_ptr,int* db_len)
{
  if (*db_len==0) 
    {
    *db_len=strlen(ps_input_vars_def);
    return;
    }
  else
    {
      memcpy(db_ptr,ps_input_vars_def, sizeof(char) * (*db_len));
      /*printf("\n test address = %p %d; \n", NULL,*db_len);*/
    }
}
