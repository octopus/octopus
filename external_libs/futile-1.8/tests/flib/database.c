#include <config.h>
#define _GNU_SOURCE
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static const char database_arr[] =
"### In this file we define some database information which is then \n"
"## inserted in the executable via the makefile\n"
"---\n"
"C60:\n"
"  JOURNAL_REF: Nature 318, 6042, 162-163 (1985)\n"
"  DESCRIPTION: Seminal Paper about C60 cage\n"
"  BIBTEX_REF: |\n"
"   @article{kroto1985,\n"
"   title={C60: buckminsterfullerene},\n"
"    author={Kroto, Harold W and Heath, James R and O'Brien, Sean C and Curl, Robert F and Smalley, Richard E and others},\n"
"   journal={Nature},\n"
"   volume={318},\n"
"   number={6042},\n"
"   pages={162--163},\n"
"   year={1985},\n"
"   publisher={London}\n"
"   }\n"
"PS_FBC:\n"
"  DESCRIPTION: Paper of the Free BC Poisson Solver\n"  ;
void FC_FUNC_(get_database, GET_DATABASE)(char* db_ptr,int* db_len)
{
  if (*db_len==0)
    {
    *db_len=strlen(database_arr);
    return;
    }
  else
    {
      memcpy(db_ptr,database_arr, sizeof(char) * (*db_len));
    }
}
