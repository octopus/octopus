/*
 Copyright (C) 2015 X. Andrade

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
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.

 $Id:$
*/

#include <stdio.h>
#include <config.h>
#include <fortran_types.h>
#include "string_f.h"
#include <string.h>

#define XML_FILE_DEBUG

fint FC_FUNC_(xml_file_init, XML_FILE_INIT)(FILE ** xml_file, STR_F_TYPE fname STR_ARG1)
{

  char *fname_c;
  TO_C_STR1(fname, fname_c);

#ifdef XML_FILE_DEBUG
  printf("Opening file \"%s\": ", fname_c);
#endif

  *xml_file = fopen(fname_c, "r");


#ifdef XML_FILE_DEBUG
  if(*xml_file != NULL){ 
    printf("success\n");
  } else {
    printf("failed\n");
  }
#endif

  free(fname_c);

  if(*xml_file == NULL){
    return 1;
  } else {
    return 0;
  }

}


void seek_tag(FILE ** xml_file, const char * tag){
  fpos_t startpos;
  char * res;
  char buffer[1000];

  fseek(*xml_file, 0, SEEK_SET);

  fgetpos(*xml_file, &startpos);

  while(fgets(buffer, sizeof(buffer), *xml_file)){

#ifdef XML_FILE_DEBUG
    //    printf("line: %s\n", buffer);
#endif
    
    /* check for the tag */
    res = strstr(buffer, tag);

    /* the string was found */
    if(res != NULL) {
    
#ifdef XML_FILE_DEBUG
    printf("Tag was found in line: %s", res);
#endif

    /* find where the tag is closed*/
    while(res[0] != '>') res++;

    /* go back to the beginning of the line */
    fsetpos(*xml_file, &startpos);

    /* and move to the end of the tag */
    fseek(*xml_file, res - buffer + 1, SEEK_CUR);

    break;
  }

    fgetpos(*xml_file, &startpos);
  }
}

fint FC_FUNC_(xml_file_read_integer_low, XML_FILE_READ_INTEGER_LOW)(FILE ** xml_file, STR_F_TYPE tag_f, fint * value STR_ARG1){
  char * tag;

  TO_C_STR1(tag_f, tag);

#ifdef XML_FILE_DEBUG
  printf("Reading tag \"%s>\" of type integer\n", tag);
#endif

  seek_tag(xml_file, tag);

  fscanf(*xml_file, "%d", value);  

#ifdef XML_FILE_DEBUG
  printf("Got value: %d\n", *value);
#endif

  free(tag);

  return 0;
}

fint FC_FUNC_(xml_file_read_double_low, XML_FILE_READ_DOUBLE_LOW)(FILE ** xml_file, STR_F_TYPE tag_f, double * value STR_ARG1){
  char * tag;

  TO_C_STR1(tag_f, tag);

#ifdef XML_FILE_DEBUG
  printf("Reading tag \"%s>\" of type double\n", tag);
#endif

  seek_tag(xml_file, tag);

  fscanf(*xml_file, "%lf", value);  

#ifdef XML_FILE_DEBUG
  printf("Got value: %lf\n", *value);
#endif

  free(tag);

  return 0;
}

void FC_FUNC_(xml_file_end, XML_FILE_END)(FILE ** xml_file)
{
  fclose(*xml_file);
}
