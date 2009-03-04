#include <stdio.h>
#include "ruby.h"
#include "spglib.h"

VALUE Example1 = Qnil;
void Init_example1(void);

VALUE method_example1(VALUE self, VALUE r_size, VALUE r_lattice, VALUE r_position, VALUE r_types);

void Init_example1(void)
{
     Example1 = rb_define_module("Example1");
     rb_define_method(Example1, "example1", method_example1, 4);
}

VALUE method_example1(VALUE self, VALUE r_size, VALUE r_lattice, VALUE r_position, VALUE r_types)
{
     int i, j, size;
     double lattice[3][3];
     
     size = NUM2INT(r_size);
     
     double position[size][3];
     int types[size];
     char symbol[21];
     char output[27];

     for (i=0; i<size; i++)
	  for (j=0; j<3; j++) {
	       position[i][j] =
		    NUM2DBL(rb_ary_entry(rb_ary_entry(r_position, i), j));
	       types[i] = NUM2DBL(rb_ary_entry(r_types, i));
	  }

     for (i=0; i<3; i++)
	  for (j=0; j<3; j++)
	       lattice[i][j] =
		    NUM2DBL(rb_ary_entry(rb_ary_entry(r_lattice, i), j));

     i = spg_get_international(symbol, lattice, position, types, size, 1e-5);
     sprintf(output, "%s (%d)", symbol, i);

     return rb_str_new2(output);
}
