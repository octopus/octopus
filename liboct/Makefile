SRC=recipes.c ylm.c parse.c parse_exp.c symbols.c varia.c interface.c
OBJ= $(SRC:.c=.o)
CPPFLAGS=
CFLAGS=-Wall -g
LIB_DIR ?= .

all: $(OBJ)
	$(AR) cr $(LIB_DIR)/liboct.a $(OBJ)

dep:
	makedepend -- $(CFLAGS) -- $(SRC)

clean:
	rm -f *~ *.o grammar.c

grammar.c: grammar.y
parse_exp.c: grammar.c

# DO NOT DELETE

recipes.o: /usr/include/stdio.h /usr/include/features.h
recipes.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
recipes.o: /usr/lib/gcc-lib/i386-redhat-linux/2.96/include/stddef.h
recipes.o: /usr/include/bits/types.h /usr/include/libio.h
recipes.o: /usr/include/_G_config.h /usr/include/wchar.h
recipes.o: /usr/include/bits/wchar.h /usr/include/gconv.h
recipes.o: /usr/lib/gcc-lib/i386-redhat-linux/2.96/include/stdarg.h
recipes.o: /usr/include/bits/stdio_lim.h /usr/include/stdlib.h
recipes.o: /usr/include/string.h /usr/include/time.h /usr/include/bits/time.h
recipes.o: config.h
ylm.o: /usr/include/stdlib.h /usr/include/features.h /usr/include/sys/cdefs.h
ylm.o: /usr/include/gnu/stubs.h
ylm.o: /usr/lib/gcc-lib/i386-redhat-linux/2.96/include/stddef.h
ylm.o: /usr/include/math.h /usr/include/bits/huge_val.h
ylm.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
ylm.o: /usr/include/gsl/gsl_sf_legendre.h /usr/include/gsl/gsl_sf_result.h
parse.o: /usr/include/stdio.h /usr/include/features.h
parse.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
parse.o: /usr/lib/gcc-lib/i386-redhat-linux/2.96/include/stddef.h
parse.o: /usr/include/bits/types.h /usr/include/libio.h
parse.o: /usr/include/_G_config.h /usr/include/wchar.h
parse.o: /usr/include/bits/wchar.h /usr/include/gconv.h
parse.o: /usr/lib/gcc-lib/i386-redhat-linux/2.96/include/stdarg.h
parse.o: /usr/include/bits/stdio_lim.h /usr/include/stdlib.h
parse.o: /usr/include/string.h /usr/include/ctype.h /usr/include/endian.h
parse.o: /usr/include/bits/endian.h /usr/include/math.h
parse.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
parse.o: /usr/include/bits/mathcalls.h liboct.h
parse.o: /usr/include/gsl/gsl_complex.h symbols.h
parse_exp.o: /usr/include/stdio.h /usr/include/features.h
parse_exp.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
parse_exp.o: /usr/lib/gcc-lib/i386-redhat-linux/2.96/include/stddef.h
parse_exp.o: /usr/include/bits/types.h /usr/include/libio.h
parse_exp.o: /usr/include/_G_config.h /usr/include/wchar.h
parse_exp.o: /usr/include/bits/wchar.h /usr/include/gconv.h
parse_exp.o: /usr/lib/gcc-lib/i386-redhat-linux/2.96/include/stdarg.h
parse_exp.o: /usr/include/bits/stdio_lim.h /usr/include/string.h
parse_exp.o: /usr/include/stdlib.h /usr/include/ctype.h /usr/include/endian.h
parse_exp.o: /usr/include/bits/endian.h /usr/include/gsl/gsl_complex.h
parse_exp.o: /usr/include/gsl/gsl_complex_math.h liboct.h symbols.h
symbols.o: /usr/include/stdio.h /usr/include/features.h
symbols.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
symbols.o: /usr/lib/gcc-lib/i386-redhat-linux/2.96/include/stddef.h
symbols.o: /usr/include/bits/types.h /usr/include/libio.h
symbols.o: /usr/include/_G_config.h /usr/include/wchar.h
symbols.o: /usr/include/bits/wchar.h /usr/include/gconv.h
symbols.o: /usr/lib/gcc-lib/i386-redhat-linux/2.96/include/stdarg.h
symbols.o: /usr/include/bits/stdio_lim.h /usr/include/stdlib.h
symbols.o: /usr/include/string.h /usr/include/ctype.h /usr/include/endian.h
symbols.o: /usr/include/bits/endian.h /usr/include/math.h
symbols.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
symbols.o: /usr/include/bits/mathcalls.h /usr/include/gsl/gsl_complex_math.h
symbols.o: /usr/include/gsl/gsl_complex.h symbols.h
interface.o: /usr/include/stdio.h /usr/include/features.h
interface.o: /usr/include/sys/cdefs.h /usr/include/gnu/stubs.h
interface.o: /usr/lib/gcc-lib/i386-redhat-linux/2.96/include/stddef.h
interface.o: /usr/include/bits/types.h /usr/include/libio.h
interface.o: /usr/include/_G_config.h /usr/include/wchar.h
interface.o: /usr/include/bits/wchar.h /usr/include/gconv.h
interface.o: /usr/lib/gcc-lib/i386-redhat-linux/2.96/include/stdarg.h
interface.o: /usr/include/bits/stdio_lim.h /usr/include/string.h
interface.o: /usr/include/math.h /usr/include/bits/huge_val.h
interface.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
interface.o: /usr/include/gsl/gsl_sf_erf.h /usr/include/gsl/gsl_sf_result.h
interface.o: /usr/include/gsl/gsl_spline.h /usr/include/stdlib.h
interface.o: /usr/include/gsl/gsl_interp.h config.h liboct.h
interface.o: /usr/include/gsl/gsl_complex.h
