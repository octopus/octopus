SRC=recipes.c math.c
OBJ= $(SRC:.c=.o)
CPPFLAGS=
CFLAGS=-Wall

all: $(OBJ)
	$(AR) cr $(LIB_DIR)/liboct.a $(OBJ)

clean:
	rm -f *~ *.o
