SRC=recipes.c
OBJ= $(SRC:.c=.o)
CPPFLAGS=

all: $(OBJ)
	$(AR) cr $(LIB_DIR)/liboct.a $(OBJ)

clean:
	rm -f *.~ *.o
