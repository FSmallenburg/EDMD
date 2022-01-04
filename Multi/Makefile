# define the C compiler to use
CC=gcc

# define any compile-time flags
CFLAGS=-Wall -Ofast -funroll-loops -march=native


# define any libraries to link into executable:
LIBS=-lm


md: mdMulti.c  mdMulti.h
	$(CC) $(CFLAGS) mdMulti.c -o md  $(LIBS)
