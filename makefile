CC=gcc
DEPS=lib/lib.h
CFLAGS= -lm


FILE=metropolis

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(FILE): examples/$(FILE).o lib/calcs.o lib/info.o lib/moves.o lib/system.o
	$(CC) -o $(FILE) examples/$(FILE).o lib/calcs.o lib/info.o lib/moves.o lib/system.o $(CFLAGS)
