CC = gcc
CFLAGS = -Wall -O3 -fopenmp
LDFLAGS = -lgsl -lgslcblas -lm -fopenmp

USE_KROUPA ?= 0

ifeq ($(USE_KROUPA),1)
    CFLAGS += -DUSE_KROUPA_IMF
    OBJ = main.o stellar_dynamics_kroupa.o
else
    OBJ = main.o stellar_dynamics_uniform.o
endif

EXEC = stellar_simulation

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Compile stellar_dynamics differently based on USE_KROUPA
stellar_dynamics_kroupa.o: stellar_dynamics.c
	$(CC) $(CFLAGS) -DUSE_KROUPA_IMF -c $< -o $@

stellar_dynamics_uniform.o: stellar_dynamics.c
	$(CC) $(CFLAGS) -c $< -o $@

main.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o $(EXEC)