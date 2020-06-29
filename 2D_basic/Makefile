CC = g++
CFLAGS = -Wall -O2 -lm

OBJS = model.o bc.o flux_solver.o

all : run

run : main.o $(OBJS)
	$(CC) main.o $(OBJS) -o a.out $(CFLAGS)

%.o : %.cpp
	$(CC) -c $< $(CFLAGS)

clean : 
	rm -f *.out *.o data/*.dat
