CC=g++
OBJECTS=main.o lattice.o rand.o sim.o georgiglashow.o observables.o datacollection.o boundary_conditions.o leapfrog.o utilities.o metropolis_hastingsMC.o stopwatch.o 
INC= -I ./../eigen -I ./../includes -lm -g -fopenmp -O3

main.exe: $(OBJECTS)
	$(CC) $(INC) -o main.exe $(OBJECTS)

%.o : %.cpp
	$(CC) $(INC) -c -o $@ $<

clean:
	rm -f *.o *.exe
