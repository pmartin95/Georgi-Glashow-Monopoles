CC=g++
SIM_OBJECTS=sim.o georgiglashow.o observables.o datacollection.o boundary_conditions.o leapfrog.o utilities.o metropolis_hastingsMC.o statistics.o monopole.o
LATTICE_OBJECTS=lattice.o plattice.o lattice_site.o lattice_base.o
OBJECTS=main.o $(LATTICE_OBJECTS) rand.o $(SIM_OBJECTS) stopwatch.o
INC= -I ./../eigen -I ./../includes -lm -g -fopenmp -O3

main.exe: $(OBJECTS)
	$(CC) $(INC) -o main.exe $(OBJECTS)

%.o : %.cpp
	$(CC) $(INC) -c -o $@ $<

clean:
	rm -f *.o *.exe
