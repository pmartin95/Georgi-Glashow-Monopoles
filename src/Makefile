CC=g++
OBJECTS=main.o lattice.o rand.o sim.o
INC= -I ./../eigen -I ./../includes -lm -g -fopenmp -O3

main.exe: $(OBJECTS)
	$(CC) $(INC) -o main.exe $(OBJECTS)

%.o:%.cpp
	$(CC) $(INC) -c -o $@ $<

clean:
	rm *.o *.exe
