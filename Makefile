CC=g++
OBJECTS= main.o lattice.o rand.o
INC= -I ./eigen -I ./unsupported -lm

main.exe: $(OBJECTS)
	$(CC) $(INC) -o main.exe $(OBJECTS)

%.o : %.cpp
	$(CC) $(INC) -c -o $@ $<

clean:
	rm *.o *.exe
