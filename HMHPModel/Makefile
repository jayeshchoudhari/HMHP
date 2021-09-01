mainHMHP: hmhpModel.o dataIO.o
	g++ -O3 hmhpModel.o dataIO.o -o mainHMHP

hmhpModel.o: hmhpModel.cpp ./include/dataIO.h  ./include/namespace.h ./include/namespace.h ./include/utilities.h ./include/initialization.h
	g++ -O3 -c hmhpModel.cpp

dataIO.o: ./include/dataIO.cpp  ./include/dataIO.h  ./include/namespace.h
	g++ -O3 -c ./include/dataIO.cpp

clean:
	rm *.o mainHMHP