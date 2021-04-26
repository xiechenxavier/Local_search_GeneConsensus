test:SolveGene.o testGene.o
		g++ SolveGene.o testGene.o -o test

SolveGene.o:SolveGene.cpp
		g++ -c SolveGene.cpp

testGene.o:testGene.cpp
		g++ -c testGene.cpp
