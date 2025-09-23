EXE=gene_network_matpat.exe
CPP=gene_network_matpat.cpp main.cpp individual.cpp
HPP=gene_network_matpat.hpp parameters.hpp individual.hpp
CXX=g++
CXXFLAGS=-Wall -O3 -Wconversion -std=c++20

$(EXE) : $(CPP) $(HPP)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(CPP)

clean : 
	rm -rf $(EXE)
