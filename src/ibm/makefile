EXE=gene_network_matpat.exe
CPP=gene_network_matpat.cpp main.cpp
HPP=gene_network_matpat.hpp parameters.hpp
CXX=g++
CXXFLAGS=-Wall -O3 -std=c++20

$(EXE) : $(CPP) $(HPP)
	$(CXX) $(CXXFLAGS) -o $(EXE) $(CPP)

clean : 
	rm -rf $(EXE)
