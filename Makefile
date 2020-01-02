CXX=g++
CXXFLAGS=-Iinclude -std=c++11 -pthread -O3 -c
LIBS=-larmadillo -lgsl -lblas

all: mkfolder tests physics

clean:
	rm -r bin

# test code

physics: fits physics/HVQCD_mass_ratios.exe

fits: physics/fitSingletAxialVectorMesons.exe physics/fittingBackground.exe physics/fitBackground_kU1_not_fixed.exe

physics/fitSingletAxialVectorMesons.exe: bin/physics/fitSingletAxialVectorMesons.o bin/background.o bin/HolographicVQCD.o bin/schrodinger/common.o bin/schrodinger/solvspec.o bin/schrodinger/numerov.o bin/schrodinger/chebspec.o bin/schrodinger/schrodinger.o
	$(CXX) -o $@ $^ $(LIBS)

physics/HVQCD_mass_ratios.exe: bin/physics/HVQCD_mass_ratios.o bin/background.o bin/HolographicVQCD.o bin/schrodinger/common.o bin/schrodinger/solvspec.o bin/schrodinger/numerov.o bin/schrodinger/chebspec.o bin/schrodinger/schrodinger.o
	$(CXX) -o $@ $^ $(LIBS)

physics/fitBackground_kU1_not_fixed.exe: bin/physics/fitBackground_kU1_not_fixed.o bin/background.o bin/HolographicVQCD.o bin/schrodinger/common.o bin/schrodinger/solvspec.o bin/schrodinger/numerov.o bin/schrodinger/chebspec.o bin/schrodinger/schrodinger.o
	$(CXX) -o $@ $^ $(LIBS)

physics/fittingBackground.exe: bin/physics/fittingBackground.o bin/background.o bin/HolographicVQCD.o bin/schrodinger/common.o bin/schrodinger/solvspec.o bin/schrodinger/numerov.o bin/schrodinger/chebspec.o bin/schrodinger/schrodinger.o
	$(CXX) -o $@ $^ $(LIBS)

tests: bin/test_schrodinger.exe bin/test_HVQCD.exe

bin/test_HVQCD.exe: bin/HolographicVQCD.o bin/background.o bin/tests/test_HVQCD.o bin/schrodinger/common.o bin/schrodinger/solvspec.o bin/schrodinger/numerov.o bin/schrodinger/chebspec.o bin/schrodinger/schrodinger.o
	$(CXX) -o $@ $^ $(LIBS)


bin/test_schrodinger.exe: bin/tests/test_schrodinger.o bin/schrodinger/common.o bin/schrodinger/solvspec.o bin/schrodinger/chebspec.o bin/schrodinger/numerov.o bin/schrodinger/schrodinger.o
	$(CXX) -o $@ $^ $(LIBS)

# test files

bin/tests/%.o: tests/%.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS)

# src files

# C++ files
bin/%.o: src/%.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS)

# Physics files
bin/physics/%.o: physics/%.cpp
	$(CXX) -o $@ $^ $(CXXFLAGS)

mkfolder: bin bin/schrodinger bin/tests bin/physics

bin/tests:
	mkdir $@

bin/physics:
	mkdir $@

bin/schrodinger:
	mkdir $@
bin:
	mkdir $@