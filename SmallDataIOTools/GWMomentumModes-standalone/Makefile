IDIR = ./Source
CXX = g++
CXXFLAGS = -g -I. -I$(IDIR) -std=c++17 -O3 -march=native -fopenmp -Wall -Wno-sign-compare
EBASE = GWMomentumModes

ODIR = o
LIBS = -lstdc++fs

_DEPS = IntegrationMethod.hpp IntegrationMethodSetup.hpp SmallDataIO.hpp SmallDataIOReader.hpp TimeData.hpp WeylModeData.hpp
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = SmallDataIO.o SmallDataIOReader.o TimeData.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

all: $(EBASE)

$(ODIR)/$(EBASE).o: $(EBASE).cpp $(DEPS)
	@mkdir -p $(ODIR)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(ODIR)/%.o: $(IDIR)/%.cpp $(DEPS)
	@mkdir -p $(ODIR)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(EBASE): $(ODIR)/$(EBASE).o $(OBJ)
	$(CXX) -o $(addsuffix .ex, $@) $^ $(CXXFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -rf *.ex $(ODIR)
