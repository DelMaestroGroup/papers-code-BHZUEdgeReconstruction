EXENAME = ED_exe

CXX = g++
CPPFLAGS = -std=c++11
LDFLAGS = -llapack -lblas

#CPPFLAGS += -Ifunctions.h -IHamiltonian.h
#CPPFLAGS += -g3
CPPFLAGS += -O3

STRIP_COMMAND = true

$(EXENAME): clean main.o
				$(CXX) $(CPPFLAGS) -o $(EXENAME)  main.o $(LDFLAGS)
				$(STRIP_COMMAND) $(EXENAME)

all: $(EXENAME)

clean:
				rm -f $(EXENAME) *.o