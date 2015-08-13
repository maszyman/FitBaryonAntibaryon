FF	= gfortran
LD	= gfortran
CC	= g++

LDFLAGS	= `root-config --libs`
FFFLAGS = -g -c
CCFLAGS = `root-config --cflags` -g -c -std=c++11

all: fitBaryonAntibaryon

fitBaryonAntibaryon: fitBaryonAntibaryon.o
	$(CC) $^ -o $@ $(LDFLAGS)

%.o: %.f
	$(FF) $^ -o $@ $(FFFLAGS)

%.o: %.C
	$(CC) $^ -o $@ $(CCFLAGS)

clean:
	rm -f *.o fitBaryonAntibaryon
