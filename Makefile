CC=g++
#CC=clang++

CFLAGS = -Wall -O2 -march=native -pipe -std=c++11 -g -pedantic

ROOT_CFLAGS = `root-config --cflags`
ROOT_LDFLAGS = `root-config --ldflags`
ROOT_LIBS = `root-config --libs --glibs` -lSpectrum
BIN = "bin"

BOOST_LIBS = -lboost_program_options

all: californio check pileupcheck display gate psd counts
		 @echo "All source compiled!"

californio: data_interface.o functions.o
						$(CC) californio.cc data_interface.o functions.o \
							-o $(BIN)/$@ \
							-I . \
							$(CFLAGS) \
							$(ROOT_CFLAGS) $(ROOT_LDFLAGS) $(ROOT_LIBS) \
							$(BOOST_LIBS)

check: data_interface.o functions.o
						$(CC) check.cc data_interface.o functions.o \
							-o $(BIN)/$@ \
							-I . \
							$(CFLAGS) \
							$(ROOT_CFLAGS) $(ROOT_LDFLAGS) $(ROOT_LIBS)
							
pileupcheck: data_interface.o functions.o
						$(CC) pileupcheck.cc data_interface.o functions.o \
							-o $(BIN)/$@ \
							-I . \
							$(CFLAGS) \
							$(ROOT_CFLAGS) $(ROOT_LDFLAGS) $(ROOT_LIBS)
							
display: data_interface.o functions.o
						$(CC) display.cc data_interface.o functions.o \
							-o $(BIN)/$@ \
							-I . \
							$(CFLAGS) \
							$(ROOT_CFLAGS) $(ROOT_LDFLAGS) $(ROOT_LIBS)

gate: data_interface.o functions.o
						$(CC) gate.cc data_interface.o functions.o \
							-o $(BIN)/$@ \
							-I . \
							$(CFLAGS) \
							$(ROOT_CFLAGS) $(ROOT_LDFLAGS) $(ROOT_LIBS)
	
psd: data_interface.o functions.o
						$(CC) psd.cc data_interface.o functions.o \
							-o $(BIN)/$@ \
							-I . \
							$(CFLAGS) \
							$(ROOT_CFLAGS) $(ROOT_LDFLAGS) $(ROOT_LIBS)

counts: data_interface.o functions.o
						$(CC) counts.cc data_interface.o functions.o \
							-o $(BIN)/$@ \
							-I . \
							$(CFLAGS) \
							$(ROOT_CFLAGS) $(ROOT_LDFLAGS) $(ROOT_LIBS)
	

data_interface.o:
						$(CC) -c data_interface.cc \
							-I . \
							$(CFLAGS) \
							$(ROOT_CFLAGS) $(ROOT_LDFLAGS) $(ROOT_LIBS)

functions.o:
						$(CC) -c functions.cc \
							-I . \
							$(CFLAGS) \
							$(ROOT_CFLAGS) $(ROOT_LDFLAGS) $(ROOT_LIBS)
							

clean: 
				rm *.o $(BIN)/californio $(BIN)/check $(BIN)/pileupcheck $(BIN)/gate $(BIN)/display $(BIN)/psd $(BIN)/counts
