# raw makefile
CC = g++
# includes
INC=-I src/ini -I src/spc
CFLAGS=$(INC)
# deps
DEPS = ms1ms2.h iniparser.h dictionary.h spectra.h mparser.h result.h scoring.h
OBJ = main.o ms1ms2.o iniparser.o dictionary.o spectra.o mparser.o result.o scoring.o
VPATH = src/ini:src/spc

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<
%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

speptide: $(OBJ) 
	$(CC) $(CFLAGS) -o $@ $^
	rm -f $(OBJ)
clean:
	rm *.o
