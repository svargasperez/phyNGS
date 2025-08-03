.PHONY: all clean

CCMPI = mpic++
CCFLAGS = -O3 -m64  -Wall -pedantic
CCFLAGS += -fopenmp -std=c++17
C_PROG = phyNGSC
D_PROG = phyNGSD
I_PROG = incompresso
M_PROG = main
T_PROG = testing
OBJSC = bit_stream.o huffman.o tasks.o
OBJSD = bit_stream.o huffman.o tasks.o
OBJSI = bit_stream.o huffman.o tasks.o
OBJSM = bit_stream.o huffman.o tasks.o debug.o phyNGSC.o phyNGSD.o incompresso.o
OBJST = 
OBJS_TRIM = defs.o

.cpp.o:
	$(CCMPI) $(CCFLAGS) -c $< -o $@

# These targets removed because they don't contain main() functions
# $(C_PROG): phyNGSC.o $(OBJSC)
#	$(CCMPI) $(CCFLAGS) -o $(C_PROG) phyNGSC.o $(OBJSC)

# $(D_PROG): phyNGSD.o $(OBJSD)
#	$(CCMPI) $(CCFLAGS) -o $(D_PROG) phyNGSD.o $(OBJSD)

# $(I_PROG): incompresso.o $(OBJSI)
#	$(CCMPI) $(CCFLAGS) -o $(I_PROG) phyNGSD.o $(OBJSI)

$(M_PROG): main.o $(OBJSM)
	$(CCMPI) $(CCFLAGS) -o $(M_PROG) main.o $(OBJSM)

$(T_PROG): testing.o $(OBJST)
	$(CCMPI) $(CCFLAGS) -o $(T_PROG) testing.o $(OBJST)

temp_trim.cpp:
	$(CCMPI) $(CCFLAGS) -o temp_trim temp_trim.o $(OBJS_TRIM)

all: $(M_PROG) $(T_PROG)

clean:
	-rm *.o
	-rm $(M_PROG) $(T_PROG)
	-rm *.dif ctest*.ngsc dtest*.fastq