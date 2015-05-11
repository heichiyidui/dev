# The Makefile for bpnn
# test class bpnn


CC      = g++
CCFLAG  = -Wall -O 
EXEC    = yaspin
LIB     = bpnn.cpp hmm.cpp hmmvalues.cpp
SRC     = main.cpp
OBJS    = bpnn.o hmm.o hmmvalues.o 

$(EXEC)  : $(OBJS) $(SRC) 
	$(CC) $(CCFLAG) -o $(EXEC) $(SRC) $(OBJS)  

$(OBJS)  : $(LIB)
	$(CC) $(CCFLAG) -c $(LIB)  

clean    :
	rm $(OBJS)
	rm $(EXEC)


.PHONY   :clean 
