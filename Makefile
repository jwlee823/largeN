BINARY = largeN

CPPFLAGS := -g -Wall -O3
GSLFLAGS = -lgsl -lgslcblas -lm


GCC = g++

CPPSRC = main.cpp

# get the list of objective files from that of the soure files
CPPOBJ := $(patsubst %.cpp,%.o,$(notdir $(CPPSRC)))

HEADER = global.h

all : largeN

largeN : $(CPPOBJ) $(HEADER)
	${GCC} ${CPPFLAGS} $(GSLFLAGS) -o largeN $(CPPOBJ)


$(CPPOBJ) : %.o : %.cpp $(HEADER)
	$(GCC) $(CPPFLAGS) -o $@ -c $<

clean :
	rm -rf $(CPPOBJ) $(BINARY)


