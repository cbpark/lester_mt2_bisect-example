CXXFLAGS := -g -O2 -Wall -I. $(CXXFLAGS)
LDFLAGS  := -lm $(LDFLAGS)

.PHONY: all clean

all: mt2

mt2: mt2.o maos.o
	$(CXX) $(LDFLAGS) -o $@ $^

maos.o: maos.cc

clean::
	rm -f mt2 mt2.o maos.o
