CXXFLAGS := -g -O2 -Wall -I.

.PHONY: all clean

all: mt2

mt2: mt2.o
	$(CXX) $(LDFLAGS) -o $@ $^

clean::
	rm -f mt2 mt2.o
