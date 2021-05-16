#website that taught me to make this make file: http://nuclear.mutantstargoat.com/articles/make/
src = $(wildcard *.cpp)
obj = $(src:.cpp=.o)
dep = $(obj:.o=.d)  # one dependency file for each source
prefix = SKMC
suffix =
conda_path =

CXXFLAGS = -MMD -O3 -I$(conda_path)include -I$(PWD) -L$(conda_path)lib# option to generate a .d file during compilation
LDFLAGS =
CC = gcc
C++ = g++
$(prefix)$(suffix): $(obj)
	$(CC) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

-include $(dep)   # include all dep files in the makefile

.PHONY: clean
clean:
	rm -f $(obj) $(prefix)$(suffix)

.PHONY: cleandep
cleandep:
	rm -f $(dep)
