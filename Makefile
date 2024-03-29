OSTYPE := $(shell uname)
ifeq ($(OSTYPE),Linux)
CYGWIN=
else
CYGWIN=-Wl,--enable-auto-import
endif

GCC=g++
GCCFLAGS=-O2 -Wall -Wextra -ansi -pedantic -Wold-style-cast -Woverloaded-virtual -Wsign-promo  -Wctor-dtor-privacy -Wnon-virtual-dtor -Wreorder
DEFINE=-DINLINE_VARIABLE -DINLINE_CONSTRAINT_GRAPH -DINLINE_CONSTRAINT -DINLINE_CSP 

MSC=cl
MSCFLAGS=/EHa /W4 /Za /Zc:forScope /nologo /D_CRT_SECURE_NO_DEPRECATE /D"_SECURE_SCL 0" /O2i /GL
MSCDEFINE=/DINLINE_VARIABLE /DINLINE_CONSTRAINT_GRAPH /DINLINE_CONSTRAINT /DINLINE_CSP 

VALGRIND_OPTIONS=-q --leak-check=full
DIFF_OPTIONS=-y --strip-trailing-cr --suppress-common-lines

#everything is templetized
OBJECTS0=

DRIVER0=main.cpp

OSTYPE := $(shell uname)
ifeq ($(OSTYPE),Linux)
CYGWIN=
else
CYGWIN=-Wl,--enable-auto-import
endif

.PHONY: run
run: example simple exampleFA simpleFA queen-28-dfs queen-28-dfsFA queen-100-fc ms5-fc msbc5-fc msbc5-dfs ms6-fc msbc6-fc

example: 
	$(GCC) $(DRIVER0) -DEXAMPLE -DDFS $(CYGWIN) $(OBJECTS0) $(GCCFLAGS) $(DEFINE) -o $@.exe #ARC,DFS
simple:
	$(GCC) $(DRIVER0) -DSIMPLE -DDFS $(CYGWIN) $(OBJECTS0) $(GCCFLAGS) $(DEFINE) -o $@.exe #ARC,DFS
exampleFA:
	$(GCC) $(DRIVER0) -DEXAMPLE -DDFS -DFIRST_AVAILABLE $(CYGWIN) $(OBJECTS0) $(GCCFLAGS) $(DEFINE) -o $@.exe #ARC,DFS
simpleFA:
	$(GCC) $(DRIVER0) -DSIMPLE -DDFS -DFIRST_AVAILABLE $(CYGWIN) $(OBJECTS0) $(GCCFLAGS) $(DEFINE) -o $@.exe #ARC,DFS
	
queen-28-dfs:
	$(GCC) $(DRIVER0) -DQUEEN -DSIZE=28 -DDFS $(CYGWIN) $(OBJECTS0) $(GCCFLAGS) $(DEFINE) -o $@.exe #ARC,DFS
queen-28-dfsFA:
	$(GCC) $(DRIVER0) -DQUEEN -DSIZE=28 -DDFS -DFIRST_AVAILABLE $(CYGWIN) $(OBJECTS0) $(GCCFLAGS) $(DEFINE) -o $@.exe #ARC,DFS
queen-100-fc:
	$(GCC) $(DRIVER0) -DQUEEN -DSIZE=100 -DFC $(CYGWIN) $(OBJECTS0) $(GCCFLAGS) $(DEFINE) -o $@.exe #ARC,DFS
queen-100-arc:
	$(GCC) $(DRIVER0) -DQUEEN -DSIZE=100 -DARC $(CYGWIN) $(OBJECTS0) $(GCCFLAGS) $(DEFINE) -o $@.exe #ARC,DFS

ms5-fc:
	$(GCC) $(DRIVER0) -DMS   -DSIZE=5 -DFC $(CYGWIN) $(OBJECTS0) $(GCCFLAGS) $(DEFINE) -o $@.exe #ARC,DFS
msbc5-fc:
	$(GCC) $(DRIVER0) -DMSBC -DSIZE=5 -DFC $(CYGWIN) $(OBJECTS0) $(GCCFLAGS) $(DEFINE) -o $@.exe #ARC,DFS
msbc5-dfs:
	$(GCC) $(DRIVER0) -DMSBC -DSIZE=5 -DDFS -DFScount $(CYGWIN) $(OBJECTS0) $(GCCFLAGS) $(DEFINE) -o $@.exe #ARC,DFS
ms5-arc:
	$(GCC) $(DRIVER0) -DMS   -DSIZE=5 -DARC $(CYGWIN) $(OBJECTS0) $(GCCFLAGS) $(DEFINE) -o $@.exe #ARC,DFS
msbc5-arc:
	$(GCC) $(DRIVER0) -DMSBC -DSIZE=5 -DARC $(CYGWIN) $(OBJECTS0) $(GCCFLAGS) $(DEFINE) -o $@.exe #ARC,DFS

ms6-fc:
	$(GCC) $(DRIVER0) -DMS   -DSIZE=6 -DFC $(CYGWIN) $(OBJECTS0) $(GCCFLAGS) $(DEFINE) -o $@.exe #ARC,DFS
msbc6-fc:
	$(GCC) $(DRIVER0) -DMSBC -DSIZE=6 -DFC $(CYGWIN) $(OBJECTS0) $(GCCFLAGS) $(DEFINE) -o $@.exe #ARC,DFS

#MS compiler
msc-example:
	$(MSC) $(DRIVER0) -DEXAMPLE -DDFS  $(OBJECTS0) $(MSCFLAGS) $(MSCDEFINE) /Fe$@.exe #ARC,DFS
msc-simple:
	$(MSC) $(DRIVER0) -DSIMPLE -DDFS  $(OBJECTS0) $(MSCFLAGS) $(MSCDEFINE) /Fe$@.exe #ARC,DFS

msc-queen-28-dfs:
	$(MSC) $(DRIVER0) -DQUEEN -DSIZE=28 -DDFS  $(OBJECTS0) $(MSCFLAGS) $(MSCDEFINE) /Fe$@.exe #ARC,DFS
msc-queen-100-fc:
	$(MSC) $(DRIVER0) -DQUEEN -DSIZE=100 -DFC  $(OBJECTS0) $(MSCFLAGS) $(MSCDEFINE) /Fe$@.exe #ARC,DFS
msc-queen-100-arc:
	$(MSC) $(DRIVER0) -DQUEEN -DSIZE=100 -DARC  $(OBJECTS0) $(MSCFLAGS) $(MSCDEFINE) /Fe$@.exe #ARC,DFS

msc-ms5-fc:
	$(MSC) $(DRIVER0) -DMS   -DSIZE=5 -DFC  $(OBJECTS0) $(MSCFLAGS) $(MSCDEFINE) /Fe$@.exe #ARC,DFS
msc-msbc5-fc:
	$(MSC) $(DRIVER0) -DMSBC -DSIZE=5 -DFC  $(OBJECTS0) $(MSCFLAGS) $(MSCDEFINE) /Fe$@.exe #ARC,DFS
msc-ms5-arc:
	$(MSC) $(DRIVER0) -DMS   -DSIZE=5 -DARC  $(OBJECTS0) $(MSCFLAGS) $(MSCDEFINE) /Fe$@.exe #ARC,DFS
msc-msbc5-arc:
	$(MSC) $(DRIVER0) -DMSBC -DSIZE=5 -DARC  $(OBJECTS0) $(MSCFLAGS) $(MSCDEFINE) /Fe$@.exe #ARC,DFS

msc-ms6-fc:
	$(MSC) $(DRIVER0) -DMS   -DSIZE=6 -DFC  $(OBJECTS0) $(MSCFLAGS) $(MSCDEFINE) /Fe$@.exe #ARC,DFS
msc-msbc6-fc:
	$(MSC) $(DRIVER0) -DMSBC -DSIZE=6 -DFC  $(OBJECTS0) $(MSCFLAGS) $(MSCDEFINE) /Fe$@.exe #ARC,DFS

clean:
	rm -f *.exe *.obj *.o
