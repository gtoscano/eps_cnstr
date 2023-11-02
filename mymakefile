EXE = epa 

OBJS = main.o \
	misc.o \
	nlp.o

ADDLIBS =

ADDINCFLAGS =

# C++ Compiler command
CXX = g++

# C++ Compiler options
CXXFLAGS = -g    -std=c++2a
#CXXFLAGS = -O3    -std=c++2a

# additional C++ Compiler options for linking
CXXLINKFLAGS =  -Wl,--rpath -Wl,/usr/local/lib

prefix=/usr/local
exec_prefix=${prefix}

# Include directories
INCL = `PKG_CONFIG_PATH=/usr/local/lib/pkgconfig: pkg-config --cflags ipopt redis++ hiredis arrow parquet` $(ADDINCFLAGS)
#INCL = -I${prefix}/include/coin-or -I/usr/local/include/coin-or/hsl -I/usr/local/include/coin-or/mumps     -DIPOPTLIB_BUILD $(ADDINCFLAGS)

# Linker flags
LIBS = `PKG_CONFIG_PATH=/usr/local/lib/pkgconfig: pkg-config --libs ipopt fmt redis++ hiredis arrow parquet`
#LIBS = -L${exec_prefix}/lib -lipopt -L/usr/local/lib -lcoinhsl -lcoinmumps -llapack -lblas    -lm  -ldl

all: $(EXE)

.SUFFIXES: .cpp .o

$(EXE): $(OBJS)
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $(OBJS) $(ADDLIBS) $(LIBS)

clean:
	rm -rf $(EXE) $(OBJS) ipopt.out

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ $<
