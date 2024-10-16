CXX = clang++ -fno-color-diagnostics --stdlib=libc++ --std=c++11

# Adjust the include paths according to your directory structure
INCLUDES = -I. -I../../ -I../../lib -I../../algebra -I../../knotkit -I../../knot_parser -I../../rd_parser \
           -I/opt/local/include -I/opt/homebrew/opt/gmp/include

CXXFLAGS = -O2 -g -DHOME="\"`pwd`\"" -Wall -Wno-unused $(INCLUDES)

LDFLAGS = -L/opt/local/lib -L/opt/homebrew/opt/gmp/lib

LIBS = -lgmp -lz

# List the object files from the main repository with adjusted paths
LIB_OBJS = ../../lib/refcount.o \
           ../../lib/lib.o ../../lib/smallbitset.o ../../lib/bitset.o ../../lib/setcommon.o \
           ../../lib/io.o ../../lib/directed_multigraph.o

ALGEBRA_OBJS = ../../algebra/algebra.o ../../algebra/grading.o ../../algebra/polynomial.o

KNOTKIT_OBJS = ../../planar_diagram.o ../../dt_code.o ../../knot_diagram.o ../../cube.o \
               ../../steenrod_square.o ../../spanning_tree_complex.o ../../smoothing.o \
               ../../cobordism.o ../../knot_tables.o ../../sseq.o \
               ../../knot_parser/knot_parser.o ../../knot_parser/knot_scanner.o \
               ../../rd_parser/rd_parser.o ../../rd_parser/rd_scanner.o 

COMMON_OBJS = $(KNOTKIT_OBJS) $(ALGEBRA_OBJS) $(LIB_OBJS)

# Define the target to build your program
all: SparseDiff

SparseDiff: differentials.o $(COMMON_OBJS)
	$(CXX) $(LDFLAGS) -o SparseDiff differentials.o $(COMMON_OBJS) $(LIBS)

differentials.o: differentials.cpp
	$(CXX) -c $(CXXFLAGS) differentials.cpp -o differentials.o
