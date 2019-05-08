CC=clang++

# Python path does not need to be set for larcv3

LARCV3_INCDIR=$(shell PYTHONPATH="" python -c "import larcv; print(larcv.get_includes())")
LARCV3_LIBDIR=$(shell PYTHONPATH="" python -c "import larcv; print(larcv.get_lib_dir())")

# $(LARCV3_INCDIR)


# export PYTHONPATH=$PYTHONPATH_BACKUP

# export ROOT_FLAGS=$(root-config --cflags)
# export ROOT_LIBS=$(root-config --libs)

H5_INCDIR=/opt/local/include/
H5_LIBDIR=/opt/local/lib/

CFLAGS=-I. -I${LARCV3_INCDIR} -I${LARCV_INCDIR} -I${H5_INCDIR} $(shell root-config --cflags) -O3 -g

LDFLAGS=$(shell root-config --libs) ${ROOT_LIBS} -L ${LARCV_LIBDIR} \
-llarcv -L ${LARCV3_LIBDIR} -lbase -ldataformat \
-L${H5_LIBDIR} -lhdf5 -lhdf5_cpp

DEPS = larcv2_to_larcv3.h
OBJ = larcv2_to_larcv3.o main.o 
EXEC = larcv2_to_larcv3

all: $(EXEC)


%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(EXEC): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

clean:
	rm *.o
	rm $(EXEC)