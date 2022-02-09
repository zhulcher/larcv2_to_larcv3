
#################################################################
# Here are the variables you may need to change:

CC=/usr/bin/g++
H5_INCDIR=/usr/include/hdf5/serial/
H5_LIBDIR=/usr/lib/x86_64-linux-gnu/hdf5/serial/

#################################################################



# Python path does not need to be set for larcv3

LARCV3_INCDIR=$(shell PYTHONPATH="" python -c "import larcv; print(larcv.get_includes())")
LARCV3_LIBDIR=$(shell PYTHONPATH="" python -c "import larcv; print(larcv.get_lib_dir())")

LARCV3_INCDIR=/usr/local/lib/python3.6/dist-packages/larcv-3.4.1-py3.6-linux-x86_64.egg/larcv/include
pybind_incdir=/app/larcv3/src/pybind11_json/include
json_incdir=/app/larcv3/src/json/include
pybind2_incdir=/app/larcv3/src/pybind11/include


# export PYTHONPATH=$PYTHONPATH_BACKUP

# export ROOT_FLAGS=$(root-config --cflags)
# export ROOT_LIBS=$(root-config --libs)



CFLAGS=-I. -I${LARCV3_INCDIR} -I${LARCV_INCDIR} -I${H5_INCDIR} -I${pybind_incdir} -I${pybind2_incdir} -I${json_incdir} -I/usr/include/python3.6m -I/usr/include/python3.6m -I $(shell root-config --cflags) -g

LDFLAGS=$(shell root-config --libs) ${ROOT_LIBS} -L ${LARCV_LIBDIR} \
-llarcv -L ${LARCV3_LIBDIR} -llarcv3 \
-L${H5_LIBDIR} -lhdf5 -lhdf5_cpp

DEPS = larcv2_to_larcv3.h
OBJ = larcv2_to_larcv3.o
EXEC = larcv2_to_larcv3.so

all: $(EXEC)


%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(EXEC): $(OBJ)
	$(CC) -shared -o $@ $^ $(LDFLAGS)

clean:
	rm *.o
	rm $(EXEC)