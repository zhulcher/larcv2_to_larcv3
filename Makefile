
#################################################################
# Here are the variables you may need to change:


LARCV3_BASEDIR=/app/larcv3/
#H5_INCDIR=/usr/include/hdf5/serial/
#H5_LIBDIR=/usr/lib/x86_64-linux-gnu/hdf5/serial/

#################################################################

CC=g++


# Python path does not need to be set for larcv3


H5_INCDIR=$(shell `h5c++ -show` | tr ' ' '\n' | grep I/ -m1)
H5_LIBDIR=$(shell `h5c++ -show` | tr ' ' '\n' | grep L/ -m1)

LARCV3_INCDIR=$(shell PYTHONPATH="" python3 -c "import larcv; print(larcv.get_includes())")
LARCV3_LIBDIR=$(shell PYTHONPATH="" python3 -c "import larcv; print(larcv.get_lib_dir())")
LARCV_LIBDIR=${LARCV3_BASEDIR}/lib

pybind_incdir=${LARCV3_BASEDIR}/src/pybind11_json/include
json_incdir=${LARCV3_BASEDIR}/src/json/include
pybind2_incdir=${LARCV3_BASEDIR}/src/pybind11/include


# export PYTHONPATH=$PYTHONPATH_BACKUP

# export ROOT_FLAGS=$(root-config --cflags)
# export ROOT_LIBS=$(root-config --libs)


CFLAGS=-I. -I${LARCV3_INCDIR} -I${LARCV_INCDIR} ${H5_INCDIR} -I${pybind_incdir} -I${pybind2_incdir} -I${json_incdir} $(shell python3-config --includes) -I $(shell root-config --cflags) -g -fPIC


LDFLAGS=$(shell root-config --libs) ${ROOT_LIBS} -L ${LARCV_LIBDIR} \
-llarcv -L ${LARCV3_LIBDIR} -llarcv3 \
${H5_LIBDIR} -lhdf5 -lhdf5_cpp

OBJ23 = larcv2_to_larcv3.o
EXEC23 = larcv2_to_larcv3.so

OBJ32 = larcv3_to_larcv2.o
EXEC32 = larcv3_to_larcv2.so

all: $(EXEC23) $(EXEC32)


%.o: %.cpp %.h
	$(CC) -c -o $@ $< $(CFLAGS)

$(EXEC23): $(OBJ23)
	$(CC) -shared -o $@ $^ $(LDFLAGS)

$(EXEC32): $(OBJ32)
	$(CC) -shared -o $@ $^ $(LDFLAGS)

clean:
	rm *.o
	rm $(EXEC23) $(EXEC32) 