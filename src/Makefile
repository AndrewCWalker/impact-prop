CC = gcc
MPICC = mpicc
INCL = 
LIBDIRS = -L.
LIBS = -lgsl -lhdf5 -lgslcblas

impactor = ../../impact/o_r/drag
pysysdevel = ../../pysysdevel

ifneq (,$(CSPICE_DIR))
  ifeq (,$(findstring "$(CSPICE_DIR)",$(wildcard "$(CSPICE_DIR)")))
    local_cspice = true
  endif
endif
ifeq ($(local_cspice),true)
  INCL += -I$(CSPICE_DIR)/include
  LIBDIRS += -L$(CSPICE_DIR)/lib
  LIBS += $(CSPICE_DIR)/lib/cspice.a
else
  INCL += -I/n/packages/include/naif
  LIBDIRS += -L/n/packages/lib
  LIBS += /n/packages/lib/cspice.a
endif
CFLAGS = -Wall -g -pg -fPIC -DDATA_DIR="\"../data/\"" $(INCL) $(LIBDIRS)
MFLAGS = $(CFLAGS) -DINLINE -DPARALLEL=1


HDRS := defs.h structs.h test.h
HDRS += acceleration.h assoclegendre.h atmprop.h Cd.h design.h emu.h
HDRS += interpolation.h io.h misc.h msisinputs.h nrlmsise-00.h parameters.h
HDRS += propagator.h rk.h stddata.h transformation.h

SRCS := acceleration.c assoclegendre.c atmprop.c Cd.c emu.c config.c
SRCS += interpolation.c io.c misc.c msisinputs.c nrlmsise-00.c
SRCS += nrlmsise-00_data.c propagator.c rhseom.c rk.c transformation.c

OBJS := $(SRCS:.c=.o)

EXE_NAME := propagator
LIB_NAME := libpropagator.so


$(EXE_NAME): main.c $(LIB_NAME)
	$(CC) -o $@ $^ $(LIBS) $(CFLAGS)

$(LIB_NAME): $(OBJS)
	$(CC) -o $@ $^ $(LIBS) $(CFLAGS) -shared


$(EXE_NAME)_mpi: CC=$(MPICC)
$(EXE_NAME)_mpi: main.c $(OBJS)
	$(CC) -o $@ $^ $(LIBS) $(MFLAGS)

%.o: %.c $(HDRS)
	$(CC) -c $< $(CFLAGS)


structs.h: $(impactor)/structs.py
	@PYTHONPATH=$(pysysdevel) python -m sysdevel.distutils.headers $(impactor)/structs.py structs.h


proptest: $(OBJS) test.o
	@$(CC) -o $@ $^ $(LIBS) $(CFLAGS)

check: proptest
	./proptest


clean :
	@rm -f *.o *~

pristine : clean
	@rm -f $(EXE_NAME) $(LIB_NAME) structs.h
