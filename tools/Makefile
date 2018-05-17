#
# Makefile for interpolation code
# CVS:$Id: makefile,v 1.3 2000/04/19 22:05:59 pwjones Exp $
#

##COMP = gfortran
COMP = mpif90
FLAGS = -O3 -assume byterecl -convert big_endian -assume ieee_fpe_flags
#FLAGS = -O3 -assume byterecl -convert big_endian
NETCDF_DIR = /usr/local/netcdf-4.0-gfortran-openmpi.1.7.5_gcc48
NETCDF_DIR = /home/plgrid/plgjjakacki/local/netcdf-4.0
#NETCDF_DIR = /home/plgrid/plgjaniopan/local/netcdf-4.4.2_hdf5-1.8.13_mpich
#HDF5_DIR = /home/plgrid/plgjaniopan/local/hdf5-1.8.13_mpich
NETCDF_INC = -I$(NETCDF_DIR)/include
NETCDF_LIB = -L$(NETCDF_DIR)/lib
HDF5_LIB = -L$(HDF5_DIR)/lib
HDF5_INC = -I$(HDF5_DIR)/include
XIOS_HOME = /home/plgrid/plgjaniopan/local/xios-1.0
XIOS_LIB = -L$(XIOS_HOME)/lib -lxios
XIOS_INC = -I$(XIOS_HOME)/inc
#LIB  = -lcurl $(NETCDF_LIB) -lnetcdf -lnetcdff $(HDF5_LIB) -lhdf5 -lhdf5hl_fortran $(XIOS_LIB) -L./ -lstdc++
#INCL = $(NETCDF_INC) $(HDF5_INC) $(XIOS_INC)
LIB  = -lcurl $(NETCDF_LIB) -lnetcdf -lnetcdff
INCL = $(NETCDF_INC)
SRCDIR  =
EXEDIR  = .
OBJ  = \
	rst2init_ts.o \
        read_init_ts.o \
        read_temp.o \
        create_init_ts_WP.o \
        TS_boundaries_data.o

all: rst2init_ts \
     read_init_ts \
     read_temp   \
     create_init_ts_WP \
     TS_boundaries_data

rst2init_ts: $(OBJ) rst2init_ts.o
	$(COMP) $(FLAGS) rst2init_ts.o $(LIB) -o $(EXEDIR)/rst2init_ts

rst2init_ts.o: rst2init_ts.f90
	$(COMP) $(FLAGS) $(INCL) -c rst2init_ts.f90

read_init_ts: $(OBJ) read_init_ts.o
	$(COMP) $(FLAGS) read_init_ts.o $(LIB) -o $(EXEDIR)/read_init_ts

read_init_ts.o: read_init_ts.f90
	$(COMP) $(FLAGS) $(INCL) -c read_init_ts.f90

read_temp: $(OBJ) read_temp.o
	$(COMP) $(FLAGS) read_temp.o $(LIB) -o $(EXEDIR)/read_temp

read_temp.o: read_temp.f90
	$(COMP) $(FLAGS) $(INCL) -c read_temp.f90

create_init_ts_WP: $(OBJ) create_init_ts_WP.o
	$(COMP) $(FLAGS) create_init_ts_WP.o $(LIB) -o $(EXEDIR)/create_init_ts_WP

create_init_ts_WP.o: create_init_ts_WP.f90
	$(COMP) $(FLAGS) $(INCL) -c create_init_ts_WP.f90

TS_boundaries_data: $(OBJ) TS_boundaries_data.o
	$(COMP) $(FLAGS) TS_boundaries_data.o $(LIB) -o $(EXEDIR)/TS_boundaries_data

TS_boundaries_data.o: TS_boundaries_data.f90
	$(COMP) $(FLAGS) $(INCL) -c TS_boundaries_data.f90

clean: 
	/bin/rm *.o rst2init_ts read_init_ts create_init_ts_WP TS_boundaries_data
