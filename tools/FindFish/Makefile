#
# Makefile for interpolation code
# CVS:$Id: makefile,v 1.3 2000/04/19 22:05:59 pwjones Exp $
#

##COMP = gfortran
COMP = mpif90
COMP = ifort
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
        TS_boundaries_data.o \
        create_rest_mask.o \
        SSH_boundaries_data.o \
        create_slvl_mask.o \
        cnvts3d.o \
        fix_kmt_ff.o \
        read_ssh.o

all: rst2init_ts \
     read_init_ts \
     read_temp   \
     create_init_ts_WP \
     TS_boundaries_data \
     create_rest_mask  \
     SSH_boundaries_data \
     create_slvl_mask \
     cnvts3d \
     fix_kmt_ff \
     read_ssh

read_ssh: $(OBJ) read_ssh.o
	$(COMP) $(FLAGS) read_ssh.o $(LIB) -o $(EXEDIR)/read_ssh

read_ssh.o: read_ssh.f90
	$(COMP) $(FLAGS) $(INCL) -c read_ssh.f90

cnvts3d: $(OBJ) cnvts3d.o
	$(COMP) $(FLAGS) cnvts3d.o $(LIB) -o $(EXEDIR)/cnvts3d

cnvts3d.o: cnvts3d.f90
	$(COMP) $(FLAGS) $(INCL) -c cnvts3d.f90

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

SSH_boundaries_data: $(OBJ) SSH_boundaries_data.o
	$(COMP) $(FLAGS) SSH_boundaries_data.o $(LIB) -o $(EXEDIR)/SSH_boundaries_data

SSH_boundaries_data.o: SSH_boundaries_data.f90
	$(COMP) $(FLAGS) $(INCL) -c SSH_boundaries_data.f90

create_rest_mask: $(OBJ) create_rest_mask.o
	$(COMP) $(FLAGS) create_rest_mask.o $(LIB) -o $(EXEDIR)/create_rest_mask

create_rest_mask.o: create_rest_mask.f90
	$(COMP) $(FLAGS) $(INCL) -c create_rest_mask.f90

create_slvl_mask: $(OBJ) create_slvl_mask.o
	$(COMP) $(FLAGS) create_slvl_mask.o $(LIB) -o $(EXEDIR)/create_slvl_mask

create_slvl_mask.o: create_slvl_mask.f90
	$(COMP) $(FLAGS) $(INCL) -c create_slvl_mask.f90

fix_kmt_ff: $(OBJ) fix_kmt_ff.o
	$(COMP) $(FLAGS) fix_kmt_ff.o $(LIB) -o $(EXEDIR)/fix_kmt_ff

fix_kmt_ff.o: fix_kmt_ff.f90
	$(COMP) $(FLAGS) $(INCL) -c fix_kmt_ff.f90

clean: 
	/bin/rm *.o rst2init_ts read_init_ts create_init_ts_WP SSH_boundaries_data TS_boundaries_data create_rest_mask fix_kmt_ff cnvts3d read_temp read_ssh
