#====================================================================
#  Makefile for "md_fdps.x"
#    2017/8/23  K.Kawai Tokumasu-Lab.
#====================================================================

#--- PATH for FDPS src directry
#------ for paolo in Tokumasu Lab.
#PS_DIR = /WORK4/kawai/FDPS/FDPS_3.0_20161228
#------ for IFS
#PS_DIR = /home/NIa/NIa030/FDPS_3.0_20161228
#------ for virtual box Ubuntu
PS_DIR = $(HOME)/FDPS/FDPS_3.0_20161228

PS_PATH  = -I$(PS_DIR)/src/
PS_PATH += -I$(PS_DIR)/src/particle_mesh/
#LIB_PM   =   $(PS_DIR)/src/particle_mesh/libpm.a
LIB_PM   =   $(PS_DIR)/src/particle_mesh/libpm_debug.a

#--- PATH for FFTw library
FFTW_DIR      = $(HOME)/local/fftw-3.3.6

INCLUDE_FFTW  = -I$(FFTW_DIR)/include/
LIB_FFTW      = -L$(FFTW_DIR)/lib/ -lfftw3f_mpi -lfftw3f
LIB_FFTW_STATIC  = $(FFTW_DIR)/lib/libfftw3f_mpi.a
LIB_FFTW_STATIC += $(FFTW_DIR)/lib/libfftw3f.a

#--- add MD extention
INCLUDE_MD_EXT = -I./md_ext_include/

#--- dinamic header
ENUM_MODEL = ./src/md_fdps_enum_model.hpp

#--- compiler settings
#CC = g++
CC = mpicxx
#CFLAGS = -I./src/ -O3 -ffast-math -funroll-loops -lm -std=c++11 -g3
CFLAGS = -I./src/ -O0 -Wall -lm -std=c++11 -g3 -p

#--- enable parallelization in FDPS
#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

#--- C++ compile target
PROGRAM = md_fdps.x
SRC_CPP = ./src/md_fdps_main.cpp
OBJ_CPP = $(SRC_CPP:%.cpp=%.o)
OBJS    = $(OBJ_CPP)

#--- dependency check
DEPS = $(SRC_CPP:%.cpp=%.d)


#=======================
#  main MD program
#=======================
$(PROGRAM): $(SRC_CPP)
	$(CC) $(PS_PATH) $(CFLAGS) $(INCLUDE_FFTW) $(INCLUDE_MD_EXT) -MMD $(LIB_FFTW) $< $(LIB_PM) -o $@

#=======================
#  unit tests
#=======================
test_vec: ./unit_test/vec.cpp
	$(CC) $(PS_PATH) $(CFLAGS) $^ $(INCLUDE_MD_EXT) -o test_vec.x
	./test_vec.x

test_blz_dist: ./unit_test/boltzmann_dist.cpp
	$(CC) $(PS_PATH) $(CFLAGS) $^ $(INCLUDE_MD_EXT) -o test_blz_dist.x
	./test_blz_dist.x

test_comm: ./unit_test/comm_tool.cpp
	$(CC) $(PS_PATH) $(CFLAGS) $^ $(INCLUDE_MD_EXT) -DDEBUG_COMM_TOOL -o test_comm.x
	mpirun -n 2 ./test_comm.x

test_intra_pair: ./unit_test/intra_pair_manager.cpp
	$(CC) $(PS_PATH) $(CFLAGS) $^ $(INCLUDE_MD_EXT) -o test_intra_pair.x
	./test_intra_pair.x

test_model: ./unit_test/loading_model.cpp
	$(CC) $(PS_PATH) $(CFLAGS) $^ $(INCLUDE_FFTW) $(INCLUDE_MD_EXT) -DTEST_MOL_INSTALL -o test_model.x
	./test_model.x  AA_wat

test_condition: ./unit_test/loading_condition.cpp
	$(CC) $(PS_PATH) $(CFLAGS) $^ $(INCLUDE_FFTW) $(INCLUDE_MD_EXT) -o test_condition.x
	mpirun -n 2 ./test_condition.x

test_init: ./unit_test/initialize.cpp
	$(CC) $(PS_PATH) $(CFLAGS) $^ $(INCLUDE_FFTW) $(INCLUDE_MD_EXT) -DTEST_MOL_INSTALL -o test_init.x
	mpirun -n 2 ./test_init.x


#--- force & potential test
test_force_LJ: ./unit_test/force_LJ.cpp
	$(CC) $(PS_PATH) $(CFLAGS) $^ $(INCLUDE_FFTW) $(INCLUDE_MD_EXT) $(LIB_PM) $(LIB_FFTW_STATIC) -o test_force_LJ.x
	mpirun -n 2 ./test_force_LJ.x

test_force_coulomb: ./unit_test/force_coulomb.cpp
	$(CC) $(PS_PATH) $(CFLAGS) $^ $(INCLUDE_FFTW) $(INCLUDE_MD_EXT) $(LIB_PM) $(LIB_FFTW_STATIC) -o test_force_coulomb.x
	mpirun -n 2 ./test_force_coulomb.x

test_force_bond: ./unit_test/force_bond.cpp
	$(CC) $(PS_PATH) $(CFLAGS) $^ $(INCLUDE_FFTW) $(INCLUDE_MD_EXT) $(LIB_PM) $(LIB_FFTW_STATIC) -o test_force_bond.x
	mpirun -n 2 ./test_force_bond.x

test_force_angle: ./unit_test/force_angle.cpp
	$(CC) $(PS_PATH) $(CFLAGS) $^ $(INCLUDE_FFTW) $(INCLUDE_MD_EXT) $(LIB_PM) $(LIB_FFTW_STATIC) -o test_force_angle.x
	mpirun -n 2 ./test_force_angle.x


.PHONY: clean
clean:
	rm -rf *.x $(DEPS) $(OBJS)

-include $(DEPS)
