#====================================================================
#  Makefile for "md_fdps.x"
#    2017/8/23  K.Kawai Tokumasu-Lab.
#====================================================================

#--- PATH for FDPS src directry
#------ default
PS_DIR = $(FDPS_ROOT)
#------ absolute path
#PS_DIR = $(HOME)/FDPS/FDPS_3.0_20161228
#PS_DIR = $(HOME)/FDPS/FDPS_3.0_gcc54

#--- PATH for FFTw library
#------ default
FFTW_DIR = $(FFTW_ROOT)
#------ absolute path
#FFTW_DIR = $(HOME)/local/fftw-3.3.6
#FFTW_DIR = $(HOME)/local/fftw-3.3.6-ompi-2.1.1-gcc-5.4

PS_PATH  = -I$(PS_DIR)/src/
PS_PATH += -I$(PS_DIR)/src/particle_mesh/
#LIB_PM   =   $(PS_DIR)/src/particle_mesh/libpm.a
LIB_PM   =   $(PS_DIR)/src/particle_mesh/libpm_debug.a

INCLUDE_FFTW  = -I$(FFTW_DIR)/include/
LIB_FFTW      = -L$(FFTW_DIR)/lib/ -lfftw3f_mpi -lfftw3f
LIB_FFTW_STATIC  = $(FFTW_DIR)/lib/libfftw3f_mpi.a
LIB_FFTW_STATIC += $(FFTW_DIR)/lib/libfftw3f.a

#--- add MD extention
INCLUDE_MD_EXT = -I./md_ext_include/

#--- compiler settings
CC = mpicxx
#CFLAGS = -I./src/ -O3 -ffast-math -funroll-loops -lm -std=c++11 -g3
CFLAGS = -I./src/ -O0 -Wall -lm -std=c++11 -g3 -p

#--- enable parallelization in FDPS
#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

#--- C++ compile target
PROGRAM = md_fdps.x
SRC_CPP = ./src/md_main.cpp
OBJ_CPP = $(SRC_CPP:%.cpp=%.o)
OBJS    = $(OBJ_CPP)

#--- dependency check
DEPS = $(SRC_CPP:%.cpp=%.d)


#=======================
#  main MD program
#=======================
md_main: $(SRC_CPP)
	$(CC) $(PS_PATH) $(CFLAGS) $(INCLUDE_FFTW) $(INCLUDE_MD_EXT) -MMD $< $(LIB_PM) $(LIB_FFTW) -o $(PROGRAM)

#=======================
#  unit tests
#=======================
test_vec: ./unit_test/test_vec.cpp
	$(CC) $(PS_PATH) $(CFLAGS) $^ $(INCLUDE_MD_EXT) -o test_vec.x
	./test_vec.x

test_blz_dist: ./unit_test/test_boltzmann_dist.cpp
	$(CC) $(PS_PATH) $(CFLAGS) $^ $(INCLUDE_MD_EXT) -o test_blz_dist.x
	./test_blz_dist.x

test_comm: ./unit_test/test_comm_tool.cpp
	$(CC) $(PS_PATH) $(CFLAGS) $^ $(INCLUDE_MD_EXT) -DDEBUG_COMM_TOOL -o test_comm.x
	mpirun -n 2 ./test_comm.x

test_intra_pair: ./unit_test/test_intra_pair_manager.cpp
	$(CC) $(PS_PATH) $(CFLAGS) $^ $(INCLUDE_FFTW) $(INCLUDE_MD_EXT)  -o test_intra_pair.x
	./test_intra_pair.x

test_model: ./unit_test/test_loading_model.cpp
	$(CC) $(PS_PATH) $(CFLAGS) $^ $(INCLUDE_FFTW) $(INCLUDE_MD_EXT) -DTEST_MOL_INSTALL -o test_model.x
	./test_model.x  AA_wat_SPC_Fw

test_condition: ./unit_test/test_loading_condition.cpp
	$(CC) $(PS_PATH) $(CFLAGS) $^ $(INCLUDE_FFTW) $(INCLUDE_MD_EXT) -o test_condition.x
	mpirun -n 2 ./test_condition.x

test_init: ./unit_test/test_initialize.cpp
	$(CC) $(PS_PATH) $(CFLAGS) $^ $(INCLUDE_FFTW) $(INCLUDE_MD_EXT) -DTEST_MOL_INSTALL -o test_init.x
	mpirun -n 2 ./test_init.x


#--- force & potential test
test_force_LJ: ./unit_test/test_force_LJ.cpp
	$(CC) $(PS_PATH) $(CFLAGS) $^ $(INCLUDE_FFTW) $(INCLUDE_MD_EXT) $(LIB_PM) $(LIB_FFTW_STATIC) -o test_force_LJ.x
	mpirun -n 2 ./test_force_LJ.x

test_force_coulomb: ./unit_test/test_force_coulomb.cpp
	$(CC) $(PS_PATH) $(CFLAGS) $^ $(INCLUDE_FFTW) $(INCLUDE_MD_EXT) $(LIB_PM) $(LIB_FFTW_STATIC) -o test_force_coulomb.x
	mpirun -n 2 ./test_force_coulomb.x

test_force_bond: ./unit_test/test_force_bond.cpp
	$(CC) $(PS_PATH) $(CFLAGS) $^ $(INCLUDE_FFTW) $(INCLUDE_MD_EXT) $(LIB_PM) $(LIB_FFTW_STATIC) -o test_force_bond.x
	mpirun -n 2 ./test_force_bond.x

test_force_angle: ./unit_test/test_force_angle.cpp
	$(CC) $(PS_PATH) $(CFLAGS) $^ $(INCLUDE_FFTW) $(INCLUDE_MD_EXT) $(LIB_PM) $(LIB_FFTW_STATIC) -o test_force_angle.x
	mpirun -n 2 ./test_force_angle.x

test_force_dihedral: ./unit_test/test_force_dihedral.cpp
	$(CC) $(PS_PATH) $(CFLAGS) $^ $(INCLUDE_FFTW) $(INCLUDE_MD_EXT) $(LIB_PM) $(LIB_FFTW_STATIC) -o test_force_dihedral.x
	mpirun -n 2 ./test_force_dihedral.x

test_force_improper: ./unit_test/test_force_improper.cpp
	$(CC) $(PS_PATH) $(CFLAGS) $^ $(INCLUDE_FFTW) $(INCLUDE_MD_EXT) $(LIB_PM) $(LIB_FFTW_STATIC) -o test_force_improper.x
	mpirun -n 2 ./test_force_improper.x


.PHONY: clean
clean:
	rm -rf *.x $(DEPS) $(OBJS)

.PHONY: clean_all
clean_all:
	make clean
	rm -rf ./posdata ./pdb ./resume *.dat *.out ./src/enum_model.hpp

-include $(DEPS)
