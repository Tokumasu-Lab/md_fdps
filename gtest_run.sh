#!/bin/bash -e

#--- setting for debug of the test codes.
export GTEST_CATCH_EXCEPTIONS=0

usage_exit() {
	echo "Usage: $0 [-j n_core] [-m n_mpi] [-o n_omp]"
	echo "    n_core: # of core for make. used as 'make -j (n_core)' internally."
	echo "    n_mpi : # of process for MPI."
	echo "    n_omp : # of threads for OpenMP."
	exit 1
}

#--- default settings
EXE_DIR="./test_bin"
CXX_NUM=1
MPI_NUM=2
OMP_NUM=1

#--- get option
while getopts j:m:o:h OPT
do
	case $OPT in
		j) CXX_NUM=$OPTARG
			;;
		m) MPI_NUM=$OPTARG
			;;
		o) OMP_NUM=$OPTARG
			;;
		h) usage_exit
	esac
done

#--- compile unit test code
MAKE_ARG=""
if [ "$CXX_NUM" -gt 2 ]; then
	MAKE_ARG=" -j $CXX_NUM"
fi
make gtest $MAKE_ARG

#--- PS::Vector3<T> extension
${EXE_DIR}/gtest_vec_ext

#--- COMM_TOOL::
${EXE_DIR}/gtest_comm_tool_SerDes
mpirun -np ${MPI_NUM} ${EXE_DIR}/gtest_comm_tool_broadcast
mpirun -np ${MPI_NUM} ${EXE_DIR}/gtest_comm_tool_gather
mpirun -np ${MPI_NUM} ${EXE_DIR}/gtest_comm_tool_scatter
mpirun -np ${MPI_NUM} ${EXE_DIR}/gtest_comm_tool_allGather
mpirun -np ${MPI_NUM} ${EXE_DIR}/gtest_comm_tool_allToAll

#--- MD_EXT::boltzmann_dist
${EXE_DIR}/gtest_blz_dist

#--- MD_EXT::CellIndex
${EXE_DIR}/gtest_cell_index

#--- static array for FullParticle
${EXE_DIR}/gtest_fixed_vector
${EXE_DIR}/gtest_basic_connect

#--- logspace array
${EXE_DIR}/gtest_logspace_array

#--- IntraPair::
mpirun -np ${MPI_NUM} -x OMP_NUM_THREADS=${OMP_NUM} ${EXE_DIR}/gtest_intra_pair

#--- force test (fixed MPI process because there are few particle only)
mpirun -np 2 -x OMP_NUM_THREADS=${OMP_NUM} ${EXE_DIR}/gtest_force_LJ
mpirun -np 2 -x OMP_NUM_THREADS=${OMP_NUM} ${EXE_DIR}/gtest_force_coulomb
mpirun -np 2 -x OMP_NUM_THREADS=${OMP_NUM} ${EXE_DIR}/gtest_force_bond
mpirun -np 2 -x OMP_NUM_THREADS=${OMP_NUM} ${EXE_DIR}/gtest_force_angle
mpirun -np 2 -x OMP_NUM_THREADS=${OMP_NUM} ${EXE_DIR}/gtest_force_dihedral
mpirun -np 2 -x OMP_NUM_THREADS=${OMP_NUM} ${EXE_DIR}/gtest_force_improper

mpirun -np 2 -x OMP_NUM_THREADS=${OMP_NUM} ${EXE_DIR}/gtest_force_mask

#--- file I/O test
mpirun -np ${MPI_NUM} -x OMP_NUM_THREADS=${OMP_NUM} ${EXE_DIR}/gtest_fileIO

#--- analysis module test
mpirun -np ${MPI_NUM} -x OMP_NUM_THREADS=${OMP_NUM} ${EXE_DIR}/gtest_analysis_RDF
mpirun -np ${MPI_NUM} -x OMP_NUM_THREADS=${OMP_NUM} ${EXE_DIR}/gtest_analysis_MSD
mpirun -np ${MPI_NUM} -x OMP_NUM_THREADS=${OMP_NUM} ${EXE_DIR}/gtest_analysis_clustering
