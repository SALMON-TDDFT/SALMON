# gnu

TARGET = salmon.cpu
FC = mpif90
CC = mpicc
FFLAGS = -O3 -fopenmp -Wall -cpp -ffree-form -ffree-line-length-none
CFLAGS = -O3 -fopenmp -Wall
LIBLAPACK = -llapack -lblas
#LIBLAPACK = -lmkl_intel_thread -lmkl_intel_lp64 -lmkl_core -lpthread -ldl -liomp5 -lm
MODULE_SWITCH = -J
COMM_SET =

LIBXC_LIB =
LIBXC_INC =
# LIBXC_LIB = -L<I<libxc_install_dir>/lib -lxcf90
# LIBXC_INC = -DSALMON_USE_LIBXC -I<libxc_install_dir>/include/ 

ifneq (,$(wildcard make.body))
include make.body
else 
include gnu_makefiles/make.body
endif
