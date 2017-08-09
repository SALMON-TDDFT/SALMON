# gnu

TARGET = salmon.cpu
FC = mpif90
CC = mpicc
FFLAGS = -O3 -fopenmp -Wall -cpp -ffree-form -ffree-line-length-none
CFLAGS = -O3 -fopenmp -Wall
LIBLAPACK = -llapack -lblas
MODULE_SWITCH = -J
COMM_SET =

include $(dir $(abspath $(lastword $(MAKEFILE_LIST))))/make.body
