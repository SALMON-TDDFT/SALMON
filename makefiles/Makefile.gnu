# gnu

SALMON = $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/../)

TARGET = $(SALMON)/bin/salmon.cpu
FC = mpif90
CC = mpicc
FFLAGS = -O3 -fopenmp -Wall -cpp -ffree-form -ffree-line-length-none
CFLAGS = -O3 -fopenmp -Wall
LIBLAPACK = -llapack -lblas
#LIBLAPACK = -lmkl_intel_thread -lmkl_intel_lp64 -lmkl_core -lpthread -ldl -liomp5 -lm
MODULE_SWITCH = -J
COMM_SET =

include $(SALMON)/makefiles/make.body
