# This is a makefile for SALMON program.
# please select archtecture by deleting "#"

#ARCH = gnu
ARCH = intel
#ARCH = intel-avx
#ARCH = intel-avx2
#ARCH = fujitsu
#ARCH = intel-knl
#ARCH = intel-knc

USE_SCALAPACK = yes

#### explanation for environmental values #########
# FC: compiler                                    #
# FFLAGS: compiler option                         # 
# FILE_MATHLIB: math libraries to be used         #
# LIBLAPACK: options for math libraries           #
# LIBSCALAPACK: options for scalapack libraries   #
###################################################

ifeq ($(ARCH), gnu)
    TARGET = salmon.cpu
    FC = mpifc
    CC = mpicc
    FFLAGS = -O3 -fopenmp -Wall -cpp -ffree-form -ffree-line-length-none
    CFLAGS = -O3 -fopenmp -Wall
    FILE_MATHLIB = lapack
    LIBSCALAPACK = -lmkl_intel_thread -lmkl_intel_lp64 -lmkl_core -lpthread -ldl -liomp5 \
        -lmkl_blacs_intelmpi_lp64 -lmkl_scalapack_lp64 -lm
    MODULE_SWITCH = -J
endif

ifeq ($(ARCH), intel)
    TARGET = salmon.cpu
    FC = mpiifort
    CC = mpiicc
    FFLAGS = -O3 -qopenmp -ansi-alias -fno-alias -fpp -nogen-interface -std03 -warn all
    CFLAGS = -O3 -qopenmp -ansi-alias -fno-alias -Wall -restrict
    FILE_MATHLIB = lapack
    LIBSCALAPACK = -mkl=cluster
    MODULE_SWITCH = -module
endif

ifeq ($(ARCH), intel-avx)
    TARGET = salmon.cpu
    FC = mpiifort
    CC = mpiicc
    FLAGS = -xAVX -qopenmp -ansi-alias -fno-alias \
            -DARTED_STENCIL_OPTIMIZED \
            -DARTED_STENCIL_WITH_C \
            -DARTED_EXPLICIT_VECTORIZATION \
            -DARTED_REDUCE_FOR_MANYCORE 
    FFLAGS = $(FLAGS) -O3 -fpp -nogen-interface -std90 -warn all -diag-disable 6187,6477,6916,7025,7416,7893
    CFLAGS = $(FLAGS) -O3 -Wall -diag-disable=10388 -restrict
    FILE_MATHLIB = lapack
    LIBSCALAPACK = -mkl=cluster
    SIMD_SET = AVX
    MODULE_SWITCH = -module
endif

ifeq ($(ARCH), intel-avx2)
    TARGET = salmon.cpu
    FC = mpiifort
    CC = mpiicc
    FLAGS = -xCORE-AVX2 -qopenmp -ansi-alias -fno-alias \
            -DARTED_REDUCE_FOR_MANYCORE 
    FFLAGS = $(FLAGS) -O3 -fpp -nogen-interface -std90 -warn all -diag-disable 6187,6477,6916,7025,7416,7893
    CFLAGS = $(FLAGS) -O3 -Wall -restrict
    FILE_MATHLIB = lapack
    LIBSCALAPACK = -mkl=cluster
    SIMD_SET = AVX
    MODULE_SWITCH = -module
endif

ifeq ($(ARCH), fujitsu)
    TARGET = salmon.cpu
    FC = mpifrtpx
    CC = mpifccpx
    FFLAGS = -O3 -Kfast,openmp,simd=1 -Cpp -Kocl,nooptmsg
    CFLAGS = -O3 -Kfast,openmp,simd=1 -Kocl,nooptmsg
    FILE_MATHLIB = lapack
    LIBSCALAPACK = -SCALAPACK -SSL2BLAMP
    MODULE_SWITCH = -M
endif

ifeq ($(ARCH), intel-knl)
    TARGET = salmon.mic
    FC = mpiifort
    CC = mpiicc
    FLAGS = -xMIC-AVX512 -qopenmp -qopt-ra-region-strategy=block -ansi-alias -fno-alias \
            -DARTED_STENCIL_OPTIMIZED \
            -DARTED_STENCIL_WITH_C \
            -DARTED_EXPLICIT_VECTORIZATION \
            -DARTED_REDUCE_FOR_MANYCORE \
            -DARTED_ENABLE_SOFTWARE_PREFETCH
    FFLAGS = $(FLAGS) -O3 -fpp -nogen-interface -std03 -warn all -diag-disable 6187,6477,6916,7025,7416
    CFLAGS = $(FLAGS) -O3 -Wall -diag-disable=10388 -restrict
    FILE_MATHLIB = lapack
    LIBSCALAPACK = -mkl=cluster
    SIMD_SET = IMCI
    MODULE_SWITCH = -module
endif

ifeq ($(ARCH), intel-knc)
    TARGET = salmon.mic
    FC = mpiifort
    CC = mpiicc
    FLAGS = -mmic -qopenmp -qopt-assume-safe-padding -qopt-streaming-stores always -qopt-gather-scatter-unroll=4 \
        -qopt-ra-region-strategy=block -ansi-alias -fno-alias \
        -DARTED_STENCIL_OPTIMIZED \
        -DARTED_STENCIL_WITH_C \
        -DARTED_EXPLICIT_VECTORIZATION \
        -DARTED_REDUCE_FOR_MANYCORE \
        -DARTED_ENABLE_SOFTWARE_PREFETCH
    FFLAGS = $(FLAGS) -O3 -fpp -nogen-interface -std90 -warn all -diag-disable 6187,6477,6916,7025,7416,7893
    CFLAGS = $(FLAGS) -O3 -Wall -restrict
    FILE_MATHLIB = lapack
    LIBSCALAPACK = -mkl=cluster
    SIMD_SET = IMCI
    MODULE_SWITCH = -module
endif

####################################################
##### please do not modify following sentences #####
####################################################

SRC_main = main.f90 end_parallel.f90 read_input.f90 setup_parallel.f90

SRC_core = exc_cor.f90

SRC_GCEED = gceed.f90 read_input_gceed.f90 

SRC_SCF_GCEED = real_space_dft.f90 subspace_diag.f90 simple_mixing.f90 copy_density.f90 \
                prep_ini.f90 gram_schmidt.f90 calc_rho_in.f90 subdgemm_$(FILE_MATHLIB).f90 \
                calc_occupation.f90 read_input_scf.f90 sendrecv_copy.f90 \
                deallocate_sendrecv_groupob.f90 structure_opt.f90 ybcg.f90 \
                rmmdiis_eigen.f90 rmmdiis.f90 calc_dos.f90 calc_pdos.f90 

ifeq ($(USE_SCALAPACK), yes)
   SRC_SCF_LAPACK_GCEED = eigen_subdiag_scalapack.f90
else
   SRC_SCF_LAPACK_GCEED = eigen_subdiag_$(FILE_MATHLIB).f90
endif

SRC_RT_GCEED = real_time_dft.f90 taylor.f90 taylor_coe.f90 WriteDensity.f90 read_rt.f90 \
               time_evolution_step.f90 hpsi_groupob.f90 gradient_ex.f90 \
               calcEstatic.f90 total_energy_groupob.f90 xc_fast.f90 \
               dip.f90 read_input_rt.f90 \
               projection.f90 add_polynomial.f90 calcVbox.f90 

SRC_COMMON_GCEED = calcELF.f90 psl.f90 hartree.f90 ylm.f90 xc.f90 OUT_IN_data.f90 \
                   laplacianh.f90 nabla.f90 copyVlocal.f90 \
                   inner_product3.f90 inner_product4.f90 \
                   setlg.f90 setmg.f90 setng.f90 init_wf.f90 calc_force.f90 \
                   calc_gradient_fast.f90 \
                   calc_gradient_fast_c.f90 check_numcpu.f90 writepsi.f90 \
                   calcJxyz.f90 storevpp.f90 calcuV.f90 calcJxyz2nd.f90 \
                   calcVpsl.f90 bisection.f90 calc_Mps3rd.f90 \
                   setbN.f90 setcN.f90 hartree_cg.f90 hartree_boundary.f90 \
                   calc_myob.f90 check_corrkob.f90 calc_allob.f90 calc_ob_num.f90 \
                   calc_pmax.f90 calc_iroot.f90 calc_iquotient.f90 set_isstaend.f90 \
                   conv_p0.f90 conv_p.f90 set_ispin.f90 read_copy_pot.f90 \
                   set_gridcoo.f90 calc_force_c.f90 conv_core_exc_cor.f90

MOD_GCEED = scf_data.f90 allocate_mat.f90 new_world.f90 \
            init_sendrecv.f90 sendrecv.f90 laplacian2.f90 \
            gradient2.f90 hpsi2.f90 \
            deallocate_mat.f90 copy_psi_mesh.f90 inner_product.f90 \
            rmmdiis_eigen_$(FILE_MATHLIB).f90 share_mesh_1d_old.f90 \
            read_pslfile.f90 total_energy.f90 calc_density.f90 \
            change_order.f90 allocate_sendrecv_groupob.f90 allocate_psl.f90 \
            gradient.f90 sendrecvh.f90 calc_invA_$(FILE_MATHLIB).f90 \
            writebox_rt.f90 readbox_rt.f90 sendrecv_groupob.f90 \
            sendrecv_groupob_ngp.f90 

SRC_ARTED0 = main.f90 common/Exc_Cor.f90 common/Hartree.f90 common/hpsi.f90 \
            common/ion_force.f90 common/preprocessor.f90 common/psi_rho.f90 \
            common/reentrance.f90 common/total_energy.f90 common/Ylm_dYlm.f90 \
            FDTD/beam.f90 FDTD/FDTD.f90 GS/CG.f90 GS/Density_Update.f90 \
            GS/diag.f90 GS/Fermi_Dirac_distribution.f90 GS/Gram_Schmidt.f90 \
            GS/Occupation_Redistribution.f90 GS/sp_energy.f90 \
            GS/write_GS_data.f90 preparation/fd_coef.f90 preparation/init.f90 \
            preparation/init_wf.f90 preparation/input_ps.f90 \
            preparation/prep_ps.f90 RT/current.f90 RT/dt_evolve.f90 \
            RT/Fourier_tr.f90 RT/hamiltonian.f90 RT/init_Ac.f90 \
            RT/k_shift_wf.f90 

C_SRC_ARTED0 = modules/env_variables_internal.c

ifdef SIMD_SET
    SRC_ARTED = $(SRC_ARTED0) 
    C_SRC_ARTED = $(C_SRC_ARTED0) stencil/C/$(SIMD_SET)/current.c stencil/C/$(SIMD_SET)/hpsi.c stencil/C/$(SIMD_SET)/total_energy.c 
else
    SRC_ARTED = $(SRC_ARTED0) stencil/F90/current.f90 stencil/F90/hpsi.f90 stencil/F90/total_energy.f90 
    C_SRC_ARTED = $(C_SRC_ARTED0)
endif

MOD_ARTED = modules/backup_routines.f90 modules/communication.f90 \
            modules/env_variables.f90 modules/global_variables.f90 \
            modules/misc_routines.f90 modules/nvtx.f90 modules/timer.f90 \
            control/inputfile.f90 modules/opt_variables.f90 \
            modules/performance_analyzer.f90 control/control_ms.f90 \
            control/control_sc.f90 





OBJDIR = obj

OBJ_GCEED = $(addprefix $(OBJDIR)/GCEED/,$(SRC_GCEED:.f90=.o))
OBJ_SCF_GCEED = $(addprefix $(OBJDIR)/GCEED/scf/,$(SRC_SCF_GCEED:.f90=.o))
OBJ_SCF_LAPACK_GCEED = $(addprefix $(OBJDIR)/GCEED/scf/,$(SRC_SCF_LAPACK_GCEED:.f90=.o))
OBJ_RT_GCEED = $(addprefix $(OBJDIR)/GCEED/rt/,$(SRC_RT_GCEED:.f90=.o))
OBJ_COMMON_GCEED = $(addprefix $(OBJDIR)/GCEED/common/,$(SRC_COMMON_GCEED:.f90=.o))
OBJM_GCEED = $(addprefix $(OBJDIR)/GCEED/modules/,$(MOD_GCEED:.f90=.o))
OBJS_GCEED = $(OBJM_GCEED) $(OBJ_SCF_LAPACK_GCEED) $(OBJ_SCF_GCEED) $(OBJ_RT_GCEED) $(OBJ_COMMON_GCEED) $(OBJ_GCEED)  

OBJ_ARTED = $(addprefix $(OBJDIR)/ARTED/,$(SRC_ARTED:.f90=.o))
C_OBJ_ARTED = $(addprefix $(OBJDIR)/ARTED/,$(C_SRC_ARTED:.c=.o))
OBJM_ARTED= $(addprefix $(OBJDIR)/ARTED/,$(MOD_ARTED:.f90=.o))
OBJS_ARTED= $(OBJM_ARTED) $(OBJ_ARTED) $(C_OBJ_ARTED)

OBJ_main = $(addprefix $(OBJDIR)/main/,$(SRC_main:.f90=.o))
OBJ_core = $(addprefix $(OBJDIR)/src/core/,$(SRC_core:.f90=.o))

OBJS = $(OBJS_GCEED) $(OBJS_ARTED) $(OBJ_core) $(OBJ_main)

.SUFFIXES:
.SUFFIXES: .F .F90 .o

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ -I $(OBJDIR) $(OBJS) $(LIBSCALAPACK)

$(OBJDIR)/%.o : %.f90
	@if [ ! -d $(dir $@) ]; then mkdir -p $(dir $@); fi
	$(FC) $(FFLAGS) $(MODULE_SWITCH) $(OBJDIR) -o $@ -c $<

$(OBJDIR)/%.o : %.c
	@if [ ! -d $(dir $@) ]; then mkdir -p $(dir $@); fi
	$(CC) $(CFLAGS) -o $@ -c $<

$(OBJS_GCEED): $(addprefix GCEED/modules/,$(MOD_GCEED))
$(OBJS_ARTED): $(addprefix ARTED/,$(MOD_ARTED))

clean: 
	rm -f $(TARGET)
	rm -rf $(OBJDIR)/*
