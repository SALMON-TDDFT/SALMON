# Input variables of SALMON

- [&calculation](#calculation)
- [&control](#control)
- [&units](#units)
- [&parallel](#parallel)
- [&system](#system)
- [&pseudo](#pseudo)
- [&functional](#functional)
- [&rgrid](#rgrid)
- [&kgrid](#kgrid)
- [&tgrid](#tgrid)
- [&propagation](#propagation)
- [&scf](#scf)
- [&emfield](#emfield)
- [&linear_response](#linear_response)
- [&analysis](#analysis)
- [&hartree](#hartree)
- [&ewald](#ewald)

## &calculation
<dl>
<dt>calc_mode; <code>Character</code></dt>
<dd>Choice of Calculation modes. <code>'GS'</code>,
<code>'RTLR'</code>,<code>'RTPulse'</code>, <code>'GS_RTLR'</code>, 
and <code>'GS_RTPulse'</code>
can be chosen.
</dd>

<dt>use_ehrenfest_md; <code>Character</code></dt>
<dd>
Enable(<code>'y'</code>)/disable(<code>'n'</code>) 
Ehrenfest dynamics.
Default is <code>'n'</code>.
</dd>

<dt>use_ms_maxwell; <code>Character</code></dt>
<dd>
Enable(<code>'y'</code>)/disable(<code>'n'</code>) 
Multi-scale Maxwell-Kohn-Sham coupling. 
Default is <code>'n'</code> 
</dd>

<dt>use_force; <code>Character</code></dt>
<dd>
Enable(<code>'y'</code>)/disable(<code>'n'</code>) 
force calculation.
Default is <code>'n'</code>.
</dd>

<dt>use_geometry_opt; <code>Character</code></dt>
<dd>
Enable(<code>'y'</code>)/disable(<code>'n'</code>) 
geometry optimization.
Default is <code>'n'</code>.
</dd>

</dl>

## &control
<dl>

<dt>restart_option; <code>Character</code></dt>
<dd>Flag for restart. <code>'new'</code> or 
<code>'restart'</code> can be chosen.
<code>'new'</code> is default.
</dd>

<dt>backup_frequency; <code>Integer</code></dt>
<dd>Frequency of backup during the time-propagation. 
If <code>0</code> is set, the backup is not performed.
Default is <code>0</code>.
</dd>


<dt>time_shutdown; <code>Real(8)</code></dt>
<dd>Timer for automatic shutdown. The unit is always second.
If negative time is chosen, the automatic shutdown is not performed.
Default is <code>-1</code> sec.
</dd>

<dt>sysname; <code>Character</code></dt>
<dd>Name of calculation. This is used for a prefix of output files.
Default is <code>default</code>.
</dd>


</dl>


## &units
<dl>

<dt>unit_time; <code>Character</code></dt>
<dd>Unit of time for input variables. 
Atomic unit <code>'a.u.'</code> and 
femtosecond <code>'fs'</code> can be chosen.
Default is <code>'a.u.'</code>.
</dd>

<dt>unit_length; <code>Character</code></dt>
<dd>Unit of length for input variables. 
Atomic unit <code>'a.u.'</code> and 
Aungstrom <code>angstrom</code> can be chosen.
Default is <code>'a.u.'</code>.
</dd>

<dt>unit_energy; <code>Character</code></dt>
<dd>Unit of energy for input variables. 
Atomic unit <code>'a.u.'</code> and 
electron-volt <code>'ev'</code> can be chosen.
Default is <code>'a.u.'</code>.
</dd>

<dt>unit_energy; <code>Character</code></dt>
<dd>Unit of charge for input variables. 
Atomic unit <code>'a.u.'</code> is only available.
</dd>

</dl>

## &parallel
<dl>

<dt>domain_parallel; <code>Character</code></dt>
<dd>
Enable(<code>'y'</code>)/disable(<code>'n'</code>) 
domain parallelization.
Default is <code>'n'</code>.
</dd>

<dt>nproc_ob; <code>Integer</code></dt>
<dd>Number of MPI parallelization for orbitals.
Default is <code>0</code>.
</dd>

<dt>nproc_domain(3); <code>Integer</code></dt>
<dd>Number of MPI parallelization for each direction in 
real-space.
Default is <code>(0,0,0)</code>.
</dd>

<dt>nproc_domain_s(3); <code>Integer</code></dt>
<dd>Number of MPI parallelization for each direction in 
real-space.
Default is <code>(0,0,0)</code>.
</dd>

<dt>num_datafiles_in; <code>Integer</code></dt>
<dd>Number of input files.
Default is <code>0</code>.
</dd>

<dt>num_datafiles_out; <code>Integer</code></dt>
<dd>Number of output files.
Default is <code>0</code>.
</dd>


</dl>

## &system 
<dl>

<dt>iperiodic; <code>Integer</code></dt>
<dd>Dimension for periodic boundary condition.
<code>0</code> is for solated systems, and 
<code>3</code> is for solids.
Default is <code>0</code>.
</dd>

<dt>ispin; <code>Integer</code></dt>
<dd>spin??
</dd>

<dt>al(3); <code>Real(8)</code></dt>
<dd>Lattice constants. Unit of the length can be 
chosen by <code>&units/unit_length</code>.
</dd>

<dt>isym; <code>Integer</code></dt>
<dd>Number of symmetries that can be used for 
reduction of k-points.
Default is <code>0</code>.
</dd>

<dt>crystal_structure; <code>Character</code></dt>
<dd>Name of symmetry that can be used for the resuction of # of k-points.
Default is <code>'none'</code>.
</dd>

<dt>nstate; <code>Integer</code></dt>
<dd>Number of states/bands.
</dd>

<dt>nelec; <code>Integer</code></dt>
<dd>Number of valence electrons.
</dd>

<dt>natom; <code>Integer</code></dt>
<dd>Number of atoms in a calculation cell.
</dd>

<dt>file_atom; <code>Character</code></dt>
<dd>File name of atomic positions.
</dd>

</dl>

## &pseudo
<dl>

<dt>pseudodir; <code>Character</code></dt>
<dd>Directry name for pseudopotential files.
</dd>

<dt>Lmax_ps(:); <code>Integer</code></dt>
<dd>Maximum angular momentum of pseudopotential projectors.
</dd>

<dt>Lloc_ps(:); <code>Integer</code></dt>
<dd>Angular momentum of pseudopotential thant will be treated as local.
</dd>

<dt>iZatom(:); <code>Integer</code></dt>
<dd>Atomic number.
</dd>

<dt>ps_format(:); <code>Character</code></dt>
<dd>Formats for pseudopotentials. 
Yabana-Bertsch format (_rps.dat) <code>KY</code>, 
ABINIT format (.pspnc) <code>ABINIT</code>, 
FHI format (.cpi) <code>FHI</code>, and
ABINITFHI format (.fhi) <code>ABINITFHI</code> can be chosen.
Default is <code>KY</code>.
</dd>

<dt>psmask_option(:); <code>Character</code></dt>
<dd>Enable(<code>'y'</code>)/disable(<code>'n'</code>) 
Fourier filtering for pseudopotentials. 
Default is <code>'n'</code>.
</dd>

<dt>alpha_mask(:); <code>Real(8)</code></dt>
<dd>Parameter for the Fourier filtering for pseudopotential.
Default is <code>'0.8'</code>.
</dd>

<dt>gamma_mask(:); <code>Real(8)</code></dt>
<dd>Parameter for the Fourier filtering for pseudopotential.
Default is <code>'1.8'</code>.
</dd>

<dt>eta_maskk(:); <code>Real(8)</code></dt>
<dd>Parameter for the Fourier filtering for pseudopotential.
Default is <code>'15.0'</code>.
</dd>

</dl>

## &functional
<dl>

<dt>eta_maskk(:); <code>Real(8)</code></dt>
<dd>Parameter for the Fourier filtering for pseudopotential.
Default is <code>'15.0'</code>.
</dd>
Exchange-correlation functionals.
At the moment, the following functionals are avelable.
<ul>
<li>
<code>'PZ'</code>: Perdew-Zunger LDA :Phys. Rev. B 23, 5048 (1981).
</li>

<li>
<code>'PZM'</code>: Perdew-Zunger LDA with modification to
improve sooth connection between high density form and low density one.
:J. P. Perdew and Alex Zunger, Phys. Rev. B 23, 5048 (1981).
</li>


<li>
<code>'TBmBJ'</code>: Tran-Blaha meta-GGA exchange with Perdew-Wang correlation.
:Fabien Tran and Peter Blaha, Phys. Rev. Lett. 102, 226401 (2008). 
John P. Perdew and Yue Wang, Phys. Rev. B 45, 13244 (1992).
</li>
</ul>

<dt>cval(:); <code>Real(8)</code></dt>
<dd>Mixing parameter in Tran-Blaha meta-GGA exchange potential. If <code>cval</code> is set to a minus value, the mixing-parameter computed
by the formula in the original paper [Phys. Rev. Lett. 102, 226401 (2008)].
Default is <code>'1.0'</code>.
</dd>
</dl>

## &rgrid
<dl>

<dt>dl(3); <code>Real(8)</code></dt>
<dd>Spacing of real-space grids. Unit of length can be chosen by
<code>&units/unit_length</code>.
This valiable cannot be set with 
<code>&rgrid/num_rgrid</code>.
</dd>


<dt>num_rgrid(3); <code>Integer</code></dt>
<dd>Number of real-space grids.
This valiable cannot be set with 
<code>&rgrid/dl</code>.
</dd>

</dl>

## &kgrid
<dl>

<dt>num_kgrid(3); <code>Integer</code></dt>
<dd>Number of grids discretizing
the Brillouin zone.
</dd>

<dt>file_kw; <code>Character</code></dt>
<dd>
Name of a file for flexible k-point sampling.
This file will be read if <code>num_kgrid</code> are
all negative.
</dd>


</dl>

## &tgrid
<dl>

<dt>nt; <code>Integer</code></dt>
<dd>
Number of total time steps for real-time propagation.
</dd>

<dt>dt; <code>Real(8)</code></dt>
<dd>
Time step. Unit of time can be chosen by <code>
&units/unit_time
</code>.
</dd>


</dl>

## &propagation
<dl>

<dt>n_hamil; <code>Integer</code></dt>
<dd>
Order of Taylor expansion of a propagation operator.
Default is <code>4</code>.
</dd>

<dt>propagator; <code>Character</code></dt>
<dd>
Choice of Propagator.
<code>middlepoint</code> is an propagator
with the Hamiltoinan at midpoint of two-times.

<code>middlepoint</code> is Enforced time-reversal symmetry propagator.
Midpoint Order of Taylor expansion of a propagation operator 
[M.A.L. Marques, A. Castro, G.F. Bertsch, and A. Rubio, 
Comput. Phys. Commun., 151 60 (2003)].
Default is <code>middlepoint</code>.
</dd>

</dl>

## &scf
<dl>

<dt>ncg; <code>Integer</code></dt>
<dd>
Number of interation of Conjugate-Gradient method for each scf-cycle.
Default is <code>5</code>.
</dd>

<dt>nmemory_mb; <code>Integer</code></dt>
<dd>
Number of stored densities at previous scf-cycles for 
the modified-Broyden method.
Default is <code>8</code>.
</dd>

<dt>nmemory_mb; <code>Integer</code></dt>
<dd>
Parameter of the modified-Broyden method.
Default is <code>0.75</code>.
</dd>

<dt>fsset_option; <code>Character</code></dt>
<dd>
Probably, we should remove this function
since we can replace it with occupaion smoothing with temepratre.
</dd>

<dt>nfsset_start; <code>Integer</code></dt>
<dd>
Probably, we should remove this function
since we can replace it with occupaion smoothing with temepratre.
</dd>

<dt>nfsset_every; <code>Integer</code></dt>
<dd>
Probably, we should remove this function
since we can replace it with occupaion smoothing with temepratre.
</dd>

<dt>nscf; <code>Integer</code></dt>
<dd>
Number of maximum scf cycle.
</dd>

<dt>ngeometry_opt; <code>Integer</code></dt>
<dd>
Number of iteration of geometry optimization.
</dd>

<dt>subspace_diagonalization; <code>Character</code></dt>
<dd>
<dd>Enable(<code>'y'</code>)/disable(<code>'n'</code>) 
subspace diagonalization during scf cycle.
</dd>

<dt>cmixing; <code>Character</code></dt>
<dd>
Methods for densiy/potential mixing for scf cycle.
Default is <code>broyden</code>.
</dd>

<dt>rmixrate; <code>Real(8)</code></dt>
<dd>
Mixingratio. Default is <code>0.5</code>
</dd>

<dt>convergence; <code>Character</code></dt>
<dd>
Choice of quantity that will be used for convergence check in scf calculation.
<code>'rho'</code> and <code>'vh'</code> can be chosen.
Default is <code>'rho'</code>.
</dd>

<dt>threshold; <code>Real(8)</code></dt>
<dd>
Threshold for convergence check.
Default is <code>1d-6</code>
</dd>

</dl>

## &emfield
<dl>

<dt>trans_longi; <code>Character</code></dt>
<dd>
Geometry of solid-state calculations.
Transverse <code>'tr'</code> and longitudinal <code>'lo'</code>
can be chosen.
Default is <code>'tr'</code>.
</dd>

<dt>ae_shape1/ae_shape2; <code>Character</code></dt>
<dd>
Shape of the first/second pulse.
<ul>
<li>
<code>'impulse'</code>: Impulsive fields.
</li>

<li>
<code>'Asin2cos'</code>:
</li>

<li>
<code>'Asin4cos'</code>:
</li>

</ul>
</dd>


<dt>amplitude1/amplitude2; <code>Real(8)</code></dt>
<dd>
Amplitude of electric fields for the first/second pulse.
This valiable has the dimension of electric field, energy/(length*charge).
</dd>

<dt>rlaser_int1/rlaser_int2; <code>Real(8)</code></dt>
<dd>
Peak laser intensity (W/cm^2) the first/second pulse.
</dd>

<dt>pulse_tw1/pulse_tw2; <code>Real(8)</code></dt>
<dd>
Duration of the first/second pulse. Unit of time can be chosend 
by <code>&units/unit_time</code>.
</dd>

<dt>omega1/omega2; <code>Real(8)</code></dt>
<dd>
Mean photon energy of the first/second pulse. Unit of energy can be chosend 
by <code>&units/unit_energy</code>.
</dd>

<dt>epdir_re1(3)/epdir_re2(3); <code>Real(8)</code></dt>
<dd>
Real part of polarization vector the first/second pulse.
</dd>

<dt>epdir_im1(3)/epdir_im2(3); <code>Real(8)</code></dt>
<dd>
Imaginary part of polarization vector the first/second pulse.
</dd>

<dt>phi_cep1/phi_cep2; <code>Real(8)</code></dt>
<dd>
Carrier emvelope phase of the first/second pulse.
</dd>

<dt>t1_t2; <code>Real(8)</code></dt>
<dd>
Time-delay between the first and the second pulses.
Unit of time can be chosen by <code>&units/unit_time</code>.
</dd>

<dt>quadrupole; <code>Character</code></dt>
<dd>
Quadrupole potential can be employed if 
<code>quadrupole</code> is set to <code>'y'</code>.
Default is <code>'n'</code>.
</dd>

<dt>quadrupole_pot; <code>Character</code></dt>
<dd>
Form of a quadrupole potential.
</dd>

<dt>alocal_laser; <code>Character</code></dt>
<dd>
The pulse is applied to a specific domain.
Default is <code>'n'</code>.
</dd>

<dt>rlaserbound_sta(3)/rlaserbound_end(3); <code>Real(8)</code></dt>
<dd>
The edge of the domain where the pulse is applied.
These parameters are effective only when <code>alocal_laser</code> is <code>'y'</code>.
Default is <code>-1d7/1d7</code> in atomic unit.
Unit of length can be chosen by <code>&units/unit_length</code>.
</dd>

</dl>

## &linear_response
<dl>

<dt>e_impulse; <code>Real(8)</code></dt>
<dd>
Momentum of impulsive perturbation.
This valiable has the dimention of momentum, energy*time/length.
Defalt value is <code>5d-5</code> a.u.
</dd>
</dl>

## &multiscale
<dl>


<dt>fdtddim; <code>Character</code></dt>
<dd>
Dimension of FDTD calculation for multi-scale Maxwell-Kohn-Sham method.
Defalt value is <code>'1D'</code>.
</dd>

<dt>twod_shape; <code>Character</code></dt>
<dd>
Boundary condision of the second dimension for FDTD calculation with 
multi-scale Maxwell-Kohn-Sham method.
Defalt value is <code>'periodic'</code>.
</dd>

<dt>nx_m/ny_m/nz_m; <code>Integer</code></dt>
<dd>
Number of macroscopic grid points inside materials for (x/y/z)-direction.
</dd>

<dt>hx_m/hy_m/hz_m; <code>Real(8)</code></dt>
<dd>
Spacing of macroscopic grid points inside materials for (x/y/z)-direction.
Unit of length can be chosen by <code>&units/unit_length</code>.
</dd>

<dt>nksplit; <code>Integer</code></dt>
<dd>
Number of MPI processers that take care electron dynamics at each single macroscopic 
point. 
</dd>

<dt>nxysplit; <code>Integer</code></dt>
<dd>
Number of macroscopic points that will be taken care by a single MPI processer.
</dd>

<dt>nxvacl_m/nxvacr_m; <code>Integer</code></dt>
<dd>
Number of macroscopic grid points for vacumm region.
<code>nxvacl_m</code> gives the number for negative x-direction in front of material,
while
<code>nxvacr_m</code> gives the number for positive x-direction behind the material.
</dd>

</dl>

## &analysis
<dl>
<dt>projection_option; <code>Character</code></dt>
<dd>
Methods of projection.
<ul>
<li>
<code>'no'</code>: no projection.
</li>

<li>
<code>'gs'</code>: projection to eigenstates of ground-state Hamiltonian.
</li>

<li>
<code>'rt'</code>: projection to eigenstates of instantaneous Hamiltonian.
</li>
</ul>
</dd>

<dt>nenergy; <code>Integer</code></dt>
<dd>
Number of energy grids for analysis.
</dd>

<dt>de; <code>Real(8)</code></dt>
<dd>
Energy spacing for analysis.
Unit of energy can be chosen by <code>&units/unit_energy</code>
</dd>

<dt>out_psi; <code>Character</code></dt>
<dd>
If <code>'y'</code>, wavefunctions are output.
Default is <code>'n'</code>.
</dd>

<dt>out_dos; <code>Character</code></dt>
<dd>
If <code>'y'</code>, density of state is output.
Default is <code>'n'</code>.
</dd>

<dt>out_pdos; <code>Character</code></dt>
<dd>
If <code>'y'</code>, projected density of state is output.
Default is <code>'n'</code>.
</dd>

<dt>out_dns; <code>Character</code></dt>
<dd>
If <code>'y'</code>, density is output.
Default is <code>'n'</code>.
</dd>

<dt>out_dns_rt/out_dns_rt_step; <code>Character/Integer</code></dt>
<dd>
If <code>'y'</code>, density during real-time time-propagation is output
every <code>outdns_rt_step</code> time steps.
Default is <code>'n'</code>.
</dd>

<dt>out_elf_rt/out_elf_rt_step; <code>Character/Integer</code></dt>
<dd>
If <code>'y'</code>, electron-localization function 
during real-time time-propagation is output
every <code>out_elf_rt_step</code> time steps.
Default is <code>'n'</code>.
</dd>

<dt>format3d; <code>Character</code></dt>
<dd>
Format for three dimensional data.
<code>'avs'</code>, <code>'cube'</code>, and <code>'vtk'</code>
can be chosen.
Default is <code>'avs'</code>.
</dd>

</dl>

## &hartree
<dl>

<dt>meo; <code>Integer</code></dt>
<dd>
Order of multi-pole expansion for calculation of Hartree potential.
Default is <code>3</code>.
</dd>

<dt>num_pole_xyz(3); <code>Integer</code></dt>
<dd>
Number of multi-poles.
</dd>

</dl>

## &ewald
<dl>
<dt>newald; <code>Integer</code></dt>
<dd>
Parameter for Ewald method. 
Short-range part of Ewald sum is calculated within <code>newald</code>th
nearlist neighbor cells.
Default is <code>4</code>.
</dd>

<dt>aewald; <code>Real(8)</code></dt>
<dd>
Range separation parameter for Ewald method. 
Default is <code>0.5</code>.
</dd>

</dl>
