# Input variables of SALMON

- [&calculation](#calculation)
- [&control](#control)
- [&units](#units)
- [&parallel](#parallel)
- [&system](#system)
- [&atomic_red_coor](#atomic_red_coor)
- [&atomic_coor](#atomic_coor)
- [&pseudo](#pseudo)
- [&functional](#functional)
- [&rgrid](#rgrid)
- [&kgrid](#kgrid)
- [&tgrid](#tgrid)
- [&propagation](#propagation)
- [&scf](#scf)
- [&emfield](#emfield)
- [&analysis](#analysis)
- [&hartree](#hartree)
- [&ewald](#ewald)

## &calculation
<dl>
<dt>calc_mode; <code>Character</code>; 0d/3d</dt>
<dd>Choice of Calculation modes. <code>'GS'</code>,
<code>'RTLR'</code>,<code>'RTPulse'</code>, <code>'GS_RTLR'</code>, 
and <code>'GS_RTPulse'</code>
can be chosen.
</dd>

<dt>use_ehrenfest_md; <code>Character</code>; 0d/3d</dt>
<dd>
Enable(<code>'y'</code>)/disable(<code>'n'</code>) 
Ehrenfest dynamics.
Default is <code>'n'</code>.
</dd>

<dt>use_ms_maxwell; <code>Character</code>; 3d</dt>
<dd>
Enable(<code>'y'</code>)/disable(<code>'n'</code>) 
Multi-scale Maxwell-Kohn-Sham coupling. 
Default is <code>'n'</code> 
</dd>

<dt>use_force; <code>Character</code>; 0d</dt>
<dd>
Enable(<code>'y'</code>)/disable(<code>'n'</code>) 
force calculation.
Default is <code>'n'</code>.
</dd>

<dt>use_geometry_opt; <code>Character</code>; 0d</dt>
<dd>
Enable(<code>'y'</code>)/disable(<code>'n'</code>) 
geometry optimization.
Default is <code>'n'</code>.
</dd>

</dl>

## &control
<dl>

<dt>restart_option; <code>Character</code>; 3d</dt>
<dd>Flag for restart. <code>'new'</code> or 
<code>'restart'</code> can be chosen.
<code>'new'</code> is default.
</dd>

<dt>backup_frequency; <code>Integer</code>; 3d</dt>
<dd>Frequency of backup during the time-propagation. 
If <code>0</code> is set, the backup is not performed.
Default is <code>0</code>.
</dd>


<dt>time_shutdown; <code>Real(8)</code>; 3d</dt>
<dd>Timer for automatic shutdown. The unit is always second.
If negative time is chosen, the automatic shutdown is not performed.
Default is <code>-1</code> sec.
</dd>

<dt>sysname; <code>Character</code>; 0d/3d</dt>
<dd>Name of calculation. This is used for a prefix of output files.
Default is <code>default</code>.
</dd>


</dl>


## &units
<dl>

<dt>unit_time; <code>Character</code>; 0d/3d</dt>
<dd>Unit of time for input variables. 
Atomic unit <code>'a.u.'</code> and 
femtosecond <code>'fs'</code> can be chosen.
Default is <code>'a.u.'</code>.
</dd>

<dt>unit_length; <code>Character</code>; 0d/3d</dt>
<dd>Unit of length for input variables. 
Atomic unit <code>'a.u.'</code> and 
Aungstrom <code>angstrom</code> can be chosen.
Default is <code>'a.u.'</code>.
</dd>

<dt>unit_energy; <code>Character</code>; 0d/3d</dt>
<dd>Unit of energy for input variables. 
Atomic unit <code>'a.u.'</code> and 
electron-volt <code>'ev'</code> can be chosen.
Default is <code>'a.u.'</code>.
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

<dt>nproc_ob; <code>Integer</code>; 0d</dt>
<dd>Number of MPI parallelization for orbitals.
Default is <code>0</code>.
</dd>

<dt>nproc_domain(3); <code>Integer</code>; 0d</dt>
<dd>Number of MPI parallelization for each direction in 
real-space.
Default is <code>(0,0,0)</code>.
</dd>

<dt>nproc_domain_s(3); <code>Integer</code>; 0d</dt>
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

<dt>iperiodic; <code>Integer</code>; 0d/3d</dt>
<dd>Dimension for periodic boundary condition.
<code>0</code> is for isolated systems, and 
<code>3</code> is for solids.
Default is <code>0</code>.
</dd>

<dt>ispin; <code>Integer</code>; 0d</dt>
<dd>Variable for classification of closed shell systems and open shell systems.
<code>0</code> is for closed shell systems, and
<code>1</code> is for open shell systems.
Default is <code>0</code>
</dd>

<dt>al(3); <code>Real(8)</code>; 0d/3d</dt>
<dd>Lattice constants. Unit of the length can be 
chosen by <code>&units/unit_length</code>.
</dd>

<dt>isym; <code>Integer</code>; 3d</dt>
<dd>Number of symmetries that can be used for 
reduction of k-points.
Default is <code>0</code>.
</dd>

<dt>crystal_structure; <code>Character</code>; 3d</dt>
<dd>Name of symmetry that can be used for the resuction of # of k-points.
Default is <code>'none'</code>.
</dd>

<dt>nstate; <code>Integer</code>; 0d/3d</dt>
<dd>Number of states/bands.
</dd>

<dt>nstate_spin(2); <code>Integer</code>; 0d</dt>
<dd>Number of states/bands can be specified independently
by <code>nstate_spin(1)/nstate_spin(2)</code>.
This option is incompatible with <code>nstate</code>
</dd>

<dt>nelec; <code>Integer</code>; 0d/3d</dt>
<dd>Number of valence electrons.
</dd>

<dt>nelec_spin(2); <code>Integer</code>; 0d</dt>
<dd>Number of up/down-spin electrons can be specified independently
by <code>nelec_spin(1)/nelec_spin(2)</code>.
This option is incompatible with <code>nelec</code>
</dd>

<dt>temperature; <code>Real(8)</code>; 3d</dt>
<dd>Temperature of electrons.
Unit of the energy can be chosen <code>&units/unit_energy</code>.
</dd>


<dt>nelem; <code>Integer</code>; 0d/3d</dt>
<dd>Number of elements that will be used in calculations.
</dd>

<dt>natom; <code>Integer</code>; 0d/3d</dt>
<dd>Number of atoms in a calculation cell.
</dd>

<dt>file_atom_red_coor; <code>Character</code></dt>
<dd>File name of atomic positions. In this file, 
the atomic coordinates can be written in reduced coordinates.
This option is incompatible with 
<code>&system/file_atom_coor</code>,
<code>&atomic_coor</code>, and 
<code>&atomic_red_coor</code>.
</dd>

<dt>file_atom_coor; <code>Character</code>; 0d</dt>
<dd>File name of atomic positions. In this file, 
the atomic coordinates can be written in Cartesian cooridnates.
The unit of the length can be chosen by 
<code>&units/unit_length</code>.
This option is incompatible with 
<code>&system/file_atom_red_coor</code>,
<code>&atomic_coor</code>, and 
<code>&atomic_red_coor</code>.
</dd>


</dl>

## &atomic_red_coor
In &atomic_red_coor, positions of atoms can be written in reduced coordinates
as follows: <br>
` 'Si'	0.00	0.00	0.00	1 ` <br>
`	'Si'	0.25	0.25	0.25	1 ` <br>
` ... ` <br>
Here, the information of atoms is ordered in row. For example, the first row gives
the information of the first atom. The number of rows must be equal to 
`&system/natom`.
The first coloum can be any caracters and does not affect calculations.
The second, third and fourth columns are reduced coordinates for
the first, second and third directions, respectively. 
The fifth column is a serial number of the spieces, which is used in 
`&pseudo`.
This option is incompatible with 
<code>&system/file_atom_red_coor</code>,
<code>&system/file_atom_coor</code>, and
<code>&atomic_coor</code>.


## &atomic_coor
In &atomic_coor, positions of atoms can be written in Cartesian coordinates.
The structure is same as &atomic_red_coor.
The unit of the length can be chosen by 
<code>&units/unit_length</code>.
This option is incompatible with 
<code>&system/file_atom_red_coor</code>,
<code>&system/file_atom_coor</code>, and
<code>&atomic_red_coor</code>.


## &pseudo
<dl>

<dt>pseudo_file(:); <code>Character</code>; 0d/3d</dt>
<dd>Name of pseudopotential files.
</dd>

<dt>Lmax_ps(:); <code>Integer</code>; 0d/3d</dt>
<dd>Maximum angular momentum of pseudopotential projectors.
</dd>

<dt>Lloc_ps(:); <code>Integer</code>; 0d/3d</dt>
<dd>Angular momentum of pseudopotential thant will be treated as local.
</dd>

<dt>iZatom(:); <code>Integer</code>; 0d/3d</dt>
<dd>Atomic number.
</dd>


<dt>psmask_option(:); <code>Character</code>; 3d</dt>
<dd>Enable(<code>'y'</code>)/disable(<code>'n'</code>) 
Fourier filtering for pseudopotentials. 
Default is <code>'n'</code>.
</dd>

<dt>alpha_mask(:); <code>Real(8)</code>; 3d</dt>
<dd>Parameter for the Fourier filtering for pseudopotential.
Default is <code>'0.8'</code>.
</dd>

<dt>gamma_mask(:); <code>Real(8)</code>; 3d</dt>
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

<dt>xc; <code>Character</code>; 3d</dt>
<dd>
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
</dd>

<dt>cval(:); <code>Real(8)</code>; 3d</dt>
<dd>Mixing parameter in Tran-Blaha meta-GGA exchange potential. If <code>cval</code> is set to a minus value, the mixing-parameter computed
by the formula in the original paper [Phys. Rev. Lett. 102, 226401 (2008)].
Default is <code>'1.0'</code>.
</dd>
</dl>

## &rgrid
<dl>

<dt>dl(3); <code>Real(8)</code>; 0d</dt>
<dd>Spacing of real-space grids. Unit of length can be chosen by
<code>&units/unit_length</code>.
This valiable cannot be set with 
<code>&rgrid/num_rgrid</code>.
</dd>


<dt>num_rgrid(3); <code>Integer</code>; 3d</dt>
<dd>Number of real-space grids.
This valiable cannot be set with 
<code>&rgrid/dl</code>.
</dd>

</dl>

## &kgrid
<dl>

<dt>num_kgrid(3); <code>Integer</code>; 3d</dt>
<dd>Number of grids discretizing
the Brillouin zone.
</dd>

<dt>file_kw; <code>Character</code>; 3d</dt>
<dd>
Name of a file for flexible k-point sampling.
This file will be read if <code>num_kgrid</code> are
all negative.
</dd>


</dl>

## &tgrid
<dl>

<dt>nt; <code>Integer</code>; 0d/3d</dt>
<dd>
Number of total time steps for real-time propagation.
</dd>

<dt>dt; <code>Real(8)</code>; 0d/3d</dt>
<dd>
Time step. Unit of time can be chosen by <code>
&units/unit_time
</code>.
</dd>


</dl>

## &propagation
<dl>

<dt>n_hamil; <code>Integer</code>; 0d</dt>
<dd>
Order of Taylor expansion of a propagation operator.
Default is <code>4</code>.
</dd>

<dt>propagator; <code>Character</code>; 3d</dt>
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

<dt>amin_routine; <code>Character</code>; 0d</dt>
<dd>
Minimization routine for the ground state calculation. 
<code>'cg'</code>, <code>'diis'</code>, and <code>'cg-diis'</code> can be chosen.
Default is <code>'cg'</code>.
</dd>

<dt>ncg; <code>Integer</code>; 0d/3d</dt>
<dd>
Number of interation of Conjugate-Gradient method for each scf-cycle.
Default is <code>5</code>.
</dd>

<dt>amixing; <code>Character</code>; 0d</dt> 
<dd>
Methods for density/potential mixing for scf cycle. <code>simple</code> and <code>broyden</code> can be chosen.
Default is <code>broyden</code>.
</dd>

<dt>rmixrate; <code>Real(8)</code>; 0d</dt>
<dd>
Mixing ratio for simple mixing. Default is <code>0.5</code>.
</dd>

<dt>nmemory_mb; <code>Integer</code>; 0d/3d</dt>
<dd>
Number of stored densities at previous scf-cycles for 
the modified-Broyden method. Default is <code>8</code>. 
If <code>&system/iperiodic</code> is <code>0</code>, <code>nmemory_mb</code> must be less than 21.
</dd>

<dt>alpha_mb; <code>Real(8)</code>; 0d/3d</dt>
<dd>
Parameter of the modified-Broyden method.
Default is <code>0.75</code>.
</dd>

<dt>fsset_option; <code>Character</code>; 3d</dt>
<dd>
Probably, we should remove this function
since we can replace it with occupaion smoothing with temepratre.
</dd>

<dt>nfsset_start; <code>Integer</code>; 3d</dt>
<dd>
Probably, we should remove this function
since we can replace it with occupaion smoothing with temepratre.
</dd>

<dt>nfsset_every; <code>Integer</code>; 3d</dt>
<dd>
Probably, we should remove this function
since we can replace it with occupaion smoothing with temepratre.
</dd>

<dt>nscf; <code>Integer</code>; 0d/3d</dt>
<dd>
Number of maximum scf cycle.
</dd>

<dt>ngeometry_opt; <code>Integer</code>; 0d</dt>
<dd>
Number of iteration of geometry optimization.
</dd>

<dt>subspace_diagonalization; <code>Character</code>; 0d</dt>
<dd>
<dd>Enable(<code>'y'</code>)/disable(<code>'n'</code>) 
subspace diagonalization during scf cycle.
</dd>

<dt>convergence; <code>Character</code>; 0d</dt>
<dd>
Choice of quantity that is used for convergence check in a scf calculation. 
Default is <code>'rho_dng'</code>. The following can be chosen.
<ul>
<li>
<code>'rho'</code>: Convergence is checked by ||rho(i)-rho(i-1)||<sup>2</sup>, where i is an iteration number of the scf calculation.
</li>

<li>
<code>'rho_dng'</code>: Convergence is checked by ||rho(i)-rho(i-1)||<sup>2</sup>/(number of grids). "dng" means "devided by number of grids".
</li>

<li>
<code>'pot'</code>: Convergence is checked by ||Vlocal(i)-Vlocal(i-1)||<sup>2</sup>, where Vlocal is Vh + Vxc + Vps_local.
</li>

<li>
<code>'pot_dng'</code>: Convergence is checked by ||Vlocal(i)-Vlocal(i-1)||<sup>2</sup>/(number of grids).
</li>
</ul>
</dd>

<dt>threshold; <code>Real(8)</code>; 0d</dt>
<dd>
Threshold for convergence check that is used when either <code>'rho'</code> or <code>'rho_dng'</code> is specified.
Default is <code>1d-17</code> (= 6.75d-17Å<sup>-3</sup>)
</dd>

<dt>threshold_pot; <code>Real(8)</code>; 0d</dt>
<dd>
Threshold for convergence check that is used when either <code>'pot'</code> or <code>'pot_dng'</code> is specified. <code>threshold_pot</code> must be set when either <code>'pot'</code> or <code>'pot_dng'</code> is specified.
Default is <code>-1d0</code> (1 a.u.= 1.10d2 Å<sup>3</sup>eV<sup>2</sup>)
</dd>

</dl>

## &emfield
<dl>

<dt>trans_longi; <code>Character</code>; 3d</dt>
<dd>
Geometry of solid-state calculations.
Transverse <code>'tr'</code> and longitudinal <code>'lo'</code>
can be chosen.
Default is <code>'tr'</code>.
</dd>

<dt>ae_shape1/ae_shape2; <code>Character</code>; 0d/3d</dt>
<dd>
Shape of the first/second pulse.
<ul>
<li>
<code>'impulse'</code>: Impulsive fields.
</li>

<li>
<code>'Acos2'</code>,
<code>'Acos3'</code>,
<code>'Acos4'</code>,
<code>'Acos6'</code>, and
<code>Acos8'</code>: Envelope of cos<sup>2</sup>,cos<sup>3</sup>,cos<sup>4</sup>
cos<sup>6</sup>, and cos<sup>8</sup> for a vector potential.
</li>

<li>
<code>'Esin2cos'</code>:
</li>

<li>
<code>'Asin2cos'</code>:
</li>

<li>
<code>'Asin4cos'</code>:
</li>

</ul>
</dd>

<dt>e_impulse; <code>Real(8)</code>; 0d/3d</dt>
<dd>
Momentum of impulsive perturbation.
This valiable has the dimention of momentum, energy*time/length.
Defalt value is <code>1d-2</code> a.u.
</dd>

<dt>amplitude1/amplitude2; <code>Real(8)</code>; 0d/3d</dt>
<dd>
Amplitude of electric fields for the first/second pulse.
This valiable has the dimension of electric field, energy/(length*charge).
</dd>

<dt>rlaser_int1/rlaser_int2; <code>Real(8)</code>; 0d/3d</dt>
<dd>
Peak laser intensity (W/cm<sup>2</sup>) the first/second pulse.
</dd>

<dt>pulse_tw1/pulse_tw2; <code>Real(8)</code>; 0d/3d</dt>
<dd>
Duration of the first/second pulse. Unit of time can be chosend 
by <code>&units/unit_time</code>.
</dd>

<dt>omega1/omega2; <code>Real(8)</code>; 0d/3d</dt>
<dd>
Mean photon energy of the first/second pulse. Unit of energy can be chosend 
by <code>&units/unit_energy</code>.
</dd>

<dt>epdir_re1(3)/epdir_re2(3); <code>Real(8)</code>; 0d/3d</dt>
<dd>
Real part of polarization vector the first/second pulse.
</dd>

<dt>epdir_im1(3)/epdir_im2(3); <code>Real(8)</code>; 0d/3d</dt>
<dd>
Imaginary part of polarization vector the first/second pulse.
</dd>

<dt>phi_cep1/phi_cep2; <code>Real(8)</code>; 0d/3d</dt>
<dd>
Carrier emvelope phase of the first/second pulse.
Default is <code>0d0/0d0</code>.
</dd>

<dt>t1_t2; <code>Real(8)</code>; 0d/3d</dt>
<dd>
Time-delay between the first and the second pulses.
Unit of time can be chosen by <code>&units/unit_time</code>.
</dd>

<dt>quadrupole; <code>Character</code>; 0d</dt>
<dd>
Quadrupole potential can be employed if 
<code>quadrupole</code> is set to <code>'y'</code>.
Default is <code>'n'</code>.
</dd>

<dt>quadrupole_pot; <code>Character</code>; 0d</dt>
<dd>
Form of a quadrupole potential.
</dd>

<dt>alocal_laser; <code>Character</code>; 0d</dt>
<dd>
The pulse is applied to a specific domain.
Default is <code>'n'</code>.
</dd>

<dt>rlaserbound_sta(3)/rlaserbound_end(3); <code>Real(8)</code>; 0d</dt>
<dd>
The edge of the domain where the pulse is applied.
These parameters are effective only when <code>alocal_laser</code> is <code>'y'</code>.
Default is <code>-1d7/1d7</code> in atomic unit.
Unit of length can be chosen by <code>&units/unit_length</code>.
</dd>

</dl>

## &multiscale
<dl>

<dt>fdtddim; <code>Character</code>; 3d</dt>
<dd>
Dimension of FDTD calculation for multi-scale Maxwell-Kohn-Sham method.
Defalt value is <code>'1D'</code>.
</dd>

<dt>twod_shape; <code>Character</code>; 3d</dt>
<dd>
Boundary condision of the second dimension for FDTD calculation with 
multi-scale Maxwell-Kohn-Sham method.
Defalt value is <code>'periodic'</code>.
</dd>

<dt>nx_m/ny_m/nz_m; <code>Integer</code>; 3d</dt>
<dd>
Number of macroscopic grid points inside materials for (x/y/z)-direction.
</dd>

<dt>hx_m/hy_m/hz_m; <code>Real(8)</code>; 3d</dt>
<dd>
Spacing of macroscopic grid points inside materials for (x/y/z)-direction.
Unit of length can be chosen by <code>&units/unit_length</code>.
</dd>

<dt>nksplit; <code>Integer</code>; 3d</dt>
<dd>
Number of MPI processers that take care electron dynamics at each single macroscopic 
point. 
</dd>

<dt>nxysplit; <code>Integer</code>; 3d</dt>
<dd>
Number of macroscopic points that will be taken care by a single MPI processer.
</dd>

<dt>nxvacl_m/nxvacr_m; <code>Integer</code>; 3d</dt>
<dd>
Number of macroscopic grid points for vacumm region.
<code>nxvacl_m</code> gives the number for negative x-direction in front of material,
while
<code>nxvacr_m</code> gives the number for positive x-direction behind the material.
</dd>

</dl>

## &analysis
<dl>
<dt>projection_option; <code>Character</code>; 3d</dt>
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

<dt>nenergy; <code>Integer</code>; 0d/3d</dt>
<dd>
Number of energy grids for analysis.
</dd>

<dt>de; <code>Real(8)</code>; 0d/3d</dt>
<dd>
Energy spacing for analysis.
Unit of energy can be chosen by <code>&units/unit_energy</code>
</dd>

<dt>out_psi; <code>Character</code>; 0d</dt>
<dd>
If <code>'y'</code>, wavefunctions are output.
Default is <code>'n'</code>.
</dd>

<dt>out_dos; <code>Character</code>; 0d/3d</dt>
<dd>
If <code>'y'</code>, density of state is output.
Default is <code>'n'</code>.
</dd>

<dt>out_pdos; <code>Character</code>; 0d</dt>
<dd>
If <code>'y'</code>, projected density of state is output.
Default is <code>'n'</code>.
</dd>

<dt>out_dns; <code>Character</code>; 0d/3d</dt>
<dd>
If <code>'y'</code>, density is output.
Default is <code>'n'</code>.
</dd>

<dt>out_elf; <code>Character</code>; 0d</dt>
<dd>
If <code>'y'</code>, electron localization function is output.
Default is <code>'n'</code>.
</dd>

<dt>out_dns_rt/out_dns_rt_step; <code>Character/Integer</code>; 0d/3d</dt>
<dd>
If <code>'y'</code>, density during real-time time-propagation is output
every <code>outdns_rt_step</code> time steps.
Default is <code>'n'</code>.
</dd>

<dt>out_elf_rt/out_elf_rt_step; <code>Character/Integer</code>; 0d</dt>
<dd>
If <code>'y'</code>, electron localization function 
during real-time time-propagation is output
every <code>out_elf_rt_step</code> time steps.
Default is <code>'n'</code>.
</dd>

<dt>out_estatic_rt/out_estatic_rt_step; <code>Character/Integer</code>; 0d</dt>
<dd>
If <code>'y'</code>, static electric field
during real-time time-propagation is output
every <code>out_estatic_rt_step</code> time steps.
Default is <code>'n'</code>.
</dd>

<dt>format3d; <code>Character</code>; 0d/3d</dt>
<dd>
Format for three dimensional data.
<code>'avs'</code>, <code>'cube'</code>, and <code>'vtk'</code>
can be chosen.
Default is <code>'cube'</code>.
</dd>

<dt>numfiles_out_3d; <code>Integer</code>; 0d</dt>
<dd>
Number of separated files for three dimensional data.
Effective only when <code>format3d</code> is <code>'avs'</code>.
<code>numfiles_out_3d</code> must be less than or equal to number of processes.
Default is <code>1</code>.
</dd>

</dl>

## &hartree
<dl>

<dt>meo; <code>Integer</code>; 0d</dt>
<dd>
A variable to determine how to put multipoles in the Hartree potential calculation. Default is <code>3</code>.
<ul>
<li>
<code>1</code>: A single pole is put at the center.
</li>

<li>
<code>2</code>: Multipoles are put at the center of atoms.
</li>

<li>
<code>3</code>: Multipoles are put at the center of mass of electrons in prepared cuboids.
</li>
</ul>
</dd>

<dt>num_pole_xyz(3); <code>Integer</code>; 0d</dt>
<dd>
Number of multipoles when <code>meo</code> is <code>3</code>. Default is <code>0,0,0</code>. When default is set, number of multipoles is calculated automatically.
</dd>

</dl>

## &ewald
<dl>
<dt>newald; <code>Integer</code>; 3d</dt>
<dd>
Parameter for Ewald method. 
Short-range part of Ewald sum is calculated within <code>newald</code>th
nearlist neighbor cells.
Default is <code>4</code>.
</dd>

<dt>aewald; <code>Real(8)</code>; 3d</dt>
<dd>
Range separation parameter for Ewald method. 
Default is <code>0.5</code>.
</dd>

</dl>




