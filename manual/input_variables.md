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
>- [&opt](#opt)
>- [&md](#md)

## &calculation
<dl>
<dt>calc_mode; <code>Character</code>; 0d/3d</dt>
<dd>Choice of Calculation modes. <code>'GS'</code>,
<code>'RTLR'</code>,<code>'RTPulse'</code>, <code>'GS_RTLR'</code>, 
and <code>'GS_RTPulse'</code>
can be chosen.
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

><dt>use_ehrenfest_md; <code>Character</code>; 0d/3d</dt>
><dd>
>Enable(<code>'y'</code>)/disable(<code>'n'</code>) 
>Ehrenfest molecular dynamics.
>Default is <code>'n'</code>.
></dd>

><dt>use_geometry_opt; <code>Character</code>; 0d/3d</dt>
><dd>
>Enable(<code>'y'</code>)/disable(<code>'n'</code>) 
>geometry optimization.
>Default is <code>'n'</code>.
></dd>

</dl>

## &control
<dl>

<dt>sysname; <code>Character</code>; 0d/3d</dt>
<dd>Name of calculation. This is used for a prefix of output files.
Default is <code>default</code>.
</dd>

<dt>directory; <code>Character</code>; 0d/3d</dt>
<dd>Name of a default directory, where the basic results will be written down.
Default is the current directoy, <code>./</code>.
</dd>

><dt>restart_option; <code>Character</code>; 3d</dt>
><dd>Flag for restart. <code>'new'</code> or 
><code>'restart'</code> can be chosen.
><code>'new'</code> is default.
></dd>

><dt>backup_frequency; <code>Integer</code>; 3d</dt>
><dd>Frequency of backup during the time-propagation. 
>If <code>0</code> is set, the backup is not performed.
>Default is <code>0</code>.
></dd>

><dt>time_shutdown; <code>Real(8)</code>; 3d</dt>
><dd>Timer for automatic shutdown. The unit is always second.
>If negative time is chosen, the automatic shutdown is not performed.
>Default is <code>-1</code> sec.
></dd>

><dt>dump_filename; <code>Character</code>; 3d</dt>
><dd>Name of a filename for the restart calculation.
></dd>

><dt>read_gs_wfn_k; <code>Character</code>; 3d</dt>
><dd>
>Read ground state wave function as initial guess (pre-calculated "gs_wfn_k" directory printed by <code>calc_mode=GS</code>) if this is <code>y</code>.
>Default is <code>n</code>. (The option is available for restart of ground state SCF calculation, geometrical optimization, etc. But the data is automatically read when "calc_mode=RT")
></dd>

><dt>write_gs_wfn_k; <code>Character</code>; 3d</dt>
><dd>
>Write ground state wave function into "gs_wfn_k" directory if this is <code>y</code>. (Useful for <code>calc_mode=GS_RT</code>. But the data is always written in the case of <code>calc_mode=GS</code> calculation.) 
>Default is <code>n</code>.
></dd>

><dt>modify_gs_wfn_k; <code>Character</code>; 3d</dt>
><dd>
>Option to modify initial guess wave function (pre-calculated "gs_wfn_k" directory) in combination with <code>read_gs_wfn_k = y</code>.
>If <code>copy_1stk_to_all</code> is set, the first k-point data file, wfn_gs_k0000001.wfn (supposed to be obtained by gamma-point calculation), is copied to all other points.
>Default is <code>n</code>.
></dd>

><dt>write_rt_wfn_k; <code>Character</code>; 3d</dt>
><dd>
>Write RT wave function into "rt_wfn_k" directory if this is <code>y</code>.
>Default is <code>n</code>.
></dd>

><dt>read_rt_wfn_k; <code>Character</code>; 3d</dt>
><dd>
>Read RT wave function (pre-calculated "rt_wfn_k" directory printed by <code>calc_mode=RT</code> with <code>write_rt_wfn_k=y</code>) if this is <code>y</code>. This is used for restarting <code>calc_mode=RT</code> (supporting <code>use_ehrenfest_md=y</code>, too: Coordinates and velocities of atoms for restarting are also included), then, "gs_wfn_k" directory is also necessary (actually used only for some analysis options). Note that, currently, field is taken over after restarting only if <code>ae_shape1=AcosX</code> type is used.
>Default is <code>n</code>.
></dd>

><dt>read_gs_wfn_k_ms; <code>Character</code>; 3d</dt>
><dd>
>Read ground state wave function as initial state for multiscale calculation. This should be used together with <code>use_ms_maxwell='y'</code>, <code>calc_mode='RT'</code> and <code>set_ini_coor_vel='y'</code>. The ground state wave function data must be pre-calculated for each configuration(or atomic coordinate) and be put the obtained directories 'gs_wfn_k' into specific directories: <code>directory</code>/multiscale/MXXXXXX/ where is the index number of the macro-grid point of the material region usually starting from '000001' up to the number of macro-grid point (the same place of 'ini_coor_vel.dat' used by the option <code>set_ini_coor_vel</code>).
>Default is <code>n</code>.
></dd>

><dt>read_rt_wfn_k_ms; <code>Character</code>; 3d</dt>
><dd> Read wave function and field information as initial state in multiscale calculation. These are generated in pre-calculation by using <code>write_rt_wfn_k_ms='y'</code>. If you give incident pulse from input file option, the field is added to the read initial data.
>Default is <code>n</code>.
></dd>

><dt>write_rt_wfn_k_ms; <code>Character</code>; 3d</dt>
><dd> Write wave function and field information at the last step in multiscale calculation. It is used for restarting by using <code>read_rt_wfn_k_ms='y'</code>. 
>Default is <code>n</code>.
></dd>


</dl>


## &units
<dl>

<dt>unit_system; <code>Character</code>; 0d/3d</dt>
<dd>Unit for input variables. 
If <code>'au'</code> or <code>'a.u.'</code>, atomic unit system is used. 
If <code>'A_eV_fs'</code>, Angstrom-eV-fs unit system is used. 
</dd>

</dl>

## &parallel
<dl>

<dt>nproc_ob/nproc_domain(3)/nproc_domain_s(3); <code>Integer</code>; 0d</dt>
<dd> Followings are explanation of each variable.
<ul>
<li>
<code>nproc_ob</code>: Number of MPI parallelization for orbitals that related to the wavefunction calculation.
</li>
<li>
<code>nproc_domain(3)'</code>: Number of MPI parallelization for each direction in real-space that related to the wavefunction calculation. 
</li>
<li>
<code>nproc_domain_s(3)'</code>: Number of MPI parallelization for each direction in real-space that related to the electron density calculation. 
</li>
</ul>
Defaults are <code>0</code> for <code>nproc_ob</code>, <code>(0/0/0)</code> for <code>nproc_domain</code>, and <code>(0/0/0)</code> for <code>nproc_domain_s</code>. If users use the defauls, automatic proccess assignment is done. Users can also specify <code>nproc_ob</code>, <code>nproc_domain</code>, and <code>nproc_domain_s</code> manually. In that case, followings must be satisfied. 
<ul>
<li>
<code>nproc_ob</code> * <code>nproc_domain(1)</code> * <code>nproc_domain(2)</code>* <code>nproc_domain(3)</code>=total number of processors
</li>
<li>
<code>nproc_domain_s(1)</code> * <code>nproc_domain_s(2)</code>* <code>nproc_domain_s(3)</code>=total number of processors
</li>
<li>
<code>nproc_domain_s(1)</code> is a multiple of <code>nproc_domain(1)</code>
</li>
<li>
<code>nproc_domain_s(2)</code> is a multiple of <code>nproc_domain(2)</code>
</li>
<li>
<code>nproc_domain_s(3)</code> is a multiple of <code>nproc_domain(3)</code>
</li>
</ul>
<dt>num_datafiles_in/num_datafiles_out; <code>Integer</code>; 0d</dt>
<dd>Number of input/output files for wavefunction.
Defaults are <code>1</code>. If <code>num_datafiles_in</code>/<code>num_datafiles_out</code> are 1, wave functions are read from/ written in a regular intermediate file. If <code>num_datafiles_in</code>/<code>num_datafiles_out</code> are larger than or equal to 2, the wave functions are read from/ written in separated intermediate data files, and number of files are equal to <code>num_datafiles_in</code>/<code>num_datafiles_out</code>. These variables must be equal to nth power of 2. (n: 0 or positive integer)
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

><dt>temperature; <code>Real(8)</code>; 3d</dt>
><dd>Temperature of electrons.
>Unit of the energy can be chosen <code>&units/unit_energy</code>.
></dd>


<dt>nelem; <code>Integer</code>; 0d/3d</dt>
<dd>Number of elements that will be used in calculations.
</dd>

<dt>natom; <code>Integer</code>; 0d/3d</dt>
<dd>Number of atoms in a calculation cell.
</dd>

><dt>file_atom_red_coor; <code>Character</code></dt>
><dd>File name of atomic positions. In this file, 
>the atomic coordinates can be written in reduced coordinates.
>This option is incompatible with 
><code>&system/file_atom_coor</code>,
><code>&atomic_coor</code>, and 
><code>&atomic_red_coor</code>.
></dd>

><dt>file_atom_coor; <code>Character</code>; 0d</dt>
><dd>File name of atomic positions. In this file, 
>the atomic coordinates can be written in Cartesian cooridnates.
>The unit of the length can be chosen by 
><code>&units/unit_length</code>.
>This option is incompatible with 
><code>&system/file_atom_red_coor</code>,
><code>&atomic_coor</code>, and 
><code>&atomic_red_coor</code>.
></dd>


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

<dt>lmax_ps(:); <code>Integer</code>; 0d/3d</dt>
<dd>Maximum angular momentum of pseudopotential projectors.
</dd>

<dt>lloc_ps(:); <code>Integer</code>; 0d/3d</dt>
<dd>Angular momentum of pseudopotential thant will be treated as local.
</dd>

<dt>izatom(:); <code>Integer</code>; 0d/3d</dt>
<dd>Atomic number.
</dd>


><dt>psmask_option(:); <code>Character</code>; 3d</dt>
><dd>Enable(<code>'y'</code>)/disable(<code>'n'</code>) 
>Fourier filtering for pseudopotentials. 
>Default is <code>'n'</code>.
></dd>

><dt>alpha_mask(:); <code>Real(8)</code>; 3d</dt>
><dd>Parameter for the Fourier filtering for pseudopotential.
>Default is <code>'0.8'</code>.
></dd>

><dt>gamma_mask(:); <code>Real(8)</code>; 3d</dt>
><dd>Parameter for the Fourier filtering for pseudopotential.
>Default is <code>'1.8'</code>.
></dd>

><dt>eta_maskk(:); <code>Real(8)</code></dt>
><dd>Parameter for the Fourier filtering for pseudopotential.
>Default is <code>'15.0'</code>.
></dd>

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

><dt>no_update_func; <code>character(1)</code>; 3d</dt>
><dd>Option not to update functional (or Hamiltonian) in RT time step, i.e., keep ground state Hamiltonian during time-evolution.
>Default is <code>'n'</code>.
></dd>

</dl>

## &rgrid
<dl>

<dt>dl(3); <code>Real(8)</code>; 0d/3d</dt>
<dd>Spacing of real-space grids. Unit of length can be chosen by
<code>&units/unit_length</code>.
This valiable cannot be set with 
<code>&rgrid/num_rgrid</code>.
If <code>&system/iperiodic</code> is set to <code>3</code>,
the actual grid spacing is automatically refined in calculations
so that the size of the simulation box
<code>&system/al(3)</code> becomes divisible by the spacing.
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
smaller than 1.
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

<dt>propagator; <code>Character</code>; 3d</dt>
<dd>
Choice of Propagator.
<code>middlepoint</code> is an propagator
with the Hamiltoinan at midpoint of two-times.

<code>etrs</code> is enforced time-reversal symmetry propagator.
[M.A.L. Marques, A. Castro, G.F. Bertsch, and A. Rubio, 
Comput. Phys. Commun., 151 60 (2003)].
Default is <code>middlepoint</code>.
</dd>

><dt>n_hamil; <code>Integer</code>; 0d</dt>
><dd>
>Order of Taylor expansion of a propagation operator.
>Default is <code>4</code>.
></dd>

</dl>

## &scf
<dl>

<dt>nscf; <code>Integer</code>; 0d/3d</dt>
<dd>
Number of maximum scf cycle.
</dd>

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

><dt>fsset_option; <code>Character</code>; 3d</dt>
><dd>
>Probably, we should remove this function
>since we can replace it with occupaion smoothing with temepratre.
></dd>

><dt>nfsset_start; <code>Integer</code>; 3d</dt>
><dd>
>Probably, we should remove this function
>since we can replace it with occupaion smoothing with temepratre.
></dd>

><dt>nfsset_every; <code>Integer</code>; 3d</dt>
><dd>
>Probably, we should remove this function
>since we can replace it with occupaion smoothing with temepratre.
></dd>

<dt>subspace_diagonalization; <code>Character</code>; 0d</dt>
<dd>
<dd>Enable(<code>'y'</code>)/disable(<code>'n'</code>) 
subspace diagonalization during scf cycle.
</dd>

<dt>convergence; <code>Character</code>; 0d/3d</dt>
<dd>
Choice of quantity that is used for convergence check in a scf calculation. 
Default is <code>'rho_dne'</code>. 
<ul>
<li>
<code>'rho_dne'</code>: Convergence is checked by sum_ix|rho(ix,iter)-rho(ix,iter-1)|dx/N, where iter is an iteration number of the scf calculation and N is <code>&system/nelec</code>, the number of the valence electrons.
</li>
</ul>

For isolated systems, the followings can also be chosen.
<ul>
<li>
<code>'norm_rho'</code>: Convergence is checked by the square of the norm of difference of density, ||rho_iter(ix)-rho_iter-1(ix)||<sup>2</sup>=sum_ix|rho(ix,iter)-rho(ix,iter-1)|<sup>2</sup>. 
</li>

<li>
<code>'norm_rho_dng'</code>: Convergence is checked by ||rho_iter(ix)-rho_iter-1(ix)||<sup>2</sup>/(number of grids). "dng" means "devided by number of grids".
</li>

<li>
<code>'norm_pot'</code>: Convergence is checked by ||Vlocal_iter(ix)-Vlocal_iter-1(ix)||<sup>2</sup>, where Vlocal is Vh + Vxc + Vps_local.
</li>

<li>
<code>'pot_dng'</code>: Convergence is checked by ||Vlocal_iter(ix)-Vlocal_iter-1(ix)||<sup>2</sup>/(number of grids).
</li>
</ul>
</dd>

<dt>threshold; <code>Real(8)</code>; 0d/3d</dt>
<dd>
Threshold for convergence check that is used when <code>'rho_dne'</code> is specified.
Default is <code>1d-17</code>. 
</dd>

<dt>threshold_norm_rho; <code>Real(8)</code>; 0d</dt>
<dd>
Threshold for convergence check that is used when either <code>'norm_rho'</code> or <code>'norm_rho_dng'</code> is specified. <code>threshold_norm_rho</code> must be set when either <code>'norm_rho'</code> or <code>'norm_rho_dng'</code> is specified.
Default is <code>-1d0</code> a.u. (1 a.u.= 45.54Å<sup>-6</sup>)
</dd>

<dt>threshold_norm_pot; <code>Real(8)</code>; 0d</dt>
<dd>
Threshold for convergence check that is used when either <code>'norm_pot'</code> or <code>'norm_pot_dng'</code> is specified. <code>threshold_norm_pot</code> must be set when either <code>'norm_pot'</code> or <code>'norm_pot_dng'</code> is specified.
Default is <code>-1d0</code> a.u. (1 a.u.= 33.72x10<sup>4</sup> Å<sup>-6</sup>eV<sup>2</sup>)
</dd>

<dt>omp_loop; <code>Character</code>; 3d</dt>
<dd>
Loop for OpenMP parallelization if periodic boundary system is used. 
<li>
<code>k</code>: parallelization by k-point loop (Default).
</li>
<li>
<code>b</code>: parallelization mainly by band orbital loop (sometimzes space grid loop too). This works efficiently if the number of k-point to be calculated in a node is small (e.x. the case of single k-point for each node)
</li>
</dd>

><dt>skip_gsortho; <code>Character</code>; 3d</dt>
><dd>
>Flag to skip Gram-Schmidt orthogonalization in CG loop if periodic boundary system is used. If this is skipped the more iteration number is necessary to get convergence but each iteration step gets faster. If <code>omp_loop=b</code>, this flag is always applied.
>Default is <code>n</code>
></dd>

><dt>ngeometry_opt; <code>Integer</code>; 0d</dt>
><dd>
>Number of iteration of geometry optimization.
></dd>

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
<code>'Acos2'</code>: Envelope of cos<sup>2</sup> for a vector potential.
</li>
<li>
<code>'Ecos2'</code>: Envelope of cos<sup>2</sup> for a scalar potential.
</li>
</ul>
If <code>&system/iperiodic</code> is <code>3</code>, following can be also chosen,
<ul>
<li>
<code>'Acos3'</code>,
<code>'Acos4'</code>,
<code>'Acos6'</code>, and
<code>'Acos8'</code>: Envelopes of cos<sup>3</sup>,cos<sup>4</sup>
cos<sup>6</sup>, and cos<sup>8</sup> for vector potentials.
</li>
</ul>
and <code>'Esin2sin'</code>, <code>'Asin2cos'</code>, <code>'Asin2cw'</code>, 
<code>'input'</code>, and <code>'none'</code> can be also chosen but explanation is skipped.

</ul>
</dd>

<dt>e_impulse; <code>Real(8)</code>; 0d/3d</dt>
<dd>
Momentum of impulsive perturbation.
This valiable has the dimention of momentum, energy*time/length.
Default value is <code>1d-2</code> a.u.
</dd>

<dt>amplitude1/amplitude2; <code>Real(8)</code>; 0d/3d</dt>
<dd>
Amplitude of electric fields for the first/second pulse.
This valiable has the dimension of electric field, energy/(length*charge).
</dd>

<dt>rlaser_int_wcm2_1/rlaser_int_wcm2_2; <code>Real(8)</code>; 0d/3d</dt>
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

<dt>t1_delay; <code>Real(8)</code>; 3d</dt>
<dd>
Time-delay of the first pulses.
Unit of time can be chosen by <code>&units/unit_time</code>.
Default is <code>0d0</code>.
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

><dt>alocal_laser; <code>Character</code>; 0d</dt>
><dd>
>The pulse is applied to a specific domain.
>Default is <code>'n'</code>.
></dd>

><dt>rlaserbound_sta(3)/rlaserbound_end(3); <code>Real(8)</code>; 0d</dt>
><dd>
>The edge of the domain where the pulse is applied.
>These parameters are effective only when <code>alocal_laser</code> is <code>'y'</code>.
>Default is <code>-1d7/1d7</code> in atomic unit.
>Unit of length can be chosen by <code>&units/unit_length</code>.
></dd>

</dl>

## &multiscale
<dl>

<dt>fdtddim; <code>Character</code>; 3d</dt>
<dd>
Dimension of FDTD calculation for multi-scale Maxwell-Kohn-Sham method.
<code>1d,2d,3d</code> is available in the present version.
Default value is <code>'1D'</code>. 
</dd>

<dt>twod_shape; <code>Character</code>; 3d</dt>
<dd>
Boundary condision of the second dimension for FDTD calculation with 
multi-scale Maxwell-Kohn-Sham method.
Default value is <code>'periodic'</code>.
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

<dt>nxvacl_m/nxvacr_m; <code>Integer</code>; 3d</dt>
<dd>
Number of macroscopic grid points for vacumm region.
<code>nxvacl_m</code> gives the number for negative x-direction in front of material,
while
<code>nxvacr_m</code> gives the number for positive x-direction behind the material.
</dd>

<dt>nx_origin_m/ny_origin_m/nz_origin_m; <code>Integer</code>; 3d</dt>
<dd>
Origin coordinat of the grid points.
Default value is <code>'1'</code>.
</dd>

><dt>set_ini_coor_vel; <code>Character</code>; 3d</dt>
><dd>
>Set initial atomic coordinates and velocities for each macro-grid point. This must be given with specific directories and files: 
>Prepare <code>directory</code>/multiscale/MXXXXXX/ini_coor_vel.dat, where 'XXXXXX' is the index number of the macro-grid point of the material region usually starting from '000001' up to the number of macro-grid point. The format of the file 'ini_coor_vel.dat' is just Rx, Ry, Rz, Vx, Vy, Vz (with space separation) for each atom (i.e. for each line), where the unit of the coordinates, Rx, Ry, Rz, is angstrom or a.u. speficied by <code>unit_system</code> but that of velocities is always a.u.. This option should be used together with <code>read_gs_wfn_k_ms</code> which is the option to read the ground state wave function for each macro-grid point. 
>Default value is <code>'y'</code>.
></dd>

><dt>nmacro_write_group; <code>Integer</code>; 3d</dt>
><dd>
>If the number of macroscopic grids are very large, computers can be unstable by writing all information of all macroscopic grid points at the same time. To avoid that, the writings are divided by specifying this option. Writings will be done by each <code>nmacro_write_group</code> macroscopic grid points. (this number must be aliquot part of the total number of macroscopic grid points)
>Default value is <code>'-1'</code>.
></dd>

><dt>file_macropoint; <code>Character</code>; 3d</dt>
><dd>
>If this variable is specified, the coordinates of the macropoints are set from the file.
>Default value is <code>''</code>.
></dd>

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

><dt>projection_decomp; <code>Character</code>; 3d</dt>
><dd>
>If <code>'atom'</code> combined with <code>projection_option='gs'</code>, 
>the number of excited electron is decomposed into each atom 
>(this is printed in <code>SYSname</code>_nex_atom.data)
>Default is <code>'n'</code>.
></dd>

<dt>out_projection_step; <code>Integer</code>; 3d</dt>
<dd>
Interval time step of projection analysis 
if <code>projection_option</code> is not <code>'no'</code>.
Default is <code>100</code>.
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

<dt>out_psi; <code>Character</code>; 0d/3d</dt>
<dd>
If <code>'y'</code>, wavefunctions are output.
For periodic system (<code>iperiodic=3</code>), it works only for ground state calculation. The converged wave functions of all orbitals with all k-points are printed in gs_wfn_cube or gs_wfn_vtk directory. The format is speficied by <code>format3d</code>. 
Default is <code>'n'</code>.
</dd>

<dt>out_dos; <code>Character</code>; 0d/3d</dt>
<dd>
If <code>'y'</code>, density of state is output.
Default is <code>'n'</code>.
</dd>

<dt>out_dos_start; <code>Real(8)</code>; 0d/3d</dt>
<dd>
Start point (energy) of the density of state spectra.
If this value is lower than a specific value near the lowest energy level, 
this value is overwrited by that value. 
Default value is <code>-1.d10</code> eV.
</dd>

<dt>out_dos_end; <code>Real(8)</code>; 0d/3d</dt>
<dd>
End point (energy) of the density of state spectra.
If this value is higher than a specific value near the highest energy level, 
this value is overwrited by that value. 
Default value is <code>1.d10</code> eV.
</dd>

<dt>iout_dos_nenergy; <code>Integer</code>; 0d/3d</dt>
<dd>
Number of energies which are calculated in DOS part. 
Default is <code>601</code>.
</dd>

<dt>out_dos_smearing; <code>Real(8)</code>; 0d/3d</dt>
<dd>
Smearing width of the density of state spectra.
Default is <code>0.1</code> eV.
</dd>

<dt>out_dos_method; <code>Character</code>; 0d/3d</dt>
<dd>
Choise of smearing method for the density of state spectra.
<code>gaussian</code> and <code>lorentzian</code> function are available.
Default is <code>gaussian</code>.
</dd>

<dt>out_dos_fshift; <code>Character</code>; 0d/3d</dt>
<dd>
If <code>'y'</code>, the electron energy is shifted to fix the Fermi energy as zero point.
For <code>&system/iperiodic</code> is <code>0</code>, <code> out_dos_fshift</code> is not used 
if <code>&system/nstate</code> is equal to <code>&system/nelec</code>/2.
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

<dt>out_dns_rt/out_dns_rt_step; <code>Character/Integer</code>; 0d/3d</dt>
<dd>
If <code>'y'</code>, density during real-time time-propagation is output
every <code>outdns_rt_step</code> time steps.
Default is <code>'n'</code>.
</dd>

><dt>out_dns_trans/out_dns_trans_energy; <code>Character/Real(8)</code>; 3d</dt>
><dd>
>If <code>'y'</code>, transition in different density from the ground state at specified field frequency omega(given by <code>out_dns_trans_energy</code>) is calculated by drho(r,omega)=FT(rho(r,t)-rho_gs(r))/T.
>Default is <code>'n'/1.55eV</code>.
></dd>

<dt>out_elf; <code>Character</code>; 0d</dt>
<dd>
If <code>'y'</code>, electron localization function is output.
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

><dt>out_rvf_rt/out_rvf_rt_step; <code>Character/Integer</code>; 3d</dt>
><dd>
>If <code>'y'</code>, coordinates[A], velocities[au], forces[au] on atoms
>during real-time time-propagation are printed in <code>SYSname</code>_trj.xyz
>every <code>out_rvf_rt_step</code> time steps.
>If <code>use_ehrenfest_md='y'</code>, 
>the printing option is automatically turned on.
>Defaults are <code>'n'/10</code>.
></dd>

><dt>out_tm; <code>Character</code>; 3d</dt>
><dd>
>If <code>'y'</code>, trandition moments between occupied and virtual orbitals are printed into <code>SYSname</code>_tm.data after the ground state calculation.
>Defaults are <code>'n'/10</code>.
></dd>

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

><dt>timer_process; <code>Character</code>; 0d</dt>
><dd>
>Basically, elapsed times are written in the output file. 
>But if <code>timer_process</code> is <code>'y'</code>, 
>files of elapsed times for every process are also generated. 
>This variable is effective only for the real-time caululation.
>Default is <code>'n'</code>.
></dd>

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

>## &opt
><dl>
>
><dt>cg_alpha_up; <code>Real(8)</code>; 3d</dt>
><dd>
>Parameter for up-rate of step length in line search in conjugated gradient method.
>Default is <code>1.3</code>.
></dd>
>
><dt>cg_alpha_down; <code>Real(8)</code>; 3d</dt>
><dd>
>Parameter for down-rate of step length in line search in conjugated gradient method.
>Default is <code>0.5</code>.
></dd>
>
><dt>cg_alpha_ini; <code>Real(8)</code>; 3d</dt>
><dd>
>Parameter for initial step length in line search in conjugated gradient method. (currently not available)
>Default is <code>0.8</code>.
></dd>
>
><dt>convrg_scf_ene; <code>Real(8)</code>; 3d</dt>
><dd>
>Convergence threshold of SCF calculation in energy difference. If negative number no threshold (SCF loop is up to <code>Nscf</code>). The other SCF thresholds such as <code>threshold</code> in <code>&scf</code> are also applied (if you do not want to use it, set very small number). 
>Default is <code>-1.0</code>.
></dd>
>
><dt>convrg_scf_force; <code>Real(8)</code>; 3d</dt>
><dd>
>Convergence threshold of SCF calculation in force (average over atoms) difference. If negative number no threshold (SCF loop is up to <code>Nscf</code>). The other SCF thresholds such as <code>threshold</code> in <code>&scf</code> are also applied (if you do not want to use it, set very small number). 
>Default is <code>-1.0</code>.
></dd>
>
><dt>convrg_opt_fmax; <code>Real(8)</code>; 3d</dt>
><dd>
>Convergence threshold of optimization in maximum force.
>Default is <code>1d-3</code>.
></dd>
>
><dt>convrg_opt_ene; <code>Real(8)</code>; 3d</dt>
><dd>
>Convergence threshold of optimization in energy difference. (currently not available)
>Default is <code>1d-6</code>.
></dd>
>
></dl>


>## &md
><dl>
>
><dt>ensemble; <code>Character</code>; 3d</dt>
><dd>
>Ensemble in MD option: "NVE" or "NVT"
>Default is <code>"NVE"</code>.
></dd>
>
><dt>thermostat; <code>Character</code>; 3d</dt>
><dd>
>Thermostat in "NVT" option: (currently only <code>nose-hoover</code>)
>Default is <code>"nose-hoover"</code>.
></dd>
>
><dt>step_velocity_scaling; <code>Integer</code>; 3d</dt>
><dd>
>Time step interval for velocity-scaling. Velocity-scaling is applied if this is set to positive.
>Default is <code>-1</code>.
></dd>
>
><dt>step_update_ps/step_update_ps2; <code>Integer/Integer</code>; 3d</dt>
><dd>
>Time step interval for updating pseudopotential (Larger number makes calculation time reduce greatly, but gets inaccurate) in case of <code>use_ehrenfest_md=y</code>. <code>step_update_ps</code> is for full update and <code>step_update_ps2</code> is for update without changing grid points array.
>Default is <code>10/1</code>.
></dd>
>
><dt>temperature0_ion; <code>Real(8)</code>; 3d</dt>
><dd>
>Setting temperature in NVT ensemble and velocity scaling and for generating initial velocities.
>Default is <code>298.15</code>.
></dd>
>
><dt>set_ini_velocity; <code>Character</code>; 3d</dt>
><dd>
>Initial velocities are set.
>Default is <code>n</code>.
><ul>
><li>
><code>y</code>: Generate initial velocity with Maxwell-Bortzman distribution.
></li>
><li>
><code>r</code>: Read initial velocity from file specified by keyword of <code>file_ini_velocity</code>. This is, for example, used for restarting MD from the previous run. The last atomic coordinates and velocities are printed in <code>SYSname</code>_trj.xyz in that run. (atomic coordinate also should be copied from that and put in the next input file for restart)
></li>
></ul>
></dd>
>
><dt>file_ini_velocity; <code>Character</code>; 3d</dt>
><dd>
>File name for initial velocities. This is read when <code>set_ini_velocity</code> is <code>'r'</code>. The format is simply vx(iatom) vy(iatom) vz(iatom) in each line. The order of atoms must be the same as the given coordinates in the main input file. In case of using nose-hoover thermostat, a thermostat variable should be put at the last line (all atomic unit). 
>Default is <code>none</code>.
></dd>
>
><dt>file_set_shake; <code>Character</code>; 3d</dt>
><dd>
>Setting file for SHAKE method in ground-state MD is read. (now not supported yet)
>Default is <code>none</code>.
></dd>
>
><dt>thermostat_tau; <code>Real(8)</code>; 3d</dt>
><dd>
>Parameter in Nose-Hoover method: controlling time constant for temperature.
>Default is <code>41.34[au] or 1.0[fs]</code>.
></dd>
>
><dt>stop_system_momt; <code>Character</code>; 3d</dt>
><dd>
>Center of mass is stopped every time step.
>Default is <code>n</code>.
></dd>
>
></dl>

>***

>**Following variables are moved from the isolated part. Some of them may be added to common input, be combined to it, and be removed.**

>## &group_fundamental
><dl>
>
><dt>iditerybcg; <code>Integer</code>; 0d</dt>
><dd>
>Iterations for which ybcg is calculated if <code>&scf/amin_routine</code> is <code>'cg-diis'</code>.
>Default is <code>20</code>.
></dd>
>
><dt>iditer_nosubspace_diag; <code>Integer</code>; 0d</dt>
><dd>
>Iterations for which subspace diagonalization is not done if <code>&scf/subspace_diagonalization</code> is <code>'y'</code>.
>Default is <code>10</code>.
></dd>
>
><dt>ntmg; <code>Integer</code>; 0d</dt>
><dd>
>Number of multigrid calculation for gs. At the moment, there is a malfunction in this variable, and recovery is needed.
>Default is <code>1</code>.
></dd>
>
><dt>idisnum(2); <code>Integer</code>; 0d</dt>
><dd>
>Label numbers for two atoms which are measured the distance. 
>Default is <code>(/1,2/)</code>.
></dd>
>
><dt>iwrite_projection; <code>Integer</code>; 0d</dt>
><dd>
>A variable for projection. 
>Default is <code>0</code>.
></dd>
>
><dt>itwproj; <code>Integer</code>; 0d</dt>
><dd>
>The projection is calculated every <code>itwproj</code> time steps. 
>Default is <code>-1</code>.
></dd>
>
><dt>iwrite_projnum; <code>Integer</code>; 0d</dt>
><dd>
>There is a malfunction in this variable.
></dd>
>
><dt>itcalc_ene; <code>Integer</code>; 0d</dt>
><dd>
>Total energy is calculated every <code>itcalc_ene</code> time steps. There may be a malfunction in this variable.
>Default is <code>1</code>.
></dd>
>
></dl>


>## &group_parallel
><dl>
>
><dt>isequential; <code>Integer</code>; 0d</dt>
><dd>
>A variable to determine the way of assignment of processes.
>Default is <code>2</code>.
></dd>
>
><dt>imesh_s_all; <code>Integer</code>; 0d</dt>
><dd>
>A variable to determine how to use processes if total number of processes 
>and number of processes for Hartree/Exc calculation differ. 
>There may be a malfunction in this variable.
>Default is <code>1</code>.
></dd>
>
><dt>iflag_comm_rho; <code>Integer</code>; 0d</dt>
><dd>
>This variable may be removed. 
></dd>
>
></dl>


>## &group_hartree
><dl>
>
><dt>hconv; <code>Real(8)</code>; 0d</dt>
><dd>
>A convergence value for the Hartree-cg calculation. 
>The convergence is checked by ||tVh(i)-tVh(i-1)||<sup>2</sup>/(number of grids).
>Default is <code>1d-15</code> a.u. (= 1.10d-13 Å<sup>3</sup>eV<sup>2</sup>)
></dd>
>
><dt>lmax_meo; <code>Integer</code>; 0d</dt>
><dd>
>A maximum angular momentum for multipole expansion in the Hartree-cg calculation. 
>Default is <code>4</code>.
></dd>
>
></dl>


>## &group_file
><dl>
>
><dt>ic; <code>Integer</code>; 0d</dt>
><dd>
>A variable to check whether reentrance is done or not in the ground state calculation. 
>Default is <code>0</code> 
></dd>
>
><dt>oc; <code>Integer</code>; 0d</dt>
><dd>
>A variable to check whether intermediate files are generated in the ground state calculation. 
>Default is <code>1</code>.
></dd>
>
><dt>ic_rt; <code>Integer</code>; 0d</dt>
><dd>
>A variable to check whether reentrance is done or not in the time propagation calculation. 
>Default is <code>0</code> 
></dd>
>
><dt>oc_rt; <code>Integer</code>; 0d</dt>
><dd>
>A variable to check whether intermediate files are generated in the time propagation calculation. 
>Default is <code>0</code>.
></dd>
>
></dl>


>## &group_others
><dl>
>
><dt>iparaway_ob; <code>Integer</code>; 0d</dt>
><dd>
>A variable to determine the way of division for orbitals. 
><code>1</code> is block division, and <code>2</code> is cyclic division.
>Default is <code>2</code> 
></dd>
>
><dt>iscf_order; <code>Integer</code>; 0d</dt>
><dd>
>A variable to determine the order of the calculation for the ground state one. 
>Default is <code>1</code>.
></dd>
>
><dt>iswitch_orbital_mesh; <code>Integer</code>; 0d</dt>
><dd>
>A variable to apply descending order for orbitals in the ground state calculation.
>Default is <code>0</code> 
></dd>
>
><dt>iflag_psicube; <code>Integer</code>; 0d</dt>
><dd>
>A variable to generate cube files for wave functions. This variable will be removed.
></dd>
>
><dt>lambda1_diis/lambda2_diis; <code>Real(8)</code>; 0d</dt>
><dd>
>Parameters for the diis calculation.
>Default is <code>0.5/0.3</code>.
></dd>
>
><dt>file_ini; <code>Character</code>; 0d</dt>
><dd>
>A input file to align wavefunctions. 
>Default is <code>'file_ini'</code>.
></dd>
>
><dt>num_projection; <code>Interger</code>; 0d</dt>
><dd>
>Number of orbitals to write projections.
>Default is <code>1</code>.
></dd>
>
><dt>iwrite_projection_ob(200); <code>Interger</code>; 0d</dt>
><dd>
>Orbital number to be written as projections.
>Default is <code>(1/2/3/.../200)</code>.
></dd>
>
><dt>iwrite_projection_k(200); <code>Interger</code>; 0d</dt>
><dd>
>This variable will be removed.
></dd>
>
><dt>filename_pot; <code>Character</code>; 0d</dt>
><dd>
>Name of file to be written local potentials. 
>Default is <code>'pot'</code>.
></dd>
>
><dt>iwrite_external; <code>Integer</code>; 0d</dt>
><dd>
>A variable to generate file to be written local potentials. 
>Default is <code>0</code>.
></dd>
>
><dt>iflag_dip2; <code>Integer</code>; 0d</dt>
><dd>
>A variable to determine whether dipole moments are calculated in divided area. 
>Default is <code>0</code>.
></dd>
>
><dt>iflag_intelectron; <code>Integer</code>; 0d</dt>
><dd>
>A variable related to the quadrupole caluclation.
>Default is <code>0</code>.
></dd>
>
><dt>num_dip2; <code>Integer</code>; 0d</dt>
><dd>
>Number of area where dipole moments are calculated.
>Default is <code>1</code>.
></dd>
>
><dt>dip2boundary(100); <code>Real(8)</code>; 0d</dt>
><dd>
>Boundary position of area where dipole moments are calculated.
>Default is <code>0</code> a.u.
></dd>
>
><dt>dip2center(100); <code>Real(8)</code>; 0d</dt>
><dd>
>Origin in the dipole moment calculation. 
>Default is <code>0</code> a.u.
></dd>
>
><dt>iflag_fourier_omega; <code>integer</code>; 0d</dt>
><dd>
>A variable to determine whether Fourier transformation of 3d data for difference of density is calclated. 
>Default is <code>0</code>.
></dd>
>
><dt>num_fourier_omega; <code>Integer</code>; 0d</dt>
><dd>
>Number of energies for which the Fourier transformation is calclated. 
>Default is <code>1</code>.
></dd>
>
><dt>fourier_omega(200); <code>Real(8)</code>; 0d</dt>
><dd>
>Energies for which the Fourier transformation is calclated. 
>Default is <code>0</code> a.u.
></dd>
>
><dt>itotntime2; <code>Integer</code>; 0d</dt>
><dd>
>Number of time steps in the reentrance for real-time calculation.
>There may be a malfunction in this variable.
>Default is <code>0</code>.
></dd>
>
><dt>iwdenoption; <code>Integer</code>; 0d</dt>
><dd>
>A variable to determine whether 3d output is generated in real-time calculation. 
>This variable will be removed.
></dd>
>
><dt>iwdenstep; <code>Integer</code>; 0d</dt>
><dd>
>3d output is generated every <code>iwdenstep</code> time steps.
>This variable will be removed.
></dd>
>
><dt>iflag_estatic; <code>Integer</code>; 0d</dt>
><dd>
>A variable to determine whether 3d output for the static electric field is generated in real-time calculation. 
>This variable will be removed.
></dd>
>
></dl>
