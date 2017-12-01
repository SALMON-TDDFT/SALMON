We here summarize namelists that appear in this Tutorial. A thorough list of the namelist variables may be found in the downloaded file in 'SALMON/manual/input_variables.md'.

# &units

Mandatory: none

```
&units
  unit_system='A_eV_fs'
/
```

This namelist specifies the unit system to be used in the input file. Options are 'A_eV_fs' for Angstrom, eV, and fs, and 'a.u.' or 'au' for atomic units. If you do not specify it, atomic unit will be used as default.

For isolated systems (specified by `iperiodic = 0` in `&system`), the unit of 1/eV is used for the output files of DOS and PDOS if `unit_system = 'A_eV_fs'` is specified, while atomic unit is used if not. For other output files, the Angstrom/eV/fs units are used irrespective of the namelist value.

For periodic systems (specified by `iperiodic =3` in `&system`), the unit system specified by this namelist variable is used for most output files. See the first few lines of output files to confirm the unit system adopted in the file.

# &calculation

Mandatory: calc_mode

```
&calculation
  calc_mode = 'GS'
/
```

The value of the `calc_mode` should be one of `'GS'`, `'RT'`, and `'GS-RT'`. For isolated systems (specified by `iperiodic = 3` in ``&system`), the ground state (`'GS'`) and the real time (`'RT'`) calculations should be done separately and sequentially. For periodic systems (specified by `iperiodic = 3` in`` &system`), both ground state and real time calculations should be carried out as a single task (`calc_mode = 'GS_RT'`).

For Maxwell + TDDFT multi-scale calculation, add the following namelist.

`use_ms_maxwell = 'y'`

# &control

Mandatory: none

```
&control
  sysname = 'C2H2'
/
```

'C2H2' defined by `sysname = 'C2H2'` will be used in the filenames of output files. If you do not specify it, the file name will start with 'default'.

# &functional

```
&functional
  xc ='PZ'
/
```

`xc ='PZ'` indicates that (adiabatic) local density approximation is adopted (Perdew-Zunger: Phys. Rev. B23, 5048 (1981)). This is the default choice.

For isolated systems (specified by `iperiodic = 0` in `&system`), only the default choice of 'PZ' is available at present.

For periodic systems (specified by `iperiodic = 3` in `&system`), the following functionals may be available in addition to 'PZ':

`xc = 'PZM'`

Perdew-Zunger LDA with modification to improve sooth connection between high density form and low density one. :J. P. Perdew and Alex Zunger, Phys. Rev. B 23, 5048 (1981).

`xc = 'TBmBJ' cval = 1.0`

Tran-Blaha meta-GGA exchange with Perdew-Wang correlation. :Fabien Tran and Peter Blaha, Phys. Rev. Lett. 102, 226401 (2009). John P. Perdew and Yue Wang, Phys. Rev. B 45, 13244 (1992). This potential is known to provide a reasonable description for the bandage of various insulators. For this choice, the additional mixing parameter 'cval' may be specified. If cval is set to a minus value, the mixing-parameter will be computed following the formula in the original paper [Phys. Rev. Lett. 102, 226401 (2009)]. The default value for this parameter is 1.0.

# &system

Mandatory: iperiodic, al, nstate, nelem, natom

**For an isolated molecule (Tutorial-1, 2, 3)**:

```
&system
  iperiodic = 0
  al = 16d0, 16d0, 16d0
  nstate = 5
  nelem = 2
  natom = 4
  nelec = 10
/
```

`iperiodic = 0` indicates that the isolated boundary condition will be used in the calculation. `al = 16d0, 16d0, 16d0` specifies the lengths of three sides of the rectangular parallelepiped where the grid points are prepared. `nstate = 8` indicates the number of Kohn-Sham orbitals to be solved. `nelec = 10` indicate the number of valence electrons in the system. Since the present code assumes that the system is spin saturated, `nstate` should be equal to or larger than `nelec/2`. `nelem = 2` and `natom = 4` indicate the number of elements and the number of atoms in the system, respectively.

**For a periodic system (Tutorial-4, 5)**:

```
&system
  iperiodic = 3
  al = 10.26d0,10.26d0,10.26d0
  nstate = 32
  nelec = 32
  nelem = 1
  natom = 8
/
```

`iperiodic = 3` indicates that three dimensional periodic boundary condition (bulk crystal) is assumed. `al = 10.26d0, 10.26d0, 10.26d0` specifies the lattice constans of the unit cell. `nstate = 32` indicates the number of Kohn-Sham orbitals to be solved. `nelec = 32` indicate the number of valence electrons in the system. `nelem = 1` and `natom = 8` indicate the number of elements and the number of atoms in the system, respectively.

**For Maxwell - TDDFT multi scale calculation (Tutorial-6)**:

```
&system
  iperiodic = 3
  al = 10.26d0,10.26d0,10.26d0
  isym = 8
  crystal_structure = 'diamond'
  nstate = 32
  nelec = 32
  nelem = 1
  natom = 8
/
```

The difference from the above case is the variables, `isym = 8` and `crystal_structure = 'diamond'`, which indicates that the spatial symmetry of the unit cell is used in the calculation. Although the use of the symmetry substantially reduces the computational cost, it should be used very carefully. At present, the spatial symmetry has been implemented only for the case of the diamond structure.

# &pseudo

Mandatory: pseudo_file, izatom

**For C2H2 molecule**:

```
&pseudo
  izatom(1)=6
  izatom(2)=1
  pseudo_file(1)='C_rps.dat'
  pseudo_file(2)='H_rps.dat'
  lmax_ps(1)=1
  lmax_ps(2)=0
  lloc_ps(1)=1
  lloc_ps(2)=0
/
```

Parameters related to atomic species and pseudopotentials. `izatom(1) = 6` specifies the atomic number of the element #1\. `pseudo_file(1) = 'C_rps.dat'` indicates the filename of the pseudopotential of element #1\. `lmax_ps(1) = 1` and `lloc_ps(1) = 1` specify the maximum angular momentum of the pseudopotential projector and the angular momentum of the pseudopotential that will be treated as local, respectively.

**For crystalline Si**:

```
&pseudo
  izatom(1)=14
  pseudo_file(1) = './Si_rps.dat'
  lloc_ps(1)=2
/
```

`izatom(1) = 14` indicates the atomic number of the element #1\. `pseudo_file(1) = 'Si_rps.dat'` indicates the pseudopotential filename of element #1\. `lloc_ps(1) = 2` indicate the angular momentum of the pseudopotential that will be treated as local.

# &atomic_coor

Mandatory: atomic_coor or atomic_red_coor (they may be provided as a separate file)

**For C2H2 molecule**:

```
&atomic_coor
  'C' 0.000000 0.000000 0.599672 1
  'H' 0.000000 0.000000 1.662257 2
  'C' 0.000000 0.000000 -0.599672 1
  'H' 0.000000 0.000000 -1.662257 2
/
```

Cartesian coordinates of atoms. The first column indicates the element. Next three columns specify Cartesian coordinates of the atoms. The number in the last column labels the element.

# &atomic_red_coor

Mandatory: atomic_coor or atomic_red_coor (they may be provided as a separate file)

**For a crystalline silicon**:

```
&atomic_red_coor
  'Si' .0 .0 .0 1
  'Si' .25 .25 .25 1
  'Si' .5 .0 .5 1
  'Si' .0 .5 .5 1
  'Si' .5 .5 .0 1
  'Si' .75 .25 .75 1
  'Si' .25 .75 .75 1
  'Si' .75 .75 .25 1
/
```

Cartesian coordinates of atoms are specified in a reduced coordinate system. First column indicates the element, next three columns specify reduced Cartesian coordinates of the atoms, and the last column labels the element.

# &rgrid

Mandatory: dl or num_rgrid

This namelist provides grid spacing of Cartesian coordinate system. `dl(3)` specify the grid spacing in three Cartesian coordinates. This is adopted for C2H2 calculation (Tutorial-1).

```
&rgrid
dl = 0.25d0, 0.25d0, 0.25d0
/
```

`num_rgrid(3)` specify the number of grid points in each Cartesian direction. This is adopted for crystalline Is calculation (Tutorial-4, 5, 6).

```
&rgrid
  num_rgrid = 12,12,12
/
```

# &kgrid

Mandatory: none

This namelist provides grid spacing of k-space for periodic systems.

```
&kgrid
num_kgrid = 4,4,4
/
```

# &scf

Mandatory: nscf

This namelists specify parameters related to the self-consistent field calculation.

```
&scf
  ncg = 4
  nscf = 1000
  convergence = 'norm_rho_dng'
  threshold_norm_rho = 1.d-15
/
```

`ncg = 4` is the number of conjugate-gradient iterations in solving the Kohn-Sham equation. Usually this value should be 4 or 5\. `nscf = 1000` is the number of scf iterations. For isolated systems specified by `&system/iperiodic = 0`, the scf loop in the ground state calculation ends before the number of the scf iterations reaches `nscf`, if a convergence criterion is satisfied. There are several options to examine the convergence. If the value of `norm_rho_dng` is specified, the convergence is examined by the squared difference of the electron density,

# &hartree

Mandatory: none

```
&hartree
  meo = 3
  num_pole_xyz = 2,2,2
/
```

`meo` specifies the order of multipole expansion of electron density that is used to prepare boundary condition for the Hartree potential.

- `meo=1`: Monopole expansion (spherical boundary condition).
- `meo=2`: Multipole expansions around each atom.
- `meo=3`: Multipole expansion around the center of mass of electrons in cubits that are defined by `num_pole_xyz`.

`num_pole_xyz(3)` defines the division of space when `meo = 3` is specified.

A default for `meo` is `3`, and defaults for `num_pole_xyz` are `(0,0,0)`. When default is set for `num_pole_xyz`, the division of space is carried out using a prescribed method.

# &tgrid

Mandatory: dt, Nt

```
&tgrid
  dt=1.25d-3
  nt=5000
/
```

`dt=1.25d-3` specifies the time step of the time evolution calculation. `nt=5000` specifies the number of time steps in the calculation.

# &propagation

This namelist specifies the numerical method for time evolution calculations of electron orbitals.

```
&propagation
  propagator='etrs'
/
```

`propagator = 'etrs'` indicates the use of enforced time-reversal symmetry propagator. [M.A.L. Marques, A. Castro, G.F. Bertsch, and A. Rubio, Comput. Phys. Commun., 151 60 (2003)](https://doi.org/10.1016/S0010-4655(02)00686-0).

```
&propagation
  propagator='middlepoint'
/
```

`propagation='middlepoint'` indicates that Hamiltonian at midpoint of two-times is used.

The default is _middlepoint_.

# &emfield

This namelist specifies the pulse shape of an electric filed applied to the system in time evolution calculations. We explain below separating two cases, [#Linear response calculations](#Linear_response_calculations "wikilink") and [#Pulsed electric field calculations](#Pulsed_electric_field_calculations "wikilink").

## Linear response calculations

A weak impulsive field is applied at _t=0_. For this case, `ae_shape1 = 'impulse'` should be described.

Mandatory: ae_shape1

```
&emfield
  ae_shape1 = 'impulse'
  epdir_re1 = 0.d0,0.d0,1.d0
/
```

`epdir_rel(3)` specify a unit vector that indicates the direction of the impulse.

For a periodic system specified by `iperiodic = 3`, one may add `trans_longi`. It has the value, `'tr'`(transverse) or `'lo'`(longitudinal), that specifies the treatment of the polarization in the time evolution calculation. The default is `'tr'`.

```
&emfield
  trans_longi = 'tr'
  ae_shape1 = 'impulse'
  epdir_re1 = 0.,0.,1.
/
```

The magnitude of the impulse of the pulse may be explicitly specified by, for example, `e_impulse = 1d-2`. The default is '1d-2' in atomic unit.

## Pulsed electric field calculations

A Pulsed electric field of finite time duration is applied. For this case, `as_shape1` should be specified. It indicates the shape of the envelope of the pulse. The options include 'Acos2' and 'Ecos2' (See below for other options).

Mandatory: ae_shape1, epdir_re1, {rlaser_int1 or amplitude1}, omega1, pulse_tw1, phi_cep1

```
&emfield
  ae_shape1 = 'Ecos2'
  epdir_re1 = 0.d0,0.d0,1.d0
  rlaser_int_wcm2_1 = 1.d8
  omega1=9.28d0
  pulse_tw1=6.d0
  phi_cep1=0.75d0
/
```

`ae_shape1 = 'Ecos2'` specifies the envelope of the pulsed electric field, 'Ecos2' for the cos\^2 envelope for the electric field. If 'Acos2' is specified, this gives cos\^2 envelope for the vector potential. Note that 'phi_cep1' must be 0.75 (or 0.25) if one employs 'Ecos2' pulse shape, since otherwise the time integral of the electric field does not vanish. There is no such restriction for the 'Acos2' pulse shape.

`epdir_re1 = 0.d0,0.d0,1.d0` specifies the real part of the unit polarization vector of the pulsed electric field. If only the real part is specified, it describes a linearly polarized pulse. Using both real ('epdir_re1') and imaginary ('epdir_im1') parts of the polarization vector, circularly (and general ellipsoidary) polarized pulses may be described.

`laser_int_wcm2_1 = 1.d8` specifies the maximum intensity of the applied electric field in unit of W/cm\^2\. It is also possible to specify the maximum intensity of the pulse by `amplitude1`.

`omega1=9.26d0` specifies the average photon energy (frequency multiplied with hbar).

`pulse_tw1=6.d0` specifies the pulse duration. Note that it is not the FWHM but a full duration of the cos\^2 envelope.

`phi_cep1=0.75d0` specifies the carrier envelope phase of the pulse. As noted above, 'phi_cep1' must be 0.75 (or 0.25) if one employs 'Ecos2' pulse shape, since otherwise the time integral of the electric field does not vanish. There is no such restriction for the 'Acos2' pulse shape.

It is possible to use two pulses simultaneously to simulate pump-probe experiments, adding information for two pulses. To specify the second pulse, change from 1 to 2 in the namelist variables, like `ae_shape2`. The time delay between two pulses is specified by the variable 't1_t2'.

For a periodic system specified by `iperiodic = 3`, one may add `trans_longi`. It has the value, `'tr'`(transverse) or `'lo'`(longitudinal), that specifies the treatment of the polarization in the time evolution calculation. The default is `'tr'`. For a periodic system, it is also specify 'Acos3', 'Acos4', 'Acos6', 'Acos8' for `ae_shape1`.

# &analysis

Mandatory: none

The following namelists specify whether the output files are created or not after the calculation. In the ground state calculation of isolated systems (Tutorial-1):

```
&analysis
  out_psi = 'y'
  out_dns = 'y'
  out_dos = 'y'
  out_pdos = 'y'
  out_elf = 'y'
/
```

In the time evolution calculation of isolated systems (Tutorial-3):

```
&analysis
  out_dns_rt = 'y'
  out_elf_rt = 'y'
  out_estatic_rt = 'y'
/
```

In the following namelists, variables related to time-frequency Fourier analysis are specified.

```
&analysis
  nenergy=1000
  de=0.001
/
```

`nenergy=1000` specifies the number of energy steps, and `de=0.001` specifies the energy spacing in the time-frequency Fourier transformation.

# &multiscale

This namelist specifies information necessary for Maxwell - TDDFT multiscale calculations.

```
&multiscale
  fdtddim = '1D'
  twod_shape = 'periodic'
  nx_m = 4
  ny_m = 1
  hX_m = 250d0
  nxvacl_m = -2000
  nxvacr_m = 256
/
```

`fdtddim` specifies the spatial dimension of the macro system. `fdtddim='1D'` indicates that one-dimensional equation is solved for the macroscopic vector potential.

`nx_m = 4` specifies the number of the macroscopic grid points in for x-direction in the spatial region where the material exists.

`hx_m = 250d0` specifies the grid spacing of the macroscopic grid in x-direction.

`nxvacl_m = -2000` and `nxvacr_m = 256` indicate the number of grid points in the vacuum region, `nxvacl_m` for the left and `nxvacr_m` for the right from the surface of the material.

# &parallel

When you execute a job with MPI parallelization, you are not required to specify any parameters that describe the assignment of the parallelization; the assignment is carried out automatically. You may also specify the parameters explicitly as below.

Mandatory: none

```
&parallel
  nproc_ob = 1
  nproc_domain = 1,1,1
  nproc_domain_s = 1,1,1
/ 
```

- `nproc_ob` specifies the number of MPI parallelization to divide the electron orbitals. The default value is 0 (automatic parallelization).
- `nproc_domain(3)`specifies the number of MPI parallelization to divide the spatial grids of the electron orbitals in three Cartesian directions. The default values are (0/0/0) (automatic parallelization).
- `nproc_domain_s(3)'` specifies the number of MPI parallelization to divide the spatial grids related to the electron density in three Cartesian directions. The default values are (0/0/0) (automatic parallelization).

The following conditions must be satisfied.

- The total number of processors must be equal to both `nproc_ob * nproc_domain(1) * nproc_domain(2) * nproc_domain(3)` and also `nproc_domain_s(1) * nproc_domain_s(2) * nproc_domain_s(3)`.
- `nproc_domain_s(1)` is a multiple of `nproc_domain(1)`, and the same relations to the second and third components.
