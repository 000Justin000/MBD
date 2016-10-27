Overview
-----------------------------------------------------------
The many-body dispersion (MBD) energy is computed using the coupled 
Quantum Harmonic Oscillator (QHO) model Hamiltonian, that maps the 
valence-electron response of finite-gap systems onto a set of atom-centered 
three-dimensional QHOs. In this version the MBD method is developed to
be coupled with semi-local DFT functionals by using the so-called
range-separated self-consistent screening (rsSCS) procedure. 
The short-range part of the dipole interaction is used to self-consistently 
screen the atomic polarizablities. Such short-ranged screened polarizabilities 
are consequently used as input in the MBD energy expression.

The MBD energy contains two crucial contributions that are missing in 
simple pairwise dispersion approaches: 
(1) many-body energy (Axilrod-Teller and higher-order),
(2) long-range Coulomb response (screening) that effectively 
modifies the polarizabilities of interacting species.

(This code is modified on top of the effort by Alexandre Tkatchenko)

Download
-----------------------------------------------------------
https://github.com/000Justin000/MBD

Contact
-----------------------------------------------------------
Junteng Jia (jj585@cornell.edu)


Compilation
-----------------------------------------------------------

To build the MBD code, modify the Makefile and change the FC,
LD and MKL_ROOT flag to to specify fortran complier and 
(SCA)LAPACK/BLAS library. 

example :- 

FC=mpiifort
MKL_ROOT=/opt/intel/composer_xe_2013.2.146/mkl/lib/intel64

After this, you should be able to generate the `pmbd.x` 
binary by running:

make


Use
-----------------------------------------------------------

In order to use MBD, you need a geometry file similar to XYZ format.
The fifth column in the XYZ file specifies the relative Hirshfeld
volume (V_eff / V_free) corresponding to a given atom in the molecule
or solid. See Eq.(8) in Ref. [1] below.


example:-
-------------------------geometry.in-----------------------------
12

C -1.04163594 -1.42199238  0.00000000 0.8062
C -1.45013466 -0.85487018  1.21057904 0.8062
C -1.45013466 -0.85487018 -1.21057904 0.8062
C -2.26713561  0.27938790  1.21059830 0.8062
C -2.67555716  0.84658359  0.00000000 0.8062
C -2.26713561  0.27938790 -1.21059830 0.8062
H -1.13113521 -1.29793512 -2.15594541 0.6017
H -2.58624335  0.72214218 -2.15607894 0.6017
H -3.31360362  1.73238019  0.00000000 0.6016
H -2.58624335  0.72214218  2.15607894 0.6017
H -1.13113521 -1.29793512  2.15594541 0.6017
H -0.40357262 -2.30778356  0.00000000 0.6016
lattice_vector   20.0       0.0      0.0
lattice_vector    0.0      20.0      0.0
lattice_vector    0.0       0.0     20.0
-----------------------------------------------------------------
if string "lattice_vector" is present in "geometry.xyz" file then
the code uses periodic boundary conditions.


and a "setting.in" file which contains the parameters for the MBD code

example:-
--------------------------setting.in-----------------------------
xc                    1
mbd_cfdm_dip_cutoff   100.d0
mbd_supercell_cutoff  25.d0
mbd_scs_dip_cutoff    120.0
mbd_scs_vacuum_axis   .false. .false. .false.
nprow                 2
npcol                 3
-----------------------------------------------------------------

Keyword information:--------------------

xc (value)
Specify Exchange-correction type 
xc   1   #(for PBE type)
xc   2   #(for PBE0 type)
xc   3   #(for HSE type)

mbd_cfdm_dip_cutoff (value)
Radius used to integrate dipole field in periodic MBD calculation (internal default=100.0 Angstrom)

mbd_scs_dip_cutoff (value)
Radius used to integrate dipole field in periodic SCS	 calculation (internal default=100.0 Angstrom)
NOTE: Our numerical tests suggest that C6 coefficient should be converged by tuning the mbd_scs_dip_cutoff parameter.
This is especially important when studying low-dimensional systems. For example for graphene the C6 coefficient
is of the order 140 hartree·bohr6 which requires a 700 Angstrom cutoff.

mbd_supercell_cutoff (value)
Radius used to construct the supercell in periodic MBD calcula-
tions. (default = 25 Angstrom) NOTE: The convergence wrt the value of
mbd_supercell_cutoff must be carefully tested for periodic systems with low
symmetry.

mbd_scs_vacuum_axis (flag) (flag) (flag)
This keyword specifies directions to be treated as vacuum in the case of low dimensional systems.
Default: No vacuum, i.e. mbd_scs_vacuum_axis .false. .false. .false.
For example in the case of periodic slab along XY directions the keyword
mbd_scs_vacuum_axis .false. .false. .true. is needed to specify Z as
the vacuum direction in MBD/SCS calculations.

nprow (value)
number of processes along the row direction for the BLACS grid

npcol (value)
number of processes along the column direction for the BLACS grid

(nprow * npcol shoule equal the total number of processes used to run the program)

You can run following command to compute the MBD energy: 
mpiexec -np 6 `pwd`/pmbd.x geometry.in setting.in


References and citations
-----------------------------------------------------------

The basic references for the MBD method are:

[1] Accurate Molecular Van Der Waals Interactions from Ground-State Electron Density and Free-Atom Reference Data.
    A. Tkatchenko and M. Scheffler,
    Phys. Rev. Lett. 102, 073005 (2009).

[2] Accurate and Efficient Method for Many-Body van der Waals Interactions.
   A. Tkatchenko R. A. DiStasio Jr., R. Car, and M. Scheffler
   Phys. Rev. Lett. 108, 236402 (2012)

[3] Long-range correlation energy calculated from coupled atomic response functions.
   A. Ambrosetti, A. M. Reilly, R. A. DiStasio Jr., and A. Tkatchenko
   J. Chem. Phys., to be published. 

A concise review of the MBD methodology is available at:
http://th.fhi-berlin.mpg.de/site/uploads/Publications/psik-ver1-3_20121212.pdf,
to be published in J. Phys.: Condens. Matter.
