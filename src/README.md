This is a repository for testing time-dependent calculation of electron-hole excitation with 1D/2D Hubbard model.

# File structure

- File `Hubbard.f90` stores basic parameters and properties of the system;
- File `1_mean_field.f90` carries out Hartree-Fock mean field calculation to obtain the 1-particle orbitals, and saves them in file `HFOrbitals.dat`;
- File `2_casida.f90` carries out static (and later dynamic) W Casida equation calculation, and saves them in file `Static.dat`.