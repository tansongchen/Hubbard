# Hubbard

This library gives 4 kinds of method to calculate the optical spectrum, i. e. the low-lying excited states that are reachable by optical excitations.

- `src/0_Hubbard.f90` defines the Hubbard model system
- `src/1_mean_field.f90` calculates the spectrum with Hartree-Fock method
- `src/2_casida.f90` calculates the spectrum with the Casida form of Bethe-Salpeter equation
- `src/3_TDBSE.f90` calculates the spectrum with the time-dependent Bethe-Salpeter equation (main work)
- `src/4_full_CI.f90` calculates the spectrum with full CI for comparison
