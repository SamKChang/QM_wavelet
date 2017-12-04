
Data set GDB10-GDB17 with ~10k for each subset
==============================================

Thermochemical properties for ~70k small organic molecules at the DFT/B3LYP(2df,p) level of theory.

Molecules
---------

There are in total 8 datasests, each consists of ~10k molecules randomly chosen from the complete 
GDB_n dataset.

Please be aware that some molecule Ids are missing due to either it's hard to converge or dissociated
after optimization.

Format
------

Each molecule is stored in its own file, ending in ".xyz".
The format is an ad hoc extension of the XYZ format [3].

Line       Content
----       -------
1          Number of atoms na
2          Properties 1-17 (see below)
3,...,na+2 Element type, coordinate (x,y,z) (Angstrom), and Mulliken partial charge (e) of atom

The properties stored in the second line of each file:

I.  Property  Unit         Description
--  --------  -----------  --------------
 1  tag       -            "gdb9"; string constant to ease extraction via grep
 2  index     -            Consecutive, 1-based integer identifier of molecule
 3  A         GHz          Rotational constant A
 4  B         GHz          Rotational constant B
 5  C         GHz          Rotational constant C
 6  mu        Debye        Dipole moment
 7  alpha     Bohr^3       Isotropic polarizability
 8  homo      Hartree      Energy of Highest occupied molecular orbital (HOMO)
 9  lumo      Hartree      Energy of Lowest occupied molecular orbital (LUMO)
10  gap       Hartree      Gap, difference between LUMO and HOMO
11  r2        Bohr^2       Electronic spatial extent
12  zpve      Hartree      Zero point vibrational energy
13  U0        Hartree      Internal energy at 0 K
14  U         Hartree      Internal energy at 298.15 K
15  H         Hartree      Enthalpy at 298.15 K
16  G         Hartree      Free energy at 298.15 K
17  Cv        cal/(mol K)  Heat capacity at 298.15 K

I. = Property index (properties are given in this order)
For the 6095 isomers, properties 12-16 were calculated at the G4MP2 level of theory.
All other calculations were done at the DFT/B3LYP/6-31G(2df,p) level of theory.


