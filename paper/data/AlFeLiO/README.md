# Al-Fe-Li-O synthetic XRD dataset

* `alfeli.npy`

Contains all 231 patterns, each has a length of 650. The q range is 15-79.9 nm⁻¹ with 0.1 nm⁻¹ increment.

* `sol.csv`

Contains the activation of the ground truth phases for each of the 231 patterns.
There are 6 phases, list in order: `Al2O3_R-3cH`, `Li2O_Fm-3m`, `Fe2O3_R-3cH`, `LiAl5O8_P4332`, `LiAlO2_R-3mH`, `LiFeO2_R-3mH`

* `sticks.csv`

Input file for creating `CrystalPhase` object.
