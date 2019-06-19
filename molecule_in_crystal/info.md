For molecular crystals, it is possible to calculate the polarizability tensors of its molecular units 
with a KRR approach, and use them as a baseline for the calculation of the full crystal polarizability.

:older_man: _The whole is different than the sum of its parts._

In this directory you have the necessary files to do so. The steps to follow are:

1) Break down the crystal trajectory into different trajectories for every molecule in the primitive cell:
(you need to adapt it to your system :no_mouth:): `split_crystal_into_mol.py`.
2) Align these individual trajectories to a reference frame (the same you used to train the molecule on), 
and store the corresponding rotation matrices: `align-molecules_rot.py`.
3) Predict the polarizability tensors for every molecule with `KRR_polar.py`.
4) Rotate these molecular polarizabilities back to their original orientation within the cell:
`transform_polar.py`
5) Sum up all molecular polarizability tensors: `sum_molecular_polar.py`.
You can then use this quantity as a baseline (to specify in the control file) for the full crystal.


