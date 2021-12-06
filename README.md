# KonrathWu_CellReports2021

For current release version, we write the code to run on our cluster only. A newer version of automatic code has been developed and will be avaiable soon. To run this code, currently, you need to install Rosetta and Slurmq and modify the path to your Rosetta Bin and your Slurmq and run the following codes sequentially:

1. python relax.py --pdb fhsr_mon.pdb

2. python glycan_auto.py --pdb gs_relax.pdb

3. python model_pngs.py --pdb gs_relax.pdb --chain A --start 59 --end 416 (position in chain)

4. python model_glycan.py --pdb gs_relax.pdb --chain Z --start 18 --end 359

5. python st_select.py --pdb gs_relax.pdb --chain Z --start 110 --end 214
