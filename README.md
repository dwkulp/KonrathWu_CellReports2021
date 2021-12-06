# KonrathWu_CellReports2021

python relax.py --pdb fhsr_mon.pdb

python glycan_auto.py --pdb gs_relax.pdb

python model_pngs.py --pdb gs_relax.pdb --chain A --start 59 --end 416 (position in chain)

python model_glycan.py --pdb gs_relax.pdb --chain Z --start 18 --end 359

python st_select.py --pdb gs_relax.pdb --chain Z --start 110 --end 214
