#st_selection.py
import os
import optparse
from pyrosetta import *
init()

def main(posi):

    pdb_name_t = "./pp/pp_N"+str(posi)+"T.pdb"
    pdb_name_s = "./pp/pp_N"+str(posi)+"S.pdb"
    gly_name_t = "./gly_gly/mjg_N"+str(posi)+"T_0001.pdb"
    gly_name_s = "./gly_gly/mjg_N"+str(posi)+"S_0001.pdb"
    pose_t = pose_from_pdb(pdb_name_t)
    pose_s = pose_from_pdb(pdb_name_s)
    
    scorefxn = get_fa_scorefxn()

    #move lower energy pp pdbs
    if(scorefxn(pose_t)<=scorefxn(pose_s)):
        cmd="cp "+str(pdb_name_t)+" ./pp/pp_select"
        os.system(cmd)

    else:
        cmd="cp "+str(pdb_name_s)+" ./pp/pp_select"
        os.system(cmd)

    #move lower energy gly_pp pdbs
    if(scorefxn(pose_t)<=scorefxn(pose_s)):
        cmd="cp "+str(gly_name_t)+" ./gly_gly/gg_select"
        os.system(cmd)

    else:
        cmd="cp "+str(gly_name_s)+" ./gly_gly/gg_select"
        os.system(cmd)


#Set Up Option Parser
parser = optparse.OptionParser()
parser.add_option('--posi', dest = 'posi',
    default = '1',    # default example position
    help = 'position for sorting' )

(options,args) = parser.parse_args()
# PDB file option
posi = options.posi
main(posi)
