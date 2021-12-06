#native_gly.py

import optparse
import os
from pyrosetta import *
init()

def main(pose):

    pngs_list = ""
    pngs_array = []
    for i in range(1, pose.total_residue()+1):
        if(is_sequeon(i) and is_sequeon(i+1)):
            pass

        if(is_sequeon(i) and not is_sequeon(i+1)):
            pose_to_pdb = pose.pdb_info().pose2pdb(i)
            sp = pose_to_pdb.split(" ")
            chain_n = str(sp[0])+str(sp[1])
            pngs_array.append(chain_n)

    pngs_list = str(pngs_array[0])
    for i in range(1,len(pngs_array)):
        pngs_list += ","+pngs_array[i]

    sb = "/home/dwkulp/software/Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease -s "+str(pdb)+" -include_sugars -beta -write_pdb_link_records -maintain_links -parser:protocol /home/ywu/wistar/software/tools/Glycan_Scanner/core/glycans_tree_v2020_jw.xml -parser:script_vars positions="+str(pngs_list)+" glycosylation=man9"
    cmd = "/wistar/kulp/software/slurmq --sbatch \""+sb+"\""
    
    #print(cmd)
    os.system(cmd)


def is_sequeon(posi):
    if(pose.residue(posi).name()=='ASN') and (pose.residue(posi+1).name()!='PRO') and (pose.pdb_info().chain(posi)==pose.pdb_info().chain(posi+2)):
        if(pose.residue(posi+2).name()=='SER' or pose.residue(posi+2).name()=='THR'):
            return True
    return False

#Set Up Option Parser
parser = optparse.OptionParser()
parser.add_option('--pdb', dest = 'pdb',
    default = 'gs_relax.pdb',    # default example position
    help = 'the input struture' )
(options,args) = parser.parse_args()
# PDB file option
pdb = options.pdb

pose = Pose()
pose = pose_from_pdb(pdb)

main(pose)
