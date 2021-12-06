#filter.py

import os
import optparse
from pyrosetta import *
init()

from pyrosetta.rosetta.core.pose import *

def main(pose,chain):

    start = c_start(pose,chain)
    end = c_end(pose,chain,start)

    #move lower energy pp pdbs
    os.system("rm -rf ./pp/pp_filter")
    os.mkdir("./pp/pp_filter")

    #move lower energy gly_pp pdbs
    os.system("rm -rf ./gly_gly/gg_filter")
    os.mkdir("./gly_gly/gg_filter")

    for i in range(start,end+1):
        sb = "/home/ywu/wistar/software/source/Python-3.7.3/bin/python3.7 /home/ywu/wistar/software/tools/Glycan_Scanner/core/unit_filter.py --pdb gs_relax.pdb --chain " + str(chain) + " --posi " + str(i)
        cmd = "/wistar/kulp/software/slurmq --sbatch \""+sb+"\""
        os.system(cmd)



def c_start(pose, chain):
    # find the start residue for the chain
    res_num=0
    count = 0
    while res_num==0:
        count += 1
        res_num = pose.pdb_info().pdb2pose(chain,count)
    return count

def c_end(pose, chain, start):
    # find the end residue for the chain
    c_num = get_chain_id_from_chain(chain,pose)
    c_end = chain_end_res(pose, c_num)
    res_num=pose.pdb_info().pdb2pose(chain,start)
    count = start
    while res_num<c_end:
        count += 1
        res_num = pose.pdb_info().pdb2pose(chain,count)
    return count

#Set Up Option Parser
parser = optparse.OptionParser()
parser.add_option('--pdb', dest = 'pdb',
    default = 'gs_relax.pdb',    # default example position
    help = 'the input struture' )
parser.add_option('--chain', dest = 'chain',
    default = 'A',    # default example position
    help = 'the mutated chain' )

(options,args) = parser.parse_args()

# PDB file option
pdb = options.pdb
chain = options.chain

pose = Pose()
# -load the data from pdb_file into the pose
pose_from_file(pose, pdb)

main(pose, chain)
