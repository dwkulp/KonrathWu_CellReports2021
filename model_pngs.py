import os
import optparse
from pyrosetta import *
init()

from pyrosetta.rosetta.core.pose import *

def main(pose,chain,start,end):

    if(start==0 and end ==0):

        start = c_start(pose,chain)
        end = c_end(pose,chain,start)

    try:
        os.rmdir("./pp")
        os.rmdir("./gly_pp")

    except OSError:
        pass

    os.mkdir("./pp")
    os.mkdir("./gly_pp")

    #protein pngs
    os.chdir("./pp")
    for i in range(start,end+1):
        sb = "/home/ywu/wistar/software/source/Python-3.7.3/bin/python3.7 /home/ywu/wistar/software/tools/Glycan_Scanner/core/pngs.py --pdb gs_relax.pdb --chain "+str(chain)+" --posi "+str(i)
        cmd = "/wistar/kulp/software/slurmq --sbatch \""+sb+"\""
        os.system(cmd)
    os.chdir("..")

    #glycan pngs
    os.chdir("./gly_pp")
    for i in range(start,end+1):
        sb = "/home/ywu/wistar/software/source/Python-3.7.3/bin/python3.7 /home/ywu/wistar/software/tools/Glycan_Scanner/core/glycan_pngs.py --pdb gs_glycan.pdb --chain "+str(chain)+" --posi "+str(i)
        cmd = "/wistar/kulp/software/slurmq --sbatch \""+sb+"\""
        os.system(cmd)
    os.chdir("..")


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
parser.add_option('--start', dest = 'start',
    default = '0',    # default example position
    help = 'the input start position' )
parser.add_option('--end', dest = 'end',
    default = '0',    # default example position
    help = 'the input end position' )

(options,args) = parser.parse_args()

# PDB file option
pdb = options.pdb
chain = options.chain
start = int(options.start)
end = int(options.end)
pose = Pose()
# -load the data from pdb_file into the pose
pose_from_file(pose, pdb)

main(pose,chain,start,end)
