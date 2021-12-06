#model_glycan.py

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
        os.rmdir("./gly_gly")

    except OSError:
        pass

    os.mkdir("./gly_gly")

    #model glycans on NxT
    os.chdir("./gly_gly")
    for i in range(start,end+1):
        posi = str(i)+str(chain)
        sb = "/home/dwkulp/software/Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease -s ../gly_pp/mjg_N"+str(i)+"T.pdb -include_sugars -beta -write_pdb_link_records -maintain_links -parser:protocol /home/ywu/wistar/software/tools/Glycan_Scanner/core/glycans_tree_v2020_jw.xml -parser:script_vars positions="+posi+" glycosylation=man9"
        cmd = "/wistar/kulp/software/slurmq --sbatch \""+sb+"\""
        os.system(cmd)

    #model glycans on NxS
    for i in range(start,end+1):
        posi = str(i)+str(chain)
        sb = "/home/dwkulp/software/Rosetta/main/source/bin/rosetta_scripts.linuxgccrelease -s ../gly_pp/mjg_N"+str(i)+"S.pdb -include_sugars -beta -write_pdb_link_records -maintain_links -parser:protocol /home/ywu/wistar/software/tools/Glycan_Scanner/core/glycans_tree_v2020_jw.xml -parser:script_vars positions="+posi+" glycosylation=man9"
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
    end = chain_end_res(pose, c_num)
    res_num=pose.pdb_info().pdb2pose(chain,start)
    count = start
    while res_num<end:
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
