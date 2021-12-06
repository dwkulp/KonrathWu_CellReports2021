import optparse
import mutate_module
from pyrosetta import *
init()

def main(i,pdb,chain):

    pose = pose_from_pdb("../"+pdb)
    posi = pose.pdb_info().pdb2pose(chain,i)
    
    testPose = Pose()
    testPose.assign(pose)
    mutate_module.model(testPose,posi,'N')
    mutate_module.model(testPose,posi+2,'T')
    pdb_name = 'pp_'+'N'+str(i)+'T.pdb'
    testPose.dump_pdb(pdb_name)

    testPose.assign(pose)
    mutate_module.model(testPose,posi,'N')
    mutate_module.model(testPose,posi+2,'S')
    pdb_name = 'pp_'+'N'+str(i)+'S.pdb'
    testPose.dump_pdb(pdb_name)

#Set Up Option Parser
parser = optparse.OptionParser()
parser.add_option('--pdb', dest = 'pdb',
    default = 'nterm_m.0001.pdb',    # default example position
    help = 'the input struture' )
parser.add_option('--chain', dest = 'chain',
    default = 'A',    # default example position
    help = 'the mutated chain' )
parser.add_option('--posi', dest = 'i',
    default = '1',    # default example position
    help = 'the mutated position' )


(options,args) = parser.parse_args()

# PDB file option
i = int(options.i)
pdb = options.pdb
chain = options.chain

main(i,pdb,chain)
