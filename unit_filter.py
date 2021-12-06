#unit_filter.py

import os
import optparse
from pyrosetta import *
init()

from pyrosetta.rosetta.core.pose import *
from pyrosetta.rosetta.core.scoring.sasa import *
from pyrosetta.rosetta.protocols.geometry import *

def main(pose,chain,posi):

    pdb_name_t = "./pp/pp_N"+str(posi)+"T.pdb"
    pdb_name_s = "./pp/pp_N"+str(posi)+"S.pdb"
    gly_name_t = "./gly_gly/mjg_N"+str(posi)+"T_0001.pdb"
    gly_name_s = "./gly_gly/mjg_N"+str(posi)+"S_0001.pdb"

    pose_t = pose_from_pdb(pdb_name_t)
    pose_s = pose_from_pdb(pdb_name_s)

    scorefxn = get_fa_scorefxn()

    #move lower energy pp pdbs
    if(scorefxn(pose_t)<=scorefxn(pose_s)):
        if(pass_filter(pose,chain,posi)):
            cmd="cp "+str(pdb_name_t)+" ./pp/pp_filter"
            os.system(cmd)

    else:
        if(pass_filter(pose,chain,posi)):
            cmd="cp "+str(pdb_name_s)+" ./pp/pp_filter"
            os.system(cmd)

    #move lower energy gly_pp pdbs
    if(scorefxn(pose_t)<=scorefxn(pose_s)):
        if(pass_filter(pose,chain,posi)):
            cmd="cp "+str(gly_name_t)+" ./gly_gly/gg_filter"
            os.system(cmd)

    else:
        if(pass_filter(pose,chain,posi)):
            cmd="cp "+str(gly_name_s)+" ./gly_gly/gg_filter"
            os.system(cmd)


def pass_filter(pose, chain, posi):
    if(dist_all(pose,chain,posi) and non_cys(pose,chain,posi) and non_glycan(pose,chain,posi)):
        return True
    return False

def exposed(pose, posi):
    if(is_res_exposed(pose, posi)):
        print("position is exposed")
        return True
    print("position is not exposed")
    return False

def dist_all(pose, chain, posi):
    dist = True
    c_num = get_chain_id_from_chain(chain,pose)
    for i in get_chains(pose):
        if (i != c_num):
            target = get_chain_from_chain_id(i, pose)
            if(not distant(pose,target,posi)):
                print("distance not meet")
                dist = False
    print("distance is meet")
    return dist

def distant(pose, target, posi):
    c_posi = posi
    # Get range
    start = c_start(pose,target)
    end = c_end(pose,target,start)
    for i in range(start,end+1):
        dist = (center_of_mass(pose,c_posi,c_posi)-center_of_mass(pose,i,i)).norm()
        if (dist<6):
            return False
    return True

def non_cys(pose, chain, posi):
    pdb_pos = pose.pdb_info().pdb2pose(chain,posi)
    c_aa = str(pose.residue(posi).aa())

    if(c_aa == 'AA.aa_cys'):
        print("current not cys")
        return False
    print("current is cys")
    return True

def non_glycan(pose, chain, posi):
    pdb_pos = pose.pdb_info().pdb2pose(chain,posi)
    c_aa = str(pose.residue(posi).aa())
    c_aa_plus = str(pose.residue(posi+2).aa())
    c_aa_minus = str(pose.residue(posi-2).aa())

    if((c_aa == 'AA.aa_asn' and (c_aa_plus == 'AA.aa_ser' or c_aa_plus == 'AA.aa_thr')) 
        or (c_aa_minus == 'AA.aa_asn' and (c_aa == 'AA.aa_ser' or c_aa == 'AA.aa_thr'))):
            print("this is a glycan")
            return False
    print("this is not a glycan")
    return True

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
    help = 'pose for range selection' )
parser.add_option('--chain', dest = 'chain',
    default = 'A',    # default example position
    help = 'chain for sorting' )
parser.add_option('--posi', dest = 'posi',
    default = '1',    # default example position
    help = 'position for sorting' )

(options,args) = parser.parse_args()
# PDB file option

pdb = options.pdb
chain = options.chain
posi = int(options.posi)

pose = Pose()
# -load the data from pdb_file into the pose
pose_from_file(pose, pdb)

main(pose,chain,posi)
