import optparse
import os

def main(pdb):

    sb = "/home/dwkulp/software/Rosetta/main/source/bin/relax.linuxgccrelease -in:file:s "+str(pdb)+" -out:file:o gs_relax.pdb -in:file:fullatom -relax:quick -relax:constrain_relax_to_start_coords -relax:sc_cst_maxdist 0.5"
    cmd = "/wistar/kulp/software/slurmq --sbatch \""+sb+"\""
    os.system(cmd)


#Set Up Option Parser
parser = optparse.OptionParser()
parser.add_option('--pdb', dest = 'pdb',
    default = 'loveu3000.pdb',    # default example position
    help = 'the input struture' )
(options,args) = parser.parse_args()
# PDB file option
pdb = options.pdb
main(pdb)