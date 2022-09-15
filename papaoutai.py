from pysam import VariantFile
import argparse, numpy

"""
Reconstruction of a paternal genomes using two children and a mother
"""

#-----------------------------------------------------------------------------#
# 0. User's Parameters 
#-----------------------------------------------------------------------------#

parser = argparse.ArgumentParser()
parser.add_argument("--vcf", required=True, help= "Vcf should contain the id\
    of the father and of the children")
parser.add_argument("--m", required=True, help="Maternal id")
parser.add_argument("--c1", required=True, help="Child1 id")
parser.add_argument("--c2", required=True, help="Child2 id")
parser.add_argument("--f", required=False, help="Father id for validation")
args = parser.parse_args()

#-----------------------------------------------------------------------------#
# 1. Prepare two table
# - path_table -> hold the probabilities of each possible paternal genotypes
# for each combination of mother, child1 and child2 genotypes
# - gp_table -> hold the probability of the genotype of the mother, child1
# and child2 given by glimpse
#-----------------------------------------------------------------------------#

path_table_gp_way = numpy.array([
            [4/5, 1/5, 0],
            [0, 1, 0],
            [0, 1, 0],
            [0, 1/5, 4/5],

            [4/5, 1/5, 0],
            [1/2, 1/2, 0],
            [0, 1, 0],
            [1/2, 1/2, 0],
            [1/3, 1/3, 1/3],
            [0, 1/2, 1/2], 
            [0, 1, 0],
            [0, 1/2, 1/2],
            [0, 1/5, 4/5],

            [0, 1/5, 4/5],
            [0, 1, 0],
            [0, 1, 0],
            [4/5, 1/5, 0]
],
dtype=float)


#-----------------------------------------------------------------------------#
# 2. Open the vcf containg mother, child1 and child2
#-----------------------------------------------------------------------------#

class GPvector:
    def __init__(self, p_GPm, p_GPc1, p_GPc2):
        self.vector = numpy.array([0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0], dtype=float)
        self.GPm, self.GPc1, self.GPc2 = p_GPm, p_GPc1, p_GPc2
        self.go_first_branch, self.go_third_branch = True, True

        # populate vector
        self.maternal_branches_validation()
        self.maternal_branches_modification()
        self.multiply_GPs()
    
    def maternal_branches_validation(self):
        # first branch
        if self.GPc1[:2] == (0, 0) or self.GPc2[:2] == (0, 0):
            self.go_first_branch = False

        # last branch
        if self.GPc1[1:] == (0, 0) or self.GPc2[1:] == (0, 0):
            self.go_third_branch = False
    
    def maternal_branches_modification(self):
        if not self.go_first_branch and not self.go_third_branch:
            self.GPm = (0, 1, 0)
        elif not self.go_first_branch:
            self.GPm = (0, self.GPm[1]/sum(self.GPm[1:]), self.GPm[2]/sum(self.GPm[1:]))
        elif not self.go_third_branch:
            self.GPm = (self.GPm[0]/sum(self.GPm[:2]), self.GPm[1]/sum(self.GPm[:2]), 0)
    
    def multiply_GPs(self):
        v=0
        # first branch
        if self.go_first_branch:
            for c1, c2 in zip((0,0,1,1), (0,1,0,1)):
                self.vector[v] = self.GPm[0]*(self.GPc1[c1]/sum(self.GPc1[:2]))*(self.GPc2[c2]/sum(self.GPc2[:2]))
                v+=1
        else:
            v+=4
        # second branch
        for c1, c2 in zip((0,0,0,1,1,1,2,2,2), (0,1,2,0,1,2,0,1,2)):
            self.vector[v] = self.GPm[1]*self.GPc1[c1]*self.GPc2[c2]
            v+=1
        # third branch
        if self.go_third_branch:
            for c1, c2 in zip((1,1,2,2), (1,2,1,2)):
                self.vector[v] = self.GPm[2]*(self.GPc1[c1]/sum(self.GPc1[1:]))*(self.GPc2[c2]/sum(self.GPc2[1:]))
                v+=1
    
    def __str__(self):
        return str(self.vector)

M, C1, C2 = args.m, args.c1, args.c2
vcf_reader = VariantFile(args.vcf, 'r')
vcf_reader.subset_samples([M, C1, C2])

for rec in vcf_reader:

    GPm, GPc1, GPc2 = rec.samples[M]['GP'], rec.samples[C1]['GP'], rec.samples[C2]['GP']

    # do not take in account impossible cases
    if not (((GPc1[:2] == (0, 0) or GPc2[:2] == (0, 0)) and GPm[1:] == (0, 0)) or \
        (GPc1[1:] == (0, 0) or GPc2[1:] == (0, 0) and GPm[:2] == (0,0))):

        gp_vec = GPvector(GPm, GPc1, GPc2)
        GPf_0, GPf_1, GPf_2= numpy.dot(gp_vec.vector, path_table_gp_way)
        print(GPf_0, GPf_1, GPf_2)

        if sum([GPf_0, GPf_1, GPf_2]) > 1.01 or sum([GPf_0, GPf_1, GPf_2]) < 0.99:
            raise Exception("No 0", sum([GPf_0, GPf_1, GPf_2]))

vcf_reader.close()
