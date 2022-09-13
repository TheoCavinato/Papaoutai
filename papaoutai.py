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

gp_table = numpy.array(
[
    [
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0]
    ],
    [
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0]
    ],
    [
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0]
    ]
],
dtype=float)
print(gp_table[2, 2, 2])

path_table = numpy.array(
[
    [ # M = 0/0
        [ # C1 = 0/0
            [4/5, 1/5, 0],
            [0, 1, 0],
            [0, 0, 0]
        ],
        [ # C1 = 0/1
            [0, 1, 0],
            [0, 1/5, 4/5],
            [0, 0, 0]
        ],
        [ # C1 = 1/1
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]
        ]
    ],
    [ # M = 0/1
        [ # C1 = 0/0
            [4/5, 1/5, 0],
            [1/2, 1/2, 0],
            [0, 1, 0]
        ],
        [ # C1 = 0/1
            [1/2, 1/2, 0],
            [1/3, 1/3, 1/3],
            [0, 1/2, 1/2]
        ],
        [ # C1 = 1/1
            [0, 1, 0],
            [0, 1/2, 1/2],
            [0, 1/5, 4/5]
        ]
    ],
    [ # M = 1/1
        [ # C1 = 0/0
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]
        ],
        [ # C1 = 0/1
            [0, 0, 0],
            [4/5, 1/5, 0],
            [0, 1, 0]
        ],
        [ # C1 = 1/1
            [0, 0, 0],
            [0, 1, 0],
            [0, 1/5, 4/5]
        ]
    ],
],
dtype=float)

#-----------------------------------------------------------------------------#
# 2. Open the vcf containg mother, child1 and child2
#-----------------------------------------------------------------------------#

answer_table = numpy.array(
[
    [
        [
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]
        ],
        [
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]
        ],
        [
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]
        ]
    ],
    [
        [
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]
        ],
        [
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]
        ],
        [
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]
        ]
    ],
    [
        [
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]
        ],
        [
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]
        ],
        [
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]
        ]
    ],
],
dtype=float)



M, C1, C2 = args.m, args.c1, args.c2
vcf_reader = VariantFile(args.vcf, 'r')
vcf_reader.subset_samples([M, C1, C2])
for rec in vcf_reader:

    # update gp_table
    GPm, GPc1, GPc2 = rec.samples[M]['GP'], rec.samples[C1]['GP'], rec.samples[C2]['GP']
    for a in range(3):
        for b in range(3):
            for c in range(3):
                gp_table[a, b, c] = GPm[a] * GPc1[b] * GPc2[c]
    
    # combine path table with GT_table
    for a in range(3):
        for b in range(3):
            for c in range(3):
                for d in range(3):
                    answer_table[a, b, c, d] = gp_table[a, b, c] * \
                        path_table[a, b, c, d]
    
    # sum the columns corresponding to each genotypes
    GPf_0, GPf_1, GPf_2 = numpy.sum(answer_table, axis=(2, 0, 1))
    #print(rec)
    #print("Path")
    #print(path_table)
    #print("GP")
    #print(gp_table)
    #print("Answer")
    #print(answer_table)
    #print("Genotype probabilities of the father")
    print(sum([GPf_0, GPf_1, GPf_2]))
    

vcf_reader.close()