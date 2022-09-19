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
parser.add_argument("--threshold", required=False, help="Probability threshold to\
    use a reconstructed heterzoygous site", default=0.9)
args = parser.parse_args()

#-----------------------------------------------------------------------------#
# 1. Prepare table holding the probabilities of each possible paternal genotypes
# for each combination of mother, child1 and child2 genotypes
#-----------------------------------------------------------------------------#

path_table_gp_way = numpy.array([
            [4/5, 1/5, 0],
            [0, 1, 0],
            [0, 0, 0],
            [0, 1, 0],
            [0, 1/5, 4/5],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],

            [4/5, 1/5, 0],
            [1/2, 1/2, 0],
            [0, 1, 0],
            [1/2, 1/2, 0],
            [1/3, 1/3, 1/3],
            [0, 1/2, 1/2], 
            [0, 1, 0],
            [0, 1/2, 1/2],
            [0, 1/5, 4/5],

            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],

            [0, 0, 0],
            [4/5, 1/5, 0],
            [0, 1, 0],

            [0, 0, 0],
            [0, 1, 0],
            [0, 1/5, 4/5]
],
dtype=float)

m_path = [i//9 for i in range(27)]
c1_path = [(i//3)%3 for i in range(27)]
c2_path = [i%3 for i in range(27)]

#-----------------------------------------------------------------------------#
# 2. Open the vcf containg mother, child1 and child2
# and reconstruct the sites
#-----------------------------------------------------------------------------#

M, C1, C2 = args.m, args.c1, args.c2
vcf_reader = VariantFile(args.vcf, 'r')
vcf_reader.subset_samples([M, C1, C2])

imputation = []
impute_none = 0
for rec in vcf_reader:

    GPm, GPc1, GPc2 = rec.samples[M]['GP'], rec.samples[C1]['GP'], rec.samples[C2]['GP']

    # do not take in account impossible cases AND only consider biallelic sites
    if not (((GPc1[:2] == (0, 0) or GPc2[:2] == (0, 0)) and GPm[1:] == (0, 0)) or \
        (GPc1[1:] == (0, 0) or GPc2[1:] == (0, 0) and GPm[:2] == (0,0))) \
            and len(rec.alts)==1:

        gp_vec = [GPm[m]*GPc1[c1]*GPc2[c2] for m, c1, c2 in zip(m_path, c1_path, c2_path)]
        GPf_0, GPf_1, GPf_2= numpy.dot(gp_vec, path_table_gp_way)
        normGPf_0, normGPf_1, normGPf_2 = GPf_0/(GPf_0+GPf_1+GPf_2), GPf_1/(GPf_0+GPf_1+GPf_2), GPf_2/(GPf_0+GPf_1+GPf_2)
        genos_prob = [normGPf_0, normGPf_1, normGPf_2 ]

        if GPf_0 != GPf_1 != GPf_2 != 1/3:
            imputation.append(genos_prob.index(max(genos_prob)))
        else:
            imputation.append(None)
            impute_none+=1
    else:
        imputation.append(None)
        impute_none+=1

vcf_reader.close()

#-----------------------------------------------------------------------------#
# 3. Write the imputed genotypes of the father
#-----------------------------------------------------------------------------#

geno_to_alleles = {
    0: (0,0),
    1: (0,1),
    2: (1,1),
    None: (None, None)
}

# write header with father
vcf_reader = VariantFile(args.vcf, 'r')
vcf_reader.subset_samples([args.c1]) # I change this in dividual to FATHER at each SNPs

vcf_header = vcf_reader.header.__str__().split('\n')
print('\n'.join(vcf_header[:-2]))
print('\t'.join(['#CHROM', "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "FATHER"]))

writing_None = 0
itr=0
for rec, imput in zip(vcf_reader, imputation):

    # put paternal GT
    rec.samples[args.c1]['GT'] = geno_to_alleles[imput]
    print(str(rec)[:-1])

vcf_reader.close()
