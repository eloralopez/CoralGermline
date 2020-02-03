#!/usr/bin/env python

#DESCRIPTION: Defines each putative mutations as either a transition or transversion, either  denovo or LoH, and whether it's been seen before in the population or not. Outputs the information into melted data table 

import sys

polymorphic = sys.argv[1] #example: /Users/eloralopez/Documents/SomaticMutations/OfuAug/AH09datatable20181029.txt
#popdata = sys.argv[2] #example: /Users/eloralopez/Documents/SomaticMutations/RTE_new/Genotypes_trimmed.txt

melted = sys.argv[2] #Example: /Users/eloralopez/Documents/SomaticMutations/OfuAug/meltedAH09datatable20181029.txt
intersection = sys.argv[3]
output = sys.argv[4] #example: /Users/eloralopez/Documents/SomaticMutations/OfuAug/meltedAH09datatable20181029_withtransitionandpopinfo.txt

def most_common(lst):
    return max(set(lst), key=lst.count)
def least_common(lst):
    return min(set(lst), key=lst.count)
AtoG = 0
AtoGstr = (str(AtoG))

AtoC = 0
AtoCstr = (str(AtoC))

AtoT = 0
AtoTstr = (str(AtoT))

GtoA = 0
GtoAstr = (str(GtoA))

GtoC = 0
GtoCstr = (str(GtoC))

GtoT = 0
GtoTstr = (str(GtoT))

CtoT = 0
CtoTstr = (str(CtoT))

CtoG = 0
CtoGstr = (str(CtoG))

CtoA = 0
CtoAstr = (str(CtoA))

TtoC = 0
TtoCstr = (str(TtoC))

TtoA = 0
TtoAstr = (str(TtoA))

TtoG = 0
TtoGstr = (str(TtoG))  
het2homo = 0
homo2het = 0 

with open(polymorphic, "r") as vcf_file:
    with open(output, "w") as out_file:
        header = ["chrom.pos", "sample", "ref", "alt", "genotype", "DeNovo_LoH", "TiTv", "WhattoWhat", "MutantAlleleDepth", "NormalAlleleDepth", "TrueorFalse"]
        header_output = "\t".join(header)
        out_file.write(header_output + '\n')

        for line in vcf_file:
            if line.startswith('Chrom.pos'): #skips the header line

                continue

            else:
                split = line.rstrip().split('\t')
                concatenated = split[0]
                TrueorFalse = split[15]
                ref = split[1]
    
                alt = split[2]
                genotypes= split[3:10] #only include the parental genotypes, not the sperm genotypes
                geno_list = []
                for genotype in genotypes:
                    
                    split = genotype.split(',')
                    
                    the_genotype = split[0]
                    refdepth = split[1]
                    altdepth= split[2]
                #
                    geno_list.append(the_genotype)
                    # the_genotype = '\t'.join(geno_list)
                # makelist = list(the_genotype)
                # print(geno_list)
                TiTv = "Transversion" #"Transversion" if mutation is a transverstion at locus, "Transition" otherwise
                DeNovo_or_LoH = "DeNovo" #"DeNovo" if mutation is homo2het, "LoH" if it is het2homo
                PopulationPoly = False #False if not seen in population before, True if seen before
                WhattoWhat = "placeholder"
                ToRef_or_ToAlt = "placeholder"
                coding_or_not = "NotCoding"
                majgeno = (most_common(geno_list))
                mingeno = (least_common(geno_list))

                if (majgeno == "C/T" or majgeno == "T/C") and mingeno == "C/C" in line:
                    TtoC +=1
                    het2homo+=1
                    TiTv = "Transition"
                    DeNovo_or_LoH = "LoH"
                    WhattoWhat = "TtoC"
                    # if ref = "C":
                    #     ToRef_or_ToAlt = "ToRef"
                    # else:
                    #     ToRef_or_ToAlt = "ToAlt"

                elif (mingeno == "T/C" or mingeno == "C/T") and majgeno == "T/T" in line:
                    TtoC +=1
                    homo2het +=1
                    TiTv = "Transition"
                    WhattoWhat = "TtoC"

                elif (majgeno == "A/T" or majgeno == "T/A") and mingeno == "A/A" in line:
                    TtoA +=1
                    het2homo+=1
                    DeNovo_or_LoH = "LoH"
                    WhattoWhat = "TtoA"

                elif (mingeno == "T/A" or mingeno == "A/T") and majgeno == "T/T" in line:
                    TtoA +=1
                    homo2het +=1
                    WhattoWhat = "TtoA"

                elif (majgeno == "G/T" or majgeno == "T/G") and mingeno == "G/G" in line:
                    TtoG +=1
                    het2homo+=1
                    DeNovo_or_LoH = "LoH"
                    WhattoWhat = "TtoG"

                elif (mingeno == "T/G" or mingeno == "G/T") and majgeno == "T/T" in line:
                    TtoG +=1
                    homo2het +=1
                    WhattoWhat = "TtoG"

                elif (majgeno == "C/T" or majgeno == "T/C") and mingeno == "T/T" in line:
                    CtoT +=1
                    het2homo+=1
                    TiTv = "Transition"
                    DeNovo_or_LoH = "LoH"
                    WhattoWhat = "CtoT"

                elif (mingeno == "C/T" or mingeno == "T/C") and majgeno == "C/C" in line:
                    CtoT +=1
                    homo2het +=1
                    TiTv = "Transition"
                    WhattoWhat = "CtoT"

                elif (majgeno == "C/G" or majgeno == "G/C") and mingeno == "G/G" in line:
                    CtoG +=1
                    het2homo+=1
                    DeNovo_or_LoH = "LoH"
                    WhattoWhat = "CtoG"

                elif (mingeno == "C/G" or mingeno == "G/C") and majgeno == "C/C" in line:
                    CtoG +=1
                    homo2het +=1
                    WhattoWhat = "CtoG"

                elif (majgeno == "C/A" or majgeno == "A/C") and mingeno == "A/A" in line:
                    CtoA +=1
                    het2homo+=1
                    DeNovo_or_LoH = "LoH"
                    WhattoWhat = "CtoA"

                elif (mingeno == "C/A" or mingeno == "A/C") and majgeno == "C/C" in line:
                    CtoA +=1
                    homo2het +=1
                    WhattoWhat = "CtoA"

                elif (majgeno == "A/G" or majgeno == "G/A") and mingeno == "A/A" in line:
                    GtoA +=1
                    het2homo+=1
                    TiTv = "Transition"
                    DeNovo_or_LoH = "LoH"
                    WhattoWhat = "GtoA"

                elif (mingeno == "A/G" or mingeno == "G/A") and majgeno == "G/G" in line:
                    GtoA +=1
                    homo2het +=1
                    TiTv = "Transition"
                    WhattoWhat = "GtoA"

                elif (majgeno == "C/G" or majgeno == "G/C") and mingeno == "C/C" in line:
                    GtoC +=1
                    het2homo+=1
                    DeNovo_or_LoH = "LoH"
                    WhattoWhat = "GtoC"

                elif (mingeno == "C/G" or mingeno == "G/C") and majgeno == "G/G" in line:
                    GtoC +=1
                    homo2het +=1
                    WhattoWhat = "GtoC"

                elif (majgeno == "T/G" or majgeno == "G/T") and mingeno == "T/T" in line:
                    GtoT +=1
                    het2homo+=1
                    DeNovo_or_LoH = "LoH"
                    WhattoWhat = "GtoT"

                elif (mingeno == "T/G" or mingeno == "G/T") and majgeno == "G/G" in line:
                    GtoT +=1
                    homo2het +=1
                    WhattoWhat = "GtoT"

                if (majgeno == "A/G" or majgeno == "G/A") and mingeno == "G/G" in line:
                    AtoG +=1
                    het2homo+=1
                    TiTv = "Transition"
                    DeNovo_or_LoH = "LoH"
                    WhattoWhat = "AtoG"

                elif (mingeno == "A/G" or mingeno == "G/A") and majgeno == "A/A" in line:
                    AtoG +=1
                    homo2het +=1
                    TiTv = "Transition"
                    WhattoWhat = "AtoG"

                elif (majgeno == "A/C" or majgeno == "C/A") and mingeno == "C/C" in line:
                    AtoC +=1
                    het2homo+=1
                    DeNovo_or_LoH = "LoH"
                    WhattoWhat = "AtoC"

                elif (mingeno == "A/C" or mingeno == "C/A") and majgeno == "A/A" in line:
                    AtoC +=1
                    homo2het +=1
                    WhattoWhat = "AtoC"

                elif (majgeno == "A/T" or majgeno == "T/A") and mingeno == "T/T" in line:
                    AtoT +=1
                    het2homo+=1
                    DeNovo_or_LoH = "LoH"
                    WhattoWhat = "AtoT"

                elif (mingeno == "A/T" or mingeno == "T/A") and majgeno == "A/A" in line:
                    AtoT +=1
                    homo2het +=1
                    WhattoWhat = "AtoT"


                # with open(popdata, "r") as popdata_file:
                #     for line in popdata_file:
                #         split = line.rstrip().split('\t')
                #         popdata_chrom = split[0]
                #         popdata_pos = split[1]
                #         popconcatenated = popdata_chrom + "." + popdata_pos
                #
                #         if concatenated == popconcatenated:
                #             PopulationPoly = True
                #             break
                what = WhattoWhat.split('to')
                normalallele = what[0]
                mutantallele = what[1]
                
                if mutantallele == ref:
                    mutant_allele_depth = refdepth
                    normal_allele_depth = altdepth
                    #print("MUT=REF")
                else:
                    mutant_allele_depth = altdepth
                    normal_allele_depth = refdepth
                    #print("MUT=ALT")
                    #print(concatenated, mutant_allele_depth, normal_allele_depth, refdepth, altdepth)
                
                #print(WhattoWhat,normalallele,mutantallele)
                
                with open(intersection, "r") as intersection_file:
                    intersectionconc_list = []
                    period_list = []
                    for line in intersection_file:
                        split = line.rstrip().split('\t')
                        intersection_chrom = split[0]
                        intersection_pos = split[1]
                        intersection_concatenated = intersection_chrom + "." + intersection_pos
                        # print(intersection_concatenated)
                        period = split[2]
                        
                        intersectionconc_list.append(intersection_concatenated)
                        period_list.append(period)
                        
                    intersectionconc_and_period = zip(intersectionconc_list, period_list)
                    intersectiondict = dict(intersectionconc_and_period)
                    # print(intersectiondict)
                    if concatenated in intersectiondict:
                        coding_or_not = "IsCoding"
                        # if concatenated == intersection_concatenated:
#                             coding_or_not = "IsCoding"
#                             break
                        # else:
#                             coding_or_not = "NotCoding"
                    with open(melted, "r") as melted_file:

                        for line in melted_file:

                            if line.startswith('Chrom.pos'):

                                continue

                            else:
                                meltsplit = line.rstrip().split('\t')
                                meltconcatenated=meltsplit[0]
                                meltsample=meltsplit[1]
                                meltgenotypes=meltsplit[2]
                                #print(meltgenotypes)
                                for meltgenotype in meltgenotypes:

                                    split = genotype.split(',')

                                    melt_the_genotype = split[0]
                                    #refdepth = split[1]
                                    #altdepth= split[2]
                                #print(split, melt_the_genotype, refdepth, altdepth)

                                if meltconcatenated == concatenated:
                                    AtoGstr = (str(AtoG))
                                    AtoCstr = (str(AtoC))

                                    AtoTstr = (str(AtoT))

                                    GtoAstr = (str(GtoA))

                                    GtoCstr = (str(GtoC))

                                    GtoTstr = (str(GtoT))

                                    CtoTstr = (str(CtoT))

                                    CtoGstr = (str(CtoG))

                                    CtoAstr = (str(CtoA))

                                    TtoCstr = (str(TtoC))

                                    TtoAstr = (str(TtoA))

                                    TtoGstr = (str(TtoG))
                                    # if mutantallele == ref:
    #                                     mutant_allele_depth = refdepth
    #                                     normal_allele_depth = altdepth
    #                                     #print("MUT=REF")
    #                                 else:
    #                                     mutant_allele_depth = altdepth
    #                                     normal_allele_depth = refdepth
    #                                     #print("MUT=ALT")
    #                                     print(meltconcatenated, mutant_allele_depth, normal_allele_depth, refdepth, altdepth)
                                    outlist = [meltconcatenated, meltsample, ref, alt, meltgenotypes, DeNovo_or_LoH, TiTv, WhattoWhat, mutant_allele_depth, normal_allele_depth, TrueorFalse, coding_or_not]
                                    outstring = '\t'.join(outlist)
                                    print(outstring)
                                    out_file.write(outstring+'\n')



