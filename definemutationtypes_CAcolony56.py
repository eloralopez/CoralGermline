#!/usr/bin/env python

#DESCRIPTION: Defines each putative mutations as either a transition or transversion, either  denovo or LoH, and whether it's been seen before in the population or not. Outputs the information into melted data table 
# python3 definemutationtypeS_CAcolony56.py *ii*.txt melted*ii*txt output*muts.txt
import sys

datatable = sys.argv[1] #example: headCAcolony56_ii_20200317.txt
#popdata = sys.argv[2] #example: /Users/eloralopez/Documents/SomaticMutations/RTE_new/Genotypes_trimmed.txt

melted = sys.argv[2] #Example: /Users/eloralopez/Documents/SomaticMutations/OfuAug/meltedAH09datatable20181029.txt
numberofparents = int(sys.argv[3]) 
# output = sys.argv[4] #example: /Users/eloralopez/Documents/SomaticMutations/OfuAug/meltedAH09datatable20181029_withtransitionandpopinfo.txt

def most_common(lst):
    return max(set(lst), key=lst.count)
def least_common(lst):
    return min(set(lst), key=lst.count)
    
def definemutations(polymorphic, typeofmutation, specifiedgenotypes):
    with open(polymorphic, "r") as mutations_table:
        # with open(output, "w") as out_file:
        #     header = ["chrom.pos", "sample", "ref", "alt", "genotype", "DeNovo_LoH", "TiTv", "WhattoWhat", "TrueorFalse"]
        #     header_output = "\t".join(header)
        #     out_file.write(header_output + '\n')

            for line in mutations_table:
                if line.startswith('chrom.pos'): #skips the header line

                    continue

                else:
                    split = line.rstrip().split('\t')
                    concatenated = split[0]
                    TrueorFalse = split[((numberofparents*2)+3)]
                    TypeofMutation = split[((numberofparents*2)+4)]
                    MutantSampleID = split[((numberofparents*2)+5)]
                    ref = split[1]
    
                    alt = split[2]
                    if TypeofMutation ==typeofmutation:
                        
                        genotypes= split[3:specifiedgenotypes]
                    
                    # print(genotypes)
                    
                    # if TypeofMutation == "SomaticMutation":
#
#                         genotypes= split[3:((numberofparents)+3)] #only include the parental genotypes, not the sperm genotypes
#
#                     elif TypeofMutation == "UniqueGermlineMutation":
#                         genotypes= split[3:((numberofparents)+3)] + split[spermsample+3]
#                         print(genotypes)
                        
                    # elif TypeofMutation == "GlobalGermlineMutation":
                    #     genotypes = split[3:((numberofparents)+4)]   
                    
                    #print(genotypes)
                        geno_list = []
                        for genotype in genotypes:

                            genosplit = genotype.split(',')

                            the_genotype = genosplit[0]
                            the_depth = genosplit[1]
                            #print(genotype)
                        #
                            geno_list.append(the_genotype)
                            # the_genotype = '\t'.join(geno_list)
                        # makelist = list(the_genotype)
                        # print(geno_list)
                        TiTv = "placeholder" #"Transversion" if mutation is a transverstion at locus, "Transition" otherwise
                        DeNovo_or_LoH = "placeholder" #"DeNovo" if mutation is homo2het, "LoH" if it is het2homo
                        WhattoWhat = "placeholder"
                        ToRef_or_ToAlt = "placeholder"
                    #
                        majgeno = (most_common(geno_list))
                        mingeno = (least_common(geno_list))
                        #print(geno_list)
                        if (majgeno == "C/T" or majgeno == "T/C") and mingeno == "C/C" in line:
                            # TtoC +=1
        #                     het2homo+=1
                            TiTv = "Transition"
                            DeNovo_or_LoH = "LoH"
                            WhattoWhat = "TtoC"
                            # if ref == "C":
                            #     ToRef_or_ToAlt = "ToRef"
                            # else:
                            #     ToRef_or_ToAlt = "ToAlt"

                        elif (mingeno == "T/C" or mingeno == "C/T") and majgeno == "T/T" in line:
                            # TtoC +=1
        #                     homo2het +=1
                            TiTv = "Transition"
                            WhattoWhat = "TtoC"
                            DeNovo_or_LoH = "DeNovo"

                        elif (majgeno == "A/T" or majgeno == "T/A") and mingeno == "A/A" in line:
                            # TtoA +=1
                            # het2homo+=1
                            TiTv = "Transversion"
                            DeNovo_or_LoH = "LoH"
                            WhattoWhat = "TtoA"

                        elif (mingeno == "T/A" or mingeno == "A/T") and majgeno == "T/T" in line:
                            # TtoA +=1
        #                     homo2het +=1
                            TiTv = "Transversion"
                            DeNovo_or_LoH = "DeNovo"
                            WhattoWhat = "TtoA"

                        elif (majgeno == "G/T" or majgeno == "T/G") and mingeno == "G/G" in line:
                            # TtoG +=1
        #                     het2homo+=1
                            TiTv = "Transversion"
                            DeNovo_or_LoH = "LoH"
                            WhattoWhat = "TtoG"

                        elif (mingeno == "T/G" or mingeno == "G/T") and majgeno == "T/T" in line:
                            # TtoG +=1
        #                     homo2het +=1
                            TiTv = "Transversion"
                            DeNovo_or_LoH = "DeNovo"
                            WhattoWhat = "TtoG"

                        elif (majgeno == "C/T" or majgeno == "T/C") and mingeno == "T/T" in line:
                            # CtoT +=1
        #                     het2homo+=1
                            TiTv = "Transition"
                            DeNovo_or_LoH = "LoH"
                            WhattoWhat = "CtoT"

                        elif (mingeno == "C/T" or mingeno == "T/C") and majgeno == "C/C" in line:
                            # CtoT +=1
        #                     homo2het +=1
                            TiTv = "Transition"
                            DeNovo_or_LoH = "DeNovo"
                            WhattoWhat = "CtoT"

                        elif (majgeno == "C/G" or majgeno == "G/C") and mingeno == "G/G" in line:
                            # CtoG +=1
        #                     het2homo+=1
                            TiTv = "Transversion"
                            DeNovo_or_LoH = "LoH"
                            WhattoWhat = "CtoG"

                        elif (mingeno == "C/G" or mingeno == "G/C") and majgeno == "C/C" in line:
                            # CtoG +=1
        #                     homo2het +=1
                            TiTv = "Transversion"
                            DeNovo_or_LoH = "DeNovo"
                            WhattoWhat = "CtoG"

                        elif (majgeno == "C/A" or majgeno == "A/C") and mingeno == "A/A" in line:
                            # CtoA +=1
        #                     het2homo+=1
                            TiTv = "Transversion"
                            DeNovo_or_LoH = "LoH"
                            WhattoWhat = "CtoA"

                        elif (mingeno == "C/A" or mingeno == "A/C") and majgeno == "C/C" in line:
                            # CtoA +=1
        #                     homo2het +=1
                            TiTv = "Transversion"
                            DeNovo_or_LoH = "DeNovo"
                            WhattoWhat = "CtoA"

                        elif (majgeno == "A/G" or majgeno == "G/A") and mingeno == "A/A" in line:
                            # GtoA +=1
        #                     het2homo+=1
                            TiTv = "Transition"
                            DeNovo_or_LoH = "LoH"
                            WhattoWhat = "GtoA"

                        elif (mingeno == "A/G" or mingeno == "G/A") and majgeno == "G/G" in line:
                            # GtoA +=1
        #                     homo2het +=1
                            TiTv = "Transition"
                            DeNovo_or_LoH = "DeNovo"
                            WhattoWhat = "GtoA"

                        elif (majgeno == "C/G" or majgeno == "G/C") and mingeno == "C/C" in line:
                            # GtoC +=1
        #                     het2homo+=1
                            TiTv = "Transversion"
                            DeNovo_or_LoH = "LoH"
                            WhattoWhat = "GtoC"

                        elif (mingeno == "C/G" or mingeno == "G/C") and majgeno == "G/G" in line:
                            # GtoC +=1
        #                     homo2het +=1
                            TiTv = "Transversion"
                            DeNovo_or_LoH = "DeNovo"
                            WhattoWhat = "GtoC"

                        elif (majgeno == "T/G" or majgeno == "G/T") and mingeno == "T/T" in line:
                            # GtoT +=1
        #                     het2homo+=1
                            TiTv = "Transversion"
                            DeNovo_or_LoH = "LoH"
                            WhattoWhat = "GtoT"

                        elif (mingeno == "T/G" or mingeno == "G/T") and majgeno == "G/G" in line:
                            # GtoT +=1
        #                     homo2het +=1
                            TiTv = "Transversion"
                            DeNovo_or_LoH = "DeNovo"
                            WhattoWhat = "GtoT"

                        elif (majgeno == "A/G" or majgeno == "G/A") and mingeno == "G/G" in line:
                            # AtoG +=1
        #                     het2homo+=1
                            TiTv = "Transition"
                            DeNovo_or_LoH = "LoH"
                            WhattoWhat = "AtoG"

                        elif (mingeno == "A/G" or mingeno == "G/A") and majgeno == "A/A" in line:
                            # AtoG +=1
        #                     homo2het +=1
                            TiTv = "Transition"
                            DeNovo_or_LoH = "DeNovo"
                            WhattoWhat = "AtoG"

                        elif (majgeno == "A/C" or majgeno == "C/A") and mingeno == "C/C" in line:
                            # AtoC +=1
        #                     het2homo+=1
                            TiTv = "Transversion"
                            DeNovo_or_LoH = "LoH"
                            WhattoWhat = "AtoC"

                        elif (mingeno == "A/C" or mingeno == "C/A") and majgeno == "A/A" in line:
                            # AtoC +=1
        #                     homo2het +=1
                            TiTv = "Transversion"
                            DeNovo_or_LoH = "DeNovo"
                            WhattoWhat = "AtoC"

                        elif (majgeno == "A/T" or majgeno == "T/A") and mingeno == "T/T" in line:
                            # AtoT +=1
                            # het2homo+=1
                            TiTv = "Transversion"
                            DeNovo_or_LoH = "LoH"
                            WhattoWhat = "AtoT"

                        elif (mingeno == "A/T" or mingeno == "T/A") and majgeno == "A/A" in line:
                            # AtoT +=1
        #                     homo2het +=1
                            TiTv = "Transversion"
                            DeNovo_or_LoH = "DeNovo"
                            WhattoWhat = "AtoT"
                            
                        if TrueorFalse == "False":
                            
                            if DeNovo_or_LoH == "DeNovo":
                                spermrep1 = split[numberofparents+3:numberofparents+4]
                                spermrep2 = split[numberofparents+4:numberofparents+5]
                                #print(spermrep1, spermrep2)
                                for genotype in spermrep1:
                                    gsplit = genotype.split(',')
                                    refdepth = int(gsplit[1])
                                    altdepth = int(gsplit[2])
                                for genotype in spermrep2:
                                    gsplit2 = genotype.split(',')
                                    refdepth2 = int(gsplit2[1])
                                    altdepth2 = int(gsplit2[2])
                                        
                                if refdepth !=0 and altdepth!=0 and refdepth2 !=0 and altdepth2 != 0:
                                    TrueorFalse ="True"
                                    #print("actually true")
                                    #print(refdepth, altdepth, refdepth2, altdepth2)

                        with open(melted, "r") as melted_file:

                            for line in melted_file:

                                if line.startswith('chrom.pos'):

                                    continue

                                else:
                                    meltsplit = line.rstrip().split('\t')
                                    meltconcatenated=meltsplit[0]
                                    meltsample=meltsplit[1]
                                    meltgenotypes=meltsplit[2]
                                    #print(meltconcatenated, meltsample, meltgenotypes)
                                    for meltgenotype in meltgenotypes:

                                        split = genotype.split(',')

                                        melt_the_genotype = split[0]
                                        melt_the_depth = split[1]
                    #
                    #
                                    if meltconcatenated ==concatenated:

                                        outlist = [meltconcatenated, meltsample, ref, alt, meltgenotypes, DeNovo_or_LoH, TiTv, WhattoWhat, TrueorFalse, TypeofMutation, MutantSampleID]
                                        outstring = '\t'.join(outlist)
                                        print(outstring)

somatic = definemutations(datatable, "SomaticMutation", (numberofparents+3))

uglm1 = definemutations(datatable, "UniqueGermlineMutation", ((numberofparents*2)+3))
#

# uglm2 = definemutations(datatable,(numberofparents,genotype_2ndsperm))
#
# uglm3 = definemutations(datatable,)
#
# if numberofparents ==8:
#     uglm4 = definemutations(datatable,)

gglm = definemutations(datatable, "GlobalGermlineMutation", (numberofparents+4))    
                    
                    
                    
                    
                    


