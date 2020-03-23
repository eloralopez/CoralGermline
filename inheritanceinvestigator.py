#DESCRIPTION:this script calls a multi-sample VCF that contains genotypes for both parents and sperm samples from a single colony. The "fileinfo" function defined at the top of the script generates a list of all of the sites where the two replicate (you designate which two) samples in a single file have the SAME GENOTYPE. Then, the  list of matches is compared to the other parent genotypes in the file and keeps only those where the matching replicates have a UNIQUE genotype not seen in any of the other parent samples. Then the corresponding sperm genotype is checked and deemed "inherited" if the sperm genotype matches the unique parent genotype and "not inherited" if it does not match. 

#USAGE: python3 inheritanceinvestigator.py x.vcf

import sys
import re

input1 = sys.argv[1]
replicate1 = int(sys.argv[2]) #column number for the first replicate of the first parent ie CAP22-1
replicate2 = int(sys.argv[3]) #column number for the second replicate of the first parent ie CAP22-2

# replicate3 = int(sys.argv[4]) #THIS IS THE SAMPLE YOU WILL LOOK FOR UNIQUE MUTATIONS IN
# replicate4 = int(sys.argv[5]) #THIS IS THE SAMPLE YOU WILL LOOK FOR UNIQUE MUTATIONS IN
numberofparents = int(sys.argv[4]) #the number of parent samples in the VCF.
doublenumber = numberofparents*2
# spermrep1 = int(sys.argv[7]) # the position of the desired sperm sample in the list of sample columns
# spermrep2 = int(sys.argv[8]) # the position of the desired sperm sample in the list of sample columns


replicate3= replicate1 + 2

replicate4 = replicate2 + 2

replicate5 = replicate3 + 2

replicate6 = replicate4 + 2

replicate7 = "NULL"

replicate8 = "NULL"

spermrep1 = replicate1 + numberofparents

spermrep2 = replicate2 + numberofparents

spermrep3 = replicate3 + numberofparents

spermrep4 = replicate4 + numberofparents

spermrep5 = replicate5 + numberofparents

spermrep6 = replicate6 + numberofparents


if numberofparents == 8:
    

    replicate7 = replicate5 + 2

    replicate8 = replicate6 + 2

    spermrep7 = replicate7 + numberofparents

    spermrep8 = replicate8 + numberofparents

#print(replicate1, replicate2, replicate3, replicate4, replicate5, replicate6, replicate7, replicate8)
#print(spermrep1, spermrep2, spermrep3, spermrep4, spermrep5, spermrep6, spermrep7, spermrep8)

def fileinfo(inputfile, rep1, rep2):
    listofoutlist = [] #starts an empty list to which the matches between replicate libraries will be added later

    with open(inputfile, 'r') as f: #opens whichever file is specified

        for line in f:
            if line.startswith('##'): #ignores all the header lines
                continue

            if line.startswith('#CHROM'): #ignores all the header lines

                line = line.strip()

                items = line.split('\t')

                samplenames = items[9:]
                samplenames = '\t'.join(samplenames)
                #print(samplenames)
                header = ["chrom.pos", "ref", "alt", samplenames, samplenames]
                #print(header)
                header = '\t'.join(header)
                #print(header)
                #print(samplenames_out)
                listofoutlist.append(header)
                #header = "Chrom.pos","ref", "alt", "sample1","sample2","sample3","sample4","sample5","sample6","sample7", "sample8", "sample9", "sample10", "sample11", "sample12", "TrueorFalse"
                #header_string = '\t'.join(header)
                #print(header_string)

            else:

                line = line.strip()

                items = line.split('\t')

                contig = items[0] #the #CHROM column in the VCF

                position = items[1] #the POS column in the VCF

                concatenated = contig + "." + position #creates a string that contains both the chromosome and the position information

                ref_allele = items[3] #the REF column in the VCF


                alt_allele = items[4] #the ALT column in the VCF

                genotypes = items[9:] #each of the samplename columns in the VCF (e.g. AHAS57-1)

                geno_list=[] #empty list to which genotypes will be added later

                genoplusdepthlist=[] #empty list to which genotype depths will be added later

                for genotype in genotypes:

                    genos = genotype.split(':') #splits the genotype column up by ":"
                    alleles = re.split('; |, |\||\/',genos[0])#genos[0].split('/') #splits the first part of the genotype column into its two respective alleles
                    #print(alleles)
                    GQscore = genos[3]

                    alleledepths = genos[1].split(',') #depth per allele at locus
                    #FROM GATK: AD is the unfiltered allele depth, i.e. the number of reads that support each of the reported alleles. All reads at the position (including reads that did not pass the variant caller’s filters) are included in this number, except reads that were considered uninformative. Reads are considered uninformative when they do not provide enough statistical evidence to support one allele over another.

                    refdepth = alleledepths[0]

                    altdepth = alleledepths[1]

#                    totaldepth_atlocus = genos[2] #totaldepth at locus

                    #FROM GATK: DP is the filtered depth, at the sample level. This gives you the number of filtered reads that support each of the reported alleles. You can check the variant caller’s documentation to see which filters are applied by default. Only reads that passed the variant caller’s filters are included in this number. However, unlike the AD calculation, uninformative reads are included in DP.

                    if alleles[0] == ".": #ignore the sites where the depth is './.'
                        continue
                    if int(alleles[0]) == 0:

                        a1 = ref_allele #gives the corresponding character (A,C, G, or T) for the allele

        #
                    else:

                        a1 = alt_allele #gives the corresponding character (A,C, G, or T) for the allele

                    if int(alleles[1]) == 0:

                        a2 = ref_allele #gives the corresponding character (A,C, G, or T) for the allele

                    else:

                        a2 = alt_allele #gives the corresponding character (A,C, G, or T) for the allele
                    genostring = a1 + '/' + a2 #gives the genotype (e.g. A/T)
                    genoPLUSdepth = genostring + "," + refdepth + ","  + altdepth + "," + GQscore #gives a string of the genotype plus its corresponding depths plus its GQ

                    genoplusdepthlist.append(genoPLUSdepth)
                    genoPLUSdepth = '\t'.join(genoplusdepthlist)

                    geno_list.append(genostring)
                    genostring = '\t'.join(geno_list)

                if geno_list[rep1] == geno_list[rep2]: #outputs just the sites where the two replicate libraries are the same genotype
                #if all(x==geno_list[0] for x in geno_list):
                    outlist = [concatenated, ref_allele, alt_allele, genostring, genoPLUSdepth] # collect each of these for each line where the genotypes are the same
                    outstring = '\t'.join(outlist)

                    listofoutlist.append(outstring) #append each line to the list started up at the top
                    #print(concatenated, outlist)
                #else:
                    #continue #skips all sites that have different genotype calls betweeen the two technical replicates
    return listofoutlist
    print(listoutoutlist)
    return dictionary #returns the full list of genotype matches in the file

    f.close()

genos1=fileinfo(input1, replicate1, replicate2) #CAP22 #calls the function for a given set of replicate libraries (whichever 2 columns you want)   #FLAG
#print(genos1)
genos2=fileinfo(input1, replicate3, replicate4) #CAP23

genos3=fileinfo(input1, replicate5, replicate6) #CAP24

if numberofparents == 8:

    genos4=fileinfo(input1, replicate7, replicate8) #CAP26
 #calls the function for CAP24
#print(genos2)
spermgenos=fileinfo(input1, spermrep1, spermrep2) #CAS22

spermgenos2=fileinfo(input1, spermrep3, spermrep4) #CAS23

spermgenos3=fileinfo(input1, spermrep5, spermrep6) #CAS24

if numberofparents == 8:

    spermgenos4=fileinfo(input1, spermrep7, spermrep8) #CAS26

#print(spermgenos)
#print(genos2)
# genos3=fileinfo(input1, 4, 5)
# #
# genos4=fileinfo(input1, 6, 7)

conc1list = []
geno1list = []

def makedictionary(fileinfo_output, replicatenumber):
    conc1list = []
    geno1list = []
    for line1 in fileinfo_output: #goes line by line in the first input file's list of matches
            #print(line1)
            if line1.startswith('chrom.pos'):
                continue
            else:
                line1 = line1.strip()

                items1 = line1.split('\t')
                #print(items1)
                conc1 = items1[0]

                genotypes1 = items1[3:]
                geno1 = genotypes1[replicatenumber] #FLAG
                #print(geno1)

                conc1list.append(conc1)
                geno1list.append(geno1)


    conc1_and_geno1 = zip(conc1list, geno1list)

    dictionary = dict(conc1_and_geno1)
    return dictionary

dictOfWords = makedictionary(genos1, replicate1)
dictOfWords2 = makedictionary(genos2, replicate3)
#dictOfWords3 = makedictionary(genos3, replicate5)
#dictOfWords4 = makedictionary(genos4, replicate7)


#print(dictOfWords)

spermdict = makedictionary(spermgenos, spermrep1)
spermdict2 = makedictionary(spermgenos2, spermrep3)
spermdict3 = makedictionary(spermgenos3, spermrep5)
if numberofparents == 8:
    spermdict4 = makedictionary(spermgenos4, spermrep7)

#print(spermdict)
#print(dictOfWords)
# conc3list = []
# geno3list = []
# for line3 in genos3: #goes line by line in the first input file's list of matches
#         #print(line3)
#         line3 = line3.strip()
#
#         items3 = line3.split('\t')
#         #print(items3)
#         conc3 = items3[0]
#
#         genotypes3 = items3[3:]
#         geno3 = genotypes3[1]
#         #print(geno3)
#
#         conc3list.append(conc3)
#         geno3list.append(geno3)
#
#
# conc3_and_geno3 = zip(conc3list, geno3list)
#
# dictOfWords3 = dict(conc3_and_geno3)
# conc4list = []
# geno4list = []
# for line4 in genos4: #goes line by line in the first input file's list of matches
#         #print(line4)
#         line4 = line4.strip()
#
#         items4 = line4.split('\t')
#         #print(items4)
#         conc4 = items4[0]
#
#         genotypes4 = items4[3:]
#         geno4 = genotypes4[1]
#         #print(geno4)
#
#         conc4list.append(conc4)
#         geno4list.append(geno4)
#
#
# conc4_and_geno4 = zip(conc4list, geno4list)
#
# dictOfWords4 = dict(conc4_and_geno4)
#print(dictOfWords)
#
conc2list = []
geno2list = []

#to find CAP23 mutants:
def findmutations(genos2, dictOfWords, replicatesample, replicate5, replicate6, replicate7, replicate8, spermrep3, spermrep4, spermdict2):
    samplenames = []
    for line2 in genos2: #now goes line by line in the list of matches from the second set of replicates. this set is the MUTANT set.
            #print(line2)
            if line2.startswith('chrom.pos'):
                #print(line2)
                line2 = line2.strip()

                items2 = line2.split('\t')
                total = (numberofparents * 2) + 3 #this gives the number of columns from the first line to output. It is the sum of the number of parent samples, the number of sperm samples (Which is equal to the number of parent samples), and the three columns that denote chrom.pos, ref, and alt
                samples = items2[0:total]
                    #print(header)
                for ele in samples:
                    sub = ele.split(',')
                    samplenames.append(sub)
                
                # header = items2[0:total]
                # header = '\t'.join(header)
                # finalheader = [header, "TypeofMutation", "TrueorFalse"]
                # finalheader = '\t'.join(finalheader)
                #
                # print(finalheader) # this is the headerline for the file
            else:
                line2 = line2.strip()

                items2 = line2.split('\t')
                #print(items2)
                conc2 = items2[0]
                ref = items2[1]
                #print(ref)
                alt = items2[2]

                genotypes2= items2[3:]
                #print(genotypes2)
                geno2= genotypes2[replicatesample]

                #deal with somatic mutations for the second set of parent replicates:
                if conc2 in dictOfWords:
                    geno1 = dictOfWords[conc2] #outputs just the geno1 values that match conc2 (that is, just the sites that are present in the lists from genos1 and genos2)
                    #print(geno1)
                    if numberofparents == 6:
                        if geno1 != geno2 and geno2 != genotypes2[replicate5] and geno2 != genotypes2[replicate6]:# and geno2 != genotypes2[replicate7] and geno2 != genotypes2[replicate8]: #this indicates that the mutation found in replicate3 and replicate4 (now referred to as geno2) is UNIQUE among the parent samples; that genotype is never seen at that site in any other parent sample
                            genolist = genotypes2[7:]
                            Type_of_Mutation = "SomaticMutation"
                            #print(*samplenames[replicate3+3])
                            Mutant_Sample_ID = samplenames[replicatesample+3]
                            if genotypes2[replicatesample] == genotypes2[spermrep3] == genotypes2[spermrep4]: #if the corresponding sperm pools are the same genotype as the mutant parent, inheritance is true. If not, inheritance is false.
                                #print("yaaa")
                                Match = True
                            else:
                                #perhaps add extra conditions to deal with the low minor allele frequency sperm genotypes?
                                Match = False
                            writeout = [conc2, ref, alt, genotypes2[12],genotypes2[13], genotypes2[14], genotypes2[15], genotypes2[16], genotypes2[17], genotypes2[18], genotypes2[19], genotypes2[20], genotypes2[21], genotypes2[22], genotypes2[23], str(Match), str(Type_of_Mutation), *Mutant_Sample_ID]
                            #here is your writeout of all of the somatic mutations unique to a particular branch
                            writeout_string = '\t'.join(writeout)

                            #print("somatic mutation")
                            print(writeout_string)
                    elif numberofparents == 8:
                        if geno1 != geno2 and geno2 != genotypes2[replicate5] and geno2 != genotypes2[replicate6] and geno2 != genotypes2[replicate7] and geno2 != genotypes2[replicate8]: #this indicates that the mutation found in replicate3 and replicate4 (now referred to as geno2) is UNIQUE among the parent samples; that genotype is never seen at that site in any other parent sample
                            genolist = genotypes2[7:]
                            Type_of_Mutation = "SomaticMutation"
                            Mutant_Sample_ID = samplenames[replicatesample+3]
                            if genotypes2[replicatesample] == genotypes2[spermrep3] == genotypes2[spermrep4]: #if the corresponding sperm pools are the same genotype as the mutant parent, inheritance is true. If not, inheritance is false.
                                #print("yaaa")
                                Match = True
                            else:
                                #perhaps add extra conditions to deal with the low minor allele frequency sperm genotypes?
                                Match = False
                            genotypes2string = genotypes2[doublenumber:(doublenumber*2)]    
                            writeout = [conc2, ref, alt, *genotypes2string, str(Match), str(Type_of_Mutation), *Mutant_Sample_ID]    
                            #writeout = [conc2, ref, alt, genotypes2[16],genotypes2[17], genotypes2[18], genotypes2[19], genotypes2[20], genotypes2[21], genotypes2[22], genotypes2[23], genotypes2[24], genotypes2[25], genotypes2[26], genotypes2[27], genotypes2[28], genotypes2[29], genotypes2[30], genotypes2[31], str(Match), str(Type_of_Mutation), *Mutant_Sample_ID]
                            #here is your writeout of all of the somatic mutations unique to a particular branch
                            writeout_string = '\t'.join(writeout)

                            print(writeout_string)
                if conc2 in spermdict2:

                    spermgeno = spermdict2[conc2]
                    #print(spermgeno, geno2)
                    #genotypes_minusspermgeno = genotypes2[0:8]+ genotypes2[12:16]
                    #x= genotypes2.remove(genotypes2[spermrep3])

                    #genotypes_minusspermgeno = x.remove(x[spermrep4])
                    #print(genotypes_minusspermgeno)
                    justparentgenotypes = genotypes2[0:numberofparents]
                   
                    justspermgenotypes = genotypes2[numberofparents:doublenumber]
                    #print(justspermgenotypes)
                    genotypes_minusspermgeno = genotypes2[0:doublenumber]
                    del genotypes_minusspermgeno[spermrep3:(spermrep4+1)]
                    #print(genotypes2)
                    #print(genotypes_minusspermgeno)
                    #print(conc2, genotypes2)
                    #print(genotypes2[spermrep3])
                    if not any(ele in spermgeno for ele in genotypes_minusspermgeno): # if the sperm genotype does not match the genotype of any other parent or sperm sample:
                         Type_of_Mutation = "UniqueGermlineMutation"
                         Match = "NA"
                         Mutant_Sample_ID = samplenames[spermrep3+3]
                         #print(samplenames[spermrep3+3])
                         #print(spermgeno, genotypes_minusspermgeno[2], conc2)
                         genolist = genotypes2[7:]
                         genotypes2string = genotypes2[doublenumber:(doublenumber*2)]
                         writeout_uniqueglm = [conc2, ref, alt, *genotypes2string, str(Match), str(Type_of_Mutation), *Mutant_Sample_ID]
                         #writeout_uniqueglm = [conc2, ref, alt, genotypes2[16],genotypes2[17], genotypes2[18], genotypes2[19], genotypes2[20], genotypes2[21], genotypes2[22], genotypes2[23], genotypes2[24], genotypes2[25], genotypes2[26], genotypes2[27], genotypes2[28], genotypes2[29], genotypes2[30], genotypes2[31], str(Match), str(Type_of_Mutation), *Mutant_Sample_ID]
    #                     #here is your writeout of all of the somatic mutations unique to a particular branch
                         writeout_uniqueglm_string = '\t'.join(writeout_uniqueglm)
                         print(writeout_uniqueglm_string)

                    # if all(x==justspermgenotypes[0] for x in justspermgenotypes) and not any(ele in spermgeno for ele in justparentgenotypes): # if the sperm genotype matches all other sperm genotypes but none of the parent genotypes:
                    #     Type_of_Mutation = "GlobalGermlineMutation"
                    #     Match = "NA"
                    #     genolist = genotypes2[7:]
                    #     writeout_globalglm = [conc2, ref, alt, genotypes2[16],genotypes2[17], genotypes2[18], genotypes2[19], genotypes2[20], genotypes2[21], genotypes2[22], genotypes2[23], genotypes2[24], genotypes2[25], genotypes2[26], genotypes2[27], genotypes2[28], genotypes2[29], genotypes2[30], genotypes2[31], str(Match), str(Type_of_Mutation)]
                    #     #here is your writeout of all of the somatic mutations unique to a particular branch
                    #
                    #     writeout_globalglm_string = '\t'.join(writeout_globalglm)
                    #     #print("Global germ line mutation")
                    #     print(writeout_globalglm_string)
# #                     print(conc2)
               # print(genotypes_minusspermgeno)
                # if spermgeno !=



#headerline and global glms (keeping this out of the findmutations function because you only need to do it once:)
for line2 in genos2: #now goes line by line in the list of matches from the second set of replicates. this set is the MUTANT set.
        #print(line2)
        if line2.startswith('chrom.pos'):
            #print(line2)
            line2 = line2.strip()

            items2 = line2.split('\t')
            total = (numberofparents * 2) + 3 #this gives the number of columns from the first line to output. It is the sum of the number of parent samples, the number of sperm samples (Which is equal to the number of parent samples), and the three columns that denote chrom.pos, ref, and alt
            header = items2[0:total]
            header = '\t'.join(header)
            finalheader = [header, "TrueorFalse","TypeofMutation", "MutantSampleID"]
            finalheader = '\t'.join(finalheader)

            print(finalheader) # this is the headerline for the file
        else:
            line2 = line2.strip()

            items2 = line2.split('\t')
            #print(items2)
            conc2 = items2[0]
            ref = items2[1]
            #print(ref)
            alt = items2[2]

            genotypes2= items2[3:]
            #print(genotypes2)
            geno2= genotypes2[replicate3]

            if conc2 in spermdict2:

                spermgeno = spermdict2[conc2]
                #print(spermgeno, geno2)
                #genotypes_minusspermgeno = genotypes2[0:10]+ genotypes2[12:16]
                #print(genotypes_minusspermgeno)
                justparentgenotypes = genotypes2[0:numberofparents]
                #print(justparentgenotypes)
                justspermgenotypes = genotypes2[numberofparents:doublenumber]

                if all(x==justspermgenotypes[0] for x in justspermgenotypes) and not any(ele in spermgeno for ele in justparentgenotypes): # if the sperm genotype matches all other sperm genotypes but none of the parent genotypes:
                    Type_of_Mutation = "GlobalGermlineMutation"
                    Match = "NA"
                    Mutant_Sample_ID = "NA"
                    genolist = genotypes2[7:]
                    genotypes2string = genotypes2[doublenumber:(doublenumber*2)]
                    writeout_globalglm = [conc2, ref, alt, *genotypes2string, str(Match), str(Type_of_Mutation), str(Mutant_Sample_ID)]
                    #writeout_globalglm = [conc2, ref, alt, genotypes2[16],genotypes2[17], genotypes2[18], genotypes2[19], genotypes2[20], genotypes2[21], genotypes2[22], genotypes2[23], genotypes2[24], genotypes2[25], genotypes2[26], genotypes2[27], genotypes2[28], genotypes2[29], genotypes2[30], genotypes2[31], str(Match), str(Type_of_Mutation), str(Mutant_Sample_ID)]
                    #here is your writeout of all of the somatic mutations unique to a particular branch

                    writeout_globalglm_string = '\t'.join(writeout_globalglm)
                    print(writeout_globalglm_string)

#if numberofparents == 8:                   
sample1 = findmutations(genos1, dictOfWords2, replicate1, replicate5, replicate6, replicate7, replicate8, spermrep1, spermrep2, spermdict)                
sample2 = findmutations(genos2, dictOfWords, replicate3, replicate5, replicate6, replicate7, replicate8, spermrep3, spermrep4, spermdict2)
sample3 = findmutations(genos3, dictOfWords, replicate5, replicate3, replicate4, replicate7, replicate8, spermrep5, spermrep6, spermdict3)
if numberofparents == 8:                   

    sample4 = findmutations(genos4, dictOfWords, replicate7, replicate3, replicate4, replicate6, replicate7, spermrep7, spermrep8, spermdict4)


                    
                    
                    
                
            
                


