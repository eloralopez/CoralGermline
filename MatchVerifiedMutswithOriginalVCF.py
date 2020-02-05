#DESCRIPTION: find positions of the verified mutations in the annotated VCF, then output each the mutations file with a new annotation column added
#USAGE: python3 Match.py muts.txt annotated.vcf

import sys

mutations_file = sys.argv[1]

annotated_vcf = sys.argv[2]

chrom_poslist = []
genotypelist = []
linelist = []
with open(mutations_file, 'r') as mutations:

    for line in mutations:

        # if line.startswith('chrom.pos'):
#             continue
#
#         else:

        line = line.strip()

        items = line.split('\t')

        chrom_pos = items[0]

        genotype = items[4]

        chrom_poslist.append(chrom_pos)
        genotypelist.append(genotype)
        linelist.append(line)
#     #mutations.close()
        chrom_pos_and_line = zip (chrom_poslist, linelist)

        mutations_dictionary = dict(chrom_pos_and_line)
        #print(mutations_dictionary)

with open(annotated_vcf, 'r') as vcf:

    for line2 in vcf: 
        if line2.startswith('#'): #ignores all the header lines

            continue

        else:

            line2 = line2.strip()

            items2 = line2.split('\t')

            contig = items2[0] #the #CHROM column in the VCF

            position = items2[1] #the POS column in the VCF

            concatenated = contig + "." + position #creates a string that contains both the chromosome and the position information
    
            formatfield = items2[7]
    
            things = formatfield.split(';')
            # print(things)
            ann = things[14]
            details = ann.split('|')

            mutation_type = details[1]
            #print(mutation_type)
            mutation_strength = details[2]
            # print(mutation_strength)

    
        if concatenated in mutations_dictionary:
            
            line = mutations_dictionary[concatenated]


            writeout = [line, mutation_type, mutation_strength]
    
            writeout_string = '\t'.join(writeout)

            print(writeout_string)
        
                # else:
                #     continue
                    
                    










        #
#
#
#
#
            
            