

import sys

gff3_file = sys.argv[1]



with open(gff3_file, 'r') as gff3:
    length_list = []
    for l in gff3:


            l = l.strip()

            i = l.split('\t')
            
            exon_start = int(i[3])

            exon_end = int(i[4])
            
            exon_length = exon_end - exon_start
            
            length_list.append(exon_length)
            
            #print(exon_length)
    print(sum(length_list))        


