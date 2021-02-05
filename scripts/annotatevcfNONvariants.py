

import sys
vcf_file = sys.argv[1]

gff3_file = sys.argv[2]

with open(vcf_file, 'r') as vcf:
    for line in vcf:
                if line.startswith('#'): #ignores all the header lines
                    continue

                else:
                    line = line.strip()

                    items = line.split('\t')

                    pos = int(items[1])

                    with open(gff3_file, 'r') as gff3:

                        for l in gff3:


                                l = l.strip()

                                i = l.split('\t')
                                
                                exon_start = int(i[3])

                                exon_end = int(i[4])

                                if exon_start <= pos <= exon_end:
                                    codeornot = "coding"
                                    print("coding",pos, exon_start, exon_end)


