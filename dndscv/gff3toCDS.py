import sys
import random
import string

def randomString(stringLength=10):
    """Generate a random string of fixed length """
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))

gff3_file = sys.argv[1]

with open(gff3_file, 'r') as gff3:
    
    strand_type = "placeholder"
    gene_name_list = []
    feature_start_list = []
    feature_end_list = []
    start_end_list = []
    
    for line in gff3: 


            line = line.strip()
            
            items = line.split('\t')
            
            chromosome = items[0]
            
            feature_type = items[2]
            
            feature_start = items[3]
            
            feature_end = items[4]
            
            length = int(feature_end) - int(feature_start)
            
            strand = items[6]
            
            
            if strand == "+":

                strand_type = "1"
            else:
                strand_type = "-1"
            
            if feature_type == "exon":
            
                gene_name = randomString()
                gene_id = gene_name
                cds_id = gene_name
                
                chromosome_start = feature_start
                chromosome_end = feature_end
                
                feature_start_adj = 1
                
                feature_end_adj = int(length)
                
                gene_name_list.append(gene_name)
                feature_start_list.append(feature_start)
                feature_end_list.append(feature_end)
                
                start_end = gene_name, feature_start, feature_end
                start_end_string = '\t'.join(start_end)
                start_end_list.append(start_end_string)
                #print(start_end_list)
                
            elif feature_type == "CDS":
                # if int(feature_start) > int(dict_start_end)
                for numbers in start_end_list:
                    number = numbers.split('\t')
                    name = number[0]
                    chromosome_start = int(number[1])
                    chromosome_end = int(number[2])
                    if int(feature_start) >= chromosome_start and int(feature_start) <= chromosome_end:
                        gene_id = name
                        cds_id = gene_id+randomString()
                        feature_start_adj = (int(feature_start) - chromosome_start) + 1
                        feature_end_adj = int(feature_start_adj) + int(length)
                        
                    #print(start, end)
                writeout = gene_id,gene_name,cds_id, chromosome, str(chromosome_start), str(chromosome_end), str(feature_start_adj), str(feature_end_adj), str(length), strand_type
                writeout_string = '\t'.join(writeout)
                print(writeout_string)
    gff3.close()
    
genename_and_featurestart = zip(gene_name_list, feature_start_list)

genename_and_featureend = zip(gene_name_list, feature_end_list)

start_and_end = zip(feature_start_list, feature_end_list)

dict_featurestart = dict(genename_and_featurestart)

dict_featureend = dict(genename_and_featureend)

dict_start_end = dict(start_and_end)

start_end_list = feature_start_list, feature_end_list 

# print(start_end_list)
                #print(dict_featurestart)
                    
# with open(gff3_file, 'r') as gff3:
#
#     strand_type = "placeholder"
#     gene_name_list = []
#     feature_start_list = []
#     feature_end_list = []
#
#     for line in gff3:
#
#
#             line = line.strip()
#
#             items = line.split('\t')
#
#             chromosome = items[0]
#
#             feature_type = items[2]
#
#             feature_start = items[3]
#
#             feature_end = items[4]
#
#             length = int(feature_end) - int(feature_start)
#
#             strand = items[6]
#
#
#             if strand == "+":
#
#                 strand_type = "1"
#             else:
#                 strand_type = "-1"
#
#             for line1 in start_end_list:
#                 line1 = line1.strip()
#
#                 items1 = line1.split('\t')
#                 print(items1)
                # chrom_start = items1[0]
#                 chrom_end = items1[1]
#                 print(chrom.start, chrom.end)
            # if feature_type == "CDS":# or feature_type == "exon":
#
#                 if int(feature_start) > dict_start_end[] in dict_featurestart:
#
#                     print(feature_start)
#
#                     gene_name = dict_featurestart[feature_start]
#
#                     print(gene_name)
#
#                 if feature_end in dict_featureend:
#
#                     gene_name = dict_featureend[featureend]

                #print("gene.id","gene.name","cds.id", chromosome, "chr.start", "chr.end", feature_start, feature_end, length, strand_type)