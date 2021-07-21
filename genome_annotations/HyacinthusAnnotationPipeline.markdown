**Repeat Annotation:**
```
RepeatModeler-2.0.1/BuildDatabase -name repeatsV1 -engine ncbi   Ahya_genome/jasmine-sta1387-mb-hirise-eyx9u_12-15-2020__final_assembly.fasta
RepeatModeler-2.0.1/RepeatModeler -pa 8 -engine ncbi -database repeatsV1 2>&1
```
**First MAKER Round:**

MAKER rnd1_opts.ctl file; other .ctl files were left at default settings
```
#-----Genome (these are always required)
genome=Ahya_genome/jasmine-sta1387-mb-hirise-eyx9u_12-15-2020__final_assembly.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est= transcript_data/ahyacinthus_july2014/ahya.fasta,transcript_data/84259_FinalCoralDinoformapping.fasta #set of ESTs or assembled mRNA-seq in fasta format
altest=transcript_data/moya.transcripts.fasta,transcript_data/amil_v2.01.transcripts.fasta,transcript_data/aten_july2014/aten.fasta #EST/cDNA sequence file in fasta format from an alternate organism
est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=protein_data/adigi.peptides.fasta,protein_data/Amil.all.maker.proteins.fasta  #protein sequence file in fasta format (i.e. from mutiple organisms)
protein_gff=  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org= #select a model organism for RepBase masking in RepeatMasker
rmlib=repeat_annotation/repeat_modeler.v1/consensi.fa.classified #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff= #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm= #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species= #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
run_evm=0 #run EvidenceModeler, 1 = yes, 0 = no
est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
snoscan_meth= #-O-methylation site fileto have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no
allow_overlap=0 #allowed gene overlap fraction (value from 0 to 1, blank for default)

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=16 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1 #skip genome contigs below this length (under 10kb are often useless)
pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)
split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
min_intron=20 #minimum intron length (used for alignment polishing)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes
tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files
```
Running MAKER
```
maker/bin/maker -c 16  -base wgenome-rnd1 rnd1_opts.ctl rnd1_bopts.ctl rnd1_exe.ctl -fix_nucleotides
```
Extract transcript and protein fastas from MAKER output
```
maker/bin/fasta_merge -d  wgenome-rnd1_master_datastore_index.log 
```
Extract gff files from MAKER output
```
maker/bin/gff3_merge -n -s -d wgenome-rnd1.maker.output/wgenome-rnd1_master_datastore_index.log > wgenome-rnd1.all.maker.noseq.gff
awk '{ if ($2 ~ "repeat") print $0 }' wgenome-rnd1.maker.output/wgenome-rnd1.all.maker.noseq.gff > wgenome-rnd1.all.maker.repeats.gff
awk '{ if ($2 == "est2genome") print $0 }' wgenome-rnd1.maker.output/wgenome-rnd1.all.maker.noseq.gff > wgenome-rnd1.all.maker.est2genome.gff
awk '{ if ($2 == "protein2genome") print $0 }' wgenome-rnd1.maker.output/wgenome-rnd1.all.maker.noseq.gff  > wgenome-rnd1.all.maker.protein2genome.gff
```

**Training _ab initio_ Gene Prediction Models (Round 1):**

SNAP training
```
maker/bin/maker2zff -x 0.25 -l 50 -d wgenome-rnd1.maker.output/wgenome-rnd1_master_datastore_index.log 
ls | awk '{split($1,a,".");print "mv", $1,"wgenome-rnd1."a[2]}' | bash
SNAP/fathom wgenome-rnd1.ann wgenome-rnd1.dna -gene-stats > gene-stats.log
SNAP/fathom wgenome-rnd1.ann wgenome-rnd1.dna -validate > validate.log
SNAP/fathom wgenome-rnd1.ann wgenome-rnd1.dna -categorize 1000 > categorize.log
SNAP/fathom wgenome-rnd1.ann wgenome-rnd1.dna -export 1000 -plus > uni-plus.log 2>&1
mkdir params
cd params
SNAP/forge ../export.ann ../export.dna > ../forge.log
cd ..
SNAP/hmm-assembler.pl wgenome-rnd1 params/ > wgenome-rnd1.hmm
```

Augustus training
```
awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }'  wgenome-rnd1.maker.output/wgenome-rnd1.all.maker.noseq.gff |  \
awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }' | \
bedtools getfasta -fi Ahya_genome/jasmine-sta1387-mb-hirise-eyx9u_12-15-2020__final_assembly.fasta  -bed - -fo wgenome-rnd1.all.maker.transcripts1000.fasta
busco/bin/busco -i wgenome-rnd1.all.maker.transcripts1000.fasta -o wgenome_rnd1 -l eukaryota_odb10/ -m genome -c 8  --augustus_parameters='--progress=true' \
--augustus_species human -f --augustus --long --force
```

**Second MAKER Round:**

MAKER rnd2_opts.ctl file; other .ctl files were left at default settings
```
#-----Genome (these are always required)
genome=Ahya_genome/jasmine-sta1387-mb-hirise-eyx9u_12-15-2020__final_assembly.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est= #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff=wgenome-rnd1.maker.output/wgenome-rnd1.all.maker.est2genome.gff #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=  #protein sequence file in fasta format (i.e. from mutiple organisms)
protein_gff=wgenome-rnd1.maker.output/wgenome-rnd1.all.maker.protein2genome.gff  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org= #select a model organism for RepBase masking in RepeatMasker
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff=wgenome-rnd1.maker.output/wgenome-rnd1.all.maker.repeats.gff #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm=gene_prediction_training/snap/wg1/wgenome-rnd1.hmm #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species=BUSCO_wgenome_rnd1 #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
run_evm=0 #run EvidenceModeler, 1 = yes, 0 = no
est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
snoscan_meth= #-O-methylation site fileto have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no
allow_overlap=0 #allowed gene overlap fraction (value from 0 to 1, blank for default)

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=24 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1 #skip genome contigs below this length (under 10kb are often useless)
pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)
split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
min_intron=20 #minimum intron length (used for alignment polishing)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes
tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files
```

Run MAKER
```
maker/bin/maker -c 16  -base wgenome-rnd2 rnd2_opts.ctl rnd2_bopts.ctl rnd2_exe.ctl -fix_nucleotides
```
Extract transcript and protein fastas from MAKER output
```
maker/bin/fasta_merge -d  wgenome-rnd2_master_datastore_index.log 
```

Extract gff files from MAKER output
```
maker/bin/gff3_merge -n -s -d wgenome-rnd2.maker.output/wgenome-rnd2_master_datastore_index.log > wgenome-rnd2.all.maker.noseq.gff
awk '{ if ($2 ~ "repeat") print $0 }' wgenome-rnd2.maker.output/wgenome-rnd2.all.maker.noseq.gff > wgenome-rnd2.all.maker.repeats.gff
awk '{ if ($2 == "est_gff:est2genome") print $0 }' wgenome-rnd2.maker.output/wgenome-rnd2.all.maker.noseq.gff > wgenome-rnd2.all.maker.est2genome.gff
awk '{ if ($2 == "protein_gff:protein2genome") print $0 }' wgenome-rnd2.maker.output/wgenome-rnd2.all.maker.noseq.gff  > wgenome-rnd2.all.maker.protein2genome.gff
```
**Training _ab initio_ Gene Prediction Models (round 2):**

SNAP training
```
maker/bin/maker2zff -x 0.25 -l 50 -d wgenome-rnd2.maker.output/wgenome-rnd2_master_datastore_index.log 
ls | awk '{split($1,a,".");print "mv", $1,"wgenome-rnd2."a[2]}' | bash
SNAP/fathom wgenome-rnd2.ann wgenome-rnd2.dna -gene-stats > gene-stats.log
SNAP/fathom wgenome-rnd2.ann wgenome-rnd2.dna -validate > validate.log
SNAP/fathom wgenome-rnd2.ann wgenome-rnd2.dna -categorize 1000 > categorize.log
SNAP/fathom wgenome-rnd2.ann wgenome-rnd2.dna -export 1000 -plus > uni-plus.log 2>&1
mkdir params
cd params
SNAP/forge ../export.ann ../export.dna > ../forge.log
cd ..
SNAP/hmm-assembler.pl wgenome-rnd2 params/ > wgenome-rnd2.hmm
```

Augustus training
```
awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }'  wgenome-rnd2.maker.output/wgenome-rnd2.all.maker.noseq.gff |  \
awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }' | \
bedtools getfasta -fi Ahya_genome/jasmine-sta1387-mb-hirise-eyx9u_12-15-2020__final_assembly.fasta  -bed - -fo wgenome-rnd2.all.maker.transcripts1000.fasta
busco/bin/busco -i wgenome-rnd2.all.maker.transcripts1000.fasta -o wgenome_rnd2 -l eukaryota_odb10/ -m genome -c 8  --augustus_parameters='--progress=true' \
--augustus_species human -f --augustus --long --force
```

**Third MAKER Round:**

MAKER rnd3_opts.ctl file; other .ctl files were left at default settings
```
#-----Genome (these are always required)
genome=Ahya_genome/jasmine-sta1387-mb-hirise-eyx9u_12-15-2020__final_assembly.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est= #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff=wgenome-rnd2.maker.output/wgenome-rnd2.all.maker.est2genome.gff #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=  #protein sequence file in fasta format (i.e. from mutiple organisms)
protein_gff=wgenome-rnd2.maker.output/wgenome-rnd2.all.maker.protein2genome.gff  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org= #select a model organism for RepBase masking in RepeatMasker
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff=wgenome-rnd2.maker.output/wgenome-rnd1.all.maker.repeats.gff #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm=gene_prediction_training/snap/wg2/wgenome-rnd2.hmm #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species=BUSCO_wgenome_rnd2 #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
run_evm=0 #run EvidenceModeler, 1 = yes, 0 = no
est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
snoscan_meth= #-O-methylation site fileto have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no
allow_overlap=0 #allowed gene overlap fraction (value from 0 to 1, blank for default)

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=24 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1 #skip genome contigs below this length (under 10kb are often useless)
pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)
split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
min_intron=20 #minimum intron length (used for alignment polishing)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes
tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files
```
Running MAKER
```
maker/bin/maker -c 16  -base wgenome-rnd3 rnd3_opts.ctl rnd3_bopts.ctl rnd3_exe.ctl -fix_nucleotides
```
Extract transcript and protein fastas from MAKER output
```
maker/bin/fasta_merge -d  wgenome-rnd3_master_datastore_index.log 
```
Extract gff files from MAKER output
```
maker/bin/gff3_merge -n -s -d wgenome-rnd3.maker.output/wgenome-rnd3_master_datastore_index.log > wgenome-rnd3.all.maker.noseq.gff
awk '{ if ($2 ~ "repeat") print $0 }' wgenome-rnd3.maker.output/wgenome-rnd3.all.maker.noseq.gff > wgenome-rnd3.all.maker.repeats.gff
awk '{ if ($2 == "est_gff:est2genome") print $0 }' wgenome-rnd3.maker.output/wgenome-rnd2.all.maker.noseq.gff > wgenome-rnd3.all.maker.est2genome.gff
awk '{ if ($2 == "protein_gff:protein2genome") print $0 }' wgenome-rnd3.maker.output/wgenome-rnd3.all.maker.noseq.gff  > wgenome-rnd3.all.maker.protein2genome.gff
```
