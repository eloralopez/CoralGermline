
setwd("~/Documents/GitHub/CoralGermline/")
library("reshape2")
library(ggplot2)
library(stringr)
library(sciplot)
library(sinaplot)
library(ggforce)
library(gridExtra)
library(dplyr)
library(patchwork)
library(ggpubr)
sessionInfo()

noreplicatesmutationcount<-c(35290, 86641, 31446, 80126, 94972, 92588, 95348, 78499, 47429, 97193, 90480)

withreplicatesmutationcount<- c( 10584, 16468, 10930, 18738, 25898, 22058, 24558, 12250, 14706, 14022, 18714)
plot(noreplicatesmutationcount,withreplicatesmutationcount)
proportion = withreplicatesmutationcount/noreplicatesmutationcount
percentdiff = (noreplicatesmutationcount - withreplicatesmutationcount)/ noreplicatesmutationcount
mean(percentdiff)
mean(proportion)
df = rbind(noreplicatesmutationcount,withreplicatesmutationcount)
barplot(df,beside=TRUE,ylim=c(0,100000), ylab="Number of first-pass putative mutations")
CAcolony60<-read.delim("CAcolony60_ii_output20191126.txt")
mdf<-melt(CAcolony60, id.vars="Chrom.pos", measure.vars=c("sample1","sample2","sample3","sample4", "sample5", "sample6", "sample7", "sample8", "sample9", "sample10", "sample11", "sample12"), value.name="genotype",variable.name="sample")
write.table(mdf, file="meltedCAcolony60_ii_output20191126.txt",sep="\t",quote=FALSE, row.name=FALSE)

CAcolony60_CAP23<-read.delim("CAcolony60_ii_output20191126CAP23.txt")
mdf<-melt(CAcolony60_CAP23, id.vars="Chrom.pos", measure.vars=c("sample1","sample2","sample3","sample4", "sample5", "sample6", "sample7", "sample8", "sample9", "sample10", "sample11", "sample12"), value.name="genotype",variable.name="sample")
write.table(mdf, file="meltedCAcolony60_ii_output20191126CAP23.txt",sep="\t",quote=FALSE, row.name=FALSE)

CAcolony60_CAP24<-read.delim("CAcolony60_ii_output20191126CAP24.txt")
mdf<-melt(CAcolony60_CAP24, id.vars="Chrom.pos", measure.vars=c("sample1","sample2","sample3","sample4", "sample5", "sample6", "sample7", "sample8", "sample9", "sample10", "sample11", "sample12"), value.name="genotype",variable.name="sample")
write.table(mdf, file="meltedCAcolony60_ii_output20191126CAP24.txt",sep="\t",quote=FALSE, row.name=FALSE)

CAcolony56_CAP12<-read.delim("CAcolony56_ii_output20191205CAP12.txt")
mdf<-melt(CAcolony56_CAP12, id.vars="Chrom.pos", measure.vars=c("sample1","sample2","sample3","sample4", "sample5", "sample6", "sample7", "sample8", "sample9"), value.name="genotype",variable.name="sample")
write.table(mdf, file="meltedCAcolony56_ii_output20191205CAP12.txt",sep="\t",quote=FALSE, row.name=FALSE)

CAcolony56_CAP6<-read.delim("CAcolony56_ii_output20191205CAP6.txt")
mdf<-melt(CAcolony56_CAP6, id.vars="Chrom.pos", measure.vars=c("sample1","sample2","sample3","sample4", "sample5", "sample6", "sample7", "sample8", "sample9"), value.name="genotype",variable.name="sample")
write.table(mdf, file="meltedCAcolony56_ii_output20191205CAP6.txt",sep="\t",quote=FALSE, row.name=FALSE)

CAcolony56_CAP8<-read.delim("CAcolony56_ii_output20191205CAP8.txt")
mdf<-melt(CAcolony56_CAP8, id.vars="Chrom.pos", measure.vars=c("sample1","sample2","sample3","sample4", "sample5", "sample6", "sample7", "sample8", "sample9"), value.name="genotype",variable.name="sample")
write.table(mdf, file="meltedCAcolony56_ii_output20191205CAP8.txt",sep="\t",quote=FALSE, row.name=FALSE)


CAcolony65_CAP11<-read.delim("CAcolony65_ii_output20200107CAP11.txt")
mdf<-melt(CAcolony65_CAP11, id.vars="Chrom.pos", measure.vars=c("sample1","sample2","sample3","sample4", "sample5", "sample6", "sample7", "sample8"), value.name="genotype",variable.name="sample")
write.table(mdf, file="meltedCAcolony65_ii_output20200107CAP11.txt",sep="\t",quote=FALSE, row.name=FALSE)

CAcolony65_CAP9<-read.delim("CAcolony65_ii_output20200107CAP9.txt")
mdf<-melt(CAcolony65_CAP9, id.vars="Chrom.pos", measure.vars=c("sample1","sample2","sample3","sample4", "sample5", "sample6", "sample7", "sample8"), value.name="genotype",variable.name="sample")
write.table(mdf, file="meltedCAcolony65_ii_output20200107CAP9.txt",sep="\t",quote=FALSE, row.name=FALSE)

par(mfrow=c(1,1))

##to look at all data combined together:
files<-list.files(path="~/Documents/GitHub/CoralGermline/WithSpermReplicates", pattern="*.ann.txt", full.names=T, recursive=FALSE) #path to all the files you want to include in the analysis
metadata= NULL
for (i in 1:length(files)) { 
  file =files[i]
  data<-read.delim(file) #read in each file in "files"
  data<-data.frame(data) # transform the data from each file into a dataframe
  base<-basename(file)
  colony<-strsplit(base, "\\_")[[1]][1]
  print(colony)
  len<-nrow(data) 
  colonyrep<-rep(colony, len)
  withcolony<-data.frame(data, colonyrep) #combines the colonyname column with each dataframe
  metadata <- rbind(metadata, withcolony) #adds each dataframe to the overall metatadata, so that the information from all of the files are now in "metadata"
}

genoanddepth<-(metadata$genotype) #names the column
split<-str_split_fixed(genoanddepth, ",", 4) #split the genotype, depths, and GQ score into their own separate strings

genotypes<-split[,1] #defines the genotype as the first string in "split"
position<-metadata$chrom.pos #names the column
positionsplit<-str_split_fixed(position, "[.]", 2) #split the chromosome number and the position on the chromosome into their own separate strings

chr<-positionsplit[,1] #defines the chromosome as the first string in "positionsplit"
pos<-positionsplit[,2] #defines the position as the second string in "positionsplit"
#totaldepth<-as.numeric(split[,2])
refdepth<-as.numeric(split[,2]) #number of reference reads
altdepth<-as.numeric(split[,3]) #number of alternate reads
what<-metadata$WhattoWhat
allelesplit<-str_split_fixed(what, "to", 2) #splits the normal and mutant alleles into different strings
normalallele<-allelesplit[,1]
mutantallele<-allelesplit[,2]
totaldepth<-refdepth+altdepth #sum of all reference and alternate reads
GQscore<-as.numeric(split[,4])
mutationtype<-metadata$MutationType #eg intergenic_region 3_prime_UTR_variant 5_prime_UTR_premature_start_codon_gain_variant 5_prime_UTR_variant
mutationstrength<-metadata$MutationStrength #eg HIGH LOW MODERATE MODIFIER
mutant_alleledepth = rep("A", nrow(metadata))
for (i in 1:nrow(metadata)){
  if (mutantallele[i] == metadata$alt[i]) {
    mutant_alleledepth[i] = altdepth[i]
    #print(mutant_alleledepth[i])
  } else {
    mutant_alleledepth[i] = refdepth[i]
    #print(mutant_alleledepth[i])
  }  #print("ALT", mutant_alleledepth, refdepth)
}
#print(mutant_alleledepth[1:10])
#print(refdepth[1:10])
#print(altdepth[1:10])

metadatadf<-data.frame("chrom.pos" = metadata$chrom.pos, "chrom"=chr, "pos"=pos,	"sample"= metadata$sample, "ref" = metadata$ref, "alt" = 	metadata$alt, "normal_allele"= normalallele, "mutant_allele" = mutantallele, "mutant_allele_depth" = as.numeric(mutant_alleledepth), "genotype"= genotypes, "totaldepth"=totaldepth, 	"refdepth"=refdepth, "altdepth"=altdepth, "GQscore"= GQscore,	"GoH_or_LoH"=metadata$DeNovo_LoH, "Ti/Tv"=metadata$TiTv, 	"WhattoWhat" = metadata$WhattoWhat,"TrueorFalse" =metadata$TrueorFalse, "MutationClass"=metadata$TypeofMutation, "MutantParent1"=metadata$MutantParent1,"MutantParent2"=metadata$MutantParent2, "MutantSperm1"=metadata$MutantSperm1, "MutantSperm2"=metadata$MutantSperm2,"MutationType"=mutationtype, "MutationStrength"=mutationstrength, "ColonyName"=metadata$colonyrep, stringsAsFactors=FALSE)#  ColonyName"=metadata$colonyrep)
i <- sapply(metadatadf, is.factor) #gives TRUE or FALSE for whether each column is a factor
metadatadf[i] <- lapply(metadatadf[i], as.character) #sets columns to be characters
metadatadf[which(metadatadf$sample==metadatadf$MutantParent1),]
DepthMeansdf<-aggregate(totaldepth~chrom.pos, metadatadf, 			FUN=mean) #calculate mean depth per locus

DepthMinsdf<-aggregate(totaldepth~chrom.pos, metadatadf, 			FUN=min) #colculate min depth per locus

GQaverage<-aggregate(GQscore~chrom.pos, metadatadf, FUN=mean) #calculate average GQ perlocus

GQmin<-aggregate(GQscore~chrom.pos, metadatadf, FUN=min) #calculate min GQ per locus



metadatadf.00<-merge(metadatadf, DepthMeansdf[, c("chrom.pos", 	"totaldepth")], by="chrom.pos")
metadatadf.0<-merge(metadatadf.00, GQaverage[,c("chrom.pos","GQscore")], by="chrom.pos")
metadatadf.0<-merge(metadatadf.0, DepthMinsdf[,c("chrom.pos","totaldepth")], by="chrom.pos")
metadatadf.1<-merge(metadatadf.0, GQmin[,c("chrom.pos","GQscore")], by="chrom.pos")
metadatadf.0<-subset(metadatadf.1, startsWith(metadatadf.1$chrom, 'chr')==TRUE)
DeNovos<-subset(metadatadf.0, GoH_or_LoH=="DeNovo")
#somatic denovos:
withoutmutparents<-subset(DeNovos,  MutationClass == "SomaticMutation" & (sample != MutantParent1 & sample != MutantParent2)) #subset the samples that are not the mutantparents
withoutsperm<-subset(withoutmutparents, startsWith(withoutmutparents$sample, 'CAS')==FALSE) #subsets the non-sperm samples
withoutsperm_freq<-as.data.frame(table(withoutsperm$chrom.pos)) #gives how many nonmutant parent samples are present per locus
trueDeNovos_somatic<-subset(withoutsperm, refdepth == 0 | altdepth ==0) #subsets just the nonmutant parent samples that have either refdepth=0 or altdepth=0. this eliminates the nonmutants for which the putatively mutant allele is present in "nonmutants" at low frequency
trueDeNovos_somatic_freq<-as.data.frame(table(trueDeNovos_somatic$chrom.pos)) #gives  how many nonmutant parent samples are present per locus once they have been subsetted
withoutsperm_freq_matched<- withoutsperm_freq[match(trueDeNovos_somatic_freq$Var1, withoutsperm_freq$Var1),] #eliminates the loci that are not present at all in trueDeNovos_somatic_freq
total <- merge(trueDeNovos_somatic_freq, withoutsperm_freq_matched,by="Var1") #now the number of samples from trueDeNovos_somatic_freq is Freq.x and withoutsperm_freq_matched is Freq.y
allcorrect<-subset(total, Freq.x == Freq.y) #subsets just the loci where Freq.x and Freq.y are equal

#nondups<-subset(trueDeNovos_somatic ,duplicated(chrom.pos)==FALSE | duplicated(chrom.pos, fromLast=TRUE)==FALSE)

mutantparent1<-subset(DeNovos, sample==MutantParent1)#| sample==MutantParent2) #pulls out just the mutantparent1s
correspondingmutantparent1<-mutantparent1[match(allcorrect$Var1,mutantparent1$chrom.pos),] #subsets just the mutantparent1 loci that appear in allcorrect

mutantparent2<-subset(DeNovos, sample==MutantParent2)#| sample==MutantParent2)
correspondingmutantparent2<-mutantparent2[match(allcorrect$Var1,mutantparent2$chrom.pos),]

mutantsperm1<-subset(DeNovos, sample==MutantSperm1)# | sample==MutantSperm2)
correspondingmutantsperm1<-mutantsperm1[match(allcorrect$Var1,mutantsperm1$chrom.pos),]

mutantsperm2<-subset(DeNovos, sample==MutantSperm2)
correspondingmutantsperm2<-mutantsperm2[match(allcorrect$Var1,mutantsperm2$chrom.pos),]

trueDenovos_somatic_mpandms<-rbind(correspondingmutantparent1, correspondingmutantparent2, correspondingmutantsperm1, correspondingmutantsperm2) #this is the dataset that has the mutantparent and mutantsperm samples at just the loci that appear in allcorrect
as.data.frame(table(trueDenovos_somatic_mpandms$chrom.pos))

#uglm denovos:
withoutmutsperm_uglmmetadatadf<- subset(DeNovos, MutationClass == "UniqueGermlineMutation" & (sample != MutantSperm1 & sample != MutantSperm2)) #subsets the non-MutantSperms
withoutparents<- subset(withoutmutsperm_uglmmetadatadf, startsWith(withoutmutsperm_uglmmetadatadf$sample, 'CAP')==FALSE) #subsets the non-parent samples
withoutparents_freq<-as.data.frame(table(withoutparents$chrom.pos))
trueDeNovos_uglm<-subset(withoutparents, refdepth == 0 | altdepth ==0) #subsets just the nonmutant sperm samples that have either refdepth=0 or altdepth=0. this eliminates the nonmutants for which the putatively mutant allele is present in "nonmutants" at low frequency
trueDeNovos_uglm_freq<-as.data.frame(table(trueDeNovos_uglm$chrom.pos)) #gives how many nonmutant sperm samples are present per locus once they have been subsetted
withoutparents_freq_matched<- withoutparents_freq[match(trueDeNovos_uglm_freq$Var1, withoutparents_freq$Var1),] #eliminates the loci that are not present at all in trueDeNovos_uglm_freq
total_uglm <- merge(trueDeNovos_uglm_freq, withoutparents_freq_matched,by="Var1") #now the number of samples from trueDeNovos_uglm_freq is Freq.x and withoutparents_freq_matched is Freq.y
allcorrect_uglm<-subset(total_uglm, Freq.x == Freq.y) #subsets just the loci where Freq.x and Freq.y are equal

corresponding_uglmmutantparent1<-mutantparent1[match(allcorrect_uglm$Var1,mutantparent1$chrom.pos),] ##subsets just the mutantparent1 loci that appear in allcorrect_uglm

corresponding_uglmmutantparent2<-mutantparent2[match(allcorrect_uglm$Var1,mutantparent2$chrom.pos),]

corresponding_uglmmutantsperm1<-mutantsperm1[match(allcorrect_uglm$Var1,mutantsperm1$chrom.pos),]

corresponding_uglmmutantsperm2<-mutantsperm2[match(allcorrect_uglm$Var1,mutantsperm2$chrom.pos),]
trueDenovos_uglm_mpandms<-rbind(corresponding_uglmmutantparent1, corresponding_uglmmutantparent2, corresponding_uglmmutantsperm1, corresponding_uglmmutantsperm2) #this is the dataset that has the mutantparent and mutantsperm samples at just the loci that appear in allcorrect_uglm

#somatic LoH:
LoH<-subset(metadatadf.0, GoH_or_LoH =="LoH")#subsets to just the LoH
trueLoHp1<-subset(LoH, MutationClass == "SomaticMutation" & sample==MutantParent1)#| sample==MutantParent2)) #subsets to just somatic mutations and just MutantParent1
trueLoHp1.1<-subset(LoH, MutationClass == "SomaticMutation" & sample==MutantParent2)#subsets to just somatic mutations and just MutantParent2
trueLoHp2<-subset(trueLoHp1, refdepth =="0" | altdepth=="0") #subsets to just the MutantParent1s that have zero minor allele frequency
trueLoHp2.1<-subset(trueLoHp1.1, refdepth =="0" | altdepth=="0") #subsets to just the MutantParent2s that have zero minor allele frequency
trueLoHmutantparents<- rbind(trueLoHp2, trueLoHp2.1) #combines trueLoHp2 and trueLoHp2.1
removeuniques_mp<-subset(trueLoHmutantparents,duplicated(chrom.pos)==TRUE | duplicated(chrom.pos, fromLast=TRUE)==TRUE) #subsets just the loci that appear in both MutantParent1 and MutantParent2 from trueLoHp2 and trueLoHp2.1


sperm1<-subset(LoH, MutationClass == "SomaticMutation" & sample ==MutantSperm1) #subsets just the LoH somatic mutations that are MutantSperm1
correspondingsperm1<-sperm1[match(removeuniques_mp$chrom.pos,sperm1$chrom.pos),] #subsets just the sperm1 loci that match loci in removeuniques_mp
sperm2<-subset(LoH, MutationClass == "SomaticMutation" & sample ==MutantSperm2) #subsets just the LoH somatic mutations that are MutantSperm2
correspondingsperm2<-sperm2[match(removeuniques_mp$chrom.pos,sperm2$chrom.pos),]#subsets just the sperm2 loci that match loci in removeuniques_mp
trueLoH_somatic_mpandms<-rbind(removeuniques_mp, unique(correspondingsperm1), unique(correspondingsperm2)) #correspondingsperm1 and 2 have duplicates of each line, use unique to remove duplicates, then rbind with removeuniques_mp for the dataset that has the mutantparent and mutantsperm samples at just the loci that are true somatic LoH
#as.data.frame(table(trueLoH_somatic_mpandms$chrom.pos))
#nrow(subset(trueDenovos_somatic_mpandms, TrueorFalse=="True"& sample==MutantSperm1))

#uglm true LoH:
trueLoHsperm1<- subset(LoH, MutationClass == "UniqueGermlineMutation" & sample==MutantSperm1) #subset just the LoH uglms that are MutantSperm1
trueLoHsperm1.1<- subset(LoH, MutationClass == "UniqueGermlineMutation" & sample==MutantSperm2) #subset just the LoH uglms that are MutantSperm2

trueLoHsperm2<- subset(trueLoHsperm1, refdepth =="0" | altdepth=="0") #subsets to just the trueLoHsperm1s that have zero minor allele frequency
trueLoHsperm2.1<- subset(trueLoHsperm1.1, refdepth =="0" | altdepth=="0") #subsets to just the trueLoHsperm1.1s that have zero minor allele frequency
trueLoHmutantsperm<- rbind(trueLoHsperm2, trueLoHsperm2.1) #combines trueLoHsperm2 and trueLoHsperm2.1
removeuniques_ms<-subset(trueLoHmutantsperm,duplicated(chrom.pos)==TRUE | duplicated(chrom.pos, fromLast=TRUE)==TRUE) #subsets just the loci that appear in both MutantSperm1 and MutantSperm2 from trueLoHsperm2 and trueLoHsperm2.1

mp1<-subset(LoH, MutationClass == "UniqueGermlineMutation" & sample==MutantParent1) #subsets just the LoH uglms that are MutantParent1
correspondingmp1<- mp1[match(removeuniques_ms$chrom.pos, mp1$chrom.pos),] #subsets just the mp1 loci that match loci in removeuniques_ms
mp2<-subset(LoH, MutationClass == "UniqueGermlineMutation" & sample==MutantParent2) #subsets just the LoH uglms that are MutantParent2
correspondingmp2<- mp2[match(removeuniques_ms$chrom.pos, mp2$chrom.pos),] #subsets just the mp2 loci that match loci in removeuniques_ms
trueLoH_uglm_mpandms<-rbind(removeuniques_ms, unique(correspondingmp1), unique(correspondingmp2)) #correspondingmp1 and 2 have duplicates of each line, use unique to remove duplicates, then rbind with removeuniques_ms for the dataset that has the mutantparent and mutantsperm samples at just the loci that are true uglm LoH

#trueLoHp1<-subset(trueLoHp, sample=="mutparent1" | sample=="mutparent2")
#trueLoHp2<-trueLoHp1[trueLoHp1$chrom.pos %in% names(which(table(trueLoHp1$chrom.pos) > 1)), ]
#sperm1<-subset(metadatadf.0, sample ==MutantSperm1)
#sperm2<-subset(metadatadf.0, sample ==MutantSperm2)
#trueLoHsperm<-sperm[match(trueLoHp2$chrom.pos, sperm$chrom.pos),]
#trueLoHsperm1<-subset(sperm1, GoH_or_LoH =="LoH")
#trueLoHsperm2<-subset(sperm2, GoH_or_LoH =="LoH")
#trueLoHp2<-subset(trueLoHp, sample=="mutparent2")

globalglm<-subset(metadatadf.0, MutationClass=="GlobalGermlineMutation")
#gglm true denovos:
justparents<-subset(globalglm, startsWith(globalglm$sample, 'CAP')==TRUE) #subsets to just the parent samples
justsperm<-subset(globalglm, startsWith(globalglm$sample, 'CAS')==TRUE)
justparents_freq<-as.data.frame(table(justparents$chrom.pos)) #gives how many parent samples are present per locus
trueDeNovos_gglm<-subset(justparents, refdepth == 0 | altdepth ==0) #subsets just the parent samples that have either refdepth=0 or altdepth=0. this eliminates the nonmutants for which the putatively mutant allele is present in "nonmutants" at low frequency
trueDeNovos_gglm_freq<-as.data.frame(table(trueDeNovos_gglm$chrom.pos)) #gives  how many nonmutant parent samples are present per locus once they have been subsetted
justparents_freq_matched<- justparents_freq[match(trueDeNovos_gglm_freq$Var1, justparents_freq$Var1),] #eliminates the loci that are not present at all in trueDeNovos_somatic_freq
total_gglm <- merge(trueDeNovos_gglm_freq, justparents_freq_matched,by="Var1") #now the number of samples from trueDeNovos_gglm_freq is Freq.x and justparents_freq_matched is Freq.y
allcorrect_gglm<-subset(total_gglm, Freq.x == Freq.y) #subsets just the loci where Freq.x and Freq.y are equal
ms1<-subset(globalglm, sample=="MutantSperm1")# #pulls out just the mutantsperm1s
correspondingms1<-justsperm[match(allcorrect_gglm$Var1,justsperm$chrom.pos),] #subsets just a single sample at the locus that appears in allcorrect_gglm

#gglm true LoH:
justsperm<-subset(globalglm, startsWith(globalglm$sample, 'CAS')==TRUE)
trueLOHgglm = rep("A", nrow(justsperm))
for (i in 1:nrow(justsperm)){
  #print(i)
  if (justsperm$refdepth[i] == 0 | justsperm$altdepth[i] ==0) {
    trueLOHgglm[i] = justsperm$chrom.pos[i]
    #print(justsperm$refdepth[i])
  } else {
    trueLOHgglm[i] = "nottrueLOH"
  }
}
tlgdf<-data.frame("chrom.pos"=trueLOHgglm)

tlgdf_freq<-as.data.frame(table(tlgdf))
trueLOHgglm_freq<-data.frame("Var1"=tlgdf_freq$tlgdf, "Freq"=tlgdf_freq$Freq)
justsperm_freq<-as.data.frame(table(justsperm$chrom.pos))
justsperm_freq_matched<- justsperm_freq[match(trueLOHgglm_freq$Var1, justsperm_freq$Var1),]
total_LOH_gglm <- merge(trueLOHgglm_freq, justsperm_freq_matched,by="Var1") #now the number of samples from trueDeNovos_gglm_freq is Freq.x and justparents_freq_matched is Freq.y
allcorrect_LOH_gglm<-subset(total_LOH_gglm, Freq.x == Freq.y) #subsets just the loci where Freq.x and Freq.y are equal

correspondingms1_LOH<-justsperm[match(allcorrect_LOH_gglm$Var1,justsperm$chrom.pos),] #subsets just a single sample at the locus that appears in allcorrect_gglm

#CA65somatic:
CA65somatic<-subset(metadatadf.0, ColonyName=="CAcolony65" & MutationClass=="SomaticMutation")
CA65somaticmpandms<-subset(CA65somatic, sample==MutantParent1 | sample==MutantParent2 | sample==MutantSperm1 | sample==MutantSperm2)

#CA65uglm:
CA65uglm<-subset(metadatadf.0, ColonyName=="CAcolony65" & MutationClass=="UniqueGermlineMutation")
CA65uglmmpandms<-subset(CA65uglm, sample==MutantParent1 | sample==MutantParent2 | sample==MutantSperm1 | sample==MutantSperm2)

#metadatadf<-rbind( trueDeNovos_somatic, trueDeNovos_uglmmetadatadf, trueLoHp2, trueLoHsperm2,globalglm)#, trueLoHp2)
somatic_metadatadf.0<-subset(metadatadf.0, MutationClass == "SomaticMutation")
uniquesomatic_metadatadf.0<-somatic_metadatadf.0[match(unique(somatic_metadatadf.0$chrom.pos), somatic_metadatadf.0$chrom.pos),]

somaticmetadatadf<-rbind(trueDenovos_somatic_mpandms, trueLoH_somatic_mpandms)#, CA65somaticmpandms) #combines all of the filtered somatic datasets into one
somaticmetadatadf_denovo<-subset(somaticmetadatadf, GoH_or_LoH=="DeNovo")
somaticmetadatadf_loh<-subset(somaticmetadatadf, GoH_or_LoH=="LoH")
uniquesomaticmetadatadf<-somaticmetadatadf[match(unique(somaticmetadatadf$chrom.pos), 					somaticmetadatadf$chrom.pos),]
uniquesomaticmetadatadf_denovo<-subset(uniquesomaticmetadatadf, GoH_or_LoH=="DeNovo")
uniquesomaticmetadatadf_loh<-subset(uniquesomaticmetadatadf, GoH_or_LoH=="LoH")

mutantsperms1<-subset(somaticmetadatadf,sample==MutantSperm1)
mutantsperms2<-subset(somaticmetadatadf,sample==MutantSperm2)
mutantsperm1df<-data.frame("chrom.pos"=mutantsperms1$chrom.pos, "mutantalleledepth"=mutantsperms1$mutant_allele_depth, "GoH_or_LoH"=mutantsperms1$GoH_or_LoH)
mutantsperm2df<-data.frame("chrom.pos"=mutantsperms2$chrom.pos, "mutantalleledepth"=mutantsperms2$mutant_allele_depth, "GoH_or_LoH"=mutantsperms2$GoH_or_LoH)

Inheritance = rep("A", nrow(mutantsperm1df))
Inherited = rep("inherited",nrow(mutantsperm1df))
NotInherited = rep("notinherited",nrow(mutantsperm1df))
Inheritancedf<-as.data.frame(Inheritance)
#Inheriteddf<-as.data.frame(Inherited)
#NotInheriteddf<-as.data.frame(NotInherited)

for (i in 1:nrow(mutantsperm1df)){
  if (mutantsperm1df$GoH_or_LoH[i] =="DeNovo" & mutantsperm1df$mutantalleledepth[i] > 0 & mutantsperm2df$mutantalleledepth[i] > 0) {
    Inheritance[i] = Inherited[i]
    #print("inherited")
  } else if (mutantsperm1df$GoH_or_LoH[i] =="DeNovo" & mutantsperm1df$mutantalleledepth[i] == 0 | mutantsperm2df$mutantalleledepth[i] == 0) {
    Inheritance[i] = NotInherited[i]
    #print("notinherited")
  } else if (mutantsperm1df$GoH_or_LoH[i] =="LoH" & mutantsperm1df$mutantalleledepth[i] == 1 & mutantsperm2df$mutantalleledepth[i] == 1) {
    Inheritance[i] = Inherited[i]
  } else if (mutantsperm1df$GoH_or_LoH[i] =="LoH" & mutantsperm1df$mutantalleledepth[i] != 1 | mutantsperm2df$mutantalleledepth[i] != 1) {
    Inheritance[i] = NotInherited[i]
  } else (Inheritance[i] =="na")
}
 
  
inheriteddenovosomaticmetadatadf<-subset(somaticmetadatadf, TrueorFalse=="True" & GoH_or_LoH=="DeNovo") #this will be definition exlude all CA65 samples

pairwisesomatic<-subset(metadatadf.0, MutationClass=="SomaticMutation")
pairwisesomaticmpandms<-subset(pairwisesomatic, sample==MutantParent1 | sample==MutantParent2 | sample==MutantSperm1 | sample==MutantSperm2)
somatic6and8<-subset(pairwisesomaticmpandms, startsWith(pairwisesomaticmpandms$sample,'CAP6')==TRUE |
                       startsWith(pairwisesomaticmpandms$sample,'CAP8')==TRUE |
                       startsWith(pairwisesomaticmpandms$sample,'CAS8')==TRUE |
                       startsWith(pairwisesomaticmpandms$sample,'CAS6')==TRUE)

somatic22and23<-subset(pairwisesomaticmpandms, startsWith(pairwisesomaticmpandms$sample,'CAP22')==TRUE |
                         startsWith(pairwisesomaticmpandms$sample,'CAP23')==TRUE |
                         startsWith(pairwisesomaticmpandms$sample,'CAS22')==TRUE |
                         startsWith(pairwisesomaticmpandms$sample,'CAS23')==TRUE)
somatic22and24<-subset(pairwisesomaticmpandms, startsWith(pairwisesomaticmpandms$sample,'CAP22')==TRUE |
                         startsWith(pairwisesomaticmpandms$sample,'CAP24')==TRUE |
                         startsWith(pairwisesomaticmpandms$sample,'CAS22')==TRUE |
                         startsWith(pairwisesomaticmpandms$sample,'CAS24')==TRUE)
somatic22and26<-subset(pairwisesomaticmpandms, startsWith(pairwisesomaticmpandms$sample,'CAP22')==TRUE |
                         startsWith(pairwisesomaticmpandms$sample,'CAP26')==TRUE |
                         startsWith(pairwisesomaticmpandms$sample,'CAS22')==TRUE |
                         startsWith(pairwisesomaticmpandms$sample,'CAS26')==TRUE)
somatic23and26<-subset(pairwisesomaticmpandms, startsWith(pairwisesomaticmpandms$sample,'CAP23')==TRUE |
                         startsWith(pairwisesomaticmpandms$sample,'CAP26')==TRUE |
                         startsWith(pairwisesomaticmpandms$sample,'CAS23')==TRUE |
                         startsWith(pairwisesomaticmpandms$sample,'CAS26')==TRUE)
somatic24and26<-subset(pairwisesomaticmpandms, startsWith(pairwisesomaticmpandms$sample,'CAP24')==TRUE |
                         startsWith(pairwisesomaticmpandms$sample,'CAP26')==TRUE |
                         startsWith(pairwisesomaticmpandms$sample,'CAS24')==TRUE |
                         startsWith(pairwisesomaticmpandms$sample,'CAS26')==TRUE)
somatic23and24<-subset(pairwisesomaticmpandms, startsWith(pairwisesomaticmpandms$sample,'CAP23')==TRUE |
                         startsWith(pairwisesomaticmpandms$sample,'CAP24')==TRUE |
                         startsWith(pairwisesomaticmpandms$sample,'CAS23')==TRUE |
                         startsWith(pairwisesomaticmpandms$sample,'CAS24')==TRUE)
somaticCA56<-subset(somaticmetadatadf, startsWith(somaticmetadatadf$sample,'CAP12')==TRUE | 
  startsWith(somaticmetadatadf$sample,'CAS12')==TRUE |
  startsWith(somaticmetadatadf$sample,'CAP6')==TRUE | 
  startsWith(somaticmetadatadf$sample,'CAS6')==TRUE |
  startsWith(somaticmetadatadf$sample,'CAP8')==TRUE |
  startsWith(somaticmetadatadf$sample,'CAS8')==TRUE)
uniquesomaticCA56<-somaticCA56[match(unique(somaticCA56$chrom.pos), 					somaticCA56$chrom.pos),]
uniquesomaticCA56_denovo<-subset(uniquesomaticCA56, GoH_or_LoH=="DeNovo")
uniquesomaticCA56_loh<-subset(uniquesomaticCA56, GoH_or_LoH=="LoH")
somaticCA56_denovo<-subset(somaticCA56, GoH_or_LoH=="DeNovo")
somaticCA60<-subset(somaticmetadatadf, startsWith(somaticmetadatadf$sample,'CAP22')==TRUE | 
                      startsWith(somaticmetadatadf$sample,'CAS22')==TRUE |
                      startsWith(somaticmetadatadf$sample,'CAP23')==TRUE | 
                      startsWith(somaticmetadatadf$sample,'CAS23')==TRUE |
                      startsWith(somaticmetadatadf$sample,'CAP24')==TRUE |
                      startsWith(somaticmetadatadf$sample,'CAS24')==TRUE |
                      startsWith(somaticmetadatadf$sample,'CAP26')==TRUE |
                      startsWith(somaticmetadatadf$sample,'CAS26')==TRUE )
uniquesomaticCA60<-somaticCA60[match(unique(somaticCA60$chrom.pos), 					somaticCA60$chrom.pos),]
uniquesomaticCA60_denovo<-subset(uniquesomaticCA60, GoH_or_LoH=="DeNovo")
uniquesomaticCA60_loh<-subset(uniquesomaticCA60, GoH_or_LoH=="LoH")
somaticCA60_denovo<-subset(somaticCA60, GoH_or_LoH=="DeNovo")
somaticCA65<-subset(somaticmetadatadf, startsWith(somaticmetadatadf$sample,'CAP10')==TRUE | 
                      startsWith(somaticmetadatadf$sample,'CAS10')==TRUE |
                      startsWith(somaticmetadatadf$sample,'CAP11')==TRUE | 
                      startsWith(somaticmetadatadf$sample,'CAS11')==TRUE |
                      startsWith(somaticmetadatadf$sample,'CAP9')==TRUE |
                      startsWith(somaticmetadatadf$sample,'CAS9')==TRUE)
uniquesomaticCA65<-somaticCA65[match(unique(somaticCA65$chrom.pos), 					somaticCA65$chrom.pos),]
uniquesomaticCA65_denovo<-subset(uniquesomaticCA65, GoH_or_LoH=="DeNovo")
uniquesomaticCA65_loh<-subset(uniquesomaticCA65, GoH_or_LoH=="LoH")


uglmmetadatadf<-rbind(trueDenovos_uglm_mpandms, trueLoH_uglm_mpandms)#, CA65uglmmpandms)
unique_uglmmetadatadf<-uglmmetadatadf[match(unique(uglmmetadatadf$chrom.pos), 					uglmmetadatadf$chrom.pos),]
#USE THESE 5/17/2020
dndscvdf_uglmmetadatadf<- data.frame("sampleID"= unique_uglmmetadatadf$sample,"chr"= unique_uglmmetadatadf$chrom,"pos"= unique_uglmmetadatadf$pos, "ref" = unique_uglmmetadatadf$ref, "alt"= unique_uglmmetadatadf$alt) # do not use ref and alt!
#write.table(dndscvdf_uglmmetadatadf, file="~/Documents/GitHub/CoralGermline/dndscv/unique_uglmmetadatadf_onlychrs20200615_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)

uglmCA56<-subset(uglmmetadatadf, startsWith(uglmmetadatadf$sample,'CAP12')==TRUE | 
                   startsWith(uglmmetadatadf$sample,'CAS12')==TRUE |
                   startsWith(uglmmetadatadf$sample,'CAP6')==TRUE | 
                   startsWith(uglmmetadatadf$sample,'CAS6')==TRUE |
                   startsWith(uglmmetadatadf$sample,'CAP8')==TRUE |
                   startsWith(uglmmetadatadf$sample,'CAS8')==TRUE)
uniqueuglmCA56<-uglmCA56[match(unique(uglmCA56$chrom.pos), 					uglmCA56$chrom.pos),]
uglmCA56_denovo<-subset(uglmCA56, GoH_or_LoH=="DeNovo")
uglmCA60<-subset(uglmmetadatadf, startsWith(uglmmetadatadf$sample,'CAP22')==TRUE | 
                   startsWith(uglmmetadatadf$sample,'CAS22')==TRUE |
                   startsWith(uglmmetadatadf$sample,'CAP23')==TRUE | 
                   startsWith(uglmmetadatadf$sample,'CAS23')==TRUE |
                   startsWith(uglmmetadatadf$sample,'CAP24')==TRUE |
                   startsWith(uglmmetadatadf$sample,'CAS24')==TRUE |
                   startsWith(uglmmetadatadf$sample,'CAP26')==TRUE |
                   startsWith(uglmmetadatadf$sample,'CAS26')==TRUE )
uniqueuglmCA60<-uglmCA60[match(unique(uglmCA60$chrom.pos), 					uglmCA60$chrom.pos),]

uglmCA65<-subset(uglmmetadatadf, startsWith(uglmmetadatadf$sample,'CAP10')==TRUE | 
                   startsWith(uglmmetadatadf$sample,'CAS10')==TRUE |
                   startsWith(uglmmetadatadf$sample,'CAP11')==TRUE | 
                   startsWith(uglmmetadatadf$sample,'CAS11')==TRUE |
                   startsWith(uglmmetadatadf$sample,'CAP9')==TRUE |
                   startsWith(uglmmetadatadf$sample,'CAS9')==TRUE) 	
uniqueuglmCA65<-uglmCA65[match(unique(uglmCA65$chrom.pos), 					uglmCA65$chrom.pos),]

gglmmetadatadf<- rbind(correspondingms1, correspondingms1_LOH)
unique_gglmmetadatadf<-gglmmetadatadf[match(unique(gglmmetadatadf$chrom.pos), 					gglmmetadatadf$chrom.pos),]
#USE THESE 5/17/2020
dndscvdf_gglmmetadatadf<- data.frame("sampleID"= unique_gglmmetadatadf$sample,"chr"= unique_gglmmetadatadf$chrom,"pos"= unique_gglmmetadatadf$pos, "ref" = unique_gglmmetadatadf$ref, "alt"= unique_gglmmetadatadf$alt) # do not use ref and alt!
#write.table(dndscvdf_gglmmetadatadf, file="~/Documents/GitHub/CoralGermline/dndscv/unique_gglmmetadatadf_onlychrs20200615_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)

gglmCA56<-subset(gglmmetadatadf, startsWith(gglmmetadatadf$sample,'CAP12')==TRUE | 
                   startsWith(gglmmetadatadf$sample,'CAS12')==TRUE |
                   startsWith(gglmmetadatadf$sample,'CAP6')==TRUE | 
                   startsWith(gglmmetadatadf$sample,'CAS6')==TRUE |
                   startsWith(gglmmetadatadf$sample,'CAP8')==TRUE |
                   startsWith(gglmmetadatadf$sample,'CAS8')==TRUE)
gglmCA56_denovo<-subset(gglmCA56, GoH_or_LoH=="DeNovo")
gglmCA60<-subset(gglmmetadatadf, startsWith(gglmmetadatadf$sample,'CAP22')==TRUE | 
                   startsWith(gglmmetadatadf$sample,'CAS22')==TRUE |
                   startsWith(gglmmetadatadf$sample,'CAP23')==TRUE | 
                   startsWith(gglmmetadatadf$sample,'CAS23')==TRUE |
                   startsWith(gglmmetadatadf$sample,'CAP24')==TRUE |
                   startsWith(gglmmetadatadf$sample,'CAS24')==TRUE |
                   startsWith(gglmmetadatadf$sample,'CAP26')==TRUE |
                   startsWith(gglmmetadatadf$sample,'CAS26')==TRUE )
gglmCA65<-subset(gglmmetadatadf, startsWith(gglmmetadatadf$sample,'CAP10')==TRUE | 
                   startsWith(gglmmetadatadf$sample,'CAS10')==TRUE |
                   startsWith(gglmmetadatadf$sample,'CAP11')==TRUE | 
                   startsWith(gglmmetadatadf$sample,'CAS11')==TRUE |
                   startsWith(gglmmetadatadf$sample,'CAP9')==TRUE |
                   startsWith(gglmmetadatadf$sample,'CAS9')==TRUE) 	

colony<-c(rep("CA56" , 4) , rep("CA60" , 4) , rep("CA65" , 4))
class<-factor(rep(c("Non-inherited Mosaic", "Inherited Mosaic", "Unique Germline" , "Global Germline") , 3),levels = c("Non-inherited Mosaic", "Inherited Mosaic","Unique Germline","Global Germline"))
value<-c(nrow(a_notinherited), nrow(a_inherited), nrow(uniqueuglmCA56), nrow(gglmCA56), nrow(b_notinherited), nrow(b_inherited), nrow(uniqueuglmCA60), nrow(gglmCA65), nrow(c_notinherited), nrow(c_inherited), nrow(uniqueuglmCA65), nrow(gglmCA65))

colonyprop<-c(rep("CAcolony56" , 3) , rep("CAcolony60" , 3) , rep("CAcolony65" , 3))
classprop<-factor(rep(c("Inherited Mosaic", "Unique Germline" , "Global Germline") , 3),levels = c("Inherited Mosaic","Unique Germline","Global Germline"))
proportions<-c(nrow(a_inherited)/(nrow(a_inherited)+nrow(uniqueuglmCA56)+nrow(gglmCA56)), nrow(uniqueuglmCA56)/(nrow(a_inherited)+nrow(uniqueuglmCA56)+nrow(gglmCA56)), nrow(gglmCA56)/(nrow(a_inherited)+nrow(uniqueuglmCA56)+nrow(gglmCA56)), nrow(b_inherited)/(nrow(b_inherited)+nrow(uniqueuglmCA60)+nrow(gglmCA60)), nrow(uniqueuglmCA60)/(nrow(b_inherited)+nrow(uniqueuglmCA60)+nrow(gglmCA60)), nrow(gglmCA60)/(nrow(b_inherited)+nrow(uniqueuglmCA60)+nrow(gglmCA60)), nrow(c_inherited)/(nrow(c_inherited)+nrow(uniqueuglmCA65)+nrow(gglmCA65)), nrow(uniqueuglmCA65)/(nrow(c_inherited)+nrow(uniqueuglmCA65)+nrow(gglmCA65)), nrow(gglmCA65)/(nrow(c_inherited)+nrow(uniqueuglmCA65)+nrow(gglmCA65)))
prop100<-100*proportions
mutationdataprop<-data.frame(colonyprop, classprop, prop100)

#value<-c(nrow(a_notinherited), nrow(a_inherited), nrow(b_notinherited), nrow(b_inherited),nrow(c_notinherited), nrow(c_inherited))
mutationdata<-data.frame(colony, class, value)
comparisonplot<-ggplot(mutationdata, aes(x=colony, y=value, fill=class)) + 
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  scale_fill_brewer(palette = "Paired") +
  ylab("Number of mutations") + xlab("")+
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25))+
  theme(legend.text=element_text(size=25))
  #+ 
  #scale_x_discrete(limits=mutationdata$class)
comparisonplotprop<-ggplot(mutationdataprop, aes(x=colonyprop, y=prop100, fill=classprop)) + 
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  scale_fill_brewer(palette = "Paired") +
  ylab("Percent of sperm mutations") + xlab("")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15))

metadatadf<- rbind(somaticmetadatadf, uglmmetadatadf, gglmmetadatadf)
#metadatadf<-rbind( DeNovos, trueLoHp1)#2,trueLoHsperm)#, trueLoHp2)
#write.table(metadatadf, file="CAcolony60_CAP22-23-24muts_20191125.txt",sep="\t",quote=FALSE, row.name=FALSE)

uniquemetadatadf<- metadatadf[match(unique(metadatadf$chrom.pos), 					metadatadf$chrom.pos),]
inheriteddenovo<-subset(metadatadf, SpermAverage>0.0 & ParentAverage>0.1 & ParentAverage<1)
inheritedloh<-subset(metadatadf, SpermAverage==1 & ParentAverage==1)
inherited<-rbind(inheriteddenovo, inheritedloh)
unique_inherited<-inherited[match(unique(inherited$chrom.pos), 					inherited$chrom.pos),]
unique_notinherited<-subset(uniquemetadatadf,TrueorFalse=="False")
unique_somatic<-subset(uniquemetadatadf, MutationClass=="SomaticMutation")
unique_uglmmetadatadf<-subset(uniquemetadatadf, MutationClass=="UniqueGermlineMutation")
unique_globalglm<-subset(uniquemetadatadf, MutationClass=="GlobalGermlineMutation")

dndscvdf_somatic<- data.frame("sampleID"= uniquesomaticmetadatadf$sample,"chr"= uniquesomaticmetadatadf$chrom,"pos"= uniquesomaticmetadatadf$pos, "ref" = uniquesomaticmetadatadf$ref, "alt"= uniquesomaticmetadatadf$alt) # do not use ref and alt!
#write.table(dndscvdf_somatic, file="~/Documents/GitHub/CoralGermline/dndscv/uniquesomaticmetadatadf_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)

dndscvdf_somatic_denovo<- data.frame("sampleID"= uniquesomaticmetadatadf_denovo$sample,"chr"= uniquesomaticmetadatadf_denovo$chrom,"pos"= uniquesomaticmetadatadf_denovo$pos, "ref" = uniquesomaticmetadatadf_denovo$ref, "alt"= uniquesomaticmetadatadf_denovo$alt) # do not use ref and alt!
#write.table(dndscvdf_somatic_denovo, file="~/Documents/GitHub/CoralGermline/dndscv/uniquesomaticmetadatadf_denovo_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)

dndscvdf_somatic_loh<- data.frame("sampleID"= uniquesomaticmetadatadf_loh$sample,"chr"= uniquesomaticmetadatadf_loh$chrom,"pos"= uniquesomaticmetadatadf_loh$pos, "ref" = uniquesomaticmetadatadf_loh$ref, "alt"= uniquesomaticmetadatadf_loh$alt) # do not use ref and alt!
#write.table(dndscvdf_somatic_loh, file="~/Documents/GitHub/CoralGermline/dndscv/uniquesomaticmetadatadf_loh_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)

dndscvdf_inherited<- data.frame("sampleID"= unique_inherited$sample,"chr"= unique_inherited$chrom,"pos"= unique_inherited$pos, "ref" = unique_inherited$ref, "alt"= unique_inherited$alt) # do not use ref and alt!
dndscvdf_notinherited<- data.frame("sampleID"= unique_notinherited$sample,"chr"= unique_notinherited$chrom,"pos"= unique_notinherited$pos, "ref" = unique_notinherited$ref, "alt"= unique_notinherited$alt) # do not use ref and alt!
#dndscvdf_somatic<- data.frame("sampleID"= unique_somatic$sample,"chr"= unique_somatic$chrom,"pos"= unique_somatic$pos, "ref" = unique_somatic$ref, "alt"= unique_somatic$alt) # do not use ref and alt!
dndscvdf_uglmmetadatadf<- data.frame("sampleID"= unique_uglmmetadatadf$sample,"chr"= unique_uglmmetadatadf$chrom,"pos"= unique_uglmmetadatadf$pos, "ref" = unique_uglmmetadatadf$ref, "alt"= unique_uglmmetadatadf$alt) # do not use ref and alt!
dndscvdf_globalglm<- data.frame("sampleID"= unique_globalglm$sample,"chr"= unique_globalglm$chrom,"pos"= unique_globalglm$pos, "ref" = unique_globalglm$ref, "alt"= unique_globalglm$alt) # do not use ref and alt!

#write.table(dndscvdf, file="~/Documents/GitHub/CoralGermline/dndscv/CAcolony65CAP11_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)
#write.table(dndscvdf_inherited, file="~/Documents/GitHub/CoralGermline/dndscv/CAcolony56and60and65_inherited_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)
#write.table(dndscvdf_notinherited, file="~/Documents/GitHub/CoralGermline/dndscv/CAcolony56and60and65_notinherited_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)
#write.table(dndscvdf_somatic, file="~/Documents/GitHub/CoralGermline/dndscv/CAcolony56and60and65_somatic_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)
#write.table(dndscvdf_uglmmetadatadf, file="~/Documents/GitHub/CoralGermline/dndscv/CAcolony56and60and65_uglmmetadatadf_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)
#write.table(dndscvdf_globalglm, file="~/Documents/GitHub/CoralGermline/dndscv/CAcolony56and60and65_globalglm_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)


#dndscvdf_mutparent1<- subset(dndscvdf,sampleID=="mutparent1")
#write.table(uniquemetadatadf, file="CAcolony60CAP22uniquemuts.txt",sep="\t",quote=FALSE, row.name=FALSE)

#x<-c(0,2,3,"a","d")
#y<-c(3, 2, 0, "a", "d")
#xy<-data.frame(x,y)
# for (i in 1:nrow(xy)){
#   if (xy$x[i] == xy$y[i]) {
#     jk = 2
#     print("yes")
#   } else {
#     print("no")
#     jk = 5
#   }  
# }
# print(mutant)


#for (i in 1:nrow(metadatadf.0)){
#  if (metadatadf.0$mutant_allele[i] == metadatadf.0$ref[i]) {
#    mutant_alleledepth = refdepth

#  } else mutant_alleledepth = altdepth
#}
scatterplot_func<-function(somaticmetadatadf) {
  #########comparison for all somatic mutations:##########
  mutantparentdf<-subset(somaticmetadatadf, sample==MutantParent1)# & TrueorFalse=="True") #for CAP22-1 and CAP22-2
  attach(mutantparentdf)
  mutantparentdf <- mutantparentdf[order(chrom.pos),]
  detach(mutantparentdf)
  
  #mutantparentdf_true<- mutantparentdf_true[match(uniquemetadatadf$chrom.pos, mutantparentdf_true$chrom.pos),]
  
  mutantparents2<-subset(somaticmetadatadf,sample==MutantParent2)# & TrueorFalse=="True")
  attach(mutantparents2)
  mutantparents2 <- mutantparents2[order(chrom.pos),]
  detach(mutantparents2)
  
  #mutantparents2_true<- mutantparents2_true[match(uniquemetadatadf$chrom.pos, mutantparents2_true$chrom.pos),]
  
  props1<-mutantparentdf$mutant_allele_depth/mutantparentdf$totaldepth.x
  props2<-mutantparents2$mutant_allele_depth/mutantparents2$totaldepth.x
  
  parentsaverage<-(props1+props2)/2

  spermdf1<-unique(subset(somaticmetadatadf, sample==MutantSperm1)) #removes duplicates that arise from the sperm singleton # & TrueorFalse=="True")
  attach(spermdf1)
  spermdf1 <- spermdf1[order(chrom.pos),]
  detach(spermdf1)
  
  spermdf2<-unique(subset(somaticmetadatadf, sample==MutantSperm2))# & TrueorFalse=="True")
  attach(spermdf2)
  spermdf2 <- spermdf2[order(chrom.pos),]
  detach(spermdf2)
  #spermdf_true<- spermdf_true[match(uniquemetadatadf$chrom.pos, spermdf_true$chrom.pos),]
  
  spermprops1<- spermdf1$mutant_allele_depth/spermdf1$totaldepth.x
  spermprops2<- spermdf2$mutant_allele_depth/spermdf2$totaldepth.x
  
  spermaverage<-(spermprops1+spermprops2)/2
  dfx<-cbind(spermdf1$chrom.pos, spermprops1,spermprops2,spermaverage)
  Inheritance = rep("A", nrow(spermdf1))
  LaxInheritance = rep("A",nrow(spermdf1))
  Inherited = rep("inherited",nrow(spermdf1))
  NotInherited = rep("notinherited",nrow(spermdf1))
  
  for (i in 1:nrow(spermdf1)){
    if (spermdf1$GoH_or_LoH[i] =="DeNovo" & spermaverage[i] >0) {
      Inheritance[i] = Inherited[i]
    } else if (spermdf1$GoH_or_LoH[i] =="DeNovo" & spermaverage[i] == 0) {
      Inheritance[i] = NotInherited[i]
    } else if (spermdf1$GoH_or_LoH[i] =="LoH" & spermdf1$mutant_allele_depth[i]/spermdf1$totaldepth.x[i] == 1 & spermdf2$mutant_allele_depth[i]/spermdf2$totaldepth.x[i] == 1) {
      Inheritance[i] = Inherited[i]
    } else if (spermdf1$GoH_or_LoH[i] =="LoH" & (spermdf1$mutant_allele_depth[i]/spermdf1$totaldepth.x[i] != 1 | spermdf2$mutant_allele_depth[i]/spermdf2$totaldepth.x[i] != 1)) {
      Inheritance[i] = NotInherited[i]
    } else (Inheritance[i] =="na")
  }
  
  df1<-data.frame(chrom=mutantparentdf$chrom, pos=mutantparentdf$pos, sample=mutantparentdf$sample, ref=mutantparentdf$ref, alt=mutantparentdf$alt, ParentAverage=parentsaverage, SpermAverage=spermaverage, TrueorFalse = as.factor(mutantparentdf$TrueorFalse), GoHorLoH = as.factor(mutantparentdf$GoH_or_LoH), MutationType = as.factor(mutantparentdf$MutationType), Inheritance = as.factor(Inheritance), LaxInheritance = as.factor(LaxInheritance), s1=spermdf1$mutant_allele_depth, s2=spermdf2$mutant_allele_depth)
  df<-subset(df1, ParentAverage>=0.1)
  df_denovo<-subset(df, ParentAverage<1)
  otherframe<-data.frame(chrom.pos=mutantparentdf$chrom.pos, sampleID=mutantparentdf$sample, ref=mutantparentdf$ref, alt=mutantparentdf$alt,  chrom=mutantparentdf$chrom, pos=mutantparentdf$pos, sample=mutantparentdf$sample, ParentAverage=parentsaverage, SpermAverage=spermaverage, s1=spermdf1$mutant_allele_depth, s2=spermdf2$mutant_allele_depth, TrueorFalse = as.factor(mutantparentdf$TrueorFalse), GoHorLoH = as.factor(mutantparentdf$GoH_or_LoH), MutationType = as.factor(mutantparentdf$MutationType),WhattoWhat=mutantparentdf$WhattoWhat)
  dfk<-subset(df_denovo, SpermAverage>0)
  #merged$TrueorFalse<-as.factor(merged$TrueorFalse)
  
  #parentspermcomparison<-lm(df_true$SpermAverage~df_true$ParentAverage)
  pall<-ggplot(df, aes(x=ParentAverage, y=SpermAverage, color=Inheritance)) + geom_point(size=2.5) +

    theme_bw() +
    ylim(0, 1) + xlim(0, 1) +
    scale_color_manual(values = c("red","black")) +
    ylab("Average Variant Allele Frequency in the Mutant Sperm Pool") + xlab("Average Variant Allele Frequency in the Mutant Parent")+
    theme(axis.text=element_text(size=25),
          axis.title=element_text(size=25))+
    theme(legend.text=element_text(size=25))
  
  p<-ggplot(df_denovo, aes(x=ParentAverage, y=SpermAverage, color=Inheritance)) + geom_point(size=2.5) +
    #geom_smooth(method=lm, aes(group=1)) +
    theme_bw() +
    #stat_cor(label.x=0.2,label.y=0.45, aes(group=1)) +
    ylim(0, .75) + xlim(0, 0.6) +
    scale_color_manual(values = c("red","black")) +
    ylab("Average Variant Allele Frequency in the Mutant Sperm Pool") + xlab("Average Variant Allele Frequency in the Mutant Parent")+
    theme(axis.text=element_text(size=25),
          axis.title=element_text(size=25))+
    theme(legend.text=element_text(size=25))
  
  pk<-ggplot(dfk, aes(x=ParentAverage, y=SpermAverage, color=Inheritance)) + geom_point(size=2.5) +
    geom_smooth(method=lm, aes(group=1)) +
    theme_bw() +
    stat_cor(label.x=0.2,label.y=0.45, aes(group=1), size=8) +
    ylim(0, .75) + xlim(0, .6) +
    scale_color_manual(values = c("red","black")) +
    ylab("Average Variant Allele Frequency in the Mutant Sperm Pool") + xlab("Average Variant Allele Frequency in the Mutant Parent")+
    theme(axis.text=element_text(size=20),
          axis.title=element_text(size=20))+
    geom_abline(intercept=0,slope=1)+
    theme(legend.text=element_text(size=20))
  
    #geom_text(x = 25, y = 300, label = lm_eqn(df_true), parse = TRUE)
  #return(p)
  #return(parentspermcomparison)
  #abline(parentspermcomparison, lwd=2)
  #cor<-cor.test(df_true$SpermAverage, df_true$ParentAverage)
  #return(summary(parentspermcomparison)$r.squared)
  newlist<-list(p, cor)
  #parentspermcomparison
  #return(newlist)
  #return(p | pk)
  #return(df)
  #return(otherframe)
  #return(pall)
  return(df)
  #dndscvdf<- data.frame("sampleID"= mutantparentdf$sample,"chr"= mutantparentdf$chrom,"pos"= mutantparentdf$pos, "ref" = mutantparentdf$ref, "alt"= uniquemetadatadf$alt) # do not use ref and alt!
  
}
somaticmetadatadfCODING<-subset(somaticmetadatadf, MutationType=="missense_variant" | MutationType=="synonymous_variant"| MutationType=="synonymous_variant")
allcolonies<-scatterplot_func(somaticmetadatadf) #shows significant positive correlation btwn parent avg and sperm avg
allcoloniesdndscv<-data.frame("sampleID"= allcolonies$sample,"chr"= allcolonies$chrom,"pos"= allcolonies$pos, "ref" = allcolonies$ref, "alt"= allcolonies$alt) # do not use ref and alt!

#write.table(allcoloniesdndscv, file="~/Documents/GitHub/CoralGermline/dndscv/allcolonies_onlychrs20200615_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)

a_dndscv<-data.frame("sampleID"= a_$sample,"chr"= a_$chrom,"pos"= a_$pos, "ref" = a_$ref, "alt"= a_$alt) # do not use ref and alt!
#write.table(a_dndscv, file="~/Documents/GitHub/CoralGermline/dndscv/CAcolony56somatic_onlychrs20200616_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)
b_dndscv<-data.frame("sampleID"= b_$sample,"chr"= b_$chrom,"pos"= b_$pos, "ref" = b_$ref, "alt"= b_$alt) # do not use ref and alt!
#write.table(b_dndscv, file="~/Documents/GitHub/CoralGermline/dndscv/CAcolony60somatic_onlychrs20200616_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)
c_dndscv<-data.frame("sampleID"= c_$sample,"chr"= c_$chrom,"pos"= c_$pos, "ref" = c_$ref, "alt"= c_$alt) # do not use ref and alt!
#write.table(c_dndscv, file="~/Documents/GitHub/CoralGermline/dndscv/CAcolony65somatic_onlychrs20200616_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)

allcoloniesCODING<-subset(allcolonies, MutationType=="missense_variant" | MutationType=="synonymous_variant"| MutationType=="synonymous_variant")
plotsomaticmetadatadfCODING<-scatterplot_func(somaticmetadatadfCODING)

allc2<-scatterplot_func(somaticmetadatadf)
allcolonies
s6and8<-scatterplot_func(somatic6and8)
s22and23<-scatterplot_func(somatic22and23)
s22and24
s22and26
s23and26
s23and24
s24and26

inheriteddenovoallcolonies<-subset(allcolonies, GoHorLoH=="DeNovo" & (s1>0 & s2>0))
inheritedlohallcolonies<-subset(allcolonies, GoHorLoH=="LoH" & SpermAverage==1 & ParentAverage==1)
inheritedallcolonies<-rbind(inheriteddenovoallcolonies, inheritedlohallcolonies)
inheritedallcoloniesCODING<-subset(inheritedallcolonies, MutationType=="missense_variant" | MutationType=="synonymous_variant"| MutationType=="synonymous_variant")

#uniqueinheritedallcolonies<-inheritedallcolonies[match(unique(inheritedallcolonies$chrom.pos), 					inheritedallcolonies$chrom.pos),]
inheritedcountallcolonies<-nrow(inheriteddenovoallcolonies)+nrow(inheritedlohallcolonies)
notinheriteddenovoallcolonies<-subset(allcolonies, GoHorLoH=="DeNovo" & (s1==0 | s2==0))

notinheritedlohallcolonies<-subset(allcolonies, GoHorLoH=="LoH" & SpermAverage<1)
notinheritedallcolonies<-rbind(notinheriteddenovoallcolonies, notinheritedlohallcolonies)
inheritedallcoloniesdndscv<-data.frame("sampleID"= inheritedallcolonies$sample,"chr"= inheritedallcolonies$chrom,"pos"= inheritedallcolonies$pos, "ref" = inheritedallcolonies$ref, "alt"= inheritedallcolonies$alt) # do not use ref and alt!
#write.table(inheritedallcoloniesdndscv, file="~/Documents/GitHub/CoralGermline/dndscv/inheritedallcolonies_onlychrs20200615_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)

inheritedlohallcoloniesdndscv<-data.frame("sampleID"= inheritedlohallcolonies$sample,"chr"= inheritedlohallcolonies$chrom,"pos"= inheritedlohallcolonies$pos, "ref" = inheritedlohallcolonies$ref, "alt"= inheritedlohallcolonies$alt) # do not use ref and alt!
#write.table(inheritedlohallcoloniesdndscv, file="~/Documents/GitHub/CoralGermline/dndscv/inheritedlohallcolonies_onlychrs20200615_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)

inheriteddenovoallcoloniesdndscv<-data.frame("sampleID"= inheriteddenovoallcolonies$sample,"chr"= inheriteddenovoallcolonies$chrom,"pos"= inheriteddenovoallcolonies$pos, "ref" = inheriteddenovoallcolonies$ref, "alt"= inheriteddenovoallcolonies$alt) # do not use ref and alt!
#write.table(inheriteddenovoallcoloniesdndscv, file="~/Documents/GitHub/CoralGermline/dndscv/inheriteddenovoallcolonies_onlychrs20200615_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)



notinheritedallcoloniesdndscv<-data.frame("sampleID"= notinheritedallcolonies$sample,"chr"= notinheritedallcolonies$chrom,"pos"= notinheritedallcolonies$pos, "ref" = notinheritedallcolonies$ref, "alt"= notinheritedallcolonies$alt) # do not use ref and alt!
#write.table(notinheritedallcoloniesdndscv, file="~/Documents/GitHub/CoralGermline/dndscv/notinheritedallcolonies_onlychrs20200615_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)

notinheritedlohallcoloniesdndscv<-data.frame("sampleID"= notinheritedlohallcolonies$sample,"chr"= notinheritedlohallcolonies$chrom,"pos"= notinheritedlohallcolonies$pos, "ref" = notinheritedlohallcolonies$ref, "alt"= notinheritedlohallcolonies$alt) # do not use ref and alt!
#write.table(notinheritedlohallcoloniesdndscv, file="~/Documents/GitHub/CoralGermline/dndscv/notinheritedlohallcolonies_onlychrs20200615_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)

notinheriteddenovoallcoloniesdndscv<-data.frame("sampleID"= notinheriteddenovoallcolonies$sample,"chr"= notinheriteddenovoallcolonies$chrom,"pos"= notinheriteddenovoallcolonies$pos, "ref" = notinheriteddenovoallcolonies$ref, "alt"= notinheriteddenovoallcolonies$alt) # do not use ref and alt!
#write.table(notinheriteddenovoallcoloniesdndscv, file="~/Documents/GitHub/CoralGermline/dndscv/notinheriteddenovoallcolonies_onlychrs20200615_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)

  
a_<-scatterplot_func(somaticCA56)
a_CODING<-subset(a_, MutationType=="missense_variant" | MutationType=="synonymous_variant"| MutationType=="synonymous_variant")
a_inherited<-subset(a_, Inheritance=="inherited")
a_notinherited<-subset(a_, Inheritance=="notinherited")
a_loh<-subset(a_, GoHorLoH=="LoH")
a_denovo<-subset(a_, GoHorLoH=="DeNovo")
stackedframe<-
stackedbarplot<-ggplot(a_, x=GoHorLoH, y=)

inheritedcount_func<-function(a_){
  

  inheriteddenovoa_<-subset(a_, GoHorLoH=="DeNovo" & (s1>0 & s2>0))
  inheritedloha_<-subset(a_, GoHorLoH=="LoH" & SpermAverage==1 & ParentAverage==1)
  notinheriteddenovoa_<-subset(a_, GoHorLoH=="DeNovo" & (s1==0 | s2==0))
  notinheritedloha_<-subset(a_, GoHorLoH=="LoH" & SpermAverage<1)
  inheritedcounta_<-nrow(inheriteddenovoa_)+nrow(inheritedloha_)
  inheritedpropa_<- inheritedcounta_/nrow(a_)
  return(inheritedpropa_)
}  
CAP6<-subset(a_, sample=="CAP6-1_S47") # when scatterplot_func returns "otherframe"
CAP8<-subset(a_, sample=="CAP8-1_S40")
CAP12<-subset(a_,sample=="CAP12-1_S46")
CAP6ic<-nrow(subset(CAP6, Inheritance=="inherited"))/nrow(CAP6)
CAP8ic<-nrow(subset(CAP8, Inheritance=="inherited"))/nrow(CAP8)
CAP12ic<-nrow(subset(CAP12, Inheritance=="inherited"))/nrow(CAP12)

b_<-scatterplot_func(somaticCA60)
b_CODING<-subset(b_, MutationType=="missense_variant" | MutationType=="synonymous_variant"| MutationType=="synonymous_variant")
b_loh<-subset(b_, GoHorLoH=="LoH")
b_denovo<-subset(b_, GoHorLoH=="DeNovo")

b_inherited<-subset(b_, Inheritance=="inherited")
b_notinherited<-subset(b_, Inheritance=="notinherited")
inheriteddenovob_<-subset(b_, SpermAverage>0 &ParentAverage>=0.1 &ParentAverage<1)
inheritedlohb_<-subset(b_, SpermAverage==1 & ParentAverage==1)
notinheriteddenovob_<-subset(b_, SpermAverage==0)
notinheritedlohb_<-subset(b_, ParentAverage==1 & SpermAverage<1)						  	
inheritedcountb_<-nrow(inheriteddenovob_)+nrow(inheritedlohb_)
CAP22<-subset(b_, sample=="CAP22-1_S43")
CAP23<-subset(b_, sample=="CAP23-1_S45")
CAP24<-subset(b_,sample=="CAP24-1_S42")
CAP26<-subset(b_,sample=="CAP26")
CAP22ic<-nrow(subset(CAP22, Inheritance=="inherited"))/nrow(CAP22)
CAP23ic<-nrow(subset(CAP23, Inheritance=="inherited"))/nrow(CAP23)
CAP24ic<-nrow(subset(CAP24, Inheritance=="inherited"))/nrow(CAP24)

c_<-scatterplot_func(somaticCA65)
c_CODING<-subset(c_, MutationType=="missense_variant" | MutationType=="synonymous_variant"| MutationType=="synonymous_variant")
c_loh<-subset(c_, GoHorLoH=="LoH")
c_denovo<-subset(c_, GoHorLoH=="DeNovo")

c_inherited<-subset(c_, Inheritance=="inherited")
c_notinherited<-subset(c_, Inheritance=="notinherited")

CAP10<-subset(c_, sample=="CAP10-1_S44")
CAP11<-subset(c_, sample=="CAP11-1_S41")
CAP9<-subset(c_,sample=="CAP9-1_S48")
CAP10ic<-nrow(subset(CAP10, Inheritance=="inherited"))/nrow(CAP10)
CAP11ic<-nrow(subset(CAP11, Inheritance=="inherited"))/nrow(CAP11)
CAP9ic<-nrow(subset(CAP9, Inheritance=="inherited"))/nrow(CAP9)
unfsom56<-subset(metadatadf.0, MutationClass=="SomaticMutation" & ColonyName=="CAcolony56")
unfsom56mpandms<-subset(unfsom56, sample==MutantParent1 | sample==MutantParent2 | sample==MutantSperm1 | sample==MutantSperm2)
unfsom60<-subset(metadatadf.0, MutationClass=="SomaticMutation" & ColonyName=="CAcolony60")
unfsom60mpandms<-subset(unfsom60, sample==MutantParent1 | sample==MutantParent2 | sample==MutantSperm1 | sample==MutantSperm2)

d_<-scatterplot_func(unfsom56mpandms)#subset(metadatadf.0, MutationClass=="SomaticMutation" & ColonyName=="CAcolony56")
e_<-scatterplot_func(unfsom60mpandms)
inheriteddenovoc_<-subset(c_, SpermAverage>0 &ParentAverage>=0.1 &ParentAverage<1)
inheritedlohc_<-subset(c_, SpermAverage==1 & ParentAverage==1)
notinheriteddenovoc_<-subset(c_, SpermAverage==0)
notinheritedlohc_<-subset(c_, ParentAverage==1 & SpermAverage<1)	
inheritedcountc_<-nrow(inheriteddenovoc_)+nrow(inheritedlohc_)
#barplot for inherited vs not inherited
colony<-c(rep("CA56" , 2) , rep("CA60" , 2) , rep("CA65" , 2))
class<-factor(rep(c("Not Inherited" , "Inherited") , 3),levels = c("Not Inherited","Inherited"))
value<-c(nrow(a_notinherited), nrow(a_inherited), nrow(b_notinherited), nrow(b_inherited),sum(CA65n), sum(CA65i))
valueifallwerenotinherited<-c((nrow(a_notinherited)+nrow(a_inherited)), 0, (nrow(b_notinherited)+ nrow(b_inherited)),0,(nrow(c_notinherited)+ nrow(c_inherited)),0)
valueifallwereinherited<-c(0,(nrow(a_notinherited)+nrow(a_inherited)), 0, (nrow(b_notinherited)+ nrow(b_inherited)),0,(nrow(c_notinherited)+ nrow(c_inherited)))

mutationdata<-data.frame(colony, class, value)

#per sample comparison of mutation types:
sample<-c("CAP6","CAP8","CAP12","CAP22","CAP23","CAP24","CAP26","CAP10","CAP11","CAP9")
colonies<-c(rep(1 , 4) , rep(2 , 4) , rep(3 , 4))

types<-factor(rep(c("Not Inherited Somatic" , "Inherited Somatic","Unique Germline","Global Germline") , 10),levels = c("Not Inherited","Inherited"))
numbers<-  
#GOHvsLoHinheritance:
colonies<-c(rep("CAcolony56" , 4) , rep("CAcolony60" , 4) , rep("CAcolony65" , 4))
inh<-rep(c("Inherited","Not inherited"),6)
counts<- c(nrow(subset(a_denovo,Inheritance=="inherited")), nrow(subset(a_denovo,Inheritance=="notinherited")),nrow(subset(a_loh,Inheritance=="inherited")),nrow(subset(a_loh,Inheritance=="notinherited")),
  nrow(subset(b_denovo,Inheritance=="inherited")), nrow(subset(b_denovo,Inheritance=="notinherited")),nrow(subset(b_loh,Inheritance=="inherited")),nrow(subset(b_loh,Inheritance=="notinherited")),
  nrow(subset(c_denovo,Inheritance=="inherited")), nrow(subset(c_denovo,Inheritance=="notinherited")),nrow(subset(c_loh,Inheritance=="inherited")),nrow(subset(c_loh,Inheritance=="notinherited")))

GorL<-rep(c("GoH","GoH","LoH","LoH"),3)
stackedframe<-data.frame("Colony"=colonies, "Inheritance"=inh, "Count"=counts, "GoHorLoH"=GorL)  
GOHframe<-subset(stackedframe, GoHorLoH=="GoH")
LOHframe<-subset(stackedframe, GoHorLoH=="LoH")
barwidth=0.35
stackedplot<-ggplot() + 
  geom_bar(data=GOHframe,
           mapping=aes(x=Colony, y=Count, fill=Inheritance),
           position="stack", 
           stat="identity",
           width = barwidth) +
  geom_bar(data=LOHframe,
             mapping=aes(x=Colony+ barwidth + 0.01, y=Count, fill=Inheritance),
             position="stack", 
             stat="identity",
             width = barwidth)#  


mutationdataifallwerenotinherited<-data.frame(colony,class, valueifallwerenotinherited)
mutationdataifallwereinherited<-data.frame(colony,class, valueifallwereinherited)
comparisonplotinh<-ggplot(mutationdata, aes(x=colony, y=value, fill=class)) + 
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  scale_fill_manual(values = c("black","red")) +
  ylab("Number of somatic mutations") + xlab("")+
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25))+
  theme(legend.text=element_text(size=25))
comparisonplotinhifallnotinh<-ggplot(mutationdataifallwerenotinherited, aes(x=colony, y=valueifallwerenotinherited, fill=class)) + 
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  scale_fill_manual(values = c("black","red")) +
  ylab("Number of somatic mutations") + xlab("")+
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25))+
  theme(legend.text=element_text(size=25))
comparisonplotinhifallinh<-ggplot(mutationdataifallwereinherited, aes(x=colony, y=valueifallwereinherited, fill=class)) + 
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  scale_fill_manual(values = c("black","red")) +
  ylab("Number of somatic mutations") + xlab("")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15))+
  theme(legend.text=element_text(size=15))
#boxplot for inherited proportion:
colony<-c(rep("CA56" , 3) , rep("CA60" , 3) , rep("CA65" , 2))
proportion<- 100*c(CAP6ic,CAP8ic, CAP12ic, CAP22ic, CAP23ic, CAP24ic, CAP11ic, CAP9ic)  
df_inheritedprops<-data.frame(colony,proportion)
box <- ggplot(df_inheritedprops, aes(x = colony, y = proportion))
box<- box + geom_boxplot()
box2 <- box + labs(x="", y="Percent of Somatic Mutations inherited by Sperm") + geom_boxplot(fill="red",
  position = position_dodge(0.9)
) +
  
  theme_bw()+
  #stat_compare_means(label.x=1,label.y=30,size=10)+
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25))

fit = lm(proportion ~ colony, df_inheritedprops)
fit
anova(fit)
oneway<-aov(proportion~colony)
summary(oneway)
plot(oneway)
pairwise.wilcox.test(df_inheritedprops$proportion, df_inheritedprops$colony,
                     p.adjust.method = "BH")
kruskal.test(proportion ~ colony, data = df_inheritedprops)

totala_<-inheritedcounta_+ nrow(uglmCA56)+ nrow(gglmCA56)
totalb_<- inheritedcountb_+ nrow(uglmCA60)+ nrow(gglmCA60)
totalc_<- inheritedcountc_+ nrow(uglmCA65)+ nrow(gglmCA65)
inheritedprop_func<-function(inheritedcounta_, a_){
  inheritedpropa_<- inheritedcounta_/nrow(a_)
  return(inheritedpropa_)
}
aprop<-inheritedprop_func(inheritedcounta_,a_)
bprop<-inheritedprop_func(inheritedcountb_,b_)
cprop<-inheritedprop_func(inheritedcountc_,c_)
  
#show what proportion of mutants in the sperm pool came from each type of mutation:
colony<-c(rep("CA56" , 3) , rep("CA60" , 3) , rep("CA65" , 3))
class<-factor(rep(c("Inherited Somatic", "Unique Germline", "Global Germline") , 3),levels = c("Inherited Somatic","Unique Germline","Global Germline"))
value<-c(inheritedcounta_*100/totala_, nrow(uglmCA56)*100/totala_, nrow(gglmCA56)*100/totala_, inheritedcountb_*100/totalb_, nrow(uglmCA60)*100/totalb_, nrow(gglmCA60)*100/totalb_, inheritedcountc_*100/totalc_, nrow(uglmCA65)*100/totalc_, nrow(gglmCA65)*100/totalc_)
mutationdata<-data.frame(colony, class, value)

proportionsplot<-ggplot(mutationdata, aes(x=colony, y=value, fill=class)) + 
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  scale_fill_brewer(palette = "Paired") +
  ylab("Percent of mutations seen in the sperm") + xlab("")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20))
(allcolonies | a) / (b | c)

unfilteredsomatic<-subset(metadatadf.0, MutationClass=="SomaticMutation")
allcolonies_unfiltered<-scatterplot_func(unfilteredsomatic)

allcolonies | allcolonies_unfiltered

CA56i<-c(nrow(subset(CAP6, Inheritance=="inherited")), nrow(subset(CAP8, Inheritance=="inherited")), nrow(subset(CAP12, Inheritance=="inherited")))
CA56ip<-CA56i/(CA56i+CA56uglm+CA56gglm)
CA60i<-c(nrow(subset(CAP22, Inheritance=="inherited")),nrow(subset(CAP23, Inheritance=="inherited")),nrow(subset(CAP24, Inheritance=="inherited")))
CA60ip<-CA60i/(CA60i+CA60uglm+CA60gglm)
CA65i<-c(nrow(subset(CAP11, Inheritance=="inherited")),nrow(subset(CAP9, Inheritance=="inherited")))
CA65ip<-CA65i/(CA65i+CA65uglm+CA65gglm)
CA56n<-c(nrow(subset(CAP6, Inheritance=="notinherited")), nrow(subset(CAP8, Inheritance=="notinherited")), nrow(subset(CAP12, Inheritance=="notinherited")))
CA56np<-CA56n/(CA56n+CA56uglm+CA56gglm)
CA60n<-c(nrow(subset(CAP22, Inheritance=="notinherited")),nrow(subset(CAP23, Inheritance=="notinherited")),nrow(subset(CAP24, Inheritance=="notinherited")))
CA60np<-CA60n/(CA60n+CA60uglm+CA60gglm)
CA65n<-c(nrow(subset(CAP11, Inheritance=="notinherited")),nrow(subset(CAP9, Inheritance=="notinherited")))
CA65np<-CA65n/(CA65n+CA65uglm+CA65gglm)
CA56uglm<-c(nrow(subset(uglmmetadatadf, startsWith(uglmmetadatadf$sample,'CAP12')==TRUE )) ,
nrow(subset(uglmmetadatadf, startsWith(uglmmetadatadf$sample,'CAP6')==TRUE ) ),
     nrow(subset(uglmmetadatadf, startsWith(uglmmetadatadf$sample,'CAP8')==TRUE )))
CA56uglmp<-CA56uglm/(CA56i+CA56uglm+CA56gglm)
CA60uglm<-c( nrow(subset(uglmmetadatadf, startsWith(uglmmetadatadf$sample,'CAP22')==TRUE ) ),
nrow(subset(uglmmetadatadf, startsWith(uglmmetadatadf$sample,'CAP23')==TRUE ) ),
nrow(subset(uglmmetadatadf, startsWith(uglmmetadatadf$sample,'CAP24')==TRUE )))#,
nrow(subset(uglmmetadatadf, startsWith(uglmmetadatadf$sample,'CAP26')==TRUE )))
CA60uglmp<-CA60uglm/(CA60i+CA60uglm+CA60gglm)
CA65uglm<- c(
nrow(subset(uglmmetadatadf, startsWith(uglmmetadatadf$sample,'CAP11')==TRUE ) ),
nrow(subset(uglmmetadatadf, startsWith(uglmmetadatadf$sample,'CAP9')==TRUE )))
CA65uglmp<-CA65uglm/(CA65i+CA65uglm+CA65gglm)
CA56gglm<-c(nrow(subset(gglmmetadatadf, startsWith(gglmmetadatadf$sample,'CAS12')==TRUE )) ,
            nrow(subset(gglmmetadatadf, startsWith(gglmmetadatadf$sample,'CAS6')==TRUE ) ),
            nrow(subset(gglmmetadatadf, startsWith(gglmmetadatadf$sample,'CAS8')==TRUE )))
CA56gglmp<-CA56gglm/(CA56i+CA56uglm+CA56gglm)
CA60gglm<-c( nrow(subset(gglmmetadatadf, startsWith(gglmmetadatadf$sample,'CAS22')==TRUE ) ),
             nrow(subset(gglmmetadatadf, startsWith(gglmmetadatadf$sample,'CAS23')==TRUE ) ),
             nrow(subset(gglmmetadatadf, startsWith(gglmmetadatadf$sample,'CAS24')==TRUE )))#,
             #nrow(subset(gglmmetadatadf, startsWith(gglmmetadatadf$sample,'CAS26')==TRUE )))
CA60gglmp<-CA60gglm/(CA60i+CA60uglm+CA60gglm)
CA65gglm<- c(
             nrow(subset(gglmmetadatadf, startsWith(gglmmetadatadf$sample,'CAS11')==TRUE ) ),
             nrow(subset(gglmmetadatadf, startsWith(gglmmetadatadf$sample,'CAS9')==TRUE )))
CA65gglmp<-CA65gglm/(CA65i+CA65uglm+CA65gglm)
sumswithnon<-c(sum(CA56n) ,sum(CA56i), sum(CA56uglm), sum(CA56gglm), sum(CA60n), sum(CA60i), sum(CA60uglm), sum(CA60gglm), sum(CA65n), sum(CA65i), sum(CA65uglm), sum(CA65gglm))
colonysums<-c(rep("CA56" , 4) , rep("CA60" , 4) , rep("CA65" , 4))
classsums<-factor(rep(c("Not Inherited Mosaic", "Inherited Mosaic", "Unique Germline" , "Global Germline") , 3),levels = c("Not Inherited Mosaic","Inherited Mosaic","Unique Germline","Global Germline"))
mutationsums<-data.frame(colonysums, classsums, sumswithnon)
z<- ggplot(mutationsums, aes(x=colonysums, y=sumswithnon,fill=classsums)) +
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  scale_fill_brewer(palette = "Paired") +
  ylab("Number of mutations") + xlab("")+
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25))   +
  theme(legend.text=element_text(size=25))
means<-100*c(mean(CA56ip), mean(CA56uglmp), mean(CA56gglmp), mean(CA60ip), mean(CA60uglmp), mean(CA60gglmp), mean(CA65ip), mean(CA65uglmp), mean(CA65gglmp))
se<-function(x) sd(x)/sqrt(length(x))
ses<-100*c(se(CA56ip), se(CA56uglmp), se(CA56gglmp), se(CA60ip), se(CA60uglmp), se(CA60gglmp), se(CA65ip), se(CA65uglmp), se(CA65gglmp))
colonyprop<-c(rep("CAcolony56" , 3) , rep("CAcolony60" , 3) , rep("CAcolony65" , 3))
classprop<-factor(rep(c("Inherited Mosaic", "Unique Germline" , "Global Germline") , 3),levels = c("Inherited Mosaic","Unique Germline","Global Germline"))
mutationprop<-data.frame(colonyprop, classprop, means,ses)
z<- ggplot(mutationprop, aes(x=colonyprop, y=means)) +
  geom_bar(stat="identity") +
  #ylim(0, 0.35)+
  theme_bw()+
  #ylim(0,0.6)+
  ylab("Inherited mutation percent") +
  
  theme(axis.text.x = element_text(angle = 90, size=18), axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+
  geom_errorbar(aes(ymin=means-(2*ses), ymax=means+(2*ses)))
zprop<-ggplot(mutationprop, aes(x=colonyprop, y=means, fill=classprop)) + 
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  scale_fill_manual(values = c("#1F78B4" ,"#B2DF8A", "#33A02C")) +
  ylab("Percent of mutations in each sperm pool") + xlab("")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15))       +
  geom_errorbar(position="dodge",aes(ymin=means-(2*ses), ymax=means+(2*ses)))
VAFvsmuttype_func<-function(somaticmetadatadf) {
  #########comparison for all somatic mutations:##########
  mutantparentdf<-subset(somaticmetadatadf, sample==MutantParent1)# & TrueorFalse=="True") #for CAP22-1 and CAP22-2
  #mutantparentdf_true<- mutantparentdf_true[match(uniquemetadatadf$chrom.pos, mutantparentdf_true$chrom.pos),]
  
  mutantparents2<-subset(somaticmetadatadf,sample==MutantParent2)# & TrueorFalse=="True")
  #mutantparents2_true<- mutantparents2_true[match(uniquemetadatadf$chrom.pos, mutantparents2_true$chrom.pos),]
  
  props1<-mutantparentdf$mutant_allele_depth/mutantparentdf$totaldepth.x
  props2<-mutantparents2$mutant_allele_depth/mutantparents2$totaldepth.x
  
  parentsaverage<-(props1+props2)/2
  
  spermdf1<-unique(subset(somaticmetadatadf, sample==MutantSperm1)) #removes duplicates that arise from the sperm singleton # & TrueorFalse=="True")
  spermdf2<-unique(subset(somaticmetadatadf, sample==MutantSperm2))# & TrueorFalse=="True")
  
  #spermdf_true<- spermdf_true[match(uniquemetadatadf$chrom.pos, spermdf_true$chrom.pos),]
  
  spermprops1<- spermdf1$mutant_allele_depth/spermdf1$totaldepth.x
  spermprops2<- spermdf2$mutant_allele_depth/spermdf2$totaldepth.x
  
  spermaverage<-(spermprops1+spermprops2)/2
  df<-data.frame(ParentAverage=parentsaverage, SpermAverage=spermaverage, TrueorFalse = as.factor(mutantparentdf$TrueorFalse), GoHorLoH = as.factor(mutantparentdf$GoH_or_LoH), MutationType = as.factor(mutantparentdf$MutationType))
  #df_true<-subset(df, GoHorLoH=="DeNovo")
  df_sub_missense<-subset(df, MutationType=="missense_variant")# | MutationType=="synonymous_variant" | MutationType=="nonsense_variant")
  df_sub_syn<-subset(df, MutationType=="synonymous_variant")
  df_sub_non<-subset(df, MutationType=="stop_gained")
  missense_count<-nrow(df_sub_missense)
  nonsense_count<-nrow(df_sub_non)
  synonymous_count<-nrow(df_sub_syn)
  missense_VAF<-mean(df_sub_missense$ParentAverage)
  syn_VAF<-mean(df_sub_syn$ParentAverage)
  non_VAF<-mean(df_sub_non$ParentAverage)
  ParentVAF<-c(syn_VAF, missense_VAF, non_VAF)
  MutType<-c("synonymous","missense","nonsense")
  standarddeviation<-c(sd(df_sub_syn$ParentAverage),sd(df_sub_missense$ParentAverage),sd(df_sub_non$ParentAverage))
  counts<-c(synonymous_count,missense_count, nonsense_count)
  allcounts<-synonymous_count+missense_count+ nonsense_count
  props<-c(synonymous_count/allcounts,missense_count/allcounts, nonsense_count/allcounts)
  pmdf<-data.frame(ParentVAF,MutType,standarddeviation,counts,props)
  pmdf$MutType<- factor(pmdf$MutType,levels=c("synonymous","missense","nonsense"))#,"April","May","June"))
  p<-ggplot(pmdf, aes(x=MutType, y=ParentVAF)) + 
    geom_bar(position="dodge", stat="identity") +
    theme_bw() +
    geom_errorbar(aes(ymin=ParentVAF-standarddeviation, ymax=ParentVAF+standarddeviation), width=.2,
                  position=position_dodge(.9)) +
    ylim(0,1) +
    theme(axis.text.x = element_text(size=15), axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))
  
  #j<-hist(df_sub_syn$ParentAverage)
  pcounts<-ggplot(pmdf, aes(x=MutType, y=counts)) + 
    geom_bar(position="dodge", stat="identity") +
    theme_bw() +
    ylim(0,100) + 
    theme(axis.text.x = element_text(size=18), axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))
  pprops<-ggplot(pmdf, aes(x=MutType, y=props)) + 
    geom_bar(position="dodge", stat="identity") +
    theme_bw() +
    ylim(0,1)
    
  #return(p)
  #return(pcounts)
  return(pcounts)

}  

all<-VAFvsmuttype_func(somaticmetadatadf) + ggtitle("All somatic mutations")
denovo<-VAFvsmuttype_func(somaticmetadatadf_denovo) + ggtitle("GoH somatic mutations")
loh<-VAFvsmuttype_func(somaticmetadatadf_loh) + ggtitle("LoH somatic mutations")
inhe<-VAFvsmuttype_func(inheritedallcolonies) + ggtitle("Inherited somatic mutations")
notinhe<- 
  
pprops<-(all | denovo | loh)
pcounts<-(all | denovo| loh)
p<-(all | denovo | loh)
pcounts/pprops/p/dndscvplot
CA60coding<-VAFvsmuttype_func(somaticCA60)
CA65coding<-VAFvsmuttype_func(somaticCA65)
CA56coding<-VAFvsmuttype_func(somaticCA56)
pcounts/p
allsites<-c(1,168551991,1) #for 56, 60, 65, respectively
nsamples<-c(3, 4, 3)
codingcount<-c(sum(CA56coding),sum(CA60coding), sum(CA65coding))
somaticcount<-c(nrow(uniquesomaticCA56), nrow(uniquesomaticCA60), nrow(uniquesomaticCA65))
somaticpernucpersample<-somaticcount/allsites/nsamples
codingpernucpersample<-codingcount/allsites/nsamples

somaticpernucpersample<-nrow()
mutantparentdf_true<-subset(somaticmetadatadf, sample==MutantParent1 & TrueorFalse=="True") #for CAP22-1 and CAP22-2
#mutantparentdf_true<- mutantparentdf_true[match(uniquemetadatadf$chrom.pos, mutantparentdf_true$chrom.pos),]

mutantparents2_true<-subset(somaticmetadatadf,sample==MutantParent2 & TrueorFalse=="True")
#mutantparents2_true<- mutantparents2_true[match(uniquemetadatadf$chrom.pos, mutantparents2_true$chrom.pos),]

props1_true<-mutantparentdf_true$mutant_allele_depth/mutantparentdf_true$totaldepth.x
props2_true<-mutantparents2_true$mutant_allele_depth/mutantparents2_true$totaldepth.x

newdf<- data.frame("chrom.pos"=mutantparentdf_true$chrom.pos, "GOHorLOH" = mutantparentdf_true$GoH_or_LoH, "Parent1prop" = props1_true, "Parent2prop" = props2_true)

parentsaverage_true<-(props1_true+props2_true)/2
parentsaverage_true<-na.omit(parentsaverage_true)


spermdf_true1<-subset(somaticmetadatadf, sample==MutantSperm1 & TrueorFalse=="True")
spermdf_true2<-subset(somaticmetadatadf, sample==MutantSperm2 & TrueorFalse=="True")

#spermdf_true<- spermdf_true[match(uniquemetadatadf$chrom.pos, spermdf_true$chrom.pos),]

spermprops_true1<- spermdf_true1$mutant_allele_depth/spermdf_true1$totaldepth.x
spermprops_true2<- spermdf_true2$mutant_allele_depth/spermdf_true2$totaldepth.x

spermaverage_true<-(spermprops_true1+spermprops_true2)/2
#spermprops_true<-na.omit(spermprops_true)

parentspermcomparison_true<-lm(spermaverage_true~parentsaverage_true)
plot(parentsaverage_true, spermaverage_true, col=ifelse(mutantparentdf_true$GoH_or_LoH=="DeNovo","red","black"))
abline(parentspermcomparison_true, lwd=2)
abline(0,1,lwd=4,col="red")
cor.test(parentsaverage_true, spermaverage_true)
summary(parentspermcomparison_true)$r.squared
##alt allele correlation when inheritance isFALSE:
mutantparentdf_false<-subset(metadatadf, sample=="mutparent1" & TrueorFalse=="False") #for CAP22-1 and CAP22-2
#mutantparentdf_false<- mutantparentdf_false[match(uniquemetadatadf$chrom.pos, mutantparentdf_false$chrom.pos),]

mutantparents2_false<-subset(metadatadf,sample=="mutparent2" & TrueorFalse=="False")
#mutantparents2_false<- mutantparents2_false[match(uniquemetadatadf$chrom.pos, mutantparents2_false$chrom.pos),]

props1_false<-mutantparentdf_false$mutant_allele_depth/mutantparentdf_false$totaldepth.x
props2_false<-mutantparents2_false$mutant_allele_depth/mutantparents2_false$totaldepth.x
#plot(props1, props2)

#propslm<-lm(props2~props1)

parentsaverage_false<-(props1_false+props2_false)/2
parentsaverage_false<-na.omit(parentsaverage_false)


spermdf_false<-subset(metadatadf, sample=="mutsperm" & TrueorFalse=="False")
#spermdf_false<- spermdf_false[match(uniquemetadatadf$chrom.pos, spermdf_false$chrom.pos),]

spermprops_false<- spermdf_false$mutant_allele_depth/spermdf_false$totaldepth.x
spermprops_false<-na.omit(spermprops_false)

parentspermcomparison_false<-lm(spermprops_false~parentsaverage_false)
##PLOT ALL OF THE MUTS!
#plot(parentsaverage_false, spermprops_false, col=ifelse(mutantparentdf_false$GoH_or_LoH=="DeNovo","green","white"))
plot(parentsaverage_false, spermprops_false, col=ifelse(mutantparentdf_false_GOH$MutationType=="missense_variant","black","green"))

#points(parentsaverage_true,spermprops_true,col=ifelse(mutantparentdf_true$GoH_or_LoH=="DeNovo","red","white"))
points(parentsaverage_true,spermprops_true,col=ifelse(mutantparentdf_false_GOH$MutationType=="missense_variant","black","green"))
###


abline(parentspermcomparison_false)
abline(parentspermcomparison_true,col="blue")
abline(0,1,lwd=3,col="red")
summary(parentspermcomparison_false)$r.squared
cor.test(parentsaverage_false,spermprops_false)

##comparing GOH and LOH when inheritance is false:
#LOH:
mutantparentdf_false_LOH<-subset(mutantparentdf_false, GoH_or_LoH=="LoH")
mutantparents2_false_LOH<-subset(mutantparents2_false, GoH_or_LoH=="LoH")

props1_false_LOH<-mutantparentdf_false_LOH$mutant_allele_depth/mutantparentdf_false_LOH$totaldepth.x
props2_false_LOH<-mutantparents2_false_LOH$mutant_allele_depth/mutantparents2_false_LOH$totaldepth.x

parentsaverage_false_LOH<-(props1_false_LOH+props2_false_LOH)/2
parentsaverage_false_LOH<-na.omit(parentsaverage_false_LOH)

spermdf_false_LOH<- subset(spermdf_false, GoH_or_LoH=="LoH")

spermprops_false_LOH<-spermdf_false_LOH$mutant_allele_depth/spermdf_false_LOH$totaldepth.x
spermprops_false_LOH<-na.omit(spermprops_false_LOH)

parentspermcomparison_false_LOH<-lm(spermprops_false_LOH~parentsaverage_false_LOH)
plot(parentsaverage_false_LOH, spermprops_false_LOH)
abline(parentspermcomparison_false_LOH)
summary(parentspermcomparison_false_LOH)$r.squared

#GOH:
mutantparentdf_false_GOH<-subset(mutantparentdf_false, GoH_or_LoH=="DeNovo")
mutantparents2_false_GOH<-subset(mutantparents2_false, GoH_or_LoH=="DeNovo")

props1_false_GOH<-mutantparentdf_false_GOH$mutant_allele_depth/mutantparentdf_false_GOH$totaldepth.x
props2_false_GOH<-mutantparents2_false_GOH$mutant_allele_depth/mutantparents2_false_GOH$totaldepth.x

parentsaverage_false_GOH<-(props1_false_GOH+props2_false_GOH)/2
parentsaverage_false_GOH<-na.omit(parentsaverage_false_GOH)

spermdf_false_GOH<- subset(spermdf_false, GoH_or_LoH=="DeNovo")

spermprops_false_GOH<-spermdf_false_GOH$mutant_allele_depth/spermdf_false_GOH$totaldepth.x
spermprops_false_GOH<-na.omit(spermprops_false_GOH)

parentspermcomparison_false_GOH<-lm(spermprops_false_GOH~parentsaverage_false_GOH)
plot(parentsaverage_false_GOH, spermprops_false_GOH,xlab="Proportion of alt allele in the parent",ylab="Proportion of alt allele in sperm",xlim=c(0,.5),ylim=c(0,0.5))


summary(parentspermcomparison_false_GOH)$r.squared  
cor.test(parentsaverage_false_LOH, spermprops_false_LOH)
cor.test(parentsaverage_false_GOH, spermprops_false_GOH)

##comparing GOH and LOH when inheritance is true:
#LOH:
mutantparentdf_true_LOH<-subset(mutantparentdf_true, GoH_or_LoH=="LoH")
mutantparents2_true_LOH<-subset(mutantparents2_true, GoH_or_LoH=="LoH")

props1_true_LOH<-mutantparentdf_true_LOH$mutant_allele_depth/mutantparentdf_true_LOH$totaldepth.x
props2_true_LOH<-mutantparents2_true_LOH$mutant_allele_depth/mutantparents2_true_LOH$totaldepth.x

parentsaverage_true_LOH<-(props1_true_LOH+props2_true_LOH)/2
parentsaverage_true_LOH<-na.omit(parentsaverage_true_LOH)


#spermdf_true<-subset(metadatadf.0, sample=="mutsperm" & TrueorFalse=="True")
#spermdf_true<- spermdf_true[match(uniquemetadatadf$chrom.pos, spermdf_true$chrom.pos),]

spermdf_true_LOH<- subset(spermdf_true, GoH_or_LoH=="LoH")

spermprops_true_LOH<-spermdf_true_LOH$mutant_allele_depth/spermdf_true_LOH$totaldepth.x
#spermprops_true_LOH<-na.omit(spermprops_true_LOH)

parentspermcomparison_true_LOH<-lm(spermprops_true_LOH~parentsaverage_true_LOH)
plot(parentsaverage_true_LOH, spermprops_true_LOH, col="pink")
abline(parentspermcomparison_true_LOH,col="pink")
summary(parentspermcomparison_true_LOH)$r.squared

#GOH:
mutantparentdf_true_GOH<-subset(mutantparentdf_true, GoH_or_LoH=="DeNovo")
mutantparents2_true_GOH<-subset(mutantparents2_true, GoH_or_LoH=="DeNovo")

props1_true_GOH<-mutantparentdf_true_GOH$mutant_allele_depth/mutantparentdf_true_GOH$totaldepth.x
props2_true_GOH<-mutantparents2_true_GOH$mutant_allele_depth/mutantparents2_true_GOH$totaldepth.x

parentsaverage_true_GOH<-(props1_true_GOH+props2_true_GOH)/2
parentsaverage_true_GOH<-na.omit(parentsaverage_true_GOH)


#spermdf_true<-subset(metadatadf.0, sample=="mutsperm" & TrueorFalse=="True")
#spermdf_true<- spermdf_true[match(uniquemetadatadf$chrom.pos, spermdf_true$chrom.pos),]

spermdf_true_GOH<- subset(spermdf_true, GoH_or_LoH=="DeNovo")

spermprops_true_GOH<-spermdf_true_GOH$mutant_allele_depth/spermdf_true_GOH$totaldepth.x
spermprops_true_GOH<-na.omit(spermprops_true_GOH)

parentspermcomparison_true_GOH<-lm(spermprops_true_GOH~parentsaverage_true_GOH)
plot(parentsaverage_true_GOH, spermprops_true_GOH,xlim=c(0,1),ylim=c(0,1))
points(parentsaverage_true_LOH, spermprops_true_LOH, col="blue")
abline(parentspermcomparison_true_GOH)
abline(parentspermcomparison_true_LOH,col="blue")
summary(parentspermcomparison_true_GOH)$r.squared  
cor.test(parentsaverage_true_LOH, spermprops_true_LOH)
cor.test(parentsaverage_true_GOH, spermprops_true_GOH)

plot(parentsaverage_false_GOH, spermprops_false_GOH, col=ifelse(mutantparentdf_false_GOH$MutationType=="missense_variant","red","black"), xlab="Proportion of mutant allele in the parent",ylab="Proportion of mutant allele in sperm",xlim=c(0,1),ylim=c(0,1))
points(parentsaverage_false_LOH, spermprops_false_LOH, col=ifelse(mutantparentdf_false_LOH$MutationType=="missense_variant","red","black"))
points(parentsaverage_true_GOH, spermprops_true_GOH, col=ifelse(mutantparentdf_true_GOH$MutationType=="missense_variant","red","black"))
points(parentsaverage_true_LOH,spermprops_true_LOH,col=ifelse(mutantparentdf_true_LOH$MutationType=="missense_variant","red","black"))
abline(parentspermcomparison_false_GOH,col="green")
abline(parentspermcomparison_false_LOH,col="blue")
abline(parentspermcomparison_true_GOH,col="red")
abline(parentspermcomparison_true_LOH,col="pink")
abline(0,1, lwd=3)
legend("bottomright",
       legend = c("Not inherited+GOH","Not inherited+LOH","Inherited+GOH","Inherited+LOH"),
       col=c("green", "blue","red","pink"),
       pch=c(16,16))

summary(parentspermcomparison_true_GOH)$r.squared  
summary(parentspermcomparison_true_LOH)$r.squared 
summary(parentspermcomparison_false_GOH)$r.squared 
summary(parentspermcomparison_false_LOH)$r.squared 
cor.test(parentsaverage_true_GOH, spermprops_true_GOH)
cor.test(parentsaverage_true_LOH, spermprops_true_LOH)
cor.test(parentsaverage_false_GOH, spermprops_false_GOH)
cor.test(parentsaverage_false_LOH, spermprops_false_LOH)
par(mfrow=c(1,1))
residtrueGOH<-resid(parentspermcomparison_true_GOH)
plot((mutantparents2_true_GOH$totaldepth.x + mutantparentdf_true_GOH$totaldepth.x)/2, residtrueGOH, xlim=c(0,500), ylim=c(-0.6,0.8), col="red",xlab="Total Parent Read Depth", ylab="Parent-Sperm Residuals", main = "Inherited +GOH")
points((mutantparents2_true_LOH$totaldepth.x + mutantparentdf_true_LOH$totaldepth.x)/2, residtrueLOH, col="pink")
points((mutantparents2_false_GOH$totaldepth.x + mutantparentdf_false_GOH$totaldepth.x)/2,residfalseGOH, col="green")
plot((mutantparents2_false_LOH$totaldepth.x + mutantparentdf_false_LOH$totaldepth.x)/2,residefalseLOH, col="blue",xlim=c(0,500), ylim=c(-0.6,0.5), xlab="Total Parent Read Depth", ylab="Parent-Sperm Residuals")
points(mutantparents2_true_GOH$totaldepth.x, residtrueGOH,col="red")
#abline(lm(residtrueGOH~mutantparents2_true_GOH$totaldepth.x))

residtrueLOH<-resid(parentspermcomparison_true_LOH)
plot(mutantparents2_true_LOH$totaldepth.x, residtrueLOH, xlim=c(0,500), ylim=c(-0.6,0.8), col="pink", main="Inherited+LOH",xlab="Total Parent Read Depth", ylab="Parent-Sperm Residuals")

residfalseGOH<-resid(parentspermcomparison_false_GOH)
plot(mutantparents2_false_GOH$totaldepth.x,residfalseGOH, col="green", xlim=c(0,500), ylim=c(-0.6,0.8), xlab="Total Parent Read Depth", ylab="Parent-Sperm Residuals", main="Not Inherited + GOH")

residefalseLOH<-resid(parentspermcomparison_false_LOH)
plot(mutantparents2_false_LOH$totaldepth.x,residefalseLOH, col="blue" , xlim=c(0,500), ylim=c(-0.6,0.8), xlab="Total Parent Read Depth", ylab="Parent-Sperm Residuals", main="Not Inherited + LOH")

cor.test(parentsaverage_true_LOH, spermprops_true_LOH)
cor.test(parentsaverage_false_GOH, spermprops_false_GOH)
cor.test(parentsaverage_false_LOH, spermprops_false_LOH)


#######
######
######
falseGOH<-data.frame(parentsaverage_false_GOH, spermprops_false_GOH)
falseGOH0.5<-falseGOH[falseGOH$parentsaverage_false_GOH <=0.5,]
falseGOH0.5.1<-falseGOH0.5[falseGOH0.5$spermprops_false_GOH <=0.5, ]
plot(falseGOH0.5.1, ylim=c(0,0.5),xlim=c(0,0.5))
falseGOHlm<-lm(falseGOH0.5.1$spermprops_false_GOH~falseGOH0.5.1$parentsaverage_false_GOH)
abline(falseGOHlm)
cor.test(falseGOH0.5.1$parentsaverage_false_GOH, falseGOH0.5.1$spermprops_false_GOH)

falseLOH<-data.frame(parentsaverage_false_LOH, spermprops_false_LOH)
falseLOH0.5<-falseLOH[falseLOH$parentsaverage_false_LOH <=0.5,]
falseLOH0.5.1<-falseLOH0.5[falseLOH0.5$spermprops_false_LOH <=0.5, ]
plot(falseLOH0.5.1, ylim=c(0,0.5),xlim=c(0,0.5))
falseLOHlm<-lm(falseLOH0.5.1$spermprops_false_LOH~falseLOH0.5.1$parentsaverage_false_LOH)
abline(falseLOHlm)
cor.test(falseLOH0.5.1$parentsaverage_false_LOH, falseLOH0.5.1$spermprops_false_LOH)

trueGOH<-data.frame(parentsaverage_true_GOH, spermprops_true_GOH)
trueGOH0.5<-trueGOH[trueGOH$parentsaverage_true_GOH <=0.5,]
trueGOH0.5.1<-trueGOH0.5[trueGOH0.5$spermprops_true_GOH <=0.5, ]
plot(trueGOH0.5.1, ylim=c(0,0.5),xlim=c(0,0.5))
trueGOHlm<-lm(trueGOH0.5.1$spermprops_true_GOH~trueGOH0.5.1$parentsaverage_true_GOH)
abline(trueGOHlm)
cor.test(trueGOH0.5.1$parentsaverage_true_GOH, trueGOH0.5.1$spermprops_true_GOH)

trueLOH<-data.frame(parentsaverage_true_LOH, spermprops_true_LOH)
trueLOH0.5<-trueLOH[trueLOH$parentsaverage_true_LOH <=0.5,]
trueLOH0.5.1<-trueLOH0.5[trueLOH0.5$spermprops_true_LOH <=0.5, ]
plot(trueLOH0.5.1, ylim=c(0,0.5),xlim=c(0,0.5))
trueLOHlm<-lm(trueLOH0.5.1$spermprops_true_LOH~trueLOH0.5.1$parentsaverage_true_LOH)
abline(trueLOHlm)
cor.test(trueLOH0.5.1$parentsaverage_true_LOH, trueLOH0.5.1$spermprops_true_LOH)

plot(falseGOH0.5.1, xlim=c(0,0.5),ylim=c(0,0.5))

points(falseLOH0.5.1, col="blue")
points(trueGOH0.5.1,col="red")
points(trueLOH0.5.1,col="pink")
abline(falseGOHlm, lwd=4)
abline(falseLOHlm, col="blue",lwd=4)
abline(trueGOHlm, col="red", lwd=4)
abline(trueLOHlm, col="pink",lwd=4)
legend("bottomright",
       legend = c("Not inherited+GOH","Not inherited+LOH","Inherited+GOH","Inherited+LOH"),
       col=c("black", "blue","red","pink"),
       pch=c(16,16))

summary(falseGOH0.5.1$parentsaverage_false_GOH)
newxfalseGOH<-seq(min(falseGOH0.5.1), max(falseGOH0.5.1),by = 0.01)#length.out=length(falseGOH0.5.1$parentsaverage_false_GOH))
conffalseGOH<-predict(falseGOHlm, newdata=data.frame(parentsaverage_false_GOH=newxfalseGOH), interval="confidence", level=0.95)
lines(newxfalseGOH, conffalseGOH[,2])
lines(newxfalseGOH, conffalseGOH[,3])

#plot(trueGOH,col="red")
#newxtrueGOH<-seq(min(trueGOH), max(trueGOH),length.out=length(trueGOH0.5.1$parentsaverage_true_GOH))
#conftrueGOH<-predict(trueGOHlm, newdata=data.frame(newxtrueGOH), interval="confidence")
#lines(newxtrueGOH, conftrueGOH[,2], col="black", lty=2)
#lines(newxtrueGOH, conftrueGOH[,3], col="black",lty=2)
#plot(spermprops_true_GOH~parentsaverage_true_GOH)
#lines(lm(spermprops_true_GOH~parentsaverage_true_GOH), col = Pal()[1], lwd = 2, lty = "solid",
type = "l", n = 100, conf.level = 0.95, args.cband = NULL,
pred.level = NA, args.pband = NULL)
uniquemetadatadf<- metadatadf[match(unique(metadatadf$chrom.pos), 					metadatadf$chrom.pos),]

par(mfrow=c(1,1))

#write.table(uniquemetadatadf, file="CAcolony60_CAP23-24muts_unique_20191125.txt",sep="\t",quote=FALSE, row.name=FALSE)

inh<-as.data.frame(Inheritance)
attach(somaticmetadatadf)
somaticmetadatadf <- somaticmetadatadf[order(chrom.pos),]
detach(somaticmetadatadf)
withinh<-cbind(mutantparentdf,inh)
TRUEset<-subset(withinh, Inheritance=="inherited")
print(nrow(TRUEset))
TRUEsetLoH<-subset(TRUEset, GoH_or_LoH =="LoH")
TRUEsetGoH<-subset(TRUEset, GoH_or_LoH =="DeNovo")

TRUEsetdepth<-TRUEset$totaldepth.y
TRUEsetmindepth<-TRUEset$totaldepth
TRUEsetGQ<-TRUEset$GQscore.y
TRUEsetminGQ<-TRUEset$GQscore
FALSEset<-subset(withinh, Inheritance=="notinherited")
print(nrow(FALSEset))
FALSEsetLoH<-subset(FALSEset, GoH_or_LoH =="LoH")
FALSEsetGoH<-subset(FALSEset, GoH_or_LoH =="DeNovo")

FALSEsetGQ<-FALSEset$GQscore.y
FALSEsetminGQ<-FALSEset$GQscore

FALSEsetdepth<-FALSEset$totaldepth.y
FALSEsetmindepth<-FALSEset$totaldepth

#plot minimum depth fortrue vs FALSE
x<- c(TRUEsetdepth, FALSEsetdepth)
groups<-c(rep("Inherited",nrow(TRUEset)),rep("Not inherited",nrow(FALSEset)))
CA56depths<-somaticCA56$totaldepth
CA60depths<-somaticCA60$totaldepth
CA65depths<-somaticCA65$totaldepth
y<-c(CA56depths, CA60depths, CA65depths)
groupsy<-c(rep("CA56",length(CA56depths)),rep("CA60",length(CA60depths)), rep("CA65",length(CA65depths)))
dfy<-data.frame(groupsy, y)
py<-ggplot(dfy,aes(groupsy,y)) + geom_boxplot() 
df<-data.frame(groups, x)
p<- ggplot(df, aes(groups,x))
depthsplots<- p +geom_boxplot() + geom_sina(aes(color=groups),size=1 ) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "", y = "Average Read Depth")+
  ylim(0,500)# here is your depths plot, FIGURE S1
depthsplots
wilcox.test(TRUEsetdepth,FALSEsetdepth) #use wilcoxon instead of t test
wilcox.test(CA56depths,CA65depths)
hist(FALSEsetdepth)
mean(TRUEsetdepth)
median(TRUEsetdepth)
mean(FALSEsetdepth)
median(FALSEsetdepth)
#plot minimum GQ score for true vs false
x<- c(TRUEsetGQ, FALSEsetGQ)
groups<-c(rep("Inherited",nrow(TRUEset)),rep("Not inherited",nrow(FALSEset)))
df<-data.frame(groups, x)
p<- ggplot(df, aes(groups,x))
GQplots<- p +geom_boxplot() + geom_sina(aes(color=groups),size=1 ) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "", y = "Average GQ at locus") # here is your depths plot, FIGURE S1
wilcox.test(TRUEsetGQ,FALSEsetGQ) #use wilcoxon instead of t test
hist(FALSEsetminGQ)
hist(TRUEsetminGQ)
median(FALSEsetGQ)
median(TRUEsetGQ)
GQplots
#plot min depthvs min gq
lminherited<-lm(TRUEsetmindepth~TRUEsetminGQ)

lmnotinherited<-lm(FALSEsetmindepth~FALSEsetminGQ)

#plot(TRUEsetdepth, TRUEsetGQ)
plot(FALSEsetdepth, FALSEsetGQ, ylim=c(30,100),xlab="Read Depth",ylab="GQ Score", xlim=c(0,500)) #only this plot!!
points(TRUEsetdepth, TRUEsetGQ,col="red")


abline(lmnotinherited)
abline(lminherited, col="red")

#plot where each mutation occurs
TRUEscaffolds<-TRUEset$chrom

FALSEscaffolds<-FALSEset$chrom

TRUEtable<-table(TRUEscaffolds)
head(TRUEtable,14)
FALSEtable<-table(FALSEscaffolds)
head(FALSEtable,14)
barplot(head(TRUEtable,14))

barplot(head(FALSEtable,14))

TOTALscaffolds<-uniquemetadatadf$chrom
TOTALtable<-table(TOTALscaffolds)
barplot(head(TOTALtable,14))
barplot(head(TRUEtable,14), col="blue")
barplot(head(FALSEtable,14),col="red", beside=TRUE)

table<-table(FALSEscaffolds,TRUEscaffolds)
#FOR PLOTS OF EACH INDIV SAMPLE:
files<-list.files(path="~/Documents/GitHub/Germline", pattern="*muts.txt", full.names=T, recursive=FALSE)

par(mfrow=c(3,2))

lapply(files,function(x) {
  metadata<-read.delim(x)
  genoanddepth<-(metadata$genotype)
  split<-str_split_fixed(genoanddepth, ",", 4)
  genotypes<-split[,1]
  #totaldepth<-as.numeric(split[,2])
  refdepth<-as.numeric(split[,2])
  altdepth<-as.numeric(split[,3])
  totaldepth<-refdepth+altdepth
  GQscore<-as.numeric(split[,4])
  
  position<-metadata$chrom.pos
  positionsplit<-str_split_fixed(position, "[.]", 2)
  
  chr<-positionsplit[,1]
  pos<-positionsplit[,2]
  
  metadatadf<-data.frame("chrom.pos" = metadata$chrom.pos, "chrom"=chr, "pos"=pos,	"sample"= metadata$sample, "ref" = metadata$ref, "alt" = 	metadata$alt, "genotype"= genotypes, "totaldepth"=totaldepth, 	"refdepth"=refdepth, "altdepth"=altdepth, "GQscore"= GQscore,	"GoH_or_LoH"=metadata$DeNovo_LoH, "Ti/Tv"=metadata$TiTv, 	"WhattoWhat" = metadata$WhattoWhat,"TrueorFalse" =metadata$TrueorFalse)#  ColonyName"=metadata$colonyrep)
  
  DepthMeansdf<-aggregate(totaldepth~chrom.pos, metadatadf, 			FUN=mean)
  
  DepthMinsdf<-aggregate(totaldepth~chrom.pos, metadatadf, 			FUN=min)
  
  GQaverage<-aggregate(GQscore~chrom.pos, metadatadf, FUN=mean)
  
  GQmin<-aggregate(GQscore~chrom.pos, metadatadf, FUN=min)
  
  metadatadf.00<-merge(metadatadf, DepthMeansdf[, c("chrom.pos", 	"totaldepth")], by="chrom.pos")
  metadatadf.0<-merge(metadatadf.00, GQaverage[,c("chrom.pos","GQscore")], by="chrom.pos")
  metadatadf.0<-merge(metadatadf.0, DepthMinsdf[,c("chrom.pos","totaldepth")], by="chrom.pos")
  metadatadf.0<-merge(metadatadf.0, GQmin[,c("chrom.pos","GQscore")], by="chrom.pos")
  
  DeNovos<-subset(metadatadf.0, GoH_or_LoH=="DeNovo")
  sample3<-subset(DeNovos, sample== "sample3")
  trueDenovos_sample3<-subset(sample3, refdepth =="0" | altdepth=="0")
  
  sample4<-subset(DeNovos, sample== "sample4")
  trueDenovos_sample4<-subset(sample4, refdepth =="0" | altdepth=="0")
  
  sample5<-subset(DeNovos, sample== "sample5")
  trueDenovos_sample5<-subset(sample5, refdepth =="0" | altdepth=="0")
  
  sample6<-subset(DeNovos, sample== "sample6")
  trueDenovos_sample6<-subset(sample6, refdepth =="0" | altdepth=="0")
  
  sample7<-subset(DeNovos, sample== "sample7")
  trueDenovos_sample7<-subset(sample7, refdepth =="0" | altdepth=="0")
  
  sample8<-subset(DeNovos, sample== "sample8")
  trueDenovos_sample8<-subset(sample8, refdepth =="0" | altdepth=="0")
  
  truedenovos3_8<-rbind(trueDenovos_sample3, trueDenovos_sample4, trueDenovos_sample5, trueDenovos_sample6, trueDenovos_sample7, trueDenovos_sample8)
  
  LoH<-subset(metadatadf.0, GoH_or_LoH =="LoH")#
  trueLoH<-subset(LoH, refdepth =="0" | altdepth=="0")
  metadatadf<-rbind( DeNovos, trueLoH)
  #write.table(metadatadf, file="CAcolony60_CAP22-23-24muts_20191125.txt",sep="\t",quote=FALSE, row.name=FALSE)
  
  uniquemetadatadf<- metadatadf[match(unique(metadatadf$chrom.pos), 					metadatadf$chrom.pos),]
  #write.table(uniquemetadatadf, file="CAcolony60_CAP23-24muts_unique_20191125.txt",sep="\t",quote=FALSE, row.name=FALSE)
  
  ##alt allele correlation when inheritance is TRUE:
  mutantparentdf_true<-subset(metadatadf.0, sample=="mutparent1" & TrueorFalse=="True") #for CAP22-1 and CAP22-2
  mutantparentdf_true<- mutantparentdf_true[match(uniquemetadatadf$chrom.pos, mutantparentdf_true$chrom.pos),]
  
  mutantparents2_true<-subset(metadatadf.0,sample=="mutparent2" & TrueorFalse=="True")
  mutantparents2_true<- mutantparents2_true[match(uniquemetadatadf$chrom.pos, mutantparents2_true$chrom.pos),]
  
  props1_true<-mutantparentdf_true$altdepth/mutantparentdf_true$totaldepth.x
  props2_true<-mutantparents2_true$altdepth/mutantparents2_true$totaldepth.x
  #plot(props1, props2)
  
  #propslm<-lm(props2~props1)
  
  parentsaverage_true<-(props1_true+props2_true)/2
  parentsaverage_true<-na.omit(parentsaverage_true)
  
  
  spermdf_true<-subset(metadatadf.0, sample=="mutsperm" & TrueorFalse=="True")
  spermdf_true<- spermdf_true[match(uniquemetadatadf$chrom.pos, spermdf_true$chrom.pos),]
  
  spermprops_true<- spermdf_true$altdepth/spermdf_true$totaldepth.x
  spermprops_true<-na.omit(spermprops_true)
  
  parentspermcomparison_true<-lm(spermprops_true~parentsaverage_true)
  #plot(parentsaverage_true, spermprops_true)
  #abline(parentspermcomparison_true)
  
  
  ##alt allele correlation when inheritance isFALSE:
  mutantparentdf_false<-subset(metadatadf.0, sample=="mutparent1" & TrueorFalse=="False") #for CAP22-1 and CAP22-2
  mutantparentdf_false<- mutantparentdf_false[match(uniquemetadatadf$chrom.pos, mutantparentdf_false$chrom.pos),]
  
  mutantparents2_false<-subset(metadatadf.0,sample=="mutparent2" & TrueorFalse=="False")
  mutantparents2_false<- mutantparents2_false[match(uniquemetadatadf$chrom.pos, mutantparents2_false$chrom.pos),]
  
  props1_false<-mutantparentdf_false$altdepth/mutantparentdf_false$totaldepth.x
  props2_false<-mutantparents2_false$altdepth/mutantparents2_false$totaldepth.x
  #plot(props1, props2)
  
  #propslm<-lm(props2~props1)
  
  parentsaverage_false<-(props1_false+props2_false)/2
  parentsaverage_false<-na.omit(parentsaverage_false)
  
  
  spermdf_false<-subset(metadatadf.0, sample=="mutsperm" & TrueorFalse=="False")
  spermdf_false<- spermdf_false[match(uniquemetadatadf$chrom.pos, spermdf_false$chrom.pos),]
  
  spermprops_false<- spermdf_false$altdepth/spermdf_false$totaldepth.x
  spermprops_false<-na.omit(spermprops_false)
  
  parentspermcomparison_false<-lm(spermprops_false~parentsaverage_false)
  #plot(parentsaverage_false, spermprops_false)
  #points(parentsaverage_true,spermprops_true,col="blue")
  #abline(parentspermcomparison_false)
  #abline(parentspermcomparison_true,col="blue")
  summary(parentspermcomparison_false)$r.squared
  cor.test(parentsaverage_false,spermprops_false)
  
  ##comparing GOH and LOH when inheritance is false:
  #LOH:
  mutantparentdf_false_LOH<-subset(mutantparentdf_false, GoH_or_LoH=="LoH")
  mutantparents2_false_LOH<-subset(mutantparents2_false, GoH_or_LoH=="LoH")
  
  props1_false_LOH<-mutantparentdf_false_LOH$altdepth/mutantparentdf_false_LOH$totaldepth.x
  props2_false_LOH<-mutantparents2_false_LOH$altdepth/mutantparents2_false_LOH$totaldepth.x
  
  parentsaverage_false_LOH<-(props1_false_LOH+props2_false_LOH)/2
  parentsaverage_false_LOH<-na.omit(parentsaverage_false_LOH)
  
  spermdf_false_LOH<- subset(spermdf_false, GoH_or_LoH=="LoH")
  
  spermprops_false_LOH<-spermdf_false_LOH$altdepth/spermdf_false_LOH$totaldepth.x
  spermprops_false_LOH<-na.omit(spermprops_false_LOH)
  
  parentspermcomparison_false_LOH<-lm(spermprops_false_LOH~parentsaverage_false_LOH)
  #plot(parentsaverage_false_LOH, spermprops_false_LOH)
  #abline(parentspermcomparison_false_LOH)
  summary(parentspermcomparison_false_LOH)$r.squared
  
  #GOH:
  mutantparentdf_false_GOH<-subset(mutantparentdf_false, GoH_or_LoH=="DeNovo")
  mutantparents2_false_GOH<-subset(mutantparents2_false, GoH_or_LoH=="DeNovo")
  
  props1_false_GOH<-mutantparentdf_false_GOH$altdepth/mutantparentdf_false_GOH$totaldepth.x
  props2_false_GOH<-mutantparents2_false_GOH$altdepth/mutantparents2_false_GOH$totaldepth.x
  
  parentsaverage_false_GOH<-(props1_false_GOH+props2_false_GOH)/2
  parentsaverage_false_GOH<-na.omit(parentsaverage_false_GOH)
  
  spermdf_false_GOH<- subset(spermdf_false, GoH_or_LoH=="DeNovo")
  
  spermprops_false_GOH<-spermdf_false_GOH$altdepth/spermdf_false_GOH$totaldepth.x
  spermprops_false_GOH<-na.omit(spermprops_false_GOH)
  
  parentspermcomparison_false_GOH<-lm(spermprops_false_GOH~parentsaverage_false_GOH)
  #plot(parentsaverage_false_GOH, spermprops_false_GOH,xlab="Proportion of alt allele in the parent",ylab="Proportion of alt allele in sperm",xlim=c(0,.5),ylim=c(0,0.5))
  
  
  summary(parentspermcomparison_false_GOH)$r.squared  
  cor.test(parentsaverage_false_LOH, spermprops_false_LOH)
  cor.test(parentsaverage_false_GOH, spermprops_false_GOH)
  
  ##comparing GOH and LOH when inheritance is true:
  #LOH:
  mutantparentdf_true_LOH<-subset(mutantparentdf_true, GoH_or_LoH=="LoH")
  mutantparents2_true_LOH<-subset(mutantparents2_true, GoH_or_LoH=="LoH")
  
  props1_true_LOH<-mutantparentdf_true_LOH$altdepth/mutantparentdf_true_LOH$totaldepth.x
  props2_true_LOH<-mutantparents2_true_LOH$altdepth/mutantparents2_true_LOH$totaldepth.x
  
  parentsaverage_true_LOH<-(props1_true_LOH+props2_true_LOH)/2
  parentsaverage_true_LOH<-na.omit(parentsaverage_true_LOH)
  
  
  spermdf_true<-subset(metadatadf.0, sample=="mutsperm" & TrueorFalse=="True")
  spermdf_true<- spermdf_true[match(uniquemetadatadf$chrom.pos, spermdf_true$chrom.pos),]
  
  spermdf_true_LOH<- subset(spermdf_true, GoH_or_LoH=="LoH")
  
  spermprops_true_LOH<-spermdf_true_LOH$altdepth/spermdf_true_LOH$totaldepth.x
  spermprops_true_LOH<-na.omit(spermprops_true_LOH)
  
  parentspermcomparison_true_LOH<-lm(spermprops_true_LOH~parentsaverage_true_LOH)
  #plot(parentsaverage_true_LOH, spermprops_true_LOH)
  #abline(parentspermcomparison_true_LOH)
  summary(parentspermcomparison_true_LOH)$r.squared
  
  #GOH:
  mutantparentdf_true_GOH<-subset(mutantparentdf_true, GoH_or_LoH=="DeNovo")
  mutantparents2_true_GOH<-subset(mutantparents2_true, GoH_or_LoH=="DeNovo")
  
  props1_true_GOH<-mutantparentdf_true_GOH$altdepth/mutantparentdf_true_GOH$totaldepth.x
  props2_true_GOH<-mutantparents2_true_GOH$altdepth/mutantparents2_true_GOH$totaldepth.x
  
  parentsaverage_true_GOH<-(props1_true_GOH+props2_true_GOH)/2
  parentsaverage_true_GOH<-na.omit(parentsaverage_true_GOH)
  
  
  spermdf_true<-subset(metadatadf.0, sample=="mutsperm" & TrueorFalse=="True")
  spermdf_true<- spermdf_true[match(uniquemetadatadf$chrom.pos, spermdf_true$chrom.pos),]
  
  spermdf_true_GOH<- subset(spermdf_true, GoH_or_LoH=="DeNovo")
  
  spermprops_true_GOH<-spermdf_true_GOH$altdepth/spermdf_true_GOH$totaldepth.x
  spermprops_true_GOH<-na.omit(spermprops_true_GOH)
  
  parentspermcomparison_true_GOH<-lm(spermprops_true_GOH~parentsaverage_true_GOH)
  #plot(parentsaverage_true_GOH, spermprops_true_GOH)
  #points(parentsaverage_true_LOH, spermprops_true_LOH, col="blue")
  #abline(parentspermcomparison_true_GOH)
  #abline(parentspermcomparison_true_LOH,col="blue")
  summary(parentspermcomparison_true_GOH)$r.squared  
  cor.test(parentsaverage_true_LOH, spermprops_true_LOH)
  cor.test(parentsaverage_true_GOH, spermprops_true_GOH)
  #plot(parentsaverage_false_GOH, spermprops_false_GOH,xlab="Proportion of alt allele in the parent",ylab="Proportion of alt allele in sperm")#,xlim=c(0,.5),ylim=c(0,0.5))
  #points(parentsaverage_false_LOH, spermprops_false_LOH, col="blue")
  #points(parentsaverage_true_GOH, spermprops_true_GOH, col="red")
  #points(parentsaverage_true_LOH,spermprops_true_LOH,col="pink")
  #abline(parentspermcomparison_false_GOH)
  #abline(parentspermcomparison_false_LOH,col="blue")
  #abline(parentspermcomparison_true_GOH,col="red")
  #abline(parentspermcomparison_true_LOH,col="pink")
  
  falseGOH<-data.frame(parentsaverage_false_GOH, spermprops_false_GOH)
  falseGOH0.5<-falseGOH[falseGOH$parentsaverage_false_GOH <=0.5,]
  falseGOH0.5.1<-falseGOH0.5[falseGOH0.5$spermprops_false_GOH <=0.5, ]
  #plot(falseGOH0.5.1, ylim=c(0,0.5),xlim=c(0,0.5))
  falseGOHlm<-lm(falseGOH0.5.1$spermprops_false_GOH~falseGOH0.5.1$parentsaverage_false_GOH)
  #abline(falseGOHlm)
  cor.test(falseGOH0.5.1$parentsaverage_false_GOH, falseGOH0.5.1$spermprops_false_GOH)
  
  falseLOH<-data.frame(parentsaverage_false_LOH, spermprops_false_LOH)
  falseLOH0.5<-falseLOH[falseLOH$parentsaverage_false_LOH <=0.5,]
  falseLOH0.5.1<-falseLOH0.5[falseLOH0.5$spermprops_false_LOH <=0.5, ]
  #plot(falseLOH0.5.1, ylim=c(0,0.5),xlim=c(0,0.5))
  falseLOHlm<-lm(falseLOH0.5.1$spermprops_false_LOH~falseLOH0.5.1$parentsaverage_false_LOH)
  #abline(falseLOHlm)
  cor.test(falseLOH0.5.1$parentsaverage_false_LOH, falseLOH0.5.1$spermprops_false_LOH)
  
  trueGOH<-data.frame(parentsaverage_true_GOH, spermprops_true_GOH)
  trueGOH0.5<-trueGOH[trueGOH$parentsaverage_true_GOH <=0.5,]
  trueGOH0.5.1<-trueGOH0.5[trueGOH0.5$spermprops_true_GOH <=0.5, ]
  #plot(trueGOH0.5.1, ylim=c(0,0.5),xlim=c(0,0.5))
  trueGOHlm<-lm(trueGOH0.5.1$spermprops_true_GOH~trueGOH0.5.1$parentsaverage_true_GOH)
  #abline(trueGOHlm)
  cor.test(trueGOH0.5.1$parentsaverage_true_GOH, trueGOH0.5.1$spermprops_true_GOH)
  
  trueLOH<-data.frame(parentsaverage_true_LOH, spermprops_true_LOH)
  trueLOH0.5<-trueLOH[trueLOH$parentsaverage_true_LOH <=0.5,]
  trueLOH0.5.1<-trueLOH0.5[trueLOH0.5$spermprops_true_LOH <=0.5, ]
  #plot(trueLOH0.5.1, ylim=c(0,0.5),xlim=c(0,0.5))
  trueLOHlm<-lm(trueLOH0.5.1$spermprops_true_LOH~trueLOH0.5.1$parentsaverage_true_LOH)
  #abline(trueLOHlm)
  cor.test(trueLOH0.5.1$parentsaverage_true_LOH, trueLOH0.5.1$spermprops_true_LOH)
  
  #plot(falseGOH0.5.1, xlim=c(0,0.5),ylim=c(0,0.5))
  points(falseLOH0.5.1, col="blue")
  points(trueGOH0.5.1,col="red")
  points(trueLOH0.5.1,col="pink")
  abline(falseGOHlm, lwd=4)
  abline(falseLOHlm, col="blue",lwd=4)
  abline(trueGOHlm, col="red", lwd=4)
  abline(trueLOHlm, col="pink",lwd=4)
  legend("bottomright",
         legend = c("Not inherited+GOH","Not inherited+LOH","Inherited+GOH","Inherited+LOH"),
         col=c("black", "blue","red","pink"),
         pch=c(16,16))
})
lapply(files,function(x) {
  metadata<-read.delim(x)
  genoanddepth<-(metadata$genotype)
  split<-str_split_fixed(genoanddepth, ",", 4)
  genotypes<-split[,1]
  #totaldepth<-as.numeric(split[,2])
  refdepth<-as.numeric(split[,2])
  altdepth<-as.numeric(split[,3])
  totaldepth<-refdepth+altdepth
  GQscore<-as.numeric(split[,4])
  
  position<-metadata$chrom.pos
  positionsplit<-str_split_fixed(position, "[.]", 2)
  
  chr<-positionsplit[,1]
  pos<-positionsplit[,2]
  
  metadatadf<-data.frame("chrom.pos" = metadata$chrom.pos, "chrom"=chr, "pos"=pos,	"sample"= metadata$sample, "ref" = metadata$ref, "alt" = 	metadata$alt, "genotype"= genotypes, "totaldepth"=totaldepth, 	"refdepth"=refdepth, "altdepth"=altdepth, "GQscore"= GQscore,	"GoH_or_LoH"=metadata$DeNovo_LoH, "Ti/Tv"=metadata$TiTv, 	"WhattoWhat" = metadata$WhattoWhat,"TrueorFalse" =metadata$TrueorFalse)#  ColonyName"=metadata$colonyrep)
  
  DepthMeansdf<-aggregate(totaldepth~chrom.pos, metadatadf, 			FUN=mean)
  
  DepthMinsdf<-aggregate(totaldepth~chrom.pos, metadatadf, 			FUN=min)
  
  GQaverage<-aggregate(GQscore~chrom.pos, metadatadf, FUN=mean)
  
  GQmin<-aggregate(GQscore~chrom.pos, metadatadf, FUN=min)
  
  metadatadf.00<-merge(metadatadf, DepthMeansdf[, c("chrom.pos", 	"totaldepth")], by="chrom.pos")
  metadatadf.0<-merge(metadatadf.00, GQaverage[,c("chrom.pos","GQscore")], by="chrom.pos")
  metadatadf.0<-merge(metadatadf.0, DepthMinsdf[,c("chrom.pos","totaldepth")], by="chrom.pos")
  metadatadf.0<-merge(metadatadf.0, GQmin[,c("chrom.pos","GQscore")], by="chrom.pos")
  
  DeNovos<-subset(metadatadf.0, GoH_or_LoH=="DeNovo")
  sample3<-subset(DeNovos, sample== "sample3")
  trueDenovos_sample3<-subset(sample3, refdepth =="0" | altdepth=="0")
  
  sample4<-subset(DeNovos, sample== "sample4")
  trueDenovos_sample4<-subset(sample4, refdepth =="0" | altdepth=="0")
  
  sample5<-subset(DeNovos, sample== "sample5")
  trueDenovos_sample5<-subset(sample5, refdepth =="0" | altdepth=="0")
  
  sample6<-subset(DeNovos, sample== "sample6")
  trueDenovos_sample6<-subset(sample6, refdepth =="0" | altdepth=="0")
  
  sample7<-subset(DeNovos, sample== "sample7")
  trueDenovos_sample7<-subset(sample7, refdepth =="0" | altdepth=="0")
  
  sample8<-subset(DeNovos, sample== "sample8")
  trueDenovos_sample8<-subset(sample8, refdepth =="0" | altdepth=="0")
  
  truedenovos3_8<-rbind(trueDenovos_sample3, trueDenovos_sample4, trueDenovos_sample5, trueDenovos_sample6, trueDenovos_sample7, trueDenovos_sample8)
  
  LoH<-subset(metadatadf.0, GoH_or_LoH =="LoH")#
  trueLoH<-subset(LoH, refdepth =="0" | altdepth=="0")
  metadatadf<-rbind( DeNovos, trueLoH)
  #write.table(metadatadf, file="CAcolony60_CAP22-23-24muts_20191125.txt",sep="\t",quote=FALSE, row.name=FALSE)
  
  uniquemetadatadf<- metadatadf[match(unique(metadatadf$chrom.pos), 					metadatadf$chrom.pos),]
  #write.table(uniquemetadatadf, file="CAcolony60_CAP23-24muts_unique_20191125.txt",sep="\t",quote=FALSE, row.name=FALSE)
  TRUEset<-subset(uniquemetadatadf, TrueorFalse=="True")
  print(nrow(TRUEset))
  TRUEsetLoH<-subset(TRUEset, GoH_or_LoH =="LoH")
  TRUEsetGoH<-subset(TRUEset, GoH_or_LoH =="DeNovo")
  
  TRUEsetdepth<-TRUEset$totaldepth.y
  TRUEsetmindepth<-TRUEset$totaldepth
  TRUEsetGQ<-TRUEset$GQscore.y
  TRUEsetminGQ<-TRUEset$GQscore
  FALSEset<-subset(uniquemetadatadf, TrueorFalse=="False")
  print(nrow(FALSEset))
  FALSEsetLoH<-subset(FALSEset, GoH_or_LoH =="LoH")
  FALSEsetGoH<-subset(FALSEset, GoH_or_LoH =="DeNovo")
  
  FALSEsetGQ<-FALSEset$GQscore.y
  FALSEsetminGQ<-FALSEset$GQscore
  
  FALSEsetdepth<-FALSEset$totaldepth.y
  FALSEsetmindepth<-FALSEset$totaldepth
  plot(FALSEsetdepth, FALSEsetGQ, ylim=c(30,100),xlab="Read Depth",ylab="GQ Score", xlim=c(0,500)) #only this plot!!
  points(TRUEsetdepth, TRUEsetGQ,col="red")
})
#plot minimum depth fortrue vs FALSE
x<- c(TRUEsetmindepth, FALSEsetmindepth)
groups<-c(rep("Inherited",nrow(TRUEset)),rep("Not inherited",nrow(FALSEset)))
df<-data.frame(groups, x)
p<- ggplot(df, aes(groups,x))
depthsplots<- p +geom_boxplot() + geom_sina(aes(color=groups),size=1 ) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "", y = "Read Depth") # here is your depths plot, FIGURE S1
#depthsplots
wilcox.test(TRUEsetmindepth,FALSEsetmindepth) #use wilcoxon instead of t test
#hist(FALSEsetmindepth)
mean(TRUEsetmindepth)
mean(FALSEsetmindepth)
#plot minimum GQ score for true vs false
x<- c(TRUEsetminGQ, FALSEsetminGQ)
groups<-c(rep("Inherited",nrow(TRUEset)),rep("Not inherited",nrow(FALSEset)))
df<-data.frame(groups, x)
p<- ggplot(df, aes(groups,x))
GQplots<- p +geom_boxplot() + geom_sina(aes(color=groups),size=1 ) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "", y = "Average GQ at locus") # here is your depths plot, FIGURE S1
wilcox.test(TRUEsetminGQ,FALSEsetminGQ) #use wilcoxon instead of t test
#hist(FALSEsetminGQ)
#hist(TRUEsetminGQ)
mean(FALSEsetminGQ)
mean(TRUEsetminGQ)
#GQplots
#plot min depthvs min gq
lminherited<-lm(TRUEsetmindepth~TRUEsetminGQ)

lmnotinherited<-lm(FALSEsetmindepth~FALSEsetminGQ)

plot(TRUEsetdepth, TRUEsetGQ)
#plot(FALSEsetmindepth, FALSEsetminGQ, ylim=c(0,100))
points(FALSEsetdepth, FALSEsetGQ,col="red")
#abline(lmnotinherited)
#abline(lminherited, col="red")

#plot where each mutation occurs
TRUEscaffolds<-TRUEset$chrom

FALSEscaffolds<-FALSEset$chrom

TRUEtable<-table(TRUEscaffolds)
head(TRUEtable,14)
FALSEtable<-table(FALSEscaffolds)
head(FALSEtable,14)
#barplot(head(TRUEtable,14))

#barplot(head(FALSEtable,14))

TOTALscaffolds<-uniquemetadatadf$chrom
TOTALtable<-table(TOTALscaffolds)
barplot(head(TOTALtable,14))
})


mat<-matrix(c((138/(138+804)), (139/(139+546)), (147/(147+1159)), (760/(760+1363)), (699/(699+793)), (484/(484+838)), (375/(375+861)), (328/(328+899)), (804/(138+804)), (546/(139+546)), (1159/(147+1159)), (1363/(760+1363)), (793/(699+793)), (838/(484+838)), (861/(375+861)), (899/(328+899))), ncol=8, byrow=TRUE)
rownames(mat)<-c("Inherited","Not Inherited")    
#rownames(mat)<-c("CAP22")
mat<-as.table(mat)
barplot(mat,col=c("red","black"), ylab="Proportion of somatic mutations")      

(804/(138+804)), (546/(139+546)), (1159/(147+1159)), (1363/(760+1363)), (793/(699+793)), (838/(484+838)), (861/(375+861)), (899/(328+899))

##ANOVA AND BOXPLOTS: proportion of inherited som muts ##
runs<-c((138/(138+804)), (139/(139+546)), (147/(147+1159)), (760/(760+1363)), (699/(699+793)), (484/(484+838)), (375/(375+861)), (328/(328+899)))#, (804/(138+804)), (546/(139+546)), (1159/(147+1159)), (1363/(760+1363)), (793/(699+793)), (838/(484+838)), (861/(375+861)), (899/(328+899)))
group<-c("CA60", "CA60", "CA60",  "CA65", "CA65", "CA56", "CA56", "CA56")  

inheritedproportion<- inheritedcounta_/(nrow(inheriteddenovoa_)+nrow(inheritedloha_))
colony<-c("CA56", "CA56", "CA56", "CA60", "CA60", "CA60",  "CA65", "CA65")
withinRunStats = function(x) c(sum = sum(x), mean = mean(x), var = var(x), n = length(x))
tapply(runs, group, withinRunStats)
data = data.frame(y = runs, group = factor(group))
data
fit = lm(runs ~ group, data)
fit
anova(fit)
oneway<-aov(runs~group)
summary(oneway)
plot(oneway)
df<-data.frame(runs,group)
degreesOfFreedom = anova(fit)[, "Df"]
names(degreesOfFreedom) = c("treatment", "error")
degreesOfFreedom
anova(fit)["Residuals", "Mean Sq"]
anova(fit)["group", "Mean Sq"]
inheritedpropplot<-boxplot(runs~group, xlab="Colony Name", ylab="Proportion of Somatic Mutations Inherited by Sperm", las=1, col=c("light blue","light gray","light green"))
a<-ggplot(df,aes(y=runs, x=group))
a2<- a+ geom_boxplot()+theme_bw() + ylim(0, 0.5) + labs(x = "", y = "Proportion of Somatic Mutations Inherited by Sperm")


stripchart(runs ~ group,
           vertical = TRUE, 
           pch = 21, col = "red", bg = "bisque",
           add = TRUE)
## ##

coding<-read.delim("~/Documents/GitHub/CoralGermline/codingsommut_types.txt")

data("ToothGrowth")
ToothGrowth$dose <- as.factor(ToothGrowth$dose)
head(ToothGrowth)
e <- ggplot(coding, aes(x = type, y = percent))
e<- e + geom_boxplot()
e2 <- e + labs(x="Coding Mutation Type", y="Percent of Total Coding Mutations") + geom_boxplot(
  aes(fill = colonynumber),
  position = position_dodge(0.9)
) +
  scale_fill_manual(values = c("#999999", "#E69F00","red"))
e2 


##plotting positions of LOH muts##
#first,splitintochrs:
CA56splitloh<-split(a_loh, a_loh$chrom)
CAP6split<-split(CAP6, CAP6$chrom)
chr1posloh<-CA56splitloh$chr1$pos
yaxloh<-rep(1, length(chr1posloh))
plot(chr1posloh, yaxloh)

CA56splitdenovo<-split(a_denovo, a_denovo$chrom)
chr1posdenovo<-CA56splitdenovo$chr1$pos
yaxdenovo<-rep(2, length(chr1posdenovo))
plot(chr1posdenovo, yaxdenovo)
difference_func<-function(chr1posloh) {
  difference<-rep("A",length(chr1posloh))
  for (i in 1:length(chr1posloh)){
    
    chr1posloh=as.numeric(as.character(chr1posloh))
    difference[i] = chr1posloh[i+1] - chr1posloh[i]
    
  }
  return(difference)
} 

splitfunc<-function(CA56splitloh) {
  chr1pos<-CA56splitloh$chr1$pos
  yax1<-rep(1, length(chr1pos))
  
  chr2pos<-CA56splitloh$chr2$pos
  yax2<-rep(2, length(chr2pos))
  
  chr3pos<-CA56splitloh$chr3$pos
  yax3<-rep(3, length(chr3pos))
  
  chr4pos<-CA56splitloh$chr4$pos
  yax4<-rep(4, length(chr4pos))
  
  chr5pos<-CA56splitloh$chr5$pos
  yax5<-rep(5, length(chr5pos))
  
  chr6pos<-CA56splitloh$chr6$pos
  yax6<-rep(6, length(chr6pos))
  
  chr7pos<-CA56splitloh$chr7$pos
  yax7<-rep(7, length(chr7pos))
  
  chr8pos<-CA56splitloh$chr8$pos
  yax8<-rep(8, length(chr8pos))
  
  chr9pos<-CA56splitloh$chr9$pos
  yax9<-rep(9, length(chr9pos))
  
  chr10pos<-CA56splitloh$chr10$pos
  yax10<-rep(10, length(chr10pos))
  
  chr11pos<-CA56splitloh$chr11$pos
  yax11<-rep(11, length(chr11pos))
  
  chr12pos<-CA56splitloh$chr12$pos
  yax12<-rep(12, length(chr12pos))
  
  chr13pos<-CA56splitloh$chr13$pos
  yax13<-rep(13, length(chr13pos))
  
  chr14pos<-CA56splitloh$chr14$pos
  yax14<-rep(14, length(chr14pos))
  
  difference1<- difference_func(chr1pos)
  chr1frame<-data.frame("pos"=as.numeric(as.character(chr1pos)), "yax"=yax1, "difference"=difference1)
  
  difference2<- difference_func(chr2pos)
  chr2frame<-data.frame("pos"=as.numeric(as.character(chr2pos)), "yax"=yax2, "difference"=difference2)
  
  difference3<- difference_func(chr3pos)
  chr3frame<-data.frame("pos"=as.numeric(as.character(chr3pos)), "yax"=yax3, "difference"=difference3)
  
  difference4<- difference_func(chr4pos)
  chr4frame<-data.frame("pos"=as.numeric(as.character(chr4pos)), "yax"=yax4, "difference"=difference4)
  
  difference5<- difference_func(chr5pos)
  chr5frame<-data.frame("pos"=as.numeric(as.character(chr5pos)), "yax"=yax5, "difference"=difference5)
  
  difference6<- difference_func(chr6pos)
  chr6frame<-data.frame("pos"=as.numeric(as.character(chr6pos)), "yax"=yax6, "difference"=difference6)
  
  difference7<- difference_func(chr7pos)
  chr7frame<-data.frame("pos"=as.numeric(as.character(chr7pos)), "yax"=yax7, "difference"=difference7)
  
  difference8<- difference_func(chr8pos)
  chr8frame<-data.frame("pos"=as.numeric(as.character(chr8pos)), "yax"=yax8, "difference"=difference8)
  
  difference9<- difference_func(chr9pos)
  chr9frame<-data.frame("pos"=as.numeric(as.character(chr9pos)), "yax"=yax9, "difference"=difference9)
  
  difference10<- difference_func(chr10pos)
  chr10frame<-data.frame("pos"=as.numeric(as.character(chr10pos)), "yax"=yax10, "difference"=difference10)
  
  difference11<- difference_func(chr11pos)
  chr11frame<-data.frame("pos"=as.numeric(as.character(chr11pos)), "yax"=yax11, "difference"=difference11)
  
  difference12<- difference_func(chr12pos)
  chr12frame<-data.frame("pos"=as.numeric(as.character(chr12pos)), "yax"=yax12, "difference"=difference12)
  
  difference13<- difference_func(chr13pos)
  chr13frame<-data.frame("pos"=as.numeric(as.character(chr13pos)), "yax"=yax13, "difference"=difference13)
  
  difference14<- difference_func(chr14pos)
  chr14frame<-data.frame("pos"=as.numeric(as.character(chr14pos)), "yax"=yax14, "difference"=difference14)
  
  chrplot<-ggplot(chr1frame, aes(x=pos, y=yax )) + geom_point(size=1) +
    
    theme_bw() +
    geom_point(data=chr2frame, aes(x=pos, y=yax), size=1) +
    geom_point(data=chr3frame, aes(x=pos, y=yax), size=1) +
    geom_point(data=chr4frame, aes(x=pos, y=yax), size=1) +
    geom_point(data=chr5frame, aes(x=pos, y=yax), size=1) +
    geom_point(data=chr6frame, aes(x=pos, y=yax), size=1) +
    geom_point(data=chr7frame, aes(x=pos, y=yax), size=1) +
    geom_point(data=chr8frame, aes(x=pos, y=yax), size=1) +
    geom_point(data=chr9frame, aes(x=pos, y=yax), size=1) +
    geom_point(data=chr10frame, aes(x=pos, y=yax), size=1) +
    geom_point(data=chr11frame, aes(x=pos, y=yax), size=1) +
    geom_point(data=chr12frame, aes(x=pos, y=yax), size=1) +
    geom_point(data=chr13frame, aes(x=pos, y=yax), size=1) +
    geom_point(data=chr14frame, aes(x=pos, y=yax), size=1) 
  #return(list(chr1frame,chr2frame,chr3frame,chr4frame,chr5frame,chr6frame,chr7frame,chr8frame,chr9frame,chr10frame,chr11frame,chr12frame,chr13frame,chr14frame))
  return(chrplot)
}
CA56splitlohframes<-splitfunc(CA56splitloh)
CA56splitdenovoframes<-splitfunc(CA56splitdenovo)
CAP6splitplot<-splitfunc(CAP6split)

lohdifference<-difference_func(chr1posloh)
denovodifference<-difference_func(chr1posdenovo)
CA56chr1lohframe<-data.frame("yax"=yaxloh, "pos"=as.numeric(as.character(chr1posloh)), "difference"=lohdifference)
CA56chr1denovoframe<-data.frame("yax"=yaxdenovo, "pos"=as.numeric(as.character(chr1posdenovo)),"difference"=denovodifference)

CA56chr1plot<-ggplot(CA56chr1lohframe, aes(x=pos, y=yax )) + geom_point(size=1) +
  
  theme_bw() +
  geom_point(data=CA56chr1denovoframe, aes(x=pos, y=yax),color="red")
  #scale_color_manual(values = c("red","black")) +
  #ylab("Average Variant Allele Frequency in the Mutant Sperm Pool") + xlab("Average Variant Allele Frequency in the Mutant Parent")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15))+
  theme(legend.text=element_text(size=15))#+
  geom_point(data=CA56chr1denovoframe, x=pos, y=yax,color="red")
