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
library(seqinr)
sessionInfo()
#denominators:
CA56denom<-79943111-(895*14)#11/17  #80718433-(892*14) #9/14
#CA60denom<-   51791529-(895*14)#11/17    #52564929-(892*14) #9/14
CA65denom<-   47410341-(895*14)#11/17    #48319810-(892*14) #9/14
CA60denom<- 76882242-(895*14)#12/3
CA56denom_coding<- 14123208
CA60denom_coding<- 3823690+9622397 #12/3
CA65denom_coding<- 8462566

#DENOMINATORS MAPPED TO AHYA (don't use these, use the ones in DownstreamAnalysis_clean.R)
#CA56denom<-113503344-(14*949)
#CA60denom<-104072836-(14*949)
#CA65denom<-56924223-(14*949)
#CA60_denom_coding<-
#surface areas:
CA56sa<-36.25
CA60sa<-253.02
CA65sa<-704.15

threedfs_func<- function(files) {
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
  
  metadatadf<-data.frame("chrom.pos" = metadata$chrom.pos, "chrom"=chr, "pos"=pos,	"sample"= metadata$sample, "ref" = metadata$ref, "alt" = 	metadata$alt, "normal_allele"= normalallele,
                         "mutant_allele" = mutantallele, "mutant_allele_depth" = as.numeric(mutant_alleledepth), "genotype"= genotypes, "totaldepth"=totaldepth, 	"refdepth"=refdepth,
                         "altdepth"=altdepth, "GQscore"= GQscore,	"GoH_or_LoH"=metadata$DeNovo_LoH, "Ti/Tv"=metadata$TiTv, 	"WhattoWhat" = metadata$WhattoWhat,
                         "TrueorFalse" =metadata$TrueorFalse, "MutationClass"=metadata$TypeofMutation, "MutantParent1"=metadata$MutantParent1,"MutantParent2"=metadata$MutantParent2, 
                         "MutantSperm1"=metadata$MutantSperm1, "MutantSperm2"=metadata$MutantSperm2, "ColonyName"=metadata$colonyrep, stringsAsFactors=FALSE)#  ColonyName"=metadata$colonyrep)
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
  metadatadf.0<-merge(metadatadf.0, GQmin[,c("chrom.pos","GQscore")], by="chrom.pos")
  #metadatadf.0<-subset(metadatadf.1, startsWith(metadatadf.1$chrom, 'chr')==TRUE)
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
  #trueLoH_somatic_mpandms<-subset(trueLoH_somatic_mpandms, chrom.pos != "chr9.10238326")
  #trueLoH_somatic_mpandms<-subset(trueLoH_somatic_mpandms, chrom.pos != "chr13.20768582")
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
  
  somaticmetadatadf<-rbind(trueDenovos_somatic_mpandms, trueLoH_somatic_mpandms)#, CA65somaticmpandms) #combines all of the filtered somatic datasets into one
  somaticmetadatadf_denovo<-subset(somaticmetadatadf, GoH_or_LoH=="DeNovo")
  somaticmetadatadf_loh<-subset(somaticmetadatadf, GoH_or_LoH=="LoH")
  uniquesomaticmetadatadf<-somaticmetadatadf[match(unique(somaticmetadatadf$chrom.pos), 					somaticmetadatadf$chrom.pos),]
  uniquesomaticmetadatadf_denovo<-subset(uniquesomaticmetadatadf, GoH_or_LoH=="DeNovo")
  uniquesomaticmetadatadf_loh<-subset(uniquesomaticmetadatadf, GoH_or_LoH=="LoH")
  #gglm true LoH:
  globalglm<-subset(metadatadf.0, MutationClass=="GlobalGermlineMutation")
  
  justsperm<-subset(globalglm, startsWith(globalglm$sample, 'CAS')==TRUE)
  
  if (nrow(justsperm) <1) {
    gglmmetadatadf<- data.frame(matrix(ncol = 30, nrow = 0))
    colnames(gglmmetadatadf)<-colnames(somaticmetadatadf)
  } else {
    #justsperm<-subset(globalglm, startsWith(globalglm$sample, 'CAS')==TRUE)
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
    gglmmetadatadf<- rbind(correspondingms1, correspondingms1_LOH)
  }
  
  #CA65somatic:
  CA65somatic<-subset(metadatadf.0, ColonyName=="CAcolony65" & MutationClass=="SomaticMutation")
  CA65somaticmpandms<-subset(CA65somatic, sample==MutantParent1 | sample==MutantParent2 | sample==MutantSperm1 | sample==MutantSperm2)
  
  #CA65uglm:
  CA65uglm<-subset(metadatadf.0, ColonyName=="CAcolony65" & MutationClass=="UniqueGermlineMutation")
  CA65uglmmpandms<-subset(CA65uglm, sample==MutantParent1 | sample==MutantParent2 | sample==MutantSperm1 | sample==MutantSperm2)
  
  #metadatadf<-rbind( trueDeNovos_somatic, trueDeNovos_uglmmetadatadf, trueLoHp2, trueLoHsperm2,globalglm)#, trueLoHp2)
  somatic_metadatadf.0<-subset(metadatadf.0, MutationClass == "SomaticMutation")
  uniquesomatic_metadatadf.0<-somatic_metadatadf.0[match(unique(somatic_metadatadf.0$chrom.pos), somatic_metadatadf.0$chrom.pos),]
  
  
  
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
  
  
  threedfs<-list("somatic"=somaticmetadatadf, "uglm"=uglmmetadatadf, "gglm"=gglmmetadatadf, "metadatadf.0"=metadatadf.0)
  return(threedfs)
} 
scatterplot_func<-function(somaticmetadatadf) {
  #########comparison for all somatic mutations:##########
  somaticmetadatadf<-na.omit(somaticmetadatadf)
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
  
  df1<-data.frame(chrom=mutantparentdf$chrom, pos=mutantparentdf$pos, sample=mutantparentdf$sample, ref=mutantparentdf$ref, alt=mutantparentdf$alt, ParentAverage=parentsaverage, SpermAverage=spermaverage, TrueorFalse = as.factor(mutantparentdf$TrueorFalse), GoHorLoH = as.factor(mutantparentdf$GoH_or_LoH), Inheritance = as.factor(Inheritance), LaxInheritance = as.factor(LaxInheritance), s1=spermdf1$mutant_allele_depth, s2=spermdf2$mutant_allele_depth)
  df<-subset(df1, ParentAverage>=0.1)
  
  df_denovo<-subset(df, ParentAverage<1)
  otherframe<-data.frame(chrom.pos=mutantparentdf$chrom.pos, sampleID=mutantparentdf$sample, ref=mutantparentdf$ref, alt=mutantparentdf$alt,  chrom=mutantparentdf$chrom, pos=mutantparentdf$pos, sample=mutantparentdf$sample, ParentAverage=parentsaverage, SpermAverage=spermaverage, s1=spermdf1$mutant_allele_depth, s2=spermdf2$mutant_allele_depth, TrueorFalse = as.factor(mutantparentdf$TrueorFalse), GoHorLoH = as.factor(mutantparentdf$GoH_or_LoH),WhattoWhat=mutantparentdf$WhattoWhat)
  dfk<-subset(df_denovo, SpermAverage>0)
  #merged$TrueorFalse<-as.factor(merged$TrueorFalse)
  
  #parentspermcomparison<-lm(df_true$SpermAverage~df_true$ParentAverage)
  pall<-ggplot(df, aes(x=ParentAverage, y=SpermAverage, color=Inheritance)) + geom_point(size=2.5) +
    
    theme_bw() +
    ylim(0, 1) + xlim(0, 1) +
    #scale_color_manual(labels= c("Parent and Sperm","Parent Only"), values = c("red","black")) +
    scale_color_manual(labels= c("Parent and Sperm","Parent Only"),values = c("#999999","red")) +
    labs(color="SNV Type\n")+
    ylab("Average Variant Allele Frequency in the Mutant Sperm Pool") + xlab("Average Variant Allele Frequency in the Mutant Parent")+
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=15))+
    theme(legend.text=element_text(size=15))
  
  
  p<-ggplot(df_denovo, aes(x=ParentAverage, y=SpermAverage, color=Inheritance)) + geom_point(size=2.5) +
    #geom_smooth(method=lm, aes(group=1)) +
    theme_bw() +
    #stat_cor(label.x=0.2,label.y=0.45, aes(group=1)) +
    ylim(0, .4) + xlim(0, 0.5) +
    #scale_color_manual(values = c("red","black")) +
    scale_color_manual(values = c("#999999","lightgoldenrod1")) +
    ylab("Average Variant Allele Frequency in the Mutant Sperm Pool") + xlab("Average Variant Allele Frequency in the Mutant Parent")+
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=15))+
    theme(legend.text=element_text(size=15),legend.title=element_blank())
  
  pk<-ggplot(dfk, aes(x=ParentAverage, y=SpermAverage, color=Inheritance)) + geom_point(size=2.5) +
    geom_smooth(method=lm, aes(group=1)) +
    theme_bw() +
    stat_cor(label.x=0.05,label.y=0.5, aes(group=1, label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
    stat_regline_equation(label.x = .05, label.y = 0.55)+
    ylim(0, 1) + xlim(0, 1) +
    #scale_color_manual(values = c("red","black")) +
    scale_color_manual(labels= c("Parent and Sperm","Parent Only"),values = c("#999999","red")) +
    ylab("Average Variant Allele Frequency in the Mutant Sperm Pool") + xlab("Average Variant Allele Frequency in the Mutant Parent")+
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=15))+
    geom_abline(intercept=0,slope=1)+
    theme(legend.text=element_blank(),legend.title=element_blank())
  pk_a<-ggplot(dfk, aes(x=ParentAverage, y=SpermAverage, color=Inheritance)) + geom_point(size=2.5) +
    geom_smooth(method=lm, aes(group=1)) +
    theme_bw() +
    stat_cor(label.x=0.05,label.y=0.35, aes(group=1, label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
    stat_regline_equation(label.x = .05, label.y = 0.3)+
    ylim(0, 0.4) + xlim(0, 0.4) +
    #scale_color_manual(values = c("red","black")) +
    scale_color_manual(labels= c("Parent and Sperm","Parent Only"),values = c("#999999","red")) +
    ylab("Average Variant Allele Frequency in the Mutant Sperm Pool") + xlab("Average Variant Allele Frequency in the Mutant Parent")+
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=15))+
    geom_abline(intercept=0,slope=1)+
    theme(legend.text=element_blank(),legend.title=element_blank())
  
  pk2<-ggplot(dfk, aes(x=ParentAverage, y=SpermAverage, color=Inheritance)) + geom_point(size=2.5) +
    geom_smooth(method=lm, aes(group=1)) +
    theme_bw() +
    stat_cor(label.x=0.5,label.y=0.5, aes(group=1)) +
    stat_regline_equation(label.x = .5, label.y = 0.3)+
    ylim(0, 1) + xlim(0.5, 1) +
    scale_color_manual(values = c("red","black")) +
    ylab("Average Variant Allele Frequency in the Mutant Sperm Pool") + xlab("Average Variant Allele Frequency in the Mutant Parent")+
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=15))+
    geom_abline(intercept=0,slope=1)+
    theme(legend.text=element_blank(),legend.title=element_blank())
  
  #geom_text(x = 25, y = 300, label = lm_eqn(df_true), parse = TRUE)
  #return(p)
  #return(parentspermcomparison)
  #abline(parentspermcomparison, lwd=2)
  #cor<-cor.test(df_true$SpermAverage, df_true$ParentAverage)
  #return(summary(parentspermcomparison)$r.squared)
  newlist<-list(p, cor)
  #parentspermcomparison
  #return(newlist)
  #return(pall | pk ) #this is supplementary figure 1
  #return(pk | pk_a) #for defense ppt
  return(df) #use this for the rest of the analyses
  #return(otherframe)
  #return(pall)
  #return(p)
  #dndscvdf<- data.frame("sampleID"= mutantparentdf$sample,"chr"= mutantparentdf$chrom,"pos"= mutantparentdf$pos, "ref" = mutantparentdf$ref, "alt"= uniquemetadatadf$alt) # do not use ref and alt!
  
}
Ahya_fasta<-read.fasta("~/Documents/GitHub/CoralGermline/mappedtoahya/jasmine-sta1387-mb-hirise-eyx9u_12-15-2020__final_assembly.fasta",forceDNAtolower = FALSE)
spectrum<- function(dataframe, dataframecount) {
  #Ahya_fasta<-read.fasta("~/Documents/GitHub/CoralGermline/dndscv/Amil.v2.01.chrs.fasta",forceDNAtolower = FALSE)
  WhattoWhat<-rep("Placeholder",nrow(dataframe))
  TiorTv<-rep("Placeholder",nrow(dataframe))
  
  ref<-as.character(dataframe$ref)
  alt<-as.character(dataframe$alt)
  for (i in 1:nrow(dataframe)) {
    if (ref[i] == "C" & alt[i] == "A") {
      WhattoWhat[i] <- 'CtoA'
      TiorTv[i]<-"Tv"
      #gh[i] <- 'CtoA'
      #print("ya")
      #print(i)
      #print(AtG[i])
    } else if (ref[i] == "A" & alt[i] == "C")  {
      WhattoWhat[i] <- "AtoC"
      TiorTv[i]<-"Tv"
    } else if (ref[i] == "A" & alt[i] == "G") {
      WhattoWhat[i] <- "AtoG"
      TiorTv[i]<-"Ti"
      
    } else if (ref[i] == "G" & alt[i] == "A") {
      WhattoWhat[i] <- 'GtoA'
      TiorTv[i]<-"Ti"
    } else if (ref[i] == "G" & alt[i] == "T") {
      WhattoWhat[i] <- "GtoT"
      TiorTv[i]<-"Tv"
    } else if (ref[i] == "T" & alt[i] == "G") {
      WhattoWhat[i] <- "TtoG"
      TiorTv[i]<-"Tv"
    } else if (ref[i] == "T" & alt[i] == "C") {
      WhattoWhat[i] <- "TtoC"
      TiorTv[i]<-"Ti"
    } else if (ref[i] == "C" & alt[i] == "T") {
      WhattoWhat[i] <-"CtoT"
      TiorTv[i]<-"Ti"
    } else if (ref[i] == "C" & alt[i] == "G") {
      WhattoWhat[i] <- "CtoG"
      TiorTv[i]<-"Tv"
    } else if (ref[i] == "G" & alt[i] == "C") {
      WhattoWhat[i] <-"GtoC"
      TiorTv[i]<-"Tv"
    } else if (ref[i] == "T" & alt[i] == "A") {
      WhattoWhat[i] <- "TtoA"
      TiorTv[i]<-"Tv"
    } else if (ref[i] == "A" & alt[i] == "T") {
      WhattoWhat[i] <-"AtoT"
      TiorTv[i]<-"Tv"
    }
  }
  
  together<-cbind(dataframe, WhattoWhat,TiorTv)				 
  togethercount<-nrow(together)
  
  verATGClist<-together[ which(together$WhattoWhat=="AtoG" | together$WhattoWhat=="TtoC"),]
  verATGCcount<-nrow(verATGClist)
  verATGCprop<-verATGCcount/togethercount
  SEverATGCprop<-se(verATGCprop)
  
  verGCATlist<-together[ which(together$WhattoWhat=="GtoA" | together$WhattoWhat=="CtoT"),]
  verGCATcount<-nrow(verGCATlist)
  verGCATprop<-verGCATcount/togethercount
  
  CpGdf<-data.frame(Date=as.Date(character()))
  CpHdf<-data.frame(Date=as.Date(character()))
  all<-NULL
  for (i in 1:nrow(together)) {
    if (together$WhattoWhat[i] == "CtoT") {
      character<-as.character(together$pos[i])
      integer<-as.integer(character)
      next_position<-integer + 1
      chr_num<-as.character(together$chrom[i])
      
      chroms<-Ahya_fasta[[chr_num]]
      ref_nuc<- chroms[next_position]
      ref_nuc<-toupper(ref_nuc)
      allnucs<-data.frame(ref_nuc)
      all<-rbind(all,allnucs)
      if (ref_nuc == "G") {
        #print(ref_nuc)
        Gnucs<-data.frame(ref_nuc)
        CpGdf<-rbind(CpGdf, Gnucs)
      } else if (ref_nuc == "A"| ref_nuc == "C"| ref_nuc == "T") {
        othernucs<-data.frame(ref_nuc)
        CpHdf<-rbind(CpHdf,othernucs)
      }
    } else if (together$WhattoWhat[i] == "GtoA") {   
      character<-as.character(together$pos[i])
      integer<-as.integer(character)
      next_position<-integer + 1
      chr_num<-as.character(together$chrom[i])
      
      chroms<-Ahya_fasta[[chr_num]]
      ref_nuc<- chroms[next_position]
      ref_nuc<-toupper(ref_nuc)
      allnucs<-data.frame(ref_nuc)
      all<-rbind(all,allnucs)
      
      if (ref_nuc == "C") {
        #print(ref_nuc)
        Gnucs<-data.frame(ref_nuc)
        CpGdf<-rbind(CpGdf, Gnucs)
      } else if (ref_nuc == "A"| ref_nuc == "G"| ref_nuc == "T") {
        othernucs<-data.frame(ref_nuc)
        CpHdf<-rbind(CpHdf,othernucs)
      }
    }
  }    
  CpGcount<-nrow(CpGdf)
  CpHcount<-nrow(CpHdf)
  
  CpGprop<- CpGcount/togethercount
  CpHprop<- CpHcount/togethercount
  
  verATTAlist<-together[ which(together$WhattoWhat=="AtoT" | together$WhattoWhat=="TtoA"),]
  verATTAcount<-nrow(verATTAlist)
  verATTAprop<-verATTAcount/togethercount
  
  verACTGlist<-together[ which(together$WhattoWhat=="AtoC" | together$WhattoWhat=="TtoG"),]
  verACTGcount<-nrow(verACTGlist)
  verACTGprop<-verACTGcount/togethercount
  
  verGCCGlist<-together[ which(together$WhattoWhat=="CtoG" | together$WhattoWhat=="GtoC"),]
  verGCCGcount<-nrow(verGCCGlist)
  verGCCGprop<-verGCCGcount/togethercount
  
  verGCTAlist<-together[ which(together$WhattoWhat=="GtoT" | together$WhattoWhat=="CtoA"),]
  verGCTAcount<-nrow(verGCTAlist)
  verGCTAprop<-verGCTAcount/togethercount
  
  binomATGC<-dbinom(verATGCcount, togethercount, verATGCprop)
  sdATGC<-togethercount*verATGCprop*(1-verATGCprop)
  seATGC<- sdATGC/sqrt(togethercount)
  
  binomGCAT<-dbinom(verGCATcount, togethercount, verGCATprop)
  sdGCAT<-togethercount*verGCATprop*(1-verGCATprop)
  seGCAT<- sdGCAT/sqrt(togethercount)
  
  binomCpG<-dbinom(CpGcount, togethercount, CpGprop)
  sdCpG<-togethercount*CpGprop*(1-CpGprop)
  seCpG<-sdCpG/sqrt(togethercount)
  
  binomCpH<-dbinom(CpHcount, togethercount, CpHprop)
  sdCpH<-togethercount*CpHprop*(1-CpHprop)
  seCpH<-sdCpH/sqrt(togethercount)
  
  binomATTA<-dbinom(verATTAcount, togethercount, verATTAprop)
  sdATTA<-togethercount*verATTAprop*(1-verATTAprop)
  seATTA<- sdATTA/sqrt(togethercount)
  
  binomACTG<-dbinom(verACTGcount, togethercount, verACTGprop)
  sdACTG<-togethercount*verACTGprop*(1-verACTGprop)
  seACTG<- sdACTG/sqrt(togethercount)
  
  binomGCCG<-dbinom(verGCCGcount, togethercount, verGCCGprop)
  sdGCCG<-togethercount*verGCCGprop*(1-verGCCGprop)
  seGCCG<- sdGCCG/sqrt(togethercount)
  
  binomGCTA<-dbinom(verGCTAcount, togethercount, verGCTAprop)
  sdGCTA<-togethercount*verGCTAprop*(1-verGCTAprop)
  seGCTA<- sdGCTA/sqrt(togethercount)
  
  standarderrors<-c(seATGC, seCpG, seCpH, seATTA, seACTG, seGCCG, seGCTA)
  se_forprops<-standarderrors/togethercount
  vertypesDF<-data.frame(Types=c("A>G/T>C","CpG", "CpH", "A>T/T>A","A>C/T>G","G>C/C>G","G>T/C>A"), coralProportion=c(verATGCprop, CpGprop, CpHprop, verATTAprop,verACTGprop, verGCCGprop, verGCTAprop),se_forprops)
  #pmdf$MutType<- factor(pmdf$MutType,levels=c("synonymous","missense","nonsense"))
  vertypesDF$Types<-factor(vertypesDF$Types,levels=c("A>G/T>C","CpG", "CpH", "A>T/T>A","A>C/T>G","G>C/C>G","G>T/C>A"))
  #MUTATION SPECTRUM PLOT:
  #vertypesDFplot<-barplot(vertypesDF$coralProportion,names.arg=vertypesDF$Types, ylim=c(0,0.7), main="allmuts", las=2) #Figure 4a
  z<- ggplot(vertypesDF, aes(x=Types, y=coralProportion)) +
    #scale_fill_manual(value=c("white"))+
    geom_bar(stat="identity", fill="white", color="blue") +
    #geom_point(size=3)+
    #geom_boxplot()+
    ylim(0, 0.35)+
    theme_bw()+
    ylab("Proportion of SNVs") +
    theme(axis.text.x = element_text(angle = 90, size=18), axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+#, axis.text.x=element_text(size=25))+
    geom_errorbar(aes(ymin=coralProportion-(se_forprops), ymax=coralProportion+(se_forprops)), width=0.4)#, width=.2,
  
  Transtions<- verATGCprop + verGCATprop
  Transversions<- verATTAprop+verACTGprop+verGCCGprop+verGCTAprop
  TiTv<- Transtions/Transversions
  counts<-c(verATGCcount,CpGcount, CpHcount, verATTAcount, verACTGcount, verGCCGcount, verGCTAcount)
  #return(z) #this is fig: supplementary fig 2
  #return(list(vertypesDF, TiTv))
  #return(list(z, vertypesDF))
  #return(z)
  return(counts) #this is for supplementary stats (chisqs)
}

