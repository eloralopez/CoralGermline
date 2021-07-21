#USE ME JULY 2021#######
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
#denominators:
#CA56denom<-79943111-(895*14)#11/17  #80718433-(892*14) #9/14
#CA60denom<-   51791529-(895*14)#11/17    #52564929-(892*14) #9/14
#CA65denom<-   47410341-(895*14)#11/17    #48319810-(892*14) #9/14
#CA60denom<- 76882242-(895*14)#12/3
#CA56denom_coding<- 14123208
#CA60denom_coding<- 3823690+9622397 #12/3
#CA65denom_coding<- 8462566
#coding denom mapped to ahya 6/2021:
CA56denom_coding<- 13167782
CA60denom_coding<- 13090241
CA65denom_coding<- 9105020
#DENOMINATORS MAPPED TO AHYA
CA56denom<-118928278-(14*949)
CA60denom<-118886418-(14*949)
CA65denom<-74856149-159102+5316571-(14*949)
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
    trueLoH_somatic_mpandms<-subset(trueLoH_somatic_mpandms, chrom.pos != "chr9.10238326")
    trueLoH_somatic_mpandms<-subset(trueLoH_somatic_mpandms, chrom.pos != "chr13.20768582")
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
#threedfs_CA56<-threedfs_func(list.files(path="~/Documents/GitHub/CoralGermline/cleanpipeline", pattern="*56_dm_20201210.txt.ann.txt", full.names=T, recursive=FALSE))
#ahya ref:
#threedfs_CA56<-threedfs_func(list.files(path="~/Documents/GitHub/CoralGermline/mappedtoahya", pattern="*56_dm_20210107.txt", full.names=T, recursive=FALSE))
#ahya ref WITH ANNOTATIONS:
threedfs_CA56<-threedfs_func(list.files(path="~/Documents/GitHub/CoralGermline/datafiles", pattern="*56_dm_20210107_chrs.txt.ann.txt", full.names=T, recursive=FALSE))
somaticCA56df<-threedfs_CA56$somatic
uglmCA56df<-threedfs_CA56$uglm
uniqueuglmCA56df<-uglmCA56df[match(unique(uglmCA56df$chrom.pos), 					uglmCA56df$chrom.pos),]
gglmCA56df<-threedfs_CA56$gglm
metadatadf.0CA56<-threedfs_CA56$metadatadf.0

#threedfs_CA60<-threedfs_func(list.files(path="~/Documents/GitHub/CoralGermline/cleanpipeline", pattern="*60_dm_20201210.txt.ann.txt", full.names=T, recursive=FALSE))
#threedfs_CA60<-threedfs_func(list.files(path="~/Documents/GitHub/CoralGermline/mappedtoahya", pattern="*60_dm_20210107.txt", full.names=T, recursive=FALSE))
threedfs_CA60<-threedfs_func(list.files(path="~/Documents/GitHub/CoralGermline/datafiles", pattern="*60_dm_20210107_chrs.txt.ann.txt", full.names=T, recursive=FALSE))

somaticCA60df<-threedfs_CA60$somatic
uglmCA60df<-threedfs_CA60$uglm
uniqueuglmCA60df<-uglmCA60df[match(unique(uglmCA60df$chrom.pos), 					uglmCA60df$chrom.pos),]
gglmCA60df<-threedfs_CA60$gglm
metadatadf.0CA60<-threedfs_CA60$metadatadf.0

#threedfs_CA65<-threedfs_func(list.files(path="~/Documents/GitHub/CoralGermline/cleanpipeline", pattern="*65_dm_20201210.txt.ann.txt", full.names=T, recursive=FALSE))
#threedfs_CA65<-threedfs_func(list.files(path="~/Documents/GitHub/CoralGermline/mappedtoahya", pattern="*65_dm_20210107.txt", full.names=T, recursive=FALSE))
threedfs_CA65<-threedfs_func(list.files(path="~/Documents/GitHub/CoralGermline/datafiles", pattern="*65_dm_20210326_chrs.txt.ann.txt", full.names=T, recursive=FALSE))

somaticCA65df<-threedfs_CA65$somatic
uglmCA65df<-threedfs_CA65$uglm
uniqueuglmCA65df<-uglmCA65df[match(unique(uglmCA65df$chrom.pos), 					uglmCA65df$chrom.pos),]
gglmCA65df<-threedfs_CA65$gglm
metadatadf.0CA65<-threedfs_CA65$metadatadf.0

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
    scale_color_manual(labels= c("Parent and Sperm","Parent Only"), values = c("red","black")) +
    labs(color="SNV Type\n")+
    ylab("Average Variant Allele Frequency in the Mutant Sperm Pool") + xlab("Average Variant Allele Frequency in the Mutant Parent")+
    theme(axis.text=element_text(size=25),
          axis.title=element_text(size=25))+
    theme(legend.text=element_text(size=25),legend.title=element_text(size=25))

  
  p<-ggplot(df_denovo, aes(x=ParentAverage, y=SpermAverage, color=Inheritance)) + geom_point(size=2.5) +
    #geom_smooth(method=lm, aes(group=1)) +
    theme_bw() +
    #stat_cor(label.x=0.2,label.y=0.45, aes(group=1)) +
    ylim(0, .4) + xlim(0, 0.5) +
    scale_color_manual(values = c("red","black")) +
    ylab("Average Variant Allele Frequency in the Mutant Sperm Pool") + xlab("Average Variant Allele Frequency in the Mutant Parent")+
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=15))+
    theme(legend.text=element_text(size=15),legend.title=element_blank())

  pk<-ggplot(dfk, aes(x=ParentAverage, y=SpermAverage, color=Inheritance)) + geom_point(size=2.5) +
    geom_smooth(method=lm, aes(group=1)) +
    theme_bw() +
    #stat_cor(label.x=0.05,label.y=0.35, aes(group=1)) +
    #stat_regline_equation(label.x = .05, label.y = 0.4)+
    ylim(0, 0.5) + xlim(0, .5) +
    scale_color_manual(values = c("red","black")) +
    ylab("Average Variant Allele Frequency in the Mutant Sperm Pool") + xlab("Average Variant Allele Frequency in the Mutant Parent")+
    theme(axis.text=element_text(size=25),
          axis.title=element_text(size=25))+
    geom_abline(intercept=0,slope=1)+
    theme(legend.position="none")
  
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
  #return(pall +labs(tag="a.")+theme(plot.tag=element_text(face="bold",size=30)) | pk +labs(tag="b.")+
  #         theme(plot.tag=element_text(face="bold",size=30))) #this is supplementary figure 1
  
  return(df) #use this for the rest of the analyses
  #return(otherframe)
  #return(pall)
  #return(p)
  #dndscvdf<- data.frame("sampleID"= mutantparentdf$sample,"chr"= mutantparentdf$chrom,"pos"= mutantparentdf$pos, "ref" = mutantparentdf$ref, "alt"= uniquemetadatadf$alt) # do not use ref and alt!
  
}
allcolonies<-scatterplot_func(rbind(somaticCA56df,somaticCA60df,somaticCA65df))
allcolonies_inherited_GOH<-subset(allcolonies, Inheritance=="inherited" & GoHorLoH=="DeNovo")
allcolonies_notinherited_GOH<-subset(allcolonies, Inheritance=="notinherited" & GoHorLoH=="DeNovo")
VAFs<-data.frame("PVAF"=c(allcolonies_inherited_GOH$ParentAverage, allcolonies_notinherited_GOH$ParentAverage), 
           "shared"=c(rep("Parent and Sperm",nrow(allcolonies_inherited_GOH)),rep("Parent Only",nrow(allcolonies_notinherited_GOH))))
VAFmeans<-c(mean(subset(VAFs, shared=="Parent and Sperm")$PVAF),mean(subset(VAFs, shared=="Parent Only")$PVAF))
mu<-data.frame(VAFmeans, shared=c("Parent and Sperm","Parent Only"))
SuppFig<-ggplot(data=VAFs, aes(x=PVAF*100, fill=shared, ..density..)) +
  geom_density(alpha=0.5) +
  #geom_histogram()+
  scale_fill_manual(values = c("#999999","goldenrod")) +
  #xlim(3.9,8.5) +
  ylab("Probability Density Function")+ xlab("Average Variant Allele Frequency in the Mutant Parent (as %)")+
  geom_vline(data=mu, aes(xintercept=VAFmeans*100, color=shared),
             
             linetype="dashed",size=2)+
  scale_color_manual(values = c("#999999","goldenrod")) +
  theme_bw() +
  theme(axis.text=element_text(size=25), 
        axis.text.x=element_text(size=25, angle = 0),
        axis.title=element_text(size=25))+
  theme(legend.text=element_text(size=25), legend.title = element_blank())

SuppFig_hist<-ggplot(data=VAFs, aes(x=PVAF, fill=shared)) +
  #geom_density(alpha=0.3) +
  geom_histogram(binwidth=0.01, closed = TRUE)+
  scale_fill_manual(values = c("#999999","lightgoldenrod1")) +
  #xlim(3.9,8.5) +
  geom_vline(data=mu, aes(xintercept=VAFmeans, color=shared),
             linetype="dashed",size=2)+
  scale_color_manual(values = c("#999999","lightgoldenrod1")) +
  theme_bw() 
SuppFig+scale_y_continuous(expand = c(0,0.01))+theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
SuppFig+ geom_line(data = VAFs, aes(x = x, y = y), color = "red")
VAFs<-data.frame("PVAF"=c(allcolonies_inherited_GOH$ParentAverage, allcolonies_notinherited_GOH$ParentAverage), 
                 "shared"=c(rep("Parent and Sperm",nrow(allcolonies_inherited_GOH)),rep("Parent Only",nrow(allcolonies_notinherited_GOH))))
plot(density(subset(VAFs, shared=="Parent Only")$PVAF))
ggplot(data=VAFs, aes(x=PVAF, fill=shared)) +
  geom_density(alpha=0.5) +
  #xlim(0.0,0.5) +
  #geom_histogram(binwidth=0.05, alpha=0.5, nrow(POsition="identity") +
  theme_minimal()
PO15<-nrow(subset(VAFs, shared=="Parent Only" & (PVAF <0.15) & (PVAF > 0)))
     PO20<-nrow(subset(VAFs, shared=="Parent Only" & (PVAF >=0.15) & (PVAF < .2)))
          PO25<-nrow(subset(VAFs, shared=="Parent Only" & (PVAF >=0.2) & (PVAF < .25)))
               PO30<-nrow(subset(VAFs, shared=="Parent Only" & (PVAF >=0.25) & (PVAF < .3)))
                    PO35<-nrow(subset(VAFs, shared=="Parent Only" & (PVAF >=0.3) & (PVAF < .35)))
                         PO40<-nrow(subset(VAFs, shared=="Parent Only" & (PVAF >=0.35) & (PVAF < .40)))
                              PO45<-nrow(subset(VAFs, shared=="Parent Only" & (PVAF >=0.40) & (PVAF < .45)))
                                   POcounts<-c(PO15, PO20, PO25, PO30, PO35, PO40, PO45)
                                                                                                                PS15<-nrow(subset(VAFs, shared=="Parent and Sperm" & (PVAF <0.15) & (PVAF > 0)))
                                                                                                                PS20<-nrow(subset(VAFs, shared=="Parent and Sperm" & (PVAF >=0.15) & (PVAF < .2)))
                                                                                                                PS25<-nrow(subset(VAFs, shared=="Parent and Sperm" & (PVAF >=0.2) & (PVAF < .25)))
                                                                                                                PS30<-nrow(subset(VAFs, shared=="Parent and Sperm" & (PVAF >=0.25) & (PVAF < .3)))
                                                                                                                PS35<-nrow(subset(VAFs, shared=="Parent and Sperm" & (PVAF >=0.3) & (PVAF < .35)))
                                                                                                                PS40<-nrow(subset(VAFs, shared=="Parent and Sperm" & (PVAF >=0.35) & (PVAF < .40)))
                                                                                                                PS45<-nrow(subset(VAFs, shared=="Parent and Sperm" & (PVAF >=0.40) & (PVAF < .45)))
                                                                                                                
PScounts<-c(PS15, PS20, PS25, PS30, PS35, PS40, PS45)
PO_v_PS<-data.frame(POcounts,PScounts)
chisq.test(PO_v_PS)
CA56_<-scatterplot_func(somaticCA56df)
CA56_CODING<-subset(CA56_, MutationType=="missense_variant" | MutationType=="synonymous_variant"| MutationType=="synonymous_variant")
CA56_missense<-subset(CA56_CODING,MutationType=="missense_variant")
CA56_syn<-subset(CA56_CODING,MutationType=="synonymous_variant")
CA56_inherited<-subset(CA56_, Inheritance=="inherited")
CA56_notinherited<-subset(CA56_, Inheritance=="notinherited")
CA56_loh<-subset(CA56_, GoHorLoH=="LoH")
CA56_denovo<-subset(CA56_, GoHorLoH=="DeNovo")

inheriteddenovoCA56_<-subset(CA56_, GoHorLoH=="DeNovo" & (s1>0 & s2>0))
inheritedlohCA56_<-subset(CA56_, GoHorLoH=="LoH" & SpermAverage==1 & ParentAverage==1)
notinheriteddenovoCA56_<-subset(CA56_, GoHorLoH=="DeNovo" & (s1==0 | s2==0))
notinheritedlohCA56_<-subset(CA56_, GoHorLoH=="LoH" & SpermAverage<1)
inheritedcountCA56_<-nrow(inheriteddenovoCA56_)+nrow(inheritedlohCA56_)
inheritedpropCA56_<- inheritedcountCA56_/nrow(CA56_)

CAP6<-subset(CA56_, sample=="CAP6-1_S47") # when scatterplot_func returns "otherframe"
CAP8<-subset(CA56_, sample=="CAP8-1_S40")
CAP12<-subset(CA56_,sample=="CAP12-1_S46")
CAP6i<-subset(CAP6, Inheritance=="inherited")
CAP8i<-subset(CAP8, Inheritance=="inherited")
CAP6n<-subset(CAP6, Inheritance=="notinherited")
CAP8n<-subset(CAP8, Inheritance=="notinherited")

CA60_<-scatterplot_func(somaticCA60df)
CA60_CODING<-subset(CA60_, MutationType=="missense_variant" | MutationType=="synonymous_variant"| MutationType=="synonymous_variant")
CA60_loh<-subset(CA60_, GoHorLoH=="LoH")
CA60_denovo<-subset(CA60_, GoHorLoH=="DeNovo")

CA60_inherited<-subset(CA60_, Inheritance=="inherited")
CA60_notinherited<-subset(CA60_, Inheritance=="notinherited")
inheriteddenovoCA60_<-subset(CA60_, SpermAverage>0 &ParentAverage>=0.1 &ParentAverage<1)
inheritedlohCA60_<-subset(CA60_, SpermAverage==1 & ParentAverage==1)
notinheriteddenovoCA60_<-subset(CA60_, SpermAverage==0)
notinheritedlohCA60_<-subset(CA60_, ParentAverage==1 & SpermAverage<1)						  	
inheritedcountCA60_<-nrow(inheriteddenovoCA60_)+nrow(inheritedlohCA60_)
CAP22<-subset(CA60_, sample=="CAP22-1_S43")
CAP23<-subset(CA60_, sample=="CAP23-1_S45")
CAP24<-subset(CA60_,sample=="CAP24-1_S42")
CAP26<-subset(CA60_,sample=="CAP26-1_S39")
CAP22i<-subset(CAP22, Inheritance=="inherited")
CAP23i<-subset(CAP23, Inheritance=="inherited")
CAP26i<-subset(CAP26, Inheritance=="inherited")
CAP22n<-subset(CAP22, Inheritance=="notinherited")
CAP23n<-subset(CAP23, Inheritance=="notinherited")
CAP26n<-subset(CAP26, Inheritance=="notinherited")


CA65_<-scatterplot_func(somaticCA65df)
CA65_CODING<-subset(CA65_, MutationType=="missense_variant" | MutationType=="synonymous_variant"| MutationType=="synonymous_variant")
CA65_loh<-subset(CA65_, GoHorLoH=="LoH")
CA65_denovo<-subset(CA65_, GoHorLoH=="DeNovo")

CA65_inherited<-subset(CA65_, Inheritance=="inherited")
CA65_notinherited<-subset(CA65_, Inheritance=="notinherited")
inheriteddenovoCA65_<-subset(CA65_, SpermAverage>0 &ParentAverage>=0.1 &ParentAverage<1)
inheritedlohCA65_<-subset(CA65_, SpermAverage==1 & ParentAverage==1)
notinheriteddenovoCA65_<-subset(CA65_, SpermAverage==0)
notinheritedlohCA65_<-subset(CA65_, ParentAverage==1 & SpermAverage<1)
inheritedcountCA65_<-nrow(inheriteddenovoCA65_)+nrow(inheritedlohCA65_)

CAP10<-subset(CA65_, sample=="CAP10-1_S44")
CAP11<-subset(CA65_, sample=="CAP11-1_S41")
CAP9<-subset(CA65_,sample=="CAP9-1_S48")
CAP11i<-subset(CAP11, Inheritance=="inherited")
CAP9i<-subset(CAP9, Inheritance=="inherited")
CAP11n<-subset(CAP11, Inheritance=="notinherited")
CAP9n<-subset(CAP9, Inheritance=="notinherited")

GOH_func<-function(CAP) {
  CAP_GOH<-subset(CAP, Inheritance=="inherited")
  return(CAP_GOH)
}
LOH_func<-function(CAP) {
  CAP_LOH<-subset(CAP, Inheritance=="inherited")
  return(CAP_LOH)
}

CAP11i_GOH<-GOH_func(CAP11i)
CAP11i_LOH<-LOH_func(CAP11i)


totalCA56_<-inheritedcountCA56_+ nrow(uniqueuglmCA56df)+ nrow(gglmCA56df)
totalCA60_<- inheritedcountCA60_+ nrow(uniqueuglmCA60df)+ nrow(gglmCA60df)
totalCA65_<- inheritedcountCA65_+ nrow(uniqueuglmCA65df)+ nrow(gglmCA65df)

CA56i<-c(nrow(subset(CAP6, Inheritance=="inherited")), nrow(subset(CAP8, Inheritance=="inherited")))
#CA56ip<-CA56i/(CA56i+CA56uglm+CA56gglm)
CA60i<-c(nrow(subset(CAP22, Inheritance=="inherited")),nrow(subset(CAP23, Inheritance=="inherited")),nrow(subset(CAP26, Inheritance=="inherited")))
#CA60ip<-CA60i/(CA60i+CA60uglm+CA60gglm)
CA65i<-c(nrow(subset(CAP11, Inheritance=="inherited")),nrow(subset(CAP9, Inheritance=="inherited")))
#CA65ip<-CA65i/(CA65i+CA65uglm+CA65gglm)
CA56n<-c(nrow(subset(CAP6, Inheritance=="notinherited")), nrow(subset(CAP8, Inheritance=="notinherited")))
#CA56np<-CA56n/(CA56n+CA56uglm+CA56gglm)
CA60n<-c(nrow(subset(CAP22, Inheritance=="notinherited")),nrow(subset(CAP23, Inheritance=="notinherited")),nrow(subset(CAP26, Inheritance=="notinherited")))
#CA60np<-CA60n/(CA60n+CA60uglm+CA60gglm)
CA65n<-c(nrow(subset(CAP11, Inheritance=="notinherited")),nrow(subset(CAP9, Inheritance=="notinherited")))
#CA65np<-CA65n/(CA65n+CA65uglm+CA65gglm)
CAPn<-c(CA56n, CA60n, CA65n)
CAPi<-c(CA56i, CA60i, CA65i)
CA56uglm<-c(nrow(subset(uniqueuglmCA56df, startsWith(uniqueuglmCA56df$MutantParent1,'CAP6')==TRUE ) ),
            nrow(subset(uniqueuglmCA56df, startsWith(uniqueuglmCA56df$MutantParent1,'CAP8')==TRUE )))
CAS6uglm<-nrow(subset(uniqueuglmCA56df, startsWith(uniqueuglmCA56df$MutantParent1,'CAP6')==TRUE ) )
nrow(CAS6uglm)
CAS6uglm<-subset(uniqueuglmCA56df, startsWith(uniqueuglmCA56df$MutantParent1,'CAP6')==TRUE ) 

CAS8uglm<-subset(uniqueuglmCA56df, startsWith(uniqueuglmCA56df$MutantParent1,'CAP8')==TRUE )
#CA56uglmp<-CA56uglm/(CA56i+CA56uglm+CA56gglm)
CA60uglm<-c( nrow(subset(uniqueuglmCA60df, startsWith(uniqueuglmCA60df$MutantParent1,'CAP22')==TRUE ) ),
             nrow(subset(uniqueuglmCA60df, startsWith(uniqueuglmCA60df$MutantParent1,'CAP23')==TRUE ) ),
             nrow(subset(uniqueuglmCA60df, startsWith(uniqueuglmCA60df$MutantParent1,'CAP26')==TRUE )))#,
#nrow(subset(uniqueuglmCA60, startsWith(uniqueuglmCA60$MutantParent1,'CAP26')==TRUE ))
CAS22uglm<-subset(uniqueuglmCA60df, startsWith(uniqueuglmCA60df$MutantParent1,'CAP22')==TRUE ) 
CAS23uglm<-subset(uniqueuglmCA60df, startsWith(uniqueuglmCA60df$MutantParent1,'CAP23')==TRUE ) 
CAS26uglm<-subset(uniqueuglmCA60df, startsWith(uniqueuglmCA60df$MutantParent1,'CAP26')==TRUE )
#CA60uglmp<-CA60uglm/(CA60i+CA60uglm+CA60gglm)
CA65uglm<- c(
  nrow(subset(uniqueuglmCA65df, startsWith(uniqueuglmCA65df$MutantParent1,'CAP11')==TRUE ) ),
  nrow(subset(uniqueuglmCA65df, startsWith(uniqueuglmCA65df$MutantParent1,'CAP9')==TRUE )))
CAS9uglm<-subset(uniqueuglmCA65df, startsWith(uniqueuglmCA65df$MutantParent1,'CAP9')==TRUE ) 
CAS11uglm<-subset(uniqueuglmCA65df, startsWith(uniqueuglmCA65df$MutantParent1,'CAP11')==TRUE ) 
#CA65uglmp<-CA65uglm/(CA65i+CA65uglm+CA65gglm)
CA56gglm<-c(
  nrow(subset(gglmCA56df, startsWith(gglmCA56df$sample,'CAS6')==TRUE ) ),
  nrow(subset(gglmCA56df, startsWith(gglmCA56df$sample,'CAS8')==TRUE )))
#CA56gglmp<-CA56gglm/(CA56i+CA56uglm+CA56gglm)
CA60gglm<-c(nrow(subset(gglmCA60df, startsWith(gglmCA60df$sample,'CAS22')==TRUE ) ),
             nrow(subset(gglmCA60df, startsWith(gglmCA60df$sample,'CAS23')==TRUE ) ),
             nrow(subset(gglmCA60df, startsWith(gglmCA60df$sample,'CAS26')==TRUE )))#,
#nrow(subset(gglmCA60, startsWith(gglmCA60$sample,'CAS26')==TRUE )))
#CA60gglmp<-CA60gglm/(CA60i+CA60uglm+CA60gglm)
CA65gglm<- c(
  nrow(subset(gglmCA65df, startsWith(gglmCA65df$sample,'CAS11')==TRUE ) ),
  nrow(subset(gglmCA65df, startsWith(gglmCA65df$sample,'CAS9')==TRUE )))

##GOH and LOH proportions
CAP6inh<-nrow(subset(CAP6, Inheritance=="inherited"))
CAP6inhLOH<-nrow(subset(CAP6, Inheritance=="inherited" & GoHorLoH=="LoH"))
CAP6inhGOH<-CAP6inh-CAP6inhLOH
CAP6notinh<-nrow(CAP6)-CAP6inh
CAP6notinhLOH<-nrow(subset(CAP6, Inheritance=="notinherited" & GoHorLoH=="LoH"))
CAP6notinhGOH<-CAP6notinh-CAP6notinhLOH
CAP8inh<-nrow(subset(CAP8, Inheritance=="inherited"))
CAP8inhLOH<-nrow(subset(CAP8, Inheritance=="inherited" & GoHorLoH=="LoH"))
CAP8inhGOH<-CAP8inh-CAP8inhLOH
CAP8notinh<-nrow(CAP8)-CAP8inh
CAP8notinhLOH<-nrow(subset(CAP8, Inheritance=="notinherited" & GoHorLoH=="LoH"))
CAP8notinhGOH<-CAP8notinh-CAP8notinhLOH
CAP22inh<-nrow(subset(CAP22, Inheritance=="inherited"))
CAP22inhLOH<-nrow(subset(CAP22, Inheritance=="inherited" & GoHorLoH=="LoH"))
CAP22inhGOH<-CAP22inh-CAP22inhLOH
CAP22notinh<-nrow(CAP22)-CAP22inh
CAP22notinhLOH<-nrow(subset(CAP22, Inheritance=="notinherited" & GoHorLoH=="LoH"))
CAP22notinhGOH<-CAP22notinh-CAP22notinhLOH

CAP23inh<-nrow(subset(CAP23, Inheritance=="inherited"))
CAP23inhLOH<-nrow(subset(CAP23, Inheritance=="inherited" & GoHorLoH=="LoH"))
CAP23inhGOH<-CAP23inh-CAP23inhLOH
CAP23notinh<-nrow(CAP23)-CAP23inh
CAP23notinhLOH<-nrow(subset(CAP23, Inheritance=="notinherited" & GoHorLoH=="LoH"))
CAP23notinhGOH<-CAP23notinh-CAP23notinhLOH

CAP26inh<-nrow(subset(CAP26, Inheritance=="inherited"))
CAP26inhLOH<-nrow(subset(CAP26, Inheritance=="inherited" & GoHorLoH=="LoH"))
CAP26inhGOH<-CAP26inh-CAP26inhLOH
CAP26notinh<-nrow(CAP26)-CAP26inh
CAP26notinhLOH<-nrow(subset(CAP26, Inheritance=="notinherited" & GoHorLoH=="LoH"))
CAP26notinhGOH<-CAP26notinh-CAP26notinhLOH

CAP9inh<-nrow(subset(CAP9, Inheritance=="inherited"))
CAP9inhLOH<-nrow(subset(CAP9, Inheritance=="inherited" & GoHorLoH=="LoH"))
CAP9inhGOH<-CAP9inh-CAP9inhLOH
CAP9notinh<-nrow(CAP9)-CAP9inh
CAP9notinhLOH<-nrow(subset(CAP9, Inheritance=="notinherited" & GoHorLoH=="LoH"))
CAP9notinhGOH<-CAP9notinh-CAP9notinhLOH

CAP11inh<-nrow(subset(CAP11, Inheritance=="inherited"))
CAP11inhLOH<-nrow(subset(CAP11, Inheritance=="inherited" & GoHorLoH=="LoH"))
CAP11inhGOH<-CAP11inh-CAP11inhLOH
CAP11notinh<-nrow(CAP11)-CAP11inh
CAP11notinhLOH<-nrow(subset(CAP11, Inheritance=="notinherited" & GoHorLoH=="LoH"))
CAP11notinhGOH<-CAP11notinh-CAP11notinhLOH

allGOH<-sum(CAP6inhGOH, CAP6notinhGOH,CAP8inhGOH, CAP8notinhGOH,CAP22inhGOH, CAP22notinhGOH,CAP23inhGOH, CAP23notinhGOH,CAP26inhGOH, CAP26notinhGOH,CAP9inhGOH, CAP9notinhGOH,CAP11inhGOH, CAP11notinhGOH)
allLOH<-sum(CAP6inhLOH, CAP6notinhLOH,CAP8inhLOH, CAP8notinhLOH,CAP22inhLOH, CAP22notinhLOH,CAP23inhLOH, CAP23notinhLOH,CAP26inhLOH, CAP26notinhLOH,CAP9inhLOH, CAP9notinhLOH,CAP11inhLOH, CAP11notinhLOH)  
inhGOH<-sum(CAP6inhGOH, CAP8inhGOH, CAP22inhGOH, CAP23inhGOH,CAP26inhGOH, CAP9inhGOH,CAP11inhGOH)
inhGOHcounts<-c(CAP6inhGOH, CAP8inhGOH, CAP22inhGOH, CAP23inhGOH,CAP26inhGOH, CAP9inhGOH,CAP11inhGOH)
inhLOH<-sum(CAP6inhLOH, CAP8inhLOH, CAP22inhLOH, CAP23inhLOH,CAP26inhLOH, CAP9inhLOH,CAP11inhLOH)
inhLOHcounts<-c(CAP6inhLOH, CAP8inhLOH, CAP22inhLOH, CAP23inhLOH,CAP26inhLOH, CAP9inhLOH,CAP11inhLOH)

notinhGOH<-sum(CAP6notinhGOH, CAP8notinhGOH, CAP22notinhGOH, CAP23notinhGOH,CAP26notinhGOH, CAP9notinhGOH,CAP11notinhGOH)
notinhGOHcounts<-c(CAP6notinhGOH, CAP8notinhGOH, CAP22notinhGOH, CAP23notinhGOH,CAP26notinhGOH, CAP9notinhGOH,CAP11notinhGOH)

notinhLOH<-sum(CAP6notinhLOH, CAP8notinhLOH, CAP22notinhLOH, CAP23notinhLOH,CAP26notinhLOH, CAP9notinhLOH,CAP11notinhLOH)
notinhLOHcounts<-c(CAP6notinhLOH, CAP8notinhLOH, CAP22notinhLOH, CAP23notinhLOH,CAP26notinhLOH, CAP9notinhLOH,CAP11notinhLOH)

uglmGOH<-sum(nrow(subset(uniqueuglmCA56df, GoH_or_LoH=="DeNovo")), nrow(subset(uniqueuglmCA60df, GoH_or_LoH=="DeNovo")), nrow(subset(uniqueuglmCA65df, GoH_or_LoH=="DeNovo")))
uglmGOHcounts<-c(nrow(subset(uniqueuglmCA56df, GoH_or_LoH=="DeNovo")), nrow(subset(uniqueuglmCA60df, GoH_or_LoH=="DeNovo")), nrow(subset(uniqueuglmCA65df, GoH_or_LoH=="DeNovo")))

uglmLOH<-sum(nrow(subset(uniqueuglmCA56df, GoH_or_LoH=="LoH")), nrow(subset(uniqueuglmCA60df, GoH_or_LoH=="LoH")), nrow(subset(uniqueuglmCA65df, GoH_or_LoH=="LoH")))
gglmGOH<-sum(nrow(subset(gglmCA56df, GoH_or_LoH=="DeNovo")), nrow(subset(gglmCA60df, GoH_or_LoH=="DeNovo")), nrow(subset(gglmCA65df, GoH_or_LoH=="DeNovo")))
gglmGOHcounts<-c(nrow(subset(gglmCA56df, GoH_or_LoH=="DeNovo")), nrow(subset(gglmCA60df, GoH_or_LoH=="DeNovo")), nrow(subset(gglmCA65df, GoH_or_LoH=="DeNovo")))

gglmLOH<-sum(nrow(subset(gglmCA56df, GoH_or_LoH=="LoH")), nrow(subset(gglmCA60df, GoH_or_LoH=="LoH")), nrow(subset(gglmCA65df, GoH_or_LoH=="LoH")))
gglmLOHcounts<-c(nrow(subset(gglmCA56df, GoH_or_LoH=="LoH")), nrow(subset(gglmCA60df, GoH_or_LoH=="LoH")), nrow(subset(gglmCA65df, GoH_or_LoH=="LoH")))

subset<-c(rep("all",2),rep("inherited",2),rep("not inherited",2))
typez<-rep(c("LoH","GoH"),3)
numbers<-c(allLOH, allGOH, inhLOH,inhGOH,notinhLOH,notinhGOH)
percents<-c(allLOH/(allLOH+allGOH), allGOH/(allLOH+allGOH), inhLOH/(inhLOH+inhGOH),inhGOH/(inhLOH+inhGOH),notinhLOH/(notinhLOH+notinhGOH),notinhGOH/(notinhLOH+notinhGOH))
#sepercents<-c(se(percents)
set<-data.frame(subset,typez,numbers)
z<- ggplot(set, aes(x=subset, y=numbers,fill=typez)) +
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  scale_fill_brewer(palette = "Paired") +
  ylab("Number of mutations") + xlab("")+
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25))   +
  theme(legend.text=element_text(size=25))

GOHallsom<- c((CAP6inhGOH+ CAP6notinhGOH)/nrow(CAP6), (CAP8inhGOH+ CAP8notinhGOH)/nrow(CAP8), (CAP22inhGOH+ CAP22notinhGOH)/nrow(CAP22), (CAP23inhGOH+ CAP23notinhGOH)/nrow(CAP23), (CAP26inhGOH+ CAP26notinhGOH)/nrow(CAP26), (CAP9inhGOH+ CAP9notinhGOH)/nrow(CAP9), (CAP11inhGOH+ CAP11notinhGOH)/nrow(CAP11))
meanGOHallsom<- mean(GOHallsom)
seGOHallsom<-se( GOHallsom)
LOHallsom<-c((CAP6inhLOH+ CAP6notinhLOH)/nrow(CAP6), (CAP8inhLOH+ CAP8notinhLOH)/nrow(CAP8), (CAP22inhLOH+ CAP22notinhLOH)/nrow(CAP22), (CAP23inhLOH+ CAP23notinhLOH)/nrow(CAP23), (CAP26inhLOH+ CAP26notinhLOH)/nrow(CAP26), (CAP9inhLOH+ CAP9notinhLOH)/nrow(CAP9), (CAP11inhLOH+ CAP11notinhLOH)/nrow(CAP11))
meanLOHallsom<-mean(LOHallsom)
seLOHallsom<-se( LOHallsom)
GOHinh<-c				  (CAP6inhGOH/CAP6inh, (CAP8inhGOH)/CAP8inh, (CAP22inhGOH)/CAP22inh, (CAP23inhGOH)/CAP23inh, 
                (CAP26inhGOH)/CAP26inh, (CAP9inhGOH)/CAP9inh, (CAP11inhGOH)/CAP11inh)
meanGOHinh<-mean(GOHinh)
seGOHinh<-se(GOHinh)
LOHinh<- c(CAP6inhLOH/CAP6inh, (CAP8inhLOH)/CAP8inh, (CAP22inhLOH)/CAP22inh, (CAP23inhLOH)/CAP23inh, 
           (CAP26inhLOH)/CAP26inh, (CAP9inhLOH)/CAP9inh, (CAP11inhLOH)/CAP11inh)
meanLOHinh<-mean(LOHinh)
GOHnotinh<-c				  (CAP6notinhGOH/CAP6notinh, (CAP8notinhGOH)/CAP8notinh, (CAP22notinhGOH)/CAP22notinh, (CAP23notinhGOH)/CAP23notinh, 
                   (CAP26notinhGOH)/CAP26notinh, (CAP9notinhGOH)/CAP9notinh, (CAP11notinhGOH)/CAP11notinh)
LOHnotinh<-c				  (CAP6notinhLOH/CAP6notinh, (CAP8notinhLOH)/CAP8notinh, (CAP22notinhLOH)/CAP22notinh, (CAP23notinhLOH)/CAP23notinh, 
                   (CAP26notinhLOH)/CAP26notinh, (CAP9notinhLOH)/CAP9notinh, (CAP11notinhLOH)/CAP11notinh)
GOHuglm<-c(nrow(subset(uniqueuglmCA56df, startsWith(uniqueuglmCA56df$sample,'CAP6')==TRUE & GoH_or_LoH=="DeNovo")), nrow(subset(uniqueuglmCA56df, startsWith(uniqueuglmCA56df$sample,'CAP8')==TRUE & GoH_or_LoH=="DeNovo")),
           nrow(subset(uniqueuglmCA60df, startsWith(uniqueuglmCA60df$sample,'CAP22')==TRUE & GoH_or_LoH=="DeNovo")), nrow(subset(uniqueuglmCA60df, startsWith(uniqueuglmCA60df$sample,'CAP23')==TRUE & GoH_or_LoH=="DeNovo")), nrow(subset(uniqueuglmCA60df, startsWith(uniqueuglmCA60df$sample,'CAP26')==TRUE & GoH_or_LoH=="DeNovo")),
           nrow(subset(uniqueuglmCA65df, startsWith(uniqueuglmCA65df$sample,'CAP9')==TRUE & GoH_or_LoH=="DeNovo")), nrow(subset(uniqueuglmCA65df, startsWith(uniqueuglmCA65df$sample,'CAP11')==TRUE & GoH_or_LoH=="DeNovo")))
LOHuglm<-c(nrow(subset(uniqueuglmCA56df, startsWith(uniqueuglmCA56df$sample,'CAS6')==TRUE & GoH_or_LoH=="LoH")), nrow(subset(uniqueuglmCA56df, startsWith(uniqueuglmCA56df$sample,'CAS8')==TRUE & GoH_or_LoH=="LoH")),
           nrow(subset(uniqueuglmCA60df, startsWith(uniqueuglmCA60df$sample,'CAS22')==TRUE & GoH_or_LoH=="LoH")), nrow(subset(uniqueuglmCA60df, startsWith(uniqueuglmCA60df$sample,'CAS23')==TRUE & GoH_or_LoH=="LoH")), nrow(subset(uniqueuglmCA60df, startsWith(uniqueuglmCA60df$sample,'CAS26')==TRUE & GoH_or_LoH=="LoH")),
           nrow(subset(uniqueuglmCA65df, startsWith(uniqueuglmCA65df$sample,'CAS9')==TRUE & GoH_or_LoH=="LoH")), nrow(subset(uniqueuglmCA65df, startsWith(uniqueuglmCA65df$sample,'CAS11')==TRUE & GoH_or_LoH=="LoH")))
uglmsums<-GOHuglm + LOHuglm
GOHuglmprop<- GOHuglm/ uglmsums
LOHuglmprop<- LOHuglm/ uglmsums
meanGOHuglmprop<- mean(GOHuglmprop)
meanLOHuglmprop<- mean(LOHuglmprop)

GOHgglm<-c( nrow(subset(gglmCA56df, GoH_or_LoH=="DeNovo")), nrow(subset(gglmCA60df, GoH_or_LoH=="DeNovo")), nrow(subset(gglmCA65df, GoH_or_LoH=="DeNovo")))
LOHgglm<-c( nrow(subset(gglmCA56df, GoH_or_LoH=="LoH")), nrow(subset(gglmCA60df, GoH_or_LoH=="LoH")), nrow(subset(gglmCA65df, GoH_or_LoH=="LoH")))

gglmsums<- GOHgglm + LOHgglm

GOHgglmprop<- GOHgglm/ gglmsums
GOHgglmprop<-na.omit(GOHgglmprop)
LOHgglmprop<- LOHgglm/ gglmsums
LOHgglmprop<-na.omit(LOHgglmprop)

meanGOHgglmprop<- mean(GOHgglmprop)
meanLOHgglmprop<- mean(LOHgglmprop)

subsetz<-factor(c(rep("all somatic",2), rep("inherited",2),rep("not inherited",2), rep("uglm",2), rep("gglm",2)),level=c("all somatic","inherited","not inherited","uglm","gglm"))
means<-c(meanGOHallsom, meanLOHallsom, meanGOHinh, meanLOHinh, mean(GOHnotinh), mean(LOHnotinh), meanGOHuglmprop, meanLOHuglmprop, meanGOHgglmprop, meanLOHgglmprop)
allse<-c(seGOHallsom, seLOHallsom, seGOHinh, se(LOHinh),se(GOHnotinh), se(LOHnotinh), se(GOHuglmprop), se(LOHuglmprop), se(GOHgglmprop), se(LOHgglmprop))
type<-rep(c("GoH","LoH"),5)
#human sperm from Wang et al. 2012 Cell:
GOHcounts<-c(36,
             30,
             35,
             33,
             27,
             34,
             25,
             31)
LOHcounts<-c(10,
             12,
             10,
             13,
             6,
             10,
             16,
             13)
GOHprop<-GOHcounts/(GOHcounts+LOHcounts)
LOHprop<-1-GOHprop
meanGOHprop<-mean(GOHprop)
meanLOHprop<-mean(LOHprop)
seGOH<-se(GOHprop)
seLOH<-se(LOHprop)

df<-data.frame(subsetz, type, means, allse )
df$meanspercent<-100*df$means
df$allsepercent<-100*allse
df$shared<-factor(c(rep("All Parent",2),rep("P+S",2),rep("PO",2),rep("SSPO",2),rep("ASP",2)),
                  levels=c("All Parent","PO","P+S","SSPO","ASP"))
df$newcol<-factor(paste(df$shared, df$typez2))
df_withoutallsom<-df[-(1:2),]
df_eachdatapoint<-data.frame("percents"=c(GOHallsom, LOHallsom, GOHinh, LOHinh, GOHnotinh, LOHnotinh, GOHuglmprop, LOHuglmprop, GOHgglmprop, LOHgglmprop),
                             "shared"=factor(c(rep("All Parent",14),rep("P+S",14),rep("PO",14),rep("SSPO",14),rep("ASP",6)),
                                             levels=c("All Parent","PO","P+S","SSPO","ASP")),
                             "type"=c(rep(c(rep("GoH",7),rep("LoH",7)),4),rep("GoH",3),rep("LoH",3)))
z2<-ggplot(df, aes(x=shared, y=meanspercent)) +
  scale_color_manual(values=c("red","blue"))+
  geom_pointrange(data=df, size=1.5, position=position_dodge(width=0.6), aes(color=type,ymin=meanspercent-(allsepercent), ymax=meanspercent+(allsepercent))) +
  
  geom_jitter(data=df_eachdatapoint,size=3,alpha=0.5,position = position_jitterdodge(
    jitter.width = 0.1,
    jitter.height = 0,
    dodge.width = 0.6,
    seed = NA
  ),
  aes(shared, percents*100,colour=type
      ))+

  theme_bw() +
  
  #scale_fill_brewer(palette = "Paired") +
  #scale_fill_manual(values = c("blue", "red"))+
  ylab("Percent of SNVs") + xlab("")+
  guides(fill=guide_legend(title="Type"))+
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25))   +
  theme(legend.text=element_text(size=25),legend.title=element_text(size=25),
        axis.text.x = element_text(angle = 0, colour=c("black","goldenrod","#999999","purple","seagreen3")))
        #legend.title=element_blank())+
  #geom_errorbar(position="dodge",aes(ymin=meanspercent-(allsepercent), ymax=meanspercent+(allsepercent)))

z2+theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  

##this is figure 2d

##chisq for fig 2d:
notinheritedcounts<- c(df$meanspercent[5:6])
inheritedcounts<-  c(df$meanspercent[3:4])
sspocounts<-c(df$meanspercent[7:8])
aspcounts<-c(df$meanspercent[9:10])
n_v_i<-data.frame(notinheritedcounts,inheritedcounts)
n_v_sspo<-data.frame(notinheritedcounts,sspocounts)
i_v_sspo<-data.frame(inheritedcounts,sspocounts)
i_v_asp<-data.frame(inheritedcounts,aspcounts)
sspo_v_asp<-data.frame(sspocounts,aspcounts)
n_v_asp<-data.frame(notinheritedcounts,aspcounts)
chisq.test(n_v_i)
chisq.test(n_v_sspo)
chisq.test(n_v_asp)

chisq.test(i_v_sspo)
chisq.test(i_v_asp)
chisq.test(sspo_v_asp)

df

##test differences in mean GOH by grouping:
t.test(GOHinh, GOHnotinh, alternative=c("two.sided"),paired=TRUE)
ggqqplot(GOHnotinh)
ggdensity(GOHinh)
shapiro.test(GOHinh)
#use wilcox test instead!!
wilcox.test(GOHinh, GOHnotinh, paired = TRUE, alternative = "two.sided")

t.test(GOHinh, GOHuglm, alternative=c("two.sided"),paired=TRUE)
t.test(GOHnotinh, GOHuglm, alternative=c("two.sided"),paired=TRUE)
t.test(GOHgglm, GOHuglm, alternative=c("two.sided"),paired=TRUE)
t.test(GOHnotinh, GOHallsom, alternative=c("two.sided"),paired=TRUE)
t.test(GOHinh, GOHallsom, alternative=c("two.sided"),paired=TRUE)

###coding###
allsom_mis_counts<-c(nrow(subset(CAP6, MutationType=="missense_variant"))/CA56denom_coding, nrow(subset(CAP8, MutationType=="missense_variant"))/CA56denom_coding, 
                     nrow(subset(CAP22, MutationType=="missense_variant"))/CA60denom_coding, nrow(subset(CAP23, MutationType=="missense_variant"))/CA60denom_coding, nrow(subset(CAP26, MutationType=="missense_variant"))/CA60denom_coding,
                     nrow(subset(CAP11, MutationType=="missense_variant"))/CA65denom_coding,nrow(subset(CAP9, MutationType=="missense_variant"))/CA65denom_coding)
allsom_syn_counts<-c(nrow(subset(CAP6, MutationType=="synonymous_variant"))/CA56denom_coding, nrow(subset(CAP8, MutationType=="synonymous_variant"))/CA56denom_coding, 
                     nrow(subset(CAP22, MutationType=="synonymous_variant"))/CA60denom_coding, nrow(subset(CAP23, MutationType=="synonymous_variant"))/CA60denom_coding, nrow(subset(CAP26, MutationType=="synonymous_variant"))/CA60denom_coding,
                     nrow(subset(CAP11, MutationType=="synonymous_variant"))/CA65denom_coding,nrow(subset(CAP9, MutationType=="synonymous_variant"))/CA65denom_coding)

n_mis_counts<-c(nrow(subset(CAP6n, MutationType=="missense_variant"))/CA56denom_coding, nrow(subset(CAP8n, MutationType=="missense_variant"))/CA56denom_coding, 
                nrow(subset(CAP22n, MutationType=="missense_variant"))/CA60denom_coding, nrow(subset(CAP23n, MutationType=="missense_variant"))/CA60denom_coding, nrow(subset(CAP26n, MutationType=="missense_variant"))/CA60denom_coding,
                nrow(subset(CAP11n, MutationType=="missense_variant"))/CA65denom_coding,nrow(subset(CAP9n, MutationType=="missense_variant"))/CA65denom_coding)
n_syn_counts<-c(nrow(subset(CAP6n, MutationType=="synonymous_variant"))/CA56denom_coding, nrow(subset(CAP8n, MutationType=="synonymous_variant"))/CA56denom_coding, 
                nrow(subset(CAP22n, MutationType=="synonymous_variant"))/CA60denom_coding, nrow(subset(CAP23n, MutationType=="synonymous_variant"))/CA60denom_coding, nrow(subset(CAP26n, MutationType=="synonymous_variant"))/CA60denom_coding,
                nrow(subset(CAP11n, MutationType=="synonymous_variant"))/CA65denom_coding,nrow(subset(CAP9n, MutationType=="synonymous_variant"))/CA65denom_coding)
i_mis_counts<-c(nrow(subset(CAP6i, MutationType=="missense_variant"))/CA56denom_coding, nrow(subset(CAP8i, MutationType=="missense_variant"))/CA56denom_coding, 
                nrow(subset(CAP22i, MutationType=="missense_variant"))/CA60denom_coding, nrow(subset(CAP23i, MutationType=="missense_variant"))/CA60denom_coding, nrow(subset(CAP26i, MutationType=="missense_variant"))/CA60denom_coding,
                nrow(subset(CAP11i, MutationType=="missense_variant"))/CA65denom_coding,nrow(subset(CAP9i, MutationType=="missense_variant"))/CA65denom_coding)
i_syn_counts<-c(nrow(subset(CAP6i, MutationType=="synonymous_variant"))/CA56denom_coding, nrow(subset(CAP8i, MutationType=="synonymous_variant"))/CA56denom_coding, 
                nrow(subset(CAP22i, MutationType=="synonymous_variant"))/CA60denom_coding, nrow(subset(CAP23i, MutationType=="synonymous_variant"))/CA60denom_coding, nrow(subset(CAP26i, MutationType=="synonymous_variant"))/CA60denom_coding,
                nrow(subset(CAP11i, MutationType=="synonymous_variant"))/CA65denom_coding,nrow(subset(CAP9i, MutationType=="synonymous_variant"))/CA65denom_coding)
uglm_mis_counts<-c(nrow(subset(CAS6uglm, MutationType=="missense_variant"))/CA56denom_coding, nrow(subset(CAS8uglm, MutationType=="missense_variant"))/CA56denom_coding, 
                   nrow(subset(CAS22uglm, MutationType=="missense_variant"))/CA60denom_coding, nrow(subset(CAS23uglm, MutationType=="missense_variant"))/CA60denom_coding, nrow(subset(CAS26uglm, MutationType=="missense_variant"))/CA60denom_coding,
                   nrow(subset(CAS11uglm, MutationType=="missense_variant"))/CA65denom_coding,nrow(subset(CAS9uglm, MutationType=="missense_variant"))/CA65denom_coding)
uglm_syn_counts<-c(nrow(subset(CAS6uglm, MutationType=="synonymous_variant"))/CA56denom_coding, nrow(subset(CAS8uglm, MutationType=="synonymous_variant"))/CA56denom_coding, 
                   nrow(subset(CAS22uglm, MutationType=="synonymous_variant"))/CA60denom_coding, nrow(subset(CAS23uglm, MutationType=="synonymous_variant"))/CA60denom_coding, nrow(subset(CAS26uglm, MutationType=="synonymous_variant"))/CA60denom_coding,
                   nrow(subset(CAS11uglm, MutationType=="synonymous_variant"))/CA65denom_coding,nrow(subset(CAS9uglm, MutationType=="synonymous_variant"))/CA65denom_coding)
allsomratio<-allsom_mis_counts/(allsom_syn_counts+allsom_mis_counts)
nratio<-n_mis_counts/(n_syn_counts+n_mis_counts)
iratio<-i_mis_counts/(i_syn_counts+i_mis_counts)
uglmratio<-uglm_mis_counts/(uglm_syn_counts+uglm_mis_counts)
ratiomeans<-c(mean(allsomratio),mean(nratio), mean(iratio), mean(uglmratio))
ratioses<-c(se(allsomratio), se(nratio), se(iratio), se(uglmratio))
ratios_forgeompoint<-data.frame(ratiomeans, ratioses, types=c("All Somatic","PO","P+S","SSPO"))
ratiopoint1<-ggplot(ratios_forgeompoint, aes(x = types, y = ratiomeans*100, color=types))+
  geom_point(size=3)+
  scale_color_manual(values=c("black","red","blue","green"))+
  geom_errorbar(aes(ymin=ratiomeans*100-ratioses*100, ymax=ratiomeans*100+ratioses*100))+
  ylab("% of coding mutations that are missense")+xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(), axis.title = element_text(size=25),axis.text = element_text(size=25),legend.text = element_text(size=25))
ratios<-c(allsomratio,nratio, iratio, uglmratio)
ratiotypes<-c(rep("All Somatic",7), rep("PO",7), rep("P+S",7),rep("SSPO",7))
ratiosamples<-rep(c("CAP6","CAP8","CAP22","CAP23","CAP26","CAP11","CAP9"),4)
ratio_df<-data.frame(ratios,ratiotypes,ratiosamples)
ratiopoint2<-ggplot(ratio_df, aes(x = ratiosamples, y = ratios*100, color=ratiotypes))+
  geom_jitter(height=0, width=0.15,size=3)+
  scale_color_manual(values=c("black","red","blue","green"))+
  geom_hline(yintercept=ratiomeans[1]*100,color="black",linetype="dashed")+geom_hline(yintercept=ratiomeans[2]*100,color="blue",linetype="dashed")+geom_hline(yintercept=ratiomeans[3]*100,color="red",linetype="dashed")+geom_hline(yintercept=ratiomeans[4]*100,color="green",linetype="dashed")+
  ylab("% of coding mutations that are missense")+xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(), axis.title = element_text(size=25),axis.text = element_text(size=25),legend.text = element_text(size=25))
#use this 20210720:
ratiopoint<-ggplot(ratios_forgeompoint, aes(x = types, y = ratiomeans*100))+
  scale_color_manual(values=c("black","red","blue","green"))+
  
  geom_pointrange(size=1.5,position=position_dodge(width=0.6), aes(color=types, ymin=100*(ratiomeans-ratioses), ymax=100*(ratiomeans+ratioses))) +
  geom_jitter(data=ratio_df,alpha=0.5,size=3,aes(ratiotypes, 100*ratios, color=ratiotypes),
              position = position_jitter(0.2)
  )+
  ylab("% of coding mutations that are missense")+xlab("")+
  theme_bw()+
  theme(legend.title=element_blank(), axis.title = element_text(size=25),axis.text = element_text(size=25),legend.text = element_text(size=25))
counts<-c(allsom_mis_counts,n_mis_counts, i_mis_counts, uglm_mis_counts,
          allsom_syn_counts,n_syn_counts, i_syn_counts, uglm_syn_counts)
mean(allsom_mis_counts+allsom_syn_counts)
mean(n_mis_counts+n_syn_counts)
mean(i_mis_counts+i_syn_counts)
mean(uglm_mis_counts+uglm_syn_counts)

muttypes<-c(rep("missense",28),rep("synonymous",28))
counts_df<-data.frame(counts,muttypes,types=rep(ratiotypes,2),samples=rep(ratiosamples,2))
counts_df$second<-factor(paste0(counts_df$types , counts_df$muttypes))
countsplot<-ggplot(counts_df, aes(x = samples, y = counts, color=types, shape=muttypes, group=types))+
  geom_point(size=3, position=position_dodge(width=0.6))+ylab("# of mutations")+
  scale_color_manual(values=c("black","red","blue","green"))+
  ylab("# of mutations per bp in coding region")+ xlab("")+
  theme_bw()+
  geom_hline(yintercept=mean(allsom_mis_counts+allsom_syn_counts),color="black",linetype="dashed")+geom_hline(yintercept=mean(n_mis_counts+n_syn_counts),color="blue",linetype="dashed")+geom_hline(yintercept=mean(i_mis_counts+i_syn_counts),color="red",linetype="dashed")+geom_hline(yintercept=mean(uglm_mis_counts+uglm_syn_counts),color="green",linetype="dashed")+
  theme(legend.title=element_blank(), axis.title = element_text(size=25),axis.text = element_text(size=25),legend.text = element_text(size=25))

 #suppfig6 20210720: 
(countsplot+ labs(tag = "a.")+theme(plot.tag=element_text(face="bold",size=30))) | (ratiopoint+ labs(tag = "b.")+theme(plot.tag=element_text(face="bold",size=30))) 
wilcox.test(iratio,uglmratio)
wilcox.test(iratio, uglmratio, paired = TRUE, alternative = "two.sided")
wilcox.test(nratio, uglmratio, paired = TRUE, alternative = "two.sided")
wilcox.test(allsomratio, uglmratio, paired = TRUE, alternative = "two.sided")

t.test(nratio,uglmratio)
t.test(allsomratio,uglmratio)

box<-ggplot(ratio_df, aes(x = ratiotypes, y = ratios, color=ratiotypes))+
  geom_boxplot()+
  scale_color_manual(values=c("red","blue","green"))+
  theme_bw()
box2 <- box + labs(x="", y="Percent of Somatic Mutations inherited by Sperm") + geom_boxplot(fill="red",
                                                                                             position = position_dodge(0.9)
) +
  
  theme_bw()+
  #stat_compare_means(label.x=1,label.y=30,size=10)+
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25))
Ahya_fasta$chr1[1:10]

