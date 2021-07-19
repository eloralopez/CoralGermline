library("seqinr")
library(patchwork)
Amil_fasta<-read.fasta("~/Documents/GitHub/CoralGermline/dndscv/Amil.v2.01.chrs.fasta",forceDNAtolower = FALSE)
files<-list.files(path="~/Documents/GitHub/CoralGermline/annotatedfiles", pattern="*ann.txt", full.names=T, recursive=FALSE) #path to all the files you want to include in the analysis

CA56_samples<-list.files(path="~/Documents/GitHub/CoralGermline/annotatedfiles", pattern="*relabeled56.txt.ann.txt", full.names=T, recursive=FALSE)
files<-list.files(path="~/Documents/GitHub/CoralGermline/annotatedfiles", pattern="*relabeled60.txt.ann.txt", full.names=T, recursive=FALSE)
CA<-list.files(path="~/Documents/GitHub/CoralGermline/annotatedfiles", pattern="*relabeled65.txt.ann.txt", full.names=T, recursive=FALSE)

# bigfunction<- function(files) {
# #files<-list.files(path="~/Documents/GitHub/CoralGermline/annotatedfiles", pattern="*ann.txt", full.names=T, recursive=FALSE) #path to all the files you want to include in the analysis
#   metadata= NULL
#   for (i in 1:length(files)) { 
#     file =files[i]
#     data<-read.delim(file) #read in each file in "files"
#     data<-data.frame(data) # transform the data from each file into a dataframe
#     base<-basename(file)
#     colony<-strsplit(base, "\\_")[[1]][2]
#     len<-nrow(data) 
#     colonyrep<-rep(colony, len)
#     withcolony<-data.frame(data, colonyrep) #combines the colonyname column with each dataframe
#     metadata <- rbind(metadata, withcolony) #adds each dataframe to the overall metatadata, so that the information from all of the files are now in "metadata"
#   }
#   
#   genoanddepth<-(metadata$genotype) #names the column
#   split<-str_split_fixed(genoanddepth, ",", 4) #split the genotype, depths, and GQ score into their own separate strings
#   
#   genotypes<-split[,1] #defines the genotype as the first string in "split"
#   position<-metadata$chrom.pos #names the column
#   positionsplit<-str_split_fixed(position, "[.]", 2) #split the chromosome number and the position on the chromosome into their own separate strings
#   
#   chr<-positionsplit[,1] #defines the chromosome as the first string in "positionsplit"
#   pos<-positionsplit[,2] #defines the position as the second string in "positionsplit"
#   #totaldepth<-as.numeric(split[,2])
#   refdepth<-as.numeric(split[,2])
#   altdepth<-as.numeric(split[,3])
#   what<-metadata$WhattoWhat
#   allelesplit<-str_split_fixed(what, "to", 2)
#   normalallele<-allelesplit[,1]
#   mutantallele<-allelesplit[,2]
#   totaldepth<-refdepth+altdepth
#   GQscore<-as.numeric(split[,4])
#   mutationtype<-metadata$MutationType
#   mutationstrength<-metadata$MutationStrength
#   
#   mutant_alleledepth = rep("A", nrow(metadata))
#   for (i in 1:nrow(metadata)){
#     if (mutantallele[i] == metadata$alt[i]) {
#       mutant_alleledepth[i] = altdepth[i]
#       #print(mutant_alleledepth[i])
#     } else {
#       mutant_alleledepth[i] = refdepth[i]
#       #print(mutant_alleledepth[i])
#     }  #print("ALT", mutant_alleledepth, refdepth)
#   }
#   #print(mutant_alleledepth[1:10])
#   #print(refdepth[1:10])
#   #print(altdepth[1:10])
#   
#   metadatadf<-data.frame("chrom.pos" = metadata$chrom.pos, "chrom"=chr, "pos"=pos,	"sample"= metadata$sample, "ref" = metadata$ref, "alt" = 	metadata$alt, "normal_allele"= normalallele, "mutant_allele" = mutantallele, "mutant_allele_depth" = as.numeric(mutant_alleledepth), "genotype"= genotypes, "totaldepth"=totaldepth, 	"refdepth"=refdepth, "altdepth"=altdepth, "GQscore"= GQscore,	"GoH_or_LoH"=metadata$DeNovo_LoH, "Ti/Tv"=metadata$TiTv, 	"WhattoWhat" = metadata$WhattoWhat,"TrueorFalse" =metadata$TrueorFalse, "MutationType"=mutationtype, "MutationStrength"=mutationstrength)#  ColonyName"=metadata$colonyrep)
#   
#   DepthMeansdf<-aggregate(totaldepth~chrom.pos, metadatadf, 			FUN=mean)
#   
#   DepthMinsdf<-aggregate(totaldepth~chrom.pos, metadatadf, 			FUN=min)
#   
#   GQaverage<-aggregate(GQscore~chrom.pos, metadatadf, FUN=mean)
#   
#   GQmin<-aggregate(GQscore~chrom.pos, metadatadf, FUN=min)
#   
#   
#   
#   metadatadf.00<-merge(metadatadf, DepthMeansdf[, c("chrom.pos", 	"totaldepth")], by="chrom.pos")
#   metadatadf.0<-merge(metadatadf.00, GQaverage[,c("chrom.pos","GQscore")], by="chrom.pos")
#   metadatadf.0<-merge(metadatadf.0, DepthMinsdf[,c("chrom.pos","totaldepth")], by="chrom.pos")
#   metadatadf.0<-merge(metadatadf.0, GQmin[,c("chrom.pos","GQscore")], by="chrom.pos")
#   
#   DeNovos<-subset(metadatadf.0, GoH_or_LoH=="DeNovo")
#   # sample3<-subset(DeNovos, sample== "sample3")
#   # trueDenovos_sample3<-subset(sample3, refdepth =="0" | altdepth=="0")
#   # 
#   # sample4<-subset(DeNovos, sample== "sample4")
#   # trueDenovos_sample4<-subset(sample4, refdepth =="0" | altdepth=="0")
#   # 
#   # sample5<-subset(DeNovos, sample== "sample5")
#   # trueDenovos_sample5<-subset(sample5, refdepth =="0" | altdepth=="0")
#   # 
#   # sample6<-subset(DeNovos, sample== "sample6")
#   # trueDenovos_sample6<-subset(sample6, refdepth =="0" | altdepth=="0")
#   # 
#   # sample7<-subset(DeNovos, sample== "sample7")
#   # trueDenovos_sample7<-subset(sample7, refdepth =="0" | altdepth=="0")
#   # 
#   # sample8<-subset(DeNovos, sample== "sample8")
#   # trueDenovos_sample8<-subset(sample8, refdepth =="0" | altdepth=="0")
#   # 
#   # truedenovos3_8<-rbind(trueDenovos_sample3, trueDenovos_sample4, trueDenovos_sample5, trueDenovos_sample6, trueDenovos_sample7, trueDenovos_sample8)
#   
#   LoH<-subset(metadatadf.0, GoH_or_LoH =="LoH")#
#   trueLoHp<-subset(LoH,refdepth =="0" | altdepth=="0")
#   trueLoHp1<-subset(trueLoHp, sample=="mutparent1" | sample=="mutparent2")
#   trueLoHp2<-trueLoHp1[trueLoHp1$chrom.pos %in% names(which(table(trueLoHp1$chrom.pos) > 1)), ]
#   sperm<-subset(metadatadf.0, sample =="mutsperm")
#   trueLoHsperm<-sperm[match(trueLoHp2$chrom.pos, sperm$chrom.pos),]
#   trueLoHsperm<-unique(trueLoHsperm)
#   
#   #trueLoHp2<-subset(trueLoHp, sample=="mutparent2")
#   metadatadf<-rbind( DeNovos, trueLoHp2,trueLoHsperm)#, trueLoHp2)
#   ##write.table(metadatadf, file="CAcolony60_CAP22-23-24muts_20191125.txt",sep="\t",quote=FALSE, row.name=FALSE)
#   
#   uniquemetadatadf<- metadatadf[match(unique(metadatadf$chrom.pos), 					metadatadf$chrom.pos),]
#   
#   uniquemetadatadfcount<-nrow(uniquemetadatadf)
#   
#   missense<-subset(uniquemetadatadf, MutationType=="missense_variant")
#   missensecount<-nrow(missense)
#   
#   syn<-subset(uniquemetadatadf, MutationType=="synonymous_variant")
#   syncount<-nrow(syn)
#   #par(mfrow=c(1,3))
#   
#   dataframe<-uniquemetadatadf
#   dataframecount<-uniquemetadatadfcount
#   
dataframe<-CA56_
dataframecount<-nrow(CA56_)
spectrum<- function(dataframe, dataframecount) {
  #Amil_fasta<-read.fasta("~/Documents/GitHub/CoralGermline/dndscv/Amil.v2.01.chrs.fasta",forceDNAtolower = FALSE)
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
      
      chroms<-Amil_fasta[[chr_num]]
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
      
      chroms<-Amil_fasta[[chr_num]]
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
    geom_bar(stat="identity") +
    ylim(0, 0.35)+
    theme_bw()+
    ylab("Proportion of SNVs") +
    theme(axis.text.x = element_text(angle = 90, size=18), axis.text.y=element_text(size=15),axis.title.y=element_text(size=15))+#, axis.text.x=element_text(size=25))+
    geom_errorbar(aes(ymin=coralProportion-(se_forprops), ymax=coralProportion+(se_forprops)))#, width=.2,
  
  Transtions<- verATGCprop + verGCATprop
  Transversions<- verATTAprop+verACTGprop+verGCCGprop+verGCTAprop
  TiTv<- Transtions/Transversions
  counts<-c(verATGCcount,CpGcount, CpHcount, verATTAcount, verACTGcount, verGCCGcount, verGCTAcount)
  return(z)
  #return(list(vertypesDF, TiTv))
  #return(list(z, vertypesDF))
  #return(z)
  #return(counts)
}
#10/19/2020:
CA56_somatic_all<-spectrum(a_,nrow(a_))
CA60_somatic_all<-spectrum(b_,nrow(b_))
CA65_somatic_all<-spectrum(c_,nrow(c_))
allcol<-rbind(a_,b_,c_)
allcol_somatic_all<-spectrum(allcol,nrow(allcol)) + ggtitle("All Somatic Mutations")
allcol_somatic_all_dndscv<-data.frame("sampleID"= allcol$sample,"chr"= allcol$chrom,"pos"= allcol$pos, "ref" = allcol$ref, "alt"= allcol$alt)
##write.table(allcol_somatic_all_dndscv, file="~/Documents/GitHub/CoralGermline/dndscv/somaticallcolonies_20201019_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)

allcol_somatic_GOH<-spectrum(subset(allcol, GoHorLoH=="DeNovo"),nrow(subset(allcol, GoHorLoH=="DeNovo")))+ggtitle("GoH SNVs found in Parents") 
allcol_somatic_GOH_dndscv<-data.frame("sampleID"= subset(allcol, GoHorLoH=="DeNovo")$sample,"chr"= subset(allcol, GoHorLoH=="DeNovo")$chrom,"pos"= subset(allcol, GoHorLoH=="DeNovo")$pos, "ref" = subset(allcol, GoHorLoH=="DeNovo")$ref, "alt"= subset(allcol, GoHorLoH=="DeNovo")$alt)
##write.table(allcol_somatic_GOH_dndscv, file="~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesGOH_20201019_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)

allcol_somatic_LOH<-spectrum(subset(allcol, GoHorLoH=="LoH"),nrow(subset(allcol, GoHorLoH=="LoH")))+ggtitle("LoH SNVs found in Parents") 
allcol_somatic_LOH_dndscv<-data.frame("sampleID"= subset(allcol, GoHorLoH=="LoH")$sample,"chr"= subset(allcol, GoHorLoH=="LoH")$chrom,"pos"= subset(allcol, GoHorLoH=="LoH")$pos, "ref" = subset(allcol, GoHorLoH=="LoH")$ref, "alt"= subset(allcol, GoHorLoH=="LoH")$alt)
##write.table(allcol_somatic_LOH_dndscv, file="~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesLOH_20201019_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)  

CA56_somatic_inherited<-spectrum(subset(a_, Inheritance=="inherited"), nrow(subset(a_, Inheritance=="inherited")))
CA60_somatic_inherited<-spectrum(subset(b_, Inheritance=="inherited"), nrow(subset(b_, Inheritance=="inherited")))
CA65_somatic_inherited<-spectrum(subset(c_, Inheritance=="inherited"), nrow(subset(c_, Inheritance=="inherited")))
allcol_somatic_inherited<-spectrum(subset(allcol, Inheritance=="inherited"), nrow(subset(allcol, Inheritance=="inherited"))) +ggtitle("SNVs shared by Parent and Sperm (n=)") 
allcol_somatic_inherited_dndscv<-data.frame("sampleID"= subset(allcol, Inheritance=="inherited")$sample,"chr"= subset(allcol, Inheritance=="inherited")$chrom,"pos"= subset(allcol, Inheritance=="inherited")$pos, "ref" = subset(allcol, Inheritance=="inherited")$ref, "alt"= subset(allcol, Inheritance=="inherited")$alt)
##write.table(allcol_somatic_inherited_dndscv, file="~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesinherited_20201019_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)  
allcol_somatic_inherited_GOH_dndscv<-data.frame("sampleID"= subset(allcol, Inheritance=="inherited" & GoHorLoH=="DeNovo")$sample,"chr"= subset(allcol, Inheritance=="inherited" & GoHorLoH=="DeNovo")$chrom,"pos"= subset(allcol, Inheritance=="inherited" & GoHorLoH=="DeNovo")$pos, "ref" = subset(allcol, Inheritance=="inherited" & GoHorLoH=="DeNovo")$ref, "alt"= subset(allcol, Inheritance=="inherited" & GoHorLoH=="DeNovo")$alt)
allcol_somatic_inherited_LOH_dndscv<-data.frame("sampleID"= subset(allcol, Inheritance=="inherited" & GoHorLoH=="LoH")$sample,"chr"= subset(allcol, Inheritance=="inherited" & GoHorLoH=="LoH")$chrom,"pos"= subset(allcol, Inheritance=="inherited" & GoHorLoH=="LoH")$pos, "ref" = subset(allcol, Inheritance=="inherited" & GoHorLoH=="LoH")$ref, "alt"= subset(allcol, Inheritance=="inherited" & GoHorLoH=="LoH")$alt)
##write.table(allcol_somatic_inherited_GOH_dndscv, file="~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesinherited_GOH_20201019_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)  
##write.table(allcol_somatic_inherited_LOH_dndscv, file="~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesinherited_LOH_20201019_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)  

CA56_somatic_notinherited<-spectrum(subset(a_, Inheritance=="notinherited"), nrow(subset(a_, Inheritance=="notinherited")))
CA60_somatic_notinherited<-spectrum(subset(b_, Inheritance=="notinherited"), nrow(subset(b_, Inheritance=="notinherited")))
CA65_somatic_notinherited<-spectrum(subset(c_, Inheritance=="notinherited"), nrow(subset(c_, Inheritance=="notinherited")))
allcol_somatic_notinherited<-spectrum(subset(allcol, Inheritance=="notinherited"), nrow(subset(allcol, Inheritance=="notinherited")))+ ggtitle("SNVs in Parent Only") 
allcol_somatic_notinherited_dndscv<-data.frame("sampleID"= subset(allcol, Inheritance=="notinherited")$sample,"chr"= subset(allcol, Inheritance=="notinherited")$chrom,"pos"= subset(allcol, Inheritance=="notinherited")$pos, "ref" = subset(allcol, Inheritance=="notinherited")$ref, "alt"= subset(allcol, Inheritance=="notinherited")$alt)
##write.table(allcol_somatic_notinherited_dndscv, file="~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesnotinherited_20201019_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)  
allcol_somatic_notinherited_GOH_dndscv<-data.frame("sampleID"= subset(allcol, Inheritance=="notinherited" & GoHorLoH=="DeNovo")$sample,"chr"= subset(allcol, Inheritance=="notinherited" & GoHorLoH=="DeNovo")$chrom,"pos"= subset(allcol, Inheritance=="notinherited" & GoHorLoH=="DeNovo")$pos, "ref" = subset(allcol, Inheritance=="notinherited" & GoHorLoH=="DeNovo")$ref, "alt"= subset(allcol, Inheritance=="notinherited" & GoHorLoH=="DeNovo")$alt)
allcol_somatic_notinherited_LOH_dndscv<-data.frame("sampleID"= subset(allcol, Inheritance=="notinherited" & GoHorLoH=="LoH")$sample,"chr"= subset(allcol, Inheritance=="notinherited" & GoHorLoH=="LoH")$chrom,"pos"= subset(allcol, Inheritance=="notinherited" & GoHorLoH=="LoH")$pos, "ref" = subset(allcol, Inheritance=="notinherited" & GoHorLoH=="LoH")$ref, "alt"= subset(allcol, Inheritance=="notinherited" & GoHorLoH=="LoH")$alt)
##write.table(allcol_somatic_notinherited_GOH_dndscv, file="~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesnotinherited_GOH_20201019_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)  
##write.table(allcol_somatic_notinherited_LOH_dndscv, file="~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesnotinherited_LOH_20201019_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)  

CA56_uglm_spec<-spectrum(uniqueuglmCA56,nrow(uniqueuglmCA56))
CA60_uglm_spec<-spectrum(uniqueuglmCA60,nrow(uniqueuglmCA60))
CA65_uglm_spec<-spectrum(uniqueuglmCA65,nrow(uniqueuglmCA65))
allcol_uglm<-rbind(uniqueuglmCA56, uniqueuglmCA60, uniqueuglmCA65)
allcol_uglm_GOH<- subset(allcol_uglm, GoH_or_LoH=="DeNovo")
allcol_uglm_GOH_dndscv<-data.frame("sampleID"= allcol_uglm_GOH$sample,"chr"= allcol_uglm_GOH$chrom,"pos"= allcol_uglm_GOH$pos, "ref" = allcol_uglm_GOH$ref, "alt"= allcol_uglm_GOH$alt)
##write.table(allcol_uglm_GOH_dndscv,file="~/Documents/GitHub/CoralGermline/dndscv/uglmallcolonies_GOH_20201019_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)

allcol_uglm_LOH<- subset(allcol_uglm, GoH_or_LoH=="LoH")
allcol_uglm_LOH_dndscv<-data.frame("sampleID"= allcol_uglm_LOH$sample,"chr"= allcol_uglm_LOH$chrom,"pos"= allcol_uglm_LOH$pos, "ref" = allcol_uglm_LOH$ref, "alt"= allcol_uglm_LOH$alt)
##write.table(allcol_uglm_LOH_dndscv,file="~/Documents/GitHub/CoralGermline/dndscv/uglmallcolonies_LOH_20201019_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)


allcol_uglm_spec<-spectrum(allcol_uglm, nrow(allcol_uglm)) + ggtitle("SNVs in Single Sperm Pool Only") 
allcol_uglm_dndscv<-data.frame("sampleID"= allcol_uglm$sample,"chr"= allcol_uglm$chrom,"pos"= allcol_uglm$pos, "ref" = allcol_uglm$ref, "alt"= allcol_uglm$alt)
##write.table(allcol_uglm_dndscv,file="~/Documents/GitHub/CoralGermline/dndscv/uglmallcolonies_20201019_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)

allcol_gglm<-rbind(gglmCA56, gglmCA60, gglmCA65)
allcol_gglm_spec<-spectrum(allcol_gglm, nrow(allcol_gglm)) + ggtitle("SNVs in All Sperm Pools") 
allcol_gglm_dndscv<-data.frame("sampleID"= allcol_gglm$sample,"chr"= allcol_gglm$chrom,"pos"= allcol_gglm$pos, "ref" = allcol_gglm$ref, "alt"= allcol_gglm$alt)
##write.table(allcol_gglm_dndscv,file="~/Documents/GitHub/CoralGermline/dndscv/gglmallcolonies_20201019_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)

#this is figure 6:
(allcol_somatic_all | allcol_somatic_inherited | allcol_somatic_notinherited) / (allcol_uglm_spec | allcol_gglm_spec)
allcol_somatic_all | allcol_somatic_GOH | allcol_somatic_LOH

suppfig_Spec<-(allcol_somatic_inherited | allcol_somatic_notinherited) / (allcol_uglm_spec | allcol_gglm_spec)

##end figure 6##


#11/17/2020:
CA56_somatic_all<-spectrum(CA56_,nrow(CA56_))
CA60_somatic_all<-spectrum(CA60_,nrow(CA60_))
CA65_somatic_all<-spectrum(CA65_,nrow(CA65_))
allcol<-rbind(CA56_,CA60_,CA65_)
all_som_GOH_inherited<-subset(allcol,GoHorLoH=="DeNovo"&Inheritance=="inherited"&ParentAverage>0.5)
all_som_GOH_notinherited<-subset(allcol,GoHorLoH=="DeNovo"& Inheritance=="notinherited")

allcol_somatic_all<-spectrum(allcol,nrow(allcol)) + ggtitle("All Somatic Mutations")
allcol_somatic_all_dndscv<-data.frame("sampleID"= allcol$sample,"chr"= allcol$chrom,"pos"= allcol$pos, "ref" = allcol$ref, "alt"= allcol$alt)
#write.table(allcol_somatic_all_dndscv, file="~/Documents/GitHub/CoralGermline/dndscv/somaticallcolonies_20201210_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)

allcol_somatic_GOH<-spectrum(subset(allcol, GoHorLoH=="DeNovo"),nrow(subset(allcol, GoHorLoH=="DeNovo")))+ggtitle("GoH SNVs found in Parents") 
allcol_somatic_GOH_dndscv<-data.frame("sampleID"= subset(allcol, GoHorLoH=="DeNovo")$sample,"chr"= subset(allcol, GoHorLoH=="DeNovo")$chrom,"pos"= subset(allcol, GoHorLoH=="DeNovo")$pos, "ref" = subset(allcol, GoHorLoH=="DeNovo")$ref, "alt"= subset(allcol, GoHorLoH=="DeNovo")$alt)
#write.table(allcol_somatic_GOH_dndscv, file="~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesGOH_20201210_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)

allcol_somatic_LOH<-spectrum(subset(allcol, GoHorLoH=="LoH"),nrow(subset(allcol, GoHorLoH=="LoH")))+ggtitle("LoH SNVs found in Parents") 
allcol_somatic_LOH_dndscv<-data.frame("sampleID"= subset(allcol, GoHorLoH=="LoH")$sample,"chr"= subset(allcol, GoHorLoH=="LoH")$chrom,"pos"= subset(allcol, GoHorLoH=="LoH")$pos, "ref" = subset(allcol, GoHorLoH=="LoH")$ref, "alt"= subset(allcol, GoHorLoH=="LoH")$alt)
#write.table(allcol_somatic_LOH_dndscv, file="~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesLOH_20201210_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)  
all_GOH<-subset(allcol, GoHorLoH=="DeNovo")
all_LOH<-subset(allcol, GoHorLoH=="LoH")  
CA56_somatic_inherited<-spectrum(subset(CA56_, Inheritance=="inherited"), nrow(subset(CA56_, Inheritance=="inherited")))
CA56_somatic_inherited_dndscv<-data.frame("sampleID"= subset(CA56_, Inheritance=="inherited")$sample,"chr"= subset(CA56_, Inheritance=="inherited")$chrom,"pos"= subset(CA56_, Inheritance=="inherited")$pos, "ref" = subset(CA56_, Inheritance=="inherited")$ref, "alt"= subset(CA56_, Inheritance=="inherited")$alt)
#write.table(CA56_somatic_inherited_dndscv, file="~/Documents/GitHub/CoralGermline/dndscv/somaticCA56inherited_20201210_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)  

CA60_somatic_inherited<-spectrum(subset(CA60_, Inheritance=="inherited"), nrow(subset(CA60_, Inheritance=="inherited")))
CA60_somatic_inherited_dndscv<-data.frame("sampleID"= subset(CA60, Inheritance=="inherited")$sample,"chr"= subset(CA60, Inheritance=="inherited")$chrom,"pos"= subset(CA60, Inheritance=="inherited")$pos, "ref" = subset(CA60, Inheritance=="inherited")$ref, "alt"= subset(CA60, Inheritance=="inherited")$alt)
#write.table(CA60_somatic_inherited_dndscv, file="~/Documents/GitHub/CoralGermline/dndscv/somaticCA60inherited_20201203_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)  
CA65_somatic_inherited<-spectrum(subset(CA65_, Inheritance=="inherited"), nrow(subset(CA65_, Inheritance=="inherited")))
CA65_somatic_inherited_dndscv<-data.frame("sampleID"= subset(CA65, Inheritance=="inherited")$sample,"chr"= subset(CA65, Inheritance=="inherited")$chrom,"pos"= subset(CA65, Inheritance=="inherited")$pos, "ref" = subset(CA65, Inheritance=="inherited")$ref, "alt"= subset(CA65, Inheritance=="inherited")$alt)
#write.table(CA65_somatic_inherited_dndscv, file="~/Documents/GitHub/CoralGermline/dndscv/somaticCA65oniesinherited_20201203_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)  
allcol_somatic_inherited<-spectrum(subset(allcol, Inheritance=="inherited"), nrow(subset(allcol, Inheritance=="inherited"))) +ggtitle("Shared Parent and Sperm SNVs") 
allcol_somatic_inherited_dndscv<-data.frame("sampleID"= subset(allcol, Inheritance=="inherited")$sample,"chr"= subset(allcol, Inheritance=="inherited")$chrom,"pos"= subset(allcol, Inheritance=="inherited")$pos, "ref" = subset(allcol, Inheritance=="inherited")$ref, "alt"= subset(allcol, Inheritance=="inherited")$alt)
#write.table(allcol_somatic_inherited_dndscv, file="~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesinherited_20201210_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)  
allcol_somatic_inherited_GOH_dndscv<-data.frame("sampleID"= subset(allcol, Inheritance=="inherited" & GoHorLoH=="DeNovo")$sample,"chr"= subset(allcol, Inheritance=="inherited" & GoHorLoH=="DeNovo")$chrom,"pos"= subset(allcol, Inheritance=="inherited" & GoHorLoH=="DeNovo")$pos, "ref" = subset(allcol, Inheritance=="inherited" & GoHorLoH=="DeNovo")$ref, "alt"= subset(allcol, Inheritance=="inherited" & GoHorLoH=="DeNovo")$alt)
allcol_somatic_inherited_LOH_dndscv<-data.frame("sampleID"= subset(allcol, Inheritance=="inherited" & GoHorLoH=="LoH")$sample,"chr"= subset(allcol, Inheritance=="inherited" & GoHorLoH=="LoH")$chrom,"pos"= subset(allcol, Inheritance=="inherited" & GoHorLoH=="LoH")$pos, "ref" = subset(allcol, Inheritance=="inherited" & GoHorLoH=="LoH")$ref, "alt"= subset(allcol, Inheritance=="inherited" & GoHorLoH=="LoH")$alt)
#write.table(allcol_somatic_inherited_GOH_dndscv, file="~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesinherited_GOH_20201210_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)  
#write.table(allcol_somatic_inherited_LOH_dndscv, file="~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesinherited_LOH_20201210_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)  

CA56_somatic_notinherited<-spectrum(subset(CA56_, Inheritance=="notinherited"), nrow(subset(CA56_, Inheritance=="notinherited")))
CA60_somatic_notinherited<-spectrum(subset(CA60_, Inheritance=="notinherited"), nrow(subset(CA60_, Inheritance=="notinherited")))
CA65_somatic_notinherited<-spectrum(subset(CA65_, Inheritance=="notinherited"), nrow(subset(CA65_, Inheritance=="notinherited")))
allcol_somatic_notinherited<-spectrum(subset(allcol, Inheritance=="notinherited"), nrow(subset(allcol, Inheritance=="notinherited")))+ ggtitle("Parent Only SNVs") 
allcol_somatic_notinherited_dndscv<-data.frame("sampleID"= subset(allcol, Inheritance=="notinherited")$sample,"chr"= subset(allcol, Inheritance=="notinherited")$chrom,"pos"= subset(allcol, Inheritance=="notinherited")$pos, "ref" = subset(allcol, Inheritance=="notinherited")$ref, "alt"= subset(allcol, Inheritance=="notinherited")$alt)
#write.table(allcol_somatic_notinherited_dndscv, file="~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesnotinherited_20201210_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)  
allcol_somatic_notinherited_GOH_dndscv<-data.frame("sampleID"= subset(allcol, Inheritance=="notinherited" & GoHorLoH=="DeNovo")$sample,"chr"= subset(allcol, Inheritance=="notinherited" & GoHorLoH=="DeNovo")$chrom,"pos"= subset(allcol, Inheritance=="notinherited" & GoHorLoH=="DeNovo")$pos, "ref" = subset(allcol, Inheritance=="notinherited" & GoHorLoH=="DeNovo")$ref, "alt"= subset(allcol, Inheritance=="notinherited" & GoHorLoH=="DeNovo")$alt)
allcol_somatic_notinherited_LOH_dndscv<-data.frame("sampleID"= subset(allcol, Inheritance=="notinherited" & GoHorLoH=="LoH")$sample,"chr"= subset(allcol, Inheritance=="notinherited" & GoHorLoH=="LoH")$chrom,"pos"= subset(allcol, Inheritance=="notinherited" & GoHorLoH=="LoH")$pos, "ref" = subset(allcol, Inheritance=="notinherited" & GoHorLoH=="LoH")$ref, "alt"= subset(allcol, Inheritance=="notinherited" & GoHorLoH=="LoH")$alt)
#write.table(allcol_somatic_notinherited_GOH_dndscv, file="~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesnotinherited_GOH_20201210_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)  
#write.table(allcol_somatic_notinherited_LOH_dndscv, file="~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesnotinherited_LOH_20201210_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)  

CA56_uglm_spec<-spectrum(uniqueuglmCA56df,nrow(uniqueuglmCA56df0))
CA60_uglm_spec<-spectrum(uniqueuglmCA60df,nrow(uniqueuglmCA60df0))
CA65_uglm_spec<-spectrum(uniqueuglmCA65df,nrow(uniqueuglmCA65df))
allcol_uglm<-rbind(uniqueuglmCA56df, uniqueuglmCA60df, uniqueuglmCA65df)
allcol_uglm_GOH<- subset(allcol_uglm, GoH_or_LoH=="DeNovo")
allcol_uglm_GOH_dndscv<-data.frame("sampleID"= allcol_uglm_GOH$sample,"chr"= allcol_uglm_GOH$chrom,"pos"= allcol_uglm_GOH$pos, "ref" = allcol_uglm_GOH$ref, "alt"= allcol_uglm_GOH$alt)
#write.table(allcol_uglm_GOH_dndscv,file="~/Documents/GitHub/CoralGermline/dndscv/uglmallcolonies_GOH_20201210_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)

allcol_uglm_LOH<- subset(allcol_uglm, GoH_or_LoH=="LoH")
allcol_uglm_LOH_dndscv<-data.frame("sampleID"= allcol_uglm_LOH$sample,"chr"= allcol_uglm_LOH$chrom,"pos"= allcol_uglm_LOH$pos, "ref" = allcol_uglm_LOH$ref, "alt"= allcol_uglm_LOH$alt)
#write.table(allcol_uglm_LOH_dndscv,file="~/Documents/GitHub/CoralGermline/dndscv/uglmallcolonies_LOH_20201210_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)


allcol_uglm_spec<-spectrum(allcol_uglm, nrow(allcol_uglm)) + ggtitle("Single Sperm Pool Only SNVs") 
allcol_uglm_dndscv<-data.frame("sampleID"= allcol_uglm$sample,"chr"= allcol_uglm$chrom,"pos"= allcol_uglm$pos, "ref" = allcol_uglm$ref, "alt"= allcol_uglm$alt)
#write.table(allcol_uglm_dndscv,file="~/Documents/GitHub/CoralGermline/dndscv/uglmallcolonies_20201210_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)

allcol_gglm<-rbind(gglmCA56df, gglmCA60df, gglmCA65df)
allcol_gglm_spec<-spectrum(allcol_gglm, nrow(allcol_gglm)) + ggtitle("All Sperm Pool SNVs") 
allcol_gglm_dndscv<-data.frame("sampleID"= allcol_gglm$sample,"chr"= allcol_gglm$chrom,"pos"= allcol_gglm$pos, "ref" = allcol_gglm$ref, "alt"= allcol_gglm$alt)
#write.table(allcol_gglm_dndscv,file="~/Documents/GitHub/CoralGermline/dndscv/gglmallcolonies_20201210_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)
allcol_gglm_LOH<- subset(allcol_gglm, GoH_or_LoH=="LoH")
allcol_gglm_GOH<- subset(allcol_gglm, GoH_or_LoH=="DeNovo")

#this is figure 6:
(allcol_somatic_all | allcol_somatic_inherited | allcol_somatic_notinherited) / (allcol_uglm_spec | allcol_gglm_spec)
allcol_somatic_all | allcol_somatic_GOH | allcol_somatic_LOH
##end figure 6##
suppfig2<-(allcol_somatic_inherited | allcol_somatic_notinherited) / (allcol_uglm_spec | allcol_gglm_spec) / (allcol_somatic_GOH | allcol_somatic_LOH)

##chisquares:## return(counts) in the spectrum function!
#If you give chisq.test a matrix of counts, it's going to assume that it's a matrix of observed counts and will perform a chi-square test of independence. If you want to perform a chi-square goodness of fit test, which seems to be what you want to do, then you would need to use the command

#chisq.test(x,p)
#where x is the vector of observed counts c(11, 44233), and p is the vector of probabilities from the null hypothesis that were used to calculated the expected counts.
allcol_si_counts<-spectrum(subset(allcol, Inheritance=="inherited"), nrow(subset(allcol, Inheritance=="inherited")))
allcol_sn_counts<-spectrum(subset(allcol, Inheritance=="notinherited"), nrow(subset(allcol, Inheritance=="notinherited")))
allcol_uglm_counts<-spectrum(allcol_uglm, nrow(allcol_uglm))
allcol_gglm_counts<-spectrum(allcol_gglm, nrow(allcolgglm))
allcol_GOH_counts<-spectrum(subset(allcol, GoHorLoH=="DeNovo"),nrow(subset(allcol, GoHorLoH=="DeNovo")))
allcol_LOH_counts<-spectrum(subset(allcol, GoHorLoH=="LoH"),nrow(subset(allcol, GoHorLoH=="LoH")))
i_LOH<-spectrum(subset(allcol, Inheritance=="inherited" & GoHorLoH=="LoH"), nrow(subset(allcol, Inheritance=="inherited" & GoHorLoH=="LoH")))
i_GOH<-spectrum(subset(allcol, Inheritance=="inherited" & GoHorLoH=="DeNovo"), nrow(subset(allcol, Inheritance=="inherited" & GoHorLoH=="DeNovo")))

n_LOH<-spectrum(subset(allcol, Inheritance=="notinherited" & GoHorLoH=="LoH"), nrow(subset(allcol, Inheritance=="notinherited" & GoHorLoH=="LoH")))
n_GOH<-spectrum(subset(allcol, Inheritance=="notinherited" & GoHorLoH=="DeNovo"), nrow(subset(allcol, Inheritance=="notinherited" & GoHorLoH=="DeNovo")))

uglm_LOH<-spectrum(subset(allcol_uglm, GoH_or_LoH=="LoH"), nrow(subset(allcol_uglm,GoH_or_LoH=="LoH")))
uglm_GOH<-spectrum(subset(allcol_uglm, GoH_or_LoH=="DeNovo"), nrow(subset(allcol_uglm, GoH_or_LoH=="DeNovo")))
gglm_LOH<-spectrum(subset(allcol_gglm, GoH_or_LoH=="LoH"), nrow(subset(allcol_gglm,GoH_or_LoH=="LoH")))
gglm_GOH<-spectrum(subset(allcol_gglm, GoH_or_LoH=="DeNovo"), nrow(subset(allcol_gglm, GoH_or_LoH=="DeNovo")))

si_v_sn<-data.frame(allcol_si_counts,allcol_sn_counts)
uglm_v_gglm<-data.frame(allcol_uglm_counts,allcol_gglm_counts)
GOH_v_LOH<-data.frame(allcol_GOH_counts, allcol_LOH_counts)
chisq.test(si_v_sn)
sum(si_v_sn)
chisq.test(uglm_v_gglm)
sum(uglm_v_gglm)
sum(uglm_v_gglm)
chisq.test(GOH_v_LOH)
sum(GOH_v_LOH)

i_GOH_LOH<-data.frame(i_LOH, i_GOH)
chisq.test(i_GOH_LOH)
n_GOH_LOH<-data.frame(n_LOH, n_GOH)
chisq.test(n_GOH_LOH)
uglm_GOH_LOH<-data.frame(uglm_LOH, uglm_GOH)
chisq.test(uglm_GOH_LOH)
gglm_GOH_LOH<-data.frame(gglm_LOH, gglm_GOH)
chisq.test(gglm_GOH_LOH)

uniquesomaticCA56_spec<-spectrum(uniquesomaticCA56,nrow(uniquesomaticCA56))
uniquesomaticCA60_spec<-spectrum(uniquesomaticCA56,nrow(uniquesomaticCA60))
nonuniquesomaticCA65_spec<-spectrum(uniquesomaticCA65,nrow(uniquesomaticCA65))
#somatic_allcolonies_spec<-spectrum(somaticmetadatadf,nrow(somaticmetadatadf))
uniquesomatic_allcolonies_spec<-spectrum(uniquesomaticmetadatadf,nrow(uniquesomaticmetadatadf))

uniquesomaticCA56_spec  | uniquesomaticCA60_spec | nonuniquesomaticCA65_spec | uniquesomatic_allcolonies_spec
somatic_allcolonies_spec | uniquesomatic_allcolonies_spec
unique_inherited_spectrum<-spectrum(unique_inherited, nrow(unique_inherited))
unique_notinherited_spectrum<-spectrum(unique_notinherited, nrow(unique_notinherited))

unique_somatic_spectrum<-spectrum(unique_somatic, nrow(unique_somatic))
unique_uniqueglm_spectrum<-spectrum(unique_uniqueglm, nrow(unique_uniqueglm))
unique_globalglm_spectrum<-spectrum(unique_globalglm, nrow(unique_globalglm))
som<-unique_somatic_spectrum + ggtitle("Somatic Mutations") 
uglm<-unique_uniqueglm_spectrum + ggtitle("Unique Germline Mutations") 
gglm<-unique_globalglm_spectrum  + ggtitle("Global Germline Mutations")
(unique_inherited_spectrum + unique_notinherited_spectrum + unique_somatic_spectrum)# / (unique_uniqueglm_spectrum + unique_globalglm_spectrum)
unique_somatic_spectrum + unique_uniqueglm_spectrum + unique_globalglm_spectrum
som + uglm +gglm

inherited<-spectrum(inheritedallcolonies,nrow(inheritedallcolonies)) +ggtitle("Inherited Somatic Mutations")
notinherited<-spectrum(notinheritedallcolonies, nrow(notinheritedallcolonies)) +ggtitle("Non-inherited Somatic Mutations")
som + inherited + notinherited

uglm<-spectrum(unique_uniqueglm, nrow(unique_uniqueglm))
gglm<-spectrum(unique_globalglm, nrow(unique_globalglm))
ivsni<-data.frame(inherited,notinherited)
ivsgglm<-data.frame(inherited,gglm)
ivsuglm<-data.frame(inherited,uglm)
gglmvsuglm<-data.frame(gglm, uglm)
somvsuglm<-data.frame(unique_somatic_spectrum, uglm)
somvsgglm<-data.frame(unique_somatic_spectrum, gglm)
CA56vs65<-data.frame(CA56_somatic_all,CA65_somatic_all)
chisq.test(CA56vs65)

CA56vs65<-data.frame(CA60_somatic_notinherited,CA56_somatic_notinherited)
chisq.test(CA56vs65)

chi<- chisq.test(ivsni)  
chi_ivsgglm<-chisq.test(ivsgglm)
chi_ivsuglm<-chisq.test(ivsuglm)
chi_uglmvsgglm<-chisq.test(gglmvsuglm)
chi_somvsuglm<-chisq.test(somvsuglm)
chi_somvsgglm<-chisq.test(somvsgglm)

inherited | notinherited | uglm | gglm   

CA56<-bigfunction(CA56samples)
CA60<-bigfunction(CA60samples)
CA65<-bigfunction(CA65samples)

CA56 | CA60 | CA65

s1<-spectrum(missense, missensecount)
s2<-spectrum(syn,syncount)
s3<-spectrum(uniquemetadatadf,uniquemetadatadfcount)
s4<-spectrum(inherited,inheritedcount)
s5<-spectrum(not_inherited,not_inheritedcount)
s6<-spectrum(GOH, GOHcount)
s7<-spectrum(LOH, LOHcount)

library(patchwork)

bound<-rbind(s1,s2,s3)
subset<-c(rep("missense",nrow(s1)), rep("syn",nrow(s2)),rep("unique",nrow(s3)))
alltogether<-data.frame(bound,subset)


p1<- ggplot(alltogether, aes(x=Types, y=coralProportion)) + 
  geom_point(aes(colour=subset), size=3, position = position_dodge(width=0.3)) + 
  ylim(0, 0.4) +
  theme_bw()  
  #geom_jitter(width = 0.1, height = 0)
  #position_dodge(width = NULL, preserve = c("total", "single"))
p1
p1 + geom_jitter(width = 0.1, height = 0)

  

inherited<-subset(uniquemetadatadf, TrueorFalse=="True")
inheritedcount<-nrow(inherited)

not_inherited<-subset(uniquemetadatadf, TrueorFalse=="False")
not_inheritedcount<-nrow(not_inherited)

GOH<-subset(uniquemetadatadf, GoH_or_LoH =="DeNovo")
GOHcount<-nrow(GOH)
LOH<-subset(uniquemetadatadf, GoH_or_LoH == "LoH")
LOHcount<-nrow(LOH)

##Transiton/Transversion ratio ###
Transtions<- verATGCprop + verGCATprop
Transversions<- verATTAprop+verACTGprop+verGCCGprop+verGCTAprop
TiTv<- Transtions/Transversions

