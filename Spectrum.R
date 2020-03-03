library("seqinr")
Amil_fasta<-read.fasta("~/Documents/GitHub/CoralGermline/dndscv/Amil.v2.01.chrs.fasta",forceDNAtolower = FALSE)
files<-list.files(path="~/Documents/GitHub/CoralGermline/annotatedfiles", pattern="*ann.txt", full.names=T, recursive=FALSE) #path to all the files you want to include in the analysis

CA56_samples<-list.files(path="~/Documents/GitHub/CoralGermline/annotatedfiles", pattern="*relabeled56.txt.ann.txt", full.names=T, recursive=FALSE)
files<-list.files(path="~/Documents/GitHub/CoralGermline/annotatedfiles", pattern="*relabeled60.txt.ann.txt", full.names=T, recursive=FALSE)
CA<-list.files(path="~/Documents/GitHub/CoralGermline/annotatedfiles", pattern="*relabeled65.txt.ann.txt", full.names=T, recursive=FALSE)

bigfunction<- function(files) {
#files<-list.files(path="~/Documents/GitHub/CoralGermline/annotatedfiles", pattern="*ann.txt", full.names=T, recursive=FALSE) #path to all the files you want to include in the analysis
  metadata= NULL
  for (i in 1:length(files)) { 
    file =files[i]
    data<-read.delim(file) #read in each file in "files"
    data<-data.frame(data) # transform the data from each file into a dataframe
    base<-basename(file)
    colony<-strsplit(base, "\\_")[[1]][2]
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
  refdepth<-as.numeric(split[,2])
  altdepth<-as.numeric(split[,3])
  what<-metadata$WhattoWhat
  allelesplit<-str_split_fixed(what, "to", 2)
  normalallele<-allelesplit[,1]
  mutantallele<-allelesplit[,2]
  totaldepth<-refdepth+altdepth
  GQscore<-as.numeric(split[,4])
  mutationtype<-metadata$MutationType
  mutationstrength<-metadata$MutationStrength
  
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
  
  metadatadf<-data.frame("chrom.pos" = metadata$chrom.pos, "chrom"=chr, "pos"=pos,	"sample"= metadata$sample, "ref" = metadata$ref, "alt" = 	metadata$alt, "normal_allele"= normalallele, "mutant_allele" = mutantallele, "mutant_allele_depth" = as.numeric(mutant_alleledepth), "genotype"= genotypes, "totaldepth"=totaldepth, 	"refdepth"=refdepth, "altdepth"=altdepth, "GQscore"= GQscore,	"GoH_or_LoH"=metadata$DeNovo_LoH, "Ti/Tv"=metadata$TiTv, 	"WhattoWhat" = metadata$WhattoWhat,"TrueorFalse" =metadata$TrueorFalse, "MutationType"=mutationtype, "MutationStrength"=mutationstrength)#  ColonyName"=metadata$colonyrep)
  
  DepthMeansdf<-aggregate(totaldepth~chrom.pos, metadatadf, 			FUN=mean)
  
  DepthMinsdf<-aggregate(totaldepth~chrom.pos, metadatadf, 			FUN=min)
  
  GQaverage<-aggregate(GQscore~chrom.pos, metadatadf, FUN=mean)
  
  GQmin<-aggregate(GQscore~chrom.pos, metadatadf, FUN=min)
  
  
  
  metadatadf.00<-merge(metadatadf, DepthMeansdf[, c("chrom.pos", 	"totaldepth")], by="chrom.pos")
  metadatadf.0<-merge(metadatadf.00, GQaverage[,c("chrom.pos","GQscore")], by="chrom.pos")
  metadatadf.0<-merge(metadatadf.0, DepthMinsdf[,c("chrom.pos","totaldepth")], by="chrom.pos")
  metadatadf.0<-merge(metadatadf.0, GQmin[,c("chrom.pos","GQscore")], by="chrom.pos")
  
  DeNovos<-subset(metadatadf.0, GoH_or_LoH=="DeNovo")
  # sample3<-subset(DeNovos, sample== "sample3")
  # trueDenovos_sample3<-subset(sample3, refdepth =="0" | altdepth=="0")
  # 
  # sample4<-subset(DeNovos, sample== "sample4")
  # trueDenovos_sample4<-subset(sample4, refdepth =="0" | altdepth=="0")
  # 
  # sample5<-subset(DeNovos, sample== "sample5")
  # trueDenovos_sample5<-subset(sample5, refdepth =="0" | altdepth=="0")
  # 
  # sample6<-subset(DeNovos, sample== "sample6")
  # trueDenovos_sample6<-subset(sample6, refdepth =="0" | altdepth=="0")
  # 
  # sample7<-subset(DeNovos, sample== "sample7")
  # trueDenovos_sample7<-subset(sample7, refdepth =="0" | altdepth=="0")
  # 
  # sample8<-subset(DeNovos, sample== "sample8")
  # trueDenovos_sample8<-subset(sample8, refdepth =="0" | altdepth=="0")
  # 
  # truedenovos3_8<-rbind(trueDenovos_sample3, trueDenovos_sample4, trueDenovos_sample5, trueDenovos_sample6, trueDenovos_sample7, trueDenovos_sample8)
  
  LoH<-subset(metadatadf.0, GoH_or_LoH =="LoH")#
  trueLoHp<-subset(LoH,refdepth =="0" | altdepth=="0")
  trueLoHp1<-subset(trueLoHp, sample=="mutparent1" | sample=="mutparent2")
  trueLoHp2<-trueLoHp1[trueLoHp1$chrom.pos %in% names(which(table(trueLoHp1$chrom.pos) > 1)), ]
  sperm<-subset(metadatadf.0, sample =="mutsperm")
  trueLoHsperm<-sperm[match(trueLoHp2$chrom.pos, sperm$chrom.pos),]
  trueLoHsperm<-unique(trueLoHsperm)
  
  #trueLoHp2<-subset(trueLoHp, sample=="mutparent2")
  metadatadf<-rbind( DeNovos, trueLoHp2,trueLoHsperm)#, trueLoHp2)
  #write.table(metadatadf, file="CAcolony60_CAP22-23-24muts_20191125.txt",sep="\t",quote=FALSE, row.name=FALSE)
  
  uniquemetadatadf<- metadatadf[match(unique(metadatadf$chrom.pos), 					metadatadf$chrom.pos),]
  
  uniquemetadatadfcount<-nrow(uniquemetadatadf)
  
  missense<-subset(uniquemetadatadf, MutationType=="missense_variant")
  missensecount<-nrow(missense)
  
  syn<-subset(uniquemetadatadf, MutationType=="synonymous_variant")
  syncount<-nrow(syn)
  #par(mfrow=c(1,3))
  
  dataframe<-uniquemetadatadf
  dataframecount<-uniquemetadatadfcount
#spectrum<- function(dataframe, dataframecount) {
  verATGClist<-dataframe[ which(dataframe$WhattoWhat=="AtoG" | dataframe$WhattoWhat=="TtoC"),]
  verATGCcount<-nrow(verATGClist)
  verATGCprop<-verATGCcount/dataframecount
  SEverATGCprop<-se(verATGCprop)
  
  verGCATlist<-dataframe[ which(dataframe$WhattoWhat=="GtoA" | dataframe$WhattoWhat=="CtoT"),]
  verGCATcount<-nrow(verGCATlist)
  verGCATprop<-verGCATcount/dataframecount
  
  CpGdf<-NULL
  CpHdf<-NULL
  all<-NULL
  for (i in 1:nrow(dataframe)) {
    if (dataframe$WhattoWhat[i] == "CtoT") {
      character<-as.character(dataframe$pos[i])
      integer<-as.integer(character)
      next_position<-integer + 1
      chr_num<-as.character(dataframe$chrom[i])
      
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
    } else if (dataframe$WhattoWhat[i] == "GtoA") {   
      character<-as.character(dataframe$pos[i])
      integer<-as.integer(character)
      next_position<-integer + 1
      chr_num<-as.character(dataframe$chrom[i])
      
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
  
  CpGprop<- CpGcount/dataframecount
  CpHprop<- CpHcount/dataframecount
  
  verATTAlist<-dataframe[ which(dataframe$WhattoWhat=="AtoT" | dataframe$WhattoWhat=="TtoA"),]
  verATTAcount<-nrow(verATTAlist)
  verATTAprop<-verATTAcount/dataframecount
  
  verACTGlist<-dataframe[ which(dataframe$WhattoWhat=="AtoC" | dataframe$WhattoWhat=="TtoG"),]
  verACTGcount<-nrow(verACTGlist)
  verACTGprop<-verACTGcount/dataframecount
  
  verGCCGlist<-dataframe[ which(dataframe$WhattoWhat=="CtoG" | dataframe$WhattoWhat=="GtoC"),]
  verGCCGcount<-nrow(verGCCGlist)
  verGCCGprop<-verGCCGcount/dataframecount
  
  verGCTAlist<-dataframe[ which(dataframe$WhattoWhat=="GtoT" | dataframe$WhattoWhat=="CtoA"),]
  verGCTAcount<-nrow(verGCTAlist)
  verGCTAprop<-verGCTAcount/dataframecount
  vertypesDF<-data.frame(Types=c("A>G/T>C","CpG", "CpH", "A>T/T>A","A>C/T>G","G>C/C>G","G>T/C>A"), coralProportion=c(verATGCprop, CpGprop, CpHprop, verATTAprop,verACTGprop, verGCCGprop, verGCTAprop))
  
  #MUTATION SPECTRUM PLOT:
  vertypesDFplot<-barplot(vertypesDF$coralProportion,names.arg=vertypesDF$Types, ylim=c(0,0.7), main="allmuts", las=2) #Figure 4a
  z<- ggplot(vertypesDF, aes(x=Types, y=coralProportion)) +
    geom_bar(stat="identity") +
    ylim(0, 0.4)
  Transtions<- verATGCprop + verGCATprop
  Transversions<- verATTAprop+verACTGprop+verGCCGprop+verGCTAprop
  TiTv<- Transtions/Transversions
  #return(z)
  return(list(vertypesDF, TiTv))
  #return(list(z, vertypesDF))
}

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

