setwd("~/Documents/GitHub/CoralGermline/")
library("reshape2")
library(ggplot2)
library(stringr)
library(sciplot)
library(sinaplot)
library(ggforce)
library(gridExtra)
library(dplyr)

CAcolony60<-read.delim("CAcolony60_ii_20200324.txt",check.names = FALSE)
header<-colnames(CAcolony60)
samplenames<-header[4:19]
vas<-c(samplenames)
CAcolony60mdf<-melt(CAcolony60, id.vars="chrom.pos", measure.vars=c(samplenames), value.name="genotype",variable.name="sample")
write.table(CAcolony60mdf, file="meltedCAcolony60_ii_20200324.txt",sep="\t",quote=FALSE, row.name=FALSE)

CAcolony56<-read.delim("CAcolony56_ii_20200324.txt", check.names = FALSE)
header<-colnames(CAcolony56)
samplenames<-header[4:15]
CAcolony56mdf<-melt(CAcolony56, id.vars="chrom.pos", measure.vars=c(samplenames), value.name="genotype",variable.name="sample")
write.table(CAcolony56mdf, file="meltedCAcolony56_ii_20200324.txt",sep="\t",quote=FALSE, row.name=FALSE)

CAcolony56<-read.delim("CAcolony56_noCAS12-2_ii_20200413.txt", check.names = FALSE)
header<-colnames(CAcolony56)
samplenames<-header[4:14]
#CAcolony56noNA<-subset(CAcolony56, TypeofMutation!= "NA")
CAcolony56mdf<-melt(CAcolony56noNA, id.vars="chrom.pos", measure.vars=c(samplenames), value.name="genotype",variable.name="sample")
write.table(CAcolony56mdf, file="meltedCAcolony56_ii_20200413.txt",sep="\t",quote=FALSE, row.name=FALSE)

CAcolony65<-read.delim("CAcolony65_ii_20200331.txt", check.names = FALSE)
header<-colnames(CAcolony65)
samplenames<-header[4:15]
CAcolony65mdf<-melt(CAcolony65, id.vars="chrom.pos", measure.vars=c(samplenames), value.name="genotype",variable.name="sample")
write.table(CAcolony65mdf, file="meltedCAcolony65_ii_20200331.txt",sep="\t",quote=FALSE, row.name=FALSE)

CAcolony65<-read.delim("CAcolony65_noCAS10-2_ii_20200414_2.txt", check.names = FALSE)
header<-colnames(CAcolony65)
samplenames<-header[4:14]
CAcolony65mdf<-melt(CAcolony65, id.vars="chrom.pos", measure.vars=c(samplenames), value.name="genotype",variable.name="sample")
write.table(CAcolony65mdf, file="meltedCAcolony65_noCAS10-2_ii_20200414_2.txt",sep="\t",quote=FALSE, row.name=FALSE)


somaticCA56<-subset(CAcolony56, TypeofMutation == "SomaticMutation")
somaticCAP6<-subset(somaticCA56, MutantSampleID=="CAP6-1_S47")
somaticCAP8<-subset(somaticCA56, MutantSampleID=="CAP8-1_S40")
somaticCAP12<-subset(somaticCA56, MutantSampleID=="CAP12-1_S46")
somatic_countpersample<-c(nrow(somaticCAP6), nrow(somaticCAP8),nrow(somaticCAP12))
uniqueCA56<-subset(CAcolony56, TypeofMutation == "UniqueGermlineMutation")
uniqueCAS6<-subset(uniqueCA56, MutantSampleID=="CAS6-2_S9")
uniqueCAS8<-subset(uniqueCA56, MutantSampleID=="CAS8-2_S2")
uniqueCAS12<-subset(uniqueCA56, MutantSampleID=="CAS12-2_S8")
write.table(uniqueCAS8, file="uniqueCAS8.txt", sep="\t",quote=FALSE, row.name=FALSE)
unique_countpersample<-c(nrow(uniqueCAS6),nrow(uniqueCAS8), nrow(uniqueCAS12))
globalCA56<-subset(CAcolony56, TypeofMutation == "GlobalGermlineMutation")

counts<-c(somatic_countpersample, unique_countpersample, nrow(globalCA56))

somaticCA60<-subset(CAcolony60, TypeofMutation == "SomaticMutation")
somaticCAP22<-subset(somaticCA60, MutantSampleID=="CAP22-1_S43")
somaticCAP23<-subset(somaticCA60, MutantSampleID=="CAP23-1_S45")
somaticCAP24<-subset(somaticCA60, MutantSampleID=="CAP24-1_S42")
somaticCAP26<-subset(somaticCA60, MutantSampleID=="CAP26-1_S39")

somatic_countpersample<-c(nrow(somaticCAP22), nrow(somaticCAP23),nrow(somaticCAP24), nrow(somaticCAP26))

uniqueCA60<-subset(CAcolony60, TypeofMutation == "UniqueGermlineMutation")
uniqueCAS22<-subset(uniqueCA60, MutantSampleID=="CAS22-2_S5")
uniqueCAS23<-subset(uniqueCA60, MutantSampleID=="CAS23-2_S7")
uniqueCAS24<-subset(uniqueCA60, MutantSampleID=="CAS24-2_S4")
uniqueCAS26<-subset(somaticCA60, MutantSampleID=="CAS26-2_S1")

#write.table(uniqueCAS12, file="uniqueCAS12.txt", sep="\t",quote=FALSE, row.name=FALSE)
unique_countpersample<-c(nrow(uniqueCAS22),nrow(uniqueCAS23), nrow(uniqueCAS24), nrow(uniqueCAS26))
globalCA60<-subset(CAcolony60, TypeofMutation == "GlobalGermlineMutation")

counts<-c(somatic_countpersample, unique_countpersample, nrow(globalCA60))
