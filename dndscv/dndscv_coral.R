#dndscv
setwd("~/Documents/GitHub/CoralGermline/dndscv")
library(dndscv)
?buildref()
buildref("AmilCDStable20200225.txt", "~/Downloads/Amil_v2.01/Amil.v2.01.chrs.fasta", "AmilRefCDS20200226.rda", numcode = 1, excludechrs= NULL)
filepath<-"CAcolony56_allmuts_dndscv.txt"
dnds_func<-function(filepath) {
  dataset<-read.delim(filepath)
  base<-basename(filepath)
  colony<-strsplit(base, "\\_")[[1]][1]

  dndsout = dndscv(dataset, refdb="AmilRefCDS20200226.rda", cv=NULL)
  
  global<-dndsout$globaldnds
  sel_loc<-dndsout$sel_loc
  annot_muts<-dndsout$annotmuts
  codon_sub<- annot_muts$codonsub
  codon_sub<-factor(codon_sub)
  return(codon_sub)
  dataglobal<-data.frame(global)
  mle_score<-global$mle
  wall<-subset(global,name=="wall")
  wall<-data.frame(wall)
  return(wall)
}  
# use rbind to put all the unique mutations from one colony into a single file,then run dnds_func
CAP12<-read.delim("CAcolony56CAP12_dndscv.txt")
CAP6<-read.delim("CAcolony56CAP8_dndscv.txt")
CAP8<-read.delim("CAcolony56CAP6_dndscv.txt")
CAP22<-read.delim("CAcolony60CAP22_dndscv.txt")
CAP23<-read.delim("CAcolony60CAP23_dndscv.txt")
CAP24<-read.delim("CAcolony60CAP24_dndscv.txt")
CAP11<-read.delim("CAcolony65CAP11_dndscv.txt")
CAP9<-read.delim("CAcolony65CAP9_dndscv.txt")

CAP22_i<-read.delim("CAcolony60CAP22_inherited_dndscv.txt")
CAP23_i<-read.delim("CAcolony60CAP23_inherited_dndscv.txt")
CAP24_i<-read.delim("CAcolony60CAP24_inherited_dndscv.txt")

CAP22_n<-read.delim("CAcolony60CAP22_notinherited_dndscv.txt")
CAP23_n<-read.delim("CAcolony60CAP23_notinherited_dndscv.txt")
CAP24_n<-read.delim("CAcolony60CAP24_notinherited_dndscv.txt")

CAP12_i<-read.delim("CAcolony56CAP12_inherited_dndscv.txt")
CAP6_i<-read.delim("CAcolony56CAP6_inherited_dndscv.txt")
CAP8_i<-read.delim("CAcolony56CAP8_inherited_dndscv.txt")

CAP12_n<-read.delim("CAcolony56CAP12_notinherited_dndscv.txt")
CAP6_n<-read.delim("CAcolony56CAP6_notinherited_dndscv.txt")
CAP8_n<-read.delim("CAcolony56CAP8_notinherited_dndscv.txt")

CAP11_i<-read.delim("CAcolony65CAP11_inherited_dndscv.txt")
CAP9_i<-read.delim("CAcolony65CAP9_inherited_dndscv.txt")

CAP11_n<-read.delim("CAcolony65CAP11_notinherited_dndscv.txt")
CAP9_n<-read.delim("CAcolony65CAP9_notinherited_dndscv.txt")

CA56<-rbind(CAP12,CAP6,CAP8)
CA56_i<-rbind(CAP12_i,CAP6_i, CAP8_i)
CA56_n<-rbind(CAP12_n,CAP6_n, CAP8_n)

CA60<-rbind(CAP22, CAP23,CAP24)
CA65<-rbind(CAP11, CAP9)
CA65_i<-rbind(CAP11_i,CAP9_i)
CA65_n<-rbind(CAP11_n, CAP9_n)

CA60_inherited<-rbind(CAP22_i, CAP23_i,CAP24_i)
CA60_notinherited<-rbind(CAP22_n, CAP23_n, CAP24_n)

write.table(CA56, file="~/Documents/GitHub/CoralGermline/dndscv/CAcolony56_allmuts_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)
write.table(CA56_i, file="~/Documents/GitHub/CoralGermline/dndscv/CAcolony56_allmuts_inherited_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)
write.table(CA56_n, file="~/Documents/GitHub/CoralGermline/dndscv/CAcolony56_allmuts_notinherited_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)

write.table(CA60, file="~/Documents/GitHub/CoralGermline/dndscv/CAcolony60_allmuts_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)
write.table(CA65, file="~/Documents/GitHub/CoralGermline/dndscv/CAcolony65_allmuts_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)
write.table(CA65_i, file="~/Documents/GitHub/CoralGermline/dndscv/CAcolony65_allmuts_inherited_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)
write.table(CA65_n, file="~/Documents/GitHub/CoralGermline/dndscv/CAcolony65_allmuts_notinherited_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)

write.table(CA60_inherited, file="~/Documents/GitHub/CoralGermline/dndscv/CAcolony60_allmuts_inherited_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)
write.table(CA60_notinherited, file="~/Documents/GitHub/CoralGermline/dndscv/CAcolony60_allmuts_notinherited_dndscv.txt",sep="\t",quote=FALSE, row.name=FALSE)

CA56_dnds<-dnds_func("CAcolony56_allmuts_dndscv.txt")
CA56_i_dnds<-dnds_func("CAcolony56_allmuts_inherited_dndscv.txt")
CA56_n_dnds<-dnds_func("CAcolony56_allmuts_notinherited_dndscv.txt")

CA60_dnds<-dnds_func("CAcolony60_allmuts_dndscv.txt")
CA65_dnds<-dnds_func("CAcolony65_allmuts_dndscv.txt")
CA65_i_dnds<-dnds_func("CAcolony65_allmuts_inherited_dndscv.txt")
CA65_n_dnds<-dnds_func("CAcolony65_allmuts_notinherited_dndscv.txt")

CA60_inherited_dnds<-dnds_func("CAcolony65_allmuts_inherited_dndscv.txt")
CA60_notinherited_dnds<-dnds_func("CAcolony65_allmuts_notinherited_dndscv.txt")

wall<-rbind(CA56_dnds, CA56_i_dnds, CA56_n_dnds, CA60_dnds, CA60_inherited_dnds, CA60_notinherited_dnds, CA65_dnds, CA65_i_dnds, CA65_n_dnds)
colony<-c("CA56", "CA56_i","CA56_n","CA60", "CA60_i","CA60_n", "CA65", "CA65_i", "CA65_n")#, "CA60", "CA60","CA60", "CA65", "CA65","CA65_i","CA65_i","CA65_i")
type<-c("all","i","n","all","i","n","all","i","n")
wall<-data.frame(wall, type)
mle<-wall$mle
cilow<-wall$cilow
cihigh<-wall$cihigh
barplot(mle)
boxplot(mle~colony, xlab="Colony Name", ylab="global dN/dS", las=1, col=c("light blue","light gray","light green"))
dndsplot<-ggplot(wall, 
       aes(x = colony, 
           y = mle, 
           group = 1
           )) +
  geom_point(aes(colour=type),size = 3) +
  theme_bw() +
  labs(x = "", y = "dN/dS") +
  ylim(0, 2) +
  geom_errorbar(aes(ymin = cilow, 
                    ymax = cihigh), 
                width = .1) +
  geom_hline(yintercept =1, linetype="dashed", color = "red") +
  scale_colour_manual(values = c("#999999", "#E69F00","red"))


#for each individual sample:
CA56CAP12<-dnds_func("CAcolony56CAP12_dndscv.txt")
CA56CAP6<-dnds_func("CAcolony56CAP6_dndscv.txt")
CA56CAP8<-dnds_func("CAcolony56CAP8_dndscv.txt")
CA60CAP22<-dnds_func("CAcolony60CAP22_dndscv.txt")
CA60CAP23<-dnds_func("CAcolony60CAP23_dndscv.txt")
CA60CAP24<-dnds_func("CAcolony60CAP24_dndscv.txt")
CA65CAP11<-dnds_func("CAcolony65CAP11_dndscv.txt")
CA65CAP9<-dnds_func("CAcolony65CAP9_dndscv.txt")

wall<-rbind(CA56CAP12, CA56CAP6, CA56CAP8,CA60CAP22, CA60CAP23, CA60CAP24, CA65CAP11, CA65CAP9)

colony<-c("CA56", "CA56", "CA56", "CA60", "CA60","CA60", "CA65", "CA65")
mle<-wall$mle
cilow<-wall$cilow
cihigh<-wall$cihigh
barplot(mle)
boxplot(mle~colony, xlab="Colony Name", ylab="global dN/dS", las=1, col=c("light blue","light gray","light green"))
ggplot(wall, 
       aes(x = colony, 
           y = mle, 
           group = 1)) +
  geom_point(size = 3) +

  geom_errorbar(aes(ymin = cilow, 
                    ymax = cihigh), 
                width = .1)

#trinucleotide context:
dnds_func<-function(filepath) {
  dataset<-read.delim(filepath)
  base<-basename(filepath)
  colony<-strsplit(base, "\\_")[[1]][1]
  
  dndsout = dndscv(dataset, refdb="AmilRefCDS20200226.rda", cv=NULL)
  
  return()
}
library(patchwork)
dndsplot/a2 

}  