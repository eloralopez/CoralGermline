#ahya ref:
threedfs_CA56<-threedfs_func(list.files(path="~/Documents/CoralGermlineManuscript", pattern="*56_dm_onelib_20210409_1_chrs.txt", full.names=T, recursive=FALSE))

somaticCA56df<-threedfs_CA56$somatic
uglmCA56df<-threedfs_CA56$uglm
uniqueuglmCA56df<-uglmCA56df[match(unique(uglmCA56df$chrom.pos), 					uglmCA56df$chrom.pos),]
gglmCA56df<-threedfs_CA56$gglm
metadatadf.0CA56<-threedfs_CA56$metadatadf.0

#threedfs_CA60<-threedfs_func(list.files(path="~/Documents/GitHub/CoralGermline/cleanpipeline", pattern="*60_dm_20201210.txt.ann.txt", full.names=T, recursive=FALSE))
threedfs_CA60<-threedfs_func(list.files(path="~/Documents/CoralGermlineManuscript", pattern="*60_dm_onelib_20210409_1_chrs.txt", full.names=T, recursive=FALSE))

somaticCA60df<-threedfs_CA60$somatic
uglmCA60df<-threedfs_CA60$uglm
uniqueuglmCA60df<-uglmCA60df[match(unique(uglmCA60df$chrom.pos), 					uglmCA60df$chrom.pos),]
gglmCA60df<-threedfs_CA60$gglm
metadatadf.0CA60<-threedfs_CA60$metadatadf.0

#threedfs_CA65<-threedfs_func(list.files(path="~/Documents/GitHub/CoralGermline/cleanpipeline", pattern="*65_dm_20201210.txt.ann.txt", full.names=T, recursive=FALSE))
threedfs_CA65<-threedfs_func(list.files(path="~/Documents/CoralGermlineManuscript", pattern="*65_dm_onelib_20210409_1_chrs.txt", full.names=T, recursive=FALSE))

somaticCA65df<-threedfs_CA65$somatic
uglmCA65df<-threedfs_CA65$uglm
uniqueuglmCA65df<-uglmCA65df[match(unique(uglmCA65df$chrom.pos), 					uglmCA65df$chrom.pos),]
gglmCA65df<-threedfs_CA65$gglm
metadatadf.0CA65<-threedfs_CA65$metadatadf.0

twolibs<-read.delim("~/Documents/CoralGermlineManuscript/SupplementaryTable2_20210625.txt")
twolibs_CA56_somatic<-subset(twolibs, ColonyName=="CA56" & MutationClass=="SomaticMutation")
onelib_CA56_somatic<-subset(everything, ColonyName=="CA56" & MutationClass =="SomaticMutation")
CA56_matched<-na.omit(onelib_CA56_somatic[match(onelib_CA56_somatic$chrom.pos, twolibs_CA56_somatic$chrom.pos),])
uniqueonelib_CA56_somatic<-onelib_CA56_somatic[match(unique(onelib_CA56_somatic$chrom.pos), onelib_CA56_somatic$chrom.pos),]
uniqueCA56_matched<-CA56_matched[match(unique(CA56_matched$chrom.pos), CA56_matched$chrom.pos),]

twolibs_CA60_somatic<-subset(twolibs, ColonyName=="CA60" & MutationClass=="SomaticMutation")
onelib_CA60_somatic<-subset(everything, ColonyName=="CA60" & MutationClass =="SomaticMutation")
CA60_matched<-na.omit(onelib_CA60_somatic[match(onelib_CA60_somatic$chrom.pos, twolibs_CA60_somatic$chrom.pos),])
uniqueonelib_CA60_somatic<-onelib_CA60_somatic[match(unique(onelib_CA60_somatic$chrom.pos), onelib_CA60_somatic$chrom.pos),]
uniqueCA60_matched<-CA60_matched[match(unique(CA60_matched$chrom.pos), CA60_matched$chrom.pos),]

twolibs_CA65_somatic<-subset(twolibs, ColonyName=="CA65" & MutationClass=="SomaticMutation")
onelib_CA65_somatic<-subset(everything, ColonyName=="CA65" & MutationClass =="SomaticMutation")
CA65_matched<-na.omit(onelib_CA65_somatic[match(onelib_CA65_somatic$chrom.pos, twolibs_CA65_somatic$chrom.pos),])
uniqueonelib_CA65_somatic<-onelib_CA65_somatic[match(unique(onelib_CA65_somatic$chrom.pos), onelib_CA65_somatic$chrom.pos),]
uniqueCA65_matched<-CA65_matched[match(unique(CA65_matched$chrom.pos), CA65_matched$chrom.pos),]

oneandtwo<-rbind(uniqueonelib_CA56_somatic, uniqueCA56_matched, uniqueonelib_CA60_somatic, uniqueCA60_matched, uniqueonelib_CA65_somatic, uniqueCA65_matched)
oneandtwo$Replicates<-factor(c(rep("without",nrow(uniqueonelib_CA56_somatic)),rep("with",nrow(uniqueCA56_matched)),
                 rep("without",nrow(uniqueonelib_CA60_somatic)),rep("with",nrow(uniqueCA60_matched)),
                 rep("without",nrow(uniqueonelib_CA65_somatic)),rep("with",nrow(uniqueCA65_matched))),levels=c("without","with"))
oneandtwo$lib2<-c(rep("one2",nrow(uniqueonelib_CA56_somatic)),rep("matched2",nrow(uniqueCA56_matched)),
                 rep("one2",nrow(uniqueonelib_CA60_somatic)),rep("matched2",nrow(uniqueCA60_matched)),
                 rep("one2",nrow(uniqueonelib_CA65_somatic)),rep("matched2",nrow(uniqueCA65_matched)))
###this is your supp fig boxplots:
avggq<-ggplot(oneandtwo,aes(x=ColonyName))+
  
  geom_sina(aes(color=Replicates, y=GQscore.y),size=1)+
  geom_boxplot(fill = NA, aes(color=lib2, y=GQscore.y),size=1 )+
  scale_color_manual(values = c("without"="#56B4E9", "with"="#D55E00", "one2" ="black","matched2" = "black"))+
  ylab("Average GQ score at locus")+
  theme_bw()+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25),legend.title = element_text(size=25),legend.text = element_text(size=25))

avgdp<-ggplot(oneandtwo,aes(x=ColonyName))+
  
  geom_sina(aes(color=Replicates, y=totaldepth.y),size=1)+
  geom_boxplot(fill = NA, aes(color=lib2, y=totaldepth.y),size=1 )+
  scale_color_manual(values = c("without"="#56B4E9", "with"="#D55E00", "one2" ="black","matched2" = "black"))+
  ylab("Average read depth at locus")+
  theme_bw()+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25),legend.title = element_text(size=25),legend.text = element_text(size=25))

mingq<-ggplot(oneandtwo,aes(x=ColonyName))+
  
  geom_sina(aes(color=Replicates, y=GQscore),size=1)+
  geom_boxplot(fill = NA, aes(color=lib2, y=GQscore),size=1 )+
  scale_color_manual(values = c("without"="#56B4E9", "with"="#D55E00", "one2" ="black","matched2" = "black"))+
  ylab("Lowest GQ score at locus")+
  theme_bw()+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25),legend.title = element_text(size=25),legend.text = element_text(size=25))
#mingq  
"#D55E00","#E69F00","#0072B2","#56B4E9"
mindp<-ggplot(oneandtwo,aes(x=ColonyName))+
  
  geom_sina(aes(color=Replicates, y=totaldepth),size=1)+
  geom_boxplot(fill = NA, aes(color=lib2, y=totaldepth),size=1 )+
  scale_color_manual(values = c("without"="#56B4E9", "with"="#D55E00", "one2" ="black","matched2" = "black"))+
  ylab("Lowest read depth at locus")+
  theme_bw()+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=25),legend.title = element_text(size=25),legend.text = element_text(size=25))
mindp
#suppfig4 20210720
(avggq +labs(tag = "a.")+theme(plot.tag=element_text(face = "bold",size=30))| avgdp +labs(tag = "b.")+theme(plot.tag=element_text(face = "bold",size=30)))/ (mingq +labs(tag = "c.")+theme(plot.tag=element_text(face = "bold",size=30))| mindp +labs(tag = "d.")+theme(plot.tag=element_text(face = "bold",size=30)))
#height = 1251, width = 2500
#CA56_matched_inverse<-na.omit(twolibs[match(twolibs$chrom.pos, metadatadf.0CA56$chrom.pos),])
twolibs_CA60<-subset(twolibs, ColonyName=="CA60")
CA60_matched<-na.omit(somaticCA60df[match(somaticCA60df$chrom.pos, twolibs_CA60$chrom.pos),])

twolibs_CA65<-subset(twolibs, ColonyName=="CA65")
CA65_matched<-na.omit(somaticCA65df[match(somaticCA65df$chrom.pos, twolibs_CA65$chrom.pos),])

hist(CA56_matched$GQscore)
hist(metadatadf.0CA56$GQscore)
hist(metadatadf.0CA60$GQscore)
hist(metadatadf.0CA65$GQscore)

###########################################################
allcolonies<-scatterplot_func(rbind(somaticCA56df,somaticCA60df,somaticCA65df))
allcolonies_twolib<-scatterplot_func(rbind(somaticCA56df_twolib, somaticCA60df_twolib, somaticCA65df_twolib))
allcolonies_inherited_GOH<-subset(allcolonies, Inheritance=="inherited" & GoHorLoH=="DeNovo")
allcolonies_notinherited_GOH<-subset(allcolonies, Inheritance=="notinherited" & GoHorLoH=="DeNovo")
VAFs<-data.frame("PVAF"=c(allcolonies_inherited_GOH$ParentAverage, allcolonies_notinherited_GOH$ParentAverage), 
                 "shared"=c(rep("Parent and Sperm",nrow(allcolonies_inherited_GOH)),rep("Parent Only",nrow(allcolonies_notinherited_GOH))))

allcolonies_twolib_inherited_GOH<-subset(allcolonies_twolib, Inheritance=="inherited" & GoHorLoH=="DeNovo")
allcolonies_twolib_notinherited_GOH<-subset(allcolonies_twolib, Inheritance=="notinherited" & GoHorLoH=="DeNovo")
VAFs_twolib<-data.frame("PVAF"=c(allcolonies_twolib_inherited_GOH$ParentAverage, allcolonies_twolib_notinherited_GOH$ParentAverage), 
                 "shared"=c(rep("Parent and Sperm",nrow(allcolonies_twolib_inherited_GOH)),rep("Parent Only",nrow(allcolonies_twolib_notinherited_GOH))))

vaf_one<-ggplot(data=VAFs, aes(x=PVAF, fill=shared)) +
  geom_density(alpha=0.5) +
  xlim(0.1,0.65) + ylim(0,23)+
  #scale_fill_manual(values=c("pink","green"))+
  theme_minimal()
vaf_two<-ggplot(data=VAFs_twolib, aes(x=PVAF, fill=shared)) +
  geom_density(alpha=0.5) +
  xlim(0.1,0.65) +ylim(0,23)+
  theme_minimal()
vaf_one | vaf_two
VAFs<-data.frame("PVAF"=c(allcolonies_inherited_GOH$ParentAverage, allcolonies_notinherited_GOH$ParentAverage), 
                 "shared"=c(rep("Parent and Sperm",nrow(allcolonies_inherited_GOH)),rep("Parent Only",nrow(allcolonies_notinherited_GOH))))
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
CA56_inherited<-subset(CA56_, Inheritance=="inherited")
CA56_notinherited<-subset(CA56_, Inheritance=="notinherited")
CA56_loh<-subset(CA56_, GoHorLoH=="LoH")
CA56_denovo<-subset(CA56_, GoHorLoH=="DeNovo")

inheriteddenovoCA56_<-subset(CA56_, SpermAverage>0 &ParentAverage>=0.1 &ParentAverage<1)
inheritedlohCA56_<-subset(CA56_, GoHorLoH=="LoH" & SpermAverage==1 & ParentAverage==1)
notinheriteddenovoCA56_<-subset(CA56_, GoHorLoH=="DeNovo" & SpermAverage==0)
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

CA56gglmGOH<-nrow(subset(gglmCA56df, startsWith(gglmCA56df$sample,'CAS6')==TRUE & GoH_or_LoH == "DeNovo" ) )
CA56gglmLOH<-nrow(subset(gglmCA56df, startsWith(gglmCA56df$sample,'CAS6')==TRUE & GoH_or_LoH == "LoH" ) )

CA60gglmGOH<-nrow(subset(gglmCA60df, startsWith(gglmCA60df$sample,'CAS22')==TRUE & GoH_or_LoH == "DeNovo" ) )
CA60gglmLOH<-nrow(subset(gglmCA60df, startsWith(gglmCA60df$sample,'CAS22')==TRUE & GoH_or_LoH == "LoH" ) )

CA65gglmGOH<-nrow(subset(gglmCA65df, startsWith(gglmCA65df$sample,'CAS11')==TRUE & GoH_or_LoH == "DeNovo" ) )
CA65gglmLOH<-nrow(subset(gglmCA65df, startsWith(gglmCA65df$sample,'CAS11')==TRUE & GoH_or_LoH == "LoH" ) )

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
typez2<-rep(c("GoH","LoH"),5)
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

df<-data.frame(subsetz, typez2, means, allse )
df$meanspercent<-100*df$means
df$allsepercent<-100*allse
df$shared<-factor(c(rep("All Parent",2),rep("P+S",2),rep("PO",2),rep("SSPO",2),rep("ASP",2)),
                  levels=c("All Parent","PO","P+S","SSPO","ASP"))
df_withoutallsom<-df[-(1:2),]
z2<- ggplot(df_withoutallsom, aes(x=shared, y=meanspercent,fill=typez2)) +
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  #scale_fill_brewer(palette = "Paired") +
  ylab("Percent of SNVs") + xlab("")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15))   +
  theme(legend.text=element_text(size=15),
        axis.text.x = element_text(angle = 90),
        legend.title=element_blank())+
  geom_errorbar(position="dodge",aes(ymin=meanspercent-(allsepercent), ymax=meanspercent+(allsepercent)))
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

onelib<-c(CA56n_GOH,CA56n_LOH, CA60n_GOH, CA60n_LOH, CA65n_GOH, CA65n_LOH, CA56i_GOH, CA56i_LOH, CA60i_GOH, CA60i_LOH, 
          CA65i_GOH, CA65i_LOH, CA56uglm_GOH, CA56uglm_LOH, CA60uglm_GOH, CA60uglm_LOH, CA65uglm_GOH, CA65uglm_LOH)
twolib<-c(47,
          61,
          192,
          243,
          66,
          19,
          61,
          160,
          88,
          112,
          47, 22, 129,60, 67,
          35, 17,
          22, 24,
          34,
          34, 17,
          9,          11, 60, 53, 14, 11, 64,65,65,51,41,47,26,29,44,43,18,28,19,22)
counts<-c(onelib,twolib)
lib<-c(rep("one", length(onelib)), rep("two",length(twolib)))
colony<-rep(c(rep("CA56",4),rep("CA60",6), rep("CA65",4)),6)
sample<-rep(c(rep(c("CAP6", "CAP8"),2) , rep(c("CAP22", "CAP23", "CAP26"),2) , rep(c("CAP11", "CAP9"),2)),6)
gl<-rep(c("GOH","GOH", "LOH", "LOH", "GOH", "GOH","GOH","LOH","LOH","LOH","GOH","GOH","LOH","LOH"),6)
percent<-(twolib/onelib)*100
inhe<-rep(c(rep("not shared",14),rep("shared",14),rep("uglm",14)),2)
countframe<-data.frame(counts, lib, colony,sample,gl,inhe, percent)
countframe$new<-factor(paste0(gl,inhe))
countframe$new2<-factor(paste0(countframe$new,sample))
colony
justshared<-subset(countframe, inhe=="shared")
jsplot<-ggplot(justshared, aes(x=new2, y=counts, fill=lib))+
  geom_bar(position="dodge", stat="identity")+
  theme_bw()+
  scale_y_continuous(trans='log10')

cfplot<-ggplot(countframe, aes(x=new2, y=counts, fill=lib))+
  geom_bar(position="dodge", stat="identity")+
  theme_bw()+
  scale_y_continuous(trans='log10')+
  theme(axis.text.x=element_text(size=8, angle = 90), axis.text.y=element_text(size=8, angle = 0),
        axis.title=element_text(size=8))+
  facet_wrap(~sample, strip.position = "bottom", scales = "free_x")+
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(size = 15))
cfplot
nsG1<-mean(subset(countframe, inhe=="not shared" & gl=="GOH" & lib=="one")$counts)
nsL1<-mean(subset(countframe, inhe=="not shared" & gl=="LOH" & lib=="one")$counts)
sG1<-mean(subset(countframe, inhe=="shared" & gl=="GOH" & lib=="one")$counts)
sL1<-mean(subset(countframe, inhe=="shared" & gl=="LOH" & lib=="one")$counts)
uG1<-mean(subset(countframe, inhe=="uglm" & gl=="GOH" & lib=="one")$counts)
uL1<-mean(subset(countframe, inhe=="uglm" & gl=="LOH" & lib=="one")$counts)
nsG2<-mean(subset(countframe, inhe=="not shared" & gl=="GOH" & lib=="two")$counts)
nsL2<-mean(subset(countframe, inhe=="not shared" & gl=="LOH" & lib=="two")$counts)
sG2<-mean(subset(countframe, inhe=="shared" & gl=="GOH" & lib=="two")$counts)
sL2<-mean(subset(countframe, inhe=="shared" & gl=="LOH" & lib=="two")$counts)
uG2<-mean(subset(countframe, inhe=="uglm" & gl=="GOH" & lib=="two")$counts)
uL2<-mean(subset(countframe, inhe=="uglm" & gl=="LOH" & lib=="two")$counts)

nsG1edp<-(subset(countframe, inhe=="not shared" & gl=="GOH" & lib=="one")$counts)
nsL1edp<-(subset(countframe, inhe=="not shared" & gl=="LOH" & lib=="one")$counts)
sG1edp<-(subset(countframe, inhe=="shared" & gl=="GOH" & lib=="one")$counts)
sL1edp<-(subset(countframe, inhe=="shared" & gl=="LOH" & lib=="one")$counts)
uG1edp<-(subset(countframe, inhe=="uglm" & gl=="GOH" & lib=="one")$counts)
uL1edp<-(subset(countframe, inhe=="uglm" & gl=="LOH" & lib=="one")$counts)
nsG2edp<-(subset(countframe, inhe=="not shared" & gl=="GOH" & lib=="two")$counts)
nsL2edp<-(subset(countframe, inhe=="not shared" & gl=="LOH" & lib=="two")$counts)
sG2edp<-(subset(countframe, inhe=="shared" & gl=="GOH" & lib=="two")$counts)
sL2edp<-(subset(countframe, inhe=="shared" & gl=="LOH" & lib=="two")$counts)
uG2edp<-(subset(countframe, inhe=="uglm" & gl=="GOH" & lib=="two")$counts)
uL2edp<-(subset(countframe, inhe=="uglm" & gl=="LOH" & lib=="two")$counts)

nsG1se<-se(subset(countframe, inhe=="not shared" & gl=="GOH" & lib=="one")$counts)
nsL1se<-se(subset(countframe, inhe=="not shared" & gl=="LOH" & lib=="one")$counts)
sG1se<-se(subset(countframe, inhe=="shared" & gl=="GOH" & lib=="one")$counts)
sL1se<-se(subset(countframe, inhe=="shared" & gl=="LOH" & lib=="one")$counts)
uG1se<-se(subset(countframe, inhe=="uglm" & gl=="GOH" & lib=="one")$counts)
uL1se<-se(subset(countframe, inhe=="uglm" & gl=="LOH" & lib=="one")$counts)
nsG2se<-se(subset(countframe, inhe=="not shared" & gl=="GOH" & lib=="two")$counts)
nsL2se<-se(subset(countframe, inhe=="not shared" & gl=="LOH" & lib=="two")$counts)
sG2se<-se(subset(countframe, inhe=="shared" & gl=="GOH" & lib=="two")$counts)
sL2se<-se(subset(countframe, inhe=="shared" & gl=="LOH" & lib=="two")$counts)
uG2se<-se(subset(countframe, inhe=="uglm" & gl=="GOH" & lib=="two")$counts)
uL2se<-se(subset(countframe, inhe=="uglm" & gl=="LOH" & lib=="two")$counts)

means<-c(nsG1, nsL1, sG1, sL1, uG1, uL1, 
         nsG2, nsL2, sG2, sL2, uG2, uL2) 
ses<-c(nsG1se, nsL1se, sG1se, sL1se, uG1se, uL1se, 
       nsG2se, nsL2se, sG2se, sL2se, uG2se, uL2se)  
i<-factor(rep(c("PO, ", "PO, ", "P+S, ", "P+S, ", "SSPO, ","SSPO, "),2),levels = c("PO, ","P+S, ","SSPO, "))
gl<-rep(c("GOH","LOH"),6)
Replicates<-factor(c(rep("without",6),rep("with",6)),levels=c("without","with"))
f<-data.frame(means,ses, i ,gl,Replicates)
f_edp<-data.frame("counts"=c(nsG1edp, nsL1edp, sG1edp, sL1edp, uG1edp, uL1edp, nsG2edp, nsL2edp, sG2edp, sL2edp, uG2edp, uL2edp),
                  "i"=factor(rep(c(rep("PO, ",14),rep("P+S, ",14), rep("SSPO, ",14)),2),levels = c("PO, ","P+S, ","SSPO, ")),
                  "gl"=rep(c(rep("GOH",7),rep("LOH",7)),6),
                  "Replicates"=factor(c(rep("without",42),rep("with",42)),levels=c("without","with"))
)
f$ne<-factor(paste0(i,gl))
f_edp$ne<-factor(paste0(i,gl))

#supp fig 20210625:
ggplot(f, aes(x=ne, y=means, fill=Replicates))+
  geom_bar(position="dodge", stat="identity")+
  theme_bw()+
  ylab("Mutation count") + 
  xlab("")+
  #theme(legend.title = element_text("Libraries per sample"))+
  geom_errorbar(position=position_dodge(width=0.9), aes(ymin=means-ses, ymax=means+ses))+
  scale_fill_manual(values=c("pink","gray"))+
  theme(axis.text.x=element_text(size=25, angle = 90), axis.text.y=element_text(size=25, angle = 0),
        axis.title=element_text(size=25),
        legend.text = element_text(size=25), legend.title = element_text(size=25))+
  scale_y_continuous(trans='log10')
#suppfig2 20210720:
ggplot(f, aes(x=ne, y=means))+
  scale_color_manual(values=c("#D55E00","#56B4E9"))+
  geom_pointrange(data=f, size=1.5, position=position_dodge(width=0.6),aes(color=Replicates, ymin=means-ses, ymax=means+ses))+
  theme_bw()+
  geom_jitter(data=f_edp,size=3,alpha=0.5,position = position_jitterdodge(
    jitter.width = 0.1,
    jitter.height = 0,
    dodge.width = 0.6,
    seed = NA
  ),
  aes(ne, counts,colour=Replicates
  ))+
  ylab("Mutation count") + 
  xlab("")+
  theme(axis.text.x=element_text(size=25, angle = 90), axis.text.y=element_text(size=25, angle = 0),
        axis.title=element_text(size=25),
        legend.text = element_text(size=25), legend.title = element_text(size=25))+
  scale_y_continuous(trans='log10')
  ylab("Mutation count") + 
  xlab("")+
  #theme(legend.title = element_text("Libraries per sample"))+
  geom_errorbar(position=position_dodge(width=0.9), aes(ymin=means-ses, ymax=means+ses))+
  scale_fill_manual(values=c("pink","gray"))+
  theme(axis.text.x=element_text(size=25, angle = 90), axis.text.y=element_text(size=25, angle = 0),
        axis.title=element_text(size=25),
        legend.text = element_text(size=25), legend.title = element_text(size=25))+
  scale_y_continuous(trans='log10')
meanpercentplot2<-ggplot(percentframe, aes(x=percentinhe, y=percentmeans))+
  scale_color_manual(values=c("red","coral1","blue","royalblue1"))+
  geom_pointrange(data=percentframe, size=1.5, position=position_dodge(width=0.6),aes(color=percentgl, ymin=percentmeans-percentses, ymax=percentmeans+percentses))+
  geom_jitter(data=percentframe_eachdatapoint,size=3,position = position_jitterdodge(
    jitter.width = 0.1,
    jitter.height = 0,
    dodge.width = 0.6,
    seed = NA
  ),
  aes(percentinhe, percents,colour=percentgl2
  ))+
  theme_bw()+
  #scale_fill_manual(values = c("red","blue"))+
  ylab("Mutation count with replicates as a percent of mutation count without replicates") + 
  xlab("")+
  theme(legend.title = element_blank())+
  #geom_errorbar(position=position_dodge(width=0.9), aes(ymin=percentmeans-percentses, ymax=percentmeans+percentses))+
  
  theme(axis.text.x=element_text(size=25, angle = 0), axis.text.y=element_text(size=25, angle = 0),
        axis.title=element_text(size=25),legend.text=element_text(size=25))

ggplot(countframe, aes(x=sample, y=counts, group=new, fill=lib))+
    geom_bar(position="dodge", stat="identity")+
    theme_bw()#+
  scale_y_continuous(trans='log10')

  cfpercentplot<-ggplot(subset(countframe,lib=="one"), aes(x=inhe, y=percent, fill=gl))+
    geom_bar(position="dodge", stat="identity")+
    theme_bw()+
    theme(axis.text.x=element_text(size=8, angle = 90), axis.text.y=element_text(size=8, angle = 0),
          axis.title=element_text(size=8))+
    scale_fill_manual(values = c("red","blue"))+
    ylab("Mutation count with replicates as a percent of mutation count without replicates") + 
    xlab("")+
    theme(legend.title = element_blank())+
    facet_wrap(~sample, strip.position = "bottom", scales = "free_x")+
    theme(panel.spacing = unit(0, "lines"), 
          strip.background = element_blank(),
          strip.placement = "outside",
          strip.text.x = element_text(size = 15))
cfpercentplot
meanGOHnotshared<-mean(subset(countframe, lib == "one" & new== "GOHnot shared")$percent)
meanLOHnotshared<-mean(subset(countframe, lib == "one" & new== "LOHnot shared")$percent)
meanGOHshared<-mean(subset(countframe, lib == "one" & new== "GOHshared")$percent)
meanLOHshared<-mean(subset(countframe, lib == "one" & new== "LOHshared")$percent)
meanGOHuglm<-mean(subset(countframe, lib == "one" & new== "GOHuglm")$percent)
meanLOHuglm<-mean(subset(countframe, lib == "one" & new== "LOHuglm")$percent)
meanGOHgglm<-mean(100*c(9,5,5)/c(CA56gglmGOH, CA60gglmGOH, CA65gglmGOH))
meanLOHgglm<-mean(100*c(7,0,16)/c(CA56gglmLOH, CA60gglmLOH, CA65gglmLOH))

GOHnotshared<-(subset(countframe, lib == "one" & new== "GOHnot shared")$percent)
LOHnotshared<-(subset(countframe, lib == "one" & new== "LOHnot shared")$percent)
GOHshared<-(subset(countframe, lib == "one" & new== "GOHshared")$percent)
LOHshared<-(subset(countframe, lib == "one" & new== "LOHshared")$percent)
GOHuglm<-(subset(countframe, lib == "one" & new== "GOHuglm")$percent)
LOHuglm<-(subset(countframe, lib == "one" & new== "LOHuglm")$percent)
GOHgglm<-(100*c(9,5,5)/c(CA56gglmGOH, CA60gglmGOH, CA65gglmGOH))
LOHgglm<-(100*c(7,0,16)/c(CA56gglmLOH, CA60gglmLOH, CA65gglmLOH))

seGOHnotshared<-se(subset(countframe, lib == "one" & new== "GOHnot shared")$percent)
seLOHnotshared<-se(subset(countframe, lib == "one" & new== "LOHnot shared")$percent)
seGOHshared<-se(subset(countframe, lib == "one" & new== "GOHshared")$percent)
seLOHshared<-se(subset(countframe, lib == "one" & new== "LOHshared")$percent)
seGOHuglm<-se(subset(countframe, lib == "one" & new== "GOHuglm")$percent)
seLOHuglm<-se(subset(countframe, lib == "one" & new== "LOHuglm")$percent)
seGOHgglm<-se(c(CA56gglmGOH, CA60gglmGOH, CA65gglmGOH))
seLOHgglm<-se(c(CA56gglmLOH, CA60gglmLOH, CA65gglmLOH))
seGOHgglm<-se(100*c(9,5,5)/c(CA56gglmGOH, CA60gglmGOH, CA65gglmGOH))
seLOHgglm<-se(100*c(7,0,16)/c(CA56gglmLOH, CA60gglmLOH, CA65gglmLOH))

percentmeans<-c(meanGOHnotshared, meanLOHnotshared, meanGOHshared, meanLOHshared, meanGOHuglm, meanLOHuglm, meanGOHgglm, meanLOHgglm)
percentses<-c(seGOHnotshared, seLOHnotshared, seGOHshared, seLOHshared, seGOHuglm, seLOHuglm, seGOHgglm, seLOHgglm)
percentinhe<-factor(c("PO", "PO", "P+S", "P+S", "SSPO","SSPO", "ASP","ASP"),levels = c("PO","P+S","SSPO","ASP"))
percentgl<-c("GoH","LoH","GoH","LoH","GoH","LoH", "GoH","LoH")
percentframe<-data.frame(percentmeans, percentses, percentinhe, percentgl)
"shared"=factor(c(rep("All Parent",14),rep("P+S",14),rep("PO",14),rep("SSPO",14),rep("ASP",6)),
                levels=c("All Parent","PO","P+S","SSPO","ASP")),
"typez"=c(rep(c(rep("GoH2",7),rep("LoH2",7)),4),rep("GoH2",3),rep("LoH2",3)))

percentframe_eachdatapoint<-data.frame("percents"=c(GOHnotshared, LOHnotshared, GOHshared, LOHshared, GOHuglm, LOHuglm, GOHgglm, LOHgglm),
                                       "percentinhe"=factor(c(rep("PO",14),rep("P+S",14),rep("SSPO",14),rep("ASP",6)),levels = c("PO","P+S","SSPO","ASP")),
                                       "percentgl"=c(rep(c(rep("GoH",7),rep("LoH",7)),3),rep("GoH",3),rep("LoH",3)))

meanpercentplot<-ggplot(percentframe, aes(x=percentinhe, y=percentmeans, fill=percentgl))+
  geom_bar(position="dodge", stat="identity")+
  theme_bw()+
  scale_fill_manual(values = c("red","blue"))+
  ylab("Mutation count with replicates as a percent of mutation count without replicates") + 
  xlab("")+
  theme(legend.title = element_blank())+
  geom_errorbar(position=position_dodge(width=0.9), aes(ymin=percentmeans-percentses, ymax=percentmeans+percentses))+

  theme(axis.text.x=element_text(size=25, angle = 0), axis.text.y=element_text(size=25, angle = 0),
        axis.title=element_text(size=25),legend.text=element_text(size=25))
meanpercentplot
#suppfig3 20210720:
meanpercentplot2<-ggplot(percentframe, aes(x=percentinhe, y=percentmeans))+
  scale_color_manual(values=c("red","blue"))+
  geom_pointrange(data=percentframe, size=1.5, position=position_dodge(width=0.6),aes(color=percentgl, ymin=percentmeans-percentses, ymax=percentmeans+percentses))+
  geom_jitter(data=percentframe_eachdatapoint,size=3,alpha=0.5,position = position_jitterdodge(
    jitter.width = 0.1,
    jitter.height = 0,
    dodge.width = 0.6,
    seed = NA
  ),
  aes(percentinhe, percents,colour=percentgl
  ))+
  theme_bw()+
  #scale_fill_manual(values = c("red","blue"))+
  ylab("Mutation count with replicates as a percent of mutation count without replicates") + 
  xlab("")+
  theme(legend.title = element_blank())+
  #geom_errorbar(position=position_dodge(width=0.9), aes(ymin=percentmeans-percentses, ymax=percentmeans+percentses))+
  
  theme(axis.text.x=element_text(size=25, angle = 0), axis.text.y=element_text(size=25, angle = 0),
        axis.title=element_text(size=25),legend.text=element_text(size=25))
