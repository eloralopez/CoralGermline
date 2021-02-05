#dndscv
setwd("~/Documents/GitHub/CoralGermline/dndscv")
library(dndscv)
?buildref()
#buildref("AmilCDStable20200225.txt", "~/Downloads/Amil_v2.01/Amil.v2.01.chrs.fasta", "AmilRefCDS20200226.rda", numcode = 1, excludechrs= NULL)
filepath<-"CAcolony56_allmuts_dndscv.txt"
filepath<-"~/Documents/GitHub/CoralGermline/dndscv/"
dnds_func<-function(filepath) {
  dataset<-read.delim(filepath)
  base<-basename(filepath)
  colony<-strsplit(base, "\\_")[[1]][1]

  dndsout = dndscv(dataset, refdb="AmilRefCDS20200226.rda", cv=NULL)
  #return(dndsout)
  global<-dndsout$globaldnds
  sel_loc<-dndsout$sel_loc
  annot_muts<-dndsout$annotmuts
  codon_sub<- annot_muts$codonsub
  codon_sub<-factor(codon_sub)
  #return(codon_sub)
  dataglobal<-data.frame(global)
  mle_score<-global$mle
  wall<-subset(global,name=="wall")
  wall<-data.frame(wall)
  
  wmis<-subset(global,name=="wmis")
  wmis<-data.frame(wmis)
  
  wnon<-subset(global,name=="wnon")
  wnon<-data.frame(wnon)
  wlist<-list("wall"=wall, "wmis"=wmis, "wnon"=wnon)
  #return(wall)
  return(wlist)
}  

dndsout_somatic_inherited<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/CAcolony56and60and65_inherited_dndscv.txt")
dndsout_somatic_notinherited<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/CAcolony56and60and65_notinherited_dndscv.txt")
dndsout_somatic<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/CAcolony56and60and65_somatic_dndscv.txt")
dndsout_uniqueglm<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/CAcolony56and60and65_uniqueglm_dndscv.txt")
dndsout_globalglm<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/CAcolony56and60and65_globalglm_dndscv.txt")

#USE THESE (5/17/2020)
dndsout_somatic<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/uniquesomaticmetadatadf_dndscv.txt")
dndsout_somatic_denovo<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/uniquesomaticmetadatadf_denovo_dndscv.txt")
dndsout_somatic_loh<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/uniquesomaticmetadatadf_loh_dndscv.txt")
dndsout_somatic_inherited<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/inheritedallcolonies_dndscv.txt")
dndsout_somatic_notinherited<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/notinheritedallcolonies_dndscv.txt")

dndsout_somatic_inherited_loh<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/inheritedlohallcolonies_dndscv.txt")
dndsout_somatic_notinherited_loh<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/notinheritedlohallcolonies_dndscv.txt")

dndsout_somatic_inherited_denovo<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/inheriteddenovoallcolonies_dndscv.txt")
dndsout_somatic_notinherited_denovo<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/notinheriteddenovoallcolonies_dndscv.txt")

dndsout_uglm<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/unique_uglmmetadatadf_dndscv.txt")
dndsout_gglm<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/unique_gglmmetadatadf_dndscv.txt")

#FOR 6/15/2020
dndsout_somatic<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/allcolonies_onlychrs20200615_dndscv.txt")
dndsout_somatic_56<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/CAcolony56somatic_onlychrs20200616_dndscv.txt")
dndsout_somatic_60<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/CAcolony60somatic_onlychrs20200616_dndscv.txt")
dndsout_somatic_65<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/CAcolony65somatic_onlychrs20200616_dndscv.txt")

dndsout_somatic_denovo<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/uniquesomaticmetadatadf_denovo_onlychrs20200615_dndscv.txt")
dndsout_somatic_loh<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/uniquesomaticmetadatadf_loh_onlychrs20200615_dndscv.txt")
dndsout_somatic_inherited<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/inheritedallcolonies_onlychrs20200615_dndscv.txt")
dndsout_somatic_notinherited<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/notinheritedallcolonies_onlychrs20200615_dndscv.txt")

dndsout_somatic_inherited_loh<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/inheritedlohallcolonies_onlychrs20200615_dndscv.txt")
dndsout_somatic_notinherited_loh<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/notinheritedlohallcolonies_onlychrs20200615_dndscv.txt")

dndsout_somatic_inherited_denovo<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/inheriteddenovoallcolonies_onlychrs20200615_dndscv.txt")
dndsout_somatic_notinherited_denovo<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/notinheriteddenovoallcolonies_onlychrs20200615_dndscv.txt")

dndsout_uglm<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/unique_uglmmetadatadf_onlychrs20200615_dndscv.txt")
dndsout_gglm<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/unique_gglmmetadatadf_onlychrs20200615_dndscv.txt")

#FOR 8/11/2020:
dndsout08_somatic<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcolonies_20200811_dndscv.txt")
dndsout08_somatic_GOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesGOH_20200811_dndscv.txt")
dndsout08_somatic_LOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesLOH_20200811_dndscv.txt")
dndsout08_somatic_inherited<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesinherited_20200811_dndscv.txt")
dndsout08_somatic_inherited_GOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesinherited_GOH_20200811_dndscv.txt")
dndsout08_somatic_inherited_LOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesinherited_LOH_20200811_dndscv.txt")

dndsout08_somatic_notinherited<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesnotinherited_20200811_dndscv.txt")
dndsout08_somatic_notinherited_GOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesnotinherited_GOH_20200811_dndscv.txt")
dndsout08_somatic_notinherited_LOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesnotinherited_LOH_20200811_dndscv.txt")

dndsout08_uglm<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/uglmallcolonies_20200811_dndscv.txt")
dndsout08_uglm_GOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/uglmallcolonies_GOH_20200811_dndscv.txt")
dndsout08_uglm_LOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/uglmallcolonies_LOH_20200811_dndscv.txt")
dndsout08_gglm<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/gglmallcolonies_20200811_dndscv.txt")

#FOR 10/19/2020:
dndsout_somatic<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcolonies_20201019_dndscv.txt")
dndsout_somatic_GOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesGOH_20201019_dndscv.txt")
dndsout_somatic_LOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesLOH_20201019_dndscv.txt")
dndsout_somatic_inherited<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesinherited_20201019_dndscv.txt")
dndsout_somatic_inherited_GOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesinherited_GOH_20201019_dndscv.txt")
dndsout_somatic_inherited_LOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesinherited_LOH_20201019_dndscv.txt")

dndsout_somatic_notinherited<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesnotinherited_20201019_dndscv.txt")
dndsout_somatic_notinherited_GOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesnotinherited_GOH_20201019_dndscv.txt")
dndsout_somatic_notinherited_LOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesnotinherited_LOH_20201019_dndscv.txt")

dndsout_uglm<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/uglmallcolonies_20201019_dndscv.txt")
dndsout_uglm_GOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/uglmallcolonies_GOH_20201019_dndscv.txt")
dndsout_uglm_LOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/uglmallcolonies_LOH_20201019_dndscv.txt")
dndsout_gglm<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/gglmallcolonies_20201019_dndscv.txt")

#FOR 11/17/2020:
dndsout_somatic<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcolonies_20201117_dndscv.txt")
dndsout_somatic_GOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesGOH_20201117_dndscv.txt")
dndsout_somatic_LOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesLOH_20201117_dndscv.txt")
dndsout_somatic_inherited<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesinherited_20201117_dndscv.txt")
dndsout_somatic_inherited_GOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesinherited_GOH_20201117_dndscv.txt")
dndsout_somatic_inherited_LOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesinherited_LOH_20201117_dndscv.txt")

dndsout_somatic_notinherited<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesnotinherited_20201117_dndscv.txt")
dndsout_somatic_notinherited_GOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesnotinherited_GOH_20201117_dndscv.txt")
dndsout_somatic_notinherited_LOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesnotinherited_LOH_20201117_dndscv.txt")

dndsout_uglm<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/uglmallcolonies_20201117_dndscv.txt")
dndsout_uglm_GOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/uglmallcolonies_GOH_20201117_dndscv.txt")
dndsout_uglm_LOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/uglmallcolonies_LOH_20201117_dndscv.txt")
dndsout_gglm<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/gglmallcolonies_20201117_dndscv.txt")

#12/3
dndsout_somatic<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcolonies_20201203_dndscv.txt")
dndsout_somatic_GOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesGOH_20201203_dndscv.txt")
dndsout_somatic_LOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesLOH_20201203_dndscv.txt")
dndsout_somatic_inherited<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesinherited_20201203_dndscv.txt")
dndsout_somatic_inherited_GOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesinherited_GOH_20201203_dndscv.txt")
dndsout_somatic_inherited_LOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesinherited_LOH_20201203_dndscv.txt")

dndsout_somatic_notinherited<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesnotinherited_20201203_dndscv.txt")
dndsout_somatic_notinherited_GOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesnotinherited_GOH_20201203_dndscv.txt")
dndsout_somatic_notinherited_LOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesnotinherited_LOH_20201203_dndscv.txt")

dndsout_uglm<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/uglmallcolonies_20201203_dndscv.txt")
dndsout_uglm_GOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/uglmallcolonies_GOH_20201203_dndscv.txt")
dndsout_uglm_LOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/uglmallcolonies_LOH_20201203_dndscv.txt")
dndsout_gglm<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/gglmallcolonies_20201203_dndscv.txt")

dndsout_GOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/GOH_all_20201203_dndscv.txt")
dndsout_LOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/LOH_all_20201203_dndscv.txt")

#12/10
dndsout_somatic<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcolonies_20201210_dndscv.txt")
dndsout_somatic_GOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesGOH_20201210_dndscv.txt")
dndsout_somatic_LOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesLOH_20201210_dndscv.txt")
dndsout_somatic_inherited<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesinherited_20201210_dndscv.txt")
dndsout_somatic_inherited_GOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesinherited_GOH_20201210_dndscv.txt")
dndsout_somatic_inherited_LOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesinherited_LOH_20201210_dndscv.txt")

dndsout_somatic_notinherited<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesnotinherited_20201210_dndscv.txt")
dndsout_somatic_notinherited_GOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesnotinherited_GOH_20201210_dndscv.txt")
dndsout_somatic_notinherited_LOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesnotinherited_LOH_20201210_dndscv.txt")

dndsout_uglm<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/uglmallcolonies_20201210_dndscv.txt")
dndsout_uglm_GOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/uglmallcolonies_GOH_20201210_dndscv.txt")
dndsout_uglm_LOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/uglmallcolonies_LOH_20201210_dndscv.txt")
dndsout_gglm<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/gglmallcolonies_20201210_dndscv.txt")

dndsout_GOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/GOH_all_20201210_dndscv.txt")
dndsout_LOH<-dnds_func("~/Documents/GitHub/CoralGermline/dndscv/LOH_all_20201210_dndscv.txt")

wall<-rbind(dndsout_somatic_inherited$wall, dndsout_somatic_notinherited$wall, dndsout_uglm$wall)
wmis<-rbind(dndsout_somatic_inherited$wmis, dndsout_somatic_notinherited$wmis, dndsout_uglm$wmis)
wnon<-rbind(dndsout_somatic_inherited$wnon, dndsout_somatic_notinherited$wnon, dndsout_uglm$wnon)

wall<-rbind(dndsout_uglm_GOH$wall, dndsout_uglm_LOH$wall, dndsout_somatic_inherited_GOH$wall, 
            dndsout_somatic_notinherited_GOH$wall, dndsout_somatic_notinherited_LOH$wall)#$wall, dndsout_somatic_notinherited$wall, dndsout_somatic_notinherited_loh$wall, dndsout_somatic_notinherited_denovo)
cat<-c("SSPO, GOH","SSPO, LOH","P+S, GOH",
       "PO, GOH","PO, LOH")
wall<-rbind(dndsout_somatic_GOH$wall, dndsout_somatic_inherited_GOH$wall, dndsout_somatic_notinherited_GOH$wall, dndsout_uglm_GOH$wall)
wmis<-rbind(dndsout_somatic$wmis, dndsout_uglm$wmis, dndsout_uglm_GOH$wmis, dndsout_uglm_LOH$wmis, dndsout_somatic_inherited_GOH$wmis, 
            dndsout_somatic_notinherited_GOH$wmis, dndsout_somatic_notinherited_LOH$wmis)
wnon<-rbind(dndsout_somatic$wnon, dndsout_uglm$wnon, dndsout_uglm_GOH$wnon, dndsout_uglm_LOH$wnon, dndsout_somatic_inherited_GOH$wnon, 
            dndsout_somatic_notinherited_GOH$wnon, dndsout_somatic_notinherited_LOH$wnon)
colony<-c("Somatic", "Not Inherited Somatic", "Not inherited LoH","Not inherited denovo")#, "CA60","CA60", "CA65", "CA65")
colony<-c("all","denovo","loh", "inherited","notinherited", "uglm", "not inherited_loh","not inherited_denovo")
colony<-c("All Somatic","Unique GLM", "uglm GOH", "uglm LOH" , "Inh GOH","NotINh GOH","NotInh LOH")
wall<-rbind(dndsout_somatic_inherited, dndsout_somatic_notinherited, dndsout_uglm)
colony<-c("All","Inherited Somatic", "Not Inherited Somatic","Germline")
mle<-wall$mle
cilow<-wall$cilow
cihigh<-wall$cihigh

wbind<-rbind(wall,wmis,wnon)
mle<-wbind$mle
cilow<-wbind$cilow
cihigh<-wbind$cihigh
w_typez<-wbind$name
#category<-rep(c("All Somatic","Unique GLM", "uglm GOH", "uglm LOH" , "Inh GOH","NotINh GOH","NotInh LOH"),3)
category<-rep(c("Inherited","Not inherited","Unique GLM"),3)
category_wall<-c("Inherited","Not inherited","Unique GLM")
mle_wall<-wall$mle
barplot(mle)
boxplot(mle~colony, xlab="Colony Name", ylab="global dN/dS", las=1, col=c("light blue","light gray","light green"))
dndscvplot_wbind<-ggplot(wbind, 
       aes(x = category, 
           y = mle, 
           group = w_typez, colour=w_typez)) +
  geom_point(size = 3, position=position_dodge(width=0.3)) +
  
  geom_errorbar(aes(ymin = cilow, 
                    ymax = cihigh), 
                width = .1, position=position_dodge(width=0.3)) +
  theme_bw() +coord_cartesian(ylim=c(0,8))+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=25),
        axis.text.x = element_text(angle = 90))+
  labs(y="dN/dS", x="") + 
  geom_hline(yintercept=1, linetype="dashed", 
             color = "red", size=1)
#+
dndscvplot_wall<-ggplot(wall, 
                   aes(x = cat, 
                       y = mle, 
                       group = 1)) +
  geom_point(size = 3, position=position_dodge(width=0.3)) +
  
  geom_errorbar(aes(ymin = cilow, 
                    ymax = cihigh), 
                width = 1, position=position_dodge(width=0.3)) +
  theme_bw() +coord_cartesian(ylim=c(0,9))+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=25),
        axis.text.x = element_text(angle = 90))+
  labs(y="dN/dS", x="") + 
  geom_hline(yintercept=1, linetype="dashed", 
             color = "red", size=1)
  
sel_cv = dndsout_somatic$sel_cv
print(head(sel_cv), digits = 3)  
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

dndsout10 = dndscv(read.delim("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesGOH_20201019_dndscv.txt"), refdb="AmilRefCDS20200226.rda", cv=NULL)
dndsout08 = dndscv(read.delim("~/Documents/GitHub/CoralGermline/dndscv/somaticallcoloniesGOH_20200811_dndscv.txt"), refdb="AmilRefCDS20200226.rda", cv=NULL)
