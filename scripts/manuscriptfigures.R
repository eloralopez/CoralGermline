#mean freqs per sample, grouped by colony
#not inherited somatic:
freq_func<-function(CA56n, CA56denom, CA60n, CA60denom, CA65n, CA65denom, CA56i, CA60i, CA65i, CA56uglm,CA60uglm,CA65uglm) {
  CA56nfreq<-CA56n/CA56denom
  CA56nfreqse<-se(CA56nfreq)
  CA60nfreq<-CA60n/CA60denom
  CA60nfreqse<-se(CA60nfreq)
  CA65nfreq<-CA65n/CA65denom
  CA65nfreqse<-se(CA65nfreq)
  
  nfreqs<-c(CA56nfreq,CA60nfreq, CA65nfreq)
  nfreqses<- c(CA56nfreqse, CA60nfreqse, CA65nfreqse)
  nfreqmeans<-c(mean(CA56nfreq), mean(CA60nfreq), mean(CA65nfreq))
  
  #inherited somatic:
  CA56ifreq<-CA56i/CA56denom
  CA56ifreqse<-se(CA56ifreq)
  CA60ifreq<-CA60i/CA60denom
  CA60ifreqse<-se(CA60ifreq)
  CA65ifreq<-CA65i/CA65denom
  CA65ifreqse<-se(CA65ifreq)
  
  ifreqs<-c(CA56ifreq,CA60ifreq, CA65ifreq)
  ifreqses<- c(CA56ifreqse, CA60ifreqse, CA65ifreqse)
  ifreqmeans<-c(mean(CA56ifreq), mean(CA60ifreq), mean(CA65ifreq))
  
  CA56total<-CA56ifreq + CA56nfreq
  CA60total<-CA60ifreq + CA60nfreq
  CA65total<-CA65ifreq + CA65nfreq
  
  CA56totalse<-se(CA56total)
  CA60totalse<-se(CA60total)
  CA65totalse<-se(CA65total)
  
  CA56mean<-mean(CA56total)
  CA60mean<-mean(CA60total)
  CA65mean<-mean(CA65total)
  
  totalfreqs<-nfreqs+ifreqs
  totalfreqses<-se(totalfreqs)
  totalfreqsmean<-mean(totalfreqs)
  
  #return(list(CA56mean, CA60mean,CA65mean, CA56totalse, CA60totalse, CA65totalse))
  
  #uglm:
  CA56uglmfreq<-CA56uglm/CA56denom
  CA56uglmfreqse<-se(CA56uglmfreq)
  CA60uglmfreq<-CA60uglm/CA60denom
  CA60uglmfreqse<-se(CA60uglmfreq)
  CA65uglmfreq<-CA65uglm/CA65denom
  CA65uglmfreqse<-se(CA65uglmfreq)
  
  uglmfreqs<-c(CA56uglmfreq,CA60uglmfreq, CA65uglmfreq)
  uglmfreqses<- c(CA56uglmfreqse, CA60uglmfreqse, CA65uglmfreqse)
  uglmfreqmeans<-c(mean(CA56uglmfreq), mean(CA60uglmfreq), mean(CA65uglmfreq))
  
  meanfreqpersample<-c(nfreqmeans, ifreqmeans, uglmfreqmeans)
  meanfreqpersamplese<-c(nfreqses, ifreqses,uglmfreqses)
  types<-c(rep("notinherited",3), rep("inherited",3),rep("uglm",3))
  colony<-rep(c("CA56","CA60","CA65"),3)
  area<-
  df<-data.frame("meanpersamplefreqs"=meanfreqpersample, "se"=meanfreqpersamplese, "type"=types,"colony"=colony)
  meanfreqpersampleplot<-ggplot(df, 
                                  aes(x = type, 
                                      y = meanpersamplefreqs, 
                                      group = colony, colour=colony)) +
    geom_point(size = 3, position=position_dodge(width=0.3)) +
    #geom_point(position="dodge", stat="identity")
    geom_errorbar(aes(ymin=meanpersamplefreqs-(2*se), ymax=meanpersamplefreqs+(2*se)),
                  width=0.1, position=position_dodge(width=0.3)) +
    theme_bw() +
    theme(axis.text=element_text(size=15),
          axis.title=element_text(size=15),
          axis.text.x = element_text(angle = 90))+
    labs(y="# mutations bp per sample", x="") ##Figure 1
  persamplemeans<-c(mean(nfreqs), mean(ifreqs),mean(uglmfreqs))
  persampleses<-c(se(nfreqs), se(ifreqs),se(uglmfreqs))
  persampletypes<-c("notinherited","inherited","uglm")
  persampledf<-data.frame(persamplemeans,persampleses,persampletypes)
  persampledf_notmeans<-data.frame("persample"=c(nfreqs,ifreqs,uglmfreqs),"persampletypes"=c(rep("notinherited",7),rep("inherited",7),rep("uglm",7)))
  return(list("means"=persampledf,"eachdatapoint"=persampledf_notmeans))
}
meanfreqpersampleplot<-freq_func(CA56n, CA56denom, CA60n, CA60denom, CA65n, CA65denom, CA56i, CA60i, CA65i, CA56uglm,CA60uglm,CA65uglm)
meanfreqpersampleplot$Type<-rep("LOH+GOH",3)

##Inherited vs Not inherited per colony
class<-factor(rep(c("Not Inherited" , "Inherited") , 3),levels = c("Not Inherited","Inherited"))


nfreqspermillion<-nfreqs*1000000
ifreqspermillion<-ifreqs*1000000
freqspermillion<-c(nfreqspermillion, ifreqspermillion)
n_percent<-nfreqs*100/(nfreqs+ifreqs)
i_percent<-ifreqs*100/(nfreqs+ifreqs)
mean(n_percent)
se(n_percent)
mean(i_percent)
se(i_percent)
n_perc<-(100*CAPn)/(CAPn+CAPi)
i_perc<-(100*CAPi)/(CAPn+CAPi)
##We identified post-embryonic SNVs in every parent branch, and we found that on average 35.6% ± 4.1% (1 s.e.m.) SNVs were shared between a parent branch and its respective sperm pool, but significantly more, 64.4% ± 4.1% (1 s.e.m.) post-embryonic SNVs found in a branch were not in the sperm (two-tailed, unpaired t-test, t12 = 4.9092, x, P = 0.00036).##
t.test(n_percent, i_percent,alternative=c("two.sided"),paired=FALSE)
colonies<-rep(c(rep("Colony 56",2),rep("Colony 60",3),rep("Colony 65",2)),2)
inh<-factor(c(rep("Not Inherited",7),rep("Inherited",7)),levels = c("Not Inherited","Inherited"))

samples<-factor(c("CAP6","CAP8","CAP22","CAP23","CAP26","CAP9","CAP11"),
                levels = c("CAP6","CAP8","CAP22","CAP23","CAP26","CAP9","CAP11"))
compdf<-data.frame(freqspermillion, colonies, inh,samples)
compdf$shared<-c(rep("PO",length(nfreqspermillion)),rep("P+S",length(ifreqspermillion)))
compdf$counts<-c(CA56n, CA60n, CA65n,CA56i, CA60i, CA65i)
freqspmNI<-subset(compdf, inh=="Not Inherited")
freqspmI<-subset(compdf, inh=="Inherited")
freqspmNImean<-mean(freqspmNI$freqspermillion)
freqspmImean<-mean(freqspmI$freqspermillion)
proportions<-ifreqspermillion/(nfreqspermillion+ifreqspermillion)

mean(proportions) #average  inherited
se(proportions) #SE inherited
1-proportions
mean(1-proportions)
se(1-proportions)
comparisonplotinh2<-ggplot(compdf, aes(x=samples, y=freqspermillion/1000000, color=shared)) + 
 geom_point(position="identity", size=5) +#, stat="identity",position_dodge(width = 1)) +
  theme_bw() +
  scale_color_manual(values = c("#999999","goldenrod")) +
  ylab("# SNVs per bp") + xlab("")+
  theme(axis.text=element_text(size=25), 
    axis.text.x=element_text(size=25, angle = 90),
        axis.title=element_text(size=25))+
  theme(legend.text=element_text(size=25), legend.title = element_blank())+
  facet_wrap(~colonies, strip.position = "top", scales = "free_x")+
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside",
        strip.text.x = element_text(size = 25))
cpi2<-comparisonplotinh2+ scale_y_continuous(expand = c(0,0.00000002))+theme(panel.grid.major = element_blank(),
                            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
cpi2
colony<-c(rep("CA56",2),rep("CA60",3),rep("CA65",2))

df_inheritedprops<-data.frame(colony,proportions)
box <- ggplot(df_inheritedprops, aes(x = colony, y = proportions))
box<- box + geom_boxplot()
box2 <- box + labs(x="", y="Percent of Somatic Mutations inherited by Sperm") + geom_boxplot(fill="red",
                                                                                             position = position_dodge(0.9)
) +
  
  theme_bw()+
  #stat_compare_means(label.x=1,label.y=30,size=10)+
  theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25))
model  <- lm(proportions ~ colony, data = df_inheritedprops)
ggqqplot(residuals(model))
shapiro.test(residuals(model))
ggqqplot(df_inheritedprops, "proportions", facet.by = "colony")
##patchwork workaround: 
meanfreqpersampleplot / comparisonplotinh2 & theme(strip.placement = NULL)

CA56nfreqspermillionsum<-nfreqspermillion[1] +nfreqspermillion[2]
                                                               
CA60nfreqspermillionsum<-nfreqspermillion[3] +nfreqspermillion[4] + nfreqspermillion[5]

CA65nfreqspermillionsum<-nfreqspermillion[6] +nfreqspermillion[7]

CA56ifreqspermillionsum<-ifreqspermillion[1] +ifreqspermillion[2]

CA60ifreqspermillionsum<-ifreqspermillion[3] +ifreqspermillion[4] + ifreqspermillion[5]

CA65ifreqspermillionsum<-ifreqspermillion[6] +ifreqspermillion[7]

CA56uglmpermillionsum <-sum(CA56uglmfreq) *1000000
CA60uglmpermillionsum <-sum(CA60uglmfreq)*1000000
CA65uglmpermillionsum <-sum(CA65uglmfreq)*1000000

CA56gglmpermillion <- (sum(CA56gglm)/CA56denom) *1000000
CA60gglmpermillion <- (sum(CA60gglm)/CA60denom) *1000000
CA65gglmpermillion <- (sum(CA65gglm)/CA65denom) *1000000

colony<-c(rep("CA56" , 3) , rep("CA60" , 3) , rep("CA65" , 3))
class<-factor(rep(c("Inherited Somatic", "Unique Germline", "Global Germline") , 3),levels = c("Inherited Somatic","Unique Germline","Global Germline"))
value<-c(CA56ifreqspermillionsum, CA56uglmpermillionsum, CA56gglmpermillion, 
         CA60ifreqspermillionsum, CA60uglmpermillionsum, CA60gglmpermillion,
         CA65ifreqspermillionsum, CA65uglmpermillionsum, CA65gglmpermillion)
gglmfreqs<-c(rep(CA56gglm[1]/CA56denom,2),rep(CA60gglm[1]/CA60denom,3),rep(CA65gglm[1]/CA65denom,2))
valuemean<-c(mean(ifreqs),mean(uglmfreqs),mean(c(CA56gglm[1]/CA56denom,CA60gglm[1]/CA60denom,CA65gglm[1]/CA65denom)))*1000000
valuemeanse<-c(se(ifreqs),se(uglmfreqs),se(c(CA56gglm[1]/CA56denom,CA60gglm[1]/CA60denom,CA65gglm[1]/CA65denom)))*1000000
meandata<-data.frame(valuemean,valuemeanse,"shared"=factor(c("Parent and Sperm", "Single Sperm Pool Only", "All Sperm Pools"),
                                                           levels = c("Parent and Sperm","Single Sperm Pool Only","All Sperm Pools")))
ifreqpercents<-ifreqs*100/(ifreqs+uglmfreqs+gglmfreqs)
uglmfreqpercents<-uglmfreqs*100/(ifreqs+uglmfreqs+gglmfreqs)
gglmfreqpercents<-gglmfreqs*100/(ifreqs+uglmfreqs+gglmfreqs)
valuemeanpercent<-c(mean(ifreqpercents),mean(uglmfreqpercents),mean(gglmfreqpercents))
valuemeanpercentse<-c(se(ifreqpercents),se(uglmfreqpercents),se(gglmfreqpercents))
meanpercentdata<-data.frame(valuemeanpercent,valuemeanpercentse,"shared"=factor(c("Parent and Sperm", "Single Sperm Pool Only", "All Sperm Pools"),
                                                                  levels = c("Parent and Sperm","Single Sperm Pool Only","All Sperm Pools")))
percentdata<-data.frame("percents"=c(ifreqpercents, uglmfreqpercents, gglmfreqpercents),"shared"=factor(c(rep("Parent and Sperm",7), rep("Single Sperm Pool Only",7), rep("All Sperm Pools",7)),
                                                                                             levels = c("Parent and Sperm","Single Sperm Pool Only","All Sperm Pools")),
                        "colz"=c(rep("1",7),rep("2",7),rep("3",7)))
percentdataplot<-ggplot(percentdata, aes(shared, percents)) +
  geom_jitter(size=3,aes(color=colz,group=shared),
    position = position_jitter(0.2)
  )
#use 20210719
  meanpercentdataplot2<-ggplot(meanpercentdata, aes(shared,valuemeanpercent))+
    geom_pointrange(size=1.5,data = meanpercentdata,aes(shared,valuemeanpercent,group=shared, color=shared, ymin=valuemeanpercent-valuemeanpercentse, ymax=valuemeanpercent+valuemeanpercentse))+
    geom_jitter(data=percentdata,size=3,alpha=0.5,aes(shared, percents, color=shared,group=shared),
                position = position_jitter(0.2)
    )+
  scale_color_manual(values = c("#999999","purple","seagreen3")) +
    theme_bw()+
  ylab("% of SNV type per sperm pool") + xlab("")+
    theme(axis.text.x = element_text(size=25,colour=c("#999999","purple","seagreen3")),
          axis.title=element_text(size=25),
          axis.text.y = element_text(size=25))+
    theme(legend.position = "none") +
    scale_x_discrete(labels=c("Parent and Sperm" = "P+S", "Single Sperm Pool Only" = "SSPO",
                              "All Sperm Pools" = "ASP"))
  #geom_errorbar(aes(ymin=valuemeanpercent-valuemeanpercentse, ymax=valuemeanpercent+valuemeanpercentse))
mutationdata<-data.frame(colony, class, value)
mutationdata$shared<-factor(rep(c("Parent and Sperm", "Single Sperm Pool Only", "All Sperm Pools") , 3),levels = c("Parent and Sperm","Single Sperm Pool Only","All Sperm Pools"))

inheritedpermillionplot<-ggplot(mutationdata, aes(x=colony, y=value, fill=shared)) + 
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  scale_fill_manual(values = c("black", "purple", "seagreen3")) +
  ylab("# SNVs per Mbp") + xlab("")+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15))+
  theme(legend.text=element_text(size=15), legend.title = element_blank()) 

meanpercentdataplot<-ggplot(meanpercentdata, aes(x=shared, y=valuemeanpercent, fill=shared)) + 
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  scale_fill_manual(values = c("#999999", "purple", "seagreen3")) +
  ylab("% of SNV type per sperm pool") + xlab("")+
  geom_pointrange(aes(ymin=valuemeanpercent-valuemeanpercentse, ymax=valuemeanpercent+valuemeanpercentse))+
  theme(axis.text.x = element_text(size=25,colour=c("#999999","purple","seagreen3")),
        axis.title=element_text(size=25),
        axis.text.y = element_text(size=25))+
  theme(legend.position = "none") +
  scale_x_discrete(labels=c("Parent and Sperm" = "P+S", "Single Sperm Pool Only" = "SSPO",
                            "All Sperm Pools" = "ASP"))
  
#find only coding muts:
coding_func<-function(df) {
  df_sub_missense<-subset(df, MutationType=="missense_variant")# | MutationType=="synonymous_variant" | MutationType=="nonsense_variant")
  df_sub_syn<-subset(df, MutationType=="synonymous_variant")
  df_sub_non<-subset(df, MutationType=="stop_gained")
  allcoding<-rbind(df_sub_missense, df_sub_non, df_sub_syn)
  return(allcoding)
}

#notinheritedcodingfreq<-c(nrow(coding_func(CAP6n))/CA56denom_coding, nrow(coding_func(CAP8n))/CA56denom_coding,
                     # nrow(coding_func(CAP22n))/CA60denom_coding, nrow(coding_func(CAP23n))/CA60denom_coding, nrow(coding_func(CAP24n))/CA60denom_coding,
                      #nrow(coding_func(CAP11n))/CA65denom_coding, nrow(coding_func(CAP9n))/CA65denom_coding)

#inheritedcodingfreq<-c(nrow(coding_func(CAP6i))/CA56denom_coding, nrow(coding_func(CAP8i))/CA56denom_coding,
                     # nrow(coding_func(CAP22i))/CA60denom_coding, nrow(coding_func(CAP23i))/CA60denom_coding, nrow(coding_func(CAP24i))/CA60denom_coding,
                     # nrow(coding_func(CAP11i))/CA65denom_coding, nrow(coding_func(CAP9i))/CA65denom_coding)


#uglmcodingfreq<-c(nrow(coding_func(CAS6uglm))/CA56denom_coding, nrow(coding_func(CAS8uglm))/CA56denom_coding,
              #  nrow(coding_func(CAS22uglm))/CA60denom_coding, nrow(coding_func(CAS23uglm))/CA60denom_coding, nrow(coding_func(CAS24uglm))/CA60denom_coding,
               # nrow(coding_func(CAS11uglm))/CA65denom_coding, nrow(coding_func(CAS9uglm))/CA65denom_coding)

CA56n_coding<-c(nrow(coding_func(CAP6n)), nrow(coding_func(CAP8n)))
CA56i_coding<-c(nrow(coding_func(CAP6i)), nrow(coding_func(CAP8i)))
CA56uglm_coding<-c(nrow(coding_func(CAS6uglm)), nrow(coding_func(CAS8uglm)))

CA60n_coding<-c(nrow(coding_func(CAP22n)), nrow(coding_func(CAP23n)), nrow(coding_func(CAP26n)))
CA60i_coding<-c(nrow(coding_func(CAP22i)), nrow(coding_func(CAP23i)), nrow(coding_func(CAP26i)))
CA60uglm_coding<-c(nrow(coding_func(CAS22uglm)), nrow(coding_func(CAS23uglm)), nrow(coding_func(CAS26uglm)))

CA65n_coding<-c(nrow(coding_func(CAP11n)), nrow(coding_func(CAP9n)))
CA65i_coding<-c(nrow(coding_func(CAP11i)), nrow(coding_func(CAP9i)))
CA65uglm_coding<-c(nrow(coding_func(CAS11uglm)), nrow(coding_func(CAS9uglm)))

##lists of coding mutations:
all_inherited_coding<-rbind(coding_func(CAP6i), coding_func(CAP8i), coding_func(CAP22i), coding_func(CAP23i), coding_func(CAP26i), coding_func(CAP11i), coding_func(CAP9i))
all_noninherited_coding<-rbind(coding_func(CAP6n), coding_func(CAP8n), coding_func(CAP22n), coding_func(CAP23n), coding_func(CAP26n), coding_func(CAP11n), coding_func(CAP9n))
all_uglm_coding<-rbind(coding_func(CAS6uglm), coding_func(CAS8uglm), coding_func(CAS22uglm), coding_func(CAS23uglm), coding_func(CAS26uglm), coding_func(CAS11uglm), coding_func(CAS9uglm))  
#write.table(all_inherited_coding, file="~/Documents/GitHub/CoralGermline/suppfiles/all_inherited_coding.txt",sep="\t",quote=FALSE, row.name=FALSE)  
#write.table(all_noninherited_coding, file="~/Documents/GitHub/CoralGermline/suppfiles/all_noninherited_coding.txt",sep="\t",quote=FALSE, row.name=FALSE)  
#write.table(all_uglm_coding, file="~/Documents/GitHub/CoralGermline/suppfiles/all_uglm_coding.txt",sep="\t",quote=FALSE, row.name=FALSE)  

meanfreqpersampleplot_coding<-freq_func(CA56n_coding, CA56denom_coding, CA60n_coding, CA60denom_coding, CA65n_coding, CA65denom_coding, 
          CA56i_coding, CA60i_coding, CA65i_coding, CA56uglm_coding, CA60uglm_coding,CA65uglm_coding)

CA56n_coding_GOH<-c(nrow(coding_func(subset(CAP6n, GoHorLoH == "DeNovo"))), nrow(coding_func(subset(CAP8n, GoHorLoH == "DeNovo"))))
CA56i_coding_GOH<-c(nrow(coding_func(subset(CAP6i, GoHorLoH == "DeNovo"))), nrow(coding_func(subset(CAP8i, GoHorLoH == "DeNovo"))))
CA56uglm_coding_GOH<-c(nrow(coding_func(subset(CAS6uglm, GoH_or_LoH == "DeNovo"))), nrow(coding_func(subset(CAS8uglm, GoH_or_LoH == "DeNovo"))))

CA60n_coding_GOH<-c(nrow(coding_func(subset(CAP22n, GoHorLoH == "DeNovo"))), nrow(coding_func(subset(CAP23n, GoHorLoH == "DeNovo"))), nrow(coding_func(subset(CAP26n, GoHorLoH == "DeNovo"))))
CA60i_coding_GOH<-c(nrow(coding_func(subset(CAP22i, GoHorLoH == "DeNovo"))), nrow(coding_func(subset(CAP23i, GoHorLoH == "DeNovo"))), nrow(coding_func(subset(CAP26i, GoHorLoH == "DeNovo"))))
CA60uglm_coding_GOH<-c(nrow(coding_func(subset(CAS22uglm, GoH_or_LoH == "DeNovo"))), nrow(coding_func(subset(CAS23uglm, GoH_or_LoH == "DeNovo"))), nrow(coding_func(subset(CAS26uglm, GoH_or_LoH == "DeNovo"))))

CA65n_coding_GOH<-c(nrow(coding_func(subset(CAP11n, GoHorLoH == "DeNovo"))), nrow(coding_func(subset(CAP9n, GoHorLoH == "DeNovo"))))
CA65i_coding_GOH<-c(nrow(coding_func(subset(CAP11i, GoHorLoH == "DeNovo"))), nrow(coding_func(subset(CAP9i, GoHorLoH == "DeNovo"))))
CA65uglm_coding_GOH<-c(nrow(coding_func(subset(CAS11uglm, GoH_or_LoH == "DeNovo"))), nrow(coding_func(subset(CAS9uglm, GoH_or_LoH == "DeNovo"))))
meanfreqpersampleplot_coding_GOH<-freq_func(CA56n_coding_GOH, CA56denom_coding, CA60n_coding_GOH, CA60denom_coding, CA65n_coding_GOH, CA65denom_coding, 
                                        CA56i_coding_GOH, CA60i_coding_GOH, CA65i_coding_GOH, CA56uglm_coding_GOH, CA60uglm_coding_GOH,CA65uglm_coding_GOH)

CA56n_coding_LOH<-c(nrow(coding_func(subset(CAP6n, GoHorLoH == "LoH"))), nrow(coding_func(subset(CAP8n, GoHorLoH == "LoH"))))
CA56i_coding_LOH<-c(nrow(coding_func(subset(CAP6i, GoHorLoH == "LoH"))), nrow(coding_func(subset(CAP8i, GoHorLoH == "LoH"))))
CA56uglm_coding_LOH<-c(nrow(coding_func(subset(CAS6uglm, GoH_or_LoH == "LoH"))), nrow(coding_func(subset(CAS8uglm, GoH_or_LoH == "LoH"))))

CA60n_coding_LOH<-c(nrow(coding_func(subset(CAP22n, GoHorLoH == "LoH"))), nrow(coding_func(subset(CAP23n, GoHorLoH == "LoH"))), nrow(coding_func(subset(CAP26n, GoHorLoH == "LoH"))))
CA60i_coding_LOH<-c(nrow(coding_func(subset(CAP22i, GoHorLoH == "LoH"))), nrow(coding_func(subset(CAP23i, GoHorLoH == "LoH"))), nrow(coding_func(subset(CAP26i, GoHorLoH == "LoH"))))
CA60uglm_coding_LOH<-c(nrow(coding_func(subset(CAS22uglm, GoH_or_LoH == "LoH"))), nrow(coding_func(subset(CAS23uglm, GoH_or_LoH == "LoH"))), nrow(coding_func(subset(CAS26uglm, GoH_or_LoH == "LoH"))))

CA65n_coding_LOH<-c(nrow(coding_func(subset(CAP11n, GoHorLoH == "LoH"))), nrow(coding_func(subset(CAP9n, GoHorLoH == "LoH"))))
CA65i_coding_LOH<-c(nrow(coding_func(subset(CAP11i, GoHorLoH == "LoH"))), nrow(coding_func(subset(CAP9i, GoHorLoH == "LoH"))))
CA65uglm_coding_LOH<-c(nrow(coding_func(subset(CAS11uglm, GoH_or_LoH == "LoH"))), nrow(coding_func(subset(CAS9uglm, GoH_or_LoH == "LoH"))))

meanfreqpersampleplot_coding_LOH<-freq_func(CA56n_coding_LOH, CA56denom_coding, CA60n_coding_LOH, CA60denom_coding, CA65n_coding_LOH, CA65denom_coding, 
                                            CA56i_coding_LOH, CA60i_coding_LOH, CA65i_coding_LOH, CA56uglm_coding_LOH, CA60uglm_coding_LOH,CA65uglm_coding_LOH)
meanfreqpersampleplot_coding_GOH$means$Type <- "GoH"
meanfreqpersampleplot_coding_LOH$means$Type <- "LoH"

meanfreqpersampleplot_coding_combined_means <- rbind(meanfreqpersampleplot_coding_GOH$means, meanfreqpersampleplot_coding_LOH$means)
meanfreqpersampleplot_coding_combined_means$secondgroup <-factor(paste0(meanfreqpersampleplot_coding_combined_means$colony , meanfreqpersampleplot_coding_combined_means$Type, meanfreqpersampleplot_coding_combined_means$type))
#eachdatapoint:
meanfreqpersampleplot_coding_GOH$eachdatapoint$Type <- "GoH"
meanfreqpersampleplot_coding_LOH$eachdatapoint$Type <- "LoH"

meanfreqpersampleplot_coding_combined_eachdatapoint <- rbind(meanfreqpersampleplot_coding_GOH$eachdatapoint, meanfreqpersampleplot_coding_LOH$eachdatapoint)
meanfreqpersampleplot_coding_combined_eachdatapoint$secondgroup <-factor(paste0(meanfreqpersampleplot_coding_combined_eachdatapoint$colony , meanfreqpersampleplot_coding_combined_eachdatapoint$Type, meanfreqpersampleplot_coding_combined_eachdatapoint$Type))

coding_combined_plot<-ggplot(meanfreqpersampleplot_coding_combined_means,
                            aes(x = persampletypes,
                                y = meanpersamplefreqs,
                                group = secondgroup, color=Type, shape=colony)) +
geom_point(size = 3, position=position_dodge(width=0.6), aes(shape=colony)) +
#geom_point(position="dodge", stat="identity")
geom_errorbar(aes(ymin=meanpersamplefreqs-(2*se), ymax=meanpersamplefreqs+(2*se)),
              width=0.1, position=position_dodge(width=0.6)) +
theme_bw() +
theme(axis.text=element_text(size=15),
      axis.title=element_text(size=15),
      axis.text.x = element_text(angle = 90))+
labs(y="# mutations bp per sample", x="")

CA56n_GOH<-c(nrow(subset(CAP6n, GoHorLoH=="DeNovo")), nrow(subset(CAP8n, GoHorLoH=="DeNovo")))
CA56i_GOH<-c(nrow(subset(CAP6i, GoHorLoH=="DeNovo")), nrow(subset(CAP8i, GoHorLoH=="DeNovo")))
CA56uglm_GOH<-c(nrow(subset(CAS6uglm, GoH_or_LoH=="DeNovo")), nrow(subset(CAS8uglm, GoH_or_LoH=="DeNovo")))
                
CA60n_GOH<-c(nrow(subset(CAP22n, GoHorLoH=="DeNovo")), nrow(subset(CAP23n, GoHorLoH=="DeNovo")), nrow(subset(CAP26n, GoHorLoH=="DeNovo")))
CA60i_GOH<-c(nrow(subset(CAP22i, GoHorLoH=="DeNovo")), nrow(subset(CAP23i, GoHorLoH=="DeNovo")), nrow(subset(CAP26n, GoHorLoH=="DeNovo")))
CA60uglm_GOH<-c(nrow(subset(CAS22uglm, GoH_or_LoH=="DeNovo")), nrow(subset(CAS23uglm, GoH_or_LoH=="DeNovo")), nrow(subset(CAS26uglm, GoH_or_LoH=="DeNovo")))
                                                                                                
CA65n_GOH<-c(nrow(subset(CAP11n, GoHorLoH=="DeNovo")), nrow(subset(CAP9n, GoHorLoH=="DeNovo")))
CA65i_GOH<-c(nrow(subset(CAP11i, GoHorLoH=="DeNovo")), nrow(subset(CAP9i, GoHorLoH=="DeNovo")))
CA65uglm_GOH<-c(nrow(subset(CAS11uglm, GoH_or_LoH=="DeNovo")), nrow(subset(CAS9uglm, GoH_or_LoH=="DeNovo")))

meanfreqpersampleplot_GOH<-freq_func(CA56n_GOH, CA56denom, CA60n_GOH, CA60denom, CA65n_GOH, CA65denom, 
                                        CA56i_GOH, CA60i_GOH, CA65i_GOH, CA56uglm_GOH, CA60uglm_GOH,CA65uglm_GOH)

CA56n_LOH<-c(nrow(subset(CAP6n, GoHorLoH=="LoH")), nrow(subset(CAP8n, GoHorLoH=="LoH")))
CA56i_LOH<-c(nrow(subset(CAP6i, GoHorLoH=="LoH")), nrow(subset(CAP8i, GoHorLoH=="LoH")))
CA56uglm_LOH<-c(nrow(subset(CAS6uglm, GoH_or_LoH=="LoH")), nrow(subset(CAS8uglm, GoH_or_LoH=="LoH")))

CA60n_LOH<-c(nrow(subset(CAP22n, GoHorLoH=="LoH")), nrow(subset(CAP23n, GoHorLoH=="LoH")), nrow(subset(CAP26n, GoHorLoH=="LoH")))
CA60i_LOH<-c(nrow(subset(CAP22i, GoHorLoH=="LoH")), nrow(subset(CAP23i, GoHorLoH=="LoH")), nrow(subset(CAP26n, GoHorLoH=="LoH")))
CA60uglm_LOH<-c(nrow(subset(CAS22uglm, GoH_or_LoH=="LoH")), nrow(subset(CAS23uglm, GoH_or_LoH=="LoH")), nrow(subset(CAS26uglm, GoH_or_LoH=="LoH")))

CA65n_LOH<-c(nrow(subset(CAP11n, GoHorLoH=="LoH")), nrow(subset(CAP9n, GoHorLoH=="LoH")))
CA65i_LOH<-c(nrow(subset(CAP11i, GoHorLoH=="LoH")), nrow(subset(CAP9i, GoHorLoH=="LoH")))
CA65uglm_LOH<-c(nrow(subset(CAS11uglm, GoH_or_LoH=="LoH")), nrow(subset(CAS9uglm, GoH_or_LoH=="LoH")))

meanfreqpersampleplot_LOH<-freq_func(CA56n_LOH, CA56denom, CA60n_LOH, CA60denom, CA65n_LOH, CA65denom, 
                                     CA56i_LOH, CA60i_LOH, CA65i_LOH, CA56uglm_LOH, CA60uglm_LOH,CA65uglm_LOH)

meanfreqpersampleplot_GOH$means$Type <- "GoH"
meanfreqpersampleplot_LOH$means$Type <- "LoH"

meanfreqpersampleplot_combined_means <- rbind(meanfreqpersampleplot_GOH$means, meanfreqpersampleplot_LOH$means)
meanfreqpersampleplot_combined_means$secondgroup <-factor(paste0(meanfreqpersampleplot_combined_means$colony , meanfreqpersampleplot_combined_means$Type, meanfreqpersampleplot_combined_means$type))
#eachdatapoint:
meanfreqpersampleplot_GOH$eachdatapoint$Type <- "GoH"
meanfreqpersampleplot_LOH$eachdatapoint$Type <- "LoH"

meanfreqpersampleplot_combined_eachdatapoint <- rbind(meanfreqpersampleplot_GOH$eachdatapoint, meanfreqpersampleplot_LOH$eachdatapoint)
meanfreqpersampleplot_combined_eachdatapoint$secondgroup <-factor(paste0(meanfreqpersampleplot_combined_eachdatapoint$colony , meanfreqpersampleplot_combined_eachdatapoint$Type, meanfreqpersampleplot_combined_eachdatapoint$Type))

# combined_plot<-ggplot(meanfreqpersampleplot_combined, 
#                              aes(x = type, 
#                                  y = meanpersamplefreqs, 
#                                  group = secondgroup, colour=Type,shape=colony)) +
#   geom_point(size = 3, position=position_dodge(width=0.6), aes(shape=colony)) +
#   #geom_point(position="dodge", stat="identity")
#   geom_errorbar(aes(ymin=meanpersamplefreqs-(2*se), ymax=meanpersamplefreqs+(2*se)),
#                 width=0.1, position=position_dodge(width=0.6)) +
#   theme_bw() +
#   theme(axis.text=element_text(size=15),
#         axis.title=element_text(size=15),
#         axis.text.x = element_text(angle = 90))+
#   labs(y="# mutations bp per sample", x="")

meanfreqpersampleplot_combined_means$Region<- "Full genome"
meanfreqpersampleplot_coding_combined_means$Region<- "Coding region only"
meanfreqpersampleplot_combined_more_means<-rbind(meanfreqpersampleplot_combined_means, meanfreqpersampleplot_coding_combined_means)
meanfreqpersampleplot_combined_more_means$thirdgroup <-factor(paste0(meanfreqpersampleplot_combined_more_means$secondgroup, meanfreqpersampleplot_combined_more_means$Region))
meanfreqpersampleplot_combined_more_means$fourthgroup <-factor(paste0(meanfreqpersampleplot_combined_more_means$group, meanfreqpersampleplot_combined_more_means$Region))
#meanfreqpersampleplot_combined_more_means$meansamplefreqspermillion<-meanfreqpersampleplot_combined_more_means$meanpersamplefreqs*1000000
#meanfreqpersampleplot_combined_more_means$meansamplefreqspermillion_se<-meanfreqpersampleplot_combined_more_means$se*1000000
#data.frame(meanfreqpersampleplot_combined_more_means$secondgroup[1:18], rate_props)
#rate_props<-meanfreqpersampleplot_combined$meanpersamplefreqs/meanfreqpersampleplot_coding_combined$meanpersamplefreqs
meanfreqpersampleplot_combined_more_means$shared<-rep(c("PO", "P+S",
                                              "SSPO"),4)
#eachdatapoint
meanfreqpersampleplot_combined_eachdatapoint$Region<- "Full genome"
meanfreqpersampleplot_coding_combined_eachdatapoint$Region<- "Coding region only"
meanfreqpersampleplot_combined_more_eachdatapoint<-rbind(meanfreqpersampleplot_combined_eachdatapoint, meanfreqpersampleplot_coding_combined_eachdatapoint)
meanfreqpersampleplot_combined_more_eachdatapoint$thirdgroup <-factor(paste0(meanfreqpersampleplot_combined_more_eachdatapoint$secondgroup, meanfreqpersampleplot_combined_more_eachdatapoint$Region))
meanfreqpersampleplot_combined_more_eachdatapoint$fourthgroup <-factor(paste0(meanfreqpersampleplot_combined_more_eachdatapoint$group, meanfreqpersampleplot_combined_more_eachdatapoint$Region))
#meanfreqpersampleplot_combined_more_eachdatapoint$eachdatapointamplefreqspermillion<-meanfreqpersampleplot_combined_more_eachdatapoint$meanpersamplefreqs*1000000
#meanfreqpersampleplot_combined_more_eachdatapoint$eachdatapointamplefreqspermillion_se<-meanfreqpersampleplot_combined_more_eachdatapoint$se*1000000
#data.frame(meanfreqpersampleplot_combined_more_eachdatapoint$secondgroup[1:18], rate_props)
#rate_props<-meanfreqpersampleplot_combined$meanpersamplefreqs/meanfreqpersampleplot_coding_combined$meanpersamplefreqs
meanfreqpersampleplot_combined_more_eachdatapoint$shared<-rep(c(rep("PO",7), rep("P+S",7),
                                                        rep("SSPO",7)),4)
my_breaks <- function(x) { if (max(x) < 1) seq(0, 1,0.25) else seq(0, 4,1) }

combined_plot_faceted<-ggplot(meanfreqpersampleplot_combined_more, 
                      aes(x = shared, 
                          y = meansamplefreqspermillion, 
                          group = secondgroup, colour=Type,shape=Region)) +
  geom_point(size = 3, position=position_dodge(width=0.6), aes(shape=Region)) +
  #geom_point(position="dodge", stat="identity")
  geom_errorbar(aes(ymin=meansamplefreqspermillion-(meansamplefreqspermillion_se), ymax=meansamplefreqspermillion+(meansamplefreqspermillion_se)),
                width=0.1, position=position_dodge(width=0.6)) +
  theme_bw() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle = 90),
        legend.text=element_text(size=),legend.title=element_text(size=15))+
  labs(y="# SNVs per Mbp per sample", x="")#+
  facet_wrap(~colony, strip.position = "bottom", scales = "free_y")+
  theme(panel.spacing = unit(0, "lines"), 
        #strip.background = element_blank(),
        #strip.placement = "outside",
        strip.text.x = element_text(size = 15))+
  scale_y_continuous(breaks = my_breaks)
  #figure 2b:
  combined_plot_faceted<-ggplot(meanfreqpersampleplot_combined_more, 
                                aes(x = shared, 
                                    y = persamplemeans, 
                                    group = secondgroup, colour=Type,shape=Region)) +
    geom_point(size = 4.5, position=position_dodge(width=0.6), aes(shape=Region)) +
    #geom_point(position="dodge", stat="identity")
    geom_errorbar(aes(ymin=persamplemeans-(persampleses), ymax=persamplemeans+(persampleses)),
                  width=0.1, position=position_dodge(width=0.6)) +
    theme_bw() +
    theme(axis.text=element_text(size=25),
          axis.title=element_text(size=25),
          axis.text.x = element_text(angle = 0,colour=c("#999999","goldenrod","purple")),
          legend.text=element_text(size=25),legend.title=element_text(size=25))+
    labs(y="# SNVs per bp per sample", x="")  
  combined_plot_faceted2<-ggplot(meanfreqpersampleplot_combined_more_means, 
                                 aes(x = shared, 
                                     y = persamplemeans, 
                                     group = thirdgroup, shape=Region)) +
    scale_color_manual(values=c("red","blue"))+
    
    geom_pointrange(size=1.5,position=position_dodge(width=0.6), aes(color=Type, shape=Region,ymin=persamplemeans-(persampleses), ymax=persamplemeans+(persampleses))) +
    geom_jitter(data=meanfreqpersampleplot_combined_more_eachdatapoint,alpha=0.5,size=3,position = position_jitterdodge(
      jitter.width = 0.1,
      jitter.height = 0,
      dodge.width = 0.6,
      seed = NA
    ),
                aes(shared, persample, 
                    group = thirdgroup, color=Type,shape=Region))+
    theme_bw() +
    theme(axis.text=element_text(size=25),
          axis.title=element_text(size=25),
          axis.text.x = element_text(angle = 0,colour=c("#999999","goldenrod","purple")),
          legend.text=element_text(size=25),legend.title=element_text(size=25))+
    labs(y="# SNVs per bp per sample", x="")  
  
  
  meanfreqpersampleplot_combined_more$permillion<-meanfreqpersampleplot_combined_more$persamplemeans*1000000
#for ahya:
meanfreqpersampleplot_combined$shared<-rep(c(rep("PO",3), rep("P+S",3),
                                                  rep("SSPO",3)),2)
justfull<-ggplot(meanfreqpersampleplot_combined, 
                              aes(x = shared, 
                                  y = 1000000*meanpersamplefreqs, 
                                  group = secondgroup, colour=Type,shape=colony)) +
  geom_point(size = 3, position=position_dodge(width=0.6), aes(shape=colony)) +
  #geom_point(position="dodge", stat="identity")
  geom_errorbar(aes(ymin=1000000*meanpersamplefreqs-1000000*(se), ymax=1000000*meanpersamplefreqs+(1000000*se)),
                width=0.1, position=position_dodge(width=0.6)) +
  theme_bw() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle = 90),
        legend.text=element_text(size=15),legend.title=element_text(size=15))+
  labs(y="# SNVs per Mbp per sample", x="")#+
 # facet_wrap(~Region, strip.position = "bottom", scales = "free_y")+
  #theme(panel.spacing = unit(0, "lines"), 
        #strip.background = element_blank(),
        #strip.placement = "outside",
      #  strip.text.x = element_text(size = 15))+
  #scale_y_continuous(breaks = my_breaks)
#meanfreqpersampleplot_combined$meanpersamplefreqs / meanfreqpersampleplot_coding$meanpersamplefreqs
ratios<-function(genomecounts,genomedenom,codingcounts,codingdenom) {
  ro<-(genomecounts/genomedenom)/(codingcounts/codingdenom) 
  sero<-se(ro)
  mean<-mean(ro)
 # return(list("ratios"=ro,"mean"=mean,"SE"=sero)) 
  meanro<-mean(genomecounts/genomedenom)/mean(codingcounts/codingdenom) 
  return(meanro)
  #return(sero)
  }
l1<-ratios(CA56n_GOH,CA56denom,CA56n_coding_GOH,CA56denom_coding)
l2<-ratios(CA60n_GOH,CA60denom,CA60n_coding_GOH,CA60denom_coding)
ratios_n_GOH<-c(ratios(CA56n_GOH,CA56denom,CA56n_coding_GOH,CA56denom_coding), ratios(CA60n_GOH,CA60denom,CA60n_coding_GOH,CA60denom_coding), ratios(CA65n_GOH,CA65denom,CA65n_coding_GOH,CA65denom_coding))
ratios_i_GOH<-c(ratios(CA56i_GOH,CA56denom,CA56i_coding_GOH,CA56denom_coding), ratios(CA60i_GOH,CA60denom,CA60i_coding_GOH,CA60denom_coding), ratios(CA65i_GOH,CA65denom,CA65i_coding_GOH,CA65denom_coding))
ratios_uglm_GOH<-c(ratios(CA56uglm_GOH,CA56denom,CA56uglm_coding_GOH,CA56denom_coding), ratios(CA60uglm_GOH,CA60denom,CA60uglm_coding_GOH,CA60denom_coding), ratios(CA65uglm_GOH,CA65denom,CA65uglm_coding_GOH,CA65denom_coding))
ratios_n_LOH<-c( ratios(CA56n_LOH,CA56denom,CA56n_coding_LOH,CA56denom_coding), ratios(CA60n_LOH,CA60denom,CA60n_coding_LOH,CA60denom_coding), ratios(CA65n_LOH,CA65denom,CA65n_coding_LOH,CA65denom_coding))
ratios_i_LOH<-c( ratios(CA56i_LOH,CA56denom,CA56i_coding_LOH,CA56denom_coding), ratios(CA60i_LOH,CA60denom,CA60i_coding_LOH,CA60denom_coding), ratios(CA65i_LOH,CA65denom,CA65i_coding_LOH,CA65denom_coding))
ratios_uglm_LOH<-c(ratios(CA56uglm_LOH,CA56denom,CA56uglm_coding_LOH,CA56denom_coding), ratios(CA60uglm_LOH,CA60denom,CA60uglm_coding_LOH,CA60denom_coding), ratios(CA65uglm_LOH,CA65denom,CA65uglm_coding_LOH,CA65denom_coding))
ratios_n<-c(ratios(CA56n,CA56denom,CA56n_coding,CA56denom_coding), ratios(CA60n,CA60denom,CA60n_coding,CA60denom_coding), ratios(CA65n,CA65denom,CA65n_coding,CA65denom_coding))
ratios_i<-c(ratios(CA56i,CA56denom,CA56i_coding,CA56denom_coding), ratios(CA60i,CA60denom,CA60i_coding,CA60denom_coding), ratios(CA65i,CA65denom,CA65i_coding,CA65denom_coding))
ratios_uglm<-c(ratios(CA56uglm,CA56denom,CA56uglm_coding,CA56denom_coding), ratios(CA60uglm,CA60denom,CA60uglm_coding,CA60denom_coding), ratios(CA65uglm,CA65denom,CA65uglm_coding,CA65denom_coding))
ratios_LOH<-c(ratios_n_LOH,ratios_i_LOH,ratios_uglm_LOH)
ratios_GOH<-c(ratios_n_GOH,ratios_i_GOH,ratios_uglm_GOH)
mean(ratios_LOH)
mean(ratios_GOH)
se(ratios_LOH)
se(ratios_GOH)
fr<-data.frame("mean"=c(mean(ratios_LOH),mean(ratios_GOH)),"SE"=c(se(ratios_LOH),se(ratios_GOH)),Type=c("LOH","GOH"))
frplot<-ggplot(fr, 
               aes(x = Type, 
                   y = mean 
                   )) +
  geom_point(size = 3) +
  #geom_point(position="dodge", stat="identity")
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE),
                width=0.1, position=position_dodge(width=0.6)) +
  theme_bw()
ratios_means<-c(mean(ratios_n_GOH),mean(ratios_i_GOH),mean(ratios_uglm_GOH), mean(ratios_n_LOH),mean(ratios_i_LOH),mean(ratios_uglm_LOH),mean(ratios_n),mean(ratios_i),mean(ratios_uglm))
ratios_ses<-c(se(ratios_n_GOH),se(ratios_i_GOH),se(ratios_uglm_GOH), se(ratios_n_LOH),se(ratios_i_LOH),se(ratios_uglm_LOH), se(ratios_n), se(ratios_i),se(ratios_uglm))

ratios_frame<-data.frame("mean"=ratios_means, "SE"=ratios_ses, "Type"=c(rep("GOH",3),rep("LOH",3),rep("both",3)),
           "shared"=rep(c("PO","P+S","SSPO"),3))
ratiosplot<-ggplot(ratios_frame, 
                   aes(x = shared, 
                       y = mean, 
                       group = Type, colour=Type)) +
  geom_point(size = 3, position=position_dodge(width=0.6)) +
  #geom_point(position="dodge", stat="identity")
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE),
                width=0.1, position=position_dodge(width=0.6)) +
  theme_bw()

persampleaveragesplot<-
##this is figure2:
(meanfreqpersampleplot | meanfreqpersampleplot_coding) /(comparisonplotinh2 & theme(strip.placement = NULL) / inheritedpermillionplot)
##this is figure2:

(comparisonplotinh2 & theme(strip.placement = NULL)) /(combined_plot_faceted & theme(strip.placement = NULL)/ inheritedpermillionplot)

  
  comparisonplotinh2 & theme(strip.placement = NULL) | meanfreqpersampleplot
  
Figure2<-  (comparisonplotinh2 & theme(strip.placement = NULL) | persampleaveragesplot) /(meanpercentdataplot | z2)
Figure2_fordiss<- (comparisonplotinh2 & theme(strip.placement = NULL) ) /(meanpercentdataplot | z2)

#Figure 2, 6/2021:
fig2b<-combined_plot_faceted+ theme(axis.line = element_line(colour = "black")) 
#Figure 2, 20210720
(comparisonplotinh2 +labs(tag="a.")+theme(plot.tag=element_text(face="bold",size=30)) & theme(strip.placement = NULL) | meanpercentdataplot2+labs(tag="b.")+theme(plot.tag=element_text(face="bold",size=30)) ) / (combined_plot_faceted2 +labs(tag="c.")+theme(plot.tag=element_text(face="bold",size=30)) | z2+labs(tag="d.")+theme(plot.tag=element_text(face="bold",size=30)))
#for ahya:
(comparisonplotinh2 & theme(strip.placement = NULL) | justfull+labs(tag="b.")+theme(plot.tag=element_text(face="bold",size=30))) /(meanpercentdataplot | z2)
justfull+labs(tag="b.")#+theme(plot.tag=element_text(face="bold",size=30))
  #png("myplot.png")
  #print(Figure2)
  #dev.off()  
##supp table 2:

CA56_inherited$chrom.pos<-paste(CA56_inherited$chrom, ".",CA56_inherited$pos, sep="")
#c56i<-c56[match(CA56_inherited$chrom.pos, c56$chrom.pos),]

CA56_notinherited$chrom.pos<-paste(CA56_notinherited$chrom, ".",CA56_notinherited$pos, sep="")

#c56i<-c56[match(CA56_inherited$chrom.pos, c56$chrom.pos),]
m =NULL
matching_func<-function(CA60_inherited,metadatadf.0CA60) {
  m =NULL
  CA60_inherited$chrom.pos<-paste(CA60_inherited$chrom, ".",CA60_inherited$pos, sep="")
  for (i in 1:nrow(CA60_inherited)) {
    subz<-subset(metadatadf.0CA60,chrom.pos==CA60_inherited$chrom.pos[i])
    m<-rbind(m,subz)
  }
  m$Inheritance<-rep(CA60_inherited$Inheritance[i], nrow(m))
  return(m)
}  
suppCA56n<-matching_func(CA56_notinherited, metadatadf.0CA56)  
suppCA56i<-matching_func(CA56_inherited, metadatadf.0CA56)  

suppCA60n<-matching_func(CA60_notinherited, metadatadf.0CA60)  
suppCA60i<-matching_func(CA60_inherited, metadatadf.0CA60)  

suppCA65n<-matching_func(CA65_notinherited, metadatadf.0CA65)  
suppCA65i<-matching_func(CA65_inherited, metadatadf.0CA65)  
matching_func_glm<-function(CA60_inherited,metadatadf.0CA60) {
  m =NULL
  for (i in 1:nrow(CA60_inherited)) {
    subz<-subset(metadatadf.0CA60,chrom.pos==CA60_inherited$chrom.pos[i])
    m<-rbind(m,subz)
  }
  m$Inheritance<-rep("N/A", nrow(m))
  return(m)
}  
suppCA56uglm<-matching_func_glm(uniqueuglmCA56df, metadatadf.0CA56)
suppCA60uglm<-matching_func_glm(uniqueuglmCA60df, metadatadf.0CA60)
suppCA65uglm<-matching_func_glm(uniqueuglmCA65df, metadatadf.0CA65)
suppCA56gglm<-matching_func_glm(gglmCA56df, metadatadf.0CA56)
suppCA60gglm<-matching_func_glm(gglmCA60df, metadatadf.0CA60)
suppCA65gglm<-matching_func_glm(gglmCA65df, metadatadf.0CA65)

supptable2_uglm$Inheritance<-rep("N/A",nrow(supptable2_uglm))
supptable2_gglm<-allcol_gglm
supptable2_gglm$Inheritance<-rep("N/A",nrow(supptable2_gglm))

everything<-rbind(suppCA56n, suppCA56i,suppCA60n,suppCA60i,suppCA65n,suppCA65i,suppCA56uglm, suppCA60uglm, suppCA65uglm, suppCA56gglm,  suppCA60gglm, suppCA65gglm)
#suppCA56n, suppCA56i,suppCA60n,suppCA60i,suppCA65n,suppCA65i,suppCA56uglm, suppCA60uglm, suppCA65uglm, suppCA56gglm,  suppCA60gglm, suppCA65gglm
unique_mut_count<-nrow(suppCA56n)/10 + nrow(suppCA56i)/10 +nrow(suppCA60n)/12 +nrow(suppCA60i)/12 +nrow(suppCA65n)/10 +nrow(suppCA65i)/10 + 
nrow(suppCA56uglm)/10 + nrow(suppCA60uglm)/12 + nrow(suppCA65uglm)/10 + nrow(suppCA56gglm)/10 +  nrow(suppCA60gglm)/12 + nrow(suppCA65gglm)/10
everything_forsupp<-subset(everything, select = -c(mutant_allele_depth,TrueorFalse,totaldepth.y,GQscore.y,totaldepth,GQscore) )
#write.table(everything_forsupp, file="~/Documents/GitHub/CoralGermline/suppfiles/SupplementaryTable2_mappedtoahya.txt",sep="\t",quote=FALSE, row.name=FALSE)
#write.table(everything_forsupp, file="~/Documents/Dissertation/Chapter2_SupplementaryTableS2.txt",sep="\t",quote=FALSE, row.name=FALSE)
#6/25 supplementary file:
#write.table(everything_forsupp, file="~/Documents/CoralGermlineManuscript/SupplementaryTable2_20210625.txt",sep="\t",quote=FALSE, row.name=FALSE)

uniqueeverything<-everything[match(unique(everything$chrom.pos), everything$chrom.pos),]

CAP11notinhLOH
CAP22notinhLOH
CAP23notinhLOH
CAP26notinhLOH
CAP6notinhLOH
CAP8notinhLOH
CAP9notinhLOH
nrow(subset(CAS11gglm, GoH_or_LoH == "LoH"))
nrow(subset(CAS22uglm, GoH_or_LoH == "LoH"))
nrow(subset(CAS23uglm, GoH_or_LoH == "LoH"))
nrow(subset(CAS26uglm, GoH_or_LoH == "LoH"))
nrow(subset(CAS6uglm, GoH_or_LoH == "LoH"))
nrow(subset(CAS8uglm, GoH_or_LoH == "LoH"))
nrow(subset(CAS9uglm, GoH_or_LoH == "LoH"))

nrow(subset(gglmCA56df, GoH_or_LoH == "LoH"))
nrow(subset(CAS22uglm, GoH_or_LoH == "LoH"))
nrow(subset(CAS23uglm, GoH_or_LoH == "LoH"))