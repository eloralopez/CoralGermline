##how to get avg depths from unix: grep -A1 'mean' *sample_summary | grep -v 'mean' | awk '{print $3}'
CAPavgdepths<-c(14.25,
                17.12,
                19.7,
                32.95,
                34.68,
                33.25,
                24.01,
                39.49,
                40.51,
                40.75,
                45.38,
                30.9,
                34.22,
                41.19,
                49.18,
                23.1,
                45.88,
                23.95)
CASavgdepths<-c(44.07,
                80.68,
                
                
                54.97,
                42.03,
                54.93,
                18.79,
                
                
                54.41,
                20.09,
                33.91,
                56.57,
                59.49,
                59.38,
                55.23,
                93.85)
mean(CAPavgdepths)
se(CAPavgdepths)
mean(CASavgdepths)
se(CASavgdepths)
CA_ahya_samplenames<-c("CAP10-1_S44",
                       "CAP10-2_S60",
                       "CAP11-1_S41",
                       "CAP11-2_S57",
                       "CAP12-1_S46",
                       "CAP12-2_S62",
                       "CAP22-1_S43",
                       "CAP22-2_S59",
                       "CAP23-1_S45",
                       "CAP23-2_S61",
                       "CAP26-1_S39",
                       "CAP26-2_S55",
                       "CAP6-1_S47",
                       "CAP6-2_S63",
                       "CAP8-1_S40",
                       "CAP8-2_S56",
                       "CAP9-1_S48",
                       "CAP9-2_S64",
                       "CAS11-2_S3",
                       "CAS11_S67",
                       "CAS22-2_S5",
                       "CAS22_S69",
                       "CAS23-2_S7",
                       "CAS23_S71", 
                       "CAS26-2_S1",
                       "CAS26_S65",
                       "CAS6-2_S9",
                       "CAS6_S73",
                       "CAS8-2_S2",
                       "CAS8_S66",
                       "CAS9-2_S10",
                       "CAS9_S74")
CA_ahya_depths<-c(17.51,
                  21.06,
                  24.34,
                  41.09,
                  42.92,
                  41.12,
                  29.53,
                  48.87,
                  50.13,
                  50.47,
                  56.24,
                  38.1,
                  42.37,
                  51.1,
                  61.04,
                  28.4,
                  57.21,
                  29.66,
                  101.91,
                  55.13,
                  52.19,
                  68.65,
                  23.03,
                  68.68,
                  24.59,
                  67.67,
                  70.81,
                  41.97,
                  74.32,
                  74.45,
                  119.22,
                  69.5)
CA_ahya_df<-data.frame(CA_ahya_samplenames, CA_ahya_depths)
CA_ahya_df$CA_ahya_samplenames<-as.character(CA_ahya_df$CA_ahya_samplenames)
CA_ahya_df_parents<-subset(CA_ahya_df, startsWith(CA_ahya_df$CA_ahya_samplenames, 'CAP')==TRUE)
CA_ahya_df_sperm<-subset(CA_ahya_df, startsWith(CA_ahya_df$CA_ahya_samplenames, 'CAS')==TRUE)
mean(CA_ahya_df_parents$CA_ahya_depths)
se(CA_ahya_df_parents$CA_ahya_depths)
mean(CA_ahya_df_sperm$CA_ahya_depths)
se(CA_ahya_df_sperm$CA_ahya_depths)
AHBavgdepths<-c(33.73,
                
                26.70,
                
                30.39,
                
                25.44,
                
                25.87,
                
                23.60,
                
                31.78,
                
                35.19,
                
                25.56,
                
                51.68,
                
                32.58,
                
                27.80,
                
                26.88,
                
                25.08,
                
                30.44,
                
                33.27,
                
                37.17,
                
                0.00,
                
                31.13,
                
                31.30,
                
                31.23,
                
                13.88,
                
                31.86,
                
                28.38,
                
                8.22,
                
                30.00,
                
                31.84,
                
                35.87,
                
                29.70,
                
                7.72,
                
                7.96,
                
                25.76,
                
                26.63,
                
                33.57,
                
                29.46,
                
                24.72,
                
                25.38,
                
                24.90,
                
                8.85,
                
                28.12,
                
                28.42,
                
                26.78,
                
                24.31,
                
                34.47,
                
                35.55,
                
                17.11,
                
                27.49,
                
                29.78,
                
                25.97
)

AHBsamplenames<-c(    "AHB125_S36",
                      "AHB126-1_S37",
                      "AHB126-2_S61",
                      "AHB127_S38",
                      "AHB128-1_S39",
                      "AHB128-2_S62",
                      "AHB129_S40",
                      "AHB145_S15",
                      "AHB146-1_S16",
                      "AHB146-2_S71",
                      "AHB147_S17",
                      "AHB148-1_S18",
                      "AHB148-2_S72",
                      "AHB149_S19",
                      "AHB151_S41",
                      "AHB152-1_S42",
                      "AHB152-2_S63",
                      "AHB153_S43",
                      "AHB154-1_S44",
                      "AHB154-2_S64",
                      "AHB155_S45",
                      "AHB176_S54",
                      "AHB177-1_S5",
                      "AHB177-2_S65",
                      "AHB178_S6",
                      "AHB179-1_S7",
                      "AHB179-2_S66",
                      "AHB180_S8",
                      "AHB195_S2",
                      "AHB196-1_S3",
                      "AHB196-2_S79",
                      "AHB197_S4",
                      "AHB198-1_S22",
                      "AHB198-2_S80",
                      "AHB199_S35",
                      "AHB70_S20",
                      "AHB71-1_S56",
                      "AHB71-2_S81",
                      "AHB72_S24",
                      "AHB73-1_S13",
                      "AHB73-2_S82",
                      "AHB74_S57",
                      "AHB90_S9",
                      "AHB91-1_S10",
                      "AHB91-2_S67",
                      "AHB92_S55",
                      "AHB93-1_S11",
                      "AHB93-2_S68",
                      "AHB94_S12"
)

AHPavgdepths<- c(43.37,
                 
                 43.78,
                 
                 21.54,
                 
                 27.10,
                 
                 21.76,
                 
                 28.50,
                 
                 28.12,
                 
                 9.33,
                 
                 15.83,
                 
                 17.84,
                 
                 4.15,
                 
                 3.26,
                 
                 31.72,
                 
                 58.64,
                 
                 37.37,
                 
                 21.80,
                 
                 45.95,
                 
                 20.46,
                 
                 4.98,
                 
                 4.24,
                 
                 30.44,
                 
                 30.95,
                 
                 22.28,
                 
                 48.64,
                 
                 23.62,
                 
                 26.49,
                 
                 30.59,
                 
                 20.58,
                 
                 37.68,
                 
                 38.70,
                 
                 36.59,
                 
                 37.05,
                 
                 33.34,
                 
                 40.48,
                 
                 32.63)
AHPsamplenames<-c(	"AHP01-1_S11",
                   "AHP01-2_S27",
                   "AHP02_S33",
                   "AHP03-1_S12",
                   "AHP03-2_S28",
                   "AHP04_S34",
                   "AHP05_S49",
                   "AHP06-1_S13",
                   "AHP06-2_S29",
                   "AHP07_S50",
                   "AHP08-1_S14",
                   "AHP08-2_S30",
                   "AHP09_S75",
                   "AHP10_S76",
                   "AHP11_S77",
                   "AHP12-1_S15",
                   "AHP12-2_S31",
                   "AHP13_S78",
                   "AHP14-1_S16",
                   "AHP14-2_S32",
                   "AHP15_S79",
                   "AHP16-1_S35",
                   "AHP16-2_S51",
                   "AHP17_S80",
                   "AHP18-1_S36",
                   "AHP18-2_S52",
                   "AHP19_S81",
                   "AHP20_S82",
                   "AHP21-1_S37",
                   "AHP21-2_S53",
                   "AHP22-1_S38",
                   "AHP22-2_S54",
                   "AHP23_S83",
                   "AHP24_S84",
                   "AHP25_S85")

AHASavgdepths<-c(9.89,
                 
                 11.53,
                 
                 13.45,
                 
                 25.21,
                 
                 17.53,
                 
                 18.92,
                 
                 0.19,
                 
                 25.36,
                 
                 7.08,
                 
                 8.39,
                 
                 21.08,
                 
                 11.79,
                 
                 12.78,
                 
                 17.27,
                
                 21.53,
                 
                 13.06,
                 
                 31.90,
                 
                 22.94,
                 
                 17.69,
                 
                 28.57,
                 
                 23.53,
                 
                 11.88,
                 
                 14.09,
                 
                 17.58,
                 
                 24.14,
                 
                 9.50,
                 
                 6.98,
                 
                 20.90,
                 
                 24.28,
                 
                 28.15,
                 
                 28.07,
                 
                 17.86,
                 
                 18.05,
                 
                 18.96,
                 
                 18.86
)
AHASsamplenames<-c(	"AHAS41_S49",
                    "AHAS42-1_S50",
                    "AHAS42-2_S73",
                    "AHAS43_S25",
                    "AHAS44-1_S26",
                    "AHAS44-2_S74",
                    "AHAS45_S21",
                    "AHAS46_S27",
                    "AHAS47-1_S51",
                    "AHAS47-2_S75",
                    "AHAS48_S52",
                    "AHAS49-1_S1",
                    "AHAS49-2_S76",
                    "AHAS50_S28",
                    "AHAS51_S29",
                    "AHAS52-1_S30",
                    "AHAS52-2_S69",
                    "AHAS53_S31",
                    "AHAS54-1_S46",
                    "AHAS54-2_S70",
                    "AHAS55_S47",
                    "AHAS56_S32",
                    "AHAS57-1_S53",
                    "AHAS57-2_S77",
                    "AHAS58_S33",
                    "AHAS59-1_S23",
                    "AHAS59-2_S78",
                    "AHAS60_S34",
                    "AHAS61_S48",
                    "AHAS62-1_S58",
                    "AHAS62-2_S83",
                    "AHAS63_S59",
                    "AHAS64-1_S60",
                    "AHAS64-2_S84",
                    "AHAS65_S14")

AHEavgdepths<-c(6.35,
                
                4.47,
                
                14.27,
                
                13.15,
                
               
                
                6.23,
                
                4.27,
                
                20.38,
                
                20.17,
                
                7.65,
                
                7.14,
                
                4.49,
                
                4.40,
                
                20.03,
                
                28.00,
                
                15.89,
                
                15.37)
AHEsamplenames<-c(	"AHE01-1_S1",
                   "AHE01-2_S17",
                   "AHE02-1_S2",
                   "AHE02-2_S18",
                   "AHE05-1_S3",
                   "AHE05-2_S19",
                   "AHE06-1_S4",
                   "AHE06-2_S20",
                   "AHE07-1_S5",
                   "AHE07-2_S21",
                   "AHE08-1_S6",
                   "AHE08-2_S22",
                   "AHE09-1_S7",
                   "AHE09-2_S23",
                   "AHE10-1_S8",
                   "AHE10-2_S24",
                   "AHE11-1_S9",
                   "AHE11-2_S25",
                   "AHE12-1_S10",
                   "AHE12-2_S26")
depths_func<-function(samplenames, avgdepths) {
  upperbounds<-qpois(0.9999, avgdepths,lower.tail=T)  
  lowerbounds<-qpois(0.0001,avgdepths,lower.tail=T)  
  depthdf<-data.frame("samplename"=samplenames, "avgdepth"=avgdepths,"upper"=upperbounds,"lower"=lowerbounds)
  #return(depthdf)
  for (i in 1:length(samplenames)) {
    #sc<-print(paste0("vc.getGenotype(", samplenames[i],".getDP()<", upperbounds[i], 
    #"&& vc.getGenotype(", samplenames[i],").getDP()>",lowerbounds[i], "&&"))
    pc<-cat("vc.getGenotype(", "\"\'", samplenames[i],"\'\"", ")",".getDP()<", upperbounds[i], 
            " && vc.getGenotype(", "\"\'", samplenames[i],"\'\"", ").getDP()>",lowerbounds[i], " &&", "\n", sep = "")
    pc
  } 
  return(pc)
}  
CA_ahya<-depths_func(CA_ahya_samplenames, CA_ahya_depths)
AHB<-depths_func(AHBsamplenames,AHBavgdepths)
AHP<-depths_func(AHPsamplenames,AHPavgdepths)
AHAS<-depths_func(AHASsamplenames,AHASavgdepths)
AHE<-depths_func(AHEsamplenames,AHEavgdepths)
AHEdf<-data.frame(AHEavgdepths, AHEsamplenames)
AHE01avg<-mean(AHEdf$AHEavgdepths[1:4])
AHE07avg<-mean(AHEdf$AHEavgdepths[9:12])
AHE09avg<-mean(AHEdf$AHEavgdepths[13:16])
AHE11avg<-mean(AHEdf$AHEavgdepths[17:20])

AHASFOURSAMPLES<-data.frame(AHASavgdepths, AHASsamplenames)[c(2:3,5:6,9:10,12:13,16:17,19:20,23:24,26:27,30:31,33:34),]
AHAS41avg<-mean(AHASFOURSAMPLES$AHASavgdepths[1:4])
AHAS46avg<-mean(AHASFOURSAMPLES$AHASavgdepths[5:8])
AHAS51avg<-mean(AHASFOURSAMPLES$AHASavgdepths[9:12])
AHAS56avg<-mean(AHASFOURSAMPLES$AHASavgdepths[13:16])
AHAS61avg<-mean(AHASFOURSAMPLES$AHASavgdepths[17:20])

AHAS4<-depths_func(as.character(AHASFOURSAMPLES[,2]), AHASFOURSAMPLES[,1])
AHBFOURSAMPLES<-data.frame(AHBavgdepths, AHBsamplenames)[c(2:3,5:6,9:10,12:13,16:17,19:20,23:24,26:27,30:31,33:34,37:38,40:41,44:45,47:48),]
AHB125avg<-mean(AHBFOURSAMPLES$AHBavgdepths[1:4])
AHB145avg<-mean(AHBFOURSAMPLES$AHBavgdepths[5:8])
AHB151avg<-mean(AHBFOURSAMPLES$AHBavgdepths[9:12])
AHB176avg<-mean(AHBFOURSAMPLES$AHBavgdepths[13:16])
AHB195avg<-mean(AHBFOURSAMPLES$AHBavgdepths[17:20])
AHB70avg<-mean(AHBFOURSAMPLES$AHBavgdepths[21:24])
AHB90avg<-mean(AHBFOURSAMPLES$AHBavgdepths[25:28])

AHB4<-depths_func(as.character(AHBFOURSAMPLES[,2]), AHBFOURSAMPLES[,1])
AHPFOURSAMPLES<-data.frame(AHPavgdepths, AHPsamplenames)[c(1:2,4:5,8:9,11:12,22:23,25:26,29:30,31:32),]
AHP01avg<-mean(AHPFOURSAMPLES$AHPavgdepths[1:4])
AHP06avg<-mean(AHPFOURSAMPLES$AHPavgdepths[5:8])
AHP16avg<-mean(AHPFOURSAMPLES$AHPavgdepths[9:12])
AHP21avg<-mean(AHPFOURSAMPLES$AHPavgdepths[13:16])

AHP4<-depths_func(as.character(AHPFOURSAMPLES[,2]), AHPFOURSAMPLES[,1])

mean(data.frame(AHPavgdepths, AHPsamplenames)[c(1:14,22:35),]$AHPavgdepths)
se(data.frame(AHPavgdepths, AHPsamplenames)[c(1:14,22:35),]$AHPavgdepths)

mean(data.frame(AHASavgdepths, AHASsamplenames)[c(2:3,5:6,8:35),]$AHASavgdepths)
se(data.frame(AHASavgdepths, AHASsamplenames)[c(2:3,5:6,8:35),]$AHASavgdepths)
for (i in 1:length(AHBsamplenames)) {
  #sc<-print(paste0("vc.getGenotype(", AHBsamplenames[i],".getDP()<", AHBupperbounds[i], 
               #"&& vc.getGenotype(", AHBsamplenames[i],").getDP()>",AHBlowerbounds[i], "&&"))
  pc<-cat("vc.getGenotype(", "\"\'", AHBsamplenames[i],"\'\"", ")",".getDP()<", AHBupperbounds[i], 
                   " && vc.getGenotype(", "\"\'", AHBsamplenames[i],"\'\"", ").getDP()>",AHBlowerbounds[i], " &&", "\n", sep = "")
  pc
}  
print("\" )
x<-cat('"',''')

x<-cat("\"\'")
y<-cat("\'\"")
z<-cat(")")


combined_plot_faceted2<-ggplot(meanfreqpersampleplot_combined_more, 
                              aes(x = shared, 
                                  y = meansamplefreqspermillion, 
                                  group = secondgroup, colour=Type,shape=corno)) +
  geom_point(size = 3, position=position_dodge(width=0.6), aes(shape=corno)) +
  #geom_point(position="dodge", stat="identity")
  geom_errorbar(aes(ymin=meansamplefreqspermillion-(meansamplefreqspermillion_se), ymax=meansamplefreqspermillion+(meansamplefreqspermillion_se)),
                width=0.1, position=position_dodge(width=0.6)) +
  theme_bw() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle = 90),
        legend.text=element_text(size=15),legend.title=element_text(size=15))+
  labs(y="# SNVs per Mbp per sample", x="")
#  facet_wrap(~corno, strip.position = "bottom", scales = "free_y")+
#  theme(panel.spacing = unit(0, "lines"), 
        #strip.background = element_blank(),
        #strip.placement = "outside",
#        strip.text.x = element_text(size = 15))+
#  scale_y_continuous(breaks = my_breaks)
colonyaverages<-c(mean(meanfreqpersampleplot_combined_more$meansamplefreqspermillion[1:3]),
                  mean(meanfreqpersampleplot_combined_more$meansamplefreqspermillion[4:6]),
                  mean(meanfreqpersampleplot_combined_more$meansamplefreqspermillion[7:9]),
                  mean(meanfreqpersampleplot_combined_more$meansamplefreqspermillion[10:12]),
                  mean(meanfreqpersampleplot_combined_more$meansamplefreqspermillion[13:15]),
                  mean(meanfreqpersampleplot_combined_more$meansamplefreqspermillion[16:18]),
                  
                  mean(meanfreqpersampleplot_combined_more$meansamplefreqspermillion[18:21]),
                  mean(meanfreqpersampleplot_combined_more$meansamplefreqspermillion[22:24]),
                  mean(meanfreqpersampleplot_combined_more$meansamplefreqspermillion[25:27]),
                  mean(meanfreqpersampleplot_combined_more$meansamplefreqspermillion[28:30]),
                  mean(meanfreqpersampleplot_combined_more$meansamplefreqspermillion[31:33]),
                  mean(meanfreqpersampleplot_combined_more$meansamplefreqspermillion[34:36]))
persampleaverages<-c(                   mean(c(CA56n_GOH/CA56denom, CA60n_GOH/CA60denom,CA65n_GOH/CA65denom)),
                                        mean(c(CA56i_GOH/CA56denom, CA60i_GOH/CA60denom,CA65i_GOH/CA65denom)),
                                        mean(c(CA56uglm_GOH/CA56denom, CA60uglm_GOH/CA60denom,CA65uglm_GOH/CA65denom)),
                                        mean(c(CA56n_LOH/CA56denom, CA60n_LOH/CA60denom,CA65n_LOH/CA65denom)),
                                        mean(c(CA56i_LOH/CA56denom, CA60i_LOH/CA60denom,CA65i_LOH/CA65denom)),
                                        mean(c(CA56uglm_LOH/CA56denom, CA60uglm_LOH/CA60denom,CA65uglm_LOH/CA65denom)),
                   mean(c(CA56n_coding_GOH/CA56denom_coding, CA60n_coding_GOH/CA60denom_coding,CA65n_coding_GOH/CA65denom_coding)),
                   mean(c(CA56i_coding_GOH/CA56denom_coding, CA60i_coding_GOH/CA60denom_coding,CA65i_coding_GOH/CA65denom_coding)),
                   mean(c(CA56uglm_coding_GOH/CA56denom_coding, CA60uglm_coding_GOH/CA60denom_coding,CA65uglm_coding_GOH/CA65denom_coding)),
                   mean(c(CA56n_coding_LOH/CA56denom_coding, CA60n_coding_LOH/CA60denom_coding,CA65n_coding_LOH/CA65denom_coding)),
                   mean(c(CA56i_coding_LOH/CA56denom_coding, CA60i_coding_LOH/CA60denom_coding,CA65i_coding_LOH/CA65denom_coding)),
                   mean(c(CA56uglm_coding_LOH/CA56denom_coding, CA60uglm_coding_LOH/CA60denom_coding,CA65uglm_coding_LOH/CA65denom_coding)))
persampleses<-c(                   se(c(CA56n_GOH/CA56denom, CA60n_GOH/CA60denom,CA65n_GOH/CA65denom)),
                                        se(c(CA56i_GOH/CA56denom, CA60i_GOH/CA60denom,CA65i_GOH/CA65denom)),
                                        se(c(CA56uglm_GOH/CA56denom, CA60uglm_GOH/CA60denom,CA65uglm_GOH/CA65denom)),
                                        se(c(CA56n_LOH/CA56denom, CA60n_LOH/CA60denom,CA65n_LOH/CA65denom)),
                                        se(c(CA56i_LOH/CA56denom, CA60i_LOH/CA60denom,CA65i_LOH/CA65denom)),
                                        se(c(CA56uglm_LOH/CA56denom, CA60uglm_LOH/CA60denom,CA65uglm_LOH/CA65denom)),
                                        se(c(CA56n_coding_GOH/CA56denom_coding, CA60n_coding_GOH/CA60denom_coding,CA65n_coding_GOH/CA65denom_coding)),
                                        se(c(CA56i_coding_GOH/CA56denom_coding, CA60i_coding_GOH/CA60denom_coding,CA65i_coding_GOH/CA65denom_coding)),
                                        se(c(CA56uglm_coding_GOH/CA56denom_coding, CA60uglm_coding_GOH/CA60denom_coding,CA65uglm_coding_GOH/CA65denom_coding)),
                                        se(c(CA56n_coding_LOH/CA56denom_coding, CA60n_coding_LOH/CA60denom_coding,CA65n_coding_LOH/CA65denom_coding)),
                                        se(c(CA56i_coding_LOH/CA56denom_coding, CA60i_coding_LOH/CA60denom_coding,CA65i_coding_LOH/CA65denom_coding)),
                                        se(c(CA56uglm_coding_LOH/CA56denom_coding, CA60uglm_coding_LOH/CA60denom_coding,CA65uglm_coding_LOH/CA65denom_coding)))
n_GOH<-c(CA56n_GOH/CA56denom, CA60n_GOH/CA60denom,CA65n_GOH/CA65denom)
n_GOH_coding<-c(CA56n_coding_GOH/CA56denom_coding, CA60n_coding_GOH/CA60denom_coding,CA65n_coding_GOH/CA65denom_coding)
i_GOH<-c(CA56i_GOH/CA56denom, CA60i_GOH/CA60denom,CA65i_GOH/CA65denom)
i_GOH_coding<-c(CA56i_coding_GOH/CA56denom_coding, CA60i_coding_GOH/CA60denom_coding,CA65i_coding_GOH/CA65denom_coding)
uglm_GOH<-c(CA56uglm_GOH/CA56denom, CA60uglm_GOH/CA60denom,CA65uglm_GOH/CA65denom)
uglm_GOH_coding<-c(CA56uglm_coding_GOH/CA56denom_coding, CA60uglm_coding_GOH/CA60denom_coding,CA65uglm_coding_GOH/CA65denom_coding)
n_LOH<-c(CA56n_LOH/CA56denom, CA60n_LOH/CA60denom,CA65n_LOH/CA65denom)
n_LOH_coding<-c(CA56n_coding_LOH/CA56denom_coding, CA60n_coding_LOH/CA60denom_coding,CA65n_coding_LOH/CA65denom_coding)
i_LOH<-c(CA56i_LOH/CA56denom, CA60i_LOH/CA60denom,CA65i_LOH/CA65denom)
i_LOH_coding<-c(CA56i_coding_LOH/CA56denom_coding, CA60i_coding_LOH/CA60denom_coding,CA65i_coding_LOH/CA65denom_coding)
uglm_LOH<-c(CA56uglm_LOH/CA56denom, CA60uglm_LOH/CA60denom,CA65uglm_LOH/CA65denom)
uglm_LOH_coding<-c(CA56uglm_coding_LOH/CA56denom_coding, CA60uglm_coding_LOH/CA60denom_coding,CA65uglm_coding_LOH/CA65denom_coding)

vars<-c(var(n_GOH), var(i_GOH), var(uglm_GOH), var(n_LOH), var(i_LOH), var(uglm_LOH))
vars_coding<-c(var(n_GOH_coding), var(i_GOH_coding), var(uglm_GOH_coding), var(n_LOH_coding), var(i_LOH_coding), var(uglm_LOH_coding))
##The rate of SNVs per Mbp was significantly higher across the full callable genome (see Methods for calculation of callable genome size) than the rate of SNVs in the callable coding regions of the genome all SNV types except for LOH shared by parent and sperm (see Supplementary for all t-test results) (Fig. 2b). 
t.test(n_GOH*1000000, n_GOH_coding*1000000, alternative=c("two.sided"),paired=TRUE)
t.test(i_GOH*1000000, i_GOH_coding*1000000, alternative=c("two.sided"),paired=TRUE)
t.test(uglm_GOH*1000000, uglm_GOH_coding*1000000, alternative=c("two.sided"),paired=TRUE)
t.test(n_LOH*1000000, n_LOH_coding*1000000, alternative=c("two.sided"),paired=TRUE)
t.test(i_LOH*1000000, i_LOH_coding*1000000, alternative=c("two.sided"),paired=TRUE)
t.test(uglm_LOH*1000000, uglm_LOH_coding*1000000, alternative=c("two.sided"),paired=TRUE)

meanpercolony<-data.frame("mean"=persampleaverages*1000000,"se"=persampleses*1000000,"Type"=c(rep("GOH",3),rep("LOH",3)),
                          "Region"=c(rep("Full genome",6),rep("Coding only",6)),"shared"=rep(c("PO","P+S","SSPO"),2))
persampleaveragesplot<-ggplot(meanpercolony, 
                               aes(x = shared, 
                                   y = mean, 
                                   group = Type, colour=Type,shape=Region)) +
  geom_point(size = 3, position=position_dodge(width=0.6), aes(shape=Region)) +
  #geom_point(position="dodge", stat="identity")
  geom_errorbar(aes(ymin=mean-(se), ymax=mean+(se)),
                width=0.1, position=position_dodge(width=0.6)) +
  theme_bw() +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        axis.text.x = element_text(angle = 90),
        legend.text=element_text(size=15),legend.title=element_text(size=15))+
  labs(y="# SNVs per Mbp per sample", x="")
rats<-meanpercolony$mean[1:6]/meanpercolony$mean[7:12]
fullgenomerates<- c(CA56n_GOH/CA56denom, CA60n_GOH/CA60denom,CA65n_GOH/CA65denom,
                  CA56i_GOH/CA56denom, CA60i_GOH/CA60denom,CA65i_GOH/CA65denom,
                  CA56uglm_GOH/CA56denom, CA60uglm_GOH/CA60denom,CA65uglm_GOH/CA65denom,
                  CA56n_LOH/CA56denom, CA60n_LOH/CA60denom,CA65n_LOH/CA65denom,
                  CA56i_LOH/CA56denom, CA60i_LOH/CA60denom,CA65i_LOH/CA65denom,
                  CA56uglm_LOH/CA56denom, CA60uglm_LOH/CA60denom,CA65uglm_LOH/CA65denom)
codingrates<-c(CA56n_coding_GOH/CA56denom_coding, CA60n_coding_GOH/CA60denom_coding,CA65n_coding_GOH/CA65denom_coding,
            CA56i_coding_GOH/CA56denom_coding, CA60i_coding_GOH/CA60denom_coding,CA65i_coding_GOH/CA65denom_coding,
            CA56uglm_coding_GOH/CA56denom_coding, CA60uglm_coding_GOH/CA60denom_coding,CA65uglm_coding_GOH/CA65denom_coding,
            CA56n_coding_LOH/CA56denom_coding, CA60n_coding_LOH/CA60denom_coding,CA65n_coding_LOH/CA65denom_coding,
            CA56i_coding_LOH/CA56denom_coding, CA60i_coding_LOH/CA60denom_coding,CA65i_coding_LOH/CA65denom_coding,
            CA56uglm_coding_LOH/CA56denom_coding, CA60uglm_coding_LOH/CA60denom_coding,CA65uglm_coding_LOH/CA65denom_coding)

rateratios<-fullgenomerates/codingrates
rateratios<-meanpercolony$mean[1:6]/meanpercolony$mean[7:12]

data.frame(AHASavgdepths, AHASsamplenames)
fourpopdepth<-c(37.62,
                29.74,
                33.86,
                28.33,
                28.8,
                26.25,
                35.47,
                39.75,
                28.81,
                58.49,
                36.74,
                31.37,
                30.31,
                28.27,
                34.77,
                38.01,
                30.75,
                35.6,
                35.79,
                35.69,
                15.38,
                35.48,
                31.54,
                9.08,
                33.38,
                35.4,
                39.98,
                33.27,
                8.55,
                8.83,
                28.81,
                29.69,
                37.53,
                32.99,
                27.56,
                28.34,
                27.82,
                9.77,
                31.43,
                31.67,
                29.94,
                27.1,
                38.54,
                39.74,
                19.01,
                30.7,
                33.24,
                28.98,
                11.64,
                13.39,
                15.66,
                30.03,
                20.86,
                22.48,
                0.22,
                29.13,
                8.18,
                9.71,
                24.11,
                13.42,
                14.57,
                19.76,
                23.86,
                1.64,
                14.44,
                35.43,
                25.7,
                19.61,
                31.74,
                26.12,
                2.72,
                15.67,
                19.62,
                26.95,
                10.52,
                7.72,
                23.26,
                27.11,
                31.52,
                31.41,
                19.84,
                20.11,
                21.15,
                21.01,
                7.02,
                4.93,
                15.85,
                14.59,
                1.91,
                1.71,
                0.71,
                0.77,
                6.85,
                4.68,
                22.6,
                22.35,
                8.47,
                7.9,
                4.94,
                4.84,
                22.15,
                30.96,
                17.52,
                16.93,
                52.85,
                53.32,
                25.88,
                32.55,
                26.09,
                34.44,
                34.05,
                10.38,
                17.66,
                19.96,
                4.61,
                3.63,
                35.6,
                66.39,
                46.36,
                26.88,
                57.19,
                25.17,
                5.95,
                5.07,
                38.35,
                27.5,
                32.68,
                25.48,
                43.93,
                42.02,
                37.93,
                37.01)
fourpopname<-c("AHB125_S36",
               "AHB126-1_S37",
               "AHB126-2_S61",
               "AHB127_S38",
               "AHB128-1_S39",
               "AHB128-2_S62",
               "AHB129_S40",
               "AHB145_S15",
               "AHB146-1_S16",
               "AHB146-2_S71",
               "AHB147_S17",
               "AHB148-1_S18",
               "AHB148-2_S72",
               "AHB149_S19",
               "AHB151_S41",
               "AHB152-1_S42",
               "AHB153_S43",
               "AHB154-1_S44",
               "AHB154-2_S64",
               "AHB155_S45",
               "AHB176_S54",
               "AHB177-1_S5",
               "AHB177-2_S65",
               "AHB178_S6",
               "AHB179-1_S7",
               "AHB179-2_S66",
               "AHB180_S8",
               "AHB195_S2",
               "AHB196-1_S3",
               "AHB196-2_S79",
               "AHB197_S4",
               "AHB198-1_S22",
               "AHB198-2_S80",
               "AHB199_S35",
               "AHB70_S20",
               "AHB71-1_S56",
               "AHB71-2_S81",
               "AHB72_S24",
               "AHB73-1_S13",
               "AHB73-2_S82",
               "AHB74_S57",
               "AHB90_S9",
               "AHB91-1_S10",
               "AHB91-2_S67",
               "AHB92_S55",
               "AHB93-1_S11",
               "AHB93-2_S68",
               "AHB94_S12",
               "AHAS41_S49",
               "AHAS42-1_S50",
               "AHAS42-2_S73",
               "AHAS43_S25",
               "AHAS44-1_S26",
               "AHAS44-2_S74",
               "AHAS45_S21",
               "AHAS46_S27",
               "AHAS47-1_S51",
               "AHAS47-2_S75",
               "AHAS48_S52",
               "AHAS49-1_S1",
               "AHAS49-2_S76",
               "AHAS50_S28",
               "AHAS51_S29",
               "AHAS51_S29",
               "AHAS52-1_S30",
               "AHAS52-2_S69",
               "AHAS53_S31",
               "AHAS54-1_S46",
               "AHAS54-2_S70",
               "AHAS55_S47",
               "AHAS56_S32",
               "AHAS57-1_S53",
               "AHAS57-2_S77",
               "AHAS58_S33",
               "AHAS59-1_S23",
               "AHAS59-2_S78",
               "AHAS60_S34",
               "AHAS61_S48",
               "AHAS62-1_S58",
               "AHAS62-2_S83",
               "AHAS63_S59",
               "AHAS64-1_S60",
               "AHAS64-2_S84",
               "AHAS65_S14",
               "AHE01-1_S1",
               "AHE01-2_S17",
               "AHE02-1_S2",
               "AHE02-2_S18",
               "AHE05-1_S3",
               "AHE05-2_S19",
               "AHE06-1_S4",
               "AHE06-2_S20",
               "AHE07-1_S5",
               "AHE07-2_S21",
               "AHE08-1_S6",
               "AHE08-2_S22",
               "AHE09-1_S7",
               "AHE09-2_S23",
               "AHE10-1_S8",
               "AHE10-2_S24",
               "AHE11-1_S9",
               "AHE11-2_S25",
               "AHE12-1_S10",
               "AHE12-2_S26",
               "AHP01-1_S11",
               "AHP01-2_S27",
               "AHP02_S33",
               "AHP03-1_S12",
               "AHP03-2_S28",
               "AHP04_S34",
               "AHP05_S49",
               "AHP06-1_S13",
               "AHP06-2_S29",
               "AHP07_S50",
               "AHP08-1_S14",
               "AHP08-2_S30",
               "AHP09_S75",
               "AHP10_S76",
               "AHP11_S77",
               "AHP12-1_S15",
               "AHP12-2_S31",
               "AHP13_S78",
               "AHP14-1_S16",
               "AHP14-2_S32",
               "AHP16-1_S35",
               "AHP16-2_S51",
               "AHP18-2_S52",
               "AHP20_S82",
               "AHP21-2_S53",
               "AHP22-2_S54",
               "AHP23_S83",
               "AHP25_S85")
depths_func(fourpopname, fourpopdepth)

names<-c("AHB152-2_S63","AHP18-1_S36","AHP21-1_S37","AHP22-1_S38","AHP17_S80","AHP19_S81", "AHP24_S84")
depths<-c(42.46,29.1, 42.76, 41.49, 61.09, 37.99, 46.05)  
depths_func(names,depths)
