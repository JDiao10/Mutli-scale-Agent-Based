library(dplyr)
library(ggplot2)
library(ggdist)
library(R.matlab)
library(latex2exp)




# Fast-TIV
NoMTFluD = readMat("MT100_Flu_phi.mat")
#NoMTFluDAll = readMat("MT100_Flu_phiNoAll.mat")
#NoMTFluDAda = readMat("MT100_Flu_phiNoAda.mat")
NoMTFlu = NoMTFluD$phi.no
NoMTFluAll = NoMTFluD$phi.no.all
NoMTFluAda = NoMTFluD$phi.no.Ada
NoMTFlu95 = NoMTFlu[NoMTFlu>=quantile(NoMTFlu,0.025)&NoMTFlu<=quantile(NoMTFlu,0.975)]
NoMTFluAll95 = NoMTFluAll[NoMTFluAll>=quantile(NoMTFluAll,0.025)&NoMTFluAll<=quantile(NoMTFluAll,0.975)]
NoMTFluAda95 = NoMTFluAda[NoMTFluAda>=quantile(NoMTFluAda,0.025)&NoMTFluAda<=quantile(NoMTFluAda,0.975)]
MCNum=94
Color = c(rep("Isolate all",MCNum),rep("Isolate high",MCNum),rep("$Adaptive_{0}$",MCNum))
dataFull_No = data.frame(x = c(rep("Isolate all",MCNum),rep("Isolate high",MCNum),rep("$Adaptive_{0}$",MCNum)),
                         values = c(NoMTFluAll95,NoMTFlu95,NoMTFluAda95),Groups = Color)
pFull_No<-ggplot(dataFull_No,aes(x=x,y=values,fill=Groups))+geom_violin(scale = 'width')+
  stat_summary(fun = median,geom ="point",size = 2,color="black")+
  scale_x_discrete(limits=c(limits=c("Isolate all", "Isolate high",'$Adaptive_{0}$')),labels=TeX)+
  scale_fill_manual(values=c("#b8e186","#276419","#762a83"))+
  xlab("Strategy")+
  ylab(TeX(r'(Relative effectiveness ($\varphi_{s}$))'))+
  ggtitle("Multiple tests")+
  theme_bw()+
  theme(axis.text = element_text(size=18),axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.line=element_line(colour='black'),
        panel.border=element_blank(),
        legend.position = 'none',plot.title = element_text(size=18,hjust=0.5),
        plot.margin = unit(c(1.1, 0.2, 0.1, 0.5), "cm"))+
  ylim(0,1)
ggsave("All_violin_No_FluMT_v1.pdf",plot=pFull_No,width=8,height=6)



NoSTFluD = readMat("ST100_Flu_phiNo.mat")
NoSTFluDAll = readMat("ST100_Flu_phiNoAll.mat")
NoSTFluDAda = readMat("ST100_Flu_phiNoAda.mat")
NoSTFlu = NoSTFluD$phi
NoSTFluAll = NoSTFluDAll$phi
NoSTFluAda = NoSTFluDAda$phi
NoSTFlu95 = NoSTFlu[NoSTFlu>=quantile(NoSTFlu,0.025)&NoSTFlu<=quantile(NoSTFlu,0.975)]
NoSTFluAll95 = NoSTFluAll[NoSTFluAll>=quantile(NoSTFluAll,0.025)&NoSTFluAll<=quantile(NoSTFluAll,0.975)]
NoSTFluAda95 = NoSTFluAda[NoSTFluAda>=quantile(NoSTFluAda,0.025)&NoSTFluAda<=quantile(NoSTFluAda,0.975)]
MCNum=94
Color = c(rep("Isolate all",MCNum),rep("Isolate high",MCNum),rep("$Adaptive_{0}$",MCNum))
dataFull_NoST = data.frame(x = c(rep("Isolate all",MCNum),rep("Isolate high",MCNum),rep("$Adaptive_{0}$",MCNum)),
                         values = c(NoSTFluAll95,NoSTFlu95,NoSTFluAda95),Groups = Color)
pFull_NoST<-ggplot(dataFull_NoST,aes(x=x,y=values,fill=Groups))+geom_violin(scale = 'width')+
  stat_summary(fun = median,geom ="point",size = 2,color="black")+
  scale_x_discrete(limits=c(limits=c("Isolate all", "Isolate high","$Adaptive_{0}$")),labels=TeX)+
  scale_fill_manual(values=c("#b8e186","#276419","#762a83"))+
  xlab("Strategy")+
  ylab(TeX(r'(Relative effectiveness ($varphi_{s}$))'))+
  ggtitle("Single test")+
  theme_bw()+
  theme(axis.text = element_text(size=18),axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.line=element_line(colour='black'),
        panel.border=element_blank(),
        legend.position = 'none',plot.title = element_text(size=18,hjust=0.5),
        plot.margin = unit(c(1.1, 0.2, 0.1, 0.5), "cm"))+
  ylim(0,1)
ggsave("All_violin_No_FluST_v1.pdf",plot=pFull_NoST,width=8,height=6)

NoFluD_sym = readMat("MT100_Flu_phi_sym.mat")
NoFlu_sym = NoFluD_sym$phi.no
NoFluAll_sym = NoFluD_sym$phi.no.all
NoFluAda_sym = NoFluD_sym$phi.no.Ada
NoFlu95_sym = NoFlu_sym[NoFlu_sym>=quantile(NoFlu_sym,0.025)&NoFlu_sym<=quantile(NoFlu_sym,0.975)]
NoFluAll95_sym = NoFluAll_sym[NoFluAll_sym>=quantile(NoFluAll_sym,0.025)&NoFluAll_sym<=quantile(NoFluAll_sym,0.975)]
NoFluAda95_sym = NoFluAda_sym[NoFluAda_sym>=quantile(NoFluAda_sym,0.025)&NoFluAda_sym<=quantile(NoFluAda_sym,0.975)]
MCNum=94
Color = c(rep("Isolate all",MCNum),rep("Isolate high",MCNum),rep("$Adaptive_{0}$",MCNum))
dataFull_No_sym = data.frame(x = c(rep("Isolate all",MCNum),rep("Isolate high",MCNum),rep("$Adaptive_{0}$",MCNum)),
                           values = c(NoFluAll95_sym,NoFlu95_sym,NoFluAda95_sym),Groups = Color)
pFull_No_sym<-ggplot(dataFull_No_sym,aes(x=x,y=values,fill=Groups))+geom_violin(scale = 'width')+
  stat_summary(fun = median,geom ="point",size = 2,color="black")+
  scale_x_discrete(limits=c(limits=c("Isolate all", "Isolate high","$Adaptive_{0}$")),labels=TeX)+
  scale_fill_manual(values=c("#b8e186","#276419","#762a83"))+
  xlab("Strategy")+
  ylab(TeX(r'(Relative effectiveness ($varphi_{s}$))'))+
  ggtitle("")+
  theme_bw()+
  theme(axis.text = element_text(size=18),axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.line=element_line(colour='black'),
        panel.border=element_blank(),
        legend.position = 'none',plot.title = element_text(size=18,hjust=0.5),
        plot.margin = unit(c(1.1, 0.2, 0.1, 0.5), "cm"))+
  ylim(0,1)
ggsave("All_violin_No_Flu_sym.pdf",plot=pFull_No_sym,width=8,height=6)


IsolationST = readMat("EffectivenessST.mat")
IsolationMT = readMat("EffectivenessMT.mat")
IsolationMT_sym = readMat("EffectivenessMT_asy.mat")
ST_H = IsolationST$e.HFinal
ST_B = IsolationST$e.BFinal
ST_SF = IsolationST$e.SFFinal
MT_H = IsolationMT$e.HFinal
MT_B = IsolationMT$e.BFinal
MT_SF = IsolationMT$e.SFFinal
MT_Hs = IsolationMT_sym$e.HFinal
MT_Bs = IsolationMT_sym$e.BFinal
MT_SFs = IsolationMT_sym$e.SFFinal


ST_H95 = ST_H[ST_H>=quantile(ST_H,0.025)&ST_H<=quantile(ST_H,0.975)]
ST_B95 = ST_B[ST_B>=quantile(ST_B,0.025)&ST_B<=quantile(ST_B,0.975)]
ST_SF95 = ST_SF[ST_SF>=quantile(ST_SF,0.025)&ST_SF<=quantile(ST_SF,0.975)]
MT_B95 = MT_B[MT_B>=quantile(MT_B,0.025)&MT_B<=quantile(MT_B,0.975)]
MT_H95 = MT_H[MT_H>=quantile(MT_H,0.025)&MT_H<=quantile(MT_H,0.975)]
MT_SF95 = MT_SF[MT_SF>=quantile(MT_SF,0.025)&MT_SF<=quantile(MT_SF,0.975)]
MT_B95s = MT_Bs[MT_Bs>=quantile(MT_Bs,0.025)&MT_Bs<=quantile(MT_Bs,0.975)]
MT_H95s = MT_Hs[MT_Hs>=quantile(MT_Hs,0.025)&MT_Hs<=quantile(MT_Hs,0.975)]
MT_SF95s = MT_SFs[MT_SFs>=quantile(MT_SFs,0.025)&MT_SFs<=quantile(MT_SFs,0.975)]

#ColorG4=["#984ea3","#4daf4a"];
MCNum = 94
Color = c(rep("Isolate all",MCNum),rep("Isolate high",MCNum),rep("$Adaptive_{0}$",MCNum))
dataFullIsolation = data.frame(x = c(rep("Isolate all",MCNum),rep("Isolate high",MCNum),rep("$Adaptive_{0}$",MCNum)),
                               values = c(MT_B95,MT_H95,MT_SF95),Groups = Color)
pIsolation_Flu<-ggplot(dataFullIsolation,aes(x=x,y=values,fill=Groups))+geom_violin()+
  stat_summary(fun = median,geom ="point",size = 2,color="black")+
  scale_x_discrete(limits=c(limits=c("Isolate all", "Isolate high","$Adaptive_{0}$")),labels=TeX)+
  scale_fill_manual(values=c("#b8e186","#276419","#762a83"))+
  xlab("Strategy")+
  ylab(TeX("Efficiency ($\\E_{s}$)"))+
  ggtitle("Multiple tests")+
  theme_bw()+
  theme(axis.text = element_text(size=18,family='sans'),axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.line=element_line(colour='black'),
        panel.border=element_blank(),
        legend.position = 'none',plot.title = element_text(size=18,hjust=0.5),
        plot.margin = unit(c(1.1, 0.2, 0.1, 0.5), "cm"))+
  ylim(0,8)
ggsave("Violin_IsolationMT_v1.pdf",plot=pIsolation_Flu,width=8,height=6)



MCNum = 94
Color = c(rep("Isolate all",MCNum),rep("Isolate high",MCNum),rep("$Adaptive_{0}$",MCNum))
dataFullIsolationST = data.frame(x = c(rep("Isolate all",MCNum),rep("Isolate high",MCNum),rep("$Adaptive_{0}$",MCNum)),
                               values = c(ST_B95,ST_H95,ST_SF95),Groups = Color)
pIsolation_FluST<-ggplot(dataFullIsolationST,aes(x=x,y=values,fill=Groups))+geom_violin()+
  stat_summary(fun = median,geom ="point",size = 2,color="black")+
  scale_x_discrete(limits=c(limits=c("Isolate all", "Isolate high","$Adaptive_{0}$")),labels=TeX)+
  scale_fill_manual(values=c("#b8e186","#276419","#762a83"))+
  xlab("Strategy")+
  ylab(TeX("Efficiency ($E_{s}$)"))+
  ggtitle("Single test")+
  theme_bw()+
  theme(axis.text = element_text(size=18,family='sans'),axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.line=element_line(colour='black'),
        panel.border=element_blank(),
        legend.position = 'none',plot.title = element_text(size=18,hjust=0.5),
        plot.margin = unit(c(1.1, 0.2, 0.1, 0.5), "cm"))+
  ylim(0,2)
ggsave("Violin_IsolationST_v1.pdf",plot=pIsolation_FluST,width=8,height=6)



MCNum = 94
Color = c(rep("Isolate all",MCNum),rep("Isolate high",MCNum),rep("$Adaptive_{0}$",MCNum))
dataFullIsolations = data.frame(x = c(rep("Isolate all",MCNum),rep("Isolate high",MCNum),rep("$Adaptive_{0}$",MCNum)),
                               values = c(MT_B95s,MT_H95s,MT_SF95s),Groups = Color)
pIsolation_Flus<-ggplot(dataFullIsolations,aes(x=x,y=values,fill=Groups))+geom_violin()+
  stat_summary(fun = median,geom ="point",size = 2,color="black")+
  scale_x_discrete(limits=c(limits=c("Isolate all", "Isolate high","$Adaptive_{0}$")),labels=TeX)+
  scale_fill_manual(values=c("#b8e186","#276419","#762a83"))+
  xlab("Strategy")+
  ylab(TeX("Efficiency ($E_{s}$)"))+
  ggtitle("")+
  theme_bw()+
  theme(axis.text = element_text(size=18,family='sans'),axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.line=element_line(colour='black'),
        panel.border=element_blank(),
        legend.position = 'none',plot.title = element_text(size=18,hjust=0.5),
        plot.margin = unit(c(1.1, 0.2, 0.2, 0.5), "cm"))+
  ylim(0,2)
ggsave("Violin_IsolationMT_asym_v1.pdf",plot=pIsolation_Flus,width=8,height=6)




HFDiffTime = readMat("MT_Flu_phi_DiifTime100.mat")
HF = readMat("MT100_Flu_phi.mat")
MTF__35 = HFDiffTime$phi..35
MTF__21 = HFDiffTime$phi..21
MTF__7 = HFDiffTime$phi..7
MTF_7 = HFDiffTime$phi.7
MTF_21 = HFDiffTime$phi.21
MTF_35 = HFDiffTime$phi.35
MTF = HF$phi
MCNum = 94

MTF__3595 = MTF__35[MTF__35>=quantile(MTF__35,0.025)&MTF__35<=quantile(MTF__35,0.975)]
MTF__2195 = MTF__21[MTF__21>=quantile(MTF__21,0.025)&MTF__21<=quantile(MTF__21,0.975)]
MTF__795 = MTF__7[MTF__7>=quantile(MTF__7,0.025)&MTF__7<=quantile(MTF__7,0.975)]
MTF_3595 = MTF_35[MTF_35>=quantile(MTF_35,0.025)&MTF_35<=quantile(MTF_35,0.975)]
MTF_2195 = MTF_21[MTF_21>=quantile(MTF_21,0.025)&MTF_21<=quantile(MTF_21,0.975)]
MTF_795 = MTF_7[MTF_7>=quantile(MTF_7,0.025)&MTF_7<=quantile(MTF_7,0.975)]
MTF95 = MTF[MTF>=quantile(MTF,0.025)&MTF<=quantile(MTF,0.975)]
Color = c(rep("-5",MCNum),rep("-3",MCNum),rep("-1",MCNum),rep("0",MCNum),rep("+1",MCNum),rep("+3",MCNum),rep("+5",MCNum))
dataFullDiff = data.frame(x = c(rep("-5",MCNum),rep("-3",MCNum),rep("-1",MCNum),rep("0",MCNum),rep("+1",MCNum),rep("+3",MCNum),rep("+5",MCNum)),
                          values = c(1-MTF__3595,1-MTF__2195,1-MTF__795,1-MTF95,1-MTF_795,1-MTF_2195,1-MTF_3595),Groups = Color)
pF<-ggplot(dataFullDiff,aes(x=x,y=values,fill=Groups))+geom_violin()+
  stat_summary(fun = median,geom ="point",size = 2,color="black")+
  scale_x_discrete(limits=c("-5", "-3", "-1","0","+1","+3","+5"))+
  scale_fill_manual(values=c("#762a83","#9970ab","#c2a5cf","#b8e186","#7fbc41","#4d9221","#276419"),
                    limits=c("-5", "-3", "-1","0","+1","+3","+5"))+
  xlab(TeX(r'(Offset from the peak day (weeks))'))+
  ylab(TeX(r'(Relative cost ($\psi$))'))+
  theme_bw()+
  theme(axis.text = element_text(size=18,family='sans'),axis.title.x=element_text(size=18,vjust=-1),
        axis.title.y=element_text(size=18),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.line=element_line(colour='black'),
        panel.border=element_blank(),
        legend.position = 'none',
        plot.margin = unit(c(1.1, 0.2, 0.2, 0.5), "cm"))+
  ylim(-0.5, 1.8)
ggsave("ViolinDiffTime.pdf",plot=pF,width=6,height=6)


HFDiffTimeST = readMat("ST_Flu_phi_DiifTime100.mat")
SF = readMat("ST100_Flu_phi.mat")
STF__35 = HFDiffTimeST$phi..35
STF__21 = HFDiffTimeST$phi..21
STF__7 = HFDiffTimeST$phi..7
STF_7 = HFDiffTimeST$phi.7
STF_21 = HFDiffTimeST$phi.21
STF_35 = HFDiffTimeST$phi.35
STF = SF$phi
MCNum = 94
STF__3595 = STF__35[STF__35>=quantile(STF__35,0.025)&STF__35<=quantile(STF__35,0.975)]
STF__2195 = STF__21[STF__21>=quantile(STF__21,0.025)&STF__21<=quantile(STF__21,0.975)]
STF__795 = STF__7[STF__7>=quantile(STF__7,0.025)&STF__7<=quantile(STF__7,0.975)]
STF_3595 = STF_35[STF_35>=quantile(STF_35,0.025)&STF_35<=quantile(STF_35,0.975)]
STF_2195 = STF_21[STF_21>=quantile(STF_21,0.025)&STF_21<=quantile(STF_21,0.975)]
STF_795 = STF_7[STF_7>=quantile(STF_7,0.025)&STF_7<=quantile(STF_7,0.975)]
STF95 = STF[STF>=quantile(STF,0.025)&STF<=quantile(STF,0.975)]
Color = c(rep("-5",MCNum),rep("-3",MCNum),rep("-1",MCNum),rep("0",MCNum),rep("+1",MCNum),rep("+3",MCNum),rep("+5",MCNum))
dataFullDiff = data.frame(x = c(rep("-5",MCNum),rep("-3",MCNum),rep("-1",MCNum),rep("0",MCNum),rep("+1",MCNum),rep("+3",MCNum),rep("+5",MCNum)),
                          values = c(1-STF__3595,1-STF__2195,1-STF__795,1-STF95,1-STF_795,1-STF_2195,1-STF_3595),Groups = Color)
pFST<-ggplot(dataFullDiff,aes(x=x,y=values,fill=Groups))+geom_violin(scale = 'width')+
  stat_summary(fun = median,geom ="point",size = 2,color="black")+
  scale_x_discrete(limits=c("-5", "-3", "-1","0","+1","+3","+5"))+
  scale_fill_manual(values=c("#762a83","#9970ab","#c2a5cf","#b8e186","#7fbc41","#4d9221","#276419"),
                    limits=c("-5", "-3", "-1","0","+1","+3","+5"))+
  xlab(TeX(r'(Offset from the peak day (weeks))'))+
  ylab(TeX(r'(Relative cost ($\psi$))'))+
  theme_bw()+
  theme(axis.text = element_text(size=18,family='sans'),axis.title.x=element_text(size=18,vjust=-1),
        axis.title.y=element_text(size=18),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.line=element_line(colour='black'),
        panel.border=element_blank(),
        legend.position = 'none',plot.margin = unit(c(1.1, 0.2, 0.2, 0.5), "cm"))+
  ylim(-0.5, 1.8)
ggsave("ViolinDiffTimeST.pdf",plot=pFST,width=6,height=6)


EffDiffTime = readMat("EffectivenessMT_DiffTime.mat")
EMTF__35 = EffDiffTime$e..35
EMTF__21 = EffDiffTime$e..21
EMTF__7 = EffDiffTime$e..7
EMTF__0 = EffDiffTime$e..0
EMTF_7 = EffDiffTime$e.7
EMTF_21 = EffDiffTime$e.21
EMTF_35 = EffDiffTime$e.35
EMTF__3595 = EMTF__35[EMTF__35>=quantile(EMTF__35,0.025)&EMTF__35<=quantile(EMTF__35,0.975)]
EMTF__2195 = EMTF__21[EMTF__21>=quantile(EMTF__21,0.025)&EMTF__21<=quantile(EMTF__21,0.975)]
EMTF__795 = EMTF__7[EMTF__7>=quantile(EMTF__7,0.025)&EMTF__7<=quantile(EMTF__7,0.975)]
EMTF__095 = EMTF__0[EMTF__0>=quantile(EMTF__0,0.025)&EMTF__0<=quantile(EMTF__0,0.975)]
EMTF_3595 = EMTF_35[EMTF_35>=quantile(EMTF_35,0.025)&EMTF_35<=quantile(EMTF_35,0.975)]
EMTF_2195 = EMTF_21[EMTF_21>=quantile(EMTF_21,0.025)&EMTF_21<=quantile(EMTF_21,0.975)]
EMTF_795 = EMTF_7[EMTF_7>=quantile(EMTF_7,0.025)&EMTF_7<=quantile(EMTF_7,0.975)]
MCNum=94
Color = c(rep("-5",MCNum),rep("-3",MCNum),rep("-1",MCNum),rep("0",MCNum),rep("+1",MCNum),rep("+3",MCNum),rep("+5",MCNum))
dataFullEffDiff = data.frame(x = c(rep("-5",MCNum),rep("-3",MCNum),rep("-1",MCNum),rep("0",MCNum),rep("+1",MCNum),rep("+3",MCNum),rep("+5",MCNum)),
                          values = c(EMTF__3595,EMTF__2195,EMTF__795,EMTF__095,EMTF_795,EMTF_2195,EMTF_3595),Groups = Color)
pFEff<-ggplot(dataFullEffDiff,aes(x=x,y=values,fill=Groups))+geom_violin()+
  stat_summary(fun = median,geom ="point",size = 2,color="black")+
  scale_x_discrete(limits=c("-5", "-3", "-1","0","+1","+3","+5"))+
  scale_fill_manual(values=c("#762a83","#9970ab","#c2a5cf","#b8e186","#7fbc41","#4d9221","#276419"),
                    limits=c("-5", "-3", "-1","0","+1","+3","+5"))+
  xlab(TeX(r'(Offset from the peak day (weeks))'))+
  ylab(TeX("Efficiency ($E_{Adaptive_n}$)"))+
  theme_bw()+
  theme(axis.text = element_text(size=18,family='sans'),axis.title.x=element_text(size=18,vjust=-1),
        axis.title.y=element_text(size=18),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.line=element_line(colour='black'),
        panel.border=element_blank(),
        legend.position = 'none',plot.margin = unit(c(1.1, 0.2, 0.2, 0.5), "cm"))+
  ylim(0,7)
ggsave("ViolinEffDiffTimeMT.pdf",plot=pFEff,width=6,height=6)


EffDiffTime = readMat("EffectivenessST_DiffTime.mat")
ESTF__35 = EffDiffTime$e..35
ESTF__21 = EffDiffTime$e..21
ESTF__7 = EffDiffTime$e..7
ESTF__0 = EffDiffTime$e..0
ESTF_7 = EffDiffTime$e.7
ESTF_21 = EffDiffTime$e.21
ESTF_35 = EffDiffTime$e.35
ESTF__3595 = ESTF__35[ESTF__35>=quantile(ESTF__35,0.025)&ESTF__35<=quantile(ESTF__35,0.975)]
ESTF__2195 = ESTF__21[ESTF__21>=quantile(ESTF__21,0.025)&ESTF__21<=quantile(ESTF__21,0.975)]
ESTF__795 = ESTF__7[ESTF__7>=quantile(ESTF__7,0.025)&ESTF__7<=quantile(ESTF__7,0.975)]
ESTF__095 = ESTF__0[ESTF__0>=quantile(ESTF__0,0.025)&ESTF__0<=quantile(ESTF__0,0.975)]
ESTF_3595 = ESTF_35[ESTF_35>=quantile(ESTF_35,0.025)&ESTF_35<=quantile(ESTF_35,0.975)]
ESTF_2195 = ESTF_21[ESTF_21>=quantile(ESTF_21,0.025)&ESTF_21<=quantile(ESTF_21,0.975)]
ESTF_795 = ESTF_7[ESTF_7>=quantile(ESTF_7,0.025)&ESTF_7<=quantile(ESTF_7,0.975)]
MCNum=94
Color = c(rep("-5",MCNum),rep("-3",MCNum),rep("-1",MCNum),rep("0",MCNum),rep("+1",MCNum),rep("+3",MCNum),rep("+5",MCNum))
dataFullEffDiffS = data.frame(x = c(rep("-5",MCNum),rep("-3",MCNum),rep("-1",MCNum),rep("0",MCNum),rep("+1",MCNum),rep("+3",MCNum),rep("+5",MCNum)),
                             values = c(ESTF__3595,ESTF__2195,ESTF__795,ESTF__095,ESTF_795,ESTF_2195,ESTF_3595),Groups = Color)
pFEffS<-ggplot(dataFullEffDiffS,aes(x=x,y=values,fill=Groups))+geom_violin()+
  stat_summary(fun = median,geom ="point",size = 2,color="black")+
  scale_x_discrete(limits=c("-5", "-3", "-1","0","+1","+3","+5"))+
  scale_fill_manual(values=c("#762a83","#9970ab","#c2a5cf","#b8e186","#7fbc41","#4d9221","#276419"),
                    limits=c("-5", "-3", "-1","0","+1","+3","+5"))+
  xlab(TeX(r'(Offset from the peak day (weeks))'))+
  ylab(TeX("Efficiency ($E_{Adaptive_n}$)"))+
  theme_bw()+
  theme(axis.text = element_text(size=18,family='sans'),axis.title.x=element_text(size=18,vjust=-1),
        axis.title.y=element_text(size=18),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.line=element_line(colour='black'),
        panel.border=element_blank(),
        legend.position = 'none',,plot.margin = unit(c(1.1, 0.2, 0.2, 0.5), "cm"))+
  ylim(0,2)
ggsave("ViolinEffDiffTimeST.pdf",plot=pFEffS,width=6,height=6)



# Slow-TIV

# Supplemental materials figures



HFsym= readMat("MT100_Flu_phi_sym.mat")
MTFsym = HFsym$phi
SFS_No = readMat("MT100_Flu_phi_symNo.mat")
STFsym = SFS_No$phi
Color = c(rep("Fast",MCNum),rep("Fast",MCNum))
dataFullasym = data.frame(x = c(rep("No",MCNum),rep("all",MCNum)),
                          values = c(STFsym,MTFsym),Groups = Color)
pFull_Fluasym<-ggplot(dataFullasym,aes(x=x,y=values,fill=Groups))+geom_violin()+
  stat_summary(fun = median,geom ="point",size = 2,color="black")+
  scale_x_discrete(limits=c(limits=c("No", "all")))+
  scale_fill_manual(values=c("#0072BD"))+
  xlab("")+
  ylab(TeX("$varphi_{No}$ and $varphi$"))+
  theme_bw()+
  theme(axis.text = element_text(size=18),axis.title.x=element_text(size=22),
        axis.title.y=element_text(size=22),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),panel.border=element_rect(colour='black',size=1),
        legend.position = 'none',
        plot.margin = unit(c(1.1, 0.2, 0.1, 0.2), "cm"))+
  ylim(-0.2,1.3)
ggsave("All_violin_Fluasym.pdf",plot=pFull_Fluasym,width=8,height=6)


# Slow TIV

NoMTCD = readMat("MT100_COVID_phi.mat")
NoMTC = NoMTCD$phi.no
NoMTCAll = NoMTCD$phi.no.all
NoMTCAda = NoMTCD$phi.no.Ada
NoMTC95 = NoMTC[NoMTC>=quantile(NoMTC,0.025)&NoMTC<=quantile(NoMTC,0.975)]
NoMTCAll95 = NoMTCAll[NoMTCAll>=quantile(NoMTCAll,0.025)&NoMTCAll<=quantile(NoMTCAll,0.975)]
NoMTCAda95 = NoMTCAda[NoMTCAda>=quantile(NoMTCAda,0.025)&NoMTCAda<=quantile(NoMTCAda,0.975)]
MCNum=94
Color = c(rep("Isolate all",MCNum),rep("Isolate high",MCNum),rep("Adaptive",MCNum))
dataCl_No = data.frame(x = c(rep("Isolate all",MCNum),rep("Isolate high",MCNum),rep("Adaptive",MCNum)),
                         values = c(NoMTCAll95,NoMTC95,NoMTCAda95),Groups = Color)
pCl_No<-ggplot(dataCl_No,aes(x=x,y=values,fill=Groups))+geom_violin(scale = 'width')+
  stat_summary(fun = median,geom ="point",size = 2,color="black")+
  scale_x_discrete(limits=c(limits=c("Isolate all", "Isolate high","Adaptive_{0}")),labels=TeX)+
  scale_fill_manual(values=c("#b8e186","#276419","#762a83"))+
  xlab("Strategy")+
  ylab(TeX(r'(Relative effectiveness ($varphi_{s}$))'))+
  ggtitle("Multiple tests")+
  theme_bw()+
  theme(axis.text = element_text(size=18),axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.line=element_line(colour='black'),
        panel.border=element_blank(),
        legend.position = 'none',plot.title = element_text(size=18,hjust=0.5),
        plot.margin = unit(c(1.1, 0.2, 0.1, 0.5), "cm"))+
  ylim(0,1)
ggsave("All_violin_No_CMT.pdf",plot=pCl_No,width=8,height=6)

NoSTCD = readMat("ST100_COVID_phi.mat")
NoSTC = NoSTCD$phi.no
NoSTCAll = NoSTCD$phi.no.all
NoSTCAda = NoSTCD$phi.no.Ada
NoSTC95 = NoSTC[NoSTC>=quantile(NoSTC,0.025)&NoSTC<=quantile(NoSTC,0.975)]
NoSTCAll95 = NoSTCAll[NoSTCAll>=quantile(NoSTCAll,0.025)&NoSTCAll<=quantile(NoSTCAll,0.975)]
NoSTCAda95 = NoSTCAda[NoSTCAda>=quantile(NoSTCAda,0.025)&NoSTCAda<=quantile(NoSTCAda,0.975)]
MCNum=94
Color = c(rep("Isolate all",MCNum),rep("Isolate high",MCNum),rep("Adaptive",MCNum))
dataCl_NoS = data.frame(x = c(rep("Isolate all",MCNum),rep("Isolate high",MCNum),rep("Adaptive",MCNum)),
                       values = c(NoSTCAll95,NoSTC95,NoSTCAda95),Groups = Color)
pCl_NoS<-ggplot(dataCl_NoS,aes(x=x,y=values,fill=Groups))+geom_violin(scale = 'width')+
  stat_summary(fun = median,geom ="point",size = 2,color="black")+
  scale_x_discrete(limits=c(limits=c("Isolate all", "Isolate high","Adaptive_{0}")),labels=TeX)+
  scale_fill_manual(values=c("#b8e186","#276419","#762a83"))+
  xlab("Strategy")+
  ylab(TeX(r'(Relative effectiveness ($varphi_{s}$))'))+
  ggtitle("Single test")+
  theme_bw()+
  theme(axis.text = element_text(size=18),axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.line=element_line(colour='black'),
        panel.border=element_blank(),
        legend.position = 'none',plot.title = element_text(size=18,hjust=0.5),
        plot.margin = unit(c(1.1, 0.2, 0.1, 0.5), "cm"))+
  ylim(0,1)
ggsave("All_violin_No_CST.pdf",plot=pCl_NoS,width=8,height=6)

IsolationSTC = readMat("EffectivenessSTCovid.mat")
IsolationMTC = readMat("EffectivenessMTCovid.mat")
ST_HC = IsolationSTC$e.HFinal
ST_BC = IsolationSTC$e.BFinal
ST_SFC = IsolationSTC$e.SFFinal
MT_HC = IsolationMTC$e.HFinal
MT_BC = IsolationMTC$e.BFinal
MT_SFC = IsolationMTC$e.SFFinal
ST_H95C = ST_HC[ST_HC>=quantile(ST_HC,0.025)&ST_HC<=quantile(ST_HC,0.975)]
ST_B95C = ST_BC[ST_BC>=quantile(ST_BC,0.025)&ST_BC<=quantile(ST_BC,0.975)]
ST_SF95C = ST_SFC[ST_SFC>=quantile(ST_SFC,0.025)&ST_SFC<=quantile(ST_SFC,0.975)]
MT_B95C = MT_BC[MT_BC>=quantile(MT_BC,0.025)&MT_BC<=quantile(MT_BC,0.975)]
MT_H95C = MT_HC[MT_HC>=quantile(MT_HC,0.025)&MT_HC<=quantile(MT_HC,0.975)]
MT_SF95C = MT_SFC[MT_SFC>=quantile(MT_SFC,0.025)&MT_SFC<=quantile(MT_SFC,0.975)]


#ColorG4=["#984ea3","#4daf4a"];
MCNum = 94
Color = c(rep("Isolate all",MCNum),rep("Isolate high",MCNum),rep("Adaptive",MCNum))
dataFullIsolationC = data.frame(x = c(rep("Isolate all",MCNum),rep("Isolate high",MCNum),rep("Adaptive",MCNum)),
                                values = c(MT_B95C,MT_H95C,MT_SF95C),Groups = Color)
dataFullIsolationCS = data.frame(x=c(rep("Isolate all",MCNum),rep("Isolate high",MCNum),rep("Adaptive",MCNum)),
                                 values = c(ST_B95C,ST_H95C,ST_SF95C),Groups = Color)
                                
pIsolation_C<-ggplot(dataFullIsolationC,aes(x=x,y=values,fill=Groups))+geom_violin()+
  stat_summary(fun = median,geom ="point",size = 2,color="black")+
  scale_x_discrete(limits=c(limits=c("Isolate all", "Isolate high","Adaptive_{0}")))+
  scale_fill_manual(values=c("Isolate all"="#276419","Isolate high"="#762a83","Adaptive"="#b8e186"))+
  xlab("Strategy")+
  ylab(TeX(r'(Efficiency ($E_s$))'))+
  ggtitle("Mulitple tests")+
  theme_bw()+
  theme(axis.text = element_text(size=18),axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.line=element_line(colour='black'),
        panel.border=element_blank(),
        legend.position = 'none',plot.title = element_text(size=18,hjust=0.5),
        plot.margin = unit(c(1.1, 0.2, 0.1, 0.5), "cm"))+
  ylim(0,4)
ggsave("Violin_Isolation_C.pdf",plot=pIsolation_C,width=8,height=6)

pIsolation_CS<-ggplot(dataFullIsolationCS,aes(x=x,y=values,fill=Groups))+geom_violin()+
  stat_summary(fun = median,geom ="point",size = 2,color="black")+
  scale_x_discrete(limits=c(limits=c("Isolate all", "Isolate high","Adaptive_{0}")))+
  scale_fill_manual(values=c("Isolate all"="#276419","Isolate high"="#762a83","Adaptive"="#b8e186"))+
  xlab("Strategy")+
  ylab(TeX(r'(Efficiency ($E_s$))'))+
  ggtitle("Single test")+
  theme_bw()+
  theme(axis.text = element_text(size=18),axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),axis.line=element_line(colour='black'),
        panel.border=element_blank(),
        legend.position = 'none',plot.title = element_text(size=18,hjust=0.5),
        plot.margin = unit(c(1.1, 0.2, 0.1, 0.5), "cm"))+
  ylim(0,0.5)
ggsave("Violin_Isolation_CS.pdf",plot=pIsolation_CS,width=8,height=6)



IsolationPerST = readMat("EfficiencyPerDayST.mat")
IsolationPerMT = readMat("EfficiencyPerDayMT.mat")
#solationPerMT_sym = readMat("EffectivenessMT_asy.mat")
STPer_H = IsolationPerST$ePerTime.H
STPer_B = IsolationPerST$ePerTime.B
STPer_SF = IsolationPerST$ePerTime.SF
MTPer_H = IsolationPerMT$ePerTime.H
MTPer_B = IsolationPerMT$ePerTime.B
MTPer_SF = IsolationPerMT$ePerTime.SF
#MT_Hs = IsolationMT_sym$e.HFinal
#MT_Bs = IsolationMT_sym$e.BFinal
#MT_SFs = IsolationMT_sym$e.SFFinal


STPer_H95 = STPer_H[STPer_H>=quantile(STPer_H,0.025)&STPer_H<=quantile(STPer_H,0.975)]
STPer_B95 = STPer_B[STPer_B>=quantile(STPer_B,0.025)&STPer_B<=quantile(STPer_B,0.975)]
STPer_SF95 = STPer_SF[STPer_SF>=quantile(STPer_SF,0.025)&STPer_SF<=quantile(STPer_SF,0.975)]
MTPer_B95 = MTPer_B[MTPer_B>=quantile(MTPer_B,0.025)&MTPer_B<=quantile(MTPer_B,0.975)]
MTPer_H95 = MTPer_H[MTPer_H>=quantile(MTPer_H,0.025)&MTPer_H<=quantile(MTPer_H,0.975)]
MTPer_SF95 = MTPer_SF[MTPer_SF>=quantile(MTPer_SF,0.025)&MTPer_SF<=quantile(MTPer_SF,0.975)]
#MTPer_B95s = MT_Bs[MT_Bs>=quantile(MT_Bs,0.025)&MT_Bs<=quantile(MT_Bs,0.975)]
#MTPer_H95s = MT_Hs[MT_Hs>=quantile(MT_Hs,0.025)&MT_Hs<=quantile(MT_Hs,0.975)]
#MT_SF95s = MT_SFs[MT_SFs>=quantile(MT_SFs,0.025)&MT_SFs<=quantile(MT_SFs,0.975)]
