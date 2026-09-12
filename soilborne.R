
########################NWS soilbore plant pathogen R code
library(Rmisc)
library(ggplot2)
library(ggsignif)
library(stringr)
library(vegan)
library(readxl)
library(tidyr)
library(Matrix)
library(reshape2)
library(tidyverse)
library(ieggr)
library(ggpubr)
library(microeco)
library(ggprism)
library(ggalluvial)
library(dplyr)
library(ggpmisc)


library(lme4)
library(car)
library(lmerTest)
library(MuMIn)

######https://mp.weixin.qq.com/s/cLOs61U3sq_42OrGDZ9JJw
###########安装 FUNGuildR
###devtools::install_github('brendanf/FUNGuildR')
#下载数据库（执行一次保存在本地即可，以后直接在本地调用）
#fung <- get_funguild_db()
#saveRDS(fung, 'funguild.rds')
#setwd("/Volumes/SHUN HAN/data from Xue/soil-borne pathogens")
setwd("/Users/data from Xue/soil-borne pathogens/raw data")
library(FUNGuildR)
fung <- readRDS('funguild.rds')  #加载数据库


Fungi_otu = read.delim(file = "zotutab_raw_ITS_2009to2020withreseq0608.txt",
                       sep="\t",row.names = 1,header=T,check.names=F,na.strings = "")
Fungi_taxon=read.delim(file = "taxonomy_ITS.txt",
                       sep="\t",row.names = 1,header=T,check.names=F,na.strings = "")

Fungi_taxon1=Fungi_taxon[match(row.names(Fungi_otu),row.names(Fungi_taxon)),]
otu=cbind(Fungi_otu,Fungi_taxon1[,1])
names(otu)[553] <- c('taxonomy')


otu=read.delim(file = "/Users/data from Xue/soil-borne pathogens/metageonomices/all_fungal_contigs.taxonomy copy.txt",
                       sep="\t",row.names = 1,header=T,check.names=F,na.strings = "")


otu=read.table(pipe("pbpaste"),header=T,row.names=1)
fung <- readRDS('/Users/data from Xue/soil-borne pathogens/funguild.rds')  #加载数据库

otu=read.delim(file = "/Users/data from Xue/soil-borne pathogens/raw data/taxonomy_fungi_18S.txt",sep="\t",header=T,check.names=F,na.strings = "")
fung_guilds <- funguild_assign(otu, db = fung, tax_col = 'taxonomy')

write.table(fung_guilds, 'all_fungal_contigs.taxonomy_funguilds_Eukdetect.txt', sep = '\t', row.names = F, quote = FALSE)

write.table(fung_guilds, 'all_fungal_taxonomy_funguilds_18S.txt', sep = '\t', row.names = F, quote = FALSE)

###
min(colSums(Fungi_otu)) 
Fungi_rare = as.data.frame(t(rrarefy(t(Fungi_otu), min(colSums(Fungi_otu)))))
Fungi_rare = Fungi_rare[which(rowSums(Fungi_rare) != 0 ),]
colSums(Fungi_rare)
write.table(Fungi_rare, "even_selected_fungi_otu.txt", row.names=T, col.names=T,sep="\t",
            quote = FALSE)

Fungi_otu = Fungi_rare


######################
soilbore = read.delim(file = "selected_fungi_otu_funguilds_alltreatment.txt",
                      sep="\t",row.names = 1,header=T,check.names=F)
treatment=read.table("clipboard",header=T,row.names=1)   #uparse
# treatment=read.table(pipe("pbpaste"),header=T,row.names=1)  

warm_precip=treatment[!(treatment$Clip %in% c("Cliping")), ]
warm_only=warm_precip[(warm_precip$Precipitation %in% c("normal")), ]

soilbore1=t(soilbore)
soilbore_warm=soilbore1[match(row.names(warm_only),row.names(soilbore1)),]
soilbore_warm1=t(soilbore_warm)
soilbore_warm1= soilbore_warm1[which(rowSums(soilbore_warm1) != 0 ),]
min(colSums(soilbore_warm1)) 

Fungi_rare = as.data.frame(t(rrarefy(t(soilbore_warm1), min(colSums(soilbore_warm1)))))
Fungi_rare = Fungi_rare[which(rowSums(Fungi_rare) != 0 ),]
colSums(Fungi_rare)
write.table(Fungi_rare, "even_selected_fungi_soilborne_otu_onlywarm.txt", row.names=T, col.names=T,sep="\t",
            quote = FALSE)


Fungi_taxon=read.delim(file = "classification.txt",
                       sep="\t",row.names = 1,header=T,check.names=F,na.strings = "")
Fungi_rare=read.delim(file = "even_selected_fungi_soilborne_otu_onlywarm.txt",
                       sep="\t",row.names = 1,header=T,check.names=F,na.strings = "")

Fungi_taxon1=Fungi_taxon[match(row.names(Fungi_rare),row.names(Fungi_taxon)),]
write.table(Fungi_taxon1, "even_selected_fungi_soilborne_otu_taxon_onlywarm.txt", row.names=T, col.names=T,sep="\t",
            quote = FALSE)


#######
######  match OTUid with representative fasta reads
otu_table <- read.delim(file = "even_selected_fungi_soilborne_otu_onlywarm.txt",
                        sep="\t",header=T,check.names=F,na.strings = "")
names(otu_table)[1]="taxonID"

sequences <- readDNAStringSet("zotus_ITS_2009to2020withreseq0608.fasta")  # 读取fasta格式的序列文件

selected_otuid <- unique(otu_table$taxonID)
selected_otu_table <- otu_table %>% filter(taxonID %in% selected_otuid)  # 使用dplyr筛选出选定的OTU
taxonid_to_seqname <- setNames(strsplit(names(sequences), " "), names(sequences))
selected_sequences <- DNAStringSet() 

for (taxonID in selected_otuid) {
  seq_name <- taxonid_to_seqname[[taxonID]]
  otu_sequence <- sequences[seq_name]
  selected_sequences <- append(selected_sequences, otu_sequence)
}

writeXStringSet(selected_sequences, "selected_sequences_soilborne_onlywarm.fasta") # 将选定的序列写入fasta文件


##################### two-sided repeated measures ANOVA 
res.aov <- anova_test(
  data = data, dv = RA, wid = block,
  within = c(Warm,year1))

get_anova_table(res.aov)










setwd("/Users/data from Xue/soil-borne pathogens")


#####################
library(ggplot2)
library(gghalves)
data=read.delim(file = "data_analysis_pathogen.txt",sep="\t",row.names = 1,header=T,check.names=F)
data$RAbiomass=data[,4]*data[,5]

dat2 <- summarySE(data, measurevar = "richness", groupvar = c("Warm","block"))
dat3 <- summarySE(data, measurevar = "FlC3", groupvar = c("Warm","block"))

dat2$richness1=dat2$sd/dat2$richness
dat3$c3=dat3$sd/dat3$FlC3

data4=cbind(dat2[,c(1,2,8)],dat3[,8])


colnames(data4)[4] <- "stability"
anova_model <- aov(stability ~ Warm, data = data4)
summary(anova_model)


cor.test(data4$richness1,data4$stability,method = "pearson")

ggplot(data4, aes(x = richness1, y = stability)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  theme_classic() +
  labs(x = "Richness",y = "Stability")







richness = ggplot(data4,aes(x=richness1,y=data4[,4],  colour= Warm))+
  geom_point(size=5,alpha=0.8,aes(shape=Warm,colour=Warm))+ 
  #scale_shape_manual(values = c(16,1))+   
  theme_bw()+   
  theme(legend.position="none")+
  geom_smooth(method = "lm",formula =y ~ poly(x, 2))+
  stat_poly_eq(aes(label = paste( ..rr.label.., stat(p.value.label), sep = '~`,`~')),
               formula = y ~ poly(x, 2), na.rm=T,parse = TRUE, label.x.npc = 0,  label.y.npc = 800, size = 5)

ggsave("Figure relationship plant biomass and pathogen RA.pdf",richness,width=8,height=6,units="cm") 

ggsave("Figure dispersion over year.pdf",richness,width=8,height=6,units="cm") 

RA <- ggplot(dat2,aes(x=year1,y=stability,  colour= Warm))+
  geom_point(size=5,alpha=0.8,aes(shape=Warm,colour=Warm))+ scale_color_viridis_d()+
  #scale_shape_manual(values = c(16,1))+   
  theme_bw()+   
  theme(legend.position="none")+
  geom_smooth(method = "lm",formula =y ~ poly(x, 1))+
  stat_poly_eq(aes(label = paste( ..rr.label.., stat(p.value.label), sep = '~`,`~')),
               formula = y ~ poly(x, 1), na.rm=T,parse = TRUE, label.x.npc = 0,  label.y.npc = 800, size = 5)

ggsave("Figure relationship plant richness and pathogen .pdf",RA,width=8,height=6,units="cm")


data2016=data[!(data$year %in% c("2015","2016","2017","2018","2019","2020")), ]


env_c=read.delim('clipboard',row.names=1, head=T)  
env_w=read.delim('clipboard',row.names=1, head=T)  
cohend(env_c$RA,env_w$RA)



fm1=lmer(scale(FlC3) ~ scale(Alternaria) + (1|year)+ (1|block),data=data) 
summary(fm1)
AIC(fm1)

fm1=lmer(scale(Richenss_Euk) ~ warm + (1|year)+ (1|block),data=data) 
summary(fm1)
AIC(fm1)

data_c=data[(data$Warm %in% c("unwarming")), ]
data_w=data[(data$Warm %in% c("warming")), ]

fm1=lmer(scale(ratio) ~ scale(year1) +  (1|block),data=data_c) 
summary(fm1)
AIC(fm1)

fm1=lmer(scale(ratio) ~ scale(year1) +  (1|block),data=data_w) 
summary(fm1)
AIC(fm1)

car::Anova(fm1,type=2)
r.squaredGLMM(fm1)

cohend(data_c$RA,data_w$RA)

fm1=lmer(scale(Totfungi) ~ Warm +  (1|year)+(1|block),data=data) 
summary(fm1)
AIC(fm1)


###############  effect size
dat1=read.table("clipboard",header=T,row.names=1)   #uparse
#dat1=read.table(pipe("pbpaste"),header=T,row.names=1)  

dat1$index<- factor(dat1$index, levels=c("richness","pd","shannon","pielou","simpson","GeoChip","Metagenomic"), ordered=TRUE)

effect.size = ggplot(data = dat1, aes(x = index, y =effect,fill = index))+ ylim(-0.9,0.2)+
  scale_color_viridis_d()+ 
  geom_bar(stat ="identity",width = 0.6,position = position_dodge())+
  geom_errorbar(aes(ymin=effect-se, ymax=effect+se),  size=0.2, 
                width=0.2,position = position_dodge(0.5))+
  theme_bw()+theme(legend.position="none")+xlab("Multiple facets of soil-borne pathogens diversity")+ylab("Effect size")
ggsave("Figure effect size diversity_pathogen1.pdf",effect.size ,width=8,height=6,units="cm") 




############### NMDS
nmds= ggplot (data, aes(x = data[,55], y = data[,56],color = as.factor(year)))+ #divided by site and treatment
  geom_point(size=1,alpha=1, aes(shape=Warm)) + scale_color_viridis_d()+
  scale_shape_manual(values = c(16, 17, 15, 3, 12, 8))+
  #scale_shape_manual(values = c(15, 16, 17, 18, 19, 21, 22, 23, 24, 25))+
  xlab(paste("NMDS1")) +ylab(paste("NMDS2"))+
  theme_bw()+
  theme(legend.position="none")+
  #annotate('text', label = 'unifrac.Wt', x =0, y =-0.3, size =5.5)+
  theme(axis.text.x = element_text(size = 12),  # 设置 x 轴坐标字体大小
        axis.text.y = element_text(size = 12),  # 设置 y 轴坐标字体大小
        axis.title.x = element_text(size = 14),  # 设置 x 轴标题字体大小
        axis.title.y = element_text(size = 14))  # 设置 y 轴标题字体大小

ggsave("Figure community_sorense.pdf",nmds,width=8,height=6,units="cm") 





############################### Fungal composition phylum
##############################
otu = read.delim(file = "even_selected_fungi_plant_pathogen_onlywarm.txt",row.names = 1,header=T,check.names=F)
taxon= read.delim(file = "even_selected_fungi_plant_pathogen_taxon_onlywarm.txt",row.names = 1,header=T,check.names=F)

taxon1=taxon[match(row.names(otu),row.names(taxon)),]

otu_abun <- otu/colSums(otu)
p<- aggregate(otu_abun, by=list(taxon1$genus), sum) #门水平
rownames(p) <- p$Group.1
p <- p[,-1]
top <- names(head(sort(rowSums(p), decreasing = T), 65)) #基于相对丰度进行排序
####phylum level
top_10<- c("Phaeosphaeria", "Alternaria", "Setophoma" ,"Bipolaris","Pseudopithomyces" ,"Paraphoma", "Marasmius","Rhizopus","Neosetophoma")
####class level
#top_5<- c("Alphaproteobacteria", "Thermoleophilia", "Actinobacteria", "Bacilli" , "Gammaproteobacteria")

taxon1$genus <- as.character(taxon1$genus)
taxon1$genus[!(taxon1$genus)%in%top_10] <- "Others"
p_top <- aggregate(otu_abun, by=list(taxon1$genus), sum) #获取top10
rownames(p_top) <- p_top[,1]
p_top <- p_top[,-1]
p_top <- p_top[order(rowSums(p_top)),] #按每行和排序
q_top <- t(p_top)
q_top <- as.data.frame(q_top)
q_top$sample <- rownames(q_top)
q <- melt(q_top,ID="names")
colnames(q)[names(q)=="variable"]<-"Taxa"

#write.table(p_top, "fungal_p_topclasf_taxon_RA_1.txt", sep="\t",row.names=T, quote = FALSE)
#phylum=p_top
#class=p_top

#datcom=rbind(phylum,class[c(1,5),])
#write.table(datcom, "Fungi_relative.abundance_phylum.txt", sep="\t",row.names=T, quote = FALSE)
############################
###########################
colors<-c("darkolivegreen3","gold","dodgerblue4","darkseagreen",
          "chartreuse4","darkorange","burlywood2","brown3","#984EA3","cyan3","steelblue")
#colors<-c("grey50","darkolivegreen3","gold","dodgerblue4","darkseagreen",
#         "chartreuse4","darkorange","burlywood2","brown3","#984EA3","cyan3","steelblue")

composition=ggplot(dat1, aes( x = treat, y = mean*100, fill = index))+
  geom_bar(position = "fill", stat = "identity")+
  theme_bw()+
  #scale_fill_manual(values=colors)+ 
  scale_y_continuous(expand = c(0,0))+
  labs(x="",y="Relative Abundance (%)",fill="Genera")+
  theme(text=element_text(size=12),
        axis.text.y=element_text(size=12,color = "black"),
        axis.text.x=element_text(size=12,color = "black",angle = 45, hjust = 0.5, vjust = 0.5),
        legend.title=element_text(size=12), 
        legend.text=element_text(size=12))+ 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(color = 'black', fill = 'transparent')) +
  guides(fill=guide_legend(keywidth = 1.3, keyheight = 2))

ggsave("Figure Fungal phylum composition.pdf",composition,width=12,height=8,units="cm")



###############  soil-borne plant pathogens (%)
dat1=read.table("clipboard",header=T,row.names=1)   #uparse
#dat1=read.table(pipe("pbpaste"),header=T,row.names=1)  

dat1$index<- factor(dat1$index, levels=c("Phaeosphaeria","Alternaria","Setophoma","Bipolaris","Pseudopithomyces","Paraphoma","Marasmius","Rhizopus","Neosetophoma","others"), ordered=TRUE)

percentage = ggplot(data = dat1, aes(x = index, y =mean,fill = index))+ ylim(0,25)+
  scale_color_viridis_c()+ 
  geom_bar(stat ="identity",width = 0.6,position = position_dodge())+
  theme_bw()+theme(legend.position="none")+xlab("")+ylab("Soil-borne pathogens (%)")+
  theme(axis.text.x=element_text(size=12,color = "black",angle = 45, hjust = 0.5, vjust = 0.5))

ggsave("Figure Soil-borne pathogens percentage_pathogen.pdf",percentage ,width=8,height=8,units="cm") 




fm1=lmer(scale(Neosetophoma) ~ Warm + (1|year)+ (1|block),data=data) 
summary(fm1)
AIC(fm1)

car::Anova(fm1,type=2)
r.squaredGLMM(fm1)



library(Hmisc) 
#注意要将数据框转换为矩阵
cor.mat <- rcorr(as.matrix(dat[,c(18:39)]), type = "pearson")
cor.mat$r
cor.mat$P


########################
######################## correlating diversity with functions_LMM from Linwei
# function to get the R and p value from the linear mixed model
corenvs<-function(divtest,envs){
  sapply(colnames(envs),function(x){
    message(x)
    if(sum(abs(divtest-envs[,x]),na.rm = T)<0.01 | length(unique(envs[,x])) < 7 | length(unique(divtest)) < 7) { 
      result=list(r=NA,pvalue=NA)
    }else{
      div<-data.frame(divtest=divtest,env=envs[,x],dat)
      div<-div[(!is.na(div$divtest)) & (!is.na(div$env)),]
      if (length(unique(div$year))<2) {
        #fm1<-lmer(env~divtest+(1|block),data=div)
      } else {
        fm1<-lmer(env~divtest+(1|block),data=div)
      }
      presult<-car::Anova(fm1,type=2)
      coefs<-coef(summary(fm1))[ , "Estimate"]
      pvalue=presult[,3]
      r2<-r.squaredGLMM(fm1)
      r=ifelse(coefs["divtest"]>0,(r2[1,1])^0.5,-(r2[1,1])^0.5)
      names(r)<-NULL
      result=list(r=r,pvalue=pvalue)
    }
    result
  },simplify = F)
}

dat=read.delim(file = "data_analysis_pathogen.txt",sep="\t",row.names = 1,header=T,check.names=F)
divs<-dat[,c(43,44)]  ##richness of bacterial and fungal groups
envs<-dat[,c(29:32,37:42)]   ##ecosystem functions



data_c=dat[(dat$Warm %in% c("unwarming")), ]
data_w=dat[(dat$Warm %in% c("warming")), ]

divs<-data_c[,c(42,3)]  ##richness of bacterial and fungal groups
envs<-data_c[,c(20:41)]   ##ecosystem functions

divs<-data_w[,c(42,3)]  ##richness of bacterial and fungal groups
envs<-data_w[,c(20:41)]   ##ecosystem functions



library(lme4)
library(MuMIn)
divs1<-sapply(colnames(divs),function(y){
  message(y)
  divtest=divs[,y]
  corenvs(divtest=divtest,envs = envs)
},simplify = F)

test<-sapply(divs1, function(x){
  sapply(x, function(y){
    y$r
  })
})

#write.csv(t(test),"LMM_R_richness and environemntal_variables_all_RA.csv")

test2<-sapply(divs1, function(x){
  sapply(x, function(y){
    y$pvalue
  })
})

#write.csv(t(test2),"LMM_Pvalue_richness and environemntal_variables_warm_RA.csv")

##ADJUST P

pvalues<-t(test2)

library(Hmisc)
adjust.p<-p.adjust (as.matrix(pvalues), method="fdr")
ad.pmatrix<-matrix(adjust.p,nrow=nrow(pvalues))
row.names(ad.pmatrix)=row.names(pvalues)
colnames(ad.pmatrix)=colnames(pvalues)
#write.csv(ad.pmatrix,"LMM_adjust P_taxa richness and functions.csv")

library(corrplot)
corrplot(as.matrix(t(test)),is.corr = FALSE,win.asp = .5, method= "square",tl.col="black",tl.cex = 0.6,tl.srt=45,p.mat =as.matrix(ad.pmatrix),sig.level = c(.001, .01, .05), pch.cex = .8,
         insig = "label_sig", pch.col = "white",cl.cex=0.6,cl.ratio=0.35,cl.align.text="l",cl.offset=0.2)

library(corrplot)

heatmap= corrplot(as.matrix(t(test)),is.corr = FALSE,win.asp = .5, method= "square",tl.col="black",tl.cex = 0.6,tl.srt=45,p.mat =as.matrix(pvalues),sig.level = c(.001, .01, .05,0.1), pch.cex = .8,
                  insig = "label_sig", pch.col = "white",cl.cex=0.6,cl.ratio=0.35,cl.align.text="l",cl.offset=0.2)

setwd("C:/Users/Jiajie/Desktop")
ggsave("Figure nematode trophic group relative abundance.pdf",heatmap,width=18,height=13,units="cm")








################################
################################
library(linkET)
library(dplyr)
library(ggplot2)
library(cols4all)
library(qgraph)
library(corrplot)


envs=read.table(pipe("pbpaste"),head=T,row.names=1,sep = "\t")  
lmm=read.table(pipe("pbpaste"),head=T,row.names=1,sep = "\t")  

set_corrplot_style()
mycols=c( "#d8b365", "#5ab4ac", "#01665e", "skyblue")
qcorrplot(correlate(envs), type = "upper", diag = FALSE) +  
  geom_square(alpha = 1,color = "white") + 
  scale_fill_gradient(high="#01665e", low="skyblue") +
  geom_mark(size = 3,only_mark= T,sig_level= c(0.05, 0.01, 0.001), sig_thres= 0.05, colour= 'white')+
  geom_couple(aes(colour = cd, size = rd,linetype=pd), data = lmm , curvature = nice_curvature()) +
  scale_size_manual(values = c(0.2, 0.4, 0.6,0.8,1.2)) + #连线粗细 
  scale_colour_manual(values = mycols) + #连线配色 
  guides(size = guide_legend(title = "LMM's r",override.aes = list(colour = "grey35"),order = 2),
         colour = guide_legend(title = "LMM's p", override.aes = list(size = 3),order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))



###demo example data
{
  data("varechem", package = "vegan")
  data("varespec", package = "vegan")
  
  mantel <- mantel_test(varespec, varechem, spec_select = list(Spec01 = 1:7,Spec02 = 8:18, Spec03 = 19:37,Spec04 = 38:44)) %>% mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.4, Inf),                                      labels = c("< 0.2", "0.2 - 0.4", ">= 0.4")), pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
                                                                                                                                                                                                                                                                              labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
  
  #write.csv(mantel, 'mental&Correlation1.csv', quote = FALSE)
  
  ## 现学现用
  set_corrplot_style()
  qcorrplot(correlate(varechem), type = "upper", diag = FALSE) +
    geom_square() +
    geom_couple(aes(colour = pd, size = rd), data = mantel, curvature = nice_curvature()) +
    scale_size_manual(values = c(0.5, 1, 2)) +
    scale_colour_manual(values = color_pal(3)) +
    guides(size = guide_legend(title = "Mantel's r",override.aes = list(colour = "grey35"),order = 2),
           colour = guide_legend(title = "Mantel's p", override.aes = list(size = 3),order = 1),
           fill = guide_colorbar(title = "Pearson's r", order = 3))
  }





#############################SEM
library(lavaan)
library(piecewiseSEM)
library(nlme)
library(lme4)
data=read.delim(file = "data_analysis_pathogen.txt",sep="\t",row.names = 1,header=T,check.names=F)


#final model
model4<-psem(
  lmer(richness ~ warm + pH  + annual_moisture + FlTotl+(1|year) + (1|block),data=data),
  
  lmer(temperature_annual ~ warm  + (1|year)+ (1|block),data=data),
  lmer(annual_moisture ~ warm  +temperature_annual +  (1|year)+ (1|block),data=data),
  lmer(pH ~   warm+NO3 + (1|year)+ (1|block),data=data),
  lmer(NO3 ~   temperature_annual+ (1|year)+ (1|block),data=data),
  #lmer(plant_richness ~  annual_moisture + (1|year)+ (1|block),data=data),
  lmer(FlTotl ~  warm+annual_moisture   + (1|year)+ (1|block),data=data))

summary(model4, .progressBar = F)
plot(model4)


#final model  # remove plant
model4<-psem(
  lmer(richness ~ warm + pH  + annual_moisture + FlTotl+(1|year) + (1|block),data=data),
  
  lmer(temperature_annual ~ warm  + (1|year)+ (1|block),data=data),
  lmer(annual_moisture ~ warm  +temperature_annual +  (1|year)+ (1|block),data=data),
  lmer(pH ~   warm+NO3 + (1|year)+ (1|block),data=data),
  lmer(NO3 ~   temperature_annual+ (1|year)+ (1|block),data=data),
  #lmer(plant_richness ~  annual_moisture + (1|year)+ (1|block),data=data),
  lmer(FlTotl ~  warm+annual_moisture   + (1|year)+ (1|block),data=data))

summary(model4, .progressBar = F)
plot(model4)


### 计算直接效应 和 间接效应 
#提取计算数据
kk1 <- summary(model4, .progressBar = F)
kk2 <- kk1$coefficients[,-(3:6)]
kk3 <- as.data.frame(kk2)
kk4 <- kk3[!grepl("~~", kk3$Predictor), ]
kk4

library(dplyr)
# 结果整理为便于计算的数据框
result1 <- kk4 %>%
  dplyr::relocate(
    from   = Predictor,
    to     = Response,
    weight = Std.Estimate,
    p      = P.Value
  )
#查看整理后结果
result1  

#开始计算
#===========================直接效应============================
# 示例数据
data <-result1 

# 定义计算直接效应的函数
calculateDirectEffects <- function(data, factors) {
  # 初始化直接效应向量
  direct_effects <- numeric(nrow(data))
  
  # 遍历数据集中的每一行
  for (i in 1:nrow(data)) {
    # 检查路径是否在因子列表中，如果是则直接效应即为weight值，否则为0
    if (paste(data$from[i], data$to[i], sep = "_") %in% factors) {
      direct_effects[i] <- data$weight[i]
    } else {
      direct_effects[i] <- 0
    }
  }
  
  # 过滤 Direct_Effect 等于 0 的行
  non_zero_rows <- direct_effects != 0
  direct_effects <- direct_effects[non_zero_rows]
  from_to <- paste(data$from[non_zero_rows], "→", data$to[non_zero_rows])
  
  # 创建包含直接效应的数据框
  direct_effects_df <- data.frame(
    from_to = from_to,
    Direct_Effect = direct_effects
  )
  
  return(direct_effects_df)
}


# 定义因子列表
# 自动化生成直接效应的因子列表的函数
generateDirectEffectFactors <- function(data) {
  # 从数据集中获取唯一的 from 和 to 组合，符合直接效应的条件
  unique_direct_combinations <- unique(paste(data$from, data$to, sep = "_"))
  
  # 提取起点和终点
  from_factors <- unique(data$from)
  to_factors <- unique(data$to)
  
  # 拆分因子并转换为向量
  direct_factors <- paste(rep(from_factors, each = length(to_factors)), rep(to_factors, length(from_factors)), sep = "_")
  
  return(direct_factors)
}

# 调用函数生成直接效应的因子列表
direct_factors <- generateDirectEffectFactors(data)

# 调用函数计算直接效应
direct_effects_result <- calculateDirectEffects(data, direct_factors)

print(direct_effects_result)





#==========================================间接效应==========================================
# 示例数据
data <- result1  

# 从数据中生成路径列表
unique_from <- unique(data$from)
unique_to <- unique(data$to)
paths <- expand.grid(from = unique_from, to = unique_to)
paths <- split(paths, seq(nrow(paths)))
paths <- lapply(paths, function(x) c(as.character(x$from), as.character(x$to), "richness"))

# 定义计算间接效应的函数
calculateIndirectEffects <- function(data, from_factor, through_factor, to_factor) {
  through_weight <- data$weight[data$from == from_factor & data$to == through_factor]
  to_weight <- data$weight[data$from == through_factor & data$to == to_factor]
  
  if (length(through_weight) == 0 | length(to_weight) == 0) {
    message(paste("Skipping invalid path:", from_factor, "→", through_factor, "→", to_factor))
    return(NULL)
  }
  
  indirect_effect <- through_weight * to_weight
  
  return(indirect_effect)
}

# 创建空的数据框来存储间接效应结果
indirect_effects_df <- data.frame(from_to = character(), Indirect_Effect = numeric())

# 循环计算间接效应并存储结果
for (path in paths) {
  from_factor <- path[1]
  through_factor <- path[2]
  to_factor <- path[3]
  
  indirect_effect <- calculateIndirectEffects(data, from_factor, through_factor, to_factor)
  if (!is.null(indirect_effect)) {
    from_to <- paste(from_factor, through_factor, to_factor, sep = " → ")
    indirect_effects_df <- rbind(indirect_effects_df, data.frame(from_to = from_to, Indirect_Effect = indirect_effect))
  }
}

# 打印结果
print(indirect_effects_df)





#然后将间接效应的和求出来
# 创建示例数据框
indirect_effects_df <- indirect_effects_df

# 定义自动化函数
calculate_total_indirect_effect <- function(data) {
  # 提取开头和结尾相同的值，并对其"Indirect_Effect"列求和
  total_indirect_effect <- data %>%
    mutate(
      start_pattern = sub(" → .*", "", from_to),
      end_pattern = sub(".* → ", "", from_to)
    ) %>%
    group_by(start_pattern, end_pattern) %>%
    summarise(total_indirect_effect = sum(Indirect_Effect), .groups = "drop") %>%
    ungroup() %>%
    arrange(start_pattern, end_pattern) %>%
    mutate(from_to = paste0(start_pattern, " → ", end_pattern)) %>%
    select(from_to, total_indirect_effect)
  
  return(total_indirect_effect)
}

# 调用自动化函数并打印结果
result_total_indirect_effect <- calculate_total_indirect_effect(indirect_effects_df)
print(result_total_indirect_effect)




#计算总效应 
#提取三个主要结果
print(direct_effects_result)
print(indirect_effects_df)
print(result_total_indirect_effect)

#合并三个结果
# 直接效应数据框
direct_effects_result <- direct_effects_result

# 间接效应数据框
indirect_effects_df <- indirect_effects_df
# 总间接效应数据框
total_indirect_effect_df <- as.data.frame(result_total_indirect_effect)

# 合并直接效应和间接效应数据框
total_effects_df <- bind_rows(direct_effects_result, indirect_effects_df, total_indirect_effect_df)

# 打印结果
print(total_effects_df)



# 创建数据框
data <- total_effects_df

# 使用dplyr包中的group_by和summarise函数，将相同因子的数据转移到相同行并求和
data_processed <- data %>%
  group_by(from_to) %>%
  summarise(Direct_Effect = sum(Direct_Effect, na.rm = TRUE),
            Indirect_Effect = sum(Indirect_Effect, na.rm = TRUE),
            total_indirect_effect = sum(total_indirect_effect, na.rm = TRUE))

# 打印处理后的数据框
print(data_processed)



# 使用dplyr包中的mutate函数，将Direct_Effect和total_indirect_effect列求和
data_processed2 <- data_processed %>%
  mutate(Total_Effect = Direct_Effect + total_indirect_effect)

# 打印处理后的数据框
print(data_processed2)

# 过滤不需要的路径，我们的目标因子是Richness，所以我们不需要那些中间路径的值
df_filtered1 <- data_processed2 %>%
  filter(grepl("richness$", from_to))

# 输出结果
df_filtered1

# 将其转换为数据框
data_frame <- as.data.frame(df_filtered1)

# 使用dplyr的过滤函数来保留只有一条箭头路径的数据
filtered_data1 <- data_frame %>%
  filter(Indirect_Effect == 0)

# 打印筛选后的数据
print(filtered_data1)




#整理绘图数据
filtered_data2 <- filtered_data1[,-(3)]
filtered_data2


# 假设您的数据框名为data_frame，将其转换为数据框（如果不是数据框的话）
data_frame <- as.data.frame(filtered_data2)

# 使用正则表达式提取箭头前的字母部分，并创建新列from_to_clean
data_frame$from_to_clean <- sub(" .*", "", data_frame$from_to)

# 打印处理后的数据框，只包含箭头前的字母部分
kk4 <- print(data_frame[, c("from_to_clean", "Direct_Effect", "total_indirect_effect", "Total_Effect")])
kk4 


###  plot 
data=read.delim('clipboard',row.names=1, head=T)  #
#  data=read.table(pipe("pbpaste"),header=T,row.names=1)  
# 使用ggplot2绘制柱形图
# 定义X轴的顺序
data$from<- factor(data$from_to_clean,levels = c("warm","temperature_annual","annual_moisture","pH","NO3.N","FlTotl"))
# 设置因子变量，并指定顺序
data$from <- factor(data$from,levels = c("warm","temperature_annual","annual_moisture","pH","NO3","FlTotl"))

library(ggplot2)

abc<- ggplot(data, aes(x= from, y= Total_Effect, fill= from)) +
  geom_bar(stat="identity", position=position_dodge(), color="white", width=0.8) +
  theme_bw()+ ylim(-0.4,0.4)+
  #scale_fill_manual(values=c("#5ab4ac", "#d8b365"))+
  scale_fill_manual(values = c("#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e", "skyblue"))+
  labs( y = 'Standardized effect size based on SEM')+
  theme(legend.justification=c(1,1),legend.position="none")+
  theme(panel.grid = element_blank(),panel.background = element_rect(color='black',fill='transparent'))+
  theme(axis.text.x=element_text(size=12,color = "black",angle = 30, hjust = 0.5, vjust = 0.5))

ggsave("Figure Standardized effect size based on SEM.pdf",abc,width=8,height=8,units="cm")


# setwd("J:/Nematode analysis_2021/SEM modeling")









data_c=data[(data$Warm %in% c("unwarming")), ]
data_w=data[(data$Warm %in% c("warming")), ]


library(simba)

#比较回归系数，p<0.05 表示二者存在区别，反之无差别，详情 ?diffslope
set.seed(123)
diffslope(data_c$richness, data_c$plant_richness, data_w$richness, data_w$plant_richness)

#比较截距，p<0.05 表示二者存在区别，反之无差别，详情 ?diffic
set.seed(123)
diffslope(data_c$richness, data_c$FlTotl, data_w$richness, data_w$FlTotl)



set.seed(123)
diffslope(data_c$RA, data_c$FlTotl, data_w$RA, data_w$FlTotl)


RICH=lm(formula = FlTotl ~ RA , data = data_c)
summary(RICH)

RICH=lm(formula = FlTotl ~ RA , data = data_w)
summary(RICH)

slope <- coef(RICH)[2]


RICH=lm(formula = FlTotl ~ richness , data = data_c)
summary(RICH)

RICH=lm(formula = FlTotl ~ richness , data = data_w)
summary(RICH)




RICH=lm(formula = plant_richness ~ richness , data = data_c)
summary(RICH)

RICH=lm(formula = plant_richness ~ richness , data = data_w)
summary(RICH)




data=read.delim(file = "data_analysis_pathogen.txt",sep="\t",row.names = 1,header=T,check.names=F)
fungi=summarySE(data,measurevar="richness",groupvar=c("Warm","year"))
fungi


c=read.table("clipboard",head=T, row.names=1)

w=read.table("clipboard",head=T, row.names=1)

cohend(c$annual_moisture,w$annual_moisture,paired=TRUE)




###############  effect size
dat1=read.table("clipboard",header=T,row.names=1)   #uparse
#dat1=read.table(pipe("pbpaste"),header=T,row.names=1)  

dat1$index<- factor(dat1$index, levels=c("temperature_annual","annual_moisture","FlTotl","FlC3","FlC4","plant_richness","RA","richness"), ordered=TRUE)

value = ggplot(data = dat1, aes(x = treatment, y =ava,color = treatment))+ 
  scale_color_viridis_d()+ 
  geom_bar(aes(fill = treatment),stat ="identity",width = 0.6,position = position_dodge())+
  scale_fill_viridis_d()+ 
  geom_errorbar(aes(ymin=ava-se, ymax=ava+se),  size=0.2, 
                width=0.2,position = position_dodge(0.5))+
  facet_wrap(index~.,scales="free",ncol=4)+
  theme_bw()+theme(legend.position="none")+ylab("value")

ggsave("Figure data_value.pdf",value ,width=15,height=9,units="cm") 




library(ggplot2)

value= ggplot(dat1, aes(x=phylum, y=ra100))+
  geom_col()+
  geom_bar(stat="identity", fill="steelblue")+
  theme_bw()+theme(legend.position="none")+ylab("value")+
  theme(axis.text.x = element_text(angle=60, hjust = 1))

ggsave("Figure data_value-1.pdf",value ,width=7,height=6,units="cm") 






bacteria=summarySE(data,measurevar="RA",groupvar=c("year","Warm"))

bacteria$stability=bacteria$sd/bacteria$RA




bacteria1=summarySE(bacteria,measurevar="stability",groupvar=c("Warm"))
bacteria1




c=read.table("clipboard",head=T, row.names=1)







setwd("/Users/data from Xue/soil-borne pathogens")
################ CCA 
otu <- read.delim('even_nematodes.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
sampledata <- data.frame(otu)
env <- read.delim('cca_env.txt', row.names = 1, sep = '\t',     stringsAsFactors  =FALSE, check.names = FALSE)
group<- read.delim('cca_group.txt', row.names = 1, sep = '\t')


sampledata <- t(sampledata)
#??OTU???ݽ???hellingerת??
sampledata <- decostand(sampledata,method = "hellinger")
group <- as.list(data$Warm)
#????????????????ɫ
col <- c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")
#?Ƚ???DCA????
dca <- decorana(veg = sampledata)
dca1 <- max(dca$rproj[,1])
dca2 <- max(dca$rproj[,2])
dca3 <- max(dca$rproj[,3])
dca4 <- max(dca$rproj[,4])
GL <- data.frame(DCA1 = c(dca1), DCA2 = c(dca2), DCA3 = c(dca3), DCA4 = c(dca4))
GL
rownames(GL) <- c("Gradient length")
write.csv(GL, file = "dca.csv")
#?ٽ???CCA????  
cca <- cca(sampledata, data[,21:32], scale = TRUE)
cca_time <- cca(sampledata, data[,21],data[,22:32], scale = TRUE)
cca_plant <- cca(sampledata, data[,29:32],data[,21:28], scale = TRUE)
cca_soil <- cca(sampledata, data[,22:28],data[,c(21,29:32)], scale = TRUE)




ccascore <- scores(cca)
ccascore$sites
cca$CCA$biplot
ccascore$species
#??ȡ??Ҫ??Ϣ
write.csv(ccascore$sites, file = "cca.sample.csv")
write.csv(cca$CCA$biplot, file = "cca.env.csv")
write.csv(ccascore$species, file = "cca.species.csv")
#Ϊ??ͼ׼??????
CCAE <- as.data.frame(cca$CCA$biplot[,1:2])
CCAE=CCAE*2

CCAS1 <- ccascore$sites[,1]
CCAS2 <- ccascore$sites[,2]

plotdata <- data.frame(rownames(ccascore$sites), CCAS1, CCAS2, data$year,data$Warm)
colnames(plotdata) <- c("sample","CCAS1","CCAS2","year","Warm")

cca1 <- round(cca$CCA$eig[1]/sum(cca$CCA$eig)*100,2)
cca2 <- round(cca$CCA$eig[2]/sum(cca$CCA$eig)*100,2)



#????CCAͼ
P <- ggplot(plotdata, aes(CCAS1, CCAS2)) +scale_color_viridis_d()+
  geom_point(size=2, aes(shape=data$Warm, color = factor(data$year))) +
  xlab(paste("CCA1 ( ",cca1,"%"," )", sep = "")) + 
  ylab(paste("CCA2 ( ",cca2,"%"," )", sep = "")) +
  geom_segment(data = CCAE, aes(x = 0, y = 0, xend = CCAE[,1], yend = CCAE[,2]),
               colour = "black", size = 0.8,
               arrow = arrow(angle = 30, length = unit(0.4, "cm"))) +
  geom_text_repel(data = CCAE, segment.colour = "black",
                  aes(x = CCAE[,1], y = CCAE[,2], 
                      label = rownames(CCAE)),size=8) +
  geom_vline(aes(xintercept = 0), linetype = "dotted") +
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid = element_blank(),
        axis.title = element_text(color = "black", size = 15),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(colour = "black", size = 15),
        axis.title.y = element_text(colour="black", size = 15),
        axis.text = element_text(colour = "black", size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 15), legend.key = element_blank(),
        plot.title = element_text(size = 15, colour = "black", 
                                  face = "bold", hjust = 0.5))

P

ggsave("Figure CCA.pdf", P, width =16, height = 12, units="cm")


envfit(cca,data[,21:32],permu=999)

permutest(cca,permu=999)

summary(cca)





######################### temporal scaling 
######################### species turnover
otu = read.delim(file = "even_selected_fungi_plant_pathogen_onlywarm.txt",row.names = 1,header=T,check.names=F)
data=read.delim(file = "data_analysis_pathogen.txt",sep="\t",row.names = 1,header=T,check.names=F)
otu=t(otu)

#B1
select_time=subset(data, Warm %in% c("unwarming"))
select_times=subset(select_time, block %in% c("B1"))
times = select_times[, 11, drop = FALSE]  
rownames(times) = rownames(select_times) 
filtered_OTU = otu[match(rownames(times), rownames(otu)), ]
comm=filtered_OTU
grp=times
str.test=STR(comm=comm,grp=times,fit.method = "OLS",mean.method="am",STR.model = c("power", "log"))
setwd("C:/Users/Jiajie/Desktop")
save.file(str.test$details,prefix = "NWS",filename = "STR.Details_B1")

#B2
select_time=subset(data, Warm %in% c("unwarming"))
select_times=subset(select_time, block %in% c("B2"))
times = select_times[, 11, drop = FALSE]  
rownames(times) = rownames(select_times) 
filtered_OTU = otu[match(rownames(times), rownames(otu)), ]
comm=filtered_OTU
grp=times
str.test=STR(comm=comm,grp=times,fit.method = "OLS",mean.method="am",STR.model = c("power", "log"))
setwd("C:/Users/Jiajie/Desktop")
save.file(str.test$details,prefix = "NWS",filename = "STR.Details_B2")

#B3
select_time=subset(data, Warm %in% c("unwarming"))
select_times=subset(select_time, block %in% c("B3"))
times = select_times[, 11, drop = FALSE]  
rownames(times) = rownames(select_times) 
filtered_OTU = otu[match(rownames(times), rownames(otu)), ]
comm=filtered_OTU
grp=times
str.test=STR(comm=comm,grp=times,fit.method = "OLS",mean.method="am",STR.model = c("power", "log"))
setwd("C:/Users/Jiajie/Desktop")
save.file(str.test$details,prefix = "NWS",filename = "STR.Details_B3")

#B4
select_time=subset(data, Warm %in% c("unwarming"))
select_times=subset(select_time, block %in% c("B4"))
times = select_times[, 11, drop = FALSE]  
rownames(times) = rownames(select_times) 
filtered_OTU = otu[match(rownames(times), rownames(otu)), ]
comm=filtered_OTU
grp=times
str.test=STR(comm=comm,grp=times,fit.method = "OLS",mean.method="am",STR.model = c("power", "log"))
setwd("C:/Users/Jiajie/Desktop")
save.file(str.test$details,prefix = "NWS",filename = "STR.Details_B4")



turnover=read.table("clipboard",header=T,row.names=NULL)   #uparse
###########
fm1<-lmer(log(richness) ~ log(time.scale) + (1|Block),data=turnover)
summary(fm1)
AIC(fm1)
car::Anova(fm1,type=2)
r.squaredGLMM(fm1)

###################
###################
plot1=read.table(pipe("pbpaste"),header=T,row.names=1)  # like NEE paper

turnover= ggplot(data = plot1, aes(x = year, y = X16S,color = Treatment,group = Treatment))+
  scale_color_viridis_d() + ylim(0.5,4)+
  geom_point(size=1.6,alpha=0.8)+theme_bw()+ theme(legend.position="none")+
  geom_smooth(method="lm", formula = y ~ poly(x, 1),se = F, fullrange = F)+	
  stat_poly_eq(aes(label = paste( ..rr.label.., stat(p.value.label), sep = '~`,`~')),
               formula = y ~ poly(x, 1), parse = T)+
  xlab("ln [time(year)]")+	ylab("ln [species richness]")

ggsave("Figure phylogenetic diversity turnover.pdf",turnover,width=8,height=6,units="cm")


setwd("C:/Users/Jiajie/Desktop")	
ggsave("Figure species turnover and STR.pdf",turnover,width=8,height=6,units="cm")	







######################
######################
library(piecewiseSEM)
library(MASS)
library(nlme)
library(vegan)
        
cdata=read.table(pipe("pbpaste"),header=T,row.names=1)     # like NEE paper
wdata=read.table(pipe("pbpaste"),header=T,row.names=1)     # like NEE paper

#####################control treatment

mainplot=cdata$year

cmodlist = list(lme(X16S ~ year, random = ~ 1|block/mainplot,data=cdata))

cmod=rsquared(cmodlist)

lme.16S=lme(X16S ~ year, random = ~ 1|block/mainplot,data=cdata)
c16Sslope=summary(lme.16S)$coefficients$fix[2]

cmod$Marginal
cmod$Conditional


#####################warm treatment

mainplot=wdata$year

wmodlist = list(
  lme(X16S ~ year, random = ~ 1|block/mainplot,data=wdata))

wmod=rsquared(wmodlist)
lme.16S=lme(X16S ~ year, random = ~ 1|block/mainplot,data=wdata)

w16Sslope=summary(lme.16S)$coefficients$fix[2]

wmod$Marginal
wmod$Conditional


########################################## 1000 time permutations
c16S=cdata$X16S
cblock=cdata$block
cyear=cdata$year

w16S=wdata$X16S
wblock=wdata$block
wyear=wdata$year


perm.c16Sslope = c(rep(NA,1000))
perm.w16Sslope = c(rep(NA,1000))

for (i in 1:1000)
{
  control.16S=sample(c16S)
  cperm16S=lme(control.16S ~ cyear, random = ~ 1|cblock/mainplot)
  perm.c16Sslope[i]=summary(cperm16S)$coefficients$fix[2]
  
  
  warm.16S=sample(w16S)
  wperm16S=lme(warm.16S ~ wyear, random = ~ 1|wblock/mainplot)
  perm.w16Sslope[i]=summary(wperm16S)$coefficients$fix[2]
  
}


###################one sample permutation test
control16S.test=t.test(perm.c16Sslope,mu=c16Sslope,alternative="two.sided")
control16S.t=control16S.test$statistic
control16S.perm.p=control16S.test$p.value

warm16S.test=t.test(perm.w16Sslope,mu=w16Sslope,alternative="two.sided")
warm16S.t=warm16S.test$statistic
warm16S.perm.p=warm16S.test$p.value


#####################control vs. warm permutation test
Cwarm16S.test=t.test(perm.c16Sslope-c16Sslope,perm.w16Sslope-w16Sslope)
Cwarm16S.t=Cwarm16S.test$statistic
Cwarm16S.p=Cwarm16S.test$p.value







library(ggplot2)
library(ggpmisc)


ggplot(data, aes(x = year, y = Virulence, color = Warm)) +
  geom_point(size = 2.8, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  stat_poly_eq(
    formula = y ~ x,
    aes(label = paste(..rr.label.., ..p.value.label.., sep = "~~~")),
    parse = TRUE,
    size = 4,
    label.x.npc = c(0.05, 0.05),
    label.y.npc = c(0.95, 0.85)) +
  scale_color_manual(values = c("unwarming" = "#1f77b4", "warming" = "#d62728")) +
  theme_bw(base_size = 14) +
  labs(x = "Pathogen diversity", y = "Plant diversity")



fm1=lmer(scale(Virulence) ~ Warm + (1|year)+ (1|block),data=data) 
summary(fm1)
AIC(fm1)
car::Anova(fm1,type=2)
r.squaredGLMM(fm1)





fm1=lmer(scale(plant_richness) ~ scale(richness) + (1|year)+ (1|block),data=data_c) 
summary(fm1)
AIC(fm1)
car::Anova(fm1,type=2)
r.squaredGLMM(fm1)


fm1=lmer(scale(plant_richness) ~ scale(richness) + (1|block),data=data_c) 
summary(fm1)
AIC(fm1)
car::Anova(fm1,type=2)
r.squaredGLMM(fm1)




fm2=lmer(scale(FlC4) ~ scale(RA) + (1|year)+ (1|block),data=data_w) 
summary(fm2)
AIC(fm2)
car::Anova(fm2,type=2)
r.squaredGLMM(fm2)






library(ggplot2)
library(dplyr)
library(tidyr)

cor_df <- data.frame(
  Variable = c("Plant biomass", "C3 biomass", "C4 biomass", "Plant richness"),
  control   = c(-0.253, 0.032, -0.192, 0.132),
  warming   = c( 0.052, 0.061, -0.041, -0.057)
)

p_df <- data.frame(
  Variable = c("Plant biomass", "C3 biomass", "C4 biomass", "Plant richness"),
  control   = c(0.0349, 0.763, 0.125, 0.302),
  warming   = c(0.668, 0.404, 0.753, 0.646)
)

plot_df <- cor_df %>%
  pivot_longer(cols = c(control, warming), names_to = "Treatment", values_to = "r") %>%
  left_join(
    p_df %>% pivot_longer(cols = c(control, warming), names_to = "Treatment", values_to = "p"),
    by = c("Variable", "Treatment")
  ) %>%
  mutate(
    label = ifelse(p < 0.05,
                   paste0(sprintf("%.2f", r), "*"),
                   sprintf("%.2f", r)),
    Variable = factor(Variable, levels = rev(c("Plant biomass", "C3 biomass", "C4 biomass", "Plant richness"))),
    Treatment = factor(Treatment, levels = c("control", "warming"))
  )

p1 <- ggplot(plot_df, aes(Treatment, Variable, fill = r)) +
  geom_tile(color = "white", linewidth = 0.9) +
  geom_text(aes(label = label), size = 4) +
  scale_fill_gradient2(
    low = "#D65F5F",
    mid = "white",
    high = "#4C78A8",
    midpoint = 0,
    limits = c(-0.3, 0.5),
    name = expression(italic(r))
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(color = "black"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

p1




library(ggplot2)
library(dplyr)
library(tidyr)

cor_df <- data.frame(
  Variable = c("Plant biomass", "C3 biomass", "C4 biomass", "Plant richness"),
  control = c(0.428, 0.236, 0.187, 0.062),
  warming = c(0.026, 0.031, 0.018, 0.164)
)

p_df <- data.frame(
  Variable = c("Plant biomass", "C3 biomass", "C4 biomass", "Plant richness"),
  control = c(0.0004, 0.0296, 0.129, 0.663),
  warming = c(0.837, 0.696, 0.893, 0.193))
  

plot_df <- cor_df %>%
  pivot_longer(cols = c(control, warming), names_to = "Treatment", values_to = "r") %>%
  left_join(
    p_df %>% pivot_longer(cols = c(control, warming), names_to = "Treatment", values_to = "p"),
    by = c("Variable", "Treatment")
  ) %>%
  mutate(
    label = ifelse(p < 0.05,
                   paste0(sprintf("%.2f", r), "*"),
                   sprintf("%.2f", r)),
    Variable = factor(Variable, levels = rev(c("Plant biomass", "C3 biomass", "C4 biomass", "Plant richness"))),
    Treatment = factor(Treatment, levels = c("control", "warming"))
  )

p2 <- ggplot(plot_df, aes(Treatment, Variable, fill = r)) +
  geom_tile(color = "white", linewidth = 0.9) +
  geom_text(aes(label = label), size = 4) +
  scale_fill_gradient2(
    low = "#D65F5F",
    mid = "white",
    high = "#4C78A8",
    midpoint = 0,
    limits = c(-0.3, 0.5),
    name = expression(italic(r))
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(color = "black"),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

p2

library(ggpubr)
figures<-ggarrange(p1,p2, ncol=2,nrow=1,align= "v")  #align = "v"??????????
figures
ggsave("Figure pathogen richness and plant traits.pdf",figures,width=16,height=8,units="cm")






# install.packages(c("sf", "dplyr", "ggplot2", "rnaturalearth", "rnaturalearthdata"))

library(sf)
library(dplyr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)

#========================
# 1. 读取 RESOLVE 数据
#========================

shp_path <- "/Users/shun/Downloads/Ecoregions2017/Ecoregions2017.shp"   # 改成你的shp路径
eco <- st_read(shp_path, quiet = TRUE)

#========================
# 2. 提取草地相关区域
#   这里只是数据筛选时用到 BIOME_NAME
#   画图时不会显示 biome
#========================

grass_biomes <- c(
  "Tropical & Subtropical Grasslands, Savannas & Shrublands",
  "Temperate Grasslands, Savannas & Shrublands",
  "Flooded Grasslands & Savannas",
  "Montane Grasslands & Shrublands"
)

# 合并成一个整体图层

grass <- eco %>%
  filter(BIOME_NAME %in% grass_biomes) %>%
  st_make_valid() %>%
  summarise(geometry = st_union(geometry))


#========================
# 3. 世界底图
#========================

world <- ne_countries(scale = "medium", returnclass = "sf")

#========================
# 4. Robinson 投影
#========================

robin_crs <- "+proj=robin +lon_0=0 +datum=WGS84 +units=m +no_defs"

world_r <- st_transform(world, robin_crs)
grass_r <- st_transform(grass, robin_crs)

#========================
# 5. 简单全球草原分布图
#========================

p <- ggplot() +
  geom_sf(data = world_r, fill = "grey92", color = "grey70", linewidth = 0.15) +
  geom_sf(data = grass_r, fill = "#8CCB6E", color = NA) +
  coord_sf(crs = st_crs(robin_crs), expand = FALSE) +
  labs(
    title = "Global grassland distribution"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.background = element_rect(fill = "#EAF4FB", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "#D7E6F5", linewidth = 0.25),
    axis.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18)
  )

p









# =========================================================
# Calculate temporal stability of RA at the block × treatment level
# Data structure:
# - 4 blocks
# - each block has one warming plot and one unwarming plot
# - 12 years of repeated observations
# Temporal stability = mean / sd across years
# =========================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# ---------------------------------
# 0. Read data
# ---------------------------------
# Replace with your real file path if needed
# df <- read.csv("your_data.csv", stringsAsFactors = FALSE)

# If your data are already loaded as df, skip this step.

df=data
# ---------------------------------
# 1. Check whether each block × treatment has 12 years
# ---------------------------------
year_check <- df %>%
  group_by(block, Warm) %>%
  summarise(
    n_year = n_distinct(year),
    years = paste(sort(unique(year)), collapse = ", "),
    .groups = "drop"
  )

print(year_check)


# ---------------------------------
# 2. Calculate temporal stability of RA
# ---------------------------------
# Temporal stability = mean(RA across years) / sd(RA across years)

ra_stability <- df %>%
  group_by(block, Warm) %>%
  summarise(
    n_year = sum(!is.na(richness)),
    mean_RA = mean(richness, na.rm = TRUE),
    sd_RA = sd(richness, na.rm = TRUE),
    cv_RA = sd_RA / mean_RA,
    stability_RA = mean_RA / sd_RA,
    .groups = "drop"
  )

print(ra_stability)


# ---------------------------------
# 3. Reshape to wide format for paired comparison
# ---------------------------------
ra_stability_wide <- ra_stability %>%
  select(block, Warm, mean_RA, sd_RA, stability_RA) %>%
  pivot_wider(
    names_from = Warm,
    values_from = c(mean_RA, sd_RA, stability_RA)
  )

print(ra_stability_wide)


# ---------------------------------
# 4. Paired tests: warming vs unwarming
# ---------------------------------
# Because there are only 4 blocks, Wilcoxon paired test is recommended.
# Paired t-test can also be reported as reference.

# Temporal stability
wilcox_stability <- wilcox.test(
  ra_stability_wide$stability_RA_warming,
  ra_stability_wide$stability_RA_unwarming,
  paired = TRUE
)

ttest_stability <- t.test(
  ra_stability_wide$stability_RA_warming,
  ra_stability_wide$stability_RA_unwarming,
  paired = TRUE
)

print(wilcox_stability)
print(ttest_stability)

# Mean RA
wilcox_mean <- wilcox.test(
  ra_stability_wide$mean_RA_warming,
  ra_stability_wide$mean_RA_unwarming,
  paired = TRUE
)

ttest_mean <- t.test(
  ra_stability_wide$mean_RA_warming,
  ra_stability_wide$mean_RA_unwarming,
  paired = TRUE
)

print(wilcox_mean)
print(ttest_mean)

# SD of RA
wilcox_sd <- wilcox.test(
  ra_stability_wide$sd_RA_warming,
  ra_stability_wide$sd_RA_unwarming,
  paired = TRUE
)

ttest_sd <- t.test(
  ra_stability_wide$sd_RA_warming,
  ra_stability_wide$sd_RA_unwarming,
  paired = TRUE
)

print(wilcox_sd)
print(ttest_sd)


# ---------------------------------
# 5. Plot yearly trajectories of RA within each block
# ---------------------------------
p1 <- ggplot(df, aes(x = year, y = RA, color = Warm, group = Warm)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ block) +
  theme_bw(base_size = 12) +
  labs(
    x = "Year",
    y = "RA",
    color = "Treatment",
    title = "Yearly trajectories of RA in each block"
  )

print(p1)


# ---------------------------------
# 6. Plot temporal stability by block and treatment
# ---------------------------------
p2 <- ggplot(ra_stability, aes(x = block, y = stability_RA, fill = Warm)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  theme_bw(base_size = 12) +
  labs(
    x = "Block",
    y = "Temporal stability of RA",
    fill = "Treatment",
    title = "Temporal stability of RA by block and treatment"
  )

print(p2)


# ---------------------------------
# 7. Optional: plot mean RA and SD separately
# ---------------------------------
p3 <- ggplot(ra_stability, aes(x = block, y = mean_RA, fill = Warm)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  theme_bw(base_size = 12) +
  labs(
    x = "Block",
    y = "Mean RA",
    fill = "Treatment",
    title = "Mean RA across 12 years"
  )

p4 <- ggplot(ra_stability, aes(x = block, y = sd_RA, fill = Warm)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  theme_bw(base_size = 12) +
  labs(
    x = "Block",
    y = "SD of RA",
    fill = "Treatment",
    title = "Interannual variation in RA"
  )

print(p3)
print(p4)


# ---------------------------------
# 8. Save results to files
# ---------------------------------
write.csv(year_check, "RA_year_check_by_block_treatment.csv", row.names = FALSE)
write.csv(ra_stability, "RA_temporal_stability_by_block_treatment.csv", row.names = FALSE)
write.csv(ra_stability_wide, "RA_temporal_stability_wide.csv", row.names = FALSE)

ggsave("RA_yearly_trajectories_by_block.pdf", p1, width = 8, height = 6)
ggsave("RA_temporal_stability_by_block.pdf", p2, width = 7, height = 5)
ggsave("RA_mean_by_block.pdf", p3, width = 7, height = 5)
ggsave("RA_sd_by_block.pdf", p4, width = 7, height = 5)














#计算每个 taxon 的 temporal prevalence 
library(dplyr) 
library(tidyr)
otu = read.delim(file = "even_selected_fungi_plant_pathogen_onlywarm.txt",row.names = 1,header=T,check.names=F)
data=read.delim(file = "data_analysis_pathogen.txt",sep="\t",row.names = 1,header=T,check.names=F)



comm=t(otu)


meta=data[,c(11,13)]
meta <- meta %>%
  as.data.frame() %>%
  tibble::rownames_to_column("SampleID")


library(dplyr)
library(tidyr)
library(tibble)

# 1. OTU表转成长表，并合并metadata
dat_long <- comm %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  left_join(meta[, c("SampleID", "year", "Warm")], by = "SampleID") %>%
  pivot_longer(
    cols = -c(SampleID, year, Warm),
    names_to = "Taxon",
    values_to = "Abundance"
  )

# 2. 在每个 treatment × year 内，判断某 taxon 是否出现
year_presence <- dat_long %>%
  dplyr::group_by(Warm, year, Taxon) %>%
  dplyr::summarise(
    Present = as.integer(any(Abundance > 0))
  ) %>%
  dplyr::ungroup()

head(year_presence)


# 3. 计算按 treatment 分组的 temporal prevalence
temporal_prev_treat <- year_presence %>%
  dplyr::group_by(Warm, Taxon) %>%
  dplyr::summarise(
    n_year_present = sum(Present),
    total_years = dplyr::n_distinct(year),
    temporal_prevalence = n_year_present / total_years
  ) %>%
  dplyr::ungroup()

head(temporal_prev_treat)

table(temporal_prev_treat$Warm)


library(ggplot2)

ggplot(temporal_prev_treat, aes(x = Warm, y = temporal_prevalence, fill = Warm)) +
  geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.12, size = 0.9, alpha = 0.35) +
  theme_bw(base_size = 13) +
  labs(x = NULL, y = "Temporal prevalence")

wilcox.test(temporal_prevalence ~ Warm, data = temporal_prev_treat)


temporal_prev_treat %>%
  group_by(Warm) %>%
  summarise(
    mean_temporal_prevalence = mean(temporal_prevalence, na.rm = TRUE),
    sd_temporal_prevalence = sd(temporal_prevalence, na.rm = TRUE),
    median_temporal_prevalence = median(temporal_prevalence, na.rm = TRUE),
    .groups = "drop"
  )


library(ggplot2)
library(dplyr)

# 先做 Wilcoxon 检验
wt <- wilcox.test(temporal_prevalence ~ Warm, data = temporal_prev_treat)

# 生成标题文本
p_lab <- if (wt$p.value < 0.001) {
  "p < 0.001"
} else {
  paste0("p = ", signif(wt$p.value, 3))
}

title_text <- paste0("W = ", round(wt$statistic), ", ", p_lab)

# 设置分组顺序
temporal_prev_treat$Warm <- factor(
  temporal_prev_treat$Warm,
  levels = c("unwarming", "warming")
)

# 画图
p <- ggplot(temporal_prev_treat, aes(x = Warm, y = temporal_prevalence, fill = Warm)) +
  geom_violin(
    width = 0.75,
    alpha = 0.28,
    color = "black",
    linewidth = 0.7,
    trim = FALSE
  ) +
  geom_boxplot(
    width = 0.16,
    alpha = 0.75,
    outlier.shape = NA,
    color = "black",
    linewidth = 0.7
  ) +
  geom_jitter(
    width = 0.14,
    height = 0,
    size = 1.5,
    alpha = 0.22,
    color = "black"
  ) +
  scale_fill_manual(values = c("unwarming" = "#BFD3E6", "warming" = "#E7C0C0")) +
  coord_cartesian(ylim = c(0, 1.05)) +
  labs(
    x = NULL,
    y = "Temporal prevalence",
    title = title_text
  ) +
  theme_bw(base_size = 15) +
  theme(
    panel.grid.major = element_line(color = "grey85", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    axis.ticks.length = unit(0.2, "cm"),
    legend.position = "none",
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )

p

ggsave("Figure temporal_prevalence_violin_boxplot.pdf",p,width=12,height=10,units="cm")  




ggplot(temporal_prev_treat, aes(x = temporal_prevalence, color = Warm, fill = Warm)) +
  geom_histogram(position = "identity", alpha = 0.25, bins = 13) +
  scale_fill_manual(values = c("unwarming" = "#BFD3E6", "warming" = "#E7C0C0")) +
  scale_color_manual(values = c("unwarming" = "#7EA6C2", "warming" = "#C98E8E")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.title = element_blank()
  ) +
  labs(
    x = "Temporal prevalence",
    y = "Count"
  )




library(vegan)
library(dplyr)
library(tibble)

otu2=t(otu)


# 计算 Bray-Curtis 距离
dist_mat <- vegdist(otu2, method = "bray")
dist_mat <- vegdist(otu2, method = "jaccard", binary = TRUE)
# Sørensen distance
dist_mat <- vegdist(otu2, method = "bray", binary = TRUE)

meta=data
# 计算 dispersion
bd <- betadisper(dist_mat, group = meta$Warm)

# 查看结果
bd
anova(bd)
permutest(bd, permutations = 999)

dispersion_df <- data.frame(
  SampleID = names(bd$distances),
  Warm = meta$Warm,
  dispersion = bd$distances
)

wilcox.test(dispersion ~ Warm, data = dispersion_df)


dat3 <- summarySE(dispersion_df , measurevar = "dispersion", groupvar = c("Warm"))
dat3


head(dispersion_df)



library(vegan)

# 1. 先计算距离矩阵
# 例子：Bray-Curtis
dist_mat <- vegdist(t(otu), method = "bray")

# 2. adonis2
adon_res <- adonis2(
  dist_mat ~ Warm + year * block,
  data = data,
  permutations = 999
)

adon_res


adon_res <- adonis2(
  dist_mat ~ Warm + year * block,
  data = data,
  permutations = 999,
  by = "terms"
)

adon_res






library(dplyr)

path_fit <- dat_long %>%
  group_by(Warm, year, Taxon) %>%
  summarise(
    mean_abundance = mean(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Taxon, Warm, year) %>%
  group_by(Taxon, Warm) %>%
  mutate(
    fitness_rt = log((lead(mean_abundance) + 1e-6) / (mean_abundance + 1e-6))
  ) %>%
  summarise(
    mean_fitness = mean(fitness_rt, na.rm = TRUE),
    .groups = "drop"
  )

head(path_fit)




library(vegan)
library(dplyr)
library(ggplot2)

# 确保样本顺序一致
meta2 <- meta[rownames(otu2), , drop = FALSE]

# Bray-Curtis distance
dist_mat <- vegdist(otu2, method = "bray")

# 转成矩阵，方便提取子矩阵
dist_m <- as.matrix(dist_mat)

# 年内、处理内 pairwise distance 的均值
years <- sort(unique(meta2$year))
treats <- unique(meta2$Warm)

res_list <- list()

for (yy in years) {
  for (tr in treats) {
    
    sam <- rownames(meta2)[meta2$year == yy & meta2$Warm == tr]
    
    # 至少两个样本才能计算组内距离
    if (length(sam) >= 2) {
      subm <- dist_m[sam, sam, drop = FALSE]
      vals <- subm[upper.tri(subm)]
      
      res_list[[paste(yy, tr, sep = "_")]] <- data.frame(
        year = yy,
        Warm = tr,
        mean_bc = mean(vals, na.rm = TRUE),
        median_bc = median(vals, na.rm = TRUE),
        sd_bc = sd(vals, na.rm = TRUE),
        n_pairs = length(vals)
      )
    }
  }
}

within_year_bc <- bind_rows(res_list)

within_year_bc


wilcox.test(mean_bc ~ Warm, data = within_year_bc)
t.test(mean_bc ~ Warm, data = within_year_bc)


bc_summary <- within_year_bc %>%
  group_by(Warm) %>%
  summarise(
    mean_bc = mean(mean_bc, na.rm = TRUE),
    sd_bc = sd(mean_bc, na.rm = TRUE),
    n_year = n(),
    se_bc = sd_bc / sqrt(n_year),
    .groups = "drop"
  )

bc_summary



library(dplyr)
library(ggplot2)
within_year_bc=read.table(pipe("pbpaste"),header=T,row.names=1)  
# 汇总
bc_summary <- within_year_bc %>%
  group_by(Warm) %>%
  summarise(
    mean_bc_over_years = mean(mean_bc, na.rm = TRUE),
    sd_bc = sd(mean_bc, na.rm = TRUE),
    n = n(),
    se_bc = sd_bc / sqrt(n),
    .groups = "drop"
  )


p <- ggplot(bc_summary, aes(x = Warm, y = mean_bc_over_years, fill = Warm)) +
  geom_col(width = 0.65, alpha = 0.95) +
  geom_errorbar(
    aes(ymin = mean_bc_over_years - se_bc,
        ymax = mean_bc_over_years + se_bc),
    width = 0.12,
    linewidth = 0.8
  ) +
  scale_fill_manual(values = c("unwarming" = "#89A5C9", "warming" = "#D9B7B0")) +
  theme_bw(base_size = 17) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 14)
  ) +
  labs(y = "Taxonomic dissimilarity (Bray–Curtis)")

ggsave("/Volumes/SHUN HAN/data from Xue/soil-borne pathogens/iCAMP_pathogen/bin36_new icamp/Taxonomic dissimilarity.pdf", p,width=5,height=10,units="cm") 




















library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

df <- data.frame(
  ID  = c("bin1","bin2","bin3","bin4","bin5","bin6","bin7","bin8"),
  HoS = c(0.07942446, 0.004329,   0.13228558, 0.00679232, 0.13617782, 0.42807263, 0.17502787, 0.03994132),
  HeS = c(0.02033572, 0.02131202, 0.00721248, 0,          0,          0,          0.034699,   0.00906848),
  DL  = c(0.54043165, 0.56843157, 0.17017544, 0.7739854,  0.25660024, 0.15095438, 0.25947603, 0.74548243),
  HD  = c(0.02359712, 0.00333,    0.05151072, 0.01969774, 0.0284449,  0.01716713, 0.0780379,  0.00853504),
  DR  = c(0.33621103, 0.4025974,  0.63881579, 0.19952454, 0.5787704,  0.40380587, 0.52299331, 0.19697273)
)


plot_df <- df %>%
  pivot_longer(
    cols = -ID,
    names_to = "Process",
    values_to = "Value"
  ) %>%
  mutate(
    Process = factor(Process, levels = c("DL", "DR", "HD", "HeS", "HoS"))
  ) %>%
  group_by(ID) %>%
  mutate(
    Percent = Value / sum(Value) * 100
  ) %>%
  ungroup()

plot_df




donut_fun <- function(dat, bin_name){
  subdat <- dat %>% filter(ID == bin_name)
  
  subdat <- subdat %>%
    mutate(
      legend_lab = paste0(Process, " (", round(Percent, 1), "%)")
    )
  
  ggplot(subdat, aes(x = 2, y = Percent, fill = Process)) +
    geom_col(color = "black", linewidth = 1, width = 1) +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +
    theme_void(base_size = 16) +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    ) +
    annotate("text", x = 0, y = 0, label = bin_name, size = 8, fontface = "bold") +
    scale_fill_manual(
      values = c(
        "DL"  = "#D8D06B",  # 黄绿色
        "DR"  = "#A8A9AD",  # 灰色（Drift）
        "HD"  = "#BFD5E8",  # 浅蓝色
        "HeS" = "#889EAF",  # 蓝灰色
        "HoS" = "#BE9A70"   # 棕色
      ),
      breaks = subdat$Process,
      labels = subdat$legend_lab
    )
}

p1 <- donut_fun(plot_df, "bin1")
p1

bins <- unique(plot_df$ID)

plot_list <- lapply(bins, function(x) donut_fun(plot_df, x))
names(plot_list) <- bins

p_all <- wrap_plots(plot_list, ncol = 2)
p_all

ggsave("/Volumes/SHUN HAN/data from Xue/soil-borne pathogens/iCAMP_pathogen/bin36_new icamp/all_bins_donut.pdf", p_all, width = 14, height = 20)








#==================================================
# ieggr time-decay analysis   结果和下面是一致的
#==================================================
setwd("/Users/data from Xue/soil-borne pathogens")
otu = read.delim(file = "even_selected_fungi_plant_pathogen_onlywarm.txt",row.names = 1,header=T,check.names=F)
taxon= read.delim(file = "even_selected_fungi_plant_pathogen_taxon_onlywarm.txt",row.names = 1,header=T,check.names=F)
data=read.delim(file = "data_analysis_pathogen.txt",sep="\t",row.names = 1,header=T,check.names=F)
library(ape)
tree <- read.tree("/Users/data from Xue/soil-borne pathogens/plant_pathogen_tree_onlywarm.nwk.txt")

meta2=data
otu2=t(otu)

common_taxa <- intersect(colnames(otu2), tree$tip.label)
length(common_taxa)
head(common_taxa)

otu2_sub <- otu2[, common_taxa, drop = FALSE]
tree_sub <- keep.tip(tree, common_taxa)
otu2_sub <- otu2_sub[, tree_sub$tip.label, drop = FALSE]
is.rooted(tree_sub)
library(phytools)
tree_sub <- midpoint.root(tree_sub)
is.rooted(tree_sub)



time_mat <- data.frame(year = meta2$year)
rownames(time_mat) <- rownames(meta2)

treat_mat <- data.frame(Warm = meta2$Warm)
rownames(treat_mat) <- rownames(meta2)

ind_mat <- data.frame(block = meta2$plot)
rownames(ind_mat) <- rownames(meta2)

perm_strata <- data.frame(block = meta2$block)
rownames(perm_strata) <- rownames(meta2)

td2 <- tdecay(
  time = time_mat,
  comm = otu2,
  dist.method = "bray",
  abundance.weighted = TRUE,
  treat = treat_mat,
  individuals = ind_mat,
  compare.yn = TRUE,
  perm.test = TRUE,
  boot.strap = TRUE,
  rand = 999,
  Dmax = 1,
  perm.strata.within = perm_strata,
  perm.strata.between = perm_strata,
  boot.strata = perm_strata
)

data_bc=td2$data
#write.csv(data_bc,'/Volumes/SHUN HAN/data from Xue/soil-borne pathogens/data_bc_TDR.csv', quote = FALSE)

TDR= ggplot(data = data_bc, aes(x = log.time.diff, y = log.Sim,color = treatment,group=treatment))+
  geom_point(aes(color=treatment),size=1,alpha=1)+
  xlab("Ln[Time(year)]")+ylab("Ln(Community similarity)")+
  geom_smooth(method = 'lm', formula = y~x, se = TRUE, show.legend = FALSE)+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black',fill = 'transparent'))+theme_bw()+
  theme(legend.key = element_rect(fill = 'transparent'), legend.title = element_blank(),legend.background = element_blank())+
  theme(legend.position="none")



#==================================================
# Manual time-decay analysis   结果和上面是一致的
#==================================================

# 需要的包
library(vegan)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(ggplot2)

#==================================================
# 1. 检查数据
#==================================================

# otu2: 行 = sample, 列 = taxon
# meta2: 行名 = sample, 包含 year / Warm / plot / block

# 确保样本顺序一致
meta2 <- meta2[rownames(otu2), , drop = FALSE]

# 添加 SampleID 列，方便后面 join
meta2$SampleID <- rownames(meta2)

# 简单检查
stopifnot(all(rownames(otu2) == rownames(meta2)))

# 查看基本信息
dim(otu2)
head(rownames(otu2))
colnames(meta2)

#==================================================
# 2. 计算群落 dissimilarity
#==================================================

# 方案 A：Bray-Curtis（丰度加权）
dist_mat <- vegdist(otu2, method = "bray", binary = F)

# 如果你想更接近 presence-absence / Sørensen，可改成：
# dist_mat <- vegdist(otu2, method = "bray", binary = TRUE)

# 转成矩阵，后面方便按样本名提取
dist_m <- as.matrix(dist_mat)

#==================================================
# 3. 构建同一 plot 内、不同年份的样本对
#==================================================

pair_list <- list()

all_plots <- unique(meta2$plot)

for (p in all_plots) {
  
  sub <- meta2 %>%
    filter(plot == p) %>%
    arrange(year)
  
  # 至少两个年份才能做样本对
  if (nrow(sub) < 2) next
  
  # plot内所有两两组合
  cmb <- t(combn(sub$SampleID, 2))
  cmb <- as.data.frame(cmb, stringsAsFactors = FALSE)
  colnames(cmb) <- c("s1", "s2")
  
  pair_list[[as.character(p)]] <- cmb
}

pair_df <- bind_rows(pair_list)

# 检查
head(pair_df)
dim(pair_df)

#==================================================
# 4. 给每个样本对补充元数据
#==================================================

pair_df <- pair_df %>%
  left_join(
    meta2[, c("SampleID", "year", "Warm", "plot", "block")],
    by = c("s1" = "SampleID")
  ) %>%
  rename(
    year1  = year,
    Warm1  = Warm,
    plot1  = plot,
    block1 = block
  ) %>%
  left_join(
    meta2[, c("SampleID", "year", "Warm", "plot", "block")],
    by = c("s2" = "SampleID")
  ) %>%
  rename(
    year2  = year,
    Warm2  = Warm,
    plot2  = plot,
    block2 = block
  )

# 检查
head(pair_df)

#==================================================
# 5. 只保留同一 treatment 内的比较
#==================================================

pair_df <- pair_df %>%
  filter(Warm1 == Warm2) %>%
  mutate(
    Warm = Warm1,
    dT   = abs(year2 - year1)
  ) %>%
  filter(dT > 0)

# 检查
table(pair_df$Warm)
summary(pair_df$dT)

#==================================================
# 6. 提取样本对 dissimilarity，并换算 similarity
#==================================================

pair_df$dissimilarity <- mapply(
  function(a, b) dist_m[a, b],
  pair_df$s1,
  pair_df$s2)

# similarity = 1 - dissimilarity
pair_df$similarity <- 1 - pair_df$dissimilarity

# 避免 log(0)
eps <- 1e-6

pair_df <- pair_df %>%
  mutate(
    similarity = pmax(similarity, eps),
    lnS = log(similarity),
    lnT = log(dT))

# 检查
head(pair_df)
summary(pair_df$similarity)
summary(pair_df$lnS)
summary(pair_df$lnT)

#==================================================
# 7. 分 treatment 拟合 time-decay 模型
# 论文思路：ln(S) ~ ln(T) + (1 + ln(T) | plot)
#==================================================

# warming
td_warm <- pair_df %>%filter(Warm == "warming")
m_warm <- lmer(lnS ~ lnT + (1 + lnT | plot1),data = td_warm)
# unwarming
td_ctrl <- pair_df %>%filter(Warm == "unwarming")
m_ctrl <- lmer(lnS ~ lnT + (1 + lnT | plot1),data = td_ctrl)
# 查看模型
summary(m_warm)
summary(m_ctrl)
R2_warm <- MuMIn::r.squaredGLMM(m_warm)
R2_ctrl <- MuMIn::r.squaredGLMM(m_ctrl)
#==================================================
# 8. 提取 slope、SE、P
# 论文里 v = - slope
#==================================================
coef_warm <- summary(m_warm)$coefficients
coef_ctrl <- summary(m_ctrl)$coefficients

v_warm <- -coef_warm["lnT", "Estimate"]
se_warm <- coef_warm["lnT", "Std. Error"]
p_warm <- coef_warm["lnT", "Pr(>|t|)"]

v_ctrl <- -coef_ctrl["lnT", "Estimate"]
se_ctrl <- coef_ctrl["lnT", "Std. Error"]
p_ctrl <- coef_ctrl["lnT", "Pr(>|t|)"]

cat("warming:   v =", v_warm, "SE =", se_warm, "P =", p_warm, "\n")
cat("unwarming: v =", v_ctrl, "SE =", se_ctrl, "P =", p_ctrl, "\n")
#==================================================
# 9. 检验两处理 slope 是否不同
# lnT:Warmwarming 这一项就是关键
#==================================================
pair_df$Warm <- factor(pair_df$Warm, levels = c("unwarming", "warming"))
m_int <- lmer(lnS ~ lnT * Warm + (1 + lnT | plot1), data = pair_df)
summary(m_int)
anova(m_int)
#==================================================
# 10. 生成拟合线数据
#==================================================
new_warm <- data.frame(
  lnT = seq(min(td_warm$lnT), max(td_warm$lnT), length.out = 100))
new_ctrl <- data.frame(
  lnT = seq(min(td_ctrl$lnT), max(td_ctrl$lnT), length.out = 100))
new_warm$pred <- predict(m_warm, newdata = new_warm, re.form = NA)
new_ctrl$pred <- predict(m_ctrl, newdata = new_ctrl, re.form = NA)

#==================================================
# 11. 生成标注文字
#==================================================

lab_ctrl <- paste0(
  "v = ", round(v_ctrl, 3),
  " \u00B1 ", round(se_ctrl, 3),
  ", P = ", signif(p_ctrl, 2))

lab_warm <- paste0(
  "v = ", round(v_warm, 3),
  " \u00B1 ", round(se_warm, 3),
  ", P = ", signif(p_warm, 2))

#==================================================
# 12. 作图
#==================================================

p <- ggplot() +
  geom_point(data = td_warm, aes(x = lnT, y = lnS),
    color = "#C00000", size = 2.2, alpha = 0.8) +
  geom_point(data = td_ctrl,aes(x = lnT, y = lnS),color = "#1F5AA6",size = 2.2, alpha = 0.8,shape = 17) +
  geom_line(data = new_warm,aes(x = lnT, y = pred),color = "#C00000",linewidth = 1) +
  geom_line(data = new_ctrl,aes(x = lnT, y = pred),color = "#1F5AA6",linewidth = 1) +
  annotate( "text",x = max(td_ctrl$lnT) - 0.1, y = max(td_ctrl$lnS),label = lab_ctrl,color = "#1F5AA6",hjust = 1,size = 5) +
  annotate("text",x = max(td_warm$lnT) - 0.1,y = min(td_warm$lnS) + 0.08,
    label = lab_warm,color = "#C00000",hjust = 1,size = 5) +
  labs(x = "ln[Time (yr)]",y = "ln[Community similarity]")+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black',fill = 'transparent'))+theme_bw()+
  theme(legend.key = element_rect(fill = 'transparent'), legend.title = element_blank(),legend.background = element_blank())+
  theme(legend.position="none")
p
#ggsave("/Users/data from Xue/soil-borne pathogens/TDR_Sørensen.pdf", p,width=12,height=10,units="cm") 












library(phyloseq)

ps <- phyloseq(
  otu_table(as.matrix(otu2_sub), taxa_are_rows = FALSE),
  sample_data(meta2),
  phy_tree(tree_sub)
)

#wu_dist <- phyloseq::distance(ps, method = "unifrac")
wu_dist <- phyloseq::distance(ps, method = "wunifrac")
wu_m <- as.matrix(wu_dist)

library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)

#-----------------------------------
# 1. 确保 meta2 有 SampleID
#-----------------------------------
meta2$SampleID <- rownames(meta2)

# 检查
stopifnot(all(rownames(wu_m) == rownames(meta2)))

#-----------------------------------
# 2. 构建同一 plot 内、不同年份的样本对
#-----------------------------------
pair_list <- list()

all_plots <- unique(meta2$plot)

for (p in all_plots) {
  
  sub <- meta2 %>%
    filter(plot == p) %>%
    arrange(year)
  
  if (nrow(sub) < 2) next
  
  cmb <- t(combn(sub$SampleID, 2))
  cmb <- as.data.frame(cmb, stringsAsFactors = FALSE)
  colnames(cmb) <- c("s1", "s2")
  
  pair_list[[as.character(p)]] <- cmb
}

pair_df <- bind_rows(pair_list)

#-----------------------------------
# 3. 合并样本属性
#-----------------------------------
pair_df <- pair_df %>%
  left_join(
    meta2[, c("SampleID", "year", "Warm", "plot", "block")],
    by = c("s1" = "SampleID")
  ) %>%
  rename(
    year1  = year,
    Warm1  = Warm,
    plot1  = plot,
    block1 = block
  ) %>%
  left_join(
    meta2[, c("SampleID", "year", "Warm", "plot", "block")],
    by = c("s2" = "SampleID")
  ) %>%
  rename(
    year2  = year,
    Warm2  = Warm,
    plot2  = plot,
    block2 = block
  )

#-----------------------------------
# 4. 只保留同一 treatment 内的比较
#-----------------------------------
pair_df <- pair_df %>%
  filter(Warm1 == Warm2) %>%
  mutate(
    Warm = Warm1,
    dT   = abs(year2 - year1)
  ) %>%
  filter(dT > 0)

# 检查
table(pair_df$Warm)
summary(pair_df$dT)

#-----------------------------------
# 5. 提取 weighted UniFrac dissimilarity
#-----------------------------------
pair_df$dissimilarity <- mapply(
  function(a, b) wu_m[a, b],
  pair_df$s1,
  pair_df$s2
)

# weighted UniFrac similarity = 1 - dissimilarity
pair_df$similarity <- 1 - pair_df$dissimilarity

# 避免 log(0)
eps <- 1e-6

pair_df <- pair_df %>%
  mutate(
    similarity = pmax(similarity, eps),
    lnS = log(similarity),
    lnT = log(dT)
  )

head(pair_df)
summary(pair_df$similarity)

#-----------------------------------
# 6. 分别拟合 warming 和 unwarming 的 time-decay
# 论文思路：ln(S) ~ ln(T) + (1 + ln(T) | plot)
#-----------------------------------
td_warm <- pair_df %>%
  filter(Warm == "warming")

td_ctrl <- pair_df %>%
  filter(Warm == "unwarming")

m_warm <- lmer(
  lnS ~ lnT + (1 + lnT | plot1),
  data = td_warm
)

m_ctrl <- lmer(
  lnS ~ lnT + (1 + lnT | plot1),
  data = td_ctrl
)

summary(m_warm)
summary(m_ctrl)
R2_warm <- MuMIn::r.squaredGLMM(m_warm)
R2_ctrl <- MuMIn::r.squaredGLMM(m_ctrl)
#-----------------------------------
# 7. 提取 turnover rate v
# 论文写法：ln(S) = c - v ln(T)
# 所以 v = - slope
#-----------------------------------
coef_warm <- summary(m_warm)$coefficients
coef_ctrl <- summary(m_ctrl)$coefficients

v_warm <- -coef_warm["lnT", "Estimate"]
se_warm <-  coef_warm["lnT", "Std. Error"]
p_warm <-   coef_warm["lnT", "Pr(>|t|)"]

v_ctrl <- -coef_ctrl["lnT", "Estimate"]
se_ctrl <-  coef_ctrl["lnT", "Std. Error"]
p_ctrl <-   coef_ctrl["lnT", "Pr(>|t|)"]

cat("warming:   v =", v_warm, "SE =", se_warm, "P =", p_warm, "\n")
cat("unwarming: v =", v_ctrl, "SE =", se_ctrl, "P =", p_ctrl, "\n")

#-----------------------------------
# 8. 检验两处理 slope 是否不同
# 看交互项 lnT:Warmwarming
#-----------------------------------
pair_df$Warm <- factor(pair_df$Warm, levels = c("unwarming", "warming"))

m_int <- lmer(
  lnS ~ lnT * Warm + (1 + lnT | plot1),
  data = pair_df
)

summary(m_int)
anova(m_int)

#-----------------------------------
# 9. 生成预测线
#-----------------------------------
new_warm <- data.frame(
  lnT = seq(min(td_warm$lnT), max(td_warm$lnT), length.out = 100)
)

new_ctrl <- data.frame(
  lnT = seq(min(td_ctrl$lnT), max(td_ctrl$lnT), length.out = 100)
)

new_warm$pred <- predict(m_warm, newdata = new_warm, re.form = NA)
new_ctrl$pred <- predict(m_ctrl, newdata = new_ctrl, re.form = NA)

#-----------------------------------
# 10. 图中文字
#-----------------------------------
lab_ctrl <- paste0(
  "v = ", round(v_ctrl, 3),
  " \u00B1 ", round(se_ctrl, 3),
  ", P = ", signif(p_ctrl, 2)
)

lab_warm <- paste0(
  "v = ", round(v_warm, 3),
  " \u00B1 ", round(se_warm, 3),
  ", P = ", signif(p_warm, 2)
)

#-----------------------------------
# 11. 作图
#-----------------------------------
p <- ggplot() +
  geom_point(data = td_warm, aes(x = lnT, y = lnS),
             color = "#C00000", size = 2.2, alpha = 0.8) +
  geom_point(data = td_ctrl,aes(x = lnT, y = lnS),color = "#1F5AA6",size = 2.2, alpha = 0.8,shape = 17) +
  geom_line(data = new_warm,aes(x = lnT, y = pred),color = "#C00000",linewidth = 1) +
  geom_line(data = new_ctrl,aes(x = lnT, y = pred),color = "#1F5AA6",linewidth = 1) +
  annotate( "text",x = max(td_ctrl$lnT) - 0.1, y = max(td_ctrl$lnS),label = lab_ctrl,color = "#1F5AA6",hjust = 1,size = 5) +
  annotate("text",x = max(td_warm$lnT) - 0.1,y = min(td_warm$lnS) + 0.08,
           label = lab_warm,color = "#C00000",hjust = 1,size = 5) +
  labs(x = "ln[Time (yr)]",y = "ln[Community similarity]")+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black',fill = 'transparent'))+theme_bw()+
  theme(legend.key = element_rect(fill = 'transparent'), legend.title = element_blank(),legend.background = element_blank())+
  theme(legend.position="none")
p
ggsave("/Users/data from Xue/soil-borne pathogens/TDR_weighted UniFrac.pdf", p,width=12,height=10,units="cm") 





library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)


sample_depth <- rowSums(otu2)
min_depth <- min(sample_depth)
min_depth

# 设置稀释步长
depth_seq <- seq(1, min_depth, by = 20)
if (tail(depth_seq, 1) != min_depth) {
  depth_seq <- c(depth_seq, min_depth)
}

# 计算每个 sample 在不同深度下的 rarefied richness
rare_list <- lapply(rownames(otu2), function(sam) {
  x <- as.numeric(otu2[sam, ])
  data.frame(
    SampleID = sam,
    Depth = depth_seq,
    Richness = sapply(depth_seq, function(d) rarefy(x, sample = d))
  )
})

rare_df <- bind_rows(rare_list)

head(rare_df)


p_rare <- ggplot(rare_df, aes(x = Depth, y = Richness, group = SampleID)) +
  geom_line(color = "grey55", linewidth = 0.6, alpha = 0.8) +
  geom_vline(xintercept = 387, linetype = 2, color = "red", linewidth = 0.8) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank()
  ) +
  labs(
    x = "Sequencing depth (reads)",
    y = "Observed richness"
  )

p_rare










##### GeoChip
library(vegan)
library(dplyr)
dat <- read.table(pipe("pbpaste"), header = TRUE, sep = "\t", check.names = FALSE, row.names=1)
dim(dat)
head(dat[, 1:12])
# 2. 提取 metadata 和 geochip probe matrix
# 前10列是 metadata，第11列以后是 probe 数据
#--------------------------------------------------
meta <- dat[, 1:9]
comm <- dat[, 10:70]

# 用 ID 作为样本名
rownames(meta) <- dat$ID
rownames(comm) <- dat$ID

# 修正 metadata 列名，避免空列名/非法列名
colnames(meta) <- make.names(colnames(meta), unique = TRUE)

# 检查列名
colnames(meta)

#==================================================
# 3. 整理 metadata
#==================================================
meta$Warm  <- as.character(meta$Warm)
meta$block <- as.factor(meta$block)
meta$plot  <- as.character(meta$plot)
meta$year <- as.factor(meta$year)

#==================================================
# 4. community matrix 数值化
#==================================================
comm[] <- lapply(comm, as.numeric)
comm[is.na(comm)] <- 0

# 去掉所有样本都为0的 probe
comm <- comm[, colSums(comm) > 0, drop = FALSE]

#==================================================
# 5. 数据变换
# GeoChip 强度推荐 log1p
#==================================================
comm_use <- log1p(comm)

#==================================================
# 6. 计算群落距离矩阵
#==================================================
dist_mat <- vegdist(comm_use, method = "bray")

#==================================================
# 7. 检查组内 dispersion
#==================================================
bd <- betadisper(dist_mat, group = meta$Warm)
bd_anova <- anova(bd)
bd_perm  <- permutest(bd, permutations = 999)

bd_anova
bd_perm

#==================================================
# 8. Adonis / PERMANOVA
# 你的目标公式：dissimilarity ~ warming + year × block
# 在 R 中写成：dist_mat ~ Warm + year * block
#==================================================

adon_res <- adonis2(
  dist_mat ~ Warm + year * block,
  data = meta,
  permutations = 999,
  by = "terms"
)

adon_res

#==================================================
# 9. 如果你还想限制置换结构
# 注意：如果公式里已经放了 block，再加 strata = block
# 有时会比较保守，但也可以试
#==================================================
adon_res_strata <- adonis2(
  dist_mat ~ Warm + year * block,
  data = meta,
  permutations = 999,
  by = "terms",
  strata = meta$block
)

adon_res_strata

#==================================================
# 10. ANOSIM
# ANOSIM 不能像 adonis2 一样直接放 year * block 公式
# 它只能检验一个 grouping，因此通常只检验 Warm
#==================================================
anosim_res <- anosim(
  x = dist_mat,
  grouping = meta$Warm,
  permutations = 999,
  strata = meta$block)

summary(anosim_res)

#==================================================
# 11. MRPP
# MRPP 也主要是单分组检验
#==================================================
mrpp_res <- mrpp(
  dat = dist_mat,
  grouping = meta$Warm,
  permutations = 999,
  strata = meta$block)

mrpp_res









library(vegan)
library(dplyr)

#==================================================
# 1. 读入数据（Mac：直接从剪贴板）
# Windows 可改成 "clipboard"
#==================================================
dat <- read.table(
  pipe("pbpaste"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

#==================================================
# 2. 提取 metadata 和 probe matrix
# 前10列是 metadata：
# ID, Warm, block, Precipitation, treatment, Clip, plot, warm, precip, year
# 第11列之后是 GeoChip 探针数据
#==================================================
meta <- dat[, 1:9]
comm <- dat[, 10:ncol(dat)]

# 用 ID 作为样本名
rownames(meta) <- dat$ID
rownames(comm) <- dat$ID

# 修正 metadata 列名，避免非法列名
colnames(meta) <- make.names(colnames(meta), unique = TRUE)

# 检查列名
colnames(meta)

#==================================================
# 3. 整理 metadata
#==================================================
meta$Warm  <- as.factor(meta$Warm)
meta$block <- as.factor(meta$block)
meta$plot  <- as.factor(meta$plot)
meta$year  <- as.factor(meta$year)

# 查看分组
table(meta$Warm)
table(meta$block)
table(meta$year)

#==================================================
# 4. community matrix 数值化
#==================================================
comm[] <- lapply(comm, as.numeric)

# 缺失值设为0
comm[is.na(comm)] <- 0

# 去掉所有样本都为0的探针
comm <- comm[, colSums(comm) > 0, drop = FALSE]

#==================================================
# 5. 数据变换
# GeoChip 强度建议用 log1p，保留0值
#==================================================
comm_use <- log1p(comm)

#==================================================
# 6. Bray-Curtis 距离矩阵
#==================================================
dist_mat <- vegdist(comm_use, method = "bray")

#==================================================
# 7. 组内 dispersion 检验
#==================================================
bd <- betadisper(dist_mat, group = meta$Warm)

bd_anova <- anova(bd)
bd_perm  <- permutest(bd, permutations = 9999)

bd_anova
bd_perm

#==================================================
# 8. Adonis / PERMANOVA
# 公式：dissimilarity ~ Warm + year * block
#==================================================
adon_res <- adonis2(
  dist_mat ~ Warm + year * block,
  data = meta,
  permutations = 999,
  by = "terms"
)

adon_res

#==================================================
# 9. 可选：更保守版本，限制置换在 block 内
#==================================================
adon_res_block <- adonis2(
  dist_mat ~ Warm + year * block,
  data = meta,
  permutations = 999,
  by = "terms",
  strata = meta$block
)

adon_res_block

#==================================================
# 10. ANOSIM
# ANOSIM 只能对一个 grouping 做检验
# 这里检验 Warm，并把置换限制在 year 内
#==================================================
anosim_res <- anosim(
  x = dist_mat,
  grouping = meta$Warm,
  permutations = 999,
  strata = meta$year
)

summary(anosim_res)

# 如果需要图：
# plot(anosim_res)

#==================================================
# 11. MRPP
# 同样检验 Warm，并把置换限制在 year 内
#==================================================
mrpp_res <- mrpp(
  dat = dist_mat,
  grouping = meta$Warm,
  permutations = 999,
  strata = meta$year
)

mrpp_res

#==================================================
# 12. 提取简洁结果表
#==================================================
result_summary <- data.frame(
  Method = c("Adonis", "Adonis_block", "ANOSIM", "MRPP"),
  Statistic = c(
    adon_res$F[1],
    adon_res_block$F[1],
    anosim_res$statistic,
    mrpp_res$delta
  ),
  P = c(
    adon_res$`Pr(>F)`[1],
    adon_res_block$`Pr(>F)`[1],
    anosim_res$signif,
    mrpp_res$Pvalue
  )
)

result_summary
write.csv(result_summary, "community_difference_summary.csv", row.names = FALSE)

#==================================================
# 13. 保存详细结果
#==================================================
capture.output(adon_res, file = "adonis_Warm_yearXblock.txt")
capture.output(adon_res_block, file = "adonis_Warm_yearXblock_blockstrata.txt")
capture.output(summary(anosim_res), file = "anosim_Warm_strataYear.txt")
capture.output(mrpp_res, file = "mrpp_Warm_strataYear.txt")
capture.output(bd_anova, file = "betadisper_anova.txt")
capture.output(bd_perm, file = "betadisper_permutest.txt")




library(dplyr)
library(readr)
library(stringr)
library(tidyr)
##############################################################
###############################   NST results
#==================================================
setwd("/Users/data from Xue/soil-borne pathogens")
data=read.delim(file = "data_analysis_pathogen.txt",sep="\t",header=T,check.names=F)
#

NST=read.csv(file = "/Users/data from Xue/soil-borne pathogens/iCAMP_pathogen/NWS.borne.tNST.pairwise_sorensen.Warm.csv",row.names = 1,header=T,check.names=F)
#
NST=read.csv(file = "/Users/data from Xue/soil-borne pathogens/iCAMP_pathogen/NWS.borne.pNST.pairwise.Warm.csv",row.names = 1,header=T,check.names=F)

normalize_id <- function(x) {
  x <- as.character(x)
  # X17plot10N -> X2017a10N
  x <- gsub("^X17plot([0-9]+)([NS])$", "X2017a\\1\\2", x)
  x <- gsub("^X18plot([0-9]+)([NS])$", "X2018a\\1\\2", x)
  x <- gsub("^X19plot([0-9]+)([NS])$", "X2019a\\1\\2", x)
  x <- gsub("^X20plot([0-9]+)([NS])$", "X2020a\\1\\2", x)
  return(x)
}

meta <- data %>% mutate(id_norm = normalize_id(id))
nst <- NST %>%mutate(name1_norm = normalize_id(name1),name2_norm = normalize_id(name2))
meta_small = meta %>% select(id,id_norm,year,Warm,warm,block,plot,Precipitation,precip,clip,Clip)
table(meta_small$Warm, useNA = "ifany")
table(meta_small$year, meta_small$Warm)

# 4. 给 name1 匹配分组信息
nst_meta <- nst %>%
  left_join(meta_small %>%rename_with(~ paste0(.x, "1"), everything()), by = c("name1_norm" = "id_norm1")) %>%
  left_join(meta_small %>%rename_with(~ paste0(.x, "2"), everything()), by = c("name2_norm" = "id_norm2"))

nst_same_year_warm <- nst_meta %>% mutate( same_year = year1 == year2, same_Warm = Warm1 == Warm2,Warm_group = Warm1) %>%
                      filter(same_year,same_Warm,Warm_group %in% c("warming", "unwarming"))
# 查看每年每组 pair 数量
table(nst_same_year_warm$year1, nst_same_year_warm$Warm_group)

# 5. 按 year 和 Warm 分组筛选 pairwise comparisons
#tNST base on sorensen distance
nst_same_year_warm <- nst_meta %>%
  mutate( same_year = year1 == year2, same_Warm = Warm1 == Warm2, Warm_group = Warm1) %>%
  filter(same_year,same_Warm, Warm_group %in% c("warming", "unwarming")) %>%
  select(year = year1,Warm_group,group,name1,name2,block1,block2,plot1,plot2,Warm1,Warm2,warm1,warm2,
    `C.ij.sorensen`,`E.ij.sorensen`,`ST.ij.sorensen`,`NST.ij.sorensen`,`MST.ij.sorensen`,`SES.ij.sorensen`,`RC.ij.sorensen`)

#pNST
nst_same_year_warm <- nst_meta %>%
  mutate( same_year = year1 == year2, same_Warm = Warm1 == Warm2, Warm_group = Warm1) %>%
  filter(same_year,same_Warm, Warm_group %in% c("warming", "unwarming")) %>%
  select(year = year1,Warm_group,group,name1,name2,block1,block2,plot1,plot2,Warm1,Warm2,warm1,warm2,
         `C.ij.bMNTD`,`E.ij.bMNTD`,`ST.ij.bMNTD`,`NST.ij.bMNTD`,`MST.ij.bMNTD`,`SES.ij.bMNTD`,`RC.ij.bMNTD`)


# 查看结果
head(nst_same_year_warm)
table(nst_same_year_warm$year, nst_same_year_warm$Warm_group)



nst_year_summary <- nst_same_year_warm %>%
  dplyr::group_by(year, Warm_group) %>%
  dplyr::summarise(
    n_pairs = dplyr::n(),
    
    mean_NST = mean(`NST.ij.bMNTD`, na.rm = TRUE),
    sd_NST   = sd(`NST.ij.bMNTD`, na.rm = TRUE),
    se_NST   = sd_NST / sqrt(n_pairs),
    
    mean_ST  = mean(`ST.ij.bMNTD`, na.rm = TRUE),
    sd_ST    = sd(`ST.ij.bMNTD`, na.rm = TRUE),
    se_ST    = sd_ST / sqrt(n_pairs),
    
    mean_MST = mean(`MST.ij.bMNTD`, na.rm = TRUE),
    sd_MST    = sd(`MST.ij.bMNTD`, na.rm = TRUE),
    se_MST    = sd_MST / sqrt(n_pairs),
    
    .groups = "drop"
  )

nst_year_summary





nst_year_test <- nst_same_year_warm %>%
  group_by(year) %>%
  group_modify(~ {
    
    dat_i <- .x
    
    if (length(unique(dat_i$Warm_group)) < 2) {
      return(data.frame(
        n_warming = sum(dat_i$Warm_group == "warming"),
        n_unwarming = sum(dat_i$Warm_group == "unwarming"),
        mean_warming = NA,
        mean_unwarming = NA,
        diff_warming_minus_unwarming = NA,
        t_value = NA,
        p_value = NA
      ))
    }
    
    tt <- t.test(`MST.ij.bMNTD` ~ Warm_group, data = dat_i)
    
    mean_warming <- mean(dat_i$`MST.ij.bMNTD`[dat_i$Warm_group == "warming"], na.rm = TRUE)
    mean_unwarming <- mean(dat_i$`MST.ij.bMNTD`[dat_i$Warm_group == "unwarming"], na.rm = TRUE)
    
    data.frame(
      n_warming = sum(dat_i$Warm_group == "warming"),
      n_unwarming = sum(dat_i$Warm_group == "unwarming"),
      mean_warming = mean_warming,
      mean_unwarming = mean_unwarming,
      diff_warming_minus_unwarming = mean_warming - mean_unwarming,
      t_value = unname(tt$statistic),
      p_value = tt$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH"),
    Significant = case_when(
      p_adj_BH < 0.001 ~ "***",
      p_adj_BH < 0.01  ~ "**",
      p_adj_BH < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

nst_year_test



library(lme4)
library(lmerTest)

nst_same_year_warm <- nst_same_year_warm %>%
  mutate(
    year = factor(year),
    block_pair = paste(block1, block2, sep = "_"),
    Warm_group = factor(Warm_group, levels = c("unwarming", "warming")))

m_nst2 <- lmer(
  `MST.ij.bMNTD` ~ Warm_group + (1 | year) + (1 | block_pair),
  data = nst_same_year_warm, REML = TRUE)

summary(m_nst2)
anova(m_nst2)



library(ggplot2)

p_nst_year <- ggplot(
  nst_year_summary,
  aes(x = year, y = mean_MST, color = Warm_group, group = Warm_group)
) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.6) +
  geom_errorbar(
    aes(
      ymin = mean_MST - se_MST,
      ymax = mean_MST + se_MST
    ),
    width = 0.15,
    linewidth = 0.5
  ) +
  scale_color_manual(
    values = c(
      "unwarming" = "#4C78A8",
      "warming" = "#D95F02"
    )
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  labs(
    x = NULL,
    y = "MST"
  )

p_nst_year


############### NST plot Overall community stochasticity under warming and control conditions.
dataplot=read.table(pipe("pbpaste"),header=T,row.names=1)

library(ggplot2)

df <- data.frame(
  Distance = c("Sorensen", "Sorensen", "bMNTD", "bMNTD"),
  Warm = c("control", "warming", "control", "warming"),
  mean = c(66.20, 57.21, 68.36, 65.29),
  se = c(2.4, 2.8, 2.1, 1.9)
)

df$Distance <- factor(df$Distance, levels = c("Sorensen", "bMNTD"))
df$Warm <- factor(df$Warm, levels = c("control", "warming"))

p= ggplot(df, aes(x = Warm, y = mean, fill = Warm)) + ylim(0,80)+
  geom_col(width = 0.65, color = "black") +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0.2,
    size = 0.6
  ) +
  facet_wrap(~Distance, scales = "free_x") +
  scale_fill_manual(values = c("control" = "#4DBBD5", "warming" = "#E64B35")) +
  labs(
    x = NULL,
    y = "Mean ± SE",
    fill = NULL
  ) + ylab("Stochasticity (%)")+
  theme_bw() + theme(legend.position="none")
p

# ggsave("/Users/data from Xue/soil-borne pathogens/Figure NST.pdf", p,width=10,height=8,units="cm") 










otu_table=otu


otu_table[] <- lapply(otu_table, as.numeric)
otu_table[is.na(otu_table)] <- 0

samples <- intersect(colnames(otu_table), rownames(meta))

otu <- otu_table[, samples, drop = FALSE]
meta <- meta[samples, , drop = FALSE]

# 删除全零 zOTU 和全零样品
otu <- otu[rowSums(otu) > 0, , drop = FALSE]

keep <- colSums(otu) > 0
otu <- otu[, keep, drop = FALSE]
meta <- meta[colnames(otu), , drop = FALSE]

stopifnot(all(colnames(otu) == rownames(meta)))

# 整理处理名称
meta$Warm <- tolower(trimws(meta$Warm))

meta$Warm <- case_when(
  meta$Warm %in% c("control", "unwarming", "unwarmed", "ambient") ~ "Control",
  meta$Warm %in% c("warming", "warmed", "warm") ~ "Warming",
  TRUE ~ NA_character_
)

meta$Warm <- factor(meta$Warm, levels = c("Control", "Warming"))
meta$year <- factor(meta$year)
meta$block <- factor(meta$block)
meta$plot <- factor(meta$plot)

# 删除处理信息缺失样品
meta <- meta[!is.na(meta$Warm), , drop = FALSE]
otu <- otu[, rownames(meta), drop = FALSE]

# community matrix：行是样品，列是 zOTU
comm <- t(otu)

#====================================================
# 3. 每年分别计算 community dispersion
#====================================================

calc_dispersion <- function(comm, meta, method = "bray") {
  
  result <- lapply(levels(meta$year), function(yr) {
    
    ids <- rownames(meta)[meta$year == yr]
    
    comm_y <- comm[ids, , drop = FALSE]
    meta_y <- meta[ids, , drop = FALSE]
    
    # 删除该年份全零 zOTU
    comm_y <- comm_y[, colSums(comm_y) > 0, drop = FALSE]
    
    # 每个处理至少需要两个样品
    if (length(table(meta_y$Warm)) < 2 || any(table(meta_y$Warm) < 2)) {
      return(NULL)
    }
    
    # Bray-Curtis 或 Sørensen
    if (method == "bray") {
      dis <- vegdist(comm_y, method = "bray")
      distance_name <- "Bray-Curtis"
    }
    
    if (method == "sorensen") {
      dis <- vegdist(comm_y, method = "bray", binary = TRUE)
      distance_name <- "Sorensen"
    }
    
    bd <- betadisper(
      dis,
      group = meta_y$Warm,
      type = "median",
      bias.adjust = TRUE
    )

    data.frame(
      SampleID = names(bd$distances),
      Dispersion = as.numeric(bd$distances),
      Distance = distance_name,
      stringsAsFactors = FALSE
    ) %>%
      left_join(
        meta_y %>% rownames_to_column("SampleID"),
        by = "SampleID"
      )
  })
  
  bind_rows(result)
}


# 分别计算两种距离
disp_bray <- calc_dispersion(comm, meta, method = "bray")
disp_sorensen <- calc_dispersion(comm, meta, method = "sorensen")

dispersion <- bind_rows(disp_bray, disp_sorensen)



dat2 <- summarySE(disp_sorensen, measurevar = "Dispersion", groupvar = c("Warm"))


fm1=lmer(scale(Sum) ~ warm + (1|year)+ (1|block),data=data) 
summary(fm1)
AIC(fm1)













############### based on metageonomics

setwd("/Users/data from Xue/soil-borne pathogens/metageonomices")
library(FUNGuildR)
fung <- readRDS('/Users/data from Xue/soil-borne pathogens/funguild.rds')  #加载数据库


Fungi_otu = read.delim(file = "all_samples_fungal_species_new_est_reads_96.txt",sep="\t",row.names = 1,header=T,check.names=F,na.strings = "")
Fungi_otu=Fungi_otu[,c(4:99)]
Fungi_taxon=read.delim(file = "dataset2_taxonomy_format.txt",sep="\t",row.names = 1,header=T,check.names=F,na.strings = "")

Fungi_taxon1=Fungi_taxon[match(row.names(Fungi_otu),row.names(Fungi_taxon)),]
otu=cbind(Fungi_otu,Fungi_taxon1)
names(otu)[97] <- c('taxonomy')


fung_guilds <- funguild_assign(otu, db = fung, tax_col = 'taxonomy')

write.table(fung_guilds, 'all_samples_fungal_species_new_est_reads_96_funguilds.txt', sep = '\t', row.names = F, quote = FALSE)



###### validation
richness <- ggplot(data,aes(x=richness,y=richness_meta_no_Bracken))+
  geom_point(size=2,alpha=0.8)+ scale_color_viridis_d()+
  #scale_shape_manual(values = c(16,1))+   
  theme_bw()+   
  theme(legend.position="none")+
  geom_smooth(method = "lm",formula =y ~ poly(x, 1))+
  stat_poly_eq(aes(label = paste( ..rr.label.., stat(p.value.label), sep = '~`,`~')),
               formula = y ~ poly(x, 1), na.rm=T,parse = TRUE, label.x.npc = 0,  label.y.npc = 800, size = 5)

RA <- ggplot(data,aes(x=RA,y=relativeabundance_Bracken))+
  geom_point(size=2,alpha=0.8)+ scale_color_viridis_d()+
  #scale_shape_manual(values = c(16,1))+   
  theme_bw()+   
  theme(legend.position="none")+
  geom_smooth(method = "lm",formula =y ~ poly(x, 1))+
  stat_poly_eq(aes(label = paste( ..rr.label.., stat(p.value.label), sep = '~`,`~')),
               formula = y ~ poly(x, 1), na.rm=T,parse = TRUE, label.x.npc = 0,  label.y.npc = 800, size = 5)


NMDS <- ggplot(data,aes(x=bray.wt_2D_NMDS1,y=NMDS1_Bracken))+
  geom_point(size=2,alpha=0.8)+ scale_color_viridis_d()+
  #scale_shape_manual(values = c(16,1))+   
  theme_bw()+   
  theme(legend.position="none")+
  geom_smooth(method = "lm",formula =y ~ poly(x, 1))+
  stat_poly_eq(aes(label = paste( ..rr.label.., stat(p.value.label), sep = '~`,`~')),
               formula = y ~ poly(x, 1), na.rm=T,parse = TRUE, label.x.npc = 0,  label.y.npc = 800, size = 5)


library(ggpubr)
figures<-ggarrange(RA,richness,NMDS, ncol=3,nrow=1,align= "v")  #align = "v"??????????
figures
ggsave("Figure ITS and metagenomic_validation.pdf",figures,width=18,height=6,units="cm")






library(vegan)
library(dplyr)

#==================================================
# 1. 读入数据（Mac：直接从剪贴板）
# Windows 可改成 "clipboard"
#==================================================
dat <- read.table(
  pipe("pbpaste"),
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

#==================================================
# 2. 提取 metadata 和 probe matrix
# 前10列是 metadata：
# ID, Warm, block, Precipitation, treatment, Clip, plot, warm, precip, year
# 第11列之后是 GeoChip 探针数据
#==================================================
meta <- dat[, 1:12]
comm <- dat[, 13:ncol(dat)]

# 用 ID 作为样本名
rownames(meta) <- dat$ID
rownames(comm) <- dat$ID

# 修正 metadata 列名，避免非法列名
colnames(meta) <- make.names(colnames(meta), unique = TRUE)

# 检查列名
colnames(meta)

#==================================================
# 3. 整理 metadata
#==================================================
meta$Warm  <- as.factor(meta$Warm)
meta$block <- as.factor(meta$block)
meta$plot  <- as.factor(meta$plot)
meta$year  <- as.factor(meta$year)

# 查看分组
table(meta$Warm)
table(meta$block)
table(meta$year)

#==================================================
# 4. community matrix 数值化
#==================================================
comm[] <- lapply(comm, as.numeric)

# 缺失值设为0
comm[is.na(comm)] <- 0

# 去掉所有样本都为0的探针
comm <- comm[, colSums(comm) > 0, drop = FALSE]

#==================================================
# 6. Bray-Curtis 距离矩阵
#==================================================
dist_mat <- vegdist(comm, method = "bray")

#==================================================
# 7. 组内 dispersion 检验
#==================================================
bd <- betadisper(dist_mat, group = meta$Warm)

bd_anova <- anova(bd)
bd_perm  <- permutest(bd, permutations = 9999)

bd_anova
bd_perm

#==================================================
# 8. Adonis / PERMANOVA
# 公式：dissimilarity ~ Warm + year * block
#==================================================
adon_res <- adonis2(
  dist_mat ~ Warm + year * block,
  data = meta,
  permutations = 999,
  by = "terms"
)

adon_res

#==================================================
# 9. 可选：更保守版本，限制置换在 block 内
#==================================================
adon_res_block <- adonis2(
  dist_mat ~ Warm + year * block,
  data = meta,
  permutations = 999,
  by = "terms",
  strata = meta$block
)

adon_res_block

#==================================================
# 10. ANOSIM
# ANOSIM 只能对一个 grouping 做检验
# 这里检验 Warm，并把置换限制在 year 内
#==================================================
anosim_res <- anosim(
  x = dist_mat,
  grouping = meta$Warm,
  permutations = 999,
  strata = meta$year
)

summary(anosim_res)

# 如果需要图：
# plot(anosim_res)

#==================================================
# 11. MRPP
# 同样检验 Warm，并把置换限制在 year 内
#==================================================
mrpp_res <- mrpp(
  dat = dist_mat,
  grouping = meta$Warm,
  permutations = 999,
  strata = meta$year
)

mrpp_res

#==================================================
# 12. 提取简洁结果表
#==================================================
result_summary <- data.frame(
  Method = c("Adonis", "Adonis_block", "ANOSIM", "MRPP"),
  Statistic = c(
    adon_res$F[1],
    adon_res_block$F[1],
    anosim_res$statistic,
    mrpp_res$delta
  ),
  P = c(
    adon_res$`Pr(>F)`[1],
    adon_res_block$`Pr(>F)`[1],
    anosim_res$signif,
    mrpp_res$Pvalue
  )
)

result_summary
write.csv(result_summary, "community_difference_summary.csv", row.names = FALSE)

#==================================================
# 13. 保存详细结果
#==================================================
capture.output(adon_res, file = "adonis_Warm_yearXblock.txt")
capture.output(adon_res_block, file = "adonis_Warm_yearXblock_blockstrata.txt")
capture.output(summary(anosim_res), file = "anosim_Warm_strataYear.txt")
capture.output(mrpp_res, file = "mrpp_Warm_strataYear.txt")
capture.output(bd_anova, file = "betadisper_anova.txt")
capture.output(bd_perm, file = "betadisper_permutest.txt")











###############
##### CCA replot
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(scales)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(scico)
library(grid)
#====================================================
# 1. 重新读取原始数据
#====================================================

sample_df <- read.csv(
  "cca.sample.csv",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

env_df <- read.csv(
  "cca.env.csv",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

colnames(sample_df)[1] <- "SampleID"
colnames(env_df)[1] <- "Variable"

#====================================================
# 2. 整理样品信息
#====================================================

sample_df <- sample_df %>%
  mutate(
    year = factor(
      as.character(year),
      levels = as.character(2009:2020)
    ),
    
    Warm = case_when(
      Warm == "unwarming" ~ "Control",
      Warm == "warming"   ~ "Warming",
      TRUE                ~ Warm
    ),
    
    Warm = factor(
      Warm,
      levels = c("Control", "Warming")
    )
  )

# 检查
table(sample_df$Warm, useNA = "always")
table(sample_df$year, useNA = "always")

year_cols <- c(
  "2009" = "#3F4A6D",
  "2010" = "#50678A",
  "2011" = "#6083A1",
  "2012" = "#6A9EAD",
  "2013" = "#5AA5A4",
  "2014" = "#4E9D8D",
  "2015" = "#6FAE91",
  "2016" = "#92B486",
  "2017" = "#B7AD73",
  "2018" = "#D09A67",
  "2019" = "#D87A62",
  "2020" = "#BC5C5B"
)

env_df <- env_df %>%
  mutate(
    label = case_when(
      Variable == "year1"              ~ "Year",
      Variable == "temperature_annual" ~ "Soil temperature",
      Variable == "annual_moisture"    ~ "Soil moisture",
      Variable == "NO3"                ~ "NO3-N",
      Variable == "NH4"                ~ "NH4-N",
      Variable == "TN"                 ~ "TN",
      Variable == "TC"                 ~ "TC",
      Variable == "pH"                 ~ "pH",
      Variable == "FlTotl"             ~ "Plant biomass",
      Variable == "FlC4"               ~ "C4 plant biomass",
      Variable == "FlC3"               ~ "C3 plant biomass",
      Variable == "plant_richness"     ~ "Plant richness",
      TRUE                             ~ Variable
    )
  )

arrow_scale <- 2.5
env_df_plot <- env_df %>%mutate(
    CCA1_plot = CCA1 * arrow_scale,
    CCA2_plot = CCA2 * arrow_scale)





p_cca <- ggplot() +
  geom_hline(yintercept = 0,linetype = "dashed",linewidth = 0.3,colour = "#BDBDBD") +
  geom_vline(xintercept = 0,linetype = "dashed",linewidth = 0.3,colour = "#BDBDBD") +
  geom_point(data = sample_df,aes(x = CCA1,y = CCA2, colour = year,shape = Warm),size = 2.7, alpha = 0.9) +
  geom_segment(data = env_df_plot,aes(x = 0, y = 0, xend = CCA1_plot, yend = CCA2_plot),
               arrow = arrow(length = unit(0.14, "cm"),type = "closed"),linewidth = 0.45,colour = "#3B3B3B") +
  geom_text_repel(data = env_df_plot,aes(x = CCA1_plot,y = CCA2_plot,label = label),size = 3.0,
    colour = "#262626",box.padding = 0.25,
    point.padding = 0.12,segment.colour = "#AAAAAA",
    min.segment.length = 0,max.overlaps = Inf,seed = 123) +
  scale_colour_manual(values = year_cols,drop = FALSE, name = "Year") +
  scale_shape_manual(values = c("Control" = 16,"Warming" = 17),name = "Treatment") +
  labs(x = "CCA1 (17.84%)",y = "CCA2 (12.52%)") +
  annotate("text",x = -Inf,y = Inf,label = "F = 1.16\np = 0.008",hjust = -0.08,vjust = 1.15,size = 3.5) +
  guides(colour = guide_legend(ncol = 1,byrow = TRUE,order = 2,override.aes = list(size = 2.5)),
    shape = guide_legend(order = 1,override.aes = list(size = 2.8))) +
  theme_classic(base_size = 11) +
  theme(axis.title = element_text(size = 11, colour = "#222222"),
    axis.text = element_text(size = 9,colour = "#333333"),
    axis.line = element_line(linewidth = 0.45,colour = "#333333"),
    axis.ticks = element_line(linewidth = 0.35,colour = "#333333"),
    legend.position = "right",
    legend.title = element_text(size = 9.5),
    legend.text = element_text(size = 8.5),
    legend.key.height = unit(0.33, "cm"),
    plot.margin = margin(7, 12, 7, 7))+ theme_bw()

p_cca

ggsave("Figure CCA.pdf", p_cca, width =16, height = 12, units="cm")












library(ggplot2)
library(dplyr)
library(tidyr)

#====================================================
# 1. 输入数据
#====================================================

cor_df <- data.frame(
  Variable = c(
    "FlTotl",
    "FlC4",
    "FlC3",
    "Plant richness",
    "ER",
    "GPP",
    "NEE",
    "Autotrophic",
    "Heterotrophic",
    "Respiration"
  ),
  
  richness = c(
    0.2866222,
    0.08930915,
    0.21537932,
    0.234377,
    0.18979064,
    0.04707784,
    -0.0834843,
    0.192762,
    0.2552552,
    0.3169059
  ),
  
  pd = c(
    0.2866222,
    0.08930915,
    0.21537932,
    0.234377,
    0.18979064,
    0.04707784,
    -0.0834843,
    0.192762,
    0.2552552,
    0.3169059
  )
)

#====================================================
# 2. 转成长表
#====================================================

plot_df <- cor_df %>%
  pivot_longer(
    cols = c(richness, pd),
    names_to = "Metric",
    values_to = "r"
  ) %>%
  mutate(
    
    label = sprintf("%.2f", r),
    
    Variable = case_when(
      Variable == "FlTotl" ~ "Plant biomass",
      Variable == "FlC4" ~ "C4 biomass",
      Variable == "FlC3" ~ "C3 biomass",
      Variable == "Plant richness" ~ "Plant richness",
      Variable == "ER" ~ "ER",
      Variable == "GPP" ~ "GPP",
      Variable == "NEE" ~ "NEE",
      Variable == "Autotrophic" ~ "Autotrophic respiration",
      Variable == "Heterotrophic" ~ "Heterotrophic respiration",
      Variable == "Respiration" ~ "Soil respiration",
      TRUE ~ Variable
    ),
    
    # 控制上下顺序
    Variable = factor(
      Variable,
      levels = rev(c(
        "Plant biomass",
        "C4 biomass",
        "C3 biomass",
        "Plant richness",
        "ER",
        "GPP",
        "NEE",
        "Autotrophic respiration",
        "Heterotrophic respiration",
        "Soil respiration"
      ))
    ),
    
    # 列顺序
    Metric = factor(
      Metric,
      levels = c("richness", "pd"),
      labels = c("Richness", "PD")
    )
  )


#====================================================
# 3. Heatmap
#====================================================

p2 <- ggplot(
  plot_df,
  aes(
    x = Metric,
    y = Variable,
    fill = r
  )
) +
  
  geom_tile(
    color = "white",
    linewidth = 1.2
  ) +
  
  geom_text(
    aes(label = label),
    size = 4.2,
    color = "black"
  ) +
  
  scale_fill_gradient2(
    low = "#C95D4A",      # negative = red
    mid = "#F7F7F7",      # zero = light grey
    high = "#3F78B5",     # positive = blue
    midpoint = 0,
    limits = c(-0.4, 0.4),
    name = expression(italic(r))
  ) +
  
  coord_fixed() +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    axis.title = element_blank(),
    
    axis.text.x = element_text(
      face = "bold",
      size = 12,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 11,
      color = "black"
    ),
    
    axis.ticks = element_blank(),
    
    axis.line = element_blank(),
    
    legend.title = element_text(
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    legend.key.height = unit(
      1.6,
      "cm"
    ),
    
    panel.background = element_blank()
  )

p2









library(vegan)
library(ggplot2)
library(dplyr)

#====================================================
# 1. Read OTU table
#====================================================

otu <- read.delim(
  "even_selected_fungi_plant_pathogen_onlywarm.txt",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

# zOTU × sample
dim(otu)

#====================================================
# 2. Read metadata
#====================================================

meta <- read.delim(
  "data_analysis_pathogen.txt",
  header = TRUE,
  check.names = FALSE
)

dim(meta)

head(meta[, c("id", "year", "Warm", "block")])


#====================================================
# 3. Check sample names
#====================================================

all(colnames(otu) %in% meta$id)

# 应该是 TRUE


#====================================================
# 4. Arrange metadata in exactly the same order
#====================================================

meta2 <- meta[
  match(colnames(otu), meta$id),
]

# 再检查
identical(
  colnames(otu),
  meta2$id
)

# 应该 TRUE


#====================================================
# 5. Transpose OTU table
# sample × zOTU
#====================================================

comm <- t(otu)

dim(comm)

# 应该：
# 96 samples × 556 zOTUs


#====================================================
# 6. Remove zero-sum taxa
#====================================================

comm <- comm[
  ,
  colSums(comm) > 0
]

dim(comm)



#====================================================
# 7. DCA
#====================================================

dca <- decorana(comm)

dca
summary(dca)
#====================================================
# 8. Site scores
#====================================================

site_scores <- scores(
  dca,
  display = "sites",
  choices = c(1, 2)
)

dca_df <- as.data.frame(site_scores)

dca_df$SampleID <- rownames(dca_df)
#====================================================
# 9. Add metadata
#====================================================

dca_df <- dca_df %>%
  left_join(
    meta %>%
      select(
        id,
        year,
        Warm,
        block
      ),
    by = c("SampleID" = "id")
  )

head(dca_df)
dca_df <- dca_df %>%
  mutate(
    Treatment = case_when(
      Warm == "warming"   ~ "Warming",
      Warm == "unwarming" ~ "Control",
      TRUE                ~ Warm
    ),
    
    Treatment = factor(
      Treatment,
      levels = c(
        "Control",
        "Warming"
      )
    ),
    
    year = factor(
      year,
      levels = 2009:2020
    )
  )



year_cols <- c(
  "2009" = "#3F4A6D",
  "2010" = "#50678A",
  "2011" = "#6083A1",
  "2012" = "#6A9EAD",
  "2013" = "#5AA5A4",
  "2014" = "#4E9D8D",
  "2015" = "#6FAE91",
  "2016" = "#92B486",
  "2017" = "#B7AD73",
  "2018" = "#D09A67",
  "2019" = "#D87A62",
  "2020" = "#BC5C5B"
)

p_dca_year <- ggplot(
  dca_df,
  aes(
    x = DCA1,
    y = DCA2,
    color = year,
    shape = Treatment
  )
) +
  
  geom_point(
    size = 2.8,
    alpha = 0.9
  ) +
  
  scale_color_manual(
    values = year_cols,
    name = "Year"
  ) +
  
  scale_shape_manual(
    values = c(
      "Control" = 16,
      "Warming" = 17
    ),
    name = "Treatment"
  ) +
  
  labs(
    x = "DCA1",
    y = "DCA2"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  guides(
    color = guide_legend(
      ncol =1,
      byrow = TRUE,
      order = 2
    ),
    
    shape = guide_legend(
      order = 1
    )
  ) +
  
  theme(
    axis.text = element_text(
      color = "black"
    ),
    
    axis.title = element_text(
      color = "black"
    ),
    
    legend.position = "right"
  )

p_dca_year

centroid_df <- dca_df %>%
  group_by(
    year,
    Treatment
  ) %>%
  summarise(
    DCA1 = mean(DCA1),
    DCA2 = mean(DCA2),
    .groups = "drop"
  ) %>%
  mutate(
    year_num = as.numeric(
      as.character(year)
    )
  )



p_traj <- ggplot() +
  
  # 原始样品点
  geom_point(
    data = dca_df,
    aes(
      x = DCA1,
      y = DCA2,
      color = Treatment
    ),
    size = 1.8,
    alpha = 0.25
  ) +
  
  # 年均轨迹
  geom_path(
    data = centroid_df,
    aes(
      x = DCA1,
      y = DCA2,
      color = Treatment,
      group = Treatment
    ),
    linewidth = 1.0,
    arrow = arrow(
      length = unit(0.15, "cm"),
      type = "open"
    )
  ) +
  
  # 年均点
  geom_point(
    data = centroid_df,
    aes(
      x = DCA1,
      y = DCA2,
      fill = year,
      color = Treatment
    ),
    shape = 21,
    size = 3.5,
    stroke = 0.8
  ) +
  
  scale_color_manual(
    values = c(
      "Control" = "#4C78A8",
      "Warming" = "#D96B5F"
    ),
    name = "Treatment"
  ) +
  
  scale_fill_manual(
    values = year_cols,
    name = "Year"
  ) +
  
  labs(
    x = "DCA1",
    y = "DCA2"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    axis.text = element_text(
      color = "black"
    ),
    
    axis.title = element_text(
      color = "black"
    ),
    
    legend.position = "right"
  ) + theme_bw()

p_traj
ggsave("Figure DCA.pdf", p_traj, width =16, height = 12, units="cm")













######### geochip proble level plot
library(tidyverse)

# 读取数据
dat <- read.csv("/Users/data from Xue/soil-borne pathogens/ResponseRatio/GeoChip_probe_overall_response_2009_2020.csv",header = TRUE,check.names = FALSE)
# 如果列名中存在转义符号，统一清理
names(dat) <- gsub("\\\\_", "_", names(dat))


#==========================================================
# 2. 设置基因的排列顺序
#    以下顺序对应图中从下到上
#==========================================================

gene_order <- c(
  "Cyanide_hydratase_Fungi",
  "Scytalone_dehydratase_Fungi",
  "Catalase_KatG_fungi",
  "Metalloprotease_fungi",
  "Superoxide_dismutase_fungi",
  "Conidial_laccase",
  "Patulin_6MSAS"
)


#==========================================================
# 3. 设置不同基因的颜色
#==========================================================

gene_colors <- c(
  "Cyanide_hydratase_Fungi"     = "#4E79A7",
  "Scytalone_dehydratase_Fungi" = "#76B7B2",
  "Catalase_KatG_fungi"         = "#59A14F",
  "Metalloprotease_fungi"       = "#F28E2B",
  "Superoxide_dismutase_fungi"  = "#B07AA1",
  "Conidial_laccase"            = "#E15759",
  "Patulin_6MSAS"               = "#EDC948"
)


#==========================================================
# 4. 整理绘图数据
#==========================================================

plot_dat <- dat %>%
  filter(
    !is.na(gene),
    !is.na(mean_lnRR),
    !is.na(se_lnRR)
  ) %>%
  mutate(
    gene = factor(gene, levels = gene_order),
    
    # 计算 mean ± SE
    lower_se = mean_lnRR - se_lnRR,
    upper_se = mean_lnRR + se_lnRR,
    
    # 显著性类别
    significant = Response != "NS"
  ) %>%
  filter(!is.na(gene)) %>%
  
  # 先按基因排列，再在每个基因内按 mean_lnRR 排列
  arrange(gene, mean_lnRR) %>%
  
  # 不同基因之间添加两个单位的间隔
  mutate(
    y = row_number() + (as.integer(gene) - 1) * 2
  )


#==========================================================
# 5. 计算每个基因的标签位置和数量
#==========================================================

gene_labels <- plot_dat %>%
  group_by(gene) %>%
  summarise(
    y = mean(range(y)),
    
    # 显著受抑制的探针数量
    n_suppressed = sum(Response == "Suppressed", na.rm = TRUE),
    
    # 显著增强的探针数量
    n_stimulated = sum(Response == "Stimulated", na.rm = TRUE),
    
    # 总探针数量
    n_total = n(),
    
    .groups = "drop"
  ) %>%
  mutate(
    clean_gene = as.character(gene) %>%
      str_replace_all("_Fungi|_fungi", "") %>%
      str_replace_all("_", " "),
    
    # 标签显示：显著受抑制探针数/总探针数
    label = paste0(
      clean_gene,
      " (", n_suppressed, "/", n_total, ")"
    )
  )


#==========================================================
# 6. 计算总体统计信息
#==========================================================

n_suppressed <- sum(
  plot_dat$Response == "Suppressed",
  na.rm = TRUE
)

n_stimulated <- sum(
  plot_dat$Response == "Stimulated",
  na.rm = TRUE
)

n_total <- nrow(plot_dat)

suppressed_pct <- round(
  100 * n_suppressed / n_total,
  digits = 1
)


#==========================================================
# 7. 自动确定横坐标和标签位置
#==========================================================

data_min <- min(plot_dat$lower_se, na.rm = TRUE)
data_max <- max(plot_dat$upper_se, na.rm = TRUE)

x_range <- data_max - data_min

# 右侧基因名称的位置
label_x <- data_max + 0.05 * x_range

# 给标签留出足够空间
x_lower <- data_min - 0.05 * x_range
x_upper <- data_max + 0.48 * x_range


#==========================================================
# 8. 绘图
#==========================================================

p <- ggplot(
  plot_dat,
  aes(
    y = y,
    colour = gene
  )
) +
  
  #--------------------------------------------------------
# mean lnRR ± SE 横向误差线
#--------------------------------------------------------
geom_segment(
  aes(
    x = lower_se,
    xend = upper_se,
    yend = y
  ),
  linewidth = 0.55,
  alpha = 0.30,
  lineend = "round"
) +
  
  #--------------------------------------------------------
# 连接同一基因内按 mean_lnRR 排序后的点
# 形成参考图中的“S形”轨迹
#--------------------------------------------------------
geom_path(
  aes(
    x = mean_lnRR,
    group = gene
  ),
  linewidth = 2.8,
  alpha = 0.65,
  lineend = "round"
) +
  
  #--------------------------------------------------------
# 所有探针的 mean lnRR
#--------------------------------------------------------
geom_point(
  aes(x = mean_lnRR),
  size = 2.2,
  alpha = 0.90
) +
  
  #--------------------------------------------------------
# 用黑色空心圆突出显著变化的探针
#--------------------------------------------------------
geom_point(
  data = plot_dat %>%
    filter(significant),
  aes(
    x = mean_lnRR,
    y = y
  ),
  inherit.aes = FALSE,
  shape = 21,
  fill = "white",
  colour = "black",
  stroke = 0.85,
  size = 3
) +
  
  #--------------------------------------------------------
# lnRR = 0 参考线
#--------------------------------------------------------
geom_vline(
  xintercept = 0,
  linetype = "dashed",
  colour = "#E64B35",
  linewidth = 0.75
) +
  
  #--------------------------------------------------------
# 添加每个基因的名称和统计数量
#--------------------------------------------------------
geom_text(
  data = gene_labels,
  aes(
    x = label_x,
    y = y,
    label = label,
    colour = gene
  ),
  inherit.aes = FALSE,
  hjust = 0,
  size = 3.7,
  fontface = "italic"
) +
  
  #--------------------------------------------------------
# 添加顶部总体统计信息
#--------------------------------------------------------
annotate(
  geom = "text",
  x = data_min,
  y = max(plot_dat$y) + 6,
  hjust = 0,
  vjust = 1,
  size = 3.8,
  label = paste0(
    "Significantly suppressed / total =\n",
    n_suppressed, "/", n_total,
    " (", suppressed_pct, "%)"
  )
) +
  
  #--------------------------------------------------------
# 设置颜色
#--------------------------------------------------------
scale_colour_manual(
  values = gene_colors,
  drop = FALSE
) +
  
  #--------------------------------------------------------
# 坐标轴
#--------------------------------------------------------
scale_x_continuous(
  breaks = scales::pretty_breaks(n = 6)
) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.14))
  ) +
  
  coord_cartesian(
    xlim = c(x_lower, x_upper),
    clip = "off"
  ) +
  
  labs(
    x = "Effect of warming (mean lnRR ± SE)",
    y = NULL
  ) +
  
  #--------------------------------------------------------
# 图形主题
#--------------------------------------------------------
theme_classic(base_size = 12) +
  
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    axis.title.x = element_text(
      size = 12,
      colour = "black",
      margin = margin(t = 8)
    ),
    
    axis.text.x = element_text(
      size = 10,
      colour = "black"
    ),
    
    legend.position = "none",
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.7
    ),
    
    # 右侧留出空间显示基因名称
    plot.margin = margin(
      t = 10,
      r = 175,
      b = 10,
      l = 10
    )
  ) +theme_bw()

p


#==========================================================
# 9. 保存图片
#==========================================================

ggsave("Figure response ratios.pdf", p, width =18, height = 20, units="cm")








