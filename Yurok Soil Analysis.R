### title: Yurok Soil Analysis
### author: Lisa Eash
### date started: 20240930
### date updated: 20241015

#Hi Soils team!#
#Comment 2#

####Load packages####
    library(reshape)
####Read in data####
    df23.12<-read.csv("./YurokSoil_Ward_Jan24.csv")
####Clean and merge datasets####
    #Remove columns
    head(df23.10)
    unique(df23.10$Sample.ID)
    df23.10<-df23.10[,-c(1,2,4)]
    head(df23.12)
    df23.12<-df23.12[,-c(1,2,4)]
    head(dfb23.10)
    dfb23.10<-dfb23.10[,-c(1:13)]
    head(dfb23.12)
    dfb23.12<-dfb23.12[,-c(1:13)]
    #Merge df
    dfs<-rbind(data.frame(Time="Pre",df23.10),
               data.frame(Time="Post",df23.12))
    dfb<-rbind(data.frame(Time="Pre",dfb23.10),
               data.frame(Time="Post",dfb23.12))
    #Clean column names
    soil.col<-c("Time","Sample.ID","B.Depth","E.Depth","SM","Sand",
                "Silt","Clay","Text","SOC","Dry.Weight","BD")
    colnames(dfs)<-soil.col
    bio.col<-c("Time","Sample.ID","Biomass","Div.Index","Bacteria.perc","Bact.biomass",
               "Actino.perc","Actino.biomass","Gramneg.perc","Gramneg.biomass",
               "Rhizobia.perc","Rhizobia.biomass","Fungi.perc","Fungi.biomass",
               "AMF.perc","AMF.biomass","Sapro.perc","Sapro.biomass","Proto.perc",
               "Proto.biomass","Grampos.biomass","Grampos.perc","Undiff.perc","Undiff.biomass",
               "Fungi.bact.ratio","Predator.prey.ratio","Grampos.gramneg.ratio","Sat.unsat.ratio",
               "Mono.poly.ratio","Pre16.1w7c.cy17","Pre18.1w7c.cy19","B.Depth","E.Depth","Ace.gkg")
    colnames(dfb)<-bio.col
    head(dfb)
    #Merge depths
    str(dfs)
    dfs$Depth<-NA
    dfs$Depth<-ifelse(dfs$E.Depth==20,"Deep","Shallow")
    dfs$SOC<-as.numeric(dfs$SOC)
    dfs<-ddply(dfs,c("Time","Depth","Sample.ID"),summarize,
               SOC=mean(SOC,na.rm=TRUE),
               BD=mean(BD,na.rm=TRUE),
               Clay=mean(Clay,na.rm=TRUE),
               Silt=mean(Silt,na.rm=TRUE),
               Sand=mean(Sand,na.rm=TRUE),
               SM=mean(SM,na.rm=TRUE),
               Dry.weight=mean(Dry.Weight,na.rm=TRUE))
    dfs<-dfs[!is.na(dfs$Time),]
    #Define treatments
    dfs$Treatment<-ifelse(dfs$Sample.ID<200,"T","C")
    dfb$Treatment<-ifelse(dfb$Sample.ID<200,"T","C")

#### Soil property analysis ####
#Remove SOC outlier
dfs<-dfs[!dfs$SOC>15,]
#Calculate SOC stock and SOC change over time
tail(dfs)
dfs$SOC_stock<-dfs$SOC/100*dfs$BD*10/ #cm shallow depth
  1000000* # g to t
  100000000 #cm2 to ha
plot(dfs$SOC)
df_wide<-reshape(dfs[,c("Treatment","Time","Depth","Sample.ID","SOC","SOC_stock")],
                 idvar = c("Treatment","Sample.ID","Depth"),
                 timevar = "Time", direction = "wide")
df_wide$SOC_stock_diff<-df_wide$SOC_stock.Post-df_wide$SOC_stock.Pre
plot(df_wide$SOC_stock_diff)
df_wide$SOC_diff<-df_wide$SOC.Post-df_wide$SOC.Pre
plot(df_wide$SOC_diff)
#Create summary table with mean and sd
sum_soil<-plyr::ddply(dfs, c("Treatment","Time","Depth"),summarise,
                      SOC_mean=mean(SOC, na.rm=TRUE),
                      SOC_se=sd(SOC, na.rm=TRUE)/sqrt(length(!is.na(SOC))),
                      SOC_stock_mean=mean(SOC_stock,na.rm=TRUE),
                      SOC_stock_se=sd(SOC_stock,na.rm=TRUE)/sqrt(length(!is.na(SOC_stock))),
                      BD_mean=mean(BD, na.rm=TRUE),
                      BD_se=sd(BD, na.rm=TRUE)/sqrt(length(!is.na(BD))),
                      SM_mean=mean(SM, na.rm=TRUE),
                      SM_se=sd(SM, na.rm=TRUE)/sqrt(length(!is.na(SM))),
                      Clay_mean=mean(Clay,na.rm=TRUE),
                      Clay_se=sd(Clay, na.rm=TRUE)/sqrt(length(!is.na(Clay))),
                      Silt_mean=mean(Silt, na.rm=TRUE),
                      Silt_se=sd(Silt, na.rm=TRUE)/sqrt(length(!is.na(Silt))),
                      Sand_mean=mean(Sand, na.rm=TRUE),
                      Sand_se=sd(Sand, na.rm=TRUE)/sqrt(length(!is.na(Sand))))
sum_soil<-sum_soil[order(sum_soil$Depth,sum_soil$Time,sum_soil$Treatment),]
sum_soil
sum_soil<-sum_soil[!is.na(sum_soil$Treatment),]
sum_soc<-plyr::ddply(df_wide,c("Treatment","Depth"),summarise,
                     SOCstock_mean=mean(SOC_stock_diff,na.rm=TRUE),
                     SOCstock_sd=sd(SOC_stock_diff,na.rm=TRUE),
                     SOC_mean=mean(SOC_diff,na.rm=TRUE),
                     SOC_sd=mean(SOC_diff,na.rm=TRUE))
#Assess variables for normality
norm_test<-function(var){
  res<-data.frame(p.val=NA,transf=NA)
  res[1,]<-shapiro.test(var)$p.value
  res[2,]<-shapiro.test(log(var))$p.value
  res[3,]<-shapiro.test(exp(var))$p.value
  res[4,]<-shapiro.test(sqrt(var))$p.value
  res[5,]<-shapiro.test(1/(var))$p.value
  res[,2]<-c("None","log","exp","sqrt","recip")
  max_p<-as.numeric(max(res[,1]))
  transf<-res[res$p.val==max_p,]$transf
  row<-data.frame(max_p,transf)
  return(row)
}
trans_df<-data.frame(Var=NA,P=NA,Transf=NA)
SOC_norm<-norm_test(dfs$SOC)
trans_df[1,]<-c("SOC",SOC_norm)
BD_norm<-norm_test(dfs$BD)
trans_df[2,]<-c("BD",BD_norm)
Clay_norm<-norm_test(dfs$Clay)
trans_df[3,]<-c("Clay",Clay_norm)
Sand_norm<-norm_test(dfs$Sand)
trans_df[4,]<-c("Sand",Sand_norm)
Silt_norm<-norm_test(dfs$Silt)
trans_df[5,]<-c("Silt",Silt_norm)
SM_norm<-norm_test(dfs[dfs$Time=="Post",]$SM)
trans_df[6,]<-c("SM",SM_norm)
SOCdiff_norm<-norm_test(df_wide$SOC_diff)
trans_df[7,]<-c("SOCdiff",SOCdiff_norm)
trans_df
####ANOVA tests####
library(emmeans)
#SOC
ANOVA.SOC<-aov((1/SOC) ~ Treatment*Depth, data=dfs[dfs$Time=="Post",])
summary(ANOVA.SOC)
emmeans(ANOVA.SOC, ~Treatment|Depth)
#SOC avg by depth
head(df23.10)
dfs$inches<-NA
dfs$inches<-ifelse(dfs$Depth=="Shallow",4,16)
avgSOC<-ddply(dfs,c("Time","Treatment","Sample.ID"),summarise,
              SOC=weighted.mean(SOC,inches))
ANOVA.SOCavg<-aov((SOC) ~ Treatment, data=avgSOC[avgSOC$Time=="Post",])
summary(ANOVA.SOCavg)
emmeans(ANOVA.SOCavg, ~Treatment)
#SOC diff
ANOVA.SOCdiff<-aov(SOC_stock_diff ~ Treatment, data=df_wide)
summary(ANOVA.SOCdiff)
emmeans(ANOVA.SOCdiff, ~Treatment)
#BD
ANOVA.BD<-aov(BD ~ Time*Treatment, data=dfs)
summary(ANOVA.BD)
emmeans(ANOVA.BD, ~Time|Treatment)
#Clay
ANOVA.Clay<-aov(Clay ~ Treatment*Depth, data=dfs[dfs$Time=="Post",])
summary(ANOVA.Clay)
emmeans(ANOVA.Clay, ~Time|Treatment)
#Avg by depth
avgClay<-ddply(dfs,c("Time","Treatment","Sample.ID"),summarise,
               Clay=weighted.mean(Clay,inches,na.rm=TRUE))
ANOVA.Clayavg<-aov((Clay) ~ Treatment*Time, data=avgClay)
summary(ANOVA.Clayavg)
emmeans(ANOVA.Clayavg, ~Treatment)
#Sand
ANOVA.Sand<-aov(Sand ~ Time*Treatment + Treatment*Depth, data=dfs)
summary(ANOVA.Sand)
emmeans(ANOVA.Sand, ~Time|Treatment)
#Silt
ANOVA.Silt<-aov(Silt^2 ~ Time*Treatment + Time*Depth, data=dfs)
summary(ANOVA.Silt)
emmeans(ANOVA.Silt, ~Treatment|Time)
#SM
ANOVA.SM<-aov(SM ~ Treatment*Time, data=dfs)
summary(ANOVA.SM)
emmeans(ANOVA.SM, ~Time|Treatment)
####Plots####
#SOC
ggplot(sum_soil, aes(x=Depth, y=SOC_mean, fill=Treatment)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  facet_wrap(~factor(Time,levels=c("Pre","Post")))+
  geom_errorbar(aes(ymin=SOC_mean-SOC_se, ymax=SOC_mean+SOC_se), width=.2,position=position_dodge(.9))
ggplot(sum_soil[sum_soil$Time=="Post" & sum_soil$Depth=="Shallow",], 
       aes(x=Depth, y=SOC_stock_mean, fill=Treatment)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  #facet_wrap(~Depth)+
  geom_errorbar(aes(ymin=SOC_stock_mean-SOC_stock_sd, ymax=SOC_stock_mean+SOC_stock_sd), width=.2,position=position_dodge(.9))
#BD
ggplot(sum_soil, aes(x=factor(Time,levels=c("Pre","Post")), y=BD_mean, fill=Treatment)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=BD_mean-BD_se, ymax=BD_mean+BD_se), width=.2,position=position_dodge(.9))
#SM
ggplot(sum_soil, aes(x=factor(Time,levels=c("Pre","Post")), y=SM_mean, fill=Treatment)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=SM_mean-SM_se, ymax=SM_mean+SM_se), width=.2,position=position_dodge(.9))
#Clay
ggplot(sum_soil, aes(x=Depth, y=Clay_mean, fill=factor(Time,levels=c("Pre","Post")))) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  facet_wrap(~Treatment)+
  geom_errorbar(aes(ymin=Clay_mean-Clay_se, ymax=Clay_mean+Clay_se), width=.2,position=position_dodge(.9))
#Avg Clay
plot_clay<-ddply(avgClay,c("Time","Treatment"),summarise,
                 Clay_mean=mean(Clay,na.rm=TRUE),
                 Clay_se=sd(Clay,na.rm=TRUE)/sqrt(length(!is.na(Clay))))
ggplot(plot_clay, aes(x=factor(Time,levels=c("Pre","Post")), y=Clay_mean, fill=Treatment)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=Clay_mean-Clay_se, ymax=Clay_mean+Clay_se), width=.2,position=position_dodge(.9))
#Sand
ggplot(sum_soil, aes(x=Depth, y=Sand_mean, fill=factor(Time,levels=c("Pre","Post")))) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  facet_wrap(~Treatment)+
  geom_errorbar(aes(ymin=Sand_mean-Sand_se, ymax=Sand_mean+Sand_se), width=.2,position=position_dodge(.9))
#Silt
ggdistr(dfs$Silt)
ggdensity(dfs, x = "Silt", fill = "lightgray", title = "CONT")
shapiro.test((dfs$Silt)^2)
ggplot(sum_soil, aes(x=Depth, y=Silt_mean, fill=factor(Time,levels=c("Pre","Post")))) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  facet_wrap(~Treatment)+
  geom_errorbar(aes(ymin=Silt_mean-Silt_se, ymax=Silt_mean+Silt_se), width=.2,position=position_dodge(.9))
####Soil Biology Analysis####
head(dfb)
####Summary table####
dfb$Pre18.1w7c.cy19<-as.numeric(dfb$Pre18.1w7c.cy19)
dfb$Pre1816Ratio<-dfb$Pre16.1w7c.cy17/dfb$Pre18.1w7c.cy19
temp<-dfb[dfb$E.Depth %in% c("2","4"),]
sum_bio<-plyr::ddply(dfb, c("Treatment","Time"),summarise,
                     Totbio_mean=mean(Biomass, na.rm=TRUE),
                     Totbio_se=sd(Biomass, na.rm=TRUE)/sqrt(length(!is.na(Biomass))),
                     Fungbio_mean=mean(Fungi.biomass,na.rm=TRUE),
                     Fungbio_se=sd(Fungi.biomass,na.rm=TRUE)/sqrt(length(!is.na(Fungi.biomass))),
                     Bactbio_mean=mean(Bact.biomass,na.rm=TRUE),
                     Bactbio_se=sd(Bact.biomass,na.rm=TRUE)/sqrt(length(!is.na(Bact.biomass))),
                     Div_mean=mean(Div.Index, na.rm=TRUE),
                     Div_se=sd(Div.Index, na.rm=TRUE)/sqrt(length(!is.na(Div.Index))),
                     FBRat_mean=mean(Fungi.bact.ratio, na.rm=TRUE),
                     FBRat_se=sd(Fungi.bact.ratio, na.rm=TRUE)/sqrt(length(!is.na(Fungi.bact.ratio))),
                     GPGRRat_mean=mean(Grampos.gramneg.ratio,na.rm=TRUE),
                     GPGRRat_se=sd(Grampos.gramneg.ratio, na.rm=TRUE)/sqrt(length(!is.na(Grampos.gramneg.ratio))),
                     Pre1816Rat_mean=mean(Pre1816Ratio, na.rm=TRUE),
                     Pre1816Rat_se=sd(Pre1816Ratio, na.rm=TRUE)/sqrt(length(!is.na(Pre1816Ratio))),
                     Actino_mean=mean(Actino.biomass, na.rm=TRUE),
                     Actino_se=sd(Actino.biomass, na.rm=TRUE)/sqrt(length(!is.na(Actino.biomass))),
                     Rhizobia_mean=mean(Rhizobia.biomass, na.rm=TRUE),
                     Rhizobia_se=sd(Rhizobia.biomass, na.rm=TRUE)/sqrt(length(!is.na(Rhizobia.biomass))),
                     AMF_mean=mean(AMF.biomass, na.rm=TRUE),
                     AMF_se=sd(AMF.biomass, na.rm=TRUE)/sqrt(length(!is.na(AMF.biomass))),
                     Sapro_mean=mean(Sapro.biomass, na.rm=TRUE),
                     Sapro_se=sd(Sapro.biomass, na.rm=TRUE)/sqrt(length(!is.na(Sapro.biomass))),
                     Proto_mean=mean(Proto.biomass, na.rm=TRUE),
                     Proto_se=sd(Proto.biomass, na.rm=TRUE)/sqrt(length(!is.na(Proto.biomass))))
plot_bio<-reshape(sum_bio, direction='long', 
                  varying=c('Totbio_mean','Totbio_se','Fungbio_mean',
                            'Fungbio_se','Bactbio_mean','Bactbio_se',
                            'Div_mean','Div_se','FBRat_mean','FBRat_se',
                            'GPGRRat_mean','GPGRRat_se','Pre1816Rat_mean',
                            'Pre1816Rat_se','Actino_mean','Actino_se','Rhizobia_mean','Rhizobia_se',
                            'AMF_mean','AMF_se',
                            'Sapro_mean','Sapro_se','Proto_mean','Proto_se'), 
                  timevar='var',
                  times=c('Totbio', 'Fungbio','Bactbio','Div','FBRat',
                          'GPGRRat','Pre1618Rat','Actino','Rhizobia','AMF',
                          'Sapro','Proto'),
                  v.names=c('mean', 'se'),
                  idvar=c('Treatment','Time'))

####Bar charts####
bact_fung<-dfb[,!names(dfb) %in% c("Undiff.perc","Undiff.biomass","Predator.prey.ratio",
                                   "Sat.unsat.ratio","Mono.poly.ratio","Pre16.1w7c.cy17",
                                   "Pre18.1w7c.cy19","B.Depth","Ace.gkg")]
head(bact_fung)
bact_fung<-bact_fung %>% 
  group_by(Time,Treatment, E.Depth) %>%
  dplyr::summarise_all( mean) %>%
  as.data.frame()
bact_fung_sd<-bact_fung %>% 
  group_by(Time,Treatment, E.Depth) %>%
  dplyr::summarise_all( sd) %>%
  as.data.frame()
bact_fung <- bact_fung %>% 
  pivot_longer(
    cols = `Biomass`:`Grampos.gramneg.ratio`, 
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  as.data.frame()
plot_bio$E.Depth<-as.factor(plot_bio$E.Depth)
head(plot_bio)
p_bf<-ggplot(plot_bio[
  plot_bio$Time=="Post" &
    plot_bio$var %in% 
    c(
      #"Fungbio",
      #"Bactbio"
      #"Bactbio"
      "Actino",
      #"Rhizobia"
      "AMF",
      "Sapro"
      #"Proto"
      #"Pre1618Rat"
      #"Bacteria.perc","Fungi.perc",
      #"Gramneg.perc","Grampos.perc",
      #"Rhizobia.perc",
      #"AMF.perc",
      #"Proto.perc",
      #"Sapro.perc"
    ),],
  aes(x=var, y=mean
      , fill=Treatment
  )) +
  geom_bar(position=position_dodge(), stat="identity", colour='black')+
  #facet_wrap(var)+
  ylab("Biomass (ng/g)")+
  xlab("Microbial Taxa")+
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2,position=position_dodge(.9))
p_bf
#Diversity
div_df<-dfb[dfb$E.Depth %in% c("2","4"),]
div_sum<-ddply(div_df,c("Time","Treatment"),summarise,
               Div_mean=mean(Div.Index, na.rm=TRUE),
               Div_se=sd(Div.Index, na.rm=TRUE)/sqrt(length(!is.na(Div.Index))))
p_bf<-ggplot(div_sum,
             aes(x=factor(Time,levels=c("Pre","Post")), y=Div_mean
                 , colour=Treatment
             )) +
  geom_point(position=position_dodge(.6))+
  #facet_wrap(Treatment)+
  ylab("Diversity Index")+
  xlab("Sampling Time")+
  geom_errorbar(aes(ymin=Div_mean-Div_se, ymax=Div_mean+Div_se), width=.2,position=position_dodge(.6))
p_bf
####Transformations####
trans_df<-data.frame(Var=NA,P=NA,Transf=NA)
Bio_norm<-shapiro.test(dfb$Biomass) #Biomass normal
Fungbio_norm<-shapiro.test(dfb$Fungi.biomass) #Fung bio normal
Bactbio_norm<-shapiro.test(log(dfb$Bact.biomass)) #Bact bio log
FB_norm<-shapiro.test(sqrt((dfb$Fungi.bact.ratio/max(dfb$Fungi.bact.ratio)))) #Bact bio log
####ANOVAS####
ANOVA.Bio<-aov(Biomass ~ Treatment*E.Depth, data=dfb[dfb$Time=="Post",])
summary(ANOVA.Bio)
emmeans(ANOVA.Bio, ~Time|Treatment)
ANOVA.FBio<-aov(Fungi.biomass ~ E.Depth*Treatment, data=dfb[dfb$Time=="Post",])
summary(ANOVA.FBio)
emmeans(ANOVA.Bio, ~Time|Treatment)
ANOVA.BBio<-aov(Bact.biomass ~ Time*Treatment + Treatment*E.Depth, data=dfb)
summary(ANOVA.BBio)
emmeans(ANOVA.BBio, ~Treatment|Time)
ANOVA.Div<-aov(asin((dfb$Div.Index/max(dfb$Div.Index))) ~ Time*Treatment + Time*E.Depth, data=dfb)
summary(ANOVA.Div)
emmeans(ANOVA.Div, ~Treatment|Time)
range(dfb$Div.Index)
shapiro.test(asin((dfb$Div.Index/max(dfb$Div.Index))))
ggdensity(dfb, x = "Div.Index", fill = "lightgray", title = "CONT")

####NMDS####
#Restructure df for NMDS
NMDS_df<-dfb[,c("Time","Treatment","E.Depth","Sample.ID","Bacteria.perc",
                "Actino.perc","Gramneg.perc","Rhizobia.perc",
                "Fungi.perc","AMF.perc","Sapro.perc","Proto.perc",
                "Grampos.perc")]
head(NMDS_df)
NMDS_df<-NMDS_df[NMDS_df$Time=="Post",]
#NMDS_shallow<-NMDS_shallow %>% 
#  group_by(Time,Treatment, Sample.ID) %>%
#  dplyr::summarise_all( mean) %>%
#  as.data.frame()
NMDS_df$Sample.Lab<-paste0(NMDS_df$Time,"_",NMDS_df$Treatment,"_",
                           NMDS_df$E.Depth,"_",NMDS_df$Sample.ID)
row.names(NMDS_df)<-NMDS_df$Sample.Lab
Treatment<-NMDS_df[,names(NMDS_df) %in% c("Time","Treatment","E.Depth","Sample.ID")]
NMDS_df<-NMDS_df[,!names(NMDS_df) %in% c("Sample.Lab","Time","Treatment","Depth","E.Depth","Sample.ID")]
NMDS_colmax<-sapply(NMDS_df,max)
NMDS_relcolmax <- data.frame(sweep(NMDS_df,2,NMDS_colmax,'/'))
#Restructure dfs for NMDS
head(dfs)
env<-dfs[,c("Time","Depth","Treatment","Sample.ID","SOC","Clay","Silt")]
head(env)
env$Sample.Lab<-paste0(env$Time,"_",env$Treatment,"_",
                       env$Depth,"_",env$Sample.ID)
row.names(env)<-env$Sample.Lab
env<-env[,!names(env) %in% c("Time","Treatment","Depth","Sample.ID","Sample.Lab")]
env[is.na(env$Silt),]
nrow(env)
nrow(NMDS_df)
#Standardize and calculate distance
library(vegan)
head(NMDS_relcolmax)
spp.bcd <- vegdist(NMDS_relcolmax)
spp.mds<-metaMDS(NMDS_relcolmax, trace = FALSE,autotransform=F,maxit=1000,k=3,engine="monoMDS")
ordiplot(spp.mds
         #,xlim=c(-.75,.75)
         #,ylim=c(-.4,.3)
)
spp.mds 
stressplot(spp.mds, spp.bcd)
orditorp(spp.mds, display = "species",  scaling = scl,
         col = "red", pch = 2, cex = 0.8, xlim=c(-1.5,1.5), ylim=c(-.5,.5))
ordiellipse(spp.mds, Treatment$Treatment, kind = c("se"), conf=0.95,
            draw = c("polygon"),
            #col = colvec,
            alpha = 127, label = TRUE)
ordiellipse(spp.mds, Treatment$Treatment, kind = c("se"), conf=0.95,
            draw = c("polygon"),
            #col = colvec,
            alpha = 127, label = TRUE)
str(spp.mds
)
data.scores <- as.data.frame(scores(spp.mds)$sites)  #Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores$Site.Label <- rownames(data.scores)  # create a column of site names, from the rownames of data.scores
Treatment$Site.Label<-rownames(Treatment)
data<-merge(data.scores,Treatment,by=c("Site.Label"))
species.scores <- as.data.frame(scores(spp.mds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
species.scores<-species.scores[species.scores$species %in% c("Rhizobia.perc","Proto.perc","Sapro.perc","Gramneg.perc","AMF.perc"),]
species.scores$species<-c("Gramneg","Rhizobia","AMF","Sapro","Proto")
head(species.scores)  #look at the data
head(data)
data_temp<-data[data$Time=="Pre",]
grp.a <- data_temp[data_temp$Treatment == "T", ][chull(data_temp[data_temp$Treatment == 
                                                                   "T", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp.b <- data_temp[data_temp$Treatment == "C", ][chull(data_temp[data_temp$Treatment == 
                                                                   "C", c("NMDS1", "NMDS2")]), ]  # hull values for grp B

hull.data_temp <- rbind(grp.a, grp.b)  #combine grp.a and grp.b
hull.data_temp
data_temp$E.Depth<-as.factor(data_temp$E.Depth)
ggplot() + 
  geom_polygon(data=hull.data_temp,aes(x=NMDS1,y=NMDS2,fill=Treatment,group=Treatment),alpha=0.30) + 
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=0.5) +  # add the species labels
  geom_point(data=data_temp,aes(x=NMDS1,y=NMDS2,shape=E.Depth,colour=Treatment),size=2) + # add the point markers
  #geom_text(data=data,aes(x=NMDS1,y=NMDS2,label=Treatment),size=6,vjust=0) +  # add the site labels
  #scale_colour_manual(values=c("A" = "red", "B" = "blue")) +
  coord_equal() +
  theme_bw()
data.scores$grp <- grp  #  add the grp variable created earlier
head(data.scores)  #look at the data
head(Treatment)
Treatment$Trt_Time<-paste0(Treatment$Treatment,
                           "_",Treatment$Time)
text(spp.mds, display = "species", cex=0.7, col="red")
#PERMANOVA
permanovamacro<-adonis2(spp.bcd ~ Treatment*Time, data=Treatment, perm=999)
(permanovamacro)
#indicator species
library(indicspecies)
Treatment$Trt_Time<-paste0(Treatment$Treatment,"_",Treatment$Time)
isa<-multipatt(NMDS_relcolmax,Treatment$Treatment,duleg=T) #duleg=T skips calcluating indicators for combinations of groups
#isa.5 is on cluster analysis with 5 groups retained
isa_pre<-summary(isa) #shortest output
#each species assigned to the one group that it's best for
summary(isa,alpha=.2) #shows A (specificity) and B (fidelity), essentially how indicator value was calculated
summary(isa.5, indvalcomp=TRUE, alpha=1) #longest output, shows indicator values in "stat" column, alpha controls which displayed


str(spp.mds)
#Environmental Data for dbRDA
df<-read.csv("AllData_v3.csv")
df$Block<-as.factor(df$Block)
df$TrtNo<-as.factor(df$TrtNo)
#Remove control trt#
dfnoctrl<-df[df$TrtNo!="1",]
str(dfnoctrl)
Env<-df[,c("Pot.ID","LAP",
           "TAP",
           "NAG","BG","PHOS",
           "CB","RootBiomass","AbvBiomass",
           "TotalBiomass","POXC")]
plotnames<-Env[,1]
rownames(Env)<-plotnames 
Env<-Env[,-1]
Env<-subset(Env, !is.na(LAP))
Env<-subset(Env, !is.na(TAP))
Env<-subset(Env, !is.na(NAG))
Env<-subset(Env, !is.na(BG))
Env<-subset(Env, !is.na(PHOS))
Env<-subset(Env, !is.na(CB))
Env<-subset(Env, !is.na(POXC))
env.data.z<-Env
env.data.z$LAP <- (env.data.z$LAP - mean(env.data.z$LAP))/sd(env.data.z$LAP)
env.data.z$TAP <- (env.data.z$TAP - mean(env.data.z$TAP))/sd(env.data.z$TAP)
env.data.z$NAG <- (env.data.z$NAG - mean(env.data.z$NAG))/sd(env.data.z$NAG)
env.data.z$BG <- (env.data.z$BG - mean(env.data.z$BG))/sd(env.data.z$BG)
env.data.z$PHOS <- (env.data.z$PHOS - mean(env.data.z$PHOS))/sd(env.data.z$PHOS)
env.data.z$CB <- (env.data.z$CB - mean(env.data.z$CB))/sd(env.data.z$CB)
env.data.z$POXC <- (env.data.z$POXC - mean(env.data.z$POXC))/sd(env.data.z$POXC)
env.data.z$RootBiomass <- (env.data.z$RootBiomass - mean(env.data.z$RootBiomass))/sd(env.data.z$RootBiomass)
env.data.z$AbvBiomass <- (env.data.z$AbvBiomass - mean(env.data.z$AbvBiomass))/sd(env.data.z$AbvBiomass)
env.data.z$TotalBiomass <- (env.data.z$TotalBiomass - mean(env.data.z$TotalBiomass))/sd(env.data.z$TotalBiomass)
#Exclude samples 307 and 308 and 405
row_names_df_to_remove<-c("G307","G308","G405")
env.data.z<-env.data.z[!(row.names(env.data.z) %in% row_names_df_to_remove),]
