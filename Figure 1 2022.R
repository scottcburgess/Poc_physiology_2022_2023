### This code produces Figure #, and associated analyses, in: 
### Edmunds PJ, Burgess SC. TITLE HERE.
# Code written by Scott Burgess, sburgess@bio.fsu.edu

# R version 4.2.3 (2023-03-15) -- "Shortstop Beagle"

# Load required libraries
library('dplyr')
library('glmmTMB')
library('AICcmodavg')
library('emmeans')
library('car')
library('DHARMa')

# set the colors for plotting
cols <- cbind.data.frame(Species = factor(c("P. tuahiniensis",
                                            "P. meandrina",
                                            "P. verrucosa",
                                            "P. grandis",
                                            "P. cf. effusa",
                                            "P. acuta")),
                         Color = c("#D55E00",
                                   "#0072B2",
                                   "#E69F00",
                                   "#56B4E9",
                                   "#009E73",
                                   "#e63946"))


# Import data
dat_2022 <- read.csv("2022_Experiment.csv")

# Ensure Mean.Light is a numerica variable
dat_2022$Mean.Light <- as.numeric(dat_2022$Mean.Light)

dat_2022 %>% 
  group_by(Species) %>% 
  summarise(n=n()) %>% 
  mutate(freq=n/sum(n))

tmp <- read.csv("/Users/scottburgess/My Drive/Data/Corals/Pocillipora/2022 Pm vs hap 10 experiment/2022 Pm vs hap 10 experiment Species ID.csv")
nrow(tmp)
tmp %>% 
  group_by(Species) %>% 
  summarise(n=n()) %>% 
  mutate(freq=n/sum(n))


# Function to compare models using AICc for the 2022 data
Analyze2022 <- function(dat,yname){
  dat$yvar <- dat[,names(dat)==yname,]
  
  m1 <- glmmTMB(yvar ~ poly(Mean.Light,2)*Temperature*Species + (1|Tank), data=dat)
  m2 <- update(m1,.~.-Temperature:poly(Mean.Light,2):Species)
  m3 <- update(m2,.~.-Temperature:poly(Mean.Light,2))
  m4 <- update(m2,.~.-Temperature:Species)
  m5 <- update(m2,.~.-poly(Mean.Light,2):Species)
  m6 <- glmmTMB(yvar ~ Temperature+poly(Mean.Light,2)+Species+Temperature:poly(Mean.Light,2) + (1|Tank), data=dat)
  m7 <- glmmTMB(yvar ~ Temperature+poly(Mean.Light,2)+Species+Species:poly(Mean.Light,2) + (1|Tank), data=dat)
  m8 <- glmmTMB(yvar ~ Temperature+poly(Mean.Light,2)+Species+Temperature:Species + (1|Tank), data=dat)
  m9 <- glmmTMB(yvar ~ Temperature+poly(Mean.Light,2)+Species + (1|Tank), data=dat)
  aic_table <- aictab(list(m1=m1,m2=m2,m3=m3,m4=m4,m5=m5,m6=m6,m7=m7,m8=m8,m9=m9))
  Anova_table <- Anova(get(aic_table$Modnames[1]),type="III")
  
  pred <- expand.grid(Temperature=unique(dat$Temperature),
                      Mean.Light=seq(min(dat$Mean.Light),max(dat$Mean.Light),length.out=50),
                      Species=unique(dat$Species),
                      Tank=NA)
  p <- predict(get(aic_table$Modnames[1]),newdat=pred,se.fit=T)
  pred$fit <- p$fit
  pred$lwr <- p$fit - 2*p$se.fit
  pred$upr <- p$fit + 2*p$se.fit
  
  DHARMaOutput <- simulateResiduals(fittedModel = get(aic_table$Modnames[1]), plot = F)
  
  list(pred=pred,
       DHARMaOutput=DHARMaOutput,
       aic_table=aic_table)
}


# Get outputs for results and plotting 
# 2022 growth (P. meandrina, P. tuahiniensis)
y <- dat_2022 %>% filter(Species %in% c("P. meandrina", "P. tuahiniensis"))
with(y,table(Tank,Light, Temperature))
y_mean <- aggregate(cbind(Growth.mg.cm2.d,Respiration.n.mol.cm2.min,Photosynthesis.nmol.cm2.min) ~ Tank+Light+Mean.Light+Temperature+Species,FUN=mean,data=y)
y_mean$Respiration.n.mol.cm2.min <- abs(y_mean$Respiration.n.mol.cm2.min)
dat_2022_growth_output <- Analyze2022(dat=y,yname="Growth.mg.cm2.d")
dat_2022_respiration_output <- Analyze2022(dat=y,yname="Respiration.n.mol.cm2.min")
dat_2022_respiration_output$pred$fit <- abs(dat_2022_respiration_output$pred$fit)
dat_2022_respiration_output$pred$lwr <- abs(dat_2022_respiration_output$pred$lwr)
dat_2022_respiration_output$pred$upr <- abs(dat_2022_respiration_output$pred$upr)
dat_2022_photosynthesis_output <- Analyze2022(dat=y,yname="Photosynthesis.nmol.cm2.min")
# plot(dat_2022_growth_output$DHARMaOutput)
# plot(dat_2022_respiration_output$DHARMaOutput)
# plot(dat_2022_photosynthesis_output$DHARMaOutput)


# Growth - hypothesis tests of factors
# m8 had lowest AICc
m8 <- glmmTMB(Growth.mg.cm2.d ~ Temperature + poly(Mean.Light,2) + Species + Temperature:Species + (1|Tank), data=y)
m9 <- glmmTMB(Growth.mg.cm2.d ~ Temperature + poly(Mean.Light,2) + Species + (1|Tank), data=y)
m10 <- glmmTMB(Growth.mg.cm2.d ~ Temperature + Species + Temperature:Species + (1|Tank), data=y)
anova(m8,m9, test="Chisq")
anova(m8,m10, test="Chisq")
emmeans(m8,pairwise~Temperature:Species)
0.0513 + c(-2,2)*0.0432
0.1093 + c(-2,2)*0.0421

# Respiration - hypothesis tests of factors
# m6 had lowest AICc
m6 <- glmmTMB(Respiration.n.mol.cm2.min ~ Temperature + poly(Mean.Light,2) + Species + Temperature:poly(Mean.Light,2) + (1|Tank), data=y)
m9 <- glmmTMB(Respiration.n.mol.cm2.min ~ Temperature + poly(Mean.Light,2) + Species + (1|Tank), data=y)
m10 <- glmmTMB(Respiration.n.mol.cm2.min ~ Temperature + poly(Mean.Light,2) + Temperature:poly(Mean.Light,2) + (1|Tank), data=y)
anova(m6,m9, test="Chisq")
anova(m6,m10, test="Chisq")

# Photosynethesis - hypothesis tests of factors
# m9 had lowest AICc
m9 <- glmmTMB(Photosynthesis.nmol.cm2.min ~ Temperature + poly(Mean.Light,2) + Species + (1|Tank), data=y)
m10 <- glmmTMB(Photosynthesis.nmol.cm2.min ~ Temperature + poly(Mean.Light,2) + (1|Tank), data=y)
m11 <- glmmTMB(Photosynthesis.nmol.cm2.min ~ Temperature + Species + (1|Tank), data=y)
m12 <- glmmTMB(Photosynthesis.nmol.cm2.min ~ poly(Mean.Light,2) + Species + (1|Tank), data=y)
anova(m9,m10, test="Chisq") # Species
anova(m9,m11, test="Chisq") # Light
anova(m9,m12, test="Chisq") # Temperature
emmeans(m10,pairwise~Temperature)
confint(m10)

# Function to make plot for 2022 data 
Make.plot2022 <- function(dat,preds,ymin,ymax,response,Sp){
  plot(c(220,825), c(ymin,ymax),
       type="n",bty="l",xaxt="n",yaxt="n",xlab="",ylab="")
  preds$cl <- cols[match(preds$Species,cols$Species),2]
  # preds$symbol <- ifelse(preds$Temperature=="cold",19,17)
    
  # preds$Temperature <- factor(preds$Temperature, levels=c("cold","hot"))
  # preds <- preds %>% arrange(Light,Temperature,Species)
  dat$response <- dat[,names(dat)==response]
  dat$cl <- cols[match(dat$Species,cols$Species),2]
  
  with(preds[preds$Temperature=="cold" & preds$Species==Sp,],
       polygon(c(Mean.Light,rev(Mean.Light)),
                 c(lwr,rev(upr)),
                   border=NA,
                   col=adjustcolor(alpha.f=0.6,cl)))
  with(dat[dat$Temperature=="cold" & dat$Species==Sp,],
       points(Mean.Light,response,
              pch=19,
              col=adjustcolor(alpha.f=0.6,cl)))
  with(preds[preds$Temperature=="cold" & preds$Species==Sp,],
       lines(Mean.Light,fit,
             lwd=2,
             col=cl))
  
  with(preds[preds$Temperature=="hot" & preds$Species==Sp,],
       polygon(c(Mean.Light,rev(Mean.Light)),
               c(lwr,rev(upr)),
               border=NA,
               col=adjustcolor(alpha.f=0.4,cl)))
  with(dat[dat$Temperature=="hot" & dat$Species==Sp,],
       points(Mean.Light,response,
              pch=17,
              col=adjustcolor(alpha.f=0.6,cl)))
  with(preds[preds$Temperature=="hot" & preds$Species==Sp,],
       lines(Mean.Light,fit,
             lwd=2,
             lty=2,
             col=cl))
}



# Make plot for 2022 data
quartz(width=5,height=6)
par(mfrow=c(3,2),mar=c(2,2,2,1),oma=c(3,5,1,1))

Make.plot2022(dat=y_mean,preds=dat_2022_growth_output$pred,ymin=0,ymax=0.6,response="Growth.mg.cm2.d",Sp="P. meandrina")
mtext(side=3,"a)",adj=0)
mtext(side=3,"P. meandrina",adj=0.5,line=1.5)
mtext(side=2,"Growth",adj=0.5,line=4.5)
mtext(side=2,expression(paste("(mg ", cm^-2,d^-1,")",sep="")),adj=0.5,line=2.5)
axis(side=2,at=seq(0,3,0.1),las=1,cex.axis=1.2)
axis(side=1,at=seq(200,1000,100))

# aggregate(Mean.Temp ~ Temperature, FUN=mean,data=dat_2022)
legend('bottomright',bty="n",pch=c(19,17),lty=c(1,2),cex=1.2,
       legend=expression(paste(26*degree,"C"), paste(29*degree,"C")),
       col=adjustcolor(alpha.f=0.6,col=c("#0072B2")))

Make.plot2022(dat=y_mean,preds=dat_2022_growth_output$pred,ymin=0,ymax=0.6,response="Growth.mg.cm2.d",Sp="P. tuahiniensis")
mtext(side=3,"b)",adj=0)
mtext(side=3,"P. tuahiniensis",adj=0.5,line=1.5)
axis(side=2,at=seq(0,3,0.1),las=1,cex.axis=1.2)
axis(side=1,at=seq(200,1000,100))

legend('bottomright',bty="n",pch=c(19,17),lty=c(1,2),cex=1.2,
       legend=expression(paste(26*degree,"C"), paste(29*degree,"C")),
       col=adjustcolor(alpha.f=0.6,col=c("#D55E00")))

Make.plot2022(dat=y_mean,preds=dat_2022_respiration_output$pred,ymin=0,ymax=20,response="Respiration.n.mol.cm2.min",Sp="P. meandrina")
mtext(side=3,"c)",adj=0)
mtext(side=2,"Respiration",,adj=0.5,line=4.5)
mtext(side=2,expression(paste("(nmol ", O[2]," ", cm^-2," ", min^-1,")")),adj=0.5,line=2.5)
axis(side=2,at=seq(0,20,2),las=1,cex.axis=1.2)
axis(side=1,at=seq(200,1000,100))

Make.plot2022(dat=y_mean,preds=dat_2022_respiration_output$pred,ymin=0,ymax=20,response="Respiration.n.mol.cm2.min",Sp="P. tuahiniensis")
mtext(side=3,"d)",adj=0)
axis(side=2,at=seq(0,20,2),las=1,cex.axis=1.2)
axis(side=1,at=seq(200,1000,100),cex.axis=1.2)

Make.plot2022(dat=y_mean,preds=dat_2022_photosynthesis_output$pred,ymin=0,ymax=20,response="Photosynthesis.nmol.cm2.min",Sp="P. meandrina")
mtext(side=3,"e)",adj=0)
mtext(side=2,"Photosynthesis",adj=0.5,line=4.5)
mtext(side=2,expression(paste("(nmol ", O[2]," ", cm^-2," ", min^-1,")")),adj=0.5,line=2.5)
axis(side=2,at=seq(0,20,2),las=1,cex.axis=1.2)
axis(side=1,at=seq(200,1000,100),cex.axis=1.2)

Make.plot2022(dat=y_mean,preds=dat_2022_photosynthesis_output$pred,ymin=0,ymax=20,response="Photosynthesis.nmol.cm2.min",Sp="P. tuahiniensis")
mtext(side=3,"f)",adj=0)
axis(side=2,at=seq(0,20,2),las=1,cex.axis=1.2)
axis(side=1,at=seq(200,1000,100),cex.axis=1.2)

mtext(side=1,expression(paste("Light (",mu,"mol photons ", m^-2, s^-1,")")),adj=0.5,outer=T,line=1)

