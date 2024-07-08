### This code produces Figure 2, and associated analyses, in: 
### Edmunds PJ, Burgess SC. Physiological niches of morphologically 
### cryptic coral species (Pocillopora spp.) in Moorea, French Polynesia
### Code written by Scott Burgess, sburgess@bio.fsu.edu
### Code finalized July 2024

# Load required libraries

library('dplyr')
library('glmmTMB')
library('AICcmodavg')
library('DHARMa')

# sessionInfo()
# R version 4.4.0 (2024-04-24)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sonoma 14.5
# [1] DHARMa_0.4.6     AICcmodavg_2.3-3
# [3] glmmTMB_1.1.9    dplyr_1.1.4   

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


growth_2023 <- read.csv("2023_Growth.csv")
yield_2023 <- read.csv("2023_Yield.csv")

foo <- growth_2023
foo$targetID <- ifelse(grepl("A",growth_2023$Coral.ID),"P. grandis","P. meandrina")
with(foo, table(Species,targetID))


# Function to compare models for 2023 data
Analyze <- function(dat,x,y){
  dat$x <- dat[,names(dat)==x]
  dat$y <- dat[,names(dat)==y,]
  m1 <- glmmTMB(y ~ poly(x,2) + (1|Tank),data=dat)
  m2 <- glmmTMB(y ~ poly(x,1) + (1|Tank), data=dat)
  m3 <- glmmTMB(y ~ 1 + (1|Tank), data=dat)
  aic_table <- aictab(list(m1=m1,m2=m2))
  LLtest1 <- anova(m1,m2,test="Chisq")
  LLtest2 <- anova(m1,m3,test="Chisq")
  LLtest3 <- anova(m2,m3,test="Chisq")
  
  DHARMaOutput <- simulateResiduals(fittedModel = m1, plot = F)
  
  x.vec <- seq(min(dat$x),max(dat$x),length.out=100)
  pred <- expand.grid(x=x.vec,Tank=NA)
  # if(aic_table$Delta_AICc[2]<2) p_use <- get(aic_table$Modnames[aic_table$Modnames==aic_table[which(aic_table$K==min(aic_table$K)),1]])
  # if(aic_table$Delta_AICc[2]>2) p_use <- get(aic_table$Modnames[1])
  p_use <- m1
  p <- predict(p_use,pred,se.fit=T)
  pred$fit <- p$fit
  pred$lwr <- p$fit - 2*p$se.fit
  pred$upr <- p$fit + 2*p$se.fit

  list(pred=pred,
       DHARMaOutput=DHARMaOutput,
       LLtest1=LLtest1,
       LLtest2=LLtest2,
       LLtest3=LLtest3,
       aic_table=aic_table)
}


# growth_2023_temp (P. verrucosa)
growth_2023_temp <- growth_2023 %>% filter(Species %in% "P. verrucosa", Experiment == "Fixed light varied Temp")
growth_2023_temp_others <- growth_2023 %>% filter(Species != "P. verrucosa", Experiment == "Fixed light varied Temp")
with(growth_2023_temp, table(Tank,Temperature))
growth_2023_temp_mean <- aggregate(Growth.mg.cm2.d ~ Tank+Temperature+Species,FUN=mean,data=growth_2023_temp)
growth_2023_temp_output <- Analyze(dat=growth_2023_temp,x="Temperature",y="Growth.mg.cm2.d")
# plot(growth_2023_temp_output$DHARMaOutput)
growth_2023_temp_others <- growth_2023 %>% filter(Species != "P. verrucosa", Experiment == "Fixed light varied Temp")
growth_2023_temp_others_mean <- aggregate(Growth.mg.cm2.d ~ Tank+Temperature+Species,FUN=mean,data=growth_2023_temp_others)
growth_2023_temp_others_output <- Analyze(dat=growth_2023_temp_others,x="Temperature",y="Growth.mg.cm2.d")

# growth_2023_light (P. verrucosa)
growth_2023_light <- growth_2023 %>% filter(Species %in% "P. verrucosa", Experiment == "Fixed Temp Varied light")
with(growth_2023_light, table(Tank,Light,Temperature))
growth_2023_light$Light <- as.numeric(growth_2023_light$Light)
growth_2023_light_mean <- aggregate(Growth.mg.cm2.d ~ Tank+Light+Species,FUN=mean,data=growth_2023_light)
growth_2023_light_mean$Light <- as.numeric(growth_2023_light_mean$Light)
growth_2023_light_output <- Analyze(dat=growth_2023_light,x="Light",y="Growth.mg.cm2.d")
# plot(growth_2023_light_output$DHARMaOutput)
growth_2023_light_others <- growth_2023 %>% filter(Species != "P. verrucosa", Experiment == "Fixed Temp Varied light")
growth_2023_light_others$Light <- as.numeric(growth_2023_light_others$Light)
growth_2023_light_others_mean <- aggregate(Growth.mg.cm2.d ~ Tank+Light+Species,FUN=mean,data=growth_2023_light_others)
growth_2023_light_others_mean$Light <- as.numeric(growth_2023_light_others_mean$Light)
growth_2023_light_others_output <- Analyze(dat=growth_2023_light_others,x="Light",y="Growth.mg.cm2.d")

# yield_2023_temp (P. verrucosa)
yield_2023_temp <- yield_2023 %>% filter(Species %in% "P. verrucosa", Experiment == "Fixed light varied Temp")
yield_2023_temp_mean <- aggregate(Fv.Fm_decimal ~ Tank+Temperature+Species,FUN=mean,data=yield_2023_temp)
yield_2023_temp_output <- Analyze(dat=yield_2023_temp,x="Temperature",y="Fv.Fm_decimal")
# plot(yield_2023_temp_output$DHARMaOutput)
yield_2023_temp_others <- yield_2023 %>% filter(Species != "P. verrucosa", Experiment == "Fixed light varied Temp")
yield_2023_temp_others_mean <- aggregate(Fv.Fm_decimal ~ Tank+Temperature+Species,FUN=mean,data=yield_2023_temp_others)
yield_2023_temp_others_output <- Analyze(dat=yield_2023_temp_others,x="Temperature",y="Fv.Fm_decimal")

# yield_2023_light (P. verrucosa)
yield_2023_light <- yield_2023 %>% filter(Species %in% "P. verrucosa", Experiment == "Fixed Temp Varied light")
yield_2023_light$Light <- as.numeric(yield_2023_light$Light)
yield_2023_light_mean <- aggregate(Fv.Fm_decimal ~ Tank+Light+Species,FUN=mean,data=yield_2023_light)
yield_2023_light_mean$Light <- as.numeric(yield_2023_light_mean$Light)
yield_2023_light_output <- Analyze(dat=yield_2023_light,x="Light",y="Fv.Fm_decimal")
# plot(yield_2023_light_output$DHARMaOutput)
yield_2023_light_others <- yield_2023 %>% filter(Species != "P. verrucosa", Experiment == "Fixed Temp Varied light")
yield_2023_light_others$Light <- as.numeric(yield_2023_light_others$Light)
yield_2023_light_others_mean <- aggregate(Fv.Fm_decimal ~ Tank+Light+Species,FUN=mean,data=yield_2023_light_others)
yield_2023_light_others_mean$Light <- as.numeric(yield_2023_light_others_mean$Light)
yield_2023_light_others_output <- Analyze(dat=yield_2023_light_others,x="Light",y="Fv.Fm_decimal")



# Function to make plot for 2023 data
Make.plot <- function(dat,x,y,preds,xmin,xmax,ymin,ymax){
  plot(c(xmin,xmax), c(ymin,ymax),
       type="n",bty="l",xaxt="n",yaxt="n",xlab="",ylab="")
  cl <- cols[match(dat$Species,cols$Species),2]
  
  with(preds,
       polygon(c(x,rev(x)),
               c(lwr,rev(upr)),
               border=NA,
               col=adjustcolor(alpha.f=0.6,cl)))
  
  x1 <- dat[,names(dat)==x]
  y1 <- dat[,names(dat)==y,]
  points(x1,y1,pch=19,cex=1.2, col=adjustcolor(alpha.f=0.6,cl))
  
  with(preds,lines(x,fit,col=cl))
}

Add.plot <- function(dat,x,y,preds){
  cl <- "grey"
  
  with(preds,
       polygon(c(x,rev(x)),
               c(lwr,rev(upr)),
               border=NA,
               col=adjustcolor(alpha.f=0.6,cl)))
  
  x1 <- dat[,names(dat)==x]
  y1 <- dat[,names(dat)==y,]
  points(x1,y1,pch=19,cex=1.2, col=adjustcolor(alpha.f=0.6,cl))
  
  with(preds,lines(x,fit,col=cl))
}


# Make plot for 2023 data
temp.vec <- seq(25,31,1)
light.vec <- seq(125,1160,200)

quartz(width=5,height=4)
par(mfrow=c(2,2),mar=c(4,4,1,1),oma=c(0,0.5,0.5,1))

Make.plot(dat=growth_2023_temp_mean,
          x="Temperature",
          y="Growth.mg.cm2.d",
          preds=growth_2023_temp_output$pred,
          xmin=25,xmax=31,ymin=0,ymax=2.1)
axis(side=1,at=temp.vec)
axis(side=2,at=seq(0,3,0.5),las=1)
mtext(side=1,expression(paste("Temperature ",degree*C)),line=2.5)
mtext(side=2,expression(paste("Growth (mg ", cm^-2,d^-1,")")),adj=0.8,line=2.5)
mtext(side=3,"a)",adj=0)
pval<-round(growth_2023_temp_output$LLtest2$'Pr(>Chisq)'[2],3)
legend('bottomleft',eval(paste("p =",pval)),bty="n")
# lo <- predict(loess(Growth.mg.cm2.d~Temperature,span=1,data=growth_2023_temp_mean))
# with(growth_2023_temp_mean, lines(Temperature,lo, col='red', lwd=1.5))

# Add.plot(dat=growth_2023_temp_others_mean,
#          x="Temperature",
#          y="Growth.mg.cm2.d",
#          preds=growth_2023_temp_others_output$pred)
# pval<-round(growth_2023_temp_others_output$LLtest$'Pr(>Chisq)'[2],3)
# legend('bottomright',eval(paste("p =",pval)),bty="n",col="grey")


Make.plot(dat=growth_2023_light_mean,
          x="Light",
          y="Growth.mg.cm2.d",
          preds=growth_2023_light_output$pred,
          xmin=120,xmax=1230,ymin=0,ymax=2.1)
axis(side=1,at=light.vec)
axis(side=2,at=seq(0,3,0.5),las=1)
mtext(side=1,expression(paste("Light (",mu,"mol photons ", m^-2, s^-1,")")),line=2.5)
mtext(side=2,expression(paste("Growth (mg ", cm^-2,d^-1,")")),adj=0.8,line=2.5)
mtext(side=3,"b)",adj=0)
pval<-round(growth_2023_light_output$LLtest3$'Pr(>Chisq)'[2],3)
legend('bottomleft',eval(paste("p =",pval)),bty="n")

# Add.plot(dat=growth_2023_light_others_mean,
#          x="Light",
#          y="Growth.mg.cm2.d",
#          preds=growth_2023_light_others_output$pred)
# pval<-round(growth_2023_light_others_output$LLtest$'Pr(>Chisq)'[2],3)
# legend('bottomright',eval(paste("p =",pval)),bty="n",col="grey")

Make.plot(dat=yield_2023_temp_mean,
          x="Temperature",
          y="Fv.Fm_decimal",
          preds=yield_2023_temp_output$pred,
          xmin=25,xmax=31,ymin=0.4,ymax=0.7)
axis(side=1,at=temp.vec)
axis(side=2,at=seq(0,1,0.1),las=1)
mtext(side=1,expression(paste("Temperature ",degree*C)),line=2.5)
mtext(side=2,"Yield (Fv/Fm)",adj=0.8,line=2.5)
mtext(side=3,"c)",adj=0)
pval<-round(yield_2023_temp_output$LLtest2$'Pr(>Chisq)'[2],3)
legend('bottomleft',eval(paste("p =",pval)),bty="n")

# Add.plot(dat=yield_2023_temp_others_mean,
#          x="Temperature",
#          y="Fv.Fm_decimal",
#          preds=yield_2023_temp_others_output$pred)
# pval<-round(yield_2023_temp_others_output$LLtest$'Pr(>Chisq)'[2],3)
# legend('bottomright',eval(paste("p =",pval)),bty="n",col="grey")


Make.plot(dat=yield_2023_light_mean,
          x="Light",
          y="Fv.Fm_decimal",
          preds=yield_2023_light_output$pred,
          xmin=120,xmax=1200,ymin=0.4,ymax=0.7)
axis(side=1,at=light.vec)
axis(side=2,at=seq(0,1,0.1),las=1)
mtext(side=1,expression(paste("Light (",mu,"mol photons ", m^-2, s^-1,")")),line=2.5)
mtext(side=2,"Yield (Fv/Fm)",adj=0.8,line=2.5)
mtext(side=3,"d)",adj=0)
pval<-round(yield_2023_light_output$LLtest3$'Pr(>Chisq)'[2],4)
pval<-format(pval,scientific=F)
legend('bottomleft',eval(paste("p <",pval)),bty="n")

# Add.plot(dat=yield_2023_light_others_mean,
#          x="Light",
#          y="Fv.Fm_decimal",
#          preds=yield_2023_light_others_output$pred)
# pval<-round(yield_2023_light_others_output$LLtest$'Pr(>Chisq)'[2],3)
# legend('bottomright',eval(paste("p =",pval)),bty="n",col="grey")

