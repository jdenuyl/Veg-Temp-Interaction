#Submission Code
#Primary Author: James Den Uyl
#Date: 21/9/2023

#Setup####

setwd("~/Documents/RFiles")

dat1<-read.csv("chap1.6.csv") 
dat1$logvalue<-log(dat1$mean.EFFLUX)


#Lmer####

mod1 <- lmer(logvalue ~ Warming*Removal*Elevation + (1|Year) + (1|Site/Plot.E), data=dat1) 
anova(mod1)
summary(mod1)

#interactions
emm <- emmeans(mod1, ~ Warming * Elevation)
emm <- pairs(emm)
emm <- summary(emm)
emm

emm <- emmeans(mod1, ~ Warming * Removal)
emm <- pairs(emm)
emm <- summary(emm)
emm

##

mod1<-lmer(mean.EFFLUX~Total_Cover_Scaled*Warming+ (1|Year) + (1|Site),dat1)
anova(mod1)


#Stepwise model identification####

imputed_data <- mice(dat1, m = 5, maxit = 50, meth = "pmm")
dat <- complete(imputed_data)
dat$logvalue<-log(dat$mean.EFFLUX)

datnew <- dat %>%
  group_by(Plot.ID.E) %>%
  summarise(avg.logvalue = mean(logvalue, na.rm = TRUE)) %>%
  left_join(dat, by = "Plot.ID.E")

datnew <- datnew %>%
  distinct(avg.logvalue, .keep_all = TRUE)

formula<- logvalue ~ mean.temp + Soil_Percent_C + VWC +  
  Inorganic_N + Soil_CN +
  Species_Richness + Total_Cover_Scaled +
  Soil_Percent_C:mean.temp + VWC:mean.temp +  
  Inorganic_N:mean.temp + Soil_CN:mean.temp +
  Species_Richness:mean.temp + Total_Cover_Scaled:mean.temp+
  mean.temp:Total_Cover_Scaled

models <- regsubsets(formula, 
                     data = datnew, nvmax = 5, 
                     method = "forward",
                     force.in =c(1))

res.sum <- summary(models)
data.frame(
  Adj.R2 = which.max(res.sum$adjr2),
  CP = which.min(res.sum$cp),
  BIC = which.min(res.sum$bic)
)

#Functions
get_model_formula <- function(id, object, outcome){
  models <- summary(object)$which[id,-1]
  predictors <- names(which(models == TRUE))
  predictors <- paste(predictors, collapse = "+")
  as.formula(paste0(outcome, "~", predictors))
}
get_cv_error <- function(model.formula, data){
  set.seed(1)
  train.control <- trainControl(method = "cv", number = 5)
  cv <- train(model.formula, data = data, method = "lm",
              trControl = train.control)
  cv$results$RMSE
}

model.ids <- 1:4
cv.errors <-  map(model.ids, get_model_formula, models, "logvalue") %>%
  map(get_cv_error, data = datnew) %>%
  unlist()
cv.errors
which.min(cv.errors)

get_model_formula(3, models, "logvalue")

mod1<-lm(logvalue ~ 
           mean.temp:Total_Cover_Scaled + mean.temp + VWC + Total_Cover_Scaled
         , datnew)
summary(mod1)

plot_model(mod1, sort.est = TRUE,
                           transform = NULL,
                           show.values = TRUE, 
                           value.offset = .3, 
                           type = "std",
                           axis.title="Coefficient",
                           axis.labels=c("Soil Moisture",
                                         "Temperature x Plant Cover", 
                                         "Soil Temperature",
                                         "Plant Cover"),
                           title='',
                           vline.color = "grey80")



#



#Figure 1####


dat1$warmingg <- ifelse(dat1$Warming == 'W',"Warming","Ambient")
my_comparisons <- list(c("Warming","Ambient"))

n <- ggboxplot(dat1, x='warmingg', y = "mean.EFFLUX",
               color = "warmingg", palette = c("#748ffc", '#fa5252'),
               facet.by = "Elevation",
               fill='warmingg',
               alpha=.4,
               outlier.shape = NA)+
  ylab(expression("Soil Respiration (µmol CO"[2]*" m"^-2*" s"^-1*")"))+xlab("")

BE<-n + stat_compare_means(method = "t.test",comparisons = my_comparisons,
                           label = "p.signif",
                           label.y = 2.9)+ ylim(0,3.5)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none",
        legend.title=element_blank())+
  ggtitle('A')



BE

#

dat1$removall <- ifelse(dat1$Removal == 'R',"Removal","Control")

i <- ggboxplot(dat1, x='warmingg', y = "mean.EFFLUX",
               color = "warmingg", palette = c("#748ffc", '#fa5252'),
               facet.by = "removall",
               fill='warmingg',
               alpha=.4,
               outlier.shape = NA)+
  ylab("Soil Respiration (µmol CO2 m-2 s-1)")+xlab("")

GE<-i + stat_compare_means(method = "t.test",comparisons = my_comparisons,
                           label = "p.signif",
                           label.y = 2.9)+ ylim(0,3.5)+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.position = "right",
        axis.ticks.x=element_blank(),
        legend.title=element_blank())+
  ggtitle('B')



GE

#

ga<-grid.arrange(BE, GE,
                 nrow=1,
                 widths=c(2.5,3))


#Figure 2####
my_comparisons <- list(c("High","Low"))

#Boxplot
r <- ggplot(dat1, aes(x=Elevation, y=mean.EFFLUX)) + 
  geom_violin()+ geom_boxplot(width=0.1)
p<-r+xlab("")+
  ylab(expression(paste("Soil respiration (", mu, "mol CO"[2], " m"^{-2}, " s"^{-1}, ")", sep="")))+
  theme(axis.text=element_text(colour="black"),
        axis.line=element_line(color="black"),
        plot.title = element_text(hjust = 0),
        panel.border = element_rect(colour = "black", 
                                    fill=NA, size =0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank())+
  ggtitle("A")+
  scale_x_discrete(labels=rep(c("High","Low")))+
  stat_compare_means(method = "t.test",comparisons = my_comparisons,
                     label = "p.signif")
p

#Coefficient Plot
coefficients <- plot_model(mod1, sort.est = TRUE,
                           transform = NULL,
                           show.values = TRUE, 
                           value.offset = .3, 
                           type = "std",
                           axis.title="Coefficient",
                           axis.labels=c("Soil Moisture",
                                         "Temperature x Plant Cover", 
                                         "Soil Temperature",
                                         "Plant Cover"),
                           title='',
                           vline.color = "grey80")

c <-coefficients+ theme(axis.text=element_text(colour="black"),
                        axis.line=element_line(color="black"),
                        plot.title = element_text(hjust = -.2),
                        panel.grid.major = element_blank(),
                        panel.grid.minor = element_blank(),
                        panel.border = element_blank(),
                        panel.background = element_blank()
)+
  ggtitle("B") 

c

#Combine
grid.arrange(p, c, nrow = 1,
             widths = c(1.3,3)
)




#Figure S3 ####

#Figure Enironmental variation#####

dat2 <- dat1[complete.cases(dat1$VWC), ]

f <- ggboxplot(dat1, x='Elevation', y = "VWC",
               color = "Elevation", palette = c("black", 'grey30'),
               facet.by = "Site", short.panel.labs = FALSE)+
  ylab("Soil Moisture (VWC)")+xlab("")+ coord_flip()+  ggtitle("A") 


labels <- c("Canada", "New Zealand", "Sweden", "United States")
names(labels) <- c("CY", "NZ","SA","US")

VWC<-f + 
  facet_grid(.~Site,labeller = labeller(Site = labels))+ 
  theme(legend.position = "none",
        strip.text.x = element_blank())+
  ylim(0,75)+  ggtitle("D") +
  stat_compare_means(method = "t.test",comparisons = my_comparisons,
                     label = "p.signif",
                     label.y = 64)


VWC

#
dat1$Elevation <- factor(dat1$Elevation, levels=c("Low", "High"))
dat1$Site.F <- factor(dat1$Site, levels=c("CY", "NZ", 'SA','US'))
my_comparisons <- list(c("High","Low"))


f <- ggboxplot(dat1, x='Elevation', y = "mean.EFFLUX",
               color = "Elevation", palette = c("black", 'grey30'),
               facet.by = "Site.F", short.panel.labs = FALSE)+
  ylab(expression(paste("Soil Respiration ("~µmol~CO[2]~m^-2~s^-1*")")))+xlab("")+ coord_flip()

labels <- c("Canada", "New Zealand", "Sweden", "United States")
names(labels) <- c("CY", "NZ","SA","US")

SR<-f + stat_compare_means(method = "t.test",comparisons = my_comparisons,
                           label = "p.signif",
                           label.y = 5.8)+ 
  facet_grid(.~Site.F,labeller = labeller(Site.F = labels))+ 
  theme(legend.position = "none")+
  ggtitle("A") +ylim(0,6.4)


SR

#
dat1$Elevation <- factor(dat1$Elevation, levels=c("Low", "High"))
dat1$Site.F <- factor(dat1$Site, levels=c("CY", "NZ", 'SA','US'))
my_comparisons <- list(c("High","Low"))


f <- ggboxplot(dat1, x='Elevation', y = "Total_Cover_Scaled",
               color = "Elevation", palette = c("black", 'grey30'),
               facet.by = "Site.F", short.panel.labs = FALSE)+
  ylab("Plant Cover (%)")+xlab("")+ coord_flip()

labels <- c("Canada", "New Zealand", "Sweden", "United States")
names(labels) <- c("CY", "NZ","SA","US")

Plant<-f + stat_compare_means(method = "t.test",comparisons = my_comparisons,
                              label = "p.signif",
                              label.y = 150)+ 
  facet_grid(.~Site.F,labeller = labeller(Site.F = labels))+ 
  theme(legend.position = "none",
        strip.text.x = element_blank())+
  ggtitle("C") +ylim(0,190)


Plant

#
dat1$Elevation <- factor(dat1$Elevation, levels=c("Low", "High"))
dat1$Site.F <- factor(dat1$Site, levels=c("CY", "NZ", 'SA','US'))
my_comparisons <- list(c("High","Low"))


f <- ggboxplot(dat1, x='Elevation', y = "mean.temp",
               color = "Elevation", palette = c("black", 'grey30'),
               facet.by = "Site.F", short.panel.labs = FALSE)+
  ylab("Soil Temperature (C)")+xlab("")+ coord_flip()

labels <- c("Canada", "New Zealand", "Sweden", "United States")
names(labels) <- c("CY", "NZ","SA","US")

Temp<-f + stat_compare_means(method = "t.test",comparisons = my_comparisons,
                             label = "p.signif",
                             label.y = 22)+ 
  facet_grid(.~Site.F,labeller = labeller(Site.F = labels))+ 
  theme(legend.position = "none",
        strip.text.x = element_blank())+
  ggtitle("B") +ylim(0,27)


Temp

#Combine
ga<-grid.arrange(SR,Temp, Plant, VWC, nrow = 4)


#

#Figure S4 ####


dat1<-read.csv("chap1.6.csv") # individual data (This is the one in my Chapter 1)
dat1$logvalue<-log(dat1$mean.EFFLUX)
dat1 <- subset(dat1, Total_Plant_Cover_All<200) #Just for imaging

dat1$warmingg <- ifelse(dat1$Warming == 'W',"Warming","Ambient")
my_comparisons <- list(c("Warming","Ambient"))

WE<-ggplot(dat1, aes(x=Total_Plant_Cover_All, y=logvalue, color=warmingg)) +
  geom_point() + 
  geom_smooth(method=lm)+
  scale_color_brewer(palette="Dark2")+
  theme_bw()+
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="right",
        legend.title = element_blank()
        
  )+
  ylab('Soil Respiration (log of µmol CO2 m-2 s-1)')+
  xlab('Plant Cover (%)')

#

dat1$warmingg <- ifelse(dat1$Warming == 'W',"Warming","Ambient")
my_comparisons <- list(c("Warming","Ambient"))

dat1$TPCcat <- ifelse(dat1$Total_Plant_Cover_All > 70,"High Plant Cover", "Low Plant Cover")


n <- ggboxplot(dat1, x='warmingg', y = "mean.EFFLUX",
               color = "warmingg", palette = c("#748ffc", '#fa5252'),
               facet.by = "TPCcat",
               fill='warmingg',
               alpha=.4,
               outlier.shape = NA)+
  ylab("Soil Respiration (µmol CO2 m-2 s-1)")+xlab("")

BE<-n + stat_compare_means(method = "t.test",comparisons = my_comparisons,
                           label = "p.signif",
                           label.y =3.1)+ ylim(0,3.5)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title=element_blank())


BE

#Combine
ga<-grid.arrange(WE, BE,
                 nrow=1,
                 widths=c(2,1))


#



#


#Averages####

dat1high <- subset(dat1, Elevation=='High')
dat1low <- subset(dat1, Elevation=='Low')
dat1highw <- subset(dat1high, Warming=='W')
dat1higha <- subset(dat1high, Warming=='A')
dat1loww <- subset(dat1low, Warming=='W')
dat1lowa <- subset(dat1low, Warming=='A')
dat1r <- subset(dat1, Removal=='R')
dat1r <- subset(dat1r, !is.na(mean.EFFLUX))
dat1c <- subset(dat1, Removal=='C')
dat1c <- subset(dat1c, !is.na(mean.EFFLUX))
dat1w <- subset(dat1, Warming=='W')
dat1w <- subset(dat1w, !is.na(mean.EFFLUX))
dat1a <- subset(dat1, Warming=='A')
dat1a <- subset(dat1a, !is.na(mean.EFFLUX))
US <- subset(dat1, Site=='US')
NZ <- subset(dat1, Site=='NZ')
SA <- subset(dat1, Site=='SA')
CY <- subset(dat1, Site=='CY')
UShigh <- subset(US, Elevation=='High')
USlow <- subset(US, Elevation=='Low')
NZhigh <- subset(NZ, Elevation=='High')
NZlow <- subset(NZ, Elevation=='Low')
SAhigh <- subset(SA, Elevation=='High')
SAlow <- subset(SA, Elevation=='Low')
CYhigh <- subset(CY, Elevation=='High')
CYlow <- subset(CY, Elevation=='Low')
dat1wr <- subset(dat1w, Removal=='R')
dat1wr <- subset(dat1wr, !is.na(mean.EFFLUX))
dat1wc <- subset(dat1w, Removal=='C')
dat1wc <- subset(dat1wc, !is.na(mean.EFFLUX))
dat1ar <- subset(dat1a, Removal=='R')
dat1ar <- subset(dat1ar, !is.na(mean.EFFLUX))
dat1ac <- subset(dat1a, Removal=='C')
dat1ac <- subset(dat1ac, !is.na(mean.EFFLUX))
USr <- subset(US, Removal=='R')
USc <- subset(US, Removal=='C')
NZr <- subset(NZ, Removal=='R')
NZc <- subset(NZ, Removal=='C')
SAr <- subset(SA, Removal=='R')
SAc <- subset(SA, Removal=='C')
CYr <- subset(CY, Removal=='R')
CYc <- subset(CY, Removal=='C')
USw <- subset(US, Warming=='W')
USa <- subset(US, Warming=='A')
NZw <- subset(NZ, Warming=='W')
NZa <- subset(NZ, Warming=='A')
SAw <- subset(SA, Warming=='W')
SAa <- subset(SA, Warming=='A')
CYw <- subset(CY, Warming=='W')
CYa <- subset(CY, Warming=='A')
CYhighw <- subset(CYhigh, Warming=='W')
CYhigha <- subset(CYhigh, Warming=='A')
CYloww <- subset(CYlow, Warming=='W')
CYlowa <- subset(CYlow, Warming=='A')
CYloww <- subset(CYloww, mean.EFFLUX>0)
CYlowa <- subset(CYlowa, mean.EFFLUX>0)
dat1w <- na.omit(dat1w$mean.EFFLUX)
dat1a <- na.omit(dat1a$mean.EFFLUX)


mean(dat1w$mean.EFFLUX)
mean(dat1a$mean.EFFLUX)
mean(dat1r$mean.EFFLUX)
mean(dat1c$mean.EFFLUX)
mean(dat1wr$mean.EFFLUX)
mean(dat1wc$mean.EFFLUX)
mean(dat1ar$mean.EFFLUX)
mean(dat1ac$mean.EFFLUX)


#