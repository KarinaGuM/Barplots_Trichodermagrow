# karina.gutierrez@cinvestav.mx
# Karina Gutierrez Moreno
# PhD Student on Integrative Biology (Cinvestav Irapuato Mexico)
# Supervisor: Martin Heil
# Statistical analysis of Greenhouse and field experiments (2018-2019)
# June 26th 2020 
# Edited: August 7th 2020


# Library
# Install if necessary
#########################################################
### A) Installing and loading required packages
#########################################################

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}

if (!require("Rmisc")) {
  install.packages("Rmisc", dependencies = TRUE)
  library(Rmisc)
}

if (!require("plyr")) {
  install.packages("plyr", dependencies = TRUE)
  library(plyr)
}

if (!require("dplyr")) {
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}

if (!require("ggpubr")) {
  install.packages("ggpubr", dependencies = TRUE)
  library(ggpubr)
}

# SPORES PRODUCTION ON CULTURE MEDIA

# Set on your working directory
setwd("~/R/Heatmaps_SPME_VOCs_2020")
spores <- read.table("150421_Sporesperml.txt", sep="\t", header=T)
head(spores)

data1 <- spores     # Replicate original data
data1$Hour <- factor(data1$Hour,    # Change ordering manually
                  levels = c("0h", "72h", "168h"))
# https://statisticsglobe.com/ordered-bars-in-ggplot2-barchart-in-r

p <- ggplot(data1, aes(x=Strain, y=Spores, fill=Hour)) +
  geom_bar(position="dodge", stat="identity") + # Quit argument to show leyend
  theme_bw() +
  geom_point() +
  facet_wrap(~Medium) +
  scale_fill_brewer(palette="Blues") +
  labs(y = "Spores/ml", fill= "Hours after inoculation") +
  theme(text = element_text(size=20),
        plot.title = element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=22, color="black"),
        axis.text.y = element_text(size=22, color="black")) +
  theme(strip.text.x = element_text(size = 22, face = "bold"))
p

## How to add standard error bars to the plot?
## VISIT: http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/

# summarySE provides the standard deviation, standard error of the mean, 
# and a (default 95%) confidence interval
tgc <- summarySE(data1, measurevar="Spores", groupvars=c("Medium","Strain","Hour"), na.rm=TRUE)
tgc

## na.rm=TRUE argument to the  end of the argument list should clear it if theres NAs**

## Then --> 
ggplot(tgc, aes(x=Strain, y=Spores, fill=Hour)) +
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  scale_fill_brewer(palette="Blues") +
  ylab("Spores/ml") +
  xlab("Trichoderma strain") +
  geom_point() +
  facet_wrap(~Medium) +
  theme(text = element_text(size=20),
        plot.title = element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=22, color="black"),
        axis.text.y = element_text(size=22, color="black")) +
  theme(strip.text.x = element_text(size = 22, face = "bold")) +
  geom_errorbar(aes(ymin=Spores-se, ymax=Spores+se), # dont forget to change variable here
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))

# STAISTICAL ANALYSIS
cfu.aov <- aov(Spores ~ Medium * Strain * Hour, data=data1)
anova(cfu.aov)

model.tables(cfu.aov, type="means", se=T)

# Tukey post-hoc analysis
PDA <- data1[as.factor(data1$Medium)=="PDA",]
summary(PDA)

cfu.aov1 <- aov(Spores ~ Strain * Hour, data=PDA)
anova(cfu.aov1)

TukeyHSD(cfu.aov1, which="Hour")
TukeyHSD(cfu.aov1, which="Strain:Hour")

# Tukey post-hoc analysis
SEA <- data1[as.factor(data1$Medium)=="SEA",]
summary(SEA)

cfu.aov2 <- aov(Spores ~ Strain * Hour, data=SEA)
anova(cfu.aov2)

TukeyHSD(cfu.aov2, which="Hour")
TukeyHSD(cfu.aov2, which="Strain:Hour")


####################################################################
# RT-qPCR Analysis results
# ng DNA/g of soil
dna <- read.table("150421_ngpergrsoil.txt", sep="\t", header=T)
head(dna)

data2 <- dna     # Replicate original data
data2$Hour <- factor(data2$Hour,    # Change ordering manually
                     levels = c("0h", "72h", "168h"))
# https://statisticsglobe.com/ordered-bars-in-ggplot2-barchart-in-r

p <- ggplot(data2, aes(x=Strain, y=ngDNA, fill=Hour)) +
  geom_bar(position="dodge", stat="identity") + # Quit argument to show leyend
  theme_bw() +
  geom_point() +
  facet_wrap(~Medium) +
  scale_fill_brewer(palette="Blues") +
  labs(y = "ng of DNA/g of soil", fill= "Hours after inoculation") +
  theme(text = element_text(size=20),
        plot.title = element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=22, color="black"),
        axis.text.y = element_text(size=22, color="black")) +
  theme(strip.text.x = element_text(size = 22, face = "bold"))
p

## How to add standard error bars to the plot?
## VISIT: http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/

# summarySE provides the standard deviation, standard error of the mean, 
# and a (default 95%) confidence interval
tgc <- summarySE(data2, measurevar="ngDNA", groupvars=c("Medium","Strain","Hour"), na.rm=TRUE)
tgc

## na.rm=TRUE argument to the  end of the argument list should clear it if theres NAs**

## Then --> 
ggplot(tgc, aes(x=Strain, y=ngDNA, fill=Hour)) +
  geom_bar(position="dodge", stat="identity") +
  theme_bw() +
  scale_fill_brewer(palette="Blues") +
  ylab("ng of DNA/g of soil") +
  xlab("Trichoderma strain") +
  geom_point() +
  facet_wrap(~Medium) +
  theme(text = element_text(size=20),
        plot.title = element_text(size=20),
        axis.title=element_text(size=20,face="bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size=22, color="black"),
        axis.text.y = element_text(size=22, color="black")) +
  theme(strip.text.x = element_text(size = 22, face = "bold")) +
  geom_errorbar(aes(ymin=ngDNA-se, ymax=ngDNA+se), # dont forget to change variable here
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))

# STATISTICAL ANALYSIS
cfu.aov <- aov(ngDNA ~ Medium * Strain * Hour, data=data2)
anova(cfu.aov)

model.tables(cfu.aov, type="means", se=T)

# Tukey post-hoc analysis
Soil <- data2[as.factor(data2$Medium)=="Natural_soil",]
summary(Soil)

cfu.aov1 <- aov(ngDNA ~ Strain * Hour, data=Soil)
anova(cfu.aov1)

TukeyHSD(cfu.aov1, which="Hour")
TukeyHSD(cfu.aov1, which="Strain:Hour")

# Tukey post-hoc analysis
SEA <- data1[as.factor(data1$Medium)=="SEA",]
summary(SEA)

cfu.aov2 <- aov(Spores ~ Strain * Hour, data=SEA)
anova(cfu.aov2)

TukeyHSD(cfu.aov2, which="Hour")
TukeyHSD(cfu.aov2, which="Strain:Hour")



