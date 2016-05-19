require(ggplot2)
require(reshape2)
require(colorspace)
require(Cairo)
require(scales)
require(grid) 
CairoWin()

rm(list=ls())
setwd('C:\\Users\\Joss\\Dropbox\\Final_Project\\Computation\\Data')

#simple high contrast palette
pal <- c('#000000', '#444444', '#000099', '#3344AA', '#3377DD', '#00FFBB', '#00FF00', '#999900', '#990088', '#FF0000', '#FF9900', '#FF00BB', '#FFFF22', '#00AA00', '#ff66ff', '#ecc6ec')
pal <- sample(pal)

main_theme <- theme(axis.text.x = element_text(size = 14),
                    axis.text.y = element_text(size = 14),
                    axis.title.y = element_text(size = 20),
                    axis.title.x = element_text(size = 18),
                    plot.title = element_text(size=22, vjust=1),
                    legend.text=element_text(size=18),
                    legend.position="bottom",
                    legend.margin=unit(0,"cm"),
                    plot.margin=unit(c(0,0,0,0), "cm"))

gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}

mround <- function(x,base){ 
  base*round(x/base) 
} 

data <- read.csv('../Data/summaries/Summary.csv')

model_scores <- ggplot(data, aes(x=Model_name, fill=factor(Rank))) +
  geom_bar(position="Dodge")
model_scores

model_scores <- ggplot(data[data$Rank ==1,], aes(x=Number.of.Data.Points, fill=factor(Model_name))) +
  geom_histogram(position = "fill", binwidth=10) +
  stat_bin(aes(label=..count.., y = ((..count.. > 0) * (cumsum(..count..))/sum(..count..))), vjust=2, geom="text", position="identity", binwidth=10) +
  xlim(c(0,100)) +
  ylab("Frequency") +
  xlab("Number of Observations")
model_scores

model_analysis <- ggplot(data, aes(x=log10(Number.of.Data.Points), y = AIC, colour=Model_name)) +
  geom_jitter(position = position_jitter(width = 0.5, height = 0.0), alpha=0.5) +
  geom_smooth(method="loess", alpha = 0.1) +
  theme_bw()
model_analysis

View(data)

best_fits <- data[data$Rank==1,]

best_fits_mean <- aggregate(cbind(Max.response, Est.Tpk, Est.Tmin, Est.Tmax, E, B0, E_D)~Species:Trait*Kingdom*Phylum*Class*Order*Family*Genus*Acceptor, data=best_fits,  FUN=gm_mean)

best_fits_growth <- best_fits_mean[best_fits_mean$Trait == 'Specific Growth Rate',]
best_fits_growth$Trange <- best_fits_growth$Est.Tpk - best_fits_growth$Est.Tmin

require(plyr)
count(data.frame(best_fits_growth$Phylum))

best_fits_bacteria <-best_fits_growth[best_fits_growth$Kingdom == 'Bacteria',]
best_fits_archaea <-best_fits_growth[best_fits_growth$Kingdom == 'Archaea',]

best_fits_Proteobacteria <-best_fits_growth[best_fits_growth$Phylum == 'Proteobacteria',]
best_fits_Tenericutes <-best_fits_growth[best_fits_growth$Phylum == 'Tenericutes',]
best_fits_Firmicutes <-best_fits_growth[best_fits_growth$Phylum == 'Firmicutes',]
best_fits_Cyanophyta <-best_fits_growth[best_fits_growth$Phylum == 'Cyanophyta',]
best_fits_Euryarchaeota <-best_fits_growth[best_fits_growth$Phylum == 'Euryarchaeota',]
best_fits_Crenarchaeota <-best_fits_growth[best_fits_growth$Phylum == 'Crenarchaeota',]

best_fits_growth[best_fits_growth$Acceptor==""] <- NA
best_fits_growth[best_fits_growth$Acceptor=="Na",] <- NA
acceptors <- best_fits_growth[complete.cases(best_fits_growth$Acceptor),]

hist <- ggplot(best_fits_growth, aes(x=Est.Tpk)) +
  geom_histogram(binwidth=5)
hist

a <- ggplot(best_fits_growth, aes(x=Est.Tpk, y=E_D, colour=Acceptor, label=Species)) +
  geom_point() +
  scale_y_continuous(limits = c(0, 100))
a

a1 <- ggplot(best_fits_growth, aes(x=Est.Tpk, y=E, colour=Acceptor)) +
  geom_smooth(method='loess', alpha = 0) +
  scale_y_continuous(limits = c(0, 6)) +
  geom_point() +
  theme_bw()
a1

f <- ggplot(best_fits_growth, aes(x=Est.Tpk, y=B0, colour=Kingdom)) +
  geom_smooth(method='lm') +
  geom_point() +
  scale_y_continuous(limits = c(-70, 50))
f

c <- ggplot(best_fits_growth, aes(x=Est.Tpk, y=Trange, colour=Acceptor)) +
  geom_smooth(method='lm', alpha=0) +
  geom_point()
c

mean_cheat <- best_fits_growth
mean_cheat$Est.Tpk <- mround(mean_cheat$Est.Tpk, 5)
  
d <- ggplot(best_fits_growth, aes(x=Est.Tpk, y=Max.response)) +  
  geom_point(alpha=0.5, size=4, colour='chartreuse3') +
  theme_bw() +
  coord_trans(y = "log10") +    
  xlab('Optimum Temperature (K)') +
  ylab(expression(paste("Specific Growth Rate (Day"^"-1", ")"))) +
  facet_grid(Kingdom~.) +
  scale_y_continuous(breaks = c(1,10,100)) +
  theme(text=element_text(size=16), strip.text.y = element_text(size = 16), text = element_text(size = 16)) +
  stat_summary(data=mean_cheat, aes(x=Est.Tpk, y=Max.response, linewith=2), fun.y = "mean", colour = "red", geom = "line")
d

d1 <- ggplot(best_fits_growth, aes(x=Est.Tpk-273.15, y=(Max.response/(24*60)))) +  
  geom_point(alpha=0.7, size=2, colour='chartreuse3') +
  theme_bw() +
  xlab('Optimum Temperature (C)') +
  ylab(expression(paste("Specific Growth Rate (Day"^"-1", ")"))) +
  theme(text=element_text(size=16), strip.text.y = element_text(size = 16), text = element_text(size = 16))
d1

ggsave(file="../Results/figs/Match.png", type="cairo-png", width = 10, height = 10)

e <- ggplot(best_fits_bacteria, aes(x=Est.Tpk, y=Max.response, colour=Phylum)) +
  #geom_smooth(method='lm', alpha=0) +
  geom_point() +
  theme_bw() +
  scale_y_sqrt() +
  scale_x_sqrt() + 
  scale_colour_manual(values=pal)
e 

f <- ggplot(best_fits_archaea, aes(x=Est.Tpk, y=Max.response, colour=Phylum)) +
  scale_y_log10() + 
  geom_smooth(method='lm') + 
  geom_point() +
  theme_bw() +
  scale_colour_manual(values=pal)
f 

unique(data$Acceptor)
unique(data$Donor)

accept_phyla <- ggplot(acceptors, aes(x=Est.Tpk, y=Max.response)) +
  scale_y_continuous(limits = c(0, 80)) +
  #geom_smooth(method='glm', alpha=0.2) +
  geom_point() +
  facet_grid(Acceptor~Phylum) +
  #scale_colour_manual(values=pal) +
  theme_bw()
accept_phyla

all_phyla <- ggplot(best_fits_growth, aes(x=Est.Tpk-273.15, y=Max.response)) +
  #scale_y_log10() +
  #scale_x_log10() +
  #geom_smooth(method='glm', alpha=0.2) +
  geom_point() +
  facet_wrap( ~ Phylum, nrow = 3) +
  xlab("Temperature (C)") +
  ylab("rMax") +
  #scale_y_continuous(limits = c(0, 35)) +
  #scale_colour_manual(values=pal) +
  theme_bw()
all_phyla

accept <- ggplot(acceptors[acceptors$Acceptor  %in% c('Oxygen'),], aes(x=Est.Tpk, y=Max.response, colour=Phylum)) +
  scale_y_continuous(limits = c(0, 80)) +
  geom_smooth(method='glm', alpha=0.2) +
  geom_point() +
  #scale_colour_manual(values=pal) +
  theme_bw()
accept

accept_phyla <- ggplot(acceptors[acceptors$Acceptor  %in% c('Oxygen'),], aes(x=Est.Tpk, y=Max.response, colour=Phylum)) +
  scale_y_continuous(limits = c(0, 80)) +
  geom_smooth(method='glm', alpha=0.2) +
  geom_point() +
  geom_text(aes(label=Species, vjust=-2)) +
  #scale_colour_manual(values=pal) +
  theme_bw()
accept_phyla

acceptB0 <- ggplot(acceptors, aes(x=Est.Tpk, y=B0, colour=Acceptor)) +
  geom_point() +
  scale_y_continuous(limits = c(-70, 50)) +
  scale_colour_manual(values=pal)
acceptB0

proteoplus <- table(best_fits_Proteobacteria$Class)
proteo_plt_data <- best_fits_Proteobacteria[best_fits_Proteobacteria$Class %in% names(proteoplus[proteoplus > 5]), ]

firmiplus <- table(best_fits_Firmicutes$Order)
firmi_plt_data <- best_fits_Firmicutes[best_fits_Firmicutes$Order %in% names(firmiplus[firmiplus >= 5]), ]

Euryplus <- table(best_fits_Euryarchaeota$Order)
Eur_plt_data <- best_fits_Euryarchaeota[best_fits_Euryarchaeota$Order %in% names(Euryplus[Euryplus >= 5]), ]

Cyanoplus <- table(best_fits_Cyanophyta$Order)
Cyano_plt_data <- best_fits_Cyanophyta[best_fits_Cyanophyta$Order %in% names(Cyanoplus[Cyanoplus > 5]), ]

Crenarchplus <- table(best_fits_Crenarchaeota$Order)
Crenarch_plt_data <- best_fits_Crenarchaeota[best_fits_Crenarchaeota$Order %in% names(Crenarchplus[Crenarchplus > 5]), ]

firmi_lm <- lm(Max.response ~ Est.Tpk:Class, data=best_fits_Firmicutes)
summary(firmi_lm)

Eur_lm <- lm(Max.response ~ Est.Tpk:Class, data=Eur_plt_data)
summary(Eur_lm)

proteo <- ggplot(proteo_plt_data, aes(x=Est.Tpk, y=Max.response, colour = Class)) +
  scale_y_sqrt() +
  scale_x_sqrt() + 
  geom_smooth(method='lm', alpha=0.15) +
  geom_point(aes(size=2)) +
  guides(size=FALSE,colour=guide_legend(nrow=2,byrow=TRUE)) +
  xlab('Optimal Temperature (K)') +
  ylab(expression(paste("Specific Growth Rate (Day"^"-1", ")"))) +
  ggtitle('Proteobacteria') +
  theme_bw() +
  main_theme

firmi <- ggplot(firmi_plt_data, aes(x=Est.Tpk, y=Max.response, colour = Order, label = Species)) +
  scale_x_sqrt() + 
  scale_y_sqrt(lim=c(0, 85)) + 
  geom_smooth(method='lm', alpha = 0.25) +
  geom_point(aes(size=2)) +
  guides(size=FALSE,colour=guide_legend(nrow=2,byrow=TRUE)) +
  xlab('Optimal Temperature (K)') +
  ylab(expression(paste("Specific Growth Rate (Day"^"-1", ")"))) +
  ggtitle('Firmicutes') +
  theme_bw() +
  main_theme  +
  theme(axis.title.y=element_blank(),
        legend.title=element_blank())

Eur_plt <- ggplot(Eur_plt_data, aes(x=Est.Tpk, y=Max.response, colour = Order)) +
  scale_x_sqrt() + 
  scale_y_sqrt(lim=c(0, 40)) + 
  geom_smooth(method='lm', alpha=0.25) +
  geom_point(aes(size=2)) +
  guides(size=FALSE, colour=guide_legend(nrow=2)) +
  xlab('Optimal Temperature (K)') +
  ylab(expression(paste("Specific Growth Rate (Day"^"-1", ")"))) +
  ggtitle('Euryarchaeota') +
  theme_bw() +
  main_theme +
  theme(legend.title=element_blank())

Cyano <- ggplot(Cyano_plt_data, aes(x=Est.Tpk, y=Max.response, colour = Order, label = Species)) +
  scale_x_sqrt() +
  scale_y_sqrt() + 
  geom_smooth(method='glm', alpha=0.15) +
  geom_point(aes(size=2)) +
  guides(size=FALSE) +
  xlab('Optimal Temperature (K)') +
  ylab(expression(paste("Specific Growth Rate (Day"^"-1", ")"))) +
  ggtitle('Cyanophyta') +
  theme_bw() +
  main_theme +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        legend.title=element_blank())


Crenarch <- ggplot(Crenarch_plt_data, aes(x=Est.Tpk, y=Max.response, colour = Order, label = Species)) +
  scale_x_sqrt() + 
  scale_y_sqrt() + 
  geom_smooth(method='glm', alpha=0.15) +
  geom_point(aes(size=2)) +
  guides(size=FALSE) +
  xlab('Optimal Temperature (K)') +
  ylab(expression(paste("Specific Growth Rate (Day"^"-1", ")"))) +
  ggtitle('Crenarchaeota') +
  theme_bw() +
  main_theme +
  theme(axis.title.x=element_blank(),
        legend.title=element_blank())


proteo
firmi
Eur_plt
Cyano
Crenarch

require(gridExtra)
all_orders <- arrangeGrob(Crenarch,Cyano,Eur_plt,firmi)
all_orders
ggsave(file="../Results/figs/all_orders.png", all_orders, type="cairo-png", width = 16, height = 10)

test2 <- glm(Max.response ~ Est.Tpk:Acceptor, data=acceptors)
summary(test2)

test <- glm(Max.response ~ Est.Tpk:Phylum, data=best_fits_growth)
summary(test)
View(data)

mp <- NULL
mapWorld <- borders("world", colour="gray", fill="gray50") # create a layer of borders
mp <- ggplot(data) + mapWorld
mp <- mp + 
      geom_point(aes(x=Longitude, y=Latitude, colour=Phylum) , size=2) +
      geom_jitter(aes(x=3, y=3))
mp

twoplus <- table(best_fits_growth$Phylum)
boxplotdata <- best_fits_growth[best_fits_growth$Phylum %in% names(twoplus[twoplus >= 10]), ]
boxplotdata$Phylum <- gsub("Deinococcus-Thermus", "Deinococcus\nThermus   ", boxplotdata$Phylum)

bp1 <- ggplot(boxplotdata[boxplotdata$Kingdom %in% c('Bacteria', 'Archaea'),], aes(x=Phylum, y=Max.response)) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(aes(fill=Kingdom)) +
  scale_x_discrete(limits=c('Crenarchaeota', 'Euryarchaeota', 'Cyanophyta', 'Deinococcus\nThermus   ', 'Proteobacteria', 'Tenericutes')) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12)) +
  coord_trans(y = "log10") +
  scale_y_continuous(breaks = c(1,10,100)) +
  scale_fill_manual(values=c('cyan3','chartreuse2')) +
  ylab(expression(paste("Specific Growth Rate (Day"^"-1", ")")))
bp1

ggsave(file="../Results/figs/Bacteria_Phyla_Opts.png", type="cairo-png", width = 10, height = 4)


threeplus <- table(best_fits_growth$Acceptor)
boxplotdata2 <- best_fits_growth[best_fits_growth$Acceptor %in% names(threeplus[threeplus >= 10]), ]

bp2 <- ggplot(boxplotdata2[boxplotdata2$Kingdom %in% c('Bacteria', 'Archaea'),], aes(x=Acceptor, y=Max.response)) +
  geom_boxplot() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        legend.text=element_text(size=12),
        legend.title=element_text(size=12)) +
  coord_trans(y = "log10") +
  scale_y_continuous(breaks = c(1,10,100)) +
  ylab(expression(paste("Specific Growth Rate (Day"^"-1", ")")))
bp2
