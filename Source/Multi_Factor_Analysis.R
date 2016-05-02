#install.packages("ggplot2")
#install.packages("reshape2")
#install.packages("colorspace")
#install.packages("scales")
#install.packages("plyr")
#install.packages("grid")
#install.packages("gridExtra")

#------Package requirements------

require(ggplot2)
require(reshape2)
require(plyr)
require(Cairo)
require(scales)
require(grid) 
require(gridExtra)

CairoWin() #On windows cairo plotting has antialiasing while native plotting does not

#------Plotting Utilities------

#Simple black and white theme that works nicely with gridExtra
main_theme <- theme_bw() + 
              theme(axis.text.x = element_text(size = 14),
                    axis.text.y = element_text(size = 14),
                    axis.title.y = element_text(size = 16),
                    axis.title.x = element_text(size = 16),
                    plot.title = element_text(size=18, vjust=1),
                    legend.text=element_text(size=14),
                    legend.position="bottom",
                    legend.margin=unit(0,"cm"),
                    plot.margin=unit(c(0.1,0.1,0.1,0.1), "cm"))

#simple large high contrast palette + order randomisation 
HC_pal <- c('#000000', '#444444', '#000099', '#3344AA', '#3377DD', '#00FFBB', '#00FF00', '#999900', '#990088', '#FF0000', '#FF9900', '#FF00BB', '#FFFF22', '#00AA00', '#ff66ff', '#ecc6ec')
HC_pal <- sample(pal)

#------Declare Global Functions------

# Geometric mean that handles NAs in the dataset well

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

#Round a column to arbitary number (useful for binning temperature values)

mround <- function(x,base){ 
  base*round(x/base) 
} 

# Extract the legend from a plot

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

#------Read and Setup Data------
setwd('C:\\Users\\Joss\\Dropbox\\Final_Project\\Computation\\Data')

data <- read.csv('../Data/summaries/summary_TS_edit.csv')

#Replace blanks in respiration type column with NA
data$Respiration.Type[data$Respiration.Type == ''] <- NA

#Summary contains all fitted models so select only the best ones
best_fits <- data[data$Rank==1,]

#Get growth curves only 
best_fits_growth <- best_fits[best_fits$Trait == 'Specific Growth Rate',]

#Geometric means for strains with more than one entry
best_fits_mean_growth <- aggregate(cbind(Max.response, Est.Tpk, Est.Tmin, Est.Tmax)~Species:Trait*Kingdom*Phylum*Class*Order*Family*Genus*Acceptor, data=best_fits_growth,  FUN=gm_mean)

#Take logs
log_data <- best_fits_mean_growth
log_data$Est.Tpk <- log(log_data$Est.Tpk)
log_data$Max.response <- log(log_data$Max.response)

#Separate by kingdom
best_fits_bacteria <-best_fits_mean[best_fits_mean_growth$Kingdom == 'Bacteria',]
best_fits_archaea <-best_fits_mean[best_fits_mean_growth$Kingdom == 'Archaea',]

best_fits_bacteria_log <-log_data[log_data$Kingdom == 'Bacteria',]
best_fits_archaea_log <-log_data[log_data$Kingdom == 'Archaea',]

#------Quick Analysis of Which models are best based on the number of data points in the curve------

model_scores <- ggplot(best_fits, aes(x=Number.of.Data.Points, fill=factor(Model_name))) +
  geom_histogram(position = "fill", binwidth=10) +
  geom_text(position = "stack") +
  xlim(c(0,100)) +
  ylab("Frequency") +
  xlab("Number of Observations")
model_scores

#------Create Plots and Log Plots of All Topt and TPK values------

#Just run everything from here to the arrange histograms section, none of it is called immediately

#All species

all_reponses <- ggplot(best_fits_growth, aes(x=Est.Tpk, y=Max.response, colour=Kingdom)) +  
  geom_point(alpha=0.5, size=4) +
  main_theme +
  xlab('Optimum Temperature (K)') +
  ylab(expression(paste("Specific Growth Rate (Day"^"-1", ")")))  +
  theme(legend.position = "none") 
#all_reponses

all_reponses_log <- ggplot(log_data, aes(x=Est.Tpk, y=Max.response, colour=Kingdom)) +  
  geom_point(alpha=0.5, size=4) +
  main_theme +
  xlab('Log(Optimum Temperature (K))') +
  ylab(expression(paste("Log(Specific Growth Rate (Day"^"-1", "))")))  +
  theme(legend.position = "none") 
#all_reponses_log 

#Bacteria Only

bacteria_reponses <- ggplot(best_fits_bacteria, aes(x=Est.Tpk, y=Max.response)) +  
  geom_point(alpha=0.5, size=4, colour='#00BFC4') +
  main_theme +
  xlab('Optimum Temperature (K)') +
  ylab(expression(paste("Specific Growth Rate (Day"^"-1", ")")))  +
  theme(legend.position = "none") 
#bacteria_reponses

bacteria_reponses_log <- ggplot(best_fits_bacteria_log, aes(x=Est.Tpk, y=Max.response)) +  
  geom_point(alpha=0.5, size=4, colour='#00BFC4') +
  main_theme +
  xlab('Log(Optimum Temperature (K))') +
  ylab(expression(paste("Log(Specific Growth Rate (Day"^"-1", "))")))  +
  theme(legend.position = "none") 
#bacteria_reponses_log

#Archaea only

archaea_reponses <- ggplot(best_fits_archaea, aes(x=Est.Tpk, y=Max.response)) +  
  geom_point(alpha=0.5, size=4, colour='#F8766D') +
  main_theme +
  xlab('Optimum Temperature (K)') +
  ylab(expression(paste("Specific Growth Rate (Day"^"-1", ")"))) 
#archaea_reponses

archaea_reponses_log <- ggplot(best_fits_archaea_log, aes(x=Est.Tpk, y=Max.response)) +  
  geom_point(alpha=0.5, size=4, colour='#F8766D') +
  main_theme +
  xlab('Log(Optimum Temperature (K))') +
  ylab(expression(paste("Log(Specific Growth Rate (Day"^"-1", "))"))) +
  theme(legend.position = "none") 
#archaea_reponses_log


#------Create Histograms of Topt and Tpk values------
#Histograms of max response frequences and tpk frequencies

hist_all_TPK <- ggplot(best_fits_growth, aes(x=Est.Tpk, fill=Kingdom)) +
  main_theme +
  xlab(' ') +
  ylab(expression(paste("Count"^"1"))) +
  ylim(c(0, 100)) +
  theme(legend.position = "none") +
  geom_histogram()
#hist_all_TPK

hist_all_Resp <- ggplot(best_fits_growth, aes(x=Max.response, fill=Kingdom)) +
  coord_flip() +
  main_theme +
  xlab(' ') +
  ylab('Count') + 
  theme(legend.position = "none") +
  geom_histogram()
#hist_all_Resp

#logs of same thing

hist_all_TPK_log <- ggplot(log_data, aes(x=Est.Tpk, fill=Kingdom)) +
  main_theme +
  xlab(' ') +
  ylab(expression(paste("Count"^"1"))) +
  theme(legend.position = "none") +
  geom_histogram()
#hist_all_TPK_log

hist_all_Resp_log <- ggplot(log_data, aes(x=Max.response, fill=Kingdom)) +
  coord_flip() +
  main_theme +
  xlab(' ') +
  ylab('Count') + 
  theme(legend.position = "none") +
  geom_histogram()
#hist_all_Resp_log

#Again for Bacteria

hist_bacteria_TPK <- ggplot(best_fits_bacteria, aes(x=Est.Tpk)) +
  main_theme +
  xlab(' ') +
  ylab(expression(paste("Count"^"1"))) +
  ylim(c(0, 100)) +
  theme(legend.position = "none") +
  geom_histogram(fill='#00BFC4')
#hist_bacteria_TPK

hist_bacteria_Resp <- ggplot(best_fits_bacteria, aes(x=Max.response)) +
  coord_flip() +
  main_theme +
  xlab(' ') +
  ylab('Count') + 
  theme(legend.position = "none") +
  geom_histogram(fill='#00BFC4')
#hist_bacteria_Resp

hist_bacteria_TPK_log <- ggplot(best_fits_bacteria_log, aes(x=Est.Tpk)) +
  main_theme +
  xlab(' ') +
  ylab(expression(paste("Count"^"1"))) +
  theme(legend.position = "none") +
  geom_histogram(fill='#00BFC4')
#hist_bacteria_TPK_log

hist_bacteria_Resp_log <- ggplot(best_fits_bacteria_log, aes(x=Max.response)) +
  coord_flip() +
  main_theme +
  xlab(' ') +
  ylab('Count') + 
  theme(legend.position = "none") +
  geom_histogram(fill='#00BFC4')
#hist_bacteria_Resp_log

#And the same again for archaea

hist_archaea_TPK <- ggplot(best_fits_archaea, aes(x=Est.Tpk)) +
  main_theme +
  xlab(' ') +
  ylab(expression(paste("Count"^"1"))) +
  scale_y_continuous(breaks=c(0,5,10)) +
  theme(legend.position = "none") +
  geom_histogram(fill='#F8766D')
#hist_archaea_TPK

hist_archaea_Resp <- ggplot(best_fits_archaea, aes(x=Max.response)) +
  coord_flip() +
  main_theme +
  xlab(' ') +
  ylab('Count') + 
  theme(legend.position = "none") +
  geom_histogram(fill='#F8766D')
#hist_archaea_Resp

hist_archaea_TPK_log <- ggplot(best_fits_archaea_log, aes(x=Est.Tpk)) +
  main_theme +
  xlab(' ') +
  ylab(expression(paste("Count"^"1"))) +
  theme(legend.position = "none") +
  geom_histogram(fill='#F8766D')
#hist_archaea_TPK_log

hist_archaea_Resp_log <- ggplot(best_fits_archaea_log, aes(x=Max.response)) +
  coord_flip() +
  main_theme +
  xlab(' ') +
  ylab('Count') + 
  theme(legend.position = "none") +
  geom_histogram(fill='#F8766D')
#hist_archaea_Resp_log

#------Setup Grid Utilities------

#Create a plot which has a nicely formatted legend
legend_plot <- ggplot(best_fits_growth, aes(x=Est.Tpk, y=Max.response, colour=Kingdom)) +
  main_theme +
  geom_point(size=4)

#Extract said legend
legend <- g_legend(legend_plot)

#------Arrange the histograms, scatters and legends into a single plot------
#Some axes don't quite line up - they will if the limits are set manually

#All species
grid.arrange(hist_all_TPK, legend, all_reponses, hist_all_Resp, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

#Log all species
grid.arrange(hist_all_TPK_log, legend, all_reponses_log, hist_all_Resp_log, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

#Bacteria only
grid.arrange(hist_bacteria_TPK, legend, bacteria_reponses, hist_bacteria_Resp, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

#Log bacteria Only
grid.arrange(hist_bacteria_TPK_log, legend, bacteria_reponses_log, hist_bacteria_Resp_log, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

#Archaea only
grid.arrange(hist_archaea_TPK, legend, archaea_reponses, hist_archaea_Resp, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

#Log Archaea only
grid.arrange(hist_archaea_TPK_log, legend, archaea_reponses_log, hist_archaea_Resp_log, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))

#------Multi Variable Analysis------

overall_lm <- lm(Max.response ~ Est.Tpk*Kingdom, data=best_fits_mean_growth)
summary(overall_lm)

bacteria_lm <- lm(Max.response ~ Est.Tpk, data=best_fits_bacteria)
summary(bacteria_lm)

archaea_lm <- lm(Max.response ~ Est.Tpk, data=best_fits_archaea)
summary(archaea_lm)
