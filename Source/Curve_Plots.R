require(ggplot2)
require(Cairo)
require(grid) 
CairoWin()

setwd('C:\\Users\\Joss\\Dropbox\\Final_Project\\Computation\\Data')

pal <- c('#000000', '#444444', '#000099', '#3344AA', '#3377DD', '#00FFBB', '#00FF00', '#999900', '#990088', '#FF0000', '#FF9900', '#FF00BB', '#FFFF22', '#00AA00', '#ff66ff')
pal <- sample(pal)

data <- read.csv('../Data/summaries/methanogens_summary.csv')

best_fits <- data[data$Observed == "False",]
best_fits_growth <- best_fits[which(best_fits$Trait == 'Specific Growth Rate'),]
best_fits_growth <- best_fits[which(best_fits_growth$R_Squared > 0.9),]
best_fits_growth$Trange <- best_fits_growth$Est.Tpk - best_fits_growth$Est.Tmin

unique(best_fits_growth$Species)

View(best_fits_growth)

a <- ggplot(data = best_fits_growth, aes(x=Temperature, y=Response, group = Species, colour=Family)) +
  geom_line() +
  scale_y_sqrt() +
  scale_x_sqrt() + 
  theme_bw() +
  ylab("Observed Growth Rate") +
  theme_bw() +
  xlab("Temperature (K)") +
  ylab(expression(paste("Specific Growth Rate (day"^"-1",")"))) +
  theme(axis.text=element_text(size=24), 
        text = element_text(size=30),
        legend.key.height=unit(3,"line"),
        legend.key.width=unit(3,"line"))
a


ggsave(file="../Results/figs/all_methanogens.png", type="cairo-png", width = 20, height = 40/3)

example <- data[data$Species == "Methanogenium frittonii",]

b <- ggplot(data = example[example$Observed == "False",], aes(x=Temperature, y=Response, group = Species)) +
  geom_line(colour="#4292c6", size=3) +
  geom_point(data = example[example$Observed == "True",], colour="#31a354", size=7, alpha = 1) +
  theme_bw() +
  xlab("Temperature (K)") +
  ylab(expression(paste("Specific Growth Rate (day"^"-1",")"))) +
  theme(axis.text=element_text(size=24), text = element_text(size=30)) +
  guides(fill=FALSE) +
  scale_y_continuous(limits = c(0, 18))
b


ggsave(file="../Results/figs/M_frittoni.png", type="cairo-png", width = 20, height = 40/3)
         