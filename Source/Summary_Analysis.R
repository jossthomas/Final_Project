require(ggplot2)
require(reshape2)
require(colorspace)
require(Cairo)
require(grid) 
CairoWin()

setwd('C:\\Users\\Joss\\Dropbox\\Final_Project\\Computation\\Data')

data <- read.csv('../Data/summaries/Summary.csv')

#simple high contrast palette
pal <- c('#000000', '#444444', '#000099', '#3344AA', '#3377DD', '#00FFBB', '#00FF00', '#999900', '#990088', '#FF0000', '#FF9900', '#FF00BB', '#FFFF22', '#00AA00', '#ff66ff')
pal <- sample(pal)

View(data)

best_fits <- data[data$Plotted == "True",]
best_fits <- best_fits[best_fits$R_Squared > 0.75,]
best_fits_growth <- best_fits[which(best_fits$Trait == 'Specific Growth Rate'),]
best_fits_growth$Trange <- best_fits_growth$Est.Tpk - best_fits_growth$Est.Tmin

best_fits_bacteria <-best_fits_growth[best_fits_growth$Kingdom == 'Bacteria',]
best_fits_archaea <-best_fits_growth[best_fits_growth$Kingdom == 'Archaea',]

best_fits_Proteobacteria <-best_fits_growth[best_fits_growth$Phylum == 'Proteobacteria',]
best_fits_Tenericutes <-best_fits_growth[best_fits_growth$Phylum == 'Tenericutes',]
best_fits_Firmicutes <-best_fits_growth[best_fits_growth$Phylum == 'Firmicutes',]

hist <- ggplot(best_fits, aes(x=Est.Tpk)) +
  geom_histogram()
hist

a <- ggplot(best_fits_growth, aes(x=Est.Tpk, y=E_D, colour=Kingdom, label=Species)) +
  geom_point() +
  geom_text()
a

a1 <- ggplot(best_fits_growth, aes(x=Est.Tpk, y=E, colour=Kingdom)) +
  geom_smooth(method='lm') +
  scale_y_continuous(limits = c(0, 6)) +
  geom_point()
a1

f <- ggplot(best_fits_growth, aes(x=Est.Tpk, y=B0, colour=Kingdom)) +
  geom_smooth(method='lm') +
  geom_point() +
  scale_y_continuous(limits = c(-70, 50))
f


c <- ggplot(best_fits_growth, aes(x=Est.Tpk, y=Trange, colour=Kingdom)) +
  geom_smooth(method='lm') +
  geom_point()
c


d <- ggplot(best_fits_growth, aes(x=Est.Tpk, y=Max.response, colour=Kingdom)) +  
  geom_point() +
  geom_smooth(method='lm') +
  theme_bw() +
  xlab('Optimum Temperature') +
  ylab('Specific Growth Rate')
d

e <- ggplot(best_fits_bacteria, aes(x=Est.Tpk, y=Max.response, colour=Phylum)) +
  #geom_smooth(method='lm', alpha=0) +
  geom_point() +
  theme_bw() +
  scale_colour_manual(values=pal)
e 

f <- ggplot(best_fits_archaea, aes(x=Est.Tpk, y=Max.response, colour=Phylum)) +
  #geom_smooth(method='lm') +
  scale_y_continuous(limits = c(0, 100)) +  
  geom_point() +
  theme_bw() +
  scale_colour_manual(values=pal)
f 

best_fits_growth[best_fits_growth==""] <- NA
acceptors <- best_fits_growth[complete.cases(best_fits_growth$Acceptor),]

unique(data$Acceptor)
unique(data$Donor)

accept <- ggplot(acceptors[acceptors$Acceptor  %in% c('Sulfate'),], aes(x=Est.Tpk, y=Max.response)) +
  scale_y_continuous(limits = c(0, 100)) +
  geom_smooth(method='glm', alpha=0.2) +
  geom_point() +
  scale_colour_manual(values=pal)
accept

proteo <- ggplot(best_fits_Proteobacteria, aes(x=Est.Tpk, y=Max.response, colour = Acceptor)) +
  geom_point()
proteo

tener <- ggplot(best_fits_Tenericutes, aes(x=Est.Tpk, y=Max.response, colour = Genus)) +
  geom_point()
tener

firmi <- ggplot(best_fits_Firmicutes, aes(x=Est.Tpk, y=Max.response, colour = Class, label = Species)) +
  geom_point()
firmi

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
