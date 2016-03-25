require(ggplot2)
require(reshape2)

setwd('C:\\Users\\Joss\\Dropbox\\Final_Project\\Computation\\Data')

data <- read.csv('../Data/Summary.csv')

View(data)

best_fits <- data[data$Plotted == "TRUE",]
best_fits <- best_fits[best_fits$R_Squared > 0.75,]
best_fits_growth <- best_fits[which(best_fits$Trait == 'Specific Growth Rate'),]
best_fits_growth$Trange <- best_fits_growth$Est.Tpk - best_fits_growth$Est.Tmin

hist <- ggplot(best_fits, aes(x=Est.Tpk)) +
  geom_histogram()
hist

a <- ggplot(best_fits_growth, aes(x=Est.Tpk, y=E_D, colour=Kingdom)) +
     geom_point() +
     scale_y_continuous(limits = c(0, 100))
a

a1 <- ggplot(best_fits_growth, aes(x=Est.Tpk, y=E, colour=Kingdom)) +
  geom_smooth(method='lm') +
  scale_y_continuous(limits = c(0, 6)) +
  geom_point()
a1


b <- ggplot(best_fits_growth, aes(x=Est.Tpk, y=Est.Tpk-Est.Tmin, colour=Kingdom)) +
  geom_smooth(method='lm') +
  geom_point()
b

f <- ggplot(best_fits_growth, aes(x=Est.Tpk, y=B0, colour=Kingdom)) +
  geom_smooth(method='lm') +
  geom_point() +
  scale_y_continuous(limits = c(-70, 50))
f


c <- ggplot(best_fits_growth, aes(x=Est.Tpk, y=Trange, colour=Kingdom)) +
  geom_smooth(method='lm') +
  geom_point()
c


d <- ggplot(best_fits_growth, aes(x=Est.Tpk, y=Max.response, colour=Phylum)) +
  #geom_smooth(method='lm') +
  scale_y_continuous(limits = c(0, 100)) +  
  geom_point()
d

require(stringr)

best_fits_growth[best_fits_growth==""] <- NA
acceptors <- best_fits_growth[complete.cases(best_fits_growth$Acceptor),]


accept <- ggplot(acceptors[acceptors$Acceptor  %in% c('Acetaldehyde', 'Acetyl-Phosphate', 'Pyruvate', 'Pyruvate ', 'Formate'),], aes(x=Est.Tpk, y=Max.response, colour=Acceptor)) +
  scale_y_continuous(limits = c(0, 100)) +
  geom_smooth(method='lm', alpha=0.2) +
  geom_point()
accept

test2 <- glm(Max.response ~ Est.Tpk:Acceptor:Phylum, data=acceptors)
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