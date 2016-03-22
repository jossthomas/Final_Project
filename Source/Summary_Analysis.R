require(ggplot2)
require(reshape2)

setwd('C:\\Users\\Joss\\Dropbox\\Final_Project\\Computation\\Source')

data <- read.csv('../Data/Summary.csv')

View(data)

best_fits <- data[data$Plotted == "True",]
best_fits <- best_fits[best_fits$R_Squared > 0.75,]

a <- ggplot(best_fits, aes(x=T_pk, y=E_D)) +
     geom_point() +
     scale_y_continuous(limits = c(0, 100))
a

a1 <- ggplot(best_fits, aes(x=T_pk, y=E)) +
  geom_point() +
  scale_y_continuous(limits = c(0, 6))
a1


b <- ggplot(best_fits, aes(x=T_pk, y=Est.Tpk-Est.Tmin)) +
  geom_point()
b

f <- ggplot(best_fits, aes(x=T_pk, y=B0)) +
  geom_point() +
  scale_y_continuous(limits = c(-70, 50))
f

best_fits_growth <- best_fits[which(best_fits$Trait == 'Specific Growth Rate'|best_fits$Trait == 'Population Growth Rate'|best_fits$Trait == 'Specific Growth'),]
best_fits_growth <- best_fits_growth[best_fits_growth$Max.response < 200,]
best_fits_growth$Trange <- best_fits_growth$Est.Tpk - best_fits_growth$Est.Tmin


c <- ggplot(best_fits_growth, aes(x=Est.Tpk, y=Trange)) +
  geom_point()
c


d <- ggplot(best_fits_growth, aes(x=Est.Tpk, y=Max.response)) +
  geom_point()
d

View(data)

mp <- NULL
mapWorld <- borders("world", colour="gray50", fill="gray50") # create a layer of borders
mp <- ggplot(data) + mapWorld
mp <- mp + 
      geom_point(aes(x=Longitude, y=Latitude, colour=Phylum, alpha =0.1) , size=3) +
      geom_jitter(aes(x=1, y=1))
mp

test <- glm(Max.response~Est.Tpk+Phylum, data=best_fits_growth)
summary(test)
