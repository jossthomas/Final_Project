
rm(list=ls())

#---------------------||  Packages  ||---------------------

#install.packages('ggplot2')
#install.packages('scales')
#install.packages('Cairo')

require(ggplot2)
require(scales)
require(Cairo)
require(reshape)

#---------------------||  Functions and Themes  ||---------------------

main_theme <- theme_bw() + 
              theme(axis.text.x = element_text(size = 14),
                    axis.text.y = element_text(size = 14),
                    axis.title.y = element_text(size = 16),
                    axis.title.x = element_text(size = 16),
                    plot.title = element_text(size=18, vjust=1),
                    legend.text=element_text(size=14),
                    strip.text.x = element_text(size = 16))

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

check_overlaps = function(ES, EG){
  overlaps <- EG
  
  ES.upper <- ES$E[,1]
  ES.mean <- ES$E[,2]
  ES.lower <- ES$E[,3]
  
  overlaps$ES.max <- unlist(as.list(ES.upper))
  overlaps$ES.mean <- unlist(as.list(ES.mean))
  overlaps$ES.min <-unlist(as.list(ES.lower))

  overlaps$ES_GT_0 <- overlaps$ES.min > 0
  overlaps$EG_GT_0 <- overlaps$EG.min > 0

  overlaps$ES_EQ_EG <- ((overlaps$ES.min < overlaps$EG.min) & (overlaps$ES.max > overlaps$EG.min)) | 
                 ((overlaps$ES.min < overlaps$EG.max) & (overlaps$ES.max > overlaps$EG.max)) | 
                 ((overlaps$ES.min > overlaps$EG.min) & (overlaps$ES.max < overlaps$EG.max))
  overlaps$EG_GT_ES <- overlaps$ES.max < overlaps$EG.min

  overlaps$best_hypothesis <- 'Biochemical\nAdaptation'
  
  hotter_better <- (apply( t(t(overlaps[,c("ES_GT_0","EG_GT_0","ES_EQ_EG", "EG_GT_ES")]) == c(TRUE, TRUE, TRUE, FALSE)), 1, all))
  hotter_much_better <- (apply(t(t(overlaps[,c("ES_GT_0","EG_GT_0","ES_EQ_EG", "EG_GT_ES")]) == c(TRUE, TRUE, FALSE, TRUE)), 1, all))
  weak_adaption <- (apply( t(t(overlaps[,c("ES_GT_0","EG_GT_0","ES_EQ_EG", "EG_GT_ES")]) == c(TRUE, TRUE, FALSE, FALSE)), 1, all))
  
  overlaps$best_hypothesis[hotter_better] <- 'Hotter Is\nBetter'
  overlaps$best_hypothesis[hotter_much_better] <- 'Hotter Is\nBetter'
  overlaps$best_hypothesis[weak_adaption] <- 'Weak Biochemical\nAdaptation'

  return(overlaps)
}

#Bootstrap functions, these are really slow

resample_all <- function(x){
  boot_nums = c()
  len_x = length(x)
  for (i in 1:5000){
    boot_nums <- c(boot_nums, mean(sample(x, len_x, replace=TRUE)))
  }
  return( boot_nums )
} 

bootstrap_upper <- function(x){
  sample <- resample_all(x)
  quant = as.single(quantile(sample, 0.975))[1]
  return(quant)
}

bootstrap_mean <- function(x){
  sample <- resample_all(x)
  mn <- as.single(mean(sample))[1]
  return(mn)
}

bootstrap_lower <- function(x){
  sample <- resample_all(x)
  quant = as.single(quantile(sample, 0.025))[1]
  return(quant)
}

create_plot_tab <- function(x){
  SGroup <- x$Group
  SGroupNum <- as.integer(x$Group)
  
  Kingdom <- x$ConKingdom
  df <- data.frame(Kingdom, SGroup, SGroupNum)
  df <- rbind(df, df)
  
  EG <- replicate(length(SGroup), "EG")
  ES <- replicate(length(SGroup), "ES")
  
  df$Evar <- c(ES, EG)
  
  df$E.mean <- c(x$EG.mean, x$ES.mean)
  df$E.min <- c(x$EG.min, x$ES.min)
  df$E.max <- c(x$EG.max, x$ES.max)
  
  return(df)
}

#---------------------||  Load and Sort Data  ||---------------------

setwd('C://Users/Joss/Dropbox/Final_Project/Computation/Source')

data <- read.csv('../Data/summaries/summary_activation.csv')

best_fits <- data[data$Rank==1,]
best_fits <- best_fits[best_fits$E < 15,] #remove a poorly fitted outlier, not a good e estimate
best_fits_growth <- best_fits[best_fits$Trait == 'Specific Growth Rate',]
best_fits_mean_growth <- aggregate(cbind(Max.response, Est.Tpk, E)~Species:Trait*ConKingdom*ConPhylum*ConClass*ConOrder*ConFamily*ConGenus*TempPref, data=best_fits_growth,  FUN=gm_mean)

View(best_fits_mean_growth[best_fits_mean_growth$E == min(best_fits_mean_growth$E),])
table(best_fits$Model_name)

C1 <- best_fits_mean_growth[(best_fits_mean_growth$ConPhylum == 'Thermotogae'),]
View(C1)
C1P <- ggplot(C1, aes(x=Est.Tpk-273, y=Max.response, colour=ConClass, label=Species)) +
  geom_point() +
  geom_text(vjust=1)
C1P 
#---------------------||  Find Groups With 5+ Members ||---------------------

Phylum_activations
PhylumTab <- table(best_fits_mean_growth$ConPhylum)
PhylumData <- best_fits_mean_growth[best_fits_mean_growth$ConPhylum %in% names(PhylumTab[PhylumTab > 4]), ]

ClassTab <- table(best_fits_mean_growth[best_fits_mean_growth$ConClass != "",]$ConClass)
ClassData <- best_fits_mean_growth[best_fits_mean_growth$ConClass %in% names(ClassTab[ClassTab  > 4]), ]

OrderTab <- table(best_fits_mean_growth[best_fits_mean_growth$ConOrder != "",]$ConOrder)
OrderData <- best_fits_mean_growth[best_fits_mean_growth$ConOrder %in% names(OrderTab[OrderTab  > 4]), ]

FamilyTab <- table(best_fits_mean_growth[best_fits_mean_growth$ConFamily != "",]$ConFamily)
FamilyData <- best_fits_mean_growth[best_fits_mean_growth$ConFamily %in% names(FamilyTab[FamilyTab  > 4]), ]

GenusTab <- table(best_fits_mean_growth[best_fits_mean_growth$ConGenus != "",]$ConGenus)
GenusTab
GenusData <- best_fits_mean_growth[best_fits_mean_growth$ConGenus %in% names(GenusTab[GenusTab  > 4]), ]

Metabolismtab <- table(best_fits_mean_growth$Best_Guess)
MetabolismData <- best_fits_mean_growth[best_fits_mean_growth$Best_Guess %in% names(Metabolismtab[Metabolismtab > 4]), ]

TempTab <- table(best_fits_mean_growth$TempPref)
TempData <- best_fits_mean_growth[best_fits_mean_growth$TempPref %in% names(TempTab[TempTab > 4]), ]

KingdomLen <- aggregate(Est.Tpk~ConKingdom, FUN= length, data = best_fits_mean_growth)
PhylumLen <- aggregate(Est.Tpk~ConPhylum, FUN= length, data = PhylumData)
ClassLen <- aggregate(Est.Tpk~ConClass, FUN= length, data = ClassData)
OrderLen <- aggregate(Est.Tpk~ConOrder, FUN= length, data = OrderData)
FamilyLen <- aggregate(Est.Tpk~ConFamily, FUN= length, data = FamilyData)
GenusLen <- aggregate(Est.Tpk~ConGenus, FUN= length, data = GenusData)


#---------------------||  Separate EG for Bacteria and Archaea Respiratory Pathways ||---------------------

best_fits_mean_bac <- best_fits_growth[best_fits_growth$ConKingdom == 'Bacteria',]
best_fits_mean_arch <- best_fits_growth[best_fits_growth$ConKingdom == 'Archaea',]

agg_mean_bac <- aggregate(cbind(Max.response, Est.Tpk, E)~Species:Trait*ConKingdom*ConPhylum*ConClass*ConOrder*ConFamily*ConGenus*Best_Guess*TempPref, data=best_fits_mean_bac,  FUN=gm_mean)
agg_mean_arch <- aggregate(cbind(Max.response, Est.Tpk, E)~Species:Trait*ConKingdom*ConPhylum*ConClass*ConOrder*ConFamily*ConGenus*Best_Guess*TempPref, data=best_fits_mean_arch,  FUN=gm_mean)

Metabolismtab_bac <- table(agg_mean_bac$Best_Guess)
Metabolismtab_arch <- table(agg_mean_arch$Best_Guess)
Metabolismtab_bac
MetabolismData_bac <- best_fits_growth[best_fits_growth$Best_Guess %in% names(Metabolismtab_bac[Metabolismtab_bac > 4]), ]
MetabolismData_arch <- best_fits_growth[best_fits_growth$Best_Guess %in% names(Metabolismtab_arch[Metabolismtab_arch > 4]), ]
metabolism_join <- rbind(MetabolismData_bac, MetabolismData_arch)
metabolism_final <- aggregate(cbind(Max.response, Est.Tpk, E, Metabolism_Seperate_Activation.min, Metabolism_Seperate_Activation, Metabolism_Seperate_Activation.max)~Species:Trait*ConKingdom*ConPhylum*ConClass*ConOrder*ConFamily*ConGenus*Best_Guess*TempPref, data=metabolism_join,  FUN=gm_mean)

View(metabolism_final)
#---------------------||  Find ED and Limits ||---------------------

#Bootstrap the data - very slow!
Kingdom_E_SE <- aggregate(E~ConKingdom, FUN=function(x) c(upper=bootstrap_upper(x), mean=bootstrap_mean(x), lower=bootstrap_lower(x)), data=best_fits_mean_growth)
Phylum_E_SE <- aggregate(E~ConPhylum*ConKingdom, FUN=function(x) c(upper=bootstrap_upper(x), mean=bootstrap_mean(x), lower=bootstrap_lower(x)), data=PhylumData)
Class_E_SE <- aggregate(E~ConClass*ConKingdom, FUN=function(x) c(upper=bootstrap_upper(x), mean=bootstrap_mean(x), lower=bootstrap_lower(x)), data=ClassData)
Order_E_SE <- aggregate(E~ConOrder*ConKingdom, FUN=function(x) c(upper=bootstrap_upper(x), mean=bootstrap_mean(x), lower=bootstrap_lower(x)), data=OrderData)
Family_E_SE <- aggregate(E~ConFamily*ConKingdom, FUN=function(x) c(upper=bootstrap_upper(x), mean=bootstrap_mean(x), lower=bootstrap_lower(x)), data=FamilyData)
Genus_E_SE <- aggregate(E~ConGenus*ConKingdom, FUN=function(x) c(upper=bootstrap_upper(x), mean=bootstrap_mean(x), lower=bootstrap_lower(x)), data=GenusData)
Metabolism_E_all_SE <- aggregate(E~Best_Guess, FUN=function(x) c(upper=bootstrap_upper(x), mean=bootstrap_mean(x), lower=bootstrap_lower(x)), data=MetabolismData)
Metabolism_E_grouped_SE <- aggregate(E~Best_Guess:ConKingdom, FUN=function(x) c(upper=bootstrap_upper(x), mean=bootstrap_mean(x), lower=bootstrap_lower(x)), data=metabolism_final)
Temp_E_SE <- aggregate(E~TempPref:ConKingdom, FUN=function(x) c(upper=bootstrap_upper(x), mean=bootstrap_mean(x), lower=bootstrap_lower(x)), data=TempData)
Kingdom_E_SE
Phylum_E_SE
Class_E_SE
Order_E_SE
Family_E_SE
Genus_E_SE
Metabolism_E_all_SE
Metabolism_E_grouped_SE
Temp_E_SE
#---------------------||  Find EG ||---------------------
FamilyData_act <- best_fits_growth[best_fits_growth$ConFamily %in% names(FamilyTab[FamilyTab  > 5]), ]
GenusData_act <- best_fits_growth[best_fits_growth$ConGenus %in% names(GenusTab[GenusTab  > 5]), ]

Kingdom_activations <- aggregate(cbind(ConKingdom.E.min, ConKingdom_Activation, ConKingdom.E.max)~ConKingdom, FUN=mean, data=best_fits_growth)
Phylum_activations <- aggregate(cbind(ConPhylum.E.min, ConPhylum_Activation, ConPhylum.E.max)~ConPhylum*ConKingdom, FUN=mean, data=best_fits_growth)
Class_activations <- aggregate(cbind(ConClass.E.min, ConClass_Activation, ConClass.E.max)~ConClass*ConKingdom, FUN=mean, data=best_fits_growth)
Order_activations <- aggregate(cbind(ConOrder.E.min, ConOrder_Activation, ConOrder.E.max)~ConOrder*ConKingdom, FUN=mean, data=best_fits_growth)
Family_activations <- aggregate(cbind(ConFamily.E.min, ConFamily_Activation, ConFamily.E.max)~ConFamily*ConKingdom, FUN=mean, data=FamilyData_act)
Genus_activations <- aggregate(cbind(ConGenus.E.min, ConGenus_Activation, ConGenus.E.max)~ConGenus*ConKingdom, FUN=mean, data=GenusData_act)
Temp_activations <- aggregate(cbind(TempPref_Activation.min, TempPref_Activation, TempPref_Activation.max)~TempPref:ConKingdom, FUN=mean, data=best_fits_growth)

Metabolism_activations <- aggregate(cbind(Metabolism_Seperate_Activation.min, Metabolism_Seperate_Activation, Metabolism_Seperate_Activation.max)~Best_Guess:ConKingdom, FUN=mean, data=best_fits_growth)
Metabolism_activations_bac <- Metabolism_activations[Metabolism_activations$ConKingdom=='Bacteria',]
Metabolism_activations_arch <- Metabolism_activations[Metabolism_activations$ConKingdom=='Archaea',]
Metabolism_activations_filter_bac <- Metabolism_activations_bac[Metabolism_activations_bac$Best_Guess %in% names(Metabolismtab_bac[Metabolismtab_bac > 4]), ]
Metabolism_activations_filter_arch <- Metabolism_activations_arch[Metabolism_activations_arch$Best_Guess %in% names(Metabolismtab_bac[Metabolismtab_arch > 4]), ]
Metabolism_activations_sep <- rbind(Metabolism_activations_filter_bac, Metabolism_activations_filter_arch)
Metabolism_activations_sep <- Metabolism_activations_sep[order(Metabolism_activations_sep$ConKingdom),]

colnames_new_kingdom <- c('Group', 'ConKingdom', 'EG.min', 'EG.mean', 'EG.max')
colnames_new <- c('Group', 'EG.min', 'EG.mean', 'EG.max')

colnames(Kingdom_activations) <- colnames_new
colnames(Phylum_activations) <- colnames_new_kingdom
colnames(Class_activations) <- colnames_new_kingdom
colnames(Order_activations) <- colnames_new_kingdom
colnames(Family_activations) <- colnames_new_kingdom
colnames(Genus_activations) <- colnames_new_kingdom
colnames(Temp_activations) <- colnames_new_kingdom
colnames(Metabolism_activations_sep) <- colnames_new_kingdom
Kingdom_activations$ConKingdom <- Kingdom_activations$Group

#---------------------||  Find Overlaps ||---------------------
Kingdom_Overlaps <- check_overlaps(Kingdom_E_SE, Kingdom_activations)
Phylum_Overlaps <- check_overlaps(Phylum_E_SE, Phylum_activations)
Class_Overlaps <- check_overlaps(Class_E_SE, Class_activations)
Order_E_SE <- Order_E_SE[Order_E_SE$ConOrder!='Nostocales',]
Order_Overlaps <- check_overlaps(Order_E_SE, Order_activations)
Family_Overlaps <- check_overlaps(Family_E_SE, Family_activations)

Genus_Overlaps <- check_overlaps(Genus_E_SE, Genus_activations)
GenusTab
Metabolism_sep_Overlaps <- check_overlaps(Metabolism_E_grouped_SE, Metabolism_activations_sep)
Temp_Overlaps <- check_overlaps(Temp_E_SE, Temp_activations)

Kingdom_Overlaps$level <- 'Kingdom'
Phylum_Overlaps$level <- 'Phylum'
Class_Overlaps$level <- 'Class'
Order_Overlaps$level <- 'Order'
Family_Overlaps$level <- 'Family'
Genus_Overlaps$level <- 'Genus'
Metabolism_sep_Overlaps$level <- 'Metabolism'
Temp_Overlaps$level <- 'Temperature'

Kingdom_Overlaps
Phylum_Overlaps
Class_Overlaps
Order_Overlaps
Family_Overlaps
Genus_Overlaps
Metabolism_sep_Overlaps
Temp_Overlaps


#---------------------||  Count At Each Level  ||---------------------
Kingdom_counts <- as.data.frame(table(Kingdom_Overlaps[,c('level','best_hypothesis', 'ConKingdom')]))
Phylum_counts <- as.data.frame(table(Phylum_Overlaps[,c('level','best_hypothesis', 'ConKingdom')]))
Class_counts <- as.data.frame(table(Class_Overlaps[,c('level','best_hypothesis', 'ConKingdom')]))
Order_counts <- as.data.frame(table(Order_Overlaps[,c('level','best_hypothesis', 'ConKingdom')]))
Family_counts <- as.data.frame(table(Family_Overlaps[,c('level','best_hypothesis', 'ConKingdom')]))
Genus_counts <- as.data.frame(table(Genus_Overlaps[,c('level','best_hypothesis', 'ConKingdom')]))
Metabolism_counts <- as.data.frame(table(Metabolism_Overlaps[,c('level','best_hypothesis')]))
Metabolism_sep_counts <- as.data.frame(table(Metabolism_sep_Overlaps[,c('level','best_hypothesis', 'ConKingdom')]))
Temp_counts <-as.data.frame(table(Temp_Overlaps[,c('level','best_hypothesis', 'ConKingdom')]))

table(Kingdom_Overlaps$EG.mean > 0)
table(Phylum_Overlaps$EG.mean > 0)
table(Class_Overlaps$EG.mean > 0)
table(Order_Overlaps$EG.mean > 0)
table(Family_Overlaps$EG.mean > 0)
table(Genus_Overlaps$EG.mean > 0)


total <- merge_all(list(Kingdom_counts,Phylum_counts, Class_counts, Order_counts, Family_counts, Genus_counts, Metabolism_sep_counts))
total 

#---------------------||  Count At Each Level > 15 Datapoints  ||---------------------

Kingdom_Overlaps_N <- Kingdom_Overlaps[order(Kingdom_Overlaps$Group),]
Phylum_Overlaps_N <- Phylum_Overlaps[order(Phylum_Overlaps$Group),]
Class_Overlaps_N <- Class_Overlaps[order(Class_Overlaps$Group),]
Order_Overlaps_N <- Order_Overlaps[order(Order_Overlaps$Group),]
Family_Overlaps_N <- Family_Overlaps[order(Family_Overlaps$Group),]
Genus_Overlaps_N <- Genus_Overlaps[order(Genus_Overlaps$Group),]

Kingdom_Overlaps_N$Nvar <- KingdomLen$Est.Tpk
Phylum_Overlaps_N$Nvar <- PhylumLen$Est.Tpk
Class_Overlaps_N$Nvar <- ClassLen$Est.Tpk
Order_Overlaps_N$Nvar <- OrderLen$Est.Tpk
Family_Overlaps_N$Nvar <- FamilyLen$Est.Tpk
Genus_Overlaps_N$Nvar <- GenusLen$Est.Tpk

Kingdom_Overlaps_N <- Kingdom_Overlaps_N[Kingdom_Overlaps_N$Nvar > 19,]
Phylum_Overlaps_N <- Phylum_Overlaps_N[Phylum_Overlaps_N$Nvar > 19,]
Class_Overlaps_N <- Class_Overlaps_N[Class_Overlaps_N$Nvar > 19,]
Order_Overlaps_N <- Order_Overlaps_N[Order_Overlaps_N$Nvar > 19,]
Family_Overlaps_N <- Family_Overlaps_N[Family_Overlaps_N$Nvar > 19,]
Genus_Overlaps_N <- Genus_Overlaps_N[Genus_Overlaps_N$Nvar > 19,]

Kingdom_counts_N <- as.data.frame(table(Kingdom_Overlaps_N[,c('level','best_hypothesis', 'ConKingdom')]))
Phylum_counts_N <- as.data.frame(table(Phylum_Overlaps_N[,c('level','best_hypothesis', 'ConKingdom')]))
Class_counts_N <- as.data.frame(table(Class_Overlaps_N[,c('level','best_hypothesis', 'ConKingdom')]))
Order_counts_N <- as.data.frame(table(Order_Overlaps_N[,c('level','best_hypothesis', 'ConKingdom')]))
Family_counts_N <- as.data.frame(table(Family_Overlaps_N[,c('level','best_hypothesis', 'ConKingdom')]))
Genus_counts_N <- as.data.frame(table(Genus_Overlaps_N[,c('level','best_hypothesis', 'ConKingdom')]))


total_N <- merge_all(list(Kingdom_counts_N,Phylum_counts_N, Class_counts_N, Order_counts_N, Family_counts_N, Genus_counts_N))
total_N 

#---------------------||  Barplot of Percentages at Each level  ||---------------------
reorder(total$best_hypothesis, c('Hotter Is\nBetter', 'Weak Biochemical\nAdaptation', 'Biochemical\nAdaptation'))

hypotheses_bp <- ggplot(total, aes(x=level, y=Freq, fill=factor(best_hypothesis))) +
  geom_bar(stat="identity", position = "fill") +
  ylab('Frequency') +
  xlab('Analysis Level') +
  main_theme +
  scale_fill_manual('', 
                    breaks=c('Hotter Is\nBetter', 'Weak Biochemical\nAdaptation', 'Biochemical\nAdaptation'),
                    values=c('Hotter Is\nBetter'='#f8766d', 'Weak Biochemical\nAdaptation'='#619cff', 'Biochemical\nAdaptation'='#00ba38'),
                    labels=c('Hotter Is\nBetter', 'Weak Thermal\nConstraint', 'Biochemical\nAdaptation')) +  
  theme(legend.position = 'bottom')
hypotheses_bp 

hypotheses_bp_kingdom <- ggplot(total, aes(x=level, y=Freq, fill=best_hypothesis)) +
  geom_bar(stat="identity", position = "fill") +
  ylab('Frequency') +
  xlab('Analysis Level') +
  scale_fill_manual('', 
                    breaks=c('Hotter Is\nBetter', 'Weak Biochemical\nAdaptation', 'Biochemical\nAdaptation'),
                    values=c('Hotter Is\nBetter'='#f8766d', 'Weak Biochemical\nAdaptation'='#619cff', 'Biochemical\nAdaptation'='#00ba38'),
                    labels=c('Hotter Is\nBetter', 'Weak Thermal\nConstraint', 'Biochemical\nAdaptation')) +  
  main_theme +
  facet_wrap(~ConKingdom, ncol=1) +
  theme(legend.position = 'bottom')
hypotheses_bp_kingdom

ggsave(file = '../Results/Maxima_fits/Fit_Stats/best_hypotheses.png', hypotheses_bp, type="cairo-png", width = 12, height = 6)
ggsave(file = '../Results/Maxima_fits/Fit_Stats/best_hypothese_kingdom.png', hypotheses_bp_kingdom, type="cairo-png", width = 12, height = 7)

hypotheses_bp_N <- ggplot(total_N, aes(x=level, y=Freq, fill=factor(best_hypothesis))) +
  geom_bar(stat="identity", position = "fill") +
  ylab('Frequency') +
  xlab('Analysis Level') +
  main_theme +
  scale_fill_manual('', 
                    breaks=c('Hotter Is\nBetter', 'Weak Biochemical\nAdaptation', 'Biochemical\nAdaptation'),
                    values=c('Hotter Is\nBetter'='#f8766d', 'Weak Biochemical\nAdaptation'='#619cff', 'Biochemical\nAdaptation'='#00ba38'),
                    labels=c('Hotter Is\nBetter', 'Weak Thermal\nConstraint', 'Biochemical\nAdaptation')) +  
  theme(legend.position = 'bottom')
hypotheses_bp_N

#---------------------||  Overlap Plots  ||---------------------
Sys.setlocale("LC_CTYPE", "Latvian") #Needed to make E appear

dodge = position_dodge(width=0.18)
dodgewide = position_dodge(width=0.3)

kingdom_plt_back<- data.frame(SGroup=best_fits_mean_growth$ConKingdom, SGroupNum=as.integer(best_fits_mean_growth$ConKingdom), E=best_fits_mean_growth$E)
phylum_plt_back<- data.frame(SGroup=PhylumData$ConPhylum, SGroupNum=as.numeric(PhylumData$ConPhylum), E=PhylumData$E, Kingdom=PhylumData$ConKingdom)
class_plt_back<- data.frame(SGroup=ClassData$ConClass, SGroupNum=as.integer(ClassData$ConClass), E=ClassData$E, Kingdom=ClassData$ConKingdom)
order_plt_back<- data.frame(SGroup=OrderData$ConOrder, SGroupNum=as.integer(OrderData$ConOrder), E=OrderData$E, Kingdom=OrderData$ConKingdom)
family_plt_back<- data.frame(SGroup=FamilyData$ConFamily, SGroupNum=as.integer(FamilyData$ConFamily), E=FamilyData$E, Kingdom=FamilyData$ConKingdom)
genus_plt_back<- data.frame(SGroup=GenusData$ConGenus, SGroupNum=as.integer(GenusData$ConGenus), E=GenusData$E, Kingdom=GenusData$ConKingdom)
metabolism_plt_back <- data.frame(SGroup=metabolism_final$Best_Guess, SGroupNum=as.integer(metabolism_final$Best_Guess), E=metabolism_final$E, Kingdom=metabolism_final$ConKingdom)                         
Temp_plt_back <- data.frame(SGroup=TempData$TempPref, E=TempData$E, Kingdom=TempData$ConKingdom)  

kingdom_plt <- create_plot_tab(Kingdom_Overlaps)
phylum_plt <- create_plot_tab(Phylum_Overlaps)
class_plt <- create_plot_tab(Class_Overlaps)
order_plt <- create_plot_tab(Order_Overlaps)
family_plt <- create_plot_tab(Family_Overlaps)
genus_plt <- create_plot_tab(Genus_Overlaps)
metabolism_sep_plt <- create_plot_tab(Metabolism_sep_Overlaps)
Temp_plt <- create_plot_tab(Temp_Overlaps)

Kingdom_Overlap_Plot <- ggplot(kingdom_plt, aes(x = SGroup)) +
  geom_point(data=kingdom_plt_back, aes(y=E, colour='3'), alpha = 0.7, position = position_jitter(width = 0.3, height = 0.0)) +
  geom_errorbar(aes(ymax = E.max, ymin = E.min, group=Evar), width=0.15, size=1, position=dodge) +
  geom_point(aes(y=E.mean, colour = Evar, group=Evar), size = 4, position=dodge) +
  ylab('E') +
  xlab('Kingdom') +
  scale_colour_manual(guide='legend',
                      values =c('ES'='navyblue',
                                'EG'='dodgerblue',
                                '3'='chartreuse3'),
                      labels = c(expression('E'['S']),
                                 quote('\U0112'[S]),
                                 expression('E'['G']))
  ) +
  main_theme +
  theme(legend.title=element_blank(),
        legend.position='bottom')


Phylum_Overlap_Plot <- ggplot(phylum_plt, aes(x = factor(SGroup))) +
  geom_point(data=phylum_plt_back, aes(y=E, colour='3'), alpha = 0.7, position = position_jitter(width = 0.3, height = 0.0)) +
  geom_errorbar(aes(ymax = E.max, ymin = E.min, group=Evar), width=0.15, size=1, position=dodge) +
  geom_point(aes(y=E.mean, colour = Evar, group=Evar), size = 4, position=dodge) +
  ylab('E') +
  xlab('Phylum') +
  coord_cartesian(ylim=c(-3, 8)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  scale_colour_manual(guide='legend',
                      values =c('ES'='navyblue',
                                'EG'='dodgerblue',
                                '3'='chartreuse3'),
                      labels = c(expression('E'['S']),
                                 quote('\U0112'[S]),
                                 expression('E'['G']))
  ) +
  main_theme +
  theme(legend.title=element_blank(),
        legend.position='bottom',
        axis.text.x = element_text(angle = 45, hjust = 1))
Phylum_Overlap_Plot

Class_Overlap_Plot <- ggplot(class_plt, aes(x = SGroup)) +
  geom_point(data=class_plt_back, aes(y=E, colour='3'), alpha = 0.7, position = position_jitter(width = 0.3, height = 0.0)) +
  geom_errorbar(aes(ymax = E.max, ymin = E.min, group=Evar), width=0.15, size=1, position=dodge) +
  geom_point(aes(y=E.mean, colour = Evar, group=Evar), size = 4, position=dodge) +
  ylab('E') +
  xlab('Class') +
  coord_cartesian(ylim=c(-5,10)) +
  facet_wrap(~Kingdom, scales="free_x", ncol=1) +
  scale_colour_manual(guide='legend',
                      values =c('ES'='navyblue',
                                'EG'='dodgerblue',
                                '3'='chartreuse3'),
                      labels = c(expression('E'['S']),
                                 quote('\U0112'[S]),
                                 expression('E'['G']))
  ) +
  main_theme +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  theme(legend.title=element_blank(),
        legend.position='bottom',
        axis.text.x = element_text(angle = 45, hjust = 1))
Class_Overlap_Plot

Order_Overlap_Plot <- ggplot(order_plt, aes(x = SGroup)) +
  geom_point(data=order_plt_back, aes(y=E, colour='3'), alpha = 0.7, position = position_jitter(width = 0.3, height = 0.0)) +
  geom_errorbar(aes(ymax = E.max, ymin = E.min, group=Evar), width=0.15, size=1, position=dodge) +
  geom_point(aes(y=E.mean, colour = Evar, group=Evar), size = 4, position=dodge) +
  ylab('E') +
  xlab('Order') +
  coord_cartesian(ylim=c(-5,10)) +
  facet_wrap(~Kingdom, scales="free", ncol=1) +
  scale_colour_manual(guide='legend',
                      values =c('ES'='navyblue',
                                'EG'='dodgerblue',
                                '3'='chartreuse3'),
                      labels = c(expression('E'['S']),
                                 quote('\U0112'[S]),
                                 expression('E'['G']))
  ) +
  main_theme +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  theme(legend.title=element_blank(),
        legend.position='bottom',
        axis.text.x = element_text(angle = 45, hjust = 1))

Family_Overlap_Plot <- ggplot(family_plt, aes(x = SGroup)) +
  geom_point(data=family_plt_back, aes(y=E, colour='3'), alpha = 0.7, position = position_jitter(width = 0.3, height = 0.0)) +
  geom_errorbar(aes(ymax = E.max, ymin = E.min, group=Evar), width=0.15, size=1, position=dodge) +
  geom_point(aes(y=E.mean, colour = Evar, group=Evar), size = 4, position=dodge) +
  ylab('E') +
  xlab('Family') +
  coord_cartesian(ylim=c(-5,10)) +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  facet_wrap(~Kingdom, scales="free", ncol=1) +
  scale_colour_manual(guide='legend',
                      values =c('ES'='navyblue',
                                'EG'='dodgerblue',
                                '3'='chartreuse3'),
                      labels = c(expression('E'['S']),
                                 quote('\U0112'[S]),
                                 expression('E'['G']))) +
  main_theme +
  theme(legend.title=element_blank(),
        legend.position='bottom',
        axis.text.x = element_text(angle = 45, hjust = 1))


Genus_Overlap_Plot <- ggplot(genus_plt, aes(x = SGroup)) +
  geom_point(data=genus_plt_back, aes(y=E, colour='3'), alpha = 0.7, position = position_jitter(width = 0.3, height = 0.0)) +
  geom_errorbar(aes(ymax = E.max, ymin = E.min, group=Evar), width=0.15, size=1, position=dodge) +
  geom_point(aes(y=E.mean, colour = Evar, group=Evar), size = 4, position=dodge) +
  ylab('E') +
  xlab('Genus') +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  coord_cartesian(ylim=c(-5,10)) +
  facet_wrap(~Kingdom, scales="free_x", ncol=1) +
  scale_colour_manual(guide='legend',
                      values =c('ES'='navyblue',
                                'EG'='dodgerblue',
                                '3'='chartreuse3'),
                      labels = c(expression('E'['S']),
                                 quote('\U0112'[S]),
                                 expression('E'['G']))
  ) +
  main_theme +
  theme(legend.title=element_blank(),
        legend.position='bottom',
        axis.text.x = element_text(angle = 45, hjust = 1))


Metabolism_sep_Overlap_Plot <- ggplot(metabolism_sep_plt, aes(x = SGroup)) +
  geom_point(data=metabolism_plt_back, aes(y=E, colour='3'), alpha = 0.7, position = position_jitter(width = 0.3, height = 0.0)) +  
  geom_errorbar(aes(ymax = E.max, ymin = E.min, group=Evar), width=0.15, size=1, position=dodge) +
  geom_point(aes(y=E.mean, colour = Evar, group=Evar), size = 4, position=dodge) +
  ylab('E') +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  scale_x_discrete('Energy Generating Reaction', 
                   breaks=c('Aerobic', 'C16', 'C19', 'D1', 'C1', 'Fermentation', 'PhotosystemII'),
                   labels=c('Oxygen\nReduction', 'Sulphur\nOxidation', 'Sulphur\nReduction', 'Methanogenesis',
                            'Sulfate\nReduction', 'Fermentation', 'Photophosphorylation')) +
  
  facet_wrap(~Kingdom, scales="free", ncol=1) +
  scale_colour_manual(guide='legend',
                      values =c('ES'='navyblue',
                                'EG'='dodgerblue',
                                '3'='chartreuse3'),
                      labels = c(expression('E'['S']),
                                 quote('\U0112'[S]),
                                 expression('E'['G']))
  ) +
  main_theme +
  theme(legend.title=element_blank(),
        legend.position='bottom',
        axis.text.x = element_text(angle = 45, hjust = 1))

Temp_overlap_plot <- ggplot(Temp_plt, aes(x = SGroup)) +
  geom_point(data=Temp_plt_back, aes(y=E, colour='3'), alpha = 0.7, position = position_jitter(width = 0.3, height = 0.0)) +  
  geom_errorbar(aes(ymax = E.max, ymin = E.min, group=Evar), width=0.15, size=1, position=dodge) +
  geom_point(aes(y=E.mean, colour = Evar, group=Evar), size = 4, position=dodge) +
  ylab('E') +
  xlab('Temperature Group') +
  geom_hline(aes(yintercept=0), linetype='dotted') +
  facet_wrap(~Kingdom, scales="free", ncol=1) +
  scale_colour_manual(guide='legend',
                      values =c('ES'='navyblue',
                                'EG'='dodgerblue',
                                '3'='chartreuse3'),
                      labels = c(expression('E'['S']),
                                 quote('\U0112'[S]),
                                 expression('E'['G']))
  ) +
  main_theme +
  theme(legend.title=element_blank(),
        legend.position='bottom')

Kingdom_Overlap_Plot
Phylum_Overlap_Plot
Class_Overlap_Plot
Order_Overlap_Plot
Family_Overlap_Plot
Genus_Overlap_Plot
Metabolism_sep_Overlap_Plot
Temp_overlap_plot

#---------------------||  Save Overlaps  ||---------------------

ggsave(file = '../Results/Maxima_fits/Fit_Stats/1_Kingdom_Overlap_Plot.png', Kingdom_Overlap_Plot, width = 6, height = 6)
ggsave(file = '../Results/Maxima_fits/Fit_Stats/1_Phylum_Overlap_Plot.png', Phylum_Overlap_Plot, width = 12, height = 6)
ggsave(file = '../Results/Maxima_fits/Fit_Stats/1_Class_Overlap_Plot.png', Class_Overlap_Plot, width = 12, height = 10)
ggsave(file = '../Results/Maxima_fits/Fit_Stats/1_Order_Overlap_Plot.png', Order_Overlap_Plot, width = 12, height = 8)
ggsave(file = '../Results/Maxima_fits/Fit_Stats/1_Family_Overlap_Plot.png', Family_Overlap_Plot, width = 12, height = 8)
ggsave(file = '../Results/Maxima_fits/Fit_Stats/1_Genus_Overlap_Plot.png', Genus_Overlap_Plot, width = 12, height = 8)
ggsave(file = '../Results/Maxima_fits/Fit_Stats/1_Metabolism_sep_Overlap_Plot.png', Metabolism_sep_Overlap_Plot, width = 8, height = 8)
ggsave(file = '../Results/Maxima_fits/Fit_Stats/1_Temp_Overlap_Plot.png', Temp_overlap_plot, width = 8, height = 8)

#--------------------|| Find Minima and Maxima of Each Grouping ||--------------------

mesothermo <- function(x){
  group <- x[,1]
  kingdom <- x[,2]
  mins <- as.vector(x$Est.Tpk[,3])
  maxs <- as.vector(x$Est.Tpk[,1])
  all_thermo <- mins > (273.15+50)
  all_meso <- maxs < (273.15+50)
  
  ret = data.frame(group, kingdom, mins, maxs)
  ret$pref <- 'Mixed'
  ret$pref[all_thermo] <- 'Thermophilic'
  ret$pref[all_meso] <- 'Mesophilic'

  return(ret)
}

Kingdom_Min_Max <- aggregate(Est.Tpk~ConKingdom+ConKingdom, FUN=function(x) c(upper=quantile(x, .95), med=median(x), lower=quantile(x, .05)), data=best_fits_mean_growth)
Phylum_Min_Max  <- aggregate(Est.Tpk~ConPhylum*ConKingdom, FUN=function(x) c(upper=quantile(x, .95), med=median(x), lower=quantile(x, .05)), data=PhylumData)
Class_Min_Max  <- aggregate(Est.Tpk~ConClass*ConKingdom, FUN=function(x) c(upper=quantile(x, .95), med=median(x), lower=quantile(x, .05)), data=ClassData)
Order_Min_Max  <- aggregate(Est.Tpk~ConOrder*ConKingdom, FUN=function(x) c(upper=quantile(x, .95), med=median(x), lower=quantile(x, .05)), data=OrderData)
Family_Min_Max <- aggregate(Est.Tpk~ConFamily*ConKingdom, FUN=function(x) c(upper=quantile(x, .95), med=median(x), lower=quantile(x, .05)), data=FamilyData)
Genus_Min_Max  <- aggregate(Est.Tpk~ConGenus*ConKingdom, FUN=function(x) c(upper=quantile(x, .95), med=median(x), lower=quantile(x, .05)), data=GenusData)
Metabolism_Min_Max  <- aggregate(Est.Tpk~Best_Guess*ConKingdom, FUN=function(x) c(upper=quantile(x, .95), med=median(x), lower=quantile(x, .05)), data=metabolism_final)

Kingdom_Min_Max$bonus <- Kingdom_Min_Max$ConKingdom
Kingdom_Min_Max <- Kingdom_Min_Max[,c('ConKingdom', 'bonus', 'Est.Tpk')]
Kingdom_Min_Max <- mesothermo(Kingdom_Min_Max)
Phylum_Min_Max <- mesothermo(Phylum_Min_Max)
Class_Min_Max <- mesothermo(Class_Min_Max)
Order_Min_Max <- mesothermo(Order_Min_Max)
Family_Min_Max <- mesothermo(Family_Min_Max)
Genus_Min_Max <- mesothermo(Genus_Min_Max)
Metabolism_Min_Max <- mesothermo(Metabolism_Min_Max)

Kingdom_Min_Max$best_hypothesis <- Kingdom_Overlaps$best_hypothesis
Phylum_Min_Max$best_hypothesis <- Phylum_Overlaps$best_hypothesis
Class_Min_Max$best_hypothesis <- Class_Overlaps$best_hypothesis
Order_Min_Max$best_hypothesis <- Order_Overlaps$best_hypothesis
Family_Min_Max$best_hypothesis <- Family_Overlaps$best_hypothesis
Genus_Min_Max$best_hypothesis <- Genus_Overlaps$best_hypothesis
Metabolism_Min_Max$best_hypothesis <- Metabolism_sep_Overlaps$best_hypothesis


Kingdom_Min_Max$level <- 'Kingdom'
Phylum_Min_Max$level <- 'Phylum'
Class_Min_Max$level <- 'Class'
Order_Min_Max$level <- 'Order'
Family_Min_Max$level <- 'Family'
Genus_Min_Max$level <- 'Genus'
Metabolism_Min_Max$level <- 'Metabolism'

Kingdom_Min_Max

all_by_pref <- rbind(Kingdom_Min_Max, Phylum_Min_Max, Class_Min_Max, Order_Min_Max, Family_Min_Max, Genus_Min_Max, Metabolism_Min_Max)
all_by_pref$level2 <- factor(all_by_pref$level, levels=c('Kingdom','Phylum', 'Class', 'Order', 'Family', 'Genus', 'Metabolism'))


by_pref <- ggplot(all_by_pref, aes(x=pref, fill=best_hypothesis)) +
  geom_bar()+
  main_theme +
  xlab('Group Temperature Preference') +
  ylab('Count of Best Hypothesis') +
  scale_y_continuous(breaks=c(0,5,10,15,20,25,30))  +
  ylab('Count of Best Hypothesis') +
  theme(legend.position='bottom', legend.title=element_blank())
by_pref

by_pref_wrap <- ggplot(all_by_pref, aes(x=pref, fill=best_hypothesis)) +
  geom_bar()+
  main_theme +
  xlab('Group Temperature Preference') +
  ylab('Count of Best Hypothesis') +
  ylab('Count of Best Hypothesis') +
  scale_y_continuous(breaks=c(0,2,4,6,8,10))  +
  theme(legend.position='bottom', 
        legend.title=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual('', 
                    breaks=c('Hotter Is\nBetter', 'Weak Biochemical\nAdaptation', 'Biochemical\nAdaptation'),
                    values=c('Hotter Is\nBetter'='#f8766d', 'Weak Biochemical\nAdaptation'='#619cff', 'Biochemical\nAdaptation'='#00ba38'),
                    labels=c('Hotter Is\nBetter', 'Weak Thermal\nConstraint', 'Biochemical\nAdaptation')) +  
  facet_grid(.~level2)
by_pref_wrap

by_pref_kingdom <- ggplot(all_by_pref, aes(x=pref, fill=best_hypothesis)) +
  geom_bar()+
  main_theme +
  xlab('Group Temperature Preference') +
  ylab('Count of Best Hypothesis') +
  scale_y_continuous(breaks=c(0,2,4,6,8))  +
  theme(legend.position='bottom', 
        legend.title=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(kingdom~level2)
by_pref_kingdom

ggsave(file = '../Results/Maxima_fits/Fit_Stats/3_Temp_Pref.png', by_pref, width = 6, height = 6)
ggsave(file = '../Results/Maxima_fits/Fit_Stats/3_Temp_Group.png', by_pref_wrap, width = 12, height = 6)
ggsave(file = '../Results/Maxima_fits/Fit_Stats/3_Temp_Kingdom.png', by_pref_kingdom, width = 12, height = 6)

#----------------------- E vs TPK ---------------------------------
cor.test(best_fits_mean_growth$E, best_fits_mean_growth$Est.Tpk)
cor.test(log(best_fits_mean_growth$E), log(best_fits_mean_growth$Est.Tpk))

ETPK <- ggplot(best_fits_mean_growth, aes(x=Est.Tpk, y=log(E))) +
  geom_point(size=3, alpha=0.5, colour='chartreuse3') +
  xlab(expression('T'['pk']*' (K)')) +
  main_theme +
  geom_smooth(method='lm')
ETPK

qqnorm(log(1+best_fits_mean_growth$E))

ggsave(file = '../Results/Maxima_fits/Fit_Stats/Evals.png', ETPK, type="cairo-png", width = 6, height = 6)

