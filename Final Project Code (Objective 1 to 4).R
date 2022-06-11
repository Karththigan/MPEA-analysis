library(tidyverse)

################IMPORTANT NOTE###############
#PLEASE SET THE WORKING DIRECTORY TO A FOLDER ON YOUR COMPUTER
#THAT CONTAINTS MPEA_dataset.csv AND compositions.csv
#############################################
setwd("C:/Users/Owner/Desktop/Final Project")

#CLEANING MPEA DATASET AND MERGING COMPOSITIONAL DATA

untidy_mp <- read.csv("MPEA_dataset.csv", check.names = FALSE)

#deleting unnecessary data from data set
untidy_mp[ , c("Unnamed: 0","REFERENCE: doi_x","Internal Reference #",
               "IDENTIFIER: Reference ID","REFERENCE: year","REFERENCE: url",
               "REFERENCE: tag","REFERENCE: comment","Original DOI","Unnamed: 22",
               "REFERENCE: title","REFERENCE: doi_y","REFERENCE: doi")] <- NULL

#Remove "Property" from all column names
colnames(untidy_mp) <- gsub('PROPERTY: ','', colnames(untidy_mp))

#Remove brackets and text within from column names to make it easier to type column names
colnames(untidy_mp) <- gsub( " *\\(.*?\\) *", "", colnames(untidy_mp))

#create tidy data set, mp, and separate phases into primary, secondary, and other phases
mp <- separate(untidy_mp,"Type of phases", into=c("Primary phase", "Secondary phase", "Tertiary phase", "Quaternary phase", "Other phases"), sep="\\+")

#combine phases after secondary phase into single "Other phase" column
mp <-unite(mp,"Other phases", "Tertiary phase", "Quaternary phase", sep="+", na.rm = T)

#import composition data and join with tidy MPEA data
compositions <- read.csv("compositions.csv")
mp <- left_join(mp, compositions, by = c("FORMULA" = "Alloy.name"))




#OBJECTIVE 1: Determine which alloys will have strength of 500 MPa at 1000K
objective1 <- subset(x = mp, subset = `Test temperature` >= 1000 & YS >= 500)[,1]

#probability that an MPEA will have YS >=500 at 1000K
p_obj1 <- length(objective1)/length(unique(mp$FORMULA))




#OBJECTIVE 2: Which element has greatest effect on MoNbTi system

#Get all alloys containing MoNbTi
elements_to_find <- grepl("Mo", mp$FORMULA) &
  grepl("Nb", mp$FORMULA) &
  grepl("Ti", mp$FORMULA)

monbti_alloys <- mp[elements_to_find,]  

#Only keep columns with phase and compositions data
monbti_alloys <- monbti_alloys[, c(1:4,5,27:56)]

#Remove element composition columns with only NA values (i.e. elements not in alloy)
monbti_alloys <- monbti_alloys[,colSums(is.na(monbti_alloys))<nrow(monbti_alloys)]
monbti_alloys[is.na(monbti_alloys)] <- 0

#remove duplicates of alloy compositions
monbti_alloys <- monbti_alloys[!duplicated(monbti_alloys$FORMULA),]
monbti_alloys <- arrange(monbti_alloys, )

##Function that calculates probability that an MoNbTi alloy containing 
#a specific additive element, "element", is multiphase
prob_multiphase <- function(element){
  total_num_alloys <- nrow(monbti_alloys)
  alloys_with_element <- monbti_alloys[grepl(element,monbti_alloys$FORMULA),]
  num_multiphase <- sum(alloys_with_element$`Single/Multiphase` == "M")
  probability_multiphase <- num_multiphase / nrow(alloys_with_element)
  return(probability_multiphase)
}

#Create data frame listing probability of having a multiphase MoNbTi alloy
#depending on additive element
multiphases_probs <- data.frame(colnames(monbti_alloys),"Probability")
colnames(multiphases_probs) <- c("Additive element", "Probability of turning MoNbTi multiphase")
multiphases_probs <- multiphases_probs[-c(1:5),]
multiphases_probs[,2]= apply(multiphases_probs, 1, prob_multiphase)





#OBJECTIVE 3: YS and Grain Size

#get alloys from MPEA dataset that have multiple grainsize values
grainsize_ys <- mp[!is.na(mp$`grain size`), c("FORMULA","YS","Test temperature", "grain size")]

#get alloys from MPEA dataset that have multiple grain size values
freq_grainsize_ys <- as.data.frame(table(grainsize_ys$FORMULA))

multi_grainsizeonly <- left_join(grainsize_ys , freq_grainsize_ys, 
                          by = c("FORMULA" = "Var1"))

#only choosing alloys with at least 8 (grain size, YS) data points
multi_grainsizeonly <- multi_grainsizeonly[multi_grainsizeonly$Freq >=8, ]
multi_grainsizeonly <- multi_grainsizeonly[order(multi_grainsizeonly[,'FORMULA']), ]

#removing alloys with NA as YS
multi_grainsizeonly <- multi_grainsizeonly[!is.na(multi_grainsizeonly$YS),]

#correcting frequency of alloys in dataframe
correctedfreq <- as.data.frame(table(multi_grainsizeonly$FORMULA))
multi_grainsizeonly <- left_join(multi_grainsizeonly,correctedfreq, by = c("FORMULA" = "Var1"))

#removing incorrect frequency column and renaming column
multi_grainsizeonly$Freq.x <- NULL
colnames(multi_grainsizeonly)[5] <- "Freq"


#produce plot of YS vs grain size
ggplot(multi_grainsizeonly, aes(x=`grain size`, y= `YS`, color= `FORMULA`)) + 
  geom_point() +
  geom_smooth(method=lm, formula = y ~ x, se=F) + 
  labs(y = "Yield Strength (MPa)", x = "Grain Size (μm)", color= "FORMULA") + 
  theme(legend.position = c(0.7,0.85)) +
  theme(legend.key.size = unit(0.2, "cm"))


#get average grain size of each alloying composition 
avg_grainsize <- left_join(mp, as.data.frame(table(mp$FORMULA)), by = c("FORMULA" = "Var1"))
avg_grainsize <- avg_grainsize[avg_grainsize$Freq >= 2, ]
avg_grainsize <- aggregate(avg_grainsize$`grain size` ~ avg_grainsize$FORMULA, avg_grainsize, mean)
colnames(avg_grainsize) <- c("Alloy", "Mean Grain Size")


#create histogram of mean grain size
ggplot(avg_grainsize, aes(x = avg_grainsize$`Mean Grain Size`)) +
  geom_histogram(binwidth = 50, color="black", fill="darkseagreen") +
  labs(y= "Number of alloys", x = "Mean Grain Size (μm)")
avg_grainsize$Alloy <- as.character(avg_grainsize$Alloy)

#create boxplot of mean grain size histogram
ggplot(avg_grainsize, aes(x = `Mean Grain Size`)) + 
  geom_boxplot(notch=FALSE, outlier.colour="red", outlier.size=3)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  coord_flip()

#get median and mean
median(avg_grainsize$`Mean Grain Size`)
mean(avg_grainsize$`Mean Grain Size`)
quantile(avg_grainsize$`Mean Grain Size`)


#create multiple regression of YS ~ temp and grain size
mult_regression <- lm(YS ~ `Test temperature` +
                        `grain size`, data=mp)
summary(mult_regression)

#This is a very poor fit, because there is too many alloys of 
#same composition that have different YS values for both different grain
#size and different temp values

#Probability that a random alloy picked from
#MPEA data set has multi-phase structure given that it is an MoNbTi alloy
p_monbtialloy <- nrow(monbti_alloys[monbti_alloys$`Single/Multiphase` == "M", ]) / nrow(mp)
p_multiphase_and_monballoy <- length(monbti_alloys[monbti_alloys$`Single/Multiphase`== "M",]) / nrow(mp)
p_multiphase_given_mobnalloy <- p_multiphase_and_monballoy / p_monbtialloy 




#OBJECTIVE 4: Yield strength and temperature relationship

#Get average YS for each alloying composition at each test temperature
avg_ys_compandtemp <- aggregate(YS~ FORMULA+`Test temperature`, mp, mean)
avg_ys_comp <- aggregate(YS~ FORMULA, mp, mean)

#make data frame with avg YS for each composition and YS for varying
#temperatures for same composition
ys_and_comp <- left_join(avg_ys_comp, avg_ys_compandtemp, by="FORMULA")

#rename columns for clarity
names(ys_and_comp)[names(ys_and_comp) == 'YS.x'] <- 'AvgYs'
names(ys_and_comp)[names(ys_and_comp) == 'YS.y'] <- 'YS'

#determine which alloys have YS values at multiple temps

#create table of alloy frequencies
frequency_of_alloys <- as.data.frame(table(ys_and_comp$FORMULA))
colnames(frequency_of_alloys)[colnames(frequency_of_alloys)  == "Var1"] <- "Formula"

multi_ysonly <- left_join(ys_and_comp, frequency_of_alloys, 
                          by = c("FORMULA" = "Formula"))

#Only choosing alloys that have at least 8 (Temp,YS) data points; any less
#does not make sense to conduct linear regression with less data points since not enough data
# For reference, only 17 alloys have at least 7 (temp, YS) data points
# 10 alloys have at least 7, and only 7 alloys >= 8
# length(unique(multi_ysonly[multi_ysonly$Freq >= 6,]$FORMULA))
# [1] 17
# length(unique(multi_ysonly[multi_ysonly$Freq >= 7,]$FORMULA))
# [1] 10
# length(unique(multi_ysonly[multi_ysonly$Freq >= 8,]$FORMULA))
# [1] 7

multi_ysonly <- multi_ysonly[multi_ysonly$Freq >= 8,]

#Create scatterplot of YS vs temp with different color for each alloy system
ggplot(multi_ysonly, aes(x=`Test temperature`, y= `YS`, color= `FORMULA`)) + 
  geom_point() +
  geom_smooth(method=lm, formula = y ~ x, se=F) + 
  labs(y = "Yield Strength (MPa)", x = "Test Temperature (Celsius)", color= "Alloy") + 
  theme(legend.position = c(0.7,0.85)) +
  theme(legend.key.size = unit(0.2, "cm"))


#Function that extracts the y-intercept and slope of each alloy's regression line
ys_temp_geteqn <- function(formula){
  test <- lm(`YS` ~ `Test temperature`, 
             data = multi_ysonly[multi_ysonly$FORMULA == formula,])
  return(unname(summary(test)$coefficients[c(1,2),1]))
}

#Create data frame with least-squares coefficients for alloys in graph
lm_coeffs <- data.frame(unique(apply(unique(multi_ysonly[,1,drop=F]), 1, ys_temp_geteqn)))
colnames(lm_coeffs) = unique(multi_ysonly$FORMULA)
lm_coeffs <- as.data.frame(t(lm_coeffs))
colnames(lm_coeffs) = c("Y-intercept", "Slope")

#Get r^2 values for each linear regression
ys_temp_geteqn_rsqrd <- function(formula){
  test <- lm(`YS` ~ `Test temperature`, 
             data = multi_ysonly[multi_ysonly$FORMULA == formula,])
  return(summary(test)$r.squared)
}

#Create column with r squared values and rename
lm_coeffs$Formula <- seq.int(nrow(lm_coeffs)) 
lm_coeffs$Formula <- row.names(lm_coeffs) #changing row names to numbers
row.names(lm_coeffs) <- seq.int(nrow(lm_coeffs))
lm_coeffs <- lm_coeffs[,c("Formula", "Slope", "Y-intercept")]
lm_coeffs[,4] <- data.frame(unique(apply(unique(multi_ysonly[,1,drop=F]), 1, ys_temp_geteqn_rsqrd)))
colnames(lm_coeffs)[4] <- "R-squared"

#Calculate correlation coefficient of each linear fit for each alloy
lm_coeffs$`R value` <- '^'(lm_coeffs$`R-squared`, 1/2)
