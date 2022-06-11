#Objective 5: Water Quality Control Charts

################IMPORTANT NOTE###############
#PLEASE SET THE WORKING DIRECTORY TO A FOLDER ON YOUR COMPUTER
#THAT CONTAINS BeachesData.csv
#############################################
setwd("C:/Users/Owner/Desktop/Final Project")

library(qcc)
library(tidyverse)
wq <- read.csv("BeachesData.csv")

#remove advisory columns (columns with text) to facilitate analysis
advisory_columns <- select(wq,contains("advisory"))
wq <- wq[,!(names(wq) %in% colnames(advisory_columns))]

#remove last three collection date rows because they have no data
wq <- wq[-c(96,97,98),]

#create rational subgroups of 5 days (so that 19 subgroups )
Sample <- rep(1:as.numeric(nrow(wq)/5), each=5)
wq <- cbind(Sample, wq)

#Create df of Marie Curtis Park East Beach E.coli levels
wq_marie <- qcc.groups(wq$data.0.eColi, wq$Sample)

#Create df of Sunnyside Beach Beach E.coli levels
wq_sunny <- qcc.groups(wq$data.1.eColi, wq$Sample)

#Create df of Hanlan's Point Beach E.coli levels
wq_hanlan <- qcc.groups(wq$data.2.eColi, wq$Sample)

#Create df of Gibraltar Point Beach E.coli levels
wq_gibraltar <- qcc.groups(wq$data.3.eColi, wq$Sample)

#Create df of Centre Island Beach E.coli levels
wq_centre <- qcc.groups(wq$data.4.eColi, wq$Sample)

#Create df of Ward's Island Beach E.coli levels
wq_ward <- qcc.groups(wq$data.5.eColi, wq$Sample)

#Create df of Cherry Beach Beach E.coli levels
wq_cherry <- qcc.groups(wq$data.6.eColi, wq$Sample)

#Create df of Woodbine Beaches E.coli levels
wq_woodbine <- qcc.groups(wq$data.7.eColi, wq$Sample)

#Create df of Kew Balmy Beach E.coli levels
wq_kew <- qcc.groups(wq$data.8.eColi, wq$Sample)

#Create df of Bluffer's Beach Park E.coli levels
wq_bluffer <- qcc.groups(wq$data.8.eColi, wq$Sample)

# Produce different types of control charts for each beach
get_controlchart <- function(wq_beach,charttype){
  if (charttype == "xbar"){
    xbar <- qcc(wq_beach, type="xbar", std.dev = "UWAVE-R",
                title = paste("X-bar chart for", deparse(substitute(wq_beach))))
    return(xbar)
  } else if (charttype == "S") {
    chart <- qcc(wq_beach, type="S",
        title = paste("S chart for", deparse(substitute(wq_beach))))
    return(chart)
  } else if (charttype == "R") {
    return(qcc(wq_beach, type="R",
               title = paste("R chart for", deparse(substitute(wq_beach)))))
  }
}

#Returns the number of out of control points from xbar chart
# #UCL is 100 because ecoli levels >100 is deemed unsafe in Toronto
num_violations <- function(wq_beach){
  chart <- qcc(wq_beach, 
               type="xbar", 
               limits = c(0, 100),
               title = paste("X-bar chart for", deparse(substitute(wq_beach))))
  chart <- chart$statistics[chart$statistics > 100]
  return(length(chart))
}

num_violations(wq_bluffer)
num_violations(wq_centre)
num_violations(wq_cherry)
num_violations(wq_gibraltar)
num_violations(wq_hanlan)
num_violations(wq_kew)
num_violations(wq_marie)
num_violations(wq_sunny)
num_violations(wq_ward)
num_violations(wq_woodbine)

get_controlchart(wq_bluffer, "S")
get_controlchart(wq_centre, "S")
get_controlchart(wq_cherry, "S")
get_controlchart(wq_gibraltar, "S")
get_controlchart(wq_hanlan, "S")
get_controlchart(wq_kew, "S")
get_controlchart(wq_marie, "S")
get_controlchart(wq_sunny, "S")
get_controlchart(wq_ward, "S")
get_controlchart(wq_woodbine, "S")
