library(readxl)
library(synergyfinder)
library(ggplot2)

SFDND41 <- read_excel("SF_R_DND41.xlsx")

head(test)   
source("Reshape.R")
dose.response.mat <- mReshapeData(SFDND41,
                                 data.type = "inhibition",
                                 impute = TRUE,
                                 noise = TRUE,
                                 correction = "non")

str(dose.response.mat)
PlotDoseResponse(dose.response.mat)

#PlotDoseResponse(dose.response.mat, save.file = TRUE)

#Other reference models can be chosen by setting the method parameter 
# as 'HSA','Loewe' , "ZIP" or 'Bliss'.
synergy.score <- CalculateSynergy(data = dose.response.mat,
                                  method = 'HSA')
str(synergy.score)

PlotSynergy(synergy.score, type = "all", save.file = TRUE)
