library(MASS)
library(drc)
library(readxl)
library(ggplot2)
x <- read.csv("/Users/Hui/Drug_synergy/DND41_Isobologram.csv", header = TRUE)
x <- x[,1:7]
x[,4:7] <- x[, 4:7] * 100

combination.vib <- x[x[, 2] != 0 & x[,3] != 0, 7]

A1.vib <- c(as.numeric(unlist(x[x[,2] == 0, 4])), 
            as.numeric(unlist(x[x[,2] == 0, 5])), 
            as.numeric(unlist(x[x[,2] == 0, 6])))/100
A1.vib[A1.vib >1] <- 1
A1.conc <- rep(as.numeric(as.character(unlist(x[x[,2] == 0, 3]))), 3)

A2.vib <- c(as.numeric(unlist(x[x[,3] == 0, 4])), 
             as.numeric(unlist(x[x[,3] == 0, 5])), 
             as.numeric(unlist(x[x[,3] == 0, 6])))/100
 
A2.vib[A2.vib > 1] <- 1
A2.conc <- rep(as.numeric(as.character(unlist(x[x[,3] == 0, 2]))), 3)

#### 3 parameters with lower bound 0 for ED50
Fors.LL.3 <- drm(A1.vib ~ A1.conc, fct = LL.3())
ED(Fors.LL.3, c(30, 50, 70), interval = "delta")
##### DND4 — Kayleigh GSK4716
Dex.LL.3 <- drm(A2.vib ~ A2.conc, fct = LL.3())
ED(Dex.LL.3, c(20, 30, 50, 70, 80), interval = "delta")

##### ##### DND4 — Kayleigh Dex
Fors.dose <- ED(Fors.LL.3, (100 - combination.vib), interval = "delta")
Dex.dose <- ED(Dex.LL.3, (100 - combination.vib), interval = "delta")


isob.data <- cbind( Dex.dose = Dex.dose[,1], Fors.dose = Fors.dose[,1], comb.Dex.dose = x[x[,2] != 0 & x[,3] != 0, 2], comb.Fors.dose = x[x[,2] != 0 & x[,3] != 0, 3])

isob.data <- transform(isob.data, CI= isob.data[,4]/isob.data[,2] + isob.data[,3] /isob.data[,1])

isob.data <- transform(isob.data, normalized.comb.Fors.dose =  comb.Fors.dose / Fors.dose, normalized.comb.Dex.dose = comb.Dex.dose / Dex.dose)

isob.data <- transform(isob.data, IC = gsub("e:1:", "",rownames(isob.data)))
colnames(isob.data)
colnames(isob.data) <- gsub("Fors", "GSK4716", colnames(isob.data))
rownames(isob.data) <- gsub("e:1:", "",rownames(isob.data))

#plot
gg<-ggplot(isob.data, aes(x=normalized.comb.GSK4716.dose, y= normalized.comb.Dex.dose))
gg1 <- gg + geom_point() + labs(title="Isobologram") + coord_fixed(ratio = 1)+ coord_cartesian(xlim=c(0,1.5), ylim=c(0, 1.5)) + geom_abline(intercept = 1, slope = -1, color="blue")
gg1
#output
write.table(isob.data, file ="ForIsobologramPlot-DND41.xls", sep="\t", row.names = FALSE)