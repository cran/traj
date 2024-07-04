## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----plottrajdata, echo = FALSE-----------------------------------------------
library(traj)
data(trajdata) 
dat <- trajdata[, -c(1:2)]
wA <- which(trajdata$Group == "A")
wB <- which(trajdata$Group == "B")
wC <- which(trajdata$Group == "C")
wD <- which(trajdata$Group == "D")

plot(x = 0, y = 0, xlim = c(1, 6), ylim = c(min(dat), max(dat) + 30), type = "n", ylab = "", xlab = "")

for(k in wA){
  lines(x = 1:6, y = dat[k, ], type = "l", col = "black ")
}

for(k in wB){
  lines(x = 1:6, y = dat[k, ], type = "l", col = "blue")
}

for(k in wC){
  lines(x = 1:6, y = dat[k, ], type = "l", col = "red")
}

for(k in wD){
  lines(x = 1:6, y = dat[k, ], type = "l", col = "green")
}

legend("topright",legend = c(paste("A (n = ", 50, ")", sep = ""), paste("B (n = ", 40, ")", sep = ""), paste("C (n = ", 30, ")", sep = ""), paste("D (n = ", 10, ")", sep = "")), col = c("black", "blue", "red", "green"), lty = 1)

## ----loadtraj-----------------------------------------------------------------
library(traj)
data(trajdata) 
head(trajdata)
dat <- trajdata[, -c(1,2)]

## ----ex1.step1----------------------------------------------------------------
step1 <- Step1Measures(Data = dat, measures = 1:18) 

summary(step1)

## ----ex1.step2a---------------------------------------------------------------
step2 <- Step2Selection(trajMeasures = step1) 
summary(step2)

## ----ex1.step2b---------------------------------------------------------------
print(step2)

## ----ex1.step3a---------------------------------------------------------------
library(cluster)
set.seed(1337)
step3 <- Step3Clusters(trajSelection = step2) 

## ----ex1.step3b---------------------------------------------------------------
par(mfrow = c(1, 1))
plot(step3, which.plots = 1, ask = FALSE)

## ----ex1.step3e, echo = FALSE-------------------------------------------------
par(mfrow = c(1, 1))
plot(step3, which.plots = 2, ask = FALSE)

## ----ex1.step3f, echo = FALSE-------------------------------------------------
par(mfrow = c(1, 1))
plot(step3, which.plots = 3, ask = FALSE)

## ----ex1.step3g, echo = FALSE-------------------------------------------------
par(mfrow = c(1, 1))
plot(step3, which.plots = 4, ask = FALSE)

## ----ex1.step3h, echo = FALSE-------------------------------------------------
par(mfrow = c(1, 1))
plot(step3, which.plots = 5, ask = FALSE)

## ----ex1.step3i, echo = FALSE-------------------------------------------------
par(mfrow = c(1, 1))
plot(step3, which.plots = 6, ask = FALSE)

## ----ex1.step3n---------------------------------------------------------------
color.pal <- palette.colors(palette = "Okabe-Ito", alpha = 1)
par(mfrow = c(1, 1))
for(k in 1:4){
  w <- which(step3$partition$Cluster == k)
  dat.w <- dat[w, ]
  plot(y = 0, x = 0, ylim = c(floor(min(dat)), ceiling(max(dat))), xlim = c(1,6), xlab="", ylab="", type="n", main = paste("Cluster ", k, " (n = ", step3$partition.summary[k], ")", sep = ""))
  for(i in 1:length(w)){
    lines(y = dat.w[i, ], x = 1:6, col = color.pal[k])
  }
}

