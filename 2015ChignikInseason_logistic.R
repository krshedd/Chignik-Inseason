## Testing with 2014 data
# Read in values from Birch's 2014 logistic worksheet in excel
sd <- as.numeric(readClipboard())
sd <- sd[!is.na(sd)]

day <- as.numeric(readClipboard())
day <- day[!is.na(day)]

est <- readClipboard()
est <- est[!est==""]
est <- as.numeric(sapply(est, function(i) unlist(strsplit(x = i, split = "%"))[1]))/100

df <- data.frame(day, est, sd)
df$invVar <- 1 / (df$sd^2)

# Determine GLM logistic fit based on IRWLS method
fit <- glm(est ~ day, data = df, weights = invVar, family = binomial(logit))
fit

# Plot estimates
plot(est ~ day, data = df, pch = 16, cex = 2)
# GLM fit
lines(y = predict(fit, data.frame(day = 1:69), type = "response"), x = 1:69)
# Birch's WLS fit
lines(y = 1 / (1 + exp(-0.1911 * (1:69 - 44.735))), x = 1:69, col = "red")


## How to determine Birch's coefficients for WLS logistic regression
# Brute force approach!!!
# Define a broad scale of coefficient parameters and try to minimize the Sum of Squares
# Broad scale look
kappas <- seq(from = 0.15, to = 0.30, by = 0.001)
gammas <- seq(from = 40, to = 50, by = 0.1)

logistic <- function(kappa, gamma, day) {
  x <- 1 / (1 + exp(-kappa * (day - gamma)))
  return(x)
}

ssq <- matrix(data = NA, nrow = length(kappas), ncol = length(gammas), dimnames = list(kappas, gammas))
for(k in seq(kappas)){
  for(g in seq(gammas)){
    res <- df$est - logistic(kappa = kappas[k], gamma = gammas[g], day = df$day)
    wt.res <- df$invVar * res^2
    ssq[k, g] <- sum(wt.res)
  }
}
persp(ssq, theta = -30, phi = -15)

ind <- which(ssq == min(ssq), arr.ind = TRUE)
kappas[ind[1]]; gammas[ind[2]]

# Zoom in on minimized region to determine fine scale parameters (4 decimals for kappa, 3 for gamma)
# fine scale look
fine.kappas <- seq(from = 0.190, to = 0.192, by = 0.0001)
fine.gammas <- seq(from = 44.6, to = 44.8, by = 0.001)

logistic <- function(kappa, gamma, day) {
  x <- 1 / (1 + exp(-kappa * (day - gamma)))
  return(x)
}

ssq <- matrix(data = NA, nrow = length(fine.kappas), ncol = length(fine.gammas), dimnames = list(fine.kappas, fine.gammas))
for(k in seq(fine.kappas)){
  for(g in seq(fine.gammas)){
    res <- df$est - logistic(kappa = fine.kappas[k], gamma = fine.gammas[g], day = df$day)
    wt.res <- df$invVar * res^2
    ssq[k, g] <- sum(wt.res)
  }
}
persp(ssq, theta = -30, phi = -15)

ind <- which(ssq == min(ssq), arr.ind = TRUE)
fine.kappas[ind[1]]; fine.gammas[ind[2]]



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2015 Round 2 July 1, 2015 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(plot3D)
GeneticEstimates2015 <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/GeneticEstimates2015.txt")
GeneticEstimates2015

df.new <- data.frame("day" = c(1, 69, GeneticEstimates2015[!is.na(GeneticEstimates2015[, "Day"]), "Day"]),
                     "est" = c(0, 1, 1 - GeneticEstimates2015[!is.na(GeneticEstimates2015[, "BlackMean"]), "BlackMean"]),
                     "sd" = c(0.025, 0.025, GeneticEstimates2015[!is.na(GeneticEstimates2015[, "sd"]), "sd"]))
df.new$invVar <- 1 / (df.new$sd^2)
df.new

# Define a broad scale of coefficient parameters and try to minimize the Sum of Squares
kappas <- seq(from = 0.10, to = 0.50, by = 0.001)
gammas <- seq(from = 30, to = 60, by = 0.1)

logistic <- function(kappa, gamma, day) {
  x <- 1 / (1 + exp(-kappa * (day - gamma)))
  return(x)
}

ssq <- matrix(data = NA, nrow = length(kappas), ncol = length(gammas), dimnames = list(kappas, gammas))
for(k in seq(kappas)){
  for(g in seq(gammas)){
    res <- df.new$est - logistic(kappa = kappas[k], gamma = gammas[g], day = df.new$day)
    wt.res <- df.new$invVar * res^2
    ssq[k, g] <- sum(wt.res)
  }
}
persp(ssq, theta = -30, phi = -15)
persp3D(ssq)

ind <- which(ssq == min(ssq), arr.ind = TRUE)
kappas[ind[1]]; gammas[ind[2]]

# Zoom in on minimized region to determine fine scale parameters (4 decimals for kappa, 3 for gamma)
fine.kappas <- seq(from = 0.220, to = 0.240, by = 0.0001)
fine.gammas <- seq(from = 45, to = 50, by = 0.001)

logistic <- function(kappa, gamma, day) {
  x <- 1 / (1 + exp(-kappa * (day - gamma)))
  return(x)
}

ssq <- matrix(data = NA, nrow = length(fine.kappas), ncol = length(fine.gammas), dimnames = list(fine.kappas, fine.gammas))
for(k in seq(fine.kappas)){
  for(g in seq(fine.gammas)){
    res <- df.new$est - logistic(kappa = fine.kappas[k], gamma = fine.gammas[g], day = df.new$day)
    wt.res <- df.new$invVar * res^2
    ssq[k, g] <- sum(wt.res)
  }
}
persp(ssq, theta = -30, phi = -15)

ind <- which(ssq == min(ssq), arr.ind = TRUE)
fine.kappas[ind[1]]; fine.gammas[ind[2]]

# Plot
NewData <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Estimates objects/SCHIG15_2_Jul01_Estimates.txt")
Period = 2
NumSampled = 189
NumAnalyzed = 189
Included = 188
Month = "July"
Day = 01

GeneticEstimates2015 <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/GeneticEstimates2015.txt")
birch.lwr <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/birch.lwr.txt")
birch.upr <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/birch.upr.txt")

par(mar = c(4.2, 8, 3, 8))
par(cex = 1)
par(family = 'serif')
plot(NA, xlim = c(27, 68), ylim = c(0, 100), bty = "n", axes = FALSE, xlab = "Sample date", ylab = "Black Lake (Early Run) % of sample", cex.lab = 1.5)
axis(side = 1, at = seq(from = 27, to = 67, by = 5), labels = NA, pos = 0)
text(x = seq(from = 27, to = 67, by = 5), y = -10, labels = c("20-Jun", "25-Jun", "30-Jun", "5-Jul", "10-Jul", "15-Jul", "20-Jul", "25-Jul", "30-Jul"), srt = 45, xpd = TRUE)
axis(side = 2, pos = 27, at = seq(from = 0, to = 100, by = 25), labels = NA)
text(x = 22, y = seq(from = 0, to = 100, by = 25), labels = seq(from = 0, to = 100, by = 25), srt = 0, xpd = TRUE, cex = 1.3)
polygon(x = c(27:68, 68:27), y = c(birch.upr[-c(1:5)] * 100, rev(birch.lwr[-c(1:5)] * 100)), col = "grey70", border = FALSE)  ## col from grey50 to grey70; THD on 062615
arrows(x0 = GeneticEstimates2015[, "Day"], y0 = GeneticEstimates2015[, "5%"] * 100, x1 = GeneticEstimates2015[, "Day"], y1 = GeneticEstimates2015[, "95%"] * 100, angle = 90, code = 3, lwd = 3, col = "black", length = 0.1)  ## lwd from 3 to 2; length from 0.15 to 0.1; THD on 062615
points(x = GeneticEstimates2015[, "Day"], y = GeneticEstimates2015[, "BlackMean"] * 100, col = "black", cex = 2, pch = 21, bg='red')  ## cex from 2 to 1.5; pch from 16 to 21; THD on 062615
abline(h = 50, lwd = 2)
segments(x0 = 27, y0 = 0, x1 = 68, y1 = 0)
segments(x0 = 27, y0 = 0, x1 = 27, y1 = 100)

legend("topright", legend = c("Average\ntransition\n2010-2014", "2015", "Estimated\ntransition\n2015"), pch = c(22, 21, NA), lwd= c(NA, NA, 3), lty = c(NA, NA, 2), col = c("grey70", "black", "black"), pt.bg=c("grey70", "red"),cex = 1.2, bty = "n", pt.cex = 2)  ## pt.cex from 3 to 2; pch from 16 to 21; THD on 062615
lines(x = 27:68, y = (1 - logistic(kappa = fine.kappas[ind[1]], gamma = fine.gammas[ind[2]], day = 27:68)) * 100, lwd = 3, lty = 2, col = "black")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2015 Round 3 July 5, 2015 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(plot3D)
GeneticEstimates2015 <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/GeneticEstimates2015.txt")
GeneticEstimates2015

df.new <- data.frame("day" = c(1, 69, GeneticEstimates2015[!is.na(GeneticEstimates2015[, "Day"]), "Day"]),
                     "est" = c(0, 1, 1 - GeneticEstimates2015[!is.na(GeneticEstimates2015[, "BlackMean"]), "BlackMean"]),
                     "sd" = c(0.025, 0.025, GeneticEstimates2015[!is.na(GeneticEstimates2015[, "sd"]), "sd"]))
df.new$invVar <- 1 / (df.new$sd^2)
df.new

# Define a broad scale of coefficient parameters and try to minimize the Sum of Squares
kappas <- seq(from = 0.10, to = 0.50, by = 0.001)
gammas <- seq(from = 30, to = 60, by = 0.1)

logistic <- function(kappa, gamma, day) {
  x <- 1 / (1 + exp(-kappa * (day - gamma)))
  return(x)
}

ssq <- matrix(data = NA, nrow = length(kappas), ncol = length(gammas), dimnames = list(kappas, gammas))
for(k in seq(kappas)){
  for(g in seq(gammas)){
    res <- df.new$est - logistic(kappa = kappas[k], gamma = gammas[g], day = df.new$day)
    wt.res <- df.new$invVar * res^2
    ssq[k, g] <- sum(wt.res)
  }
}
persp(ssq, theta = -30, phi = -15)
persp3D(ssq)

ind <- which(ssq == min(ssq), arr.ind = TRUE)
kappas[ind[1]]; gammas[ind[2]]

# Zoom in on minimized region to determine fine scale parameters (4 decimals for kappa, 3 for gamma)
fine.kappas <- seq(from = 0.220, to = 0.240, by = 0.0001)
fine.gammas <- seq(from = 45, to = 50, by = 0.001)

logistic <- function(kappa, gamma, day) {
  x <- 1 / (1 + exp(-kappa * (day - gamma)))
  return(x)
}

ssq <- matrix(data = NA, nrow = length(fine.kappas), ncol = length(fine.gammas), dimnames = list(fine.kappas, fine.gammas))
for(k in seq(fine.kappas)){
  for(g in seq(fine.gammas)){
    res <- df.new$est - logistic(kappa = fine.kappas[k], gamma = fine.gammas[g], day = df.new$day)
    wt.res <- df.new$invVar * res^2
    ssq[k, g] <- sum(wt.res)
  }
}
persp(ssq, theta = -30, phi = -15)

ind <- which(ssq == min(ssq), arr.ind = TRUE)
fine.kappas[ind[1]]; fine.gammas[ind[2]]

# Plot
GeneticEstimates2015 <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/GeneticEstimates2015.txt")
birch.lwr <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/birch.lwr.txt")
birch.upr <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/birch.upr.txt")

plot.new()
par(mar = c(4.2, 8, 3, 8))
par(cex = 1)
par(family = 'serif')
plot(NA, xlim = c(27, 68), ylim = c(0, 100), bty = "n", axes = FALSE, xlab = "Sample date", ylab = "Black Lake (Early Run) % of sample", cex.lab = 1.5)
axis(side = 1, at = seq(from = 27, to = 67, by = 5), labels = NA, pos = 0)
text(x = seq(from = 27, to = 67, by = 5), y = -10, labels = c("20-Jun", "25-Jun", "30-Jun", "5-Jul", "10-Jul", "15-Jul", "20-Jul", "25-Jul", "30-Jul"), srt = 45, xpd = TRUE)
axis(side = 2, pos = 27, at = seq(from = 0, to = 100, by = 25), labels = NA)
text(x = 22, y = seq(from = 0, to = 100, by = 25), labels = seq(from = 0, to = 100, by = 25), srt = 0, xpd = TRUE, cex = 1.3)
polygon(x = c(27:68, 68:27), y = c(birch.upr[-c(1:5)] * 100, rev(birch.lwr[-c(1:5)] * 100)), col = "grey70", border = FALSE)  ## col from grey50 to grey70; THD on 062615
arrows(x0 = GeneticEstimates2015[, "Day"], y0 = GeneticEstimates2015[, "5%"] * 100, x1 = GeneticEstimates2015[, "Day"], y1 = GeneticEstimates2015[, "95%"] * 100, angle = 90, code = 3, lwd = 3, col = "black", length = 0.1)  ## lwd from 3 to 2; length from 0.15 to 0.1; THD on 062615
points(x = GeneticEstimates2015[, "Day"], y = GeneticEstimates2015[, "BlackMean"] * 100, col = "black", cex = 2, pch = 21, bg='red')  ## cex from 2 to 1.5; pch from 16 to 21; THD on 062615
abline(h = 50, lwd = 2)
segments(x0 = 27, y0 = 0, x1 = 68, y1 = 0)
segments(x0 = 27, y0 = 0, x1 = 27, y1 = 100)

legend("topright", legend = c("Average\ntransition\n2010-2014", "2015", "Estimated\ntransition\n2015"), pch = c(22, 21, NA), lwd= c(NA, NA, 3), lty = c(NA, NA, 2), col = c("grey70", "black", "black"), pt.bg=c("grey70", "red"),cex = 1.2, bty = "n", pt.cex = 2)  ## pt.cex from 3 to 2; pch from 16 to 21; THD on 062615
lines(x = 27:68, y = (1 - logistic(kappa = fine.kappas[ind[1]], gamma = fine.gammas[ind[2]], day = 27:68)) * 100, lwd = 3, lty = 2, col = "black")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2015 Round 4 July 12, 2015 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(plot3D)
GeneticEstimates2015 <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/GeneticEstimates2015.txt")
GeneticEstimates2015

df.new <- data.frame("day" = c(1, 69, GeneticEstimates2015[!is.na(GeneticEstimates2015[, "Day"]), "Day"]),
                     "est" = c(0, 1, 1 - GeneticEstimates2015[!is.na(GeneticEstimates2015[, "BlackMean"]), "BlackMean"]),
                     "sd" = c(0.025, 0.025, GeneticEstimates2015[!is.na(GeneticEstimates2015[, "sd"]), "sd"]))
df.new$invVar <- 1 / (df.new$sd^2)
df.new

# Define a broad scale of coefficient parameters and try to minimize the Sum of Squares
kappas <- seq(from = 0.10, to = 0.50, by = 0.001)
gammas <- seq(from = 30, to = 60, by = 0.1)

logistic <- function(kappa, gamma, day) {
  x <- 1 / (1 + exp(-kappa * (day - gamma)))
  return(x)
}

ssq <- matrix(data = NA, nrow = length(kappas), ncol = length(gammas), dimnames = list(kappas, gammas))
for(k in seq(kappas)){
  for(g in seq(gammas)){
    res <- df.new$est - logistic(kappa = kappas[k], gamma = gammas[g], day = df.new$day)
    wt.res <- df.new$invVar * res^2
    ssq[k, g] <- sum(wt.res)
  }
}
persp(ssq, theta = -30, phi = -15)
persp3D(ssq)

ind <- which(ssq == min(ssq), arr.ind = TRUE)
kappas[ind[1]]; gammas[ind[2]]

# Zoom in on minimized region to determine fine scale parameters (4 decimals for kappa, 3 for gamma)
fine.kappas <- seq(from = 0.270, to = 0.280, by = 0.0001)
fine.gammas <- seq(from = 53, to = 57, by = 0.001)

logistic <- function(kappa, gamma, day) {
  x <- 1 / (1 + exp(-kappa * (day - gamma)))
  return(x)
}

ssq <- matrix(data = NA, nrow = length(fine.kappas), ncol = length(fine.gammas), dimnames = list(fine.kappas, fine.gammas))
for(k in seq(fine.kappas)){
  for(g in seq(fine.gammas)){
    res <- df.new$est - logistic(kappa = fine.kappas[k], gamma = fine.gammas[g], day = df.new$day)
    wt.res <- df.new$invVar * res^2
    ssq[k, g] <- sum(wt.res)
  }
}
persp(ssq, theta = -30, phi = -15)

ind <- which(ssq == min(ssq), arr.ind = TRUE)
fine.kappas[ind[1]]; fine.gammas[ind[2]]

# Plot
GeneticEstimates2015 <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/GeneticEstimates2015.txt")
birch.lwr <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/birch.lwr.txt")
birch.upr <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/birch.upr.txt")

plot.new()
par(mar = c(4.2, 8, 3, 8))
par(cex = 1)
par(family = 'serif')
plot(NA, xlim = c(27, 68), ylim = c(0, 100), bty = "n", axes = FALSE, xlab = "Sample date", ylab = "Black Lake (Early Run) % of sample", cex.lab = 1.5)
axis(side = 1, at = seq(from = 27, to = 67, by = 5), labels = NA, pos = 0)
text(x = seq(from = 27, to = 67, by = 5), y = -10, labels = c("20-Jun", "25-Jun", "30-Jun", "5-Jul", "10-Jul", "15-Jul", "20-Jul", "25-Jul", "30-Jul"), srt = 45, xpd = TRUE)
axis(side = 2, pos = 27, at = seq(from = 0, to = 100, by = 25), labels = NA)
text(x = 22, y = seq(from = 0, to = 100, by = 25), labels = seq(from = 0, to = 100, by = 25), srt = 0, xpd = TRUE, cex = 1.3)
polygon(x = c(27:68, 68:27), y = c(birch.upr[-c(1:5)] * 100, rev(birch.lwr[-c(1:5)] * 100)), col = "grey70", border = FALSE)  ## col from grey50 to grey70; THD on 062615
arrows(x0 = GeneticEstimates2015[, "Day"], y0 = GeneticEstimates2015[, "5%"] * 100, x1 = GeneticEstimates2015[, "Day"], y1 = GeneticEstimates2015[, "95%"] * 100, angle = 90, code = 3, lwd = 3, col = "black", length = 0.1)  ## lwd from 3 to 2; length from 0.15 to 0.1; THD on 062615
points(x = GeneticEstimates2015[, "Day"], y = GeneticEstimates2015[, "BlackMean"] * 100, col = "black", cex = 2, pch = 21, bg='red')  ## cex from 2 to 1.5; pch from 16 to 21; THD on 062615
abline(h = 50, lwd = 2)
segments(x0 = 27, y0 = 0, x1 = 68, y1 = 0)
segments(x0 = 27, y0 = 0, x1 = 27, y1 = 100)

legend("topright", legend = c("Average\ntransition\n2010-2014", "2015", "Estimated\ntransition\n2015"), pch = c(22, 21, NA), lwd= c(NA, NA, 3), lty = c(NA, NA, 2), col = c("grey70", "black", "black"), pt.bg=c("grey70", "red"),cex = 1.2, bty = "n", pt.cex = 2)  ## pt.cex from 3 to 2; pch from 16 to 21; THD on 062615
lines(x = 27:68, y = (1 - logistic(kappa = fine.kappas[ind[1]], gamma = fine.gammas[ind[2]], day = 27:68)) * 100, lwd = 3, lty = 2, col = "black")


writeClipboard(as.character(round(logistic(kappa = fine.kappas[ind[1]], gamma = fine.gammas[ind[2]], day = 1:71), 4)))




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 2015 Round 5 July 18, 2015 ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(plot3D)
GeneticEstimates2015 <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/GeneticEstimates2015.txt")
GeneticEstimates2015

df.new <- data.frame("day" = c(1, 69, GeneticEstimates2015[!is.na(GeneticEstimates2015[, "Day"]), "Day"]),
                     "est" = c(0, 1, 1 - GeneticEstimates2015[!is.na(GeneticEstimates2015[, "BlackMean"]), "BlackMean"]),
                     "sd" = c(0.025, 0.025, GeneticEstimates2015[!is.na(GeneticEstimates2015[, "sd"]), "sd"]))
df.new$invVar <- 1 / (df.new$sd^2)
df.new

# Define a broad scale of coefficient parameters and try to minimize the Sum of Squares
kappas <- seq(from = 0.10, to = 0.50, by = 0.001)
gammas <- seq(from = 30, to = 60, by = 0.1)

logistic <- function(kappa, gamma, day) {
  x <- 1 / (1 + exp(-kappa * (day - gamma)))
  return(x)
}

ssq <- matrix(data = NA, nrow = length(kappas), ncol = length(gammas), dimnames = list(kappas, gammas))
for(k in seq(kappas)){
  for(g in seq(gammas)){
    res <- df.new$est - logistic(kappa = kappas[k], gamma = gammas[g], day = df.new$day)
    wt.res <- df.new$invVar * res^2
    ssq[k, g] <- sum(wt.res)
  }
}
persp(ssq, theta = -30, phi = -15)
persp3D(ssq)

ind <- which(ssq == min(ssq), arr.ind = TRUE)
kappas[ind[1]]; gammas[ind[2]]

# Zoom in on minimized region to determine fine scale parameters (4 decimals for kappa, 3 for gamma)
fine.kappas <- seq(from = 0.200, to = 0.210, by = 0.0001)
fine.gammas <- seq(from = 52, to = 56, by = 0.001)

logistic <- function(kappa, gamma, day) {
  x <- 1 / (1 + exp(-kappa * (day - gamma)))
  return(x)
}

ssq <- matrix(data = NA, nrow = length(fine.kappas), ncol = length(fine.gammas), dimnames = list(fine.kappas, fine.gammas))
for(k in seq(fine.kappas)){
  for(g in seq(fine.gammas)){
    res <- df.new$est - logistic(kappa = fine.kappas[k], gamma = fine.gammas[g], day = df.new$day)
    wt.res <- df.new$invVar * res^2
    ssq[k, g] <- sum(wt.res)
  }
}
persp(ssq, theta = -30, phi = -15)

ind <- which(ssq == min(ssq), arr.ind = TRUE)
fine.kappas[ind[1]]; fine.gammas[ind[2]]

# Plot
GeneticEstimates2015 <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/GeneticEstimates2015.txt")
birch.lwr <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/birch.lwr.txt")
birch.upr <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Objects/birch.upr.txt")

plot.new()
par(mar = c(4.2, 8, 3, 8))
par(cex = 1)
par(family = 'serif')
plot(NA, xlim = c(27, 68), ylim = c(0, 100), bty = "n", axes = FALSE, xlab = "Sample date", ylab = "Black Lake (Early Run) % of sample", cex.lab = 1.5)
axis(side = 1, at = seq(from = 27, to = 67, by = 5), labels = NA, pos = 0)
text(x = seq(from = 27, to = 67, by = 5), y = -10, labels = c("20-Jun", "25-Jun", "30-Jun", "5-Jul", "10-Jul", "15-Jul", "20-Jul", "25-Jul", "30-Jul"), srt = 45, xpd = TRUE)
axis(side = 2, pos = 27, at = seq(from = 0, to = 100, by = 25), labels = NA)
text(x = 22, y = seq(from = 0, to = 100, by = 25), labels = seq(from = 0, to = 100, by = 25), srt = 0, xpd = TRUE, cex = 1.3)
polygon(x = c(27:68, 68:27), y = c(birch.upr[-c(1:5)] * 100, rev(birch.lwr[-c(1:5)] * 100)), col = "grey70", border = FALSE)  ## col from grey50 to grey70; THD on 062615
arrows(x0 = GeneticEstimates2015[, "Day"], y0 = GeneticEstimates2015[, "5%"] * 100, x1 = GeneticEstimates2015[, "Day"], y1 = GeneticEstimates2015[, "95%"] * 100, angle = 90, code = 3, lwd = 3, col = "black", length = 0.1)  ## lwd from 3 to 2; length from 0.15 to 0.1; THD on 062615
points(x = GeneticEstimates2015[, "Day"], y = GeneticEstimates2015[, "BlackMean"] * 100, col = "black", cex = 2, pch = 21, bg='red')  ## cex from 2 to 1.5; pch from 16 to 21; THD on 062615
abline(h = 50, lwd = 2)
segments(x0 = 27, y0 = 0, x1 = 68, y1 = 0)
segments(x0 = 27, y0 = 0, x1 = 27, y1 = 100)

legend("topright", legend = c("Average\ntransition\n2010-2014", "2015", "Estimated\ntransition\n2015"), pch = c(22, 21, NA), lwd= c(NA, NA, 3), lty = c(NA, NA, 2), col = c("grey70", "black", "black"), pt.bg=c("grey70", "red"),cex = 1.2, bty = "n", pt.cex = 2)  ## pt.cex from 3 to 2; pch from 16 to 21; THD on 062615
lines(x = 27:68, y = (1 - logistic(kappa = fine.kappas[ind[1]], gamma = fine.gammas[ind[2]], day = 27:68)) * 100, lwd = 3, lty = 2, col = "black")


writeClipboard(as.character(1 - round(logistic(kappa = fine.kappas[ind[1]], gamma = fine.gammas[ind[2]], day = 1:71), 4)))
