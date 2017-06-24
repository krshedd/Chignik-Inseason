## 2017 Chignik Escapement mixture analysis - Inseason!
## This is the set-up script where create file structure, update tables from last year, and generate the inseason pdf function
## Kyle Shedd Tue Jun 20 10:08:55 2017

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Initial Setup ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ls()
rm(list = ls(all = TRUE))
search()
getwd()
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures")

# Create file structure based on 2016
dirs2016 <- list.dirs(path = "2016", recursive = FALSE, full.names = FALSE)
setwd("2017")
sapply(dirs2016[-2], dir.create)
sapply(c("Control", "Mixture", "Output"), function(folder) {dir.create(path = paste(getwd(), "BAYES", folder, sep = "/"))})

# Copy baseline
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures")
file.copy(from = "2016/BAYES/ChignikPops24Loci.bse", to = "2017/BAYES/ChignikPops24Loci.bse")
# dir.create(path = "Contacts and schedule")  # Done

## save.image("2016/2016ChignikInseason.RData")
## load("2016/2016ChignikInseason.RData")

#This sources all of the new GCL functions from our local GitHub repository to this workspace
source("H:/R Source Scripts/Functions.GCL_KS.R")
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Set up Estimates Objects with previous year's data ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures")

Black2010to2016Means <- dget(file="2016/Objects/Black2010to2016Means.txt")
Black2010to2016Uppers <- dget(file="2016/Objects/Black2010to2016Uppers.txt")
Black2010to2016Lowers <- dget(file="2016/Objects/Black2010to2016Lowers.txt")

Black2010to2017Means <- rbind(Black2010to2016Means, "2017" = rep(NA, 6))
dput(x = Black2010to2017Means, file = "2017/Objects/Black2010to2017Means.txt")
Black2010to2017Uppers <- rbind(Black2010to2016Uppers, "2017" = rep(NA, 6))
dput(x = Black2010to2017Uppers, file = "2017/Objects/Black2010to2017Uppers.txt")
Black2010to2017Lowers <- rbind(Black2010to2016Lowers, "2017" = rep(NA, 6))
dput(x = Black2010to2017Lowers, file = "2017/Objects/Black2010to2017Lowers.txt")

# This is a new object for this year
GeneticEstimates2016 <- dget(file = "2016/Objects/GeneticEstimates2016.txt")
GeneticEstimates2017 <- matrix(data = NA, nrow = 6, ncol = 5, dimnames = list(1:6, c("Day", "BlackMean", "5%", "95%", "sd")))
dput(x = GeneticEstimates2017, file = "2017/Objects/GeneticEstimates2017.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Making a formula to put out report NEW ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

require(xlsx)

Genetics.Estimates <- do.call(rbind, 
  lapply(2010:2016, function(yr) {
    yr.dat <- read.xlsx(file = "2016/Chignik inseason summary data.xlsx", sheetName = paste0(yr, "_Estimates"), stringsAsFactors = FALSE, colClasses = c("character", rep("numeric", 12)))
    yr.dat$Year <- rep(as.character(yr), nrow(yr.dat))
    yr.dat
  })
)
str(Genetics.Estimates)

Birch.WLS.RegCoeffs <- data.matrix(read.xlsx(file = "2016/Chignik inseason summary data.xlsx", sheetName = "Birch.WLS.RegCoeffs", stringsAsFactors = FALSE, header = TRUE, colIndex = 2:3))
rownames(Birch.WLS.RegCoeffs) <- 2010:2016
str(Birch.WLS.RegCoeffs)

logistic.kg.f <- function(x, kappa, gamma) {
  p = 1 / (1 + exp(-kappa * (x - gamma)))
  return(1 - p)
}

# Plot all genetic estimates
plot(Genetics.Estimates$BlackMean ~ Genetics.Estimates$Day, pch = 16, bty = "n", xlab = "Day", ylab = "Proportion Black Lake", col = as.numeric(Genetics.Estimates$Year)-2009, cex = 2, xlim = c(27, 68), axes = FALSE)  # xlim = c(22, 68) for full dataset
axis(side = 1, at = seq(from = 27, by = 5, length.out = 9), labels = format(seq(from = as.Date(x = "20-Jun", format = "%d-%b"), by = 5, length.out = 9), "%d-%b"))
axis(side = 2, at = seq(from = 0, to = 1, by = 0.1), labels = paste0(seq(from = 0, to = 100, by = 10), "%"))
legend("bottomleft", legend = unique(Genetics.Estimates$Year), fill = seq(unique(Genetics.Estimates$Year)), bty = "n")
arrows(x0 = Genetics.Estimates$Day, y0 = Genetics.Estimates$Black5CI, x1 = Genetics.Estimates$Day, y1 = Genetics.Estimates$Black95CI, angle = 90, code = 3)
points(Genetics.Estimates$BlackMean ~ Genetics.Estimates$Day, pch = 16, col = as.numeric(Genetics.Estimates$Year)-2009, cex = 2)

# Matrix of all daily stock comps from Birch WLS Logistic Regression Coefficients
model.out <- apply(Birch.WLS.RegCoeffs, 1, function(yr) {
  logistic.kg.f(x = 27:68, kappa = yr["kappa"], gamma = yr["gamma"])  # June 20 (day 27) through July 31 (day 68)
})
rownames(model.out) <- 27:68

# Plot year specific regressions
apply(model.out, 2, function(yr) {lines(x = 27:68, y = yr)})

# Plot daily averages from year specific regressions
lines(x = 27:68, y = rowMeans(model.out), lwd = 5)

# Daily average +/- 1.64 SD (90% PI)
model.intervals <- apply(model.out, 1, function(day) {c("lower" = mean(day) + qnorm(0.05) * sd(day), "upper" = mean(day) + qnorm(0.95) * sd(day))} )
model.intervals["lower", ][model.intervals["lower", ] < 0] <- 0
model.intervals["upper", ][model.intervals["upper", ] > 1] <- 1
polygon(x = c(27:68, 68:27), y = c(model.intervals["lower", ], rev(model.intervals["upper", ])), col = "grey50", border = FALSE)

# Daily average +/- 2 SE (95% CI)
model.intervals <- apply(model.out, 1, function(day) {c("lower" = mean(day) - 2 * sd(day) / sqrt(length(day)), "upper" = mean(day) + 2 * sd(day) / sqrt(length(day)))} )
polygon(x = c(27:68, 68:27), y = c(model.intervals["lower", ], rev(model.intervals["upper", ])), col = "grey80", border = FALSE)

# Daily average +/- 2 fixed SD = 0.056
# polygon(x = c(27:68, 68:27), y = c(rowMeans(model.out) - 2 * 0.056, rev(rowMeans(model.out) + 2 * 0.056)), col = "grey80", border = FALSE)

legend("topright", fill = c("grey80", "grey50"), legend = c("+/- 2SE, 95% CI", "+/- 1.64SD, 90 PI"), bty = 'n')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# E-mailed Birch Foster and he agreed that the 95% CI is what we want
avg.transition.lwr <- model.intervals["lower", ]
avg.transition.upr <- model.intervals["upper", ]

names(avg.transition.lwr) <- names(avg.transition.upr) <- format(seq(from = as.Date("6/20", format = "%m/%d"), 
                                                                     to = as.Date("7/31", format = "%m/%d"), 
                                                                     by = "day"), format = "%m/%d")
dput(x = avg.transition.lwr, file = "2017/Objects/avg.transition.lwr.txt")
dput(x = avg.transition.upr, file = "2017/Objects/avg.transition.upr.txt")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#Test ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
NewData <- dget(file = "2016/Estimates objects/SCHIG16_1_Jun27_Estimates.txt")
Period = 1
NumSampled = 190
NumAnalyzed = 190
Included = 189
Month = "June"
Day = 27
Year = 2016
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures")

ChignikInseasonReport.f <- function(NewData = XX, Period = period, NumSampled = sampled, NumAnalyzed = anlayzed, Included = 190, Month = "June", Day = 27) {
  
  if(length(NewData) == 2) {
    NewData <- NewData$Stats
  }
  
  start.wd <- getwd()
  setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures")
  
  GeneticEstimates2017 <- dget(file = "2017/Objects/GeneticEstimates2017.txt")
  day <- as.numeric(as.Date(paste(Month, Day), format = "%B %d") - as.Date("05/24", format = "%m/%d"))
  GeneticEstimates2017[Period, ] <- c(day, NewData[[1]]["Black Lake", c("mean", "5%", "95%", "sd")]) 
  
  Black2010to2017Means <- dget(file = "2017/Objects/Black2010to2017Means.txt")
  Black2010to2017Uppers <- dget(file = "2017/Objects/Black2010to2017Uppers.txt")
  Black2010to2017Lowers <- dget(file = "2017/Objects/Black2010to2017Lowers.txt")
  
  Year <- rev(rownames(Black2010to2017Means))[1]
  
  ## Adding this set of estimates to summary object
  Black2010to2017Means[Year, Period] <- NewData[[1]]["Black Lake", "mean"]
  Black2010to2017Uppers[Year, Period] <- NewData[[1]]["Black Lake", "95%"]
  Black2010to2017Lowers[Year, Period] <- NewData[[1]]["Black Lake", "5%"]
  
  wrapper <- dget(file = "2013/Objects/wrapper.txt")
  
  layoutmatrix2 <- dget(file = "2013/Objects/layoutmatrix2.txt")
  
  layout(layoutmatrix2, widths = c(0.3, 1, 0.3), heights = c(0.75, 0.25, 1, 0.1), respect = FALSE)
  #layout.show(4)
  
  pdf(paste("2017/Updates/2017ChignikInseason", Period, ".pdf", sep = ""), family = "Times", width = 8.5, height = 11, title = paste(Year, " Chignik Inseason", Period, sep = ""))
  layout(layoutmatrix2, widths = c(0.3, 1, 0.3), heights = c(0.6, 0.3, 1, 0.3), respect = FALSE)
  par(family = "serif")
  
  ## Plot 1:
  par(mar = c(0, 4, 6, 4))
  plot(1:8, type = "n", axes = FALSE, xlab = "", ylab = "")
  text(x = 4.5, y = 8, labels = "Chignik Sockeye Salmon Fishery", adj = c(0.5, 1), font = 2, cex = 2.75)
  text(x = 4.5, y = 6, labels = paste("Chignik Weir Sockeye Salmon Stock Composition Summary #", Period, sep = ""), adj = c(0.5, 1), font = 2, cex = 2)
  
  ## Update these lines
  text(x = 4.5, y = 5, labels = paste(Month, " ", Day, ", ", Year, sep = ""), adj = c(0.5, 1), font = 2, cex = 2)
  text(x = 1, y = 2, labels = wrapper(paste("Genetic stock composition estimates for sockeye salmon from the Chignik Weir for ", Month, " ", Day, ", ", Year, ". A total of ", NumSampled, " fish were sampled and ", NumAnalyzed, " were analyzed (", Included, " had adequate data to include in the analysis).", sep = ""), width = 100), adj = c(0, 0), cex = 1.5)
  
  ## Plot 2:
  par(mar = c(1, 1, 0, 1))
  plot(1:5, type = "n", axes = FALSE, xlab = "", ylab = "")
  arrows(x0 = 0, y0 = 5, x1 = 5, y1 = 5, length = 0)
  text(x = 4.5, y = 4.75, labels = "90%", adj = 0.5, font = 2, cex = 1.5)
  text(x = 3, y = 4.75, labels = "Stock", adj = 0.5, font = 2, cex = 1.5)
  text(x = c(3, 4.5), y = 4, labels = c("Composition", "Confidence Intervals"), adj = 0.5, font = 2, cex = 1.5)
  text(x = 1, y = 3.25, labels = "Reporting Group", adj = 0, font = 2, cex = 1.5)
  text(x = c(3, 4.2, 4.7), y = 3.25, labels = c("Estimate", "Lower", "Upper"), adj = 0.5, font = 2, cex = 1.5)
  arrows(x0 = 0, y0 = 3, x1 = 5, y1 = 3, length = 0)
  text(x = 1, y = c(2.5, 1.75), labels = c("Black (Early Run)", "Chignik (Late Run)"), adj = 0, cex = 1.5)
  
  ## Change the estimate object name here:
  text(x = 3, y = c(2.5, 1.75), labels = sprintf("%.1f%%", NewData[[1]][, "mean"]*100), adj = 0.5, cex = 1.5)
  text(x = 4.2, y = c(2.5, 1.75), labels = sprintf("%.1f%%", NewData[[1]][, "5%"]*100), adj = 0.5, cex = 1.5)
  text(x = 4.7, y = c(2.5, 1.75), labels = sprintf("%.1f%%", NewData[[1]][, "95%"]*100), adj = 0.5, cex = 1.5)
  
  arrows(x0 = 0, y0 = 1.5, x1 = 5, y1 = 1.5, length = 0)
  
  ## Plot 3:
  avg.transition.lwr <- dget(file = "2017/Objects/avg.transition.lwr.txt")
  avg.transition.upr <- dget(file = "2017/Objects/avg.transition.upr.txt")
  
  # The x-axis is Day where 5/25 = 1, thus c(22, 68) is 6/15 - 7/31
  par(mar = c(4.2, 8, 3, 8))
  par(cex = 1)
  par(family = 'serif')
  plot(NA, xlim = c(27, 68), ylim = c(0, 100), bty = "n", axes = FALSE, xlab = "Sample date", ylab = "Black Lake (Early Run) % of sample", cex.lab = 1.5)
  axis(side = 1, at = seq(from = 27, to = 67, by = 5), labels = NA, pos = 0)
  text(x = seq(from = 27, to = 67, by = 5), y = -10, labels = c("20-Jun", "25-Jun", "30-Jun", "5-Jul", "10-Jul", "15-Jul", "20-Jul", "25-Jul", "30-Jul"), srt = 45, xpd = TRUE)
  axis(side = 2, pos = 27, at = seq(from = 0, to = 100, by = 25), labels = NA)
  text(x = 22, y = seq(from = 0, to = 100, by = 25), labels = seq(from = 0, to = 100, by = 25), srt = 0, xpd = TRUE, cex = 1.3)
  polygon(x = c(27:68, 68:27), 
          y = c(avg.transition.upr[which(names(avg.transition.upr) == "06/20"):which(names(avg.transition.upr) == "07/31")] * 100,
                rev(avg.transition.lwr[which(names(avg.transition.upr) == "06/20"):which(names(avg.transition.upr) == "07/31")] * 100)), 
          col = "grey70", border = FALSE)  ## col from grey50 to grey70; THD on 062615
  arrows(x0 = GeneticEstimates2017[, "Day"], y0 = GeneticEstimates2017[, "5%"] * 100, x1 = GeneticEstimates2017[, "Day"], y1 = GeneticEstimates2017[, "95%"] * 100, angle = 90, code = 3, lwd = 3, col = "black", length = 0.1)  ## lwd from 3 to 2; length from 0.15 to 0.1; THD on 062615
  points(x = GeneticEstimates2017[, "Day"], y = GeneticEstimates2017[, "BlackMean"] * 100, col = "black", cex = 2, pch = 21, bg='red')  ## cex from 2 to 1.5; pch from 16 to 21; THD on 062615
  abline(h = 50, lwd = 2)
  segments(x0 = 27, y0 = 0, x1 = 68, y1 = 0)
  segments(x0 = 27, y0 = 0, x1 = 27, y1 = 100)
  legend("topright", legend = c("Average\ntransition\n2010-2016", "2017"), pch = c(22, 21), col = c("grey70", "black"), pt.bg=c("grey70", "red"),cex = 1.2, bty = "n", pt.cex = 2)  ## pt.cex from 3 to 2; pch from 16 to 21; THD on 062615
  
  
  ## Change the date here:
  mtext(wrapper(paste("Genetic Stock Composition Estimates of Black (Early Run) Sockeye Salmon Sampled at the Chignik Weir, ", Month, " ", Day, ", ", Year, ".", sep = ""), width = 65), side = 3, line = 1, cex = 1.25, font = 2, adj = 0.5)
  
  ## Plot 4:
  par(mar = c(0, 1, 1, 1))
  plot(1:3, type = "n", axes = FALSE, xlab = "", ylab = "")
  text(x = 1.1, y = 1.95, labels = wrapper("This project is funded by the Alaska Sustainable Salmon Fund (http://www.akssf.org/Default.aspx?id=3423). Samples were collected at the Chignik weir and analyzed at the Gene Conservation Laboratory by Commercial Fisheries Division staff. These results are in-season estimates and are therefore preliminary in nature. Final quality control will be conducted post-season."
                                           , width = 115), adj = c(0, 0), cex = 0.8)
  text(x = 1.1, y = 1.5, labels = wrapper("Page 1 of 1", width = 115), adj = c(0, 0), cex = 0.8)
  text(x = 2.88, y = 1.5, labels = paste("Reported as of: ", format(Sys.time(), "%B %d, %Y"), sep = ''), adj = c(1, 0), cex = 0.8)
  text(x = 2.88, y = 1.2, labels = paste(format(Sys.time(), "%l:%M%p")), adj = c(1, 0), cex = 0.8)
  
  dev.off()
  
  
  dput(x = GeneticEstimates2017, file = "2017/Objects/GeneticEstimates2017.txt")
  dput(Black2010to2017Means, file = "2017/Objects/Black2010to2017Means.txt")
  dput(Black2010to2017Uppers, file = "2017/Objects/Black2010to2017Uppers.txt")
  dput(Black2010to2017Lowers, file = "2017/Objects/Black2010to2017Lowers.txt")
  
  setwd(start.wd)
}



# Testing the function
ChignikInseasonReport.f(NewData = NewData, Period = 1, NumSampled = 190, NumAnalyzed = 190, Included = 189, Month = "June", Day = 27)

# It works!!!
# Go reset "Lowers" "Means" and "Uppers"
# Save function
dput(ChignikInseasonReport.f, file = "2017/Objects/ChignikInseasonReport.f.txt")
