## 2016 Chignik Escapement mixture analysis - Inseason!
## This is the set-up script where create file structure, update tables from last year, and generate the inseason pdf function
## Kyle Shedd Mon Jun 13 09:19:55 2016
## Modified by Tyler Dann Fri June 26 2015

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Initial Setup ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ls()
rm(list = ls(all = TRUE))
search()
getwd()
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures")

# Create file structure based on 2015
dirs2015 <- list.dirs(path = "2015", recursive = FALSE, full.names = FALSE)
setwd("2016")
sapply(dirs2015, dir.create)
sapply(c("Control", "Mixture", "Output"), function(folder) {dir.create(path = paste(getwd(), "BAYES", folder, sep = "/"))})

# Copy baseline
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures")
file.copy(from = "2015/BAYES/ChignikPops24Loci.bse", to = "2016/BAYES/ChignikPops24Loci.bse")
# dir.create(path = "Contacts and schedule")  # Done

## save.image("2016/2016ChignikInseason.RData")
## load("2016/2016ChignikInseason.RData")

#This sources all of the new GCL functions from our local GitHub repository to this workspace
source("H:/R Source Scripts/Functions.GCL_KS.R")
source("C:/Users/krshedd/Documents/R/Functions.GCL.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set up Estimates Objects with previous year's data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures")

Black2010to2015Means <- dget(file="2015/Objects/Black2010to2015Means.txt")
Black2010to2015Uppers <- dget(file="2015/Objects/Black2010to2015Uppers.txt")
Black2010to2015Lowers <- dget(file="2015/Objects/Black2010to2015Lowers.txt")

Black2010to2016Means <- rbind(Black2010to2015Means, "2016" = rep(NA, 6))
dput(x = Black2010to2016Means, file = "2016/Objects/Black2010to2016Means.txt")
Black2010to2016Uppers <- rbind(Black2010to2015Uppers, "2016" = rep(NA, 6))
dput(x = Black2010to2016Uppers, file = "2016/Objects/Black2010to2016Uppers.txt")
Black2010to2016Lowers <- rbind(Black2010to2015Lowers, "2016" = rep(NA, 6))
dput(x = Black2010to2016Lowers, file = "2016/Objects/Black2010to2016Lowers.txt")

# This is a new object for this year
GeneticEstimates2015 <- dget(file = "2015/Objects/GeneticEstimates2015.txt")
GeneticEstimates2016 <- matrix(data = NA, nrow = 6, ncol = 5, dimnames = list(1:6, c("Day", "BlackMean", "5%", "95%", "sd")))
dput(x = GeneticEstimates2016, file = "2016/Objects/GeneticEstimates2016.txt")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Making a formula to put out report NEW ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

birch.dat <- read.table(file = "2016/averageChigISGmodel2010to2015.txt", header = TRUE)
str(birch.dat)

birch.lwr <- birch.dat$lower[which(birch.dat$date == "15-Jun"):which(birch.dat$date == "31-Jul")]
birch.upr <- birch.dat$upper[which(birch.dat$date == "15-Jun"):which(birch.dat$date == "31-Jul")]
names(birch.lwr) <- format(seq(from = as.Date("6/15", format = "%m/%d"), to = as.Date("7/31", format = "%m/%d"), by = "day"), format = "%m/%d")
names(birch.upr) <- names(birch.lwr)
dput(x = birch.lwr, file = "2016/Objects/birch.lwr.txt")
dput(x = birch.upr, file = "2016/Objects/birch.upr.txt")


#Test ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
NewData <- dget(file = "2016/Estimates objects/SCHIG15_1_Jun27_Estimates.txt")
Period = 1
NumSampled = 190
NumAnalyzed = 190
Included = 190
Month = "June"
Day = 27

NewData <- dget(file = "2015/Estimates objects/SCHIG15_2_Jul01_Estimates.txt")
Period = 2
NumSampled = 189
NumAnalyzed = 189
Included = 188
Month = "July"
Day = 01
Year <- 2015 ## Added year
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures")

ChignikInseasonReport.f <- function(NewData = XX, Period = period, NumSampled = sampled, NumAnalyzed = anlayzed, Included = 190, Month = "June", Day = 28) {
  
  start.wd <- getwd()
  setwd("V:/Analysis/4_Westward/Sockeye/Chignik Inseason 2012-2017/Mixtures")
  
  GeneticEstimates2016 <- dget(file = "2016/Objects/GeneticEstimates2016.txt")
  day <- as.numeric(as.Date(paste(Month, Day), format = "%B %d") - as.Date("05/24", format = "%m/%d"))
  GeneticEstimates2016[Period, ] <- c(day, NewData[[1]][[1]][1, c("mean", "5%", "95%", "sd")]) 
  
  Black2010to2016Means <- dget(file = "2016/Objects/Black2010to2016Means.txt")
  Black2010to2016Uppers <- dget(file = "2016/Objects/Black2010to2016Uppers.txt")
  Black2010to2016Lowers <- dget(file = "2016/Objects/Black2010to2016Lowers.txt")
  
  Year <- rev(rownames(Black2010to2016Means))[1]
  
  ## Adding this set of estimates to summary object
  Black2010to2016Means[Year, Period] <- NewData$Stats[[1]][1, "mean"]
  Black2010to2016Uppers[Year, Period] <- NewData$Stats[[1]][1, "95%"]
  Black2010to2016Lowers[Year, Period] <- NewData$Stats[[1]][1, "5%"]
  
  wrapper <- dget(file = "2013/Objects/wrapper.txt")
  
  layoutmatrix2 <- dget(file = "2013/Objects/layoutmatrix2.txt")
  
  layout(layoutmatrix2, widths = c(0.3, 1, 0.3), heights = c(0.75, 0.25, 1, 0.1), respect = FALSE)
  #layout.show(4)
  
  pdf(paste("2016/Updates/2016ChignikInseason", Period, ".pdf", sep = ""), family = "Times", width = 8.5, height = 11, title = paste(Year, " Chignik Inseason", Period, sep = ""))
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
  text(x = 4.5, y = 4.75, labels = "90%", adj = 0.5, font = 2, cex = 1.5)  ## cex from 1.25 to 1.5; THD on 062615
  text(x = 3, y = 4.75, labels = "Stock", adj = 0.5, font = 2, cex = 1.5)  ## cex from 1.25 to 1.5; THD on 062615
  text(x = c(3, 4.5), y = 4, labels = c("Composition", "Confidence Intervals"), adj = 0.5, font = 2, cex = 1.5)  ## cex from 1.25 to 1.5; THD on 062615
  text(x = 1, y = 3.25, labels = "Reporting Group", adj = 0, font = 2, cex = 1.5)  ## cex from 1.25 to 1.5; THD on 062615
  text(x = c(3, 4.2, 4.7), y = 3.25, labels = c("Estimate", "Lower", "Upper"), adj = 0.5, font = 2, cex = 1.5)  ## cex from 1.25 to 1.5; THD on 062615
  arrows(x0 = 0, y0 = 3, x1 = 5, y1 = 3, length = 0)
  text(x = 1, y = c(2.5, 1.75), labels = c("Black (Early Run)", "Chignik (Late Run)"), adj = 0, cex = 1.5)  ## cex from 1.25 to 1.5; THD on 062615
  
  ## Change the estimate object name here:
  text(x = 3, y = c(2.5, 1.75), labels = sprintf("%.1f%%", NewData$Stats[[1]][, "mean"]*100), adj = 0.5, cex = 1.5)  ## cex from 1.25 to 1.5; THD on 062615
  text(x = 4.2, y = c(2.5, 1.75), labels = sprintf("%.1f%%", NewData$Stats[[1]][, "5%"]*100), adj = 0.5, cex = 1.5)  ## cex from 1.25 to 1.5; THD on 062615
  text(x = 4.7, y = c(2.5, 1.75), labels = sprintf("%.1f%%", NewData$Stats[[1]][, "95%"]*100), adj = 0.5, cex = 1.5)  ## cex from 1.25 to 1.5; THD on 062615
  
  arrows(x0 = 0, y0 = 1.5, x1 = 5, y1 = 1.5, length = 0)
  
  ## Plot 3:
  birch.lwr <- dget(file = "2016/Objects/birch.lwr.txt")
  birch.upr <- dget(file = "2016/Objects/birch.upr.txt")
  
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
          y = c(birch.upr[which(names(birch.upr) == "06/20"):which(names(birch.upr) == "07/31")] * 100,
                rev(birch.lwr[which(names(birch.upr) == "06/20"):which(names(birch.upr) == "07/31")] * 100)), 
          col = "grey70", border = FALSE)  ## col from grey50 to grey70; THD on 062615
  arrows(x0 = GeneticEstimates2016[, "Day"], y0 = GeneticEstimates2016[, "5%"] * 100, x1 = GeneticEstimates2016[, "Day"], y1 = GeneticEstimates2016[, "95%"] * 100, angle = 90, code = 3, lwd = 3, col = "black", length = 0.1)  ## lwd from 3 to 2; length from 0.15 to 0.1; THD on 062615
  points(x = GeneticEstimates2016[, "Day"], y = GeneticEstimates2016[, "BlackMean"] * 100, col = "black", cex = 2, pch = 21, bg='red')  ## cex from 2 to 1.5; pch from 16 to 21; THD on 062615
  abline(h = 50, lwd = 2)
  segments(x0 = 27, y0 = 0, x1 = 68, y1 = 0)
  segments(x0 = 27, y0 = 0, x1 = 27, y1 = 100)
  legend("topright", legend = c("Average\ntransition\n2010-2015", "2016"), pch = c(22, 21), col = c("grey70", "black"), pt.bg=c("grey70", "red"),cex = 1.2, bty = "n", pt.cex = 2)  ## pt.cex from 3 to 2; pch from 16 to 21; THD on 062615
  
  
  ## Change the date here:
  mtext(wrapper(paste("Genetic Stock Composition Estimates of Black (Early Run) Sockeye Salmon Sampled at the Chignik Weir, ", Month, " ", Day, ", ", Year, ".", sep = ""), width = 65), side = 3, line = 1, cex = 1.25, font = 2, adj = 0.5)
  
  ## Plot 4:
  par(mar = c(3, 1, 1, 1))
  plot(1:3, type = "n", axes = FALSE, xlab = "", ylab = "")
  text(x = 1.1, y = 1, labels = wrapper("This project is funded by the Alaska Sustainable Salmon Fund (http://www.akssf.org/Default.aspx?id=3423). Samples were collected at the Chignik weir and analyzed at the Gene Conservation Laboratory by Commercial Fisheries Division staff. These results are in-season estimates and are therefore preliminary in nature. Final quality control will be conducted post-season."
                                        , width = 115), adj = c(0, 0), cex = 0.8)
  dev.off()
  
  
  dput(x = GeneticEstimates2016, file = "2016/Objects/GeneticEstimates2016.txt")
  dput(Black2010to2016Means, file = "2016/Objects/Black2010to2016Means.txt")
  dput(Black2010to2016Uppers, file = "2016/Objects/Black2010to2016Uppers.txt")
  dput(Black2010to2016Lowers, file = "2016/Objects/Black2010to2016Lowers.txt")
  
  setwd(start.wd)
}



# Testing the function
SCHIG15_1_Jun27_Estimates_TEST <- dget(file = "V:/WORK/Sockeye/Chignik/2015 Chignik Inseason/Kyle/Estimates objects/SCHIG14_1_Jun28_Estimates.txt")

ChignikInseasonReport.f(NewData = SCHIG15_1_Jun27_Estimates_TEST, Period = 1, NumSampled = 190, NumAnalyzed = 190, Included = 189, Month = "June", Day = 27)

# It works!!!
# Go reset "Lowers" "Means" and "Uppers"
# Save function
dput(ChignikInseasonReport.f, file = "2016/Objects/ChignikInseasonReport.f.txt")
