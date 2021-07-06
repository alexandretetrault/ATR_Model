#Load Packages
library(baseline)
library(signal)
library(pls)
library(caTools)
library(dplyr)
library(Hmisc)

#### Input Code for Compounds ####

#Load ATR Spectral Data File Names
file_names <- list.files(path = "Spectra", pattern = "*.asp")

#Create String of Compounds
Compounds <- print(paste(gsub("\\(.*", "", file_names)))

#Compile Compounds into List
spectra <- list()

for (i in file_names){
  spectra[[i]] <- list()
}

spectra <- lapply(spectra, function(x){
  x$Wavenumber <- read.csv("Wavenumbers.csv", 
                           sep = ",", 
                           col.names = "Wavenumber")
  }
)

##[UNDECIDED]Round Wavenumber Values##
#spectra <- lapply(spectra, function(x){
#  x$Wavenumber <- round(x$Wavenumber)
#  return(x)
#}
#)

#Add ATR Data
spectra <- sapply(file_names,
                  USE.NAMES = TRUE,
                  simplify = FALSE,
                  function(x){
  cbind(spectra[[x]], read.csv(paste0("./Spectra/", x),
                               skip = 6,
                               sep = ",",
                               col.names = "Transmittance"
  )
  )
}
)

#Transform Transmittance to Absorbance
spectra <- lapply(spectra, function(x){
  x$Absorbance <- as.vector(
    2-log10(x$Transmittance)
  )
  return(x)
}
)

#Baseline Correct Absorbance 
spectra <- lapply(spectra, function(x){ 
  x$Absorbance <- t(
    baseline.modpolyfit(
      t(x$Absorbance),
      degree = 4
    )[[2]]
  )
  return(x)
}
)

#Min-Max Normalize Absorbance
spectra <- lapply(spectra, function(x){
  x$Normalized <- as.vector(
    scale(x$Absorbance,
          center = FALSE,
          scale = max(x$Absorbance) -
            min(x$Absorbance)
    )
  )
  return(x)
}
)
####END####

####Average the Normalized Absorbance Spectra####
spectra_avg <- matrix(nrow = length(spectra[[1]]$Wavenumber))

file_names <- NULL 

for(i in seq(0, length(spectra)-5, 5)){
  avg <- as.vector(
    (spectra[[i+1]][4] + 
       spectra[[i+2]][4] + 
       spectra[[i+3]][4] +
       spectra[[i+4]][4] +
       spectra[[i+5]][4])/5)
  spectra_avg <- cbind(spectra_avg, avg)
  file_names <- c(file_names, Compounds[i +1])
  
}
spectra_avg <- data.frame(spectra_avg[ , -1])
colnames(spectra_avg) <- file_names
###END###

####Load Isotherm Results####
isotherm.raw <- lapply(paste0("./Isotherms/", file_names, ".csv"), read.delim)

#Assign Names
names(isotherm.raw) <- file_names

#Perform Linear Regressions
isotherm.fit <- lapply(isotherm.raw, function(i) lm(i[[3]] ~ i[[2]],data = i))

#Extract Kd Values from Slopes
Kd <- data.frame(file_names, row.names = 1, Kd = NA)

for(i in file_names){
  Kd[i,1] <- round((coef(isotherm.fit[[i]])[2]), 1) 
}
####END####

####Add Compounds + Kd values to Target Variables.csv if necessary####
#Write the file names of compounds to csv
write.csv(file_names, file = "File Names.csv")

#Write Kd Values to csv
write.csv(Kd, file = "Kd.csv")

## -> Manually copy file names to Target Variables.csv and add in Kd data 

#Load Target Variables Data
Targets <- read.csv("Target Variables.csv", sep = ",", row.names = 1)
###END###


####[DEPRECATED]Plot ATR Spectra####

#Create a Reversed Wavenumber Scale for Conventional Plotting
wavenumber <- rev(range(spectra[[1]][1]))

lapply(names(spectra), function(x){
  #pdf(gsub("_.*", "", x), width = 10, height =5)
  plot(spectra[[x]][c(1,2)], 
       type = "l",
       main = gsub("_.*", "", x), 
       xlim = wavenumber
  )
  minor.tick(nx = 5)
  #dev.off()
}
)
###END###

#####[DEPRECATED]Plot Normalized Absorbance ATR Spectra####
lapply(names(spectra), function(x){
  tmp <- data.frame(spectra[[x]][1], spectra[[x]][4] )
  plot(tmp,
       type = "l",
       main = print(paste(gsub("_.*", "", x), "Normalized Absorbance")), 
       xlim = wavenumber
  )
  minor.tick(nx = 5)
}
)
###END###

####[DEPRECATED]Export Spectra for Peak Fitting####
write.csv(cbind(wavenumbers, spectra_avg$`Mar30_TPC_0-0-1_HMW`),
          row.names = FALSE,
          file = "/home/alex/Desktop/CNOM_HMW.csv"
)

write.csv(cbind(wavenumbers, spectra_avg$`Apr6_TPC_1-0-0_HMW`),
          row.names = FALSE,
          file = "/home/alex/Desktop/TNOM_HMW.csv"
)

write.csv(cbind(wavenumbers, spectra_avg$`Mar22_TPC_0-1-0_HMW`),
          row.names = FALSE,
          file = "/home/alex/Desktop/PNOM_HMW.csv"
)

###END###

#####[DEPRECATED]Calculate First Derivative with Savitzky-Golay Filter & Min-Max Normalize####
spectra <- lapply(spectra, function(x){
  deriv <- sgolayfilt(x$Absorbance,
                      p = 3,
                      n = 5,
                      m = 1
  )
  deriv <- as.vector(
    scale(
      deriv,
      center = FALSE,
      scale = max(deriv) - min(deriv)
    )
  )
  x$Deriv <- deriv
  return(x)
}
)

#Plot First Derivative ATR Spectra
lapply(names(spectra), function(x){
  tmp <- data.frame(spectra[[x]][1], spectra[[x]][5] )
  plot(tmp,
       type = "l",
       main = print(paste(gsub("_.*", "", x), "1st Derivative Absorbance")), 
       xlim = wavenumber
  )
  minor.tick(nx = 5)
}
)
####End####

####Dataframe Preparation for PLS - Full Spectrum####
spectra.df <- data.frame(row.names = names(spectra_avg))

spectra.df$Kd <- Targets$Kd

#Log-transform KD Values
spectra.df$logKd <- log10(Targets$Kd)
####END####

####Make Nested Matrix of MIR Spectra for All Compounds####
MIR <- as.matrix(spectra_avg)
MIR <- t(MIR)
spectra.df$MIR <- MIR
####END####

####Plot the ATR Spectra####
wavenumber <- rev(range(spectra[[1]][1]))
wavenumbers <- as.vector(spectra[[1]][1])

for (i in 1:ncol(spectra_avg)){
  plot(wavenumbers[, 1],
       spectra_avg[ , i],
       type = "l",
       xlim = wavenumber,
       xlab = "Wavenumber",
       ylab = "Absorbance",
       main = names(spectra_avg[i])
  )
}

col_set <- ifelse(spectra.df[,1] > 500, "red",
                  ifelse(between(spectra.df[,1], 300, 500), "yellow", "blue"))

matplot(wavenumbers,
        spectra_avg,
        type = "l",
        xlim = wavenumber,
        col = col_set,
        xlab = "Wavenumber (1/cm)",
        ylab = "Absorbance",
        main = "Overlaid ATR-FTIR Spectra"
)
minor.tick(nx = 5)
legend("top", 
       legend = c("Kd > 600", "Kd > 300", "Kd < 300"), 
       fill = c("red", "yellow", "blue")
       )
legend("top",legend = names(spectra_avg), fill = col_set)
###END###

####Calculate Second Derivative for Peak Assignment####
#Note: n = 7 was found to best capture peaks in absorbance, with find_peaks tolerance (below)
# set to m = 5.  However this might be overfitting so we will use n = 11:

MIR.deriv <- apply(spectra.df$MIR, 
                   MARGIN = 1, 
                   FUN = sgolayfilt, 
                   p = 3, n = 11, m = 2, ts = 1
                   )

spectra.df$MIR.deriv <- t(MIR.deriv)

#Peak-finding Alogrithm from https://stats.stackexchange.com/questions/22974/how-to-find-local-peaks-valleys-in-a-series-of-data
find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

#Create a List of All Peaks Calculated for 2nd Derivative Spectra
peaks.list <- apply(
  spectra.df$MIR.deriv, 
  MARGIN = 1, 
  function(i)wavenumbers[find_peaks(-i, m = 4),]
  ) 

#Retain Only Region Below 1800 cm-1
peaks.list <- lapply(peaks.list, function(i)i[which(i < 1800)])

#Extract Peaks for One Compound
#TNOM_HMW_peaks <- peaks.list$`Apr6_TPC_1-0-0_HMW`
#CNOM_HMW_peaks <- peaks.list$`Mar30_TPC_0-0-1_HMW`
#PNOM_HMW_peaks <- peaks.list$`Mar22_TPC_0-1-0_HMW`

#write.csv(c(TNOM_HMW_peaks, PNOM_HMW_peaks, CNOM_HMW_peaks), "peaks.csv")

#Plot 2nd Derivative Spectra of All Compounds
for (i in colnames(spectra_avg)){
  plot(wavenumbers[,1],
       spectra.df$MIR.deriv[i,],
       type = "l",
       xlim = c(1800, 650),
       main = i
  )
  minor.tick(nx = 4)
  abline(v = peaks.list[[i]], col = "red")
}

#Plot Found Peaks Against Absorbance Spectrum for All Compounds
for (i in colnames(spectra_avg)){
     plot(wavenumbers[,1],
          spectra_avg[[i]],
          type = "l",
          xlim = c(1800, 650),
          main = i
          )
  minor.tick(nx = 4)
  abline(v = peaks.list[[i]], col = "red")
}
####END####

####Remove CO2 Region & Split Data for Training on Spectra####
#Make the Wavenumbers Scale
wavenumbers <- as.vector(spectra[[1]][1])

#Spectral Regions of Interest
cm3600_2200 <- which(between(wavenumbers[[1]], 2200, 3600) == TRUE)
cm1900_650 <- which(between(wavenumbers[[1]], 650, 1900) == TRUE)
wavenumbers_NoCO2 <- wavenumbers[c(cm3600_2200, cm1900_650), 1]

spectra.df$MIR <- spectra.df$MIR[, c(cm3600_2200, cm1900_650)]
set.seed(123)
split <- sample.split(spectra.df$Kd, SplitRatio = 0.8)

train <- spectra.df[-(45:50),]
train <- subset(spectra.df, split == T)

test <- spectra.df[45:50,]
test <- subset(spectra.df, split == F)
####END####


####Run Model on Full Spectrum####
plsr.fit <- plsr(logKd ~ MIR,
                     ncomp = 4,
                     data = train,
                     validation = "CV",
                     scale = FALSE,
                     center = TRUE,
                     method = "simpls" 
)

#Change MIR Dimension Names for Easier Plotting
dimnames(plsr.fit[[3]])[[1]] <- wavenumbers_NoCO2

summary(plsr.fit)

plot(RMSEP(plsr.fit))
####END####

####Plot the Scores####
plot(plsr.fit, plottype = "scores", comps = 1:3, labels = "numbers")

##If only one component selected, plot u vs. t
plot(plsr.fit[[2]][,1], 
     plsr.fit[[4]][,1],
     xlab = "t",
     ylab = "u"
)
abline(a = 0, b = 1)

text(x = plsr.fit[[2]][,1], 
     y = plsr.fit[[4]][,1], 
     labels = paste(gsub("_.*","", names(plsr.fit[[2]][,1]))), 
     pos = 3
)
####END####

####Plot the Model-Fitted Coefficients####
plot(plsr.fit[[1]],
     type = "l")

##Plot Coefficients at Specific no. Components
#Choose Number of Components
Num.components <- 2

plot(wavenumbers_NoCO2,
     plsr.fit[[1]][(1422*Num.components):((1422*Num.components)+1421)],
     type = "l",
     xlim = c(3600,650),
     xlab = "Wavenumber",
     ylab = "Coefficient",
     main = paste(c("Regression Coefficients by Wavenumber", Num.components, "Components"))
     )
abline(h = 0)
minor.tick(nx = 5)

#Store Coefficients as Object
Fit.Coefficients <-plsr.fit[[1]][(1422*Num.components):((1422*Num.components)+1421)] 

#Create a List of All Peaks Calculated for Coefficients Below 1800 cm-1
#peaks.list <- wavenumbers_NoCO2[find_peaks(Fit.Coefficients, m = 3)]
Fit.deriv <- sgolayfilt(Fit.Coefficients, 
                        p = 3, 
                        n = 19, 
                        m = 2, 
                        ts = 1
) 

peaks.list <- wavenumbers_NoCO2[c(find_peaks(Fit.deriv, m = 3),find_peaks(-Fit.deriv, m = 3))]


#Retain Only Region Below 1800 cm-1
peaks.list <- peaks.list[which(peaks.list < 1800)]

###Not quite..
for (i in 1:length(Fit.deriv)){
  if(abs(Fit.deriv[i]) > 0.00005) Fit.deriv[i] <- Fit.deriv[i]
  else Fit.deriv[i] <- NA
}

#Plot Derivative
plot(wavenumbers_NoCO2,
     Fit.deriv,
     type = "l",
     xlim = c(1800, 650)
)

#Plot Coefficients
plot(wavenumbers_NoCO2,
     Fit.Coefficients,
     type = "l",
     xlim = c(1800, 650)
     )

abline(v = peaks.list, col = "red")

#plot(wavenumbers_NoCO2,
#     Fit.Coefficients,
#     type = "l",
#     xlim = c(1800,650),
#     xlab = "Wavenumber",
#     ylab = "Coefficient",
#     main = paste(c("Regression Coefficients by Wavenumber", Num.components, "Components"))
#)
#abline(v = peaks.list, col = "red")


##Working on this: Not functional yet, trying to limit peak find to "real" peaks and troughs (June 21)
temp <- find_peaks((abs(Fit.Coefficients)>0.005), m = 4)
abline(v = temp, col = "red")
####END####

####Test to see if coefficients are indeed the regression coefficients Y onto X####
#Account for mean centering by calculating the intercept to add
coef1 <- coef(plsr.fit, intercept = TRUE, ncomp = 1)[[1]]

(test$MIR)%*%plsr.fit[[1]][1:1422]+coef1 

predict(plsr.fit, ncomp = 1, newdata = test)
####END####

#####Plot Fit for Training Data####
plot(plsr.fit, ncomp = 2, line = TRUE,
     xlim = c(1.5, 3),
     ylim = c(1.5, 3))#,
labels = "numbers")
#axes = FALSE
#)
#axis(1, pos=0)
#axis(2, pos=0)

#Add labels
text(x = train[,2], 
     y = plsr.fit[[9]][36:70], 
     labels = paste(gsub("_.*","", names(plsr.fit[[2]][,2])))#, 
     #pos = 3
)
RMSEP(plsr.fit)
####END####

####Fit of Test Data####
plot(plsr.fit, ncomp = 2, line = TRUE, newdata = test,
     xlim = c(1.5,3),
     ylim = c(1.5,3)
)

RMSEP(plsr.fit, ncomp = 2, newdata = test)
####END####

####Compare Predicted vs. Actual Kd Values of Test Set as log-Untransformed####
predicted_values <- 10^predict(plsr.fit, ncomp = 2, newdata = test)

plot(test[,1],
     predicted_values,
     type = "p"
     )
abline( a = 0, b = 1)
####END####

#Obtain Predicted KD values of Training Data
predict(plsr.fit, ncomp = 2)

#Predict KD values of Test Data
predict(plsr.fit, ncomp = 2, newdata = test)

#Plot Loadings of First Three Components
plot(plsr.fit, 
     plottype = "loadings", 
     comps = 1:3, 
     legendpos = "topleft",
     labels = "numbers", 
     xlab = "wavenumber"
)
abline(h = 0)
minor.tick(nx = 5)
###END###

####Prediction of Wild Sediments####
plot(plsr.fit, ncomp = 2, line = TRUE, newdata = sediment,
     xlim = c(1.5,3),
     ylim = c(1.5,3)
)

RMSEP(plsr.fit, ncomp = 2, newdata = sediment)


####Build Dataframe of Test Prediction Errors for Later Use####
#Homemade RMSEP Function
alex.RMSEP <- function(model, data, target){
  temp <- data.frame(predict(model, newdata = data))
  temp <- temp - data[[target]]
  temp[,] <- temp[,]^2
  temp <- colSums(temp)
  temp <- temp/length(data[[target]])
  temp <- sqrt(temp)
}

RMSEP.table <- data.frame(
  alex.RMSEP(model = plsr.fit, data = test, target = "logKd")
)

colnames(RMSEP.table) <- "Full Spectrum"
####END####

####Compare homebrew RMSEP to Package RMSEP function on Full Spectrum Model####
RMSEP.table$Full.Stock <- t(data.frame(
  RMSEP(plsr.fit, estimate = "test", newdata = test, intercept = FALSE)$val)
)
####END####

####Moving Window Partial Least Squares####
##Set Total Number of Spectral Points (i.e w/o CO2)
total.points <- length(spectra.df$MIR[1,])

#Create an object to store the results
residues <- NULL

for(i in 1:(total.points - 4)){
  #Set the midpoint of the moving window of 5 spectral points
  midpoint <- wavenumbers_NoCO2[i+2]
  
  #Train the Model on Scaled Absorbance Data
  mwindow.fit <- plsr(logKd ~ MIR[ ,i:(i+4)],
                      ncomp = 3, 
                      data = train, 
                      validation = "none", 
                      scale = FALSE, 
                      center = TRUE,
                      method = "simpls"
  )
  #Extract the Residual Sum of Squares for the Fit (RSS)
  RSS <- mvrValstats(mwindow.fit, estimate = "train", intercept = "TRUE")
  RSS <- RSS$SSE

  #Append Data to Dataframe
  residues <- rbind(
    residues, 
    data.frame(midpoint, RSS))
  
}

#Plot Residue Lines by Midpoint of Spectral Window
matplot(residues$midpoint, 
        log10(residues[,-1]), 
        type = "l", 
        xlim = c(3600,650), 
        lty = 1,
        ylab = "log(RSS)",
        xlab = "Wavenumber"
)
minor.tick(nx = 5)
abline(h = 0)
title(main = "Traces des Residus, Taille de Fenêtre = 5")
####END####

####Define Windows####
#Window 1
#Define Spectral Region 1
window.1 <- which(between(wavenumbers_NoCO2, 1550, 1800) == TRUE)

#Fit PLS
#Train the Model on Windows of Absorbance Data
window1.fit <- plsr(logKd ~ MIR[ , window.1],
                        ncomp = 4, 
                        data = train, 
                        validation = "CV", 
                        scale = FALSE, 
                        center = TRUE,
                        method = "simpls"
)

#Calculate RMSEP Window 1
plot(RMSEP(window1.fit))
RMSEP(window1.fit)

#Window 2
#Define Spectral Region 2
window.2 <- which(between(wavenumbers_NoCO2, 1150, 1450) == TRUE)

#Fit PLS
#Train the Model on Windows of Absorbance Data
window2.fit <- plsr(logKd ~ MIR[ , window.2],
                    ncomp = 4, 
                    data = train, 
                    validation = "CV", 
                    scale = FALSE, 
                    center = TRUE,
                    method = "simpls"
)

#Calculate RMSEP Window 2
plot(RMSEP(window2.fit))
RMSEP(window2.fit)

#Window 3
#Define Spectral Region 2
window.3 <- which(between(wavenumbers_NoCO2, 800, 1100) == TRUE)

#Fit PLS
#Train the Model on Windows of Absorbance Data
window3.fit <- plsr(logKd ~ MIR[ , window.3],
                    ncomp = 4, 
                    data = train, 
                    validation = "CV", 
                    scale = FALSE, 
                    center = TRUE,
                    method = "simpls"
)

#Calculate RMSEP Window 3
plot(RMSEP(window3.fit))
RMSEP(window3.fit)

#Window 4
#Define Spectral Region 4
window.4 <- which(between(wavenumbers_NoCO2, 2800, 3000) == TRUE)

#Fit PLS
#Train the Model on Windows of Absorbance Data
window4.fit <- plsr(logKd ~ MIR[ , window.4],
                    ncomp = 4, 
                    data = train, 
                    validation = "CV", 
                    scale = FALSE, 
                    center = TRUE,
                    method = "simpls"
)

#Calculate RMSEP Window 4
plot(RMSEP(window4.fit))
RMSEP(window4.fit)

#Build Model on all windows combined
combined.fit <- plsr(logKd ~ 
                       MIR[ , window.1] + 
                       MIR[ , window.2] +
                       MIR[ , window.3] +
                       MIR[ , window.4],
                     ncomp = 4,
                     data = train, 
                     validation = "CV",
                     scale = FALSE, 
                     center = TRUE,
                     method = "simpls"
)
####END####

####RMSEC Data Comparison for Optimal Window Selection####
#Combine Calibration Prediction Data of all Models for Comparison
RMSEC.table <- t(data.frame(RMSEP(plsr.fit, estimate = "CV", intercept = FALSE)$val))

RMSEC.table <- data.frame(RMSEC.table)

colnames(RMSEC.table)[colnames(RMSEC.table) == "CV"] <- "Full Spectrum"

RMSEC.table$MovingWindow <- t(
  data.frame(RMSEP(combined.fit, estimate = "CV", intercept = FALSE)$val)
)

RMSEC.table$Window1 <- t(
  data.frame(RMSEP(window1.fit, estimate = "CV", intercept = FALSE)$val)
)

RMSEC.table$Window2 <- t(
  data.frame(RMSEP(window2.fit, estimate = "CV", intercept = FALSE)$val)
)

RMSEC.table$Window3 <- t(
  data.frame(RMSEP(window3.fit, estimate = "CV", intercept = FALSE)$val)
)

RMSEC.table$Window4 <- t(
  data.frame(RMSEP(window4.fit, estimate = "CV", intercept = FALSE)$val)
)
####END####

####CSMWPLS on Base Region####
##Based on RMSEC curve of Best Window, ideal no. components should be __
Num.components <- 2

##CSMWPLS on Base Region (Window _)

#Select window
window <- window.3

#Clear container
CSMWPLS <- NULL

set.seed(123)
#Set size of moving window
for(w in 4:length(window)){
  
  #Factor region by window size to determine no. windows that span it
  span <- length(window) - w + 1
  
  for(j in 1:span){
    #Train the Model on Scaled Absorbance Data
    base1.fit <- plsr(logKd ~ 
                        MIR[ , window[j]:window[j + w - 1]],
                      ncomp = Num.components, 
                      data = train, 
                      validation = "CV", 
                      scale = FALSE, 
                      center = TRUE,
                      method = "simpls"
    )
    
    #Extract the Residual Sum of Squares for the Fit (RSS)
    RSS <- RMSEP(base1.fit, estimate = "CV", intercept = "TRUE")
    RSS <- RSS$val
    
    #Extract the minimum RMSEP and its associated number of LV's
    LV <- which.min(RSS)
    
    .RMSEP <- min(RSS)
    
    temp.matrix <- rbind(.RMSEP, LV)
    
    colnames(temp.matrix) <- print(
      paste(
        wavenumbers_NoCO2[window[j]],
        "-",
        wavenumbers_NoCO2[window[j + w - 1]]
      )
    ) 
    
    #Append Data to Dataframe
    CSMWPLS <- cbind(CSMWPLS, temp.matrix )
    
    rm(temp.matrix)
  }
}

#Query for Spectral Window of lowest RMSEP 
which.min(CSMWPLS[1,])                              

CSMWPLS[, which.min(CSMWPLS[1,])]

Base.Region <- which(between(wavenumbers_NoCO2, 1039, 1048) == TRUE)

set.seed(123)
#Run PLSR on Base Region and Append for General Comparison
base1.fit <- plsr(logKd ~ MIR[ , Base.Region],
                        ncomp = 4, 
                        data = train, 
                        validation = "CV", 
                        scale = FALSE, 
                        center = TRUE,
                        method = "simpls"
)

RMSEC.table$Base.Region <- t(
  data.frame(RMSEP(base1.fit, estimate = "CV", intercept = FALSE)$val)
)
####END####

####SCSMWPLS on Base Region + Window_####

#Select next window
window <- window.2 

#Clear container
CSMWPLS <- NULL

set.seed(123)
#Set size of moving window
for(w in 4:length(window)){
  
  #Factor region by window size to determine no. windows that span it
  span <- length(window) - w + 1
  
  for(j in 1:span){
    #Train the Model on Scaled Absorbance Data
    base12.fit <- plsr(logKd ~ 
                         MIR[ , Base.Region] +
                         MIR[ , window[j]:window[j + w - 1]],
                       ncomp = Num.components, 
                       data = train, 
                       validation = "CV", 
                       scale = FALSE, 
                       center = TRUE,
                       method = "simpls"
    )
    
    #Extract the Residual Sum of Squares for the Fit (RSS)
    RSS <- RMSEP(base12.fit, estimate = "CV", intercept = "TRUE")
    RSS <- RSS$val
    
    #Extract the minimum RMSEP and its associated number of LV's
    LV <- which.min(RSS)
    
    .RMSEP <- min(RSS)
    
    temp.matrix <- rbind(.RMSEP, LV)
    
    colnames(temp.matrix) <- print(
      paste(
        wavenumbers_NoCO2[window[j]],
        "-",
        wavenumbers_NoCO2[window[j + w - 1]]
      )
    ) 
    
    #Append Data to Dataframe
    CSMWPLS <- cbind(CSMWPLS, temp.matrix )
    
    rm(temp.matrix)
  }
}

#Query for Spectral Window of lowest RMSEP 
which.min(CSMWPLS[1,])                              

CSMWPLS[, which.min(CSMWPLS[1,])]

Base.Region2 <- which(between(wavenumbers_NoCO2, 1218, 1354) == TRUE)

#Run PLSR on Base Region 1 & 2 and Append for General Comparison
base12.fit <- plsr(logKd ~ MIR[ , Base.Region] + MIR[ , Base.Region2],
                           ncomp = 4, 
                           data = train, 
                           validation = "CV", 
                           scale = FALSE, 
                           center = TRUE,
                           method = "simpls"
)

#Strip the "MIR[]" off of the returned dimension names for easier plotting later
dimnames(Base.Region1.2.fit[[3]])[[1]] <- gsub(".*]", "", dimnames(Base.Region1.2.fit[[3]])[[1]])

RMSEC.table$Base.Region12 <- t(
  data.frame(RMSEP(base12.fit, estimate = "CV", intercept = FALSE)$val)
)
####END####

####SCSMWPLS on Base Region12 + Window _####
#Select window
window <- window.1

#Clear container
CSMWPLS <- NULL

set.seed(123)
#Set size of moving window
for(w in 4:length(window)){
  
  #Factor region by window size to determine no. windows that span it
  span <- length(window) - w + 1
  
  for(j in 1:span){
    #Train the Model on Scaled Absorbance Data
    base123.fit <- plsr(logKd ~ 
                          MIR[ , Base.Region] +
                          MIR[ , Base.Region2] +
                          MIR[ , window[j]:
                                 window[j + w - 1]],
                        ncomp = Num.components, 
                        data = train, 
                        validation = "CV", 
                        scale = FALSE, 
                        center = TRUE,
                        method = "simpls"
    )
    
    #Extract the Residual Sum of Squares for the Fit (RSS)
    RSS <- RMSEP(base123.fit, estimate = "CV", intercept = "TRUE")
    RSS <- RSS$val
    
    #Extract the minimum RMSEP and its associated number of LV's
    LV <- which.min(RSS)
    
    .RMSEP <- min(RSS)
    
    temp.matrix <- rbind(.RMSEP, LV)
    
    colnames(temp.matrix) <- print(
      paste(
        wavenumbers_NoCO2[window[j]],
        "-",
        wavenumbers_NoCO2[window[j + w - 1]]
      )
    ) 
    
    #Append Data to Dataframe
    CSMWPLS <- cbind(CSMWPLS, temp.matrix )
    
    rm(temp.matrix)
  }
}

#Query for Spectral Window of lowest RMSEP 
which.min(CSMWPLS[1,])                              

CSMWPLS[, which.min(CSMWPLS[1,])]

Base.Region3 <- which(between(wavenumbers_NoCO2, 1768, 1775) == TRUE)

#Run PLSR on Base Region 123 and Append for General Comparison
base123.fit <- plsr(logKd ~ 
                      MIR[ , Base.Region] + 
                      MIR[ , Base.Region2] +
                      MIR[ , Base.Region3],
                    ncomp = 4, 
                    data = train, 
                    validation = "CV", 
                    scale = FALSE, 
                    center = TRUE,
                    method = "simpls"
)

#Strip the "MIR[]" off of the returned dimension names for easier plotting later
dimnames(Base.Region1.2.fit[[3]])[[1]] <- gsub(".*]", "", dimnames(Base.Region1.2.fit[[3]])[[1]])

RMSEC.table$Base.Region123 <- t(
  data.frame(RMSEP(base123.fit, estimate = "CV", intercept = FALSE)$val)
)
####END####

####SCSMWPLS on Base Region123 + Window _####
#Select window
window <- window.4

#Clear container
CSMWPLS <- NULL

set.seed(123)
#Set size of moving window
for(w in 4:length(window)){
  
  #Factor region by window size to determine no. windows that span it
  span <- length(window) - w + 1
  
  for(j in 1:span){
    #Train the Model on Scaled Absorbance Data
    base1234.fit <- plsr(logKd ~ 
                          MIR[ , Base.Region] +
                          MIR[ , Base.Region2] +
                          MIR[ , Base.Region3] +
                          MIR[ , window[j]:
                                 window[j + w - 1]],
                        ncomp = Num.components, 
                        data = train, 
                        validation = "CV", 
                        scale = FALSE, 
                        center = TRUE,
                        method = "simpls"
    )
    
    #Extract the Residual Sum of Squares for the Fit (RSS)
    RSS <- RMSEP(base123.fit, estimate = "CV", intercept = "TRUE")
    RSS <- RSS$val
    
    #Extract the minimum RMSEP and its associated number of LV's
    LV <- which.min(RSS)
    
    .RMSEP <- min(RSS)
    
    temp.matrix <- rbind(.RMSEP, LV)
    
    colnames(temp.matrix) <- print(
      paste(
        wavenumbers_NoCO2[window[j]],
        "-",
        wavenumbers_NoCO2[window[j + w - 1]]
      )
    ) 
    
    #Append Data to Dataframe
    CSMWPLS <- cbind(CSMWPLS, temp.matrix )
    
    rm(temp.matrix)
  }
}

#Query for Spectral Window of lowest RMSEP 
which.min(CSMWPLS[1,])                              

CSMWPLS[, which.min(CSMWPLS[1,])]

Base.Region4 <- which(between(wavenumbers_NoCO2, 2993, 2999) == TRUE)

#Run PLSR on Base Region 1234 and Append for General Comparison
base1234.fit <- plsr(logKd ~ 
                      MIR[ , Base.Region] + 
                      MIR[ , Base.Region2] +
                      MIR[ , Base.Region3] +
                      MIR[ , Base.Region4],
                    ncomp = 4, 
                    data = train, 
                    validation = "CV", 
                    scale = FALSE, 
                    center = TRUE,
                    method = "simpls"
)

#Strip the "MIR[]" off of the returned dimension names for easier plotting later
dimnames(Base.Region1.2.fit[[3]])[[1]] <- gsub(".*]", "", dimnames(Base.Region1.2.fit[[3]])[[1]])

RMSEC.table$Base.Region1234 <- t(
  data.frame(RMSEP(base1234.fit, estimate = "CV", intercept = FALSE)$val)
)
####END####

####Compare Predictive Values of Each Fitting Model####

#RMSEP.$Truncated <- alex.RMSEP(model = truncated.fit, data = test, target = "XlogP3") 

RMSEP.table$MovingWindow <- alex.RMSEP(model = combined.fit, data = test, target = "logKd") 

RMSEP.table$Window1 <- alex.RMSEP(model = window1.fit, data = test, target = "logKd")

RMSEP.table$Window2 <- alex.RMSEP(model = window2.fit, data = test, target = "logKd")

RMSEP.table$Window3 <- alex.RMSEP(model = window3.fit, data = test, target = "logKd")

RMSEP.table$Base1 <- alex.RMSEP(model = base1.fit, data = test, target = "logKd")

RMSEP.table$Base12 <- alex.RMSEP(model = base12.fit, data = test, target = "logKd")

RMSEP.table$Base123 <- alex.RMSEP(model = base123.fit, data = test, target = "logKd")

RMSEP.table$Base1234 <- alex.RMSEP(model = base1234.fit, data = test, target = "logKd")


#Plot the line of Measured vs. Best Prediction
plot(base1234.fit, ncomp = 2, asp = 1, line = TRUE, newdata = test)
plot(plsr.fit, ncomp = 2, asp = 1, line = TRUE, newdata = test)


#### ->Left off here on June 16
RMSEP(base1234.fit, ncomp = 2, newdata = test)

#Get Actual Prediction Values
predict(base1234.fit, ncomp = 2, newdata = test)

#Plot Loadings of First Two Components
dimnames(base123.fit[[3]])[[1]] <- gsub(".*]", "", dimnames(base123.fit[[3]])[[1]])

plot(base123.fit, 
     plottype = "loadings", 
     comps = 1:2, 
     #legendpos = "topleft",
     labels = "numbers", 
     xlab = "wavenumber",
     ylim = c(-1.2, 1.2),
     xlim = c(650, 3500)
)
abline(h = 0)
minor.tick(nx = 5)

plot(base123.fit, plottype = "scores", comps = 1:3)#, labels = "names")
