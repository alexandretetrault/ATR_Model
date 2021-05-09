#Load Packages
library(baseline)
library(signal)
library(pls)
library(caTools)
library(dplyr)
library(Hmisc)

#### Input Code for Compounds ####

#Load ATR Spectral Data File Names
file_names <- list.files(pattern = "*.asp")

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
spectra <- lapply(spectra, function(x){
  x$Wavenumber <- round(x$Wavenumber)
  return(x)
}
)

#Add ATR Data
spectra <- sapply(file_names,
                  USE.NAMES = TRUE,
                  simplify = FALSE,
                  function(x){
  cbind(spectra[[x]], read.csv(x,
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
      degree = 6
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

for(i in seq(0, length(spectra)-3, 3)){
  avg <- as.vector((spectra[[i+1]][4] + spectra[[i+2]][4] + spectra[[i+3]][4])/3)
  spectra_avg <- cbind(spectra_avg, avg)
  file_names <- c(file_names, Compounds[i +1])
  
}
spectra_avg <- data.frame(spectra_avg[ , -1])
colnames(spectra_avg) <- file_names

#Write the file names to csv
write.csv(file_names, file = "File Names.csv")

##Add Files + KD values to Target Variables.csv if necessary
#Load Target Variables Data
Targets <- read.csv("Target Variables.csv", sep = ",", row.names = 1)
###END###

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

col_set <- rainbow(ncol(spectra_avg))
matplot(wavenumbers,
        spectra_avg,
        type = "l",
        xlim = wavenumber,
        col = col_set
        )
minor.tick(nx = 5)
legend("top", legend = names(spectra_avg), lty=c(1,1), lwd=c(2.5,2.5), col = col_set)
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

#####[Deprecated]Plot Normalized Absorbance ATR Spectra####
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

####Export Spectra for Peak Fitting####
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
###End###


####Dataframe Preparation for PLS - Full Spectrum####
spectra.df <- data.frame(row.names = names(spectra_avg))

spectra.df$Kd <- Targets$Kd

##[UNDECIDED]Log-transform KD Values
spectra.df$logKd <- log10(Targets)
###END###


####[DEPRECATED]Make Nested Matrix of MIR Spectra for All Compounds####
#For Normalized Absorbance
MIR <- sapply(spectra, function(x){
  MIR <- as.matrix(x$Normalized) })

MIR <- t(MIR)

colnames(MIR) <- spectra[[1]][ ,1]

spectra.df$MIR <- MIR

rm(MIR)

rownames(spectra.df) <- spectra.df[ ,1]

spectra.df <- spectra.df[ ,-1]
###END###


####Make Nested Matrix of MIR Spectra for All Compounds####
MIR <- as.matrix(spectra_avg)
MIR <- t(MIR)
spectra.df$MIR <- MIR
###END###
####[UNDECIDED]Add Molecular Weight Tag####
spectra.df$MW <- c(1,0,1,1,1,1,1,1,1,1,1,0,0,1,1,1,0,0,1,1,1,1,1,0,1,1,1)
###END###


####Calculate Second Derivative for Peak Assignment####
MIR.deriv <- apply(spectra.df$MIR, 
                   MARGIN = 2, 
                   FUN = sgolayfilt, 
                   p = 2, n = 5, m = 2, ts = 1
                   )
spectra.df$MIR.deriv <- MIR.deriv

matplot(wavenumbers[c(cm3600_2200, cm1900_650), ],
        t(MIR.deriv),
        type = "l"
)
###END###




#####[DEPRECATED]For Derivative of Absorbance####
MIR.deriv <-  sapply(spectra, function(x){
  MIR.deriv <- as.matrix(x$Deriv) })
MIR.deriv <- t(MIR.deriv)

colnames(MIR.deriv) <- spectra[[1]][ ,1]

spectra.df$MIR.deriv <- MIR.deriv

rm(MIR.deriv)
###END###




####Remove CO2 Region & Split Data for Training on Spectra####
#Make the Wavenumbers Scale
wavenumbers <- as.vector(spectra[[1]][1])

#Spectral Regions of Interest
cm3600_2200 <- which(between(wavenumbers[[1]], 2200, 3600) == TRUE)
cm1900_650 <- which(between(wavenumbers[[1]], 650, 1900) == TRUE)

spectra.df$MIR <- spectra.df$MIR[, c(cm3600_2200, cm1900_650)]
set.seed(123)
split <- sample.split(spectra.df$Kd, SplitRatio = 0.7)

train <- spectra.df
train <- subset(spectra.df, split == T)
test <- subset(spectra.df, split == F)
###END###


####Run Model on Full Spectrum####
plsr.fit <- plsr(logKd ~ MIR,
                     ncomp = 5,
                     data = train,
                     validation = "CV",
                     scale = FALSE,
                     center = TRUE,
                     method = "simpls" 
)

summary(plsr.fit)

plot(RMSEP(plsr.fit))

plot(plsr.fit, plottype = "scores", comps = 1:3, labels = "numbers")

#For only one component, plot u vs. t
plot(plsr.fit[[2]][,1], 
     plsr.fit[[4]][,1],
     xlab = "t",
     ylab = "u"
)

text(x = plsr.fit[[2]][,1], 
     y = plsr.fit[[4]][,1], 
     labels = paste(gsub("_.*","", names(plsr.fit[[2]][,1]))), 
     pos = 3
)
     
##Plot the coefficients
plot(plsr.fit[[1]],
     type = "l")

#For one component only
plot(wavenumbers[c(cm3600_2200, cm1900_650), ],
     plsr.fit[[1]][(1422*3):((1422*3)+1421)],
     type = "l",
     xlim = c(3600,650),
     xlab = "Wavenumber",
     ylab = "Coefficient",
     main = "Regression Coefficients by Wavenumber"
     )
abline(h = 0)
minor.tick(nx = 5)

##Test to see if coefficients are indeed the regression coefficients Y onto X
#Account for mean centering by calculating the intercept to add
coef1 <- coef(plsr.fit, intercept = TRUE, ncomp = 1)[[1]]

(test$MIR)%*%plsr.fit[[1]][1:1422]+coef1 

predict(plsr.fit, ncomp = 3, newdata = test)

#Store Predicted Values as log-Untransformed
predicted_values <- 10^predict(plsr.fit, ncomp = 3, newdata = test)

plot(predicted_values,
     10^test[,1],
     type = "p",
     asp = 1,
     axes = FALSE
     )
abline( a = 0, b = 1)
axis(1, pos=0)
axis(2, pos=0)


#Plot fit for training data
plot(plsr.fit, ncomp = 3, line = TRUE,
     xlim = c(1, 4),
     ylim = c(1, 4))#,
     #axes = FALSE
#)
#axis(1, pos=0)
#axis(2, pos=0)

text(x = spectra.df[,1], 
     y = plsr.fit[[9]][28:54], 
     labels = paste(gsub("_.*","", names(plsr.fit[[2]][,2]))), 
     pos = 3
)

#Obtain Predicted KD values of Training Data
predict(plsr.fit, ncomp = 3)

#Predict KD values of Test Data
predict(plsr.fit, ncomp = 3, newdata = test)

plot(plsr.fit, ncomp = 3, line = TRUE, newdata = test,
     xlim = c(1,4),
     ylim = c(1,4)
     )

RMSEP(plsr.fit, ncomp = 3, newdata = test)

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


####Build Dataframe of Test Prediction Errors for Later Use####
RMSEP. <- data.frame(
  alex.RMSEP(model = plsr.fit, data = test, target = "XlogP3")
)

colnames(RMSEP.) <- "Full"

#Compare homebrew RMSEP to Package RMSEP function on Full Spectrum Model
RMSEP.$Full.Stock <- t(data.frame(
  RMSEP(plsr.fit, estimate = "test", newdata = test, intercept = FALSE)$val)
)
###END###


####Moving Window Partial Least Squares####
##Set Total Number of Spectral Points (Full or w/o CO2)

#total.points <- length(spectra.df$MIR[1,])
total.points <- length(region.df$MIR[1,])

##Set the Wavenumber Scale (Full or w/o CO2)

#wavenumbers <- as.vector(spectra[[1]][1])
wavenumbers <- data.frame(wavenumbers[c(cm3600_2200, cm1900_650), ])

#Create an object to store the results
residues <- NULL

for(i in 1:(total.points - 4)){
  #Set the midpoint of the moving window of 5 spectral points
  midpoint <- wavenumbers[(i+2), ]
  
  #Train the Model on Scaled Absorbance Data
  mwindow.fit <- plsr(XlogP3 ~ MIR[ ,i:(i+4)],
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
        ylim = c(2.05, 2.45), 
        lty = 1,
        ylab = "log(RSS)",
        xlab = "Wavenumber"
)
minor.tick(nx = 5)
title(main = "Traces des Residus, Taille de Fenêtre = 5")
####END####

####Define Windows####
#Window 1
#Define Spectral Region 1
window.1 <- which(between(wavenumbers[[1]], 1100, 1350) == TRUE)

#Fit PLS
#Train the Model on Windows of Absorbance Data
window1.fit <- plsr(XlogP3 ~ MIR[ , window.1],
                        ncomp = 10, 
                        data = train, 
                        validation = "LOO", 
                        scale = FALSE, 
                        center = TRUE,
                        method = "simpls"
)

#Calculate RMSEP Window 1
plot(RMSEP(window1.fit))
RMSEP(window1.fit)

#Window 2
#Define Spectral Region 2
window.2 <- which(between(wavenumbers[[1]], 1350, 1600) == TRUE)

#Fit PLS
#Train the Model on Windows of Absorbance Data
window2.fit <- plsr(XlogP3 ~ MIR[ , window.2],
                    ncomp = 10, 
                    data = train, 
                    validation = "LOO", 
                    scale = FALSE, 
                    center = TRUE,
                    method = "simpls"
)

#Calculate RMSEP Window 2
plot(RMSEP(window2.fit))
RMSEP(window2.fit)

#Window 3
#Define Spectral Region 2
window.3 <- which(between(wavenumbers[[1]], 2500, 3200) == TRUE)

#Fit PLS
#Train the Model on Windows of Absorbance Data
window3.fit <- plsr(XlogP3 ~ MIR[ , window.3],
                    ncomp = 10, 
                    data = train, 
                    validation = "LOO", 
                    scale = FALSE, 
                    center = TRUE,
                    method = "simpls"
)

#Calculate RMSEP Window 2
plot(RMSEP(window3.fit))
RMSEP(window3.fit)

#Build Model on all windows combined
combined.fit <- plsr(XlogP3 ~ 
                       MIR[ , window.1] + 
                       MIR[ , window.2] +
                       MIR[ , window.3],
                     ncomp = 10,
                     data = train, 
                     validation = "LOO",
                     scale = FALSE, 
                     center = TRUE,
                     method = "simpls"
)
####END####

####RMSEC Data Comparison####
#Combine Calibration Prediction Data of all Models for Comparison
RMSEC. <- t(data.frame(RMSEP(truncated.fit, estimate = "CV", intercept = FALSE)$val))

RMSEC. <- data.frame(RMSEC.)

colnames(RMSEC.)[colnames(RMSEC.) == "CV"] <- "Full"

#RMSEC.$Truncated <- t(
#  data.frame(RMSEP(truncated.fit, estimate = "CV", intercept = FALSE)$val)
#)

RMSEC.$MovingWindow <- t(
  data.frame(RMSEP(combined.fit, estimate = "CV", intercept = FALSE)$val)
)

RMSEC.$Window1 <- t(
  data.frame(RMSEP(window1.fit, estimate = "CV", intercept = FALSE)$val)
)

RMSEC.$Window2 <- t(
  data.frame(RMSEP(window2.fit, estimate = "CV", intercept = FALSE)$val)
)

RMSEC.$Window3 <- t(
  data.frame(RMSEP(window3.fit, estimate = "CV", intercept = FALSE)$val)
)
  

####CSMWPLS on Base Region####
#Based on RMSEC curve, ideal no. components should be __

##CSMWPLS on Base Region (Window 3)

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
    base1.fit <- plsr(XlogP3 ~ 
                        MIR[ , window[j]:window[j + w - 1]],
                      ncomp = 1, 
                      data = train, 
                      validation = "LOO", 
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
        wavenumbers[window[j], ],
        "-",
        wavenumbers[window[j + w - 1], ]
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

Base.Region <- which(between(wavenumbers[[1]], 3012, 3017) == TRUE)

set.seed(123)
#Run PLSR on Base Region and Append for General Comparison
base1.fit <- plsr(XlogP3 ~ MIR[ , Base.Region],
                        ncomp = 3, 
                        data = train, 
                        validation = "LOO", 
                        scale = FALSE, 
                        center = TRUE,
                        method = "simpls"
)



#for(i in length(RMSEC.[[1]])){
#  temp <- t(data.frame(
#    RMSEP(Base.Region.fit, estimate = "CV", intercept = FALSE)$val
#    ))
#  c
#}

RMSEC.$Base.Region <- t(
  data.frame(RMSEP.$val)
)
####END####

####SCSMWPLS on Base Region + Window_####

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
    base12.fit <- plsr(XlogP3 ~ 
                         MIR[ , Base.Region] +
                         MIR[ , window[j]:window[j + w - 1]],
                       ncomp = 1, 
                       data = train, 
                       validation = "LOO", 
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
        wavenumbers[window[j], ],
        "-",
        wavenumbers[window[j + w - 1], ]
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

Base.Region2 <- which(between(wavenumbers[[1]], 1159, 1165) == TRUE)

#Run PLSR on Base Region 1 & 2 and Append for General Comparison
base12.fit <- plsr(XlogP3 ~ MIR[ , Base.Region] + MIR[ , Base.Region2],
                           ncomp = 8, 
                           data = train, 
                           validation = "LOO", 
                           scale = FALSE, 
                           center = TRUE,
                           method = "simpls"
)

#Strip the "MIR[]" off of the returned dimension names for easier plotting later
dimnames(Base.Region1.2.fit[[3]])[[1]] <- gsub(".*]", "", dimnames(Base.Region1.2.fit[[3]])[[1]])

RMSEC.$Base.Region1.2 <- t(
  data.frame(RMSEP(base12.fit, estimate = "CV", intercept = FALSE)$val)
)
####END####

####SCSMWPLS on Base Region1.2 + Window _####
#Select window
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
    base123.fit <- plsr(XlogP3 ~ 
                          MIR[ , Base.Region] +
                          MIR[ , Base.Region2] +
                          MIR[ , window[j]:
                                 window[j + w - 1]],
                        ncomp = 1, 
                        data = train, 
                        validation = "LOO", 
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
        wavenumbers[window[j], ],
        "-",
        wavenumbers[window[j + w - 1], ]
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

Base.Region3 <- which(between(wavenumbers[[1]], 1411, 1418) == TRUE)

#Run PLSR on Base Region 123 and Append for General Comparison
base123.fit <- plsr(XlogP3 ~ 
                      MIR[ , Base.Region] + 
                      MIR[ , Base.Region2] +
                      MIR[ , Base.Region3],
                    ncomp = 10, 
                    data = train, 
                    validation = "LOO", 
                    scale = FALSE, 
                    center = TRUE,
                    method = "simpls"
)

#Strip the "MIR[]" off of the returned dimension names for easier plotting later
dimnames(Base.Region1.2.fit[[3]])[[1]] <- gsub(".*]", "", dimnames(Base.Region1.2.fit[[3]])[[1]])

RMSEC.$Base.Region1.2.3 <- t(
  data.frame(RMSEP(base123.fit, estimate = "CV", intercept = FALSE)$val)
)
####END####

####Compare Predictive Values of Each Fitting Model####

RMSEP.$Truncated <- alex.RMSEP(model = truncated.fit, data = test, target = "XlogP3") 

RMSEP.$MovingWindow <- alex.RMSEP(model = combined.fit, data = test, target = "XlogP3") 

RMSEP.$Window1 <- alex.RMSEP(model = window1.fit, data = test, target = "XlogP3")

RMSEP.$Window2 <- alex.RMSEP(model = window2.fit, data = test, target = "XlogP3")

RMSEP.$Window3 <- alex.RMSEP(model = window3.fit, data = test, target = "XlogP3")

RMSEP.$Base1 <- alex.RMSEP(model = base1.fit, data = test, target = "XlogP3")

RMSEP.$Base12 <- alex.RMSEP(model = base12.fit, data = test, target = "XlogP3")

RMSEP.$Base123 <- alex.RMSEP(model = base123.fit, data = test, target = "XlogP3")


#Plot the line of Measured vs. Best Prediction
plot(base123.fit, ncomp = 1, asp = 1, line = TRUE, newdata = test)
plot(truncated.fit, ncomp = 1, asp = 1, line = TRUE, newdata = test)

#And second-best fit:
plot(gas.window2.fit, ncomp = 2, asp = 1, line = TRUE, newdata = test)

#Get Actual Prediction Values
predict(base123.fit, ncomp = 1, newdata = test)

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
