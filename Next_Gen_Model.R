#Load Packages
library(baseline)
library(signal)
library(pls)
library(caTools)
library(dplyr)
library(Hmisc)
library(HotellingEllipse)
library(RColorBrewer)

####Homemade RMSEP Function####
alex.RMSEP <- function(model, data, target){
  temp <- data.frame(predict(model, newdata = data))
  temp <- temp - data[[target]]
  temp[,] <- temp[,]^2
  temp <- colSums(temp)
  temp <- temp/length(data[[target]])
  temp <- sqrt(temp)
}
####END####

####Peak-Finding Algorithm####
#https://stats.stackexchange.com/questions/22974/how-to-find-local-peaks-valleys-in-a-series-of-data
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
####END####

####All Salts Matrix-Matched####
salt_file_names <- list.files(path = "Salts", pattern = "*.asp")

salt_spectra <- list()

for (i in salt_file_names){
  salt_spectra[[i]] <- list()
}

salt_spectra <- lapply(salt_spectra, function(x){
  x$Wavenumber <- read.csv("Wavenumbers.csv", 
                           sep = ",", 
                           col.names = "Wavenumber")
}
)

#Add ATR Data
salt_spectra <- sapply(salt_file_names,
                  USE.NAMES = TRUE,
                  simplify = FALSE,
                  function(x){
                    cbind(salt_spectra[[x]], 
                          read.csv(paste0("./Salts/", x),
                                   skip = 6,
                                   sep = ",",
                                   col.names = "Transmittance"
                                   )
                          )
                  }
)

#Transform Transmittance to Absorbance
salt_spectra <- lapply(salt_spectra, function(x){
  x$Absorbance <- as.vector(
    2-log10(x$Transmittance)
  )
  return(x)
}
)

#Baseline Correct Absorbance 
salt_spectra <- lapply(salt_spectra, function(x){ 
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
salt_spectra <- lapply(salt_spectra, function(x){
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

#Plot Individual Spectra
wavenumber <- rev(range(salt_spectra[[1]][1]))
wavenumbers <- as.vector(salt_spectra[[1]][1])

lapply(names(salt_spectra), function(x) {
  plot(salt_spectra[[x]]$Wavenumber, 
       #x$Normalized, 
       salt_spectra[[x]]$Absorbance,
       type = "l", 
       xlim = wavenumber,
       main = x
       )
  }
)
####END####

####Single Salt Correction####
salt_spectra <- data.frame(
  read.csv("Wavenumbers.csv", 
                        sep = ",",
                        col.names = "Wavenumber"
                        ),
  Transmittance = read.csv("./Salts/May31_TPC_06-04-0_LMW_salt.asp",
                           skip = 6,
                           sep = ",",
                           col.names = "Transmittance"
                           )
)

#Transform Transmittance to Absorbance
salt_spectra$Absorbance <- as.vector(
    2-log10(salt_spectra$Transmittance)
  )

#Baseline Correct Absorbance 
salt_spectra$Absorbance <- t(
    baseline.modpolyfit(
      t(salt_spectra$Absorbance),
      degree = 4
    )[[2]]
)

#Min-Max Normalize Absorbance
salt_spectra$Normalized <- as.vector(
    scale(salt_spectra$Absorbance,
          center = FALSE,
          scale = max(salt_spectra$Absorbance) -
            min(salt_spectra$Absorbance)
    )
)

plot(salt_spectra$Wavenumber,
     salt_spectra$Normalized,
     type = "l",
     xlim = wavenumber)
####END####

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

#Plot
lapply(spectra, function(x) {
  plot(x$Wavenumber, 
       x$Normalized, 
       type = "l", 
       xlim = wavenumber
       )
  }
)
####END####

####[DEPRECATED]Correct LMW Using Matrix-Matched####
##If Salt Spectra Name Matches Spectra Name, Substract the Former Absorbance from the Latter
for(i in names(spectra)){
  for(j in names(salt_spectra)){
    if(paste(gsub("\\_salt.asp*", "", j)) == paste(gsub("\\(.*", "", i))){
      spectra[[i]][4] <- spectra[[i]][4] - salt_spectra[[j]][4]
    }
  }
}

#Change Negative Absorbance No.'s to Zero
for(i in grep("LMW", names(spectra))){
  for(j in 1:nrow(spectra_avg)){
    if(spectra_avg[[i]][j] < 0){
      spectra_avg[[i]][j] <- 0
    }
  }
}


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
####END####

####Correct LMW Using Matrix-Matched####
#First Plot Salt Matrix vs. Averaged Spectra
for(i in names(spectra_avg)){
  for(j in names(salt_spectra)){
    if(paste(gsub("\\_salt.asp*", "", j)) == i){
      plot(wavenumbers[,1],
           spectra_avg[[i]],
           type = "l",
           xlim = wavenumber,
           main = i)
      lines(wavenumbers[,1],
            salt_spectra[[j]][,4],
            type = "l",
            col = "red")
      minor.tick(nx = 5)
    }
  }
}

##If Salt Spectra Name Matches Spectra_Avg Name, Substract the Former Absorbance from the Latter
for(i in names(spectra_avg)){
  for(j in names(salt_spectra)){
    if(paste(gsub("\\_salt.asp*", "", j)) == i){
      spectra_avg[[i]] <- spectra_avg[[i]] - salt_spectra[[j]][,4]
    }
  }
}

#Change Negative Absorbance No.'s to Zero
for(i in grep("LMW", names(spectra_avg))){
  for(j in 1:nrow(spectra_avg)){
    if(spectra_avg[[i]][j] <= 0){
      spectra_avg[[i]][j] <- 0
    }
  }
}

#Min-Max Normalize Again
for(i in grep("LMW", names(spectra_avg))){
  spectra_avg[[i]] <- scale(spectra_avg[[i]],
                            center = FALSE,
                            scale = max(spectra_avg[[i]]) -
                              min(spectra_avg[[i]])
  )
}
####END####

####Correct LMW Averaged Spectra with Single Salt Background####
#Subtract Background from LMW Spectra
spectra_avg[grep("LMW", names(spectra_avg))] <- 
  spectra_avg[grep("LMW", names(spectra_avg))] - salt_spectra$Normalized 

#Change Negative Absorbance No.'s to Zero
for(i in grep("LMW", names(spectra_avg))){
  for(j in 1:nrow(spectra_avg)){
  if(spectra_avg[[i]][j] < 0){
    spectra_avg[[i]][j] <- 0
  }
 }
}

#Min-Max Normalize Again
for(i in grep("LMW", names(spectra_avg))){
  spectra_avg[[i]] <- scale(spectra_avg[[i]],
                            center = FALSE,
                            scale = max(spectra_avg[[i]]) -
                              min(spectra_avg[[i]])
                            )
}
####END####

####Plot the ATR Spectra####
wavenumber <- rev(range(spectra[[1]][1]))
wavenumbers <- as.vector(spectra[[1]][1])

#As Transmittance
lapply(names(spectra), function(x) {
  plot(spectra[[x]]$Wavenumber, 
       spectra[[x]]$Transmittance, 
       type = "l", 
       xlim = wavenumber,
       main = x
       )
  }
)

for (i in 1:ncol(spectra_avg)){
  plot(wavenumbers[, 1],
       spectra_avg[ , i],
       type = "l",
       xlim = wavenumber,
       xlab = "Wavenumber",
       ylab = "Absorbance",
       main = names(spectra_avg[i])
  )
  minor.tick(nx = 5)
}

#col_set <- brewer.pal(6, "RdYlGn") 
#col_set <- colorRampPalette(coul)(6)
#pie(rep(1, length(col_set)), col = col_set , main="") 

#rank <- as.factor( as.numeric( cut(spectra.df[1:40,1], 6)))

col_set <- ifelse(spectra.df[1:40,1] > 500, "darkgreen",
                  ifelse(between(spectra.df[1:40,1], 300, 500), "yellow", "red"))

matplot(wavenumbers,
        spectra_avg[,1:40],
        type = "l",
        xlim = c(3600,800),
        col = col_set,
        xlab = expression(Wavenumber ~ "(cm"^{-1}*")"),
        ylab = "Absorbance",
        main = "Overlaid ATR-FTIR Spectra"
)
minor.tick(nx = 5)
legend("topleft", 
       legend = c("Kd > 600", "Kd [300-600]", "Kd < 300"), 
       fill = c("red", "yellow", "blue"), bty = "n"
)
#legend("top",legend = names(spectra_avg), fill = col_set)
####END####

####Compute Isotherm Results####
isotherm.raw <- lapply(paste0("./Isotherms/", file_names, ".csv"),
                       read.delim, 
                       colClasses = c("NULL","numeric","numeric","NULL")
                       )

isotherm.raw <- lapply(isotherm.raw, setNames, c("X","Y"))

#Assign Names
names(isotherm.raw) <- file_names

#Perform Linear Regressions
linear.kd <- lapply(isotherm.raw, lm, formula = Y ~ X)
                    
#Perform Freundlich Fit
#log.isotherm.raw <- lapply(isotherm.raw, log10)

#freundlich.fit <- lapply(log.isotherm.raw, nls, formula = Y ~ logK+N*X, start = list(N = 1, logK = 1), algorithm = "port")

#Extract Kd Values from Slopes
Kd <- data.frame(file_names, 
                 row.names = 1, 
                 Kd = NA, 
                 SDEV = NA,
                 CI = NA, 
                 R2 = NA#,
                 #Kf = NA,
                 #Nf = NA
                 )

for(i in file_names){
  Kd[i,1] <- round((coef(linear.kd[[i]])[2]), 1)
  mod_summary <- summary(linear.kd[[i]])
  Kd[i,2] <- round(mod_summary$coefficients[2,2], digits = 2)
  temp <- confint(linear.kd[[i]])
  Kd[i,3] <- round((temp[2,2]-temp[2,1])/2, digits = 1)
  Kd[i,4] <- round(summary(linear.kd[[i]])$r.squared, digits = 2)
  #Kd[i,4] <- round(as.numeric(10^coef(freundlich.fit[[i]])[2]), digits = 1)
  #Kd[i,5] <- round(as.numeric(coef(freundlich.fit[[i]])[1]), digits = 2)
}

#Plot Isotherm Fits
for(i in names(isotherm.raw)){
  plot(isotherm.raw[[i]][[1]],
       isotherm.raw[[i]][[2]],
       main = i,
       xlab = expression(C[i] ~ "("*mu*"g/mL)"),
       ylab = expression(A[i] ~ "("*mu*"g/g)"),
       lines(isotherm.raw[[i]][[1]], predict(linear.kd[[i]]), col = "green")
  )
  r.sq <- vector("expression",2)
  r.sq[1] = substitute(expression(K[d] == MYVALUE), 
                     list(MYVALUE = round((coef(linear.kd[[i]])[2]), digits = 1)))[2]
  r.sq[2] = substitute(expression(italic(R)^2 == MYOTHERVALUE), 
                     list(MYOTHERVALUE = round(summary(linear.kd[[i]])$r.squared, digits = 2)))[2]
  legend("topleft", legend = r.sq, bty = "n")
}

#Write Kd Values to csv
write.csv(Kd, file = "Kd.csv")
####END####

####Langmuir Isotherms####
langmuir.raw <- read.delim(paste("./Removed/Mar30_TPC_0-0-1_HMW.csv"))
colnames(langmuir.raw)[2:3] <- c("X","Y")
langmuir.raw <- langmuir.raw[-4,]
plot(langmuir.raw[,c("X","Y")])
with(langmuir.raw,lines(X,4000*0.6*X/(1+0.6*X),col='red'))

langmuir.fit <- nls(formula = Y ~ Q*b*X/(1+b*X), data = langmuir.raw, start = list(Q = 4000, b = 1), algorithm = "port")
summary(langmuir.fit)
lines(langmuir.raw$X,predict(langmuir.fit),col='green')

z <- 1/langmuir.raw[,2:3]
plot(Y~X,z)
abline(lm(Y~X,z))
M <- lm(Y~X,z)

Q <- 1/coef(M)[1]

b <- coef(M)[1]/coef(M)[2]

lines(langmuir.raw$X,Q*b*langmuir.raw$X/(1+b*langmuir.raw$X), col = "darkgreen")

langmuir.raw <- read.delim(paste("./Isotherms/Apr06_TPC_1-0-0_HMW.csv"), 
                           colClasses = c("NULL", "numeric","numeric","NULL"),
                           header = TRUE
                           )
colnames(langmuir.raw) <- c("X","Y")
z <- 1/langmuir.raw
plot(Y~X,z)
abline(lm(Y~X,z))
M <- lm(Y~X,z)
Q <- 1/coef(M)[1]
b <- coef(M)[1]/coef(M)[2]

plot(langmuir.raw[,c("X","Y")])
lines(langmuir.raw$X,Q*b*langmuir.raw$X/(1+b*langmuir.raw$X), col = "darkgreen")

with(langmuir.raw,lines(X,5000*0.075*X/(1+0.075*X),col='red'))

langmuir.fit <- nls(formula = Y ~ Q*b*X/(1+b*X), data = langmuir.raw, start = list(Q = mean(langmuir.raw$Y), b = 0.075), algorithm = "port")
summary(langmuir.fit)
lines(langmuir.raw$X,predict(langmuir.fit),col='green')
####END####

####Freundlich Isotherm####
langmuir.raw <- read.delim(paste("./Isotherms/Apr06_TPC_1-0-0_HMW.csv"), 
                           colClasses = c("NULL", "numeric","numeric","NULL"),
                           header = TRUE
)
colnames(langmuir.raw) <- c("X","Y")
z <- log10(langmuir.raw)
freundlich.fit <- nls(formula = Y ~ logK+N*X, data = z, start = list(N = 1, logK = 1), algorithm = "port")
summary(freundlich.fit)
K <-as.numeric(10^coef(freundlich.fit)[2])
N <- as.numeric(coef(freundlich.fit)[1])
langmuir.raw$CN <- langmuir.raw$X^N
plot(langmuir.raw[,1:2])
lines(langmuir.raw$X,K*langmuir.raw$CN,col = "green")

langmuir.raw <- read.delim(paste("./Isotherms/Mar30_TPC_0-0-1_HMW.csv"), 
                           colClasses = c("NULL", "numeric","numeric","NULL"),
                           header = TRUE
)
colnames(langmuir.raw) <- c("X","Y")
z <- log10(langmuir.raw)
freundlich.fit <- nls(formula = Y ~ logK+N*X, data = z, start = list(N = 1, logK = 1), algorithm = "port")
summary(freundlich.fit)
K <-as.numeric(10^coef(freundlich.fit)[2])
N <- as.numeric(coef(freundlich.fit)[1])
langmuir.raw$CN <- langmuir.raw$X^N
plot(langmuir.raw[,1:2])
lines(langmuir.raw$X,K*langmuir.raw$CN,col = "red")
####END####

####Add Compounds + Kd values to Target Variables.csv if necessary####
#Write the file names of compounds to csv
write.csv(file_names, file = "File Names.csv")


## -> Manually copy file names to Target Variables.csv and add in Kd data 

#Load Target Variables Data
Targets <- read.csv("Target Variables.csv", sep = ",", row.names = 1)
####END####

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

####Remove Diamond Region & Split Data for Training on Spectra####
#Make the Wavenumbers Scale
wavenumbers <- as.vector(spectra[[1]][1])

#Spectral Regions of Interest
cm3600_2200 <- which(between(wavenumbers[[1]], 2200, 3600) == TRUE)
cm1900_650 <- which(between(wavenumbers[[1]], 650, 1900) == TRUE)
wavenumbers_NoCO2 <- wavenumbers[c(cm3600_2200, cm1900_650), 1]

spectra.df_NoSed <- spectra.df

spectra.df_NoSed$MIR <- spectra.df$MIR[, c(cm3600_2200, cm1900_650)]

#Make sediment test dataframe
sediment <- spectra.df_NoSed[(nrow(spectra.df_NoSed)-5):nrow(spectra.df_NoSed),]

#Remove sediment from train dataframe
spectra.df_NoSed<- spectra.df_NoSed[-((nrow(spectra.df_NoSed)-5):nrow(spectra.df_NoSed)),]

#spectra.df_NoSed <- cbind(spectra.df_NoSed, n = 1:nrow(spectra.df_NoSed))
#Remove LMW Compounds
#spectra.df_NoSed <- spectra.df_NoSed[-c(4:6,15,16,19,25,31,34:40),]

#Split between train and test df
set.seed(123)
#set.seed(222)
split <- sample.split(spectra.df_NoSed$Kd, SplitRatio = 0.8)

train <- subset(spectra.df_NoSed, split == T)

test <- subset(spectra.df_NoSed, split == F)

#Number rows of train dataframe
train <- cbind(train, n = 1:nrow(train))

####END####

####Calculate Second Derivative for Peak Assignment####
MIR.deriv <- apply(spectra.df$MIR, 
                   MARGIN = 1, 
                   FUN = sgolayfilt, 
                   p = 2, n = 11, m = 2, ts = 1
)

spectra.df$MIR.deriv <- t(MIR.deriv)

#Create a List of All Peaks Calculated for 2nd Derivative Spectra
peaks.list <- apply(
  spectra.df$MIR.deriv, 
  MARGIN = 1, 
  function(i)wavenumbers[find_peaks(-i, m = 4),]
) 

#Retain Only Region Below 1800 cm-1
peaks.list <- lapply(peaks.list, function(i)i[which(i < 1800)])

#Extract Peaks for End-Members
#TNOM HMW
TNOM_HMW_peaks <- peaks.list$`Apr06_TPC_1-0-0_HMW`
#Plot 
plot(wavenumbers[,1],
     spectra_avg$`Apr06_TPC_1-0-0_HMW`,
     type = "l",
     xlim = c(1800, 800),
     main = "Apr6_TPC_1-0-0_HMW"
)

minor.tick(nx = 4)
abline(v = TNOM_HMW_peaks, col = "red")
#Triage Peaks by Eye
TNOM_HMW_peaks <- c(TNOM_HMW_peaks[c(1,5,11,15,25:27,30:34,39:40)],1560,1397,1389)

#CNOM HMW
CNOM_HMW_peaks <- peaks.list$`Mar30_TPC_0-0-1_HMW`
#Plot 
plot(wavenumbers[,1],
     spectra_avg$`Mar30_TPC_0-0-1_HMW`,
     type = "l",
     xlim = c(1800, 800),
     main = "Mar30_TPC_0-0-1_HMW"
)

minor.tick(nx = 4)
abline(v = CNOM_HMW_peaks, col = "red")
#Triage Peaks by Eye
CNOM_HMW_peaks <- c(CNOM_HMW_peaks[c(4,6,8:10,17,21,25:27,29,36)],1543,1509,1410)

#C2NOM HMW
C2NOM_HMW_peaks <- peaks.list$`Jun16_TPC_0-0-1_HMW`
#Plot 
plot(wavenumbers[,1],
     spectra_avg$`Jun16_TPC_0-0-1_HMW`,
     type = "l",
     xlim = c(1800, 650),
     main = "Jun16_TPC2_0-0-1_HMW"
)

minor.tick(nx = 4)
abline(v = C2NOM_HMW_peaks, col = "red")
#Triage Peaks by Eye
C2NOM_HMW_peaks <- c(C2NOM_HMW_peaks[c(6,7,8,10,12,14,15,17,18,21,23,28,30,33)],1730,1654,1026)

#PNOM HMW
PNOM_HMW_peaks <- peaks.list$`Mar22_TPC_0-1-0_HMW`
#Plot 
plot(wavenumbers[,1],
     spectra_avg$`Mar22_TPC_0-1-0_HMW`,
     type = "l",
     xlim = c(1800, 800),
     main = "Mar22_TPC_0-1-0_HMW"
)

minor.tick(nx = 4)
abline(v = PNOM_HMW_peaks, col = "red")
#Triage Peaks by Eye
PNOM_HMW_peaks <- c(PNOM_HMW_peaks[c(5,12,13,16,23,25,29,30,35,39,40)],1637,1543,1458,1397,1047,975)

#TNOM LMW
TNOM_LMW_peaks <- peaks.list$`Apr08_TPC_1-0-0_LMW`
#Plot 
plot(wavenumbers[,1],
     spectra_avg$`Apr08_TPC_1-0-0_LMW`,
     type = "l",
     xlim = c(1800, 650),
     main = "Apr08_TPC_1-0-0_LMW"
)

minor.tick(nx = 4)
abline(v = TNOM_LMW_peaks, col = "red")
#Triage Peaks by Eye
TNOM_LMW_peaks <- c(TNOM_LMW_peaks[c(3,9,23,26,28:31,33,38,40)],1433,1370)

#CNOM LMW
CNOM_LMW_peaks <- peaks.list$`Apr12_TPC_0-0-1_LMW`
#Plot 
plot(wavenumbers[,1],
     spectra_avg$`Apr12_TPC_0-0-1_LMW`,
     type = "l",
     xlim = c(1800, 650),
     main = "Apr12_TPC_0-0-1_LMW"
)

minor.tick(nx = 4)
abline(v = CNOM_LMW_peaks, col = "coral", lty = 2)
#Triage Peaks by Eye
CNOM_LMW_peaks <- c(CNOM_LMW_peaks[c(6,9,11,14,17,20,23,27,31,36,37,39,44,45)],1108)

#PNOM LMW
PNOM_LMW_peaks <- peaks.list$`Apr09_TPC_0-1-0_LMW`
#Plot 
plot(wavenumbers[,1],
     spectra_avg$`Apr09_TPC_0-1-0_LMW`,
     type = "l",
     xlim = c(1800, 650),
     main = "Apr09_TPC_0-1-0_LMW"
)

minor.tick(nx = 4)
abline(v = PNOM_LMW_peaks, col = "red")
#Triage Peaks by Eye
PNOM_LMW_peaks <- c(PNOM_LMW_peaks[c(5,14,16:18,20,26,28:30,32,33,35,38)],1570,1508)

#Sediment
STN19_peaks <- peaks.list$`Sediment_STN19_Jun23`
#Plot 
plot(wavenumbers[,1],
     spectra_avg$`Sediment_STN19_Jun23`,
     type = "l",
     xlim = c(3600, 650),
     main = "Sediment_STN19_Jun23"
)

minor.tick(nx = 4)
abline(v = STN19_peaks, col = "red")
#Triage Peaks by Eye
STN19_peaks <- c(STN19_peaks[c(8:10,28,33)], 1115)

#Sediment
SAG05_peaks <- peaks.list$`Sediment_SAG05_Jun09`
#Plot 
plot(wavenumbers[,1],
     spectra_avg$`Sediment_SAG05_Jun09`,
     type = "l",
     xlim = c(1800, 650),
     main = "Sediment_SAG05_Jun09"
)

minor.tick(nx = 4)
abline(v = SAG05_peaks, col = "red")
#Triage Peaks by Eye
SAG05_peaks <- c(SAG05_peaks[c(8:10,12,20,37)],1115)

peaks.df <- data.frame(TNOM_HMW = TNOM_HMW_peaks[1:20], 
                       PNOM_HMW = PNOM_HMW_peaks[1:20], 
                       CNOM_HMW = CNOM_HMW_peaks[1:20],
                       C2NOM_HMW = C2NOM_HMW_peaks[1:20],
                       TNOM_LMW = TNOM_LMW_peaks[1:20], 
                       PNOM_LMW = PNOM_LMW_peaks[1:20], 
                       CNOM_LMW = CNOM_LMW_peaks[1:20],
                       STN19 = STN19_peaks[1:20],
                       SAG05 = SAG05_peaks[1:20]
)

write.csv(peaks.df, "Endmember_peaks.csv")

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

####Plot HMW Fingerprint Region####
#Graphically parameter reset
dev.off()
#or
par(mfrow=c(1,1))

par(mar=c(0.5, 0.5, 0.2, 0.2), 
    mfrow=c(1,2),
    oma = c(5, 5, 0.2, 0.2)
    )
layout(matrix(c(1,2,2), nrow = 1, ncol = 3))

plot(wavenumbers[,1],
     spectra_avg$`Apr06_TPC_1-0-0_HMW`,
     type = "l",
     #xlim = c(1800, 800),
     xlim = c(3600, 2800),
     ylim = c(0,2.3),
     xlab = expression("Wavenumber (cm"^-1*")"),
     ylab = "Absorbance",
     col = "dark green"
)
mtext("Absorbance", side = 2, line = 2.5)
mtext(expression("Wavenumber"~italic("(cm"^-1*")")), 
      side = 1, 
      outer = TRUE,
      line = 2.5
      )
lines(wavenumbers[,1],
      spectra_avg$`Mar22_TPC_0-1-0_HMW`+0.5,
      col = "dark blue"
)

lines(wavenumbers[,1],
      spectra_avg$`Mar30_TPC_0-0-1_HMW`+1,
      col = "red"
)

minor.tick(nx = 2)

abline(v = c(3275,3060,2920,2850),lty = 2, col = "darkgrey")

text(x = c(3275,3060,2920,2850), 
     y = 2.3, 
     labels = c(3275,3060,2920,2850),
     cex = 1
)

text(x = c(3275,3060,2885), 
     y = 2.2, 
     labels = c(expression(nu*"(O-H)"),
                expression(nu*"(C-H)"[italic(sp^2)]),
                expression(nu*"(C-H)"[italic(Aliph.)])
                ),
     cex = 1
)

legend("topleft", 
       legend = c("Corn", "Plankton", "Soil"), 
       col = c("red","darkblue","darkgreen"),
       lty = 1,
       bty = "n"
)

#1800-800 Range HMW
plot(wavenumbers[,1],
     spectra_avg$`Apr06_TPC_1-0-0_HMW`,
     type = "l",
     #xlim = c(1800, 800),
     xlim = c(1800, 800),
     ylim = c(0,2.3),
     xlab = expression("Wavenumber (cm"^-1*")"),
     yaxt = "n",
     ylab = "",
     col = "dark green"
)

lines(wavenumbers[,1],
      spectra_avg$`Mar22_TPC_0-1-0_HMW`+0.5,
      col = "dark blue"
)

lines(wavenumbers[,1],
      spectra_avg$`Mar30_TPC_0-0-1_HMW`+1,
      col = "red"
)

minor.tick(nx = 4)

abline(v = c(1630,1576,1560,1511,1448,1397,1314,1239,1124,1038,975), lty = 2, col = "darkgrey")

text(x = c(1630,1511,1397,1239,1038), 
     y = 2.3, 
     labels = c(1630,1511,1397,1239,1038),
     cex = 1
)

text(x = c(1653,1511,1397,1239,1038), 
     y = 2.2, 
     labels = c(expression(nu*"(C=O)"[italic("Amide I")]),
                expression(nu*"(C=C)"[italic(Arom.)]),
                expression(nu*"(C-OH)"[italic(Phenol)]),
                expression(delta*"(COH)"[italic(Polysacch.)]),
                expression(italic(Polysacch.))
     ),
     cex = 1
)

text(x = c(1576,1448,1314,1124,975), 
     y = 2.1, 
     labels = c(1576,1448,1314,1124,975),
     cex = 1
)

text(x = c(1576,1448,1314,975), 
     y = 2, 
     labels = c(expression(nu*"(C=C)"[italic(Arom.)]),
                expression(delta*"(C-H)"),
                expression(delta*"(CH"[2]*")"[italic(Polysacch.)]),
                expression(nu[s]*"(C-O-C)"[italic(Glycosidic)])
      ),
     cex = 1
)






plot(wavenumbers[,1],
     spectra_avg$`Apr06_TPC_1-0-0_HMW`,
     type = "l",
     #xlim = c(1800, 800),
     xlim = c(3600, 800),
     ylim = c(0,2.5),
     xlab = expression("Wavenumber (cm"^-1*")"),
     ylab = "Absorbance",
     col = "dark green"
)

lines(wavenumbers[,1],
      spectra_avg$`Mar22_TPC_0-1-0_HMW`+0.5,
      col = "dark blue"
)

lines(wavenumbers[,1],
      spectra_avg$`Mar30_TPC_0-0-1_HMW`+1,
      col = "red"
)

minor.tick(nx = 4)

abline(v = c(1650,1637,1595,1574,1560,1541,1509,1458,1397,1389,1305,1239,1124,1075,1038,975),lty = 2, col = "gray")

text(x = c(1653,1574), 
     y = 2.5, 
     labels = c(1653,1574),
     cex = 0.70
)
text(x = c(1653,1574), 
     y = 2.4, 
     labels = c(expression(nu*"(C=O)"["Amide I"]),
                expression(nu*"(C=C)"[Aromatic])
              
     ),
     cex = 0.70
)




text(x = c(1637,1560), 
     y = 2.3, 
     labels = c(1637,1560),
     cex = 0.70
)

text(x = c(1637,1560), 
     y = 2.2, 
     labels = c(expression(nu[as]*"(COO"^"-"*")")
                
                
     ),
     cex = 0.70
)

text(x = c(1595), 
     y = 2.1, 
     labels = c(1595),
     cex = 0.70
)

text(x = c(1595), 
     y = 2.0, 
     labels = c(expression(nu*"(C=C)"[Aromatic])
                
                
     ),
     cex = 0.70
)

text(x = c(1653,1595,1560,1509,1397,1305,1124,1038), 
     y = 2, 
     labels = c(1653,1595,1560,1509,1397,1305,1124,1038),
     cex = 0.70
     )
text(x = c(1637,1574,1541,1458,1389,1239,1075,975), 
     y = 1.9, 
     labels = c(1637,1574,1543,1458,1389,1239,1075,975), 
     cex = 0.70
     )
#abline(v = c(1731,1651,1637,1593,1560,1510,1390,1239,1124,1085,1075,1038),lty =2, col = "gray")

legend("topleft", 
       legend = c("Corn", "Plankton", "Soil"), 
       col = c("red","darkblue","darkgreen"),
       lty = 1,
       bty = "n"
)
####END####

####Plot LMW Fingerprint Region####
plot(wavenumbers[,1],
     spectra_avg$`Apr08_TPC_1-0-0_LMW`,
     type = "l",
     xlim = c(1800, 800),
     ylim = c(0,2.5),
     xlab = expression("Wavenumber (cm"^-1*")"),
     ylab = "Absorbance",
     col = "dark green"
)

lines(wavenumbers[,1],
      spectra_avg$`Apr09_TPC_0-1-0_LMW`+0.5,
      col = "dark blue"
)

lines(wavenumbers[,1],
      spectra_avg$`Apr12_TPC_0-0-1_LMW`+1.1,
      col = "red"
)

minor.tick(nx = 4)

abline(v = c(1762,1654,1636,1570,1508,1408,1370,1120,1108,1049), lty =2, col = "gray")
#abline(v = c(1762,1683,1654,1636,1570,1560,1542,1508,1458,1433,1408,1398,1371,1355,1191,1120,1108,1049,990,879,865,833,824), lty =2, col = "gray")

text(x = c(1762,1636,1508,1408), 
     y = 2.5, 
     labels = c(1762,1636,1508,1408),
     cex = 0.70
)
text(x = c(1762,1636,1508,1408), 
     y = 2.4, 
     labels = c(expression(nu*"(C=O)"[vinyl/lactone]),
                expression(nu[as]*"(COO"^"-"*")"),
                expression(nu*"(C=C)"[arom]),
                expression(nu*"(C-O)"[phenol])
                ),
     cex = 0.70
)

text(x = c(1654,1570,1370), 
     y = 2.2, 
     labels = c(1654,1570,1370),
     cex = 0.70
)

text(x = c(1654,1570,1370), 
     y = 2.1, 
     labels = c(expression(nu*"(C=O)"["Amide I"]),
                expression(nu*"(C=C)"[arom]),
                expression(delta*"(C-H)"[aliphatic])
                
     ),
     cex = 0.70
)

text(x = c(1637,1570,1541,1458,1389,1239,1075,975), 
     y = 1.9, 
     labels = c(1637,1570,1543,1458,1389,1239,1075,975), 
     cex = 0.70
)

legend("topleft", 
       legend = c("Corn", "Plankton", "Soil"), 
       col = c("red","darkblue","darkgreen"),
       lty = 1,
       bty = "n"
)
####END####

####Plot Sediment Overlaid####
#Plot 
plot(wavenumbers[,1],
     spectra_avg$`Sediment_SAG05_Jun09`,
     type = "l",
     xlim = c(3600, 650),
     ylim = c(0, 2),
     main = "Sediment_STN19_Jun23",
     xlab = expression("Wavenumber (cm"^-1*")"),
     ylab = "Absorbance",
     col = "dark green"
)

lines(wavenumbers[,1],
      spectra_avg$`Sediment_STN19_Jun23`+0.5,
      col = "dark blue"
)

minor.tick(nx = 5)

legend("topleft", 
       legend = c("SAG05", "STN19"), 
       col = c("darkgreen","darkblue"),
       lty = 1,
       bty = "n"
)

abline(v = c(3511,3349,3243,1636,1610,1424,1115,995), lty =2, col = "gray")
####END####

####Plot Fresh -> Degraded Fingerprint Region####
plot(wavenumbers[,1],
     spectra_avg$`Apr06_TPC_1-0-0_HMW`,
     type = "l",
     xlim = c(1800, 800),
     ylim = c(0,2.2),
     xlab = expression("Wavenumber (cm"^-1*")"),
     ylab = "Absorbance",
     col = "dark green"
)

lines(wavenumbers[,1],
      spectra_avg$`Mar30_TPC_0-0-1_HMW`+0.5,
      col = "darkred"
)

lines(wavenumbers[,1],
      spectra_avg$`Jun16_TPC_0-0-1_HMW`+1.2,
      col = "red"
)

minor.tick(nx = 4)

abline(v = 1023, lwd = 3, col = "green")

abline(v = 1372, lwd = 3, col = "red")

legend("topleft", 
       legend = c("Fresh Corn", "Aged Corn", "Soil"), 
       col = c("red","darkred","darkgreen"),
       lty = 1
)
####END####

####Run Model on Full Spectrum####
set.seed(123)
plsr.fit <- plsr(logKd ~ MIR,
                     ncomp = 10,
                     data = spectra.df_NoSed,
                     #data = train,
                     validation = "CV",
                     segments = 5,
                     scale = FALSE,
                     center = TRUE,
                     method = "simpls" 
)

#Change MIR Dimension Names for Easier Plotting
#dimnames(plsr.fit[[3]])[[1]] <- wavenumbers_NoCO2

summary(plsr.fit)

#RMSEP by Cross-Validation
plot(RMSEP(plsr.fit), legendpos = "topright")
minor.tick(nx = 2)

#Inspect Hotelling's T^2 for Outliers
plsr.scores <- data.frame(comp1 = plsr.fit$scores[,1], 
                          comp2 = plsr.fit$scores[,2],
                          comp3 = plsr.fit$scores[,3]
)

plsr.scores <- plsr.scores %>%
  as_tibble() %>%
  print()

T2.ellipse <- ellipseParam(data = plsr.scores, k = 3, pcx = 1, pcy = 2)

plsr.scores %>%
  ggplot(aes(x = comp1, y = comp2)) +
  #ggforce::geom_ellipse(aes(x0 = 0, y0 = 0, a = T2.ellipse$Ellipse$a.99pct, b = T2.ellipse$Ellipse$b.99pct, angle = 0), size = .5, linetype = "dashed", fill = "white") +
  ggforce::geom_ellipse(aes(x0 = 0, y0 = 0, a = T2.ellipse$Ellipse$a.95pct, b = T2.ellipse$Ellipse$b.95pct, angle = 0), size = .5, linetype = "dotted", fill = "white") +
  geom_point(aes(fill = T2.ellipse$Tsquare$value), shape = 21, size = 3, color = "black") +
  scale_fill_viridis_c(option = "viridis") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = .2) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = .2) 
#labs(title = "Scatterplot of PCA scores", subtitle = "PC1 vs. PC2", x = "PC1", y = "PC2", fill = "T2", caption = "Figure 1: Hotelling's T2 ellipse obtained\n using the ellipseParam function") +
#theme_grey()

tibble(
  T2 = purrr::pluck(T2.ellipse, "Tsquare", "value"), 
  obs = 1:nrow(plsr.scores)
) %>%
  ggplot() +
  geom_point(aes(x = obs, y = T2, fill = T2), shape = 21, size = 3, color = "black") +
  geom_segment(aes(x = obs, y = T2, xend = obs, yend = 0), size = .5) +
  scale_fill_gradient(low = "black", high = "red", guide = "none") +
  geom_hline(yintercept = purrr::pluck(T2.ellipse, "cutoff.99pct"), linetype = "dashed", color = "darkred", size = .5) +
  geom_hline(yintercept = purrr::pluck(T2.ellipse, "cutoff.95pct"), linetype = "dashed", color = "darkblue", size = .5) +
  annotate("text", x = 15, y = 12.4, label = "99% limit", color = "darkred") +
  annotate("text", x = 15, y = 8.6, label = "95% limit", color = "darkblue") +
  labs(x = "Observations", y = "Hotelling’s T-square (3 PCs)", caption = "Figure 3: Hotelling’s T2-value vs. Observations") +
  theme_bw() 

#RMSEP(plsr.fit, newdata = test)

#RMSEP(plsr.fit, ncomp = 3, newdata = sediment)

pls::R2(plsr.fit)

predplot(plsr.fit, 
         ncomp = 2, 
         line = TRUE,
         main = "Model Validation",
         xlim = c(1,3.5),
         ylim = c(1,3.5),
         pch = 21,
         bg = "gray"
         )
####END####

####Plot the Scores & Loadings####
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

loadingplot(plsr.fit, comps = 1:6)
####END####

####Test to see if coefficients are indeed the regression coefficients Y onto X####
#Account for mean centering by calculating the intercept to add
coef1 <- coef(plsr.fit, intercept = TRUE, ncomp = 2)[[1]]

(test$MIR)%*%plsr.fit[[1]][,1,2]+coef1 

predict(plsr.fit, ncomp = 2, newdata = test)
####END####

#####Fit for Training Data####
#plot(plsr.fit, 
#     ncomp = 2, 
#     xlim = c(1, 3.5),
#     ylim = c(1, 3.5),
#     line = TRUE,
     #asp = 1
#     )
#minor.tick(nx = 5, ny = 5)

plot(spectra.df_NoSed$logKd,
    #train[,2],
     predict(plsr.fit, ncomp = 2),
     xlab = "measured",
     ylab = "predicted",
     main = "Validation - Full Spectrum",
     xlim = c(1, 3.5),
     ylim = c(1, 3.5)
     )
abline(a = 0, b = 1)

#Add labels
text(x = train[,2], 
     y = plsr.fit[9]$fitted.value[,1,2], 
     labels = paste(gsub("_.*","", names(plsr.fit[[2]][,2]))), 
     pos = 3
)

#points(test[,2],
#       predict(plsr.fit, ncomp = 2, newdata = test),
#       col = "blue"
#       )

#points(sediment[,2],
#       predict(plsr.fit, ncomp = 2, newdata = sediment),
#       col = "red"
#)

#legend("topleft", 
#       legend = c("train", "test", "sediment"), 
#       fill = c("purple","blue","red")
#       )

RMSEP(plsr.fit)
####END####

####Fit of Test & Sediment Data####
#plot(plsr.fit, ncomp = 2, line = TRUE, newdata = test,
#     xlim = c(1.5,3),
#     ylim = c(1.5,3)
#)

RMSEP(plsr.fit, ncomp = 2, newdata = test)

plot(test[,2],
     predict(plsr.fit, ncomp = 2, newdata = test),
     #asp = 1,
     xlim = c(1.5, 3),
     ylim = c(1.5, 3),
)
abline(a = 0, b = 1)

RMSEP(plsr.fit, ncomp = 2, newdata = sediment)

plot(sediment[,2],
     predict(plsr.fit, ncomp = 3, newdata = sediment),
     #asp = 1,
     xlim = c(1.5, 3),
     ylim = c(1.5, 3),
)
abline(a = 0, b = 1)
####END####

####Compare Predicted vs. Actual Kd Values####
predicted_values <- 10^predict(plsr.fit, ncomp = 2, newdata = test)

plot(test[,1],
     predicted_values,
     type = "p"
     )
abline( a = 0, b = 1)

#Obtain Predicted KD values of Training Data
predict(plsr.fit, ncomp = 2)

#Predict KD values of Test Data
predict(plsr.fit, ncomp = 2, newdata = test)
####END####

####Compare homebrew RMSEP to Package RMSEP function on Full Spectrum Model####
RMSEP.table$Full.Stock <- t(data.frame(
  RMSEP(plsr.fit, estimate = "test", newdata = test, intercept = FALSE)$val)
)
####END####

####Run Model on 650-1800, 2899-2999 cm^-1 Regions Exclusively####
region_650_1800 <- which(between(wavenumbers_NoCO2, 650, 1800) == TRUE)
region_2800_2999 <- which(between(wavenumbers_NoCO2, 2801, 2999) == TRUE)

#region.df <- spectra.df_NoSed

set.seed(123)
region.fit <- plsr(logKd ~
                     MIR[ , region_2800_2999] +
                     MIR[ , region_650_1800],
                   ncomp = 10, 
                   data = spectra.df_NoSed,
                   validation = "CV",
                   segments = 5,
                   scale = FALSE,
                   center = TRUE,
                   method = "simpls"
)

summary(region.fit)

plot(RMSEP(region.fit, estimate = "CV"))

pls::R2(region.fit, estimate = "CV")

print(alex.RMSEP(model = region.fit, data = test, target = "logKd"))

print(alex.RMSEP(model = region.fit, data = sediment, target = "logKd"))

#predplot(region.fit, ncomp = 2, line = TRUE)#, labels = "numbers")

#predplot(base12345.fit, ncomp = 2, line = TRUE)#, labels = "numbers")

plot(spectra.df_NoSed$logKd,
     #train[,2],
     predict(region.fit, ncomp = 2),
     xlab = "measured",
     ylab = "predicted",
     main = "Validation - Full Spectrum",
     xlim = c(1, 3.5),
     ylim = c(1, 3.5)
)
abline(a = 0, b = 1)

text(x = spectra.df_NoSed$logKd, 
     y = region.fit[9]$fitted.value[,1,2], 
     labels = paste(gsub("_.*","", names(region.fit[[2]][,2]))), 
     pos = 3
)
####END####

####Repeated k-Fold####
set.seed(123)
train_control <- trainControl(method = "repeatedcv", 
                              number = 5, 
                              repeats = 3
                              )

foo <- data.frame(spectra.df_NoSed[,2])
foo$MIR <- spectra.df_NoSed$MIR[,c(region_2800_2999 , region_650_1800)]
foo$Kd <- spectra.df_NoSed$Kd
colnames(foo) <- c("logKd", "MIR", "Kd")

set.seed(123)
rpeat.kfold.fit <- train(form = logKd ~
                      MIR,
                   method = "simpls",
                   trControl = train_control,
                   tuneLength = 4, 
                   data = foo,
                   preProcess = "center"
                   )


rpeat.kfold.fit

plot(rpeat.kfold.fit$finalModel, line = TRUE)
sediment$MIR <- sediment$MIR[,c(region_2800_2999 , region_650_1800)]
predictions <- predict(region.fit, newdata = sediment)

RMSEP(region.fit$finalModel, newdata = sediment)

##Sediments
#Log-transformed Q^2
Q2 <- 1-(sum((predictions - sediment$logKd)^2)) / 
  sum((sediment$logKd - mean(sediment$logKd))^2)

#Q^2 on original values
Q2 <- 1-(sum((10^predictions - sediment$Kd)^2)) / 
  sum((sediment$Kd - mean(sediment$Kd))^2)

##OM
1- (sum((region.fit$finalModel$fitted.values[,,2] - spectra.df_NoSed$logKd)^2)) / 
     sum((spectra.df_NoSed$logKd - mean(spectra.df_NoSed$logKd))^2)
####END####

####Moving Window Partial Least Squares####
##Set Total Number of Spectral Points (i.e w/o CO2)
total.points <- length(train$MIR[1,])
#total.points <- length(spectra.df_NoSed$MIR[1,])

#Create an object to store the results
residues <- NULL

for(i in 1:(total.points - 9)){
  #Set the midpoint of the moving window of 10 spectral points
  midpoint <- wavenumbers_NoCO2[i+4]

  #Train the Model on Scaled Absorbance Data
  set.seed(123)
  mwindow.fit <- plsr(logKd ~ MIR[ ,i:(i+9)],
                      ncomp = 10, 
                      data = train,
                      #data = spectra.df_NoSed,
                      validation = "none",
                      #segments = 5,
                      scale = FALSE, 
                      center = TRUE,
                      method = "simpls"
  )
  #Extract the Residual Sum of Squares for the Fit (RSS)
  RSS <- mvrValstats(mwindow.fit, estimate = "train", intercept = FALSE)#, estimate = "CV", intercept = "TRUE")
  RSS <- RSS$SSE
  
  #RSS <- RMSEP(mwindow.fit)#, estimate = "CV", intercept = "TRUE")
  #RSS <- RSS$val

  #Append Data to Dataframe
  residues <- rbind(
    residues, 
    data.frame(midpoint, RSS))
  
}

#Plot Residue Lines by Midpoint of Spectral Window
matplot(residues$midpoint, 
        log10(residues[,-1]), 
        #residues[,-1],
        type = "l", 
        xlim = c(3600,650), 
        #xlim = c(1500,750),
        #ylim = c(-1,1),
        lty = 1,
        ylab = "log(RSS)",
        xlab = "Wavenumber"
)
minor.tick(nx = 5)
#minor.tick(nx = 4)
abline(h = 0)
title(main = "Moving Window PLS Residuals Trace")
legend("bottomleft", 
       legend = c("1 comp", "2 comps", "3 comps"), 
       col = c("black", "red", "green"),
       lty = 1
       )
####END####

####Define Windows & Components####
Num.components <- 6

#Window 1
#HMW Only
window.1 <- which(between(wavenumbers_NoCO2, 900, 1100) == TRUE)
#BEST
window.1 <- which(between(wavenumbers_NoCO2, 900, 999) == TRUE)
#TRIAL
window.1 <- which(between(wavenumbers_NoCO2, 850, 950) == TRUE)

#Window 2
#HMW Only
window.2 <- which(between(wavenumbers_NoCO2, 1150, 1400 ) == TRUE)
#BEST
window.2 <- which(between(wavenumbers_NoCO2, 1000, 1250) == TRUE)
#TRIAL
window.2 <- which(between(wavenumbers_NoCO2, 950, 1100) == TRUE)

#Window 3
#HMW Only
window.3 <- which(between(wavenumbers_NoCO2, 1500, 1700) == TRUE)
#BEST
window.3 <- which(between(wavenumbers_NoCO2, 1250, 1450) == TRUE)
#TRIAL
window.3 <- which(between(wavenumbers_NoCO2, 1300, 1450) == TRUE)

#Window 4
#Define Spectral Region 4
window.4 <- which(between(wavenumbers_NoCO2, 1300, 1425) == TRUE)
##window.4 <- which(between(wavenumbers_NoCO2, 1400, 1450) == TRUE)

#Window 5
#Define Spectral Region 5
window.5 <- which(between(wavenumbers_NoCO2, 1425, 1475) == TRUE)
#window.5 <- which(between(wavenumbers_NoCO2, 920, 975) == TRUE)

#Fit PLS
#Train the Model on Windows of Absorbance Data
set.seed(123)
window1.fit <- plsr(logKd ~ MIR[ , window.1],
                        ncomp = Num.components, 
                        data = train,
                        #data = spectra.df_NoSed,
                        validation = "none",
                        #validation = "CV",
                        segments = 5,
                        scale = FALSE, 
                        center = TRUE,
                        method = "simpls"
)

#Calculate RMSEP Window 1
plot(RMSEP(window1.fit))
RMSEP(window1.fit)

#Fit PLS
#Train the Model on Windows of Absorbance Data
set.seed(123)
window2.fit <- plsr(logKd ~ MIR[ , window.2],
                    ncomp = Num.components, 
                    data = train,
                    #data = spectra.df_NoSed,
                    validation = "none",
                    #validation = "CV",
                    segments = 5,
                    scale = FALSE, 
                    center = TRUE,
                    method = "simpls"
)

#Calculate RMSEP Window 2
plot(RMSEP(window2.fit))
RMSEP(window2.fit)

#Fit PLS
#Train the Model on Windows of Absorbance Data
set.seed(123)
window3.fit <- plsr(logKd ~ MIR[ , window.3],
                    ncomp = Num.components, 
                    data = train,
                    #data = spectra.df_NoSed, 
                    validation = "none",
                    #validation = "CV",
                    segments = 5,  
                    scale = FALSE, 
                    center = TRUE,
                    method = "simpls"
)

#Calculate RMSEP Window 3
plot(RMSEP(window3.fit))
RMSEP(window3.fit)

#Fit PLS
#Train the Model on Windows of Absorbance Data
set.seed(123)
window4.fit <- plsr(logKd ~ MIR[ , window.4],
                    ncomp = Num.components,
                    data = train,
                    #data = spectra.df_NoSed, 
                    #validation = "none",
                    validation = "CV",
                    segments = 5, 
                    scale = FALSE, 
                    center = TRUE,
                    method = "simpls"
)

#Calculate RMSEP Window 4
plot(RMSEP(window4.fit))
RMSEP(window4.fit)

#Fit PLS
#Train the Model on Windows of Absorbance Data
set.seed(123)
window5.fit <- plsr(logKd ~ MIR[ , window.5],
                    ncomp = Num.components,
                    data = train,
                    #data = spectra.df_NoSed, 
                    #validation = "none",
                    validation = "CV",
                    segments = 5, 
                    scale = FALSE, 
                    center = TRUE,
                    method = "simpls"
)

#Calculate RMSEP Window 5
plot(RMSEP(window5.fit))
RMSEP(window5.fit)

#Build Model on all windows combined
set.seed(123)
combined.fit <- plsr(logKd ~ 
                       MIR[ , window.1] + 
                       MIR[ , window.2] +
                       MIR[ , window.3] +
                       MIR[ , window.4] +
                       MIR[ , window.5],
                     ncomp = 10,
                     data = train,
                     #data = spectra.df_NoSed, 
                     validation = "CV",
                     segments = 5,
                     scale = FALSE, 
                     center = TRUE,
                     method = "simpls"
)
RMSEP(combined.fit)
plot(RMSEP(combined.fit))
####END####

####CSMWPLS on Base Region1####
##Based on RMSEC curve of Best Window, ideal no. components should be __
#BEST
Num.components <- 3
#TRIAL
Num.components <- 3

##CSMWPLS on Base Region (Window _)

#Select window
#HMW Only
window <- window.3
#BEST
window <- window.3
#TRIAL
window <- window.2

#Clear container
CSMWPLS <- NULL

#Set size of moving window
for(w in 10:length(window)){
#for(w in 4:length(window)){
    
  #Factor region by window size to determine no. windows that span it
  span <- length(window) - w + 1
  
  for(j in 1:span){
    #Train the Model on Scaled Absorbance Data
    set.seed(123)
    base1.fit <- plsr(logKd ~ 
                        MIR[ , window[j]:window[j + w - 1]],
                      ncomp = Num.components, 
                      data = train,
                      #data = spectra.df_NoSed,
                      #validation = "none",
                      validation = "CV",
                      segments = 5,
                      scale = FALSE, 
                      center = TRUE,
                      method = "simpls"
    )
    
    #Extract the Residual Sum of Squares for the Fit (RSS)
    RSS <- RMSEP(base1.fit, intercept = "FALSE", estimate = "CV")
    RSS <- RSS$val
    
    #RSS <- mvrValstats(base1.fit, estimate = "CV", intercept = "TRUE")
    #RSS <- RSS$SSE
    
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

#HMW Only
Base.Region1 <- which(between(wavenumbers_NoCO2, 1584, 1605) == TRUE)
#BEST
Base.Region1 <- which(between(wavenumbers_NoCO2, 1269, 1335) == TRUE)
#TRIAL
Base.Region1 <- which(between(wavenumbers_NoCO2, 997, 1029) == TRUE)

#Run PLSR on Base Region and Append for General Comparison
set.seed(123)
base1.fit <- plsr(logKd ~ MIR[ , Base.Region1],
                        ncomp = 6, 
                        #ncomp = 4, 
                        data = train,
                        #data = spectra.df_NoSed,
                        #validation = "none",
                        validation = "CV",
                        segments = 5,
                        scale = FALSE, 
                        center = TRUE,
                        method = "simpls"
)

RMSEP(base1.fit)#, estimate = "train")
####END####

####SCSMWPLS on Base Region1 + Window_####
#Select next window
#HMW Only
window <- window.1 
#BEST
window <- window.1 
#TRIAL
window <- window.1 

#Clear container
CSMWPLS <- NULL

#Set size of moving window
for(w in 10:length(window)){
#for(w in 4:length(window)){  
  #Factor region by window size to determine no. windows that span it
  span <- length(window) - w + 1
  
  for(j in 1:span){
    #Train the Model on Scaled Absorbance Data
    set.seed(123)
    base12.fit <- plsr(logKd ~ 
                         MIR[ , window[j]:window[j + w - 1]] +
                         MIR[ , Base.Region1],
                       ncomp = Num.components, 
                       data = train,
                       #data = spectra.df_NoSed,
                       validation = "CV",
                       #validation = "none",
                       segments = 5,
                       scale = FALSE, 
                       center = TRUE,
                       method = "simpls"
    )
    
    #Extract the Residual Sum of Squares for the Fit (RSS)
    RSS <- RMSEP(base12.fit, intercept = "FALSE", estimate = "CV")
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

#HMW Only
Base.Region2 <- which(between(wavenumbers_NoCO2, 903, 925) == TRUE)
#BEST
Base.Region2 <- which(between(wavenumbers_NoCO2, 900, 977) == TRUE)
#TRIAL
Base.Region2 <- which(between(wavenumbers_NoCO2, 892, 947) == TRUE)


#Run PLSR on Base Region 1 & 2 and Append for General Comparison
set.seed(123)
base12.fit <- plsr(logKd ~ MIR[ , Base.Region1] + MIR[ , Base.Region2],
                           ncomp = 10, 
                           data = train,
                           #data = spectra.df_NoSed,
                           #validation = "none",
                           validation = "CV",
                           segments = 5,
                           scale = FALSE, 
                           center = TRUE,
                           method = "simpls"
)

RMSEP(base12.fit)#, estimate = "train")
RMSEP(base1.fit)#, estimate = "train")

print(alex.RMSEP(base1.fit, data = test, target = "logKd"))
print(alex.RMSEP(base12.fit, data = test, target = "logKd"))

#Strip the "MIR[]" off of the returned dimension names for easier plotting later
dimnames(base12.fit[[3]])[[1]] <- gsub(".*]", "", dimnames(base12.fit[[3]])[[1]])
####END####

####SCSMWPLS on Base Region12 + Window _####
#Select window
#HMW Only
window <- window.2
#BEST
window <- window.2
#TRIAL
window <- window.2

#Clear container
CSMWPLS <- NULL

#Set size of moving window
for(w in 10:length(window)){
  
  #Factor region by window size to determine no. windows that span it
  span <- length(window) - w + 1
  
  for(j in 1:span){
    #Train the Model on Scaled Absorbance Data
    set.seed(123)
    base123.fit <- plsr(logKd ~ 
                          MIR[ , Base.Region1] +
                          MIR[ , Base.Region2] +
                          MIR[ , window[j]:
                                 window[j + w - 1]],
                        ncomp = Num.components, 
                        data = train,
                        #data = spectra.df_NoSed,
                        #validation = "none",
                        validation = "CV",
                        segments = 5,
                        scale = FALSE, 
                        center = TRUE,
                        method = "simpls"
    )
    
    #Extract the Residual Sum of Squares for the Fit (RSS)
    RSS <- RMSEP(base123.fit, intercept = "FALSE", estimate = "CV") 
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

#HMW Only
Base.Region3 <- which(between(wavenumbers_NoCO2, 1189, 1234) == TRUE)
#BEST
Base.Region3 <- which(between(wavenumbers_NoCO2, 1144, 1173) == TRUE)
#TRIAL
#Base.Region3 <- which(between(wavenumbers_NoCO2, 1175, 1249) == TRUE)

#Run PLSR on Base Region 123 and Append for General Comparison
set.seed(123)
base123.fit <- plsr(logKd ~ 
                      MIR[ , Base.Region1] + 
                      MIR[ , Base.Region2] +
                      MIR[ , Base.Region3],
                    ncomp = 10, 
                    data = train,
                    #data = spectra.df_NoSed,
                    #validation = "none",
                    validation = "CV",
                    segments = 5,
                    scale = FALSE, 
                    center = TRUE,
                    method = "simpls"
)

RMSEP(base123.fit)#, estimate = "train")
RMSEP(base12.fit)#, estimate = "train")
RMSEP(base1.fit)#, estimate = "train")

print(alex.RMSEP(base1.fit, data = test, target = "logKd"))
print(alex.RMSEP(base12.fit, data = test, target = "logKd"))
print(alex.RMSEP(base123.fit, data = test, target = "logKd"))

#Strip the "MIR[]" off of the returned dimension names for easier plotting later
dimnames(Base.Region1.2.fit[[3]])[[1]] <- gsub(".*]", "", dimnames(Base.Region1.2.fit[[3]])[[1]])
####END####

####SCSMWPLS on Base Region123 + Window 4####
#Select window
window <- window.1

#Clear container
CSMWPLS <- NULL

#Set size of moving window
for(w in 4:length(window)){
  
  #Factor region by window size to determine no. windows that span it
  span <- length(window) - w + 1
  
  for(j in 1:span){
    #Train the Model on Scaled Absorbance Data
    set.seed(123)
    base1234.fit <- plsr(logKd ~
                           MIR[ , window[j]:window[j + w - 1]] +
                           MIR[ , Base.Region2] +
                           MIR[ , Base.Region] +
                           MIR[ , Base.Region3],
                        ncomp = Num.components, 
                        data = train,
                        #data = spectra.df_NoSed,
                        #validation = "none",
                        validation = "CV",
                        segments = 5,
                        scale = FALSE, 
                        center = TRUE,
                        method = "simpls"
    )
    
    #Extract the Residual Sum of Squares for the Fit (RSS)
    RSS <- RMSEP(base1234.fit, intercept = "FALSE", estimate = "CV")
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

Base.Region4 <- which(between(wavenumbers_NoCO2, 814, 821) == TRUE)
#Base.Region4 <- which(between(wavenumbers_NoCO2, 2993, 2999) == TRUE)
##Base.Region4 <- which(between(wavenumbers_NoCO2, 816, 822) == TRUE)

#Run PLSR on Base Region 1234 and Append for General Comparison
set.seed(123)
base1234.fit <- plsr(logKd ~
                      MIR[ , Base.Region2] +
                      MIR[ , Base.Region3] + 
                      MIR[ , Base.Region] + 
                      MIR[ , Base.Region4],
                    ncomp = 4, 
                    data = train,
                    #data = spectra.df_NoSed,
                    #validation = "none",
                    validation = "CV",
                    segments = 5,
                    scale = FALSE, 
                    center = TRUE,
                    method = "simpls"
)

RMSEP(base1234.fit, estimate = "CV")

#Strip the "MIR[]" off of the returned dimension names for easier plotting later
dimnames(Base.Region1.2.fit[[3]])[[1]] <- gsub(".*]", "", dimnames(Base.Region1.2.fit[[3]])[[1]])
####END####

####SCSMWPLS on Base Region1234 + Window 5####
#Select window
window <- window.5

#Clear container
CSMWPLS <- NULL

#Set size of moving window
for(w in 4:length(window)){
  
  #Factor region by window size to determine no. windows that span it
  span <- length(window) - w + 1
  
  for(j in 1:span){
    #Train the Model on Scaled Absorbance Data
    set.seed(123)
    base12345.fit <- plsr(logKd ~
                           MIR[ , window[j]:window[j + w - 1]] +
                           MIR[ , Base.Region4] +
                           MIR[ , Base.Region2] +
                           MIR[ , Base.Region] +
                           MIR[ , Base.Region3],
                         ncomp = Num.components, 
                         data = train,
                         #data = spectra.df_NoSed,
                         #validation = "none",
                         validation = "CV",
                         segments = 5,
                         scale = FALSE, 
                         center = TRUE,
                         method = "simpls"
    )
    
    #Extract the Residual Sum of Squares for the Fit (RSS)
    RSS <- RMSEP(base12345.fit, intercept = FALSE, estimate = "CV")
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

Base.Region5 <- which(between(wavenumbers_NoCO2, 1196, 1240) == TRUE)
#Base.Region5 <- which(between(wavenumbers_NoCO2, 1351, 1406) == TRUE)
##Base.Region5 <- which(between(wavenumbers_NoCO2, 1455, 1467) == TRUE)

#Run PLSR on Base Region 1234 and Append for General Comparison
set.seed(123)
base12345.fit <- plsr(logKd ~
                       MIR[ , Base.Region5] +
                       MIR[ , Base.Region4] +
                       MIR[ , Base.Region2] + 
                       MIR[ , Base.Region] + 
                       MIR[ , Base.Region3],
                     ncomp = 4, 
                     data = train,
                     #data = spectra.df_NoSed,
                     #validation = "none",
                     validation = "CV",
                     segments = 5,
                     scale = FALSE, 
                     center = TRUE,
                     method = "simpls"
)

RMSEP(base12345.fit, estimate = "CV")
RMSEP(base1234.fit, estimate = "CV")

print(alex.RMSEP(base12345.fit, data = test, target = "logKd"))
print(alex.RMSEP(base1234.fit, data = test, target = "logKd"))

#Strip the "MIR[]" off of the returned dimension names for easier plotting later
dimnames(Base.Region1.2.fit[[3]])[[1]] <- gsub(".*]", "", dimnames(Base.Region1.2.fit[[3]])[[1]])
####END####

####Fully Algorithmic SCSMWPLS####
#Select window
window <- region_800_1800

#Clear container
CSMWPLS <- NULL

#Set size of moving window
for(w in 10:length(window)){
  
  #Factor region by window size to determine no. windows that span it
  span <- length(window) - w + 1
  
  for(j in 1:span){
    #Train the Model on Scaled Absorbance Data
    set.seed(123)
    mw.fit <- plsr(logKd ~ 
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
    RSS <- RMSEP(mw.fit, estimate = "CV", intercept = "TRUE")
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

#Query for Spectral Windows of lowest RMSEP 
sort(CSMWPLS[1,])[1:100]

which.min(CSMWPLS[1,])

CSMWPLS[, which.min(CSMWPLS[1,])]

CSMWPLS.1 <- CSMWPLS

df <- data.frame(matrix(ncol = 0, nrow = 2))

for(i in 9:17){
  for(j in 1:10){
    temp1 <- names(which.min(CSMWPLS.1[1, grep(paste0("^",i), names(CSMWPLS.1[1,]))]))
    temp2 <- min(CSMWPLS.1[1, grep(paste0("^",i), CSMWPLS.1(CSMWPLS.1[1,]))]) 
    df <- cbind(df, c(temp1, temp2))
    CSMWPLS.1 <- CSMWPLS.1[ ,-which(colnames(CSMWPLS.1) == temp1)]
  }
}

which.min(CSMWPLS[1, grep("^9", names(CSMWPLS[1,]))])
min(CSMWPLS[1, grep("^9", names(CSMWPLS[1,]))])

which.min(CSMWPLS[1, grep("^10", names(CSMWPLS[1,]))])
min(CSMWPLS[1, grep("^10", names(CSMWPLS[1,]))])

which.min(CSMWPLS[1, grep("^11", names(CSMWPLS[1,]))])
min(CSMWPLS[1, grep("^11", names(CSMWPLS[1,]))])

which.min(CSMWPLS[1, grep("^12", names(CSMWPLS[1,]))])
min(CSMWPLS[1, grep("^12", names(CSMWPLS[1,]))])

which.min(CSMWPLS[1, grep("^13", names(CSMWPLS[1,]))])
min(CSMWPLS[1, grep("^13", names(CSMWPLS[1,]))])

which.min(CSMWPLS[1, grep("^14", names(CSMWPLS[1,]))])
min(CSMWPLS[1, grep("^14", names(CSMWPLS[1,]))])

which.min(CSMWPLS[1, grep("^15", names(CSMWPLS[1,]))])
min(CSMWPLS[1, grep("^15", names(CSMWPLS[1,]))])

algo.window.1 <- which(between(wavenumbers_NoCO2, 926, 960) == TRUE)
algo.window.2 <- which(between(wavenumbers_NoCO2, 1026, 1044) == TRUE)
algo.window.3 <- which(between(wavenumbers_NoCO2, 1394, 1503) == TRUE)

#Run PLSR on Algorithmically-Determined Windows and Append for General Comparison
set.seed(123)
algo.fit <- plsr(logKd ~ 
                       MIR[ , algo.window.1] + 
                       MIR[ , algo.window.2] +
                       MIR[ , algo.window.3],
                     ncomp = 4, 
                     data = train, 
                     validation = "CV", 
                     scale = FALSE, 
                     center = TRUE,
                     method = "simpls"
)

#Strip the "MIR[]" off of the returned dimension names for easier plotting later
#dimnames(Base.Region1.2.fit[[3]])[[1]] <- gsub(".*]", "", dimnames(Base.Region1.2.fit[[3]])[[1]])

summary(algo.fit)

print(alex.RMSEP(algo.fit, test, "logKd"))

print(alex.RMSEP(algo.fit, sediment, "logKd"))
####END####

####RMSEC Data Comparison for Optimal Window Selection####
#Combine Calibration Prediction Data of all Models for Comparison
RMSEC.table <- t(data.frame(RMSEP(plsr.fit, estimate = "CV", intercept = FALSE, ncomp = 1:4)$val))

RMSEC.table <- data.frame(RMSEC.table)

colnames(RMSEC.table)[colnames(RMSEC.table) == "CV"] <- "Full Spectrum"

RMSEC.table$FunctionalRegion <- t(
  data.frame(RMSEP(region.fit, estimate = "CV", intercept = FALSE)$val)
)

RMSEC.table$AllWindows <- t(
  data.frame(RMSEP(combined.fit, estimate = "CV", intercept = FALSE)$val)
)

RMSEC.table$Base.Region <- t(
  data.frame(RMSEP(base1.fit, estimate = "train", intercept = FALSE)$val)
)

RMSEC.table$Base.Region12 <- t(
  data.frame(RMSEP(base12.fit, estimate = "train", intercept = FALSE)$val)
)

RMSEC.table$Base.Region123 <- t(
  data.frame(RMSEP(base123.fit, estimate = "train", intercept = FALSE)$val)
)

RMSEC.table$Base.Region1234 <- t(
  data.frame(RMSEP(base1234.fit, estimate = "CV", intercept = FALSE)$val)
)

RMSEC.table$AlgoFit <- t(
  data.frame(RMSEP(algo.fit, estimate = "CV", intercept = FALSE)$val)
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

####Compare Test Predictive Values of Each Fitting Model####
RMSEP.table <- data.frame(
  alex.RMSEP(model = plsr.fit, data = spectra.df_NoSed, target = "logKd")[1:4]
)

colnames(RMSEP.table) <- "Full Spectrum"

RMSEP.table$FunctionalRegion <- alex.RMSEP(model = region.fit, data = test, target = "logKd") 

RMSEP.table$AllWindows <- alex.RMSEP(model = combined.fit, data = test, target = "logKd") 

RMSEP.table$Base1 <- alex.RMSEP(model = base1.fit, data = test, target = "logKd")

RMSEP.table$Base12 <- alex.RMSEP(model = base12.fit, data = test, target = "logKd")

RMSEP.table$Base123 <- alex.RMSEP(model = base123.fit, data = test, target = "logKd")

RMSEP.table$Base1234 <- alex.RMSEP(model = base1234.fit, data = test, target = "logKd")

RMSEP.table$AlgoFit <- alex.RMSEP(model = algo.fit, data = test, target = "logKd")

RMSEP.table$Window1 <- alex.RMSEP(model = window1.fit, data = test, target = "logKd")

RMSEP.table$Window2 <- alex.RMSEP(model = window2.fit, data = test, target = "logKd")

RMSEP.table$Window3 <- alex.RMSEP(model = window3.fit, data = test, target = "logKd")
####END####

####Compare Sediment Predictive Values of Each Fitting Model####
RMSEP.sediment <- data.frame(
  alex.RMSEP(model = plsr.fit, data = sediment, target = "logKd")
)

colnames(RMSEP.sediment) <- "Full Spectrum"

RMSEP.sediment$FunctionalRegion <- alex.RMSEP(model = region.fit, data = sediment, target = "logKd") 

RMSEP.sediment$AllWindows <- alex.RMSEP(model = combined.fit, data = sediment, target = "logKd") 

RMSEP.sediment$Base1 <- alex.RMSEP(model = base1.fit, data = sediment, target = "logKd")

RMSEP.sediment$Base12 <- alex.RMSEP(model = base12.fit, data = sediment, target = "logKd")

RMSEP.sediment$Base123 <- alex.RMSEP(model = base123.fit, data = sediment, target = "logKd")

RMSEP.sediment$Base1234 <- alex.RMSEP(model = base1234.fit, data = sediment, target = "logKd")

RMSEP.sediment$AlgoFit <- alex.RMSEP(model = algo.fit, data = sediment, target = "logKd")

RMSEP.sediment$Window1 <- alex.RMSEP(model = window1.fit, data = sediment, target = "logKd")

RMSEP.sediment$Window2 <- alex.RMSEP(model = window2.fit, data = sediment, target = "logKd")

RMSEP.sediment$Window3 <- alex.RMSEP(model = window3.fit, data = sediment, target = "logKd")
####END####

#Explore Best Fits
####Hotelling's T^2####
base123.scores <- data.frame(comp1 = base123.fit$scores[,1], 
                             comp2 = base123.fit$scores[,2],
                             comp3 = base123.fit$scores[,3]
)

base123.scores <- base123.scores %>%
  as_tibble() %>%
  print()

T2.ellipse <- ellipseParam(data = base123.scores, k = 2, pcx = 1, pcy = 2)
#Note: Ellipse is a circle b/c scores are normalized. This was verified
#by normalizing the GitHub example data, which then produced a circle
#instead of an ellipse

base123.scores %>%
  ggplot(aes(x = comp1, y = comp2)) +
  #ggforce::geom_ellipse(aes(x0 = 0, y0 = 0, a = T2.ellipse$Ellipse$a.99pct, b = T2.ellipse$Ellipse$b.99pct, angle = 0), size = .5, linetype = "dashed", fill = "white") +
  ggforce::geom_ellipse(aes(x0 = 0, y0 = 0, a = T2.ellipse$Ellipse$a.95pct, b = T2.ellipse$Ellipse$b.95pct, angle = 0), size = .5, linetype = "dotted", fill = "white") +
  geom_point(aes(fill = T2.ellipse$Tsquare$value), shape = 21, size = 3, color = "black") +
  scale_fill_viridis_c(option = "viridis") +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", size = .2) +
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = .2) 
  #labs(title = "Scatterplot of PCA scores", subtitle = "PC1 vs. PC2", x = "PC1", y = "PC2", fill = "T2", caption = "Figure 1: Hotelling's T2 ellipse obtained\n using the ellipseParam function") +
  #theme_grey()

tibble(
    T2 = purrr::pluck(T2.ellipse, "Tsquare", "value"), 
    obs = 1:nrow(base123.scores)
) %>%
ggplot() +
geom_point(aes(x = obs, y = T2, fill = T2), shape = 21, size = 3, color = "black") +
geom_segment(aes(x = obs, y = T2, xend = obs, yend = 0), size = .5) +
scale_fill_gradient(low = "black", high = "red", guide = "none") +
geom_hline(yintercept = purrr::pluck(T2.ellipse, "cutoff.99pct"), linetype = "dashed", color = "darkred", size = .5) +
geom_hline(yintercept = purrr::pluck(T2.ellipse, "cutoff.95pct"), linetype = "dashed", color = "darkblue", size = .5) +
annotate("text", x = 15, y = 12.4, label = "99% limit", color = "darkred") +
annotate("text", x = 15, y = 8.6, label = "95% limit", color = "darkblue") +
labs(x = "Observations", y = "Hotelling’s T-square (3 PCs)", caption = "Figure 3: Hotelling’s T2-value vs. Observations") +
theme_bw()  
####END####

####RMSEP####
summary(plsr.fit)
summary(base123.fit)
RMSEP(plsr.fit, estimate = "train")
RMSEP(base123.fit, estimate = "CV")
RMSEP(base12.fit, estimate = "CV")

plot(RMSEP(plsr.fit))
plot(RMSEP(base12.fit))

summary(base12.fit)
summary(base123.fit)
summary(base1234.fit)

pls::R2(plsr.fit)
pls::R2(base123.fit, estimate = "train", newdata = test)

print(alex.RMSEP(model = plsr.fit, data = spectra.df_NoSed, target = "logKd"))
print(alex.RMSEP(model = base123.fit, data = test, target = "logKd"))
print(alex.RMSEP(model = base1234.fit, data = test, target = "logKd"))
print(alex.RMSEP(model = base12.fit, data = test, target = "logKd"))

print(alex.RMSEP(model = plsr.fit, data = sediment, target = "logKd"))
print(alex.RMSEP(model = base1234.fit, data = sediment, target = "logKd"))

####END####

####R^2 & Q^2####
RQ <- data.frame(c("Full", "MW"),
                 row.names = 1,
                 R2 = NA, 
                 Q2_test = NA, 
                 Q2_sed = NA
                 )


pls::R2(plsr.fit, ncomp = 3, estimate = "train")#, estimate = "CV", intercept = TRUE)
pls::R2(region.fit, ncomp = 3, estimate = "CV", intercept = FALSE)
pls::R2(base1234.fit, ncomp = 3, estimate = "CV", intercept = FALSE)

RQ[1,1] <- round(1-(
  sum((plsr.fit$fitted.values[,1,3] - spectra.df_NoSed$logKd)^2) / 
    sum((spectra.df_NoSed$logKd - mean(spectra.df_NoSed$logKd))^2)
), digits = 2
)

#RQ[1,2] <- round(1-(
#  sum((predict(plsr.fit, ncomp = 3, newdata = test) - test$logKd)^2) / 
#    sum((test$logKd - mean(test$logKd))^2)
#), digits = 2
#)

#RQ[1,3] <- round(1-(
#  sum((predict(plsr.fit, ncomp = 3, newdata = sediment) - sediment$logKd)^2) / 
#    sum((sediment$logKd - mean(sediment$logKd))^2)
#  
#), digits = 2
#)

RQ[2,1] <- round(
  1-
  sum((base123.fit$fitted.values[,1,2] - train$logKd)^2) / 
  sum((train$logKd - mean(train$logKd))^2),
  digits = 2
)

RQ[2,2] <- round(1-(
  sum((predict(base123.fit, ncomp = 2, newdata = test) - test$logKd)^2) / 
    sum((test$logKd - mean(test$logKd))^2)
  
), digits = 2
)

#RQ[2,3] <- round(
#  1-(
#  sum((predict(base123.fit, ncomp = 2, newdata = sediment) - sediment$logKd)^2) / 
#    sum((sediment$logKd - mean(sediment$logKd))^2)
#  ), digits = 2)

RQ[3,1] <- round(1-(
  sum((base1234.fit$fitted.values[,1,3] - train$logKd)^2) / 
    sum((train$logKd - mean(train$logKd))^2)
), digits = 2
)

RQ[3,2] <- round(1-(
  sum((predict(base1234.fit, ncomp = 3, newdata = test) - test$logKd)^2) / 
    sum((test$logKd - mean(test$logKd))^2)
  
), digits = 2
)

RQ[3,3] <- round(1-(
  sum((predict(base1234.fit, ncomp = 3, newdata = sediment) - sediment$logKd)^2) / 
    sum((sediment$logKd - mean(sediment$logKd))^2)
  
), digits = 2
)

pls::R2(region.fit, ncomp = 1:4, estimate = "CV", intercept = FALSE)

pls::R2(region.fit, ncomp = 1:4, estimate = "test", newdata = test, intercept = FALSE)


#Log-transformed Q^2
Q2 <- 1-(sum((predictions - sediment$logKd)^2)) / 
  sum((sediment$logKd - mean(sediment$logKd))^2)

#Q^2 on original values
Q2 <- 1-(sum((10^predictions - sediment$Kd)^2)) / 
  sum((sediment$Kd - mean(sediment$Kd))^2)

##OM
1- (sum((region.fit$finalModel$fitted.values[,,2] - spectra.df_NoSed$logKd)^2)) / 
  sum((spectra.df_NoSed$logKd - mean(spectra.df_NoSed$logKd))^2)
####END####

####Fit Plots####

#Full Spectra
predplot(plsr.fit,
         ncomp = 3,
         line = TRUE,
         xlim = c(1, 3.5),
         ylim = c(1, 3.5),
         xlab = "measured log(Kd)",
         ylab = "predicted log(Kd)",
         pch = 21,
         bg = "gray"
)
legend("topleft", legend = "Train", fill = "gray")

pls::R2(plsr.fit)

RMSEP(plsr.fit)

#Add points for Model Fit of Built on Full Data
points(spectra.df_NoSed[,2],
       predict(plsr.fit, ncomp = 3, estimate = "train"),
       pch = 21,
       bg = "blue"
)

#legend("topleft", 
#       legend = c("Train", "Test"), 
#       fill = c("gray","blue")
#)

#plot(spectra.df_NoSed[,2],
#     predict(plsr.fit, ncomp = 3),
#     xlab = "measured",
#     ylab = "predicted",
#     main = "Model Validation",
#     xlim = c(1, 3.5),
#     ylim = c(1, 3.5),
#     pch = 21,
#     bg = "gray"
#)
#abline(a = 0, b = 1)
#legend("topleft", legend = "Train", fill = "gray")


RMSEP(plsr.fit, newdata = test)
print(alex.RMSEP(plsr.fit, data = test, target = "logKd"))

#points(sediment[,2],
#       predict(plsr.fit, ncomp = 3, newdata = sediment),
#       pch = 21,
#       bg = "coral"
#)

#legend("topleft", 
#       legend = c("Train", "Test", "Sediment"), 
#       fill = c("gray","blue", "coral")
#)

#RMSEP(plsr.fit, newdata = sediment)
#print(alex.RMSEP(plsr.fit, data = sediment, target = "logKd"))

##Best MW Fit
predplot(base123.fit,
         ncomp = 6,
         line = TRUE,
         xlim = c(1, 3.5),
         ylim = c(1, 3.5),
         pch = 21,
         bg = "gray"
         )
legend("topleft", legend = "Train", fill = "gray")

pls::R2(base123.fit, estimate = "CV")
RMSEP(base1.fit)

#Add Data Fit Prediction on Full Dataset
points(train[,2],
       predict(base12.fit, ncomp = 6),
       pch = 21,
       bg = "blue"
)

#plot(train[,2],
#     predict(base123.fit, ncomp = 4),
#     predict(plsr.fit, ncomp = 3),
#     predict(base1234.fit, ncomp = 3),
#     xlab = "measured",
#     ylab = "predicted",
#     main = "Model Validation",
#     xlim = c(1, 3.5),
#     ylim = c(1, 3.5),
#     pch = 21,
#     bg = "gray"
#)
#abline(a = 0, b = 1)
#legend("topleft", legend = "Train", fill = "gray")

points(test[,2],
       predict(base123.fit, ncomp = 6, newdata = test),
       pch = 21,
       #bg = "pink"
       bg = "blue"
)

legend("topleft", 
       legend = c("Train", "Test"), 
       fill = c("gray","blue")
)

print(alex.RMSEP(base123.fit, data = test, target = "logKd"))

points(sediment[,2],
       predict(base123.fit, ncomp = 4, newdata = sediment),
       pch = 21,
       bg = "coral"
)

print(alex.RMSEP(base1234.fit, data = sediment, target = "logKd"))


##Add labels
text(x = train[,2], 
     y = region.fit[9]$fitted.value[,1,2], 
     labels = paste(gsub("_.*","", names(plsr.fit[[2]][,2]))), 
     pos = 3
)


#Add labels
text(x = test[,2], 
     y = predict(region.fit, ncomp = 2, newdata = test), 
     labels = paste(gsub("_.*","", row.names(test))), 
     pos = 3
)
#Add labels
text(x = sediment[,2], 
     y = predict(region.fit, ncomp = 2, newdata = sediment), 
     labels = paste(gsub(".*_","", row.names(sediment))), 
     pos = 3
)
####END####

####Score Plot####
plot(plsr.fit, plottype = "scores", comps = 1:3, labels = "numbers")
plot(base1234.fit, plottype = "scores", comps = 1:3, labels = "numbers")

#Plot Loadings of First Two Components
plot(wavenumbers_NoCO2,
     plsr.fit$loadings[,1],
     type = "l",
     xlab = "wavenumber",
     ylab = "loading",
     main = "Loadings by Component",
     xlim = c(3000, 650),
     ylim = c(-1.5, 2)
)

lines(wavenumbers_NoCO2,
      plsr.fit$loadings[,2],
      col = "red"
      )     

lines(wavenumbers_NoCO2,
      plsr.fit$loadings[,3],
      col = "green"
)     

abline(h = 0)
minor.tick(nx = 5)

legend("top", legend = c("Component 1", "Component 2"), fill = c("black", "red"))
####END####

####Coefficients Plot####
plot(base123.fit[[1]],
     type = "l"
     )

#Choose Number of Components
Num.components <- 3

#wavenumbers_region <- wavenumbers_NoCO2[c(region_2800_2999,region_650_1800)]

#wavenumbers_MW <- wavenumbers_NoCO2[c(Base.Region2,Base.Region3,Base.Region,Base.Region4)]
wavenumbers_MW <- wavenumbers_NoCO2[c(Base.Region1)]#,Base.Region2)]#,Base.Region3)]

plot(wavenumbers_NoCO2,
     plsr.fit$coefficients[,1,3],
     type = "l",
     xlim = c(3600,650),
     xlab = "Wavenumber",
     ylab = "Coefficient",
     main = paste(c("Regression Coefficients by Wavenumber", Num.components, "Components"))
)
abline(h = 0, col = "blue")
minor.tick(nx = 5)

plot(wavenumbers_NoCO2,
     plsr.fit$coefficients[,1,3],
     type = "l",
     xlim = c(1800,800),
     ylim = c(-0.025,0.025),
     xlab = expression("Wavenumber (cm"^-1*")"),
     ylab = "Coefficient",
     main = paste(c("Regression Coefficients by Wavenumber", Num.components, "Components"))
)
abline(h = 0, col = "blue")
minor.tick(nx = 4)

#-or-
par(mar=c(5.1,4.1,4.1,4.1), new = TRUE)

matplot(wavenumbers,
        spectra_avg[,1:40],
        type = "l",
        xlim = c(1800,800),
        col = col_set,
        xlab = expression(Wavenumber ~ "(cm"^{-1}*")"),
        ylab = "Absorbance",
        main = "Overlaid ATR-FTIR Spectra"
)
minor.tick(nx = 4)

legend("topleft", 
       legend = c("Kd > 600", "Kd [300-600]", "Kd < 300"), 
       fill = c("darkgreen", "yellow", "red"), bty = "n"
)

#Superimpose MW Coefficients on Full Spectrum Coefs Plot
par(new = TRUE)

plot(wavenumbers_NoCO2[Base.Region1[1]:Base.Region1[length(Base.Region1)]],
     #wavenumbers_MW[1:49],
     base123.fit$coefficients[,1,4][grep("Base.Region1", names(base123.fit$coefficients[,1,4]))],
     #base123.fit$coefficients[,1,4][1:length(Base.Region1)],
     #base1234.fit$coefficients[,1,3][1:49],
     type = "l",
     col = "red",
     lwd = 2,
     axes = FALSE,
     ylim = c(-0.25,0.25),
     xlim = c(1800,800),
     #xlim = c(1800,650),
     xlab = "",
     ylab = ""
     )

abline(h = 0, col = "blue")

axis(side = 4, 
     at = pretty(range(c(-0.25,0.25))), 
     col.axis = "red",
     col = "red",
     cex = 2
     )
mtext("Coefficients", side = 4, line = 2.5)

x <- wavenumbers_NoCO2[Base.Region1[1]:Base.Region1[length(Base.Region1)]]
y <- base123.fit$coefficients[,1,4][grep("Base.Region1", names(base123.fit$coefficients[,1,4]))] 
x_vert <- 20

polygon(x = c(x[x_vert:length(x)], x[length(x)], x[x_vert]), 
        y = c(y[x_vert:length(x)], 0, 0),
  col = "green")

polygon(x = c(x[1:x_vert], x[x_vert], x[1] ), 
        y = c(y[1:x_vert], 0, 0),
        col = "red")

lines(wavenumbers_NoCO2[Base.Region2[1]:Base.Region2[length(Base.Region2)]],
      base123.fit$coefficients[,1,4][grep("Base.Region2", names(base123.fit$coefficients[,1,4]))],
      col = "red",
      lwd = 2
)

x <- wavenumbers_NoCO2[Base.Region2[1]:Base.Region2[length(Base.Region2)]]
y <- base123.fit$coefficients[,1,4][grep("Base.Region2", names(base123.fit$coefficients[,1,4]))] 
x_vert <- 16

polygon(x = c(x[x_vert:length(x)], x[length(x)], x[x_vert]), 
        y = c(y[x_vert:length(x)], 0, 0),
        col = "red")

polygon(x = c(x[1:x_vert], x[x_vert], x[1] ), 
        y = c(y[1:x_vert], 0, 0),
        col = "green")

lines(wavenumbers_NoCO2[Base.Region3[1]:Base.Region3[length(Base.Region3)]],
      base123.fit$coefficients[,1,4][grep("Base.Region3", names(base123.fit$coefficients[,1,4]))],
      col = "red",
      lwd = 2
      )

x <- wavenumbers_NoCO2[Base.Region3[1]:Base.Region3[length(Base.Region3)]]
y <- base123.fit$coefficients[,1,4][grep("Base.Region3", names(base123.fit$coefficients[,1,4]))] 
x_vert <- 7

polygon(x = c(x[x_vert:length(x)], x[length(x)], x[x_vert]), 
        y = c(y[x_vert:length(x)], 0, 0),
        col = "red")

polygon(x = c(x[1:x_vert], x[x_vert], x[1] ), 
        y = c(y[1:x_vert], 0, 0),
        col = "green")

#Reset plot settings if need be
#dev.off()



#lines(wavenumbers_MW[(53+1):(53+26)],
#      base1234.fit$coefficients[,1,3][(53+1):(53+26)],
#      col = "red",
#      lwd = 2
#      )


#Plot MW Coefficients ONLY
plot(wavenumbers_MW,
     base1.fit$coefficients[,1,4],
     #base1234.fit$coefficients[,1,3][1:49],
     type = "l",
     col = "red",
     lwd = 2,
     #ylim = c(-0.25,0.25),
     xlim = c(1800,650),
     xlab = "",
     ylab = ""
)

abline(h = 0, col = "blue")
minor.tick(nx = 4)

#Store Coefficients as Object
#Full Spectrum
Fit.Coefficients <- plsr.fit$coefficients[,1,3]

#MW
#Fit.Coefficients <- base1234.fit$coefficients[,1,3]
Fit.Coefficients <- base123.fit$coefficients[,1,4]

#Peak-Finding Procedure
peaks.pos <- find_peaks(Fit.Coefficients, m = 2)

peaks.pos <- peaks.pos[which(Fit.Coefficients[peaks.pos]>0)]

peaks.neg <- find_peaks(-Fit.Coefficients, m = 3)

peaks.neg <- peaks.neg[which(Fit.Coefficients[peaks.neg]<0)]

peaks <- sort(c(peaks.pos, peaks.neg), decreasing = TRUE)

#Add lines: Full Spectrum Peaks
#abline(v = wavenumbers_NoCO2[peaks], col = "coral", lty =2)
abline(v = c(824,965,997,1016,1065,1193,1299,1351,1372,1513,1653), col = "coral", lty =2)

#Add lines: MW Peaks
abline(v = wavenumbers_MW[peaks], col = "darkgrey", lty = 2)
abline(v = 1334, col = "darkgrey", lty = 2)
#abline(v = c(1344,1232,1074), col = "coral", lty =2)


#?
#abline(v = 1513, col = "coral", lty =2)
#abline(v = 898, col = "coral", lty =2)
#abline(v = 1432, col = "coral", lty =2)
#abline(v = 1300, col = "coral", lty =2)
#abline(v = 1652, col = "coral", lty =2)
#abline(v = 975, col = "coral", lty =2)

#write.csv(wavenumbers_region[peaks], file = "Regression Coefficient Peaks.csv")
####END####

####Test if Coefficients = Regression Coefficients####
Fit.Coefficients <- plsr.fit$coefficients[,1,3]

coef1 <- coef(plsr.fit, intercept = TRUE, ncomp = 3)[[1]]
coef1 <- coef(base123.fit, intercept = TRUE, ncomp = 4)[[1]]

#(train$MIR[, c(region_2800_2999,region_650_1800)])%*%Fit.Coefficients+coef1 

#(test$MIR)%*%Fit.Coefficients
(test$MIR[, c(Base.Region1,Base.Region2,Base.Region3)])%*%Fit.Coefficients#+coef1 

(sediment$MIR)%*%Fit.Coefficients#+coef1 
(sediment$MIR[,711:1422])%*%Fit.Coefficients[711:1422]
(sediment$MIR[,1:710])%*%Fit.Coefficients[1:710]

predict(base123.fit, ncomp = 4, newdata = test)
###END###

matplot(wavenumbers,
        spectra_avg,
        type = "l",
        xlim = wavenumber,
        col = col_set,
        xlab = "Wavenumber (1/cm)",
        ylab = "Absorbance",
        main = "Overlaid ATR-FTIR Spectra",
        ylim = c(-0.5,1)
)
minor.tick(nx = 5)
legend("top", 
       legend = c("Kd > 600", "Kd > 300", "Kd < 300"), 
       fill = c("red", "yellow", "blue")
)

lines(wavenumbers_region,
      40%*%Fit.Coefficients)

abline(h = 0)

abline(v = wavenumbers_region[peaks], col = "red")
####END####

####"Bootstrap" Intra-Day KD Variation####
#Create Normal Distribution for Each Experimental KD Value
#Based on Mean, SD of slope returned from lm()
Boot <- NULL

set.seed(123)
for(i in 1:nrow(Targets)){
  Boot <- rbind(Boot, t(as.matrix(
    rnorm(n = 100, mean = Targets$Kd[i], sd = Targets$SDEV[i])
  )
  )
  )
}  

Boot_df <- data.frame(row.names = file_names)
Boot_df$Kd <- Boot
Boot_df$Kd <- log10(Boot_df$Kd)
Boot_df$MIR <- MIR
Boot_df$MIR <- Boot_df$MIR[, c(cm3600_2200, cm1900_650)]
Boot_df <-Boot_df[-((nrow(Boot_df)-5):nrow(Boot_df)),] 

#Fit N Models using randomly-selected values from each KD distribution
Boot_Coefs <- NULL
Boot_Interc <- NULL

set.seed(123)
for(i in 1:ncol(Boot)){
  temp <- plsr(Kd[,i] ~ MIR,
              ncomp = 3,
              data = Boot_df,
              scale = FALSE,
              center = TRUE,
              method = "simpls"
              )
  Boot_Coefs <- cbind(Boot_Coefs, temp$coefficients[,,3])
  Boot_Interc <- rbind(Boot_Interc, coef(temp, intercept = TRUE, ncomp = 3)[[1]] )
}

Boot_PredKD <- data.frame(file_names[1:40],
  Mean = NA,
  SDEV = NA,
  row.names = 1
)
temp <- NULL

#Calculated sample KDs from each bootstrap model, average them
for(i in 1:nrow(Boot_df)){
  for(j in 1:nrow(Boot_Interc)){
  value <- Boot_df[,"MIR"]%*%Boot_Coefs[,j] + Boot_Interc[j]
  temp <- cbind(temp, value)
  }
  Boot_PredKD[i,1] <- mean(temp[i,])
  Boot_PredKD[i,2] <- sd(temp[i,])
  }

plot(Boot_PredKD$Mean,
     spectra.df$logKd[1:40],
     line = TRUE,
     xlim = c(1,4),
     ylim = c(1,4)
     )

##Careful! Long runtime:
#matplot(wavenumbers_NoCO2,
#        Boot_Coefs,
#        xlim = wavenumber
#        )

####END####

#####"Bootstrap" Inter-Day KD Variation####

####Other####
#plsr 'projection' = r weights
#test code on plsr.fit
X0 <- spectra.df_NoSed[1,3]-plsr.fit$Xmeans
t.score <- X0[1,]%*%plsr.fit$projection[,1]
round(t.score, digits = 10) == round(plsr.fit$scores[1,1], digits = 10)

#r proportional to S for univariate y
#Calculate cross product .: = R
S1 <- t(spectra.df_NoSed[,3])%*%(spectra.df_NoSed[,2]-plsr.fit$Ymeans)
#In plsr object, R has been divided by normt. Solve for normt using relation
#normt = r/r_updated
normt <- S1/plsr.fit$projection[,1]

round(S1/normt, digits = 10) == round(plsr.fit$projection[,1], digits = 10)

#Work backwards and see if we can arrive at the score of the first observation
#starting with S matrix = r -> t = X0*r -> t = t/normt 

#reference t score:
t.score <- as.matrix(plsr.fit$scores[,1])

round(((spectra.df_NoSed[1,3]-plsr.fit$Xmeans)%*%S1)/normt[1,1], digits = 10) == 
  round(t.score[1,1], digits = 10)


FeS <- read.csv("JAN10_FeSDirect_DRIFT.CSV", col.names = c("X", "Y") )
#Transform Transmittance to Absorbance
FeS$Y <- 2-log10(FeS$Y)

#Baseline Correct Absorbance 
FeS$Y_corr <- t((baseline.modpolyfit(t(FeS$Y), degree = 4))$corrected) 

plot(FeS$X,
     FeS$Y,
     type = "l"
)
abline(v = c(649,1156, 827,1538,1637))

plot(FeS$X,
     FeS$Y_corr,
     type = "l"
     )  

x <- c(0.11, 0.31, 0.43, 0.91, 7.87, 9.95)
y <- c(6.72, 7.12, 8.10, 9.35, 9.96, 10.02)
####END####