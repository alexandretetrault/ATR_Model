#Load Packages
library(baseline)
library(signal)
library(pls)
library(caTools)
library(dplyr)
library(Hmisc)
library(HotellingEllipse)
#library(RColorBrewer)
library(caret)

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

####Input Matrix-Matched Salts####
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

#!STOP! -> go to Matrix_Matched.R
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
limits <- rev(range(salt_spectra[[1]][1]))
wavenumbers <- as.vector(salt_spectra[[1]][1])

lapply(names(salt_spectra), function(x) {
  plot(salt_spectra[[x]]$Wavenumber, 
       #salt_spectra[[x]]$Normalized, 
       salt_spectra[[x]]$Transmittance,
       type = "l", 
       xlim = limits,
       main = x
  )
}
)

#Average the Normalized Absorbance Spectra
salt_avg <- matrix(nrow = length(salt_spectra[[1]]$Wavenumber))

#Create String of Compounds
Salt_Compounds <- print(paste(gsub("\\(.*", "", salt_file_names)))

salt_file_names1 <- NULL 

for(i in seq(0, length(salt_spectra)-3, 3)){
  avg <- as.vector(
    (salt_spectra[[i+1]][4] + 
       salt_spectra[[i+2]][4] + 
       salt_spectra[[i+3]][4])/3)
  salt_avg <- cbind(salt_avg, avg)
  salt_file_names1 <- c(salt_file_names1, Salt_Compounds[i +1])
  
}
salt_avg <- data.frame(salt_avg[ , -1])
colnames(salt_avg) <- salt_file_names1

#Plot Averaged Spectra
for (i in 1:ncol(salt_avg)){
  plot(wavenumbers[, 1],
       salt_avg[ , i],
       type = "l",
       xlim = wavenumber,
       xlab = "Wavenumber",
       ylab = "Absorbance",
       main = names(salt_avg[i])
  )
  minor.tick(nx = 5)
}
rm(salt_file_names1)
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

#!STOP! Go to -> Matrix_Matched.R
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
       x$Absorbance,
       #x$Normalized, 
       type = "l", 
       xlim = wavenumber
       )
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
####END####

####Correct LMW Spectra Using Matrix-Matched####
#Make a New Object for Corrected Spectra
spectra_avg_corr <- spectra_avg

#Scale Salt_Avg Max Peaks to Corresponding Peaks in Spectra_Avg
for(i in 1:ncol(salt_avg)){
  temp_num <- which(salt_avg[[i]] == max(salt_avg[[i]]))
  temp_salt <- names(salt_avg)[i]
  temp_comp <- paste(gsub("\\_salt*", "", temp_salt))
  #Scale to Corresponding Spectra_Avg Peak
  salt_avg[[i]] <- salt_avg[[i]]*spectra_avg[[temp_comp]][temp_num]
}

#Plot Salt Matrix vs. Averaged Spectra
for(i in names(spectra_avg)){
  for(j in names(salt_avg)){
    if(paste(gsub("\\_salt*", "", j)) == i){
      plot(wavenumbers[,1],
           spectra_avg[[i]],
           type = "l",
           xlim = wavenumber,
           main = i)
      lines(wavenumbers[,1],
            salt_avg[[j]],
            type = "l",
            col = "red")
      minor.tick(nx = 5)
    }
  }
}

##If Salt Spectra Name Matches Spectra_Avg Name, Subtract the Former Absorbance from the Latter
for(i in names(spectra_avg)){
  for(j in names(salt_avg)){
    if(paste(gsub("\\_salt*", "", j)) == i){
      spectra_avg_corr[[i]] <- spectra_avg[[i]] - salt_avg[[j]]
    }
  }
}

#Change Negative Absorbance No.'s to Zero
for(i in grep("LMW", names(spectra_avg_corr))){
  for(j in 1:nrow(spectra_avg_corr)){
    if(spectra_avg_corr[[i]][j] <= 0){
      spectra_avg_corr[[i]][j] <- 0
    }
  }
}

#Change Negative Absorbance No.'s to Zero
for(i in grep("Sediment", names(spectra_avg_corr))){
  for(j in 1:nrow(spectra_avg_corr)){
    if(spectra_avg_corr[[i]][j] <= 0){
      spectra_avg_corr[[i]][j] <- 0
    }
  }
}

#Min-Max Normalize Again
for(i in grep("LMW", names(spectra_avg_corr))){
  spectra_avg_corr[[i]] <- scale(spectra_avg_corr[[i]],
                            center = FALSE,
                            scale = max(spectra_avg_corr[[i]]) -
                              min(spectra_avg_corr[[i]])
  )
}

#Plot Post-Subtraction
for (i in 1:ncol(spectra_avg_corr)){
  plot(wavenumbers[, 1],
       spectra_avg_corr[ , i],
       type = "l",
       xlim = wavenumber,
       xlab = "Wavenumber",
       ylab = "Absorbance",
       main = names(spectra_avg_corr[i])
  )
  minor.tick(nx = 5)
}
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

####Add Compounds + Kd values to Target Variables.csv if necessary####
#Write the file names of compounds to csv
write.csv(file_names, file = "File Names.csv")


## -> Manually copy file names to Target Variables.csv and add in Kd data 

#Load Target Variables Data
Targets <- read.csv("Target Variables.csv", sep = ",", row.names = 1)
####END####

####Dataframe Preparation for PLS - Full Spectrum####
spectra.df <- data.frame(row.names = names(spectra_avg_corr))

spectra.df$Kd <- Targets$Kd

#Log-transform KD Values
spectra.df$logKd <- log10(Targets$Kd)
####END####

####Make Nested Matrix of MIR Spectra for All Compounds####
MIR <- as.matrix(spectra_avg_corr)
MIR <- t(MIR)
spectra.df$MIR <- MIR

#Smoothing
MIR.smooth <- apply(spectra.df$MIR, 
                    MARGIN = 1, 
                    FUN = sgolayfilt, 
                    p = 1, n = 15, m = 0, ts = 1
)

spectra.df$MIR.smooth <- t(MIR.smooth)

#First Derivative
MIR.1deriv <- apply(spectra.df$MIR.smooth, 
                   MARGIN = 1, 
                   FUN = sgolayfilt, 
                   p = 1, n = 15, m = 1, ts = 1
)

spectra.df$MIR.1deriv <- t(MIR.1deriv)
####END####

####Remove Diamond Region/Fingerprint & Split Data for Training on Spectra####
#Make the Wavenumbers Scale
wavenumbers <- as.vector(spectra[[1]][1])

#Spectral Regions of Interest
cm3600_2200 <- which(between(wavenumbers[[1]], 2200, 3600) == TRUE)
cm1900_900 <- which(between(wavenumbers[[1]], 900, 1900) == TRUE)
wavenumbers_NoCO2 <- wavenumbers[c(cm3600_2200, cm1900_900), 1]

limits <- rev(range(wavenumbers_NoCO2))

spectra.df_NoSed <- spectra.df

spectra.df_NoSed$MIR <- spectra.df$MIR[, c(cm3600_2200, cm1900_900)]
spectra.df_NoSed$MIR.smooth <- spectra.df$MIR.smooth[, c(cm3600_2200, cm1900_900)]
spectra.df_NoSed$MIR.1deriv <- spectra.df$MIR.1deriv[, c(cm3600_2200, cm1900_900)]

#Make sediment test dataframe
sediment <- spectra.df_NoSed[(nrow(spectra.df_NoSed)-5):nrow(spectra.df_NoSed),]

#Remove sediment from train dataframe
spectra.df_NoSed <- spectra.df_NoSed[-((nrow(spectra.df_NoSed)-5):nrow(spectra.df_NoSed)),]

#Remove LMW Compounds
#spectra.df_NoSed <- spectra.df_NoSed[-c(4:6,15,16,19,25,31,34:40),]

#Split between train and test df

#Which Samples are Endmembers?
#endmem <- grep("^.*_.*?([1]).*", row.names(spectra.df_NoSed))

#set.seed(101)
set.seed(123)
#split <- sample.split(spectra.df_NoSed[-endmem,1], SplitRatio = 18)
split <- sample.split(spectra.df_NoSed$Kd, SplitRatio = 0.75)

#train <- rbind(spectra.df_NoSed[endmem,], subset(spectra.df_NoSed[-endmem,], split == T)) 
train <- subset(spectra.df_NoSed, split == T)

#test <- subset(spectra.df_NoSed[-endmem, ], split == F)
test <- subset(spectra.df_NoSed, split == F)

#Number rows of train dataframe/spectra.df 
train <- cbind(train, n = 1:nrow(train))
test <- cbind(test, n = 1:nrow(test))
spectra.df_NoSed <- cbind(spectra.df_NoSed, n = 1:nrow(spectra.df_NoSed))
####END####

####Plot the ATR Spectra####
wavenumbers <- as.vector(spectra[[1]][1])

#Plot smoothed Spectra of All Compounds
for (i in rownames(spectra.df)){
  plot(wavenumbers[,1],
       spectra.df$MIR.smooth[i,],
       type = "l",
       xlim = c(3600, 900),
       main = i
  )
  minor.tick(nx = 5)
}  

#Plot 1st-Deriv Spectra of All Compounds
for (i in rownames(spectra.df)){
  plot(wavenumbers[,1],
       spectra.df$MIR.1deriv[i,],
       type = "l",
       xlim = c(3600, 800),
       main = i
  )
}  

#Or Plot as Transmittance
lapply(names(spectra), function(x) {
  plot(spectra[[x]]$Wavenumber, 
       spectra[[x]]$Transmittance, 
       type = "l", 
       xlim = wavenumber,
       main = x
  )
}
)

for (i in 1:ncol(spectra_avg_corr)){
  plot(wavenumbers[, 1],
       spectra_avg[ , i],
       type = "l",
       xlim = wavenumber,
       xlab = "Wavenumber",
       ylab = "Absorbance",
       main = names(spectra_avg_corr[i])
  )
  minor.tick(nx = 5)
}

col_set <- brewer.pal(6, "RdYlGn") 
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

####RUN MODEL on Full Spectrum####
set.seed(123)
plsr.fit <- plsr(logKd ~ MIR.smooth,
                     ncomp = 10,
                     #data = spectra.df_NoSed,
                     data = train,
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
                          comp3 = plsr.fit$scores[,3],
                          comp4 = plsr.fit$scores[,4],
                          comp5 = plsr.fit$scores[,5],
                          comp6 = plsr.fit$scores[,6],
                          comp7 = plsr.fit$scores[,7],
                          comp8 = plsr.fit$scores[,8],
                          comp9 = plsr.fit$scores[,9]
)

plsr.scores <- plsr.scores %>%
  as_tibble() %>%
  print()

T2.ellipse <- ellipseParam(data = plsr.scores, k = 5, pcx = 1, pcy = 2)

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

pls::R2(plsr.fit)

predplot(plsr.fit, 
         ncomp = 9, 
         line = TRUE,
         main = "Model Validation",
         xlim = c(1,3.5),
         ylim = c(1,3.5),
         pch = 21,
         bg = "gray"
         )
####END####

####Plot the Scores & Loadings####
plot(plsr.fit, plottype = "scores", comps = 1:5, labels = "numbers")

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

loadingplot(plsr.fit, comps = 1:9, legendpos = "topleft")

matplot(wavenumbers_NoCO2, 
        plsr.fit$loadings[,1:9],
        type = "l",
        xlim = limits,
        col = 1:9,
        lty = 1:5
        )
minor.tick(nx = 5)

legend(x = "topleft",
       legend = 1:9,
       col = 1:9,
       lty = 1:5
       )

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

print(alex.RMSEP(model = plsr.fit, data = test1, target = "logKd"))
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

##Fit Plots##

####Full Spectra####
pls::R2(plsr.fit)

plot(train$logKd,
     predict(plsr.fit, ncomp = 9),
     xlab = "measured",
     ylab = "predicted",
     main = "Model Validation",
     xlim = c(1, 3.5),
     ylim = c(1, 3.5),
     pch = 21,
     bg = "gray"
)
abline(a = 0, b = 1)
minor.tick(nx = 5)
#legend("topleft", legend = "Train", fill = "gray")


RMSEP(plsr.fit, newdata = test)
print(alex.RMSEP(plsr.fit, data = test, target = "logKd"))

points(test1$logKd,
       predict(plsr.fit, ncomp = 9, newdata = test1),
       pch = 21,
       bg = "blue"
)
points(sediment[,2],
       predict(plsr.fit, ncomp = 9, newdata = sediment),
       pch = 21,
       bg = "coral"
)

#legend("topleft", 
#       legend = c("Train", "Test", "Sediment"), 
#       fill = c("gray","blue", "coral")
#)

#RMSEP(plsr.fit, newdata = sediment)
#print(alex.RMSEP(plsr.fit, data = sediment, target = "logKd"))
####END####

####Best MW Fit####
predplot(base123.fit,
         ncomp = 9,
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

points(test1[,2],
       predict(base123.fit, ncomp = 6, newdata = test1),
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
     y = plsr.fit[9]$fitted.value[,1,9], 
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
#plot(base123.fit[[1]],
#     type = "l"
#     )

#Choose Number of Components
Num.components <- 9

#wavenumbers_MW <- wavenumbers_NoCO2[c(Base.Region1 ,Base.Region2, Base.Region3)]

plot(wavenumbers_NoCO2,
     plsr.fit$coefficients[,1,Num.components],
     type = "l",
     xlim = c(3600,900),
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
        spectra_avg_corr[,1:40],
        type = "l",
        xlim = c(3600,900),
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
####END####

####Peak Analysis of Regression Coefficients####
#Full Spectrum
Fit.Coefficients <- plsr.fit$coefficients[,1,9]

#MW
#Fit.Coefficients <- base1234.fit$coefficients[,1,3]
#Fit.Coefficients <- base123.fit$coefficients[,1,4]

#Peak-Finding Procedure
peaks.pos <- find_peaks(Fit.Coefficients, m = 3)

peaks.pos <- peaks.pos[which(Fit.Coefficients[peaks.pos]>0)]

peaks.neg <- find_peaks(-Fit.Coefficients, m = 3)

peaks.neg <- peaks.neg[which(Fit.Coefficients[peaks.neg]<0)]

peaks <- sort(c(peaks.pos, peaks.neg), decreasing = TRUE)

#Add lines: Full Spectrum Peaks
#abline(v = wavenumbers_NoCO2[peaks], col = "coral", lty =2)
abline(v = c(3386,3146,2922,1720,1657,1577,1435,1388,1288,1232,1159,1109,1023,984), 
       col = "coral", lty = 2)

#Add lines: MW Peaks
#abline(v = wavenumbers_MW[peaks], col = "darkgrey", lty = 2)
#abline(v = 1334, col = "darkgrey", lty = 2)
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
Fit.Coefficients <- plsr.fit$coefficients[,1,9]

coef1 <- coef(plsr.fit, intercept = TRUE, ncomp = 9)[[1]]
#coef1 <- coef(base123.fit, intercept = TRUE, ncomp = 4)[[1]]

#(train$MIR[, c(region_2800_2999,region_650_1800)])%*%Fit.Coefficients+coef1 

(train$MIR.smooth)%*%Fit.Coefficients# + coef1
(test$MIR.smooth)%*%Fit.Coefficients# + coef1
#(test$MIR[, c(Base.Region1,Base.Region2,Base.Region3)])%*%Fit.Coefficients#+coef1 

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