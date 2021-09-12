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
       salt_spectra[[x]]$Normalized, 
       #salt_spectra[[x]]$Absorbance,
       type = "l", 
       xlim = wavenumber,
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
####END####

####Correct LMW Spectra Using Matrix-Matched####

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

##If Salt Spectra Name Matches Spectra_Avg Name, Substract the Former Absorbance from the Latter
for(i in names(spectra_avg)){
  for(j in names(salt_avg)){
    if(paste(gsub("\\_salt*", "", j)) == i){
      spectra_avg[[i]] <- spectra_avg[[i]] - salt_avg[[j]]
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

#Change Negative Absorbance No.'s to Zero
for(i in grep("Sediment", names(spectra_avg))){
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

#Plot Post-Subtraction
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
####END####
