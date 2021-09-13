####Average the Salt Transmittance Spectra####
salt_avg <- matrix(nrow = length(salt_spectra[[1]]$Wavenumber))

#Create String of Compounds
Salt_Compounds <- print(paste(gsub("\\(.*", "", salt_file_names)))

salt_file_names1 <- NULL 

for(i in seq(0, length(salt_spectra)-3, 3)){
  avg <- as.vector(
    (salt_spectra[[i+1]][2] + 
       salt_spectra[[i+2]][2] + 
       salt_spectra[[i+3]][2])/3)
  salt_avg <- cbind(salt_avg, avg)
  salt_file_names1 <- c(salt_file_names1, Salt_Compounds[i +1])
  
}
salt_avg <- data.frame(salt_avg[ , -1])
colnames(salt_avg) <- salt_file_names1


wavenumbers <- as.vector(spectra[[1]][1])
limits <- rev(range(salt_spectra[[1]][1]))

#Plot Averaged Spectra
for (i in 1:ncol(salt_avg)){
  plot(wavenumbers[, 1],
       salt_avg[ , i],
       type = "l",
       xlim = limits,
       xlab = "Wavenumber",
       ylab = "Absorbance",
       main = names(salt_avg[i])
  )
  minor.tick(nx = 5)
}
rm(salt_file_names1)
####END####

####Average the Compound Transmittance Spectra####
spectra_avg <- matrix(nrow = length(spectra[[1]]$Wavenumber))

file_names <- NULL 

for(i in seq(0, length(spectra)-5, 5)){
  avg <- as.vector(
    (spectra[[i+1]][2] + 
       spectra[[i+2]][2] + 
       spectra[[i+3]][2] +
       spectra[[i+4]][2] +
       spectra[[i+5]][2])/5)
  spectra_avg <- cbind(spectra_avg, avg)
  file_names <- c(file_names, Compounds[i +1])
  
}
spectra_avg <- data.frame(spectra_avg[ , -1])
colnames(spectra_avg) <- file_names

#Plot Averaged Spectra
for (i in 1:ncol(spectra_avg)){
  plot(wavenumbers[, 1],
       spectra_avg[ , i],
       type = "l",
       xlim = limits,
       xlab = "Wavenumber",
       ylab = "Absorbance",
       main = names(spectra_avg[i])
  )
  minor.tick(nx = 5)
}
####END####

####Transform Transmittance to Absorbance####
salt_avg <- 2-log10(salt_avg)

spectra_avg <- 2-log10(spectra_avg)
####END####

####Anchor Paired Spectra at 2160 cm^-1####
for(i in 1:ncol(salt_avg)){
  temp_salt <- names(salt_avg)[i]
  temp_comp <- paste(gsub("\\_salt*", "", temp_salt))
  if(salt_avg[[i]][987] > spectra_avg[[temp_comp]][987]){
    temp_diff <- salt_avg[[i]][987] - spectra_avg[[temp_comp]][987]
    salt_avg[[i]] <- salt_avg[[i]] - temp_diff
    rm(temp_diff)
  }
}
####END####

####Plot Salt Matrix vs. Averaged Spectra####
for(i in names(spectra_avg)){
  for(j in names(salt_avg)){
    if(paste(gsub("\\_salt*", "", j)) == i){
      plot(wavenumbers[,1],
           spectra_avg[[i]],
           type = "l",
           xlim = limits,
           ylim = c(0, 1),
           main = i)
      lines(wavenumbers[,1],
            salt_avg[[j]],
            type = "l",
            col = "red")
      minor.tick(nx = 5)
    }
  }
}
####END####

####Write Averaged Spectra to File of Offline Scatter Correction####
for(i in names(salt_avg)){
  write.csv(cbind(wavenumbers[[1]],salt_avg[[i]]),
            file = paste0("/media/alex/KINGSTON/temp/",i,".csv"),
            row.names = FALSE)
}

for(i in names(spectra_avg)){
  write.csv(cbind(wavenumbers[[1]],spectra_avg[[i]]),
            file = paste0("/media/alex/KINGSTON/temp/",i,".csv"),
            row.names = FALSE)
}
####END####

####Input Offline Corrected Spectra####
salt_base <- salt_avg
spectra_base <- spectra_avg

for(i in colnames(salt_base)){
  temp <- read.csv2(paste0("./Salts_Corr/",i,".csv"),
                                skip = 2,
                                nrows = 1797,
                                colClasses = c("NULL", "NULL", "numeric"),
                                header = FALSE
  )
  salt_base[[i]] <- rev(temp[,1])
  rm(temp)
}

for(i in colnames(spectra_base)){
  temp <- read.csv2(paste0("./Spectra_Corr/",i,".csv"),
                    skip = 2,
                    nrows = 1797,
                    colClasses = c("NULL", "NULL", "numeric"),
                    header = FALSE
  )
  spectra_base[[i]] <- rev(temp[,1])
  rm(temp)
}
####END####

####Min-Max Normalize the Spectra####
for(i in names(salt_base)){
  salt_base[[i]] <- as.vector(
    scale(salt_base[[i]],
          center = FALSE,
          scale = max(salt_base[[i]] -
                        min(salt_base[[i]])
          )
    )
  )
}


for(i in names(spectra_base)){
  spectra_base[[i]] <- as.vector(
    scale(spectra_base[[i]],
          center = FALSE,
          scale = max(spectra_base[[i]] -
                        min(spectra_base[[i]])
          )
    )
  )
}
####END####

####Plot Overlaid Corrected Spectra####
for(i in names(spectra_base)){
  for(j in names(salt_base)){
    if(paste(gsub("\\_salt*", "", j)) == i){
      plot(wavenumbers[,1],
           spectra_base[[i]],
           type = "l",
           xlim = limits,
           ylim = c(-0.2,1),
           main = i)
      lines(wavenumbers[,1],
            salt_base[[j]],
            type = "l",
            col = "red")
      minor.tick(nx = 5)
    }
  }
}
####END####

####May 27 LMW Salt Special Adjustment####
which(between(wavenumbers[[1]], 1342, 1500))
which.min(spectra_base[["May27_TPC_0-08-02_LMW"]][1181:1206])
which.max(salt_base[["May27_TPC_0-08-02_LMW_salt"]][1342:1425])

#Good
salt_base[["May27_TPC_0-08-02_LMW_salt"]][1181:1341] <-
  salt_base[["May27_TPC_0-08-02_LMW_salt"]][1181:1341] -
  (salt_base[["May27_TPC_0-08-02_LMW_salt"]][1192] -
     spectra_base[["May27_TPC_0-08-02_LMW"]][1192])
####END####

####April 09 LMW Salt Special Adjustment####
which(between(wavenumbers[[1]], 1100, 1200))
which.min(spectra_base[["May27_TPC_0-08-02_LMW"]][1181:1206])
which.min(salt_base[["Apr09_TPC_0-1-0_LMW_salt"]][1503:1555])

#Good
salt_base[["Apr09_TPC_0-1-0_LMW_salt"]][1288:1341] <-
  salt_base[["Apr09_TPC_0-1-0_LMW_salt"]][1288:1341] +
  (spectra_base[["Apr09_TPC_0-1-0_LMW"]][1341] -
     salt_base[["Apr09_TPC_0-1-0_LMW_salt"]][1341]) 

#Good
salt_base[["Apr09_TPC_0-1-0_LMW_salt"]][1503:1663] <-
  salt_base[["Apr09_TPC_0-1-0_LMW_salt"]][1503:1663] -
    (salt_base[["Apr09_TPC_0-1-0_LMW_salt"]][1516] -
      spectra_base[["Apr09_TPC_0-1-0_LMW"]][1516]) 
####END####     

plot(wavenumbers[[1]],
      spectra_base[["Apr09_TPC_0-1-0_LMW"]],
     type = "l",
     xlim = limits,
     ylim = c(-0.2,1))
minor.tick(nx = 5)
      

lines(wavenumbers[[1]],
      salt_base[["Apr09_TPC_0-1-0_LMW_salt"]],
      col = "blue")

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




####Smoothing####
for(i in 1:ncol(salt_avg_abs)){
  salt_avg_abs[[i]] <- sgolayfilt(salt_avg_abs[[i]],
                                  p = 1,
                                  n = 15,
                                  m = 0,
                                  ts = 1
  )
  
}

for(i in 1:ncol(spectra_avg_abs)){
  spectra_avg_abs[[i]] <- sgolayfilt(spectra_avg_abs[[i]],
                                     p = 1,
                                     n = 15,
                                     m = 0,
                                     ts = 1
  )
  
}
####END####

####Code Dump####

plot(wavenumbers[,1],
     salt_base[[7]],
     type = "l",
     xlim = wavenumber,
     ylim = c(-0.1,0.5))
lines(wavenumbers[,1], salt_avg_abs[[7]], col = "red")


for(i in 1:ncol(salt_base)){
  temp <- t(
    baseline.modpolyfit(
      t(salt_base[[i]][1086:1797]),
      degree = 3
    )[[2]]
  )
  temp <- as.vector(temp)
  salt_base[[i]][1086:1797] <- temp
  rm(temp)
}

for(i in 1:ncol(spectra_base)){
  spectra_base[[i]][1086:1797] <- t(
    baseline.modpolyfit(
      t(spectra_base[[i]][1086:1797]),
      degree = 3
    )[[2]]
  )
}

#Find Index of a given wavenumber/range of wavenumbers
which(between(wavenumbers[[1]], 800,850) == 
        max(salt_base[["May26_TPC_08-0-02_LMW_salt"]][between(wavenumbers[[1]], 800,850)]))

max(salt_base[["May26_TPC_08-0-02_LMW_salt"]][between(wavenumbers[[1]], 800,850)])


#Find index of max absorbance value in that range
(salt_base[["Apr08_TPC_1-0-0_LMW_salt"]][859:885] == max(salt_base[["Apr08_TPC_1-0-0_LMW_salt"]][859:885]))

#Place line at the
abline(v = (wavenumbers[[1]][861]))
####END####



#Plot Averaged Spectra
for (i in 1:ncol(salt_avg)){
  plot(wavenumbers[, 1],
       salt_avg[ , i],
       type = "l",
       xlim = limits,
       xlab = "Wavenumber",
       ylab = "Absorbance",
       main = names(salt_avg[i])
  )
  minor.tick(nx = 5)
}


####Baseline Correct Absorbance#### 
#Create Object
salt_base <- salt_avg
spectra_base <- spectra_avg

####RollingBall####
for(i in 1:ncol(salt_base)){
  temp <- t(
    baseline.rollingBall(
      t(salt_base[[i]][1086:1797]), wm = 300, ws = 300
    )[[2]]
  )
  temp <- as.vector(temp)
  salt_base[[i]][1086:1797] <- temp
  rm(temp)
}

for(i in 1:ncol(spectra_base)){
  temp <- t(
    baseline.rollingBall(
      t(spectra_base[[i]][1086:1797]), wm = 300, ws = 300
    )[[2]]
  )
  temp <- as.vector(temp)
  spectra_base[[i]][1086:1797] <- temp
  rm(temp)
}
####END####
####ALS####
for(i in 1:ncol(salt_base)){
  temp <- t(
    baseline.als(
      t(salt_base[[i]][1086:1797]), lambda = 6
    )[[2]]
  )
  temp <- as.vector(temp)
  salt_base[[i]][1086:1797] <- temp
  rm(temp)
}

for(i in 1:ncol(spectra_base)){
  temp <- t(
    baseline.als(
      t(spectra_base[[i]][1086:1797]), lambda = 6
    )[[2]]
  )
  temp <- as.vector(temp)
  spectra_base[[i]][1086:1797] <- temp
  rm(temp)
}
####END####


####ModPolyFit####
for(i in 1:ncol(salt_base)){
  temp <- t(
    baseline.modpolyfit(
      t(salt_base[[i]][1086:1797]), degree = 4
    )[[2]]
  )
  temp <- as.vector(temp)
  salt_base[[i]][1086:1797] <- temp
  rm(temp)
}

for(i in 1:ncol(spectra_base)){
  temp <- t(
    baseline.modpolyfit(
      t(spectra_base[[i]][1086:1797]), degree = 4
    )[[2]]
  )
  temp <- as.vector(temp)
  spectra_base[[i]][1086:1797] <- temp
  rm(temp)
}
####END####

