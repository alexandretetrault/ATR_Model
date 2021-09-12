####Moving Window Partial Least Squares####
##Set Total Number of Spectral Points (i.e w/o CO2)
total.points <- length(train$MIR.smooth[1,])
#total.points <- length(spectra.df_NoSed$MIR[1,])

#Create an object to store the results
residues <- NULL

for(i in 1:(total.points - 9)){
  #Set the midpoint of the moving window of 10 spectral points
  midpoint <- wavenumbers_NoCO2[i+4]
  
  #Train the Model on Scaled Absorbance Data
  set.seed(123)
  mwindow.fit <- plsr(logKd ~ MIR.smooth[ ,i:(i+9)],
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
Num.components <- 9

#Window 1
#HMW Only
window.1 <- which(between(wavenumbers_NoCO2, 900, 1100) == TRUE)
#BEST
window.1 <- which(between(wavenumbers_NoCO2, 900, 999) == TRUE)
#TRIAL
window.1 <- which(between(wavenumbers_NoCO2, 900, 1100) == TRUE)

#Window 2
#HMW Only
window.2 <- which(between(wavenumbers_NoCO2, 1150, 1400 ) == TRUE)
#BEST
window.2 <- which(between(wavenumbers_NoCO2, 1000, 1250) == TRUE)
#TRIAL
window.2 <- which(between(wavenumbers_NoCO2, 1120, 1220) == TRUE)

#Window 3
#HMW Only
window.3 <- which(between(wavenumbers_NoCO2, 1500, 1700) == TRUE)
#BEST
window.3 <- which(between(wavenumbers_NoCO2, 1250, 1450) == TRUE)
#TRIAL
window.3 <- which(between(wavenumbers_NoCO2, 1700, 1900) == TRUE)

#Window 4
#Define Spectral Region 4
window.4 <- which(between(wavenumbers_NoCO2, 1300, 1425) == TRUE)
#Trial
window.4 <- which(between(wavenumbers_NoCO2, 2800, 3000) == TRUE)

#Window 5
#Define Spectral Region 5
window.5 <- which(between(wavenumbers_NoCO2, 1425, 1475) == TRUE)
#Trial
window.5 <- which(between(wavenumbers_NoCO2, 3250, 3350) == TRUE)

#Fit PLS
#Train the Model on Windows of Absorbance Data
set.seed(123)
window1.fit <- plsr(logKd ~ MIR.smooth[ , window.1],
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
window2.fit <- plsr(logKd ~ MIR.smooth[ , window.2],
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
                    validation = "CV",
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
                    validation = "none",
                    #validation = "CV",
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
                    validation = "none",
                    #validation = "CV",
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
Num.components <- 9

##CSMWPLS on Base Region (Window _)

#Select window
#BEST
window <- window.3
#TRIAL
window <- window.4

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
                        MIR.smooth[ , window[j]:window[j + w - 1]],
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
Base.Region1 <- which(between(wavenumbers_NoCO2, 2838, 2984) == TRUE)

#Run PLSR on Base Region and Append for General Comparison
set.seed(123)
base1.fit <- plsr(logKd ~ MIR.smooth[ , Base.Region1],
                  ncomp = 9, 
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
#BEST
window <- window.1 
#TRIAL
window <- window.5 

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
                         MIR.smooth[ , window[j]:window[j + w - 1]] +
                         MIR.smooth[ , Base.Region1],
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
Base.Region2 <- which(between(wavenumbers_NoCO2, 3259, 3277) == TRUE)


#Run PLSR on Base Region 1 & 2 and Append for General Comparison
set.seed(123)
base12.fit <- plsr(logKd ~ MIR.smooth[ , Base.Region1] + 
                     MIR.smooth[ , Base.Region2],
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
window <- window.3

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
                          MIR.smooth[ , Base.Region1] +
                          MIR.smooth[ , Base.Region2] +
                          MIR.smooth[ , window[j]:
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
Base.Region3 <- which(between(wavenumbers_NoCO2, 1759, 1784) == TRUE)

#Run PLSR on Base Region 123 and Append for General Comparison
set.seed(123)
base123.fit <- plsr(logKd ~ 
                      MIR.smooth[ , Base.Region1] + 
                      MIR.smooth[ , Base.Region2] +
                      MIR.smooth[ , Base.Region3],
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
for(w in 10:length(window)){
  
  #Factor region by window size to determine no. windows that span it
  span <- length(window) - w + 1
  
  for(j in 1:span){
    #Train the Model on Scaled Absorbance Data
    set.seed(123)
    base1234.fit <- plsr(logKd ~
                           MIR.smooth[ , window[j]:window[j + w - 1]] +
                           MIR.smooth[ , Base.Region2] +
                           MIR.smooth[ , Base.Region1] +
                           MIR.smooth[ , Base.Region3],
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

Base.Region4 <- which(between(wavenumbers_NoCO2, 1041, 1059) == TRUE)

#Run PLSR on Base Region 1234 and Append for General Comparison
set.seed(123)
base1234.fit <- plsr(logKd ~
                       MIR.smooth[ , Base.Region2] +
                       MIR.smooth[ , Base.Region3] + 
                       MIR.smooth[ , Base.Region1] + 
                       MIR.smooth[ , Base.Region4],
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

RMSEP(base1234.fit, estimate = "CV")

#Strip the "MIR[]" off of the returned dimension names for easier plotting later
#dimnames(Base.Region1.2.fit[[3]])[[1]] <- gsub(".*]", "", dimnames(Base.Region1.2.fit[[3]])[[1]])

RMSEP(base1234.fit)#, estimate = "train")
RMSEP(base123.fit)#, estimate = "train")
RMSEP(base12.fit)#, estimate = "train")
RMSEP(base1.fit)#, estimate = "train")

print(alex.RMSEP(base1.fit, data = test1, target = "logKd"))
print(alex.RMSEP(base12.fit, data = test1, target = "logKd"))
print(alex.RMSEP(base123.fit, data = test1, target = "logKd"))
print(alex.RMSEP(base1234.fit, data = test1, target = "logKd"))
print(alex.RMSEP(plsr.fit, data = test1, target = "logKd"))
####END####