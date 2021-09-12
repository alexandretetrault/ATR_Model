
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
     predict(plsr.fit, ncomp = 5, newdata = sediment),
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

####Repeated k-Fold####
set.seed(123)
train_control <- trainControl(method = "repeatedcv", 
                              number = 5, 
                              repeats = 5
)

set.seed(123)
rpeat.kfold.fit <- train(form = logKd ~
                           MIR.smooth,
                         method = "simpls",
                         trControl = train_control,
                         tuneLength = 10, 
                         data = train,
                         preProcess = "center"
)

rpeat.kfold.fit

plot(rpeat.kfold.fit$finalModel, line = TRUE)

extractPrediction(rpeat.kfold.fit, testX = test, testY = test)

points(test$logKd, predict(finalModel, newdata = test))

plot(test1$logKd,
     predict(rpeat.kfold.fit,newdata = test1)
)
abline(a = 0, b = 1)


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

####RUN MODEL on Full Spectrum - 1st Derivative####
set.seed(123)
deriv.fit <- plsr(logKd ~ MIR.1deriv,
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

summary(deriv.fit)

#RMSEP by Cross-Validation
plot(RMSEP(deriv.fit), legendpos = "topright")
minor.tick(nx = 2)

#Inspect Hotelling's T^2 for Outliers
plsr.scores <- data.frame(comp1 = plsr.fit$scores[,1], 
                          comp2 = plsr.fit$scores[,2],
                          comp3 = plsr.fit$scores[,3],
                          comp4 = plsr.fit$scores[,4],
                          comp5 = plsr.fit$scores[,5]
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

predplot(deriv.fit, 
         ncomp = 4, 
         line = TRUE,
         main = "Model Validation",
         xlim = c(1,3.5),
         ylim = c(1,3.5),
         pch = 21,
         bg = "gray"
)
####END####

spectra.df_LMW <- spectra.df_NoSed[grep("LMW", row.names(spectra.df_NoSed)),]
spectra.df_HMW <- spectra.df_NoSed[grep("HMW", row.names(spectra.df_NoSed)),]

set.seed(123)

split <- sample.split(spectra.df_HMW$Kd, SplitRatio = 20)

train_HMW <- subset(spectra.df_HMW, split == T)

test_HMW <- subset(spectra.df_HMW, split == F)

train_HMW <- cbind(train_HMW, n = 1:nrow(train_HMW))
test_HMW <- cbind(test_HMW, n = 1:nrow(test_HMW))

set.seed(123)
LMW.fit <- plsr(logKd ~ MIR,
                 ncomp = 10,
                 data = spectra.df_LMW,
                 validation = "CV",
                 segments = 5,
                 scale = FALSE,
                 center = TRUE,
                 method = "simpls" 
)
summary(LMW.fit)

plot(RMSEP(LMW.fit), legendpos = "topright")

predplot(LMW.fit, 
         ncomp = 3, 
         line = TRUE,
         main = "Model Validation",
         xlim = c(1,3.5),
         ylim = c(1,3.5),
         pch = 21,
         bg = "gray"
)

plot(LMW.fit, plottype = "scores", comps = 1:3, labels = "numbers")
loadingplot(LMW.fit, comps = 1:8, legend = "topleft")

####Repeated k-Fold####
set.seed(123)
train_control <- trainControl(method = "repeatedcv", 
                              number = 5, 
                              repeats = 5
)

set.seed(123)
rpeat.kfold.fit <- train(form = logKd ~
                           MIR.smooth,
                         method = "simpls",
                         trControl = train_control,
                         tuneLength = 10, 
                         data = train_HMW,
                         preProcess = "center"
)

rpeat.kfold.fit
plot(rpeat.kfold.fit$finalModel, line = TRUE)



set.seed(123)
HMW.fit <- plsr(logKd ~ MIR.smooth,
                ncomp = 10,
                data = train_HMW,
                validation = "CV",
                segments = 5,
                scale = FALSE,
                center = TRUE,
                method = "simpls" 
)
summary(HMW.fit)
plot(RMSEP(HMW.fit), legendpos = "topright")

plot(train_HMW$logKd,
     predict(HMW.fit, ncomp = 3),
     xlab = "measured",
     ylab = "predicted",
     main = "Model Validation",
     xlim = c(1, 3.5),
     ylim = c(1, 3.5),
     pch = 21,
     bg = "gray"
)
minor.tick(nx = 5)
abline(a = 0, b = 1)

points(test_HMW$logKd,
       predict(plsr.fit, ncomp = 3, newdata = test_HMW),
       pch = 21,
       bg = "blue"
)

plot(HMW.fit, plottype = "scores", comps = 1:2, labels = "numbers")
loadingplot(HMW.fit, comps = 1:2, legend = "topleft")
