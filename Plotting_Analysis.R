####Calculate Second Derivative for Peak Assignment####
#Second derivative
MIR.2deriv <- apply(spectra.df$MIR.smooth, 
                   MARGIN = 1, 
                   FUN = sgolayfilt, 
                   p = 2, n = 15, m = 2, ts = 1
)

spectra.df$MIR.2deriv <- t(MIR.2deriv)

#Create a List of All Peaks Calculated for 2nd Derivative Spectra
peaks.list <- apply(
  spectra.df$MIR.2deriv, 
  MARGIN = 1, 
  function(i)wavenumbers[find_peaks(-i, m = 3),]
) 

#Retain Only Region Below 1800 cm-1
peaks.list <- lapply(peaks.list, function(i)i[which(i < 1800)])
####END####

####Extract Peaks for HMW End-Members####
#TNOM HMW
TNOM_HMW_peaks <- peaks.list$`Apr06_TPC_1-0-0_HMW`

#Plot Second Derivative 
plot(wavenumbers[,1],
     spectra.df["Apr06_TPC_1-0-0_HMW","MIR.2deriv"],
     type = "l",
     xlim = c(1800, 900),
     main = "Apr6_TPC_1-0-0_HMW"
)

#Plot 
plot(wavenumbers[,1],
     spectra.df["Apr06_TPC_1-0-0_HMW","MIR.smooth"],
     type = "l",
     xlim = c(1800, 900),
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
####END####

####Extract Peaks for LMW End-Members####
TNOM_LMW_peaks <- peaks.list$`Apr08_TPC_1-0-0_LMW`
#Plot 
plot(wavenumbers[,1],
     smoothed_spectra$`Apr08_TPC_1-0-0_LMW`,
     type = "l",
     xlim = c(3600, 900),
     main = "Apr08_TPC_1-0-0_LMW"
)

minor.tick(nx = 5)
abline(v = TNOM_LMW_peaks, col = "red")
#Triage Peaks by Eye
TNOM_LMW_peaks <- c(TNOM_LMW_peaks[c(3,9,23,26,28:31,33,38,40)],1433,1370)
TNOM_LMW_peaks[16]

#CNOM LMW
CNOM_LMW_peaks <- peaks.list$`Apr12_TPC_0-0-1_LMW`
#Plot 
plot(wavenumbers[,1],
     #spectra_avg$`Apr12_TPC_0-0-1_LMW`,
     smoothed_spectra$`Apr12_TPC_0-0-1_LMW`,
     type = "l",
     xlim = c(3600, 900),
     main = "Apr12_TPC_0-0-1_LMW"
)

minor.tick(nx = 5)
abline(v = CNOM_LMW_peaks, col = "coral", lty = 2)
#Triage Peaks by Eye
CNOM_LMW_peaks <- c(CNOM_LMW_peaks[c(6,9,11,14,17,20,23,27,31,36,37,39,44,45)],1108)

#PNOM LMW
PNOM_LMW_peaks <- peaks.list$`Apr09_TPC_0-1-0_LMW`
#Plot 
plot(wavenumbers[,1],
     smoothed_spectra$`Apr09_TPC_0-1-0_LMW`,
     #spectra_avg$`Apr09_TPC_0-1-0_LMW`,
     type = "l",
     xlim = c(3600, 900),
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

####Plot Regression Coefficients####
#Choose Number of Components
Num.components <- 9

#Graphical parameter reset
dev.off()
#or
par(mfrow=c(1,1))

par(mar=c(0.5, 0.5, 0.2, 0.2), 
    mfrow=c(1,2),
    oma = c(5, 5, 3, 0.5)
)
layout(matrix(c(1,2,2), nrow = 1, ncol = 3))


plot(wavenumbers_NoCO2,
     plsr.fit$coefficients[,1,Num.components],
     type = "l",
     xlim = c(3600,2300),
     xlab = NA,
     #xlab = expression("Wavenumber (cm"^-1*")"),
     ylab = "Coefficient"#,
     #main = paste(c("Regression Coefficients by Wavenumber", Num.components, "Components"))
)
abline(h = 0, col = "blue")
#minor.tick(nx = 5)

abline(v = c(3386,3146,2922), 
       col = "coral", lty = 2)

plot(wavenumbers_NoCO2,
     plsr.fit$coefficients[,1,Num.components],
     type = "l",
     xlim = c(1900,900),
     xlab = NA,
     #xlab = expression("Wavenumber (cm"^-1*")"),
     #ylab = "Coefficient"#,
     yaxt = "n",
     ylab = "",
     #main = paste(c("Regression Coefficients by Wavenumber", Num.components, "Components"))
)
abline(h = 0, col = "blue")

abline(v = c(1720,1657,1577,1435,1388,1288,1232,1159,1109,1023,984), 
       col = "coral", lty = 2)
####END####

####Plot HMW Fingerprint Region####
smoothed_spectra <- data.frame(MIR.smooth, check.names = FALSE)

#Graphical parameter reset
dev.off()
#or
par(mfrow=c(1,1))

par(mar=c(0.5, 0.5, 0.2, 0.2), 
    mfrow=c(1,2),
    oma = c(5, 5, 3, 0.2)
)
layout(matrix(c(1,2,2), nrow = 1, ncol = 3))

plot(wavenumbers[,1],
     smoothed_spectra$`Apr06_TPC_1-0-0_HMW`,
     #spectra_avg_corr$`Apr06_TPC_1-0-0_HMW`,
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
      smoothed_spectra$`Mar22_TPC_0-1-0_HMW`+0.5,
      col = "dark blue"
)

lines(wavenumbers[,1],
      smoothed_spectra$`Mar30_TPC_0-0-1_HMW`+1,
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
     smoothed_spectra$`Apr06_TPC_1-0-0_HMW`,
     type = "l",
     #xlim = c(1800, 800),
     xlim = c(1900, 900),
     ylim = c(0,2.3),
     xlab = expression("Wavenumber (cm"^-1*")"),
     yaxt = "n",
     ylab = "",
     col = "dark green"
)

lines(wavenumbers[,1],
      smoothed_spectra$`Mar22_TPC_0-1-0_HMW`+0.5,
      col = "dark blue"
)

lines(wavenumbers[,1],
      smoothed_spectra$`Mar30_TPC_0-0-1_HMW`+1,
      col = "red"
)

minor.tick(nx = 4)

abline(v = c(1729,1630,1576,1560,1511,1448,1397,1314,1239,1124,1038,975), lty = 2, col = "darkgrey")

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

#Graphical parameter reset
dev.off()
#or
par(mfrow=c(1,1))

par(mar=c(0.5, 0.5, 0.2, 0.2), 
    mfrow=c(1,2),
    oma = c(5, 5, 3, 0.2)
)
layout(matrix(c(1,2,2), nrow = 1, ncol = 3))

plot(wavenumbers[,1],
     smoothed_spectra$`Apr08_TPC_1-0-0_LMW`+1.6,
     type = "l",
     xlim = c(3600, 2300),
     ylim = c(0,3),
     xlab = expression("Wavenumber (cm"^-1*")"),
     ylab = "Absorbance",
     col = "dark green"
)

lines(wavenumbers[,1],
      smoothed_spectra$`Apr09_TPC_0-1-0_LMW`,#+0.5,
      col = "dark blue"
)

lines(wavenumbers[,1],
      smoothed_spectra$`Apr12_TPC_0-0-1_LMW`+0.8,#1.1,
      col = "red"
)

minor.tick(nx = 2)

abline(v = c(3386),lty =2, col = "gray")
abline(v = c(1762,1654,1636,1570,1508,1408,1370,1120,1108,1049), lty =2, col = "gray")
####END###

####Labels LMW####
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

#####1800-800 Range LMW####
plot(wavenumbers[,1],
     smoothed_spectra$`Apr08_TPC_1-0-0_LMW`+1.6,
     type = "l",
     #xlim = c(1800, 800),
     xlim = c(1900, 900),
     ylim = c(0,3),
     xlab = expression("Wavenumber (cm"^-1*")"),
     yaxt = "n",
     ylab = "",
     col = "dark green"
)

lines(wavenumbers[,1],
      smoothed_spectra$`Apr09_TPC_0-1-0_LMW`,#+0.5,
      col = "dark blue"
)

lines(wavenumbers[,1],
      smoothed_spectra$`Apr12_TPC_0-0-1_LMW`+0.8,#+1.1,
      col = "red"
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