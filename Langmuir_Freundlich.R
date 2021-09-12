
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