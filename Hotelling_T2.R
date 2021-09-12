#Hotelling T^2 = t^2_1/s^2_1 + t^2_2/s^2_2
#for two components (Dunn, K.)

#Test this function vs. HotellingEllipse package

#Extract scores
scores_1 <- plsr.fit$scores[,1]
scores_2 <- plsr.fit$scores[,2]


#Compute variances
var_1 <- var(scores_1)
var_2 <- var(scores_2)

#Check var() results with good ole' arithmetic
var_1_check <- (scores_1 - mean(scores_1))^2
var_1_check <- sum(var_1_check)/(length(scores_1)-1)
var_1 == var_1_check

#Create dataframe of T^2 values
Tsquared <- data.frame((scores_1^2)/var_1^2 + (scores_2^2)/var_2^2)
colnames(Tsquared) <- "T2"

#T^2 threshold defined as: 
# (A * (n - 1)/(n - A)) * stats::qf(p = 0.99, df1 = A, df2 = (n - A))
# where A is no. of components i.e. parameters

#p-value for individual observations is:
# ((n - A)/(A * (n - 1))) * MDsq)
#Where MDsq is the squared Mahalanobis distance for each value of x