# ACST 3059 Assignment Session 2 2022
# Abdullah Shahid ID: 45761000

### Loading the required packages ###
set.seed(100)
library(demography)
library(forecast)
library(splines)

### Older function in order to pull from former HMD website
hmd.mx <- function(country, username, password, label=country)
{
  path <- paste("https://former.mortality.org/hmd/", country, "/STATS/", "Mx_1x1.txt", sep = "")
  userpwd <- paste(username, ":", password, sep = "")
  txt <- RCurl::getURL(path, userpwd = userpwd)
  con <- textConnection(txt)
  mx <- try(utils::read.table(con, skip = 2, header = TRUE, na.strings = "."),TRUE)
  close(con)
  if(class(mx)=="try-error")
    stop("Connection error at www.mortality.org. Please check username, password and country label.")
  
  path <- paste("https://former.mortality.org/hmd/", country, "/STATS/", "Exposures_1x1.txt", sep = "")
  userpwd <- paste(username, ":", password, sep = "")
  txt <- RCurl::getURL(path, userpwd = userpwd)
  con <- textConnection(txt)
  pop <- try(utils::read.table(con, skip = 2, header = TRUE, na.strings = "."),TRUE)
  close(con)
  if(class(pop)=="try-error")
    stop("Exposures file not found at www.mortality.org")
  
  obj <- list(type="mortality",label=label,lambda=0)
  
  obj$year <- sort(unique(mx[, 1]))
  #obj$year <- ts(obj$year, start=min(obj$year))
  n <- length(obj$year)
  m <- length(unique(mx[, 2]))
  obj$age <- mx[1:m, 2]
  mnames <- names(mx)[-c(1, 2)]
  n.mort <- length(mnames)
  obj$rate <- obj$pop <- list()
  for (i in 1:n.mort)
  {
    obj$rate[[i]] <- matrix(mx[, i + 2], nrow = m, ncol = n)
    obj$rate[[i]][obj$rate[[i]] < 0] <- NA
    obj$pop[[i]] <- matrix(pop[, i + 2], nrow = m, ncol = n)
    obj$pop[[i]][obj$pop[[i]] < 0] <- NA
    dimnames(obj$rate[[i]]) <- dimnames(obj$pop[[i]]) <- list(obj$age, obj$year)
  }
  names(obj$pop) = names(obj$rate) <- tolower(mnames)
  
  suppressWarnings(obj$age <- as.numeric(as.character(obj$age)))
  if (is.na(obj$age[m]))
    obj$age[m] <- 2 * obj$age[m - 1] - obj$age[m - 2]
  return(structure(obj, class = "demogdata"))
}
AusMort <- hmd.mx("AUS","saclryvftsaeuvsola@rffff.net","1629763659")
View(AusMort)

###................. INTRODUCTION .....................###
# Obtaining the ranges of Variables provided

# Variables Available for analysis
str(AusMort)

# The Variable Range present in the data set
summary(AusMort)
sapply(AusMort, range)
range(AusMort[["rate"]]["total"], na.rm = TRUE)
range(AusMort[["pop"]]["total"], na.rm = TRUE)

###................... Preliminary Data Analysis .......................###
# The latest observed year in the data set is 2019
# Obtaining the Mortality plots for Australia's population for the year 2019

par(mfrow=c(1,1))
plot(x = AusMort[["age"]], y = log(AusMort[["rate"]][["total"]][,"2019"]),
     type = "l", xlab = "Age", ylab = "Log - Mortality",
     main = "Total Australian Mortality rate for the year 2019") # Mortality rate for everyone in year 2019
# Adding Male mortality
lines(x = AusMort[["age"]], y = log(AusMort[["rate"]][["male"]][,"2019"]),type = "l", col = "blue",
     xlab = "Age", ylab = "Log - Mortality",
     main = "Total Mortality by gender for the year 2019") # Mortality rate for male in year 2019
# Adding Female mortality
lines(x = AusMort[["age"]], y = log(AusMort[["rate"]][["female"]][,"2019"]),type = "l", col = "red") # Mortality rate for female in year 2019

legend("bottomright" , c("Total","Male","Female"),
       cex=0.8,col=c("black","blue","red"),lty=1)


###........................ Parametric curve fitting – Spline models .........................###

#Defining the data sets used for the parametric curve fiitnig (Training, Validation and test)
age_range = 18:109

#2017 training
# Extracting training data set from the data provided
training = cbind(18:109, AusMort[["rate"]][["total"]][age_range, "2017"], 
                 AusMort[["pop"]][["total"]][age_range, "2017"])

#2018 validation
# Extracting validation data set from the data provided
validation = cbind(18:109, AusMort[["rate"]][["total"]][age_range, "2018"],
                   AusMort[["pop"]][["total"]][age_range, "2018"])
#2019 test
# Extracting test data set from the data provided
test = cbind(18:109, AusMort[["rate"]][["total"]][age_range, "2019"],
                AusMort[["pop"]][["total"]][age_range, "2019"])


#Assigning Column names to the new data sets for easier usability.
dataset.names = c("Age", "mu_x", "Ec_x")
# Ensuring No row names are attached to the data sets
rownames(training) = NULL
rownames(validation) = NULL
rownames(test) = NULL
# Assigning Column names to the training data set
colnames(training) = dataset.names
# Assigning Column names to the validation data set
colnames(validation) = dataset.names
# Assigning Column names to the test data set
colnames(test) = dataset.names

#Creating list for all possible knots combinations using combn function and then combining them into single list:

knots = c()
for(i in 0:8){
  #Combinations of possible knots
  # Only considering adult data so knots position 5 and 15 ignored
  knots_combination = combn(seq(25,95,by=10),i, simplify = FALSE)
  
  knots = c(knots,knots_combination)
}

# Obtaining the mean square error of all models
# The models are fitted using training data (Calculation of spline basis) set by (2017)
# The MSE is calculated using validation data set (2018)
# We require accurate prediction therefore we do not add penalty for higher number of knots to loss function

Mean_SQE = c()
for(i in 1:length(knots)){
  # Fitting the basis of Natural cubic spline using 2017 data set. 
   nc.basis <- ns(training[,"Age"],knots = knots[[i]]) 
  
  # Fitting the spline using liner model in order to estimate the coefficients of co-variate
   nc.spline <- lm(training[,"mu_x"] ~ nc.basis, weights = training[,"Ec_x"]) 
  
  # A list to store MSE values for all models (all possible knots)
   # The Ages in the training and validation data set are same 
   # We assume no mortality improvement
  # Therefore the fitted values of training data set are predictions for validation data set
   Mean_SQE[i] = mean ((validation[,"mu_x"] - fitted(nc.spline))^2)
  
}

# Finding the index for the lowest Mean Squared Error and obtaining the respective Knots.
Min_Error_index = which.min(Mean_SQE)
Final_Knots_Position = knots[Min_Error_index]
min(Mean_SQE)
Final_Knots_Position

# Confirming that the chosen set of knots provide lowest Validation MSE.
nc.basis_check <- ns(training[,"Age"],knots = Final_Knots_Position[[1]]) 
nc.spline_check <- lm(training[,"mu_x"] ~ nc.basis_check, weights = training[,"Ec_x"]) 
mean ((validation[,"mu_x"] - fitted(nc.spline_check))^2) == min(Mean_SQE)
nc.spline_check

####### Fitting Smoothing Spline #######

# We create a grid of different values for the smoothing parameter spar
hyper_parameter <- expand.grid(seq(0.001,1, by = 0.001), c(0))
colnames(hyper_parameter) <- c("Tuning_Parameter", "Error_Rate")

for(i in 1:dim(hyper_parameter)[1]){
  
# fitting a smooth spline model on the training data set using different values of hyper parameter   
sm.spline.train <- smooth.spline(training[,"Age"],training[,"mu_x"], w = training[,"Ec_x"], spar = hyper_parameter$Tuning_Parameter[i])

# predicting the values Mean squared error for validation data set using the model fitted on training data set
sm.spline.validation <- predict(sm.spline.train,validation[,"Age"])
Error = mean((validation[,"mu_x"]-sm.spline.validation$y)^2)

hyper_parameter$Error_Rate[i] = Error
}

# Determining the optimal level of tuning parameter
optimal_tuning = hyper_parameter[which.min(hyper_parameter$Error_Rate), 1]
min(hyper_parameter$Error_Rate)
optimal_tuning 

# In order to obtain a more precice value of lambda we can use Cross Validation on the validation data set

sm.spline.cv = smooth.spline(validation[,"Age"],validation[,"mu_x"],validation[,"Ec_x"] ,cv = T)
sm.spline.cv$spar

# We are required to use 2017 data as training dataset and 2018 data as validation data set. 
# Therefore we will continue with our spar value obtained in part 1
# Final Smoothing spline model based on training data-set
sm.spline.final <- smooth.spline(training[,"Age"],training[,"mu_x"], w = training[,"Ec_x"], spar = optimal_tuning )
sm.spline.final

# Checking the model implementation
# Implementing Natural Cubic spline on 2019 data set
nc.basis_2019 <- ns(test[,"Age"],knots = Final_Knots_Position[[1]]) 
nc.spline_2019 <- lm(test[,"mu_x"] ~ nc.basis_2019, weights = test[,"Ec_x"]) 
nc.spline_2019
# Implementing Smoothing Cubic spline on 2019 data set
sm.spline.2019 = smooth.spline(test[,"Age"],test[,"mu_x"],test[,"Ec_x"] ,spar = optimal_tuning)
sm.spline.2019

###....................Comparing BOTH MODELS ON 2019 DATA-SET....................###
#Plotting the crude and graduated mortality rate from Natural Cubic spline Model
# We compare the MSE for log(mortality) rate and use log scale for plots for better understanding of the mortality rates
par(mfrow=c(1,2))
plot(test[,"Age"],log(test[,"mu_x"]), main= "Natural Cubic Splines 2019 (With Choosen Knots)", xlab= "Age", ylab= "log mortality", col="black")
lines(test[,"Age"],log(fitted(nc.spline_check)), col="red", lwd= 2) 
legend("topleft",c("Log mu_x (2019)","Projected NCS"), col = c("black","red"), lwd = 1, cex = 0.5)

#Plotting the crude and graduated mortality rate from Smoothing Cubic spline Model

plot(test[,"Age"],log(test[,"mu_x"]), main= "Smoothing Splines 2019 (With Choosen lambda)", xlab= "Age", ylab= "log mortality", col="black")
lines(test[,"Age"],log(fitted(sm.spline.final)), col="red", lwd= 2) 
legend("topleft",c("Log mu_x (2019)","Projected Smooth.S"), col = c("black","red"), lwd = 1,cex = 0.5)

par(mfrow=c(1,1))
## Calculating Mean Squared errors of both models for model selection

# MSE between log(crude rates) and log(projected rates)
mean((log(test[,"mu_x"])-log(fitted(nc.spline_check)))^2)
mean((log(test[,"mu_x"])-log(fitted(sm.spline.final)))^2)

# Weighted MSE between crude rates and projected rates
mean(test[,"Ec_x"]*(test[,"mu_x"]-fitted(nc.spline_check))^2)/sum(test[,"Ec_x"])
mean(test[,"Ec_x"]*(test[,"mu_x"]-fitted(sm.spline.final))^2)/sum(test[,"Ec_x"])


# From both the graph and the mean square residual it is evident 
# that smooth spline model is a better fit for the data

###................Graduation Tests on Smooth Spline Model (2018 Data set) ................... ###

# Fitting the smooth spline model to 2018 data set.
sm.spline.2018 = smooth.spline(validation[,"Age"],validation[,"mu_x"],validation[,"Ec_x"] ,spar = optimal_tuning)
sm.spline.2018

# Combining the mortality data for 2018
AusMort_2018 =  cbind(validation[, "Age"], validation[, "mu_x"],
                      validation[, "Ec_x"], fitted(sm.spline.2018))
# Calculating the number of observations
No_observations = dim(AusMort_2018)[1]

# Assigning the column names to the 2018 data set
colnames(AusMort_2018) = c("Age", "Crude.Mu_x", "Ec_x", "Graduated.Mu_x")
AusMort_2018

# Chi-squared test of fit.............................
Actual_Deaths = AusMort_2018[,"Crude.Mu_x"]*AusMort_2018[,"Ec_x"]
Expected_Deaths = AusMort_2018[,"Graduated.Mu_x"]*AusMort_2018[,"Ec_x"]

Chisq_test_stat = sum(((Actual_Deaths - Expected_Deaths)^2)/Expected_Deaths)
Chisq_test_stat
# Smoothing spline is non parametric. One df lost because of estimation of hyper parameter
Chisq_test_df = No_observations - 1
Chisq_test_df 
Chisq_critical_value = qchisq(0.95, df = Chisq_test_df)
Chisq_critical_value

## Test statistic lower than critical value and high p value. Model fits the data well

# Standardized deviations test..........................

#Calculating Standardized deviations
z_x = (Actual_Deaths - Expected_Deaths)/sqrt(Expected_Deaths)

# We select 4 intervals (-inf,-1),(-1,0),(0,1) and (1 inf)
# Expected Number of observations in each interval
# It is calculated by multiplying standard normal probability with No of observations
EO = c()
# Expected observations less than -1
EO [ 1 ] = pnorm( -1 ) * No_observations
# Expected observations between 0 and -1
EO [ 2 ] = (pnorm(0)-pnorm(-1))* No_observations
# Expected observations between 0 and 1
EO [ 3 ] = (pnorm(1)-pnorm(0))* No_observations
# Expected observations more than 1
EO [ 4 ] = (1-pnorm(1))* No_observations
EO 

# Actual Number of observations in each interval
AO = c()
# Actual observations less than -1
AO [ 1 ] = sum ( z_x <= (-1))
# Actual observations between 0 and -1
AO [ 2 ] = sum ( (z_x <= (0) )) - sum(( z_x <= (-1)))
# Actual observations between 0 and 1
AO [ 3 ] = sum ( (z_x > (0) )) - sum((z_x >= (1) ))
# Actual observations more than 1
AO [ 4 ] = sum( z_x >= (1))
AO

# The CHI statistics for each group
SD_cont = c()
# The individual contribution to test statistics
# Contribution by section 1
SD_cont [ 1 ] = (( AO [ 1 ] - EO [ 1 ]) ^2)/ EO [ 1 ]
# Contribution by section 2
SD_cont [ 2 ] = (( AO [ 2 ] - EO [ 2 ]) ^2)/ EO [ 2 ]
# Contribution by section 3
SD_cont [ 3 ] = (( AO [ 3 ] - EO [ 3 ]) ^2)/ EO [ 3 ]
# Contribution by section 4
SD_cont [ 4 ] = (( AO [ 4 ] - EO [ 4 ]) ^2)/ EO [ 4 ]
SD_cont

# Chi squared test statistic
z_x_test_stat = sum ( SD_cont )
z_x_test_stat
#Chis squared degree fo freedom n-1
z_x_df =length(SD_cont)-1

# Chi squared critical value
qchisq(0.95,z_x_df)

# The test statistic is below the critical value insufficient evidence to reject H0
#Signs test...............................
# Calculating the number of positive deviations
deviations_death = c(Actual_Deaths - Expected_Deaths) 
pos_deviations = length(which(deviations_death>0))
pos_deviations

# Calculating the Acceptance region
lower.criteria = qbinom(0.025,No_observations,0.5)
lower.criteria
upper.criteria = qbinom(0.975,No_observations,0.5)
upper.criteria

# The positive deviations fall in the Acceptance region, Accept Null hypothesis. Model fits the data well

#Cumulative deviations test.......................

# Calculating the test statistic
CD_test_stat = (sum(Actual_Deaths - Expected_Deaths))/sqrt(sum(Expected_Deaths))
CD_test_stat
CD_pvalue = pnorm(CD_test_stat)
CD_pvalue
CD_crit = qnorm(c(0.025,0.975))
CD_crit
deviations_death
# The test statistic falls within our Acceptance region, therefore we accept the null hypothesis
# It can be observed that very large negative deviation are being cancelled by large positive deviation.

# Grouping of signs test................................

# Getting total deviations, Positive deviations and negative deviations
Total_deviations = No_observations
pos_d = pos_deviations
pos_d
neg_d = Total_deviations - pos_d
neg_d

# Calculating the number of positive groups G
# Checking sign of the first deviation
if( z_x[1]>0){
  G_start = 1 }else{
  G_start = 0
  }
# Calculating total number of groups of positive deviations
for(i in 2:Total_deviations){
  if( (z_x[i-1]>0 && z_x[i] <0 ) ){
    G_start = G_start+1
  }
}
G = G_start
G

# total deviations = 92 therefore we can use normal approximation
## Mean and variance
mean_GS_test = ( pos_d * (neg_d + 1))/ (pos_d+neg_d)
mean_GS_test 
var_Gs_test = ((pos_d * neg_d) ^2)/( (pos_d + neg_d)^3)
var_Gs_test 

# Critical value k such that P(G <= k) >= 0.05
critical_val_GS = qnorm(0.05,mean_GS_test,sqrt(var_Gs_test))
critical_val_GS

# Finding discrete critical value
pnorm(19,mean_GS_test,sqrt(var_Gs_test))
pnorm(20,mean_GS_test,sqrt(var_Gs_test))

# Therefore critical value is 20. The number of G 26 is more than critical value. Insufficient evidence to reject H0

# Serial Correlation test ................................
# Two sequences of standardized deviations of length 92-1
z1 = z_x[1:(No_observations-1)]
z2 = z_x[2:(No_observations)]
mu_z1 = mean(z1)
mu_z1
mu_z2 =  mean(z2)
mu_z2  
# Calculation of test statistic of serial correlation test  
r_1 = sum((z1 - mu_z1 )*(z2 - mu_z2))/sqrt(sum((z1-mu_z1)^2) * sum((z2-mu_z2)^2))
r_1
SCT_test_stat = r_1*sqrt(No_observations)
SCT_test_stat

# Calculation of the critical value for test statistic
qnorm(0.95)

# Test statistic is lower than the critical vale. We only use SCT to check over graduation. Insufficient evidence to reject H0

###............................ Fitting Lee Carter Model to the data up to 2018 ......................###

# From the Ausmort data set it can be observed that the Data is missing for the ages 104+
# Obtaining data to fit the lee carter model
Maximum_age = 101
lc_data = extract.ages(extract.years(data = AusMort, years = 1921:2018),ages = 18:Maximum_age)
lc_data_2019 = extract.ages(extract.years(data = AusMort, years = 2019),ages = 18:Maximum_age)

# Fitting Lee Carter Model
LEE_C_total = demography::lca(lc_data,series = "total",max.age = Maximum_age)
LEE_C_male = demography::lca(lc_data,series = "male",max.age = Maximum_age)
LEE_C_female = demography::lca(lc_data,series = "female",max.age = Maximum_age)
#Graphical Representation of LC model residual
par(mfrow=c(1,1))
heatmap(LEE_C_total$residuals$y, Rowv=NA, Colv=NA, col = terrain.colors(300), main = "Heat map")
#Graphical Representation of LC model
plot(LEE_C_total$fitted, main ="LC model for adult mortality rates")

# MSE of lee carter projection for the year 2019...............................................
lc_2019 = forecast(LEE_C_total,h=1)
lc_MSE = mean ( ( lc_2019$rate$total - lc_data_2019$rate$total) ^2)
lc_MSE

# Graph to get understanding of MSE
lc_2019_resid = lc_2019$rate$total- lc_data_2019$rate$total
plot(18:101,lc_2019_resid, xlab = "Age", ylab = "Mortality Residuals", main = "Lee Carter Model residual evaluation")
abline(0,0, col="red")

# Plots for Lee Carter Model Parameters..............................................
par(mfrow = c(1, 3))
#Plotting Alpha
# Creating the overall mortality experience from 1921 to 2018
plot(LEE_C_total$fitted, main= "Alpha (x)", xlab= "Age", ylab= "A(x)", type= "l", col = "azure3")

# Plotting the alpha values for the three models segregated by gender
lines(x= LEE_C_total$age, y=LEE_C_total$ax, col= "black")
# Plotting the alpha values for female
lines(x= LEE_C_female$age, y= LEE_C_female$ax, col= "red")
# Plotting the alpha values for male
lines(x= LEE_C_male$age, y= LEE_C_male$ax, col= "blue") 
# legend for the plot
legend("bottomright", c("Male", "Female", "Total"), cex = 1, col =
         c("blue", "red", "black"), lty = 1)

#Plotting Beta for female
plot(x= LEE_C_female$age, y= LEE_C_female$bx, main= "Beta (x)", xlab= "Age", col= "red", ylab= "B(x)", type= "l")
#Plotting Beta for total
lines(x= LEE_C_total$age , y= LEE_C_total$bx, col= "black")
#Plotting Beta for male
lines(x= LEE_C_male$age, y= LEE_C_male$bx, col= "blue")
legend("bottomleft", c("Male", "Female", "Total"), cex = 1, col =
         c("blue", "red", "black"), lty = 1)

#Plotting kappa for total 
plot(x = LEE_C_total$year, y = LEE_C_total$kt, main= "Kappa (t)", xlab= "Year", ylab= "K(t)", type= "l")
#Plotting kappa for female 
lines(x= LEE_C_female$year, y= LEE_C_female$kt, col= "red")
#Plotting kappa for male
lines(x= LEE_C_male$year, y= LEE_C_male$kt, col= "blue")
legend("bottomleft", c("Male", "Female", "Total"), cex = 1, col = c("blue", "red", "black"), lty = 1)

###...................... Model comparison ....................###
# We assume there is no mortality improvement in the spline model
# Therefore we use the graduated mortality mortality rates for 2018 for all the years 2030,2040 and 2050

# obtaining the lee carter model projection upto year 2050
lc2050 = forecast(LEE_C_total,h=32)

par(mfrow=c(1,4))
# Plotting Fall in mortality rate and variation
plot(lc2050$kt.f, main="LC Total Mortality Rates", xlab="Years", ylab="Mortality Rates")

# Plotting graphs for comparison

#2030 Model Plot
plot(lc2050$age, log(lc2050[["rate"]][["total"]][,12]), xlab = "Age", ylab = "log mu_x", type = "l",
     main = "Projected mortality 2030")
lines(validation[,"Age"],log(fitted(sm.spline.2018)), col = "red")
legend("topleft", c("Lee Carter", "Smoothing spline"), cex = 1, col = c("black", "red"),lty = 1)

#2040 Model Plot
plot(lc2050$age, log(lc2050[["rate"]][["total"]][,22]), xlab = "Age", ylab = "log mu_x", type = "l",
     main = "Projected mortality 2040")
lines(validation[,"Age"],log(fitted(sm.spline.2018)), col = "red")
legend("topleft", c("Lee Carter", "Smoothing spline"), cex = 1, col = c("black", "red"), lty =1)

#2050 Model Plot
plot(lc2050$age, log(lc2050[["rate"]][["total"]][,32]), xlab = "Age", ylab = "log mu_x", type = "l",
     main = "Projected mortality 2050")
lines(validation[,"Age"],log(fitted(sm.spline.2018)), col = "red")
legend("topleft", c("Lee Carter", "Smoothing spline"), cex = 1, col = c("black", "red"), lty=1)
     
