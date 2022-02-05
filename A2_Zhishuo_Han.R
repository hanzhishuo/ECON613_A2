library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(plm)
options(scipen=999)
indpath = "D:/duke/ECON613/A2/ind"

assign(paste0("datind2009"), fread(paste0(indpath,"/","datind2009.csv")))

#exercise 1 OLS estimate

df <- datind2009 %>%
  filter(!is.na(wage) & !is.na(age))

#Calculate the correlation between Y and X.
correlation = cor(df[, "wage"], df[, "age"])

#Calculate the coefficients on this regression.

X = as.matrix(cbind(1, df$age))
Y = as.matrix(df$wage)

beta_hat = solve(t(X)%*%X)%*%t(X)%*%Y

#Calculate the standard errors of beta
#Using the standard formulas of the OLS

residual = as.matrix(df$wage - beta_hat[1] - beta_hat[2]*df$age)
n = nrow(df)
k = ncol(X)

Var = 1/(n-k) * as.numeric(t(residual)%*%residual) * solve(t(X)%*%X)
SE = sqrt(diag(Var))

#Using bootstrap with 49 and 499 replications respectively. Comment on the difference between
#the two strategies.


R = 49
nind = nrow(df)
set.seed(123)
outs = mat.or.vec(R, 2)
for (i in 1:R)
{
  samp = sample(1:nind, nind, rep=TRUE)
  dat_samp = df[samp,]
  X = as.matrix(cbind(1, dat_samp$age))
  Y = as.matrix(dat_samp$wage)
  beta_hat = solve(t(X)%*%X)%*%t(X)%*%Y

  residual = as.matrix(dat_samp$wage - beta_hat[1] - beta_hat[2]*dat_samp$age)
  n = nrow(dat_samp)
  k = ncol(X)
  Var = 1/(n-k) * as.numeric(t(residual)%*%residual) * solve(t(X)%*%X)
  outs[i,] = sqrt(diag(Var))
}

mean_est_49 = apply(outs, 2, mean)
sd_est_49 = apply(outs, 2, sd)

R = 499
nind = nrow(df)
outs = mat.or.vec(R, 2)
set.seed(123)

for (i in 1:R)
{
  samp = sample(1:nind, nind, rep=TRUE)
  dat_samp = df[samp,]
  X = as.matrix(cbind(1, dat_samp$age))
  Y = as.matrix(dat_samp$wage)
  beta_hat = solve(t(X)%*%X)%*%t(X)%*%Y
  
  residual = as.matrix(dat_samp$wage - beta_hat[1] - beta_hat[2]*dat_samp$age)
  n = nrow(dat_samp)
  k = ncol(X)
  Var = 1/(n-k) * as.numeric(t(residual)%*%residual) * solve(t(X)%*%X)
  outs[i,] = sqrt(diag(Var))
}

mean_est_499 = apply(outs, 2, mean)
sd_est_499 = apply(outs, 2, sd)

# Exercise 2 Detrend Data
#Read all individual datasets from 2005 to 2018. Append all these datasets.
setwd(indpath)
ind_file = list.files(path=indpath)
ind2 <- lapply(ind_file, function(x){
  df <- fread(file = x, stringsAsFactors =FALSE)
  df$idind <- as.character(df$idind)
  df$idmen <- as.character(df$idmen)
  df
  
})

ind <- do.call(rbind, ind2)
ind <- ind %>%
  filter(year >= 2005 & year <= 2018) %>%
  filter(!is.na(wage) & !is.na(age))

ind$ag <- as.factor(ifelse(ind$age < 18, 0, 
                    ifelse(ind$age >= 18 & ind$age <= 25, 1, 
                    ifelse(ind$age >= 26 & ind$age <= 30, 2,
                    ifelse(ind$age >= 31 & ind$age <= 35, 3,
                    ifelse(ind$age >= 36 & ind$age <= 40, 4,
                    ifelse(ind$age >= 41 & ind$age <= 45, 5,
                    ifelse(ind$age >= 46 & ind$age <= 50, 6,
                    ifelse(ind$age >= 51 & ind$age <= 55, 7,
                    ifelse(ind$age >= 56 & ind$age <= 60, 8,
                    ifelse(ind$age > 60, 9, 0)))))))))))

#Plot the wage of each age group across years. Is there a trend?
wage_age <- ind %>%
  filter(ag != 0) %>%
  group_by(ag, year)%>%
  mutate(mean_wage = mean(wage))%>%
  distinct(ag, year, .keep_all = TRUE)%>%
  ungroup()%>%
  select(ag, year, mean_wage)


p1 <- ggplot(wage_age,aes(x=year, y=mean_wage, color = ag, fill = ag)) + 
  geom_line() + facet_wrap(~ag)


#After including a time fixed effect, how do the estimated coefficients change?
  
reg_OLS = lm(wage ~ age, data=ind)
reg_OLS$coefficients

reg_fixed = plm(wage ~ age, data=ind, index = "year", model = "within")
reg_fixed$coefficients



#Exercise 3 Numerical Optimization

assign(paste0("datind2007"), fread(paste0(indpath,"/","datind2007.csv")))

#Exclude all individuals who are inactive.
ind_2007 <- datind2007%>%
  filter(empstat != "Inactive") %>%
  mutate(emp_dummy = ifelse(empstat == "Employed", 1, 0))

#Write a function that returns the likelihood of the probit of being employed.

flike = function(par, age, emp_dummy)
{
  xbeta           = par[1] + par[2]*age
  pr              = pnorm(xbeta)
  pr[pr>0.999999] = 0.999999
  pr[pr<0.000001] = 0.000001
  like           = emp_dummy*log(pr) + (1-emp_dummy)*log(1-pr)
  return(-sum(like))
}

#Optimize the model and interpret the coefficients.
reg1 <- glm(emp_dummy ~ age, data = ind_2007, family = binomial(link = "probit"))
start <- reg1$coefficients
flike(start, ind_2007$age, ind_2007$emp_dummy)
logLik(reg1)
res = optim(start,fn=flike,method="BFGS", age = ind_2007$age, emp_dummy = ind_2007$emp_dummy)

#Can you estimate the same model including wages as a determinant of labor market participation? Explain.
reg2 = glm(emp_dummy ~ age + wage, data = ind_2007, family = binomial(link = "probit"))

#Exercise 4 Discrete choice

setwd(indpath)
ind_file = list.files(path=indpath)
ind2 <- lapply(ind_file, function(x){
  df <- fread(file = x, stringsAsFactors =FALSE)
  df$idind <- as.character(df$idind)
  df$idmen <- as.character(df$idmen)
  df
  
})

ind <- do.call(rbind, ind2)

#Exclude all individuals who are inactive.
ind_0515 <- ind %>%
  filter(year >= 2005 & year <= 2015) %>%
  filter(empstat != "Inactive")

#Write and optimize the probit, logit, and the linear probability models.
# probit

ind_0515 <- ind_0515%>%
  mutate(emp_dummy = ifelse(empstat == "Employed", 1, 0)) %>%
  mutate(year_2006 = ifelse(year == 2006, 1, 0),
         year_2007 = ifelse(year == 2007, 1, 0),
         year_2008 = ifelse(year == 2008, 1, 0),
         year_2009 = ifelse(year == 2009, 1, 0),
         year_2010 = ifelse(year == 2010, 1, 0),
         year_2011 = ifelse(year == 2011, 1, 0),
         year_2012 = ifelse(year == 2012, 1, 0),
         year_2013 = ifelse(year == 2013, 1, 0),
         year_2014 = ifelse(year == 2014, 1, 0),
         year_2015 = ifelse(year == 2015, 1, 0)) 

probit_flike = function(par, age, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, emp_dummy)
{
  xbeta              = par[1] + par[2]*age + par[3]*x1 + par[4]*x2 + par[5]*x3 + par[6]*x4 + 
    par[7]*x5 + par[8]*x6 + par[9]*x7 + par[10]*x8 + par[11]*x9 + par[12]*x10
  pr                 = pnorm(xbeta)
  pr[pr>0.999999] = 0.999999
  pr[pr<0.000001] = 0.000001
  like               = emp_dummy*log(pr) + (1-emp_dummy)*log(1-pr)
  return(-sum(like))
}

#Optimize the model
reg3 <- glm(emp_dummy ~ age+year_2006+year_2007+year_2008+year_2009+year_2010+year_2011+year_2012+year_2013+year_2014+year_2015, 
            data = ind_0515, family = binomial(link = "probit"))
start_probit <- reg3$coefficients

res_probit = optim(start_probit,fn=probit_flike, method="BFGS", age = ind_0515$age, emp_dummy = ind_0515$emp_dummy,
            x1 = ind_0515$year_2006, x2 = ind_0515$year_2007, x3 = ind_0515$year_2008, x4 = ind_0515$year_2009,
            x5 = ind_0515$year_2010, x6 = ind_0515$year_2011, x7 = ind_0515$year_2012, x8 = ind_0515$year_2013,
            x9 = ind_0515$year_2014, x10 = ind_0515$year_2015, hessian = TRUE)


#logit

logit_flike = function(par, age, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, emp_dummy)
{
  xbeta           = par[1] + par[2]*age + par[3]*x1 + par[4]*x2 + par[5]*x3 + par[6]*x4 + 
    par[7]*x5 + par[8]*x6 + par[9]*x7 + par[10]*x8 + par[11]*x9 + par[12]*x10
  pr              = exp(xbeta)/(1+exp(xbeta))
  pr[pr>0.999999] = 0.999999
  pr[pr<0.000001] = 0.000001
  like           = emp_dummy*log(pr) + (1-emp_dummy)*log(1-pr)
  return(-sum(like))
}

#Optimize the model
reg4 <- glm(emp_dummy ~ age+year_2006+year_2007+year_2008+year_2009+year_2010+year_2011+year_2012+year_2013+year_2014+year_2015, 
            data = ind_0515, family = binomial(link = "logit"))
start_logit <- reg4$coefficients

res_logit = optim(start_logit,fn=logit_flike,method="BFGS", age = ind_0515$age, emp_dummy = ind_0515$emp_dummy,
                   x1 = ind_0515$year_2006, x2 = ind_0515$year_2007, x3 = ind_0515$year_2008, x4 = ind_0515$year_2009,
                   x5 = ind_0515$year_2010, x6 = ind_0515$year_2011, x7 = ind_0515$year_2012, x8 = ind_0515$year_2013,
                   x9 = ind_0515$year_2014, x10 = ind_0515$year_2015, hessian = TRUE)

#Linear Probability Model

lpm_flike = function(par, age, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, emp_dummy)
{
  pr          = par[1] + par[2]*age + par[3]*x1 + par[4]*x2 + par[5]*x3 + par[6]*x4 + 
    par[7]*x5 + par[8]*x6 + par[9]*x7 + par[10]*x8 + par[11]*x9 + par[12]*x10
  pr[pr>0.999999] = 0.999999
  pr[pr<0.000001] = 0.000001
  like           = emp_dummy*log(pr) + (1-emp_dummy)*log(1-pr)
  return(-sum(like))
}

#Optimize the model
reg5 <- lm(emp_dummy ~ age+year_2006+year_2007+year_2008+year_2009+year_2010+year_2011+year_2012+year_2013+year_2014+year_2015, 
            data = ind_0515)
start_lpm <- reg5$coefficients
res_lpm = optim(start_lpm,fn=lpm_flike,method="BFGS", age = ind_0515$age, emp_dummy = ind_0515$emp_dummy,
                  x1 = ind_0515$year_2006, x2 = ind_0515$year_2007, x3 = ind_0515$year_2008, x4 = ind_0515$year_2009,
                  x5 = ind_0515$year_2010, x6 = ind_0515$year_2011, x7 = ind_0515$year_2012, x8 = ind_0515$year_2013,
                  x9 = ind_0515$year_2014, x10 = ind_0515$year_2015)



#Interpret and compare the estimated coefficients. How significant are they
#probit standard error & t-statistic
probit_se = sqrt(diag(solve(res_probit$hessian)))
probit_t = res_probit$par/probit_se
#logit standard error & t-statistic
logit_se = sqrt(diag(solve(res_logit$hessian)))
logit_t = res_logit$par/logit_se
#lpm standard error & t-statistic
X = as.matrix(cbind(1, ind_0515$age, ind_0515$year_2006, ind_0515$year_2007, ind_0515$year_2008, 
                    ind_0515$year_2009, ind_0515$year_2010, ind_0515$year_2011, ind_0515$year_2012,
                    ind_0515$year_2013, ind_0515$year_2014, ind_0515$year_2015))
Y = as.matrix(ind_0515$emp_dummy)

pr = res_lpm$par[1] + res_lpm$par[2]*ind_0515$age + res_lpm$par[3]*ind_0515$year_2006 + 
  res_lpm$par[4]*ind_0515$year_2007 + res_lpm$par[5]*ind_0515$year_2008 + res_lpm$par[6]*ind_0515$year_2009 +
  res_lpm$par[7]*ind_0515$year_2010 + res_lpm$par[8]*ind_0515$year_2011 + res_lpm$par[9]*ind_0515$year_2012 + 
  res_lpm$par[10]*ind_0515$year_2013 + res_lpm$par[11]*ind_0515$year_2014 + res_lpm$par[12]*ind_0515$year_2015
pr[pr >0.999999] = 0.999999
pr[pr <0.000001] = 0.000001

Y_pr  = as.numeric(pr>=0.5)

residual = as.matrix(ind_0515$emp_dummy - Y_pr)
n = nrow(ind_0515)
k = ncol(X)

Var = 1/(n-k) * as.numeric(t(residual)%*%residual) * solve(t(X)%*%X)
lpm_se = sqrt(diag(Var))
lpm_t = res_lpm$par/lpm_se


#Exercise 5 Marginal Effects

#marginal effect of probit

predict_y = res_probit$par[1] + res_probit$par[2]*ind_0515$age + res_probit$par[3]*ind_0515$year_2006 +
  res_probit$par[4]*ind_0515$year_2007 + res_probit$par[5]*ind_0515$year_2008 + res_probit$par[6]*ind_0515$year_2009 + 
  res_probit$par[7]*ind_0515$year_2010 + res_probit$par[8]*ind_0515$year_2011 + res_probit$par[9]*ind_0515$year_2012 +
  res_probit$par[10]*ind_0515$year_2013 + res_probit$par[11]*ind_0515$year_2014 + res_probit$par[12]*ind_0515$year_2015
pdf_probit = mean(dnorm(predict_y))  
marginal_effect_probit = pdf_probit * res_probit$par

#marginal effect of logit

pdf_logit = mean(dlogis(predict_y))  
marginal_effect_logit = pdf_logit * res_probit$par


#Standard errors
#probit
boot <- 20
bootvals <- matrix(rep(NA, boot*length(res_probit$par)), nrow=boot)
set.seed(1111)

for(i in 1:boot){
  samp1 <- ind_0515[sample(1:dim(ind_0515)[1], replace=T, dim(ind_0515)[1]),]
  x1 <- optim(start_probit,fn=probit_flike, method="BFGS", age = samp1$age, emp_dummy = samp1$emp_dummy,
                     x1 = samp1$year_2006, x2 = samp1$year_2007, x3 = samp1$year_2008, x4 = samp1$year_2009,
                     x5 = samp1$year_2010, x6 = samp1$year_2011, x7 = samp1$year_2012, x8 = samp1$year_2013,
                     x9 = samp1$year_2014, x10 = samp1$year_2015, hessian = TRUE)
  predict1 <- x1$par[1] + x1$par[2]*samp1$age + x1$par[3]*samp1$year_2006 +
    x1$par[4]*samp1$year_2007 + x1$par[5]*samp1$year_2008 + x1$par[6]*samp1$year_2009 + 
    x1$par[7]*samp1$year_2010 + x1$par[8]*samp1$year_2011 + x1$par[9]*samp1$year_2012 +
    x1$par[10]*samp1$year_2013 + x1$par[11]*samp1$year_2014 + x1$par[12]*samp1$year_2015
  pdf1 <- mean(dnorm(predict1))
  bootvals[i,] <- pdf1 * x1$par
}

res <- cbind(marginal_effect_probit,apply(bootvals,2,sd),marginal_effect_probit/apply(bootvals,2,sd))

colnames(res) <- c("marginal.effect","standard.error","z.ratio")  

#logit
boot <- 20
bootvals2 <- matrix(rep(NA, boot*length(res_probit$par)), nrow=boot)
set.seed(1111)

for(i in 1:boot){
  samp2 <- ind_0515[sample(1:dim(ind_0515)[1], replace=T, dim(ind_0515)[1]),]
  x2 <- optim(start_logit,fn=logit_flike, method="BFGS", age = samp2$age, emp_dummy = samp2$emp_dummy,
              x1 = samp2$year_2006, x2 = samp2$year_2007, x3 = samp2$year_2008, x4 = samp2$year_2009,
              x5 = samp2$year_2010, x6 = samp2$year_2011, x7 = samp2$year_2012, x8 = samp2$year_2013,
              x9 = samp2$year_2014, x10 = samp2$year_2015, hessian = TRUE)
  predict2 <- x2$par[1] + x2$par[2]*samp2$age + x2$par[3]*samp2$year_2006 +
    x2$par[4]*samp2$year_2007 + x2$par[5]*samp2$year_2008 + x2$par[6]*samp2$year_2009 + 
    x2$par[7]*samp2$year_2010 + x2$par[8]*samp2$year_2011 + x2$par[9]*samp2$year_2012 +
    x2$par[10]*samp2$year_2013 + x2$par[11]*samp2$year_2014 + x2$par[12]*samp2$year_2015
  pdf2 <- mean(dlogis(predict2))
  bootvals2[i,] <- pdf2 * x2$par
}

res2 <- cbind(marginal_effect_logit,apply(bootvals2,2,sd),marginal_effect_logit/apply(bootvals2,2,sd))

colnames(res2) <- c("marginal.effect","standard.error","z.ratio")  








