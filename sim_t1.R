#Generate Data
inv_logit = function(x){
  exp(x)/(exp(x)+1)
}

#true blip function used to generate data
full_blip = function(l1,l2,u,psi){
  psi[1]*l1 + psi[2]*l2 + psi[3]*u +psi[4]*l1*l2
}

n=10000

sim_pats = function(N=n,blip_func=full_blip,psi=c(1,1,1,1)){
  U = rnorm(N,0,.1)
  L1 = rnorm(N,U)
  L2 = rnorm(N,U)
  Y0 = rnorm(N,U+L1+L2) 
  Y1_0 = rnorm(N,U+2*L1+1.5*L2)
  A = rbinom(N,1,inv_logit(.5*U+.25*L1+.25*L2))
  Y1 = Y1_0 + rnorm(N,A*mapply(blip_func,L1,L2,U,MoreArgs=list(psi)),.1)
  return(data.frame(cbind(L1,L2,A,Y0,Y1_0,Y1)))
}

sim_pats2 = function(N=n,blip_func=full_blip,psi=c(1,1,1,1)){
  U = rnorm(N,0,.1)
  mu_L = rnorm(N)
  L1 = rnorm(N,mu_L,.1)
  L2 = rnorm(N,mu_L,.1)
  Y0 = rnorm(N,U+L1+L2) 
  Y1_0 = rnorm(N,U+1.5*L1+2*L2)
  A = rbinom(N,1,inv_logit(.5*U+.25*L1+.25*L2))
  Y1 = Y1_0 + rnorm(N,A*mapply(blip_func,L1,L2,U,MoreArgs=list(psi)))
  return(data.frame(cbind(L1,L2,A,Y0,Y1_0,Y1)))
}

data = sim_pats2(N=10000000)

#unadjusted parallel trends does not hold
mean(data$Y1_0[data$A==1]-data$Y0[data$A==1])
mean(data$Y1_0[data$A==0]-data$Y0[data$A==0])

#true ATT
mean(data$Y1[data$A==1]-data$Y1_0[data$A==1])

#With covariates, how does effect vary with L1
#true linear projection
data$real_effect = data$Y1-data$Y1_0
true_lin_blip = lm(real_effect~bs(L1),data=data[data$A==1,])
true_lin_blip2 = lm(real_effect~L1+I(L1^2),data=data[data$A==1,])

summary(true_lin_blip2)

nsims=1000
abadie_ests_list = vector(mode='list',length=nsims)
dr_ests_list = vector(mode='list',length=nsims)
dr_eff_ests_list = vector(mode='list',length=nsims)
att_unadj_ests = rep(NA,nsims)
att_abadie_ests = rep(NA,nsims)
att_dr_ests = rep(NA,nsims)

library(nleqslv)
library(splines)

set.seed(1)

for(i in 1:nsims){

data = sim_pats2(N=n)

#Unadjusted DiD estimate
att_unadj_ests[i]=mean(data$Y1[data$A==1]-data$Y0[data$A==1]) - mean(data$Y1[data$A==0]-data$Y0[data$A==0])

#Adjusted DiD ATT estimate IPW Abadie
# treat_mod = gam(A ~ s(L1)*s(L2),family="binomial",data=data)
treat_mod = glm(A ~ L1*L2+I(L1^2)+I(L2^2),family="binomial",data=data)
# while(!treat_mod$converged){
#   data = sim_pats(N=n)
#   treat_mod = glm(A ~ bs(L1)*bs(L2),family="binomial",data=data)
# }
# summary(treat_mod)
data$A_hat = predict(treat_mod,newdata=data,type="response")

att_abadie_ests[i]=mean(((data$Y1-data$Y0)/mean(data$A))*(data$A-data$A_hat)/(1-data$A_hat))

#Adjusted DiD ATT estimate outcome regression
data$diff = data$Y1 - data$Y0
# treat_mod = gam(I(Y1-Y0) ~ s(L1)*s(L2),family="binomial",data=data)
diff_mod_A0 = lm(diff ~ L1*L2+I(L1^2)+I(L2^2),data=data[data$A==0,])
data$diff_preds = predict(diff_mod_A0,newdata=data)
# mean(data$diff[data$A==1]) - mean(data$diff_preds[data$A==1])


#Doubly robust ATT estimate no covariates
est_eq = function(psi){
  sum(((data$A_hat-data$A)/(1-data$A_hat))*(data$diff_preds-data$diff)-data$A*psi)
}

ss = nleqslv(x=0, fn=est_eq)
psi_hat = ss$x
att_dr_ests[i]=psi_hat


#Abadie heterogeneity
data$rho = (data$A-data$A_hat)/(data$A_hat*(1-data$A_hat))
# spline_mat = bs(data$L1)
# data = cbind(data,spline_mat)
# names(data)[12:14]=c('s1','s2','s3')
data$int=1
# X = as.matrix(data[,c('int','s1','s2','s3')])
X2 = as.matrix(cbind(data$int,data$L1,data$L1^2))
abadie_est_eq = function(theta){
  mean(data$A_hat*(data$rho*data$diff-X2%*%theta)^2)
}
abadie_ests_list[[i]] = optim(par=rep(0,3),fn=abadie_est_eq)$par


#Us, regression based


#Us, doubly robust not efficient

est_eq_dr = function(psi){
  V1=as.vector(((data$A-data$A_hat)/(1-data$A_hat))*(data$diff-data$diff_preds)-data$A*X2%*%psi)
  colSums(V1*X2)
}
ss = nleqslv(x=rep(0,3), fn=est_eq_dr)
psi_hat_dr = ss$x
# ss$termcd
# ss$fvec
dr_ests_list[[i]]=psi_hat_dr

#Us, doubly robust efficient
A_mod_tilde = glm(A~bs(L1),data=data,family='binomial')
data$A_hat_tilde = predict(A_mod_tilde,newdata=data,type="response")
V1_dr=as.vector(((data$A-data$A_hat)/(1-data$A_hat))*(data$diff-data$diff_preds)-data$A*X2%*%psi_hat_dr)
# V1_dr = pmax(pmin(V1_dr,quantile(V1_dr,.99)),quantile(V1_dr,.01))
data$var_eps = V1_dr^2
v1_mod = lm(var_eps~data$L1+I(data$L1^2),data=data)
data$E_var = pmax(.1,v1_mod$fitted.values)
# data$E_var = pmax(.01,data$E_var)
# v1_resids = v1_mod$residuals^2
# var_eps_mod = lm(v1_resids~bs(data$L1))
# data$E_var = var_eps_mod$fitted.values

est_eq_dr_eff = function(psi){
  V1=as.vector(((data$A-data$A_hat)/(1-data$A_hat))*(data$diff-data$diff_preds)-data$A*X2%*%psi)
  J = data$A_hat_tilde*X2
  colSums(V1*J/data$E_var)
}
ss = nleqslv(x=rep(0,3), fn=est_eq_dr_eff)
psi_hat_dr_eff = ss$x
# ss$termcd
# ss$fvec
dr_eff_ests_list[[i]] = psi_hat_dr_eff
}


sd(sapply(dr_eff_ests_list,function(x)x[1]))
sd(sapply(abadie_ests_list,function(x)x[1]))
sd(sapply(dr_ests_list,function(x)x[1]))

sd(sapply(dr_eff_ests_list,function(x)x[2]))
sd(sapply(abadie_ests_list,function(x)x[2]))
sd(sapply(dr_ests_list,function(x)x[2]))

sd(sapply(dr_eff_ests_list,function(x)x[3]))
sd(sapply(abadie_ests_list,function(x)x[3]))
sd(sapply(dr_ests_list,function(x)x[3]))

mean(sapply(dr_eff_ests_list,function(x)x[1]))
mean(sapply(abadie_ests_list,function(x)x[1]))
mean(sapply(dr_ests_list,function(x)x[1]))

mean(sapply(dr_eff_ests_list,function(x)x[2]))
mean(sapply(abadie_ests_list,function(x)x[2]))
mean(sapply(dr_ests_list,function(x)x[2]))

mean(sapply(dr_eff_ests_list,function(x)x[3]))
mean(sapply(abadie_ests_list,function(x)x[3]))
mean(sapply(dr_ests_list,function(x)x[3]))

