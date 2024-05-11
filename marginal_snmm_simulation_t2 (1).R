#Generate Data
inv_logit = function(x){
  exp(x)/(exp(x)+1)
}

#true blip function used to generate data
blip = function(m,k,l1,a,psi){
  a*sum(c(k-m+1,l1[m])*psi)
}

#number of time-points including time 0 when treatment is always 0
ntimes=3

sim_pat_coarse_ptT = function(times=ntimes,blip_func=blip,psi=c(1,.5)){
  L1 = rep(NA,times)
  L2 = rnorm(1)
  Y = rep(NA,times)
  Y_0 = rep(NA,times)
  A = rep(0,times)
  start=Inf
  U = rnorm(1)
  L1[1] = rnorm(1)
  Y[1] = rnorm(1,U+L1[1]+L2)
  Y_0[1] = Y[1]
  #everyone starts out at 0
  A[1]=0
  #make counterfactual untreated trajectory
  for(t in 2:(times)){
    L1[t] = rnorm(1,L1[t-1])
    Y_0[t] = rnorm(1,L1[t]+L2+U)
    A[t] = rbinom(1,1,inv_logit(-1+L1[t]+L2+U))
  }
  start = min(which(A!=0))
  #if no treatment, just return untreated trajectory
  if(start==Inf){
    return(list(A=A,Y=Y_0,L1=L1,L2=rep(L2,ntimes),Y_0=Y_0))
  }else{
    if(start<times){
    #otherwise blip up the outcomes using the blip function
      A[(start+1):times]=0
    }
    Y= Y_0
    for(k in start:times){
      Y[k] = rnorm(1,Y_0[k] + blip_func(start,k,L1,A[start],psi))
    }
  }
  return(list(A=A,Y=Y,L1=L1,L2=rep(L2,ntimes),Y_0=Y_0))
}
sim_pat_coarse_ptT()
#number of patients to simulate
N = 10000
#make a dataframe of all the observed values for each patient
#id,A,Y,L,time
nsims=1000
psi_hat_or_list = vector(mode='list',length=nsims)
psi_hat_dr_list = vector(mode='list',length=nsims)
psi_hat_dr_wrong_treat_list = vector(mode='list',length=nsims)
psi_hat_dr_wrong_outcome_list = vector(mode='list',length=nsims)
E_Y0_hat_or_list = vector(mode='list',length=nsims)
E_Y0_hat_dr_list = vector(mode='list',length=nsims)
E_Y0_hat_dr_wrong_treat_list = vector(mode='list',length=nsims)
E_Y0_hat_or_wrong_outcome_list = vector(mode='list',length=nsims)


for(i in 1:nsims){
  blip = function(m,k,l1,a,psi){
    a*sum(c(k-m+1,l1[m])*psi)
  }
data = replicate(N,sim_pat_coarse_ptT())
data = data.frame(cbind(A = unlist(data['A',]),Y=unlist(data['Y',]),L1=unlist(data['L1',]),L2=unlist(data['L2',]),Y_0=unlist(data['Y_0',])))
data$id = rep(1:N,each=ntimes)
data$time = rep(0:(ntimes-1),N)

###........................................................................................................................................


#go to wide format
for(var in c('L1','Y','A','Y_0')){
  for(y in 0:2){
    data[,paste0(var,y)] = unlist(lapply(split(data[,var],data$id),function(x)rep(x[y+1],length(x))))
  }
}

data$lag_Y = unlist(lapply(split(data$Y,data$id),function(x)c(NA,x[1:(length(x)-1)])))
#fit the regressions
data$diffs2 = data$Y2 - data$Y1
data$diffs1 = data$Y1 - data$Y0

data$past_A2 = data$A1==1
library(splines)

regression2_1 = lm(diffs2~bs(L12),
                   data=data[data$A2==1 & data$past_A2==0 & data$time==2,])
regression1_1 = lm(diffs1~bs(L11),
                   data=data[data$A1==1 & data$time==1,])

regression2_0 = lm(diffs2~bs(L12)*bs(L2)*bs(L11),
                   data=data[data$A2==0 & data$past_A2==0 & data$time==2,])
regression1_0 = lm(diffs1~bs(L11)*bs(L2),
                   data=data[data$A1==0 & data$time==1,])
data$preds2_0 = predict(regression2_0,newdata=data)


regression12_0 = lm(preds2_0~bs(L11)*bs(L2),
                    data=data[data$A1==0 & data$time==1,])
regression20_1 = lm(preds2_0~bs(L12),
                    data=data[data$A2==1 & data$past_A2==0 & data$time==2,])
data$preds20_1 = predict(regression20_1,newdata = data)
data$preds1_0 = predict(regression1_0,newdata=data)
regression10_1 = lm(preds1_0~bs(L11),
                    data=data[data$A1==1 & data$time==1,])
data$preds10_1 = predict(regression10_1,newdata = data)

data$preds1_20 = predict(regression12_0,newdata=data)
regression1_20_1 = lm(preds1_20~bs(L11),data=data[data$time==1 & data$A==1,])
data$preds1_20_1 = predict(regression1_20_1,newdata = data)

#treatment model
treat_mod1 = glm(A1~bs(L11)*bs(L2),data=data,family='binomial')
treat_mod2 = glm(A2~bs(L11)*bs(L2)*bs(L12),data=data[data$past_A2==0,],family='binomial')

data$A1_hat = predict(treat_mod1,newdata=data,type='response')
data$A2_hat = predict(treat_mod2,newdata=data,type='response')


# gamma_22 = function(L_tilde){
#   #sum(regression2_1$coefficients*c(1,L_tilde)) - sum(regression20_1$coefficients*L_tilde)
#   predict(regression2_1,newdata = data.frame(L12=L_tilde)) - 
#     predict(regression20_1,newdata = data.frame(L12=L_tilde))
# }
# 
# gamma_11 = function(L_tilde){
#   #sum(regression1_1$coefficients*c(1,L_tilde)) - sum(regression10_1$coefficients*L_tilde)
#   predict(regression1_1,newdata=data.frame(L11=L_tilde))-
#     predict(regression10_1,newdata=data.frame(L11=L_tilde))
# }
# 
# data$H11 = data$Y1 - data$A1*unlist(sapply(data$L11,gamma_11))
# data$diffs_blip11 = data$Y2 - data$H11
# regression12_blip = lm(diffs_blip11~bs(L11),data=data[data$A1==1 & data$time==1,])
# 
# gamma_12 = function(L_tilde){
#   # sum(regression12_blip$coefficients*c(1,L_tilde)) - sum(regression1_20_1$coefficients*c(1,L_tilde))
#   predict(regression12_blip,newdata=data.frame(L11=L_tilde))-
#     predict(regression1_20_1,newdata=data.frame(L11=L_tilde))
# }
# 
# #can then provide blip outputs for grid of inputs of interest
# 
# gamma_12(1)
# gamma_12(0)
# 
# gamma_22(1)
# gamma_22(0)

# data$real_effects2 = data$Y2-data$Y_02
# real_effects2_mod = lm(real_effects2~bs(L12),data=data[data$A2==1 & data$past_A2==0,])
# predict(real_effects2_mod,newdata=data.frame(L12=1))
# predict(real_effects2_mod,newdata=data.frame(L12=0))
# 
# mean(data$Y_02[data$A2==1 & data$past_A2==0 & data$time==1]-data$Y_01[data$A2==1 & data$past_A2==0 & data$time==1])
# mean(data$Y_02[data$A2==0 & data$past_A2==0 & data$time==1]-data$Y_01[data$A2==0 & data$past_A2==0 & data$time==1])
# 
# data$trend21 = data$Y_02-data$Y_01
# trend21_mod_A1 = lm(trend21~bs(L12)*bs(L11)*bs(L2),data=data[data$A1==0 & data$A2==1 & data$time==2,])
# trend21_mod_A0 = lm(trend21~bs(L12)*bs(L11)*bs(L2),data=data[data$A1==0 & data$A2==0 & data$time==2,])
# predict(trend21_mod_A0,data.frame(L11=1,L12=2,L2=1))
# predict(trend21_mod_A1,data.frame(L11=1,L12=2,L2=1))

# gamma_11(1)
# gamma_11(0)


#estimating equation style with parametric blip models

blip = function(m,k,l_m,psi){
  (psi[1] + psi[2]*l_m)*(m==1)*(k==1) + (psi[3]+psi[4]*l_m)*(m==1)*(k==2) + (psi[5]+psi[6]*l_m)*(m==2)*(k==2)
}

data$first_treat = apply(data[,c('A1','A2')],1,function(x)min(which(x==1)))
data$L1_treat = ifelse(data$first_treat==1,data$L11,
                                   ifelse(data$first_treat==2,data$L12,
                                          0))
mk = expand.grid(m=1:2,k=1:2)
mk = mk[mk$k>=(mk$m-1),]
mk = mk[order(mk$m),]
Hmk = data.frame(id=rep(unique(data$id),each=nrow(mk)),m=rep(mk$m,length(unique(data$id))),k=rep(mk$k,length(unique(data$id))))
data$k = data$time
data$m = data$time
Hmk = merge(Hmk,data[,c("id","Y","lag_Y","k")])
Hmk = merge(Hmk,data[,c("id","A","L1","L2","m","preds1_20","preds2_0","preds1_0")])
Hmk = merge(Hmk,data[data$time==1,c("id",'A1','A2',"L11","L12",'Y1','Y2','first_treat','L1_treat','A1_hat','A2_hat')])

Hmk$q1 = (Hmk$m==1)*(Hmk$k==1)
Hmk$q2 = Hmk$L1*(Hmk$m==1)*(Hmk$k==1)
Hmk$q3 = (Hmk$m==1)*(Hmk$k==2)
Hmk$q4 = Hmk$L1*(Hmk$m==1)*(Hmk$k==2)
Hmk$q5 = (Hmk$m==2)*(Hmk$k==2)
Hmk$q6 = Hmk$L1*(Hmk$m==2)*(Hmk$k==2)

est_eq_or = function(psi_hat,estmat=Hmk){
  estmat$blips_k = ifelse(estmat$k>=estmat$first_treat & estmat$m<=estmat$first_treat, 
                          sapply(1:nrow(estmat),function(i)blip(m=estmat$first_treat[i],k=estmat$k[i],l_m=estmat$L1_treat[i],psi=psi_hat)),
                          0)
  estmat$Hk = estmat$Y - estmat$blips_k
  estmat$blips_k_minus_1 = ifelse((estmat$k-1)>=estmat$first_treat & estmat$m<=estmat$first_treat, 
                                  sapply(1:nrow(estmat),function(i)blip(m=estmat$first_treat[i],k=estmat$k[i]-1,l_m=estmat$L1_treat[i],psi=psi_hat)),0)
  estmat$Hk_lag = estmat$lag_Y - estmat$blips_k_minus_1
  estmat$E_H_diff = ifelse(estmat$m==1 & estmat$k==1,estmat$preds1_0,
                           ifelse(estmat$m==1 & estmat$k==2,estmat$preds1_20,estmat$preds2_0))
  estmat = estmat[estmat$m<=estmat$first_treat & estmat$m<=estmat$k,]
  apply(estmat[,c('q1','q2','q3','q4','q5','q6')]*(estmat$Hk-estmat$Hk_lag-estmat$E_H_diff),2,sum)
}

library(nleqslv)
ss_or = nleqslv(x=rep(0,6),fn=est_eq_or)
psi_hat_or = ss_or$x
ss_or$termcd
ss_or$fvec
psi_hat_or_list[[i]]=psi_hat_or

Hmk$blips_k_or = ifelse(Hmk$k>=Hmk$first_treat & Hmk$m<=Hmk$first_treat, 
                        sapply(1:nrow(Hmk),function(i)blip(m=Hmk$first_treat[i],k=Hmk$k[i],l_m=Hmk$L1_treat[i],psi=psi_hat_or)),
                        0)
Hmk$Hk = Hmk$Y - Hmk$blips_k_or
E_Y0_hat_or_list[[i]] = c(mean(Hmk$Hk[Hmk$m==1 & Hmk$k==2]), mean(Hmk$Hk[Hmk$m==1 & Hmk$k==2]))


data$A2 = ifelse(data$A1==1,0,data$A2)
data$A2_hat[data$A1==1] = 0
blip22 = function(l_m,psi22){
  psi22[1]+psi22[2]*l_m 
}

blip12 = function(l_m,psi12){
  psi12[1]+psi12[2]*l_m 
}

blip11 = function(l_m,psi11){
  psi11[1]+psi11[2]*l_m 
}

data$q11_1 = 1
data$q11_2 = data$L11

data$q22_1 = 1
data$q22_2 = data$L12

data$q12_1 = 1
data$q12_2 = data$L11


est_eq_dr = function(psi_hat,estmat=data){
  psi_hat11 = psi_hat[1:2]
  psi_hat12 = psi_hat[3:4]
  psi_hat22 = psi_hat[5:6]
  epsilon_22 = ((estmat$A2-estmat$A2_hat)/(1-estmat$A2_hat))*(estmat$Y2-estmat$Y1-estmat$preds2_0)-estmat$A2*sapply(estmat$L12,blip22,psi22=psi_hat22)
  epsilon_11 = ((estmat$A1-estmat$A1_hat)/(1-estmat$A1_hat))*(estmat$Y1-estmat$Y0-estmat$preds1_0)-estmat$A1*sapply(estmat$L11,blip11,psi11=psi_hat11)
  epsilon_12 = (1-((1-estmat$A2)*(1-estmat$A1))/((1-estmat$A2_hat)*(1-estmat$A1_hat)))*(estmat$Y2-estmat$Y1)-
    ((1-estmat$A1)/(1-estmat$A1_hat))*(1-(1-estmat$A2)/(1-estmat$A2_hat))*estmat$preds2_0 -
    (1-(1-estmat$A1)/(1-estmat$A1_hat))*estmat$preds1_20 - (estmat$A2*sapply(1:nrow(estmat),function(i) blip22(l_m=estmat$L12[i],psi_hat22)) -
                                                              estmat$A1*sapply(1:nrow(estmat),function(i) blip11(l_m=estmat$L11[i],psi_hat11))) - 
    estmat$A1*sapply(1:nrow(estmat),function(i) blip12(l_m=estmat$L11[i],psi_hat12))
  q11 = data[,c('q11_1','q11_2')]
  q12 = data[,c('q12_1','q12_2')]
  q22 = data[,c('q22_1','q22_2')]
  c(apply(q12*(epsilon_12-q22*epsilon_22+q11*epsilon_11),2,sum),apply(q22*epsilon_22,2,sum),apply(q11*epsilon_11,2,sum))
}

ss_dr = nleqslv(x=rep(0,6),fn=est_eq_dr)
psi_hat_dr = ss_dr$x
ss_dr$termcd
ss_dr$fvec
psi_hat_dr_list[[i]]=psi_hat_dr

Hmk$blips_k_dr = ifelse(Hmk$k>=Hmk$first_treat & Hmk$m<=Hmk$first_treat, 
                        sapply(1:nrow(Hmk),function(i)blip(m=Hmk$first_treat[i],k=Hmk$k[i],l_m=Hmk$L1_treat[i],psi=psi_hat_dr)),
                        0)
Hmk$Hk = Hmk$Y - Hmk$blips_k_dr
E_Y0_hat_dr_list[[i]] = c(mean(Hmk$Hk[Hmk$m==1 & Hmk$k==2]), mean(Hmk$Hk[Hmk$m==1 & Hmk$k==2]))


#Now do dr estimator but with incorrect treatment model
data$A2_hat = ifelse(data$past_A2==1,0,.5)
data$A1_hat = .5
ss_dr = nleqslv(x=rep(0,6),fn=est_eq_dr)
psi_hat_dr_wrong_treat = ss_dr$x
ss_dr$termcd
ss_dr$fvec
psi_hat_dr_wrong_treat_list[[i]]=psi_hat_dr_wrong_treat

#now set the treatment model to be correct again and make the outcome model wrong
data$A1_hat = predict(treat_mod1,newdata=data,type='response')
data$A2_hat = predict(treat_mod2,newdata=data,type='response')
data$A2_hat[data$A1==1] = 0
data$preds2_0=.5
data$preds1_0=.5
data$preds1_20 = .5
ss_dr = nleqslv(x=rep(0,6),fn=est_eq_dr)
psi_hat_dr_wrong_outcome = ss_dr$x
ss_dr$termcd
ss_dr$fvec
psi_hat_dr_wrong_outcome_list[[i]]=psi_hat_dr_wrong_outcome
}

psi_hat_or_list
psi_hat_dr_list
psi_hat_dr_wrong_outcome_list
psi_hat_dr_wrong_treat_list

