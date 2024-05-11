library(haven)
library(labelled)
library(readxl)
mgi = ses <- read_excel("~/snmm/medicaid/MGI_2000_2020_Ver1.9_20201104_imputed.xlsx",sheet="Income Eligibility")
mgi$thresh = (mgi$WG_Medicaid_Eligibility_06 + mgi$WG_Medicaid_Eligibility_07)/2
names(mgi)[3] = 'stname'
names(mgi)[1] = 'year'
ses <- read_dta("~/snmm/medicaid/master_area_dataset_4MEPS.dta")
sesdic <- labelled::generate_dictionary(ses)

#variables 1-42, state identification and exposure identification
#82-96, state and county uninsured with all, <=135FPL and 400FPL
#1440, % poverty up to 2017

#year manipulation
#log of outcome
# impute missing with median value
#get the year of treatment initiation
# #get first year that each county deregulated
#    first_dereg_year = aggregate(data$inter_bra,by=list(data$county),function(x) min(which(x>0))) #*MS** minimum year for deregulation by county

###........................................................................................................................................

#SETUP
#Data subsetting to exposure, outcome, time, and ID vars

###------------------------------------------------------------------------------------------------------------------------------------------
library(tidyverse)

ses_sub <- ses[,c('stname', 'stabb', 'stfips',  #state codes
                  'coname', 'coname_stabb', 'cofips', #county
                  'stcofips', 'stcofips_num', #state+county
                  'year', 'st_acaexp', #time and A
                  'co_sahie1864_insn_138fpl', 'co_sahie1864_unipct_138fpl', #Y
                  'co_tot_pop', 'co_ruc_code13','co_per_black','co_per_hisp')] #L
ses_sub = merge(ses_sub,mgi[,c("stname","thresh","year")],all.x=T)

ses_sub$log_pop = log(ses_sub$co_tot_pop)

ses_sub <- ses_sub %>% mutate(across(everything(), as.vector))
str(ses_sub)    

log_pop2013 = ses_sub[ses_sub$year==2013,c('stcofips','log_pop')]
names(log_pop2013)[2] = 'log_pop2013'
ses_sub = merge(ses_sub,log_pop2013)
ses_sub$log_pop2013[is.na(ses_sub$log_pop2013)] = median(ses_sub$log_pop2013,na.rm=T)

co_per_black2013 = ses_sub[ses_sub$year==2013,c('stcofips','co_per_black')]
names(co_per_black2013)[2] = 'co_per_black2013'
ses_sub = merge(ses_sub,co_per_black2013)
ses_sub$co_per_black2013[is.na(ses_sub$co_per_black2013)] = median(ses_sub$co_per_black2013,na.rm=T)

co_per_hisp2013 = ses_sub[ses_sub$year==2013,c('stcofips','co_per_hisp')]
names(co_per_hisp2013)[2] = 'co_per_hisp2013'
ses_sub = merge(ses_sub,co_per_hisp2013)
ses_sub$co_per_hisp2013[is.na(ses_sub$co_per_hisp2013)] = median(ses_sub$co_per_hisp2013,na.rm=T)


ses_sub = ses_sub[order(ses_sub$stcofips,ses_sub$year),]

ses_sub$last_thresh = unlist(lapply(split(ses_sub$thresh,ses_sub$stcofips),function(x)c(NA,x[1:(length(x)-1)])))
ses_sub$last_thresh[is.na(ses_sub$last_thresh)] = median(ses_sub$last_thresh,na.rm=T)
ses_sub$thresh[is.na(ses_sub$thresh)] = median(ses_sub$thresh,na.rm=T)

L = c('log_pop2013', 'co_ruc_code13','co_per_black2013','co_per_hisp2013','last_thresh')
L_tilde = c('last_thresh','co_ruc_code13')

ses_sub$co_ruc_code13 = as.numeric(ses_sub$co_ruc_code13)

###........................................................................................................................................

#change time so that 2013=0 and 2019 = 6
ses_sub$yearms <- ses_sub$year - 2013

#get data starting from 2013 onwards
ses_sub <- ses_sub[ses_sub$yearms >=0 & ses_sub$yearms<=2,]  


#we imputed with the population median for missing Y
medimp <- median(ses_sub$co_sahie1864_unipct_138fpl, na.rm = T)

ses_sub$co_sahie1864_unipct_138fpl <- ifelse(is.na(ses_sub$co_sahie1864_unipct_138fpl), medimp, ses_sub$co_sahie1864_unipct_138fpl)
#check for nas
sum(is.na(ses_sub$co_sahie1864_unipct_138fpl))

ses_sub$co_ruc_code13 = as.numeric(ses_sub$co_ruc_code13)
ses_sub$co_ruc_code13[is.na(ses_sub$co_ruc_code13)]=median(ses_sub$co_ruc_code13,na.rm=T)

ses_sub$Y = ses_sub$co_sahie1864_unipct_138fpl

ses_sub$A = ses_sub$st_acaexp
#go to wide format
for(var in c(L,'Y','A')){
  for(y in 0:2){
    ses_sub[,paste0(var,y)] = unlist(lapply(split(ses_sub[,var],ses_sub$stcofips),function(x)rep(x[y+1],length(x))))
  }
}

ses_sub$lag_Y = unlist(lapply(split(ses_sub$Y,ses_sub$stcofips),function(x)c(NA,x[1:(length(x)-1)])))
#fit the regressions
ses_sub$diffs2 = ses_sub$Y2 - ses_sub$Y1
ses_sub$diffs1 = ses_sub$Y1 - ses_sub$Y0

ses_sub$past_A2 = ses_sub$A1==1

library(splines)
regression2_1 = lm(diffs2~last_thresh2+co_ruc_code132,
                   data=ses_sub[ses_sub$A2==1 & ses_sub$past_A2==0 & ses_sub$yearms==2,])
regression1_1 = lm(diffs1~last_thresh1+co_ruc_code131,
                   data=ses_sub[ses_sub$A1==1 & ses_sub$yearms==1,])

regression2_0 = lm(diffs2~co_per_black20132+co_per_hisp20132+log_pop20132+co_ruc_code132+last_thresh2,
                   data=ses_sub[ses_sub$A2==0 & ses_sub$yearms==2,])
regression1_0 = lm(diffs1~co_per_black20131+co_per_hisp20131+log_pop20131+co_ruc_code131+last_thresh1,
                   data=ses_sub[ses_sub$A1==0 & ses_sub$yearms==1,])
ses_sub$preds2_0 = predict(regression2_0,newdata=ses_sub)


regression12_0 = lm(preds2_0~co_per_black20131+co_per_hisp20131+log_pop20131+co_ruc_code131+last_thresh1,
                    data=ses_sub[ses_sub$A1==0 & ses_sub$yearms==1,])
regression20_1 = lm(preds2_0~last_thresh2+co_ruc_code132,
                    data=ses_sub[ses_sub$A2==1 & ses_sub$past_A2==0 & ses_sub$yearms==2,])
ses_sub$preds20_1 = predict(regression20_1,newdata = ses_sub)
ses_sub$preds1_0 = predict(regression1_0,newdata=ses_sub)
regression10_1 = lm(preds1_0~last_thresh1+co_ruc_code131,
                    data=ses_sub[ses_sub$A1==1 & ses_sub$yearms==1,])
ses_sub$preds10_1 = predict(regression10_1,newdata = ses_sub)

ses_sub$preds1_20 = predict(regression12_0,newdata=ses_sub)
regression1_20_1 = lm(preds1_20~last_thresh1+co_ruc_code131,data=ses_sub[ses_sub$yearms==1 & ses_sub$A==1,])
ses_sub$preds1_20_1 = predict(regression1_20_1,newdata = ses_sub)

#treatment model
treat_mod1 = glm(A1~co_per_black20131+co_per_hisp20131+log_pop20131+co_ruc_code131+last_thresh1,data=ses_sub,family='binomial')
treat_mod2 = glm(A2~co_per_black20132+co_per_hisp20132+log_pop20132+co_ruc_code132+last_thresh2,data=ses_sub[ses_sub$past_A2==0,],family='binomial')

treat_mod = glm(A1~log_pop20131+co_ruc_code131+last_thresh1,data=ses_sub,family='binomial')

ses_sub$A1_hat = predict(treat_mod1,newdata=ses_sub,type='response')
ses_sub$A1_hat_trim = pmin(.95,ses_sub$A1_hat)
ses_sub$A2_hat = predict(treat_mod2,newdata=ses_sub,type='response')
ses_sub$A2_hat_trim = pmin(.95,ses_sub$A2_hat)


gamma_22 = function(L_tilde1,L_tilde2){
  #sum(regression2_1$coefficients*c(1,L_tilde)) - sum(regression20_1$coefficients*L_tilde)
  predict(regression2_1,newdata = data.frame(last_thresh2=L_tilde1,co_ruc_code132=L_tilde2)) - 
    predict(regression20_1,newdata = data.frame(last_thresh2=L_tilde1,co_ruc_code132=L_tilde2))
}

gamma_11 = function(L_tilde1,L_tilde2){
  #sum(regression1_1$coefficients*c(1,L_tilde)) - sum(regression10_1$coefficients*L_tilde)
  predict(regression1_1,newdata=data.frame(last_thresh1=L_tilde1,co_ruc_code131=L_tilde2))-
    predict(regression10_1,newdata=data.frame(last_thresh1=L_tilde1,co_ruc_code131=L_tilde2))
}

ses_sub$H11 = ses_sub$Y1 - ses_sub$A1*unlist(mapply(gamma_11,ses_sub$last_thresh1,ses_sub$co_ruc_code131))
ses_sub$diffs_blip11 = ses_sub$Y2 - ses_sub$H11
regression12_blip = lm(diffs_blip11~last_thresh1+co_ruc_code131,data=ses_sub[ses_sub$A1==1 & ses_sub$yearms==1,])

gamma_12 = function(L_tilde1,L_tilde2){
  # sum(regression12_blip$coefficients*c(1,L_tilde)) - sum(regression1_20_1$coefficients*c(1,L_tilde))
  predict(regression12_blip,newdata=data.frame(last_thresh1=L_tilde1,co_ruc_code131=L_tilde2))-
    predict(regression1_20_1,newdata=data.frame(last_thresh1=L_tilde1,co_ruc_code131=L_tilde2))
}

#can then provide blip outputs for grid of inputs of interest

gamma_12(.1,9)
gamma_22(.1,9)
gamma_11(.1,9)


#estimating equation style with parametric blip models

blip = function(m,k,l_m1,l_m2,psi){
  (psi[1] + psi[2]*l_m1 +psi[3]*l_m2)*(m==1)*(k==1) + (psi[4]+psi[5]*l_m1+psi[6]*l_m2)*(m==1)*(k==2) + (psi[7]+psi[8]*l_m1+psi[9]*l_m2)*(m==2)*(k==2)
}

ses_sub$first_treat = apply(ses_sub[,c('A1','A2')],1,function(x)min(which(x==1)))
ses_sub$last_thresh_treat = ifelse(ses_sub$first_treat==1,ses_sub$last_thresh1,
                                   ifelse(ses_sub$first_treat==2,ses_sub$last_thresh2,
                                          0))
ses_sub$co_ruc13_treat = ifelse(ses_sub$first_treat==1,ses_sub$co_ruc_code131,
                                   ifelse(ses_sub$first_treat==2,ses_sub$co_ruc_code132,
                                          0))
mk = expand.grid(m=1:2,k=1:2)
mk = mk[mk$k>=(mk$m-1),]
mk = mk[order(mk$m),]
Hmk = data.frame(stcofips=rep(unique(ses_sub$stcofips),each=nrow(mk)),m=rep(mk$m,length(unique(ses_sub$stcofips))),k=rep(mk$k,length(unique(ses_sub$stcofips))))
ses_sub$k = ses_sub$yearms
ses_sub$m = ses_sub$yearms
Hmk = merge(Hmk,ses_sub[,c("stcofips","Y","lag_Y","k")])
Hmk = merge(Hmk,ses_sub[,c("stcofips","A",L,"m","preds1_20","preds2_0","preds1_0")])
Hmk = merge(Hmk,ses_sub[ses_sub$yearms==1,c("stcofips",'A1','A2',paste0(L,1),paste0(L,2),'Y1','Y2','first_treat','last_thresh_treat','co_ruc13_treat','A1_hat','A2_hat')])

#blip model
Hmk$q1 = (Hmk$m==1)*(Hmk$k==1)
Hmk$q2 = Hmk$last_thresh*(Hmk$m==1)*(Hmk$k==1)
Hmk$q3 = Hmk$co_ruc_code13*(Hmk$m==1)*(Hmk$k==1)
Hmk$q4 = (Hmk$m==1)*(Hmk$k==2)
Hmk$q5 = Hmk$last_thresh*(Hmk$m==1)*(Hmk$k==2)
Hmk$q6 = Hmk$co_ruc_code13*(Hmk$m==1)*(Hmk$k==2)
Hmk$q7 = (Hmk$m==2)*(Hmk$k==2)
Hmk$q8 = Hmk$last_thresh*(Hmk$m==2)*(Hmk$k==2)
Hmk$q9 = Hmk$co_ruc_code13*(Hmk$m==2)*(Hmk$k==2)


est_eq_or = function(psi_hat,estmat=Hmk){
  
  estmat$blips_k = ifelse(estmat$k>=estmat$first_treat & estmat$m<=estmat$first_treat, 
                          sapply(1:nrow(estmat),function(i)blip(m=estmat$first_treat[i],k=estmat$k[i],l_m1=estmat$last_thresh_treat[i],l_m2=estmat$co_ruc13_treat[i],psi=psi_hat)),
                          0)
  estmat$Hk = estmat$Y - estmat$blips_k
  estmat$blips_k_minus_1 = ifelse((estmat$k-1)>=estmat$first_treat & estmat$m<=estmat$first_treat, 
                                  sapply(1:nrow(estmat),function(i)blip(m=estmat$first_treat[i],k=estmat$k[i]-1,l_m1=estmat$last_thresh_treat[i],l_m2=estmat$co_ruc13_treat[i],psi=psi_hat)),0)
  estmat$Hk_lag = estmat$lag_Y - estmat$blips_k_minus_1
  estmat$E_H_diff = ifelse(estmat$m==1 & estmat$k==1,estmat$preds1_0,
                           ifelse(estmat$m==1 & estmat$k==2,estmat$preds1_20,estmat$preds2_0))
  estmat = estmat[estmat$m<=estmat$first_treat & estmat$m<=estmat$k,]
  apply(estmat[,c('q1','q2','q3','q4','q5','q6','q7','q8','q9')]*(estmat$Hk-estmat$Hk_lag-estmat$E_H_diff),2,sum)
}

library(nleqslv)
ss_or = nleqslv(x=rep(0,9),fn=est_eq_or)
psi_hat_or = ss_or$x
ss_or$termcd
ss_or$fvec
psi_hat_or


ses_sub$A2 = ifelse(ses_sub$A1==1,0,ses_sub$A2)
ses_sub$A2_hat[ses_sub$A1==1] = 0
ses_sub$A2_hat_trim[ses_sub$A1==1] = 0

blip22 = function(l_m1,l_m2,psi22){
  psi22[1]+psi22[2]*l_m1 + psi22[3]*l_m2 
}

blip12 = function(l_m1,l_m2,psi12){
  psi12[1]+psi12[2]*l_m1+psi12[3]*l_m2 
}

blip11 = function(l_m1,l_m2,psi11){
  psi11[1]+psi11[2]*l_m1+psi11[3]*l_m2 
}

ses_sub$q11_1 = 1
ses_sub$q11_2 = ses_sub$last_thresh1
ses_sub$q11_3 = ses_sub$co_ruc_code131


ses_sub$q22_1 = 1
ses_sub$q22_2 = ses_sub$last_thresh2
ses_sub$q22_3 = ses_sub$co_ruc_code132

ses_sub$q12_1 = 1
ses_sub$q12_2 = ses_sub$last_thresh1
ses_sub$q12_3 = ses_sub$co_ruc_code131


est_eq_dr = function(psi_hat,estmat=ses_sub){
  psi_hat11 = psi_hat[1:3]
  psi_hat12 = psi_hat[4:6]
  psi_hat22 = psi_hat[7:9]
  epsilon_22 = (estmat$A1==0)*(((estmat$A2-estmat$A2_hat)/(1-estmat$A2_hat))*(estmat$Y2-estmat$Y1-estmat$preds2_0)-estmat$A2*mapply(blip22,estmat$last_thresh2,estmat$co_ruc_code132,list(psi22=psi_hat22)))
  epsilon_11 = ((estmat$A1-estmat$A1_hat)/(1-estmat$A1_hat))*(estmat$Y1-estmat$Y0-estmat$preds1_0)-estmat$A1*mapply(blip11,estmat$last_thresh1,estmat$co_ruc_code131,list(psi11=psi_hat11))
  epsilon_12 = (1-((1-estmat$A2)*(1-estmat$A1))/((1-estmat$A2_hat)*(1-estmat$A1_hat)))*(estmat$Y2-estmat$Y1)-
    ((1-estmat$A1)/(1-estmat$A1_hat))*(1-(1-estmat$A2)/(1-estmat$A2_hat))*estmat$preds2_0 -
    (1-(1-estmat$A1)/(1-estmat$A1_hat))*estmat$preds1_20 - (estmat$A2*sapply(1:nrow(estmat),function(i) blip22(l_m1=estmat$last_thresh2[i],l_m2=estmat$co_ruc_code132[i],psi_hat22)) -
                                                              estmat$A1*sapply(1:nrow(estmat),function(i) blip11(l_m1=estmat$last_thresh1[i],l_m2=estmat$co_ruc_code131[i],psi_hat11))) - 
    estmat$A1*sapply(1:nrow(estmat),function(i) blip12(l_m1=estmat$last_thresh1[i],l_m2=estmat$co_ruc_code131[i],psi_hat12))
  q11 = ses_sub[,c('q11_1','q11_2','q11_3')]
  q12 = ses_sub[,c('q12_1','q12_2','q12_3')]
  q22 = ses_sub[,c('q22_1','q22_2','q22_3')]
  c(apply(q12*(epsilon_12-q22*epsilon_22+q11*epsilon_11),2,sum),apply(q22*epsilon_22,2,sum),apply(q11*epsilon_11,2,sum))
}

ss_dr = nleqslv(x=rep(0,9),fn=est_eq_dr)
psi_hat_dr = ss_dr$x
ss_dr$termcd
ss_dr$fvec
psi_hat_dr

est_eq_dr_trim = function(psi_hat,estmat=ses_sub){
  psi_hat11 = psi_hat[1:3]
  psi_hat12 = psi_hat[4:6]
  psi_hat22 = psi_hat[7:9]
  epsilon_22 = (estmat$A1==0)*(((estmat$A2-estmat$A2_hat_trim)/(1-estmat$A2_hat_trim))*(estmat$Y2-estmat$Y1-estmat$preds2_0)-estmat$A2*mapply(blip22,estmat$last_thresh2,estmat$co_ruc_code132,list(psi22=psi_hat22)))
  epsilon_11 = ((estmat$A1-estmat$A1_hat_trim)/(1-estmat$A1_hat_trim))*(estmat$Y1-estmat$Y0-estmat$preds1_0)-estmat$A1*mapply(blip11,estmat$last_thresh1,estmat$co_ruc_code131,list(psi11=psi_hat11))
  epsilon_12 = (1-((1-estmat$A2)*(1-estmat$A1))/((1-estmat$A2_hat_trim)*(1-estmat$A1_hat_trim)))*(estmat$Y2-estmat$Y1)-
    ((1-estmat$A1)/(1-estmat$A1_hat_trim))*(1-(1-estmat$A2)/(1-estmat$A2_hat_trim))*estmat$preds2_0 -
    (1-(1-estmat$A1)/(1-estmat$A1_hat_trim))*estmat$preds1_20 - (estmat$A2*sapply(1:nrow(estmat),function(i) blip22(l_m1=estmat$last_thresh2[i],l_m2=estmat$co_ruc_code132[i],psi_hat22)) -
                                                              estmat$A1*sapply(1:nrow(estmat),function(i) blip11(l_m1=estmat$last_thresh1[i],l_m2=estmat$co_ruc_code131[i],psi_hat11))) - 
    estmat$A1*sapply(1:nrow(estmat),function(i) blip12(l_m1=estmat$last_thresh1[i],l_m2=estmat$co_ruc_code131[i],psi_hat12))
  q11 = ses_sub[,c('q11_1','q11_2','q11_3')]
  q12 = ses_sub[,c('q12_1','q12_2','q12_3')]
  q22 = ses_sub[,c('q22_1','q22_2','q22_3')]
  c(apply(q12*(epsilon_12-q22*epsilon_22+q11*epsilon_11),2,sum),apply(q22*epsilon_22,2,sum),apply(q11*epsilon_11,2,sum))
}

ss_dr_trim = nleqslv(x=rep(0,9),fn=est_eq_dr_trim)
psi_hat_dr_trim = ss_dr_trim$x
ss_dr_trim$termcd
ss_dr_trim$fvec
psi_hat_dr_trim

#compare predictions from different estimators
#implicit blip specification via parametric outcome regressions
gamma_11(1)
gamma_12(1)
gamma_22(1)


#g-estimation with just outcome regression
psi_hat_or
psi_hat_or[1]+.1*psi_hat_or[2]+9*psi_hat_or[3]
psi_hat_or[1]+1*psi_hat_or[2]+1*psi_hat_or[3]
psi_hat_or[1]+.1*psi_hat_or[2]+1*psi_hat_or[3]


psi_hat_or[4]+.1*psi_hat_or[5]+9*psi_hat_or[6]
psi_hat_or[4]+1*psi_hat_or[5]+1*psi_hat_or[6]
psi_hat_or[4]+.1*psi_hat_or[5]+1*psi_hat_or[6]

psi_hat_or[7]+.1*psi_hat_or[8]+9*psi_hat_or[9]
psi_hat_or[7]+1*psi_hat_or[8]+1*psi_hat_or[9]
psi_hat_or[7]+.1*psi_hat_or[8]+1*psi_hat_or[9]

#doubly robust g-estimation untrimmed
psi_hat_dr
psi_hat_dr[1]+.1*psi_hat_dr[2]+9*psi_hat_dr[3]
psi_hat_dr[1]+1*psi_hat_dr[2]+1*psi_hat_dr[3]
psi_hat_dr[1]+.1*psi_hat_dr[2]+1*psi_hat_dr[3]


psi_hat_dr[4]+.1*psi_hat_dr[5]+9*psi_hat_dr[6]
psi_hat_dr[4]+1*psi_hat_dr[5]+1*psi_hat_dr[6]
psi_hat_dr[4]+.1*psi_hat_dr[5]+1*psi_hat_dr[6]

psi_hat_dr[7]+.1*psi_hat_dr[8]+9*psi_hat_dr[9]
psi_hat_dr[7]+1*psi_hat_dr[8]+1*psi_hat_dr[9]
psi_hat_dr[7]+.1*psi_hat_dr[8]+1*psi_hat_dr[9]

#doubly robust g-estimation trimmed
psi_hat_dr_trim
psi_hat_dr_trim[1]+.1*psi_hat_dr_trim[2]+9*psi_hat_dr_trim[3]
psi_hat_dr_trim[1]+1*psi_hat_dr_trim[2]+1*psi_hat_dr_trim[3]
psi_hat_dr_trim[1]+.1*psi_hat_dr_trim[2]+1*psi_hat_dr_trim[3]


psi_hat_dr_trim[4]+.1*psi_hat_dr_trim[5]+9*psi_hat_dr_trim[6]
psi_hat_dr_trim[4]+1*psi_hat_dr_trim[5]+1*psi_hat_dr_trim[6]
psi_hat_dr_trim[4]+.1*psi_hat_dr_trim[5]+1*psi_hat_dr_trim[6]

psi_hat_dr_trim[7]+.1*psi_hat_dr_trim[8]+9*psi_hat_dr_trim[9]
psi_hat_dr_trim[7]+1*psi_hat_dr_trim[8]+1*psi_hat_dr_trim[9]
psi_hat_dr_trim[7]+.1*psi_hat_dr_trim[8]+1*psi_hat_dr_trim[9]

(mean(ses_sub$Y2[ses_sub$A1==0&ses_sub$A2==1])-mean(ses_sub$Y1[ses_sub$A1==0&ses_sub$A2==1]))-
  (mean(ses_sub$Y2[ses_sub$A1==0&ses_sub$A2==0])-mean(ses_sub$Y1[ses_sub$A1==0&ses_sub$A2==0]))

#write general code
#for j in 0:K
#   for k in K:1
#     identify gamma_(k-j)k
max_times = 2
for(j in 0:max_times){
  for(k in max_times:1){
    if(j<k){
      #tilde regression of Y_k-Y_(k-1)(\infty) where first start is k-j
    }
  }
}




