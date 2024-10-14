library(dplyr)
library(rpact)
library(magrittr)
library(tidyverse)
library(MASS)
library(ltm)
library(gMCP)


stage1_data_generation <- function( rho, n, mu_control, delta, cov.matrix, 
                                    thres, sig_level,data=NULL, w_1, w_2, var_type
){
  
  # sample size
  n_S_A_1=n[["S_A_1"]]  
  n_S_B_1=n[["S_B_1"]]
  n_C_A_1=n[["C_A_1"]] 
  n_C_B_1=n[["C_B_1"]] 
  n_F_A_1=n[["stage1_A"]] 
  n_F_B_1=n[["stage1_B"]] 
  
  # generate correlated data (bivariate normal)
  df=sim_F_CPL(delta, mu_control, cov.matrix, n_S_B_1, n_S_A_1, n_C_B_1, n_C_A_1,var_type)
  
  # number of responders for binary outcomes and means for continuous outcomes 
  # (A is control, B is intervention)
  mean_S_C=df%>%
    pivot_longer(contains("V"),names_to = "variable",values_to = "values")%>%
    group_by(variable,gp,trt)%>%
    summarise(
      mean=mean(values),
      n.responder=sum(values),
      n=n()
    )
  
  mean_F=df%>%
    pivot_longer(contains("V"),names_to = "variable",values_to = "values")%>%
    group_by(variable,trt)%>%
    summarise(
      mean=mean(values),
      n.responder=sum(values),
      n=n()
    )
  
  # estimated delta and p-value in full pop and subgroup from stage 1
  data <- data.frame(
    # Stage 1 sample size 
    n_S_A_1,n_S_B_1,n_F_A_1,n_F_B_1
  )
  
  for(i in 1:length(var_type)){
    
    if(var_type[i]=="binary"){
      
      # estimated effects
      data[[(paste("est_delta_F_bi_1_E",i,sep = ""))]] = mean_F$mean[mean_F$trt=="B"& mean_F$variable==(paste("V",i,sep = ""))] - mean_F$mean[mean_F$trt=="A"& mean_F$variable==(paste("V",i,sep = ""))]
      data[[(paste("est_delta_S_bi_1_E",i,sep = ""))]] = mean_S_C$mean[mean_S_C$trt=="B" & mean_S_C$gp=="S" & mean_S_C$variable==(paste("V",i,sep = ""))] - mean_S_C$mean[mean_S_C$trt=="A"& mean_S_C$gp=="S" & mean_S_C$variable==(paste("V",i,sep = ""))]
      data[[(paste("est_delta_C_bi_1_E",i,sep = ""))]] = mean_S_C$mean[mean_S_C$trt=="B" & mean_S_C$gp=="C" & mean_S_C$variable==(paste("V",i,sep = ""))] - mean_S_C$mean[mean_S_C$trt=="A"& mean_S_C$gp=="C" & mean_S_C$variable==(paste("V",i,sep = ""))]
      
      # binary responders
      data[[(paste("n_B_S_1_E",i,sep = ""))]] = mean_S_C$n.responder[mean_S_C$gp=="S" & mean_S_C$trt=="B" & mean_S_C$variable==(paste("V",i,sep = ""))]
      data[[(paste("n_A_S_1_E",i,sep = ""))]] = mean_S_C$n.responder[mean_S_C$gp=="S" & mean_S_C$trt=="A" & mean_S_C$variable==(paste("V",i,sep = ""))]
      data[[(paste("n_B_F_1_E",i,sep = ""))]] = mean_F$n.responder[mean_F$trt=="B" & mean_F$variable==(paste("V",i,sep = ""))]
      data[[(paste("n_A_F_1_E",i,sep = ""))]] = mean_F$n.responder[mean_F$trt=="A" & mean_F$variable==(paste("V",i,sep = ""))]
      data[[(paste("n_F_resp_1_E",i,sep = ""))]] =  data[[(paste("n_B_F_1_E",i,sep = ""))]] + data[[(paste("n_A_F_1_E",i,sep = ""))]]
      
      # pvals and Z
      data[[(paste("pval_F_bi_1_E",i,sep = ""))]] = prop.test(c(data[[(paste("n_B_F_1_E",i,sep = ""))]],data[[(paste("n_A_F_1_E",i,sep = ""))]]),c(n_F_B_1,n_F_A_1),alternative="greater", correct = F)$p.value
      data[[(paste("pval_S_bi_1_E",i,sep = ""))]] = prop.test(c(data[[(paste("n_B_S_1_E",i,sep = ""))]],data[[(paste("n_A_S_1_E",i,sep = ""))]]),c(n_S_B_1,n_S_A_1),alternative="greater", correct = F)$p.value
      data[[(paste("z_F_bi_1_E",i,sep = ""))]] =  as.double(prop.test(c(data[[(paste("n_B_F_1_E",i,sep = ""))]],data[[(paste("n_A_F_1_E",i,sep = ""))]]),c(n_F_B_1,n_F_A_1),alternative="greater", correct = F)$statistic %>% sqrt)
      data[[(paste("z_S_bi_1_E",i,sep = ""))]] =  as.double(prop.test(c(data[[(paste("n_B_S_1_E",i,sep = ""))]],data[[(paste("n_A_S_1_E",i,sep = ""))]]),c(n_S_B_1,n_S_A_1),alternative="greater", correct = F)$statistic %>% sqrt)
      
      
    }else if(var_type[i]=="continuous"){
      
      # estimated effects
      data[[(paste("est_delta_F_cts_1_E",i,sep = ""))]] = mean_F$mean[mean_F$trt=="B"& mean_F$variable==(paste("V",i,sep = ""))] - mean_F$mean[mean_F$trt=="A"& mean_F$variable==(paste("V",i,sep = ""))]
      data[[(paste("est_delta_S_cts_1_E",i,sep = ""))]] = mean_S_C$mean[mean_S_C$trt=="B" & mean_S_C$gp=="S" & mean_S_C$variable==(paste("V",i,sep = ""))] - mean_S_C$mean[mean_S_C$trt=="A"& mean_S_C$gp=="S" & mean_S_C$variable==(paste("V",i,sep = ""))]
      data[[(paste("est_delta_C_cts_1_E",i,sep = ""))]] = mean_S_C$mean[mean_S_C$trt=="B" & mean_S_C$gp=="C" & mean_S_C$variable==(paste("V",i,sep = ""))] - mean_S_C$mean[mean_S_C$trt=="A"& mean_S_C$gp=="C" & mean_S_C$variable==(paste("V",i,sep = ""))]
      
      # mean
      data[[(paste("mean_B_S_1_E",i,sep = ""))]] = mean_S_C$mean[mean_S_C$trt=="B" & mean_S_C$gp=="S" & mean_S_C$variable==(paste("V",i,sep = ""))]
      data[[(paste("mean_A_S_1_E",i,sep = ""))]] = mean_S_C$mean[mean_S_C$trt=="A" & mean_S_C$gp=="S" & mean_S_C$variable==(paste("V",i,sep = ""))]
      data[[(paste("mean_B_F_1_E",i,sep = ""))]] = mean_F$mean[mean_F$trt=="B"& mean_F$variable==(paste("V",i,sep = ""))]
      data[[(paste("mean_A_F_1_E",i,sep = ""))]] = mean_F$mean[mean_F$trt=="A"& mean_F$variable==(paste("V",i,sep = ""))]
      data[[(paste("mean_F_cts_1_E",i,sep = ""))]] = (n_F_B_1*data[[(paste("mean_B_F_1_E",i,sep = ""))]]+n_F_A_1*data[[(paste("mean_A_F_1_E",i,sep = ""))]])/(n_F_A_1+n_F_B_1)
      
      # pvals and Z
      data[[(paste("pval_F_cts_1_E",i,sep = ""))]] = t.test(df[[(paste("V",i,sep = ""))]][df$trt=="B"],df[[(paste("V",i,sep = ""))]][df$trt=="A"],alternative="greater")$p.value
      data[[(paste("pval_S_cts_1_E",i,sep = ""))]] = t.test(df[[(paste("V",i,sep = ""))]][df$gp=="S" & df$trt=="B"],df[[(paste("V",i,sep = ""))]][df$gp=="S" & df$trt=="A"],alternative="greater")$p.value
      data[[(paste("z_F_cts_1_E",i,sep = ""))]] = as.double(t.test(df[[(paste("V",i,sep = ""))]][df$trt=="B"],df[[(paste("V",i,sep = ""))]][df$trt=="A"],alternative="greater")$statistic)
      data[[(paste("z_S_cts_1_E",i,sep = ""))]] = as.double(t.test(df[[(paste("V",i,sep = ""))]][df$gp=="S" & df$trt=="B"],df[[(paste("V",i,sep = ""))]][df$gp=="S" & df$trt=="A"],alternative="greater")$statistic)
      
    }
    
    
  } 
  
  
  
  return(data)
}


stage2_data_generation <- function( type=c("S only","F only","S and F"),rho, n, mu_control, delta, cov.matrix,
                                    thres, sig_level,data=NULL, w_1, w_2, var_type
){
  if(!(type%in%c("S only","F only","S and F"))){return(print("invalid type"))}
  
  if(type=="S only"){
    # Sample sizes stage 2
    n_S_A_2 <- n[["stage2_A"]] # all patients recruited are patients with subgroup characteristics
    n_S_B_2 <- n[["stage2_B"]] 
    
    df=sim_S_CPL(delta, mu_control, cov.matrix, n_S_B_2, n_S_A_2,var_type)
    mean_S_C=df%>%
      pivot_longer(contains("V"),names_to = "variable",values_to = "values")%>%
      group_by(variable,gp,trt)%>%
      summarise(
        mean=mean(values),
        n.responder=sum(values),
        n=n()
      )
    # estimated delta and p-value from stage 2
    data <- data.frame(
      # Stage 2 sample size 
      n_S_A_2,n_S_B_2, n_F_A_2=NA,n_F_B_2=NA
    )
    
    for(i in 1:length(var_type)){
      
      if(var_type[i]=="binary"){
        
        # estimated effects
        data[[(paste("est_delta_F_bi_2_E",i,sep = ""))]] = NA
        data[[(paste("est_delta_S_bi_2_E",i,sep = ""))]] = mean_S_C$mean[mean_S_C$trt=="B" & mean_S_C$gp=="S" & mean_S_C$variable==(paste("V",i,sep = ""))] - mean_S_C$mean[mean_S_C$trt=="A"& mean_S_C$gp=="S" & mean_S_C$variable==(paste("V",i,sep = ""))]
        data[[(paste("est_delta_C_bi_2_E",i,sep = ""))]] = NA
        
        # binary responders
        data[[(paste("n_B_S_2_E",i,sep = ""))]] = mean_S_C$n.responder[mean_S_C$gp=="S" & mean_S_C$trt=="B" & mean_S_C$variable==(paste("V",i,sep = ""))]
        data[[(paste("n_A_S_2_E",i,sep = ""))]] = mean_S_C$n.responder[mean_S_C$gp=="S" & mean_S_C$trt=="A" & mean_S_C$variable==(paste("V",i,sep = ""))]
        data[[(paste("n_B_F_2_E",i,sep = ""))]] = NA
        data[[(paste("n_A_F_2_E",i,sep = ""))]] = NA
        data[[(paste("n_F_resp_2_E",i,sep = ""))]] =  NA
        
        # pvals and Z
        data[[(paste("pval_F_bi_2_E",i,sep = ""))]] = NA
        data[[(paste("pval_S_bi_2_E",i,sep = ""))]] = prop.test(c(data[[(paste("n_B_S_2_E",i,sep = ""))]],data[[(paste("n_A_S_2_E",i,sep = ""))]]),c(n_S_B_2,n_S_A_2),alternative="greater", correct = F)$p.value
        data[[(paste("z_F_bi_2_E",i,sep = ""))]] =  NA
        data[[(paste("z_S_bi_2_E",i,sep = ""))]] =  as.double(prop.test(c(data[[(paste("n_B_S_2_E",i,sep = ""))]],data[[(paste("n_A_S_2_E",i,sep = ""))]]),c(n_S_B_2,n_S_A_2),alternative="greater", correct = F)$statistic %>% sqrt)
        
        
      }else if(var_type[i]=="continuous"){
        
        # estimated effects
        data[[(paste("est_delta_F_cts_2_E",i,sep = ""))]] = NA
        data[[(paste("est_delta_S_cts_2_E",i,sep = ""))]] = mean_S_C$mean[mean_S_C$trt=="B" & mean_S_C$gp=="S" & mean_S_C$variable==(paste("V",i,sep = ""))] - mean_S_C$mean[mean_S_C$trt=="A"& mean_S_C$gp=="S" & mean_S_C$variable==(paste("V",i,sep = ""))]
        data[[(paste("est_delta_C_cts_2_E",i,sep = ""))]] = NA
        
        # mean
        data[[(paste("mean_B_S_2_E",i,sep = ""))]] = mean_S_C$mean[mean_S_C$trt=="B" & mean_S_C$gp=="S" & mean_S_C$variable==(paste("V",i,sep = ""))]
        data[[(paste("mean_A_S_2_E",i,sep = ""))]] = mean_S_C$mean[mean_S_C$trt=="A" & mean_S_C$gp=="S" & mean_S_C$variable==(paste("V",i,sep = ""))]
        data[[(paste("mean_B_F_2_E",i,sep = ""))]] = NA
        data[[(paste("mean_A_F_2_E",i,sep = ""))]] = NA
        data[[(paste("mean_F_cts_2_E",i,sep = ""))]] = NA
        
        # pvals and Z
        data[[(paste("pval_F_cts_2_E",i,sep = ""))]] = NA
        data[[(paste("pval_S_cts_2_E",i,sep = ""))]] = t.test(df[[(paste("V",i,sep = ""))]][df$gp=="S" & df$trt=="B"],df[[(paste("V",i,sep = ""))]][df$gp=="S" & df$trt=="A"],alternative="greater")$p.value
        data[[(paste("z_F_cts_2_E",i,sep = ""))]] = NA
        data[[(paste("z_S_cts_2_E",i,sep = ""))]] = as.double(t.test(df[[(paste("V",i,sep = ""))]][df$gp=="S" & df$trt=="B"],df[[(paste("V",i,sep = ""))]][df$gp=="S" & df$trt=="A"],alternative="greater")$statistic)
        
      }
      
      
    } 
    
    return(data)
    
  }
  
  if(type=="F only"|type=="S and F"){
    
    # Sample sizes stage 2
    n_F_A_2 <- n[["stage2_A"]]
    n_F_B_2 <- n[["stage2_B"]]
    n_S_A_2 <- round(rho*n_F_A_2)
    n_S_B_2 <- round(rho*n_F_B_2) # assuming same prevalence as for stage 1 (take estimate from stage 1 know prior to interim analysis)
    n_C_A_2 <- n_F_A_2-n_S_A_2
    n_C_B_2 <- n_F_B_2-n_S_B_2
    
    df=sim_F_CPL(delta, mu_control, cov.matrix, n_S_B_2, n_S_A_2, n_C_B_2, n_C_A_2,var_type)
    
    mean_S_C=df%>%
      pivot_longer(contains("V"),names_to = "variable",values_to = "values")%>%
      group_by(variable,gp,trt)%>%
      summarise(
        mean=mean(values),
        n.responder=sum(values),
        n=n()
      )
    
    mean_F=df%>%
      pivot_longer(contains("V"),names_to = "variable",values_to = "values")%>%
      group_by(variable,trt)%>%
      summarise(
        mean=mean(values),
        n.responder=sum(values),
        n=n()
      )
    
    
    # responder
    n_B_S_2=mean_S_C$n.responder[mean_S_C$gp=="S" & mean_S_C$trt=="B"]
    n_A_S_2=mean_S_C$n.responder[mean_S_C$gp=="S" & mean_S_C$trt=="A"]
    n_B_F_2=mean_F$n.responder[mean_F$trt=="B"]
    n_A_F_2=mean_F$n.responder[mean_F$trt=="A"]
    n_F_resp_2=mean_F$n.responder[mean_F$trt=="B"] + mean_F$n.responder[mean_F$trt=="A"]
    
    
    if(type=="S and F"){
      
      data <- data.frame(n_S_A_2,n_S_B_2, n_F_A_2,n_F_B_2)
      
      for(i in 1:length(var_type)){
        
        if(var_type[i]=="binary"){
          
          # estimated effects
          data[[(paste("est_delta_F_bi_2_E",i,sep = ""))]] = mean_F$mean[mean_F$trt=="B"& mean_F$variable==(paste("V",i,sep = ""))] - mean_F$mean[mean_F$trt=="A"& mean_F$variable==(paste("V",i,sep = ""))]
          data[[(paste("est_delta_S_bi_2_E",i,sep = ""))]] = mean_S_C$mean[mean_S_C$trt=="B" & mean_S_C$gp=="S" & mean_S_C$variable==(paste("V",i,sep = ""))] - mean_S_C$mean[mean_S_C$trt=="A"& mean_S_C$gp=="S" & mean_S_C$variable==(paste("V",i,sep = ""))]
          data[[(paste("est_delta_C_bi_2_E",i,sep = ""))]] = mean_S_C$mean[mean_S_C$trt=="B" & mean_S_C$gp=="C" & mean_S_C$variable==(paste("V",i,sep = ""))] - mean_S_C$mean[mean_S_C$trt=="A"& mean_S_C$gp=="C" & mean_S_C$variable==(paste("V",i,sep = ""))]
          
          # binary responders
          data[[(paste("n_B_S_2_E",i,sep = ""))]] = mean_S_C$n.responder[mean_S_C$gp=="S" & mean_S_C$trt=="B" & mean_S_C$variable==(paste("V",i,sep = ""))]
          data[[(paste("n_A_S_2_E",i,sep = ""))]] = mean_S_C$n.responder[mean_S_C$gp=="S" & mean_S_C$trt=="A" & mean_S_C$variable==(paste("V",i,sep = ""))]
          data[[(paste("n_B_F_2_E",i,sep = ""))]] = mean_F$n.responder[mean_F$trt=="B" & mean_F$variable==(paste("V",i,sep = ""))]
          data[[(paste("n_A_F_2_E",i,sep = ""))]] = mean_F$n.responder[mean_F$trt=="A" & mean_F$variable==(paste("V",i,sep = ""))]
          data[[(paste("n_F_resp_2_E",i,sep = ""))]] =  data[[(paste("n_B_F_2_E",i,sep = ""))]] + data[[(paste("n_A_F_2_E",i,sep = ""))]]
          
          # pvals and Z
          data[[(paste("pval_F_bi_2_E",i,sep = ""))]] = prop.test(c(data[[(paste("n_B_F_2_E",i,sep = ""))]],data[[(paste("n_A_F_2_E",i,sep = ""))]]),c(n_F_B_2,n_F_A_2),alternative="greater", correct = F)$p.value
          data[[(paste("pval_S_bi_2_E",i,sep = ""))]] = prop.test(c(data[[(paste("n_B_S_2_E",i,sep = ""))]],data[[(paste("n_A_S_2_E",i,sep = ""))]]),c(n_S_B_2,n_S_A_2),alternative="greater", correct = F)$p.value
          data[[(paste("z_F_bi_2_E",i,sep = ""))]] =  as.double(prop.test(c(data[[(paste("n_B_F_2_E",i,sep = ""))]],data[[(paste("n_A_F_2_E",i,sep = ""))]]),c(n_F_B_2,n_F_A_2),alternative="greater", correct = F)$statistic %>% sqrt)
          data[[(paste("z_S_bi_2_E",i,sep = ""))]] =  as.double(prop.test(c(data[[(paste("n_B_S_2_E",i,sep = ""))]],data[[(paste("n_A_S_2_E",i,sep = ""))]]),c(n_S_B_2,n_S_A_2),alternative="greater", correct = F)$statistic %>% sqrt)
          
          
        }else if(var_type[i]=="continuous"){
          
          # estimated effects
          data[[(paste("est_delta_F_cts_2_E",i,sep = ""))]] = mean_F$mean[mean_F$trt=="B"& mean_F$variable==(paste("V",i,sep = ""))] - mean_F$mean[mean_F$trt=="A"& mean_F$variable==(paste("V",i,sep = ""))]
          data[[(paste("est_delta_S_cts_2_E",i,sep = ""))]] = mean_S_C$mean[mean_S_C$trt=="B" & mean_S_C$gp=="S" & mean_S_C$variable==(paste("V",i,sep = ""))] - mean_S_C$mean[mean_S_C$trt=="A"& mean_S_C$gp=="S" & mean_S_C$variable==(paste("V",i,sep = ""))]
          data[[(paste("est_delta_C_cts_2_E",i,sep = ""))]] = mean_S_C$mean[mean_S_C$trt=="B" & mean_S_C$gp=="C" & mean_S_C$variable==(paste("V",i,sep = ""))] - mean_S_C$mean[mean_S_C$trt=="A"& mean_S_C$gp=="C" & mean_S_C$variable==(paste("V",i,sep = ""))]
          
          # mean
          data[[(paste("mean_B_S_2_E",i,sep = ""))]] = mean_S_C$mean[mean_S_C$trt=="B" & mean_S_C$gp=="S" & mean_S_C$variable==(paste("V",i,sep = ""))]
          data[[(paste("mean_A_S_2_E",i,sep = ""))]] = mean_S_C$mean[mean_S_C$trt=="A" & mean_S_C$gp=="S" & mean_S_C$variable==(paste("V",i,sep = ""))]
          data[[(paste("mean_B_F_2_E",i,sep = ""))]] = mean_F$mean[mean_F$trt=="B"& mean_F$variable==(paste("V",i,sep = ""))]
          data[[(paste("mean_A_F_2_E",i,sep = ""))]] = mean_F$mean[mean_F$trt=="A"& mean_F$variable==(paste("V",i,sep = ""))]
          data[[(paste("mean_F_cts_2_E",i,sep = ""))]] = (n_F_B_2*data[[(paste("mean_B_F_2_E",i,sep = ""))]]+n_F_A_2*data[[(paste("mean_A_F_2_E",i,sep = ""))]])/(n_F_A_2+n_F_B_2)
          
          # pvals and Z
          data[[(paste("pval_F_cts_2_E",i,sep = ""))]] = t.test(df[[(paste("V",i,sep = ""))]][df$trt=="B"],df[[(paste("V",i,sep = ""))]][df$trt=="A"],alternative="greater")$p.value
          data[[(paste("pval_S_cts_2_E",i,sep = ""))]] = t.test(df[[(paste("V",i,sep = ""))]][df$gp=="S" & df$trt=="B"],df[[(paste("V",i,sep = ""))]][df$gp=="S" & df$trt=="A"],alternative="greater")$p.value
          data[[(paste("z_F_cts_2_E",i,sep = ""))]] = as.double(t.test(df[[(paste("V",i,sep = ""))]][df$trt=="B"],df[[(paste("V",i,sep = ""))]][df$trt=="A"],alternative="greater")$statistic)
          data[[(paste("z_S_cts_2_E",i,sep = ""))]] = as.double(t.test(df[[(paste("V",i,sep = ""))]][df$gp=="S" & df$trt=="B"],df[[(paste("V",i,sep = ""))]][df$gp=="S" & df$trt=="A"],alternative="greater")$statistic)
          
        }
        
        
      } 
      
      
      
    } 
    else if(type=="F only"){
      
      
      data <- data.frame(n_S_A_2=NA,n_S_B_2=NA,n_F_A_2,n_F_B_2)
      for(i in 1:length(var_type)){
        
        if(var_type[i]=="binary"){
          
          # estimated effects
          data[[(paste("est_delta_F_bi_2_E",i,sep = ""))]] = mean_F$mean[mean_F$trt=="B"& mean_F$variable==(paste("V",i,sep = ""))] - mean_F$mean[mean_F$trt=="A"& mean_F$variable==(paste("V",i,sep = ""))]
          data[[(paste("est_delta_S_bi_2_E",i,sep = ""))]] = NA
          data[[(paste("est_delta_C_bi_2_E",i,sep = ""))]] = NA
          
          # binary responders
          data[[(paste("n_B_S_2_E",i,sep = ""))]] = NA
          data[[(paste("n_A_S_2_E",i,sep = ""))]] = NA
          data[[(paste("n_B_F_2_E",i,sep = ""))]] = mean_F$n.responder[mean_F$trt=="B" & mean_F$variable==(paste("V",i,sep = ""))]
          data[[(paste("n_A_F_2_E",i,sep = ""))]] = mean_F$n.responder[mean_F$trt=="A" & mean_F$variable==(paste("V",i,sep = ""))]
          data[[(paste("n_F_resp_2_E",i,sep = ""))]] =  data[[(paste("n_B_F_2_E",i,sep = ""))]] + data[[(paste("n_A_F_2_E",i,sep = ""))]]
          
          # pvals and Z
          data[[(paste("pval_F_bi_2_E",i,sep = ""))]] = prop.test(c(data[[(paste("n_B_F_2_E",i,sep = ""))]],data[[(paste("n_A_F_2_E",i,sep = ""))]]),c(n_F_B_2,n_F_A_2),alternative="greater", correct = F)$p.value
          data[[(paste("pval_S_bi_2_E",i,sep = ""))]] = NA
          data[[(paste("z_F_bi_2_E",i,sep = ""))]] =  as.double(prop.test(c(data[[(paste("n_B_F_2_E",i,sep = ""))]],data[[(paste("n_A_F_2_E",i,sep = ""))]]),c(n_F_B_2,n_F_A_2),alternative="greater", correct = F)$statistic %>% sqrt)
          data[[(paste("z_S_bi_2_E",i,sep = ""))]] = NA
          
        }else if(var_type[i]=="continuous"){
          
          # estimated effects
          data[[(paste("est_delta_F_cts_2_E",i,sep = ""))]] = mean_F$mean[mean_F$trt=="B"& mean_F$variable==(paste("V",i,sep = ""))] - mean_F$mean[mean_F$trt=="A"& mean_F$variable==(paste("V",i,sep = ""))]
          data[[(paste("est_delta_S_cts_2_E",i,sep = ""))]] =NA
          data[[(paste("est_delta_C_cts_2_E",i,sep = ""))]] = NA
          
          # mean
          data[[(paste("mean_B_S_2_E",i,sep = ""))]] = NA
          data[[(paste("mean_A_S_2_E",i,sep = ""))]] = NA
          data[[(paste("mean_B_F_2_E",i,sep = ""))]] = mean_F$mean[mean_F$trt=="B"& mean_F$variable==(paste("V",i,sep = ""))]
          data[[(paste("mean_A_F_2_E",i,sep = ""))]] = mean_F$mean[mean_F$trt=="A"& mean_F$variable==(paste("V",i,sep = ""))]
          data[[(paste("mean_F_cts_2_E",i,sep = ""))]] = (n_F_B_2*data[[(paste("mean_B_F_2_E",i,sep = ""))]]+n_F_A_2*data[[(paste("mean_A_F_2_E",i,sep = ""))]])/(n_F_A_2+n_F_B_2)
          
          # unadjusted stage 2 pvals and Z : 
          data[[(paste("pval_F_cts_2_E",i,sep = ""))]] = t.test(df[[(paste("V",i,sep = ""))]][df$trt=="B"],df[[(paste("V",i,sep = ""))]][df$trt=="A"],alternative="greater")$p.value
          data[[(paste("pval_S_cts_2_E",i,sep = ""))]] = NA
          data[[(paste("z_F_cts_2_E",i,sep = ""))]] = as.double(t.test(df[[(paste("V",i,sep = ""))]][df$trt=="B"],df[[(paste("V",i,sep = ""))]][df$trt=="A"],alternative="greater")$statistic)
          data[[(paste("z_S_cts_2_E",i,sep = ""))]] = NA
          
        }
        
        
      } 
      
    }
    
    return(data)
  }
}

#-----------------------------------------------------------------------------------------
#--- sim_S_CPL: Generate correlated outcomes for S only 
#-----------------------------------------------------------------------------------------

sim_S_CPL<- function(delta, mu_control, cov.matrix, n_S_B, n_S_A,var_type){
  
  # S, for TRT group B
  data_S_B <- pnorm( mvtnorm::rmvnorm(n_S_B, sigma = cov.matrix) )
  # S, for control group A 
  data_S_A <- pnorm( mvtnorm::rmvnorm(n_S_A, sigma = cov.matrix) )
  
  # data transformation
  for(i in 1:length(var_type)){
    if (var_type[i]=="binary") {
      data_S_B[,i] = qbinom(data_S_B[,i], size = 1, prob = mu_control[[2*i-1]]+delta[[2*i-1]])
      data_S_A[,i] = qbinom(data_S_A[,i], size = 1, prob = mu_control[[2*i-1]])
    }else if (var_type[i]=="continuous") {
      data_S_B[,i] = qnorm(data_S_B[,i], mean = mu_control[[2*i-1]]+delta[[2*i-1]], sd = sqrt(cov.matrix[i,i]))
      data_S_A[,i] = qnorm(data_S_A[,i], mean = mu_control[[2*i-1]], sd = sqrt(cov.matrix[i,i]))
    }
  }
  # trt group in S
  data_S_B = as.data.frame(data_S_B)%>%
    mutate(
      gp = "S",
      trt = "B"
    )
  # control group in S
  data_S_A = as.data.frame(data_S_A)%>%
    mutate(
      gp = "S",
      trt = "A"
    )
  
  data=rbind(data_S_B,data_S_A)
  return(data)
  
}

#-----------------------------------------------------------------------------------------
#--- sim_F_CPL: Generate correlated outcomes for F
#-----------------------------------------------------------------------------------------

sim_F_CPL<- function(delta, mu_control, cov.matrix, n_S_B, n_S_A,n_C_B, n_C_A,var_type){
  
  # S, for TRT group B
  data_S_B <- pnorm( mvtnorm::rmvnorm(n_S_B, sigma = cov.matrix) )
  # S, for control group A 
  data_S_A <- pnorm( mvtnorm::rmvnorm(n_S_A, sigma = cov.matrix) )
  # C, for TRT group B
  data_C_B <- pnorm( mvtnorm::rmvnorm(n_C_B, sigma = cov.matrix) )
  # C, for TRT group A
  data_C_A <- pnorm( mvtnorm::rmvnorm(n_C_A, sigma = cov.matrix) )
  
  # data transformation
  for(i in 1:length(var_type)){
    if (var_type[i]=="binary") {
      data_S_B[,i] = qbinom(data_S_B[,i], size = 1, prob = mu_control[[2*i-1]]+delta[[2*i-1]])
      data_S_A[,i] = qbinom(data_S_A[,i], size = 1, prob = mu_control[[2*i-1]])
      data_C_B[,i] = qbinom(data_C_B[,i], size = 1, prob = mu_control[[2*i]]+delta[[2*i]])
      data_C_A[,i] = qbinom(data_C_A[,i], size = 1, prob = mu_control[[2*i]])
    }else if (var_type[i]=="continuous") {
      data_S_B[,i] = qnorm(data_S_B[,i], mean = mu_control[[2*i-1]]+delta[[2*i-1]], sd = sqrt(cov.matrix[i,i]))
      data_S_A[,i] = qnorm(data_S_A[,i], mean = mu_control[[2*i-1]], sd = sqrt(cov.matrix[i,i]))
      data_C_B[,i] = qnorm(data_C_B[,i], mean = mu_control[[2*i]]+delta[[2*i]], sd = sqrt(cov.matrix[i,i]))
      data_C_A[,i] = qnorm(data_C_A[,i], mean = mu_control[[2*i]], sd = sqrt(cov.matrix[i,i]))
    }
  }
  # trt group in S
  data_S_B = as.data.frame(data_S_B)%>%
    mutate(
      gp = "S",
      trt = "B"
    )
  # control group in S
  data_S_A = as.data.frame(data_S_A)%>%
    mutate(
      gp = "S",
      trt = "A"
    )
  # trt group in C
  data_C_B = as.data.frame(data_C_B)%>%
    mutate(
      gp = "C",
      trt = "B"
    )
  # control group in C
  data_C_A = as.data.frame(data_C_A)%>%
    mutate(
      gp = "C",
      trt = "A"
    )
  
  data=rbind(data_S_B,data_S_A,data_C_B,data_C_A)
  return(data)
}

#-----------------------------------------------------------------------------------------
#--- Interim decisions on group selection at least one endpoint
#-----------------------------------------------------------------------------------------

Interim_Decision <- function(data,thres,var_type){
  
  #-----------------------------------------------------------------------------------------
  #--- Step 1: group selection based on threshold
  #-----------------------------------------------------------------------------------------
  # continue to stage 2 with F if any of the endpoints in F pass the thresholds
  
  effect.F = as.vector(data%>% dplyr::select(dplyr::starts_with("est_delta_F")))
  effect.S = as.vector(data%>% dplyr::select(dplyr::starts_with("est_delta_S")))
  # check if true for each hypothesis
  check.F = if_else(effect.F>=unlist(thres),TRUE,FALSE)
  check.S = if_else(effect.S>=unlist(thres),TRUE,FALSE)
  # check if any true
  pass.F=any(check.F==TRUE)
  pass.S=any(check.S==TRUE)
  
  if (pass.F & !pass.S) {
    interim_decision="2- Continue to stage 2: F only"
  }
  
  # continue to stage 2 with S if any of the endpoints in S pass the thresholds
  # Sime
  
  if (!pass.F & pass.S) {
    interim_decision="1- Continue to stage 2: S only"
  }
  
  
  # continue to stage 2 with F if any of the endpoints in S and any of the endpoints in C pass the thresholds
  if (pass.F & pass.S) {
    interim_decision="3- Continue to stage 2: S and F"
  }
  
  #-----------------------------------------------------------------------------------------
  #--- Step 2: futility based on threshold
  #-----------------------------------------------------------------------------------------
  
  # continue to stage 2 with F if any of the endpoints in S and any of the endpoints in C pass the thresholds
  if (!pass.F & !pass.S) {
    interim_decision="4 - Stop after stage 1: futility"
  }
  return(interim_decision)
}

# test function 

adjusted_p <- function(method,alpha,H,graph.matrix,graph.weights){
  
  # H is the raw pvals
  if (method=="CTP"){
    
    # Hochberg adjusted p-value (Methods for Closed Testing with Simes Inequality)
    result <- p.adjust(H,"hochberg")
    
  }else if(method=="Graphical"){
    # mygraph
    
    graph <- new("graphMCP", m=graph.matrix, weights=graph.weights)
    result=gMCP(graph, H, test="Simes", alpha)@adjPValues
    names(result)=names(H)
    
  }else if(method=="Hierarchical"){
    
    g1 <- H[grepl("H_S",names(H))] 
    g2 <- H[grepl("H_F",names(H))] 
    pvals_g1_hommel <- p.adjust(g1,"hochberg")
    
    if(all(pvals_g1_hommel <= alpha)){
      pvals_g2_hommel <- p.adjust(g2,"hochberg")
    }else{
      pvals_g2_hommel<- rep(1,length(g2))
      
    }
    result=NULL
    for(i in 1:length(g1)){
      result[(paste("H_F_",i,sep = ""))]=pvals_g2_hommel[i]
      result[(paste("H_S_",i,sep = ""))]=pvals_g1_hommel[i]
    }
    
  }
  
  return(result)
}

