source("helper.R")
library(parallel)
library(foreach)
library(doParallel)
#library(doMPI)
library(doSNOW)

numCores <- detectCores(logical = FALSE)
cl <- parallel::makeCluster(2)
doParallel::registerDoParallel(cl)
comb <- function(x, ...) {
  mapply(rbind, x, ..., SIMPLIFY = F)
}
r_libs=c('dplyr','rpact','magrittr','tidyverse','MASS','ltm','gMCP') 
sapply(r_libs, function (x) library(x, character.only=TRUE))

#-----------------------------------------------------------------------------------------
#--- Analysis function using difference method ("CTP","Hierarchical","Graphical")
#-----------------------------------------------------------------------------------------

AED_analysis <- function(nsim, rho, sig_level, n, var_type, mu_control, delta, 
                         cov.matrix, thres, method=c("CTP","Hierarchical","Graphical"), 
                         graph.matrix, graph.weights, stage1_data=NULL, 
                         stage2_data=NULL, mydata=NULL, w_1=NULL, w_2=NULL
                         
){
  
  ##------- Simulation of adaptive enrichment design with interim monitoring and multiple endpoints
  # tryCatch({
  if (!(method%in%c("CTP","Hierarchical","Graphical"))){
    return(print("please enter a valid method"))
  }
  
  if((is.null(graph.matrix)|all(graph.matrix==0)) & (is.null(graph.weights)|all(graph.weights==0)) & method=="Graphical"){
    print("Missing graph.matrix or graph.weights")
    print("Fixed Sequence Test will be performed")
    graph.matrix=cbind(0,diag(1,length(var_type)*2 , length(var_type)*2))[,c(1:(length(var_type)*2))]
    graph.weights=c(1,rep(0,length(var_type)*2-1))
    
  }
  # calculate the mean(control group) and the effect size in F population
  for (i in 1:length(var_type)){
    mu_control[[paste("mu", i, "_F",sep = "")]]=rho*mu_control[[2*i-1]]+(1-rho)*mu_control[[2*i]]
    delta[[paste("delta", i, "_F",sep = "")]]=rho*delta[[2*i-1]]+(1-rho)*delta[[2*i]]
  }
  
  
  ##-------prepare parameters------------------------------
  ##-------------------------------------------------------
  # Sample size stage 1 per population
  # for A control group
  n$S_A_1 <- round(rho*n[["stage1_A"]])
  n$C_A_1 <- n[["stage1_A"]]-n[["S_A_1"]]
  # for B treatment group
  n$S_B_1 <- round(rho*n[["stage1_B"]])
  n$C_B_1 <- n[["stage1_B"]]-n[["S_B_1"]]
  
  # Weights p-value combination for stage 1 and stage 2
  n$all <- n[["stage1_A"]] + n[["stage1_B"]] + n[["stage2_A"]] + n[["stage2_B"]]
  info_1<- (n[["stage1_A"]] + n[["stage1_B"]]) / n[["all"]]
  
  w_1 <- sqrt(info_1)
  w_2 <- sqrt(1-info_1)
  
  # Calculate nominal one-significance level for stage 1 & 2 combination test using the group sequential design methodology
  sig_BD <- getDesignInverseNormal(informationRates = c(w_1^2,1), typeOfDesign ="asOF",
                                   alpha=sig_level)$stageLevels
  sig_1<- sig_BD[1]
  # sig_comb<- sig_BD[2]
  sig_comb <- sig_level - sig_1
  
  ##-------prepare stage 1 data placeholder----------------
  ##-------------------------------------------------------
  # A is Control group, B is trt group. n_S_A_1 is sample size of control group in S at stage 1 
  stage1_data <- data.frame(n_S_A_1=rep(NA,nsim),n_S_B_1=rep(NA,nsim),n_F_A_1=rep(NA,nsim),n_F_B_1=rep(NA,nsim))
  for(i in 1:length(var_type)){
    
    if(var_type[i]=="binary"){
      
      # estimated effects
      # i is index for endpoints
      stage1_data[[(paste("est_delta_F_bi_1_E",i,sep = ""))]] = rep(NA,nsim)
      stage1_data[[(paste("est_delta_S_bi_1_E",i,sep = ""))]] = rep(NA,nsim)
      stage1_data[[(paste("est_delta_C_bi_1_E",i,sep = ""))]] = rep(NA,nsim)
      
      # binary responders
      stage1_data[[(paste("n_B_S_1_E",i,sep = ""))]] = rep(NA,nsim)
      stage1_data[[(paste("n_A_S_1_E",i,sep = ""))]] = rep(NA,nsim)
      stage1_data[[(paste("n_B_F_1_E",i,sep = ""))]] = rep(NA,nsim)
      stage1_data[[(paste("n_A_F_1_E",i,sep = ""))]] = rep(NA,nsim)
      stage1_data[[(paste("n_F_resp_1_E",i,sep = ""))]] = rep(NA,nsim)
      
      # pvals and Z
      stage1_data[[(paste("pval_F_bi_1_E",i,sep = ""))]] = rep(NA,nsim)
      stage1_data[[(paste("pval_S_bi_1_E",i,sep = ""))]] = rep(NA,nsim)
      stage1_data[[(paste("z_F_bi_1_E",i,sep = ""))]] = rep(NA,nsim)
      stage1_data[[(paste("z_S_bi_1_E",i,sep = ""))]] = rep(NA,nsim)
      
    }else if(var_type[i]=="continuous"){
      
      # estimated effects
      stage1_data[[(paste("est_delta_F_cts_1_E",i,sep = ""))]] = rep(NA,nsim)
      stage1_data[[(paste("est_delta_S_cts_1_E",i,sep = ""))]] = rep(NA,nsim)
      stage1_data[[(paste("est_delta_C_cts_1_E",i,sep = ""))]] = rep(NA,nsim)
      
      # mean
      stage1_data[[(paste("mean_B_S_1_E",i,sep = ""))]] = rep(NA,nsim)
      stage1_data[[(paste("mean_A_S_1_E",i,sep = ""))]] = rep(NA,nsim)
      stage1_data[[(paste("mean_B_F_1_E",i,sep = ""))]] = rep(NA,nsim)
      stage1_data[[(paste("mean_A_F_1_E",i,sep = ""))]] = rep(NA,nsim)
      stage1_data[[(paste("mean_F_cts_1_E",i,sep = ""))]] = rep(NA,nsim)
      
      # pvals and Z
      stage1_data[[(paste("pval_F_cts_1_E",i,sep = ""))]] = rep(NA,nsim)
      stage1_data[[(paste("pval_S_cts_1_E",i,sep = ""))]] = rep(NA,nsim)
      stage1_data[[(paste("z_F_cts_1_E",i,sep = ""))]] = rep(NA,nsim)
      stage1_data[[(paste("z_S_cts_1_E",i,sep = ""))]] = rep(NA,nsim)
      
      
    }
    
    
  } 
  
  ##-------prepare stage 2 data placeholder----------------
  ##-------------------------------------------------------
  stage2_data <- data.frame(n_S_A_2=rep(NA,nsim),n_S_B_2=rep(NA,nsim),n_F_A_2=rep(NA,nsim),n_F_B_2=rep(NA,nsim))
  for(i in 1:length(var_type)){
    
    if(var_type[i]=="binary"){
      
      # estimated effects
      stage2_data[[(paste("est_delta_F_bi_2_E",i,sep = ""))]] = rep(NA,nsim)
      stage2_data[[(paste("est_delta_S_bi_2_E",i,sep = ""))]] = rep(NA,nsim)
      stage2_data[[(paste("est_delta_C_bi_2_E",i,sep = ""))]] = rep(NA,nsim)
      
      # binary responders
      stage2_data[[(paste("n_B_S_2_E",i,sep = ""))]] = rep(NA,nsim)
      stage2_data[[(paste("n_A_S_2_E",i,sep = ""))]] = rep(NA,nsim)
      stage2_data[[(paste("n_B_F_2_E",i,sep = ""))]] = rep(NA,nsim)
      stage2_data[[(paste("n_A_F_2_E",i,sep = ""))]] = rep(NA,nsim)
      stage2_data[[(paste("n_F_resp_2_E",i,sep = ""))]] = rep(NA,nsim)
      
      # pvals and Z
      stage2_data[[(paste("pval_F_bi_2_E",i,sep = ""))]] = rep(NA,nsim)
      stage2_data[[(paste("pval_S_bi_2_E",i,sep = ""))]] = rep(NA,nsim)
      stage2_data[[(paste("z_F_bi_2_E",i,sep = ""))]] = rep(NA,nsim)
      stage2_data[[(paste("z_S_bi_2_E",i,sep = ""))]] = rep(NA,nsim)
      
    }else if(var_type[i]=="continuous"){
      
      # estimated effects
      stage2_data[[(paste("est_delta_F_cts_2_E",i,sep = ""))]] = rep(NA,nsim)
      stage2_data[[(paste("est_delta_S_cts_2_E",i,sep = ""))]] = rep(NA,nsim)
      stage2_data[[(paste("est_delta_C_cts_2_E",i,sep = ""))]] = rep(NA,nsim)
      
      # mean
      stage2_data[[(paste("mean_B_S_2_E",i,sep = ""))]] = rep(NA,nsim)
      stage2_data[[(paste("mean_A_S_2_E",i,sep = ""))]] = rep(NA,nsim)
      stage2_data[[(paste("mean_B_F_2_E",i,sep = ""))]] = rep(NA,nsim)
      stage2_data[[(paste("mean_A_F_2_E",i,sep = ""))]] = rep(NA,nsim)
      stage2_data[[(paste("mean_F_cts_2_E",i,sep = ""))]] = rep(NA,nsim)
      
      # pvals and Z
      stage2_data[[(paste("pval_F_cts_2_E",i,sep = ""))]] = rep(NA,nsim)
      stage2_data[[(paste("pval_S_cts_2_E",i,sep = ""))]] = rep(NA,nsim)
      stage2_data[[(paste("z_F_cts_2_E",i,sep = ""))]] = rep(NA,nsim)
      stage2_data[[(paste("z_S_cts_2_E",i,sep = ""))]] = rep(NA,nsim)
      
      
    }
    
    
  }
  
  ##-------prepare result data placeholder----------------
  ##-------------------------------------------------------
  result <- data.frame(interim_decision=rep(NA,nsim),interim_test=rep(NA,nsim), 
                       signif_S_all=rep(NA,nsim), signif_F_all=rep(NA,nsim), 
                       signif_Sall_or_Fall=rep(NA,nsim),signif_S_any=rep(NA,nsim), 
                       signif_F_any=rep(NA,nsim), signif_Sany_or_Fany=rep(NA,nsim) )
  
  for(i in 1:length(var_type)){
    
    if(var_type[i]=="binary"){
      
      # significant or not
      result[[(paste("signif_S_E",i,sep = ""))]]=rep(NA,nsim)
      result[[(paste("signif_F_E",i,sep = ""))]]=rep(NA,nsim) 
      
      # estimated overall effect in trt group
      result[[(paste("effect_bi_1_E",i,sep = ""))]] = rep(NA,nsim)
      result[[(paste("effect_F_bi_overall_E",i,sep = ""))]] = rep(NA,nsim)
      result[[(paste("effect_S_bi_overall_E",i,sep = ""))]] = rep(NA,nsim)
      
      # combination test unadjusted pvalue
      result[[(paste("p_comb_F_E",i,sep = ""))]] = rep(NA,nsim)
      result[[(paste("p_comb_S_E",i,sep = ""))]] = rep(NA,nsim)
      
      # combination test adjusted pvalues
      result[[(paste("adj_p_comb_F_E",i,sep = ""))]] = rep(NA,nsim)
      result[[(paste("adj_p_comb_S_E",i,sep = ""))]] = rep(NA,nsim)
      
      # stage 2 control group
      result[[(paste("est_A_rate_F_bi_2_E",i,sep = ""))]] = rep(NA,nsim)
      result[[(paste("est_A_rate_S_bi_2_E",i,sep = ""))]] = rep(NA,nsim)
      
      
      # first stage info
      ## effect
      result[[(paste("effect_S_bi_1_E",i,sep = ""))]] = rep(NA,nsim)
      result[[(paste("effect_F_bi_1_E",i,sep = ""))]] = rep(NA,nsim)
      result[[(paste("effect_C_bi_1_E",i,sep = ""))]] = rep(NA,nsim)
      result[[(paste("n_F_resp_1_E",i,sep = ""))]] = rep(NA,nsim)
      ## pvals
      result[[(paste("pval_F_bi_1_E",i,sep = ""))]] = rep(NA,nsim)
      result[[(paste("pval_S_bi_1_E",i,sep = ""))]] = rep(NA,nsim)
      ## adjusted pvals
      result[[(paste("adj_pval_F_bi_1_E",i,sep = ""))]] = rep(NA,nsim)
      result[[(paste("adj_pval_S_bi_1_E",i,sep = ""))]] = rep(NA,nsim)
      # effect in control group
      result[[(paste("est_A_F_bi_1_E",i,sep = ""))]] = rep(NA,nsim)
      result[[(paste("est_A_S_bi_1_E",i,sep = ""))]] = rep(NA,nsim)
      
    }else if(var_type[i]=="continuous"){
      
      # significant or not
      result[[(paste("signif_S_E",i,sep = ""))]]=rep(NA,nsim)
      result[[(paste("signif_F_E",i,sep = ""))]]=rep(NA,nsim) 
      
      # estimated overall effect in trt group
      result[[(paste("effect_cts_1_E",i,sep = ""))]] = rep(NA,nsim)
      result[[(paste("effect_F_cts_overall_E",i,sep = ""))]] = rep(NA,nsim)
      result[[(paste("effect_S_cts_overall_E",i,sep = ""))]] = rep(NA,nsim)
      
      # combination test p value
      result[[(paste("p_comb_F_E",i,sep = ""))]] = rep(NA,nsim)
      result[[(paste("p_comb_S_E",i,sep = ""))]] = rep(NA,nsim)
      
      # combination test adjusted pvalues
      result[[(paste("adj_p_comb_F_E",i,sep = ""))]] = rep(NA,nsim)
      result[[(paste("adj_p_comb_S_E",i,sep = ""))]] = rep(NA,nsim)
      
      # stage 2 control group 
      result[[(paste("est_A_rate_F_cts_2_E",i,sep = ""))]] = rep(NA,nsim)
      result[[(paste("est_A_rate_S_cts_2_E",i,sep = ""))]] = rep(NA,nsim)
      
      # first stage info
      ## effect
      result[[(paste("effect_S_cts_1_E",i,sep = ""))]] = rep(NA,nsim)
      result[[(paste("effect_F_cts_1_E",i,sep = ""))]] = rep(NA,nsim)
      result[[(paste("effect_C_cts_1_E",i,sep = ""))]] = rep(NA,nsim)
      result[[(paste("mean_F_cts_1_E",i,sep = ""))]] = rep(NA,nsim)
      ## pvals
      result[[(paste("pval_F_cts_1_E",i,sep = ""))]] = rep(NA,nsim)
      result[[(paste("pval_S_cts_1_E",i,sep = ""))]] = rep(NA,nsim)
      ## adjusted pvals
      result[[(paste("adj_pval_F_cts_1_E",i,sep = ""))]] = rep(NA,nsim)
      result[[(paste("adj_pval_S_cts_1_E",i,sep = ""))]] = rep(NA,nsim)
      # effect in control group
      result[[(paste("est_A_F_cts_1_E",i,sep = ""))]] = rep(NA,nsim)
      result[[(paste("est_A_S_cts_1_E",i,sep = ""))]] = rep(NA,nsim)
      
    }
    
    
  }
  
  ##----simulation for loop: log stage 1 data to result dataset 
  ##-------------------------------------------------------
  foreach (i= 1:nsim)%do% {
    #source("helper.R",local = TRUE)
    set.seed(i+3000)
    print(i)
    
    mydata<-stage1_data_generation( rho, n, mu_control, delta, cov.matrix, 
                                    thres, sig_level,data=NULL, w_1, w_2, var_type)
    
    # update stage1_data
    stage1_data[i,1:ncol(stage1_data)] <- mydata
    for (j in 1:length(var_type)){
      
      if(var_type[j]=="binary"){
        
        # first stage info
        ## effect
        result[[(paste("effect_S_bi_1_E",j,sep = ""))]][i] = stage1_data[[(paste("est_delta_S_bi_1_E",j,sep = ""))]][i]
        result[[(paste("effect_F_bi_1_E",j,sep = ""))]][i] = stage1_data[[(paste("est_delta_F_bi_1_E",j,sep = ""))]][i]
        result[[(paste("effect_C_bi_1_E",j,sep = ""))]][i] = stage1_data[[(paste("est_delta_C_bi_1_E",j,sep = ""))]][i]
        result[[(paste("n_F_resp_1_E",j,sep = ""))]][i] = stage1_data[[(paste("n_F_resp_1_E",j,sep = ""))]][i]
        ## pvals
        result[[(paste("pval_F_bi_1_E",j,sep = ""))]][i] = stage1_data[[(paste("pval_F_bi_1_E",j,sep = ""))]][i]
        result[[(paste("pval_S_bi_1_E",j,sep = ""))]][i] = stage1_data[[(paste("pval_S_bi_1_E",j,sep = ""))]][i]
        
        # effect in control group
        result[[(paste("est_A_F_bi_1_E",j,sep = ""))]][i] = stage1_data[[(paste("n_A_F_1_E",j,sep = ""))]][i] / stage1_data$n_F_A_1[i]
        result[[(paste("est_A_S_bi_1_E",j,sep = ""))]][i] = stage1_data[[(paste("n_A_S_1_E",j,sep = ""))]][i] / stage1_data$n_F_A_1[i]
        
      }else if(var_type[j]=="continuous"){
        
        # first stage info
        ## effect
        result[[(paste("effect_S_cts_1_E",j,sep = ""))]][i] = stage1_data[[(paste("est_delta_S_cts_1_E",j,sep = ""))]][i]
        result[[(paste("effect_F_cts_1_E",j,sep = ""))]][i] = stage1_data[[(paste("est_delta_F_cts_1_E",j,sep = ""))]][i]
        result[[(paste("effect_C_cts_1_E",j,sep = ""))]][i] = stage1_data[[(paste("est_delta_C_cts_1_E",j,sep = ""))]][i]
        result[[(paste("mean_F_cts_1_E",j,sep = ""))]][i] =  stage1_data[[(paste("mean_F_cts_1_E",j,sep = ""))]][i]
        
        ## pvals
        result[[(paste("pval_F_cts_1_E",j,sep = ""))]][i] = stage1_data[[(paste("pval_F_cts_1_E",j,sep = ""))]][i]
        result[[(paste("pval_S_cts_1_E",j,sep = ""))]][i] = stage1_data[[(paste("pval_S_cts_1_E",j,sep = ""))]][i]
        
        # effect in control group
        result[[(paste("est_A_F_cts_1_E",j,sep = ""))]][i] = stage1_data[[(paste("mean_A_F_1_E",j,sep = ""))]][i]
        result[[(paste("est_A_S_cts_1_E",j,sep = ""))]][i] = stage1_data[[(paste("mean_A_S_1_E",j,sep = ""))]][i]
        
      }
      
    }
    
    
    ###------------population selection------------------
    decision<-Interim_Decision(mydata,thres, var_type)
    result$interim_decision[i] <- decision
    
    ###----------- stop for futility---------------------
    if(decision=="4 - Stop after stage 1: futility" ){
      
      for(j in 1:length(var_type)) {
        # significant or not
        result[[(paste("signif_S_E",j,sep = ""))]][i] = FALSE
        result[[(paste("signif_F_E",j,sep = ""))]][i] = FALSE 
      }
      
      #print("futility")
    }
    
    ###-------- stop for efficacy in S-------------------
    if (decision=="1- Continue to stage 2: S only"){
      H=NULL
      # only test endpoints for S
      for(j in 1:length(var_type)){
        if(var_type[j]=="binary"){
          # F :  set to 0 since we are looking at S only
          H[(paste("H_F_",j,sep = ""))]=1
          # S 
          H[(paste("H_S_",j,sep = ""))]=result[[(paste("pval_S_bi_1_E",j,sep = ""))]][i] 
        }else if (var_type[j]=="continuous"){
          # F
          H[(paste("H_F_",j,sep = ""))]=1
          # S 
          H[(paste("H_S_",j,sep = ""))]=result[[(paste("pval_S_cts_1_E",j,sep = ""))]][i] 
        }
        
      }
      
      
      p_adjusted <- adjusted_p(method,sig_1,H,graph.matrix,graph.weights)
      
      for(j in 1:length(var_type)){
        if(var_type[j]=="binary"){
          # S 
          result[[(paste("adj_pval_S_bi_1_E",j,sep = ""))]][i] = p_adjusted[(paste("H_S_",j,sep = ""))] 
          result[[(paste("signif_F_E",j,sep = ""))]][i] <- F
          result[[(paste("signif_S_E",j,sep = ""))]][i] <- (result[[(paste("adj_pval_S_bi_1_E",j,sep = ""))]][i] <= sig_1)
        }else if (var_type[j]=="continuous"){
          # S 
          result[[(paste("adj_pval_S_cts_1_E",j,sep = ""))]][i]= p_adjusted[(paste("H_S_",j,sep = ""))]
          result[[(paste("signif_F_E",j,sep = ""))]][i] = F 
          result[[(paste("signif_S_E",j,sep = ""))]][i] =  (result[[(paste("adj_pval_S_cts_1_E",j,sep = ""))]][i] <= sig_1)
          
        }
        
      }
      
      #print(p_adjusted)
      
      # efficacy result
      check=NULL
      for(j in 1:length(var_type)){
        check=c(check,result[[(paste("signif_S_E",j,sep = ""))]][i])
      }
      
      if(all(check==TRUE)){
        result$interim_test[i] <- "6 - Stop after stage 1: efficacy in S only"
        
      }  
    } 
    
    
    ###---------stop for efficacy in F------------------  
    if (decision=="2- Continue to stage 2: F only"){
      H=NULL
      # only test endpoints for F
      for(j in 1:length(var_type)){
        if(var_type[j]=="binary"){
          # F 
          H[(paste("H_F_",j,sep = ""))]=result[[(paste("pval_F_bi_1_E",j,sep = ""))]][i] 
          # S 
          H[(paste("H_S_",j,sep = ""))]=1
        }else if (var_type[j]=="continuous"){
          # F
          H[(paste("H_F_",j,sep = ""))]=result[[(paste("pval_F_cts_1_E",j,sep = ""))]][i] 
          # S 
          H[(paste("H_S_",j,sep = ""))]=1
        }
        
      }
      # only test H1 H3; set H2=H4=0
      p_adjusted <- adjusted_p(method,sig_1,H,graph.matrix,graph.weights)
      
      #print(p_adjusted)
      
      for(j in 1:length(var_type)){
        if(var_type[j]=="binary"){
          # S 
          result[[(paste("adj_pval_F_bi_1_E",j,sep = ""))]][i] = p_adjusted[(paste("H_F_",j,sep = ""))] 
          result[[(paste("signif_F_E",j,sep = ""))]][i] <- (result[[(paste("adj_pval_F_bi_1_E",j,sep = ""))]][i] <= sig_1)
          result[[(paste("signif_S_E",j,sep = ""))]][i] <- F 
        }else if (var_type[j]=="continuous"){
          # S 
          result[[(paste("adj_pval_F_cts_1_E",j,sep = ""))]][i]= p_adjusted[(paste("H_F_",j,sep = ""))]
          result[[(paste("signif_F_E",j,sep = ""))]][i] = (result[[(paste("adj_pval_F_cts_1_E",j,sep = ""))]][i] <= sig_1)
          result[[(paste("signif_S_E",j,sep = ""))]][i] = F 
          
        }
        
      }
      
      # efficacy result
      check=NULL
      for(j in 1:length(var_type)){
        check=c(check,result[[(paste("signif_F_E",j,sep = ""))]][i])
      }
      
      if(all(check==TRUE)){
        result$interim_test[i] <- "5 - Stop after stage 1: efficacy in F only"
        
      }  
      
    } 
    
    ### --------stop for efficacy in S and F------------ 
    if (decision=="3- Continue to stage 2: S and F"){
      H=NULL
      # only test endpoints for F
      for(j in 1:length(var_type)){
        if(var_type[j]=="binary"){
          # F 
          H[(paste("H_F_",j,sep = ""))]=result[[(paste("pval_F_bi_1_E",j,sep = ""))]][i] 
          # S 
          H[(paste("H_S_",j,sep = ""))]=result[[(paste("pval_S_bi_1_E",j,sep = ""))]][i] 
        }else if (var_type[j]=="continuous"){
          # F
          H[(paste("H_F_",j,sep = ""))]=result[[(paste("pval_F_cts_1_E",j,sep = ""))]][i] 
          # S 
          H[(paste("H_S_",j,sep = ""))]=result[[(paste("pval_S_cts_1_E",j,sep = ""))]][i]
        }
        
      }
      
      p_adjusted <- adjusted_p(method,sig_1,H,graph.matrix,graph.weights)
      #print(p_adjusted)
      for(j in 1:length(var_type)){
        if(var_type[j]=="binary"){
          # S 
          result[[(paste("adj_pval_S_bi_1_E",j,sep = ""))]][i] = p_adjusted[(paste("H_S_",j,sep = ""))] 
          result[[(paste("adj_pval_F_bi_1_E",j,sep = ""))]][i] = p_adjusted[(paste("H_F_",j,sep = ""))] 
          
          result[[(paste("signif_F_E",j,sep = ""))]][i] <- (result[[(paste("adj_pval_F_bi_1_E",j,sep = ""))]][i] <= sig_1)
          result[[(paste("signif_S_E",j,sep = ""))]][i] <- (result[[(paste("adj_pval_S_bi_1_E",j,sep = ""))]][i] <= sig_1)
        }else if (var_type[j]=="continuous"){
          # S 
          result[[(paste("adj_pval_S_cts_1_E",j,sep = ""))]][i]= p_adjusted[(paste("H_S_",j,sep = ""))]
          result[[(paste("adj_pval_F_cts_1_E",j,sep = ""))]][i] = p_adjusted[(paste("H_F_",j,sep = ""))]
          
          result[[(paste("signif_F_E",j,sep = ""))]][i] = (result[[(paste("adj_pval_F_cts_1_E",j,sep = ""))]][i] <= sig_1) 
          result[[(paste("signif_S_E",j,sep = ""))]][i] =  (result[[(paste("adj_pval_S_cts_1_E",j,sep = ""))]][i] <= sig_1)
          
        }
      }
      
      # NEED TO CHANGE
      check.F=NULL
      check.S=NULL
      
      for(j in 1:length(var_type)){
        check.F=c(check.F,result[[(paste("signif_F_E",j,sep = ""))]][i])
        check.S=c(check.S,result[[(paste("signif_S_E",j,sep = ""))]][i])
      }
      
      if(all(check.S==TRUE)){
        result$interim_test[i] <- "6 - Stop after stage 1: efficacy in S only"
      }    
      
      if(all(check.F==TRUE)){
        result$interim_test[i] <- "5 - Stop after stage 1: efficacy in F only"
      }  
      
      # efficacy result: three situation in S&F
      if(all(check.S==TRUE) & all(check.F==TRUE)){
        
        result$interim_test[i] <- "7 - Stop after stage 1: efficacy in S and F"
        
      }
      
    } 
    
    ###------------- update overall effect for those stop at interim-----------------
    if(result$interim_test[i] %in% c("5 - Stop after stage 1: efficacy in F only")){
      for(j in 1:length(var_type)){
        if(var_type[j]=="binary"){
          # S 
          result[[(paste("effect_bi_1_E",j,sep = ""))]][i] = stage1_data[[(paste("est_delta_F_bi_1_E",j,sep = ""))]][i] 
        }else if (var_type[j]=="continuous"){
          # S 
          result[[(paste("effect_cts_1_E",j,sep = ""))]][i]= stage1_data[[(paste("est_delta_F_cts_1_E",j,sep = ""))]][i]
        }
      }
    }
    
    if(decision=="4 - Stop after stage 1: futility" | 
       result$interim_test[i] %in% c("5 - Stop after stage 1: efficacy in F only",
                                     "6 - Stop after stage 1: efficacy in S only",
                                     "7 - Stop after stage 1: efficacy in S and F")){
      
      for(j in 1:length(var_type)){
        if(var_type[j]=="binary"){
          # S 
          result[[(paste("effect_F_bi_overall_E",j,sep = ""))]][i] = stage1_data[[(paste("est_delta_F_bi_1_E",j,sep = ""))]][i] 
          result[[(paste("effect_S_bi_overall_E",j,sep = ""))]][i] = stage1_data[[(paste("est_delta_S_bi_1_E",j,sep = ""))]][i] 
          
        }else if (var_type[j]=="continuous"){
          # S 
          result[[(paste("effect_F_cts_overall_E",j,sep = ""))]][i]= stage1_data[[(paste("est_delta_F_cts_1_E",j,sep = ""))]][i]
          result[[(paste("effect_S_cts_overall_E",j,sep = ""))]][i]= stage1_data[[(paste("est_delta_S_cts_1_E",j,sep = ""))]][i]
        }
      }
      
    }
    
    ### --------- continue to stage 2------------------------
    # Decision to continue to stage 2 with S only
    if (decision=="1- Continue to stage 2: S only" & is.na(result$interim_test[i])) {
      mydata2<- stage2_data_generation(type="S only",rho, n, mu_control, delta, cov.matrix,
                                       thres, sig_level,data=NULL, w_1, w_2, var_type 
      )
      stage2_data[i,1:ncol(stage2_data)] <- mydata2
      
      #print("S only")
      
      # H is a temp var for using MTP. For h we don't want to/need to test, set to 0
      H= NULL
      for(j in 1:length(var_type)){
        # pvals for F set to NA in the dataset, set to 0 for H
        result[[(paste("p_comb_F_E",j,sep = ""))]][i] = NA
        H[(paste("H_F_",j,sep = ""))]=1
        
        if(var_type[j]=="binary"){
          # Z stats for S
          assign((paste("z_comb_S_E",j,sep = "")), w_1*qnorm(1-stage1_data[[(paste("pval_S_bi_1_E",j,sep = ""))]][i])+w_2*qnorm(1-stage2_data[[(paste("pval_S_bi_2_E",j,sep = ""))]][i]))
          
        }else if (var_type[j]=="continuous"){
          assign((paste("z_comb_S_E",j,sep = "")), w_1*qnorm(1-stage1_data[[(paste("pval_S_cts_1_E",j,sep = ""))]][i])+w_2*qnorm(1-stage2_data[[(paste("pval_S_cts_2_E",j,sep = ""))]][i]))
          
        }
        
        # if significant in STAGE 1, set to NA
        if(result[[(paste("signif_S_E",j,sep = ""))]][i]){
          result[[(paste("p_comb_S_E",j,sep = ""))]][i]=NA
          H[(paste("H_S_",j,sep = ""))]=0
        }else{
          result[[(paste("p_comb_S_E",j,sep = ""))]][i] <- pnorm(-get((paste("z_comb_S_E",j,sep = "")))) 
          H[(paste("H_S_",j,sep = ""))]=result[[(paste("p_comb_S_E",j,sep = ""))]][i]
        }
        
      }
      
      p_adjusted <- adjusted_p(method,sig_comb,H,graph.matrix,graph.weights)
      
      #print(p_adjusted)
      
      # log p values and significance or not
      for(j in 1:length(var_type)){
        
        if(!is.na(result[[(paste("p_comb_S_E",j,sep = ""))]][i])){
          result[[(paste("adj_p_comb_S_E",j,sep = ""))]][i]=p_adjusted[(paste("H_S_",j,sep = ""))]
          # test decision
          if(result[[(paste("adj_p_comb_S_E",j,sep = ""))]][i]<=sig_comb){
            result[[(paste("signif_S_E",j,sep = ""))]][i] <- T
          }else{
            result[[(paste("signif_S_E",j,sep = ""))]][i] <- F
          }
        }
        
      }
      
      # log overall effects
      
      for(j in 1:length(var_type)){
        if(var_type[j]=="binary"){
          
          tot_B_resp=stage1_data[[(paste("n_B_S_1_E",j,sep = ""))]][i] + stage2_data[[(paste("n_B_S_2_E",j,sep = ""))]][i]
          tot_B=stage1_data$n_S_B_1[i] + stage2_data$n_S_B_2[i]
          tot_A_resp=stage1_data[[(paste("n_A_S_1_E",j,sep = ""))]][i] + stage2_data[[(paste("n_A_S_2_E",j,sep = ""))]][i]
          tot_A=stage1_data$n_S_A_1[i] + stage2_data$n_S_A_2[i]
          result[[(paste("effect_S_bi_overall_E",j,sep = ""))]][i] <- (tot_B_resp/tot_B)-(tot_A_resp/tot_A)
          result[[(paste("effect_F_bi_overall_E",j,sep = ""))]][i] <-  stage1_data[[(paste("est_delta_F_bi_1_E",j,sep = ""))]][i]
          
        }else if (var_type[j]=="continuous"){
          sum_B_1=stage1_data[[(paste("mean_B_S_1_E",j,sep = ""))]][i] * stage1_data$n_S_B_1[i]
          sum_B_2=stage2_data[[(paste("mean_B_S_2_E",j,sep = ""))]][i] * stage2_data$n_S_B_2[i]
          num.B=stage1_data$n_S_B_1[i] + stage2_data$n_S_B_2[i]
          sum_A_1=stage1_data[[(paste("mean_A_S_1_E",j,sep = ""))]][i] * stage1_data$n_S_A_1[i]
          sum_A_2=stage2_data[[(paste("mean_A_S_2_E",j,sep = ""))]][i] * stage2_data$n_S_A_2[i]
          num.A=stage1_data$n_S_A_1[i] + stage2_data$n_S_A_2[i]
          result[[(paste("effect_S_cts_overall_E",j,sep = ""))]][i] <- (sum_B_1+sum_B_2)/num.B-(sum_A_1+sum_A_2)/num.A
          result[[(paste("effect_F_cts_overall_E",j,sep = ""))]][i] <-  stage1_data[[(paste("est_delta_F_cts_1_E",j,sep = ""))]][i]
          
          
        }
      }
      
    }
    
    
    # Decision to continue to stage 2 with F only
    # Sime
    if (decision=="2- Continue to stage 2: F only"  & is.na(result$interim_test[i])){
      mydata2 <- stage2_data_generation(type="F only",rho, n, mu_control, delta, cov.matrix,
                                        thres, sig_level,data=NULL, w_1, w_2, var_type)
      
      stage2_data[i,1:ncol(stage2_data)] <- mydata2
      
      #print("F only")
      
      # H is a temp var for using MTP. For h we don't want to/need to test, set to 0
      H= NULL
      for(j in 1:length(var_type)){
        # pvals for S set to NA in the dataset, set to 0 for H
        result[[(paste("p_comb_S_E",j,sep = ""))]][i] = NA
        H[(paste("H_S_",j,sep = ""))]=1
        
        if(var_type[j]=="binary"){
          # Z stats for S
          assign((paste("z_comb_F_E",j,sep = "")), w_1*qnorm(1-stage1_data[[(paste("pval_F_bi_1_E",j,sep = ""))]][i])+w_2*qnorm(1-stage2_data[[(paste("pval_F_bi_2_E",j,sep = ""))]][i]))
          
          
        }else if (var_type[j]=="continuous"){
          
          assign((paste("z_comb_F_E",j,sep = "")), w_1*qnorm(1-stage1_data[[(paste("pval_F_cts_1_E",j,sep = ""))]][i])+w_2*qnorm(1-stage2_data[[(paste("pval_F_cts_2_E",j,sep = ""))]][i]))
          
          
        }
        
        # if significant in STAGE 1, set to NA
        if(result[[(paste("signif_F_E",j,sep = ""))]][i]){
          result[[(paste("p_comb_F_E",j,sep = ""))]][i]=NA
          H[(paste("H_F_",j,sep = ""))]=0
        }else{
          result[[(paste("p_comb_F_E",j,sep = ""))]][i] <- pnorm(-get((paste("z_comb_F_E",j,sep = "")))) 
          H[(paste("H_F_",j,sep = ""))]=result[[(paste("p_comb_F_E",j,sep = ""))]][i]
        }
        
      }
      
      p_adjusted <- adjusted_p(method,sig_comb,H,graph.matrix,graph.weights)
      
      #print(p_adjusted)
      
      # log p values and significance or not
      for(j in 1:length(var_type)){
        
        if(!is.na(result[[(paste("p_comb_F_E",j,sep = ""))]][i])){
          result[[(paste("adj_p_comb_F_E",j,sep = ""))]][i]=p_adjusted[(paste("H_F_",j,sep = ""))]
          # test decision
          if(result[[(paste("adj_p_comb_F_E",j,sep = ""))]][i]<=sig_comb){
            result[[(paste("signif_F_E",j,sep = ""))]][i] <- T
          }else{
            result[[(paste("signif_F_E",j,sep = ""))]][i] <- F
          }
        }
        
      }
      
      # log overall effects
      
      for(j in 1:length(var_type)){
        if(var_type[j]=="binary"){
          
          tot_B_resp=stage1_data[[(paste("n_B_F_1_E",j,sep = ""))]][i] + stage2_data[[(paste("n_B_F_2_E",j,sep = ""))]][i]
          tot_B=stage1_data$n_F_B_1[i] + stage2_data$n_F_B_2[i]
          tot_A_resp=stage1_data[[(paste("n_A_F_1_E",j,sep = ""))]][i] + stage2_data[[(paste("n_A_F_2_E",j,sep = ""))]][i]
          tot_A=stage1_data$n_F_A_1[i] + stage2_data$n_F_A_2[i]
          result[[(paste("effect_F_bi_overall_E",j,sep = ""))]][i] <- (tot_B_resp/tot_B)-(tot_A_resp/tot_A)
          result[[(paste("effect_S_bi_overall_E",j,sep = ""))]][i] <-  stage1_data[[(paste("est_delta_S_bi_1_E",j,sep = ""))]][i]
          
        }else if (var_type[j]=="continuous"){
          sum_B_1=stage1_data[[(paste("mean_B_F_1_E",j,sep = ""))]][i] * stage1_data$n_F_B_1[i]
          sum_B_2=stage2_data[[(paste("mean_B_F_2_E",j,sep = ""))]][i] * stage2_data$n_F_B_2[i]
          num.B=stage1_data$n_F_B_1[i] + stage2_data$n_F_B_2[i]
          sum_A_1=stage1_data[[(paste("mean_A_F_1_E",j,sep = ""))]][i] * stage1_data$n_F_A_1[i]
          sum_A_2=stage2_data[[(paste("mean_A_F_2_E",j,sep = ""))]][i] * stage2_data$n_F_A_2[i]
          num.A=stage1_data$n_F_A_1[i] + stage2_data$n_F_A_2[i]
          result[[(paste("effect_F_cts_overall_E",j,sep = ""))]][i] <- (sum_B_1+sum_B_2)/num.B-(sum_A_1+sum_A_2)/num.A
          result[[(paste("effect_S_cts_overall_E",j,sep = ""))]][i] <-  stage1_data[[(paste("est_delta_S_cts_1_E",j,sep = ""))]][i]
          
          
        }
      }
      
    }
    
    if(decision=="3- Continue to stage 2: S and F"  & is.na(result$interim_test[i])){
      
      # The rest: Decision to continue to stage 2 with F and S
      
      mydata2 <- stage2_data_generation(type="S and F",rho, n, mu_control, delta, cov.matrix,
                                        thres, sig_level,data=NULL, w_1, w_2, var_type
      )
      
      stage2_data[i,1:ncol(stage2_data)] <- mydata2
      
      # calculate combination test p-values
      #print("S and F")
      
      # H is a temp var for using MTP. For h we don't want to/need to test, set to 0
      H= NULL
      for(j in 1:length(var_type)){
        
        if(var_type[j]=="binary"){
          # Z stats for S
          assign((paste("z_comb_F_E",j,sep = "")), w_1*qnorm(1-stage1_data[[(paste("pval_F_bi_1_E",j,sep = ""))]][i])+w_2*qnorm(1-stage2_data[[(paste("pval_F_bi_2_E",j,sep = ""))]][i]))
          assign((paste("z_comb_S_E",j,sep = "")), w_1*qnorm(1-stage1_data[[(paste("pval_S_bi_1_E",j,sep = ""))]][i])+w_2*qnorm(1-stage2_data[[(paste("pval_S_bi_2_E",j,sep = ""))]][i]))
          
        }else if (var_type[j]=="continuous"){
          assign((paste("z_comb_F_E",j,sep = "")), w_1*qnorm(1-stage1_data[[(paste("pval_F_cts_1_E",j,sep = ""))]][i])+w_2*qnorm(1-stage2_data[[(paste("pval_F_cts_2_E",j,sep = ""))]][i]))
          assign((paste("z_comb_S_E",j,sep = "")), w_1*qnorm(1-stage1_data[[(paste("pval_S_cts_1_E",j,sep = ""))]][i])+w_2*qnorm(1-stage2_data[[(paste("pval_S_cts_2_E",j,sep = ""))]][i]))
          
        }
        
        # if significant in STAGE 1, set to NA
        if(result[[(paste("signif_F_E",j,sep = ""))]][i]){
          result[[(paste("p_comb_F_E",j,sep = ""))]][i]=NA
          H[(paste("H_F_",j,sep = ""))]=0
        }else{
          result[[(paste("p_comb_F_E",j,sep = ""))]][i] <- pnorm(-get((paste("z_comb_F_E",j,sep = "")))) 
          H[(paste("H_F_",j,sep = ""))]=result[[(paste("p_comb_F_E",j,sep = ""))]][i]
        }
        
        if(result[[(paste("signif_S_E",j,sep = ""))]][i]){
          result[[(paste("p_comb_S_E",j,sep = ""))]][i]=NA
          H[(paste("H_S_",j,sep = ""))]=0
        }else{
          result[[(paste("p_comb_S_E",j,sep = ""))]][i] <- pnorm(-get((paste("z_comb_S_E",j,sep = "")))) 
          H[(paste("H_S_",j,sep = ""))]=result[[(paste("p_comb_S_E",j,sep = ""))]][i]
        }
        
      }
      
      p_adjusted <- adjusted_p(method,sig_comb,H,graph.matrix,graph.weights)
      
      #print(p_adjusted)
      
      # log p values and significance or not
      for(j in 1:length(var_type)){
        
        if(!is.na(result[[(paste("p_comb_F_E",j,sep = ""))]][i])){
          result[[(paste("adj_p_comb_F_E",j,sep = ""))]][i]=p_adjusted[(paste("H_F_",j,sep = ""))]
          # test decision
          if(result[[(paste("adj_p_comb_F_E",j,sep = ""))]][i]<=sig_comb){
            result[[(paste("signif_F_E",j,sep = ""))]][i] <- T
          }else{
            result[[(paste("signif_F_E",j,sep = ""))]][i] <- F
          }
        }
        
        if(!is.na(result[[(paste("p_comb_S_E",j,sep = ""))]][i])){
          result[[(paste("adj_p_comb_S_E",j,sep = ""))]][i]=p_adjusted[(paste("H_S_",j,sep = ""))]
          # test decision
          if(result[[(paste("adj_p_comb_S_E",j,sep = ""))]][i]<=sig_comb){
            result[[(paste("signif_S_E",j,sep = ""))]][i] <- T
          }else{
            result[[(paste("signif_S_E",j,sep = ""))]][i] <- F
          }
        }
        
      }
      
      # log overall effects
      
      for(j in 1:length(var_type)){
        if(var_type[j]=="binary"){
          
          tot_B_resp=stage1_data[[(paste("n_B_F_1_E",j,sep = ""))]][i] + stage2_data[[(paste("n_B_F_2_E",j,sep = ""))]][i]
          tot_B=stage1_data$n_F_B_1[i] + stage2_data$n_F_B_2[i]
          tot_A_resp=stage1_data[[(paste("n_A_F_1_E",j,sep = ""))]][i] + stage2_data[[(paste("n_A_F_2_E",j,sep = ""))]][i]
          tot_A=stage1_data$n_F_A_1[i] + stage2_data$n_F_A_2[i]
          
          result[[(paste("effect_F_bi_overall_E",j,sep = ""))]][i] <-  (tot_B_resp/tot_B)-(tot_A_resp/tot_A)
          
          tot_B_resp=stage1_data[[(paste("n_B_S_1_E",j,sep = ""))]][i] + stage2_data[[(paste("n_B_S_2_E",j,sep = ""))]][i]
          tot_B=stage1_data$n_S_B_1[i] + stage2_data$n_S_B_2[i]
          tot_A_resp=stage1_data[[(paste("n_A_S_1_E",j,sep = ""))]][i] + stage2_data[[(paste("n_A_S_2_E",j,sep = ""))]][i]
          tot_A=stage1_data$n_S_A_1[i] + stage2_data$n_S_A_2[i]
          result[[(paste("effect_S_bi_overall_E",j,sep = ""))]][i] <- (tot_B_resp/tot_B)-(tot_A_resp/tot_A)
          
          
        }else if (var_type[j]=="continuous"){
          sum_B_1=stage1_data[[(paste("mean_B_F_1_E",j,sep = ""))]][i] * stage1_data$n_F_B_1[i]
          sum_B_2=stage2_data[[(paste("mean_B_F_2_E",j,sep = ""))]][i] * stage2_data$n_F_B_2[i]
          num.B=stage1_data$n_F_B_1[i] + stage2_data$n_F_B_2[i]
          sum_A_1=stage1_data[[(paste("mean_A_F_1_E",j,sep = ""))]][i] * stage1_data$n_F_A_1[i]
          sum_A_2=stage2_data[[(paste("mean_A_F_2_E",j,sep = ""))]][i] * stage2_data$n_F_A_2[i]
          num.A=stage1_data$n_F_A_1[i] + stage2_data$n_F_A_2[i]
          result[[(paste("effect_F_cts_overall_E",j,sep = ""))]][i] <-  (sum_B_1+sum_B_2)/num.B-(sum_A_1+sum_A_2)/num.A
          
          sum_B_1=stage1_data[[(paste("mean_B_S_1_E",j,sep = ""))]][i] * stage1_data$n_S_B_1[i]
          sum_B_2=stage2_data[[(paste("mean_B_S_2_E",j,sep = ""))]][i] * stage2_data$n_S_B_2[i]
          num.B=stage1_data$n_S_B_1[i] + stage2_data$n_S_B_2[i]
          sum_A_1=stage1_data[[(paste("mean_A_S_1_E",j,sep = ""))]][i] * stage1_data$n_S_A_1[i]
          sum_A_2=stage2_data[[(paste("mean_A_S_2_E",j,sep = ""))]][i] * stage2_data$n_S_A_2[i]
          num.A=stage1_data$n_S_A_1[i] + stage2_data$n_S_A_2[i]
          result[[(paste("effect_S_cts_overall_E",j,sep = ""))]][i] <- (sum_B_1+sum_B_2)/num.B-(sum_A_1+sum_A_2)/num.A
          
          
        }
      }
      
    }
    
    # log stage 2 sample size
    
    result$n_F_A_2[i]<-  stage2_data$n_F_A_2[i]
    result$n_S_A_2[i]<-  stage2_data$n_S_A_2[i]
    result$n_F_B_2[i]<-  stage2_data$n_F_B_2[i]
    result$n_S_B_2[i]<-  stage2_data$n_S_B_2[i]
    
    
    # log overall significance in S AND F
    # log control group info
    tmp_S=NULL
    tmp_F=NULL
    for(j in 1:length(var_type)){
      if(var_type[j]=="binary"){
        tmp_S=c(tmp_S,result[[(paste("signif_S_E",j,sep = ""))]][i])
        tmp_F=c(tmp_F,result[[(paste("signif_F_E",j,sep = ""))]][i])
        
        result[[(paste("n_A_F_2_E",j,sep = ""))]][i] <-  stage2_data[[(paste("n_A_F_2_E",j,sep = ""))]][i]
        result[[(paste("n_A_S_2_E",j,sep = ""))]][i] <-  stage2_data[[(paste("n_A_S_2_E",j,sep = ""))]][i]
        
        result[[(paste("est_A_rate_F_bi_2_E",j,sep = ""))]][i] = result[[(paste("n_A_F_2_E",j,sep = ""))]][i]/result$n_F_A_2[i]
        result[[(paste("est_A_rate_S_bi_2_E",j,sep = ""))]][i] = result[[(paste("n_A_S_2_E",j,sep = ""))]][i]/result$n_S_A_2[i]
        
      }else if (var_type[j]=="continuous"){
        tmp_S=c(tmp_S,result[[(paste("signif_S_E",j,sep = ""))]][i])
        tmp_F=c(tmp_F,result[[(paste("signif_F_E",j,sep = ""))]][i])
        
        result[[(paste("mean_A_F_2_E",j,sep = ""))]][i] <-  stage2_data[[(paste("mean_A_F_2_E",j,sep = ""))]][i]
        result[[(paste("mean_A_S_2_E",j,sep = ""))]][i] <-  stage2_data[[(paste("mean_A_S_2_E",j,sep = ""))]][i]
        
        result[[(paste("est_A_rate_F_cts_2_E",j,sep = ""))]][i] = result[[(paste("mean_A_F_2_E",j,sep = ""))]][i]
        result[[(paste("est_A_rate_S_cts_2_E",j,sep = ""))]][i] = result[[(paste("mean_A_S_2_E",j,sep = ""))]][i]
        
        
      }
    }
    
    
    # This is for power 
    # pass F population
    if(all(tmp_F==TRUE)){
      result$signif_F_all[i] <- TRUE
    }else{
      result$signif_F_all[i] <- FALSE
    }
    # pass S subpopulation
    if(all(tmp_S==TRUE)){
      result$signif_S_all[i] <- TRUE
    }else{
      result$signif_S_all[i] <- FALSE
    }
    # pass any one of them
    if(all(tmp_S==TRUE)|all(tmp_F==TRUE)){
      result$signif_Sall_or_Fall[i] <- TRUE
    }else{
      result$signif_Sall_or_Fall[i] <- FALSE
    }
    
    # this is for type I error
    # any false positive in F
    if(any(tmp_F==TRUE)){
      result$signif_F_any[i] <- TRUE
    }else{
      result$signif_F_any[i] <- FALSE
    }
    
    # any false positive in S
    if(any(tmp_S==TRUE)){
      result$signif_S_any[i] <- TRUE
    }else{
      result$signif_S_any[i] <- FALSE
    }
    
    # any false positive in S or F: global null
    if(any(tmp_S==TRUE)|any(tmp_F==TRUE)){
      result$signif_Sany_or_Fany[i] <- TRUE
    }else{
      result$signif_Sany_or_Fany[i] <- FALSE
    }
    
  }# END FOR LOOP
  
  
  result$interim_decision <- factor(result$interim_decision,
                                    levels = c("1- Continue to stage 2: S only","2- Continue to stage 2: F only",
                                               "3- Continue to stage 2: S and F","4 - Stop after stage 1: futility"))
  result$interim_test <- factor(result$interim_test, levels=c("5 - Stop after stage 1: efficacy in F only",
                                                              "6 - Stop after stage 1: efficacy in S only",
                                                              "7 - Stop after stage 1: efficacy in S and F"))
  
  result$sig_comb <- sig_comb
  
  return(result)
  #}, error=function(cond) browser())
}
