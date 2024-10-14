library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggeasy)
library(gridExtra)
source("AED_analysis_function.R")

#----------- data compiling function------------
sim_data = function(nsim, rho, sig_level, n, var_type, mu_control, delta_set, cov.matrix, thres, method=c("CTP","Hierarchical","Graphical"), graph.matrix, graph.weights, stage1_data=NULL, stage2_data=NULL, mydata=NULL, w_1=NULL, w_2=NULL){
  data <- list()
  id <- 1
  for (j in 1:length(delta_set[[1]])) {
    delta <- lapply(delta_set, function(x) x[[j]])
    data[[id]] <- AED_analysis(nsim, rho, sig_level, n, var_type, mu_control, delta, cov.matrix, thres, method, graph.matrix, graph.weights, stage1_data=NULL, stage2_data=NULL, mydata=NULL, w_1=NULL, w_2=NULL)
    id <- id+1
  }#j
  return(data)
}

#-------------- table generating function for type I error-------
out.tab.t1e=function(res_t1e,delta_set,nsim,n){
  
  out_tab_t1e <- list()
  
  for (j in 1:length(res_t1e)) {
    out_tab_t1e[[j]] <- list()
    delta <- lapply(delta_set, function(x) x[[j]])
    cat("Assumptions are ","\n")
    for (q in 1:length(delta)){
      
      cat(paste(names(delta)[q],"=",delta[q],sep = ""),"\n")
    }
    
    tmp <- res_t1e[[j]] 
    
    ### Decision frequency
    out_tab_t1e[[j]]$dec_freq_tab <- dec_freq_tab <- prop.table(table(tmp$interim_decision))
    kable(dec_freq_tab, digits=5, col.names=c("Outcome","Rel. frequency"), caption = 'Decision Freq.') %>% print
    cat('\n\n')
    
    ### sample size tab
    out_tab_t1e[[j]]$sample_size_tab <- sample_size_tab <- as.table(c(Stage1_Average=n[["stage1_A"]]+n[["stage1_B"]],
                                                                      Stage2_Average=sum(res_t1e[[j]]$n_F_B_2+res_t1e[[j]]$n_F_A_2,na.rm = TRUE)/nsim,TotalN_Average=n[["stage1_A"]]+n[["stage1_B"]]+sum(res_t1e[[j]]$n_F_B_2+res_t1e[[j]]$n_F_A_2,na.rm = TRUE)/nsim))
    kable(sample_size_tab, digits=5, col.names=c("Stage","Avereage Sample size"), caption = 'Average Sample size.') %>% print
    cat('\n\n')
    
    ### Type I error
    out_tab_t1e[[j]]$typeIerror_tab <- 
      typeIerror_tab <- c(
        "Type I error (F)"=mean(tmp$signif_F_any),
        "Type I error (S)"=mean(tmp$signif_S_any),
        "Type I error (all)"=mean(tmp$signif_Sany_or_Fany)
        
      )
    kable(typeIerror_tab,digits=5,col.names=c("Type I error"), caption='Type I error') %>% print
    cat('\n\n')
    
  } # end of for (j in 1:length(res))
  
  return(out_tab_t1e)
}

#-------------- table generating function for power-------
out.tab.pwr=function(res_pwr,delta_set,nsim,n){
  out_tab_pwr <- list()
  for (j in 1:length(res_pwr)) {
    out_tab_pwr[[j]] <- list()
    
    delta <- lapply(delta_set, function(x) x[[j]])
    cat("Assumptions are ","\n")
    for (q in 1:length(delta)){
      
      cat(paste(names(delta)[q],"=",delta[q],sep = ""),"\n")
    }
    
    tmp <- res_pwr[[j]] 
    
    ### Decision frequency
    out_tab_pwr[[j]]$dec_freq_tab <- dec_freq_tab <- prop.table(table(tmp$interim_decision))
    kable(dec_freq_tab, digits=5, col.names=c("Outcome","Rel. frequency"), caption = 'Decision Freq.') %>% print
    cat('\n\n')
    
    ### efficacy frequency
    eff_freq_tab <- table(tmp$interim_test)/nsim
    out_tab_pwr[[j]]$eff_freq_tab<- eff_freq_tab <- as.table(c(eff_freq_tab,Efficacy=sum(eff_freq_tab)))
    
    kable(eff_freq_tab, digits=5, col.names=c("Outcome","Rel. frequency"), caption = 'Efficacy Freq.') %>% print
    cat('\n\n')
    
    ### sample size tab
    out_tab_pwr[[j]]$sample_size_tab <- sample_size_tab <- as.table(c(Stage1_Average=n[["stage1_A"]]+n[["stage1_B"]],Stage2_Average=sum(res_pwr[[j]]$n_F_B_2+res_pwr[[j]]$n_F_A_2,na.rm = TRUE)/nsim,TotalN_Average=n[["stage1_A"]]+n[["stage1_B"]]+sum(res_pwr[[j]]$n_F_B_2+res_pwr[[j]]$n_F_A_2,na.rm = TRUE)/nsim))
    kable(sample_size_tab, digits=5, col.names=c("Stage","Avereage Sample size"), caption = 'Average Sample size.') %>% print
    cat('\n\n')
    
    
    ### Power
    out_tab_pwr[[j]]$power_tab <- 
      power_tab <- c("Power (S)"=mean(tmp$signif_S_all),
                     "Power (F)"=mean(tmp$signif_F_all),
                     "Power S or F"=mean(tmp$signif_Sall_or_Fall)
                     
      )
    kable(power_tab,digits=5,col.names=c("Power"), caption='Power') %>% print
    cat('\n\n')
  }
  return(out_tab_pwr)
}

#------------FWER FIGURE--------------------
FWER.visual=function(out_tab_t1e_CTP,out_tab_t1e_graphical,out_tab_t1e_hierarchical){
  FWER=data.frame(FWER=c(rep(NA,3*length(out_tab_t1e_CTP))),Method=c(rep(NA,3*length(out_tab_t1e_CTP))),Scenario=c(rep(NA,3*length(out_tab_t1e_CTP))))
  for (i in 1:length(out_tab_t1e_CTP)) {
    FWER$FWER[3*(i-1)+1]= out_tab_t1e_CTP[[i]]$typeIerror_tab["Type I error (all)"]
    FWER$FWER[3*(i-1)+2]=out_tab_t1e_graphical[[i]]$typeIerror_tab["Type I error (all)"]
    FWER$FWER[3*(i-1)+3]=out_tab_t1e_hierarchical[[i]]$typeIerror_tab["Type I error (all)"]
    
    FWER$Method[3*(i-1)+1]="CTP"
    FWER$Method[3*(i-1)+2]= "Graphical"
    FWER$Method[3*(i-1)+3]= "Hierarchical"
    
    FWER$Scenario[3*(i-1)+1]=i
    FWER$Scenario[3*(i-1)+2]= i
    FWER$Scenario[3*(i-1)+3]= i
  }
  FWER$FWER=round(FWER$FWER,digits = 4)
  ggplot(FWER, aes(x=as.character(Scenario),y=FWER, fill=Method)) +
    geom_bar(stat='identity', position='dodge')+scale_fill_brewer()+ggtitle("FWER")+easy_center_title()+xlab("")+geom_text(aes(label=FWER), position=position_dodge(width=0.9),size=3, vjust=-0.25)
  
}

#------------Power S FIGURE--------------------
power.S.visual=function(out_tab_pwr_CTP,out_tab_pwr_graphical,out_tab_pwr_hierarchical){
  S_pwr=data.frame(S.power=c(rep(NA,3*length(out_tab_pwr_CTP))),Method=c(rep(NA,3*length(out_tab_pwr_CTP))),Scenario=c(rep(NA,3*length(out_tab_pwr_CTP))),max_power=c(rep(NA,3*length(out_tab_pwr_CTP))))
  for (i in 1:length(out_tab_pwr_CTP)) {
    S_pwr$S.power[3*(i-1)+1]= out_tab_pwr_CTP[[i]]$power_tab["Power (S)"]
    S_pwr$S.power[3*(i-1)+2]=out_tab_pwr_graphical[[i]]$power_tab["Power (S)"]
    S_pwr$S.power[3*(i-1)+3]=out_tab_pwr_hierarchical[[i]]$power_tab["Power (S)"]
    max=max(S_pwr$S.power[3*(i-1)+1],S_pwr$S.power[3*(i-1)+2],S_pwr$S.power[3*(i-1)+3])
    
    S_pwr$max_power[3*(i-1)+1]= max
    S_pwr$max_power[3*(i-1)+2]=max
    S_pwr$max_power[3*(i-1)+3]=max
    
    S_pwr$Method[3*(i-1)+1]="CTP"
    S_pwr$Method[3*(i-1)+2]= "Graphical"
    S_pwr$Method[3*(i-1)+3]= "Hierarchical"
    
    S_pwr$Scenario[3*(i-1)+1]=as.integer(i)
    S_pwr$Scenario[3*(i-1)+2]= as.integer(i)
    S_pwr$Scenario[3*(i-1)+3]= as.integer(i)
  }
  S_pwr$Power_Difference=round(S_pwr$max_power-S_pwr$S.power,digits = 4)
  
  a=ggplot(S_pwr, aes(x=as.character(Scenario), y=S.power, fill=Method)) +
    geom_bar(stat='identity', position='dodge')+scale_fill_brewer()+ggtitle("Power (S subpopulation)")+easy_center_title()+geom_text(aes(label=S.power), position=position_dodge(width=0.9), size=2.5, vjust=-0.25)+xlab("Scenario")+ylab("Power")
  
  b=ggplot(S_pwr, aes(x=as.character(Scenario), y=Power_Difference , fill=Method)) +
    geom_bar(stat='identity', position='dodge')+scale_fill_brewer()+ggtitle("S: Difference from Maximum Power")+easy_center_title()+geom_text(aes(label=Power_Difference),position=position_dodge(width=0.9),size=2,  vjust=-0.25)+xlab("Scenario")+ylab("Power Difference")
  
  return(list(a,b))
}

#------------Power F FIGURE--------------------
power.F.visual=function(out_tab_pwr_CTP,out_tab_pwr_graphical,out_tab_pwr_hierarchical){
  F_pwr=data.frame(F_power=c(rep(NA,3*length(out_tab_pwr_CTP))),Method=c(rep(NA,3*length(out_tab_pwr_CTP))),Scenario=c(rep(NA,3*length(out_tab_pwr_CTP))),max_power=c(rep(NA,3*length(out_tab_pwr_CTP))))
  for (i in 1:length(out_tab_pwr_CTP)) {
    F_pwr$F_power[3*(i-1)+1]= out_tab_pwr_CTP[[i]]$power_tab["Power (F)"]
    F_pwr$F_power[3*(i-1)+2]=out_tab_pwr_graphical[[i]]$power_tab["Power (F)"]
    F_pwr$F_power[3*(i-1)+3]=out_tab_pwr_hierarchical[[i]]$power_tab["Power (F)"]
    max=max(F_pwr$F_power[3*(i-1)+1],F_pwr$F_power[3*(i-1)+2],F_pwr$F_power[3*(i-1)+3])
    
    F_pwr$max_power[3*(i-1)+1]= max
    F_pwr$max_power[3*(i-1)+2]=max
    F_pwr$max_power[3*(i-1)+3]=max
    
    F_pwr$Method[3*(i-1)+1]="CTP"
    F_pwr$Method[3*(i-1)+2]= "Graphical"
    F_pwr$Method[3*(i-1)+3]= "Hierarchical"
    
    F_pwr$Scenario[3*(i-1)+1]=as.integer(i)
    F_pwr$Scenario[3*(i-1)+2]= as.integer(i)
    F_pwr$Scenario[3*(i-1)+3]= as.integer(i)
  }
  F_pwr$Power_Difference=round(F_pwr$max_power-F_pwr$F_power,digits = 4)
  
  a=ggplot(F_pwr, aes(x=as.character(Scenario), y=F_power, fill=Method)) +
    geom_bar(stat='identity', position='dodge')+scale_fill_brewer()+ggtitle("Power (F population)")+easy_center_title()+geom_text(aes(label=F_power), position=position_dodge(width=0.9), size=2.5, vjust=-0.25)+xlab("Scenario")+ylab("Power")
  
  b=ggplot(F_pwr, aes(x=as.character(Scenario), y=Power_Difference , fill=Method)) +
    geom_bar(stat='identity', position='dodge')+scale_fill_brewer()+ggtitle("F: Difference from Maximum Power")+easy_center_title()+geom_text(aes(label=Power_Difference),position=position_dodge(width=0.9),size=2,  vjust=-0.25)+xlab("Scenario")+ylab("Power Difference")
  
  return(list(a,b))
}

#------------Power S or F FIGURE--------------------
power.any.visual=function(out_tab_pwr_CTP,out_tab_pwr_graphical,out_tab_pwr_hierarchical){
  any_pwr=data.frame(any_power=c(rep(NA,3*length(out_tab_pwr_CTP))),Method=c(rep(NA,3*length(out_tab_pwr_CTP))),Scenario=c(rep(NA,3*length(out_tab_pwr_CTP))),max_power=c(rep(NA,3*length(out_tab_pwr_CTP))))
  for (i in 1:length(out_tab_pwr_CTP)) {
    any_pwr$any_power[3*(i-1)+1]= out_tab_pwr_CTP[[i]]$power_tab["Power S or F"]
    any_pwr$any_power[3*(i-1)+2]=out_tab_pwr_graphical[[i]]$power_tab["Power S or F"]
    any_pwr$any_power[3*(i-1)+3]=out_tab_pwr_hierarchical[[i]]$power_tab["Power S or F"]
    max=max(any_pwr$any_power[3*(i-1)+1],any_pwr$any_power[3*(i-1)+2],any_pwr$any_power[3*(i-1)+3])
    
    any_pwr$max_power[3*(i-1)+1]= max
    any_pwr$max_power[3*(i-1)+2]=max
    any_pwr$max_power[3*(i-1)+3]=max
    
    any_pwr$Method[3*(i-1)+1]="CTP"
    any_pwr$Method[3*(i-1)+2]= "Graphical"
    any_pwr$Method[3*(i-1)+3]= "Hierarchical"
    
    any_pwr$Scenario[3*(i-1)+1]=as.integer(i)
    any_pwr$Scenario[3*(i-1)+2]= as.integer(i)
    any_pwr$Scenario[3*(i-1)+3]= as.integer(i)
  }
  any_pwr$Power_Difference=round(any_pwr$max_power-any_pwr$any_power,digits = 4)
  
  a=ggplot(any_pwr, aes(x=as.character(Scenario), y=any_power, fill=Method)) +
    geom_bar(stat='identity', position='dodge')+scale_fill_brewer()+ggtitle("Power S or F")+easy_center_title()+geom_text(aes(label=any_power), position=position_dodge(width=0.9), size=2.5, vjust=-0.25)+xlab("Scenario")+ylab("Power")
  
  b=ggplot(any_pwr, aes(x=as.character(Scenario), y=Power_Difference , fill=Method)) +
    geom_bar(stat='identity', position='dodge')+scale_fill_brewer()+ggtitle("S or F: Difference from Maximum Power")+easy_center_title()+geom_text(aes(label=Power_Difference),position=position_dodge(width=0.9),size=2,  vjust=-0.25)+xlab("Scenario")+ylab("Power Difference")
  
  return(list(a,b))
}

#------------efficacy FIGURE--------------------
eff_visual=function(out_tab_pwr_CTP,out_tab_pwr_graphical,out_tab_pwr_hierarchical){
  Efficacy=data.frame(efficacy=c(rep(NA,3*length(out_tab_pwr_CTP))),Method=c(rep(NA,3*length(out_tab_pwr_CTP))),Scenario=c(rep(NA,3*length(out_tab_pwr_CTP))),max_efficacy=c(rep(NA,3*length(out_tab_pwr_CTP))))
  for (i in 1:length(out_tab_pwr_CTP)) {
    Efficacy$efficacy[3*(i-1)+1]= out_tab_pwr_CTP[[i]]$eff_freq_tab["Efficacy"]
    Efficacy$efficacy[3*(i-1)+2]=out_tab_pwr_graphical[[i]]$eff_freq_tab["Efficacy"]
    Efficacy$efficacy[3*(i-1)+3]=out_tab_pwr_hierarchical[[i]]$eff_freq_tab["Efficacy"]
    max=max(Efficacy$efficacy[3*(i-1)+1],Efficacy$efficacy[3*(i-1)+2],Efficacy$efficacy[3*(i-1)+3])
    
    Efficacy$max_efficacy[3*(i-1)+1]= max
    Efficacy$max_efficacy[3*(i-1)+2]=max
    Efficacy$max_efficacy[3*(i-1)+3]=max
    
    Efficacy$Method[3*(i-1)+1]="CTP"
    Efficacy$Method[3*(i-1)+2]= "Graphical"
    Efficacy$Method[3*(i-1)+3]= "Hierarchical"
    
    Efficacy$Scenario[3*(i-1)+1]=as.integer(i)
    Efficacy$Scenario[3*(i-1)+2]= as.integer(i)
    Efficacy$Scenario[3*(i-1)+3]= as.integer(i)
  }
  Efficacy$Efficacy_Difference=round(Efficacy$max_efficacy-Efficacy$efficacy,digits = 4)
  
  a=ggplot(Efficacy, aes(x=as.character(Scenario), y=efficacy, fill=Method)) +
    geom_bar(stat='identity', position='dodge')+scale_fill_brewer()+ggtitle("Early Efficacy Stopping Probability")+easy_center_title()+geom_text(aes(label=efficacy), position=position_dodge(width=0.9), size=2, vjust=-0.25)+xlab("Scenario")+ylab("Probability")
  
  b=ggplot(Efficacy, aes(x=as.character(Scenario), y=Efficacy_Difference , fill=Method)) +
    geom_bar(stat='identity', position='dodge')+scale_fill_brewer()+ggtitle("Early Efficacy: Difference from Maximum")+easy_center_title()+geom_text(aes(label=Efficacy_Difference),position=position_dodge(width=0.9),size=2,  vjust=-0.25)+xlab("Scenario")+ylab("Efficacy Difference")
  return(list(a,b))
}

#------------Average Sample Size FIGURE--------------------
samplesize_visual=function(out_tab_CTP,out_tab_graphical,out_tab_hierarchical){
  SS=data.frame(Total_N=c(rep(NA,3*length(out_tab_CTP))),Method=c(rep(NA,3*length(out_tab_CTP))),Scenario=c(rep(NA,3*length(out_tab_CTP))),min_ss=c(rep(NA,3*length(out_tab_CTP))))
  for (i in 1:length(out_tab_CTP)) {
    SS$Total_N[3*(i-1)+1]= out_tab_CTP[[i]]$sample_size_tab["TotalN_Average"]
    SS$Total_N[3*(i-1)+2]=out_tab_graphical[[i]]$sample_size_tab["TotalN_Average"]
    SS$Total_N[3*(i-1)+3]=out_tab_hierarchical[[i]]$sample_size_tab["TotalN_Average"]
    min=min(SS$Total_N[3*(i-1)+1],SS$Total_N[3*(i-1)+2],SS$Total_N[3*(i-1)+3])
    
    SS$min_ss[3*(i-1)+1]= min
    SS$min_ss[3*(i-1)+2]=min
    SS$min_ss[3*(i-1)+3]=min
    
    SS$Method[3*(i-1)+1]="CTP"
    SS$Method[3*(i-1)+2]= "Graphical"
    SS$Method[3*(i-1)+3]= "Hierarchical"
    
    SS$Scenario[3*(i-1)+1]=as.integer(i)
    SS$Scenario[3*(i-1)+2]= as.integer(i)
    SS$Scenario[3*(i-1)+3]= as.integer(i)
  }
  SS$N_Difference=round(SS$Total_N-SS$min_ss,digits = 4)
  
  a=ggplot(SS, aes(x=as.character(Scenario), y=Total_N, fill=Method)) +
    geom_bar(stat='identity', position='dodge')+scale_fill_brewer()+ggtitle("Average Sample Size Across Stages")+easy_center_title()+geom_text(aes(label=Total_N),position=position_dodge(width=0.9),size=2, vjust=-0.25)+xlab("Scenario")+ylab("Average Sample Size")
  
  
  b=ggplot(SS, aes(x=as.character(Scenario), y=N_Difference , fill=Method)) +
    geom_bar(stat='identity', position='dodge')+scale_fill_brewer()+ggtitle("Average Sample size: Difference from Minimum")+easy_center_title()+geom_text(aes(label=N_Difference),position=position_dodge(width=0.9),size=2, vjust=-0.25)+xlab("Scenario")+ylab("Average Sample Size Difference")
  
  return(list(a,b))
}
