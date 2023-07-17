#############################################################################################################
# In this script, you can find the implementation of the algorithms presented in the Manuscript 'Structural restrictions in local causal discovery: identifying direct causes of a target variable'
# The main function is 'IDS_and_score_based_estimation_of_F_parents()'
# It takes as an input $Y$ as a target variable and $X=data.frame(X_1, \dots, X_p)$ as covariates. Moreover, 
# it takes as an input 'constraint_set_F'. We implemented the following choices of 'constraint_set_F':
# 'Additive'
# 'Location-scale'
# 'Gaussian'
# 'Gamma'
# 'Pareto'
# 'Gumbel'
# For the 'Additive' case, we also implemented Random forest as an estimation method. If interested, replace 'Additive' by 'Additive_RF' and run all simulations once more, but typically it gave us worse results than GAM
#  If interested, implement your own class by adding an estimation method for epsilon (step 1 of the algorithm) in the estimate_epsilon() function
#  We also implemented function 'F_identifiable_set_of_parents' that is faster but computes only F_identifiable set of parents, not the score nor F-plausibility of every set
#  See an example of usage below


library(copula)
library(boot)
library("evd")
library("mgcv")
library("gamlss")
library(CondIndTests)
library("twosamples")
library(tidyverse)
library(kSamples)
library(independence)
library(Pareto)
library(dHSIC)
library(rje)
library(stringr)
library(TauStar)
library(randomForest)
library(rfPermute)
library(DescTools)
library(raster)
library(ambient)
library(dplyr)
library(rgl)
library(reticulate)
library(feather)

####################################################Example
#n=500
#
#X1 = rnorm(n); X2 = rnorm(n)
#Y = X1*X2  + rnorm(n)              #Watchout here! if you use e.g. Y = X1^2 + X2^2 + rnorm(n) this case is non-identifiable and also sets {X1}, {X2} are F-plausbile!
#X3 = Y^2 + X1 + rnorm(n)
#
#IDS_and_score_based_estimation_of_F_parents(Y=Y, X=as.data.frame(cbind(X1, X2,X3)), constraint_set_F="Additive", show_intermediate_results=FALSE)

#And here is the typical result you should obtain
#>                          score independence significance distribution
#>Empty set                  -Inf            0            1            1
#>X1                       -5.598        0.001         0.27            1
#>X2                         -4.5        0.003         0.27            1
#>X1 X2                     4.604        0.999         0.01            1
#>X3                       -5.635        0.001         0.28            1
#>X1 X3                    -5.011        0.001         0.15            1
#>X2 X3                    -5.737        0.001         0.31            1
#>X1 X2 X3                 -4.094        0.001         0.06            1
#>---                         ---          ---          ---          ---
#>Score-based estimate:     X1 X2                                        
#>Identifiable causes:      X1 X2                                       
#>F-plausible causes:       X1 X2                                       
#> 
  
##################################################Example 2 
#You have to upload functions from the script 'Random_multivariate_functions_and_benchmarks' to run this example

#n=500
#
#X1 = rnorm(n); X2 = rnorm(n)
#f=random_function_2d()  # this generates a random 2D function, check the plot3d()
#Y = evaluation_of_f_2d(f, X1, X2)  + rnorm(n) #evaluation_of_f_2d(f, X1, X2) is basically just =f(X1, X2), but with a bit more ugly notation
#X3 = Y^2 + X1 + rnorm(n)
#
#plot3d(x=X1, y=X2, z=Y)
#
#IDS_and_score_based_estimation_of_F_parents(Y=Y, X=as.data.frame(cbind(X1, X2,X3)), constraint_set_F="Additive", show_intermediate_results=FALSE)
#F_identifiable_set_of_parents(Y=Y, X=as.data.frame(cbind(X1, X2, X3)), constraint_set_F="Additive")







###################  And here is the code for the functions

estimate_epsilon<-function(Y, X, constraint_set_F="Additive"){ 
  
  
  ### ### ### ### ### ### ### ###
  ############### Some helpful functions ##################
  rename_columns_to_1_to_d <-function(X){ #Just making sure that our variables are called X1, X2, ...
    if (is.null(nrow(X))) {X = as.data.frame(X); names(X)[1] <- "X1"; return(X) }
    return(X %>% rename_with(~ str_c("X", seq_along(.))))
  }
  X = rename_columns_to_1_to_d(X)
  names(Y)<- "Y"
  
  significance_maximum <- function(fit, number_of_parameters=2){
    significance=-Inf
    #m=str_count(toString(fit$terms[[3]]), "[0-9]") #How many numbers are in the formula; Chceme iba ziskat pocet covariates in fit
    vector = summary(fit)[[8]]
    m=length(vector)
    if (number_of_parameters==1) {
      for (i in 1:m) {
        if(vector[i]>significance) significance = vector[i]
      }}
    if (number_of_parameters==2) {
      for (i in 1:(m/number_of_parameters)) { number = min(vector[i], vector[ (m/number_of_parameters)+i])
      if(number>significance) significance = number
      }}
    if (number_of_parameters>2) print('Not implemented, sorry.')
    if (significance<0.001) {significance = 0.001}
    
    return(significance)
  }
  
  
  formula <- function(X, withY=TRUE, smoothing=TRUE, type ='additive'){ #Definition of formula. We want to distinguish between discrete and continous variables for GAM
    q = ncol(X); if (is.null(nrow(X))) {q=1} 
    if(withY==TRUE){ form="Y ~ "}
    if(withY==FALSE){form=" ~ "}
    if (type =='additive') {
      
      for (i in 1:q) {#Discrete variable is if it contains <=9 different values
        if (  length(unique(X[,i]))>9  ) { 
          if (smoothing==TRUE) { form=paste(form, paste0("s(X", i,  ")+")) }
          if (smoothing==FALSE) { form=paste(form, paste0("X", i,  "+")) }
        }
        if (  length(unique(X[,i]))<=9  ) {  form=paste(form, paste0("as.factor(X", i,  ")+"))  }
      }
      
    }
    
    if (type =='general') {
      form=paste(form, paste0("s("))
      for (i in 1:q) {
        form=paste(form, paste0("X", i,  ",")) 
      }
    }
    form=substr(form,1,nchar(form)-1)
    if (type == 'general') { form=paste(form, paste0(")"))}
    
    return(as.formula(form))
  }
  
  ###
  
  
  ### ### ### ### ### ### ### ### Here we start with the estimation of epsilons for each family separately.### ### ### ### ### ### ### ###
  
  if (constraint_set_F=="Additive_RF") {
    
    
    fit= randomForest(y=Y, x = X)
    
    residual = Y-predict(fit, newdata = X)
    residual=pnorm(residual)
    
    
    ntree=50; nrep=40
    fit_for_significance = rfPermute(Y~.,data=data.frame(Y, X),ntree=50, num.rep=nrep)$pval #Computes the p-values for each variable in terms of predictible power
    significance = max(fit_for_significance[1:nrow(fit_for_significance),1,2]) 
    significance=significance*(nrep+1)/100 #This is ugly. It is here only because if 'num.rep' is large it takes too much time. So I took only nrep=40 and decreased the resulting pvalue by a factor of ~2 
    
    return(data.frame(residual = residual, significance = significance ))
    
  }
  
  if (constraint_set_F=="Location-scale" || constraint_set_F=="Gaussian"){
    
    fit=  gam(list(formula(X, TRUE, type='general'),formula(X, FALSE, type='general')), data.frame(Y,X) ,family=gaulss(), link=list("identity","logb"))
    
    give_me_epsilons_from_GAULSS <-function(Y, fit){ #Note that GAM with GAULSS family does not estimate mu and sigma, but mu and 1/sigma. Thats why this is here
      residuals=c()
      for (i in 1:length(Y)) {
        residuals=c(residuals,  (Y[i] -fit$fitted.values[i,1])*(fit$fitted.values[i,2]) ) #fit$fitted.values[i,2]=1/sigma[i], fit$fitted.values[i,1]=mu[i]
      }
      return(residuals)
    }
    
    residual = give_me_epsilons_from_GAULSS(Y, fit)
    residual=pnorm(residual)
    
    significance = significance_maximum(gam(list(formula(X, TRUE),formula(X, FALSE)), data.frame(Y,X) ,family=gaulss(), link=list("identity","logb")), number_of_parameters = 2)
    return(data.frame(residual = residual, significance = significance ))
  }
  
  
  if (constraint_set_F=="Additive"){
    
    fit=  gam(list(formula(X, TRUE, type = 'general'),~1), data.frame(Y,X) ,family=gaulss(), link=list("identity","logb"))
    
    give_me_epsilons_from_GAULSS <-function(Y, fit){ #Note that GAM with GAULSS family does not estimate mu and sigma, but mu and 1/sigma. Thats why this is here
      residuals=c()
      for (i in 1:length(Y)) {
        residuals=c(residuals,  (Y[i] -fit$fitted.values[i,1])*(fit$fitted.values[i,2]) ) #fit$fitted.values[i,2]=1/sigma[i], fit$fitted.values[i,1]=mu[i]
      }
      return(residuals)
    }
    
    residual = give_me_epsilons_from_GAULSS(Y, fit)
    residual=pnorm(residual)
    
    ntree=50; nrep=40
    fit_for_significance = rfPermute(Y~.,data=data.frame(Y, X),ntree=50, num.rep=nrep)$pval #Computes the p-values for each variable in terms of predictible power
    significance = max(fit_for_significance[1:nrow(fit_for_significance),1,2]) 
    significance=significance*(nrep+1)/100 #This is ugly. It is here only because if 'num.rep' is large it takes too much time. So I took only nrep=40 and decreased the resulting pvalue by a factor of ~2 
    
    return(data.frame(residual = residual, significance = significance ))
  }
  
  
  
  if (constraint_set_F=="Pareto") {
    
    Y=-log(log(Y)) #Trick: -loglog of pareto is gumbel which is implemented in GAM package
    fit=  gam(list(formula(X, TRUE, type = 'general'), ~1), data.frame(Y,X) ,family=gumbls())
    residual=c()
    for (i in 1:n) {
      residual=c(residual,  pgumbel(Y[i], (fit$fitted.values[i,1]-0.577), 1)  )
    }
    
    
    significance = significance_maximum(gam(list(formula(X, TRUE), ~1), data.frame(Y,X) ,family=gumbls()), number_of_parameters = 1)
    
    return(data.frame(residual = residual, significance = significance ))
    
  }
  
  
  
  if (constraint_set_F=="Gumbel") {
    
    fit=  gam(list(formula(X, TRUE), formula(X, FALSE), type = 'general'), data.frame(Y,X) ,family=gumbls())
    residual=c()
    for (i in 1:n) {
      residual=c(residual,  pgumbel(Y[i], fit$fitted.values[i,1], exp(fit$fitted.values[i,2]))  )
    }
    
    significance = significance_maximum(gam(list(formula(X, TRUE), formula(X, FALSE)), data.frame(Y,X) ,family=gumbls()), number_of_parameters = 2)
    
    return(data.frame(residual = residual, significance = significance ))
    
  }
  
  if (constraint_set_F=="Gamma") {
    
    fit=  gam(list(formula(X, TRUE, type = 'general'),formula(X, FALSE, type = 'general')), data.frame(Y,X) ,family=gammals())
    give_me_generalized_residuals_from_gammals <-function(Y, fit){
      residuals=c()
      for (i in 1:length(Y)) {
        residuals=c(residuals,  pgamma( Y[i], shape = 1/exp(fitted(fit)[i,2]), scale = fitted(fit)[i,1]*exp(fitted(fit)[i,2])   ) ) 
      }
      return(residuals)
    }
    residual = give_me_generalized_residuals_from_gammals(Y, fit)
    
    significance = significance_maximum(gam(list(formula(X, TRUE),formula(X, FALSE)), data.frame(Y,X) ,family=gammals()), number_of_parameters = 2)
    
    return(data.frame(residual = residual, significance = significance ))
    
  }
  
  if (constraint_set_F=="Gamma with fixed scale") {
    
    fit=  gam(list(formula(X, TRUE, type = 'general'),~1), data.frame(Y,X) ,family=gammals())
    
    give_me_generalized_residuals_from_gammals <-function(Y, fit){
      residuals=c()
      for (i in 1:length(Y)) {
        residuals=c(residuals,  pgamma( Y[i], shape = 1/exp(fitted(fit)[i,2]), scale = fitted(fit)[i,1]*exp(fitted(fit)[i,2])   )) 
      }
      return(residuals)
    }
    
    residual = give_me_generalized_residuals_from_gammals(Y, fit)
    
    significance = significance_maximum( gam(list(formula(X, TRUE),~1), data.frame(Y,X) ,family=gammals()), number_of_parameters = 1)
    
    return(data.frame(residual = residual, significance = significance ))
    
  }
  
  
  print("Family not implemented, sorry.")
  return("Family not implemented, sorry.")
  
}


F_identifiable_set_of_parents<- function(Y, X, constraint_set_F="Additive"){
  ###########################Some random helpful functions ##########################
  
  rename_columns_to_1_to_d <-function(X){ 
    if (is.null(nrow(X))) { return(as.data.frame(X1=X)) }
    return(X %>% rename_with(~ str_c("X", seq_along(.))))
  }
  
  ###
  
  formula <- function(X, withY=TRUE, smoothing=TRUE, type ='additive'){ #Definition of formula. We want to distinguish between discrete and continous variables for GAM
    q = ncol(X); if (is.null(nrow(X))) {q=1} 
    if(withY==TRUE){ form="Y ~ "}
    if(withY==FALSE){form=" ~ "}
    if (type =='additive') {
      
      for (i in 1:q) {#Discrete variable is if it contains <=9 different values
        if (  length(unique(X[,i]))>9  ) { 
          if (smoothing==TRUE) { form=paste(form, paste0("s(X", i,  ")+")) }
          if (smoothing==FALSE) { form=paste(form, paste0("X", i,  "+")) }
        }
        if (  length(unique(X[,i]))<=9  ) {  form=paste(form, paste0("as.factor(X", i,  ")+"))  }
      }
      
    }
    
    if (type =='general') {
      form=paste(form, paste0("s("))
      for (i in 1:q) {
        form=paste(form, paste0("X", i,  ",")) 
      }
    }
    form=substr(form,1,nchar(form)-1)
    if (type == 'general') { form=paste(form, paste0(")"))}
    
    return(as.formula(form))
  }
  
  ###
  
  are_all_significant <- function(fit, number_of_parameters=2){
    all_are_significant=TRUE
    m=str_count(toString(fit$terms[[3]]), "[0-9]") #How many numbers are in the formula; Chceme iba ziskat pocet covariates in fit
    vector = summary(fit)[[8]]
    if (number_of_parameters==2) {
      for (i in 1:m) {
        if(min(vector[i], vector[m+i])>0.01) all_are_significant=FALSE
      }}
    if (number_of_parameters==1) {
      for (i in 1:m) {
        if(vector[i]>0.01) all_are_significant=FALSE
      }}
    
    return(all_are_significant)
  }
  
  significance_maximum <- function(fit, number_of_parameters=2){
    significance=-Inf
    #m=str_count(toString(fit$terms[[3]]), "[0-9]") #How many numbers are in the formula; Chceme iba ziskat pocet covariates in fit
    vector = summary(fit)[[8]]
    m=length(vector)
    if (number_of_parameters==1) {
      for (i in 1:m) {
        if(vector[i]>significance) significance = vector[i]
      }}
    if (number_of_parameters==2) {
      for (i in 1:(m/number_of_parameters)) { number = min(vector[i], vector[ (m/number_of_parameters)+i])
      if(number>significance) significance = number
      }}
    if (number_of_parameters>2) print('Not implemented, sorry.')
    if (significance<0.001) {significance = 0.001}
    
    return(significance)
  }
  
  is_sample_uniformly_distributed<-function(x){
    return(AndersonDarlingTest(x, "punif")$p.value)
  }
  
  
  p_value_of_independence<- function(epsilon, Z, accuracy = 0.0025){
    if (is.null(nrow(Z))) { numberofcolumns = 1} else {numberofcolumns = ncol(Z)}
    x=matrix(as.matrix(Z), ncol=numberofcolumns); #x=matrix(c(X1), ncol=1); 
    y=matrix(epsilon, ncol=1)
    
    #if (numberofcolumns==1) { pvalue=hoeffding.D.test(as.vector(x),as.vector(y), na.rm = TRUE, collisions = FALSE, precision = accuracy/26)$p.value  }
    if(length(epsilon)<=1000){pvalue=dhsic.test(list(x,y))$p.value  }
    else{  pvalue=multIndepTest(data.frame(x,y) , d=rep(1, ncol(x)+1), verbose = FALSE, N=as.integer(1/(2*accuracy)-1) )$fisher.pvalue}
    
    
    return(pvalue + 0.000001) #+0.000001 is there as a cutoff value so we do not return zero value (since we later take log(0))
  }
  
  ################# Here starts the code #############
  
  X = rename_columns_to_1_to_d(X)
  m=ncol(X); if (is.null(nrow(X))) {m=1}
  pvalues=0
  p=powerSet(colnames(X))
  
  for (i in 2:length(p)) {
    
    do_we_need_to_compute_this <- function(){ #if we already did not rejected some subset of p[[i]] then it is not nessesary to continue
      we_need_it=TRUE
      
      if(i>3){for (j in 1:(i-1)) {
        if (pvalues[j]>0.05 && all(p[[j]] %in% p[[i]])) {
          we_need_it=FALSE }}}
      return(we_need_it)}
    
    if(do_we_need_to_compute_this()==TRUE){
      
      set=powerSet(1:m)[[i]]
      epsilon_hat=estimate_epsilon(Y, as.data.frame(X[,set]), constraint_set_F)
      distribution_test =  constraint_set_F=="Additive" || constraint_set_F=="Location-scale" || AndersonDarlingTest(epsilon_hat$residual, "punif")$p.value>0.05 
      new_pvalue=0
      if (epsilon_hat$significance[1]<0.05) { #Significance
        if (distribution_test==TRUE)  #Distribution
          new_pvalue=p_value_of_independence(epsilon_hat$residual, X[,set]) }
      pvalues=c(pvalues, round( new_pvalue, digits = 6) )
    }
    else{ pvalues=c(pvalues, 1 )}
  }  
  
  #since here its just handling the output into a nice form
  z=Reduce(intersect, p[pvalues>0.05])
  if(identical( Reduce(intersect, p[pvalues>0.05]), integer(0))){z="EMPTY"}
  
  p[[1]]='Empty set'
  for (i in 2:length(p)) {p[[i]] = paste(p[[i]], collapse=", ")  }
  
  result=  as.data.frame( c(pvalues, "___________", paste(z,  collapse=" ") )  )
  colnames(result)<-c("p.values")
  rownames(result)<-c(p, "__________", "Intersection")
  
  
  return(result)
}




IDS_and_score_based_estimation_of_F_parents<- function(Y, X, constraint_set_F="Additive", show_intermediate_results=FALSE){
  ###########################Some random helpful functions ##########################
  
  rename_columns_to_1_to_d <-function(X){ 
    if (is.null(nrow(X))) { return(as.data.frame(X1=X)) }
    return(X %>% rename_with(~ str_c("X", seq_along(.))))
  }
  
  ###
  
  formula <- function(X, withY=TRUE, smoothing=TRUE, type ='additive'){ #Definition of formula. We want to distinguish between discrete and continous variables for GAM
    q = ncol(X); if (is.null(nrow(X))) {q=1} 
    if(withY==TRUE){ form="Y ~ "}
    if(withY==FALSE){form=" ~ "}
    if (type =='additive') {
      
      for (i in 1:q) {#Discrete variable is if it contains <=9 different values
        if (  length(unique(X[,i]))>9  ) { 
          if (smoothing==TRUE) { form=paste(form, paste0("s(X", i,  ")+")) }
          if (smoothing==FALSE) { form=paste(form, paste0("X", i,  "+")) }
        }
        if (  length(unique(X[,i]))<=9  ) {  form=paste(form, paste0("as.factor(X", i,  ")+"))  }
      }
      
    }
    
    if (type =='general') {
      form=paste(form, paste0("s("))
      for (i in 1:q) {
        form=paste(form, paste0("X", i,  ",")) 
      }
    }
    form=substr(form,1,nchar(form)-1)
    if (type == 'general') { form=paste(form, paste0(")"))}
    
    return(as.formula(form))
  }
  
  ###
  
  are_all_significant <- function(fit, number_of_parameters=2){
    all_are_significant=TRUE
    m=str_count(toString(fit$terms[[3]]), "[0-9]") #How many numbers are in the formula; Chceme iba ziskat pocet covariates in fit
    vector = summary(fit)[[8]]
    if (number_of_parameters==2) {
      for (i in 1:m) {
        if(min(vector[i], vector[m+i])>0.01) all_are_significant=FALSE
      }}
    if (number_of_parameters==1) {
      for (i in 1:m) {
        if(vector[i]>0.01) all_are_significant=FALSE
      }}
    
    return(all_are_significant)
  }
  
  significance_maximum <- function(fit, number_of_parameters=2){
    significance=-Inf
    #m=str_count(toString(fit$terms[[3]]), "[0-9]") #How many numbers are in the formula; Chceme iba ziskat pocet covariates in fit
    vector = summary(fit)[[8]]
    m=length(vector)
    if (number_of_parameters==1) {
      for (i in 1:m) {
        if(vector[i]>significance) significance = vector[i]
      }}
    if (number_of_parameters==2) {
      for (i in 1:(m/number_of_parameters)) { number = min(vector[i], vector[ (m/number_of_parameters)+i])
      if(number>significance) significance = number
      }}
    if (number_of_parameters>2) print('Not implemented, sorry.')
    if (significance<0.001) {significance = 0.001}
    
    return(significance)
  }
  
  is_sample_uniformly_distributed<-function(x){
    return(AndersonDarlingTest(x, "punif")$p.value)
  }
  
  
  p_value_of_independence<- function(epsilon, Z, accuracy = 0.0025){
    if (is.null(nrow(Z))) { numberofcolumns = 1} else {numberofcolumns = ncol(Z)}
    x=matrix(as.matrix(Z), ncol=numberofcolumns); #x=matrix(c(X1), ncol=1); 
    y=matrix(epsilon, ncol=1)
    
    #if (numberofcolumns==1) { pvalue=hoeffding.D.test(as.vector(x),as.vector(y), na.rm = TRUE, collisions = FALSE, precision = accuracy/26)$p.value  }
    if(length(epsilon)<=1000){pvalue=dhsic.test(list(x,y))$p.value  }
    else{  pvalue=multIndepTest(data.frame(x,y) , d=rep(1, ncol(x)+1), verbose = FALSE, N=as.integer(1/(2*accuracy)-1) )$fisher.pvalue}
    
    
    return(pvalue + 0.000001) #+0.000001 is there as a cutoff value so we do not return zero value (since we later take log(0))
  }
  
  ############### Here starts the code ########################
  
  lambda1 = 1
  lambda2 = 1
  lambda3 = 1; if(constraint_set_F=="Additive" || constraint_set_F=="Location-scale") lambda3 = 0
  
  X = rename_columns_to_1_to_d(X)
  m=ncol(X); if (is.null(nrow(X))) {m=1}
  n=length(Y)
  pvalues=0; score=-Inf; independence = 0;  significance = 1; distribution = 1
  p=powerSet(colnames(X))
  
  for (i in 2:length(p)) {
    
    set=powerSet(1:m)[[i]]
    epsilon_hat=estimate_epsilon(Y, as.data.frame(X[ ,set]), constraint_set_F) #Estimate epsilon
    
    new_significance = epsilon_hat$significance[1] #Significance
    new_independence = p_value_of_independence(epsilon_hat$residual, X[,set])#Independence
    
    if(constraint_set_F=="Additive" || constraint_set_F=="Location-scale"){#Distribution
      distribution_test = 1} else distribution_test= AndersonDarlingTest(epsilon_hat$residual, "punif")$p.value
    
    
    new_score =  lambda1*(log(new_independence)) + lambda2*(-log(new_significance)) + lambda3*log( distribution_test )
    
    new_pvalue = 0
    if (new_significance<0.05 & distribution_test>0.01) {   new_pvalue = round(new_independence, digits = 6) }
    
    independence = c(independence,new_independence)
    significance = c(significance,new_significance)
    distribution = c(distribution,distribution_test)
    pvalues=c(pvalues,new_pvalue)
    score = c(score, new_score)
    if(show_intermediate_results==TRUE) cat(set, ': ', new_score,', I, S, D: ',new_independence, '||| ',new_significance, '||| ',distribution_test, "\n")
  }  
  
  #This is just handling the output into a nice form
  plausible_parents=p[pvalues>0.05];
  if(length(plausible_parents)>0){for (i in 1:length(plausible_parents)) { plausible_parents[[i]] = paste(plausible_parents[[i]], collapse = " ")}}
  plausible_parents = unlist(plausible_parents); plausible_parents = paste(plausible_parents, collapse = " ; ")

  
  identifiable_parents=paste( Reduce(intersect, p[pvalues>0.05]),collapse = " ")
  if(identical( Reduce(intersect, p[pvalues>0.05]), character(0))||identical( Reduce(intersect, p[pvalues>0.05]), NULL)){identifiable_parents="EMPTY"}
  score_estimate = "   "; for (i in 1:length(powerSet(1:m)[[which.max(score)]])) {score_estimate = paste0(score_estimate, "X", powerSet(1:m)[[which.max(score)]][i], " ")}
  
  p[[1]]='Empty set'
  for (i in 2:length(p)) {p[[i]] = paste(p[[i]], collapse=" ")  }
  
  result= round( data.frame(score = score, independence = independence, significance = significance, distribution=distribution), digits = 3)
  result[nrow(result)+1,] = c('---', '---', '---', '---')
  result[nrow(result)+1,] = c(score_estimate, ' ', ' ', ' ')
  result[nrow(result)+1,] = c(identifiable_parents, ' ' , ' ', ' ')
  result[nrow(result)+1,] = c(plausible_parents, ' ' , ' ', ' ')
  
  rownames(result)<-c(p,"---", "Score-based estimate:", "Identifiable causes:","F-plausible causes:")
  return(result)
}


