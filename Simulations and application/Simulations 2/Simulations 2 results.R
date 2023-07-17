result_IDENTIFIABLE = c()
result_SCORE = c()
result_QCCD = c()
result_IGCIGauss = c()
result_Slope = c()

set.seed(0)
n=500
number_of_repetitions = 50



performance_of_algorithm <- function(result_graphs, benchmark = 'first', metric = 'correct_parents'){
  k = length(result_graphs)
 
  if (benchmark == 'first') {
  
  if( metric == 'wrong_parents'){
    count = 0
    for(i in 1:k){graph = result_graphs[[i]]
    if (identical(graph[2:4],c(0,0,0))==FALSE) {count = count+1}
    }
    return(1-count/k) }
  
  if( metric == 'correct_parents'){
    count = 0
    for(i in 1:k){graph = result_graphs[[i]]
    if (identical(graph[1],c(1))==TRUE) {count = count+1}
    }
    return(count/k) }
  }
  
  if (benchmark == 'second') {
    
    if( metric == 'wrong_parents')return(0)
    
    if( metric == 'correct_parents'){
      count = 0
      for(i in 1:k){graph = result_graphs[[i]]
      count = count + sum(graph)/3}
      return(count/k) 
    }}
    
    if (benchmark == 'third') {
      
      if( metric == 'wrong_parents'){
        count = 0
        for(i in 1:k){graph = result_graphs[[i]]; correct_graph = third_benchmark_dataset[[i]]$graph
        if (all(correct_graph - graph >= 0)==FALSE)count = count + 1        }
        return(1-count/k) 
      }}
        
      
      if( metric == 'correct_parents'){
        count = 0
        number_of_empty_graphs = 0
        for(i in 1:k){graph = result_graphs[[i]]; correct_graph = third_benchmark_dataset[[i]]$graph
        if(all(correct_graph==c(0,0,0))){number_of_empty_graphs = number_of_empty_graphs+1}
        else{count = count + sum(graph[correct_graph==1])/ sum(correct_graph[correct_graph==1]) }    }
        return(count/(k - number_of_empty_graphs)) 
        }
}








#########Our IDS and scoring algorithms############################################################
change_format_of_result_for_our_algorithm <- function(x, length){
  graph = rep(0,length)
  for (i in 1:length) { if(grepl(i, x)) graph[i] = 1 }
  return(graph)
}
for (j in 1:number_of_repetitions) {
  data = first_benchmark_dataset[[j]]
  Y = data$Y
  if(ncol(data)-1 ==3){X1=data$X1; X2=data$X2; X3=data$X3; X = data.frame(X1, X2, X3)}
  if(ncol(data)-1 ==4){X1=data$X1; X2=data$X2; X3=data$X3; X4=data$X4; X = data.frame(X1, X2, X3, X4)}
  
  graph = IDS_and_score_based_estimation_of_F_parents(Y=Y, X=X, constraint_set_F="Additive")
  print(graph);  print(j)
  graph_score=  graph[nrow(graph)-1,][1]
  graph_identifiable=  graph[nrow(graph),][1]
  
  result_SCORE[[j]] =   change_format_of_result_for_our_algorithm(graph_score, length = ncol(X))
  result_IDENTIFIABLE[[j]] =   change_format_of_result_for_our_algorithm(graph_identifiable, length = ncol(X))
  }  

#for (j in 1:number_of_repetitions) {
#  cat(result_IDENTIFIABLE[[j]] , ', ', result_SCORE[[j]], ', graph = ',third_benchmark_dataset[[j]]$graph, '\n')
#}
performance_of_algorithm(result_IDENTIFIABLE,benchmark = 'first', metric = "correct_parents")
performance_of_algorithm(result_IDENTIFIABLE,benchmark = 'first', metric = "wrong_parents")

performance_of_algorithm(result_SCORE,benchmark = 'first', metric = "correct_parents")
performance_of_algorithm(result_SCORE,benchmark = 'first', metric = "wrong_parents")




#########QCCD############################################################
result_QCCD = c()
for (j in 1:number_of_repetitions) {
  data = third_benchmark_dataset[[j]]$data
  
  Y = data$Y
  if(ncol(data)-1 ==3){X1=data$X1; X2=data$X2; X3=data$X3; X = data.frame(X1, X2, X3)}
  if(ncol(data)-1 ==4){X1=data$X1; X2=data$X2; X3=data$X3; X4=data$X4; X = data.frame(X1, X2, X3, X4)}
  
  graph=c()
  for (i in 1:(ncol(data)-1) ) {
    Z = data.frame(X[,i], Y)
    method = QCCD(Z, m=1)$cd
    graph = c(graph, method)
  }
  result_QCCD[[j]] = graph
  cat('Time remaining: ',number_of_repetitions-j, '\n')
}  

for (j in 1:number_of_repetitions) {
  cat( result_QCCD[[j]], ',  ', third_benchmark_dataset[[j]]$graph, '\n' )
}

performance_of_algorithm(result_QCCD,benchmark = 'third', metric = "correct_parents")
performance_of_algorithm(result_QCCD,benchmark = 'third', metric = "wrong_parents")
######### IGCI ############################################################
result_IGCIGauss = c()
for (j in 1:number_of_repetitions) {
  data = third_benchmark_dataset[[j]]$data

  Y = data$Y
  if(ncol(data)-1 ==3){X1=data$X1; X2=data$X2; X3=data$X3; X = data.frame(X1, X2, X3)}
  if(ncol(data)-1 ==4){X1=data$X1; X2=data$X2; X3=data$X3; X4=data$X4; X = data.frame(X1, X2, X3, X4)}
  
  graph=c()
  for (i in 1:(ncol(data)-1)) {
    Z = data.frame(X[,i], Y)
    method = IGCIWrap_Unif(Z)$cd #IGCIWrap_G(Z)$cd
    graph = c(graph, method)
  }
  result_IGCIGauss[[j]] = graph
  cat('Time remaining: ',j, '\n')
}
performance_of_algorithm(result_IGCIGauss,benchmark = 'third', metric = "correct_parents")
performance_of_algorithm(result_IGCIGauss,benchmark = 'third', metric = "wrong_parents")


#########Slope############################################################
result_Slope = c()
for (j in 1:number_of_repetitions) {
  data = first_benchmark_dataset[[j]]
  
  Y = data$Y
  if(ncol(data)-1 ==3){X1=data$X1; X2=data$X2; X3=data$X3; X = data.frame(X1, X2, X3)}
  if(ncol(data)-1 ==4){X1=data$X1; X2=data$X2; X3=data$X3; X4=data$X4; X = data.frame(X1, X2, X3, X4)}
  
  graph=c()
  for (i in 1:(ncol(data)-1) ) {
    Z = data.frame(X[,i], Y)
    method = SlopeWrap(Z)$cd
    graph = c(graph, method)
  }
  result_Slope[[j]] = graph
  cat('Time remaining: ',number_of_repetitions-j, '\n')
}  


performance_of_algorithm(result_Slope,benchmark = 'first', metric = "correct_parents")
performance_of_algorithm(result_Slope,benchmark = 'first', metric = "wrong_parents")



###############RESIT###################################
df <- data.frame(position = 1, first_benchmark_dataset[[1]])
for (i in 2:number_of_repetitions) {df <- rbind(df, data.frame(position = i, first_benchmark_dataset[[i]]))}

write_feather(df, "first_benchmark")

df <- data.frame(position = 1, second_benchmark_dataset[[1]])
for (i in 2:number_of_repetitions) {df <- rbind(df, data.frame(position = i, second_benchmark_dataset[[i]]))}

write_feather(df, "second_benchmark")

df <- data.frame(position = 1, third_benchmark_dataset[[1]]$data)
for (i in 2:number_of_repetitions) {df <- rbind(df, data.frame(position = i,  third_benchmark_dataset[[i]]$data))}
third_graphs = data.frame( rbind( third_benchmark_dataset[[1]]$graph))
for (i in 2:number_of_repetitions) {third_graphs <- rbind(third_graphs, data.frame( rbind( third_benchmark_dataset[[i]]$graph)))}

write_feather(df, "third_benchmark")
write_feather(third_graphs, "third_benchmark_graphs")



