# In this file you can find the three benchmark datasets, as well as the functions that generate random d-dimensional functions using Perlin noise approach
# If you ran the whole script (it will take a few minutes), you will receive  first_benchmark_dataset, second_benchmark_dataset, third_benchmark_dataset


library(raster)
library(ambient)
library(dplyr)
library(rgl)


random_function_1d<-function(freq = 0.15){
  s = seq(-10, 10, length.out = 1001)
  f <- long_grid(s)
  f$noise <- rep(0, length(s))
  while (all(f$noise[300:600]==0)) {f$noise <- gen_perlin(f$x, frequency = 0.15, fractal = 'rigid-multi')}
  amplitude = 20 / (max(f$noise)+1)
  return(f$noise*amplitude)
}


random_function_2d<-function(freq = 0.001){
  f = noise_perlin(c(1001, 1001), frequency= freq, fractal ='fbm', octaves = 2, lacunarity=2, gain=0.4)
  f=f^4
  amplitude = 20 / max(f)
  f=amplitude*f
  return(f) }


random_function_3d<-function(){
  f <- long_grid(seq(-10, 10, length.out = 50),seq(-10, 10, length.out = 50), seq(-10, 10, length.out = 50))
  f$noise <- 1000*(gen_perlin(f$x, f$y, f$z, frequency = 0.01, fractal = 'fbm'))^2+ 5000*(gen_perlin(f$x, f$y, f$z, frequency = 0.01, fractal = 'fbm'))^2
  f$noise = f$noise +abs(f$x*f$y*f$z)^(0.5)+f$x*f$y*f$z #This is here to normalize gen_perlin a bit. 
  return(f) }


random_function_2d_ugly_v2<-function(){
  f1 = noise_perlin(c(1001, 1001), frequency= 0.00001, fractal ='fbm', octaves = 2, lacunarity=2, gain=0.4)
  f1=f1^4
  amplitude = 20 / max(f1)
  f1=amplitude*f1
  
  f2 = noise_perlin(c(1001, 1001), frequency= 0.00001, fractal ='fbm', octaves = 2, lacunarity=2, gain=0.4)
  f2=f2^2
  amplitude = 20 / max(f2)
  f2=amplitude*f2
  
  f = f1
  for (i in 1:1001) { for (j in 1:1001) {
    if (j<500) {f[i,j] = f1[i,j]} 
    if (j>=500) {f[i,j] = f2[i,j]} 
    f[i,j] = f[i,j] + (j-500)/1000
  }}
  return(f) 
}


random_function_2d_ugly<-function(){
  f = noise_perlin(c(1001, 1001), frequency= 0.001, fractal ='fbm', octaves = 2, lacunarity=2, gain=0.4)
  f=f^4
  amplitude = 20 / max(f)
  f=amplitude*f
  
  f1 = noise_perlin(c(1001, 1001), frequency= 0.001, fractal ='fbm', octaves = 2, lacunarity=2, gain=0.4)
  f1=f1^4
  amplitude = 20 / max(f1)
  f1=abs( - amplitude*f1 )
  
  for (i in 1:1001) { for (j in 1:1001) {
    k = round(j/2)+1; 
    if (j<500)  {f[i,j] = f[i,k] + abs( rnorm(1,0,f[i,k]/10) )} 
    if (j>=500) {f[i,j] = f1[i,k]+ abs( rnorm(1,0,f1[i,k]/10) )} 
  }}
  return( f ) 
}



evaluation_of_f_1d<-function(f, X1){
  minim=min(X1)
  maxim = max(X1)
  
  result=c()
  for (i in 1:length(X1)) {
    result = c(result,  f[round( (X1[i]-minim)*1000/(maxim-minim) + 1 )   ])
  }
  return(result)
}

evaluation_of_f_2d<-function(f, X1,X2){
  
  function_evaluated_in_x_y<-function(f, x,y){
    
    xx = round( (x-min(X1))*1000/(max(X1)-min(X1)) + 1 )
    yy = round( (y-min(X2))*1000/(max(X2)-min(X2)) + 1 )
    
    return( f[xx,yy]   )
  }
  
  
  result=c()
  for (i in 1:length(X1)) {
    result = c(result, function_evaluated_in_x_y(f, X1[i], X2[i])  )
  }
  return(result)
}

evaluation_of_f_3d<-function(f,X1,X2,X3){
  
  function_evaluated_in_x_y_z<-function(f,x,y,z){
    index <- which.min( sqrt((f$x - x)^2 + (f$y - y)^2 + (f$z - z)^2) )
    return(f$noise[index])
  }
  
  result=c()
  for (i in 1:length(X1)) {
    result = c(result, function_evaluated_in_x_y_z(f, X1[i], X2[i],X3[i])  )
  }
  return(result)
}                      



#################Generating benchmark datasets##############
set.seed(0)
n=500
number_of_repetitions = 50

#################First benchmark##############
first_benchmark_dataset = c()
for (k in 1:number_of_repetitions) {
  
  f_Y=random_function_1d()
  f_2 = random_function_2d_ugly(); 
  f_3 = random_function_2d_ugly()
  f_4 = random_function_2d_ugly()
  dependence=rbinom(n, 1, 0.5)
  
  v=data.frame()
  for (j in 1:3) {
    for (i in 1:n) {
      v[i,j] = dependence[i]*runif(1, 0, 0.5)+(1-dependence[i])*runif(1, 0.5, 1)
    } }
  v[,4] = runif(n,0,1)
  
  X1 =  v[,1]
  Y  =  evaluation_of_f_1d(f_Y, X1)  + rnorm(n) 
  X2 =  evaluation_of_f_2d(f_2, Y, v[,2])  
  X3 =  evaluation_of_f_2d(f_3, Y, v[,3] ) 
  X4 =  evaluation_of_f_2d(f_4, Y, v[,4])  
  
  first_benchmark_dataset[[k]] = data.frame(X1, X2, X3, X4, Y)
}


#################Second benchmark##############
c=0.5
Sigma <- matrix(c(1,c,c,c,1,c,c,c,1),3,3)
second_benchmark_dataset = c()
for (j in 1:number_of_repetitions) {
  
  v = mvrnorm(n, rep(0, 3), Sigma)
  
  X1 = v[,1]
  X2 = v[,2]
  X3 = v[,3]
  f_Y=random_function_3d()
  Y = evaluation_of_f_3d(f_Y, X1, X2, X3)  + rnorm(n) 
  
  second_benchmark_dataset[[j]] = data.frame(X1, X2, X3, Y)
}





#################Third benchmark Pareto##############
n=500
third_benchmark_dataset = c()

generate_pareto <- function(thetas){
  result=c()
  for (i in 1:n) {
    result=c(result, rPareto(1, 1, alpha = thetas[i]));
  }
  return(result)
}

for (j in 1:number_of_repetitions) {
  graph = rbinom(3,1,0.5) #graph[i]=0 if Y-->Xi and graph[i]=1 if Xi-->Y 
  
  X=data.frame(rep(0,n),rep(0,n),rep(0,n))
  
  for (i in 1:3) { if (graph[i]==1) {X[,i] = rnorm(n)} }
  
  if (sum(graph)==0) { theta=abs(random_function_1d())+1;
  Y =generate_pareto( evaluation_of_f_1d(theta, rnorm(n)))
  }
  
  if (sum(graph)==1) { theta=abs(random_function_1d())+1; parents=X[,which(graph==1)]
  Y =generate_pareto(evaluation_of_f_1d(theta, parents))
  }
  
  if (sum(graph)==2) { theta=random_function_2d()+1; parents=X[,which(graph==1)]
  Y =generate_pareto( evaluation_of_f_2d(theta, parents[[1]],  parents[[2]]) )}
  
  if (sum(graph)==3) { theta=abs(random_function_3d())+1;  parents=X[,which(graph==1)]
  Y =generate_pareto( evaluation_of_f_3d(theta, parents[[1]],  parents[[2]],  parents[[3]]))}
  
  for (i in 1:3) if (graph[i]==0){ 
    f=random_function_2d_ugly();
    X[,i] = evaluation_of_f_2d(f, Y, runif(n))
  }
  
  
  X1 = X[,1];X2=X[,2];X3=X[,3];
  
third_benchmark_dataset[[j]]=list(data = data.frame(X1, X2, X3, Y), graph = graph)
}
