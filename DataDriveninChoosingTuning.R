library(evd)
library(MASS)
n  <- 100
p <- 100

starttime = proc.time()

nrepeat = 100 
deltan_chosen = rep(0,nrepeat)
cn_chosen = rep(0,nrepeat)

quantile = .95
cutoff = qgev(quantile,loc = - log(2*pi), scale =2, shape =0) + 4*log(p) - log(log(p)) 



Ip = diag(p)


Omega_half <- diag(p)

for (i in 1:p){
  h = sample(1:p,1)  
  Omega_half[i, h] <- .7
}

Omega_gen = Omega_half%*%t(Omega_half)



Sigma <- solve(Omega_gen)
##Ha included

Omega_gen1 <- matrix(0,p,p)

for (i in 1:p)
  for (j in 1:p)
  {
    if ( abs(i - j) < 2)
    {
      Omega_gen1[i, j] <- .8^ (abs(i - j))
    }
  }

Omega2 <- Omega_gen + Omega_gen1


###ESTIMATING PRECISON MATRIX GIVEN NETWORK STRUCTURE FUNCTION
precisionestimate_function = function(structure,samplecovariance){
  estimationmatrix = matrix(0,ncol(samplecovariance),ncol(samplecovariance))
  Ip = diag(ncol(samplecovariance))
  for (col in 1:ncol(samplecovariance)){
    index = which(structure[,col]!=0) #nonzero position in Omega
    B.index = Ip[,index]
    what1 <- solve(samplecovariance[index,index]) %*% t(B.index) %*% Ip[,col]
    what <- B.index %*% what1
    estimationmatrix[,col] = what
  }
  return(estimationmatrix)
}



# Function to calculate the maximum component-wise and compute Cn, A:sigma, B:Omega
cn_function <- function(A, B) {
    # Extract diagonal components of A and B
    sigma_ii <- diag(A)  # diagonal of A
    omega_ii <- diag(B)  # diagonal of B
    p = length(A[1,])
    # Initialize variable to store the max value for Cn calculation
    max_value <- 0
    
    # Loop over all pairs of i, j to compute the maximum value
    for (i in 1:p) {
      for (j in 1:p) {
        
        # Compute the denominator for each pair of i, j
        denominator <- sigma_ii[i] * sigma_ii[j] + 2 * A[i, j]^2  # A[i,j] = σ_ij
        
        # Compute the numerator for each pair of i, j
        numerator <-  (omega_ii[i] * sigma_ii[j] + 1)
        
        # Compute the value for each pair of i, j and update max_value
        value <- numerator / denominator
        if (value > max_value) {
          max_value <- value
        }
      }
    }
    
    # Calculate Cn using the maximum value found
    cn <- 4 * (log(p))^0.5 * max_value
    
    return(cn)
  }





##########Threshold function add on
indfunc = function(x,thre, cn){
  if(x>= thre){ x = x} else{
    x = x + cn
  }
  return(x)
}
############
indfunc_star = function(x,thre, cn){
  if(x>= thre){ x = x} else{
    x = cn
  }
  return(x)
}
###modified for choosing original Dntildle

tuning_standardize = function(matrixtun,cn,threshold){
  matrixout = matrix(0, ncol(matrixtun), ncol(matrixtun))
  for(col in 1:ncol(matrixtun)){
    index = which(abs(matrixtun[,col])>0.0000001)
    for(z in 1:length(index)){
      positionth = index[z]
      variance = (matrixtun[col,col]* matrixtun[positionth,positionth] + matrixtun[col,positionth]^2)/n
      matrixout[positionth, col] = indfunc(matrixtun[positionth, col], threshold*sqrt(variance), cn)
    }}
  return(matrixout)
}

##########
##############modified precision estimation dntildle star
##First extract a tuning sequence from data and corresponding structure
threshold_seq_func = function(mydata, structure){
  Sn_Dntild = cov(mydata)
  Omegahat = precisionestimate_function(structure,Sn_Dntild)
  matrixtun = Omegahat
  matrixdeltatun = matrix(0, ncol(matrixtun), ncol(matrixtun))
  for(col in 1:ncol(matrixtun)){
    index = which(abs(matrixtun[,col])>0.0000001)
    
    for(z in 1:length(index)){
      positionth = index[z]
      variance = (matrixtun[col,col]* matrixtun[positionth,positionth] + matrixtun[col,positionth]^2)/n
      matrixdeltatun[positionth, col] = matrixtun[positionth, col]/sqrt(variance)
    }}
  # Extract non-zero values
  non_zero_values <- matrixdeltatun[matrixdeltatun!= 0]
  threshold_sequence <- sort(abs(non_zero_values))
  
  return(threshold_sequence)
}

#apply the tuning accordingly to the values for Dntildlestar

tuning_standardize1 = function(matrixtun, cn, threshold){
  matrixdeltatun = matrix(0, ncol(matrixtun), ncol(matrixtun))
  matrixout = matrix(0, ncol(matrixtun), ncol(matrixtun))
  for(col in 1:ncol(matrixtun)){
    index = which(abs(matrixtun[,col])>0.0000001)
    
    for(z in 1:length(index)){
      positionth = index[z]
      variance = (matrixtun[col,col]* matrixtun[positionth,positionth] + matrixtun[col,positionth]^2)/n
      matrixout[positionth, col] = indfunc_star(matrixtun[positionth, col], threshold*sqrt(variance), cn)
    }}
  return(matrixout)
}

#########TEST STATISTIC FUNCTION

Dnfunction = function(mydata, precisionestimate){
  Sn <- cov(mydata)
  samplesize = dim(mydata)[1]
  center_samplecovariance <- (t(mydata) %*% mydata)/samplesize
  stdarray = matrix(0,ncol(center_samplecovariance),ncol(center_samplecovariance))
  Ip = diag(ncol(center_samplecovariance))
  for(col in 1:ncol(center_samplecovariance)){
    
    distance = center_samplecovariance%*% precisionestimate[,col] - Ip[col,]
    
    for (z in 1:ncol(center_samplecovariance)){
      if (z != col){
        var = (1/samplesize)*(precisionestimate[col,col]*center_samplecovariance[z,z])
        
        stdarray[z,col] = Ip[z,]%*%distance/sqrt(var)
      } else {
        var = (1/samplesize )*(precisionestimate[col,col]*center_samplecovariance[z,z] + 1)
        stdarray[z,col] = Ip[z,]%*%distance/sqrt(var)
      }
    }# end loop for z
  }  ##end loop for col
  return(stdarray^2)
}

######################## Calculate power or type 1 error
powerfunction = function(teststatistic_values,cutoff_point){
  
  testindicator = rep(0,length(teststatistic_values))
  for(i in 1:length(teststatistic_values)){
    if (teststatistic_values[i] > cutoff_point){
      testindicator[i] = 1
    }
  }
  return(mean(testindicator))
}

############
Dnhat_function = function(mydata,structure){
  Sn_Dnhat = cov(mydata)
  Omegahat_Dnhat = precisionestimate_function(structure,Sn_Dnhat)
  testvalue_Dnhat = max(Dnfunction(mydata, Omegahat_Dnhat))
  return(testvalue_Dnhat)
}


Dntilde_function = function(mydata,structure,cn,threshold){
  Sn_Dntild = cov(mydata)
  Omegahat_Dntild = precisionestimate_function(structure,Sn_Dntild)
  Omegahat_tun2 = tuning_standardize(Omegahat_Dntild,cn,threshold) #Tuning 1 applied
  testvalue_Dntild = max(Dnfunction(mydata, Omegahat_tun2)) 
  return(testvalue_Dntild)
}


Dntildestar_function = function(mydata,structure,cn,threshold){
  Sn_Dntild = cov(mydata)
  Omegahat_Dntild = precisionestimate_function(structure,Sn_Dntild)
  Omegahat_tun3 = tuning_standardize1(Omegahat_Dntild,cn,threshold) #Tuning 1 applied
  testvalue_Dntild_star = max(Dnfunction(mydata, Omegahat_tun3)) 
  return(testvalue_Dntild_star)
}

############## Testing 
start_time <- Sys.time()
teststatistic_include_seq = rep(0,nrepeat)
teststatistic_correct_seq = rep(0,nrepeat)
mystructure = Omega2


for (nrep in 1:nrepeat){
  X <- mvrnorm(n=n,rep(0,p),Sigma)#########Generate data
  ##Obtain the threshold tuning sequence from data and included structure
  mydata = X
  Sn = cov(mydata)
  Omegahat = precisionestimate_function(mystructure,Sn)
  cn = cn_function(Sn, Omegahat)
  cn_chosen[nrep] = cn
  
  threshold_seq = threshold_seq_func(mydata, mystructure) 
  rejection_decision = 0
  
  rejection_decision <- 0
  
  # While loop to iterate through threshold values
  i <- 1
  while (rejection_decision == 0 && i <= length(threshold_seq)) {
    
    threshold <- threshold_seq[i]
    #print(threshold)
    ##split data 20 times
    
    nsplit = 20
    Dntild_seq = rep(0,nsplit)
    Dntildstar_seq = rep(0,nsplit)
    
    
    for(k in 1:nsplit){
      n1  = n*.5
      
      trainid = sort(sample(1:n,n1))
      X_half1 = X[trainid,]
      X_half2 = X[-trainid,]
      
      Dntild_seq[k] =  Dntilde_function(X_half1, mystructure, cn, threshold) 
      Dntildstar_seq[k] =  Dntildestar_function(X_half2, mystructure, cn, threshold)
    }
    
    t_test_result <- t.test(Dntild_seq, Dntildstar_seq)
    
    # Display the result
    #print(t_test_result)
    
    # Define significance level
    alpha <- 0.05
    
    # Rejection decision: 1 for reject, 0 for not reject
    rejection_decision <- ifelse(t_test_result$p.value < alpha, 1, 0)
    
    
    # Move to the next threshold
    i <- i + 1
  }
  
  deltan_chosen[nrep] = threshold
  
  ###Come back to apply the Dntildle test on the included case for the original data
  
  teststatistic_include = Dntilde_function(mydata, mystructure, cn, threshold) 
  teststatistic_include_seq[nrep] =  teststatistic_include
  
  teststatistic_correct = Dntilde_function(mydata, Omega_gen, cn, threshold) 
  teststatistic_correct_seq[nrep] = teststatistic_correct
  
}


powerfunction = function(teststatistic_values,cutoff_point){
  
  testindicator = rep(0,length(teststatistic_values))
  for(i in 1:length(teststatistic_values)){
    if (teststatistic_values[i] > cutoff_point){
      testindicator[i] = 1
    }
  }
  
  return(mean(testindicator))
}


##rejection rate
RR_include = powerfunction(teststatistic_include_seq, cutoff)
RR_correct = powerfunction(teststatistic_correct_seq, cutoff)
totaltime = proc.time() - starttime

