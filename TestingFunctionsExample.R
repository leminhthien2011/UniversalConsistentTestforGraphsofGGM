args <- commandArgs(trailingOnly = TRUE)
n  <- as.integer(args[1])
p <- as.integer(args[2])
s0 <- as.integer(args[3])

n = 200
p =200
s0 = 4
library(evd)
library(MASS)

thresholdpower = 3
cn = 30 
starttime = proc.time()
nrepeat = 100 
quantile = .95
cutoff = qgev(quantile,loc = - log(2*pi), scale =2, shape =0) + 4*log(p) - log(log(p)) 

threshold = log(n)^(1/thresholdpower)
cn = cn/100


Ip = diag(p)

Omega_half <- diag(p)

for (i in 1:p){
  h = sample(1:p,1)  
  Omega_half[i, h] <- .8
}

Omega = Omega_half%*%t(Omega_half)


##Ha1 nested, first column of Omega_half nested

Omega_half1 = Omega_half
Omega_half1[1,] = Ip[1,]
Omega1 = Omega_half1%*%t(Omega_half1)

##Ha2 included
first_index = which(Omega_half[1,]!=0)
remain = setdiff(1:p,first_index)
top_on = sample(remain,1)
Omega_half2 = Omega_half
Omega_half2[1,top_on] = .8
Omega2 = Omega_half2%*%t(Omega_half2)





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




indfunc = function(x,thre, cn){
  if(x>= thre){ x = x} else{
    x = x + cn
  }
  return(x)
}

##modified precision estimation tuning standardize
tuning_standardize = function(matrixtun,cn){
  matrixout = matrix(0, ncol(matrixtun), ncol(matrixtun))
  for(col in 1:ncol(matrixtun)){
    
    index = which(abs(matrixtun[,col])>0.0000001)
    for(z in 1:length(index)){
      positionth = index[z]
      variance = (matrixtun[col,col]* matrixtun[positionth,positionth] + matrixtun[col,positionth]^2)/n
      
      matrixout[positionth, col] = indfunc(matrixtun[positionth, col], threshold*sqrt(variance), cn)
    }
    
  }
  return(matrixout)
}
#############
################# TEST STATISTIC FUNCTION

Dnfunction = function(center_samplecovariance, precisionestimate){
  stdarray = matrix(0,ncol(center_samplecovariance),ncol(center_samplecovariance))
  Ip = diag(ncol(center_samplecovariance))
  for(col in 1:ncol(center_samplecovariance)){
    
    distance = center_samplecovariance%*% precisionestimate[,col] - Ip[col,]
    
    for (z in 1:ncol(center_samplecovariance)){
      if (z != col){
        var = (1/n)*(precisionestimate[col,col]*center_samplecovariance[z,z])
        
        stdarray[z,col] = Ip[z,]%*%distance/sqrt(var)
      } else {
        var = (1/n)*(precisionestimate[col,col]*center_samplecovariance[z,z] + 1)
        stdarray[z,col] = Ip[z,]%*%distance/sqrt(var)
      }
    }# end loop for z
  }  ##end loop for col
  
  return(stdarray^2)
  
}


#STEP 1: Generate data
Sigma <- solve(Omega)
X <- mvrnorm(n=n,rep(0,p),Sigma)

#STEP 2: TESTING PROCEDURE 
Sn <- cov(X)
M11 <- t(X) %*% X
M1 <- M11/n


Omega_H0 = Omega #Structure want to test
################# INTUITIVE WAY, FREE TUNING
Omegahat = precisionestimate_function(Omega_H0,Sn) #estimating the precision matrix given Omega structure
teststatistic_naive = max(Dnfunction(M1, Omegahat)) # Test statistic for replication r

###Modified Test Statistic
Omegahat_tun = tuning_standardize(Omegahat,cn) #Tuning 1 applied
teststatistic_modified = max(Dnfunction(M1, Omegahat_tun)) # Test statistic3 for replication r

##Decision function

decision_function <- function(X,cutoff) {
  
  if (X > cutoff) {
    return("Reject the null hypothesis")
  } else {
    return("Fail to reject the null hypothesis")
  }
}


decision_function(teststatistic_naive,cutoff)
decision_function(teststatistic_modified,cutoff)
