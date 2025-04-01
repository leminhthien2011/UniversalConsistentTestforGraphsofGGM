#### demo code
library(Matrix)
library(MASS)
library(tidyr)
library(stringr)
library(glasso)
library(ggplot2)
library(ggsci)
library(flare)
library(igraph)
library(evd)

data = read.csv("COVID19_US2021_weeklycases.csv", header =T)
states = data[,1]
data00 = data[,-1]
dim(data00)

data0 = matrix(NA,51,52)


percentage = .8
data0 = data00/1000
dim(data0)
trainingsize = 12
sampleindex = sort(sample(1:51,trainingsize,replace=F))


data01 = data0[sampleindex,] #train data
data02 =data0[-sampleindex,] #test data
####Original
data01  = data0
data02 = data0


p = dim(data0)[2]


##Fourth with Independence
Omega1 = diag(p)
#H=diag(p)
s0=3
Omega2 = matrix(0,p,p)
for(i in 1:p){
  for(j in 1:p){
    if(abs(i-j)<s0){
      Omega2[i,j] = .6^{i-j}
    }
    
  }
}

Omega2 = ifelse(abs(Omega2)>0,1,0)
#########


###Block matrix
M1 <- Matrix(1, 4,4)
M <- as.matrix(kronecker(Diagonal(13), M1))
Omega3 = M
Omega3 = Omega3 + t(Omega3)
Omega3= ifelse(abs(Omega3)>0,1,0)
#########


#Y = beta0 + beta1*X1 + beta2*X2 + beta3*X3 + Epsilon

#gamma is the tuning parameter to tune the weight given to data driven
#information V and prior network topology H
#gamma=1: we fully rely on prior network topology H, 

#############Test statistic value of the pre-specified structure
network_test =  function(X, Omega){
  n = dim(X)[1]
  p = dim(X)[2]
  Ip =diag(p)
  stdarray = rep(0,p)
  max2network = rep(0,p)
  Sn <- cov(X)
  M11 <- t(X) %*% X
  M1 <- M11/n
  for (col in 1:p){
    ##Our Est
    index = which(Omega[,col]!=0)
    B.index = Ip[,index]
    what1 <- solve(Sn[index,index]) %*% t(B.index) %*% Ip[,col]
    
    ############
    what <- B.index %*% what1
    distance = M1%*%what - Ip[col,]
    
    for (z in 1:p){
      if (z != col){
        var = (1/n)*(what[col]*M1[z,z])
        estimation_ij = Ip[z,]%*%distance/sqrt(var)
        stdarray[z] = estimation_ij  
      } else {
        var = (1/n)*(what[col]*M1[col,col] + 1)
        estimation_ij = Ip[z,]%*%distance/sqrt(var)
        stdarray[z] = estimation_ij 
        
      }
      
    }
    
    max2network[col] = max(stdarray^2)
    
  }
  
  return(max(max2network))
}

#######Estimate the precision matrix with a given network structure


precision_estimation = function(X, Omega){
  precision_est = matrix(0,dim(X)[2],dim(X)[2]) 
  Ip = diag(dim(X)[2])
  Sn <- cov(X)
  for (col in 1:dim(X)[2]){
    ##Our Est
    index = which(Omega[,col]!=0)
    B.index = Ip[,index]
    what1 <- solve(Sn[index,index]) %*% t(B.index) %*% Ip[,col]
    ############
    what <- B.index %*% what1
    precision_est[,col] = what
  }
  return(precision_est)
}

Y1 = scale(data01, scale = FALSE)
test1 = network_test(Y1,Omega1)
test2 = network_test(Y1,Omega2)
test3 = network_test(Y1,Omega3)


c(test1,test2,test3)


##

data1 = data02
n = dim(data1)[1]
p = dim(data1)[2]
##


ID = kronecker(1:n, rep(1,p)) 
timepoint = kronecker(rep(1,n),1:p)

cases =  rep(0,length(ID)) 
for(i in 1: n){
  j1 = (i-1)*p+1
  j2 = (i-1)*p+p
  cases[j1:j2] = as.numeric(data1[i,])
}


data2 = data.frame(ID, timepoint, cases)
region.ne = c( "CT","MA","ME", "NH", "NJ" , "NY", "PA", "VT", "RI" ) # North East
region.mw = c( "IL", "IN", "IA", "KS",  "MI", "MN",  "MO", "ND", "OH", "SD","WI","NE") # Midwest
region.s = c( "AL", "AR", "DE", "FL", "GA","LA","MD","NC","KY","VA","WV","MS","TN",
              "TX","SC","DC","OK") #South
region.w = c( "AZ",  "CA", "CO","ID","OR", "WA", "NV","NM", "UT" , "WY","MT","HI","AK") #West






region.ne.num = which(states %in% region.ne)
region.mw.num = which(states %in% region.mw)
region.s.num = which(states %in% region.s)
region.w.num = which(states %in% region.w)
region.northeast = region.midwest = region.south = region.west = rep(0,length(ID)) #use south as reference

for (i in 1:length(ID)){
  if (ID[i] %in% region.ne.num ){
    region.northeast[i] = 1
  }
  if (ID[i] %in% region.mw.num ){
    region.midwest[i] = 1
  }
  if (ID[i] %in% region.w.num ){
    region.west[i] = 1
  }
  if (ID[i] %in% region.s.num ){
    region.south[i] = 1
  }
  
}

data2$region1 = region.northeast
data2$region2 = region.midwest
data2$region3 = region.south
data2$region4 = region.west




#################### Work with AR structure
Omega11 = Omega1
diag(Omega11)  = 0 #modified to match with the RAND function
gamma1= 1 #gamma=1: we fully rely on prior network topology H, 
names(data2)



precision_matrix1 = as.matrix(precision_estimation(data02, Omega1))

precision_matrix2 = as.matrix(precision_estimation(data02, Omega2))


precision_matrix3 = as.matrix(precision_estimation(data02, Omega3))


GEEfunction_precision =function (formula = formula(data), id = id, data = parent.frame(), b = NULL, tol = 1e-08, maxiter = 1000, family = gaussian, corstr = "IND", invfun = "solve", gamma=NULL, H=NULL, trace=F) 
{
  if ((invfun != "solve") && (invfun != "ginv") && (invfun != "shrink")  ) {
    stop("Unknown inverse function. Only solve or ginv or shrink.")
  }
  if (invfun == "ginv") 
    library(MASS)
  call <- match.call()
  m <- match.call(expand = FALSE)
  m$b <- m$tol <- m$maxiter <- m$link <- m$varfun <- m$corstr <- m$family <- m$invfun <-m$gamma <-m$H <-m$trace <- NULL
  
  if (is.null(m$id)) 
    m$id <- as.name("id")
  m[[1]] <- as.name("model.frame")
  m <- eval(m, parent.frame())
  Terms <- attr(m, "terms")
  y <- as.matrix(model.extract(m, "response"))
  x <- model.matrix(Terms, m)
  QR <- qr(x)
  if (QR$rank < ncol(x)) 
    stop("rank-deficient model matrix")
  if (dim(y)[2] > 1) 
    stop("Only one response (1 column)")
  
  id <- model.extract(m, id)
  if (is.null(id)) {
    stop("Id variable not found")
  }
  nobs <- nrow(x)
  np <- ncol(x)
  xnames <- dimnames(x)[[2]]
  if (is.null(xnames)) {
    xnames <- paste("x", 1:np, sep = "")
    dimnames(x) <- list(NULL, xnames)
  }
  if (is.character(family)) 
    family <- get(family)
  if (is.function(family)) 
    family <- family()
  if (!is.null(b)) {
    beta <- matrix(as.double(b), ncol = 1)
    if (nrow(beta) != np) {
      stop("Dim beta != ncol(x)")
    }
    message("user's initial regression estimate")
    print(beta)
  }
  else {
    mm <- match.call(expand = FALSE)
    mm$b <- mm$tol <- mm$maxiter <- mm$link <- mm$varfun <- mm$corstr <- mm$id <- mm$invfun <-mm$gamma <- mm$H <-mm$trace<- NULL
    mm[[1]] <- as.name("glm")
    beta <- eval(mm, parent.frame())$coef
    beta <- as.numeric(beta)
  }
  if (length(id) != length(y)) 
    stop("Id and y not same length")
  maxiter <- as.integer(maxiter)
  links <- c("identity", "log", "logit", "inverse")
  fams <- c("gaussian", "poisson", "binomial", "Gamma")
  varfuns <- c("constant", "mu", "mu(1-mu)", "mu^2")
  corstrs <- c("IND", "RAND")
  famv <- match(family$family, fams, -1)
  dist <- family$family
  linkv <- as.integer(match(c(family$link), links, -1))
  if (famv < 1) 
    stop("unknown family")
  if (famv <= 4) 
    varfunv <- famv
  else varfunv <- match(family$varfun, varfuns, -1)
  varfunv <- as.integer(varfunv)
  corstrv <- as.integer(match(corstr, corstrs, -1))
  if (linkv < 1) 
    stop("unknown link.")
  if (varfunv < 1) 
    stop("unknown varfun.")
  if (corstrv < 1) 
    stop("unknown corstr.")
  y <- as.matrix(y)
  x <- as.matrix(x)
  obs <- lapply(split(id, id), "length")
  nobs <- as.numeric(obs)
  nsub <- length(nobs)
  np <- dim(x)[[2]]
  time1 <- date()
  
  
  betadiff <- 1
  iteration <- 0
  betanew <- beta
  
  while (betadiff > tol && iteration < maxiter) {
    beta <- betanew
    
    sumg <- matrix(rep(0, np), nrow = np)
    sumc <- matrix(rep(0, np * np), nrow = np)
    arsumg <- matrix(rep(0, np), nrow = np)
    arsumc <- matrix(rep(0, np * np), nrow = np)
    gi <- matrix(rep(0, np), nrow = np)
    arsumgfirstdev <- matrix(rep(0, np * np), nrow = np)
    firstdev <- matrix(rep(0, np * np), nrow = np)
    
    
    
    
    n1 <- nobs[1]
    m0 <- diag(n1)
    UN.G <- matrix(rep(0, n1 * n1), n1)
    loc1 <- 0
    loc2 <- 0
    for (i in 1:nsub) {
      loc1 <- loc2 + 1
      loc2 <- loc1 + nobs[i] - 1
      yi <- as.matrix(y[loc1:loc2, ])
      xi <- x[loc1:loc2, ]
      ni <- nrow(yi)
      
      ui <- xi %*% beta
      
      
      
      UN.G <- UN.G + (yi - ui) %*% t(yi - ui)
    }
    UN.G <- 1/nsub * UN.G
    
    
    loc1 <- 0
    loc2 <- 0
    for (i in 1:nsub) {   
      loc1 <- loc2 + 1
      loc2 <- loc1 + nobs[i] - 1
      yi <- as.matrix(y[loc1:loc2, ])
      xi <- x[loc1:loc2, ]
      ni <- nrow(yi)
      m0 <- diag(ni)
      m1 <- matrix(rep(0, ni * ni), ni)
      ui <- xi %*% beta
      fui <- ui
      fui_dev <- diag(ni)
      vui <- diag(ni)
      
      #H1=gamma*vui %*% H %*% vui+(1-gamma)*UN.G #52*52 MATRIX
      H1 = H #precision_matrix
      m1=H1	
      wi <- t(xi) %*% fui_dev %*% vui %*% m0 %*% vui
      zi <- t(xi) %*% fui_dev %*% H1
      gi1 <- (1/nsub) * zi %*% (yi - ui)
      gi[1:np, ] <- gi1
      arsumc <- arsumc + gi %*% t(gi)
      arsumg <- arsumg + gi
      di1 <- -(1/nsub) * zi %*% fui_dev %*% xi
      firstdev[1:np, ] <- di1
      arsumgfirstdev <- arsumgfirstdev + firstdev
      
    }
    
    arcinv = solve(arsumc)  
    Q <- t(arsumg) %*% arcinv %*% arsumg
    arqif1dev <- t(arsumgfirstdev) %*% arcinv %*% arsumg
    arqif2dev <- t(arsumgfirstdev) %*% arcinv %*% arsumgfirstdev
    Godambe_matrix =  arqif2dev 
    invarqif2dev <- solve(arqif2dev)             ###########  ginv->solve: if solve does not work, you can try ginv  
    betanew <- beta - invarqif2dev %*% arqif1dev
    betadiff <- abs(sum(betanew - beta))
    iteration <- iteration + 1
  } #end while loop updating beta  
  
  time2 <- date()
  fit <- list()
  
  fit$nobs <- nobs
  fit$iteration <- iteration
  fit$coefficients <- as.vector(beta)
  names(fit$coefficients) <- xnames
  fit$linear.predictors <- as.matrix(x %*% beta)
  if (dist == "gaussian") {
    mu <- x %*% beta
    pearson <- y - mu
  }
  
  fit$fitted.value <- as.matrix(mu)
  fit$residuals <- y - mu
  fit$pearson.resi <- pearson
  fit$scale <- sum(pearson^2)/(length(y) - length(beta))
  pvalue <- 1 - pchisq(Q, np)
  betase <- sqrt(diag(invarqif2dev))
  Z <- as.vector(beta)/betase
  betapvalue <- 2 * (1 - pnorm(abs(Z)))
  parameter <- cbind(beta, betase, Z, betapvalue)
  dimnames(parameter)[[2]] <- c("estimate", "stderr", "Z", "pvalue")
  dimnames(parameter)[[1]] <- xnames
  fit$parameter <- parameter
  dimnames(invarqif2dev)[[1]] <- xnames
  dimnames(invarqif2dev)[[2]] <- xnames
  
  
  return(fit)
}



tmp1 = GEEfunction_precision(cases~region1 +region2 +region4, id = ID, data = data2, family = gaussian, corstr = "IND", gamma=gamma1, H=precision_matrix1)

tmp2 = GEEfunction_precision(cases~region1 +region2 +region4, id = ID, data = data2, family = gaussian, corstr = "IND", gamma=gamma1, H=precision_matrix2)

tmp3 = GEEfunction_precision(cases~region1 +region2 +region4, id = ID, data = data2, family = gaussian, corstr = "IND", gamma=gamma1, H=precision_matrix3)


test1 = network_test(Y1,Omega1)
##Test band
test2 = network_test(Y1,Omega2)
##Test block
test3 = network_test(Y1,Omega3)



c(test1, test2, test3)
tmp1$parameter
tmp2$parameter
tmp3$parameter

nboots = 100

sd1mat = sd2mat= sd3mat = matrix(0, nboots, 3)
subsamplesize = floor(percentage*nrow(data02))

for(b1 in 1:nboots){
  mysample =  sample(1:nrow(data02),subsamplesize,replace=F)  
  data1 = data02[mysample,]
  n = dim(data1)[1]
  p = dim(data1)[2]
  ##
  
  
  ID = kronecker(1:n, rep(1,p)) 
  timepoint = kronecker(rep(1,n),1:p)
  
  cases =  rep(0,length(ID)) 
  for(i in 1: n){
    j1 = (i-1)*p+1
    j2 = (i-1)*p+p
    cases[j1:j2] = as.numeric(data1[i,])
  }
  
  
  data2 = data.frame(ID, timepoint, cases)
  region.ne = c( "CT","MA","ME", "NH", "NJ" , "NY", "PA", "VT", "RI" ) # North East
  region.mw = c( "IL", "IN", "IA", "KS",  "MI", "MN",  "MO", "ND", "OH", "SD","WI","NE") # Midwest
  region.s = c( "AL", "AR", "DE", "FL", "GA","LA","MD","NC","KY","VA","WV","MS","TN",
                "TX","SC","DC","OK") #South
  region.w = c( "AZ",  "CA", "CO","ID","OR", "WA", "NV","NM", "UT" , "WY","MT","HI","AK") #West
  
  
  
  
  
  
  region.ne.num = which(states %in% region.ne)
  region.mw.num = which(states %in% region.mw)
  region.s.num = which(states %in% region.s)
  region.w.num = which(states %in% region.w)
  region.northeast = region.midwest = region.south = region.west = rep(0,length(ID)) #use south as reference
  
  for (i in 1:length(ID)){
    if (ID[i] %in% region.ne.num ){
      region.northeast[i] = 1
    }
    if (ID[i] %in% region.mw.num ){
      region.midwest[i] = 1
    }
    if (ID[i] %in% region.w.num ){
      region.west[i] = 1
    }
    if (ID[i] %in% region.s.num ){
      region.south[i] = 1
    }
    
  }
  
  data2$region1 = region.northeast
  data2$region2 = region.midwest
  data2$region3 = region.south
  data2$region4 = region.west
  
  
  # precision1 =  as.matrix(precision_estimation(data1, Omega1))
  # precision2 = as.matrix(precision_estimation(data1, Omega2))
  # precision3 = as.matrix(precision_estimation(data1, Omega3))
  # precision4 = as.matrix(precision_estimation(data1, Omega4))
  ##Use the big precision as a proxy for to stablize the Sd
  precision1 = precision_matrix1
  precision2 = precision_matrix2
  precision3 = precision_matrix3
  
  #################### Work with AR structure
  Omega11 = Omega1
  diag(Omega11)  = 0 #modified to match with the RAND function
  gamma1= 1 #gamma=1: we fully rely on prior network topology H, 
  result1 = GEEfunction_precision(cases~region1 +region2 +region4, id = ID, data = data2, family = gaussian, corstr = "IND", gamma=gamma1, H=precision1)
  
  ######## work with gand
  Omega21 = Omega2
  diag(Omega21)  = 0 #modified to match with the RAND function
  result2 = GEEfunction_precision(cases~region1 +region2 +region4, id = ID, data = data2, family = gaussian, corstr = "IND", gamma=gamma1, H=precision2)
  
  ######## work with Block
  Omega31 = Omega3
  diag(Omega31)  = 0 #modified to match with the RAND function
  result3 =  GEEfunction_precision(cases~region1 +region2 +region4, id = ID, data = data2, family = gaussian, corstr = "IND", gamma=gamma1, H=precision3)
  
  
  sd1mat[b1,] = result1$parameter[,2][2:4]
  sd2mat[b1,] = result2$parameter[,2][2:4]
  sd3mat[b1,] = result3$parameter[,2][2:4]
  
  
}


M1 = apply(sd1mat, 2, mean)
SD1 = apply(sd1mat, 2, sd)
M2 = apply(sd2mat, 2, mean)
SD2 = apply(sd2mat, 2, sd)
M3 = apply(sd3mat, 2, mean)
SD3 = apply(sd3mat, 2, sd)


##Plots Absolute error: Mean and Sd
seqdat= 1:3 #as.numeric(colnames(Y)[191:200])
mydata = data.frame(seqdat,M1,M2, M3)
colnames(mydata) = c("Coefficients", "Isolated", "Band(3)", "Block")

data_long = gather(mydata,  Network, Error, Isolated:Block, factor_key=TRUE)

data_long[,2] <- gsub('Block', 'Block(4)', data_long[,2])
sizetext  = 36
sizetext1 = 30

plt1 = ggplot(data_long, aes(x=Coefficients, y=Error, color=Network)) + 
  geom_line(size=.75,aes(linetype=Network))+ theme_classic()+
  geom_point(aes(shape=Network),size=5)+
  labs( x = "Coefficients", y = "Mean of standard errors")+
  theme(
    axis.text=element_text(size=sizetext),
    axis.title.x = element_text( size=sizetext, face ="bold"),
    axis.title.y = element_text(size=sizetext, face="bold"), legend.position=c(.7, .7),
    legend.title = element_text( size=sizetext, 
                                 face="bold"),legend.text=element_text(size=sizetext1)
  )+
  scale_color_tron()+theme(legend.title=element_blank())+
  scale_fill_tron()+ scale_x_continuous(breaks=c(1,2, 3),labels=c(expression(beta[1]), expression(beta[2]), expression(beta[3])))


plt1
#############

seqdat= 1:3 #as.numeric(colnames(Y)[191:200])
mydata1 = data.frame(seqdat,SD1,SD2, SD3)
colnames(mydata1) = c("Coefficients", "Isolated", "Band(3)",  "Block")

data_long1 = gather(mydata1,  Network, SDError, Isolated:Block, factor_key=TRUE)

data_long1[,2] <- gsub('Block', 'Block(4)', data_long1[,2])


plt2 = ggplot(data_long1, aes(x=Coefficients, y=SDError, color=Network)) + 
  geom_line(size=.75,aes(linetype=Network))+ theme_classic()+
  geom_point(aes(shape=Network),size=5)+
  labs( x = "Coefficients", y = "SD of standard errors")+
  theme(
    axis.text=element_text(size=sizetext),
    axis.title.x = element_text( size=sizetext, face ="bold"),
    axis.title.y = element_text(size=sizetext, face="bold"), legend.position=c(.7, .8),
    legend.title = element_text( size=sizetext, 
                                 face="bold"),legend.text=element_text(size=sizetext1)
  )+ylim(0, 0.5)+
  scale_color_tron()+theme(legend.title=element_blank())+
  scale_fill_tron()+ scale_x_continuous(breaks=c(1,2, 3),labels=c(expression(beta[1]), expression(beta[2]), expression(beta[3])))


plt2

# 
tmp1$parameter
tmp2$parameter
tmp3$parameter

c(test1,test2,test3)

###P-values
#####Check point ######

# 
####Q1. Test statistic and p values
quantile = .95
cutoff = qgev(quantile,loc = - log(2*pi), scale =2, shape =0) + 4*log(p) - log(log(p))
cutoff


p1 = pgev(test1, loc=- log(2*pi)+ 4*log(p) - log(log(p)), scale= 2, shape=0, lower.tail = FALSE)

p2 = pgev(test2, loc=- log(2*pi)+ 4*log(p) - log(log(p)), scale= 2, shape=0, lower.tail = FALSE)

p3 = pgev(test3, loc=- log(2*pi)+4*log(p) - log(log(p)), scale= 2, shape=0, lower.tail = FALSE)


rbind(c(test1,p1), c(test2,p2), c(test3,p3) )
# 