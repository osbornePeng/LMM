# This function solves Bayesian linear regression by using EM
# X: n-by-m design matrix for random effects
# Z: n-by-p design matrix for fixed effects, equal to zero when non
# and the first column of Z is intercept
# trainset_y: n-by-1 response matrix
# K is m by m cov matrix of random effects
# 2018/05 @HKUST

#  y = Zu + Xw + e; w~N(0,k*alpha^-1), e~N(0, beta^-1)

EM_for_LMM <- function(trainset_z, trainset_x, trainset_y, k) {
  
  t1 <- proc.time()
  
  # Initialize for simulation

# trainset_x <- matrix(rnorm(300*800),300,800)
# trainset_z <- matrix(rnorm(300*5),300,5)
# trainset_z <- cbind(rep(1,300),trainset_z)
# 
# m <- length(trainset_x[1,]) # The dimension of random efects
# n <- length(trainset_x[,1]) # Number of samples
# m_f <- length(trainset_z[1,]) # The dimension of fixed efects
# 
# w <- matrix(rnorm(m*1),m,1)
# u_true <- matrix(rnorm(m_f*1),m_f,1)
# e <- matrix(rnorm(n*1),n,1)
# k=diag(rep(1,m))
# 
# trainset_y=trainset_z%*%u_true+trainset_x%*%w+e
  
  # Initialization for parameters
  
  trainset_x <- as.matrix(trainset_x)
  trainset_z <- as.matrix(trainset_z)
  trainset_y <- as.matrix(trainset_y)
  
  trainset_x <- scale(trainset_x,center=T,scale=F)
  trainset_y <- scale(trainset_y,center=T,scale=F)
  trainset_z <- scale(trainset_z,center=T,scale=F)
  
  alpha <- 1
  beta <- c(2/var(trainset_y))
  
  maxite <- 200
  low_b <- 0 # The lower bound
  
  # Dimensions
  
  m <- length(trainset_x[1,]) # The dimension of random effects
  n <- length(trainset_y) # Number of samples
  
  # Centering
  
  # trainset_x <- trainset_x-
  #   matrix(t(replicate(n,colMeans(trainset_x))),n,m)
  # trainset_y <- trainset_y-replicate(1,colMeans(trainset_y))

  # Speed up matrix inversion of two variance components
  
  if (all(k %*% rep(1,m) == rep(1,m))){ # k is identity
    inv_k <- k
    vectors <- k
    values <- diag(k)
    eig_k <- list("vectors" = vectors, "values" = values)
  }
  else{
    eig_k <- eigen(k)
    inv_k <- (eig_k$vectors * t(replicate(m, (1 / eig_k$values)))) %*% t(eig_k$vectors) # m*m
  }
  v_1 <- t(trainset_x) %*% trainset_x # x'*x m*m
  v_2_sq <- (eig_k$vectors * t(replicate(m,(sqrt(eig_k$values))))) #v_2^(-0.5)=u_v*1/sqrt(lambda)
  qdqt <- eigen( t(v_2_sq) %*% v_1 %*% v_2_sq ) # qdq^t is the eigen decomposition of...
  U <- v_2_sq %*% qdqt$vectors
  D <- (qdqt$values) # m*1
  
  if (length(trainset_z)>1) { # If the fixed effects exist
    m_f <- length(trainset_z[1,]) # The dimension of fixed effects
    # trainset_z <- (trainset_z-
    #               matrix(t(replicate(n,colMeans(trainset_z))),n,m_f)) # Centering
    inv_z <- solve(t(trainset_z)%*%trainset_z)
    u <- rep(1,m_f)
    if (length(which(trainset_z[,1] == 1))==n) {
      u[1] <- mean(trainset_y) # If there is an intercept
    }
  } 
  else {
    inv_z <- 0
    u <- 0
    trainset_z <- as.matrix(rep(0,n))
  }
  
  # Update
  
  for (i in 1 : maxite) {
    # Posterior moments of w
    
    #a <- beta*v_1+inv_k*alpha # The posterior precision
    s <- (U * t(replicate(m, 1 / (beta * D + alpha)))) %*% t(U)
    mu <- beta * s %*% ( t(trainset_x) %*% (trainset_y-trainset_z %*% u) )
    
    beta_old <- beta
    alpha_old <- alpha
    
    # M step
    
    alpha_part_1 <- t(mu) %*% inv_k %*% mu
    alpha_part_2 <- sum(inv_k * s)
    alpha <- m / (alpha_part_2 + alpha_part_1)
    alpha <- c(alpha)
    
        if (length(trainset_z)>1){ # trainset_z is existing 
      u <- inv_z %*% ( t(trainset_z) %*% (trainset_y - trainset_x%*%mu) )
      u <- c(u)
    } else {
      u <- 0
    }
    
    beta_part <- sum((trainset_y - trainset_z %*% u - trainset_x %*% mu)^2)
    beta <- n / ( sum( (D + 1) / (beta_old * D + alpha) ) - alpha_part_2 + beta_part )
    #beta <- n/(sum(v_1*s)+beta_part)
    beta <- c(beta)
    
    # Lower bound
    
    low_b[i] <- 0.5 * n * log(beta) + 0.5 * m * log(alpha)-
      0.5*(sum( (beta*D + alpha) / (beta_old * D + alpha_old) ) +
             alpha_part_1 * alpha + beta_part * beta) +
      0.5 * (sum(-log( (beta_old*D + alpha_old) )))
    
    # print(i)
    # Break
    
    if (i>1) {
      if (abs((low_b[i] - low_b[i-1]) / low_b[i - 1]) < 1e-6) {
        break
      }
    }
  }
  
  t2 <- proc.time()
  total_time <- t2-t1 # Total time
  
  result <- list(alpha = alpha, beta = beta, mu = mu, u = u,
                 low_b = low_b, total_time = total_time)
  return(result)
}