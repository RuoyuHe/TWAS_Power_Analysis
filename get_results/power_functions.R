get_power <- function(theta, vartheta, alpha, m){
  # N: n1 + n3
  # vartheta: the vartheta given by power_fn
  # alpha: significance level
  # m: number of tests
  cutoff = alpha/m
  'upper = pnorm(qnorm(cutoff/2, sd = sqrt(vartheta), lower.tail = F), 
                mean = theta, sd = sqrt(vartheta), lower.tail = F)
  lower = pnorm(qnorm(cutoff/2, sd = sqrt(vartheta)), 
                mean = theta, sd = sqrt(vartheta))'
  upper = pnorm(qnorm(cutoff/2, lower.tail = F) - theta/sqrt(vartheta), lower.tail = F)
  lower = pnorm(qnorm(cutoff/2) - theta/sqrt(vartheta))
  return(upper+lower)
}

power_fn <- function(n1, n2, covX, varY, theta, vartheta){
  # n1, n2: sample sizes of SNP (Z) and GWAS (Y), respectively
  # covX: covariance matrix of the error term in stage 1
  # varY: variance of the error term in stage 2
  # theta: estimated effects of the stage 2 model
  # vartheta: covariance matrix of the effects in stage 2
  
  inflation = (1+(n2/n1 * (t(theta)%*%covX%*%theta) / varY))[1]
  vartheta_corrected = vartheta*inflation
  
  return(list(vartheta = vartheta_corrected))
}

get_chisq_power <- function(theta, vartheta, alpha, m){
  # N: n1 + n3
  # vartheta: the vartheta given by power_fn
  # alpha: significance level
  # m: number of tests
  sig_level = alpha/m
  cutoff = qchisq(sig_level, length(theta), lower.tail = FALSE)
  ncp = t(theta)%*%solve(vartheta)%*%theta
  power = pchisq(cutoff, length(theta), ncp = ncp, lower.tail = FALSE)
  return(power)
}
