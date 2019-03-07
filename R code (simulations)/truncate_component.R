truncate_component <- function(Y, theta=NA, mu=0, sigma=1)
{
# included the regressors in delta1 and delta2

  	delta1 <- (theta[Y] - mu) / sigma
  	delta2 <- (theta[Y+1] - mu) / sigma
  	tmp1 <- (dnorm(delta1) - dnorm(delta2)) / (pnorm(delta2) - pnorm(delta1))
  	EX <- mu + tmp1 * sigma

  	
  	delta1[delta1 < -1e+10] <- -1e+10
  	delta2[delta2 > 1e+10] <- 1e+10
  	tmp2 <- (delta1*dnorm(delta1) - delta2*dnorm(delta2)) / (pnorm(delta2) - pnorm(delta1))
  	EXX <- sigma^2 + mu^2 + sigma^2 * tmp2 + 2*mu*sigma*tmp1

	return(list(EX=EX, EXX=EXX))
}

