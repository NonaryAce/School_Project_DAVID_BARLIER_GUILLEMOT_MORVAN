oxford = function(nchain, data, init, prop.sd){
  
  ## Initialisation grace aux donnees
  N = dim(data)[1]
  r1 = data$r1
  n1 = data$n1
  r0 = data$r0
  n0 = data$n0
  year = data$year
  
  
  
  
  ## Calculs une bonne fois pour toute
  dim = length(init)
  
  
  
  
  ## Creation de la chaine et du taux d'acceptation
  chain = matrix(NA, nrow=nchain+1, ncol=dim) 
  colnames(chain) = c("Alpha","Beta_1","Beta_2","Sigma")
  chain[1,] = init
  acc.rates = rep(0, dim-1)
  
  chain_mu = matrix(NA, nchain+1, N)
  chain_mu[1,] = 0
  acc.rates.mu = 0
  
  chain_b = matrix(NA, nchain+1, N)
  chain_b[1,] = 0
  acc.rates.b = 0
  
  
  
  
  ## Mise a jour des parametres
  for (iter in 1:nchain){
    
    ## Current
    mu = chain_mu[iter,]
    b = chain_b[iter,]
    current = chain[iter, ]
    
    
    ## Mise a jour des mu avec M-H
    for (j in 1:120){
      
      # Noyau de proposition
      prop = rnorm(1,mu[j],prop.sd[4])
      
      # Calculs de p0
      p0_prop = 1/(1+exp(-prop))
      p0_current = 1/(1+exp(-mu[j]))
      
      # A l'echelle log
      top = ((-1/2)*(10^-(6))*prop^2) + r0[j]*log(p0_prop) + (n0[j]-r0[j])*log(1-p0_prop)
      bottom = ((-1/2)*(10^-(6))*mu[j]^2) + r0[j]*log(p0_current) + (n0[j]-r0[j])*log(1-p0_current)
      
      # Calcul de la probabilite d'acceptation
      acc.prob = exp(top - bottom)
      
      # Acceptation ou non de la proposition
      if (runif(1) < acc.prob){
        mu[j] = prop
        acc.rates.mu = acc.rates.mu + 1
      }
      
      chain_mu[iter+1,j] = mu[j] 
      
    }
    
    
    mu = chain_mu[iter,]
    
    # Mise a jour des b avec M-H
    for (k in 1:120){
      
      # Noyau de proposition
      prop = rnorm(1,b[k],prop.sd[5])
      
      # Calcul de p1
      p1_prop = 1/(1+exp(-mu[k] - current["Alpha"] - current["Beta_1"]*year[k] - 
                           current["Beta_2"]*(year[k]^2 - 22) - prop))
      p1_current = 1/(1+exp(-mu[k] - current["Alpha"] - current["Beta_1"]*year[k] -
                              current["Beta_2"]*(year[k]^2 - 22) - b[k]))
      
      # A l'echelle log
      top = ((-1/2)*(1/current["Sigma"]^2)*prop^2) + r1[k]*log(p1_prop) +
        (n1[k]-r1[k])*log(1-p1_prop)
      bottom = ((-1/2)*(1/current["Sigma"]^2)*b[k]^2) + r1[k]*log(p1_current) +
        (n1[k]-r1[k])*log(1-p1_current)
      
      # Calcul de la probabilite d'acceptation
      acc.prob = exp(top - bottom)
      
      # Acceptation ou non de la proposition
      if (runif(1) < acc.prob){
        b[k] = prop
        acc.rates.b = acc.rates.b + 1
      }
      
      chain_b[iter+1,k] = b[k]
      
    }
    
    
    b = chain_b[iter,]
    
    
    ## Mise a jour de sigma
    update_shape = 10^(-3) + N/2
    update_rate = 10^(-3) + sum(b^2)/2
    sigma = 1/sqrt(rgamma(1, shape=update_shape, rate=update_rate))
    
    chain[iter+1, "Sigma"] = (sigma)
    
    

    
    ## Mise a jour des alpha, beta1 et beta2 avec M-H
    K_P = current[1:3]
    mu = chain_mu[iter,]
    
    for (j in 1:3){
      
      
      # Noyau de proposition
      prop = K_P
      prop[j] = rnorm(1, current[j], prop.sd[j])
      
      # Calcul de p1,i
      
      p1_prop = 1/(1+exp(-mu - prop[1] - prop[2]*year - prop[3]*(year*year - 22) - b))
      p1_current = 1/(1+exp(-mu - current["Alpha"] - current["Beta_1"]*year -
                              current["Beta_2"]*(year^2 - 22) - b))
      
      # A l'echelle log
      top = (-1/2)*(10^(-6))*prop[j]^2 + sum(r1*log(p1_prop) + (n1-r1)*log((1-p1_prop)))
      bottom = (-1/2)*(10^(-6))*current[j]^2 + sum(r1*log(p1_current) + (n1-r1)*log((1-p1_current)))
      
      # Calcul de la probabilite d'acceptation
      acc.prob = exp(top - bottom)
      
      # Acceptation ou non de la proposition
      if (runif(1) < acc.prob){
        current[j] = prop[j]
        acc.rates[j] = acc.rates[j] + 1
      }
      
      chain[iter+1,j] = current[j]
      
    }
    
  }
  
  return(list(chain=chain, chain_mu=chain_mu, chain_b=chain_b, acc.rates=acc.rates/nchain,
              acc.rates.mu=acc.rates.mu/(nchain*N), acc.rates.b=acc.rates.b/(nchain*N)))
}



# Valeurs initiales
init = matrix(c(0,0,0,1), ncol=4)
colnames(init) = c("Alpha","Beta_1","Beta_2","Sigma")

# Choix de prop.sd
prop.sd = matrix(c(0.09,0.014,0.003,1.2,1.2), ncol=5)
colnames(prop.sd) = c("Alpha","Beta_1","Beta_2","Mu","b")




out = oxford(10^4, data, init=init, prop.sd=prop.sd)
out$acc.rates
out$acc.rates.mu
out$acc.rates.b




## Pour b
plot(out$chain_b[,12], type='l')
plot(density(out$chain_b[,12]), type='l')




## Pour mu
plot(out$chain_mu[,12], type='l')
plot(density(out$chain_mu[,12]), type='l')




## Pour alpha
alpha_burn = out$chain[,"Alpha"][1000:10000]

plot(alpha_burn, type='l')
plot(density(alpha_burn), type='l')
mean(alpha_burn)
sd(alpha_burn)




## Pour beta_1
beta_burn =  out$chain[,"Beta_1"][1000:10000]

plot(beta_burn, type='l')
plot(density(beta_burn), type='l')
mean(beta_burn)
sd(beta_burn)




## Pour beta_2
beta_burn =  out$chain[,"Beta_2"][1000:10000]


plot(beta_burn, type='l')
plot(density(beta_burn), type='l')
mean(beta_burn)
sd(beta_burn)




## Pour sigma
plot(out$chain[,"Sigma"][100:1000], type='l')
plot(density(out$chain[,"Sigma"][1000:10000]), type='l')
mean(out$chain[,"Sigma"][1000:10000])
sd(out$chain[,"Sigma"][1000:10000])




## Odd ratio
odd_ratio = exp(mean(out$chain[,"Alpha"]) + mean(out$chain[,"Beta_1"])*year + 
                  mean(out$chain[,"Beta_2"])*(year^2 - 22) + apply(out$chain_b,c(2),mean))
summary(odd_ratio[2:56])
summary(odd_ratio[56:120])
