select.model2 <-
function(pp, sigmas, r,  nlarge=10000, q=1/4, p=2, correction="iso"){

   #how many sigmas to evaluate?
   len.sig<- length(sigmas)
   # point density
   pprho<- pp$n/area.owin(pp$window) # densidad de puntos

   
   #define initial parameters for optimization
      #maximum and minimum realistic bvalues for the parameters to be fitted
      # mínimum sigma = 1/10 of the r step;   minimum rho = 1 cluster / study area
    lower=c(sigma2=r[2]/10, rho=1/area.owin(pp$window)) 
    # maximum sigma = 4 times the largest r value (usually this would equal to the shorter side of the plot)
    # maximum rho = actual density (i.e., each point would ba a cluster)
    upper=c(sigma2= (max(r)*4 )^2 , rho=pprho) 
      # initial values for each parameter: the mean one
    sigma2.0<- (upper["sigma2"]-lower["sigma2"])/2
    rho.0 <- (upper["rho"]-lower["rho"])/2
    #parscale<-c(max(upper),min(lower))
    parscale<-c(1,1)

   # fit inhomogeneous Poisson cluster models and compute inhomogenepous K funcitions
   HPPs = list()
   models <- list() 
   for ( i in 1:len.sig){
      progressreport(i, len.sig)
      # estimate intensity with each supplied sigma
      lambda<- density.ppp(pp, sigma=sigmas[i], at="points") 
      # fit inhomogeneous PC model
      hpc.model <- ipc.estK2(pp,lambda=lambda, correction=correction, r=r,
                           nlarge=nlarge, p=p, q=q, sigma2=sigma2.0, rho=rho.0,
                           method= "L-BFGS-B", lower=lower, upper=upper,
   control=list(parscale=parscale))
      models[[i]] <- lambda
      models[[i+len.sig]] <- hpc.model
      HPPs[[i]] <- Kinhom(pp, lambda=lambda, r=r, correction=correction, nlarge=nlarge)
       }

   # fit homogeneopus Poisson cluster model
   # pc.model<- ipc.estK2(pp, correction=correction, r=r, nlarge=nlarge,  p=p, q=q, sigma2=sigma2.0, rho=rho.0) # OJO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   pc.model<- ipc.estK2(pp, correction=correction, r=r, nlarge=nlarge,
p=p, q=q, sigma2=sigma2.0, rho=rho.0, method= "L-BFGS-B",
lower=lower, upper=upper,
control=list(parscale=c(max(upper),min(lower)))) 
   
   models[[(2*len.sig)+1]] <- pc.model

    # Compute homogeneous K function
    homo.lam<- pp$n/area.owin(pp$window )
    P <- Kest(pp, r=r, correction=correction, nlarge=nlarge)
    models[[(2*len.sig)+2]] <- homo.lam # por coherencia pero no sirve para nada

 #-------------------------------------------------
    # Define function to compute discrepancy between observed and expected K funncyions
    # This is the D(theta) function of Diggle 2003: XX)
                
        dtheta.fun <- function(Kobs, Kfit) return(sum((Kobs^q - Kfit^q)^p))
             
    # COmpute differences between observed and expected functions
    # and store the computed fK functions
    dtheta <- NULL
    Kas <- list()
    #first for the IPP's
    for ( i in 1:length(sigmas)){
           dtheta[i] <-   dtheta.fun(HPPs[[i]][[3]], pi*HPPs[[i]]$r^2)   #( HPPs[[i]][[3]] is the observed K function for  model i)
        Kas[[i]] <- HPPs[[i]][[3]] 
    }
   #then for the IPCP's
   for ( i in 1:length(sigmas)){
            dtheta[length(sigmas)+i] <-   dtheta.fun(models[[length(sigmas)+i]]$Kobs, models[[length(sigmas)+i]]$Kfit)
         Kas[[length(sigmas)+i]] <- models[[length(sigmas)+i]]$Kobs

   }
   #finally for the HPCP and HPP
   dtheta[2*length(sigmas)+1] <- dtheta.fun(pc.model$Kobs, pc.model$Kfit)
   dtheta[2*length(sigmas)+2] <- dtheta.fun(P[[3]],   pi*P$r^2)
   
   Kas[[2*length(sigmas)+1]] <- pc.model$Kobs
   Kas[[2*length(sigmas)+2]] <- P[[3]]
  
#-------------------------------------------------

# assign names to the stored models and search for the minor discrepancy
# beware that we call "heterogeneous" the "inhomogeneous" modesl and for this HPC = IPC, etc 
   nombres.modelos<- c(paste("HPP_sg_", sigmas), paste("HPC_sg_", sigmas), "PC", "P")
   names(dtheta) <- nombres.modelos
   best.dtheta<- which.min(dtheta)
   names(models) <- nombres.modelos

   # compute how many models have been fitted
   nHPP <- length(sigmas)
   nHPC  <- nHPP
   nP <- 1
   nPC <-1
   
   #vector with the number of parameters fitted for each model
      # P (homogeneous poisson)= 1 (RSS)
      # HPP (Heterogeneous, i.e. inhomogeneous Poisson= 2 (RSS, bw)
      #  PC (Poisson cluster)= 3 (RSS, sigma, rho)
      # HPC(Heterogeneous, i.e., inhomogeneous Poisson cluster)= 4 (RSS, sigma, rho, bw)
  # where "bw" is the sd (sigma)  of the Gaussian kernel selected to the estimate intensity 
  # and RSS is the meassure of discrepancy (dtheta)
   
    npar <- c(rep(2,nHPP),rep(4,nHPC), 3, 1)
    
    # Compute AICc of each model
    aics<- apply(cbind(dtheta,npar), 1,function(x) aic.function(r,x[1],x[2]))

  # return results 
  result<- list(dtheta=dtheta, best.dtheta = dtheta[best.dtheta], best.model=models[[best.dtheta]], models=models, HPPs=HPPs, sigmas=sigmas, aics=aics, Kas=Kas)
  class(result) <- c("selectedmod",class(result))
  return(result)

}
