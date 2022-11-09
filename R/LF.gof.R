
LF.gof<- function(X, rmin=NULL, rmax=NULL, na.rm=TRUE){
    
   # controls that is using an envelope object and that it has been computed with the argument savefuns=TRUE
    if(sum(!is.na(match(class(X),"envelope")))==0) stop("LF.gof only accept spatstat envelope objects")
    if(is.null(attributes(X)$simfuns)) stop ("did you compute the envelopes with the argument 'savefuns=TRUE'?") 
    
    # get data
    r <- as.data.frame(attributes(X)$simfuns)[,1]
    H1 <- X$obs
    Hi <- as.data.frame(attributes(X)$simfuns)[,-1]
    
    #Join observed and simulated functions
    H <- cbind(H1,Hi) 
    
    
    # select the range of integrable results
    if (!is.null(rmin) & !is.null(rmax)){ H <- H[r >= rmin & r <= rmax,]}
    if (!is.null(rmin) & is.null(rmax)) { H <- H[r >= rmin,]}
    if (is.null(rmin) & !is.null(rmax)) { H <- H[r <= rmax,]}
    
    # count how may NA values for each distance r
    na.count<-apply(H,1, function(x) sum(is.na(x)))
    na.count.by.r <-na.count[na.count>0] 

    s <- dim(H)[2]

    # compute summary statistics for the observed (i.e., first) pattern 
    Hmeans <- rowSums (H[,-1],  na.rm = na.rm)/(s-1)
    
    u <- sum((H[,1]-Hmeans)^2, na.rm=T)
    
    #width of the distance interval. the second value of the r vector (the first one is 0)
    delta_tk <- X$r[2]

     
     # compute summary statistics for the simulated patterns
     for ( i in 2:s) {
         Himeans <- rowSums (H[,-i],na.rm=na.rm)/(s-1) 
         u <- c(u, sum((H[,i]-Himeans)^2, na.rm=na.rm)) 
       }
     p <-  1-((rank(u)[1]-1)/(s))
     
return(list(u=u[1]*delta_tk, p=p,na.count.by.r=na.count.by.r))

# computatiopn according to Baddeley et al. 2014  https://doi.org/10.1890/13-2042.1
# dr <- X$r[2]
# H0 <- X$obs
# simfuns <- data.frame(attributes(X)$simfuns)[,-1]
# m <- ncol(simfuns)
# a <- ((m+1)/m)^2

# H1bar <- apply(simfuns, 1, mean)
# H2bar <- (m*H1bar /(m+1)) +(H0/(m+1))

# U <- a * sum(dr*(H0-H2bar)^2)
# for ( i in 1:m) {
#           U<- c(U, a*sum(dr*(simfuns[,i]-H2bar)^2, na.rm=T)) 
#       }
#     p <-  1-((rank(U)[1]-1)/(m+1))
# return( u = U[1], p = p)

}
