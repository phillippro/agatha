###The following function is adapted from bgls.py by Mortier et al. 2015
bgls <- function(t, y, err, ofac=1,fmax=1,fmin=NA,tspan=NULL){
    t <- t-min(t)
    dy <- err
    if(is.null(tspan)){
        tspan = max(t)-min(t)
    }
    if(tspan<100){
       ofac <- 10
    }
    step <- 1/(tspan*ofac)
#    fmax <- 1
    if(is.na(fmin)){
        fmin <- 1/(ofac*tspan)
    }
    f = seq(fmin,fmax,by=step)
    nout <- length(f)
    omegas <- 2*pi*f
    err2 <- err * err
    w = 1./err2
    W <- sum(w)
    bigY <- sum(w*y)  # Eq. (10)
    p <- c()
    constants <- rep(NA,length(omegas))
    exponents <- rep(NA,length(omegas))
    for(k in 1:length(omegas)){
        omega <- omegas[k]
        theta = 0.5*atan(sum(w*sin(2*omega*t))/sum(w*cos(2*omega*t)))
        x = omega*t - theta
        cosx = cos(x)
        sinx = sin(x)	 
        wcosx = w*cosx
        wsinx = w*sinx

        C = sum(wcosx)
        S = sum(wsinx)
        YCh = sum(y*wcosx)
        YSh = sum(y*wsinx)
        CCh = sum(wcosx*cosx)
        SSh = sum(wsinx*sinx)
        if(CCh != 0 & SSh != 0){
            K = (C*C*SSh + S*S*CCh - W*CCh*SSh)/(2.*CCh*SSh)
            L = (bigY*CCh*SSh - C*YCh*SSh - S*YSh*CCh)/(CCh*SSh)
            M = (YCh*YCh*SSh + YSh*YSh*CCh)/(2.*CCh*SSh)
            constants[k] <- 1./sqrt(CCh*SSh*abs(K))
        }else if(CCh == 0){
            K = (S*S - W*SSh)/(2.*SSh)
            L = (bigY*SSh - S*YSh)/(SSh)
            M = (YSh*YSh)/(2.*SSh)
            constants[k] <- 1./sqrt(SSh*abs(K))
        }else if(SSh == 0){
            K = (C*C - W*CCh)/(2.*CCh)
            L = (bigY*CCh - C*YCh)/(CCh)
            M = (YCh*YCh)/(2.*CCh)
            constants[k] <- 1./sqrt(CCh*abs(K))
        }
        if(K > 0){
            cat('K is positive. This should not happen.')
        }
        exponents[k] <- M - L^2/(4*K)
    }
    power  <- log(constants) + exponents
    inds <- sort(power,decreasing=TRUE,index.return=TRUE)$ix
    ps <- 1/f[inds[1:5]]
    power.opt <- power[inds[1:5]]
    ind <- which.max(power)
    omega.opt <- 2*pi*f[ind]
####optimized parameter
    Amin <- Bmin <- -2*(max(y)-min(y))
    Amax <- Bmax <- 2*(max(y)-min(y))
    Aini <- Bini <- (Amin+Amax)/2
    gamma.min <- min(y)
    gamma.max <- max(y)
    gamma.ini <- (gamma.min+gamma.max)/2
    phi <- 0
    start <- list(A=Aini,B=Bini,gamma=gamma.ini)
    par.low <- c(A=Amin,B=Bmin,gamma=gamma.min)
    par.up <- c(A=Amax,B=Bmax,gamma=gamma.max)
    par.fix <- list(omega=omega.opt,phi=phi)
    df <- list(t=t,y=y,dy=dy,par.fix=par.fix)
    out <- nls.lm(par = start,lower=par.low,upper=par.up,fn =gls.res,df=df,control=nls.lm.control(maxiter=500))
    par.full <- c(as.list(coef(out)),par.fix)
    yp <- gls.model(t,par.full)
    res <- y-yp
    level <- log(c(10,100,1000))
    return(list(P=1/f, power=power,res=res,ps=ps,power.opt=power.opt,sig.level=level))
}

###based on Zechmeister09.pdf or ZK09 and the lsp function in the 'lomb' library
gls <- function(t, y, err,ofac=1, norm="Cumming",fmax=1,fmin=NA,tspan=NULL){
    t <- t-min(t)
    dy <- err
    if(is.null(tspan)){
        tspan = max(t)-min(t)
    }
    step <- 1/(tspan*ofac)
#    f.max <- 1
    if(is.na(fmin)){
        fmin <- 1/(tspan*ofac)
    }
    f = seq(fmin,fmax,by=step)
    nout <- length(f)
#    nout = ofac * hifac * length(t)/2
#    xdif = max(t)-min(t)
#    f = 1./(xdif*ofac) + c(1:nout)/(ofac*xdif)
    omegas <- 2*pi*f
    err2 <- err * err
    W <- sum(1/err2)
    w <- 1/err2/W
    bigY <- sum(w*y)  # Eq. (7)
    p <- rep(NA,length(omegas))
    for(i in 1:length(omegas)){
        omega <- omegas[i]
        bigC <- sum(w*cos(omega*t))
        bigS <- sum(w*sin(omega*t))
        YY.hat <- sum(w*(y^2))
        YC.hat <- sum(w*y*cos(omega*t))
        YS.hat <- sum(w*y*sin(omega*t))
        CC.hat <- sum(w*cos(omega*t)^2)
        SS.hat <- sum(w*sin(omega*t)^2)
        CS.hat <- sum(w*sin(omega*t)*cos(omega*t))
        YY <- YY.hat-bigY*bigY
        YC <- YC.hat-bigY*bigC
        YS <- YS.hat -bigY*bigS
        CC <- CC.hat - bigC*bigC
        SS <- SS.hat - bigS*bigS
        CS <- CS.hat - bigC*bigS
        bigD <- CC*SS-CS^2
        p[i] <- (SS*YC^2+CC*YS^2-2*CS*YC*YS)/(YY*bigD)
    }
    N <- length(y)
 # An ad-hoc estimate of the number of independent frequencies (ZK_09 Eq. 24)
    #M <- (max(f)-min(f))*(max(t)-min(t))
    M <- 2*nout/ofac#ref lsp
    if(norm=='Scargle'){
        popvar <- 1#arbitrary value or input
        power <- p/popvar
    }
    if(norm=='HorneBaliunas'){
        power <- (N-1)*p/2
    }
    m <- 1
    if(norm=='Cumming'){
        power <- ((N-2-m)/2)*p/(1-max(p))
    }
    PN <- power
    PN.max <- max(PN)
    peak.freq <- f[PN==PN.max]
    peak.per <- 1/f[PN==PN.max]
    FAP <- c(0.1,1e-2,1e-3)#significance level of FAP
    level <- powerLevel(FAP,M,N,norm=norm,m=m)#power level
    pp <- M*prob(Pn=PN.max,N=N,norm=norm,m=m)
    if(pp>0.01){
        pp <- 1-(1-prob(Pn=PN.max,N=N,norm=norm,m=m))^M
    }
    ind <- which.max(pp)
    omega.opt <- omegas[ind]
#####optimal parameters
    Amin <- Bmin <- -2*(max(y)-min(y))
    Amax <- Bmax <- 2*(max(y)-min(y))
    Aini <- Bini <- (Amin+Amax)/2
    gamma.min <- min(y)
    gamma.max <- max(y)
    gamma.ini <- (gamma.min+gamma.max)/2
    phi <- 0
    start <- list(A=Aini,B=Bini,gamma=gamma.ini)
    par.low <- c(A=Amin,B=Bmin,gamma=gamma.min)
    par.up <- c(A=Amax,B=Bmax,gamma=gamma.max)
    par.fix <- list(omega=omega.opt,phi=phi)
    df <- list(t=t,y=y,dy=dy,par.fix=par.fix)
    out <- nls.lm(par = start,lower=par.low,upper=par.up,fn =gls.res,df=df,control=nls.lm.control(maxiter=500))
    par.full <- c(as.list(coef(out)),par.fix)
    yp <- gls.model(t,par.full)
    res <- y-yp
#    cat('level:',level,'\n')
    inds <- sort(power,decreasing=TRUE,index.return=TRUE)$ix
    ps <- 1/f[inds[1:5]]
    power.opt <- power[inds[1:5]]
    return(list(P=1/f, power=power, pvalue=pp, sig.level=level,res=res,par=par.full,ps=ps,power.opt=power.opt))
}

###generalized lomb-scargle periodogram with trend component
glst <- function(t, y, err,ofac=1, norm="Cumming",fmax=1,fmin=NA,tspan=NULL){
    unit <- 365.24#to make the elements of the matrix in the function of 'solve' on the same order
    t <- t-min(t)
    if(is.null(tspan)){
        tspan <- max(t)-min(t)
    }
    step <- 1/(tspan*ofac)
    if(is.na(fmin)){
        fmin <- 1/(tspan*ofac)
    }
    f <- seq(fmin,fmax,by=step)*unit
    nout <- length(f)
    t <- (t-min(t))/unit
    omegas <- 2*pi*f
    dy <- err
    err2 <- err * err
    W <- sum(1./err2)
    w <- 1/err2/W
    bigY <- sum(w*y)  # Eq. (7)
    p <- rep(NA,length(omegas))
    bigT <- sum(w*t)
    YY.hat <- sum(w*y^2)
    YT.hat <- sum(w*y*t)
    TT.hat <- sum(w*t^2)
    YY <- YY.hat - bigY*bigY
    YT <- YT.hat - bigY*bigT
    TT <- TT.hat - bigT*bigT
###optimized parameterse for the trend model
    d0 <- YT/TT
    c0 <- bigY-d0*bigT
    chi2.ref <- sum(W*w*(y-(c0+d0*t))^2)
    for(k in 1:length(f)){
        omega <- omegas[k]
        bigC <- sum(w*cos(omega*t))
        bigS <- sum(w*sin(omega*t))
        YY.hat <- sum(w*y^2)
        YC.hat <- sum(w*y*cos(omega*t))
        YS.hat <- sum(w*y*sin(omega*t))
        CC.hat <- sum(w*cos(omega*t)^2)
        SS.hat <- sum(w*sin(omega*t)^2)
        CS.hat <- sum(w*sin(omega*t)*cos(omega*t))
        ST.hat <- sum(w*sin(omega*t)*t)
        CT.hat <- sum(w*cos(omega*t)*t)
        YC <- YC.hat - bigY*bigC
        YS <- YS.hat - bigY*bigS
        CC <- CC.hat - bigC*bigC
        SS <- SS.hat - bigS*bigS
        ST <- ST.hat - bigS*bigT
        CS <- CS.hat - bigC*bigS
        bigD <- CC*SS-CS^2
        lin.mat <- matrix(c(CC.hat,CS.hat,bigC,CT.hat,CS.hat,SS.hat,bigS,ST.hat,bigC,bigS,1,bigT,CT.hat,ST.hat,bigT,TT.hat),byrow=TRUE,nrow=4)
        vec.rh <- c(YC.hat,YS.hat,bigY,YT.hat)
        pp <- solve(lin.mat,vec.rh,tol=1e-16)
        yp <- pp[1]*cos(omega*t)+pp[2]*sin(omega*t)+pp[3]+pp[4]*t
        chi2 <- sum(W*w*(y-yp)^2)
        p[k] <- (chi2.ref-chi2)/chi2.ref
    }
    m <- 2#floating mean:1; floating trend: 2
    N <- length(y)
    M <- 2*nout/ofac#ref lsp
    if(norm=='Scargle'){
        popvar <- 1#arbitrary value or input
        power <- p/popvar
    }
    if(norm=='HorneBaliunas'){
        power <- (N-1)*p/2
    }
    if(norm=='Cumming'){
        power <- ((N-2-m)/2)*p/(1-max(p))
    }
    PN <- power
    PN.max <- max(PN)
    peak.freq <- f[PN==PN.max]
    peak.per <- 1/f[PN==PN.max]
#    FAP <- c(0.317,0.046,3e-3)#significance level of FAP
    FAP <- c(0.1,0.01,0.001)
    level <- powerLevel(FAP,M,N,norm,m=m)#power level
    pp <- M*prob(Pn=PN.max,N=N,norm=norm,m=m)
    if(pp>0.01){
        pp <- 1-(1-prob(Pn=PN.max,N=N,norm=norm,m=m))^M
    }
    P <- unit/f
    ind.max <- which.max(power)
    omega.opt <- 2*pi*f[ind.max]
#####optimal parameters
    Amin <- Bmin <- -2*(max(y)-min(y))
    Amax <- Bmax <- 2*(max(y)-min(y))
    Aini <- Bini <- (Amin+Amax)/2
    gamma.min <- min(y)
    gamma.max <- max(y)
    gamma.ini <- (gamma.min+gamma.max)/2
    beta.min <- -(max(y)-min(y))/(max(t)-min(t))
    beta.max <- (max(y)-min(y))/(max(t)-min(t))
    beta.ini <- (beta.min+beta.max)/2
    phi <- 0
    start <- list(A=Aini,B=Bini,gamma=gamma.ini,beta=beta.ini)
    par.low <- c(A=Amin,B=Bmin,gamma=gamma.min,beta=beta.min)
    par.up <- c(A=Amax,B=Bmax,gamma=gamma.max,beta=beta.max)
    par.fix <- list(omega=omega.opt,phi=phi)
    df <- list(t=t,y=y,dy=dy,par.fix=par.fix)
    out <- nls.lm(par = start,lower=par.low,upper=par.up,fn =gls.res,df=df,control=nls.lm.control(maxiter=500))
    par.full <- c(as.list(coef(out)),par.fix)
    yp <- glst.model(t,par.full)
    res <- y-yp
    inds <- sort(power,decreasing=TRUE,index.return=TRUE)$ix
    ps <- 1/f[inds[1:5]]
    power.opt <- power[inds[1:5]]
    return(list(P=unit/f, power=power, pvalue=pp, sig.level=level,Popt=P[ind.max],ps=ps,res=res,power.opt=power.opt))
}

#give a power, calcuate the the p value
prob <- function(Pn,N,m,norm='Cumming'){
    if(norm=="Scargle") return(exp(-Pn))
    if(norm=="HorneBaliunas") return((1-2*Pn/(N-1))^((N-2-m)/2))
    if(norm=="Cumming") return((1+2*Pn/(N-2-m))^(-(N-2-m)/2))
}
###Inverse of prob
probInv <- function(Prob,N,m,norm='Cumming'){
    if(norm=="Scargle") return(-log(Prob))
    if(norm=="HorneBaliunas") return((N-1)/2*(1-Prob^(2/(N-2-m))))
    if(norm=="Cumming") return((N-2-m)/2*(Prob^(-2/(N-2-m))-1))
}
####
powerLevel <- function(FAPlevel,M,N,m,norm='Cumming'){
    return(probInv(1-(1-FAPlevel)^(1/M),N,norm,m=m))
}

lsp <- function (x, times = NULL, from = NULL, to = NULL, tspan=NULL, ofac = 1, alpha = 0.01) 
{
    times <- as.numeric(times)
    start <- min(times)
    end <- max(times)
    av.int <- mean(diff(times))
    o <- order(times)
    times <- times[o]
    x <- x[o]
    y <- cbind(times, x)
    t <- y[, 1]
    y <- y[, 2]
    n <- length(y)
    if(is.null(tspan))    tspan <- t[n] - t[1]
    step <- 1/(tspan * ofac)
    fmax <- to
    fmin <- from
    if(is.null(from)){
        fmin <- 1/(ofac*tspan)
    }
    freq = seq(fmin,fmax,by=step)
    n.out <- length(freq)
    if (n.out == 0) 
        stop("erroneous frequency range specified ")
    x <- t * 2 * pi
    if(sd(y)!=0){
        y <- y - mean(y)
    }
    if(var(y)==0){
        norm <- 1/2
    }else{
        norm <- 1/(2 * var(y))
    }
    w <- 2 * pi * freq
    PN <- rep(0, n.out)
    for (i in 1:n.out){
        wi <- w[i]
        tau <- 0.5 * atan2(sum(sin(wi * t)), sum(cos(wi * t)))/wi
        arg <- wi * (t - tau)
        cs <- cos(arg)
        sn <- sin(arg)
        A <- (sum(y * cs))^2
        B <- sum(cs * cs)
        C <- (sum(y * sn))^2
        D <- sum(sn * sn)
        PN[i] <- A/B + C/D
    }
    PN <- norm * PN
    PN.max <- max(PN)
    peak.freq <- freq[PN == PN.max][1]
    peak.at <- c(peak.freq, 1/peak.freq)
    effm <- 2 * n.out/ofac
    level <- -log(1 - (1 - alpha)^(1/effm))
    exPN <- exp(-PN.max)
    p <- effm * exPN
    if (p > 0.01) p <- 1 - (1 - exPN)^effm
####optimized parameter
    Amin <- Bmin <- -2*(max(y)-min(y))
    Amax <- Bmax <- 2*(max(y)-min(y))
    Aini <- Bini <- (Amin+Amax)/2
    omega <- 2*pi*peak.freq
    phi <- 0
    start <- list(A=Aini,B=Bini)
    par.low <- c(A=Amin,B=Bmin)
    par.up <- c(A=Amax,B=Bmax)
    par.fix <- list(omega=omega,phi=phi)
    df <- list(t=t,y=y,par.fix=par.fix)
    out <- nls.lm(par = start,lower=par.low,upper=par.up,fn =lsp.res,df=df,control=nls.lm.control(maxiter=500))
    par.full <- c(as.list(coef(out)),par.fix)
    yp <- lsp.model(t,par.full)
    res <- y-yp
####sort
    Ps <- 1/freq
    inds <- sort(Ps,index.return=TRUE,decreasing=TRUE)$ix
    power <- PN[inds]
    P <- Ps[inds]
###peaks
    inds <- sort(power,decreasing=TRUE,index.return=TRUE)$ix
    ps <- 1/freq[inds[1:5]]
    power.opt <- power[inds[1:5]]
###
    sp.out <- list(P = P, power = power, ps=ps, power.opt=power.opt, alpha = alpha, sig.level = level, peak = PN.max, peak.at = peak.at, p.value = p,par=par.full,res=res,power.opt=power.opt)
    return(sp.out)
}

lsp.model <- function(t,par){
    A <- par$A 
    B <- par$B
    omega <- par$omega
    phi <- par$phi
    A*cos(omega*t-phi)+B*sin(omega*t-phi)
}
lsp.res <- function(par,df){
    t <- df$t
    y <- df$y
    A <- par$A 
    B <- par$B
    phi <- df$par.fix$phi
    omega <- df$par.fix$omega
    y-(A*cos(omega*t-phi)+B*sin(omega*t-phi))
}

gls.model <- function(t,par){
    A <- par$A 
    B <- par$B
    gamma <- par$gamma
    omega <- par$omega
    phi <- par$phi
    A*cos(omega*t-phi)+B*sin(omega*t-phi)+gamma
}
gls.res <- function(par,df){
    t <- df$t
    y <- df$y
    dy <- df$dy
    A <- par$A 
    B <- par$B
    gamma <- par$gamma
    phi <- df$par.fix$phi
    omega <- df$par.fix$omega
    (y-(A*cos(omega*t-phi)+B*sin(omega*t-phi))+gamma)/abs(dy)
}

glst.model <- function(t,par){
    A <- par$A 
    B <- par$B
    gamma <- par$gamma
    beta <- par$beta
    omega <- par$omega
    phi <- par$phi
    A*cos(omega*t-phi)+B*sin(omega*t-phi)+gamma
}
glst.res <- function(par,df){
    t <- df$t
    y <- df$y
    dy <- df$dy
    A <- par$A 
    B <- par$B
    gamma <- par$gamma
    beta <- par$beta
    phi <- df$par.fix$phi
    omega <- df$par.fix$omega
    (y-(A*cos(omega*t-phi)+B*sin(omega*t-phi))+gamma+beta*t)/abs(dy)
}
