library(minpack.lm)
source('periodograms.R')
#library(lomb)
tol <- 1e-16
RV.model <- function(par,data){
####
    t <- data$data[,1]
    y <- data$data[,2]
    dy <- data$data[,3]
    Indices <- data$Indices
    Nma <- data$Nma
    NI <- data$NI
    if(!any('A'==names(par))){
        A <- B <- phi <- omega <- 0
    }else{
        phi <- data$par.fix$phi
        omega <- data$par.fix$omega
    }
    var <- names(par)
    d <- c()
    m <- c()
    for(k in 1:length(var)){
        assign(var[k],par[[var[k]]])
        if(grepl('^d',var[k])) d <- c(d,par[[var[k]]])
        if(grepl('^m',var[k])) m <- c(m,par[[var[k]]])
    }
    ind <- sort(t,index.return=TRUE)$ix
    if(NI>0){
        r <- A*cos(omega*t-phi)+B*sin(omega*t-phi)+gamma+beta*t+partI(d,Indices,NI)
    }else{
        r <- A*cos(omega*t-phi)+B*sin(omega*t-phi)+gamma+beta*t
    }
    val <- r
    if(Nma>0){
        for(k in 1:Nma){
            t2 <- c(rep(0,k),t[-(1:k)])
            t1 <- c(rep(0,k),t[-(length(t)+1-(1:k))])
            ys <- c(rep(0,k),y[-(length(t)+1-(1:k))])
            rs <- c(rep(0,k),r[-(length(t)+1-(1:k))])
            val <- val+m[k]*exp(-abs(t2-t1)/exp(logtau))*(ys-rs)
        }
    }
    return(val)
}
partI <- function(d,Is,NI){
    if(NI==1){
        d*Is[,1]
    }else{
        d%*%t(Is[,1:NI])
    }
}
wI <- function(w,Is,NI){
    unlist(lapply(1:NI,function(i) sum(w*Is[,i])))
}
me <- function(m,t,logtau,x){
    if(length(m)==1){
        m*exp(-abs(t[-1]-t[-length(t)])/exp(logtau))*x[-length(x)]
    }else{
#        unlist(lapply(1:length(m), function(i) m[i]*exp(-abs(t[-c(1:i)]-t[-(length(t)+1-c(1:i))])/exp(logtau))*x[-(length(x)+1-c(1:i))]))
    }
}
rv.res <- function(par,data){
    y <- data$data[,2]
    dy <- data$data[,3]
###add a large constant to make the logL always positive and thus avoid errors
#    tmp <- sqrt((y-RV.model(par,data))^2/(2*(dy^2+par$sj^2)) + 0.5*log(dy^2+par$sj^2)+100)
    tmp <- (y-RV.model(par,data))^2/(2*(dy^2+par$sj^2)) + 0.5*log(dy^2+par$sj^2)+100
    if(any(is.na(tmp))){
        cat('inside=',head((y-RV.model(par,data))^2/(2*(dy^2+par$sj^2))),'\n')
    }
    return(tmp)
}

######optimize all parameters using the LM algorithm
nlopt <- function(pars,type='noise'){
    var <- names(pars)
    for(k in 1:length(var)){
        assign(var[k],pars[[var[k]]])
    }
    par.opt <- NULL
    chi2 <- NULL
    if(type=='period'){
        par.fix <- list(omega=omega,phi=phi)
    }else{
        par.fix <- NULL
    }
    tmp <- data.frame(t,y,dy)
    df <- list(data=tmp,Indices=Indices,par.fix=par.fix,Nma=Nma,NI=NI)
    trace <- FALSE
    if(type!='noise'){
        start <- list(A=Aini,B=Bini,gamma=gamma.ini,beta=beta.ini)
        par.low <- c(A=Amin,B=Bmin,gamma=gamma.min,beta=beta.min)
        par.up <- c(A=Amax,B=Bmax,gamma=gamma.max,beta=beta.max)
    }else{
        start <- list(gamma=gamma.ini,beta=beta.ini)
        par.low <- c(gamma=gamma.min,beta=beta.min)
        par.up <- c(gamma=gamma.max,beta=beta.max)
    }
    if(NI>0){
        for(k in 1:NI){
            start[[paste0('d',k)]] <- dini
            par.low <- c(par.low,d=dmin)
            par.up <- c(par.up,d=dmax)
        }
    }
    start$sj <- sj.ini
    par.low <- c(par.low,sj=sj.min)
    par.up <- c(par.up,sj=sj.max)
    if(Nma>0){
        for(k in 1:Nma){
            start[[paste0('m',k)]] <- mini#[k]
            par.low <- c(par.low,m=mmin)
            par.up <- c(par.up,m=mmax)
        }
        start <- c(start,logtau=logtau.ini)
        par.low <- c(par.low,logtau=logtau.min)
        par.up <- c(par.up,logtau=logtau.max)
    }
#    cat('names(start)=',names(start),'\n')
#    cat('start=',unlist(start),'\n')
    names(par.up) <- names(par.low) <- names(start)
###########
    out <- nls.lm(par = start,lower=par.low,upper=par.up,fn = rv.res,data=df,control=nls.lm.control(maxiter=500))#ftol=1e-8
    opt.par <- pars <- as.list(coef(out))
    yp.full <- RV.model(par=pars,data = df)
    if(type=='period'){
        pars$A <- 0
        pars$B <- 0
        yp.noise <- RV.model(par=pars,data=df)
    }else{
        yp.noise <- yp.full
    }
    res <- as.numeric(y-yp.full)
    res.sig <- as.numeric(y-(yp.full-yp.noise))
    chi2 <- sum(res^2/dy^2)
    chi2.noise <- sum((y-yp.noise)^2/dy^2)
    logLmax <- sum(-res^2/(2*(opt.par$sj^2+dy^2))-0.5*log(2*pi)-0.5*log(opt.par$sj^2+dy^2))
    return(list(res=res,res.sig=res.sig,chi2=chi2,chi2.noise=chi2.noise,par=opt.par,logLmax=logLmax))
}
rv.red <- function(par,df){
#####setting up
    var <- names(df)
    for(k in 1:length(var)){
        assign(var[k],df[[var[k]]])
    }
    if(!is.null(par.fix)){
        phi <- par.fix$phi
        omega <- par.fix$omega
    }
    t <- data[,1]
    y <- data[,2]
    dy <- data[,3]
    logtau <- par$logtau
    m <- c()
    if(Nma>0){
        for(j in 1:Nma){
            m <- c(m,par[[paste0('m',j)]])
        }
    }else{
        m <- 0
    }
    sj <- par$sj
#####notations
    W <- sum(1/(sj^2+dy^2))
    w <- 1/(sj^2+dy^2)/W
    es <- c()
    if(Nma>0){
        for(i in 1:Nma){
            ei <- c(rep(0,i),exp(-abs(t[-(length(t)+1-(1:i))]-t[-(1:i)])/exp(logtau)))
            es <- rbind(es,ei)
        }
        yp <- y-m%*%(es*ys)
        wp <- 1-m%*%es
        tp <- t-m%*%(es*ts)
        if(type=='period'){
            cp <- cos(omega*t)-m%*%(es*cs)
            sp <- sin(omega*t)-m%*%(es*ss)
        }
    }else{
        yp <- y
        wp <- 1
        tp <- t
        if(type=='period'){
            cp <- cos(omega*t)
            sp <- sin(omega*t)
        }
    }
    WIp <- ip <- TIp <- YIp <- CIp <- SIp <- c()
    IIp <- array(data=NA,dim=c(NI,NI))
    if(NI>0){
        for(j in 1:NI){
            if(Nma>0){
                ip <- rbind(ip,Indices[j,]-m%*%(es*Is[,j,]))
            }else{
                ip <- rbind(ip,Indices[j,])
            }
            YIp <- c(YIp,sum(w*yp*ip[j,]))
            WIp <- c(WIp,sum(w*wp*ip[j,]))
            TIp <- c(TIp,sum(w*tp*ip[j,]))
            if(type=='period'){
                CIp <- c(CIp,sum(w*cp*ip[j,]))
                SIp <- c(SIp,sum(w*sp*ip[j,]))
            }
        }
        for(j in 1:NI){
            for(i in j:NI){
                IIp[j,i] <- IIp[i,j] <- sum(w*ip[i,]*ip[j,])
            }
        }
    }
    WWp <- sum(w*wp*wp)
    YWp <- sum(w*wp*yp)
    YTp <- sum(w*yp*tp)
    Wp <- sum(w*wp)
    WTp <- sum(w*wp*tp)
    TTp <- sum(w*tp^2)
    if(type=='period'){
        Cp <- sum(w*cp)
        Sp <- sum(w*sp)
        YCp <- sum(w*yp*cp)
        YSp <- sum(w*yp*sp)
        CCp <- sum(w*cp^2)
        SSp <- sum(w*sp^2)
        CSp <- sum(w*cp*sp)
        CTp <- sum(w*cp*tp)
        CWp <- sum(w*cp*wp)
        STp <- sum(w*sp*tp)
        SWp <- sum(w*sp*wp)
    }

    if(type=='noise'){
        if(NI>0){
            lin.mat <- matrix(c(WWp,WTp,WIp,WTp,TTp,TIp),byrow=TRUE,ncol=2+NI)
            for(j in 1:NI){
                lin.mat <- rbind(lin.mat,c(WIp[j],TIp[j],IIp[j,]))
            }
            vec.rh <- c(YWp,YTp,YIp)
        }else{
            dI <- 0
            d <- 0
            lin.mat <- matrix(c(WWp,WTp,WTp,TTp),byrow=TRUE,ncol=2)
            vec.rh <- c(YWp,YTp)
        }
###optimized parameterse for the trend model
    }else if(type=='period'){
#####calculate the optimal white noise model parameters:gamma,beta,dj
        lin.mat <- c()
        if(NI>0){
            for(j in 1:(4+NI)){
                if(j==1){
                    lin.mat <- rbind(lin.mat,c(CCp,CSp,CWp,CTp,CIp))
                }
                if(j==2){
                    lin.mat <- rbind(lin.mat,c(CSp,SSp,SWp,STp,SIp))
                }
                if(j==3){
                    lin.mat <- rbind(lin.mat,c(CWp,SWp,WWp,WTp,WIp))
                }
                if(j==4){
                    lin.mat <- rbind(lin.mat,c(CTp,STp,WTp,TTp,TIp))
                }
                if(j>4){
                    lin.mat <- rbind(lin.mat,c(CIp[j-4],SIp[j-4],WIp[j-4],TIp[j-4],IIp[j-4,]))
                }
            }
            vec.rh <- c(YCp,YSp,YWp,YTp,YIp)
        }else{
            lin.mat <- rbind(lin.mat,c(CCp,CSp,CWp,CTp))
            lin.mat <- rbind(lin.mat,c(CSp,SSp,SWp,STp))
            lin.mat <- rbind(lin.mat,c(CWp,SWp,WWp,WTp))
            lin.mat <- rbind(lin.mat,c(CTp,STp,WTp,TTp))
            vec.rh <- c(YCp,YSp,YWp,YTp)
        }
    }
    white.par <- try(solve(lin.mat,vec.rh,tol=tol))#gamma,beta,dj
    ind0 <- length(white.par)-NI-1
    r <- white.par[ind0]+white.par[ind0+1]*t
    if(NI>0){
        r <- r+white.par[(ind0+2):length(white.par)]%*%Indices[1:NI,]
    }
    if(type=='period'){
        r <- r+white.par[1]*cos(omega*t)+white.par[2]*sin(omega*t)
    }
    if(Nma>0){
        r <- as.numeric(r)
        rs <- c()
        for(i in 1:Nma){
            ri <- c(rep(0,i),r[-(length(t)+1-(1:i))])
            rs <- rbind(rs,ri)
        }
        v <- as.numeric(r+m%*%(es*(ys-rs)))
    }else{
        v <- as.numeric(r)
    }
    opt.par <- white.par
    return(list(v=v,par=opt.par))
}
rv.white <- function(par,df){
#####setting up
    var <- names(df)
    for(k in 1:length(var)){
        assign(var[k],df[[var[k]]])
    }
    if(!is.null(par.fix)){
        phi <- par.fix$phi
        omega <- par.fix$omega
    }
    t <- data[,1]
    y <- data[,2]
    dy <- data[,3]
    sj <- par$sj
#####notations
    W <- sum(1/(sj^2+dy^2))
    w <- 1/(sj^2+dy^2)/W
    if(type=='period'){
        c <- cos(omega*t)
        s <- sin(omega*t)
    }
    Y <- sum(w*y)
    YT <- sum(w*y*t)
    T <- sum(w*t)
    TT <- sum(w*t^2)
    if(NI>0){
        I <- TI <- YI <- c()
        II <- array(data=NA,dim=c(NI,NI))
        for(j in 1:NI){
            I <- c(I,sum(w*Indices[j,]))
            TI <- c(TI,sum(w*t*Indices[j,]))
            YI <- c(YI,sum(w*y*Indices[j,]))
            for(i in j:NI){
                II[j,i] <- II[i,j] <- sum(w*Indices[i,]*Indices[j,])
            }
        }
        if(type=='period'){
            CI <- SI <- c()
            for(j in 1:NI){
                CI <- c(CI,sum(w*cos(omega*t)*Indices[j,]))
                SI <- c(SI,sum(w*sin(omega*t)*Indices[j,]))
            }
        }
    }
####
#    if w is normalized, 1 is used; otherwise, W is used.
    if(type=='noise'){
        if(NI>0){
            lin.mat <- matrix(c(1,T,I,T,TT,TI),byrow=TRUE,ncol=2+NI)
            for(j in 1:NI){
                lin.mat <- rbind(lin.mat,c(I[j],TI[j],II[j,]))
            }
            vec.rh <- c(Y,YT,YI)
            white.par <- solve(lin.mat,vec.rh,tol=df$tol)#gamma,beta,dj
            indd <- (length(white.par)-NI+1):length(white.par)
        }else{
            dI <- 0
            d <- 0
            lin.mat <- matrix(c(1,T,T,TT),byrow=TRUE,ncol=2)
            vec.rh <- c(Y,YT)
            white.par <- solve(lin.mat,vec.rh,tol=tol)
        }
###optimized parameterse for the trend model
    }else if(type=='period'){
        C <- sum(w*cos(omega*t))
        S <- sum(w*sin(omega*t))
        CC <- sum(w*cos(omega*t)^2)
        SS <- sum(w*sin(omega*t)^2)
        CS <- sum(w*sin(omega*t)*cos(omega*t))
        CT <- sum(w*t*cos(omega*t))
        ST <- sum(w*t*sin(omega*t))
        YC <- sum(w*y*cos(omega*t))
        YS <- sum(w*y*sin(omega*t))
#####calculate the optimal white noise model parameters:gamma,beta,dj
        lin.mat <- c()#matrix(c(CC,CS,C,CT,CIW,T,I,T,TT,TI),byrow=TRUE,ncol=2+NI)
        if(NI>0){
            for(j in 1:(4+NI)){
                if(j==1){
                    lin.mat <- rbind(lin.mat,c(CC,CS,C,CT,CI))
                }
                if(j==2){
                    lin.mat <- rbind(lin.mat,c(CS,SS,S,ST,SI))
                }
                if(j==3){
                    lin.mat <- rbind(lin.mat,c(C,S,W,T,I))
                }
                if(j==4){
                    lin.mat <- rbind(lin.mat,c(CT,ST,T,TT,TI))
                }
                if(j>4){
                    lin.mat <- rbind(lin.mat,c(CI[j-4],SI[j-4],I[j-4],TI[j-4],II[j-4,]))
                }
            }
            vec.rh <- c(YC,YS,Y,YT,YI)
            white.par <- solve(lin.mat,vec.rh,tol=tol)#gamma,beta,dj
            indd <- (length(white.par)-NI+1):length(white.par)
        }else{
            lin.mat <- rbind(lin.mat,c(CC,CS,C,CT))
            lin.mat <- rbind(lin.mat,c(CS,SS,S,ST))
            lin.mat <- rbind(lin.mat,c(C,S,W,T))
            lin.mat <- rbind(lin.mat,c(CT,ST,T,TT))
            vec.rh <- c(YC,YS,Y,YT)
            white.par <- solve(lin.mat,vec.rh,tol=tol)#gamma,beta,dj
        }
    }
    ind0 <- length(white.par)-NI-1
    r <- white.par[ind0]+white.par[ind0+1]*t
    if(NI>0){
        r <- r+white.par[(ind0+2):length(white.par)]%*%Indices[1:NI,]
    }
    if(type=='period'){
        r <- r+white.par[1]*cos(omega*t)+white.par[2]*sin(omega*t)
    }
    if(NI>0){
        indd <- (length(white.par)-NI+1):length(white.par)
        opt.par <- c(white.par[-indd],d=white.par[indd])
    }else{
        opt.par <- white.par
    }
    return(list(v=r,par=opt.par))
}
rv.red.res <- function(par,df){
    y <- df$data[,2]
    dy <- df$data[,3]
#    if(df$Nma>0){
    v <- rv.red(par,df)$v
#    }else{
#        v <- rv.white(par,df)$v
#    }
    sj <- par$sj
#    tmp <- sqrt((y-v)^2/(2*(dy^2+sj^2)) + 0.5*log(dy^2+sj^2)+100)
    tmp <- (y-v)^2/(2*(dy^2+sj^2)) + 0.5*log(dy^2+sj^2)+100
    return(tmp)
}

######express other parameters as functions of correlated noise parameters, optimize correlated noise parameters using the LM algorithm
sopt <- function(pars,type='noise'){
    var <- names(pars)
    for(k in 1:length(var)){
        assign(var[k],pars[[var[k]]])
    }
####transpose Indices to make the computation faster
    if(!is.null(Indices) & !any(is.na(Indices))){
        if(nrow(Indices)==length(t)) Indices <- t(Indices)
    }
    if(type=='period'){
        par.fix <- list(omega=omega,phi=phi)
    }else{
        par.fix <- NULL
    }
    ts <- cs <- ys <- ss <- c()
    if(NI>0){
        Is <- array(data=NA,dim=c(Nma,NI,length(t)))
    }else{
        Is <- c()
    }
    if(Nma>0){
        for(i in 1:Nma){
            yi <- c(rep(0,i),y[-(length(t)+1-(1:i))])
            ys <- rbind(ys,yi)
            ti <- c(rep(0,i),t[-(length(t)+1-(1:i))])
            ts <- rbind(ts,ti)
            if(type=='period'){
                ci <- c(rep(0,i),cos(omega*t[-(length(t)+1-(1:i))]))
                cs <- rbind(cs,ci)
                si <- c(rep(0,i),sin(omega*t[-(length(t)+1-(1:i))]))
                ss <- rbind(ss,si)
            }
            if(NI==1){
                Is[i,,] <- c(rep(0,i*NI),Indices[1,-(length(t)+1-(1:i))])
            }else if(NI>1){
                Is[i,,] <- c(matrix(rep(0,i*NI),nrow=NI),Indices[1:NI,-(length(t)+1-(1:i))])
            }
        }
    }
    tmp <- data.frame(t,y,dy)
    df <- list(data=tmp,Indices=Indices,par.fix=par.fix,Nma=Nma,NI=NI,type=type,ts=ts,cs=cs,ss=ss,Is=Is,ys=ys)
####

    start <- list()
    par.low <- par.up <- c()
    start$sj <- sj.ini
    par.low <- c(par.low,sj=sj.min)
    par.up <- c(par.up,sj=sj.max)
    if(Nma>0){
        for(k in 1:Nma){
            start[[paste0('m',k)]] <- mini
            par.low <- c(par.low,mmin)
            par.up <- c(par.up,mmax)
        }
        start <- c(start,logtau=logtau.ini)
        par.low <- c(par.low,logtau=logtau.min)
        par.up <- c(par.up,logtau=logtau.max)
    }
    names(par.up) <- names(par.low) <- names(start)

##########numerical fit
    Ntry <- 1
    if(type=='noise'){
#    if(TRUE){
        start0 <- start
        pars <- c()
        likes <- c()
        for(kk in 1:Ntry){
            for(j in 1:5){
                start <- as.list(rnorm(length(start),unlist(start0),as.numeric((par.up-par.low)/100)))
                if(all(unlist(start)>par.low & unlist(start)<par.up)) break()
            }
            names(start) <- names(start0)
            out <- try(nls.lm(par = start,lower=par.low,upper=par.up,fn = rv.red.res,df=df,control=nls.lm.control(maxiter=500)),TRUE)#ftol=1e-16,ptol=1e-16#ftol=tol,ptol=tol
            if(class(out)!='try-error'){
                opt.par <- as.list(coef(out))
                yp <- rv.red(par=opt.par, df=df)$v
                pars <- rbind(pars,unlist(start))
                logLmax <- sum(-(y-yp)^2/(2*(dy^2+opt.par$sj^2))-0.5*log(2*pi)-0.5*log((dy^2+opt.par$sj^2)))
                likes <- c(likes,logLmax)
            }
        }
        ind <- which.max(likes)
        start <- as.list(pars[ind,])
        names(start) <- names(start0)
    }
    out <- nls.lm(par = start,lower=par.low,upper=par.up,fn = rv.red.res,df=df,control=nls.lm.control(maxiter=500))


#####retrieve parameters
    opt.par <- as.list(coef(out))
    tmp <- rv.red(par=opt.par,df = df)
    nams <- c('gamma','beta')
    if(NI>0) nams <- c(nams,paste0('d',1:NI))
    if(type=='period') nams <- c('A','B',nams)
    names(tmp$par) <- nams
    yp.full <- tmp$v
    opt.par <- c(unlist(tmp$par),opt.par)
####model prediction and chi2
    res <- as.numeric(y-yp.full)
    logLmax <- sum(-res^2/(2*(dy^2+opt.par$sj^2))-0.5*log(2*pi)-0.5*log((dy^2+opt.par$sj^2)))
    chi2 <- sum(res^2/dy^2)
####
    pars <- opt.par
    if(type=='period'){
        pars$A <- 0
        pars$B <- 0
        df$Indices <- t(df$Indices)
        yp.noise <- as.numeric(RV.model(par=pars,data = df))
    }else{
        yp.noise <- yp.full
    }
    res.sig <- as.numeric(y-(yp.full-yp.noise))
    return(list(res=res,res.sig=res.sig,chi2=chi2,logLmax=logLmax,par=opt.par))#chi2.noise=chi2.noise
}
par.optimize <- function(data,Indices,NI,Nma,opt.type='sl',type='noise',omega=0,pars=NULL){
###pars is a list containing different variables
###Indices include indices which linear correlate with y in the model
###data include t,y,dy
    opt.type0 <- opt.type
    var <- names(pars)
    for(k in 1:length(var)){
        assign(var[k],pars[[var[k]]])
    }
    opt.type <- opt.type0
    if(!exists('pars$Nma')){
        pars$Nma <- Nma
    }
    if(!exists('pars$NI')){
        pars$NI <- NI
    }
    t <- data[,1]
    y <- data[,2]
    dy <- data[,3]
    err2 <- dy^2
    W <- sum(1/err2)#new weight sum is 1.
    w <- 1/err2/W#normalized weighting
    Nrep <- 1
    if(opt.type=='rep'){
        Nrep <- 2
    }
    xs <- array(0,dim=c(1,length(t)))
    for(ii in 1:Nrep){
########linear optimization
        if(opt.type=='l'){
            tmp <- lopt(pars=pars,type=type)
        }else if(opt.type=='nl'){
            tmp <- nlopt(pars=pars,type=type)
        }else if(opt.type=='sl'){
            tmp <- sopt(pars=pars,type=type)
        }
    }
###the returned variables are chi2 and best fitted parameters
    return(tmp)
}
par.integral <- function(data,Indices,sj,m,d,type='noise',logtau=NULL,omega=NULL,Nma=0, NI=0){
    if(NI==0){
        dI <- 0
    }else{
        dI <- partI(d,Indices,NI)
    }
    t <- data[,1]
    y <- data[,2]
    dy <- data[,3]
    ind <- which(dy==0)
    dy[ind] <- 1e-6
    W <- sum(1/(dy^2+sj^2))#new weight sum is 1.
    w <- 1/(dy^2+sj^2)/W#normalized weighting
    if(Nma>0){
        Is <- c()
        vs <- c()
        ts <- c()
        trep <- c()
        for(k in 1:Nma){
            if(NI>1){
                Is <- rbind(Is,c(rep(0,k),partI(d,Indices[-(length(t)+1-(1:k)),],NI)))
            }else{
                Is <- 0
            }
            ts <- rbind(ts,c(rep(0,k),t[-(length(t)+1-(1:k))]))#shift forward by k points
            trep <- rbind(trep,c(rep(0,k),t[-(1:k)]))#unshifted but with the previous points chopped
            vs <- rbind(vs,c(rep(0,k),y[-(length(t)+1-(1:k))]))
        }
        c <- m*exp(-abs(trep-ts)/exp(logtau))
    }else{
        c <- 0
        ts <- 0
        vs <- 0
        Is <- 0
    }
    if(Nma>0){
        if(is.matrix(Is) | is.data.frame(Is)){
            yp <- y-colSums(c*vs)-dI-colSums(c*Is)
        }else if(length(Is)>1){
            yp <- y-colSums(c*vs)-dI-c*Is
        }else{
            yp <- y-colSums(c*vs)-dI
        }
        tp <- t+colSums(c*ts)
        wp <- 1+colSums(c)
    }else{
        yp <- y-dI
        tp <- t
        wp <- 1
    }
    YWp <- sum(w*wp*yp)
    Wp <- sum(w*wp)
    WWp <- sum(w*wp^2)
    WTp <- sum(w*wp*tp)
    TTp <- sum(w*tp^2)
    WTp <- sum(w*wp*tp)
    YTp <- sum(w*yp*tp)
    YWp <- sum(w*yp*wp)
    YYp <- sum(w*yp^2)
    logL0 <- log(2*pi/sqrt(WWp*TTp-WTp^2))-log(W)+W*(((YTp*WWp-WTp*YWp)^2/(WWp*TTp-WTp^2)+YWp^2-YYp*WWp)/(2*WWp))
    logLp <- logL <- logL0
    if(type=='period'){
###determine phase to eliminate CS
        if(Nma>0){
            cc <- colSums(c*cos(omega*ts))
            ss <- colSums(c*sin(omega*ts))
        }else{
            ss <- cc <- 0
        }
        CCpp <- sum(w*(cos(omega*t)+cc)^2)
        SSpp <- sum(w*(sin(omega*t)+ss)^2)
        CSpp <- sum(w*(sin(omega*t)+ss)*(cos(omega*t)+cc))
        phi <- 0.5*atan(2*CSpp/(CCpp-SSpp))
        phip <- 0.5*atan(sum(w*sin(2*omega*t))/sum(w*cos(2*omega*t)))
###define notations using the determined phase
        if(Nma>0){
            cc <- colSums(c*cos(omega*ts-phi))
            ss <- colSums(c*sin(omega*ts-phi))
        }else{
            ss <- cc <- 0
        }
        cp <- cos(omega*t-phi)+cc
        sp <- sin(omega*t-phi)+ss
        Cp <- sum(w*cp)
        Sp <- sum(w*sp)
        CCp <- sum(w*cp^2)
        SSp <- sum(w*sp^2)
        CSp <- sum(w*sp*cp)
        CTp <- sum(w*cp*tp)
        STp <- sum(w*sp*tp)
        CWp <- sum(w*cp*wp)
        SWp <- sum(w*sp*wp)
        YCp <- sum(w*yp*cp)
        YSp <- sum(w*yp*sp)
        U <- WWp-CWp^2/CCp-SWp^2/SSp
        R <- CWp*CTp/CCp+SWp*STp/SSp-WTp
        Q <- YWp-CWp*YCp/CCp-SWp*YSp/SSp
        V <- CCp*SSp*TTp*U-SSp*CTp^2*U-CCp*STp^2*U-CCp*SSp*R^2
        if(V<=0){
            cat('V=',V,'is not positive for f=',omega/(2*pi),'!\n')
        }
        X <- SSp*YCp^2*U+CCp*YSp^2*U-YYp*CCp*SSp*U+CCp*SSp*Q^2
        logL <- log((2*pi)^2/sqrt(abs(V)))-2*log(W)+W*((X+(CCp*SSp*YTp*U-YCp*CTp*SSp*U-YSp*STp*CCp*U+CCp*SSp*Q*R)^2/V)/(2*CCp*SSp*U))
    }
    return(list(logL0=logL0,logL=logL,logLp=logLp))
}

#notations depending on frequency omega
local.notation <- function(t,y,dy,Indices,NI,omega,phi){
    W <- sum(1/dy^2)
    w <- 1/dy^2/W
    C <- sum(w*cos(omega*t-phi))
    S <- sum(w*sin(omega*t-phi))
    YC <- sum(w*y*cos(omega*t-phi))
    YS <- sum(w*y*sin(omega*t-phi))
    CC <- sum(w*cos(omega*t-phi)^2)
    SS <- sum(w*sin(omega*t-phi)^2)
    CS <- sum(w*sin(omega*t-phi)*cos(omega*t-phi))
    ST <- sum(w*sin(omega*t-phi)*t)
    CT <- sum(w*cos(omega*t-phi)*t)
    if(NI>0){
        CI <- (w*cos(omega*t-phi))%*%Indices[,1:NI]
        SI <- (w*sin(omega*t-phi))%*%Indices[,1:NI]
    }else{
        CI <- SI <- rep(0,NI)
    }
    vars <- unique(c('omega','phi','C','S','CC','SS','CS','CT','ST','SI','CI','YC','YS'))
    pars <- list()
    for(k in 1:length(vars)){
        if(exists(vars[k])){
            pars[[vars[k]]] <- eval(parse(text = vars[k]))
        }
    }
    return(pars)
}

#####definition of notations
global.notation <- function(t,y,dy,Indices,Nma,NI){
#######notations
    data <- cbind(t,y,dy)
    W <- sum(1/dy^2)
    w <- 1/dy^2/W
    T <- sum(w*t)
    Y <- sum(w*y)
    if(NI>0){
        I <- w%*%Indices[,1:NI]
    }else{
        I <- 0
    }
    YY <- sum(w*y^2)
    YT <- sum(w*y*t)
    if(NI>0){
        YI <- (w*y)%*%Indices[,1:NI]
    }else{
        YI <- 0
    }
    TT <- sum(w*t^2)
    if(NI>0){
        TI <- (w*t)%*%Indices[,1:NI]
    }
    II <- array(data=NA,dim=c(NI,NI))
    if(NI>0){
        for(i in 1:NI){
            II[i,] <- (w*Indices[,i])%*%Indices[,1:NI]
        }
    }
#####parameter boundaries
    gamma.min <- min(y)
    gamma.max <- max(y)
    gamma.ini <- (gamma.min+gamma.max)/2
    if(NI>0){
        dmin <- -2*(max(y)-min(y))/(max(Indices)-min(Indices))
        dmax <- 2*(max(y)-min(y))/(max(Indices)-min(Indices))
        dini <- rep((dmin+dmax)/2,1)
    }else{
        dini <- dmin <- dmax <- 0#rep(0,NI)
    }
    beta.min <- -(max(y)-min(y))/(max(t)-min(t))
    beta.max <- (max(y)-min(y))/(max(t)-min(t))
    beta.ini <- (beta.min+beta.max)/2
    sj.min <- 0
    sj.max <- 2*sd(y)
    sj.ini <- (sj.min+sj.max)/2
    mmin <- -1
    mmax <- 1
    rr <- 1
    mini <- rep((mmin+rr*mmax)/(1+rr),1)
    logtau.min <- log(max(min(diff(t)),-10))
    logtau.max <- log(2*(max(t)-min(t)))
    logtau.ini <- (logtau.min+rr*logtau.max)/(1+rr)
    Amin <- Bmin <- -2*(max(y)-min(y))
    Amax <- Bmax <- 2*(max(y)-min(y))
    Aini <- Bini <- (Amin+Amax)/2
    vars <- unique(c('t','dy','err2','II','T','TT','TI','II','w','W','I','Y','YT','YI','data','Indices','y','logtau','phi','m','xs','Amin','Bmin','Amax','Bmax','logtau.min','logtau.max','mmin','mmax','beta.min','beta.max','dmin','dmax','gamma.max','gamma.min','Aini','Bini','gamma.ini','beta.ini','dini','mini','logtau.ini','sj.ini','sj.max','sj.min'))
    pars <- list()
    for(k in 1:length(vars)){
        if(exists(vars[k])){
            pars[[vars[k]]] <- eval(parse(text = vars[k]))
        }
    }
    return(pars)
}

####Marginalized likelihood periodogram
MLP <- function(t, y, dy, Nma=0, Inds=0,mar.type='part',sj=0,logtau=NULL,ofac=1,fmax=NULL,fmin=NULL,tspan=NULL,sampling='freq',model.type='MA',opt.par=NULL,Indices=NA,MLP.type='sub'){
#    unit <- 365.24#to make the elements of the matrix in the function of 'solve' on the same order
    unit <- 1
    if(all(Inds==0)){
        NI <- 0
#        Indices <- NA
    }else{
        NI <- length(Inds)
        if(NI>0){
            if(NI<2){
                Indices <- matrix(Indices[,Inds],ncol=1)
            }else{
                Indices <- Indices[,Inds]
            }
        }
    }
    if(NI==0) d <- 0
    if(Nma==0) m <- 0
    if(is.null(tspan)){
        tspan <- max(t)-min(t)
    }
    step <- 1/(tspan*ofac)
    if(is.null(fmin)){
        fmin <- 1/(tspan*ofac)
    }
    fnyq <- 0.5*length(y)/tspan
    if(is.null(fmax)){
        fmax <- fnyq
    }
#    cat('Number of independent frequencies=', (fnyq-fmin)*tspan,'\n')
    if(sampling=='logP'){
        logP.max <- log(1/fmin)
        logP.min <- log(1/fmax)
        f <- 1/exp(seq(logP.min,logP.max,length.out=1000*ofac))*unit
    }else{
        f <- seq(fmin,fmax,by=step)*unit
    }
    nout <- length(f)
    t <- (t-min(t))/unit#
    Ndata <- length(t)
    omegas <- 2*pi*f
    phi <- 0
#######define notations and variables
    y0 <- y
    vars <- global.notation(t,y,dy,Indices,Nma,NI)
    var <- names(vars)
    for(k in 1:length(var)){
        assign(var[k],vars[[var[k]]])
    }
    data <- cbind(t,y,dy)
############################################################
#####optimization; select the initial condition
############################################################
####fix the non-marginzed parameters at their optimal values
    if(is.null(opt.par) | MLP.type=='sub'){
        tmp <- par.optimize(data,Indices=Indices,NI=NI,Nma=Nma,opt.type='sl',type='noise',pars=vars)
        opt.par <- tmp$par
        if(MLP.type=='sub'){
            y <- tmp$res
            data[,2] <- y
            vars$y <- y
            vars$data <- data
        }
    }
    if(Nma>0 & MLP.type=='assign'){
        ind <- grep('^m',names(opt.par))
        m <- unlist(lapply(ind,function(i) opt.par[[i]]))
        logtau <- opt.par$logtau
    }else{
        m <- 0
        logtau <- 1
    }
    if(NI>0 & MLP.type=='assign'){
        ind <- grep('^d',names(opt.par))
        d <- unlist(lapply(ind,function(i) opt.par[[i]]))
    }else{
        d <- rep(0,NI)
    }
    if(sj==0 & any(sj==names(opt.par))){
        sj <- opt.par$sj
    }
########################################
########marginalized posterior
########################################
    if(!exists('m')) m <- 0
    tmp <- par.integral(data,Indices,sj=sj,m=m,d=d,logtau=logtau,Nma=Nma,NI=NI,type='noise')
    logL0 <- tmp$logL0
    logBFmax <- 0
    p <- pn <- rep(NA,length(omegas))
    logBF.noise <- rep(NA,length(omegas))
    logBF <- rep(NA,length(omegas))
    logLp <- rep(NA,length(omegas))
    for(kk in 1:length(f)){
        omega <- omegas[kk]
        tmp <- par.integral(data,Indices,sj=sj,m=m,d=d,type='period',logtau=logtau,omega=omega,Nma=Nma,NI=NI)
        logL1 <- tmp$logL0#signal dependent noise model log likelihood
        logL <- tmp$logL#likelihood for the full model for f=f[k]
        logLp[kk] <- tmp$logLp
        logBF.noise[kk] <- logL1-logL0
        logBF[kk] <- logL-logL0
    }

    P <- unit/f
    ind <- which(logBF>(max(logBF)+log(0.01)))
    Popt <- (unit/f)[ind]
    omega.opt <- 2*pi*(f[which.max(logBF)]/unit)
####optimized parameter
#opt <- par.optimize(data,Indices,NI,Nma,opt.type=opt.type,type='period',omega=omega,pars=var.new)
    vars <- global.notation(t,y0,dy,Indices,Nma,NI)
    tmp <- local.notation(t,y0,dy,Indices,NI,omega.opt,phi=0)
    var.new <- c(vars,tmp)
#    var.new$y <- y0
    data[,2] <- y0
    tmp <- par.optimize(data,Indices=Indices,NI=NI,Nma=Nma,opt.type='sl',type='period',pars=var.new)
    opt.par <- tmp$par
#    cat('names(opt.par)=',names(opt.par),'\n')
#    cat('opt.par=',unlist(opt.par),'\n')
    inds <- sort(logBF,decreasing=TRUE,index.return=TRUE)$ix
    ps <- P[inds[1:5]]
    power.opt <- logBF[inds[1:5]]
    return(list(P=unit/f, logp=power, ps=ps,power.opt=power.opt,logBF.opt=logBF[ind],logBF=logBF,power=logBF,logBF.noise=logBF.noise,Nma=Nma,NI=NI,par=opt.par,res=tmp$res,sig.level=log(c(10,100,1000))))
}

#####BFP-based model inference/selection
bfp.inf <- function(vars,Indices,Nmas=NULL,NI.inds=NULL){
    var <- names(vars)
    for(k in 1:length(var)){
        assign(var[k],vars[[var[k]]])
    }
    Ndata <- length(t)
    if(is.null(Nmas)){
        Nmas <- c(0,1,2)
    }
    if(is.null(NI.inds)){
        NI.inds <- list(0,1:3,1:5,c(1:3,6:10),c(1:3,11:18))
    }
    nis <- c()
    for(j in 1:length(NI.inds)){
        if(any(NI.inds[[j]]==0)){
            nis <- c(nis,0)
        }else{
            nis <- c(nis,length(NI.inds[[j]]))
        }
    }
    ni.min <- min(nis)
    Nma.opt <- Nmas[1]
    Inds.opt <- NI.inds[[1]]
    NI.opt <- length(Inds.opt)

    logLmaxs <- array(data=NA,dim=c(length(NI.inds),length(Nmas)))
    logBFs <- array(data=NA,dim=c(length(NI.inds),length(Nmas)))
    withProgress(message = 'Calculating log(BF) table', value = 0, {
    for(i in 1:length(Nmas)){
        nma <- Nmas[i]
        for(j in 1:length(NI.inds)){
            incProgress(1/(length(Nmas)*length(NI.inds)), detail = paste("MA:", nma,'; Proxies:',paste(NI.inds[[j]],collapse=',')))
            if(!all(NI.inds[[j]]==0)){
                ni <- length(NI.inds[[j]])
                Inds <- NI.inds[[j]]
                if(length(Inds)==1){
                    indices <- matrix(Indices[,Inds],ncol=1)
                }else{
                    indices <- Indices[,Inds]
                }
            }else{
                ni <- 0
                indices <- Indices
                Inds <- 0
            }
            Ntry <- 10*(round(ni/3)+2*nma)
            lls <- c()
            vars <- global.notation(t,y,dy,Indices=indices,Nma=nma,NI=ni)
            tmp <- par.optimize(data,Indices=indices,NI=ni,Nma=nma,opt.type='sl',type='noise',pars=vars)
            ps.opt <- c(logtau=vars$logtau.ini,m=vars$mini,d=vars$dini,sj=vars$sj.ini)
            pp <- tmp$par
            if(i==1 & j==1) pss <- pp
            lls <- c(lls,tmp$logLmax)
            Nsd <- 100
            if(Ntry>1){
                for(kk in 1:Ntry){
                    vars$NI <- ni
                    vars$Nma <- nma
                    if(nma>0){
                        for(ii in 1:10){
                            vars$logtau.ini <- rnorm(1,pp$logtau,(vars$logtau.max-vars$logtau.min)/Nsd)#
                            if(vars$logtau.ini>vars$logtau.min & vars$logtau.ini<vars$logtau.max) break()
                        }
                        for(ii in 1:10){
                            vars$mini <- rnorm(1,unlist(pp[grepl('^m',names(pp))])[1],(vars$mmax-vars$mmin)/Nsd)
                            if(vars$mini>vars$mmin & vars$mini<vars$mmax) break()
                        }
                    }
                    if(ni>0){
                        vars$dini <- runif(1,dmin,dmax)
                        for(ii in 1:10){
                            vars$dini <- rnorm(1,unlist(pp[grepl('^d',names(pp))][1]),(dmax-dmin)/Nsd)
                            if(vars$dini>vars$dmin & vars$dini<vars$dmax) break()
                        }
                    }
                    vars$sj.ini <- runif(1,sj.min,sj.max)
                    vars$NI <- ni
                    vars$Indices <- indices
                    tmp <- try(par.optimize(data,Indices=indices,NI=ni,Nma=nma,opt.type='sl',type='noise',pars=vars))
                    if(class(tmp)!='try-error'){
                        if(tmp$logLmax>max(lls)){
                                pp <- tmp$par
                            }
                            lls <- c(lls,tmp$logLmax)
                    }
                }
            }
            logLmaxs[j,i] <- max(lls)
            if(nma>0){
                pen <- 0.5*(ni-ni.min+nma+1)*log(Ndata)
            }else{
                pen <- 0.5*(ni-ni.min)*log(Ndata)
            }
            logBFs[j,i] <- logLmaxs[j,i]-logLmaxs[1,1]-pen
            cat('log(BF) for Indices=',Inds,'and','Nma=',nma,'is',format(logBFs[j,i],digit=3),'\n\n')
            penI <- log(10)
            penMA <- log(150)
            if(j>1 & i==1){
                if(logBFs[j,i]>(max(logBFs[1:(j-1),1:i])+penI)){
                    Inds.opt <- sort(Inds)
                    Nma.opt <- nma
                    pss <- pp
                }
            }
            if(i>1 & j==1){
                if(logBFs[j,i]>(max(logBFs[1:j,1:(i-1)])+penMA)){
                    Inds.opt <- sort(Inds)
                    Nma.opt <- nma
                    pss <- pp
                }
            }
            if(j>1 & i>1){
                if(logBFs[j,i]>(max(logBFs[1:(j-1),1:i])+penI) & logBFs[j,i]>(max(logBFs[1:j,(1:i-1)])+penMA)){
                    Inds.opt <- sort(Inds)
                    Nma.opt <- nma
                    pss <- pp
                }
            }
        }
    }
    })
    cat('The optimal Nma=',Nma.opt,'Inds=',Inds.opt,'\n')
    return(list(Nma=Nma.opt,Inds=Inds.opt,logBFs=logBFs))
}

BFP.comp <- function(t, y, dy, Nmas,NI.inds=NULL,Indices=NULL){
    t <- t-min(t)
    Ndata <- length(t)
    #######define notations and variables
    if(is.null(NI.inds)){
        vars <- global.notation(t,y,dy,Indices=Indices,Nma=0,NI=0)
    }else{
        vars <- global.notation(t,y,dy,Indices=Indices,Nma=0,NI=0)
    }
    var <- names(vars)
    for(k in 1:length(var)){
        assign(var[k],vars[[var[k]]])
    }
    data <- cbind(t,y,dy)
##########################################
#####optimizing the noise parameters
#########################################
    t1 <- proc.time()
    out <- bfp.inf(vars,Indices,Nmas=Nmas,NI.inds=NI.inds)
    t2 <- proc.time()
    dur <- format((t2-t1)[3],digit=3)
    cat('model comparison computation time:',dur,'s\n\n')
    return(list(logBF=out$logBF,Nma=out$Nma,Inds=out$Inds))
}

combine.data <- function(data,Ninds,Nmas){
###data is a list of matrices
    out <- c()
    idata <- data
    for(j in 1:length(data)){
        Inds <- Ninds[[j]]
        tab <- data[[j]]
        t0 <- tab[,1]%%2400000
        t <- tab[,1]-min(tab[,1])
        y <- tab[,2]
        dy <- tab[,3]
        if(ncol(tab)>3){
            Indices <- as.matrix(tab[,4:ncol(tab)])
        }else{
            Indices <- NA
        }
        if(all(Inds==0)){
            NI <- 0
#            Indices <- NA
        }else{
            NI <- length(Inds)
            if(NI>0){
                if(NI<2){
                    Indices <- matrix(Indices[,Inds],ncol=1)
                }else{
                    Indices <- Indices[,Inds]
                }
            }
        }
        vars <- global.notation(t,y,dy,Indices,Nmas[j],NI)
        if(is.numeric(Indices)) Indices <- matrix(Indices,ncol=1)
        tmp <- par.optimize(data=tab,Indices=Indices,NI=NI,Nma=Nmas[j],opt.type='sl',type='noise',pars=vars)
        val <- cbind(t0,tmp$res,dy)
        idata[[j]] <- val
        out <- rbind(out,val)
    }
    inds <- sort(out[,1],index.return=TRUE)$ix
    return(list(cdata=out[inds,],idata=idata))
}

####Bayes factor periodogram
BFP <- function(t, y, dy, Nma=0, Inds=Inds,Indices=Indices,opt.type='sl',sj=0,logtau=NULL,ofac=1,  norm="Cumming",hifac=1,fmax=NULL,fmin=NA,tspan=NULL,sampling='freq',model.type='MA',progress=TRUE){
  #  unit <- 365.24#to make the elements of the matrix in the function of 'solve' on the same order
    unit <- 1
    if(all(Inds==0)){
        NI <- 0
        Indices <- NA
    }else{
        NI <- length(Inds)
        if(NI>0){
            if(NI<2){
                Indices <- matrix(Indices[,Inds],ncol=1)
            }else{
                Indices <- Indices[,Inds]
            }
        }
    }
    if(NI==0) d <- 0
    if(Nma==0) m <- 0
    if(is.null(tspan)){
        tspan <- max(t)-min(t)
    }
    step <- 1/(tspan*ofac)
    if(is.na(fmin)){
        fmin <- 1/(tspan*ofac)
    }
    fnyq <- 0.5*length(y)/tspan
    if(is.null(fmax)){
        fmax <- hifac*fnyq
    }
    if(sampling=='logP'){
        logP.max <- log(1/fmin)
        logP.min <- log(1/fmax)
        f <- 1/exp(seq(logP.min,logP.max,length.out=1000*ofac))*unit
    }else{
        f <- seq(fmin,fmax,by=step)*unit
    }
    nout <- length(f)
    t <- (t-min(t))/unit#rescale time
    Ndata <- length(t)
    omegas <- 2*pi*f
    phi <- 0
    if(is.null(logtau)){
        logtau <- log(median(diff(t)))
    }
    #######define notations and variables
    vars <- global.notation(t,y,dy,Indices,Nma,NI)
    var <- names(vars)
####shor
    for(k in 1:length(var)){
        assign(var[k],vars[[var[k]]])
    }
    data <- cbind(t,y,dy)
##########################################
#####optimizing the noise parameters
#########################################
    if(model.type=='auto'){
        t1 <- proc.time()
        out <- bfp.inf(vars,Indices,logtau=NULL,m=NULL)
        t2 <- proc.time()
        dur <- format((t2-t1)[3],digit=3)
        cat('model comparison computation time:',dur,'s\n\n')
        NI <- out$NI
        Nma <- out$Nma
        vars <- out$vars
        Indices <- out$Indices
    }
    m <- rep(vars$mini,Nma)
    logtau <- vars$logtau.ini
    d <- rep(vars$dini,NI)
    tmp <- par.optimize(data,Indices=Indices,NI=NI,Nma=Nma,opt.type='sl',type='noise',pars=vars)
#####optimize multiple times just to get global maximum likelihood for the noise model
    for(j in 1:2){
        opt.par.out <- opt.par <- tmp$par
        cor.max <- 0
        vars$sj.ini <- sj <- opt.par$sj
        if(Nma>0){
            m <- c()
            for(j in 1:Nma){
                m <- c(m,opt.par[[paste0('m',j)]])
            }
            vars$mini <- m[1]
            vars$logtau.ini <- logtau <- opt.par$logtau
        }else{
            m <- 0
            logtau <- 1
        }
        if(NI>0){
            ind <- grep('^d',names(opt.par))
            d <- unlist(lapply(ind,function(i) opt.par[[i]]))
            vars$dini <- d[1]
        }else{
            d <- 0
        }
#####optimize again using the optimized parameter as initial conditions
        tmp <- par.optimize(data,Indices=Indices,NI=NI,Nma=Nma,opt.type='sl',type='noise',pars=vars)
        chi2.ref <- tmp$chi2
        logLmax0 <- tmp$logLmax
    }
    p <- rep(NA,length(omegas))
    logLmax <- rep(NA,length(omegas))
    Npar <- length(tmp$par)
    opt.pars <- array(NA,dim=c(length(f),Npar+2))
    if(progress){
        withProgress(message = 'Calculating BFP', value = 0, {
            for(kk in 1:length(f)){
#####for periodogram, set sj=0
                                        #        break()
                omega <- omegas[kk]
                tmp <- local.notation(t,y,dy,Indices,NI,omega,phi)
                var.new <- c(vars,tmp)
                var <- names(tmp)
                for(k in 1:length(var)){
                    assign(var[k],tmp[[var[k]]])
                }
################################################
####I optimization
####################################################
                opt <- par.optimize(data,Indices,NI,Nma,opt.type=opt.type,type='period',omega=omega,pars=var.new)
                logLmax[kk] <- opt$logLmax
                chi2 <- opt$chi2
                opt.pars[kk,] <- unlist(opt$par)
                                        #progress
                incProgress(1/length(omegas), detail = paste("Freq.:",format(f[kk],digit=3),"; Period:", format(2*pi/omegas[kk],digit=3),'d'))
#####power
                p[kk] <- (chi2.ref-chi2)/chi2.ref
            }
        })
    }else{
        for(kk in 1:length(f)){
#####for periodogram, set sj=0
                                        #        break()
            omega <- omegas[kk]
            tmp <- local.notation(t,y,dy,Indices,NI,omega,phi)
            var.new <- c(vars,tmp)
            var <- names(tmp)
            for(k in 1:length(var)){
                assign(var[k],tmp[[var[k]]])
            }
################################################
####I optimization
####################################################
            opt <- par.optimize(data,Indices,NI,Nma,opt.type=opt.type,type='period',omega=omega,pars=var.new)
            logLmax[kk] <- opt$logLmax
            chi2 <- opt$chi2
            opt.pars[kk,] <- unlist(opt$par)
#####power
            p[kk] <- (chi2.ref-chi2)/chi2.ref
        }
    }
    colnames(opt.pars) <- names(opt$par)
    logBF <- logLmax-logLmax0-log(length(y))#BIC-estimated BF; the extra free parameter n=2, which are {A, B}
    P <- unit/f
####signals
    ind.max <- which.max(logBF)
    inds <- which(logBF>(max(logBF)+log(0.01)))
    Popt <- P[inds]
    opt.par <- opt.pars[inds,]
#    cat('dim(opt.pars)=',dim(opt.pars),'\n')
#    cat('colnames(opt.pars)=',colnames(opt.pars),'\n')
    logBF.opt <- logBF[inds]
    ##sort
    if(length(inds)>1){
        ind.sort <- sort(logBF.opt,decreasing=TRUE,index.return=TRUE)$ix
        Popt <- Popt[ind.sort]
        opt.par <- opt.par[ind.sort,]
        logBF.opt <- logBF.opt[ind.sort]
    }
####calculate the residual
    par.fix <- list(omega=2*pi/Popt[1],phi=0)
    df <- list(data=cbind(t,y,dy),Indices=Indices,par.fix=par.fix,Nma=Nma,NI=NI)
    if(is.matrix(opt.par) | is.data.frame(opt.par)){
        pp <- as.list(opt.par[1,])
    }else{
        pp <- as.list(opt.par)
    }
    names(pp) <- colnames(opt.pars)
    yall <- RV.model(pp,data=df)
    ysig <- pp$A*cos(2*pi/Popt[1]*t)+pp$B*sin(2*pi/Popt[1]*t)
    ysigt <- ysig+pp$gamma+pp$beta*t
    ynoise <- yall-ysigt
    ynoiset <- yall-ysig
    res.nst <- y-yall
    res.n <- y-ynoise
    res.nt <- y-ynoiset
    res.s <- y-ysig
    res.st <- y-ysigt
    inds <- sort(logBF,decreasing=TRUE,index.return=TRUE)$ix
    ps <- P[inds[1:5]]
    power.opt <- logBF[inds[1:5]]
    return(list(data=data,logBF=logBF,P=P,Popt=Popt,logBF.opt=logBF.opt,par=as.list(opt.par),res.nst=res.nst,res=res.nst,res.n=res.n,res.s=res.s,res.st=res.st,res.nt=res.nt,sig.level=log(150), power=logBF,ps=ps,power.opt=power.opt))
}

#####moving periodogram
MP <- function(t, y, dy,Dt,nbin,fmax=1,ofac=1,fmin=1/1000,tspan=NULL,Indices=NA,per.type='MLP',...){
    n <- nbin-1
    dt <- (max(t)-min(t)-Dt)/n
    tstart <- min(t)+(0:n)*dt
    tend <- min(t)+(0:n)*dt+Dt
    tmid <- (tstart+tend)/2
    df <- 1/(Dt*ofac)
    rel.powers <- powers <- array(data=NA,dim=c((fmax-fmin)/df+1,nbin))
    ndata <- rep(NA,nbin)
    withProgress(message = 'Calculating moving periodogram', value = 0, {
        for(j in 0:n){
            cat(paste0(Dt,'d window'),j+1,'/',nbin,'\n')
            incProgress(1/nbin, detail = paste0('window ',j+1,'/',nbin))
            inds <- which(t>=min(t)+j*dt & t<min(t)+j*dt+Dt)
            if(any(is.na(Indices)) | is.null(Indices)){
                index <- Indices
            }else{
                index <- Indices[inds,]
            }
            if(per.type=='BGLS'){
                tmp <- bgls(t=t[inds],y=y[inds],err=dy[inds],fmax=fmax,ofac=ofac,fmin=fmin,tspan=Dt)
            }else if(per.type=='GLS'){
                tmp <- gls(t=t[inds],y=y[inds],err=dy[inds],fmax=fmax,ofac=ofac,fmin=fmin,tspan=Dt)
            }else if(per.type=='GLST'){
                tmp <- glst(t=t[inds],y=y[inds],err=dy[inds],fmax=fmax,ofac=ofac,fmin=fmin,tspan=Dt)
            }else if(per.type=='MLP'){
                tmp <- MLP(t=t[inds],y=y[inds],dy=dy[inds],fmax=fmax,ofac=ofac,fmin=fmin,tspan=Dt,Indices=index,...)
            }else if(per.type=='BFP'){
                tmp <- BFP(t=t[inds],y=y[inds],dy=dy[inds],fmax=fmax,ofac=ofac,fmin=fmin,tspan=Dt,Indices=index,...)
            }else if(per.type=='LS'){
                tmp <- lsp(times=t[inds],x=y[inds],ofac=ofac,from=fmin,to=fmax,tspan=Dt,alpha=c(0.1,0.01,0.001))
                ps <- 1/tmp$scanned
                pps <- tmp$power
                inds <- sort(ps,index.return=TRUE,decreasing=TRUE)$ix
                tmp[['power']] <- pps[inds]
                tmp[['P']] <- ps[inds]
            }
            index <- sort(tmp$P,index.return=TRUE)$ix
            powers[,j+1] <- tmp$power[index]
            rel.powers[,j+1] <- (tmp$power[index]-mean(tmp$power[index]))/(max(tmp$power[index])-mean(tmp$power[index]))
            if(per.type=='MLP' | per.type=='BGLS'){
                powers[,j+1] <- tmp$power[index]-max(tmp$power[index])
            }else{
                powers[,j+1] <- tmp$power[index]
            }
            ndata[j+1] <- length(inds)
        }
    })
    return(list(tmid=tmid,P=tmp$P[index],powers=powers,rel.powers=rel.powers,ndata=ndata))
}

####MP without progress
MP.norm <- function(t, y, dy,Dt,nbin,fmax=1,ofac=1,fmin=1/1000,tspan=NULL,Indices=NA,per.type='MLP',...){
    n <- nbin-1
    dt <- (max(t)-min(t)-Dt)/n
    tstart <- min(t)+(0:n)*dt
    tend <- min(t)+(0:n)*dt+Dt
    tmid <- (tstart+tend)/2
    df <- 1/(Dt*ofac)
    rel.powers <- powers <- array(data=NA,dim=c((fmax-fmin)/df+1,nbin))
    ndata <- rep(NA,nbin)
        for(j in 0:n){
            cat(paste0(Dt,'d window'),j+1,'/',nbin,'\n')
            inds <- which(t>=min(t)+j*dt & t<min(t)+j*dt+Dt)
            if(any(is.na(Indices)) | is.null(Indices)){
                index <- Indices
            }else{
                index <- Indices[inds,]
            }
            if(per.type=='BGLS'){
                tmp <- bgls(t=t[inds],y=y[inds],err=dy[inds],fmax=fmax,ofac=ofac,fmin=fmin,tspan=Dt)
            }else if(per.type=='GLS'){
                tmp <- gls(t=t[inds],y=y[inds],err=dy[inds],fmax=fmax,ofac=ofac,fmin=fmin,tspan=Dt)
            }else if(per.type=='GLST'){
                tmp <- glst(t=t[inds],y=y[inds],err=dy[inds],fmax=fmax,ofac=ofac,fmin=fmin,tspan=Dt)
            }else if(per.type=='MLP'){
                tmp <- MLP(t=t[inds],y=y[inds],dy=dy[inds],fmax=fmax,ofac=ofac,fmin=fmin,tspan=Dt,Indices=index,...)
            }else if(per.type=='BFP'){
                tmp <- BFP(t=t[inds],y=y[inds],dy=dy[inds],fmax=fmax,ofac=ofac,fmin=fmin,tspan=Dt,Indices=index,progress=FALSE,...)
            }else if(per.type=='LS'){
                tmp <- lsp(times=t[inds],x=y[inds],ofac=ofac,from=fmin,to=fmax,tspan=Dt,alpha=c(0.1,0.01,0.001))
                ps <- 1/tmp$scanned
                pps <- tmp$power
                inds <- sort(ps,index.return=TRUE,decreasing=TRUE)$ix
                tmp[['power']] <- pps[inds]
                tmp[['P']] <- ps[inds]
            }
            index <- sort(tmp$P,index.return=TRUE)$ix
            powers[,j+1] <- tmp$power[index]
            rel.powers[,j+1] <- (tmp$power[index]-mean(tmp$power[index]))/(max(tmp$power[index])-mean(tmp$power[index]))
            if(per.type=='MLP' | per.type=='BGLS'){
                powers[,j+1] <- tmp$power[index]-max(tmp$power[index])
            }else{
                powers[,j+1] <- tmp$power[index]
            }
            ndata[j+1] <- length(inds)
        }
    return(list(tmid=tmid,P=tmp$P[index],powers=powers,rel.powers=rel.powers,ndata=ndata))
}
