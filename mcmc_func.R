###In this new file, I have chagned the calcuation memthod of mu1 and correct a bug in the expression of Vr.ma; 
library(e1071)
library(kernlab)
library(glasso)
library(JPEN)
library(matrixcalc)
library(mvtnorm)
library(MASS)
library(rootSolve)
###modeL Of Keplerian Motions Of planets
###According to http://w.astro.berkeley.edu/~kclubb/pdf/RV_Derivation.pdf and Ford 2006
###The mass is from Berger, D. H.; et al. (2006). "First Results from the CHARA Array. IV. The Interferometric Radii of Low-Mass Stars". 
###Ms=0.404
###assign names to the parameter vectors of a RV model
#Naming system: Keplerian parameters: {{per},{K},{e},{omega},{Mo}}; trend pars: {a, b}; noise pars (white noise: {s}, Gaussian process red noise: {sigma.white,l}, ARMA red noise:{{phi},alpha,{w},beta}>)
##kepler solver
if(!exists('nepoch') | !exists('Nepoch')){
    Nepoch <- nepoch <- 1
    ep.par <- ''
}
kep.mt <- function(m,e){
    tol = 1e-10
    E <- rep(NA,length(m))
    Ntt <- 10000
    for(j in 1:length(m)){
        E0 <- m[j]
        for(k in 1:Ntt){
            E1 = E0-(E0-e*sin(E0)-m[j])/(sqrt((1-e*cos(E0))^2-(E0-e*sin(E0)-m[j])*(e*sin(E0))))
            if(abs(E1-E0)<tol) break()
            if(j==Ntt) cat('Keplerian solver does not converge!\n')
            E0 <- E1
        }
        E[j] <- E1%%(2*pi)
    }
    return(E)
}
kep.mt2 <- function(m,e){
    tol = 1e-8
    E0 <- m
    Ntt <- 1000
    for(k in 1:Ntt){
        E1 = E0-(E0-e*sin(E0)-m)/(sqrt((1-e*cos(E0))^2-(E0-e*sin(E0)-m)*(e*sin(E0))))
        if(all(abs(E1-E0)<tol)) break()
        if(k==Ntt) cat('Keplerian solver does not converge!\n')
        E0 <- E1
    }
    return(E1%%(2*pi))
}
##refer to Murison A Practical Method for Solving the Kepler Equation
eps3 <- function(m,e,x){
    t1 <- cos(x)
    t2 <- -1+e*t1
    t3 <- sin(x)
    t4 <- e*t3
    t5 <- -x+t4+m
    t6 <- t5/(1/2*t5*t4/t2+t2)
    return(t5/((1/2*t3-1/6*t1*t6)*e*t6+t2))
}
KeplerStart3 <- function(m,e){
    t34 <- e^2
    t35 <- e*t34
    t33 <- cos(m)
    return(m+(-1/2*t35+e+(t34+3/2*t33*t35)*t33)*sin(m))
}
tol <- 1e-12
kep.murison <- function(m,e,tol=tol){
    E <- rep(NA,length(m))
    Ntt <- 1000
    for(j in 1:length(m)){
        Mnorm <- m[j]%%(2*pi)
        E0 <- KeplerStart3(Mnorm,e)
        for(k in 1:Ntt){
            E1 <- E0-eps3(Mnorm,e,E0)
            if(k==Ntt) 'Kepler solver failed to converge!\n'
            if(abs(E1-E0)<tol) break()
            E0 <- E1
        }
        E[j] <- E1
    }
    return(E)
}
kep.murison2 <- function(m,e,tol=tol){
    Mnorm <- m%%(2*pi)
    E0 <- KeplerStart3(Mnorm,e)
    Ntt <- 1000
    for(k in 1:Ntt){
        E1 <- E0-eps3(Mnorm,e,E0)
        if(k==Ntt) 'Kepler solver failed to converge!\n'
        if(all(abs(E1-E0)<tol)) break()
        E0 <- E1
    }
    return(E1)
}
###my own solver of Keplerion equation
kep <- function(m,e,Mo){
    tol = 1e-8
    E <- rep(NA,length(m))
    Ntt <- 1000
    for(j in 1:length(m)){
        E0 <- m[j]
        for(k in 1:Ntt){
            E1 = m[j]+e*sin(E0)
            if(abs(E1-E0)<tol) break()
            if(j==Ntt) cat('Keplerian solver does not converge!\n')
            E0 <- E1
        }
        E[j] <- E1%%(2*pi)
    }
    return(E)
}
##R package: uniroot
kep.R <- function(m,e){
    E <- rep(NA,length(m))
    for(j in 1:length(m)){
        y <- function(x) x-e*sin(x)-m[j]%%(2*pi)
        E[j] <- uniroot(y,interval=c(0,2*pi))$root
    }
    return(E)
}
divide.pars <- function(pars){
    #select Keplerian pars
    if(Np>0){
        ind.kep <- 1:(Nkeppar*Np)
    }else{
        ind.kep <- c()
    }
    #select noise parameters for different aperture
#    ind.noise <- c()
    ind.trend <- sort((1:Npar)[grepl('\\da',names(pars)) | grepl('\\db',names(pars)) | grepl('\\ds',names(pars))])
    ind.noise <- (1:Npar)[-c(ind.kep,ind.trend)]
#    for(k1 in 1:Nw){
#        tmp <- grep(paste0(k1,'$'),names(pars))
#        tmp1 <- setdiff(tmp,ind.kep)
#        ind.noise = cbind(ind.noise,tmp1)
#    }
    return(list(kep=ind.kep,trend=ind.trend,noise=ind.noise))
}
assign.names <- function(par.vector,Np,p=0,q=0,ntj=Ntj,injection=FALSE){
    nams <- c()
    #names for keplerian part
    if(Np>0){
        if(prior.type=='e0'){
            nams.kep <- c('per','K','omega')
        }else{
            nams.kep <- c('per','K','e','omega','Mo')            
        }
        if(grepl('AK',kep.type)){
            nams.kep <- c(nams.kep,'logtau','ta')
        }
        if(grepl('QP',kep.type)){
            nams.kep <- c(nams.kep,'logtq','aq')
        }
        for(j in 1:Np){
                nams <- c(nams,paste(nams.kep,j,sep=''))
        }
    }
    if(!injection){
        for(i1 in 1:Nw){
            var <- names(par.data[[ids[i1]]])
            for(k in 1:length(var)){
                assign(var[k],par.data[[ids[i1]]][[var[k]]])
            }
            for(j1 in nepoch){
                if(j1==1 | length(nepoch)==1 | any(grepl('^a',ep.par))){
                    nams <- c(nams,paste0(j1,'a',i1))
                }
                if(j1==1 | length(nepoch)==1 | any(grepl('^b',ep.par))){
                    nams <- c(nams,paste0(j1,'b',i1))
                }
                if(j1==1 | length(nepoch)==1 | any(grepl('^s',ep.par))){
                    nams <- c(nams,paste0(j1,'s',i1))
                }
                if(noise.model=='TJ' | noise.model=='ARMATJ' | noise.model=='TJAR' | noise.model=='ARMATJAR'){
                    nams <- c(nams,paste0(j1,'ss',i1))
                    if(ntj==3){
                        nams <- c(nams,paste0(j1,c('sa','sb'),i1))
                    }
                }
                if(noise.model=='GP' | noise.model=='GPR'){
                    nams <- c(nams,paste0(j1,c('sigma.red','l'),i1))
                    if(gp.type=='quasi-periodic'){
                        nams <- c(nams,paste0(j1,c('Pgp','loglp'),i1))
                    }
                }
                if(noise.model=='ARMA' | noise.model=='ARMATJ' | noise.model=='ARMATJAR' | noise.model=='ARMAPSID'){
                                        #AR(p)
                    if(p>0){
                        if(j1==1 | any(grepl('ar',ep.par)) | length(nepoch)==1){
                            nams.ar <- 'phi'
                            for(j in 1:p){
                                nams <- c(nams,paste0(j1,nams.ar,j,i1))
                            }
                            nams <- c(nams,paste0(j1,'alpha',i1))
                        }
                    }
                                        #MA(q)
                    if(q>0){
                        if(j1==1 | any(grepl('ma',ep.par)) | length(nepoch)==1){
                            nams.ma <- 'w'
                            for(j in 1:q){
                                nams <- c(nams,paste0(j1,nams.ma,j,i1))
                            }
                            nams <- c(nams,paste0(j1,'beta',i1))
                        }
                    }
                }
                if(noise.model=='TJAR' | noise.model=='ARMATJAR' | noise.model=='PSID' | noise.model=='ARMAPSID'){
                    for(i0 in 1:Par){
                        nams <- c(nams,paste0(j1,'phis',i0,i1))
                        if(TJAR.sym=='asym'){
                            nams <- c(nams,paste0(j1,'psis',i0,i1))
                        }
                    }
                    nams <- c(nams,'alphas')
                    if(exists('BIS') & exists('FWHM') & TJAR.AB){
                        for(i0 in 1:Par){
                            nams <- c(nams,paste0(j1,'phia',i0,i1))
                            if(TJAR.sym=='asym'){
                                nams <- c(nams,paste0(j1,'psia',i0,i1))
                            }
                        }
                        nams <- c(nams,'alphaa')
                        for(i0 in 1:Par){
                            nams <- c(nams,paste0(j1,'phib',i0,i1))
                            if(TJAR.sym=='asym'){
                                nams <- c(nams,paste0(j1,'psib',i0,i1))
                            }
                        }
                        nams <- c(nams,'alphab')
                    }
                }
                if(noise.model=='TARMA'){
                                        #AR(p)
                    if(p>0){
                        for(j in 1:p){
                            nams <- c(nams,paste0(j1,'phiS',j,i1))
                        }
                        for(j in 1:p){
                            nams <- c(nams,paste0(j1,'phiC',j,i1))
                        }
                        nams <- c(nams,'alpha')
                    }
                                        #MA(q)
                    if(q>0){
                        for(j in 1:q){
                            nams <- c(nams,paste0(j1,'wS',j,i1))
                        }
                        for(j in 1:q){
                            nams <- c(nams,paste0(j1,'wC',j,i1))
                        }
                        nams <- c(nams,paste0(j1,'beta',i1))
                    }
                }
                if(indicator>0){
                    if(j1==1| length(nepoch)==1){
                        nams <- c(nams,paste0(j1,paste0('c',1:indicator),i1))
                        if(index.more){
                            nams <- c(nams,paste0(j1,c('c4','c5','c6','c7'),i1))
                        }
                    }else if(any(grepl('^c',ep.par))){
                        int <- as.integer(gsub('c','',ep.par[grep('^c',ep.par)]))
                        nams <- c(nams,paste0(j1,paste0('c',int),i1))
                    }
                }
                if(Nap>1){
                    for(j in 1:(Nap-1)){
                        if(j1==1| length(nepoch)==1){
                            nams <- c(nams,paste0(j1,'d',j+1,'-',j,'_',i1))
                        }
                    }
                    if(any(grepl('d[^c]',ep.par)) & j1!=1){
                        ind <- grepl('d[^c]',ep.par)
                        par1 <- ep.par[ind]
                        par2 <- gsub('-.+','',par1)
                        nap <- as.integer(gsub('d','',par2))
                        nams <- c(nams,paste0(j1,'d',nap,'-',nap-1,'_',i1))
                    }
                }
                if(Nall>0 & i1==1){
                    if(Naverage>0){
                        nams <- c(nams,'dc0')
                    }
                    if(ndRV>0){
                        if(j1==1 | length(nepoch)==1){
                            nams <- c(nams,paste0(j1,'dc',2:ndRV,'-',1:(ndRV-1),'_',i1))
                        }
                        if(any(grepl('dc',ep.par)) & j1!=1){
                            ind <- grepl('dc',ep.par)
                            par1 <- ep.par[ind]
                            par2 <- gsub('-.+','',par1)
                            nap <- as.integer(gsub('dc','',par2))
                            nams <- c(nams,paste0(j1,'dc',nap,'-',nap-1,'_',i1))
                        }
                    }
                }
                if((j1==1 | length(nepoch)==1 | any(grepl('^m[^a]',ep.par))) & Nmom>0){
                    nams <- c(nams,paste0(j1,mom.name[1:Nmom],'_',i1))
                }
            }
        }
    }
    if(FALSE){
        cat('length(nams)=',length(nams),'\n')
        cat('nams=',nams,'\n')
        cat('length(par.vector)=',length(par.vector),'\n')
    }
    names(par.vector) <- nams
    return(par.vector)
}

plot.labels <- function(Np,noise.model,p=0,q=0,ntj=Ntj){
    nams <- c()
    #names for keplerian part
    if(Np>0){
        nam.kep <- c()
        if(prior.type!='e0'){
            nam.kep <- c(nam.kep,c('period[day]','period[day]',expression(nu*"[day"^{-1}*']'),'K[m/s]','e',expression(omega*'[rad]'),expression(M[0]*"[rad]")))
        }else{
            nam.kep <- c(nam.kep,c('period[day]','period[day]',expression(nu*"[day"^{-1}*']'),'K[m/s]',expression(omega*'[rad]')))
        }
        if(grepl('AK',kep.type)){
            nam.kep <- c(nam.kep,c(expression("log(T[day])"),expression(t[a]*"[day]")))
        }
        if(grepl('QP',kep.type)){
            nam.kep <- c(nam.kep,c(expression("log("*t[q]*"/day)"),expression(a[q]*"[rad]")))
        }
        nams <- c(nams,rep(nam.kep,Np))
    }
    for(i3 in 1:Nw){
        var <- names(par.data[[ids[i3]]])
        for(k in 1:length(var)){
            assign(var[k],par.data[[ids[i3]]][[var[k]]])
        }
#    nams <- c(nams,c('a[m/s/year]','b[m/s]','s[m/s]'))
#    for(j1 in 1:Nes[j3]){
    for(j1 in nepoch){
 	    if(j1==1 | length(nepoch)==1 | any(grepl('^a',ep.par))){
	        nams <- c(nams,paste0(j1,'a',i3,'[m/s/year]'))
	    }
	    if(j1==1 | length(nepoch)==1 | any(grepl('^b',ep.par))){
	        nams <- c(nams,paste0(j1,'b',i3,'[m/s]'))
	    }
	    if(j1==1 | length(nepoch)==1 | any(grepl('^s',ep.par))){
	        nams <- c(nams,paste0(j1,'s',i3,'m/s'))
	    }
    if((noise.model=='TJ' | noise.model=='ARMATJ' | noise.model=='TJAR' | noise.model=='ARMATJAR') & indicator>0){
        nams <- c(nams,'ss[m/s/[Sindex]]')
        if(ntj==3){
            nams <- c(nams,c('sa[m/s/[BIS]]','sb[m/s/[FWHM]]'))
        }
    }
    if(noise.model=='GP' | noise.model=='GPR'){
        if(gp.type=='abs' | gp.type=='squared'){
            nams <- c(nams,c(expression(sigma[red]*'[m/s]'),'l[day]'))
        }
        if(gp.type=='quasi-periodic'){
            nams <- c(nams,c(expression(sigma[red]*'[m/s]'),'l[day]',expression(P[gp]*'[day]'),expression(l[gp])))
        }
    }
    if(noise.model=='ARMA' | noise.model=='ARMATJ' | noise.model=='ARMATJAR' | noise.model=='ARMAPSID'){
#AR(p)
        if(p>0){
            for(i2 in 1:p){
                if(j1==1 | any(grepl('ar',ep.par))| length(nepoch)==1){
                    nams <- c(nams,paste0(j1,expression(phi),i2))
                }
            }
            nams <- c(nams,c(paste0(j1,expression("log("*eta*"[day])"),i3)))
        }
#MA(q)
        if(q>0){
            if(j1==1 | any(grepl('ma',ep.par)) | length(nepoch)==1){
                for(i2 in 1:q){
                    nams <- c(nams,c(paste0(j1,expression(w),i2)))
                }
                nams <- c(nams,c(expression("log("*tau*"[day])")))
            }
        }
    }
    if(noise.model=='TJAR' | noise.model=='ARMATJAR' | noise.model=='PSID' | noise.model=='ARMAPSID'){
        if(Par>0){
            for(i2 in 1:Par){
                if(TJAR.sym=='asym'){
                    nams <- c(nams,c(expression(phi[s]*"[m/s]"),expression(psi[s]*"[m/s]")))
                }else{
                    nams <- c(nams,c(expression(phi[s]*"[m/s]")))
                }
            }
            nams <- c(nams,expression(alpha[s]*"[year"^{-1}*"]"))
            if(exists('BIS') & exists('FWHM') & TJAR.AB){
                for(i2 in 1:Par){
                    if(TJAR.sym=='asym'){
                        nams <- c(nams,c(expression(phi[a]*"[m/s]"),expression(psi[a]*"[m/s]")))
                    }else{
                        nams <- c(nams,c(expression(phi[a]*"[m/s]")))
                    }
                }
                nams <- c(nams,expression(alpha[a]*"[year"^{-1}*"]"))
                for(i2 in 1:Par){
                    if(TJAR.sym=='asym'){
                        nams <- c(nams,c(expression(phi[b]*"[m/s]"),expression(psi[b]*"[m/s]")))
                    }else{
                        nams <- c(nams,c(expression(phi[b]*"[m/s]")))
                    }
                }
                nams <- c(nams,expression(alpha[b]*"[year"^{-1}*"]"))
            }
        }
    }
    if(noise.model=='TARMA'){
#AR(p)
        if(p>0){
            nams <- c(nams,c(rep(expression(phi[S]),p),rep(expression(phi[C]),p),expression(eta*"[log(day)]")))
        }
#MA(q)
        if(q>0){
            nams <- c(nams,c(rep(expression(w[S]),q),rep(expression(w[C]),q),expression(tau*"[log(day)]")))
        }
    }
###indices
            if(indicator>0){
                cname <- c('s','b','f')
                if(j1==1| length(nepoch)==1){
                    nams <- c(nams,paste0(j1,paste0('c',cname[1:indicator]),i3,'[m/s]'))
                    if(index.more){
                        nams <- c(nams,paste0(j1,c('c4','c5','c6','c7'),i3,'[m/s]'))
                    }
                }else if(any(grepl('^c',ep.par))){
                    int <- as.integer(gsub('c','',ep.par[grep('^c',ep.par)]))
                    nams <- c(nams,paste0(j1,paste0('c',cname[int]),i3,'[m/s]'))
                }
            }
####differential RVs
    if(Nap>1){
        for(i1 in 1:(Nap-1)){
            if(j1==1 | any(grepl(paste0('d',i1+1,'-',i1),ep.par))| length(nepoch)==1){
                nams <- c(nams,paste0(j1,'d',i1+1,'-',i1,'_',i3))
            }
        }
    }
    if(Nall>0 & (j1==1 | length(nepoch)==1 | any(grepl('dc',ep.par)))){
        if(Naverage>0){
            nams <- c(nams,'dc0')
        }
        if(ndRV>0){
	    if(j1==1 | length(nepoch)==1){
                nams <- c(nams,paste0(j1,'dc',2:ndRV,'-',1:(ndRV-1),'_',i3))
	    }else if(any(grepl('dc',ep.par)) & j1!=1){
                 ind <- grepl('dc',ep.par)
                 par1 <- ep.par[ind]
                 par2 <- gsub('-.+','',par1)
                 nap <- as.integer(gsub('dc','',par2))
                 nams <- c(nams,paste0(j1,'dc',nap,'-',nap-1,'_',i3))
            }
        }
    }
    if((j1==1 | length(nepoch)==1 | any(grepl('m[^a]',ep.par))) & Nmom>0){
          nams <- c(nams,paste0(j1,mom.name[1:Nmom],'_',i3))
    }
    }
    }
    labs <- vector('expression',length(nams))
    for(k in 1:length(nams)){
        labs[k] <- substitute(lab.name,list(lab.name=nams[k]))
    }
#    label.expressions <- sapply(1:Npar,function(i){as.expression(substitute(y,list(y=as.name(nams[k]))))})
#    label.expressions <- nams
    return(labs)
}
####functions for mcmc algorithms for detection of keplerian signal in RV data
####Keplerian model with 1 planet
RV.kepler <- function(pars.kep,tt=NA,Np.kep=Np,prior.kep=prior.type,period.kep=period.par,act=activity,injection=FALSE,kep.only=FALSE,noise.only=FALSE){
    dVr.kep <- list()
    if(all(is.na(tt))){
        sim.kep <- FALSE
        Nrep <- Nw
    }else{
        sim.kep <- TRUE
        Nrep <- 1
    }
    for(j2 in 1:Nrep){
###assign parameters and data
        var <- names(par.data[[ids[j2]]])
        for(k in 1:length(var)){
            assign(var[k],par.data[[ids[j2]]][[var[k]]])
        }
###
        if(!sim.kep){
            tt <- trv
        }
        tt <- tt-tmin
        if(Np.kep>0 & !noise.only){
            per.mul = pars.kep[grep('per([[:digit:]]{1})',names(pars.kep))]
            K.mul = pars.kep[grep('K([[:digit:]]{1})',names(pars.kep))]
            if(prior.kep!='e0'){
                e.mul = pars.kep[grep('e([[:digit:]]{1})',names(pars.kep))]
                Mo.mul = pars.kep[grep('Mo([[:digit:]]{1})',names(pars.kep))]
            }
            omega.mul = pars.kep[grep('omega([[:digit:]]{1})',names(pars.kep))]
            if(grepl('AK',kep.type)){
                logtau.mul = pars.kep[grep('logtau([[:digit:]]{1})',names(pars.kep))]
                ta.mul = pars.kep[grep('ta([[:digit:]]{1})',names(pars.kep))]
            }
            if(grepl('QP',kep.type)){
                logtq.mul = pars.kep[grep('logtq([[:digit:]]{1})',names(pars.kep))]
                aq.mul = pars.kep[grep('aq([[:digit:]]{1})',names(pars.kep))]
            }
            dVr.p <- rep(0,length(tt))
            for(h in 1:Np.kep){
                if(period.kep=='nv'){
                    nv = per.mul[h]
                    P = 1/nv
                }else if(period.kep=='P'){
                    P = per.mul[h]
                }else if(period.kep=='logP'){
                    logP = per.mul[h]
                    P <- exp(logP)
                }
                K <- K.mul[h]
                omega <- omega.mul[h]
                if(prior.kep!='e0'){
                    e <- e.mul[h]
                    Mo <- Mo.mul[h]
                    m = Mo+2*pi*tt/P#mean anomaly
                    E <- kep.mt2(m,e)
                    T <- 2*atan(sqrt((1+e)/(1-e))*tan(E/2))
                }else{
                    E <- 2*pi*tt/P
                    T <- E
                    e <- 0
                }
                if(grepl('AK',kep.type)){
                    ta <- ta.mul[h]
                    logtau <- logtau.mul[h]
                    ak <- exp(-(tt-ta)^2/(2*exp(logtau)^2))
                }else{
                    ak <- 1
                } 
                if(grepl('QP',kep.type)){
                    logtq <- logtq.mul[h]
                    aq <- aq.mul[h]
                    phiq <- aq*cos(2*pi*t/exp(logtq))
                }else{
                    phiq <- 0
                }
                dVr.p <- dVr.p+K*ak*(cos(omega+T+phiq)+e*cos(omega+phiq))
            }
        }else{
            dVr.p <- rep(0,length(tt))
        }
        ind.na <- which(is.na(tt))
	if(length(ind.na)>0){
            dVr.p[ind.na] <- NA
 	}
###assign Keplerian variation to the vector
        dVr.kep[[ids[j2]]] <- dVr.p
###noise
        if(!kep.only){
            if(indicator>0){
                indexes <- matrix(Sindex,ncol=1)
                if(indicator>1){
                    indexes <- cbind(indexes,BIS,FWHM)
                }
            }
            par.can <- pars.kep[grep(paste0(j2,'$'),names(pars.kep))]
            mom.par = pars.kep[grep('^m',names(par.can))]#pars.kep
            for(j1 in nepoch){
               ind.ok <- which(tt>(epochs[j1]-tmin) & tt<(epochs[j1+1]-tmin))
               ind <- grep(paste0(j1,'a',j2),par.can)
	       if(length(ind)>0){
	           at <- par.can[ind]
	       }else{
	           at <- par.can[paste0(nepoch[1],'a',j2)]
	       }
               ind <- grep(paste0(j1,'b',j2),par.can)
	       if(length(ind)>0){
	           bt <- par.can[ind]
	       }else{
	           bt <- par.can[paste0(nepoch[1],'b',j2)]
	       }
###include trend component the vector
               dVr.kep[[ids[j2]]][ind.ok] <- dVr.kep[[ids[j2]]][ind.ok]+at*tt/365.24+bt 
               if(any(grepl('d[^c]',ep.par))){
                    ind <- grepl('d[^c]',ep.par)
                    par1 <- ep.par[ind]
                    par2 <- gsub('-.+','',par1)
                    nap <- as.integer(gsub('d','',par2))
               }
####indices
               if(indicator>0){
                   ind <- grep(paste0(j1,'c'),names(par.can))#c(c1,c2,c3)
                   if(length(ind)>0){
                       cs <- par.can[ind]
                       if(indicator>1){
                           tmp <- rowSums(t(cs*t(indexes[ind.ok,])))
                       }else{
                           tmp <- cs*indexes[ind.ok]
                       }
                       dVr.kep[[ids[j2]]][ind.ok] <- dVr.kep[[ids[j2]]][ind.ok] + tmp
                   }
               }
####differential RVs
               if(Nap>1 & indicator>0){
                   ind <- grep(paste0(j1,'d\\d'),names(par.can))#dj-j'
                   if(length(ind)>0){
                       ds <- par.can[ind]
                       dVr.kep[[ids[j2]]][ind.ok] <- dVr.kep[[ids[j2]]][ind.ok] + rowSums(t(ds*t(dRV[ind.ok,])))##vectorize the code
                   }
               }
####calibration data
               if(Nall>0){
                   ind <- grep(paste0(j1,'dc\\d'),names(par.can))#dcj-j'
                   if(length(ind)>0){
                       dcs <- par.can[ind]
                       dVr.kep[[ids[j2]]][ind.ok] <- dVr.kep[[ids[j2]]][ind.ok] + rowSums(t(dcs*t(dc[ind.ok,])))
                   }
               }
               if(Nmom>0){
                   ind <- grep(paste0(j1,'m'),names(par.can))#moments
                   if(length(ind)>0){
                       dVr.kep[[ids[j2]]][ind.ok] <- dVr.kep[[ids[j2]]][ind.ok] + rowSums(t(par.can[ind]*t(mom[ind.ok,])))
                   }
               }
           }
        }
    }
    return(dVr.kep)
}
####red noise kernel
ker.red <- function(x,y,sigma.red,l,Pgp=1,lp=1,type='abs'){
    if(type=='abs'){
        return(sigma.red^2*exp(-abs(x-y)/l))
    }
    if(type=='squared'){
        return(sigma.red^2*exp(-(x-y)^2/(2*l^2)))
    }
    if(type=='quasi-periodic'){
        return(sigma.red^2*exp(-(x-y)^2/(2*l^2)-(sin(pi*(x-y)/Pgp))^2/(2*lp^2)))
    }
}
###red noise model
rednoise.cov <- function(sigma.red,l,tt=trv,tol=1e-12,Pgp=1,lp=1,type='abs'){
#method 1: easy but not so many choices of kernels
#    ker <- laplacedot(1/l)
#    cov.red <- kernelMatrix(ker,tt)#+diag(eRV^2)#to guarantee that it is positive definite
#method 2: more time-consuming
#    cov.red <- outer(1:length(RV),1:length(RV),Vectorize(function(u,v) ker.red(sigma.red,l,tt[u],tt[v],Pgp=Pgp,lp=lp,type=gp.type)))
#method 3: 
    SE <- function(x,y) ker.red(x,y,sigma.red=sigma.red,l=l,Pgp=Pgp,lp=lp,type=gp.type)
    cov.red <- outer(tt, tt, SE)
    return(cov.red)
}
##arma noise model
arma <- function(pars,rv.kep,x=par.data[[1]]$RV,t=NA,j3=1,Np.kep=Np,ARMAtype=arma.type){
    dv.arma <- c()
    for(j1 in nepoch){
        if(any(grepl('ma',ep.par))){
            ind.ok <- which(t>epochs[j1] & t<epochs[j1+1])
        }else{
            ind.ok <- 1:length(t)
        }
        y <- x[ind.ok]
        par.can <- pars[grep(paste0(j3,'$'),names(pars))]
        if(ARMAtype!='index'){
            phi <- par.can[grep(paste0(j1,'phi([[:digit:]]{1})'),names(par.can))]
            w <- par.can[grep(paste0(j1,'w([[:digit:]]{1})'),names(par.can))]
            p.ar <- length(phi)
            q.ma <- length(w)
        }else{
            phiS <- par.can[grep('phiS([[:digit:]]{1})',names(par.can))]
            phiC <- par.can[grep('phiS([[:digit:]]{1})',names(par.can))]
            wS <- par.can[grep('wC([[:digit:]]{1})',names(par.can))]
            wC <- par.can[grep('wC([[:digit:]]{1})',names(par.can))]
            p.ar <- length(phiC)
            q.ma <- length(wC)
        }
        dVr.ar <- rep(0,length(y))
        dVr.ma <- rep(0,length(y))
        ##AR(p.ar) model
        if(p.ar>0){
            alpha <- par.can[paste0(nepoch[1],'alpha',j3)]
            if(any(grepl('ar',ep.par))){
                alpha <- par.can[paste0(j1,'alpha',j3)]
            }
            for(j in 1:p.ar){
                if(ARMAtype=='abs'){
                    dVr.ar <- dVr.ar + c(rep(0,j),phi[j]*exp((t[-(length(t)+1-(1:j))]-t[-(1:j)])/exp(alpha))*y[-(length(t)+1-(1:j))])
                }
                if(ARMAtype=='squared'){
                    dVr.ar <- dVr.ar + c(rep(0,j),phi[j]*exp(-(t[-(length(t)+1-(1:j))]-t[-(1:j)])^2/(2*exp(alpha)^2))*y[-(length(t)+1-(1:j))])
                }
                if(ARMAtype=='index'){
                    phiT <- phiS[j]*(Sindex[-(1:j)]/max(Sindex))+phiC[j]
                    dVr.ar <- dVr.ar + c(rep(0,j),phiT*exp((t[-(length(t)+1-(1:j))]-t[-(1:j)])/exp(alpha))*y[-(length(t)+1-(1:j))])
                }
                if(ARMAtype=='index.exponent'){
                    phiT <- exp(-phiS[j]*(Sindex[-(1:j)]/max(Sindex)))*phiC[j]
                    dVr.ar <- dVr.ar + c(rep(0,j),phiT*exp((t[-(length(t)+1-(1:j))]-t[-(1:j)])/exp(alpha))*y[-(length(t)+1-(1:j))])
                }
            }
        }
        ##MA(q.ma) model
        if(q.ma>0){
            beta <- par.can[paste0(j1,'beta',j3)]
            if(any(grepl('ma',j3))){
                beta <- par.can[paste0(j1,'beta',j3)]
            }
            res <- y-rv.kep[ind.ok]
            for(j in 1:q.ma){
                if(ARMAtype=='abs'){
                    dVr.ma <- dVr.ma + c(rep(0,j),w[j]*exp(-abs(t[-(length(t)+1-(1:j))]-t[-(1:j)])/exp(beta))*res[-(length(t)+1-(1:j))])
                }
                if(ARMAtype=='squared'){
                    dVr.ma <- dVr.ma + c(rep(0,j),w[j]*exp(-(t[-(length(t)+1-(1:j))]-t[-(1:j)])^2/(2*exp(beta)^2))*res[-(length(t)+1-(1:j))])
                }
                if(ARMAtype=='index'){
                    wT <- wS[j]*(Sindex[-(1:j)]/max(Sindex))+wC[j]
                    dVr.ma <- dVr.ma + c(rep(0,j),wT*exp((t[-(length(t)+1-(1:j))]-t[-(1:j)])/exp(beta))*res[-(length(t)+1-(1:j))])
                }
                if(ARMAtype=='index.exponent'){
                    wT <- exp(-wS[j]*(Sindex[-(1:j)]/max(Sindex)))*wC[j]
                    dVr.ma <- dVr.ma + c(rep(0,j),wT*exp((t[-(length(t)+1-(1:j))]-t[-(1:j)])/exp(beta))*res[-(length(t)+1-(1:j))])
                }
            }
        }
        dv.arma <- c(dv.arma,dVr.ar+dVr.ma)
    }
    return(dv.arma)
}
# ----- Define a function for plotting a matrix ----- #
myImagePlot <- function(x, ...){
     min <- min(x)
     max <- max(x)
     yLabels <- rownames(x)
     xLabels <- colnames(x)
     title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
       min <- Lst$zlim[1]
       max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
       yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
       xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
       title <- Lst$title
    }
  }
# check for null values
if( is.null(xLabels) ){
   xLabels <- c(1:ncol(x))
}
if( is.null(yLabels) ){
   yLabels <- c(1:nrow(x))
}

layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))

 # Red and green range from 0 to 1 while Blue ranges from 1 to 0
 ColorRamp <- rgb( seq(0,1,length=256),  # Red
                   seq(0,1,length=256),  # Green
                   seq(1,0,length=256))  # Blue
 ColorLevels <- seq(min, max, length=length(ColorRamp))

 # Reverse Y axis
 reverse <- nrow(x) : 1
 yLabels <- yLabels[reverse]
 x <- x[reverse,]

 # Data Map
 par(mar = c(3,5,2.5,2))
 image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
 ylab="", axes=FALSE, zlim=c(min,max))
 if( !is.null(title) ){
    title(main=title)
 }
axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
 axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
 cex.axis=0.7)

 # Color Scale
 par(mar = c(3,2.5,2.5,2))
 image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n")

 layout(1)
}
matrix.power <- function(mat,power){
    for(j1 in 1:ncol(mat)){
        ind <- which(mat[,j1]>1e-3)
        mat[ind,j1] <- mat[ind,j1]^power
    }
    return(mat)
}
#likelihood    
loglikelihood <- function(pars,noise.model){ 
    RV.kep = RV.kepler(pars.kep=pars)
    logLike <- 0
    for(k3 in 1:Nw){
        var <- names(par.data[[ids[k3]]])
        for(k in 1:length(var)){
            assign(var[k],par.data[[ids[k3]]][[var[k]]])
        }
#        if(Np>0 & !quantify){
#            eRV <- erv
#        }
        if(grepl('ARMA',noise.model)){
            RV.arma = arma(pars,j3=k3,rv.kep=RV.kep[[k3]],t=par.data[[k3]]$trv,x=par.data[[ids[k3]]]$RV)
        }
	sigma <- eRV
	if(exists('small.err') & !quantify){
            if(small.err) sigma <- erv
	}
#        for(j1 in 1:Nes[k3]){
        for(j1 in nepoch){
            if(any(grepl(paste0(j1,'s',k3),names(pars)))){
                s <- pars[paste0(j1,'s',k3)]
            }else{
                s <- pars[paste0(nepoch[1],'s',k3)]
            }
        if(noise.model=='GP'){
            sigma.red = pars[paste0(j1,'sigma.red',k3)]
            l = pars[paste0(j1,'l',k3)]
            if(gp.type=='abs' | gp.type=='squared'){
                cov.red = (cov.red.basic)^(1/l)*sigma.red^2
                diag(cov.red) <- diag(cov.red)+sigma^2+s^2
            }else if(gp.type=='quasi-periodic'){
                Pgp <- pars[paste0(j1,'Pgp',k3)]
                lp <- exp(paste0(j1,pars['loglp'],k3))
                                        #            tt <- proc.time()

                ind.ok <- ind.ok[trv[ind.ok,k3]>epochs[j1,k3] & trv[ind.ok,k3]<epochs[j1+1,k3]]
                cov.red = rednoise.cov(sigma.red,l,tt=trv[ind.ok,k3],tol=1e-10,Pgp=Pgp,lp=lp,type=gp.type)
            }
            ind.neg <- 0
            det.red <- determinant(cov.red)
            logLike = try(logLike -1/2*log(det.red)-1/2*(t(RV-RV.kep[[k3]])%*%chol2inv(chol(cov.red))%*%(RV-RV.kep[[k3]]))-length(trv)*log(sqrt(2*pi)),TRUE)
            if(class(logLike)=='try-error'){
                ind.neg <- 1
                logLike <- 0
            }
        }else if(noise.model=='GPR'){
            gpr <- sim.GP(pars[paste0(j1,'sigma.red',k3)],pars[paste0(j1,'l',k3)])
            rv.gpr <- gpr$RV
            ind.neg <- gpr$flag
            if(ind.neg==0){
                logLike = logLike +sum(dnorm(RV,mean=RV.kep[[k3]]+rv.gpr,sd=sqrt(sigma^2+s^2),log=T))
            }else{
                logLike = 0
            }
        }else if(noise.model=='ARMA'){
            tmp <- sum(dnorm(RV,mean=RV.kep[[k3]]+RV.arma,sd=sqrt(sigma^2+s^2),log=T))
            logLike = logLike +tmp
            ind.neg <- 0
        }else if(noise.model=='white'){
            logLike = logLike +sum(dnorm(RV,mean=RV.kep[[k3]],sd=sqrt(sigma^2+s^2),log=T))
            ind.neg <- 0
        }else if(noise.model=='TJ'){
            ss = pars[paste0(j1,'ss',k3)]#the slope of the jitter with respect to the BIS
            st = Sindex*ss+s#time varying jitter
            if(any(names(pars)=='sa')){
                st <- st+pars[paste0(j1,'sa',k3)]*BIS
            }
            if(any(names(pars)=='sb')){
                st <- st+pars[paste0(j1,'sb',k3)]*FWHM
            }
            logLike = logLike +sum(dnorm(RV,mean=RV.kep[[k3]],sd=sqrt(sigma^2+st^2),log=T))
            ind.neg <- 0
        }else if(noise.model=='TARMA'){
            RV.tarma = arma(pars,j3=k3,ARMAtype='index',rv.kep=RV.kep[[k3]],x=par.data[[ids[k3]]]$RV,t=par.data[[k3]]$trv)
            logLike = logLike +sum(dnorm(RV,mean=RV.kep[[k3]]+RV.tarma,sd=sqrt(sigma^2+s^2),log=T))
            ind.neg <- 0
        }else if(noise.model=='ARMATJ'){
            ss = pars['ss']#the slope of the jitter with respect to the BIS
            st = Sindex*ss+s#time varying jitter
            if(any(names(pars)==paste0(j1,'sa',k3))){
                st <- st+paste0(j1,pars['sa'],k3)*BIS
            }
            if(any(names(pars)==paste0(j1,'sb',k3))){
                st <- st+pars[paste0(j1,'sb',k3)]*FWHM
            }
            logLike = logLike +sum(dnorm(RV,mean=RV.kep[[k3]]+RV.arma,sd=sqrt(sigma^2+st^2),log=T))
            ind.neg <- 0
        }
            ind.neg <- 0
    }
    }
    return(list(loglike = logLike,ind.neg = ind.neg))
}

tjar <- function(t=trv,x=Sindex[,1],phi.sab,alpha.sab,psi.sab,par.tj=1,symmetry=TJAR.sym){
    rv.tjar <- 0
    for(i0 in 1:length(phi.sab)){
        if(symmetry=='sym'){
            rv.tjar <- rv.tjar + par.tj*c(rep(0,i0),phi.sab[i0]*exp(alpha.sab/365.24*(t[-(length(t)+1-(1:i0))]-t[-(1:i0)]))*x[-(length(t)+1-(1:i0))]) + par.tj*c(psi.sab[i0]*exp(-alpha.sab/365.24*(t[-(1:i0)]-t[-(length(t)+1-(1:i0))]))*x[-(1:i0)],rep(0,i0)) 
        }else if(symmetry=='past'){
            rv.tjar <- rv.tjar + par.tj*c(rep(0,i0),phi.sab[i0]*exp(alpha.sab/365.24*(t[-(length(t)+1-(1:i0))]-t[-(1:i0)]))*x[-(length(t)+1-(1:i0))]) 
        }else if(symmetry=='future'){
            rv.tjar <- rv.tjar + par.tj*c(phi.sab[i0]*exp(-alpha.sab/365.24*(t[-(1:i0)]-t[-(length(t)+1-(1:i0))]))*x[-(1:i0)],rep(0,i0)) 
        }
    }
    return(rv.tjar)
}
prior.func <- function(pars){
#    at = pars['a']
#    bt = pars['b']
    logprior <- 0
    if(Np>0){
        per.mul = pars[grep('per([[:digit:]]{1})',names(pars))]
        if(period.par=='logP'){
            if(Np>1){
                logPprior = 1/sum(log(Pups)-log(Plows))#logP is uniform#or P(nv)~1/nv or P(P)~1/P
            }else{
                logPprior = 1/(log(Pmax)-log(Pmin))
            }
            per.prior = logPprior
            logP = per.mul
            P = exp(logP)
        }else if(period.par=='P'){
            P = per.mul
            Pprior = 1/(Pmax-Pmin)
            per.prior = Pprior
        }else if(period.par=='nv'){
            nv = per.mul
            P = 1/nv
            nvprior = 1/(nvmax-nvmin)
            per.prior = nvprior
        }
        if(prior.type=='jeffrey'){
            K = pars[grep('K([[:digit:]]{1})',names(pars))]
            Kprior = (K+K0)^-1*(log(1+Kmax/K0*(Pmin/P)^(1/3)))^-1#according to Table 1 of Ford & Gregory 2007        
        }else{
            Kprior = 1/(Kmax-Kmin)
        }
        e = pars[grep('e([[:digit:]]{1})',names(pars))]
#        omega = pars[grep('omega',names(pars))]
#        Mo = pars[grep('Mo',names(pars))]
        if(prior.type!='e0'){
            eprior = dnorm(e,mean=0,sd=Esd)
        }else{
            eprior = 1
        }
        omegaprior = 1/(2*pi)
        if(prior.type!='e0'){
            Moprior = 1/(2*pi)
        }else{
            Moprior = 1
        }
        logprior = logprior + Np*log(per.prior)+sum(log(eprior))+Np*log(omegaprior)+Np*log(Moprior)
        if(prior.type=='jeffrey'){
            logprior <- logprior+sum(log(Kprior))
        }else{
            logprior <- logprior+Np*log(Kprior)
        }
        if(grepl('AK',kep.type)){
            logtau.prior <- 1/(logtau.max-logtau.min)
            ta.prior <- 1/(ta.max-ta.min)
            logprior = logprior+log(logtau.prior)+log(ta.prior)
        }
        if(grepl('QP',kep.type)){
            logtq.prior <- 1/(logtq.max-logtq.min)
            aq.prior <- 1/(aq.max-aq.min)
            logprior = logprior+log(aq.prior)+log(logtq.prior)
        }
    }
    for(i1 in 1:Nw){
        var <- names(par.data[[ids[i1]]])
        for(k in 1:length(var)){
            assign(var[k],par.data[[ids[i1]]][[var[k]]])
        }
        logprior <- logprior+na*log(1/(amax-amin))
        logprior <- logprior+nb*log(1/(bmax-bmin))
	 if(prior.type!='jeffrey'){
	    logprior <- logprior+nb*log(1/(smax-smin))
	 }
	 logprior <- logprior+nc/3*sum(log(1/(index.max-index.min)))
	 if(Nap>1){
        	 logprior <- logprior+nd*log(1/(dmax-dmin))
	 }
	 if(Nall>0){
        	 logprior <- logprior+ndc*log(1/(dcmax-dcmin))
	 }
         if(nmom>0 & Nmom>0){
             logprior <- logprior+nmom*log(1/(mmax-mmin))
         }
#        for(j1 in 1:Nes[i1]){
        for(j1 in nepoch){
	 ind <- grep(paste0(j1,'s',i1),names(pars))
	 if(length(ind)>0){
	   s <- pars[ind]
	   if(prior.type=='jeffrey'){
               sprior = (s+s0)^(-1)/(log(1+smax/s0))
               logprior = logprior+log(sprior)
           }
	 }else{
           s <- pars[paste0(nepoch[1],'s',j3)]
	 }
        if(noise.model=='TJ' | noise.model=='ARMATJ' | noise.model=='TJAR' | noise.model=='ARMATJAR'){
            ss <- pars[paste0(j1,'ss',i1)]
            ssprior = 1/(ssmax-ssmin)
            logprior <- logprior+log(ssprior)
        }
        if(noise.model=='TJAR' | noise.model=='ARMATJAR' | noise.model=='PSID' | noise.model=='ARMAPSID'){
            phis.prior = 1/(phis.max-phis.min)
                                        #        alphas.prior = 1/(pars['alphas']*log(alphas.max/alphas.min))#
            alphas.prior <- 1/(alphas.max-alphas.min)
            logprior = logprior+Par*log(phis.prior)+(Par>0)*log(alphas.prior)
            if(exists('BIS') & exists('FWHM') & TJAR.AB){
                phia.prior = 1/(phia.max-phia.min)
                phib.prior = 1/(phib.max-phib.min)
                alphaa.prior <- 1/(alphaa.max-alphaa.min)
                alphab.prior <- 1/(alphab.max-alphab.min)
                logprior = logprior+Par*log(phia.prior)+(Par>0)*log(alphaa.prior)+Par*log(phib.prior)+(Par>0)*log(alphab.prior)
            }
            if(TJAR.sym=='asym'){
                psis.prior = 1/(psis.max-psis.min)
                logprior = logprior+Par*log(psis.prior)
                if(exists('BIS') & exists('FWHM') & TJAR.AB){
                    psia.prior = 1/(psia.max-psia.min)
                    psib.prior = 1/(psib.max-psib.min)
                    logprior = logprior+Par*log(psia.prior)+Par*log(psib.prior)
                }
            }
        }
        if(noise.model=='GP'){
            sigma.red <- pars[paste0(j1,'sigma.red',i1)]
            l <- pars['l1']
            ##using Jeffery's prior
            if(prior.type=='jeffery'){
                redprior = (sigma.red+s0)^(-1)/(log(1+smax/s0))
            }else{
                redprior = 1/(smax-smin)
            }
            lprior = 1/(lmax-lmin)
            logprior = logprior+log(redprior)+log(lprior)
            if(gp.type=='quasi-periodic'){
                loglpprior <- 1/(loglpmax-loglpmin)
                Pgpprior <- 1/(Pgpmax-Pgpmin)
                logprior = logprior+log(loglpprior)+log(Pgpprior)
            }
        }
        if(noise.model=='ARMA' | noise.model=='ARMATJ' | noise.model=='ARMATJAR' | noise.model=='ARMAPSID'){
            phiprior = wprior = 1/(phi.max-phi.min)
            if(p>0){
#                alpha.prior = 1/(pars[paste0(j1,'alpha',i1)]*log(alpha.max/alpha.min))#1/(alpha.max-alpha.min)
                alpha.prior <- 1/(alpha.max-alpha.min)
            }else{
                alpha.prior <- 1
            }
            if(q>0){
#                beta.prior = 1/(pars[paste0(j1,'beta',i1)]*log(beta.max/beta.min))#
                beta.prior <- 1/(beta.max-beta.min)
            }else{
                beta.prior <- 1
            }
            logprior = logprior+p*log(phiprior)+q*log(wprior)+(p>0)*log(alpha.prior)+(q>0)*log(beta.prior)
        }
    }
    }
    return(logprior)
}

##########generate samples from prior distributions
prior.gen <- function(Nprior,noise.model){
    priorA <- runif(Nprior,amin,amax)
    priorB <- runif(Nprior,bmin,bmax)
    if(prior.type=='jeffrey'){
#P(s) ~ (s+s0)^(-1)/(log(1+smax/s0))
        Ntt <- 1000
        priorS <- rep(NA,2*Nprior)#one half is for priorS and the other is for redprior
        for(i0 in 1:2*Nprior){
            for(i1 in 1:Ntt){
                priorS[i0] <- runif(1,smin,smax)
                xtry <- runif(1,0,1/(s0+smin))
                if(xtry<1/(s0+priorS[i0])) break()
            } 
            if(i1==Ntt) cat('The rejection samplings reaches the maximum trial limit!\n')
        }
    }else{
        priorS <- runif(2*Nprior,smin,smax)
    }
    priorBasic <- cbind(priorA,priorB,priorS[1:Nprior])
    priorKep <- c()
    if(Np>0){
        for(i2 in 1:Np){
        if(period.par=='logP'){
            priorPer <- runif(Nprior,log(Pmin),log(Pmax))
            priorP <- exp(priorPer)
        }else if(period.par=='P'){
            priorPer <- runif(Nprior,Pmin,Pmax)
            priorP <- priorPer
        }else if(period.par=='nv'){
            priorPer <- runif(Nprior,nvmin,nvmax)
            priorP <- 1/priorPer
        }
        priorK <- rep(NA,Nprior)
        if(prior.type=='jeffrey'){
        Ntt <- 1000
        for(i0 in 1:Nprior){
            for(i1 in 1:Ntt){
                priorK[i0] <- runif(1,Kmin,Kmax)
                xtry <- runif(1,0,1/(K0+Kmin))
                if(xtry<1/(K0+priorK[i0])) break()
            } 
            if(i1==Ntt) cat('The rejection samplings reaches the maximum trial limit!\n')
        }
        }else{
            priorK <- runif(Nprior,Kmin,Kmax)
        }
        priorE <- dnorm(3*Nprior,mean=0,sd=Esd)
        ind <- which(priorE>=0 & priorE<=1)
        priorE <- priorE[ind[1:Nprior]]
        priorOmega <- runif(Nprior,0,2*pi)
        priorM0 <- runif(Nprior,0,2*pi)
        priorKep <- cbind(priorKep,priorPer,priorK,priorE,priorOmega,priorM0)
    }
    }
    priorARMA <- c()
    if(noise.model=='ARMA' | noise.model=='ARMATJ' | noise.model=='ARMATJAR' | noise.model=='ARMAPSID'){
        priorPhi <- priorW <- priorAlpha <- priorBeta <- c()
        if(p>0){
            priorPhi <- matrix(runif(p*Nprior,phi.min,phi.max),nrow=Nprior)
            priorAlpha <- runif(Nprior,alpha.min,alpha.max)
        }
        if(q>0){
            priorW <- matrix(runif(Nprior*q,wmin,wmax),nrow=Nprior)
            priorBeta <- runif(Nprior,beta.min,beta.max)
        }
        priorARMA <- cbind(priorPhi,priorAlpha,priorW,priorBeta)
    }
    priorTJ <- c()
    if(noise.model=='TJ' | noise.model=='ARMATJ' | noise.model=='TJAR' | noise.model=='ARMATJAR'){
        priorSs <- runif(Nprior,ssmin,ssmax)
        priorTJ <- priorSs
        if(Ntj==3){
            priorSa <- runif(Nprior,ssmin,ssmax)
            priorSb <- runif(Nprior,ssmin,ssmax)
            priorTJ <- cbind(priorTJ,priorSa,priorSb)
        }
    }
    priorGP <- c()
    if(noise.model=='GP' | noise.model=='GPR'){
        priorRed <- priorS[(Nprior+1):length(priorS)]
        priorL <- runif(Nprior,lmin,lmax)
        priorGP <- cbind(priorRed,priorL)
    }
    return(cbind(priorKep,priorBasic,priotTJ,priorARMA,priorGP))
}

#posterior distribution
posterior <- function(param,noise.model){
    llike <- loglikelihood(param)
    pr <- prior.func(param)
    post <- llike$loglike*tem+pr
    return(list(loglike=llike$loglike,logprior=pr,post=post,ind.neg=llike$ind.neg))
}

#Metropolis algorithm####
#####Initial proposal function
proposalfunction <- function(param,cov.adapt){
    Ntt <- 1000
    for(k in 1:Ntt){
        cov.small <- eps*diag(Npar)
        cov.all <- cov.adapt+cov.small
	if(fixP & Np>1){
#	   cov.all[(1:(Np-1))*Nkeppar-Nkeppar+1,(1:(Np-1))*Nkeppar-Nkeppar+1] <- 0
	   cov.all[(1:(Np-1))*Nkeppar-Nkeppar+1,]<- 0
	   cov.all[,(1:(Np-1))*Nkeppar-Nkeppar+1]<- 0
	}
        if(mean(eRV)/(max(RV)-min(RV))<0.01){
#            cov.small <- cov.small*
#            tmp <- det(cov.adapt+cov.small)
            if(!is.positive.definite(cov.all,tol=1e-10)){
#                cat('Covariance is not positively definite and use cov.small instead!\n')
                cov.all <- cov.small
            }
        }
        param.new <- try(mvrnorm(n=1,mu=param,Sigma=cov.all,tol=1e-10,empirical=FALSE))#this could cause some non-zero values for fix value, but it is too small to be accounted for.
#        cat('param.new[k]=',param.new['K1'],'\n')
        if(class(param.new)=='try-error'){
            param.new <- try(mvrnorm(n=1,mu=param,Sigma=cov.small,tol=1e-10,empirical=FALSE))
        }
        if(Np>0){
            M0.sim <- param.new[grep('Mo([[:digit:]]{1})',names(param.new))]
            param.new[grep('Mo([[:digit:]]{1})',names(param.new))] <- M0.sim%%(2*pi)
            omega.sim <- param.new[grep('omega([[:digit:]]{1})',names(param.new))]
            param.new[grep('omega([[:digit:]]{1})',names(param.new))] <- omega.sim%%(2*pi)
        }
        logplast <- param.new[grep(paste0('per',Np),names(param.new))]
        if(!quantify){
            logplast <- param.new[grep(paste0('per',Np),names(param.new))]
            if(any(logplast<logPmaxs & logplast>logPmins)){
                logic.per <- TRUE
            }else{
                logic.per <- FALSE
            }
            if(!fixP){
                if(all(param.new>par.min & param.new<par.max) & all(logic.per)) break()
            }else if(fixP & Np>1){
                ind.mul = grep('per([[:digit:]]{1})',names(param.new))
                ind <- ind.mul[-length(ind.mul)]#fixed periods
#		 ind <- ind.mul
                if(all(param.new[-ind]>par.min[-ind] & param.new[-ind]<par.max[-ind]) & all(logic.per)) break()
            }
        }else{
            if(!fixP){
                if(all(param.new>par.min & param.new<par.max)) break()
            }else{
                ind.mul = grep('per([[:digit:]]{1})',names(param.new))
                ind <- ind.mul[-length(ind.mul)]
#                ind <- ind.mul
                if(all(param.new[-ind]>par.min[-ind] & param.new[-ind]<par.max[-ind])) break()
            }
        }
    }
    names(param.new) <- names(startvalue)
    return(param.new)
}

#####decide whether to accept or adjust the new period values
newper <- function(pars0,pars,pers.low,pers.up){
    ind.last <- which(names(pars)==paste0('per',Np))
    if(Np>1){
        ind.mul = grep('per([[:digit:]]{1})',names(pars))
        ind.former <- ind.mul[ind.mul!=ind.last]
        per.former <- pars[ind.former]
        per0.former <- pars0[ind.former]
    }
    per.last <- pars[ind.last]
    tmp <- par.sec(per.last,pers.low,pers.up)
    logic.all <- logic1 <- tmp$logic
    pars[ind.last] <- tmp$par.new
###former period parameters
    if(length(pers.up)>1 & !fixP){
        pf.low <- pers.up[-length(pers.up)]
        pf.up <- pers.low[-1]
        if(Np>1){
            logic2 <- c()
            for(k1 in 1:length(per.former)){
                ind <- which(per0.former[k1]>pf.low & per0.former[k1]<pf.up)
                if(per.former[k1]>pf.low & per.former[k1]<pf.up){
                    logic2 <- c(logic2,TRUE)
                }else{
                    logic2 <- c(logic2,FALSE)
                }
            }
            logic.all <- all(logic1,logic2)
        }
    }
    return(list(pars=pars,logic=logic.all))
}
####
par.sec <- function(par,pars.low,pars.up){
    if(all(par>pars.low & par<pars.up)){
        log1 <- TRUE
        par.new <- par
    }else if(par< max(pars.up) & par>min(pars.low)){
        dpar.low <- min(abs(par-pars.low))
        dpar.up <- min(abs(par-pars.up))
        if(dpar.low<dpar.up){
            ind.min <- which.min(abs(par-pars.low))
            par.new <- min(pars.low[ind.min]+dpar.low,pars.up[ind.min])
        }else{
            ind.min <- which.min(abs(par-pars.up))
            par.new <- max(pars.up[ind.min]-dpar.up,pars.low[ind.min])
        }
        log1 <- TRUE
    }else{
        log1 <- FALSE
        par.new <- par
    }
    return(list(par.new=par.new,logic=log1))
}
####adaptive proposal function
###calculate covariance matrix
covariance.n0 <- function(mat){
    cov.par <- Sd*cov(mat)+Sd*eps*diag(ncol(mat))
    return(cov.par)
}

covariance.rep <- function(parms,cov1,mu1,n){
    if(n%%Nupdate==0){
        Sd <- 2.4^2/Npar
        mu2 <- (mu1*(n-1)+parms)/n#mu1 and mu2 are the mean vectors of parameters for n-1 and n iterations
        N <- n-1
        cov2 <- (N-1)*cov1/N+Sd*(N*mu1%*%t(mu1)-(N+1)*mu2%*%t(mu2)+parms%*%t(parms)+eps*diag(length(parms)))/N
    }else{
        cov2 <- cov1
    }
    if(fixP & Np>1 & any(grepl('per',names(parms)))){
        ind.mul = grep('per([[:digit:]]{1})',names(parms))
        index <- ind.mul[-length(ind.mul)]
#        index <- ind.mul
        cov2[index,] <- 0
        cov2[,index] <- 0
    }
    return(cov2)
}
###run MH algorithm
run.metropolis.MCMC <- function(startvalue,cov.start,iterations){
    chain = array(dim=c(iterations+1,Npar))
    post = loglike = rep(NA,iterations+1)
    colnames(chain) <- names(startvalue)
    chain[1,]<- startvalue
    mu1 <- chain[1,]
    post.out = posterior(chain[1,])
    post[1] <- post.pre <- post.out$post
    loglike[1] <- loglike.pre <- post.out$loglike
    logprior.pre <- post.out$logprior
    ind.neg <- ind.neg.pre <- post.out$ind.neg
    t.start <- proc.time()
    dt0 <- 0
    Ntt <- 1000
    cov.adapt <- array(data=0,dim=c(Npar,Npar))
    for(i in 1:iterations){
        inds <- divide.pars(startvalue)
        ind.kep <- inds$kep
        ind.noise <- inds$noise
        ind.trend <- inds$trend
        for(j in 1:Ntt){
            if(((Nw>1 | length(nepoch)>1)) & independence){
##############adjust covariance for independent variables
#the calculation is time-consuming, so treat Keplerian and noise parameters independently 
                if(length(nepoch)==1){
                    Ntot <- Nw+1
                }else{
                    Ntot <- 1+Nw+Nw*length(nepoch)
                }
#                Nnoise <- length(ind.noise)/(Nw*length(nepoch))#
                Nnoise <- (length(ind.noise)+length(ind.trend))/(Nw*length(nepoch))#regard trend as noise
                Nt <- cumsum(c(0,rep(length(nepoch),Nw)))
                Nc <- Nt[-length(Nt)]
                for(k2 in 1:Ntot){
                    if(k2==1){
                        ind <- ind.kep
                    }else if(any(Nc+(1:Nw)+1==k2)){
                        nw <- which(Nc+(1:Nw)+1==k2)
                        ind <- ind.trend[((nw-1)*3+1):(nw*3)]
                    }else{
                        nw <- max(which((k2-(Nc+(1:Nw)+1))>0))
                        nep <- k2-(Nc[nw]+(1:Nw[nw])+1)
                        ind <- ind.noise[Nc[nw]*Nnoise+((nep-1)*Nnoise+1):(nep*Nnoise)]
                    }
                    cat('ind.noise=',ind.noise,'\n')
                    cat('ind=',ind,'\n')
                    if(Np>0 | k2>1){
                        if(i == n0){
                            cov.adapt[ind,ind] <- covariance.n0(mat=chain[1:n0,ind])
                        }else if(i > n0){
                            cov.adapt[ind,ind] <- covariance.rep(chain[i,ind],cov.adapt[ind,ind],mu1[ind],i)
                        }else{
                            cov.adapt[ind,ind] <- cov.start[ind,ind]
                        }
                    }
                }
            }else{
                if(i == n0 & !fixP){
                    cov.adapt <- covariance.n0(mat=chain[1:n0,])
                }else if(i > n0){
                    cov.adapt <- covariance.rep(chain[i,],cov.adapt,mu1,i)
                }else{
                    cov.adapt <- cov.start
                }
            }
            proposal = proposalfunction(chain[i,],cov.adapt)
            if(fix.par!=''){
                proposal <- fix(proposal,fix.par)
            }
#            if(fixP & Np>1){
#                proposal[(0:(Np-2))*Nkeppar+1] <- 0
#            }
            proprop <- posterior(proposal)
            if(proprop$ind.neg==0) break()
        }
        post.prop <- proprop$post
	logprior.prop <- proprop$logprior
        loglike.prop <- proprop$loglike
        ind.neg.prop <- proprop$ind.neg

        post.cur <- post.pre
        loglike.cur <- loglike.pre
	logprior.cur <- logprior.pre
        ind.neg.cur <- ind.neg.pre
        if(is.na(post.prop)){
            probab = 0
        }else{
            probab = exp(post.prop-post.cur)
        }
        if(runif(1)<probab){
            chain[i+1,]=proposal
###values for next proposal estimation
            post.pre <- post.prop
            loglike.pre <- loglike.prop
	    logprior.pre <- logprior.prop
            ind.neg.pre <- ind.neg.prop
        }else{
            chain[i+1,]=chain[i,]
            post.pre <- post.cur
            loglike.pre <- loglike.cur
            logprior.pre <- logprior.cur
            ind.neg.pre <- ind.neg.cur
        }
##save values 
        post[i+1]<- logprior.pre + loglike.pre
        loglike[i+1] <- loglike.pre
        mu1 <- (mu1*i+chain[i+1,])/(i+1)
##moniter parameters
#        if((i%%Nverbose==0 | i==1 | i<=100) & verbose){
        if((i%%Nverbose==0 | i==1) & verbose){
#            cat('i=',i,';tem=',tem,';diag(cov.adapt)=',diag(cov.adapt),';')
            cat('i=',i,';tem=',tem,';')
            cat('acceptance percentage:',100*(1-mean(duplicated(chain[1:(i+1),1:Npar]))))
            if(Np>=1){
                pars.now <- chain[i+1,]
                per.mul = pars.now[grep('per([[:digit:]]{1})',names(pars.now))]
                Kmul = pars.now[grep('K([[:digit:]]{1})',names(pars.now))]
                emul = pars.now[grep('e([[:digit:]]{1})',names(pars.now))]
                if(period.par=='nv'){
                    P.cur <- 1/per.mul
                }else if(period.par=='P'){
                    P.cur <- per.mul
                }else if(period.par=='logP'){
                    P.cur <- exp(per.mul)
                }else{
                    cat('The period parameter is not recognized!\n')
                }
                cat('; period= ',P.cur,'; K=', Kmul)
#                cat('; period= ',P.cur,'; Pmin=',Pmin,'; Pmax=',Pmax,'K=', Kmul)
                if(prior.type!='e0'){
                    cat('; e=',emul,'\n')
                }
            }
            cat('; max(logpost)=',max(post[1:(i+1)]),'; maximum likelihood:',max(loglike[1:(i+1)]),'\n')
            cat('duration for the',round(i/Nverbose),'th',Nverbose,' iterations: ',format((proc.time()-t.start)[3],digit=3),'s\n\n')
        }
#        cat('maximum likelihood:',max(loglike[1:(i+1)]),'\n')
    }
#    cat('boundary of period=',exp(par.min[1]),exp(par.max[1]),'d\n')
#    cat('range of period=',exp(range(chain[,(Np-1)*Nkeppar+1])),'d\n')
    return(list(chain=chain,post=post,like=loglike,cov=cov.adapt))
}

binning.post <- function(par.val,post.val,like.val){
    par.min<- min(par.val)
    par.max <- max(par.val)
    bin <- (par.max-par.min)/Nbins
    post.bin <- c()
    like.bin <- c()
    par.bin <- c()
    ind.na <- c()
    for(i in 1:Nbins){
        ind <- which((par.val>=par.min+(i-1)*bin) & (par.val<par.min+i*bin))
        if(length(ind)==0){
            post.bin <- c(post.bin,min(post.val))
	    like.bin <- c(like.bin,min(like.val))
            ind.na <- c(ind.na,i)
        }else{   
#            post.bin <- c(post.bin,max(post.val[ind]))
#            like.bin <- c(like.bin,max(like.val[ind]))	
            post.bin <- c(post.bin,NA)
            like.bin <- c(like.bin,NA)
        }
        par.bin <- c(par.bin,par.min+(2*i-1)*bin/2)
    }
    return(list(likes=like.bin,posts=post.bin,pars=par.bin,ind.na=ind.na))
}
###a more efficient way to binning parameters and log posteriors/likelihoods
binning.post2 <- function(par.val,post.val,like.val,Nbins){
    ind <- sort(par.val,index.return=TRUE)$ix
    par.sort <- par.val[ind]
    post.sort <- post.val[ind]
    like.sort <- like.val[ind]
    p1 <- hist(par.sort,breaks=seq(min(par.sort),max(par.sort),length.out=Nbins+1),plot=FALSE)
    index <- c(0,cumsum(p1$counts))
    post.max <- rep(NA,length(p1$mids))
    like.max <- rep(NA,length(p1$mids))
    ind.na <- c()
    for(j in 1:(length(index)-1)){
        if(index[j+1]>index[j]){
            post.max[j] <- max(post.sort[(index[j]+1):index[j+1]])
            like.max[j] <- max(like.sort[(index[j]+1):index[j+1]])
        }
    }
    ind.na <- which(is.na(post.max))
    return(list(likes=like.max,posts=post.max,pars=p1$mids,ind.na=ind.na))
}
##simulate artificial RV data
sim.RV <- function(par.sim,sim.type,tt=trv[,1],Np.sim,injection=FALSE){
    rv.kep <- RV.kepler(pars.kep=par.sim,tt=tt,Np.kep=Np.sim,prior.kep='mt',period.kep='P',kep.only=injection)
    if(any(names(par.sim)=='s')){
        rv.white <- rnorm(length(tt),mean=0,sd=sqrt(eRV^2+par.sim['s']^2))
    }
    if(sim.type=='white'){
        rv <- rv.kep + rv.white
    }else if(sim.type=='GP' | sim.type=='GPR'){
        sigma.red = par.sim['sigma.red']
        l = par.sim['l']
        s = par.sim['s']
        cov.red.sim = (cov.red.basic)^(1/l)*sigma.red^2+diag(eRV^2+s^2)
        cov.red.sim[cov.red.sim<1e-8] <- 0
        if(!is.positive.definite(cov.red.sim,tol=1e-10)){
#            cat('cov.red.sim is not positive definite!\n')
            eigen.neg <- min(eigen(cov.red.sim)$val)
            cov.red.sim <- cov.red.sim + (1+1e-6)*abs(eigen.neg)*diag(length(trv))
        }
        rv.gp <- mvrnorm(1,rv.kep,cov.red.sim)
        rv <- rv.gp
    }else if(sim.type=='ARMA'){
        rv.arma <- arma(par.sim,Np.kep=Np.sim, x=rv.kep+rv.white,rv.kep=RV.kep[[k3]],t=par.data[[k3]]$trv)
        rv <- rv.arma+rv.kep+rv.white
    }else if(sim.type=='inject'){
        rv <- RV+rv.kep
    }
    return(rv)
}
sim.GP <- function(sigma.red,l){
    cov.red.sim = (cov.red.basic)^(1/l)*sigma.red^2#+1e-6*diag(length(trv))
    cov.red.sim[cov.red.sim<1e-8] <- 0
    flag <- 0
    if(!is.positive.definite(cov.red.sim,tol=1e-10)){
#        cat('cov.red.sim is not positive definite!\n')
#        eigen.neg <- min(eigen(cov.red.sim)$val)
        cov.red.sim <- diag(length(trv))
        flag <- 1
    }
    rv.gp <- mvrnorm(1,rep(0,length(trv)),cov.red.sim)
    return(list(RV=rv.gp,flag=flag))
}
kepler.tm <- function(Kp,wp,ep,Ea){
    Kp*sqrt(1-ep^2)*(cos(wp)*sqrt(1-ep^2)*cos(Ea)-sin(wp)*sin(Ea)/(1-ep*cos(Ea)))
}

kepler.ford <- function(Kp,wp,ep,Ea){
    Tp <- 2*atan(sqrt((1+ep)/(1-ep))*tan(Ea/2))
    rv <- Kp*(cos(wp+Tp)+ep*cos(wp))
    return(rv)
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

data.distr <- function(x,xlab,ylab,main='',oneside=FALSE,plotf=TRUE){
    xs <- seq(min(x),max(x),length.out=1e3)
    fitnorm <- fitdistr(x,"normal")
    p <- hist(x,plot=FALSE)
    xfit <- length(x)*mean(diff(p$mids))*dnorm(xs,fitnorm$estimate[1],fitnorm$estimate[2])
    ylim <- range(xfit,p$counts)
    if(plotf){
        plot(p,xlab=xlab,ylab=ylab,main=main,ylim=ylim)
        lines(xs,xfit,col='red')
    }
    x1=Mode(x)
    x2=mean(x)
    x3=sd(x)
    x4=skewness(x)
    x5=kurtosis(x)
    xs = sort(x)
    x1per = max(min(xs),xs[floor(length(xs)*0.01)])
    x99per = min(xs[ceiling(length(xs)*0.99)],max(xs))
#    abline(v=c(x1per,x99per),col='blue')
    if(plotf){
        if(!oneside){
            legend('topleft',legend=c(as.expression(bquote('mode ='~.(format(x1,digit=3)))),as.expression(bquote(mu~'='~.(format(x2,digit=3)))),as.expression(bquote(sigma~'='~.(format(x3,digit=3))))),bty='n')
            legend('topright',legend=c(as.expression(bquote(mu^3~'='~.(format(x4,digit=3)))),as.expression(bquote(mu^4~'='~.(format(x5,digit=3))))),bty='n')
        }else{
            legend('topleft',legend=c(as.expression(bquote('mode ='~.(format(x1,digit=3)))),as.expression(bquote(mu~'='~.(format(x2,digit=3)))),as.expression(bquote(sigma~'='~.(format(x3,digit=3)))),as.expression(bquote(mu^3~'='~.(format(x4,digit=3)))),as.expression(bquote(mu^4~'='~.(format(x5,digit=3))))),bty='n')
        }
    }
    return(c(x1per=x1per,x99per=x99per,mode=x1,mean=x2,sd=x3,skewness=x4,kurtosis=x5))
}
####function to select output files from MCMC chain
select.file <- function(sub.out){
   fout <- NA  
   if(length(sub.out)!=0){
      ind1 <- which(as.logical(sub.out[,5])==TRUE & as.numeric(sub.out[,3])>10 & as.numeric(sub.out[,3])<35)
      if(length(ind1)!=0){
        out2 <- sub.out[ind1,]
        if(length(ind1)==1){
           tmp <- out2
        }else{
#          ind2 <- which.max(as.numeric(out2[,4]))
           ind2 <- which.min(as.numeric(out2[,3]))
           tmp <- out2[ind2,]
        }
        fout <- as.character(tmp[1])
      }
  }
  if(!is.na(fout)){
    cat('Choose this mcmc results for further investigations: ', fout,'\n')
  }else{
    cat('No qualified mcmc chain for further investigations!\n')
  }
   return(fout)
}

###a function to reset -Inf element in an evidences
inf.rm <- function(E,B,flag=NA){
    ind.m <- which(!is.na(E) & E!=Inf & E!=-Inf,arr.ind=T)
    Emin <- min(E[ind.m])
    ind.b <- which(!is.na(B) & B!=Inf & B!=-Inf,arr.ind=T)
    Bmin <- min(B[ind.b])
    Edim <- dim(E)
    cat('Emin=',Emin,'\n')
    for(i in 1:Edim[3]){
        for(j in 1:Edim[2]){
            if(any(E[,j,i]==-Inf) | any(E[,j,i]==Inf) | any(is.na(E[,j,i])) | B[j,i]==Inf | B[j,i]==-Inf | is.na(B[j,i])){
                if(flag!=NA){
                    E[,j,i] <- Emin
                    B[j,i] <- 0
                }else{
                    E[,j,i] <- NA
                    B[j,i] <- NA
                }
            }
        }
    }
    return(list(E=E,B=B))
}
stemPlot <- function(x,y,pch=16,linecol=1,clinecol=1,add=FALSE,pair=FALSE,...){
    if(!add){
        plot(x,y,pch=pch,...)
    }else{
        points(x,y,pch=pch)
    }
    if(!pair){
        for (i in 1:length(x)){
            lines(c(x[i],x[i]), c(0,y[i]),col=linecol)
        }
    }else{
        ind.pair <- which((1:as.integer(length(x)/2))%%2==1)
        ind.solid <- c(2*ind.pair-1,2*ind.pair)
        ind.dashed <- c(1:length(x))[-ind.solid]
        for (i in 1:length(x)){
            if(any(i==ind.solid)){
                lines(c(x[i],x[i]), c(0,y[i]),col=linecol)
            }else{
                lines(c(x[i],x[i]), c(0,y[i]),col=linecol,lty=2)
            }
        }
    }
 #   lines(c(x[1]-2,x[length(x)]+2), c(0,0),col=clinecol,lty=3)
    abline(h=0,lty=2)
}
###functions
tpm <- function(priors, likes, lambda, h){
    prod.pre <- c(rep(0,h),likes[-((length(likes)-h+1):length(likes))])
    Ptpm <- sum(priors*likes/((1-lambda)*likes*priors+lambda*prod.pre))/sum(priors/((1-lambda)*likes*priors+lambda*prod.pre))
    return(log(Ptpm))
}
calc.alpha <- function(param,cov.adapt,post0,order='01'){
#        par.prop <- rbind(par.prop,par.tmp)
        post.prop <- posterior(param)$post
        if(order=='01'){
            if(is.na(post.prop)){
                probab = 0
            }else{
                probab = exp(post.prop-post0)
            }
        }else{
            if(is.na(post.prop)){
                probab = 1
            }else{
                probab = exp(post0-post.prop)
            }
        }
        min(1,probab)#since q1=q0, probab*q0/q1=probab; ref Chib and Jeliazkov 2001
}
####parallel computing
AMH <- function(nburn,Ns){
    if(nburn==0) nburn <- 1
    out.amh <- foreach(n = 1:Ncores, .combine='comb',.multicombine=TRUE) %dopar% {
        tmp <- run.metropolis.MCMC(startvalue,cov.start,Ns)
        list(out=cbind(tmp$chain[-(1:nburn),],tmp$post[-(1:nburn)],tmp$like[-(1:nburn)]),cov=tmp$cov)
#        cbind(tmp$chain[-(1:nburn),],tmp$post[-(1:nburn)],tmp$like[-(1:nburn)])
    }
#    return(out.amh)
    return(list(out=out.amh$out,cov=out.amh$cov))
}
####parallel chain returning covariance
AMH2 <- function(nburn,Ns){
    if(nburn==0) nburn <- 1
    out.amh <- foreach(n = 1:Ncores, .combine='comb',.multicombine=TRUE) %dopar% {
        tmp <- run.metropolis.MCMC(startvalue,cov.start,Ns)
        list(out=cbind(tmp$chain[-(1:nburn),],tmp$post[-(1:nburn)],tmp$like[-(1:nburn)]),cov=tmp$cov)
    }
    return(list(out=out.amh$out,cov=out.amh$cov))
}
####parallel chain with different temperture
AMH3 <- function(nburn,Ns){
    soft <- TRUE
    if(!soft){
        zeta <- tem.min
        Nlow <- floor(Ncores/2)+1
        Nup <- ceiling(Ncores/2)
        eta <- log(tem/tem.min)/log(Nlow)
        tems <- tempering(1:Nlow,zeta,eta)
        zeta <- tem
        eta <- log(1/tem)/log(Nup)
        tems <- c(tems,sort(tempering(Nup:1,zeta,eta))[-1])
    }else{
        zeta <- tem.min
        eta <- log(tem/tem.min)/log(Ncores)
        tems <- tempering(1:Ncores,zeta,eta)
    }
    tems <- sort(tems,decreasing=TRUE)
    if(nburn==0) nburn <- 1
    out.amh <- foreach(n = 1:Ncores, .combine='comb',.multicombine=TRUE) %dopar% {
        tem <- tems[n]
#        if(n==1){
#        if(n<=Ncores/2){
#        if(n<=3*Ncores/4){
#        if(n<=length(ps) & n<=Ncores/2){
        par.tmp <- pars
#        par.tmp <- par.end
        if(n<=(length(par.tmp)/Npar) & n<=(Ncores/2)){
            if(length(par.tmp)==Npar){
                startvalue <- par.tmp
            }else{
                startvalue <- par.tmp[n,]
            }
#            indP <- (Np-1)*Nkeppar+1
#            if(cov.start[indP,indP]<1e-3){
#                cov.start[indP,indP] <- startvalue[indP]*1e-3
#            }

            tem <- tems0[Ntrace[n]]
#            tem <- 1
        }else{
#            tem <- tems[1]
#            tem <- tems0[1]
            tem <- 1
            source('prepare_par.R',local=TRUE)
        }
        source('mcmc_func.R',local=TRUE)
        tmp <- run.metropolis.MCMC(startvalue,cov.start,Ns)
        list(out=cbind(tmp$chain[-(1:nburn),],tmp$post[-(1:nburn)],tmp$like[-(1:nburn)]),cov=tmp$cov)
    }
    return(list(out=out.amh$out,cov=out.amh$cov))
}
###tempering function; power function
tempering <- function(x,a,b){
    a*x^b
}
####MCMC stuck steps calculation
stuck <- function(arr,post,Nstuck.min){
#    index <- which(duplicated(arr[,((Np-1)*Nkeppar+1):Npar]))
    index <- which(duplicated(arr[,(Np-1)*Nkeppar+1]))
    Ns <- 1
    Nstuck <- c()
    ind <- c()
    ind.stuck <- c()
    for(k in 2:length(index)){
        if((index[k]-index[k-1])==1){
            Ns <- Ns+1
            ind <- c(ind,index[k])
        }else if(length(ind)>0){
            Nstuck <- c(Nstuck,Ns)
            ind.stuck <- c(ind.stuck,ind[length(ind)])
            Ns <- 1
            ind <- c()
        }else{
	    Ns <- 1
            ind <- c()
	}
    }
    index <- which(Nstuck>Nstuck.min)
    if(length(index)>0){
        nstuck <- Nstuck[index]
        indstuck <- ind.stuck[index]
        index1 <- which.max(post[indstuck])
        indstuck <- indstuck[index1]
        nstuck <- nstuck[index1]
    }else{
        index1 <- which.max(Nstuck)
        indstuck <- ind.stuck[index1]
        nstuck <- Nstuck[index1]
    }
    ind.max <- which.max(post)
    per.stuck <- arr[indstuck,(Np-1)*Nkeppar+1]
    per.max <- arr[ind.max,(Np-1)*Nkeppar+1]
    #shift the to the local maxima
    if(abs(per.stuck-per.max)<(0.01*per.stuck)){
        indstuck <- ind.max
    }
#else{
#        indstuck <- 1
#        nstuck <- 1
#    }
    return(list(Nstuck=nstuck,ind.stuck=indstuck))
}
# Print the hostname for each cluster member
GetClusterInfo <- function() {
  info <- Sys.info()[c("nodename", "machine")]
  paste("Node:", info[1], "with CPU type", info[2])
}

###new Plows and Pups for divided period space
period.division <- function(plows,pups,Pmax,Pmin){
    plows.new <- c()
    pups.new <- c()
    if(any(c(plows,pups)>Pmax) | any(c(plows,pups)<Pmin)){
        for(i0 in 1:length(plows)){
            if((pups[i0]<Pmin & plows[i0]<Pmin) | (plows[i0]>Pmax & pups[i0]>Pmax)){
                plows.new <- plows.new
                pups.new <- pups.new
            }else if(plows[i0]<=Pmin & pups[i0]>=Pmin & pups[i0]<=Pmax){
                plows.new <- c(plows.new,Pmin)
                pups.new <- c(pups.new,pups[i0])
            }else if(plows[i0]<=Pmin & pups[i0]>Pmax){
                plows.new <- c(plows.new,Pmin)
                pups.new <- c(pups.new,Pmax)
            }else if(plows[i0]>=Pmin & pups[i0]<=Pmax){
                plows.new <- c(plows.new,plows[i0])
                pups.new <- c(pups.new,pups[i0])
            }else if(plows[i0]>=Pmin & plows[i0]<=Pmax & pups[i0]>Pmax){
                plows.new <- c(plows.new,plows[i0])
                pups.new <- c(pups.new,Pmax)
            }
        }
    }else{
        plows.new <- plows
        pups.new <- pups
    }
    return(list(plow=plows.new,pup=pups.new))
}


###tell foreach how to combine output
comb <- function(x, ...) {  
      mapply(rbind,x,...,SIMPLIFY=FALSE)
}
####calculate timie
time.calc <- function(t1){
    dur <- as.numeric(proc.time()[3]-t1[3])
    time.consumed <- paste(floor(dur/3600),'h',floor((dur%%3600)/60),'m',dur%%60,'s',sep='')
    cat('Time consumed: ',time.consumed,'\n')
}
####weighted time binning
wtb <- function(t,x,ex,dt=1,sj=0){
#    t <- sort(t)
    t0 <- t[1]#the start time point
###note: weight w=1/(ex^2+sj^2)
    ts <- t[1]
    xs <- x[1]
    exs <- ex[1]
    tnew <- c()
    xnew <- c()
    exnew <- c()
    for(i1 in 2:length(t)){
        if((t[i1]-t0)<dt & i1<length(t)){
            ts <- c(ts,t[i1])
            xs <- c(xs,x[i1])
            exs <- c(exs,ex[i1])
        }else{
            tnew <- c(tnew,mean(ts))
            xnew <- c(xnew,xs%*%(1/(exs^2+sj^2))/sum(1/(exs^2+sj^2)))
#            xnew <- c(xnew,mean(xs))
            exnew <- c(exnew,sqrt(exs^2%*%(1/(exs^2+sj^2))^2)/sum(1/(exs^2+sj^2)))
            ts <- t[i1]
            xs <- x[i1]
            exs <- ex[i1]
            t0 <- t[i1]
        }
    }
    return(cbind(tnew,xnew,exnew))
}
#####extract variable values from file names
extract.Nsamp <- function(f){
    Ns <- c()
    for(i in 1:length(f)){
        f1 <- gsub('.*Nsamp','',f[i])
        Ns <- c(Ns,as.integer(gsub('_.+','',f1)))
    }
    return(Ns)
}
extract.tem <- function(f){
    f1 <- gsub('.*tem','',f)
    return(as.numeric(gsub('_acc\\d+','',f1)))
}
extract.Ndata <- function(f){
    if(grepl('Ndata',f)){
        f1 <- gsub('.*Ndata','',f)
    }else{
        f1 <- gsub('.*w0_N','',f)
    }
    return(as.integer(gsub('_.+','',f1)))
}
extract.Lmax <- function(f){
    f1 <- gsub('.+negLmax','',f)
    f2 <- gsub('-.+','',f1)
    f3 <- gsub('.pdf','',f2)
    return(-as.numeric(f3))
}
extract.acc <- function(f){
    f1 <- gsub('.*acc','',f)
    return(as.numeric(f1))
}
combine.index <- function(ind1,ind2,norm=FALSE){
    if(length(ind1)!=0){
        if(!is.matrix(ind1)) ind1 <- matrix(ind1,ncol=1)
        if(!is.matrix(ind2)) ind2 <- matrix(ind2,ncol=1)
        if(nrow(ind2)<nrow(ind1)){
            ind <- cbind(ind1,c(ind2,rep(NA,nrow(ind1)-nrow(ind2))))
        }else if(nrow(ind2)>nrow(ind1)){
            ind1 <- rbind(ind1,matrix(NA,nrow=nrow(ind2)-nrow(ind1),ncol=ncol(ind1)))
            ind <- cbind(ind1,ind2)
        }else{
            ind <- cbind(ind1,ind2)
        }
    }else{
        ind <- ind2
    }
    if(norm){
        ind <- scale(ind)
    }
    if(!is.matrix(ind)){
        ind <- matrix(ind,ncol=1)
    }
    return(ind)
}
####period par transformation
period.transformation <- function(pers,period.type=period.par){
    if(period.type=='logP'){
        Pers <- exp(pers)
    }else if(period.type=='nu'){
        Pers <- 1/pers
    }else if(period.type=='P'){
        Pers <- pers
    }else{
        Pers <- NA
    }
    return(Pers)
}
period.trans2 <- function(ps,period.type=period.par){
    if(period.type=='logP'){
        Pers <- log(ps)
    }else if(period.type=='nu'){
        Pers <- 1/ps
    }else if(period.type=='P'){
        Pers <- ps
    }else{
        Pers <- NA
    }
    return(Pers)
}

K2msini <- function(K,P,e,Ms){
    G <- 4*pi^2*(1.496e11)^3/(365.25*24*3600)^2#m^3*s^-2*Msun^-1
    Me2s <- 3.003e-6#Earth mass in solar unit
    Mj2s <- 1/1048
    Mpj <- K*Ms/Mj2s*(2*pi*G*Ms/(P*24*3600))^(-1/3)*(1-e^2)^0.5#
    Mpe <- Mpj*Mj2s/Me2s
    ap <- ((P/365.25)^2*Ms)^(1/3)#au
#    as <- K*P*24*3600/(2*pi)/1.496e11#m
#    Mp2 <- Ms*(as/ap)*1048#Mj; for circular orbit
#    Mp3 <- K*(P/365.25)^(1/3)*Ms^(2/3)/28.4#Mj; for circular orbit
#    return(c(Mpj,Mpe,ap))
    return(list(mj=Mpj,me=Mpe,a=ap))
}

fix <- function(x,par){
    if(length(par)>0){
        for(j in 1:length(par)){
            ind <- which(par[j]==names(startvalue))
            if(is.matrix(x)){
                x[ind,] <- 0
                x[,ind] <- 0
            }else{
                if(grepl('per',par[j])){
                    x[ind] <- startvalue[ind]
                }else{
                    x[ind] <- 0
                }
            }
        }
    }
    return(x)
}
lag.sta <- function(Ms,m1,a1,e1,m2,e2){
    mu1 <- m1/Ms
    mu2 <- m2/Ms
    alpha <- mu1+mu2
    gamma1 <- (1-e1^2)^0.5
    gamma2 <- (1-e2^2)^0.5
    f <- function(delta){
        alpha^-3*(mu1+mu2/delta^2)*(mu1*gamma1+mu2*gamma2*delta)^2-1-3^(4/3)*mu1*mu2/alpha^(4/3)
    }
    d <- uniroot.all(f,c(0,10))
    return(d^2*a1)
}
show.peaks <- function(ps,powers,levels=NULL,Nmax=5){
    if(is.null(levels)) levels <- max(max(powers)-log(150),median(powers))
    ind <- which(powers==max(powers) | (powers>(max(powers)-log(100)) & powers>max(levels)))
    if(max(powers)-min(powers)<5) ind <- which.max(powers)
    pmax <- ps[ind]
    ppmax <- powers[ind]
    j0 <- 1
    p0 <- pmax[1]
    pp0 <- ppmax[1]
    pms <- p0
    pos <- pp0
    if(length(pmax)>1){
        for(j in 2:length(pmax)){
            if(abs(pmax[j]-p0) < 0.1*p0){
                if(ppmax[j]>pp0){
                    j0 <- j
                    p0 <- pmax[j]
                    pp0 <- ppmax[j0]
                    pms[length(pms)] <- p0
                    pos[length(pos)] <- pp0
                }
		    }else{
                j0 <- j
                p0 <- pmax[j]
                pp0 <- ppmax[j0]
                pms <- c(pms,p0)
                pos <- c(pos,pp0)
            }
        }
    }else{
        pms <- pmax
        pos <- ppmax
    }
    if(length(pms)>Nmax){
      pms <- pms[1:Nmax]
      pos <- pos[1:Nmax]
    }
    return(cbind(pms,pos))
}
