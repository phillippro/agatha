nrc2 <- function(Nvar){
    nrow <- ceiling(sqrt(Nvar))
    if(Nvar<nrow*(nrow-1)){
        ncol <- nrow-1
    }else{
        ncol <- nrow
    }
    return(c(nrow,ncol))
}

nrc <- function(Nvar){
    nrow <- ceiling(Nvar/2)
    if(Nvar<2){
        ncol <- 1
    }else{
        ncol <- 2
    }
    return(c(nrow,ncol))
}

tv.per <- function(targets,ofac,data){
    for(target in targets){
        tab <- data[[target]]
        commandArgs <- function(trailingOnly=TRUE) c('NA',1000,100,ofac,'bgls','res')
        source('time_varying_periodogram.R',local=TRUE)
    }
}

calc.1Dper <- function(Nmax.plots, vars,per.par,data){
    var <- names(per.par)
    for(k in 1:length(var)){
        assign(var[k],per.par[[var[k]]])
    }
    per.data <- c()
    tits <- c()
    pars <- list()
    kk <- 1
    MLP.type <- 'sub'
    for(j1 in 1:length(vars)){
        for(j2 in 1:length(per.type)){
            if(per.type[j2]=='MLP' | per.type[j2]=='BFP'){
                for(j3 in 1:length(Nmas)){
                    pars[[kk]] <- list(var=vars[j1],per.type=per.type[j2],Inds=Inds,Nma=Nmas[j3])
                    kk <- kk+1
                }
            }else{
                pars[[kk]] <- list(var=vars[j1],per.type=per.type[j2],Inds=0,Nma=0)
                kk <- kk+1
            }
        }
    }
    Nvar <- min(length(pars),Nmax.plots)      
    sig.levels <- c()
    cnames <- c()
    pers <- c()
    ylabs <- c()
#    lapply(1:Nvar, function(i){
    for(i in 1:Nvar){
        var <- pars[[i]]$var
        Nma <- as.integer(pars[[i]]$Nma)
        Inds <- as.integer(pars[[i]]$Inds)
        per.type <- pars[[i]]$per.type
        instrument <- gsub(':.+','',var)
        ypar <- gsub('.+:','',var)
        tab <- data[[instrument]]
        Indices <- NULL
        if(ncol(tab)>3){
            Indices <- as.matrix(tab[,4:ncol(tab)])
            if(!is.matrix(Indices) & !is.data.frame(Indices)){
                Indices <- matrix(Indices,ncol=1)
            }
            Indices <- as.matrix(Indices)
            for(j in 1:ncol(Indices)){
                Indices[,j] <- as.numeric(scale(Indices[,j]))
            }
        }
        if(ypar=='RV'){
            dy <- tab[,3]
        }else{
            dy <- rep(0.1,nrow(tab))
        } 
        cat('per.type=',per.type,'\n')
        if(ypar!='Window Function'){
            y <- tab[,ypar]
            if(per.type=='GLST'){
                rv.ls <- glst(t=tab[,1]-min(tab[,1]),y=y,err=dy,ofac=ofac,fmin=frange[1],fmax=frange[2])
                ylab <- 'Power'
                name <- 'power'
            }else if(per.type=='GLS'){
                rv.ls <- gls(t=tab[,1]-min(tab[,1]),y=y,err=dy,ofac=ofac,fmin=frange[1],fmax=frange[2])
                ylab <- 'Power'
                name <- 'power'
            }else if(per.type=='BGLS'){
                rv.ls <- bgls(t=tab[,1]-min(tab[,1]),y=y,err=dy,ofac=ofac,fmin=frange[1],fmax=frange[2])
                ylab <- 'log(ML)'
                name <- 'logML'
#                rv.ls[['sig.level']] <- max(rv.ls$power)-log(c(10,100,1000))
            }else if(per.type=='BFP'){
                rv.ls <- BFP(t=tab[,1]-min(tab[,1]),y=y,dy=dy,
                             Nma=Nma,Inds=Inds,model.type='man',Indices=Indices,
                             ofac=ofac,fmin=frange[1],fmax=frange[2],tol=1e-12)
                ylab <- 'log(BF)'
                name <- 'logBF'
            }else if(per.type=='MLP'){
                rv.ls <- MLP(t=tab[,1]-min(tab[,1]),y=y,dy=dy,Nma=Nma,Inds=Inds,Indices=Indices,
                             ofac=ofac,fmin=frange[1],fmax=frange[2],MLP.type=MLP.type)
                ylab <- 'log(ML)'
                name <- 'logML'
            }else if(per.type=='LS'){
                rv.ls <- lsp(times=tab[,1]-min(tab[,1]),x=y,ofac=ofac,from=frange[1],to=frange[2],alpha=c(0.1,0.01,0.001))
                ylab <- 'Power'
                name <- 'power'
            }
#            tit <- paste('Periodogram:',per.type,'; Target:',instrument,'; Observable',ypar)
            tit <- paste0(per.type,'; ',instrument,';', ypar,';1 signal')
            if(length(rv.ls$sig.level)<3){
                sig.levels <- cbind(sig.levels,c(rv.ls$sig.level,rep(NA,3-length(rv.ls$sig.level))))
            }else{
                sig.levels <- cbind(sig.levels,rv.ls$sig.level)
            }
#            if(per.type=='MLP' | per.type=='BFP'){
#                tit <- paste('proxy columns:',paste(Inds,collapse=','),'; q =',Nma,';',tit)
#            }
        }else{
            y <- rep(1,nrow(tab))
            rv.ls <- bgls(t=tab[,1]-min(tab[,1]),y=rep(2,nrow(tab)),err=rep(0.2,nrow(tab)), ofac=ofac)
            tit <- paste0('BGLS;',instrument,';',ypar,'; 1 signal')
            ylab <- 'log(ML)'
            name <- 'logML'
        }
        ylabs <- c(ylabs,ylab)
        tits <- c(tits,tit)
        pers <- c(pers,per.type)
###plot
#        plotname <- paste("plot", i, sep="")
        if(per.type=='BFP'){
            yy  <- rv.ls$logBF
        }else if(per.type=='MLP'){
            yy  <- rv.ls$logBF#-max(rv.ls$logBF)
        }else if(per.type=='BGLS'){
            yy  <- rv.ls$power#-max(rv.ls$power)
        }else{
            yy <- rv.ls$power
        }
        if(!is.null(per.data)){
            if(nrow(per.data)>length(yy)){
                rv.ls$P <- c(rv.ls$P,rv.ls$P[length(rv.ls$P)])
                yy <- c(yy,yy[length(yy)])
            }else if(nrow(per.data)<length(yy)){
                rv.ls$P <- rv.ls$P[-length(rv.ls$P)]
                yy <- yy[-length(yy)]
            }
        }
        if(i==1) per.data <- cbind(per.data,rv.ls$P)
        per.data <- cbind(per.data,rv.ls$power)
        inds <- (ncol(per.data)-1):ncol(per.data)
        if(i==1)  cnames <- c(cnames,'P')
        cnames <- c(cnames,paste0(pers[i],'1signal:',ypar,':',name))
        if(exists('per.type.seq')){
            if(per.type==per.type.seq){
                source('additional_signals.R',local=TRUE)
            }
        }
    }
    colnames(per.data) <- cnames
    return(list(per.data=per.data,tits=tits,pers=pers,levels=sig.levels,ylabs=ylabs))
}

per1D.plot <- function(per.data,tits,pers,levels,ylabs,download=FALSE,index=NULL){
    if(is.null(index)){
        par(mfrow=c(ceiling(Nmax.plots/2),2),mar=c(5,5,3,1),cex.lab=1.5,cex.axis=1.5,cex=1,cex.main=1.0)
    }
    if(download & is.null(index)){
        par(mfrow=c(2,2),mar=c(5,5,3,1),cex.lab=1.2,cex.main=0.8,cex.axis=1.2,cex=1)
    }
    P <- per.data[,1]
#    cat('colnames(per.data)=',colnames(per.data),'\n')
    if(!is.null(index)){
        inds <- index
        titles <- rep('',inds)
    }else{
        inds <- 1:(ncol(per.data)-1)
        titles <- tits
    }
    for(i in inds){ 
        power <- per.data[,i+1]
        ylab <- ylabs[i]#gsub('.+:','',colnames(per.data)[i+1])
        per.type <- gsub('[[:digit:]]signal:.+','',colnames(per.data)[i+1])
        f1 <- gsub('signal:.+','',colnames(per.data)[i+1])
        Nsig <- gsub('[A-Z]','',f1)
#        cat('Nsig=',Nsig,'\n')
#        cat('per.type=',per.type,'\n')
        if(per.type=='BFP'){
#            ymin <- max(median(power),0)
            ymin <- median(power)

            ylim <- c(ymin,max(power)+0.1*(max(power)-ymin))
        }else{
            ymin <- min(power)
            ylim <- c(min(power),max(power)+0.1*(max(power)-min(power)))
        }
        ylim <- c(min(power),max(power)+0.15*(max(power)-min(power)))
        plot(P,power,xlab='Period [day]',ylab=ylab,xaxt='n',log='x',type='l',main=titles[i], ylim=ylim)
        magaxis(side=1)
        abline(h=levels[,i],lty=2)
        p <- show.peaks(ps=P,powers=power,levels=levels[,i])
#        legend('topright',legend=paste0(Nsig,'signal'),bty='n')
        if(!is.matrix(p)){
            pmaxs <- p[1]
            power.max <- p[2]
        }else{
            pmaxs <- p[,1]
            power.max <- p[,2]           
            if(length(pmaxs)>4){
                pmaxs <- p[1:2,1]
                power.max <- p[1:2,2]           
            }
        }
        if(length(pmaxs)>0){
            par(xpd=TRUE)
            offset <- c(0.08*(max(power)-ymin), 0.02*(max(power)-ymin))
            pmaxs <- pmaxs[1]
            power.max <- power.max[1]
            text(pmaxs,power.max+offset[1],pos=3,labels=format(pmaxs,digit=3),col='red',cex=1.0)
            arrows(pmaxs,power.max+offset[1],pmaxs,power.max+offset[2],col='red',length=0.05)
            par(xpd=FALSE)
        }
    }
}

per2D.data <- function(vars,per.par,data){
    var <- names(per.par)
    for(k in 1:length(var)){
        assign(var[k],per.par[[var[k]]])
    }
    pars <- list()
    kk <- 1
    for(j1 in 1:length(vars)){
        for(j2 in 1:length(per.type)){
            if(per.type[j2]=='MLP' | per.type[j2]=='BFP'){
                for(j3 in 1:length(Nmas)){
                    pars[[kk]] <- list(var=vars[j1],per.type=per.type[j2],Inds=Inds,Nma=Nmas[j3])
                    kk <- kk+1
                }
            }else{
                pars[[kk]] <- list(var=vars[j1],per.type=per.type[j2],Inds=0,Nma=0)
                kk <- kk+1
            }
        }
    }
    i <- 1
    var <- pars[[i]]$var
    Nma <- as.integer(pars[[i]]$Nma)
    Inds <- as.integer(pars[[i]]$Inds)
    per.type <- pars[[i]]$per.type
    instrument <- gsub(':.+','',var)
    ypar <- gsub('.+:','',var)
    tab <- data[[instrument]]
    Indices <- NA
    if(ncol(tab)>3){
        Indices <- as.matrix(tab[,4:ncol(tab)])
        if(!is.matrix(Indices) & !is.data.frame(Indices)){
            Indices <- matrix(Indices,ncol=1)
        }
        Indices <- as.matrix(Indices)
        for(j in 1:ncol(Indices)){
            Indices[,j] <- as.numeric(scale(Indices[,j]))
        }
    }
    t <- tab[,1]%%2400000#min(tab[,1])
    y <- tab[,2]
    dy <- tab[,3]
    mp <- MP(t=t,y=y,dy=dy,Dt=Dt,nbin=Nbin,ofac=ofac,fmin=frange[1],fmax=frange[2],per.type=per.type,sj=0,Nma=Nmas,Inds=Inds,tol=1e-12,Indices=Indices)
    x2 <- mp$tmid
    y2 <- mp$P
    z2 <- mp$powers
    z2.rel <- mp$rel.powers
    return(list(t=t,y=y,dy=dy,xx=x2,yy=y2,zz=z2,zz.rel=z2.rel))
}
plotMP <- function(vals,pars){
    var <- names(pars)
    for(k in 1:length(var)){
        assign(var[k],pars[[var[k]]])
    }
    t <- vals$t
    y <- vals$y
    dy <- vals$dy
    xx <- vals$xx
    yy <- vals$yy
    zz <- vals$zz
    zz.rel <- vals$zz.rel
    source('MP_plot.R',local=TRUE)
}

calcBF <- function(data,Nbasic,proxy.type,Nma.max,groups=NULL,Nproxy=NULL){
    t <- data[,1]
    y <- data[,2]
    dy <- data[,3]
    NI.max <- ncol(data)-3
    NI.inds <- list(0)
    if(NI.max>0){
        NI.inds <- list()
        Nvary <- NI.max-Nbasic
        if(proxy.type=='cum' & Nvary>0 & Nproxy>0){
            NI.inds[[1]] <- Nbasic
            Indices <- data[,4:ncol(data)]
            cors <- c()
            for(j in 1:ncol(Indices)){
                cors <- c(cors,abs(cor(Indices[,j],y)))
            }
            inds <- sort(cors,decreasing=TRUE,index.return=TRUE)$ix
            if(Nproxy>Nbasic){
                for(j in 1:(Nproxy-Nbasic)){
                    NI.inds[[j+1]] <- inds[1:j]
                }
            }
        }else if(proxy.type=='group' & Nvary>0){
            NI.inds <- lapply(1:(length(groups)+1),function(i) NI.inds[[i]] <- list())
            if(Nbasic>0){
                NI.inds[[1]] <- 1:Nbasic
            }else{
                NI.inds[[1]] <- 0
            }
            groups <- sort(as.integer(groups))
            for(j in 1:length(groups)){  
                if(j==1){
                    NI.inds[[j+1]] <- 1:groups[j]
                }else{
                    if(Nbasic>0){
                        NI.inds[[j+1]] <- c(1:Nbasic,(groups[j-1]+1):groups[j])
                    }else{
                        NI.inds[[j+1]] <- (groups[j-1]+1):groups[j]
                    }
                }
            }
        }else{
            NI.inds <- list(list(Nbasic:NI.max))
        }
    }
    if(ncol(data)>3){
        out <- BFP.comp(t, y, dy, Nmas=0:Nma.max,NI.inds=NI.inds,Indices=data[4:ncol(data)])
    }else{
        out <- BFP.comp(t, y, dy, Nmas=0:Nma.max,NI.inds=0,Indices=NA)
    }
    if(!is.matrix(out$logBF)){
        out$logBF <- matrix(out$logBF,ncol=1)
    }
    return(list(Inds=NI.inds,Inds.opt=out$Inds,Nmas=0:Nma.max,Nma.opt=out$Nma,logBF=out$logBF,NI.opt=out$NI))
}

MCMC.panel <- function(){
    id <- 'HD020794_TERRA_1AP1_ervab6ap_ccf'
    Niter <- 1.0e3
    Nbin.per <- 1
    nbin.per <- 1
    tem <- 1
    inicov <- 1e-3
    Pini <- 200#day
    noise.model <- 'ARMA05'#noise.model: white, GP(R), ARMA, TJ(Ntj=1,noise vary with RHK or SA index, the third column of HARPS data), TJ(Ntj=3,vary with FWHM, BIS, RHK), TARMA, TGP, ARMATJ(ARMA+TJ), GPTJ,TJAR(model the RV contributed by index as a AR(p)-like model), ARMATJAR(ARMA+TJAR), PSID (previous-subsequent index dependent, this model is similar to TJAR but without time-varying/index-dependent jitter), ARMAPSID(ARMA+PSID)
    period.par <- 'logp'
    Ncores <- 1
    Np <- 1
    mode <- 'data'#data,sim
    Dtye <- 'D'#Dtype: DE:differential exclusing the target aperture, D: different including all aperture, N: no dependence on differential RV
    Nw <- 1#fit models to multiple wavelength data sets simultaneously
    prior.type <- 'mt'
    calibration <- 0#
    commandArgs <- function(trailingOnly=TRUE){
        cat('args=',c(id,Niter,Nbin.per,nbin.per,tem,inicov,Pini,noise.model,period.par,Ncores,Np,mode,Dtype,Nw,prior.type,calibration)
           ,'\n')
        c(id,Niter,Nbin.per,nbin.per,tem,inicov,Pini,noise.model,period.par,Ncores,Np,mode,Dtype,Nw,prior.type,calibration)
    } 
    source('../mcmc_red.R',local=TRUE)
    return(list(folder=folder,pdf=gsub('.+/','',pdf.name)))
}