####This file is an example for making PeriodoFrame: the computation part
inds <- 1:2
leg.pos <- 'topright'
############################
####find additional signals
############################
for(jj in 2:Nsig.max){
    cat('\nfind',jj,'signal!\n')
    if(is.matrix(rv.ls$res)){
        rr <- rv.ls$res[1,]
    }else{
        rr <- rv.ls$res
    }
    if(per.type.seq=='BFP'){
        rv.ls <- BFP(tab[,1],rr,tab[,3],Nma=Nma,Inds=Inds,Indices=Indices,ofac=ofac,opt.type='sl',model.type='man',fmin=frange[1],fmax=frange[2],tol=tol)
        ylab <- 'log(BF)'
        name <- 'logBF'
    }
    if(per.type.seq=='MLP'){
        rv.ls <- MLP(tab[,1],rr,tab[,3],Nma=Nma,Inds=Inds,ofac=ofac,mar.type='part',model.type='man',fmin=frange[1],fmax=frange[2],opt.par=NULL,Indices=Indices,MLP.type=MLP.type)
        ylab <- 'log(ML)'
        name <- 'logML'
    }
    if(per.type.seq=='GLS'){
        rv.ls <- gls(tab[,1],rr,tab[,3],ofac=ofac,fmin=frange[1],fmax=frange[2])
        name <- ylab <- 'Power'
    }
    if(per.type.seq=='BGLS'){
        rv.ls <- bgls(tab[,1],rr,tab[,3],ofac=ofac,fmin=frange[1],fmax=frange[2])
        name <- 'logML'
        ylab <- 'log(ML)'
    }
    if(per.type.seq=='GLST'){
        rv.ls <- glst(tab[,1],rr,tab[,3],ofac=ofac,fmin=frange[1],fmax=frange[2])
        name <- ylab <- 'Power'
    }
    if(per.type.seq=='LS'){
        rv.ls <- lsp(times=tab[,1]-min(tab[,1]),x=rr,ofac=ofac,from=NULL,to=fmax,alpha=c(0.1,0.01,0.001))
        name <- ylab <- 'Power'
    }
    ylim <- c(min(rv.ls$power),max(rv.ls$power)+0.15*(max(rv.ls$power)-min(rv.ls$power)))
#    if(sd(rr)<sd(rv.ls$res)) break()
    if(per.type.seq=='BFP' | per.type.seq=='MLP'){
        if(per.type.seq=='BFP' & all(rv.ls$logBF.opt<log(150)))
            break()
        if(per.type.seq=='MLP' & all(rv.ls$logBF.opt<log(1000)))
            break()
    }else{
        if(max(rv.ls$power) < max(rv.ls$sig.level))
            break()
    }
####store data
    per.data <- cbind(per.data,rv.ls$power)
#    tit <- paste('Periodogram: BGLS; Target:',instrument,'; Observable',ypar)
    tit <- paste0(per.type.seq,';',instrument,';',ypar,';',jj,' signal')
    tits <- c(tits,tit)
    if(length(rv.ls$sig.level)<3){
        sig.levels <- cbind(sig.levels,c(rv.ls$sig.level,rep(NA,3-length(rv.ls$sig.level))))
    }else{
        sig.levels <- cbind(sig.levels,rv.ls$sig.level)
    }
    cnames <- c(cnames,paste0(per.type.seq,jj,'signal:',ypar,':',name))
}

