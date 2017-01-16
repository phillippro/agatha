####This file is an example for making PeriodoFrame: the computation part
source('prepare_data.R',local=TRUE)
source('periodoframe.R')#for PeriodoFrame
source('periodograms.R')#for other periodograms
############################################
#####PartIV: make PeriodoFrame
############################################
#per.types <- c('BFP','MLP','BGLS','GLST','GLS','LS')
per.types <- c('MLP')
if(NI>0){
    Inds <- 1:NI
}else{
    Inds <- 0
}
inds <- 1:2
sequential <- TRUE
fname <- paste0('periodograms_',paste(per.types,collapse='0'),'_Ndata_',Ndata,'modeltype',model.type,'_',star,'_NI',NI,'Nma',Nma,'_Nc',Nc,'_ofac',ofac,'_opt',opt.type,'_Nap',Nap)
fname <- paste0('results/',fname)
######plot
leg.pos <- 'topright'
pdf.name <- paste0(fname,'.pdf')
pdf(pdf.name,8,8)
size <- 1.0
par(mfrow=c(2,2),mar=c(4.5,4.5,1,1),cex.lab=size,cex.axis=size,cex=size)
for(kk in 1:length(per.types)){
    per.type <- per.types[kk]
    cat('periodogram: ',per.type,'\n')
    t <- tab[,1]-min(tab[,1])
    y <- tab[,2]
    dy <- tab[,3]
    cat('\nfind 1 signal!\n')

####BFP
if(per.type=='BFP'){
    t1 <- proc.time()
    per <- BFP(tab[,1],tab[,2],tab[,3],Nma=Nma,Inds=Inds,Indices=Indices,ofac=ofac,opt.type=opt.type,model.type=model.type,fmax=fmax,tol=tol,progress=FALSE)
    t2 <- proc.time()
    dur1 <- format((t2-t1)[3],digit=3)
    cat('BFP computation time:',dur1,'s\n\n')
    ylim <- c(median(per$logBF),max(per$logBF))
    plot(per$P,per$logBF,xlab='Period[d]',ylab='log(BF)',type='l',log='x',ylim=ylim,main=paste0(per.type,'; 1 signal'))
    if(np==0){
        ps <- per$Popt[inds]
        ks <- 2*rev(inds)
        if(length(inds)>1){
            text(x=per$Popt[inds],y=per$logBF.opt[inds],pos=c(2,4),labels=paste0(format(per$Popt[inds],digit=3),'d'),col='red')
        }else{
            text(x=per$Popt[inds],y=per$logBF.opt[inds],pos=4,offset=-0.1,labels=paste0(format(per$Popt[inds],digit=3),'d'),col='red')
        }
    }
    abline(h=per$sig.level,lty=2)
    abline(v=per$ps[inds],col='red',lty=3,lwd=2)
        if(length(inds)>1){
            text(x=per$ps[inds],y=per$power.opt[inds],pos=c(2,4),labels=paste0(format(per$ps[inds],digit=3),'d'),col='red')
        }else{
            text(x=per$ps[inds],y=per$power.opt[inds],pos=4,offset=-0.1,labels=paste0(format(per$ps[inds],digit=3),'d'),col='red')
        }

    abline(h=per$sig.level,lty=2)
    legend(leg.pos,legend=c('BFP',paste0(dur1,'s')),bty='n')
}

####MLP 
if(per.type=='MLP'){
    t1 <- proc.time()
    per <- MLP(tab[,1],tab[,2],tab[,3],Nma=Nma,Inds=Inds,ofac=ofac,mar.type='part',model.type=model.type,fmax=fmax,opt.par=NULL,Indices=Indices,MLP.type=MLP.type)
    t2 <- proc.time()
    dur2 <- format((t2-t1)[3],digit=3)
    cat('MLP computation time:',dur2,'s\n')
    Nma <- per$Nma
    NI <- per$NI
    plot(per$P,per$logBF,xlab='Period[d]',ylab='log(ML)',type='l',log='x',main=paste0(per.type,'; 1 signal'))#expression('log(ML/M'*L[max]*')')
    abline(h=per$sig.level,lty=2)
    abline(v=per$ps[inds],col='red',lty=3,lwd=2)
        if(length(inds)>1){
            text(x=per$ps[inds],y=per$power.opt[inds],pos=c(2,4),labels=paste0(format(per$ps[inds],digit=3),'d'),col='red')
        }else{
            text(x=per$ps[inds],y=per$power.opt[inds],pos=4,offset=-0.1,labels=paste0(format(per$ps[inds],digit=3),'d'),col='red')
        }
    legend(leg.pos,legend=c('MLP',paste0(dur2,'s')),bty='n')
}

#####BGLS
if(per.type=='BGLS'){
    t1 <- proc.time()
    per <- bgls(tab[,1],tab[,2],tab[,3],ofac=ofac,fmax=fmax)
    t2 <- proc.time()
    dur3 <- format((t2-t1)[3],digit=3)
    plot(per$P,per$power,xlab='Period[d]',ylab='Relative probability',type='l',log='x',main=paste0(per.type,'; 1 signal'))
    abline(h=per$sig.level,lty=2)
    abline(v=per$ps[inds],col='red',lty=3,lwd=2)
        if(length(inds)>1){
            text(x=per$ps[inds],y=per$power.opt[inds],pos=c(2,4),labels=paste0(format(per$ps[inds],digit=3),'d'),col='red')
        }else{
            text(x=per$ps[inds],y=per$power.opt[inds],pos=4,offset=-0.1,labels=paste0(format(per$ps[inds],digit=3),'d'),col='red')
        }

    legend(leg.pos,legend=c('BGLS',paste0(dur3,'s')),bty='n')
}

#####GLS
if(per.type=='GLS'){
    t1 <- proc.time()
    per <- gls(tab[,1],tab[,2],tab[,3],ofac=ofac,fmax=fmax)
    t2 <- proc.time()
    plot(per$P,per$power,xlab='Period[d]',ylab='Power',type='l',log='x',main=paste0(per.type,'; 1 signal'))
    abline(h=per$sig.level,lty=2)
    abline(v=per$ps[inds],col='red',lty=3,lwd=2)
        if(length(inds)>1){
            text(x=per$ps[inds],y=per$power.opt[inds],pos=c(2,4),labels=paste0(format(per$ps[inds],digit=3),'d'),col='red')
        }else{
            text(x=per$ps[inds],y=per$power.opt[inds],pos=4,offset=-0.1,labels=paste0(format(per$ps[inds],digit=3),'d'),col='red')
        }

    legend(leg.pos,legend=c('GLS',paste0(dur3,'s')),bty='n')
}

####GLST
if(per.type=='GLST'){
    t1 <- proc.time()
    per <- glst(tab[,1],tab[,2],tab[,3],ofac=ofac,fmax=fmax)
    t2 <- proc.time()
    dur <- format((t2-t1)[3],digit=3)
    plot(per$P,per$power,xlab='Period[d]',ylab='Power',type='l',log='x',main=paste0(per.type,'; 1 signal'))
    abline(h=per$sig.level,lty=2)
    abline(v=per$ps[inds],col='red',lty=3,lwd=2)
        if(length(inds)>1){
            text(x=per$ps[inds],y=per$power.opt[inds],pos=c(2,4),labels=paste0(format(per$ps[inds],digit=3),'d'),col='red')
        }else{
            text(x=per$ps[inds],y=per$power.opt[inds],pos=4,offset=-0.1,labels=paste0(format(per$ps[inds],digit=3),'d'),col='red')
        }

    legend(leg.pos,legend=c('GLST',paste0(dur,'s')),bty='n')
}

if(per.type=='LS'){
    t1 <- proc.time()
    per <- lsp(times=tab[,1]-min(tab[,1]),x=y,ofac=ofac,from=NULL,to=fmax,alpha=c(0.1,0.01,0.001))
    t2 <- proc.time()
    dur <- format((t2-t1)[3],digit=3)
    plot(per$P,per$power,xlab='Period[d]',ylab='Power',type='l',log='x',main=paste0(per.type,'; 1 signal'))
    abline(h=per$sig.level,lty=2)
    abline(v=per$ps[inds],col='red',lty=3,lwd=2)
        if(length(inds)>1){
            text(x=per$ps[inds],y=per$power.opt[inds],pos=c(2,4),labels=paste0(format(per$ps[inds],digit=3),'d'),col='red')
        }else{
            text(x=per$ps[inds],y=per$power.opt[inds],pos=4,offset=-0.1,labels=paste0(format(per$ps[inds],digit=3),'d'),col='red')
        }

    legend(leg.pos,legend=c('GLST',paste0(dur,'s')),bty='n')
}
############################
####find additional signals
############################
if(sequential){
    for(jj in 1:10){
        cat('\nfind',jj+1,'signal!\n')
###make another BFP for data subtracted by the signal
        if(is.matrix(per$res)){
            rr <- per$res[1,]
        }else{
            rr <- per$res
        }
#        cat('head(rr)=',head(rr),'\n')
        if(per.type=='BFP'){
            per <- BFP(tab[,1],rr,tab[,3],Nma=Nma,Inds=Inds,Indices=Indices,ofac=ofac,opt.type=opt.type,model.type=model.type,fmax=fmax,tol=tol,progress=FALSE)
            ylim <- c(min(median(per$logBF),0),1.1*max(per$logBF))
            ylab <- 'log(BF)'
        }
        if(per.type=='MLP'){
            per <- MLP(tab[,1],rr,tab[,3],Nma=Nma,Inds=Inds,ofac=ofac,mar.type='part',model.type=model.type,fmax=fmax,opt.par=NULL,Indices=Indices,MLP.type=MLP.type)
            ylim <- c(min(median(per$logBF),0),1.1*max(per$logBF))
            ylab <- 'log(ML)'
        }
        if(per.type=='GLS'){
            per <- gls(tab[,1],rr,tab[,3],ofac=ofac,fmax=fmax)
            ylim <- c(median(per$power),max(per$power)+0.1*(max(per$power)-min(per$power)))
#            ylim <- c(min(median(per$power),0),1.1*max(per$power))
            ylab <- 'Power'
        }
        if(per.type=='BGLS'){
            per <- bgls(tab[,1],rr,tab[,3],ofac=ofac,fmax=fmax)
            ylim <- c(median(per$power),max(per$power)+0.1*(max(per$power)-min(per$power)))
            ylab <- 'log(ML)'
        }
        if(per.type=='GLST'){
            per <- glst(tab[,1],rr,tab[,3],ofac=ofac,fmax=fmax)
            ylim <- c(median(per$power),max(per$power)+0.1*(max(per$power)-min(per$power)))
            ylab <- 'Power'
        }
        if(per.type=='LS'){
            per <- lsp(times=tab[,1]-min(tab[,1]),x=rr,ofac=ofac,from=NULL,to=fmax,alpha=c(0.1,0.01,0.001))
            ylim <- c(median(per$power),max(per$power)+0.1*(max(per$power)-min(per$power)))
            ylab <- 'Power'
        }
        if(sd(rr)<sd(per$res)) break()
        if(per.type=='MLP' | per.type=='BFP'){
            plot(per$P,per$logBF,xlab='Period[d]',ylab=ylab,type='l',log='x',ylim=ylim,main=paste0(per.type,'; ',jj+1,' signal'))
        }else{
            plot(per$P,per$power,xlab='Period[d]',ylab=ylab,type='l',log='x',ylim=ylim,main=paste0(per.type,'; ',jj+1,' signal'))
        }
        ind <- which(per$power.opt>log(150))
        abline(h=per$sig.level,lty=2)
        if(length(ind)>0){
#            ind <- ind[1:min(length(inds),length(ind))]
            ind <- 1
            abline(v=per$ps[ind],col='red',lty=3,lwd=2)
        }
        if(length(ind)>1){
            text(x=per$ps[ind],y=per$power.opt[ind],pos=c(2,4),labels=paste0(format(per$ps[ind],digit=3),'d'),col='red')
        }else{
            text(x=per$ps[ind],y=per$power.opt[ind],pos=4,offset=-0.1,labels=paste0(format(per$ps[ind],digit=3),'d'),col='red')
        }

        if(per.type=='BFP' | per.type=='MLP'){
            if(per.type=='BFP' & all(per$logBF.opt<log(150)))
                break()
            if(per.type=='MLP' & all(per$logBF.opt<log(1000)))
                break()
        }else{
            if(max(per$power) < max(per$sig.level))
                break()
        }
    }
}
}
################################################
#####PartV: save the data and plot periodograms
################################################
###save data
#obj.name <- paste0(fname,'.Robj')
#cat(obj.name,'\n')
#save(list = ls(all.names = TRUE),file=obj.name)

####output pdf
cat('output pdf:\n')
cat(pdf.name,'\n')
dev.off()

