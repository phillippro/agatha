####This file is an example for making PeriodoFrame: the computation part
source('prepare_data.R',local=TRUE)
source('periodoframe.R')#for PeriodoFrame
source('periodograms.R')#for other periodograms
############################################
#####PartIV: make PeriodoFrame
############################################
if(NI>0){
    Inds <- 1:NI
}else{
    Inds <- 0
}
lab <- FALSE
                                        #inds <- 1:2
inds <- NULL
sig <- TRUE
RML <- TRUE
nshow <- 1
fname <- paste0('periodograms_',paste(per.types,collapse='0'),'_Ndata_',Ndata,'modeltype',model.type,'_',star,'_NI',NI,'Nma',Nma,'_Nc',Nc,'_ofac',ofac,'_opt',opt.type,'_Nap',Nap,'_res',Np.max,loading)
fname <- paste0('results/',fname)
######plot
leg.pos <- 'topright'
pdf.name <- paste0(fname,'.pdf')
if(sequential){
    pdf(pdf.name,8,8)
    size <- 1.0
    par(mfrow=c(2,2),mar=c(4.5,4.5,1,1),cex.lab=size,cex.axis=size,cex=size)
    tit <- TRUE
}else{
    pdf(pdf.name,4,4)
    size <- 1.0
    par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),cex.lab=size,cex.axis=size,cex=size)
    tit <- FALSE
}
t <- tab[,1]-min(tab[,1])
y <- tab[,2]
dy <- tab[,3]
for(kk in 1:length(per.types)){
    per.type <- per.types[kk]
    if(per.type==per.type.seq){
        NN <- Np.max
    }else{
        NN <- 1
    }
    for(jj in 1:NN){
#        ofac <- 1+2*(jj-1)/NN
        cat('\nfind',jj,'signal using ',per.type,'!\n')
###make another BFP for data subtracted by the signal
        if(jj==1){
            rr <- y
        }else{
            if(is.matrix(per$res)){
                rr <- per$res[1,]
            }else{
                rr <- per$res
            }
        }
                                        #        cat('head(rr)=',head(rr),'\n')
        if(per.type=='BFP'){
            per <- BFP(tab[,1],rr,tab[,3],Nma=Nma,Inds=Inds,Indices=Indices,ofac=ofac,opt.type=opt.type,model.type=model.type,fmax=fmax,progress=FALSE,quantify=quantify)
            ylim <- range(min(median(per$power),0),1.1*max(per$power),5)
            ylab <- 'log(BF)'
        }
        if(per.type=='MLP'){
            per <- MLP(tab[,1],rr,tab[,3],Nma=Nma,Inds=Inds,ofac=ofac,mar.type='part',model.type=model.type,fmax=fmax,opt.par=NULL,Indices=Indices,MLP.type=MLP.type)
            ylim <- range(min(median(per$power),0),1.1*max(per$power))
            if(RML){
                ylab <- expression('log(ML/'*ML[max]*')')
            }else{
                ylab <- 'log(ML)'
            }
        }
        if(per.type=='GLS'){
            per <- gls(tab[,1],rr,tab[,3],ofac=ofac,fmax=fmax)
            ylim <- range(median(per$power),max(per$power)+0.1*(max(per$power)-min(per$power)),per$sig.level)
                                        #            ylim <- c(min(median(per$power),0),1.1*max(per$power))
            ylab <- 'Power'
        }
        if(per.type=='BGLS'){
            per <- bgls(tab[,1],rr,tab[,3],ofac=ofac,fmax=fmax)
            ylim <- range(median(per$power),max(per$power)+0.1*(max(per$power)-min(per$power)),per$sig.level)
            if(RML){
                ylab <- expression('log(ML/'*ML[max]*')')
            }else{
                ylab <- 'log(ML)'
            }
        }
        if(per.type=='GLST'){
            per <- glst(tab[,1],rr,tab[,3],ofac=ofac,fmax=fmax)
            ylim <- range(median(per$power),max(per$power)+0.1*(max(per$power)-min(per$power)),per$sig.level)
            ylab <- 'Power'
        }
        if(per.type=='LS'){
            per <- lsp(times=tab[,1]-min(tab[,1]),x=rr,ofac=ofac,from=NULL,to=fmax,alpha=c(0.1,0.01,0.001))
            ylim <- range(median(per$power),max(per$power)+0.1*(max(per$power)-min(per$power)),per$sig.level)
            ylab <- 'Power'
        }
                                        #            if(sd(rr)<sd(per$res)) break()
        main <- paste0(per.type,'; ',jj,' signal')
        if(per.type=='MLP' | per.type=='BFP'){
            plot(per$P,per$power,xlab='Period[d]',ylab=ylab,type='l',log='x',ylim=ylim,main=main)
        }else{
            plot(per$P,per$power,xlab='Period[d]',ylab=ylab,type='l',log='x',ylim=ylim,main=main)
        }
        if(length(inds)==0){
            inds <- which(per$power.opt>max(per$sig.level))
            if(length(inds)>0) inds <- inds[1:min(nshow,length(inds))]
        }
        abline(h=per$sig.level,lty=2)
        if(length(inds)>0){
            abline(v=per$ps[inds],col='red',lty=3,lwd=2)
            text(x=per$ps[inds],y=per$power.opt[inds],pos=c(2,4),labels=paste0(format(per$ps[inds],digit=3),'d'),col='red')
        }else{
            abline(v=per$ps[inds],col='red',lty=3,lwd=2)
            text(x=per$ps[inds],y=per$power.opt[inds],pos=4,offset=-0.1,labels=paste0(format(per$ps[inds],digit=3),'d'),col='red')
        }
    }
}
################################################
#####PartV: save the data and plot periodograms
################################################
###save data
####output pdf
cat('output pdf:\n')
cat(pdf.name,'\n')
dev.off()

