####This file is an example for making PeriodoFrame: the computation part
library(magicaxis)
source('prepare_data.R',local=TRUE)
source('periodoframe.R')#for PeriodoFrame
source('periodograms.R')#for other periodograms
############################################
#####PartIV: make PeriodoFrame
############################################
lab <- FALSE
#inds <- 3
inds <- NULL
sig <- TRUE
RML <- TRUE
nshow <- 1
plot.type <- 'sep'# or 'comb'; how to make plots; separately or combine
fname <- paste0('periodograms_proxy',plot.proxy,'_',paste(per.types,collapse='0'),'_Ndata_',Ndata,'_',target,'_Nma',Nma,'_Inds',paste(Inds,collapse='.'),'_ofac',ofac,'_Np',Np.max)
if(!file.exists('results/')) system(paste0('mkdir results'))
fname <- paste0('results/',fname)
######plot
leg.pos <- 'topright'
pdf.name <- paste0(fname,'.pdf')
if(plot.type!='sep'){
    pdf(pdf.name,8,8)
    size <- 1.0
    par(mfrow=c(2,2),mar=c(4.5,4.5,1,1),cex.lab=size,cex.axis=size,cex=size)
    tit <- FALSE
}
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
        if(per.type=='BFP'){
            dP <- 0.1
            ind.man <- 1
            if(subtract.manual){
                if(jj==ind.man){
                    fmin0 <- fmin
                    fmax0 <- fmax
#                    fmin <- 1/132
#                    fmax <- 1/130
                    fmin <- 1/25.9
                    fmax <- 1/24.5
                    dP <- 0.001
                }else if(jj>ind.man){
                    fmin <- fmin0 
                    fmax <- fmax0
                }
            }
            per <- BFP(tab[,1],rr,tab[,3],Nma=Nma,Inds=Inds,Indices=Indices,ofac=ofac,opt.type=opt.type,model.type=model.type,fmin=fmin,fmax=fmax,progress=FALSE,quantify=quantify,dP=dP)
            ylab <- 'ln(BF)'
        }
        if(per.type=='MLP'){
            per <- MLP(tab[,1],rr,tab[,3],Nma=Nma,Inds=Inds,ofac=ofac,mar.type='part',model.type=model.type,fmin=fmin,fmax=fmax,opt.par=NULL,Indices=Indices,MLP.type=MLP.type)
            if(RML){
                ylab <- expression('log(ML/'*ML[max]*')')
            }else{
                ylab <- 'log(ML)'
            }
            ylab <- 'log(ML)'
        }
        if(per.type=='GLS'){
            per <- gls(tab[,1],rr,tab[,3],ofac=ofac,fmax=fmax)
#            ylim <- c(min(median(per$power),0),1.1*max(per$power))
            ylab <- 'Power'
        }
        if(per.type=='BGLS'){
            per <- bgls(tab[,1],rr,tab[,3],ofac=ofac,fmax=fmax)
            if(RML){
                ylab <- expression('log(ML/'*ML[max]*')')
            }else{
                ylab <- 'log(ML)'
            }
        }
        if(per.type=='GLST'){
            per <- glst(tab[,1],rr,tab[,3],ofac=ofac,fmax=fmax)
            ylab <- 'Power'
        }
        if(per.type=='LS'){
            per <- lsp(times=tab[,1]-min(tab[,1]),x=rr,ofac=ofac,from=NULL,to=fmax,alpha=c(0.1,0.01,0.001))
            ylab <- 'Lomb-Scargle Power'
        }
                                        #            if(sd(rr)<sd(per$res)) break()
        if(plot.type=='sep'){
            pname <- gsub('.pdf',paste0('_',jj,'sig.pdf'),pdf.name)
            cat('output pdf:\n',pname,'\n')
            pdf(pname,4,4)
            size <- 1.2
            par(mar=c(4.1,4.1,0.5,0.5),cex.lab=size,cex.axis=size,cex=size,mgp=c(2.5,1,0))
            tit <- FALSE
        }
        if(tit){
            main <- paste0(per.type,'; ',jj,' signal')
        }else{
            main <- ''
        }
        if(plot.proxy=='window'){
            inds <- which(per$P<360)
            per$P <- per$P[inds]
            per$power <- per$power[inds]
        }

        ymin <- min(per$power)
#        ylim <- c(ymin,max(max(per$power)+0.2*(max(per$power)-ymin),per$sig.level[which(!is.na(per$sig.level))]))
        ylim <- c(ymin,max(per$power)+0.2*(max(per$power)-ymin))
        if(per.type=='MLP' | per.type=='BFP'){
            plot(per$P,per$power,xlab='Period [day]',ylab=ylab,type='l',log='x',ylim=ylim,main=main,xaxt='n')
            magaxis(side=1,tcl=-0.5)
        }else{
            plot(per$P,per$power,xlab='Period [day]',ylab=ylab,type='l',log='x',ylim=ylim,main=main,xaxt='n')
            magaxis(side=1,tcl=-0.5)
        }
#        if(length(inds)==0){
#            inds <- which(per$power.opt>max(per$sig.level))
#            if(length(inds)>0) inds <- inds[1:min(nshow,length(inds))]
                                        #        }
#        abline(v=131,col='steelblue')
        if(plot.proxy!='window'){
            abline(h=per$sig.level,lty=2)
            pmaxs <- per$ps[1]
            power.max <- per$power.opt[1]
            offset <- c(0.08*(max(per$power)-ymin), 0.02*(max(per$power)-ymin))
            text(pmaxs,power.max+offset[1],pos=3,labels=format(pmaxs,digit=4),col='red',cex=1.0)
            arrows(pmaxs,power.max+offset[1],pmaxs,power.max+offset[2],col='red',length=0.05)
        }else{
        abline(v=24.7,col='red')
        abline(v=92,col='red')
        }
#        abline(v=per$ps[1],col='red',lty=3,lwd=2)
#        text(x=per$ps[1],y=per$power.opt[1],pos=c(2,4),labels=paste0(format(per$ps[1],digit=3),'d'),col='red')
#        else{
#            abline(v=per$ps[inds],col='red',lty=3,lwd=2)
#            text(x=per$ps[inds],y=per$power.opt[inds],pos=4,offset=-0.1,labels=paste0(format(per$ps[inds],digit=3),'d'),col='red')
#        }
if(plot.type=='sep'){dev.off()}
    }
}
################################################
#####PartV: save the data and plot periodograms
################################################
###save data
####output pdf
if(plot.type!='sep'){
    cat('output pdf:\n')
    cat(pdf.name,'\n')
    dev.off()
}

