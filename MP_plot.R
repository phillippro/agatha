library(fields)
library(magicaxis)
if(scale) zz <- zz.rel
th <- median(zz)-alpha*sd(zz)
zlim <- c(th,max(zz))
zz[zz<th] <- th+1e-3
size <- 1.5
###plot
par(cex.axis=size,cex.lab=size,cex=size,mar=c(5,5,1,1),oma=c(0,0,1,2))
layout(matrix(1:4, ncol = 2,byrow=TRUE), widths = 1, heights = c(2,4), respect = FALSE)
###data plot
for(j in 1:2){
    if(j==1){
        par(mar = c(0, 5.1, 1, 0))#,oma=c(1,1,1,4))
    }else{
        par(mar= c(0, 0, 1, 4.1))#,oma=c(1,1,1,4))
    }
###plot
    if(j==1){
        plot(t,y,ylab='RV[m/s]',xaxt='n',pch=20,cex=0.5)
    }else{
        plot(t,y,ylab='RV[m/s]',xaxt='n',yaxt='n',pch=20,cex=0.5)
    }
    arrows(t,y-dy,t,y+dy,length=0.03,angle=90,code=3)
}
#######periodograsm plot
cols <- rainbow(length(y),start=0)#
###show signals; e.g. the signals in the HARPS RVs for HD41248
##find strong signals
power1D <- c()
for(j in 1:length(yy)){
    if(scale){
        power1D <- c(power1D,mean(zz[j,]))
    }else{
        power1D <- c(power1D,max(zz[j,]))
    }
}
inds <- sort(power1D,decreasing=TRUE,index.return=TRUE)$ix
sigs <- yy[inds]
powers <- power1D[inds]
tmp <- show.peaks(ps=sigs,powers=powers,levels=median(powers))
sigs <- tmp[,1]
ind.show <- which.min(sigs)
#cat('range(zz)=',range(zz),'\n')
#cat('length(xx)=',length(xx),'\n')
#cat('length(yy)=',length(yy),'\n')
#cat('dim(zz)=',dim(zz),'\n')
for(j in 1:2){
    if(j==1){
        par(mar = c(5.1, 5.1, 0, 0.5))
        image(xx,log10(yy),t(zz),xlab='',ylab='Period [day]',axes=FALSE,col=cols,xlim=c(min(t),max(t)),zlim=zlim)
        axis(side=1)
        magaxis(side=2,labels=FALSE,unlog=TRUE,tcl=-0.5)
        Ntick <- round(log10(max(yy))-log(min(yy)))
        ticks <- ceiling(log10(min(yy)))+(0:Ntick)*round(log10(max(yy))-log(min(yy)))/Ntick
        axis(side=2,at=ticks,labels=10^ticks)
        mtext(text='Time [JD-2400000]',side=1,outer=TRUE,line=-2,cex=0.9*size)
    }else{
        par(mar = c(5.1, 0, 0, 4.1))
        Nsig <- length(sigs)
        flow <- 0.9
        fup <- 1.2
        for(k in ind.show){
#            inds <- which(yy>12 & yy<30)
            if(exists('pmin.zoom')){
                pmin <- pmin.zoom
                pmax <- pmax.zoom
            }else{
                pmin <- max(min(yy),flow*sigs[k])
                pmax <- min(fup*sigs[k],max(yy))
            }
            inds <- which(yy>pmin & yy<pmax)
            if(length(inds)<5){
                pmin <- max(min(yy),flow*sigs[k])
                pmax <- min(2*fup*sigs[k],max(yy))
                inds <- which(yy> pmin& yy<pmax)
            }
            if(length(inds)<5){
                pmin <- min(yy)
                pmax <- 3*min(yy)
                inds <- which(yy>pmin & yy<pmax)
            }
            image(xx,log10(yy[inds]),t(zz[inds,]),xlab='',ylab='Orbital period [day]',axes=FALSE,col=cols,xlim=c(min(t),max(t)),zlim=zlim)
            at.labels <- axis(side=1,labels=FALSE)
            labs <- at.labels
            labs[2*(1:floor(length(at.labels)/2))] <- ''
            labs[1] <- ''
#            cat('labs=',labs,'\n')
            axis(side=1,at=at.labels,labels=labs)
            ticks <- seq(log10(pmin),log10(pmax),length.out=5)
            tick.lab <- c('',format(10^ticks[-c(1,length(ticks))],digit=2),'')
            axis(side=2,at=ticks,labels=tick.lab)
        }
    }
    if(exists('sigs') & show.signal){
        abline(h=log10(sigs),lty=3,lwd=2,col='grey')
        text(rep(min(t),length(sigs)),log10(sigs),labels=format(sigs,digit=3),pos=4,cex=size,col='grey',offset=-0.0)
####these signals show
        FP <- FALSE
        if(FP){
            sig.false <- c(18.4,290.4)
            abline(h=log10(sig.false),lty=3,lwd=2,col='cyan')
            text(rep(min(t),length(sig.false)),log10(sig.false),labels=format(sig.false,digit=3),pos=4,cex=size,xpd=NA,col='cyan',offset=-0.0)
        }
    }
###make pretty axes
    axis(side=3,labels=FALSE)
}
image.plot(t(zz),col=cols,legend.only=TRUE,zlim = zlim,legend.mar=10.1,xpd=NA)
