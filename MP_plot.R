library(fields)
library(magicaxis)
if(scale) zz <- zz.rel
#th <- median(zz)-alpha*sd(zz)
#zlim <- c(th,max(zz))
zlim <- range(zz)
xlim <- range(t)
#zz[zz<th] <- th+1e-3
size <- 1.5
###plot
s <- 5
par(cex.axis=size,cex.lab=size,cex=size,oma=c(s,0,0,0),mar=c(0,0,0,0))
layout(matrix(1:4, ncol = 2,byrow=TRUE), widths = 1, heights = c(2,4), respect = FALSE)
if(exists('ypar')){
    ylab <- ypar
    if(grepl('RV',ypar)){
        ylab <- 'RV [m/s]'
    }
}else{
    ylab <- 'RV [m/s]'
}
###data plot
for(j in 1:2){
    if(j==1){
        par(mar = c(0, s, 1, 0))#,oma=c(1,1,1,4))
    }else{
        par(mar= c(0, 0, 1, s))#,oma=c(1,1,1,4))
    }
###plot
    Ntarget <- length(per.target)
    if(Ntarget>8){
        cols <- rainbow(Ntarget)
    }else{
        cols <- c('black','red','blue','green','orange','brown','cyan','pink')
    }
    if(exists('idata')){
        for(k in 1:Ntarget){
            t1 <- idata[[k]][,1]
            y1 <- idata[[k]][,2]
            dy1 <- idata[[k]][,3]
            if(k==1){
                if(j==1){
                    plot(t1%%2400000,y1,ylab=ylab,xaxt='n',pch=20,cex=0.5,xlim=xlim,ylim=range(y),col=cols[1])
                }else{
                    plot(t1%%2400000,y1,ylab=ylab,xaxt='n',yaxt='n',pch=20,cex=0.5,xlim=xlim,ylim=range(y),col=cols[1])
                }
            }else{
                points(t1%%2400000,y1,col=cols[k],pch=20,cex=0.5)
            }
            if(mean(dy1)>0.01*abs(mean(y1))){
                arrows(t1,y1-dy1,t1,y1+dy1,length=0.03,angle=90,code=3,col=cols[k])
            }
        }
    }else{
        if(j==1){
            plot(t,y,ylab=ylab,xaxt='n',pch=20,cex=0.5)
        }else{
            plot(t,y,ylab=ylab,xaxt='n',yaxt='n',pch=20,cex=0.5)
        }
        if(mean(dy)>0.01*sd(y)){
            arrows(t,y-dy,t,y+dy,length=0.03,angle=90,code=3)
        }
    }
}
#######periodograsm plot
cols <- rainbow(length(y),start=alpha/10)#
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
#sigs <- c(0.85,3.7,8.9)##define sigs manually
ind.show <- which.min(sigs)
for(j in 1:2){
    if(j==1){
        par(mar = c(1, s, 0, 0))
        image(xx,log10(yy),t(zz),xlab='',ylab='Period [day]',axes=FALSE,col=cols,xlim=xlim,zlim=zlim)
        at.labels <- axis(side=1,xpd=TRUE)
        magaxis(side=2,unlog=TRUE,tcl=-0.5)
#        Ntick <- round(log10(max(yy))-log(min(yy)))
#        ticks <- ceiling(log10(min(yy)))+(0:Ntick)*round(log10(max(yy))-log(min(yy)))/Ntick
#        p <- try(axis(side=2,at=ticks,labels=10^ticks),TRUE)
#        magaxis(side=2,unlog=TRUE)
#        if(class(p)=='try-error') 
        mtext(text='Time [JD-2400000]',side=1,outer=TRUE,line=2.2,cex=0.9*size)
#        mtext(at=c(min(t)-(max(t)-min(t))/10,0.5*(log10(max(yy))-log10(min(yy)))),text='Period [day]',side=2,outer=TRUE,line=2.2,cex=0.9*size)
    }else{
        par(mar = c(1, 0, 0, s))
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
            image(xx,log10(yy[inds]),t(zz[inds,]),xlab='',ylab='Orbital period [day]',axes=FALSE,col=cols,xlim=xlim,zlim=zlim)
            if(FALSE){
                labs <- at.labels
                labs[2*(1:floor(length(at.labels)/2))] <- ''
                labs[1] <- ''
                axis(side=1,at=at.labels,labels=labs,xpd=TRUE)
            }else{
                axis(side=1)
            }
            magaxis(side=2,unlog=TRUE,tcl=-0.5)
            if(FALSE){
                ticks <- seq(log10(pmin),log10(pmax),length.out=5)
                tick.lab <- c('',format(10^ticks[-c(1,length(ticks))],digit=2),'')
                axis(side=2,at=ticks,labels=tick.lab)
            }
        }
    }
    if(exists('sigs') & show.signal){
#        abline(h=log10(sigs),lty=3,lwd=2,col='grey')
        text(rep(min(t),length(sigs)),log10(sigs),labels=format(sigs,digit=3),pos=4,cex=size,col='grey',offset=-0.0)
        arrows(min(t)-0.1*(max(t)-min(t)),log10(sigs),min(t),log10(sigs),length=0.03,angle=90,code=3,col=cols[k])
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
image.plot(t(zz),col=cols,legend.only=TRUE,zlim = zlim,legend.mar=10,xpd=NA)
