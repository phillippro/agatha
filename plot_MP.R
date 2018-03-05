########################################
#####part II: plot settings
########################################
#####setting whether to use relative power/scaled power or use absolute power for a given periodogram
type <- 'rel'#or 'abs'
FP <- FALSE#whether to show false positives 
###truncate the power to optimize the color dynamical range
#truncate <- FALSE
truncate <- TRUE
if(truncate){
    if(type=='abs'){
        alpha <- 8
    }else{
        alpha <- 2
    }
    rs <- 0.0
}else{
    alpha <- NA
    rs <- 0
}
fname <- paste0(fname0,'_',type)#alpha',alpha
####prepare data for plot
xx <- mp$tmid
yy <- mp$P
if(type=='abs'){
    zz <- mp$powers
}else{
    zz <- mp$rel.powers
}
###truncate the low values to make the color corresponding to peaks clear
if(truncate){
    th <- median(zz)-alpha*sd(zz)
    zlim <- c(th,max(zz))
    zz[zz<th] <- th+1e-3
}else{
    zlim <- range(zz)
}
size <- 1.2

#########################
#####part III: MP plot
#########################
###plot
pdf.name <- paste0(dir,'MP_',fname,'.pdf')
pdf(pdf.name,6,6)
par(cex.axis=size,cex.lab=size,cex=size)
layout(matrix(1:4, ncol = 2,byrow=TRUE), widths = 1, heights = c(2,4), respect = FALSE)
###data plot
for(j in 1:2){
    if(j==1){
        par(mar = c(0, 5.1, 1, 0))
    }else{
        par(mar= c(0, 0, 1, 4.1))
    }
###plot
    if(j==1){
        plot(t,y,ylab='RV[m/s]',xaxt='n',pch=20,cex=0.5)
    }else{
        plot(t,y,ylab='RV[m/s]',xaxt='n',yaxt='n',pch=20,cex=0.5)
    }
}
#######periodograsm plot
cols <- rainbow(length(y),start=rs)#
###show signals; e.g. the signals in the HARPS RVs for HD41248
#sigs <- c(13.4, 25.6)
sigs <- as.numeric(Popt)
ind.show <- which.min(sigs)
for(j in 1:2){
    if(j==1){
        par(mar = c(5.1, 5.1, 0, 0))
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
        fup <- 2
        for(k in ind.show){
            inds <- which(yy>10 & yy<30)
            image(xx,log10(yy[inds]),t(zz[inds,]),xlab='',ylab='Orbital period [day]',axes=FALSE,col=cols,xlim=c(min(t),max(t)),zlim=zlim)
            axis(side=1)
            ticks <- seq(log10(flow*sigs[k]),log10(fup*sigs[k]),length.out=5)
            tick.lab <- c('',format(10^ticks[-c(1,length(ticks))],digit=2),'')
            axis(side=2,at=ticks,labels=tick.lab)
        }
    }
    if(exists('sigs')){
        abline(h=log10(sigs),lty=3,lwd=2,col='grey')
        text(rep(min(t),length(sigs)),log10(sigs),labels=format(sigs,digit=3),pos=4,cex=size,col='grey',offset=-0.0)
####these signals show
        if(FP){
            sig.false <- c(18.4,290.4)
            abline(h=log10(sig.false),lty=3,lwd=2,col='cyan')
            text(rep(min(t),length(sig.false)),log10(sig.false),labels=format(sig.false,digit=3),pos=4,cex=size,xpd=NA,col='cyan',offset=-0.0)
        }
    }
###make pretty axes
    axis(side=3,labels=FALSE)
}
image.plot(t(zz),col=cols,legend.only=TRUE,zlim = zlim,legend.mar=8.1)
####
dev.off()
cat('MP computation time:',dur,'s\n')
cat('output pdf:\n')
cat(pdf.name,'\n')
