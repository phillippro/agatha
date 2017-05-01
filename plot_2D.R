########################################
#####part II: plot settings
########################################
#####setting whether to use relative power/scaled power or use absolute power for a given periodogram
colors <- c('red','green','black','blue','purple','orange','brown','cyan','pink')
t <- t%%50000
if(per.type=='BFP'){
#    type <- 'abs'#'abs' or 'rel'
    type <- 'abs'#'abs' or 'rel'
}else{
    type <- 'rel'#'abs' or 'rel'
}
FP <- FALSE#whether to show false positives 
###truncate the power to optimize the color dynamical range
truncate <- FALSE
#truncate <- TRUE
if(truncate){
    if(type=='abs'){
        alpha <- 10
    }else{
        alpha <- 5
    }
    rs <- 0.0
}else{
        alpha <- NA
        if(type=='abs'){
        rs <- 0.65
    }else{
        rs <- 0.4
}
}
fname <- paste0(fname0,'_',type)#alpha',alpha
####prepare data for plot
xx <- mp$tmid%%50000
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
par(cex.axis=size,cex.lab=size,cex=size,oma=c(5,0,0,0))
layout(matrix(1:4, ncol = 2,byrow=TRUE), widths = 1, heights = c(2,4), respect = FALSE)
###data plot
for(j in 1:Nwindow){
    if(j==1){
        par(mar = c(0, 6.1, 1, 0))
    }else{
        par(mar= c(0, 0, 1, 6.1))
    }
    if(j==1){
        yaxt <- 's'
    }else{
        yaxt <- 'n'
    }
###plot
    if(Nw>1){
        for(k in 1:Nw){
            ind.rm <- grep('LICK',data.names)
            if(!grepl('LICK',data.names[k])){
                trv <- par.data[[k]]$trv
                rv <- par.data[[k]]$RV2
                erv <- par.data[[k]]$eRV
                if(k==1){
                    plot(trv,rv,ylab='RV[m/s]',xaxt='n',yaxt=yaxt,pch=20,cex=0.5,col=colors[k],ylim=range(res.all),xlim=range(trv.all))
                }else{
                    points(trv,rv,col=colors[k],pch=20,cex=0.5)
                }
#                arrows(trv,rv-erv,trv,rv+erv,length=0.03,angle=90,code=3,col=colors[k])
            }
        }
    }else{
        plot(t,y,ylab='RV[m/s]',xaxt='n',yaxt=yaxt,pch=20,cex=0.5,col=colors[1])
    }
}
if(Nw>1){
    lname <- gsub('GJ699_|','',data.names[-ind.rm])
    lname <- gsub('_.+','',lname)
    lname[which(lname=='TERRA')] <- 'HARPS'
    lname[which(lname=='KECK')] <- 'HIRES'
    legend('topright',inset=c(-0.45,0),legend=lname,col=colors[(1:Nw)[-ind.rm]],pch=rep(20,Nw-1),xpd=NA)
}
#######periodogram plot
power1D <- c()
for(j in 1:length(yy)){
    if(type=='rel'){
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
###show signals; e.g. the signals in the HARPS RVs for HD41248
cols <- rainbow(length(y),start=rs)#
for(j in 1:Nwindow){
    if(j==1){
        par(mar = c(1.1, 6.1, 0, 0))
        image(xx,log10(yy),t(zz),xlab='',ylab='Period [day]',axes=FALSE,col=cols,xlim=c(min(t),max(t)),zlim=zlim)
        axis(side=1)
#        at.labels <- axis(side=1,xpd=TRUE)
#        magaxis(side=2)
        magaxis(side=2,unlog=TRUE,tcl=-0.5)
#        Ntick <- round(log10(max(yy))-log(min(yy)))
#        ticks <- ceiling(log10(min(yy)))+(0:Ntick)*round(log10(max(yy))-log(min(yy)))/Ntick
#        axis(side=2,at=ticks,labels=10^ticks)
        mtext(text='JD-2450000',side=1,outer=TRUE,line=2.2,cex=0.9*size)
    }else{
        par(mar = c(1.1, 0, 0, 6.1))
        Nsig <- length(sigs)
        flow <- 0.9
        fup <- 2
        for(k in ind.show){
#            inds <- which(yy>12 & yy<30)
            flow <- 0.7
            fup <- 1.3
            inds <- which(yy>max(1/fmax,flow*min(sigs)) & yy<min(1/fmin,fup*min(sigs)))
            image(xx,log10(yy[inds]),t(zz[inds,]),xlab='',ylab='Period [day]',axes=FALSE,col=cols,xlim=c(min(t),max(t)),zlim=zlim)
            axis(side=1)
            ticks <- seq(log10(flow*sigs[k]),log10(fup*sigs[k]),length.out=5)
            tick.lab <- c('',format(10^ticks[-c(1,length(ticks))],digit=2),'')
            axis(side=2,at=ticks,labels=tick.lab,col.axis='black',line=-0.6)
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
