library(fields)
library(magicaxis)
#library(grDevices)
library(minpack.lm)
source('periodoframe.R')
#source('periodograms.R')
source('prepare_data.R',local=TRUE)
########################################
#####part I: making MP
########################################
###parameters setting
Nwindow <- 10
#colors <- c('black','red','blue','green','orange','brown','cyan','pink')
if(NI>0){
    Inds <- 1:NI
}else{
    Inds <- 0
}
scale <- FALSE
Dt <- 300#time span of the time window
Nbin <- 2#the steps required to cover the whole data baseline by moving the window
ofac0 <- 5#over sampling factor
per.type <- 'BFP'##type of moving periodogram: 'MLP' or 'BFP' or 'gls' or 'bgls'
#fmax <- 0.1#truncate to Pmin=10d
fmax <- 1/0.5#truncate to Pmin=10d
#fmin <- 1e-3
fmin <- 1/20
#fmax <- 1
t1 <- proc.time()
mp <- MP.norm(t=t,y=y,dy=dy,Dt=Dt,nbin=Nbin,ofac=ofac0,fmin=fmin,fmax=fmax,per.type=per.type,sj=0,Nma=Nma,Inds=Inds)
t2 <- proc.time()
dur <- format((t2-t1)[3],digit=3)
cat('computation time for MP:', dur,'s\n')
fname0 <- fname <- paste0(star,'_',per.type,'_Dt',Dt,'_Nbin',Nbin,'_ofac',ofac0)
dir <- 'results/'
#cat('save data into Robj file:\n')
#cat(paste0(fname,'.Robj'),'\n')
###save data; optional
#save(list = ls(all.names = TRUE),file=paste0(dir,fname,'.Robj'))
xx <- mp$tmid
yy <- mp$P
if(!scale){
    zz <- mp$powers
}else{
    zz <- mp$rel.powers
}
per.target <- ''
alpha <- 10
show.signal <- TRUE
#source('plot_2D.R',local=TRUE)
pdf.name <- paste0(dir,'MP_',fname,'.pdf')
cat('output pdf:\n',pdf.name,'\n')
pdf(pdf.name,6,6)
source('MP_plot.R',local=TRUE)
dev.off()
