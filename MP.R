library(fields)
library(magicaxis)
#library(grDevices)
library(minpack.lm)
source('periodoframe.R')
#source('periodograms.R')
#source('prepare_data.R',local=TRUE)
########################################
#####part I: making MP
########################################
###parameters setting
Nwindow <- 2
#colors <- c('black','red','blue','green','orange','brown','cyan','pink')
if(NI>0){
    Inds <- 1:NI
}else{
    Inds <- 0
}
Dt <- 2000#time span of the time window
Nbin <- 10#the steps required to cover the whole data baseline by moving the window
ofac0 <- 10#over sampling factor
per.type <- 'BFP'##type of moving periodogram: 'MLP' or 'BFP' or 'gls' or 'bgls'
fmax <- 0.1#truncate to Pmin=10d
fmin <- 1e-3
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

source('plot_2D.R',local=TRUE)
