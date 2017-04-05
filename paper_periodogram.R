#f <- 'HD10700_TERRA_1AP1w0.dat'
#f <- 'HD020794_TERRA_1AP1_erv.dat'
#f <- 'HD010180_TERRA_1AP1_erv.dat'
#f <- 'HD222368_TERRA_1AP1_erv.dat'
#f <- 'HD222582_TERRA_1AP1_erv.dat'
#f <- 'WASP2_TERRA_1AP1_erv.dat'
#f <- 'HD115617_TERRA_1AP1_erv.dat'
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    np <- as.integer(args[1])
    NI <- as.integer(args[2])
    Nma <- as.integer(args[3])
    opt.type <- as.character(args[4])
    res.type <- as.logical(args[5])
    model.type <- as.character(args[6])
}else{
    np <- 2
    NI <- 3
    Nma <- 0
    opt.type <- 'sl'#l or nl or sl(semi-linear)
#    res.type <- 'kep'
    res.type <- 'cir'
   model.type <- 'man'
#   model.type <- 'auto'
#    if(res.type=='kep') fres <- 'keppure_priormt_talkpost_Ndata68_quantifyTRUE_1per1_Pmin1Pmax1684_Nw1ABTRUE_Nap3_modedata_HD177565_HARPS_TERRA_1AP1_ervab3ap_smallerrFALSE_1planet_ARMA_0momc3.2_Nsamp2000000_tem1_acc3.3_pretem0.566P44.5d_negLmax108'
    if(res.type=='kep') fres <- c('keppure_priormt_talklike_Ndata221_quantifyTRUE_1per1_Pmin0.5Pmax7378_Nw1ABTRUE_Nap0_modedata_HD41248_HARPS_TERRA_1AP1_ervab0ap_smallerrFALSE_2planet_ARMA01_0momc0_Nsamp2000000_tem1_acc1.1_pretem0.512P25.6d18.4d_negLmax493')
#fres <- c('keppure_priormt_talkpost_Ndata68_quantifyTRUE_1per1_Pmin1Pmax1684_Nw1ABTRUE_Nap3_modedata_HD177565_HARPS_TERRA_1AP1_ervab3ap_smallerrFALSE_1planet_ARMA_0momc3.2_Nsamp2000000_tem1_acc3.3_pretem0.566P44.5d_negLmax108','keppure_priormt_talklike_Ndata68_quantifyTRUE_1per1_Pmin0.5Pmax3369_Nw1ABTRUE_Nap3_modedata_HD177565_HARPS_TERRA_1AP1_ervab3ap_smallerrFALSE_2planet_ARMA01_0momc3.2_Nsamp4000000_tem1_acc0.89_pretem1P44.6d3.1d_negLmax100','keppure_priormt_talklike_Ndata68_quantifyTRUE_1per1_Pmin0.5Pmax3369_Nw1ABTRUE_Nap3_modedata_HD177565_HARPS_TERRA_1AP1_ervab3ap_smallerrFALSE_3planet_ARMA01_0momc3.2_Nsamp4000000_tem1_acc0.18_pretem1P44.3d3.1d1.4d_negLmax83.9')
}
tol <- 1e-12
Nap <- 1#minimum is 1
Nc <- 0
NI0 <- NI
if(Nap>1){
    NI <- NI+Nap-1+Nc
    NI0 <- NI0+Nc
}
#res.type <- 'kep'
#res.type <- 'cir'
fap <- FALSE
leg.pos <- 'topright'
MLP.type <- 'sub'#
MLP.type <- 'assign'
#example: Rscript compare_periodograms.R 0 3 0 part FALSE man 
if(np>0){
    f <- paste0('PlSystem',np,'_harps.rdb')
}else{
#####
#     f <- 'WASP2_TERRA_1AP1_erv.dat'
#    f <- 'HD20794_TERRA_1AP1_erv.dat'
#    f <- 'GJ628_TERRA_1AP1_erv.dat'
#    f <- 'HD26965_TERRA_1AP1_erv.dat'
#    f <- 'HD102365_TERRA_1AP1_erv.dat'
#     f <- 'HD10700_TERRA_1AP1_erv.dat'
#    f <- 'HD189567_TERRA_1AP1_erv.dat'
#    f <- 'HD102365v1_TERRA_1AP1_erv.dat'
#    f <- 'HD085512_TERRA_1AP1_erv.dat'
#    f <- 'HD102365_ucles_new.dat'
#     f <- 'HD73350_TERRA_1AP1_erv.dat'
#     f <- 'HD085512_TERRA_1AP1_erv.dat'
#     f <- 'HD202917_TERRA_1AP1_erv.dat'
#     f <- 'HD39060_TERRA_1AP1_erv.dat'
#    f <- 'HD038858_TERRA_1AP1_erv.dat'
#     f <- 'HD172167_TERRA_1AP1_erv.dat'
#     f <- 'HD172167_elodie.dat'
#     f <- 'HD38858_KECK.vels'
#     f <- 'HD20807_TERRA_1AP1_erv.dat'
#     f <- 'HD222368_TERRA_1AP1_erv.dat'
#    f <- 'HD20766_TERRA_1AP1_erv.dat'
#     f <- 'HD222368_TERRA_1AP1_erv.dat'
#     f <- 'HD150433_TERRA_1AP1_erv.dat'
#     f <- 'HD207129_TERRA_1AP1_erv.dat'
#     f <- 'Vega.dat'
#     f <- 'HD4628_TERRA_1AP1_erv.dat'
#     f <- 'HD33262_TERRA_1AP1_erv.dat'
#     f <- 'HD358623_TERRA_1AP1_erv.dat'
#     f <- 'HD172555_TERRA_1AP1_erv.dat'
#     f <- 'wbin0.5j0_HD020794_TERRA_1AP1_ervcon.dat'
#     f <- 'HD41248_TERRA_1AP1_erv.dat'
#     f <- 'HD177565_TERRA_1AP1_erv.dat'
}
f1 <- gsub('.+j0_','',f)
star.name <- gsub('_.+','',f1)
if(grepl('wbin',f)){
    name <- 'wbin0.5j0_'
}else{
    name <- ''
}
if(grepl('.rdb',f)){
    tab <- read.table(paste0('RV_challenge/',f),skip=2)
    tab <- tab[,c(1,2,3,4,6,10)]
    ##convert from km/s to m/s
    tab[,2:3] <- tab[,2:3]*1e3
}else if(grepl('HD172167_elodie',f)){
    tab <- read.table(paste0('data/',f))
    ind <- which(tab[,2]!=0 & tab[,2]>-14)
    tab <- tab[ind,]
    index <- which(diff(tab[,1])==0)
    tab <- tab[-index,]
    tab[,2] <- tab[,2]
    tab[,3] <- 0.015*tab[,3]/mean(tab[,3])#m/s
    cat('range(tab[,3])=',range(tab[,3]),'\n')
}else{
    fname <- paste0('data/',star.name,'/',f)
    cat('fname:',fname,'\n')
    if(!file.exists(fname)){
        fname <- paste0('../dwarfs/data/aperture/',star.name,'/',f)
        if(!file.exists(fname)) fname <- paste0('../dwarfs/data/aperture/',paste0(star.name,'_HARPS'),'/',f)
    }
    tab <- read.table(fname)
    if(model.type!='auto'){
        Naps <- Nap
#Nap!=0 & 
    }else if(grep('TERRA',f)){
        Naps <- c(3,6)
    }else{
        Naps <- 1
    }
    cat('Naps=',Naps,'\n')
    for(nap in Naps){
        if(nap>1){
            for(j in 2:nap){
                fdrv <- paste0('data/',star.name,'/',name,star.name,'_TERRA_',nap,'AP',j,'-',j-1,'.dat')
                cat('fdrv:',fdrv,'\n')
                if(!file.exists(fdrv)){
                    fdrv <- paste0('../dwarfs/data/aperture/',star.name,'/',star.name,'_TERRA_',nap,'AP',j,'-',j-1,'.dat')
                }
                if(!file.exists(fdrv)) fdrv <- paste0('../dwarfs/data/aperture/',paste0(star.name,'_HARPS'),'/',star.name,'_TERRA_',nap,'AP',j,'-',j-1,'.dat')
                tmp <- read.table(fdrv)
                tab <- cbind(tab,tmp[,2])
            }
        }
    }
    if(Nc>0){
        for(j in 1:Nc){
            if(Nc>1){
                name0 <- paste0(Nc+1,'AP',j+1,'-',j)
            }else{
                name0 <- '3AP2-1'
            }
            if(grepl('con',f)){
                data.type <- 'con'
            }else{
                data.type <- ''
            }
            fmom <- paste0('data/',star.name,'/',name,'calibration_',name0,'_out3',data.type,'.dat')
            cat('fmom:',fmom,'\n')
            if(!file.exists(fmom)){
                fmom <- paste0('/Users/phillippro/Documents/projects/dwarfs/data/calibration/combined/',star.name,'/calibration_',name0,'_out3',data.type,'.dat')
            }
            if(!file.exists(fmom)) fmom <- paste0('/Users/phillippro/Documents/projects/dwarfs/data/calibration/combined/',star.name,'_HARPS/calibration_',name0,'_out3',data.type,'.dat')
            tmp <- read.table(fmom)
            tab <- cbind(tab,tmp[,2])
        }
    }
}
###sort the data
ind.sort <- sort(tab[,1],index.return=TRUE)$ix
tab <- tab[ind.sort,]
###
source('periodoframe.R')#new version; old version is in periodoframe_v3.R
source('periodograms.R')
########################
####set up
########################
target <- gsub('_.+','',f)
#update <- FALSE
rescale <- FALSE#TRUE
#opt.types <- c('full','part','wr','rw','rep')
ofac <- 1
mar.type <- 'part'
Ndata <- nrow(tab)
#model.type <- 'MA'#'MA' or 'auto'
cat('file:',f,';Ndata=',Ndata,'; Nma=',Nma,'; NI=',NI,'; Nap=',Nap,'; Nc=',Nc,'; opt.type=',opt.type,'; model.type=',model.type,'; ofac=',ofac,'; res.type=',res.type,'; tol=',tol,'; MLP.type=',MLP.type,'\n')
#rescale <- FALSE
Ps <- list()
Ps$s1 <- c(0009.8916, 0023.3678, 0033.2757, 0112.4589, 0273.2000)
Ps$s2 <- c(0003.7700, 0005.7936, 0010.6370, 0020.1644, 0075.2835)
Ps$s3 <- c(0001.1188, 0017.0110, 0026.3000, 0048.7521, 0201.5000, 0595.9800, 2315.4400)
Ps$s4 <- c()
Ps$s5 <- c(0014.6632, 0026.2000, 0034.6548, 0173.1636, 0283.1000, 0616.3200)
Ps$s6 <- c()#c(0036.8000, 0072.4784, 0105.2000, 0310.0000, 0541.5740)
Ps$s7 <- c(0038.3160, 0072.4784, 0100.9920, 0303.8000, 0541.5740)
Ps$s8 <- c()
Ps$s9 <- c()
Ps$s10 <- c(0000.8200, 0056.6771, 0296.3000)
Ps$s11 <- c(0014.6632, 0034.6548, 0096.9333, 0283.1000, 3245.2000)
Ps$s12 <- c(0003.0780, 0014.6632, 0034.6548, 0096.9333, 0268.9450, 3407.4600)
Ps$s13 <- c()
Ps$s14 <- c()
Ps$s15 <- c(0000.8792, 0003.6260, 0009.4720)
Ks <- list()
Ks$s1 <- c(1.45, 1.67, 2.05, 0.38, 0.22)
Ks$s2 <- c(2.75, 0.27, 2.85, 0.34, 1.35)
Ks$s3 <- c(0.96, 3.68, 0.38, 5.14, 0.42, 1.91, 3.87)
Ks$s4 <- c()
Ks$s5 <- c(0.65, 0.44, 0.69, 0.59, 0.41, 0.55)
Ks$s6 <- c()
Ks$s7 <- c(0.34, 2.94, 0.32, 0.16, 2.38)
Ks$s8 <- c()
Ks$s9 <- c()
Ks$s10 <- c(0.67, 0.61, 0.16)
Ks$s11 <- c(0.58, 0.62, 0.54, 0.37, 1.54)
Ks$s12 <- c(0.48, 0.58, 0.62, 0.54, 0.36, 1.64)
Ks$s13 <- c()
Ks$s14 <- c()
Ks$s15 <- c(3.44, 5.85, 5.56)
##########new system
#Ps$s0 <- c(4.9, 17.9, 182.4, 252.9)
#Ps$s0 <- c(20.9,1.4,9.1)
#Ps$s0 <- c(20.9, 1.4, 9.1)
#Ps$s0 <- c(40.5, 8.9)
#Ps$s0 <- c(612.65, 163.34, 20.004, 49.361, 91.327)
#Ps$s0 <- bf$Popt[1:2]#
Ps$s0 <- c(5.4,2.9)
#Ps$s0 <- c(14.3, 33.6, 19.6)
#Ps$s0 <- c(58.43, 2993.8)
#Ps$s0 <- c(122)
#Ps$s0 <- c(12.3, 1.7,407.15)
#Ps$s0 <- c(42.4, 1.4)
#cols <- c(rep('red',length(Ps$s0)-1),'blue')
cols <- rep('red',length(Ps$s0))
ind <- which(cols=='blue')
if(length(ind)>0){
    Ks$s0 <- c((length(Ps$s0)-length(ind)):1,rep((length(Ps$s0)-length(ind)),length(ind)))
}else{
#    Ks$s0 <- length(Ps$s0):1
    Ks$s0 <- rep(1,length(Ps$s0))
}
######reshape indices
Indices <- NULL
sortInd <- FALSE
if(ncol(tab)>3){
    Indices <- as.matrix(tab[,4:ncol(tab)])
    if(!is.matrix(Indices) & !is.data.frame(Indices)){
        Indices <- matrix(Indices,ncol=1)
    }
    for(j in 1:ncol(Indices)){
        Indices[,j] <- scale(Indices[,j])
    }
    if(sortInd){
        cors <- c()
        y <- tab[,2]
        for(j in 1:ncol(Indices)){
            cors <- c(cors,abs(cor(Indices[,j],y)))
        }
        ind.sort <- sort(cors,decreasing=TRUE,index.return=TRUE)$ix
        Indices <- Indices[,ind.sort]
    }
}
fmax <- 0.5#Pmin=0.5d
#######make BFP 
t1 <- proc.time()
cat('fmax=',fmax,'\n')
bf <- BFP(tab[,1],tab[,2],tab[,3],Nma=Nma,NI=NI,Indices=Indices,ofac=ofac,opt.type=opt.type,mar.type=mar.type,model.type=model.type,update=update,fmax=fmax,rescale=rescale,tol=tol)
t2 <- proc.time()
dur1 <- format((t2-t1)[3],digit=3)
cat('BFP computation time:',dur1,'s\n\n')

#######make MLP 
t1 <- proc.time()
ml <- MLP(tab[,1],tab[,2],tab[,3],Nma=Nma,NI=NI,ofac=ofac,mar.type='part',model.type=model.type,update=update,fmax=fmax,rescale=rescale,opt.par=NULL,Indices=Indices,MLP.type=MLP.type)
t2 <- proc.time()
dur2 <- format((t2-t1)[3],digit=3)
cat('MLP computation time:',dur2,'s\n')
Nma <- ml$Nma
NI <- ml$NI
fname <- paste0('periodograms_Ndata_',Ndata,'modeltype',model.type,'_',target,'_res',res.type,'_NI',NI,'Nma',Nma,'_Nc',Nc,'_ofac',ofac,'_opt',opt.type,'_MLP',MLP.type,'_Nap',Nap,'_paper2')
if(!file.exists('results/')) system(paste0('mkdir results'))

fname <- paste0('results/',fname)

###save data
obj.name <- paste0(fname,'.Robj')
cat(obj.name,'\n')
save(list = ls(all.names = TRUE),file=obj.name)

###plot
source('plot_periodoframe.R') 
#source('presentationPer.R') 
