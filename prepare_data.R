library(shiny)
############################################
#####PartI: set parameters
############################################
#star <- 'ChallengeDataSet2'
star <- 'GJ699'
#star <- 'GJ699'
#star <- 'HD172555'
np <- 0
NI <- 0
Nma <- 0
Np <- 0
cummulative <- TRUE
Ncum <- 0
opt.type <- 'sl'#nl(nonlinear fitting all parameters; i.e. without using the formula ) or sl(semi-linear fitting)
#model.type <- 'auto'#manually setting the number of MA components and differential RVs
model.type <- 'man'#automatically determine the optimal noise model
tol <- 1e-16#numerical fitting precision tolerance
Nap <- 1#number of aperture
Nc <- 1#1
per.types <- c('BFP')#the periodograms to calculate
#per.types <- c('BFP','MLP','BGLS','GLST','GLS','LS')
per.type.seq <- per.types[1]
sequential <- FALSE
quantify <- TRUE
Np.max <- 5
###To select the numbers of differential RVs and MA components, it is better not to include calibration data.
if(model.type=='auto'){
    Nc <- 0
}
NI0 <- NI#the number of activity indices
####include all other noise proxies as indices
if(Nap>0){
    NI <- NI+Nap-1+Nc
    NI0 <- NI0+Nc
}
leg.pos <- 'topright'
MLP.type <- 'sub'#subtract the correlated noise from the data and then marginalize the non-correlated-noise parameters
#MLP.type <- 'assign'#fix the correlated noise parameters at their optimal values and then marginalize the other parameters
############################################
#####PartII: load data
############################################
#f <- paste0('data/',star,'_HARPS_TERRA.dat')
#if(!file.exists(f)){
#f <- paste0('data/',star,'_HARPS.dat')
#}
loading <- TRUE
Nw <- 1
if(loading){
#file <- 'keppure_priormt_12eppar1_Ndata221_quantifyTRUE_1per1_Pmin1Pmax3689_Nw1ABTRUE_Nap0_modedata_HD41248_HARPS_TERRA_1AP1_ervab0ap_nosim_1planet_ARMAp0q1_0momc0_Nsamp2000000_tem1_acc2.3_pretem0.311P25.6d_negLmax506'
#    file <- 'keppure_priormt_talklike_Ndata221_quantifyTRUE_1per1_Pmin0.5Pmax7378_Nw1ABTRUE_Nap0_modedata_HD41248_HARPS_TERRA_1AP1_ervab0ap_smallerrFALSE_3planet_ARMA01_0momc3.2_Nsamp10000000_tem1_acc0.12_pretem1P13.4d26.2d290.4d_negLmax474'
#    file <- 'keppure_priormt_talklike_Ndata805_quantifyTRUE_1per1_Pmin1Pmax21315_Nw6ABFALSE_Nap0_modedata_GJ699_TERRA_1AP1_ervab3ap_smallerrFALSE_2planet_ARMA02_0momc0_Nsamp8000000_tem1_acc7.1_pretem1P437d232.5d_negLmax1794'
    file <- 'keppure_priormt_talkpost_Ndata805_quantifyTRUE_1per1_Pmin0.1Pmax21315_Nw6ABFALSE_Nap0_modedata_GJ699_TERRA_1AP1_ervab3ap_smallerrFALSE_1planet_ARMA02_0momc0_Nsamp4000000_tem1_acc9.4_pretem1P231.6d_negLmax1816'
    load(paste0('../output/',star,'/',file,'.Robj'),envir=e <- new.env())
    Nw <- e$Nw
    if(Nw>1){
        res.all <- c()
        trv.all <- c()
        RV.all <- c()
        eRV.all <- c()
        par.data <- e$par.data
        data.names <- names(par.data)
        for(j in 1:Nw){
            if(!grepl('LICK',data.names[j])){
                res.all <- c(res.all,par.data[[j]]$RV2)
                RV.all <- c(RV.all,par.data[[j]]$RV)
                eRV.all <- c(eRV.all,par.data[[j]]$eRV)
                trv.all <- c(trv.all,par.data[[j]]$trv)
            }
        }
        inds <- sort(trv.all,index.return=TRUE)$ix
        trv.all <- trv.all[inds]
        RV.all <- RV.all[inds]
        eRV.all <- eRV.all[inds]
        res.all <- res.all[inds]
        sigs <- e$Popt
        tab <- cbind(trv.all,res.all,eRV.all)
    }else{
        tab <- cbind(e$trv,e$res,e$eRV)
    }
    Np <- e$Np
}else{
    file <- paste0('../data/aperture/',star,'/',star,'_HARPS.dat')
    tab <- read.table(file,header=TRUE,check.names=FALSE)
    Nproxy <- (ncol(tab)-3)/2
    tab <- cbind(tab[,1:3],tab[,2+(1:Nproxy)*2])
    if(cummulative){
        r <- c()
        for(j in 1:Nproxy){
            r <- c(r,abs(cor(tab[,2],tab[,3+j])))
        }
        ind <- sort(r,decreasing=TRUE,index.return=TRUE)$ix
        tab <- cbind(tab[,1:3],tab[,3+ind])
    }
}
if(ncol(tab)==3){
    NI <- 0
    Indices <- NA
    Inds <- 0
}
#NI0 <- NI <- ncol(tab)-3
############################################
#####PartIII: prepare data
############################################
########################
####set up
########################
#update <- FALSE
#opt.types <- c('full','part','wr','rw','rep')
ofac <- 1
fmax <- 1#Pmin=0.5d
mar.type <- 'part'
Ndata <- nrow(tab)
#model.type <- 'MA'#'MA' or 'auto'
cat('Ndata=',Ndata,'; Nma=',Nma,'; NI=',NI,'; Nap=',Nap,'; Nc=',Nc,'; opt.type=',opt.type,'; model.type=',model.type,'; ofac=',ofac,'; tol=',tol,'\n')

######rescale indices
Indices <- NA
sortInd <- FALSE
if(ncol(tab)>3){
    Indices <- as.matrix(tab[,4:ncol(tab)])
    if(!is.matrix(Indices) & !is.data.frame(Indices)){
        Indices <- matrix(Indices,ncol=1)
    }
    for(j in 1:ncol(Indices)){
        Indices[,j] <- scale(Indices[,j])
    }
}
t <- tab[,1]
y <- tab[,2]
dy <- tab[,3]
data <- cbind(t,y,dy)
