library(shiny)
############################################
#####PartI: set parameters
############################################
#star <- 'ChallengeDataSet10'
star <- 'HD41248'
#star <- 'HD172555'
np <- 0
NI <- 0
Nma <- 0
opt.type <- 'sl'#nl(nonlinear fitting all parameters; i.e. without using the formula ) or sl(semi-linear fitting)
#model.type <- 'auto'#manually setting the number of MA components and differential RVs
model.type <- 'man'#automatically determine the optimal noise model
tol <- 1e-12#numerical fitting precision tolerance 
Nap <- 1#number of aperture 
Nc <- 0#1
per.types <- c('BFP')#the periodograms to calculate
#per.types <- c('BFP','MLP','BGLS','GLST','GLS','LS')
sequential <- FALSE
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
f <- paste0('data/',star,'_HARPS_TERRA.dat')
if(!file.exists(f)){
f <- paste0('data/',star,'_HARPS.dat')
}
tab <- read.table(f,header=TRUE)

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
cat('file:',f,';Ndata=',Ndata,'; Nma=',Nma,'; NI=',NI,'; Nap=',Nap,'; Nc=',Nc,'; opt.type=',opt.type,'; model.type=',model.type,'; ofac=',ofac,'; tol=',tol,'\n')

######rescale indices
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
}
t <- tab[,1]
y <- tab[,2]
dy <- tab[,3]
data <- cbind(t,y,dy)
