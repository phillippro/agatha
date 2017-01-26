fs <- list.files('data')
indices <- list.files('../data/activity_indices',full.name=FALSE)
for(f in fs){
    tab <- read.table(paste0('data/',f))
    target <- gsub('_TERRA.dat','',f)
    star <- gsub('_.+','',f)
    cat('\ndim(tab)=',dim(tab),'\n')
    if(ncol(tab)==6){
        colnames(tab)=c('Time','RV','eRV','BIS','FWHM','S-index')
    }else if(ncol(tab)==18){
        colnames(tab)=c('Time','RV','eRV',paste0('3AP',2:3,'-',1:2),paste0('6AP',2:6,'-',1:5),paste0('9AP',2:9,'-',1:8))#harps
    }else if(ncol(tab)==19){
        colnames(tab)=c('Time','RV','eRV','C3AP2-1',paste0('3AP',2:3,'-',1:2),paste0('6AP',2:6,'-',1:5),paste0('9AP',2:9,'-',1:8))#harps
    }else if(ncol(tab)==20){
        colnames(tab)=c('Time','RV','eRV','Halpha','S-index',paste0('3AP',2:3,'-',1:2),paste0('6AP',2:6,'-',1:5),paste0('9AP',2:9,'-',1:8))#harps
    }else if(ncol(tab)==22){
        colnames(tab)=c('Time','RV','eRV','BIS','FWHM','S-index','C3AP2-1',paste0('3AP',2:3,'-',1:2),paste0('6AP',2:6,'-',1:5),paste0('9AP',2:9,'-',1:8))#harps
    }else if(ncol(tab)==21){
        colnames(tab)=c('Time','RV','eRV','BIS','FWHM','S-index',paste0('3AP',2:3,'-',1:2),paste0('6AP',2:6,'-',1:5),paste0('9AP',2:9,'-',1:8))#harps
    }else if(ncol(tab)==7 & grepl('KECK',star)){
        colnames(tab)=c('Time','RV','eRV','S-index','H-alpha','Photon Count','Integration Time')#new keck
    }else if(ncol(tab)==6 & grepl('KECK',star)){
        colnames(tab)=c('Time','RV','eRV','S-index','Photon Count','Integration Time')#old keck
    }else if(ncol(tab)==3){
        colnames(tab)=c('Time','RV','eRV')#other
    }else{
        colnames(tab)=paste0('column',1:ncol(tab))
    }
    cat('star=',star,'\n')
#    cat('indices=',indices,'\n')
    ind <- which(grepl(star,indices))
    if(length(ind)>0){
        cat('length(ind)=',length(ind),'\n')
        files <- paste0('../data/activity_indices/',indices[ind])
        k0 <- 1
        for(k in 1:length(files)){
            file <- files[k]
            cat('file=',file,'\n')
            data <- read.table(file)
            tmp <- c()
            for(j in 1:nrow(tab)){
                index <- which.min(abs(tab[j,1]%%2400000-data[,1]%%2400000))
                tmp <- c(tmp,data[index,2])
            }
            if(!(any(grepl('Halpha',colnames(tab))) & grepl('Halpha',file))){
                if(ncol(tab)>(2+k)){
                    tab <- cbind(tab[,1:(2+k0)],tmp,tab[(3+k):ncol(tab)])
                }else{
                    tab <- cbind(tab[,1:(2+k0)],tmp)
                }
                f2 <- gsub('.+/','',file)
                f2 <- gsub('.dew','',f2)
                f3 <- gsub('.+_','',f2)
                cat('f3=',f3,'\n')
                colnames(tab)[3+k0] <- f3
                k0 <- k0+1
                fname <- paste0('data/',f)
                cat('fname:\n',fname,'\n')
            }
        }
    }

    write.table(tab,file=paste0('data/',f),row.names=FALSE,quote=FALSE)
}
