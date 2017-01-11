library(shiny)
library(magicaxis)
#library(lomb)
#library(ggplot2)
# use the below options code if you wish to increase the file input limit, in this example file input limit is increased from 5MB to 9MB
# options(shiny.maxRequestSize = 9*1024^2)
source('periodoframe.R')
source("periodograms.R")
source("mcmc_func.R",local=TRUE)
source('functions.R',local=TRUE)
Nmax.plots <- 16
count0 <- 0
instruments <- c('HARPS','SOHPIE','HARPN','AAT','KECK')

shinyServer(function(input, output, session) {
####select from list
    output$about <- renderUI({HTML(paste("
<html>
<head>
<style>
p {
    text-indent: 25px;
}
</style>
</head>
<body>
<br />
<p>Agatha is the name of my wife's most favorite crime novelist, Agatha Christie. Similar to the investigations of various crimes in the detective novels, the Agatha algorithm is to find the weak signals embedded in the sea of correlated noise. 

This web app is based on the code in GitHub: <a href='https://github.com/phillippro/agatha'>https://github.com/phillippro/agatha</a>. If you use this web app in your work, please cite 'Feng F., Tuomi M., Jones H. R. A., 2017, Agatha: disentangling periodic signals from correlated noise in a periodogram framework, to be submitted'.</p>

<p>
Agatha is based on the Bayes factor periodogram (BFP) and the marginalized likelihood periodogram (MLP). The BFP is calculated by maximizing the likelihood of a combination of sinusoids, linear functions of time and noisy proxies, and the moving average model. The Bayes factor for a given frequency is derived from the maximum likelihood by approximating the Bayes factor using the Bayes Information Criterion (BIC). The MLP is calculated by marginalizing the likelihood over the amplitudes of sinusoids and the parameters in a linear function of time. Before calculating MLP, the best-fitted noise model is subtracted from the data. 
</p>

<p>
The BFP and MLP can be compared with the other periodograms, which are the Lomb-Scargle periodogram (LS), the generalized LS (GLS), the GLS with floating trend (GLST), the Bayesian GLS (BGLS). All periodograms can be computed for the sub-dataset within a moving time window to form 2D periodograms, which are also called 'moving periodograms'. Moving periodograms are used to check the consistency of signals in time. The user should adjust the 'visualization parameters' to optimize the visualization of signals in moving periodograms. 
</p>
</body>
</html>
"))})
    
  output$files <- renderUI({
    if(is.null(input$uptype)) return()
    if(input$uptype=='list'){
        selectizeInput('target','Select data files from the list',
                  choices=gsub('.dat','',gsub('_TERRA.+','',list.files('data',
                                    full.names=FALSE))),multiple=FALSE)
    }else if(input$uptype=='upload'){
        fileInput('file', 'Choose data file', multiple=FALSE)
#      selectizeInput('ins','Instrument',choices=c('HARPS','KECK','SOPHIE','AAT','HARPSN'),selected=NULL,multiple=FALSE)
    }
  })
  
    output$uptext <- renderUI({
        if(is.null(input$uptype)) return()
        if(input$uptype=='upload'){
            helpText("The file name should be 'star_instrument.fmt' where 'fmt' could be any plain text format. The file should not contain column names. The first three columns should be observation times, observables (interpreted as RVs here) and measurement uncertainties, while the other columns are noise proxies.")
        }
    })

 output$nI.max <- renderUI({
      if(is.null(input$proxy.type)) return()
      if(input$proxy.type=='cum'){
          selectInput("ni.max",'Maximum number of noise proxies',choices = 0:NI.max(),selected = min(3,NI.max()))
      }
  }) 

  NI.max <- reactive({
      if(is.null(data())) return()
      ncol(data()[[1]])-3
  })
  
    output$per.type.seq <- renderUI({
        if(is.null(input$sequence)) return()
        if(input$sequence) selectInput("per.type.seq",'Periodogram used to find additional signals',choices=input$per.type,selected=NULL,multiple=FALSE)
    })

 output$Nsig.max <- renderUI({
        if(is.null(input$sequence)) return()
        if(input$sequence) sliderInput("Nsig.max", "Maximum number of signals", min = 2, max = 10,value=3,step=1)
    })

  output$nma <- renderUI({
    if(is.null(input$per.type)) return()
    if(any(input$per.type=='MLP' | input$per.type=='BFP')){
        selectizeInput("Nmas",'Number of MA components',choices = 0:10,selected = 1,multiple=FALSE)
    }
  })

  output$nma2 <- renderUI({
    if(is.null(input$per.type)) return()
    if(any(input$per.type=='MLP' | input$per.type=='BFP')){
        selectInput("Nmas2",'Number of MA components',choices = 0:10,selected = 0,multiple=FALSE)
    }
  })
  
  output$proxy <- renderUI({
    if(is.null(input$per.type)) return()
    if((any(grepl('MLP',input$per.type)) | any(grepl('BFP',input$per.type))) & NI.max()>0){
      checkboxInput('proxy','Include noise proxies in the model',value=FALSE)
    }
  })

  output$proxy2 <- renderUI({
    if(is.null(input$per.type2)) return()
    if((any(grepl('MLP',input$per.type)) | any(grepl('BFP',input$per.type))) & NI.max()>0){
      checkboxInput('proxy2','Include noise proxies in the model',value=FALSE)
    }
  })

  output$Inds2 <- renderUI({
    if(is.null(input$proxy2)) return()
    if(input$proxy2){
        selectInput("Inds2",'Noise proxies to be included in the model',choices = 1:NI.max(),selected = 0,multiple=TRUE)
    }
  })
  
  output$Inds <- renderUI({
    if(is.null(input$proxy)) return()
    if(input$proxy){
        selectInput("Inds",'Noise proxies to be included in the model',choices = 1:NI.max(),selected = 0,multiple=TRUE)
    }
  })

  output$proxy.text <- renderUI({
      if(is.null(NI.max())) return()
      if(NI.max()>0){
    helpText("If 'cumulative' is selected, the noise proxies would be arranged in decreasing order of the Pearson correlation coefficients between proxies and RVs. Then proxies would be compared cumulatively from the basic number up to the maximum number of proxies. If 'group' is selected, the proxies would not be rearranged, and would be compared in groups, which are determined by the basic number of proxies and group division numbers. For example, if the basic number is 4 and division numbers are 6,11,19, models with proxies of {1-4}, {1-6}, {1-4, 7-11}, {1-4, 12-19} would be compared.")
      }
  })

  output$proxy.type <- renderUI({
      if(is.null(NI.max())) return()
      if(NI.max()>0){
          radioButtons("proxy.type",'Type of proxy comparison',c("Cumulative"="cum","Group"='group'))
      }else{
          output$warn <- renderText({'No indices available for comparison!'})
          verbatimTextOutput("warn")
      }
  }) 

  output$nI.basic <- renderUI({
      if(is.null(input$proxy.type)) return()
      if(NI.max()>0){
#input$proxy.type=='cum'
          selectInput("NI0",'Basic number of noise proxies',choices = 0:NI.max(),selected = 0)#NI.max())
      }
  }) 
  
  output$nI.comp <- renderUI({
    if(is.null(input$proxy.type) | is.null(input$NI0)) return()
    if(input$proxy.type=='group' & (NI.max()>as.integer(input$NI0))){
      selectInput("NI.group",'Group division numbers',
                  choices=(as.integer(input$NI0)+1):NI.max(),multiple = TRUE )
    }
  }) 
  
  target.list <- reactive({
    f1 <- gsub('_TERRA.+','',list.files('data',full.names=FALSE))
    gsub('.dat','',f1)
  })
  
  instr <- reactive({
      gsub('.+_','',target())
  })
  
  target <- reactive({
    if(is.null(input$target) & is.null(input$file)) return()
    if(input$uptype=='list'){
      return(input$target)
    }else if(input$uptype=='upload'){
        fname <- input$file$name
        f <- gsub('.dat','',fname)
        f <- gsub('_TERRA','',f)
        return(f)
    }
  })

  # added "session" because updateSelectInput requires it
  data <- eventReactive(input$show,{
    ins <- instr()
    #input$show
    if(input$uptype=='upload'){
      if(length(instr())>0){
        tmp <- list(NA)
        df <- rep(tmp,length(instr()))
        names(df) <- target()#paste0(instr())
        for(i in 1:length(instr())){
          data.path <- input$file$datapath
          tab <- read.table(data.path)
          inds <- sort(tab[,1],index.return=TRUE)$ix
          tab <- tab[inds,]
          if(any(diff(tab[,1])==0)){
              ind <- which(diff(tab[,1])==0)
              tab <- tab[-ind,]
          }
          df[[i]] <- tab
          if(ncol(df[[i]])==6 & ins[i]=='HARPS'){
            colnames(df[[i]])=c('Time','RV','eRV','BIS','FWHM','S-index')#harps
          }else if(ncol(df[[i]])==7 & ins[i]=='KECK'){
            colnames(df[[i]])=c('Time','RV','eRV','S-index','H-alpha','Photon Count','Integration Time')#new keck
          }else if(ncol(df[[i]])==6 & ins[i]=='KECK'){
            colnames(df[[i]])=c('Time','RV','eRV','S-index','Photon Count','Integration Time')#old keck
          }else if(ncol(df[[i]])==3){
            colnames(df[[i]])=c('Time','RV','eRV')#other
          }else{
            colnames(df[[i]])=c('Time','RV','eRV',paste0('proxies',1:(ncol(df[[i]])-3)))#other
          }
        }
      }
    }else if(input$uptype=='list'){
      if(is.null(input$target)) return()
      tmp <- list(NA)
      df <- rep(tmp,length(input$target))
      names(df) <- input$target
      for(i in 1:length(input$target)){
        target <- input$target[i]
        star <- gsub('_.+','',target)
        dir  <- 'data/'
        f0 <- paste0(dir,target,'_TERRA.dat')
        if(!file.exists(f0)){
            f0 <- paste0(dir,target,'.dat')
        }
        tab <- read.table(f0)
        inds <- sort(tab[,1],index.return=TRUE)$ix
        tab <- tab[inds,]
        if(any(diff(tab[,1])==0)){
            ind <- which(diff(tab[,1])==0)
            tab <- tab[-ind,]
        }
        df[[i]] <- tab
        if(ncol(df[[i]])==6){
            colnames(df[[i]])=c('Time','RV','eRV','BIS','FWHM','S-index')
        }else if(ncol(df[[i]])==18){
          colnames(df[[i]])=c('Time','RV','eRV',paste0('3AP',2:3,'-',1:2),paste0('6AP',2:6,'-',1:5),paste0('9AP',2:9,'-',1:8))#harps
        }else if(ncol(df[[i]])==19){
          colnames(df[[i]])=c('Time','RV','eRV','C3AP2-1',paste0('3AP',2:3,'-',1:2),paste0('6AP',2:6,'-',1:5),paste0('9AP',2:9,'-',1:8))#harps
      }else if(ncol(df[[i]])==22){
          colnames(df[[i]])=c('Time','RV','eRV','BIS','FWHM','S-index','C3AP2-1',paste0('3AP',2:3,'-',1:2),paste0('6AP',2:6,'-',1:5),paste0('9AP',2:9,'-',1:8))#harps
      }else if(ncol(df[[i]])==21){
          colnames(df[[i]])=c('Time','RV','eRV','BIS','FWHM','S-index',paste0('3AP',2:3,'-',1:2),paste0('6AP',2:6,'-',1:5),paste0('9AP',2:9,'-',1:8))#harps
      }else if(ncol(df[[i]])==7 & instr()[i]=='KECK'){
          colnames(df[[i]])=c('Time','RV','eRV','S-index','H-alpha','Photon Count','Integration Time')#new keck
        }else if(ncol(df[[i]])==6 & instr()[i]=='KECK'){
          colnames(df[[i]])=c('Time','RV','eRV','S-index','Photon Count','Integration Time')#old keck
        }else{
          colnames(df[[i]])=c('Time','RV','eRV')#other
        }
      }
    }
    return(df)
})

    output$data.out <- downloadHandler(
        filename = function() {
            f1 <- gsub(" ",'_',Sys.time())
            f2 <- gsub(":",'-',f1)
            paste('data-', f2, '.txt', sep='')
        },
        content = function(file) {
            write.table(data()[[1]], file,quote=FALSE,row.names=FALSE)#FALSE,col.names=FALSE
        }
    )


    output$download.data <- renderUI({
        if(is.null(input$uptype)) return()
        if(input$uptype=='list' & !is.null(data())){
            downloadButton('data.out', 'Download data')
        }
    })

  output$tab <- renderUI({
    input$show
    if(is.null(data())) return()
    if(input$show>0 & !is.null(data())){
    isolate({
       tabs <- lapply(1:length(instr()),function(i){
          if(input$uptype=='upload'){
            output[[paste0('f',input$target[i])]] <- renderDataTable(data()[[i]])
            tabPanel(target()[i],dataTableOutput(paste0('f',input$target[i])))
          }else{
            output[[paste0('f',input$target[i])]] <- renderDataTable(data()[[i]])
            tabPanel(target()[i],dataTableOutput(paste0('f', input$target[i])))
          }
        })
       do.call(tabsetPanel, tabs)
    })
    }
  })

  ns <- reactive({
    nam <- c()
    for(i in 1:length(instr())){
      names <- colnames(data()[[i]])
      names <- names[names!='Time' & names!='eRV']
      names <- c(names,'Window Function')
      nam <- c(nam,paste(names(data())[i],names,sep=':'))
    }
    return(nam)
  })
  
  ns.wt <- reactive({
      nam <- c()
      lab <- c()
      for(i in 1:length(data())){
        labs <- names <- colnames(data()[[i]])
        labs[grep('RV',names)] <- 'RV [m/s]'
        labs[grep('Time',names)] <- 'Time [JD-2400000]'
        labs[!grepl('RV',names) & !grepl('eRV',names) & !grepl('Time',names)] <- paste('Normalized',labs[!grepl('RV',names) & !grepl('eRV',names) & !grepl('Time',names)])
        nam <- c(nam,paste(names(data())[i],names,sep=':'))
        lab <- c(lab,paste(names(data())[i],labs,sep=':'))
      }
    return(list(name=nam,label=lab))
  })
  
  output$xs <- renderUI({
    if(is.null(data())){return()}
    selectizeInput("xs", "Choose x axis", 
                  choices  = ns.wt()$name,
                  selected = ns.wt()$name[grep('Time',ns.wt()$name)[1]],multiple=FALSE)
  })
  

    output$ys <- renderUI({
        if(is.null(data())){return()}
        selectizeInput("ys", "Choose y axis", 
                       choices  = ns.wt()$name,
                       selected = ns.wt()$name[grep('RV',ns.wt()$name)[1]],multiple=FALSE)
    })
   
    scatterInput <- function(){
        i <- 1
        instrument <- gsub(':.+','',input$xs[i])
        indx <- which(input$xs==ns.wt()$name)
        x <- gsub('.+:','',ns.wt()$label[indx])
        indy <- which(input$ys==ns.wt()$name)
        y <- gsub('.+:','',ns.wt()$label[indy])
        varx <- data()[[instrument]][,gsub('.+:','',input$xs[i])]
        vary <- data()[[instrument]][,gsub('.+:','',input$ys[i])]
        if(!grepl('RV|eRV|Time','',input$xs[i])) varx <- scale(varx)
        if(!grepl('RV|eRV|Time','',input$ys[i])) vary <- scale(vary)
        plot(varx,vary,xlab=x,ylab=y,pch=20,cex=0.5)
#        legend('topright',legend=paste0('r=',format(cor(varx,vary),digit=3)),bty='n',cex=1.5)
        if(grepl('RV',input$xs[i]) | grepl('RV',input$ys[i])){
            eRV <- data()[[instrument]][,grep('eRV',colnames(data()[[instrument]]))]
            if(grepl('RV',input$ys[i])){
                arrows(varx,vary-eRV,varx,vary+eRV,length=0.03,angle=90,code=3)
            }else{
                arrows(varx-eRV,vary,varx+eRV,vary,length=0.03,angle=90,code=3)
            }
        }
    }
  
    observeEvent(input$scatter,{
        output$sca <- renderPlot({
            par(mfrow=c(length(input$xs),1),cex=1,cex.axis=1.5,cex.lab=1.5,mar=c(5,5,1,1))
            scatterInput()
        })
    })

    observeEvent(input$scatter,{
        output$scatter <- renderUI({
            height <- 400*ceiling(length(input$xs)/2)
            plotOutput("sca", width = "400px", height = height)
        })
    })

    output$download.scatter <- downloadHandler(
        filename = function() {
            f1 <- gsub(" ",'_',Sys.time())
            f2 <- gsub(":",'-',f1)
            paste0("scatter_",f2,".pdf")
        },
        content = function(file) {
            pdf(file,4,4)
            par(mar=c(5,5,1,1))
            scatterInput()
            dev.off()
        })    
    
    observeEvent(input$scatter,{
        output$download.scatter.button <- renderUI({
            downloadButton('download.scatter', 'Download scatter plot')
        })
    })
    
    output$var <- renderUI({
        if(is.null(input$per.type)) return()
        if(!any(grepl('MLP',input$per.type)) & !any(grepl('BFP',input$per.type))){
#        selectInput("yvar", "Choose observables", choices  = c('all','RVs','Indices',ns()),selected = NULL,multiple=TRUE)         
            selectInput("yvar", "Choose observables", choices  = ns(),selected = ns()[grepl('RV',ns())],multiple=TRUE)         
        }else{
            selectInput("yvar", "Choose observables", 
                        choices  = ns()[grepl('RV',ns())],
                        selected = ns()[grepl('RV',ns())],multiple=FALSE)    
        }
    })

  output$var2 <- renderUI({
      if(is.null(input$per.type)) return()
      selectInput("yvar2", "Choose observables", 
                  choices  = ns()[grepl('RV',ns())],
                    selected = ns()[grepl('RV',ns())],multiple=FALSE)    
  })

  output$helpvar <- renderUI({
    if(is.null(input$yvar)) return()
    helpText("If the BFP is selected, only 'RV' is available for selection. The meaning of variables are as follows: 'all'--the periodograms of all variables, 
            'RVs'--periodograms of RVs, 'Indices'-- periodograms of Indices,
             'Instrument:Variable'--individual variables")
  })

  periodogram.var <- reactive({
    if(is.null(input$yvar)) return() 
    vars <- input$yvar[input$yvar!='all' & input$yvar!='RVs' & input$yvar!='Indices']
    return(unique(vars))
})

  periodogram.var2 <-  reactive({
    if(is.null(input$yvar2)) return() 
    vars <- c()
    if(any(input$yvar2!='Indices' & input$yvar2!='all' & input$yvar2!='RVs')){
        vars <- c(vars,input$yvar2[input$yvar2!='all' & input$yvar2!='RVs' & input$yvar2!='Indices'])
    }
    return(unique(vars))
  })

  per.par <- reactive({
      vals <- list(ofac=input$ofac,frange=10^input$frange,per.type=input$per.type)
      if(any(input$per.type=='MLP'|input$per.type=='BFP')){
          vals <- c(vals,Nmas=input$Nmas)
          if(!is.null(input$proxy)){
              if(input$proxy){
                  vals <- c(vals,Inds=list(input$Inds))
              }else{
                  vals <- c(vals,Inds=0)
              }
          }else{
              vals <- c(vals,Inds=0)
          }
      }else{
          vals <- c(vals,Nmas=0,Inds=0)
      }
      if(input$sequence){
          vals <- c(vals,Nmas=0,Inds=0,per.type.seq=input$per.type.seq, Nsig.max=as.integer(input$Nsig.max))
      }
      return(vals)
  })

  per.par2 <- reactive({
      vals <- list(ofac=input$ofac2,frange=10^input$frange2,per.type=input$per.type2)
      if(any(input$per.type=='MLP'|input$per.type=='BFP')){
          vals <- c(vals,Nmas=input$Nmas2)
          if(!is.null(input$proxy2)){
              if(input$proxy2){
                  vals <- c(vals,Inds=list(input$Inds2))
              }else{
                  vals <- c(vals,Inds=0)
              }
          }else{
              vals <- c(vals,Inds=0)
          }
      }else{
          vals <- c(vals,Nmas=0,Inds=0)
      }
      vals <- c(vals,Dt=as.integer(input$Dt),Nbin=as.integer(input$Nbin),alpha=as.integer(input$alpha),scale=input$scale,pmin.zoom=input$range.zoom[1],pmax.zoom=input$range.zoom[2],show.signal=input$show.signal)
      return(vals)
  })
  
  model.selection <- eventReactive(input$compare,{
      instrument <- gsub(':.+','',ns()[1])
      tab <- data()[[instrument]]
      if(!is.null(input$NI0)){
          Nbasic <- as.integer(input$NI0)
      }else{
          Nbasic <- 0
      }
      if(!is.null(input$proxy.type)){
          if(input$proxy.type=='group'){
              groups <- input$NI.group
              proxy.type <- 'group'
              ni <- NI.max()
          }else{
              groups <- NULL
              proxy.type <- 'cum'
              ni <- as.integer(input$ni.max)
          }
      }else{
          groups <- NULL
          proxy.type <- 'cum'
          ni <- 0
      }
      out <- calcBF(data=tab,Nbasic=Nbasic,
                    proxy.type=proxy.type,
                    Nma.max=as.integer(input$Nma.max),
                    groups=groups,Nproxy=ni)
      col.names <- c()
      for(j in 1:length(out$Nmas)){
          if(out$Nmas[j]==0){
              col.names <- c(col.names,'white noise')
          }else{
              col.names <- c(col.names,paste0('MA(',out$Nmas[j],')'))
          }
      }
      row.names <- c()
      for(j in 1:length(out$Inds)){
          if(all(out$Inds[[j]]==0)){
              row.names <- c(row.names,'no proxy')
          }else{
              row.names <- c(row.names,paste0('proxies: ',paste(out$Inds[[j]],collapse=',')))
          }
      }
      logBF <- data.frame(round(out$logBF,digit=1))
      colnames(logBF) <- col.names
      rownames(logBF) <- row.names
      cat('colnames(logBF)=',colnames(logBF),'\n')
      cat('rownames(logBF)=',rownames(logBF),'\n')
      logBF.download <- logBF
      colnames(logBF.download) <- paste0('MA',out$Nmas)
      rnames <- c()
      for(j in 1:nrow(logBF)){
          rnames <- c(rnames,paste0('proxy',paste(out$Inds[[j]],collapse='-')))
      }
      rownames(logBF.download) <- NULL#rnames
      return(list(logBF=logBF,out=out,logBF.download=logBF.download))
  })

  observeEvent(input$compare,{
      output$BFtab <- renderUI({
          output$table <- renderTable({model.selection()$logBF},digits=1,caption = "BIC-estimated Bayes factor",
                                      caption.placement = getOption("xtable.caption.placement", "top"),
                                      caption.width = getOption("xtable.caption.width", NULL))
          tableOutput('table')
      })
  })

  output$download.logBF <- downloadHandler(
      filename = function(){
          f1 <- gsub(" ",'_',Sys.time())
          f2 <- gsub(":",'-',f1)
          paste('logBF-', f2, '.txt', sep='')
      },
      content = function(file) {
          write.table(round(model.selection()$logBF.download,digit=1), file,quote=FALSE,row.names=FALSE)
      }
  )

  output$download.logBF.table <- renderUI({
      if(is.null(model.selection())) return()
#      downloadLink('download.logBF', 'Download the Bayes Factor table')
      downloadButton('download.logBF', 'Download the Bayes Factor table')
  })

  output$optNoise <- renderUI({
      if(is.null(input$compare)) return()
      if(input$compare>0){
          output$noise.opt <- renderText({
              Nma.opt <- model.selection()$out$Nma.opt
              Inds.opt <- model.selection()$out$Inds.opt
              if(Nma.opt==0){
                  t1 <- 'white noise'
              }else{
                  t1 <- paste0('MA(',Nma.opt,')')
              }
              if(all(Inds.opt==0)){
                  t2 <- 'Optimal proxies: no proxy'
              }else{
                  t2 <- paste0('Optimal proxies: ',paste(Inds.opt,collapse=','))
              }
              text1 <- paste0('Optimal noise model: ',t1)
              text2 <- t2
              HTML(paste(text1, text2, sep = '<br/>'))
          })
          htmlOutput('noise.opt')
      }
  })

  Nper <- eventReactive(input$plot1D,{
    if(is.null(input$per.type)) return()
    Nvar <- length(periodogram.var())
    Nplots <- 0
    pars <- per.par()
    Nplots <- Nplots+length(input$per.type)
    Nplots <- max(1,Nplots)*Nvar
    return(Nplots)
  })
    
  Nper2 <- eventReactive(input$plot2D,{
    if(is.null(input$per.type2)) return()
    Nvar <- length(periodogram.var2())
    Nplots <- 0
    pars <- per.par2()
    Nplots <- Nplots+length(input$per.type2)
    Nplots <- max(1,Nplots)*Nvar
    return(Nplots)
  })

  tvper <- reactive({
    if(is.null(data())) return()
    logic <- c()
    for(j in 1:length(data())){
      trv <- data()[[j]][,1]
      if(length(trv)>100 & (max(trv)-min(trv))>1000){
        logic <- c(logic,TRUE)
      }else{
        logic <- c(logic,FALSE)
      }  
    }
    return(logic)
  })

  output$Dt <- renderUI({
#      cat('input$plot2D=',input$plot2D,'\n')
      if(!is.null(data())){
          tmin <- min(data()[[1]][,1])
          tmax <- max(data()[[1]][,1])
          selectizeInput('Dt','Moving time window [time unit]',
                  choices=c(10,20,50,100,200,300,400,500,600,700,800,900,1000,1200,1500,2000,3000,4000,5000,6000,7000,8000,9000,10000),selected=1000,multiple=FALSE)#round((tmax-tmin)/150)*100
#          sliderInput("Dt", "Moving time window", min = 100, max = ,value=min(1000,round(tmax-tmin)),step=100)
      }
  })
  
  output$Nbin <- renderUI({
      if(!is.null(data())){
          selectizeInput('Nbin','Number of moving steps',
                  choices=c(5,10,20,50,100,200,500),selected=10,multiple=FALSE)
#          sliderInput("Nbin", "Number of moving steps", min = 5, max = 500,value=10)
      }
  })

  output$alpha <- renderUI({
      if(!is.null(data())){
          sliderInput('alpha','Truncate the periodogram power to the median minus X*(standard deviation)', min = 0, max = 30,value=10)
      }
  })

  output$zoom <- renderUI({
      sliderInput('range.zoom','Zoom-in period range', min = as.integer(1/10^(input$frange2[2])), max = min(100,as.integer(1/10^(input$frange2[1]))),value=c(max(10,1/10^(input$frange2[2])),min(30,1/10^(input$frange2[1]))),step=1)
  })
  
  per1D.data <- eventReactive(input$plot1D,{
#      reactive({
#function(){
          calc.1Dper(Nmax.plots, periodogram.var(),per.par(),data())
  })
    
    output$per1D.data <- downloadHandler(
        filename = function() {
            f1 <- gsub(" ",'_',Sys.time())
            f2 <- gsub(":",'-',f1)
            paste('periodogram1D-', f2, '.txt', sep='')
        },
        content = function(file) {
            tab <- per1D.data()$per.data
            write.table(tab, file,quote=FALSE,row.names=FALSE)#FALSE,col.names=FALSE
        }
    )

    output$download.per1D.data <- renderUI({
        if(is.null(per1D.data())) return()
        downloadButton('per1D.data', 'Download data of periodograms')
    })

    output$per1D.figure <- downloadHandler(
        filename = function() {
            f1 <- gsub(" ",'_',Sys.time())
            f2 <- gsub(":",'-',f1)
            paste('periodogram1D-', f2, '.pdf', sep='')
        },
      content = function(file) {
        pdf(file,8,8)
        per1D.plot(per1D.data()$per.data,per1D.data()$tits,per1D.data()$pers,per1D.data()$levels,ylabs=per1D.data()$ylabs,download=TRUE)
        dev.off()
      })

    output$plot.single <- renderUI({
        if(is.null(input$down.type) | is.null(per1D.data())) return()
        if(input$down.type=='individual'){
            selectizeInput('per1D.name','Select periodogram',
                  choices=per1D.data()$tits,multiple=FALSE)
        }
    })

    output$per1D.single <- downloadHandler(
        filename = function() {
            f1 <- gsub(" ",'_',Sys.time())
            f2 <- gsub(":",'-',f1)
            paste('periodogram1D-individual-', f2, '.pdf', sep='')
        },
      content = function(file) {
        pdf(file,4,4)
        par(mar=c(5,5,1,1))
        ind <- which(input$per1D.name==per1D.data()$tits)
        per1D.plot(per1D.data()$per.data,per1D.data()$tits,per1D.data()$pers,per1D.data()$levels,per1D.data()$ylabs,download=TRUE,index=ind)
        dev.off()
      })

    output$download.per1D.plot <- renderUI({
        if(is.null(per1D.data())) return()
        if(input$down.type=='all'){
            downloadButton('per1D.figure', 'Download periodograms')
        }else{
            downloadButton('per1D.single', 'Download periodograms')
        }
    })

    output$help.per1D <- renderUI({
        if(is.null(per1D.data())) return()
        helpText("The column names are 'P' and 'type:Observable:power', where 'name' is the periodogram type, 'P' is period, and 'power' is the periodogram power which could be logarithmic marginalized likelihood (logML; for MLP and BGLS) or Bayes factor (logBF; for BFP) or power (for other periodograms). ")
    })
    
    output$per <- renderPlot({
        if(is.null(per1D.data())) return()
        per1D.plot(per1D.data()$per.data,per1D.data()$tits,per1D.data()$pers,per1D.data()$levels,per1D.data()$ylabs)
    })

    output$plot.1Dper <- renderUI({
        plotOutput("per", width = "750px", height = 400*ceiling(Nmax.plots/2))
    })

  MP.data <- eventReactive(input$data.update,{
      per2D.data(periodogram.var2(),per.par2(),data())
  })

    output$MP.data <- downloadHandler(
        filename = function() {
            f1 <- gsub(" ",'_',Sys.time())
            f2 <- gsub(":",'-',f1)
            paste('periodogram2D-', f2, '.txt', sep='')
        },
        content = function(file) {
            tmp <- MP.data()
            if(nrow(tmp$zz)==length(tmp$xx)){
                tab <- cbind(tmp$xx,tmp$zz)
                tab <- t(rbind(c(NA,tmp$yy),tab))
            }else{
                tab <- rbind(tmp$xx,tmp$zz)
                tab <- cbind(c(NA,tmp$yy),tab)
            }
            write.table(tab, file,quote=FALSE,row.names=FALSE,col.names=FALSE)#FALSE,col.names=FALSE
        }
    )

    output$download.MP.data <- renderUI({
        if(is.null(MP.data())) return()
        downloadButton('MP.data', 'Download data of 2D periodogram')
    })

  observeEvent(input$plot2D,{
      output$per2 <- renderPlot({
          isolate({
              plotMP(MP.data(),per.par2())
          })
      })
  })

    observeEvent(input$plot2D,{
        output$plot.2Dper <- renderUI({
            plotOutput("per2", width = "600px", height = "600px")
        })
    })
    
    output$per2D.figure <- downloadHandler(
        filename = function() {
            f1 <- gsub(" ",'_',Sys.time())
            f2 <- gsub(":",'-',f1)
            paste('periodogram2D-', f2, '.pdf', sep='')
        },
        content = function(file) {
            pdf(file,8,8)
            plotMP(MP.data(),per.par2())
            dev.off()
        })
    
    output$download.per2D.plot <- renderUI({
        if(is.null(MP.data())) return()
        downloadButton('per2D.figure', 'Download 2D periodogram')
    })

#########MCMC
  output$nbin <- renderUI({
    if(is.null(input$Nbin.per)) return()
    isolate({
      sliderInput("nbin.per", "Sampling over certain period range", min = 1, max = as.integer(input$Nbin.per), 
                  value=1,step=1)
    })
  })
  output$noise.type <- renderUI({
    if(input$noise.model!='ARMA'){
      return()
    }else{
        sliderInput("p", "Number of AR components", min = 0, max = 20, 
                    value=0,step=1)
        sliderInput("q", "Number of MA components", min = 0, max = 20, 
                    value=1,step=1)
      }
    })
  
#####show the pdf
  mcmc.file <- reactive({
    if(is.null(input$run)) return()
    isolate({
      if(input$run>0){
        MCMCpanel()
      }
    })
  })

  output$mcmcpanel <- renderImage({
    if(is.null(mcmc.file()) | input$run<1) return()
    isolate({
      if(file.exists(paste0(mcmc.file()[1],mcmc.file()[2]))){
              normalizePath(file.path(paste0('../',mcmc.file()[1]),
                                            mcmc.file()[2]))
      }
    })
    })
})
