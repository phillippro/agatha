library(shiny)
shinyUI(fluidPage(
    includeCSS("my_style.css"),
#    titlePanel(HTML("<font color='blue'>Needle</font>")),
    titlePanel(HTML("Agatha")),
    h4(HTML("<em>Disentangling periodic signals from correlated noise in a periodogram framework</em>")),
    tabsetPanel(
        tabPanel("About Agatha",uiOutput('about')
                 ),               
        tabPanel("Upload File",
                 sidebarLayout(
                     sidebarPanel(
                                        # h4('Upload Type'),
                                        #  checkboxGroupInput('uptype','Upload Type',choices=c('Select from list','Upload files'),  selected=NULL),
                         radioButtons('uptype','Upload Type',
                                      choices=c('Select from the list'='list','Upload files'='upload'),  selected='list'),
#                         uiOutput('select.type'),
                         uiOutput('files'),
                         uiOutput('uptext'),
                         actionButton('show','upload and show data'),
                         uiOutput('download.data')
                                        # added interface for uploading data from
                                        # http://shiny.rstudio.com/gallery/file-upload.html
                                        #       tags$br()
                     ),
                     mainPanel(
                         uiOutput('tab')
                     )
                 )
                 ),
        tabPanel("Scatter Plot",
                 pageWithSidebar(
                     headerPanel(''),
                     sidebarPanel(
                         uiOutput("xs"),
                         uiOutput("ys"),
                         actionButton('scatter', 'show scatter plot'),
                         uiOutput('download.scatter.button')
                     ),
                     mainPanel(
                                        #plotOutput('scatter')
                         uiOutput("scatter")
                     )
                 )
                 ),
        tabPanel("Model Comparison",
                 pageWithSidebar(
                     headerPanel(''),
                     sidebarPanel(
######noise model comparison
                         sliderInput("Nma.max",'Maximum number of MA components',min=0,max=10,value=1,step=1), 
                         uiOutput('proxy.text'),
                         uiOutput('proxy.type'),
                         uiOutput('nI.basic'),
                         uiOutput('nI.max'),
                         helpText("The number of proxies is counted from the fourth column of the data."),
                                        #uiOutput('nI.max'),
                         uiOutput('nI.comp'),
                                        # uiOutput('nI.group'),
                                        # uiOutput('mlp'),
                         actionButton('compare', 'compare noise models')
                     ),
                     mainPanel(                         
                         tags$style(type="text/css", "#file_progress { max-width: 200px; }"),
                         uiOutput('BFtab'),
                         uiOutput('optNoise'),
                         uiOutput('download.logBF.table')
                     )
                 )
                 ),
####calculate 1D periodograms
        tabPanel("1D Periodogram",
                 pageWithSidebar(
                     headerPanel(''),
                     sidebarPanel(
                         helpText("There could be errors in the calculation of the BFP if the data is small (e.g. less than 20 data points) or not well sampled. "),
                         selectInput("per.type",'Periodogram type',
                                     choices=c('BFP','MLP','GLST','BGLS','GLS','LS'),selected="MLP",multiple=TRUE),
                         uiOutput('nma'),
                         uiOutput('proxy'),
                         uiOutput('Inds'),
                         sliderInput("frange","Range of frequency in base-10 log scale",min = -5,max = 1,value = c(-3,-0.2),step=0.1),
                         sliderInput("ofac", "Oversampling factor", min = 0, max = 30, value=1,step=0.2),
                                        # "Empty inputs" - they will be updated after the data is uploaded
                         helpText("If the BFP is selected, only 'RV' is available for the following observable selection."),
                         uiOutput("var"),
#                         helpText("If the BFP is selected, the periodograms are only calculated for RVs. The meaning of variables are as follows: 'all'--all observables, 'RVs'--RVs, 'Proxies'-- noise proxies,'Instrument:Variable'--individual observables"),
                                        #uiOutput('helpvar'),
                                        #uiOutput('tv'),
                         checkboxInput('sequence','Find additional signals sequentially',value=FALSE),
                         uiOutput('per.type.seq'),
                         uiOutput('Nsig.max'),
                         actionButton('plot1D', 'plot periodograms'),
                         radioButtons('down.type','Download plots',
                                      choices=c('All plots'='all','Individual plot'='individual'),  selected='individual'),
#                         uiOutput('select.type'),
                         uiOutput('plot.single'),
                         uiOutput('download.per1D.plot'),
                         helpText("The users are encouraged to make their own periodogram figures, particularly the BFP, by downloading and using the relevant data."),
                         uiOutput('download.per1D.data')
#                         br(),
                     ),
                     mainPanel(
                         uiOutput("plot.1Dper")
#                         uiOutput("mp")
                     )
                 )
                 ),
####calculate 2D periodograms
        tabPanel("2D Periodogram",
                 pageWithSidebar(
                     headerPanel(''),
                     sidebarPanel(
                         helpText("To make 2D periodograms, the time series should contain at least 100 data points over a time span beyond 100 time units."),
                         selectInput("per.type2",'Periodogram type',
                                     choices=c('BFP','MLP','GLST','BGLS','GLS','LS'),selected="MLP",multiple=FALSE),
                         uiOutput('nma2'),
                         uiOutput('proxy2'),
                         uiOutput('Inds2'),
                         sliderInput("frange2","Range of frequency in base-10 log scale",min = -4,max = 1,value = c(-3,-1),step=0.1),
                         sliderInput("ofac2", "Oversampling factor", min = 0, max = 20, value=1,step=0.2),
#                         selectInput("yvar", "Choose observables", 
#                                     choices  = 'RVs',
#                                     selected = 'RVs',multiple=FALSE),    
                         uiOutput('var2'),
#                         helpText("'all': the periodograms of all variables; 
#                             'RVs': periodograms of RVs; 
#                            'Indices': periodograms of Indices; 
#                             'Instrument:Variable': individual variables"),
                                        #uiOutput('helpvar'),
                                       #uiOutput('tv'),
                         helpText("The time window should be adjusted to guarantee the existence of a few data points in the time window for each moving step."),
                         uiOutput("Dt"),
                         uiOutput("Nbin"),
                         helpText("The above parameters are called 'calculating parameters', which are used for calculate the moving periodogram.","The following parameters are called 'visualization parameters', and are set to optimize the visulization of signals."),
                         uiOutput("alpha"),
                         uiOutput('zoom'),
                         checkboxInput('scale','Scaling power',value=FALSE),
                         checkboxInput('show.signal','Show significant signals',value=TRUE),
                         helpText("If you change the calculating parameters, click both 'calculate' and 'plot' to show the 2D periodogram.", "If you only change the visualization parameters, only click 'plot' to show the periodogram."),
                         actionButton('data.update', 'calculate'),
                         actionButton('plot2D', 'plot'),
                         uiOutput('download.per2D.plot'),
                         helpText("The users are encouraged to make their own plot of moving periodogram by downloading and using the relevant data. The first row is the centers of time windows. The first column is the periods, and the rest data is the matrix of periodogram powers."),
                         uiOutput('download.MP.data')
                     ),
                     mainPanel(
#                         plotOutput("per2", width = "750px", height = 400)
                        uiOutput("plot.2Dper")
                     )
                 )
                 )
    )
)
        )
