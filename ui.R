library(shiny)
shinyUI(fluidPage(
    includeCSS("my_style.css"),
#    titlePanel(HTML("<font color='blue'>Needle</font>")),
    titlePanel(HTML("Agatha")),
    h4(HTML("<em>Disentangling periodic signals from correlated noise in a periodogram framework</em>")),
    tabsetPanel(
        tabPanel("About Agatha",uiOutput('about')
                 ),               
        tabPanel("Choose File",
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
                         uiOutput("scatter.target"),
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
                         uiOutput('comp.target'),
                         sliderInput("Nma.max",'Maximum number of MA components',min=0,max=10,value=1,step=1), 
                         uiOutput('proxy.type'),
                         uiOutput('proxy.text'),
                         uiOutput('nI.basic'),
                         uiOutput('nI.max'),
                         helpText("The number of proxies is counted from the fourth column of the data."),
                         uiOutput('Nman'),
                         uiOutput('nI.man'),
                         uiOutput('nI.comp'),
                                        # uiOutput('nI.group'),
                                        # uiOutput('mlp'),
                         actionButton('compare', 'compare noise models')
                     ,width = 6),
                     mainPanel(                         
                         tags$style(type="text/css", "#file_progress { max-width: 200px; }"),
                         uiOutput('BFtab'),
                         uiOutput('optNoise'),
                         uiOutput('download.logBF.table')
                     ,width=6)
                 )
                 ),
####calculate 1D periodograms
        tabPanel("1D Periodogram",
                 pageWithSidebar(
                     headerPanel(''),
                     sidebarPanel(
                         uiOutput('per.target'),
                         helpText("If more than one data sets are selected, only the MLP can be calculated for the combined data. To calculate the periodogram, the data sets are combined after subtracting the best-fit noise components.."),
                         uiOutput('per.type'),
                         helpText("There could be errors in the calculation of the BFP if the data is small (e.g. less than 20 data points) or not well sampled. "),
                         uiOutput('nma'),
                         uiOutput('Inds'),
                         sliderInput("prange","Period range in base-10 log scale",min = -1,max = 5,value = c(0.1,3),step=0.1),
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
                         uiOutput('per.target2'),
                         helpText("If more than one data sets are selected, only the MLP-based 2D periodogram can be calculated for the combined data. To calculate the periodogram, the data sets are combined after subtracting the best-fitted noise components."),
                         uiOutput('per.type2'),
                         htmlOutput('text2D'),
                         tags$br(),
                         uiOutput('nma2'),
                         uiOutput('Inds2'),
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
                         uiOutput("Dt"),
                         uiOutput('textDt'),
                         uiOutput('prange2'),
                         br(),
                         uiOutput("Nbin"),
                         helpText("The above parameters are called 'calculating parameters', which are used for calculate the moving periodogram.","The following parameters are called 'visualization parameters', and are set to optimize the visulization of signals."),
                         uiOutput("alpha"),
                         uiOutput('zoom'),
                         checkboxInput('scale','Normalize power',value=TRUE),
                         checkboxInput('show.signal','Show significant signals',value=TRUE),
                         helpText("If you change the calculating parameters, click both 'calculate' and 'plot' to show the 2D periodogram.", "If you only change the visualization parameters, only click 'plot' to show the periodogram."),
                         actionButton('data.update', 'calculate'),
                         actionButton('plot2D', 'plot'),
                         uiOutput('download.per2D.plot'),
                         helpText("The users are encouraged to make their own plot of moving periodogram by downloading and using the relevant data. The first row is the centers of time windows. The first column is the periods, and the rest data is the matrix of periodogram powers."),
                         uiOutput('download.MP.data')
                     ,width=6),
                     mainPanel(
#                         plotOutput("per2", width = "750px", height = 400)
                         uiOutput("plot.2Dper"),
                         htmlOutput("color"),width=6
                     )
                 )
                 )
    )
)
        )
