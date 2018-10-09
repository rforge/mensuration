crDistnUI = function(params.run,
                     mortRange,
                     recruitRange,
                     mRange,
                     gammaRange,
                     etaRange,
                     useDT,
                     ...
                    )
{
#---------------------------------------------------------------------------
#
#   This builds the user interface for the crDistnShiny function.
#
#   Arguments...
#     params.run = vector of current CR growth and vital rate parameters
#     recruitRange = the range for R slider
#     mortRange = the range for M slider
#     m-, gamma-, etaRange = the slider ranges for these
#     useDT = TRUE: the display uses the DataTable format (very nice); FALSE:
#             uses a regular html table
#     ... = gobbled for now
#
#   Returns...
#     It returns an S3 object of class (shiny.tag.list,list)
#
#Author...	                                     Date: 8-June-2015
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jgove@fs.fed.us/jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#
    ui = fluidPage(
#       center the title in the center section; exit on the right...
        fluidRow(
          column(5, offset=3,
                 tags$div(list(HTML('<center>'),titlePanel('Chapman-Richards Distribution...'),HTML('</center>')))
                 #titlePanel("CR Distribution")
                 ),  #column
          column(4,  
                 #wellPanel(                                  #wellPanel adds more unneeded vertical space
                   actionButton('stopIt', 'Exit Application') #  )
               ) #column
          ), #fluidRow
#       widgets and plots...
        fluidRow(
          column(3,
              wellPanel(
                 sliderInput("dc.width", "DBH class width:", min = 1, max = 5,
                            value = 1, step=1, round=TRUE),
                 sliderInput("recruitment", "Recruitment", min = recruitRange[1],
                             max = recruitRange[2], value = params.run['R']),
                 sliderInput("mortality", "Mortality rate...", min = mortRange[1], max = mortRange[2],
                             value = params.run['M']),
                 sliderInput("eta", "eta", min = etaRange[1], max=etaRange[2], value = params.run['eta']),
                 sliderInput("m", "m", min = mRange[1], max=mRange[2], value = params.run['m']),
                 sliderInput("gamma", "gamma", min = gammaRange[1], max=gammaRange[2],
                            value = params.run['gamma'])
                 ), #wellPanel
                 verbatimTextOutput('stand')
                 ), #column,
          column(5,
                    plotOutput("crPlot"),   #CR function
                    plotOutput("ssdPlot"),  #numbers density
                    uiOutput('noteTiny')    #let user know if tiny trees are included
                ), #column
          column(4, 
                   tabsetPanel(
                     tabPanel("BASD", plotOutput("basdPlot")),
                     tabPanel("MortSurv",
                              plotOutput("pMortPlot"),
                              helpText('Probability of Survival: Solid line & blue dots', tags$br(),
                                       'Probability of Mortality: Bar chart & dashed line'
                                      )
                             ), #tabPanel
                     tabPanel("StandTable", 
                          if(useDT)
                            DT::dataTableOutput('crView') #no wellPanel here
                          else
                            wellPanel( tableOutput("crView") )
                              )  #tabPanel
                   ), #tabsetPanel
                 uiOutput('warnBA')   #basal area too large warning!                 
               ) #column
        ) #fluidRow
    ) #fluidPage


    return(ui) 
} #crDistnUI

