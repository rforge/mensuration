qDistnUI = function(conv,
                    useDT,
                    ...
                   )
{
#---------------------------------------------------------------------------
#
#   This builds the user interface for the qDistnShiny function.
#
#   Arguments...
#     conv = conversion factors for units 
#     useDT = TRUE: the display uses the DataTable format (very nice); FALSE: use
#             the default html table format. 
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
          column(6, offset=3,
                 tags$div(list(HTML('<center>'),titlePanel('Q Distribution...'),HTML('</center>')))
                 #titlePanel("Q Distribution")
                 ),  #column
          column(3,  
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
                 sliderInput("dbh.max", "Maximum DBH:", min = round(10*conv$d), max = round(50*conv$d),
                            value = 20, step=1, round=TRUE), 
                 sliderInput("q.value", "q value", min = 1, max = 2, value = 1.2),
                 sliderInput("bpa", "basal area...", min = round(30*conv$ba), max = round(200*conv$ba),
                             value = 80*conv$ba)
                 ), #wellPanel
                 verbatimTextOutput('stand')
                 ), #column,
          column(5, 
                 plotOutput("qPlot"),
                 plotOutput("basdPlot")
                ), #column
          column(4,
                   if(useDT)
                     DT::dataTableOutput('qView') #no wellPanel here
                   else
                     wellPanel( tableOutput("qView") )
               ) #column
        ) #fluidRow
    ) #fluidPage


    return(ui) 
} #qDistnUI

