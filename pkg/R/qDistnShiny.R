qDistnShiny = function(units = c('English', 'metric'),
                       returnVar = NA,
                       useDT = TRUE,
                       shiftZero = FALSE,
                       ...
                       )
{
#---------------------------------------------------------------------------
#
#   This is the function that drives the Shiny q distn application. The
#   UI and Server are defined in separate files and the calls below
#   illustrate how these can be used in conjunction with defining these--
#   the server especially--outside the calling environment with use of
#   environments.
#
#**>Note that I removed the dependency on the "units" program below as it
#   was really overkill here, and it will make it simpler for people to
#   use this routine in general without having to worry about installing
#   units. JHG, 5-Oct-2018.
#
#   Arguments...
#     units = 'English' or 'metric'
#     returnVar = variable name for the results of the app; NA=no results,
#                 otherwise a character string in .GlobalEnv; see the
#                 comments in qDisnShiny for more details.
#     conv = conversion factors for units 
#     useDT = TRUE: the display uses the DataTable format (very nice); FALSE: use
#             the default html table format. 
#     shiftZero = TRUE: shift the lower DBH limit down to zero diameter;
#                 FALSE: no shift, just include small trees to the lower
#                        bound of the first dbh class
#     ... = gobbled
#
#   Returns...
#     An object of class "shiny.appobj" that can be "printed" to run the app;
#     i.e., all that is needed is to run this function, the return gets
#     "printed" automatically unless it is captured in a variable, and so
#     the app gets run.
#
#   For more information on the DT package, see: http://rstudio.github.io/DT/
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
#   this version requires Mensuration to get the class widths/midpoints...
#
    #requireNamespace(Mensuration, quietly = TRUE)
    
#
#   shiny is required, and DT is optional...
#
    require("shiny", quietly = TRUE)
    if(useDT && !require("DT", quietly = TRUE)) {
      useDT = FALSE
      cat('\n***>DT packackage not avavilable: normal html tables used!')
    }

  
#
#   units are used in ba calculations and conversions for the slider widgets,
#   where the defaul is English...
#
    units = match.arg(units)
    if(units=='metric') {
      k.ba = pi/40000
      #d.conv = Units(1,'in','cm')[,3]                  #widget bounds conversions
      d.conv = 2.54                                    #as aboveL inches to cm
      #ba.conv = Units(1,'ft^2/acre','m^2/hectare')[,3] #widget bounds conversions
      ba.conv = 0.22956841                             #as above: ft^2/acre to m^2/ha
    }
    else {
      k.ba = pi/(4*144)
      d.conv = 1
      ba.conv = 1
    }
    conv = list(d=d.conv, ba=ba.conv, k.ba=k.ba, units=units)

    cat('\nClick the "exit application" button to exit\n')

    if(!is.na(returnVar) && !is.character(returnVar))
      stop('"returnVar" must be character!')

    #set up the UI...
    ui = qDistnUI(conv, useDT, ...) #returns a list
    #and the server...
    server = qDistnServer(returnVar, conv, useDT, shiftZero, ...) #returns a closure
    
    #return the function...
    return( shinyApp(ui=ui, server=server) )
} #qDistnShiny

    

