crDistnShiny = function(R = 7,
                        M = 0.051,  #mortality
                        recruitRange = c(1, 100),   #just scales the density
                        mortRange = c(0.001, 0.5),
                        #CR parameters...
                        m = 0.4,
                        gamma = 0.053,
                        eta = 0.4,
                        mRange = c(-2, 0.78),       #should be less than one
                        gammaRange = c(0.001, 0.1),
                        etaRange = c(0,1),
                        #other...
                        d0 = 1e-8,                  #>0 for the closed-form density
                        units = c('English', 'metric'),
                        returnVar = NA,
                        useDT = TRUE,
                        maxBA = 200,
                        includeTiny = FALSE,
                        ...
                       )
{
#---------------------------------------------------------------------------
#
#   This shiny app is used to display an interactive session for the
#   Chapman-Richards distribution showing the CR growth model and the
#   associated distn based on the vital rates.
#
#   This is the main function that starts the shiny app and calls the
#   UI and Server functions. To use this, see the Returns section below;
#   essentially just running it like a normal function and allowing the
#   return value to "print" starts the app in the browser.
#
#   Note that this is in English units, if you want metric the only change
#   will be to basal area, everything else is agnostic to the units.
#**>I added metric to this for basal area, but note that the defaults
#   above give a very low starting BA for metric (since it is over the same
#   dbh range, now in cm rather than inches, making a smaller stand), so
#   one needs to play with the parameters; specifically increasing R and
#   decreasing mortality are a start to getting larger dbhes and BA.
#   JHG, 21-Jan-2016
#
#   Arguments...
#     R = recruitment numbers
#     M = mortality rate
#     recruitRange = the range for R slider
#     mortRange = the range for M slider
#     m,gamma,eta = CR growth model parameters
#     m-, gamma-, etaRange = the (min,max) slider ranges for these
#     d0 = minimum diameter in the CR distn as zero is illegal
#     units = 'English' or 'metric'
#     returnVar = variable name for the results of the app; NA=no results,
#                 otherwise a character string in .GlobalEnv; see the
#                 comments in qDisnShiny for more details.
#     useDT = TRUE: the display uses the DataTable format (very nice); FALSE:
#             uses a regular html table
#     maxBA = the maximum realistic value for a stand, for anything larger than
#             maxBA, a warning is displayed
#     includeTiny = TRUE: include all trees down to zero diameter; FALSE: only include
#                   small trees to the lower bound of the first dbh class
#     ... = passed on
#
#   Returns...
#     An object of class "shiny.appobj" that can be "printed" to run the app;
#     i.e., all that is needed is to run this function, the return gets
#     "printed" automatically unless it is captured in a variable, and so
#     the app gets run.
# 
#Author...									Date: 15-June-2015
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   this version requires Mensuration to get the class widths/midpoints...
#
    #require(Mensuration, quietly = TRUE)
    
#
#   shiny is required, and DT is optional...
#
    require("shiny", quietly = TRUE)
    if(useDT && !require("DT", quietly = TRUE)) {
      useDT = FALSE
      cat('\n***>DT packackage not avavilable: normal html tables used!\n')
    }
    if(d0 > 10e-4) {
      d0 = 10e-4
      cat('\nResetting d0 =', d0, '\n')
    }
  
#
#   units are used in ba calculations, where the default is English...
#
    units = match.arg(units)
    if(units=='metric') {
      k.ba = pi/40000                                   #kappa basal area conversion
      #d.conv = Units(1,'in','cm')[,3]                  #widget bounds conversions
      #d.conv = 2.54                                    #as aboveL inches to cm
      #ba.conv = Units(1,'ft^2/acre','m^2/hectare')[,3]  #not actually needed right now
      ba.conv = 0.22956841                             #as above: ft^2/acre to m^2/ha
    }
    else {
      k.ba = pi/(4*144)
      #d.conv = 1
      ba.conv = 1
    }
    conv = list(ba=ba.conv, k.ba=k.ba, units=units)



#
#   parameter vectors...
#
    params.cr = c(m=m, gamma=gamma, eta=eta)  #CR
    params.ssd = c(R=R, M=M)                  #ssd
    params.run = c(params.cr, params.ssd)     #all
    #params.run = c(params.run, maxDBH = (gamma/eta)^(1/(m-1)))


#
#
#   Chapman-Richards growth function in terms of dbh for SSD...
#
    growth = function(d, params=params.cr, ...) {
      y = params['eta']*d^params['m'] - params['gamma']*d
      return(invisible(y))
    }

#
#   Recruitment...
#
    recruitment = function(d, params, ...)
      return(params['R'])

#
#   mortality (harvest)...
#
    mortality = function(d, params, ...)
      return(params['M'])

#
#   the closed-form solution for the numbers density as a callable function;
#   note that this version substitutes tau into the equilibrium density, the
#   full form is in ncr below, but the two are algebraically equivalent...
#    
    ssd.cr = function(dbh, params=params.run) {
      m = params['m']
      gamma = params['gamma']
      eta = params['eta']
      R = params['R']
      M = params['M']
      ##mm1 = m-1
      #tau = -1/(gamma*(m-1))*log( (dbh^(m-1)*(eta*d0^(m-1) - gamma)) /
      #                            (d0^(m-1)*(eta*dbh^(m-1) - gamma))
      #                          )
      #N2.cr = R/growth(dbh, params)*exp(-M*tau)
      ee = M/(gamma*(m-1))
      N.cr = R/growth(dbh, params)*((dbh^(m-1)*(eta*d0^(m-1) - gamma)) /
                                    (d0^(m-1)*(eta*dbh^(m-1) - gamma)))^ee
      #if(!isTRUE(all.equal(N2.cr,N.cr))) print('woops!')
      return(N.cr)
    }


    cat('\nClick the "exit application" button to exit\n')

    if(!is.na(returnVar) && !is.character(returnVar))
      stop('"returnVar" must be character!')
	
    #set up the UI...
    ui = crDistnUI(params.run, mortRange, recruitRange, mRange, gammaRange, etaRange,
                   useDT, ...) #returns a list
    #and the server...
    server = crDistnServer(growth, recruitment, mortality, ssd.cr, params.run,
                           d0, returnVar, conv, useDT, maxBA, includeTiny, ...) #returns a closure
    
    #return the function...
    return( shinyApp(ui=ui, server=server) )
}   #crDistnShiny
