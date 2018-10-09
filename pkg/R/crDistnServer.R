crDistnServer = function(growth,         #function
                         recruitment,    #function
                         mortality,      #function
                         ssd.cr,         #function
                         params.run,
                         d0,
                         returnVar,
                         conv,
                         useDT,
                         maxBA,
                         includeTiny,
                         ...
                        )
{
#---------------------------------------------------------------------------
#
#   This builds the server interface function for the crDistnShiny function.
#
#   Arguments...
#     growth = the CR growth function
#     recruitment = the recruitment function
#     mortality = the mortality function
#     ssd.cr = the CR distribution function
#     params.run = vector of current CR growth and vital rate parameters
#     d0 = minimum diameter in the CR distn as zero is illegal
#     returnVar = variable name for the results of the app; NA=no results,
#                 otherwise a character string in .GlobalEnv
#     conv = conversion factors for units 
#     useDT = TRUE: the display uses the DataTable format (very nice); FALSE:
#             uses a regular html table
#     maxBA = the maximum realistic value for a stand, for anything larger than
#             maxBA, a warning is displayed
#     includeTiny = TRUE: include all trees down to zero diameter; FALSE: only include
#                   small trees to the lower bound of the first dbh class
#     ... = passed on to, e.g. dbhClassLimits
#
#   Returns...
#     It returns a closure. 
#
#**>Note: one can not return something along with the server function and expect
#         R to evaluate it--e.g., the df below--because the running of the
#         server function is controlled by runApp(), it is only built here, so
#         no variables defined within are available from the return of this
#         function. This is why the returnVar argument is supplied: it holds
#         the final result (when the exit button is clicked) as assigned
#         in .GlobalEnv; this was the only work-around I could come up with.
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
#   define the server function...
#
    server = function(input, output) {


#
#     qmsd^2 from the GB1 2nd raw moment...
#      
      qmsd.sq = function(params) {
        gamma = params[['gamma']]
        b = params[['dmax']]
        a = 1 - params[['m']]
        p = 1
        q = params[['M']]/(gamma*a)
        qmsd = b*b*beta(p + 2/a, q)/beta(p, q)
        return(qmsd)
      } #qmsd.sq
             
#
#     this is the GB1 normalizing factor, which equals total N...
#    
      N.c = function(params) {
        m = params[['m']]
        gamma = params[['gamma']]
        eta = params[['eta']]
        R = params[['R']]
        M = params[['M']]
        b = params[['dmax']]
        a = 1 - m
        p = 1
        q = M/(gamma*a)
        c0 = R*( eta - gamma/(d0^(m-1)) )^(M/(gamma*(m-1)))
        c1 = eta^(M/(gamma*a)-1)
        N = c0 * c1 * b^(a*p) * beta(p,q)/abs(a)
        return(N)
      } #N.c
  

#
#     this function is the GB1 transformation of the CR numbers density...
#  
      gb1.cr = function(dbh, params=params.run, alpha=0) {
        m = params[['m']]
        gamma = params[['gamma']]
        M = params[['M']]
        b = params[['dmax']]
        a = 1 - m
        p = 1 + alpha/a
        q = M/(gamma*a)
        d.gb = (abs(a)*dbh^(a*p-1) * (1-(dbh/b)^a)^(q-1))/(b^(a*p) * beta(p,q))
        return(d.gb)
      } #gb1.cr


#
#    probability of survival...
#
     pSurvive = function(dbh, params = params.run) {      
       pSurv = growth(dbh, params)*ssd.cr(dbh, params)/params['R']
       return(pSurv)
     } #pSurvive

#
#    probability of mortality...
#
     pMortality = function(dbh, params = params.run) {      
       pMort = params['M']*ssd.cr(dbh, params)/params['R']
       return(pMort)
     } #pMortality
      


#-------------------------------------------------------------------------------------      
#
#     reactive() returns a function that gets updated each time an "input" changes...
#

      CRD = reactive({ #inputs first...
        dcWidth = input$dc.width                #nominal width: for all but last
        R = input$recruitment
        M = input$mortality
        eta = input$eta
        gamma = input$gamma
        m = input$m
        #parameters...
        dmax = (gamma/eta)^(1/(m-1))          #CR  constraint
        params.run = c(m=m, gamma=gamma, eta=eta, R=R, M=M, dmax=dmax)
        #dbh class bounds: 1st and last can be different than others...
        if(includeTiny)             #first class width>dcWidth to include dbh=0
          adjustLowest = dcWidth/2
        else                        #first class width==dcWidth
          adjustLowest = 0.0
        dbhc = dbhClassLimits(dcWidth, dmax, dcWidth, forceIntegerWidth=FALSE,
                              adjustLowest = adjustLowest, ...)
        dbh.lo = dbhc$dbh.lo
        dbh.hi = dbhc$dbh.hi
        dbh = dbhc$dbh.mid
        ndclass = length(dbh)
        dmax.mid = dbh[ndclass]                    #need to reset
        dbh.hi[ndclass] = dmax                    #set to exact upper limit
        dcWidths = dbh.hi - dbh.lo                  #last (& 1st may) will be different from dcWidth
        #growth...
        #kappa = pi/(4*144)
        kappa = conv$k.ba
        cr.growth = growth(dbh, params.run)
        #stand table (midpoint method)...
        ba = kappa * dbh * dbh                   #by dbh class
        ntrees.app = ssd.cr(dbh, params.run)*dcWidths   #note: *dcwidths!
        #totTrees.app = sum(ntrees.app)
        BA.app = ba*ntrees.app             #by dbh class
        #exact qmsd & total BA from basd...
        qmsd = sqrt(qmsd.sq(params.run))
        totTrees = N.c(params.run)
        B.basd = qmsd*qmsd * totTrees * kappa
        #exact numbers and ba by dbh class...
        ntrees = rep(NA, ndclass)
        BA = ntrees
        for(i in seq_len(ndclass)) {
          ntrees[i] = integrate(gb1.cr, dbh.lo[i], dbh.hi[i], params=params.run)$value
          BA[i] = integrate(gb1.cr, dbh.lo[i], dbh.hi[i], params=params.run, alpha=2.0)$value
        }
        ntrees = totTrees * ntrees
        BA = B.basd * BA
        pMort = M*ntrees/R  #mortality density is a scaled version of numbers density==>exact
        #save the results...
        df = data.frame(lo=dbh.lo, dbh, hi=dbh.hi, ntrees, BA, dd.dt = cr.growth, pMort = pMort)
        df.approx = data.frame(lo=dbh.lo, dbh, hi=dbh.hi, ntrees.app, BA.app, dd.dt = cr.growth)
        rl = list(dcl = c(dmax = dmax,
                          dmax.mid = dmax.mid,
                          ndclass = ndclass,
                          dcWidth = dcWidth),   #nominal see dcWidths for all
                  dcWidths = dcWidths,
                  totals = c(totTrees = totTrees,
                             totBA = B.basd,
                             qmsd = qmsd),
                  params.run = params.run,
                  df = df,
                  df.approx = df.approx #midpoint method
                 )
        if(!is.na(returnVar))
          assign(returnVar, rl, envir=.GlobalEnv) #save in workspace
        return(rl)
      }) #reactive

             

#
#     the CR growth plot...
#
      output$crPlot = renderPlot({
        with(CRD()$df, {
                        plot(dbh, dd.dt, pch=19, col='gray60',
                             cex.lab=1.6, cex.axis=1.6, xlab='DBH',
                             ylab = 'Annual dbh growth')
             }) #with CRD
      }) #renderPlot

#             
#     plot of the numbers density...
#
      output$ssdPlot = renderPlot({
        with(CRD(), {
          with(df, {#note: d0:dmax limits can produce huge values at the endpoints
                    #      no good for plotting, stick with midpoint limits...
                    ndc = dcl['ndclass']
                    dd = seq(dbh[1], dbh[ndc], length.out=100)
                    #note that either ssd.cr or gb1.cr can be used below...
                    y.cr = ssd.cr(dd, params.run)*dcl['dcWidth']  #note: *dcwidth!
                    #ye = N.c(params.run) * gb1.cr(dd, params.run)*dcl['dcWidth']
                    ymax = max(ntrees, y.cr)  #so the curve will not get cut off
                    barplot(ntrees, dcWidths, 0.0, col='grey90', names.arg=dbh,
                            cex.lab=1.6, cex.axis=1.6, cex.names=1.6, xlab='DBH',
                            ylim=c(0,ymax), ylab='Number of trees'
                           )
                    #must offset by dbh.lo[1] to match barplot...
                    lines(dd - lo[1], y.cr, col='blue', lty='dashed')
                    rug(totals['qmsd'], col='red', lwd=1.25)
                }) #with df
              }) #with CRD()
      }) #renderPlot

      
#
#     basal area-dbh distn plot...
#             
      output$basdPlot = renderPlot({
        with(CRD(), {
          with(df, {#note: see above for ssdPlot...
                    ndc = dcl['ndclass']
                    dd = seq(dbh[1], dbh[ndc], length.out=100)
                    a = 1 - params.run['m']
                    b = dcl['dmax']
                    p = 1
                    q = params.run['M']/(params.run['gamma']*(1 - params.run['m']))
                    #note: *dcwidth!...
                    y = totals['totBA'] * gb1.cr(dd, params.run, alpha=2.0)*dcl['dcWidth']  
                    ymax = max(BA, y)  #so the curve will not get cut off
                    barplot(BA, dcWidths, 0.0, col='grey90', names.arg=dbh,
                            cex.lab=1.6, cex.axis=1.6, cex.names=1.6, xlab='DBH',
                            ylab='Basal area', ylim=c(0, ymax)
                           )
                    #must offset by dbh.lo[1] to match barplot...
                    lines(dd - lo[1], y, col='blue', lty='dashed')
                    rug(totals['qmsd'], col='red', lwd=1.25)
          }) #with df
        }) #with CRD()
      }) #renderPlot


#             
#     plot of the mortality density...
#
      output$pMortPlot = renderPlot({
        with(CRD(), {
          with(df, {#note: d0:dmax limits can produce huge values at the endpoints
                    #      no good for plotting, stick with midpoint limits...
                    ndc = dcl['ndclass']
                    dd = seq(dbh[1], dbh[ndc], length.out=100)
                    #note that either ssd.cr or gb1.cr can be used below; also note: *dcwidth!...
                    y.pm = ssd.cr(dd, params.run)*dcl['dcWidth']*params.run['M']/params.run['R'] 
                    y.ps = pSurvive(dd, params.run) * dcl['dcWidth']
                    ymax = max(pMort, y.pm, y.ps)  #so the curve will not get cut off
                    barplot(pMort, dcWidths, 0.0, col='grey90', names.arg=dbh,
                            cex.lab=1.6, cex.axis=1.6, cex.names=1.6, xlab='DBH',
                            ylim=c(0,1.05*ymax), ylab='Probability of Mortality & Survival'
                           )
                    #must offset by dbh.lo[1] to match barplot...
                    lines(dd - lo[1], y.pm, col='blue', lty='dashed')
                    points(hi - lo[1], pSurvive(hi, params.run) * dcl['dcWidth'], pch = 20, col='cadetblue4')
                    lines(dd - lo[1], y.ps, col='grey70')
                    rug(totals['qmsd'], col='red', lwd=1.25)
                }) #with df
              }) #with CRD()
      }) #renderPlot

              
#
#     the stand summary text info...
#  
      output$stand = renderText({
        with(CRD(), {
                paste('Units =', conv$units,
                     '\nDbh class width (nominal) =', dcl['dcWidth'],
                     '\nMax dbh =', format(dcl['dmax'], digits=4),
                     '\nMax dbh =', format(dcl['dmax.mid'], digits=4), '(midpoint)',
                     '\nNumber of Trees: N =', format(totals['totTrees'], digits=5),
                     '\nBasal area: B =',format(totals['totBA'], digits=5),
                     '\nQuadratic MSD =', format(totals['qmsd'], digits=5) #,
                     #'\nM/(gamma*(m-1)) =', format(
                     #  params.run['M']/(params.run['gamma']*(params.run['m']-1)), digits=5) 
                    )
        }) #with CRD
      }) #renderText


                  
#
#     the data frame table output -- either a regular html table or a fancy DT...
#
      if(!useDT) {
        output$crView = renderTable(CRD()$df,
                                    digits=c(0,1,1,1,1,2,2),     #xtable options
                                    include.rownames=FALSE   #xtable.print options
                                   ) #renderTable
      }
      else {
        output$crView = DT::renderDataTable(
                        formatRound(datatable(CRD()$df, options = list(dom='ltip')),  #dom==>no search
                                    3:7, c(1,1,2,2,3)
                                   ) #format
                        ) #DT::
      } #if

      
#
#     the exit button...
#
      observeEvent(input$stopIt, { stopApp()  })

      
#
#     warning if BA gets too large...
#
      output$warnBA = renderUI(
        with(CRD(),
                   if(totals['totBA'] > maxBA) 
                     tags$div(tags$h2(
                              HTML(paste(tags$span(style="color:red", 'Warning:'),
                                   " stand basal area is > ",
                                   tags$span(style="color:red", maxBA), sep = ""))
                     )) #tags$div
             ) #with
        ) #renderUI

      
#
#     remind the user if tiny trees are included and show the width of non-nominal classes...
#
      output$noteTiny = renderUI(
        with(CRD(), {hi.text = paste('Highest class width: [',df$lo[dcl['ndclass']],', ',
                                  round(dcl['dmax'],4),') = ',
                                  round(dcWidths[dcl['ndclass']], 4),
                                  sep='')
                     list(          #this list is what is returned from with() to renderUI
                          if(includeTiny) 
                            tags$div(tags$h4(
                              HTML(paste(tags$span(style="color:blue", 'Note:'),
                                   " smallest class includes trees down to ",
                                         tags$span(style="color:blue", d0), ' dbh', sep = "")
                                  ),
                              helpText(paste('Lowest class width: [',d0,', ',df$hi[1],'] = ',
                                       round(dcWidths[1],4), sep='')
                                       ),
                              helpText(hi.text)
                                            ) #tags$h4
                                    ) #tags$div
                          else
                            helpText(hi.text) 
                         ) #list
           }) #with
      ) #renderUI



      
    } #server

    return(server)
} #crDistnServer

