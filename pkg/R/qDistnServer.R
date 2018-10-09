qDistnServer = function(returnVar,
                        conv,
                        useDT,
                        shiftZero,
                        ...
                       )
{
#---------------------------------------------------------------------------
#
#   This builds the server interface function for the qDistnShiny function.
#
#   Arguments...
#     returnVar = variable name for the results of the app; NA=no results,
#                 otherwise a character string in .GlobalEnv
#     conv = conversion factors for units 
#     useDT = TRUE: the display uses the DataTable format (very nice);
#             FALSE: use the default html table format. 
#     shiftZero = TRUE: shift the lower DBH limit down to zero diameter;
#                 FALSE: no shift, just include small trees to the lower
#                        bound of the first dbh class
#     ... = gobbled
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
#            reactive() returns a function that gets updated each time an "input" changes...
#  
      dbq = reactive({
        q = input$q.value
        dcWidth = input$dc.width
        dcHalf = dcWidth/2
        ##dbh.max = ceiling(input$dbh.max)        #ceiling is not needed
        dbh.max = input$dbh.max
        dmax.mid = dbh.max - dcHalf             #at the midpoint
        if(shiftZero)             #first class width>dcWidth to include dbh=0
          #adjustLowest = dcWidth/2
          dt = dcHalf
        else                        #first class width==dcWidth
          #adjustLowest = 0.0
          dt = dcWidth
        #dbhc = dbhClassLimits(dcHalf, dmax.mid, dcWidth, forceIntegerWidth=FALSE)
        dbhc = dbhClassLimits(dt, dmax.mid, dcWidth, forceIntegerWidth=FALSE)
        #dbhc = dbhClassLimits(dcWidth, dmax.mid, dcWidth, forceIntegerWidth=FALSE,
        #                      adjustLowest = adjustLowest, ...)
        dbh.lo = dbhc$dbh.lo
        dbh.hi = dbhc$dbh.hi
        dbh = dbhc$dbh.mid
        bpa = input$bpa   
        ndclass = length(dbh)
        dmax.mid = dbh[ndclass]                   #need to reset
        ba = conv$k.ba * dbh * dbh
        q.seq = q^(0:(ndclass-1))                 #largest to smallest dbh
        nt.max = bpa/(sum(rev(ba)*q.seq))         # number in the largest class
        ntrees = rev(nt.max * q.seq)              #in all classes
        totTrees = sum(ntrees)
                  BA = ba*ntrees                            #calculated BPA/BPH 
                  df = data.frame(dbh.lo, dbh, dbh.hi, ntrees, BA) #, q.seq=rev(q.seq))
                  #if(!is.na(returnVar))
                  #  assign(returnVar, df, envir=.GlobalEnv) #save in workspace
                  totBA = sum(BA)
                  qMSD = sqrt(totBA/(totTrees*conv$k.ba))
                  #the scaled exponential distribution for comparison: y=K*exp(-a*dbh)...
                  a = log(q)/dcWidth                #a parameter
                  K = ntrees[1]/exp(-a*df$dbh[1])   #K parameter: midpoint approximation
###                  K = ntrees[1]/(exp(-a*df$dbh.lo[1]) - exp(-a*df$dbh.hi[1]))   #K parameter: exact
                  #ntrees.ed = K*exp(-a*dbh)
                  #a list for the return info...
                  rl = list(stParams = c(bpa = bpa,
                                         q = q,
                                         dcWidth = dcWidth,
                                         dbh.max = dbh.max,
                                         dmax.mid = dmax.mid,
                                         ndclass = ndclass
                                        ),
                            stSum = c(totTrees = totTrees,
                                      totBA = totBA,
                                      qMSD = qMSD
                                     ),
                            df = df,
                            edParams = c(a=a, K=K)
                            #ntrees.ed = ntrees.ed
                           ) #rl
                  if(!is.na(returnVar))
                    assign(returnVar, rl, envir=.GlobalEnv) #save in workspace
                 #exit with the list...
                  return(rl)
                }) #reactive
  

#
#            the barplot of trees...
#             
             output$qPlot = renderPlot({ with(dbq(), {
                                  with(df, {#note: 0:dbh.max limits can produce huge values at 0
                                            #no good for plotting, stick with midpoint limits...
                                            ndc = stParams['ndclass']
                                            dd = seq(dbh[1], dbh.hi[ndc], length.out=100)
                                            barplot(ntrees, stParams['dcWidth'], 0.0, col='grey90',
                                                    names.arg = dbh,
                                                    cex.lab=1.6, cex.axis=1.6, cex.names=1.6,
                                                    xlab='DBH', ylab='Number of trees'
                                                    )
                                            #add the exponential curve...
                                            a = edParams['a']
                                            K = edParams['K']
                                            lines(dd-dbh.lo[1], K*exp(-a*dd), col='blue', lty='dashed', lwd=1.4)
                                            rug(stSum['qMSD'], col='red', lwd=1.25)
                                         }) #df
                                    }) #with dbq()
                             }) #renderPlot
             

#
#            basal area-dbh distn plot...
#             
             output$basdPlot = renderPlot({ with(dbq(), {
                                  with(df, {
                                            ndc = stParams['ndclass']
                                            dmax = dbh.hi[ndc]
                                            barplot(BA, stParams['dcWidth'], 0.0, col='grey90',
                                                    names.arg = dbh, cex.lab=1.6, 
                                                    cex.axis=1.6, cex.names=1.6,
                                                    xlab='DBH', xlim=c(0.0, dmax),
                                                    ylab='Basal area', ylim=c(0,1.05*max(BA))
                                                   )
                                            #the size-biased for ba == gamma...
                                            dd = seq(0.0, dmax, length.out=100)
                                            a = edParams['a']
                                            K = edParams['K']
                                            Dbarq.sq = 2/a^2
                                            N = K/a
                                            B = Dbarq.sq * N * conv$k.ba
                                            basd = B * dgamma(dd, shape=3, scale= 1/edParams['a'])
                                            lines(dd-dbh.lo[1], basd, col='blue', lty='dashed')
                                            rug(stSum['qMSD'], col='red', lwd=1.2)
                                         }) #df
                                 }) #with dbq()
                               }) #renderPlot
              
#
#           the stand summary text info...
#  
            output$stand = renderText({  with(dbq(), {
               paste('Units =', conv$units,
                     '\nDbh class width =', stParams['dcWidth'],
                     '\nMax dbh =', format(stParams['dmax.mid'], digits=4), '(midpoint)',
                     '\nN =', format(stSum['totTrees'], digits=5),
                     '\nStand BA =',format(stSum['totBA'], digits=5),
                     '\nQuadratic MSD =', format(stSum['qMSD'], digits=3)
                    )
               }) #with
            }) #renderText
                  
#
#          the data frame table output...
#
#          the data frame table output -- either a regular html table or a fancy DT...
#
           if(!useDT) {
             output$qView <- renderTable({ with(dbq(), { df  })},
                                         digits=c(0,1,1,1,1,2),     #xtable options
                                         include.rownames=FALSE   #xtable.print options
                                         ) #renderTable
           }
           else {
             output$qView = DT::renderDataTable(formatRound(
                              datatable(dbq()$df, options = list(dom='ltip')),  #dom==>no search
                              4:5, c(2,1)) )
           }

#
#            the exit button...
#
             observeEvent(input$stopIt, { stopApp()  })

             
    } #server

    return(server)
} #qDistnServer

