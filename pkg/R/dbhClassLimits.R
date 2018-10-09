dbhClassLimits = function(dbh.low = 1,
                          dbh.up = 20,
                          dbh.width = 1,
                          forceIntegerWidth = TRUE,
                          adjustLowest = NULL,
                          classNames = c('dbh.low','dbh.mid','dbh.hi'),
                          runQuiet = TRUE,
                          ...
                         )
{
#------------------------------------------------------------------------------------------
#   dbhClassLimits
#
#   This routine will simply set up the matrix of lower, mid, and upper dbh class limits
#   used as input to groupDBH().
#
#   Arguments...
#     dbh.low = smallest dbh class midpoint
#     dbh.up = largest dbh class midpoint
#     dbh.width = dbh class width
#     classNames = a vector of length 3 for the labels of the data frame returned; always
#                  in the form lower, mid and upper.
#     forceIntegerWidth = TRUE: dbh class widths must be integer; FALSE: can be real
#     adjustLowest = NULL, NA, or <=0 means no adjustment; other positive means to
#                    subtract this amount from the lower bound of the smallest class
#     runQuiet = T: no feedback; F: a little feedback
#     ... = collects unused arguments
#
#   Note: the values in dbh.low and dbh.up are for the midpoint of the class!!! <<****
#
#   Returns...
#     a data frame with the desired intervals all set up for groupDBH()
#
#   One simple use of these two functions is to produce a dbh frequency table...
#
#   table(groupDBH(c(10,13,22,22,30), dbhClassLimits(c(5,30,2))))
#   9 13 21 29
#   1  1  2  1
#
#   In the above, to get this into usable form use...
#
#   as.data.frame(table(groupDBH(c(10,13,22,22,30), dbhClassLimits(c(5,30,2)))))
#
#   which will return the above output in two column form.
#
#   Note, however, with this approach that table() has no way of returning intervals with
#   zero counts, i.e., just the diameter classes with diameters present are returned
#   with the frequencies. Thus, if we wanted to use this information subsequently, in
#   a place where all levels are required, we could not do so unless we tell table() what
#   all of the levels (classes) are; i.e., turn the results into a factor...
#
#   lims = dbhClassLimits(c(5,30,2))
#   tmp = factor( groupDBH(c(10,13,22,22,30), lims), levels=lims[,'mid'])
#   table(tmp)
#
#   An alternative to the above is to use hist(x, breaks=c(...), plot=F) and use the
#   "mids" and "counts" fields; or use the equivalent with cut(), which returns a factor
#   with all levels...
#
#   table( cut(c(10,13,22,22,30), dbhClassLimits(c(5,30,2))[,'lower']) )
#
#   The routine testMixDist.r() provides examples of each of these.
#
#
#   J. H. Gove			jhgove@christa.unh.edu         June 2006 & Feb 2010
#   USDA Forest Service
#   Northeastern Research Station
#   P.O. Box 640
#   Durham, NH 03824		(603) 868-7667
#------------------------------------------------------------------------------------------
#
#   a few checks...
#
    if(length(dbh.low)!=1 || length(dbh.up)!=1 || length(dbh.width)!=1)
      stop('All arguments must be of length 1!')
    if(dbh.low > dbh.up)
      stop('Upper limit (',dbh.up,') must be greater than lower (',dbh.low,')\n')
    if(dbh.width <= 0)
      stop('Negative dbh class sizes (',dbh.width,') are not allowed!\n')
    if(forceIntegerWidth) {
      if(dbh.width%%1 > 0) {
        dbh.width = round(dbh.width)
        warning('dbh.width rounded to nearest integer value.\n')
      }
    }
    if((dbh.up-dbh.low) < dbh.width)
      stop('Dbh class width (',dbh.width,') is larger than the dbh range!\n')

    if(is.null(adjustLowest) || is.na(adjustLowest) || adjustLowest<0)
      adjustLowest = 0
    
#   check lower limit below zero...    
    halfWidth = dbh.width/2
    lowCheck = dbh.low - halfWidth - adjustLowest
    if(lowCheck < 0)
      stop(paste('Lowest dbh limit (',lowCheck,') < 0!'))

    n.lims = 3L   #number of columns in data frame for limits
    if(length(classNames) != n.lims || !is.character(classNames))
      stop('classNames must be a character vector of length 3!')
     
    
#
#   get the dbh midpoints...
#
    dbh = seq(dbh.low, dbh.up, by=dbh.width)

#
#   setup the dbh class limits data frame...
#
    ndbhclass = trunc((dbh.up-dbh.low)/dbh.width) + 1
    dbhClassLimits = as.data.frame(matrix(NA, nrow=ndbhclass, ncol=n.lims))
    colnames(dbhClassLimits) = classNames
    dbhClassLimits[,1] = seq(dbh.low - halfWidth, dbh.up - halfWidth, by= dbh.width)
    dbhClassLimits[1,1] = dbhClassLimits[1,1] - adjustLowest
    dbhClassLimits[,3] = seq(dbh.low + halfWidth, dbh.up + halfWidth, by= dbh.width)
    dbhClassLimits[,2] = dbh

    dbhClassLimits = structure(dbhClassLimits, createdBy='dbhClassLimits')

    if(!runQuiet) {
      cat('\nNumber of DBH classes =',ndbhclass)
      cat('\nLowest class limit =',dbhClassLimits[1,1])
      cat('\nLargest class limit =',dbhClassLimits[ndbhclass,n.lims])
      cat('\nLowest limit adjustment =', adjustLowest)
      cat('\n')
    }

    return(dbhClassLimits)
}   #dbhClassLimits
