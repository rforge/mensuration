groupDBH = function(dbh, classLimits) {
#------------------------------------------------------------------------------------------
#   groupDBH
#
#   This routine will just group diameters by the first column of the classLimits (lower limit)
#   argument. The idea was derived from the built-in R function 'findInterval,' which
#   does not seem to work correctly w/r to the rightmost.close parameter.
#
#   Note: the outer() creates a num.dbh x num.limits matrix of T,F values for the test;
#         these are then summed (T=1, F=0) across rows (limits) yielding the index
#         to the correct interval, based on the lower limit of the interval, in the
#         apply() portion.
#   Note: In the above note, I have now changed it so that anything below the lowest
#         class lower limit or above the upper limit of the highest class gets assigned NA. 
#
#   Note that diameters larger than the last limit get lumped correctly because the sum
#   over all True values is the number of rows; however, for diameters smaller than the 
#   lower limit the sum is zero (all False), so we must help it along with the ifelse()
#   check to get the index for the first class limits.
#         
#
#   Input...
#       dbh  = a vector of diameters
#       classLimits = a data frame of dbh limits with column 1 being the lower, 2 being the
#               midpoint and 3 being the upper limits, though only lower and mid are used.
#
#   J. H. Gove			jhgove@christa.unh.edu         June 2005 & Feb 2010
#   USDA Forest Service
#   Northeastern Forest Experiment Station
#   P.O. Box 640
#   Durham, NH 03824		(603) 868-7667
#------------------------------------------------------------------------------------------
#
#   a quick check first...
#
    creator = attr(classLimits, 'createdBy')
    if(is.null(creator) || creator!='dbhClassLimits')
      stop('classLimits must be a data frame from dbhClassLimits!')
      
#
#   get rid of zero diameter trees...
#
    dbh = dbh[ifelse(dbh>0, TRUE, FALSE)]
    if(length(dbh)==0)                #all zero diameter trees...exit
      return(NULL)

    classNames = colnames(classLimits)       #lower, mid and upper names
    ndbhclass = nrow(classLimits)
    
#
#   index them, then get their actual dbh classes according to the class limits passed...
#
    idx = ifelse(dbh <= classLimits[1,classNames[1]], NA,
                 apply(outer(dbh, classLimits[,classNames[1]], '>'), 1, sum)
                )
    idx = ifelse(dbh >= classLimits[ndbhclass, classNames[3]], NA, idx)
    classed.dbh = classLimits[idx, classNames[2]]

    if(any(is.na(classed.dbh)))
      warning('Some diameters fell outside the class limits passed;\n these have been discarded.')
    
    return( list(dbh.grp = classed.dbh, idx = idx) )
}   #groupDBH
