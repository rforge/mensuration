findq = function(dbh.low = 0.5,
                 dbh.up = 19.5,
                 dc.width = 1,
				 B = 80,
				 N.max = 2,
				 units = c('English', 'metric'),
                 runQuiet = FALSE,
				 ...
				)
{
#---------------------------------------------------------------------------
#
#   This is the latest version of the routine to find the q value from the
#   basal area and number of trees in the maximum dbh class. The objective
#   function is defined within. The original is in the Stocking/unevenaged/q
#   workspace
#
#   Arguments...
#     dbh.low = smallest dbh class midpoint
#     dbh.up = largest dbh class midpoint
#     dc.width = dbh class width
#     B = basal area per acre/hectare
#     N.max = number of trees in the largest class
#     units = choose either English or metric
#     runQuiet = TRUE: no report; FALSE: report
#
#   Returns...
#     the estimated q-value
#
#   Arguments for internal qObj function...
#     q = the estimated q-value for an iteration of uniroot()
#     mindbh = dbh.low
#     maxdbh = dbh.up
#     dc.width = dc.width
#     baPerUA = B
#     TreesInMax = N.max
#     dc.ba = the midpoint basal area for each class
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
#   trivia...
#
    units = match.arg(units)
    if(units=='metric') {
      k.ba = pi/40000
      uarea = 'ha'
      incm = 'cm'
    }
    else {
      k.ba = pi/(4*144)
      uarea = 'acre'
      incm = 'inch'
    }

#
#   let dbhClassLimits flag errors in input...
#
    dbhc = dbhClassLimits(dbh.low, dbh.up, dc.width, ...)
    dbh.mid = dbhc$dbh.mid
    dc.ba = k.ba * dbh.mid * dbh.mid
    ndclass = length(dbh.mid)

#
#   the objective for uniroot() -- should nearly equal zero on solutio...
#
    qObj = function(q, mindbh, maxdbh, dc.width, baPerUA, TreesInMax, dc.ba) {
                    q.seq = rev(q^(0:(ndclass-1)))             #largest to smallest dbh
                    obj = baPerUA/TreesInMax - sum(dc.ba * q.seq)
	                return(obj)
                  }
    
#
#	now search for the root...
#
   	ans = uniroot(f = qObj, lower = 1, upper = 6,
                  mindbh=dbh.low, maxdbh=dbh.up, dc.width=dc.width, baPerUA=B,
                  TreesInMax=N.max, dc.ba=dc.ba
                  )
    qValue = ans$root

#
#   results if desired...
#    
    if(!runQuiet) {
      cat('\nResults...\n')
	  cat('------------------\n')
      cat('Minimum dbh =', dbh.low)
      cat('\nMaximum dbh =', dbh.up)
      cat('\nBasal area/', uarea,' = ', B, sep='')
      cat('\nTrees in the ',dbh.up, ' ', incm,' class = ', N.max, sep='')
      cat('\n',dc.width,'-',incm,' q value = ', round(qValue, 2) ,sep='')
      cat('\n...objective at q root = ',ans$f.root)
      cat('\n...iterations =',ans$iter)
      cat('\n...approximate precision =',ans$estim.prec)
      cat('\n------------------\n')
    }

#
#   complete the stand table...
#    
    q.seq = rev(qValue^(0:(ndclass-1)))
    df = dbhc
    df$N = q.seq*N.max
    df$B = dc.ba * df$N

    
#
#	and done...
#
	return(invisible(list(qValue = qValue,
                          qTable = df
                          )
                     )
          )
}   #findq
