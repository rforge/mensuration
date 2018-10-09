standTable = function(dbh,
                      wantBA = TRUE,
                      isEnglish=TRUE,
                      runQuiet=TRUE,
                      ...
                     )
{
#------------------------------------------------------------------------------------------
#
#   This little routine will calculate a stand table from some diameters. See documentation
#   for dbhClassLimits and groupDBH as these are used to get the dbh class groupings for
#   the diameters. Basal area is also calculated if desire.
#
#   Arguments...
#     dbh = a vector of tree diameters
#     wantBA = T: calculate basal area too; F: no BA
#     isEnglish = T: English; F: metric
#     runQuiet = T: no feedback; F: a little feedback
#     ... = passed on to other routines
#
#   J. H. Gove			jhgove@christa.unh.edu         24-Feb-2010
#   USDA Forest Service
#   Northeastern Research Station
#   P.O. Box 640
#   Durham, NH 03824		(603) 868-7667
#------------------------------------------------------------------------------------------
#
#
    dbh.lims = dbhClassLimits(runQuiet=runQuiet, ...)
    gd = groupDBH(dbh, dbh.lims)
    dbh.grp = gd$dbh.grp
    dbh.fac = factor(dbh.grp, levels=dbh.lims[,2])

#
#   tree frequency...
#    
    st = dbh.lims
    st$trees = as.data.frame(table(dbh.fac))$Freq

#
#   basal area...
#
    if(wantBA) {
      if(isEnglish)
        k = pi/(4*144)
      else
        k = pi/(4*10000)
      ba = k*dbh^2
      j = data.frame(gd$idx, ba)
      ba.sum = tapply(j$ba, j$gd.idx, sum)
      st$ba = NA
      bdx = as.numeric(names(ba.sum))    #dbh classes are names in ba.sum
      st[bdx,'ba'] = ba.sum
      st$ba = ifelse(is.na(st$ba), 0, st$ba)
    }

    
    if(!runQuiet) {
      line = paste(rep('-', 40), collapse='')
      cat(line,sep='')
      if(wantBA) {
        cat('\nNote: BA is calculated from raw diameters.')
        cat('\nTotal BA =', sum(st$ba))
      }
      cat('\nTotal Trees =',sum(st$trees))
      cat('\n',line,sep='')
      cat('\n\n')
    }    

    return(st)
}   #standTable

