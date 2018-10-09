Units = function(x=c(3,5), from=c('m'), to='ft' )
{
#------------------------------------------------------------------------------------------
#
#   This little routine is a wrapper for the linux "units" command. It will take conformable
#   vectors or try to make them such in the form below...
#
#   Arguments...
#     x = a vector of quantities in the "from" units to convert from
#     from = a vector of units that appliy to x; if scalar, it will be repeated to the
#            length of x
#     to = a vector of units that we want x converted to; it will be recycled as above
#
#   Returns...
#     a data frame with the command in the first column, the original units to be converted
#     in the second, and the conversion factors (or converted units total) in the 3rd
#
#
#   Note that the conversion factors can be either factors you will multiply by later,
#   or the actual answer. Examples include...
#
#     Units(1,'m','ft')     #this produces a conversion factor from meters to feet
#
#   whereas the following (more useful) gives the converted quantities of interest...
#
#     Units(c(25,50),'m^2/hectare','ft^2/acre')      #basal area
#     Units(c(250,500),'ft^3/acre','m^3/hectare')    #volume
#
#   and this illustrates that we can have mixes of units...
#
#     Units(c(250,1,1),c('ft^3/acre','acre','hectare'),c('m^3/hectare','hectare','acre'))
#
#
#   J. H. Gove			jhgove@christa.unh.edu         Oct 2011
#   USDA Forest Service
#   Northeastern Research Station
#   P.O. Box 640
#   Durham, NH 03824		(603) 868-7667
#------------------------------------------------------------------------------------------
#
#   few checks, up to the user to get thing right for now...
#
    nx = length(x)
    nf = length(from)
    nt = length(to)
    if(nf > nx)                               #truncate extras
      from = from[nx]
    if(nt > nx)
      to = to[nx]

#
#   put the whole command together into columns of a matrix first...
#
    xFrom = paste('\'',x,from,'\'', sep='')   #we need to quote the units, don't use sQuote!
    xTo = paste('\'',to,'\'', sep='')         #same
    kk = data.frame(rep('units -1', nx))
    kk$From = xFrom                           #recycling used here if possible
    kk$To = xTo                               #and here

#
#   then splice it together...
#
    aa = apply(kk,1,paste,collapse=' ')       #collapse into one string for each entry
    suppressWarnings({
      zz = sapply(aa,system,intern=T)         #send it to be evaluated
    })
    unknown = grep('Unknown', zz)
    if(length(unknown) > 0) {                  #any bad conversions?
      cat('\n***> \"units\" conversion errors at index:',unknown,'\n')
      stop('at least one of your conversions is illegal!') #could return them here!
    }
    zz = sapply(zz,strsplit,'\\*')            #split on '*' for each item
    zz = as.numeric(sapply(zz,'[[',2))        #and save only the multiplier conversion
    
    df = data.frame(as.character(aa), stringsAsFactors=FALSE)         #pack it up
    df$x = x                                  #original units outside command
    df$conv = zz                              #with the conversion
    colnames(df) = c('systemCommand','orig','conv')

    return(df)
}   #Units

