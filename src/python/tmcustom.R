## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=TRUE,
                      conf.interval=.95, .drop=TRUE) {
    require(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This is does the summary; it's not easy to understand...
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun= function(xx, col, na.rm) {
                           c( N    = length2(xx[,col], na.rm=na.rm),
                              mean = mean   (xx[,col], na.rm=na.rm),
                              sd   = sd     (xx[,col], na.rm=na.rm)
                              )
                          },
                    measurevar,
                    na.rm
             )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean"=measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}

## Norms the data within specified groups in a data frame; it normalizes each
## subject (identified by idvar) so that they have the same mean, within each group
## specified by betweenvars.
##   data: a data frame.
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   na.rm: a boolean that indicates whether to ignore NA's
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
    require(plyr)

    # Measure var on left, idvar + between vars on right of formula.
    data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
     .fun = function(xx, col, na.rm) {
        c(subjMean = mean(xx[,col], na.rm=na.rm))
      },
      measurevar,
      na.rm
    )

    # Put the subject means with original data
    data <- merge(data, data.subjMean)

    # Get the normalized data in a new column
    measureNormedVar <- paste(measurevar, "_norm", sep="")
    data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
                               mean(data[,measurevar], na.rm=na.rm)

    # Remove this subject mean column
    data$subjMean <- NULL

    return(data)
}

##Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, un-normed mean, normed mean (with same between-group mean),
##   standard deviation, standard error of the mean, and confidence interval.
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=TRUE, conf.interval=.95, .drop=TRUE) {

  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
    FUN=is.factor, FUN.VALUE=logical(1))

  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }

  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)

  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL

  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)

  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")

  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)

  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                           FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )

  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor

  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}







## some helpful threads
## https://stat.ethz.ch/pipermail/r-help/2008-September/172641.html
## http://tolstoy.newcastle.edu.au/R/e4/help/08/02/4875.html
## http://tolstoy.newcastle.edu.au/R/e2/help/07/01/8598.html

## http://www.r-statistics.com/wp-content/uploads/2011/01/boxplot-add-label-for-outliers.r.txt

## last updated: 2013-08-21: added require2 function (originally from the installr package)



#require2 <- function (package, ask = TRUE, ...) 
#{
    #package <- as.character(substitute(package))
    #if (!suppressWarnings(require(package = package, character.only = TRUE))) {
        #install_package <- ask.user.yn.question(paste("Package ", 
            #package, " is not installed. Do you want to install it now?"))
        #if (install_package) 
            #install.packages(pkgs = package)
    #}
    #require(package = package, character.only = TRUE)
#}






#boxplot.with.outlier.label <- function(y, label_name, ..., spread_text = T, data, plot = T, range = 1.5, label.col = "blue", push_text_right = 1.5, # enlarge push_text_right in order to push the text labels further from their point
#segement_width_as_percent_of_label_dist = .45, # Change this if you want to have the line closer to the label (range should be between 0 to 1
	#jitter_if_duplicate = T, jitter_only_positive_duplicates = F)
#{	
	## change log:
	## 19.04.2011 - added support to "names" and "at" parameters.


	## jitter_if_duplicate - will jitter (Actually just add a bit of numbers) so to be able to decide on which location to plot the label when having identical variables...
	#require2(plyr) # for is.formula and ddply

	## a function to jitter data in case of ties in Y's
	#jitter.duplicate <- function(x, only_positive = F)
	#{
		#if(only_positive) {
			#ss <- x > 0
		#} else {
			#ss <- T
		#}	
		#ss_dup <- duplicated(x[ss])
		## ss <- ss & ss_dup
		#temp_length <- length(x[ss][ss_dup])	
		#x[ss][ss_dup] <- x[ss][ss_dup] + seq(from = 0.00001, to = 0.00002, length.out = temp_length)
		#x
	#}
	## jitter.duplicate(c(1:5))
	## jitter.duplicate(c(1:5,5,2))
	## duplicated(jitter.duplicate(c(1:5,5,2)))
	## jitter.duplicate(c(0,0,1:5,5,2))
	## duplicated(jitter.duplicate(c(0,0,1:5,5,2)))


	
	## handle cases where 
	#if(jitter_if_duplicate) {
		## warning("duplicate jutter of values in y is ON")
		#if(!missing(data)) {	#e.g: we DO have data
			## if(exists("y") && is.formula(y)) {		# F && NULL # F & NULL
			#y_name <- as.character(substitute(y))	# I could have also used as.list(match.call())
												## credit to Uwe Ligges and Marc Schwartz for the help
												## https://mail.google.com/mail/?shva=1#inbox/12dd7ca2f9bfbc39
			#if(length(y_name) > 1) {	# then it is a formula (for example: "~", "y", "x"
				#model_frame_y <- model.frame(y, data = data)
				#temp_y <- model_frame_y[,1]
				#temp_y  <- jitter.duplicate(temp_y, jitter_only_positive_duplicates)	# notice that the default of the function is to work only with positive values...
				## the_txt <- paste(names(model_frame_y)[1], "temp_y", sep = "<<-") # wrong...
				#the_txt <- paste("data['",names(model_frame_y)[1],"'] <- temp_y", sep = "")				
				#eval(parse(text = the_txt))	# jutter out y var so to be able to handle identical values.
			#} else {	# this isn't a formula
				#data[,y_name] <- jitter.duplicate(data[,y_name], jitter_only_positive_duplicates)
				#y <- data[,y_name]	# this will make it possible for boxplot(y, data) to work later (since it is not supposed to work with data when it's not a formula, but now it does :))
			#}		
		#} else {	# there is no "data"		 
			#if(is.formula(y)) { # if(exists("y") && is.formula(y)) {		# F && NULL # F & NULL
				#temp_y <- model.frame(y)[,1]
				#temp_y  <- jitter.duplicate(temp_y, jitter_only_positive_duplicates)	# notice that the default of the function is to work only with positive values...
				#temp_y_name <- names(model.frame(y))[1]	# we must extract the "names" before introducing a new enbironment (or there will be an error)
				#environment(y) <- new.env()
				#assign(temp_y_name, temp_y, environment(y))
					## Credit and thanks for doing this goes to Niels Richard Hansen (2 Jan 30, 2011)
					## http://r.789695.n4.nabble.com/environment-question-changing-variables-from-a-formula-through-model-frame-td3246608.html
				## warning("Your original variable (in the global environemnt) was just jittered.")	# maybe I should add a user input before doing this....
				## the_txt <- paste(names(model_frame_y)[1], "temp_y", sep = "<<-")
				## eval(parse(text = the_txt))	# jutter out y var so to be able to handle identical values.
			#} else {
				#y <- jitter.duplicate(y, jitter_only_positive_duplicates)
			#}		
		#}
	#}
	## the_txt <- paste("print(",names(model_frame_y)[1], ")")
	## eval(parse(text = the_txt))	# jutter out y var so to be able to handle identical values.
	## print(ls())

	
	## y should be a formula of the type: y~x, y~a*b
	## or it could be simply y
	#if(missing(data)) {
			#boxdata <- boxplot(y, plot = plot,range = range ,...)
		#} else {
			#boxdata <- boxplot(y, plot = plot,data = data, range = range ,...)
		#}
	#if(length(boxdata$names) == 1 && boxdata$names =="") boxdata$names <- 1	# this is for cases of type: boxplot(y) (when there is no dependent group)
	#if(length(boxdata$out) == 0 ) {
		#warning("No outliers detected for this boxplot")
		#return(invisible())
		#}
	
	#if(!missing(data)) attach(data)	# this might lead to problams I should check out for alternatives for using attach here...
	

	## creating a data.frame with information from the boxplot output about the outliers (location and group)
	#boxdata_group_name <- factor(boxdata$group)
	#levels(boxdata_group_name) <- boxdata$names[as.numeric(levels(boxdata_group_name))]	# the subseting is for cases where we have some sub groups with no outliers
	#if(!is.null(list(...)$at))	{	# if the user chose to use the "at" parameter, then we would like the function to still function (added on 19.04.2011)
		#boxdata$group <- list(...)$at[boxdata$group]		
		#}
	#boxdata_outlier_df <- data.frame(group = boxdata_group_name, y = boxdata$out, x = boxdata$group)
	

	## Let's extract the x,y variables from the formula:
	#if(is.formula(y))
	#{
		#model_frame_y <- model.frame(y)
			## old solution: (which caused problems if we used the names parameter when using a 2 way formula... (since the order of the names is different then the levels order we get from using factor)
			## y <- model_frame_y[,1]
			## x <- model_frame_y[,-1]

		#y <- model_frame_y[,1]
		#x <- model_frame_y[,-1]
		#if(!is.null(dim(x))) {	# then x is a matrix/data.frame of the type x1*x2*..and so on - and we should merge all the variations...
			#x <- apply(x,1, paste, collapse = ".")
		#}
	#} else {
		## if(missing(x)) x <- rep(1, length(y))
		#x <- rep(1, length(y))	# we do this in case y comes as a vector and without x
	#}	
	
	## and put all the variables (x, y, and outlier label name) into one data.frame
	#DATA <- data.frame(label_name, x ,y)
	
	#if(!is.null(list(...)$names))	{	# if the user chose to use the names parameter, then we would like the function to still function (added on 19.04.2011)
		#DATA$x <- factor(DATA$x, levels = unique(DATA$x))
		#levels(DATA$x) = list(...)$names	# enable us to handle when the user adds the "names" parameter # fixed on 19.04.11	# notice that DATA$x must be of the "correct" order (that's why I used split above
		## warning("Careful, the use of the 'names' parameter is experimental.  If you notice any errors please e-mail me at: tal.galili@gmail.com")
		#}

	#if(!missing(data)) detach(data)	# we don't need to have "data" attached anymore.

	## let's only keep the rows with our outliers 
	#boxplot.outlier.data <- function(xx, y_name = "y")
	#{
		#y <- xx[,y_name]
		#boxplot_range <- range(boxplot.stats(y, coef = range )$stats)
		#ss <- (y < boxplot_range[1]) | (y > boxplot_range[2])
		#return(xx[ss,])	
	#}
	#outlier_df <-ddply(DATA, .(x), boxplot.outlier.data)
	

	## create propor x/y locations to handle over-laping dots...
	#if(spread_text) {
		## credit: Greg Snow
		#require2(TeachingDemos)		
		#temp_x <- boxdata_outlier_df[,"x"]
		#temp_y1 <- boxdata_outlier_df[,"y"]
		#temp_y2 <- temp_y1
		#for(i in unique(temp_x))
		#{
			#tmp <- temp_x == i
			#temp_y2[ tmp ] <- spread.labs( temp_y2[ tmp ], 1.3*strheight('A'), maxiter=6000, stepsize = 0.05) #, min=0 )
		#}
		
	#}
	

	
	## max(strwidth(c("asa", "a"))
	## move_text_right <- max(strwidth(outlier_df[,"label_name"]))	
	
	## plotting the outlier labels :)  (I wish there was a non-loop wise way for doing this)
	#for(i in seq_len(dim(boxdata_outlier_df)[1]))
	#{
		## ss <- (outlier_df[,"x"]  %in% boxdata_outlier_df[i,]$group) & (outlier_df[,"y"] %in% boxdata_outlier_df[i,]$y)

		## if(jitter_if_duplicate) {
			## ss <- (outlier_df[,"x"]  %in% boxdata_outlier_df[i,]$group) & closest.number(outlier_df[,"y"]  boxdata_outlier_df[i,]$y)
		## } else {
		#ss <- (outlier_df[,"x"]  %in% boxdata_outlier_df[i,]$group) & (outlier_df[,"y"] %in% boxdata_outlier_df[i,]$y)
		## }

		#current_label <- outlier_df[ss,"label_name"]
		#temp_x <- boxdata_outlier_df[i,"x"]
		#temp_y <- boxdata_outlier_df[i,"y"]		
		## cbind(boxdata_outlier_df,		temp_y2)
		## outlier_df

		
		
		#if(spread_text) {
			#temp_y_new <- temp_y2[i] # not ss			
			#move_text_right <- strwidth(current_label) * push_text_right
			#text( temp_x+move_text_right, temp_y_new, current_label, col = label.col)			
			## strwidth
			#segments( temp_x+(move_text_right/6), temp_y, temp_x+(move_text_right*segement_width_as_percent_of_label_dist), temp_y_new )
		#} else {
			#text(temp_x, temp_y, current_label, pos = 4, col = label.col)
		#}		
	#}

	## outputing some of the information we collected
	#list(boxdata = boxdata, boxdata_outlier_df = boxdata_outlier_df, outlier_df=outlier_df)
#}






#########################################
#### examples to see that it works

## library(plyr)
## library(TeachingDemos)
## source("http://www.r-statistics.com/wp-content/uploads/2011/01/boxplot-with-outlier-label-r.txt") # Load the function
## set.seed(210)
## n <- 20
## y <- rnorm(n)
## x1 <- sample(letters[1:3], n,T)
## lab_y <- sample(letters, n)
## boxplot.with.outlier.label(y~x1, lab_y, push_text_right = 1.5, range = .3)
## data.frame(y, x1, lab_y)

## set.seed(10)
## x2 <- sample(letters[1:3], n,T)
## boxplot.with.outlier.label(y~x1*x2, lab_y, push_text_right = 1.5, range = .3)
## data.frame(y, x1, x2, lab_y)



# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# difference between correlations
diff.corr <- function( r1, n1, r2, n2 ){ 
    Z1 <- 0.5 * log( (1+r1)/(1-r1) ) 
    Z2 <- 0.5 * log( (1+r2)/(1-r2) ) 
    diff   <- Z1 - Z2 
    SEdiff <- sqrt( 1/(n1 - 3) + 1/(n2 - 3) ) 
    diff.Z  <- diff/SEdiff 
    p <- 2*pnorm( abs(diff.Z), lower=F) 
    cat( "Two-tailed p-value", p , "\n" ) 
  } 
