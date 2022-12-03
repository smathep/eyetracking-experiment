################### IMPORTANT NOTICE #######################
# Below you'll find  4 functions: 
# TransWMatrix  
#   calculates transition matrix between AOIs - returns transtion matrix M, 
# TransWPlot 
#   plots tranistion matrix - returns TM plot with probabilities in each cell; 
# TransWEntropy  
#   calculates transition matrix entropy for each participant - returns vector
#   of Shannon's H.
#
####  1. Use TransWMatrix before using TransWPlot
####  1. Data for functions 1 to 3 have to be in the form of 'data.frame'.
####     Each raw represents seperate fixation. 
################################################################################

### Required packages
load.libraries(c('gplots','RColorBrewer','colorspace','entropy','gmodels','car','lattice'))
#library(gplots)
#library(RColorBrewer)
#library(colorspace)
#library(entropy)
#library(gmodels)
#library(car)
#library(lattice)

spanletters <- function(min, max) {
  v <- list()
  for (i in min:max) {
    v[length(v) + 1L] <- paste(i,sep="")
    }
  return(v)
}

zeroWM <- function(min,max) {
  rows = 1
  cols = max - min + 1
  AOIs <- spanletters(min,max)
# print(sprintf("AOIs:\n"))
# print(AOIs)
# print(sprintf("zeroWM: nrow x ncol: %d x %d",rows,cols))
  M <- matrix(data=0,
              nrow=rows,
              ncol=cols,
              dimnames=list(1,AOIs)) # empty matrix
# print(M)
  return(M)
}

TransWMatrix <- function(M,data,StimulusVar,SubjectsVar,SpanVar,print=FALSE) {
  data <- data[order(data[ ,StimulusVar]), ]    
  uniqS <- sort(unique(data[ ,StimulusVar])) # unique values of stimulus var
  lui <- length(unique(data[ ,StimulusVar])) # how many stimuli 
# empty matrix
  M[,] <- 0
  # for each stimulus
  for (i in 1:lui) {
    kukudf <- data[which(data[ ,StimulusVar] == uniqS[i]), ] # choose stimulus i
    luj <- dim(kukudf)[1] # how many stimuli for subject i
    # for each subject
    if (luj > 1) {
      for (j in 1:luj) {
        span <- kukudf[j, SpanVar]
        # increase occurence of detected span
        M[1,as.character(span)] <- M[1,as.character(span)] + 1
      }
    } else {
      print('Warning: ')
      print(sprintf("%s%s", 'No word spans for stimulus: ', uniqS[i]))
    }
  }
  Munst <<- round(M, 2)
  
  # ATD: if we have row summing up to zero (no transitions), then we
  #      should set each cell to the max entropy (1/siz) instead of
  #      min entropy (0) since if there were no observed transitions
  #      we can't very well predict where they're likely to go, hence
  #      we should have max entropy (max suprise)
  siz <- ncol(M)
  Srow <- sum(M[1,])
  if(abs(Srow) < 0.00001) {
    M[1,] <- 1.0/siz
  } else {
    for (j in 1:siz) {
      M[1,j] <- M[1,j]/Srow
      if(is.nan(M[1,j])) { M[1,j] = 0.0 }
    }
  }
  M <- round(M, 2)
  M <<- round(M, 2)
  
  ####### Print results
  if (print == TRUE) {
    print('***************************')
    print('Raw data Transition Matrix' )
    print('')
    print(Munst)
    print('***************************')
    print('Normalized Transition Matrix' )
    print('')
    print(M)
    print('Results written to M and Munst matrices. Use print(M) and print(Munst)')
  }
}


############################################################################
##################### TRANSITION MATRIX PLOT2 ##############################
## Funtion arguments:
#     transMatrix - transition matrix to be plotted 
#     plotName    - name of plot for saving (put in brackets)
#     plotColors  - color scheme for the plot
#     margin      - left, bottom margin
#     annColor    - color of annotations in cells (default 'black')
#     annAlpha    - transparency of cell colors (default '0.85')
## Function returns:
#     Saved pdf plot in the specified directory and name
## Usage example:
#     TransPlot2(transMatrix=M, plotName="Loveactually15_shot.pdf", 
#                     plotColors=brewer.pal(9,"Oranges"), annColor='black')
############################################################################
TransWPlot2 <- function(transMatrix,
                       plotName,
                       plotColors=brewer.pal(9,"Oranges"),
                       margin=c(4,4),
                       annColor='black',
                       annCex=1.1,
                       annAlpha=0.85,
                       title=NULL,
                       cexR=1.2,
                       cexC=1.2,
                       cexAxis=1.5,
                       xlabel='Word Span',
                       ylabel='') {
  pdf.options(family = "NimbusSan",useDingbats=FALSE)
  pdf(plotName,width=17,height=2)
  lwheatmap(transMatrix,
         ColorRamp=plotColors,
         cellColor=annColor,
         margin=margin,
         cex.cels=annCex,
         cex.rows=cexR,
         cex.cols=cexC,
         cex.lab=cexAxis,
         title=title,
         xlab=xlabel,
         ylab=ylabel
         )
  dev.off()
# embedFonts(plotName, "pdfwrite", outfile = plotName,
#            fontpaths =
#              c("/opt/local/share/texmf-texlive/fonts/type1/urw/helvetic",
#                "/usr/share/texmf/fonts/type1/urw/helvetic",
#                "/usr/local/teTeX/share/texmf-dist/fonts/type1/urw/helvetic",
#                "/usr/share/texmf-texlive/fonts/type1/urw/helvetic",
#                "/usr/local/texlive/texmf-local/fonts/type1/urw/helvetic"))
}
############################################################################
##################### TRANSITION MATRIX PLOT ###############################
## Funtion arguments:
#     transMatrix - transition matrix to be plotted 
#     plotName - name of plot for saving (put in brackets)
#     plotColots - color scheme for the plot
#     annColor - color of annotations in cells (default 'black')
## Function returns:
#     Saved pdf plot in the specified directory and name
## Usage example:
#     TransPlot(transMatrix=M, plotName="Loveactually15_shot.pdf", 
#                     plotColors=brewer.pal(9,"Oranges"), annColor='black')
############################################################################

TransPlot <- function (transMatrix,
                       plotName,
                       plotColors,
                       annColor='black',
                       annCex=2.5,
                       margin=c(15,15),
                       cexR=1.9,
                       cexC=1.9) {
  pdf.options(family = "NimbusSan",useDingbats=FALSE)
  pdf(plotName)
# ATD: use heat function in heat.R, which is a hacked version of heatmap.2
#      which accepts cex.lab to scale axis labels
# heatmap.2(x = transMatrix,
  heat(x = transMatrix,
            Rowv = FALSE,
            Colv = FALSE,
            dendrogram = 'none',
            scale = 'none',
            col = plotColors,
            cellnote = M,
            notecol = annColor,
            notecex = annCex,
            trace = 'none',
            margins = margin,
            cexRow = cexR,
            cexCol = cexC,
            key = FALSE,
            keysize = 1.0,
            density.info = "none",
            main = '',
            ylab = 'Source AOI (from)',
            xlab = 'Destination AOI (to)',
            cex.lab = 1.7,
            lmat = rbind(c(0,3),c(2,1),c(1,4)),
            lhei = c(1,8,1), # controls heights of rows
            lwid = c(1,8) # seems to control placement of main title
  )
  dev.off()
# embedFonts(plotName, "pdfwrite", outfile = plotName,
#            fontpaths =
#              c("/opt/local/share/texmf-texlive/fonts/type1/urw/helvetic",
#                "/usr/share/texmf/fonts/type1/urw/helvetic",
#                "/usr/local/teTeX/share/texmf-dist/fonts/type1/urw/helvetic",
#                "/usr/share/texmf-texlive/fonts/type1/urw/helvetic",
#                "/usr/local/texlive/texmf-local/fonts/type1/urw/helvetic"))

}


#############################################################################
#################### TRANSITION MATRIX ENTROPY ##############################
## Function arguments:
#     data - data file
#     AOInamesVar - variable which contain AOI names (put in brackets)
#     SubjectVar - variable which contain Subjects' names (put in brackets)
#     FixOrderVar - variable indicating time order of AOI hits within each
#                   participant (put in brackets)

## Returns: 
#     TMentrop - string with entropy value of transition matrix for each person
#                in the data file
## Usage example;
#     TransEntropy(data=df, AOInamesVar="AOI", SubjectsVar="Subject", FixOrderVar="")    
############################################################################

TransEntropy <- function(M, data, StimulusVar, SubjectsVar, SpanVar) {
  data <- data[order(data[ ,StimulusVar]), ]
  uniqS <- sort(unique(data[ ,StimulusVar]))
  lu <- length(uniqS)
  TMentrop <- numeric()
  for (i in 1:lu) {
    en <- numeric()
# empty matrix
    M[,] <- 0
    kukudf <- data[which(data[ ,StimulusVar] == uniqS[i]), ] # choose stimulus i
    luj <- dim(kukudf)[1] # how many stimuli for subject i
    if (luj > 1) {
      for (j in 1:luj) {
        span <- kukudf[j, SpanVar]
        # increase occurence of detected span
        M[1,as.character(span)] <- M[1,as.character(span)] + 1
      }
    } else {
        print('Warnings: ')
        print(sprintf("%s%s", 'No word spans for stimulus: ', uniqS[i]))
    }
    Munst <<- round(M, 2)

    #### Normalize to source - probability of going from i's AOI to any j's AOI
 #  for (i in 1:siz) {
 #    Srow <- sum(M[i,])
 #    for (j in 1:siz) {
 #      M[i,j] <- M[i,j]/Srow
 #      if(is.nan(M[i,j])) { M[i,j] = 0.0 }
 #    }
 #  }
    # ATD: if we have row summing up to zero (no transitions), then we
    #      should set each cell to the max entropy (1/siz) instead of
    #      min entropy (0) since if there were no observed transitions
    #      we can't very well predict where they're likely to go, hence
    #      we should have max entropy (max suprise)
    siz <- ncol(M)
    Srow <- sum(M[1,])
    if(abs(Srow) < 0.00001) {
      M[1,] <- 1.0/siz
    } else {
      for (j in 1:siz) {
        M[1,j] <- M[1,j]/Srow
        if(is.nan(M[1,j])) { M[1,j] = 0.0 }
      }
    }
    M <- round(M, 2)
    M <<- round(M, 2)

#   ATD: max entropy will be with each row's cells = 1/col, i.e.,
#   each AOI equally likely
#   max_entrop <- numeric()
#   max entropy is log_2(n) where n is the number of possible outcomes,
#   hence number of AOIs
#   max_entrop <- log2(siz)
#   default (empirical) entropy estimation is done via ML, maximum-likelihood
#   method; if there are many zero counts and the sample size is small, this
#   is very inefficient and also strongly biased
#   en <- entropy(M)
#   the bias-corrected maximum likelihood method, applying
#   the Miller-Madow correction to the empirical entropy
    # ATD: normalize entropy
#   en <- entropy(M,method="MM")
#   en <- entropy(M,method="MM")/max_entrop
    # compute normalized entropy (\eta(X)) in bits
#   en <- entropy(M,method="MM",unit="log2")/max_entrop
#   en <- entropy(M,method="ML",unit="log2")/max_entrop
    # don't use entropy function since it assumes probability distribution
    # M does not represent probability distribution, rather each of its
    # rows does
    en <- H_t(M)
#   could get NA here...
    if(is.na(en)) {
      en = 0.0
    }
    TMentrop <- c(TMentrop, en)
  }
# TMentrop <<- TMentrop
  TMentrop <<- data.frame(uniqS, TMentrop)
  names(TMentrop) <<- c("Subject", "Entropy")
}

############################################################################
################# TRANSITION ENTROPY CALCULATION ###########################
############################################################################
H_t <- function(M) {
  # ATD: max entropy will be with each row's cells = 1/col, i.e.,
  #      each AOI equally likely
  max_entrop <- numeric()
  # max entropy is log_2(n) where n is the number of possible outcomes,
  # hence number of AOIs
  max_entrop <- log2(ncol(M))

  # get stationary distribution
# p <- H_s(M)$sdist

  # see: http://www-isl.stanford.edu/~cover/papers/paper101.pdf
  # we really should do H(M) = -\sum_i \pi_i \sum_j M_{ij} log_2(M_{ij})
  H <- 0
  for (i in 1:nrow(M)) {
    c <- 0
#   for (j in 1:ncol(M)) {
#     if(abs(M[i,j]) > 0.0001) {
#       c = c + (M[i,j] * log2(M[i,j]))
#     }
#   }
    c = entropy(M[i,],method="ML",unit="log2")
    if(is.na(c)) {
      c = 0.0
    }
#   H = H + (p[i] * c)
    H = H + c
  }

  H <- H / max_entrop

# print(sprintf("H_t = %f",H))
  return(H)
}

############################################################################
################# STATIONARY ENTROPY CALCULATION ###########################
############################################################################
H_s <- function(M) {
  # matrix whose columns contain the eigenvectors
  e <- eigen(M)
  # or is it transpose of M?
  # see: http://faculty.cas.usf.edu/jkwilde/mathcamp/Linear_Algebra_II.pdf
# e <- eigen(t(M))
  # stationary distribution \pi is unchanged by the operation of transition
  # matrix \mathbf{P} on it, and so is defined by
  # \( \pi \mathbf{P} = \pi \)
  # by comparing this definition with that of an eigenvector, the two concepts
  # are related and that
  # \( \pi = \frac{e}{\sum_{i}{e_i}} \)
  # is a normalized ($\sum_{i}{\pi_i}$) multiple of a left $\mathbf{e}$ of the
  # transition matrix $\mathbf{P}$
  #
  lam <- e$values
# print(lam)
  vec <- e$vectors
# print(vec)
  v <- vec[,1]
# v <- Re(v)
  v <- Mod(v)
# print(sprintf("v, the stationary distribution eigenvector:\n"))
# print(v)
  # this gives us \pi, the stationary distribution from the eigenvector
  p <- v * 1.0/sum(v)
# print(sprintf("pi, the stationary distribution from the eigenvecotr:\n"))
# print(p)
  # $\sum_{i}{\pi_i}$ should equal 1.0
# print(sum(p))
  # see: http://www-isl.stanford.edu/~cover/papers/paper101.pdf
  # we really should do H(M) = -\sum_i \pi_i \sum_j M_{ij} log_2(M_{ij})
  H <- 0
#  for (i in 1:nrow(M)) {
#    c <- 0
#    for (j in 1:ncol(M)) {
#      if(abs(M[i,j]) > 0.0001) {
#        c = c + (M[i,j] * log2(M[i,j]))
#      }
#    }
##   print(c)
#    H = H + (p[i] * c)
#  }
  # as per the paper: H_s(M) = -\sum_i \pi_i log_2(pi_{i})
  for (i in 1:length(p)) {
    if(abs(p[i]) > 0.0001) {
      H = H + p[i] * log2(p[i])
    }
  }
  H <- -H
  # this should now be in bits per transition
# print(sprintf("H = %f\n",H))
  max_entrop <- numeric()
  # max entropy is log_2(n) where n is the number of possible outcomes,
  # hence number of AOIs
  n <- length(p)
# print(sprintf("nrow = %d, ncol = %d, n = %d\n",nrow(M),ncol(M),n))
  max_entrop <- log2(n)
# print(sprintf("max entropy = %f\n",max_entrop))
  # see also: http://aix1.uottawa.ca/~jkhoury/markov.htm
  eta <- H / max_entrop
  if(is.complex(eta)) {
    eta = Re(eta)
  }
  # eta should now be normalized stationary entropy since log2(n) is
  # max entropy if n is the number of possible outcomes
# print(sprintf("eta = %f\n",eta))
  retList <- list("sdist" = p, "eta" = eta)
  return(retList)
}

#############################################################################
#################### STATIONARY TRANSITION MATRIX ENTROPY ###################
## Function arguments:
#     data - data file
#     AOInamesVar - variable which contain AOI names (put in brackets)
#     SubjectVar - variable which contain Subjects' names (put in brackets)
#     FixOrderVar - variable indicating time order of AOI hits within each participant (put in brackets)

## Returns: 
#     TMentrop - string with entropy value of transition matrix for each person
#                in the data file
## Usage example;
#     TransEntropy(data=df, AOInamesVar="AOI", SubjectsVar="Subject", FixOrderVar="")    
############################################################################

StationaryEntropy <- function(M, data, SubjectsVar, AOInamesVar, FixOrderVar) {
  data <- data[order(data[ ,SubjectsVar], data[ ,FixOrderVar]), ]
  uniqS <- sort(unique(data[ ,SubjectsVar]))
  lu <- length(uniqS)
  uniqAOI <- sort(unique(data[ ,AOInamesVar])) # unique values of AOI variable
# siz <- length(unique(data[ ,AOInamesVar])) # how many AOIs
  siz <- nrow(M)
  if(nrow(M) != ncol(M)) {
    print('Warnings: ')
    print(sprintf("%s%dx%d", 'Matrix not square: ',nrow(M),ncol(M)))
  }
  STentrop <- numeric()
  for (i in 1:lu) {
    en <- numeric()
# empty matrix
#   M <- matrix(data=0, nrow=siz, ncol=siz, dimnames=list(uniqAOI, uniqAOI))
    M[,] <- 0
    kukudf <- data[which(data[ ,SubjectsVar] == uniqS[i]), ] # choose subject i
    luj <- dim(kukudf)[1] # how many AOIs for subject i
    if (luj > 1) {
      j <- 0
      repeat {
        j <- j+1
        from <- kukudf[j, AOInamesVar]
        to <- kukudf[j+1, AOInamesVar]
        M[as.character(from), as.character(to)] <-
        M[as.character(from), as.character(to)] + 1 
        if (j > luj-2) break()
      } 
    } else {
        print('Warnings: ')
        print(sprintf("%s%s", 'No AOI transitions for subject: ', uniqS[i]))
    }
    Munst <<- round(M, 2)

    #### Normalize to source - probability of going from i's AOI to any j's AOI
#   for (i in 1:siz) {
#     Srow <- sum(M[i,])
#     for (j in 1:siz) {
#       M[i,j] <- M[i,j]/Srow
#       if(is.nan(M[i,j])) { M[i,j] = 0.0 }
#     }
#   }
    # ATD: if we have row summing up to zero (no transitions), then we
    #      should set each cell to the max entropy (1/siz) instead of
    #      min entropy (0) since if there were no observed transitions
    #      we can't very well predict where they're likely to go, hence
    #      we should have max entropy (max suprise)
    for (i in 1:siz) {
      Srow <- sum(M[i,])
      if(abs(Srow) < 0.00001) {
        M[i,] <- 1.0/siz
      } else {
        for (j in 1:siz) {
          M[i,j] <- M[i,j]/Srow
          if(is.nan(M[i,j])) { M[i,j] = 0.0 }
        }
      }
    }
    M <- round(M, 2)
    M <<- round(M, 2)

    sen <- H_s(M)
    eta <- sen$eta
#   could get NA here...
    if(is.na(eta)) {
      eta = 0.0
    }
    STentrop <- c(STentrop, eta)
  }
# STentrop <<- STentrop
  STentrop <<- data.frame(uniqS, STentrop)
  names(STentrop) <<- c("Subject", "SEntropy")
}

#############################################################################
################## LEVENSTEIN INDEX OF SCANPATH SIMILARITIES ################
## Function arguments:
#     data - data file
#     SubjectVar - variable which contain Subjects' names (put in brackets)
#     AOInamesVar - variable which contain AOI names (put in brackets)
#     FixOrderVar - variable indicating time order of AOI hits within each participant (put in brackets)
## Returns: 
#     LM - matrix containing similarities measure (Levenstein's H) between participants' scanpaths
#     
## Usage example;
#     LevenMatrix(data=dd, SubjectsVar="Subject", AOInamesVar="AOI", FixOrderVar="")    
############################################################################

LevenMatrix <- function(data, SubjectsVar, AOInamesVar, FixOrderVar){
  data <- data[order(data[ ,SubjectsVar], data[ ,FixOrderVar]), ]
  uniqS <- sort(unique(data[ ,SubjectsVar])) # unique values of Subject variable
  lui <- length(unique(data[ ,SubjectsVar])) # how many subjects 
  uniqAOI <- sort(unique(data[ ,AOInamesVar])) # unique values of AOI variable
  M <- matrix(data=NA, nrow=lui, ncol=lui, dimnames=list(uniqS, uniqS)) # empty matrix
  for (i in 1:lui){
    s <- as.character(uniqS)[i]
    str1 <- toString(as.character(data[which(data[,SubjectsVar] == s), AOInamesVar]))
    j <- 0
    for (j in 1:lui) {
      ss <- as.character(uniqS)[j]
      str2 <- toString(as.character(data[which(data[,SubjectsVar] == ss), AOInamesVar]))
      M[i,j] <- stringMatch(str1, str2)
    }
  }
  LM <<- M
}

