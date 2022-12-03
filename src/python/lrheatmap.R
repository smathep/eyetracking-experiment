# from: http://www.phaget4.org/R/image_matrix.html
#
# ----- Define a function for plotting a matrix ----- #
lrheatmap <- function(x,
                      cex.lab=1.0,	# size adjustment for axis labels
                      cex.rows=1.0,	# size adjustment for row labels
                      cex.cols=1.0,	# size adjustment for col labels
                      cex.cels=1.0,	# size adjustment for cell numerals
                      cex.main=1.0,
                      cellColor='black',
                      cellAlpha=0.85,
                      margin=c(8,8),	# left, bottom margins
                      ColorRamp=rgb( seq(0,1,length=256),	# Red
                                     seq(0,1,length=256),	# Green
                                     seq(1,0,length=256)),	# Blue
                      transpose=TRUE,
                      ...) {

  min <- min(x)
  max <- max(x)
  yLabels <- rownames(x)
  xLabels <- colnames(x)
  title <-c()

  # check for additional function arguments
  if( length(list(...)) ) {

    Lst <- list(...)

    # min, max (for color ramp)
    if( !is.null(Lst$zlim) ) {
      min <- Lst$zlim[1]
      max <- Lst$zlim[2]
    }

    # column labels
    if( !is.null(Lst$yLabels) )
      yLabels <- c(Lst$yLabels)

    # row labels
    if( !is.null(Lst$xLabels) )
      xLabels <- c(Lst$xLabels)

    # title
    if( !is.null(Lst$title) )
      title <- Lst$title

    # x axis label
    if( !is.null(Lst$xlab) )
      xlab <- Lst$xlab

    # y axis label
    if( !is.null(Lst$ylab) )
      ylab <- Lst$ylab
  }

  # check for null values
  if( is.null(xLabels) ) {
    xLabels <- c(1:ncol(x))
  }
  if( is.null(yLabels) ) {
    yLabels <- c(1:nrow(x))
  }

  # graphic layout -- useful if you want a colorramp legend
  # [ 1 2 ]
# layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))

  ColorRamp <- adjustcolor(ColorRamp,cellAlpha)

  # color levels
  ColorLevels <- seq(min, max, length=length(ColorRamp))

  # reverse Y axis
  reverse <- nrow(x) : 1
  yLabels <- yLabels[reverse]
  x <- x[reverse,]

  # graphics output

  # margins (bottom, left, top, right)
# margins=c(margin[1],margin[2],2,2)
  margins=c(margin[1],margin[2],4,2)
  par(mar=margins)

  # cell numbers (as text)
  cellnote <- matrix(x, ncol=ncol(x), nrow=nrow(x))
  if(transpose) {
    cellnote <- t(cellnote)
  } else {
    cellnote <- cellnote
  }
  cellnote <- format(cellnote,digits=2)

  # draw the image
  if(transpose) {
    image(1:length(xLabels), 1:length(yLabels), t(x),
          col=ColorRamp,
          asp=1.0,
          xlab="",
          ylab="",
          axes=FALSE,
          zlim=c(min,max))
  } else {
    image(1:length(xLabels), 1:length(yLabels), x,
          col=ColorRamp,
          asp=1.0,
          xlab="",
          ylab="",
          axes=FALSE,
          zlim=c(min,max))
  }

  if( !is.null(title) ) {
    title(main=title,cex.main=cex.main,line=1.8)
  }

  # row labels; las: 0: parallel, 1: horiz., 2: perp., 3: vertical
  axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, cex.axis=cex.rows,
       las = HORIZONTAL<-2, tick=0, line=-0.5)
  # col labels
# axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=cex.cols,
#      las = VERTICAL<-3, tick=0, line=-0.5)
  axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=cex.cols,
       las = PARALLEL<-1, tick=0, line=-0.5)
  if(!is.null(xlab)) mtext(xlab, cex=cex.lab, side=1, line=margins[1] - 1.5)
  if(!is.null(ylab)) mtext(ylab, cex=cex.lab, side=2, line=margins[2] - 1.5)

  # cellnotes (numbers within cells)
  text(x=c(row(cellnote)),
       y=c(col(cellnote)),
       labels=c(cellnote),
       col=cellColor,
       cex=cex.cels)

}
# ----- END plot function ----- #
