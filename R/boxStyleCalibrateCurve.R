
kmCalibrateTemplate = function(df, fo, timeInc){
  coxPlot = cph(fo, data=df, surv=T, x=T, y=T, time.inc = timeInc)
  bucketNum = as.integer(dim(df)[1]/3)
  calKM = calibrate(coxPlot, u=coxPlot$time.inc, cmethod='KM', m=bucketNum, B=100)
  return(calKM)
}

calibrateKmDataSort = function(objKm){
  ## Compute confidence limits for survival based on -log survival,
  ## constraining to be in [0,1]; d = std.error of cum hazard * z value
  getCi = function(rows, direction){ # row的构成：第一列是值， 第二列是标准误
    surv = rows[1]
    d    = rows[2] * 1.96
    if(direction == 'up'){
      rst = ifelse(surv==0, 0, pmin(1, surv*exp(d)))
    }else if(direction == 'low'){
      rst = ifelse(surv==0, 0, surv*exp(-d))
    }else{
      rst = NULL
    }
    return(rst)
  }
  
  if (class(objKm) != 'calibrate'){
    stop('the object\'s class must be \'calibrate\'')
  }
  pred = objKm[,"mean.predicted"]
  cal  = objKm[,"KM"]
  se = objKm[,"std.err"]
  calUpper = apply(cbind(cal , se), 1, getCi, 'up')
  calLower = apply(cbind(cal , se), 1, getCi, 'low')
  rst = data.frame(pred, cal, 
                   # se, 
                   calUpper, calLower)
  return(rst)
}
getXYLimit = function(plotData){
  maxXY = ceiling(max(plotData)*10)/10
  minXY = floor(min(plotData)*10)/10
  return(c(pmax(0, minXY), pmin(1, maxXY)))
}
plotCalibrateKm = function(plotData, xlab = NULL, ylab = NULL, col = 'black', lwd = 2, addErrorBar = T, add = F, xyLimit = c(0, 1), xyBreak = 0.2){
  plotError <- function(xStart, yStart, yUpper, yLower, len = 1, col, lwd) {
    len <- len * 0.05
    arrows(x0 = xStart, y0 = yStart, x1 = xStart, y1 = yUpper, col = col, angle = 90, length = len, lwd = lwd)
    arrows(x0 = xStart, y0 = yStart, x1 = xStart, y1 = yLower, col = col, angle = 90, length = len, lwd = lwd)
  }
  x = plotData[, 'pred']
  y = plotData[, 'cal']
  if (add == FALSE) {
    plot(x, y, xlab = xlab, ylab = ylab, lwd = lwd, ylim = xyLimit, xlim = xyLimit, type = "l", col = col, xaxt = "n", yaxt = "n")
    points(x, y, pch = 4, col = col, lwd = lwd)
    axisSeq = seq(xyLimit[1], xyLimit[2], xyBreak)
    # 横轴
    axis(1, axisSeq, axisSeq)
    # 纵轴
    axis(2, axisSeq, axisSeq)
    abline(0,1,lty=3,lwd=2,col = '#E2E3E5')
  }else {
    lines(x, y, lwd = lwd, type = "l", col = col)
    points(x, y, pch = 4, col = col, lwd = lwd)
  }
  if (addErrorBar){
    for (i in 1: length(x)) {
      plotError(xStart = x[i], yStart = y[i], yUpper = plotData[, 'calUpper'][i], yLower = plotData[, 'calLower'][i], col = col, lwd = lwd)
    }
  }
}