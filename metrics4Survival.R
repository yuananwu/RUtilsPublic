library('survival')
library('survminer')
library('rms')
library('pec') ##验证模型
library('timeROC')
library('CsChange')
library('dplyr')
library('riskRegression')

# 这个工具类中包含了生存分析建模过程中一些指标(c-index, roc等)的计算函数
#===============================================km法相关

# https://zyicu.cn/?p=13869
library('survivalROC')
getSurvivalCutoff = function(df, tm, status, variable, method = 'km', rocPredTm = 0){
  if (method == 'km'){
    cutoff = surv_cutpoint(df, time = tm, event = status, variables = c(variable))
    print(cutoff)
    return(cutoff[["cutpoint"]]['cutpoint'])
  }else if (method == 'roc'){
    if(rocPredTm == 0){
      stop("method 'roc' need to specify the predict time")
    }
    #确定指标的最佳cutoff值，构建生存ROC函数（survivalROC）
    rocM = survivalROC(Stime = df[, tm],    #生存时间=time
                       status = df[, status], #生存状态=event
                       marker = df[, variable], #需要分析的变量
                       predict.time = mean(df[, tm]),      #预测时间
                       method="KM")          #使用生存分析KM计算
    print(paste('auc: ', rocM$AUC, sep = ''))
    #约登指数最大对应的值为最佳cutoff
    cutoff = rocM$cut.values[which.max(rocM$TP-rocM$FP)]
    return(cutoff)
  }else{
    return("")
  }
}

simplePlotKm = function(df, label, fea = 1){
  fit <- survfit(as.formula(paste(label, fea, sep = '~')),data = df) 
  ggplot = ggsurvplot(fit, data = df, pval = T, risk.table = T)
  return(ggplot)
}

#===============================================参数模型(cox)相关
#=============1 c-index计算
# 获取c-index及其95CI
getCIdxInfo = function(model, df, rstType = 'cIdx'){
  cIdxList = survival::concordance(model,newdata=df)
  if (rstType == 'cIdx'){
    return(cIdxList$concordance)
  }else if (rstType == 'var'){
    return(cIdxList$var)
  }else{return(NULL)}
}
getCindex95CI = function(model, df, roundNum = 4){
  if(!inherits(model, 'coxph')){
    stop('model 需要是来自于survival方法的coxph构建')
  }
  cIdx = getCIdxInfo(model, df)
  se = sqrt(getCIdxInfo(model, df, 'var'))
  res = sprintf('%s (%s-%s)'
                ,str_pad(round(cIdx, roundNum), roundNum + 2, side = 'right', pad = '0')
                ,str_pad(round(cIdx - 1.96*se, roundNum), roundNum + 2, side = 'right', pad = '0')
                ,str_pad(round(cIdx + 1.96*se, roundNum), roundNum + 2, side = 'right', pad = '0'))
  return(res)
}
# 对来自同一数据集的两个评分计算cIdx的p值
# p值小于0.05代表两个有差异
getCIdxPValueFromPairedData = function(df, fo1, fo2, numBoot = 100){
  md1 = coxph(as.formula(fo1), data = df)
  md2 = coxph(as.formula(fo2), data = df)
  rst = CsChange(md1, md2, data=df, nb = numBoot)
  p = as.numeric(rst[[1]]['p'])
  return(p)
}

#=============2 time-ROC计算

tmRocTemplate = function(tm,event,pred, tmList){
  # https://cloud.tencent.com/developer/article/1675084
  timeROC(
    T = tm,
    delta = event,
    marker = pred,
    cause = 1,
    weighting="marginal",
    times = tmList,
    ROC = TRUE,
    iid = TRUE
  )
}
# 返回timeROC生成的roc带95ci的结果, 组成一个dataframe
# https://cloud.tencent.com/developer/article/1675084
tmRocRst = function(roc1){
  auc = as.data.frame(roc1$AUC)
  auc['tm'] = row.names(auc)
  
  # print(auc)
  auc95ci = as.data.frame(confint(roc1, level = 0.95)$CI_AUC)
  auc95ci['tm'] = row.names(auc95ci)
  
  rst = merge(auc, auc95ci, by = 'tm', all.x = TRUE)
  colnames(rst) = c('tm', 'auc', '95low', '95upper')
  rst[, 'auc']  = rst[, 'auc']
  rst[, '95low'] = rst[, '95low']/100
  rst[, '95upper'] = rst[, '95upper']/100
  removeStr = function(x, sep = '='){
    strsplit(x, sep)[[1]][2]
  }
  rst[, 'tm'] = apply(rst['tm'], 1, removeStr)
  return(rst)
}

# riskRegression中的score方法计算的timeROC
# !!timeROC和 riskRegression 中计算auc的有个细节差异点, timeROC中会有小于0.5的, 而riskRegression会将其变为1-x
tmRocFromRR = function(df, model, tmList){
  # md = coxph(as.formula(fo), data = df, x=TRUE, y=TRUE)
  metricsRst = riskRegression::Score(list(RiskScore = model), 
                                     # summary = c("risks","IPA","riskQuantile","ibs"),
                                     summary = NULL,
                                     formula=as.formula(paste(as.character(model$formula[2]), 1, sep = '~')), 
                                     data = df,
                                     null.model=F,
                                     times = tmList,
                                     metrics =c("auc"))
  rst = as.data.frame(metricsRst$AUC$score)
  rst[, 'tm']  = rst[, 'times']
  rst[, 'auc']  = round(rst[, 'AUC'], 4)
  rst[, '95low']  = round(rst[, 'lower'], 4)
  rst[, '95upper']  = round(rst[, 'upper'], 4)

  return(rst[, c('tm', 'auc', '95low', '95upper')])

}

displayROC = function(rows, roundNum = 3){
  tmpFunc = function(x, roundNum){
    x = as.numeric(x)
    str_pad(round(x, roundNum), roundNum + 2, side = 'right', pad = '0')
  }
  rows = sapply(rows, function(x) tmpFunc(x, roundNum = roundNum))
  rst = sprintf('%s (%s-%s)', rows[2], rows[3], rows[4])
  return(rst)
}

# 使用demo
# a = tmRocRst(tmRocTemplate(dfCal[,tm], dfCal[,event], dfCal[,mdName], c(30,36,42)))
# a[, 'display'] = apply(a, 1, displayROC)
# 
# b = tmRocFromRR(dfCal, paste(label, mdName, sep = '~'), c(30,36,42))
# b[, 'display'] = apply(b, 1, displayROC)

# timeROC比较
compareTimeROC = function(roc1, roc2){
  cmp = compare(roc1, roc2, adjusted = TRUE)
  return(cmp$p_values_AUC)
}

#绘制不同时间点的time AUC折线段
plotAUCCurve = function (object, add = FALSE, col = "black", xbreaks = 12, minY = 0.4, addErrorBar = FALSE,  xlab = '', ylab = '') {
  plotError <- function(x, y, sd, len = 1, col = "black", lwd = 1) {
    len <- len * 0.05
    arrows(x0 = x, y0 = y, x1 = x, y1 = y - sd, col = col, angle = 90, length = len, lwd = lwd, add = T)
    arrows(x0 = x, y0 = y, x1 = x, y1 = y + sd, col = col, angle = 90, length = len, lwd = lwd, add = T)
  }
  if (class(object) == "ipcwsurvivalROC") {
    AUC = object$AUC[!is.na(object$AUC)]
    times = object$times[!is.na(object$AUC)]
    se = object$inference$vect_sd_1[!is.na(object$AUC)]
    lwd = 2
    if (add == FALSE) {
      plot(times, AUC, xlab = xlab, ylab = ylab, lwd = lwd, ylim = c(minY, 1), type = "l", col = col, xaxt = "n", yaxt = "n")
      # 横轴
      axis(1, seq(min(times), max(times), xbreaks), seq(min(times), max(times), xbreaks))
      # 纵轴
      axis(2,seq(minY, 1, 0.1), seq(minY, 1, 0.1))
      abline(h = 0.5, lty = 2)
    }else {
      lines(times, AUC, lwd = lwd, type = "l", col = col)
    }
    if (addErrorBar){
      for (i in 1: length(AUC)) {
        plotError(times[i], AUC[i], 1.96*se[i], col = col, lwd = lwd)
      }
    }
  }
}



#=============3 brier计算

# 计算IBS值
# https://stats.stackexchange.com/questions/574772/difference-between-predict-and-predictsurvprob-r-functions
# 这个有助于理解预测值
getIBS = function(df, cphModel){
  survObj = Surv(df[, tm], df[, event])
  timeList = sort(unique(df[, tm]))
  predProb = predictSurvProb(cphModel, newdata = df, times = timeList)
  rst = IBS(survObj, predProb) # sp_matrix行是样本, 列是不同时间点
  # ps: IBSrange <- seq(t_IBSrange[1], t_IBSrange[2], length = p) 函数中时间点的范围确定是根据最大最小值和步长来设定的, 和常规理解会有误区, 但结果和按照timeList划分的没有差
}

getBrierScore = function(df, coxModel, tm, event, tmList){
  survObj = Surv(df[, tm], df[, event])
  rst = NULL
  for (tm in tmList){
    tmPredProb = predictSurvProb(coxModel, newdata = df, times = tm)
    bScore = Brier(survObj, tmPredProb)
    rst = rbind(rst, c(tm, bScore))
  }
  return(rst)
}
# riskRegression中的score方法计算的timeROC
brierTemplateFromRR = function(df, model, tmList){
  # md = coxph(as.formula(fo), data = df, x=TRUE, y=TRUE)
  metricsRst = riskRegression::Score(list(RiskScore =model), 
                                     # summary = c("risks","IPA","riskQuantile","ibs"),
                                     summary = NULL,
                                     formula=as.formula(paste(as.character(model$formula[2]), 1, sep = '~')), 
                                     data = df,
                                     null.model=F,
                                     times = tmList,
                                     metrics =c("brier"))
  return(metricsRst)
}
brierRstFromRR = function(df, model, tmList){
  metricsRst = brierTemplateFromRR(df, model, tmList)
  rst = as.data.frame(metricsRst$Brier$score)
  rst[, 'tm']  = rst[, 'times']
  rst[, 'Brier']  = round(rst[, 'Brier'], 4)
  rst[, '95low']  = round(rst[, 'lower'], 4)
  rst[, '95upper']  = round(rst[, 'upper'], 4)

  return(rst[, c('tm', 'Brier', '95low', '95upper')])
  
}
# 比较brier得分
compareBrierFromRR = function(df, label, f1, f2, tmList, randomSeedNum = 42, bootNum = 100){
  {
    set.seed(randomSeedNum)
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)){
      runif(1)
    }
    assign(".Random.seed", get(".Random.seed", envir = .GlobalEnv, inherits = FALSE), envir = .GlobalEnv)
  }
  # 同一数据集不同方法的brier得分的差
  twoBrierInterval = function(df, indices, label, f1, f2, tmList, randomSeedNum){
    set.seed(randomSeedNum)
    dataBoot = df[indices, ]
    md1 = coxph(as.formula(paste(label, f1, sep = '~')), data = dataBoot, x=TRUE, y=TRUE)
    md2 = coxph(as.formula(paste(label, f2, sep = '~')), data = dataBoot, x=TRUE, y=TRUE)
    b1 = brierTemplateFromRR(dataBoot, md1, tmList)$Brier$score[, c('times', 'Brier')]
    b2 = brierTemplateFromRR(dataBoot, md2, tmList)$Brier$score[, c('times', 'Brier')]
    m = data.frame(merge(b1, b2, by = 'times'))
    rst = as.vector(m[,2] -m[,3])
    c(rst)
  }
  bootRst = boot(df, twoBrierInterval, label = label, f1 = f1, f2 = f2, tmList = tmList, randomSeedNum = 42, R = bootNum)
  rst = NULL
  for (i in 1:length(tmList)){
    tm = tmList[i]
    interval = bootRst$t0[i]
    bootDetail = bootRst$t[,i]
    se = sqrt(var(bootDetail))
    z = interval/se
    p = 2 * pnorm(abs(z), lower.tail = F)
    rst = rbind(rst, c(tm, p))
  }
  rst = as.data.frame(rst)
  colnames(rst) = c('time', 'p-value')
  return(rst)
}
#绘制brier得分曲线
plotBrierCurve = function (brierData, add = FALSE, col = "black", xbreaks = 12, ybreaks = 0.05, ySeq = c(0, 0.4), addErrorBar = FALSE,  xlab = '', ylab = '') {
  plotError <- function(x, y, yUp, yLow, len = 1, col = "black", lwd = 1) {
    len <- len * 0.05
    arrows(x0 = x, y0 = y, x1 = x, y1 = yLow, col = col, angle = 90, length = len, lwd = lwd)
    arrows(x0 = x, y0 = y, x1 = x, y1 = yUp, col = col, angle = 90, length = len, lwd = lwd)
  }
  times = brierData[, 'tm']
  brier = brierData[, 'Brier']
  brierUp  = brierData[, '95upper']
  brierLow = brierData[, '95low']
  lwd = 2
  if (add == FALSE) {
    plot(times, brier, xlab = xlab, ylab = ylab, lwd = lwd, ylim = ySeq, type = "l", col = col, xaxt = "n", yaxt = "n")
    # 横轴
    axis(1, seq(min(times), max(times), xbreaks), seq(min(times), max(times), xbreaks))
    # 纵轴
    axis(2, seq(ySeq[1], ySeq[2], ybreaks), seq(ySeq[1], ySeq[2], ybreaks))
  }else {
    lines(times, brier, lwd = lwd, type = "l", col = col)
  }
  if (addErrorBar){
    for (i in 1: length(brier)) {
      plotError(x = times[i], y = brier[i], yUp = brierUp[i], yLow = brierLow[i], col = col, lwd = lwd)
    }
  }
}



# ======================================================================================

#整合cox模型的输出数据
getCoxphModelResult = function(model, cal10PointFlag){
  if (class(model) != 'coxph'){
    stop("cox model's source doesn't match 'coxph'")
  }
  rst = summary(model)
  coefList = rst$coefficients
  hr95CI = rst$conf.int
  
  rstDf = as.data.frame(row.names(coefList))
  colnames(rstDf) = 'feature'
  rstDf[, 'coef'] = as.numeric(coefList[, 'coef'])
  rstDf[, 'se']   = as.numeric(coefList[, 'se(coef)'])
  rstDf[, 'coef_95ci_lower'] = rstDf[, 'coef'] - 1.96*rstDf[, 'se']
  rstDf[, 'coef_95ci_upper'] = rstDf[, 'coef'] + 1.96*rstDf[, 'se']
  rstDf[, 'HR'] = hr95CI[, 'exp(coef)']
  rstDf[, 'HR_95ci_lower'] = hr95CI[, 'lower .95']
  rstDf[, 'HR_95ci_upper'] = hr95CI[, 'upper .95']
  rstDf[, 'z_value'] = coefList[, 'z']
  rstDf[, 'p_value'] = round(coefList[, 'Pr(>|z|)'], 6)
  
  if(cal10PointFlag){
    rstDf['coef_10_point'] = round(rstDf[,'coef']/max(abs(rstDf[,'coef']))*10, 0)
  }
  return(rstDf)
}

getCoxRstListV3 = function(inputCols, coxModel, tm, event, rocTmList, extInfo = NULL, dfList, dfDescList){
  label = paste('Surv(', tm, ',', event, ')',sep = '')
  # 输出模型结果
  coxModelRstList = NULL
  if(!is.null(extInfo)){
    coxModelRstList[['拓展信息']] = extInfo
  }
  
  coxModelRstList[['模型入组特征']] = inputCols
  coxModelRstList[['最终模型的特征']] = paste(names(coxModel$assign), collapse = ', ')
    
  coxModelRstList[['模型系数']] = getCoxphModelResult(coxModel, cal10PointFlag = T)
  coxModelRstList[['模型AIC']] = extractAIC(coxModel)[2]
  if (length(dfList) != length(dfDescList)){
    stop('input df num doesnt equal df describe')
  }
  for(idx in 1:length(dfList)){
    predName = 'coxModelPred'
    df = dfList[[idx]]
    dfDesc = dfDescList[idx]
    df[,predName] = predict(coxModel, newdata = df)
    coxModelRstList[[paste(dfDesc, 'Cindex', sep = '')]] = getCindex95CI(coxModel, df)  
  }
  for(idx in 1:length(dfList)){
    predName = 'coxModelPred'
    df = dfList[[idx]]
    dfDesc = dfDescList[idx]
    df[,predName] = predict(coxModel, newdata = df)
    
    rocRst = tmRocRst(tmRocTemplate(df[,tm], df[,event],df[,predName], rocTmList))
    rocRst[, 'display'] = apply(rocRst, 1, displayROC)
    coxModelRstList[[paste(dfDesc, 'timeRoc', sep = '')]] = rocRst[, c('tm', 'display')]
  }
  for(idx in 1:length(dfList)){
    predName = 'coxModelPred'
    df = dfList[[idx]]
    dfDesc = dfDescList[idx]
    df[,predName] = predict(coxModel, newdata = df)
    
    bScore = brierRstFromRR(df, coxModel, rocTmList)
    bScore[, 'display'] = apply(bScore, 1, displayROC)
    # print(bScore[, c('tm', 'display')], row.names = F)
    coxModelRstList[[paste(dfDesc, 'BrierScore', sep = '')]] = bScore[, c('tm', 'display')]
  }
  
  return(coxModelRstList)
}

