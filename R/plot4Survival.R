library('survival')
library('survminer')
library('rms')
# library('survcomp')
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
getCindex95CI = function(model, df){
  if(!inherits(model, 'coxph')){
    stop('model 需要是来自于survival方法的coxph构建')
  }
  getCIdx = function(model, df, rstType = 'cIdx'){
    cIdxList = survival::concordance(model,newdata=df)
    if (rstType == 'cIdx'){
      return(cIdxList$concordance)
    }else if (rstType == 'var'){
      return(cIdxList$var)
    }else{return(NULL)}
  }
  cIdx = getCIdx(model, df)
  se = sqrt(getCIdx(model, df, 'var'))
  res = paste(round(cIdx,4),' (',round(cIdx - 1.96*se,4),'-',round(cIdx + 1.96*se,4),')',sep = '')
  return(res)
}
#=============2 time-ROC计算


tmRocTemplate = function(df,tm,event,pred, tmList){
  # https://cloud.tencent.com/developer/article/1675084
  timeROC(
    T = df[,tm],
    delta = df[,event],
    marker =df[,pred],
    cause = 1,
    weighting="marginal",
    times = tmList,
    ROC = TRUE,
    iid = TRUE
  )
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


tmROCResult = function(roc1){
  # https://cloud.tencent.com/developer/article/1675084
  auc = as.data.frame(roc1$AUC)
  auc['tm'] = row.names(auc)
  
  # print(auc)
  auc95ci = as.data.frame(confint(roc1, level = 0.95)$CI_AUC)
  auc95ci['tm'] = row.names(auc95ci)
  
  rst = merge(auc, auc95ci, by = 'tm', all.x = TRUE)
  colnames(rst) = c('AUC_time', 'AUC', 'AUC_95ci_lower', 'AUC_95ci_upper')
  rst['AUC_4_display'] = paste(round(rst[, 'AUC'], 4), ' (', round(rst[, 'AUC_95ci_lower']/100,4), '-', round(rst[, 'AUC_95ci_upper']/100,4), ')', sep = '')
  
  rst2 = rst[,c('AUC_time', 'AUC_4_display')]
  colnames(rst2) = c('AUC_time', 'AUC')
  return(rst2)
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


